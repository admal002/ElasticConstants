#include "KIM_API_status.h"
#define THIS_FILE_NAME __FILE__
#define TRUEFALSE(TRUTH) merge(1,0,(TRUTH))

module mod_kim
    use mod_global
    use KIM_API_F03
    implicit none

    integer(c_int), parameter :: SizeOne        = 1
    !
    ! KIM variables
    !
    character(len=KIM_KEY_STRING_LENGTH) :: NBC_Method
    ! 0- NEIGH_RVEC_H, 1- NEIGH_PURE_H, 2- NEIGH_RVEC_F, 3- NEIGH_PURE_F,
    ! 4- MI_OPBC_H,    5- MI_OPBC_F
    integer(c_int) :: ier, idum
    contains

        function getCutoff_kim(pkim)
            implicit none
            ! Passed variables
            type(c_ptr)    :: pkim
            real(c_double) :: getCutoff_kim

            ! Local variables
            real(c_double), pointer :: cutoff;type(c_ptr) :: pcutoff

            call kim_api_getm_data(pkim, ier,"cutoff",pcutoff,1)
            if (ier.lt.KIM_STATUS_OK) then
                idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                         "kim_api_getm_data", ier)
                stop
            endif
            call c_f_pointer(pcutoff, cutoff)
            getCutoff_kim = cutoff
            return
        end function getCutoff_kim

        subroutine setup_kim(pkim,nbc)
            implicit none
            ! Passed variables
            type(c_ptr)                                     :: pkim
            integer(c_int),intent(out)                      :: nbc

            ! Local variables
            character(len=KIM_KEY_STRING_LENGTH) :: NBC_Method

            ! determine which NBC scenerio to use
            ier = kim_api_get_nbc_method(pkim, NBC_Method)
            if (ier.lt.KIM_STATUS_OK) then
                idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                         "kim_api_get_nbc_method", ier)
                stop
            endif
            if (index(NBC_Method,"NEIGH_RVEC_H").eq.1) then
                nbc = 0
            elseif (index(NBC_Method,"NEIGH_PURE_H").eq.1) then
                nbc = 1
            elseif (index(NBC_Method,"NEIGH_RVEC_F").eq.1) then
                nbc = 2
            elseif (index(NBC_Method,"NEIGH_PURE_F").eq.1) then
                nbc = 3
            elseif (index(NBC_Method,"MI_OPBC_H").eq.1) then
                nbc = 4
            elseif (index(NBC_Method,"MI_OPBC_F").eq.1) then
                nbc = 5
            else
                ier = KIM_STATUS_FAIL
                idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                         "Unknown NBC method", ier)
                stop
            endif

        end subroutine setup_kim

        subroutine allocate_kim(pkim,numParticles,numSpecies)
            implicit none
            ! Passed variables
            type(c_ptr)     :: pkim
            integer(c_int)  :: numParticles,numSpecies

            ! Allocate memory via the KIM system
            call kim_api_allocate(pkim, numParticles, numSpecies, ier)
            if (ier.lt.KIM_STATUS_OK) then
                idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                         "kim_api_allocate", ier)
                stop
            endif

        end subroutine allocate_kim

        subroutine initialize_kim(pkim,cutoff)
            use mod_atomistic, only : Atomistic
            implicit none
            ! Passed variables
            type(c_ptr)             :: pkim
            real(c_double),target   :: cutoff

            ier = kim_api_set_data(pkim,"cutoff",SizeOne, c_loc(cutoff))
            if (ier.lt.KIM_STATUS_OK) then
             idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                         "kim_api_setm_data", ier)
             stop
            endif

            ! call model's init routine
            ier = kim_api_model_init(pkim)
            if (ier.lt.KIM_STATUS_OK) then
                idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                         "kim_api_model_init", ier)
                stop
            endif

        end subroutine initialize_kim

        subroutine broadcastAtomistic_kim(pkim,aS)
            use mod_atomistic, only : Atomistic,NeighObject_type,   &
                                      process_dedr,process_d2edr2
                                      

            ! Passed variables
            type(c_ptr), intent(in) :: pkim
            type(Atomistic),pointer :: aS
            
            ! Local variables
            integer(c_int) :: pbcFlag, halfFlag,sizecoor


            pbcFlag = TRUEFALSE((aS%nbc.eq.4).or.(aS%nbc.eq.5))
            halfFlag = TRUEFALSE((aS%nbc.eq.0).or.(aS%nbc.eq.1).or.(aS%nbc.eq.4))

            sizecoor = DIM*aS%numParticles
            call kim_api_setm_data(pkim, ier, &
            "numberOfParticles",           SizeOne, c_loc(aS%numParticles),   1,&
            "numberContributingParticles", SizeOne, c_loc(aS%numContributingParticles), halfFlag,&
            "numberOfSpecies",             SizeOne, c_loc(aS%numSpecies),     1,&
            "particleSpecies",     aS%numParticles, c_loc(aS%particleSpecies),1)
            if (ier.lt.KIM_STATUS_OK) then
             idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                         "kim_api_setm_data", ier)
             stop
            endif
            call kim_api_setm_data(pkim, ier, &
            "coordinates",  DIM*aS%numParticles   , c_loc(aS%coords),         1,                 &
            "boxSideLengths",              DIM, c_loc(aS%boxSideLengths), pbcFlag,           &
            "energy",                      SizeONe, c_loc(aS%energy),         1)

            if (ier.lt.KIM_STATUS_OK) then
             idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                         "kim_api_setm_data", ier)
             stop
            endif

            ier = kim_api_set_data(pkim, "neighObject", SizeOne, c_loc(aS%neighObject))
            if (ier.lt.KIM_STATUS_OK) then
                idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                         "kim_api_set_data", ier)
                stop
            endif

            ier = kim_api_set_method(pkim, "get_neigh", SizeOne, c_funloc(get_neigh))
            if (ier.lt.KIM_STATUS_OK) then
                idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                         "kim_api_set_method", ier)
                stop
            endif

            if (aS%processdEdrFlag .eq. 1) then
                ! Allocate and store pointers to process_dedr and process_d2edr2 subroutines
                ier = kim_api_set_method(pkim, "process_dEdr", SizeOne, c_funloc(process_dedr))
                if (ier.lt.KIM_STATUS_OK) then
                   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                              "kim_api_set_method", ier)
                   stop
                endif
            endif

            if (aS%processd2Edr2Flag .eq. 1) then
                ier = kim_api_set_method(pkim, "process_d2Edr2", SizeOne, c_funloc(process_d2edr2))
                if (ier.lt.KIM_STATUS_OK) then
                   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                              "kim_api_set_method", ier)
                   stop
                endif
            endif

        end subroutine broadcastAtomistic_kim

        subroutine setCompute_kim(pkim,nm,on)
            implicit none

            ! Passed variables
            type(c_ptr), intent(in)     :: pkim
            character(len=*),intent(in) :: nm
            integer(c_int),intent(in)   :: on

            call kim_api_set_compute(pkim,trim(nm),TRUEFALSE(on .eq. 1),ier)
            if (ier.lt.KIM_STATUS_OK) then
               idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                          "kim_api_set_compute", ier)
               stop
            endif

        end subroutine setCompute_kim

        function get_neigh(pkim,mode,request,part,numnei,pnei1part, &
                                      pRij) bind(c)
        use mod_atomistic, only : neighObject_type
        implicit none

        !-- Transferred variables
        type(c_ptr),    intent(in)  :: pkim
        integer(c_int), intent(in)  :: mode
        integer(c_int), intent(in)  :: request
        integer(c_int), intent(out) :: part
        integer(c_int), intent(out) :: numnei
        type(c_ptr),    intent(out) :: pnei1part
        type(c_ptr),    intent(out) :: pRij
        integer(c_int)              :: get_neigh

        !-- Local variables
        integer(c_int), save :: iterVal = 0
        integer(c_int)  N
        integer(c_int)  partToReturn
        integer(c_int), pointer :: numberOfParticles;   type(c_ptr) :: pnParts
        type(neighObject_type), pointer :: neighObject; type(c_ptr) :: pneighObject
        integer(c_int)  ier, idum

        ! unpack number of particles
        pnParts = kim_api_get_data(pkim, "numberOfParticles", ier)
        if (ier.lt.KIM_STATUS_OK) then
         idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                     "kim_api_get_data", ier)
         stop
        endif
        call c_f_pointer(pnParts, numberOfParticles)

        ! unpack neighbor list object
        pneighObject = kim_api_get_data(pkim, "neighObject", ier)
        if (ier.lt.KIM_STATUS_OK) then
         idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                     "kim_api_get_data", ier)
         stop
        endif
        call c_f_pointer(pneighObject, neighObject)

        N = size(neighObject%neighborList, 2)

        ! check mode and request
        if (mode.eq.0) then ! iterator mode
         if (request.eq.0) then ! reset iterator
            iterVal = 0
            get_neigh = KIM_STATUS_NEIGH_ITER_INIT_OK
            return
         elseif (request.eq.1) then ! increment iterator
            iterVal = iterVal + 1
            if (iterVal.gt.N) then
               get_neigh = KIM_STATUS_NEIGH_ITER_PAST_END
               return
            else
               partToReturn = iterVal
            endif
         else
            idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                        "Invalid request in get_neigh", &
                                        KIM_STATUS_NEIGH_INVALID_REQUEST)
            get_neigh = KIM_STATUS_NEIGH_INVALID_REQUEST
            return
         endif
        elseif (mode.eq.1) then ! locator mode
         if ( (request.gt.N) .or. (request.lt.1)) then
            idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                        "Invalid part ID in get_neigh", &
                                        KIM_STATUS_PARTICLE_INVALID_ID)
            get_neigh = KIM_STATUS_PARTICLE_INVALID_ID
            return
         else
            partToReturn = request
         endif
        else ! not iterator or locator mode
         idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                     "Invalid mode in get_neigh", &
                                     KIM_STATUS_NEIGH_INVALID_MODE)
         get_neigh = KIM_STATUS_NEIGH_INVALID_MODE
         return
        endif

        ! set the returned part
        part = partToReturn

        ! set the returned number of neighbors for the returned part
        numnei = neighObject%neighborList(1,part)

        ! set the location for the returned neighbor list
        pnei1part = c_loc(neighObject%neighborList(2,part))

        ! set pointer to Rij to appropriate value
        if (allocated(neighObject%RijList)) then
        pRij = c_loc(neighObject%RijList(1,1,part))
        else
        pRij = c_null_ptr
        endif

        get_neigh = KIM_STATUS_OK
        return
        end function get_neigh

        subroutine compute_kim(pkim)
            implicit none
            !-- Transferred variables
            type(c_ptr)      :: pkim

            ! Call model compute
            ier = kim_api_model_compute(pkim)
            if (ier.lt.KIM_STATUS_OK) then
                idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                         "kim_api_model_compute", ier)
                stop
            endif
        end subroutine compute_kim

        subroutine destroy_kim(pkim)
            implicit none
            type(c_ptr),    intent(in)  :: pkim

            ier = kim_api_model_destroy(pkim)
            if (ier.lt.KIM_STATUS_OK) then
                idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                         "kim_api_model_destroy", ier)
                stop
            endif

            call kim_api_free(pkim, ier)
            if (ier.lt.KIM_STATUS_OK) then
                idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                         "kim_api_free", ier)
                stop
            endif
        end subroutine destroy_kim

        subroutine write_results_kim(species,cS,cT,fsgeType)
            use mod_fsge, only : voigtStrainGradient,voigtStrain,Fsge
            use mod_crystal, only : Crystal
            implicit none

            ! Passed variables
            character(len=*),intent(in) :: cS
            type(Crystal),intent(in)    :: cT
            type(Fsge),intent(in)       :: fsgeType
            character(len=2),intent(in) :: species

            ! Local variables
            integer(c_int)     :: runit,m,n,i
            integer,external   :: GetFreeUnit
            character(len=1)   :: cr
            character(len=30)  :: cvalstr,dvalstr
            character(len=2)   :: mstr,nstr
            real(c_double)     :: a,c
            integer(c_int)     :: cpntsD(6),cpntsC(4)
            character(len=2000):: basisAtomsCoordsFractionString,   &
                                  latticeConstantString,            &
                                  crystalString,                    &
                                  speciesString,                    &
                                  CString,                        &
                                  DString,&
                                  propertyString,&
                                  tempstr

            cr = char(10)

            crystalString = '           "'  // trim(cS) //'"'  // cr

            if (trim(cS) .eq. 'fcc') then
                a = cT%LatticeType%latVec(2,1)*2._cd
                write(basisAtomsCoordsFractionString,*) &
                '            ['                                // cr //  &
                '                0'                            // cr //  &
                '                0'                            // cr //  &
                '                0'                            // cr //  &
                '            ]'                                // cr //  &
                '            ['                                // cr //  &
                '                0'                            // cr //  &
                '                0.5'                          // cr //  &
                '                0.5'                          // cr //  &
                '            ]'                                // cr //  & 
                '            ['                                // cr //  &
                '                0.5 '                         // cr //  &
                '                0 '                           // cr //  &
                '                0.5'                          // cr //  &
                '            ]'                                // cr //  &
                '            ['                                // cr //  &
                '                0.5'                          // cr //  &
                '                0.5'                          // cr //  &
                '                0'                            // cr //  &
                '            ]'                                // cr

                write(speciesString,*)                          &
                '            "' // species //'"' //  cr  //     &
                '            "' // species //'"' //  cr  //     &
                '            "' // species //'"' //  cr  //     &
                '            "' // species //'"' //  cr

                write(propertyString,*) &
                '    "property-id" "tag:staff@noreply.openkim.org,2016-05-24:property/'
                propertyString = trim(adjustl(propertyString)) // &
                'elastic-constants-first-strain-gradient-isothermal-cubic-crystal-npt"' // cr

            elseif(trim(cS) .eq. 'bcc') then
                a = cT%LatticeType%latVec(2,1)*2._cd
                write(basisAtomsCoordsFractionString,*) &
                '            ['                                // cr //  &
                '                0'                            // cr //  &
                '                0'                            // cr //  &
                '                0'                            // cr //  &
                '            ]'                                // cr //  &
                '            ['                                // cr //  &
                '                0.5'                          // cr //  &
                '                0.5'                          // cr //  &
                '                0.5'                          // cr //  &
                '            ]'                                // cr

                write(speciesString,*)                          &
                '            "' // species //'"' //  cr  //     &
                '            "' // species //'"' //  cr 

                write(propertyString,*) &
                '    "property-id" "tag:staff@noreply.openkim.org,2016-05-24:property/'
                propertyString = trim(adjustl(propertyString)) // &
                'elastic-constants-first-strain-gradient-isothermal-cubic-crystal-npt"' // cr

            elseif(trim(cS) .eq. 'hex') then
                c = cT%LatticeType%latVec(3,3)
                write(basisAtomsCoordsFractionString,*) &
                '            ['                                // cr //  &
                '                0'                            // cr //  &
                '                0'                            // cr //  &
                '                0'                            // cr //  &
                '            ]'                                // cr

                write(speciesString,*)                          &
                '            "' // species //'"' //  cr

                write(propertyString,*) &
                '    "property-id" "tag:staff@noreply.openkim.org,2016-05-24:property/'
                propertyString = trim(adjustl(propertyString)) // &
                'elastic-constants-first-strain-gradient-isothermal-monoatomic-hexagonal-crystal-npt"' // cr
            endif

            write(tempstr,*) a
            write(latticeConstantString,*) &
            '    "a" {'                                       // cr //  &
            '        "source-value" '//trim(adjustl(tempstr)) // cr //  &
            '        "source-unit" "angstrom"'                // cr //  &
            '    }'                                           // cr
            if (trim(cS) .eq. 'hex' .or. trim(cS) .eq. 'hcp') then
                write(tempstr,*) c
                latticeConstantString = trim(latticeConstantString)     //  &
                '    "c" {'                                       // cr //  &
                '        "source-value" '//trim(adjustl(tempstr)) // cr //  &
                '        "source-unit" "angstrom"'                // cr //  &
                '    }'                                           // cr
            endif


            CString=''
            do i=1,size(fsgeType%IndComponentsC,1)
                m = fsgeType%IndComponentsC(i,1); write(mstr,'(i2)') m
                n = fsgeType%IndComponentsC(i,2); write(nstr,'(i2)') n
                cpntsC(1:2) = voigtStrain(m)
                cpntsC(3:4) = voigtStrain(n)
                write(cvalstr,*) fsgeType%Cijkl(cpntsC(1),cpntsC(2),cpntsC(3),cpntsC(4))
                write(tempstr,*)&
                '    "c' //trim(adjustL(mstr)) // trim(adjustl(nstr))   // '" {'  //  cr //  &
                '    "source-value" '   //  trim(cvalstr)                         //  cr //  &
                '    "source-unit"  "eV/angstrom^3"'                              //  cr //  &
                '    }'                                                           //  cr
                CString = trim(CString) // trim(tempstr)
            enddo

            DString=''
            do i=1,size(fsgeType%IndComponentsD,1)
                m = fsgeType%IndComponentsD(i,1); write(mstr,'(i2)') m
                n = fsgeType%IndComponentsD(i,2); write(nstr,'(i2)') n
                cpntsD(1:3) = voigtStrainGradient(m)
                cpntsD(4:6) = voigtStrainGradient(n)
                write(dvalstr,'(F20.8)') fsgeType%Dijmkln(cpntsD(1),cpntsD(2),cpntsD(3),cpntsD(4),cpntsD(5),cpntsD(6))
                write(tempstr,*) &
                '    "d-' // trim(adjustl(mstr)) // '-' // trim(adjustl(nstr)) // '" {'  //  cr //   &
                ' "source-value" '  //  trim(adjustl(dvalstr))                           //  cr //   &
                '    "source-unit"  "eV/angstrom"'                                       //  cr //   &
                '    }'                                                                  //  cr
                DString = trim(Dstring) // trim(tempstr)
            enddo



            ! Now write the results.edn file
            runit = GetFreeUnit()
            open(unit=runit,file='output/results.edn',status='replace')
            write(runit,*) &
    '{'                                            // cr //  &
    !!!!!!!!!!!!!!!!
    '   "short-name" {'                            // cr //  &
    '       "source-value" ['                      // cr //  &
             trim(crystalString)                   //        &
    '       ]'                                     // cr //  &
    '   }'                                         // cr //  &
    !!!!!!!!!!!!!!!!
    trim(latticeConstantString)                    //        &
    !!!!!!!!!!!!!!!!
    trim(propertyString)                           //        &
    !!!!!!!!!!!!!!!!
    '    "basis-atom-coordinates" {'               // cr //  &
    '        "source-value" ['                     // cr //  &
          trim(basisAtomsCoordsFractionString)           //  &
    '        ]'                                    // cr //  &
    '    }'                                        // cr //  &
    !!!!!!!!!!!!!!!!
    '    "species" {'                              // cr //  &
    '        "source-value" ['                     // cr //  &
              trim(speciesString)                        //  & 
    '        ]'                                    // cr //  &
    '    }'                                        // cr //  &
    !!!!!!!!!!!!!!!!
    '    "temperature" {'                          // cr //  &
    '    "source-value" 0'                         // cr //  &
    '    "source-unit"  "K"'                       // cr //  &
    '    }'                                        // cr //  &
    !!!!!!!!!!!!!!!!
    trim(CString)                                  //        &
    !!!!!!!!!!!!!!!!
    trim(DString)                                  //        &
    !!!!!!!!!!!!!!!!
    '    "instance-id" 1'                          // cr //  &
    '}'
        end subroutine write_results_kim

            function getSpeciesCode_kim(pkim,species)
            implicit none
            type(c_ptr),     intent(in) :: pkim
            character(len=2),intent(in) :: species
            integer(c_int)              :: getSpeciesCode_kim

            getSpeciesCode_kim = kim_api_get_species_code(pkim, species, ier)
            if (ier.lt.KIM_STATUS_OK) then
                idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                         "kim_api_get_species_code", ier)
                stop
            endif
            return
            end function getSpeciesCode_kim

        function getCompute_kim(pkim,nm)
            implicit none

            ! Passed variables
            type(c_ptr), intent(in)     :: pkim
            character(len=*),intent(in) :: nm
            integer(c_int)              :: getCompute_kim

            getCompute_kim = kim_api_get_compute(pkim,trim(nm),ier)
            if (ier.lt.KIM_STATUS_OK) then
               idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                          "kim_api_get_compute", ier)
               stop
            endif
        end function getCompute_kim

subroutine WriteDescriptorFile(pkim,species,modelname)
use, intrinsic :: iso_c_binding
implicit none

!-- Transferred variables
type(c_ptr)      :: pkim
character(len=KIM_KEY_STRING_LENGTH),intent(in) :: modelname
character(len=2),     intent(in)   :: species

!-- Local variables
integer(c_int) :: i
character(len=10000)   :: descriptorfile
character(len=103) :: divider
character(len=1)   :: cr
character(len=52)  :: type_line

! Define frequently used variables
!
cr = char(10)
divider = '#######################################################################################################'

! Write Minimol descriptor file into string descriptorfile
!
descriptorfile = &
   divider                                                                            // cr // &
   '#'                                                                                // cr // &
   '# Copyright 2015 Nikhil Chandra Admal and Giacomo Po'                             // cr // &
   '# All rights reserved.'                                                           // cr // &
   '#'                                                                                // cr // &
   '#'                                                                                // cr // &
   '# See KIM_API/standard.kim for documentation about this file'                     // cr // &
   '#'                                                                                // cr // &
   divider                                                                            // cr // &
                                                                                         cr // &
                                                                                         cr // &
   'KIM_API_Version := 1.6.0'                                                         // cr // &
   'Unit_length      := A'                                                            // cr // &
   'Unit_energy      := eV'                                                           // cr // &
   'Unit_charge      := e'                                                            // cr // &
   'Unit_temperature := K'                                                            // cr // &
   'Unit_time        := ps'                                                           // cr // &
                                                                                         cr // &
                                                                                         cr // &
   divider                                                                            // cr // &
   'PARTICLES_SPECIES:'                                                               // cr // &
   '# Symbol/name           Type                    code'                             // cr

   write(type_line,'(a2,20x,''spec'',20x,i4)') species,0
   descriptorfile = trim(descriptorfile) // type_line // cr

   descriptorfile = trim(descriptorfile) // & 
                                                                                         cr // &
                                                                                         cr // &
   divider                                                                            // cr // &
   'CONVENTIONS:'                                                                     // cr // &
   '# Name                  Type'                                                     // cr // &
                                                                                         cr // &
   'OneBasedLists           flag'                                                     // cr // &
                                                                                         cr // &
   'Neigh_BothAccess        flag'                                                     // cr // &
                                                                                         cr // &
   'NEIGH_RVEC_H            flag'                                                     // cr // &
                                                                                         cr // &
   'NEIGH_PURE_H            flag'                                                     // cr // &
                                                                                         cr // &
   'NEIGH_RVEC_F            flag'                                                     // cr // &
                                                                                         cr // &
   'NEIGH_PURE_F            flag'                                                     // cr // &
                                                                                         cr // &
   'MI_OPBC_H               flag'                                                     // cr // &
                                                                                         cr // &
   'MI_OPBC_F               flag'                                                     // cr // &
                                                                                         cr // &
                                                                                         cr // &
   divider                                                                            // cr // &
   'MODEL_INPUT:'                                                                     // cr // &
   '# Name                        Type         Unit       Shape              requirements' // cr // &
                                                                                         cr // &
                                                                                         cr // &
   'numberOfParticles             integer      none       []'                         // cr // &
                                                                                         cr // &
   'numberContributingParticles   integer      none       []'                         // cr // &
                                                                                         cr // &
   'numberOfSpecies               integer      none       []'                         // cr // &
                                                                                         cr // &
   'particleSpecies               integer      none       [numberOfParticles]'        // cr // &
                                                                                         cr // &
   'coordinates                   double       length     [numberOfParticles,3]'      // cr // &
                                                                                         cr // &
   'get_neigh                     method       none       []'                         // cr // &
                                                                                         cr // &
   'neighObject                   pointer      none       []'                         // cr // &
                                                                                         cr // &
   'boxSideLengths                double       length     [3]'                        // cr // &
                                                                                         cr
   if (analytical .eq. 1) then
       descriptorfile = trim(descriptorfile)                                                // &
   'process_dEdr                  method       none       []'                         // cr // &
                                                                                         cr // &
   'process_d2Edr2                method       none       []'                         // cr
   endif
   
   descriptorfile = trim(descriptorfile)                                              // cr // &
                                                                                         cr // &
   divider                                                                            // cr // &
   'MODEL_OUTPUT:'                                                                    // cr // &
   '# Name                        Type         Unit       Shape              requirements' // cr // &
                                                                                         cr // &
   'destroy                       method       none       []'                         // cr // &
                                                                                         cr // &
   'compute                       method       none       []'                         // cr // &
                                                                                         cr // &
   'cutoff                        double       length     []'                         // cr // &
                                                                                         cr // &
   'energy                        double       energy     []'                         // cr // &
                                                                                         cr // &
   divider                                                                            // cr // char(0)

ier = kim_api_string_init(pkim, trim(descriptorfile), modelname)
if (ier.lt.KIM_STATUS_OK) then
    idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                             "kim_api_file_init", ier)
    stop
endif

return

end subroutine WriteDescriptorFile

end module mod_kim

integer function GetFreeUnit()
    implicit none
    logical InUse

    do GetFreeUnit=7,98
        inquire(unit=GetFreeUnit,opened=InUse)
        if (.not. InUse) return
    enddo
    write(*,*) "Could not obtain a free unit handle in function GetFreeUnit"
    stop
end function GetFreeUnit


