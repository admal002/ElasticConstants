!
! CDDL HEADER START
!
! The contents of this file are subject to the terms of the Common Development
! and Distribution License Version 1.0 (the "License").
!
! You can obtain a copy of the license at
! http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
! specific language governing permissions and limitations under the License.
!
! When distributing Covered Code, include this CDDL HEADER in each file and
! include the License file in a prominent location with the name LICENSE.CDDL.
! If applicable, add the following below this CDDL HEADER, with the fields
! enclosed by brackets "[]" replaced with your own identifying information:
!
! Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
!
! CDDL HEADER END
!
!
! Copyright (c) 2015--2016, University of Californa Los Angeles
! All rights reserved.
!
! Contributors:
!    Nikhil Chandra Admal
!

!
! Release: This file is part of the kim-api.git package.
!
!-------------------------------------------------------------------------------
! Notes:
!
! 1. The test works for FCC and BCC crystal structures. It can be easily
!    generalized to any simple lattice. 
! 2. The size (nCellsPerSide) is taken to be equal to seven. This size is fixed.
!    It is made large enough to ensure that the central atoms truely represents a
!    bult atom. In other words, this depends on the cutoff of the model. Seven seems
!    to be a large number for most models. 

#include "KIM_API_status.h"
#define THIS_FILE_NAME __FILE__
#define TRUEFALSE(TRUTH) merge(1,0,(TRUTH))

!-------------------------------------------------------------------------------
! module mod_global :
!     Global variables
!
module mod_global
    use, intrinsic            :: iso_c_binding
    implicit none
    integer(c_int), parameter :: DIM=3
    integer(c_int), parameter :: cd = c_double ! used for literal constants
end module mod_global

!-------------------------------------------------------------------------------
! module mod_coordinates:
!     Contains coordinates and the middle particle ID of the sysetem
!
module mod_coordinates
    use mod_global
    implicit none
    real(c_double), pointer :: coords(:,:); type(c_ptr) :: pcoor
    integer(c_int) :: MiddlePartId
end module mod_coordinates

!-------------------------------------------------------------------------------
!
! module mod_neighborlist :
!
!    Module contains type and routines related to neighbor list calculation
!
!-------------------------------------------------------------------------------

module mod_neighborlist
    use KIM_API_F03
    use mod_global

    public  get_neigh

    type neighObject_type
    integer(c_int), pointer :: neighborList(:,:)
    real(c_double), pointer :: RijList(:,:,:)
    end type neighObject_type

contains

!-------------------------------------------------------------------------------
!
! get_neigh neighbor list access function
!
! This function implements Locator and Iterator mode
!
!-------------------------------------------------------------------------------
integer(c_int) function get_neigh(pkim,mode,request,part,numnei,pnei1part, &
                                  pRij) bind(c)
  use mod_global
  implicit none

  !-- Transferred variables
  type(c_ptr),    intent(in)  :: pkim
  integer(c_int), intent(in)  :: mode
  integer(c_int), intent(in)  :: request
  integer(c_int), intent(out) :: part
  integer(c_int), intent(out) :: numnei
  type(c_ptr),    intent(out) :: pnei1part
  type(c_ptr),    intent(out) :: pRij

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
  if (associated(neighObject%RijList)) then
    pRij = c_loc(neighObject%RijList(1,1,part))
  else
    pRij = c_null_ptr
  endif

  get_neigh = KIM_STATUS_OK
  return
end function get_neigh

end module mod_neighborlist

module mod_process
use mod_global
implicit none

real(c_double) :: Cijkl(DIM,DIM,DIM,DIM),Dijmkln(DIM,DIM,DIM,DIM,DIM,DIM)
real(c_double) :: Cijkl_numerical(DIM,DIM,DIM,DIM),Dijmkln_numerical(DIM,DIM,DIM,DIM,DIM,DIM)

contains
    subroutine setup_process_KIM_access(pkim)
    use KIM_API_F03
    implicit none
    ! Passed variables
    type(c_ptr), intent(in) :: pkim

    !-- Local variables
    integer(c_int), parameter :: SizeOne = 1
    integer(c_int)            :: ier, idum

    ! Allocate and store pointers to process_dedr and process_d2edr2 subroutines
    ier = kim_api_set_method(pkim, "process_dEdr", SizeOne, c_funloc(process_dedr))
    if (ier.lt.KIM_STATUS_OK) then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                  "kim_api_set_method", ier)
       stop
    endif

    ier = kim_api_set_method(pkim, "process_d2Edr2", SizeOne, c_funloc(process_d2edr2))
    if (ier.lt.KIM_STATUS_OK) then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                  "kim_api_set_method", ier)
       stop
    endif
    end subroutine setup_process_KIM_access

    function process_dedr(pkim,deidr,r,pdx,i,j)bind(c)
    use KIM_API_F03
    use mod_coordinates
    implicit none

    ! Passed variables
    type(c_ptr),   intent(in) :: pkim,pdx
    real(c_double),intent(in) :: deidr,r
    integer(c_int),intent(in) :: i,j
    integer(c_int)            :: process_dedr

    ! Local variables
    real(c_double), pointer :: dx(:)
    real(c_double) :: R3(DIM)

    call c_f_pointer(pdx,dx,[DIM])

    R3 = (coords(:,i) + coords(:,j))/2._cd - coords(:,MiddlePartId)

    Cijkl = Cijkl - deidr * tensor_prod(dx,dx,dx,dx)/(r**3)
    Dijmkln = Dijmkln - deidr * tensor_prod_6(dx,dx,R3,dx,dx,R3)/(r**3)
    Dijmkln = Dijmkln + deidr * operatorT( tensor_prod(dx,dx,R3,R3)/(2._cd*r) )

    process_dedr = KIM_STATUS_OK

    end function process_dedr

    function process_d2edr2(pkim,d2edr2,pr,pdx,pi,pj)bind(c)
    use KIM_API_F03
    use mod_coordinates
    implicit none

    ! Passed variables
    type(c_ptr),   intent(in) :: pkim,pr,pdx
    real(c_double),intent(in) :: d2edr2
    type(c_ptr),   intent(in) :: pi,pj
    integer(c_int)            :: process_d2edr2

    ! Local variables
    integer(c_int), pointer :: i(:);
    integer(c_int), pointer :: j(:);
    real(c_double), pointer :: dx(:,:)
    real(c_double), pointer :: r(:)
    real(c_double) :: R1(DIM),R2(DIM),R3(DIM),R4(DIM)

    call c_f_pointer(pi, i, [2])
    call c_f_pointer(pj, j, [2])
    call c_f_pointer(pdx,dx,[DIM,2])
    call c_f_pointer(pr,r,[2])


    R1 = dx(:,1)
    R2 = dx(:,2)
    R3 = (coords(:,i(1)) + coords(:,j(1)))/2._cd - coords(:,MiddlePartId)
    R4 = (coords(:,i(2)) + coords(:,j(2)))/2._cd - coords(:,MiddlePartId)

    Cijkl = Cijkl + d2edr2*tensor_prod(R1,R1,R2,R2)/(r(1)*r(2))
    Dijmkln= Dijmkln + d2edr2*tensor_prod_6(R1,R1,R3,R2,R2,R4)/(r(1)*r(2))


    process_d2edr2 = KIM_STATUS_OK
    end function process_d2edr2

    function operatorT(tensor4)
    implicit none

    ! Passed variables
    real(c_double) :: operatorT(DIM,DIM,DIM,DIM,DIM,DIM)
    real(c_double),intent(in) :: tensor4(DIM,DIM,DIM,DIM)

    ! Local variables
    integer(c_int) :: i,j,m,k,l,n

    operatorT = 0._cd
    do n=1,DIM
       do l=1,DIM
           do k=1,DIM
               do m=1,DIM
                   do j=1,DIM
                       do i=1,DIM
                           if (i .eq. k) then
                               operatorT(i,j,m,k,l,n) = operatorT(i,j,m,k,l,n) + tensor4(l,m,j,n) &
                                                                           + tensor4(j,n,l,m) &
                                                                           + tensor4(j,l,m,n) &
                                                                           + tensor4(m,n,j,l)
                           endif

                           if (i .eq. l) then
                               operatorT(i,j,m,k,l,n) = operatorT(i,j,m,k,l,n) + tensor4(j,k,m,n) &
                                                                           + tensor4(j,n,k,m) &
                                                                           + tensor4(k,m,j,n) &
                                                                           + tensor4(m,n,j,k)
                           endif

                           if (i .eq. n) then
                               operatorT(i,j,m,k,l,n) = operatorT(i,j,m,k,l,n) - tensor4(j,k,l,m) &
                                                                           - tensor4(j,l,k,m) &
                                                                           - tensor4(k,m,j,l) &
                                                                           - tensor4(l,m,j,k)
                           endif

                           if (j .eq. k) then
                               operatorT(i,j,m,k,l,n) = operatorT(i,j,m,k,l,n) + tensor4(l,m,i,n) &
                                                                           + tensor4(i,n,l,m) &
                                                                           + tensor4(i,l,m,n) &
                                                                           + tensor4(m,n,i,l)
                           endif

                           if (j .eq. l) then
                               operatorT(i,j,m,k,l,n) = operatorT(i,j,m,k,l,n) + tensor4(i,k,m,n) &
                                                                           + tensor4(i,n,k,m) &
                                                                           + tensor4(k,m,i,n) &
                                                                           + tensor4(m,n,i,k)
                           endif

                           if (j .eq. n) then
                               operatorT(i,j,m,k,l,n) = operatorT(i,j,m,k,l,n) - tensor4(i,k,l,m) &
                                                                           - tensor4(i,l,k,m) &
                                                                           - tensor4(k,m,i,l) &
                                                                           - tensor4(l,m,i,k)
                           endif

                           if (m .eq. k) then
                               operatorT(i,j,m,k,l,n) = operatorT(i,j,m,k,l,n) - tensor4(i,l,j,n) &
                                                                           - tensor4(i,n,j,l) &
                                                                           - tensor4(j,l,i,n) &
                                                                           - tensor4(j,n,i,l)
                           endif
                           if (m .eq. l) then
                               operatorT(i,j,m,k,l,n) = operatorT(i,j,m,k,l,n) - tensor4(i,k,j,n) &
                                                                           - tensor4(i,n,j,k) &
                                                                           - tensor4(j,k,i,n) &
                                                                           - tensor4(j,n,i,k)
                           endif
                           if (m .eq. n) then
                               operatorT(i,j,m,k,l,n) = operatorT(i,j,m,k,l,n) + tensor4(i,k,j,l) &
                                                                           + tensor4(i,l,j,k) &
                                                                           + tensor4(j,k,i,l) &
                                                                           + tensor4(j,l,i,k)
                           endif

                       enddo
                   enddo
               enddo
           enddo
       enddo
    enddo

    operatorT = 0.5_cd * operatorT
    end function operatorT

    function tensor_prod(vec1,vec2,vec3,vec4)
    implicit none

    ! Passed variables
    real(c_double),intent(in) :: vec1(DIM),vec2(DIM),vec3(DIM),vec4(DIM)
    real(c_double)            :: tensor_prod(DIM,DIM,DIM,DIM)

    ! Local variables
    integer(c_int) :: p,q,r,s
    do s=1,DIM
        do r=1,DIM
            do q=1,DIM
                do p=1,DIM
                    tensor_prod(p,q,r,s) = vec1(p)*vec2(q)*vec3(r)*vec4(s)
                enddo
            enddo
        enddo
    enddo
    return
    end function tensor_prod

    function tensor_prod_6(vec1,vec2,vec3,vec4,vec5,vec6)
    implicit none
    real(c_double),intent(in) :: vec1(DIM),vec2(DIM),vec3(DIM),vec4(DIM),vec5(DIM),vec6(DIM)
    real(c_double) :: tensor_prod_6(DIM,DIM,DIM,DIM,DIM,DIM)

    integer(c_int) :: p,q,r,s,t,u
    do u=1,DIM
        do t=1,DIM
            do s=1,DIM
                do r=1,DIM
                    do q=1,DIM
                        do p=1,DIM
                            tensor_prod_6(p,q,r,s,t,u) = vec1(p)*vec2(q)*vec3(r)*vec4(s)*vec5(t)*vec6(u)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    return
    end function tensor_prod_6

subroutine write_results
    implicit none

    ! Local variables
    integer(c_int)             :: cunit,dunit,runit
    integer,external           :: GetFreeUnit
    character(len=1)           :: cr

    cr = char(10)

    cunit = GetFreeUnit()
    open(unit=cunit,file='output/Cijkl.out',status='replace')
    dunit = GetFreeUnit()
    open(unit=dunit,file='output/Dijmkln.out',status='replace')
    
    write(cunit, *) "c11 = ",Cijkl(1,1,1,1),Cijkl_numerical(1,1,1,1)
    write(cunit, *) "c12 = ",Cijkl(1,1,2,2),Cijkl_numerical(1,1,2,2)
    write(cunit, *) "c44 = ",Cijkl(1,2,1,2),Cijkl_numerical(1,2,1,2)

    write(dunit,*) "D111111 = ",Dijmkln(1,1,1,1,1,1),Dijmkln_numerical(1,1,1,1,1,1)
    write(dunit,*) "D111122 = ",Dijmkln(1,1,1,1,2,2),Dijmkln_numerical(1,1,1,1,2,2)
    write(dunit,*) "D111221 = ",Dijmkln(1,1,1,2,2,1),Dijmkln_numerical(1,1,1,2,2,1)
    write(dunit,*) "D112112 = ",Dijmkln(1,1,2,1,1,2),Dijmkln_numerical(1,1,2,1,1,2)
    write(dunit,*) "D112121 = ",Dijmkln(1,1,2,1,2,1),Dijmkln_numerical(1,1,2,1,2,1)
    write(dunit,*) "D112233 = ",Dijmkln(1,1,2,2,3,3),Dijmkln_numerical(1,1,2,2,3,3)
    write(dunit,*) "D112332 = ",Dijmkln(1,1,2,3,3,2),Dijmkln_numerical(1,1,2,3,3,2)
    write(dunit,*) "D121121 = ",Dijmkln(1,2,1,1,2,1),Dijmkln_numerical(1,2,1,1,2,1)
    write(dunit,*) "D121233 = ",Dijmkln(1,2,1,2,3,3),Dijmkln_numerical(1,2,1,2,3,3)
    write(dunit,*) "D123123 = ",Dijmkln(1,2,3,1,2,3),Dijmkln_numerical(1,2,3,1,2,3)

    close(cunit)
    close(dunit)

    ! Now write the results.edn file
    runit = GetFreeUnit()
    open(unit=runit,file='output/results.edn',status='replace')
    write(runit,*) &
    '{'                               // cr //  &
    '   "short-name" {'               // cr //  &
    '       "source-value" ['         // cr //  &
    '           "fcc"'                // cr //  &
    '       ]'                        // cr //  &
    '   }'                            // cr //  &
    '    "a" {'                                    // cr //  &
    '        "si-unit" "m"'                        // cr //  &
    '        "source-unit" "angstrom"'             // cr //  &
    '        "si-value" 5.24859e-10'               // cr //  &
    '        "source-value" 5.23232'               // cr //  &
    '    }'                                        // cr //  &
    '    "wyckoff-species" {'                      // cr //  &
    '        "source-value" ['                     // cr //  &
    '            "Ar" '                            // cr //  &
    '        ]'                                    // cr //  &
    '    }'                                        // cr //  &
    '    "wyckoff-multiplicity-and-letter" {'      // cr //  &
    '        "source-value" ['                     // cr //  &
    '            "4a"'                             // cr //  &
    '        ]'                                    // cr //  &
    '    }'                                        // cr //  &
    '    "property-id" "tag:staff@noreply.openkim.org,2014-04-15:property/cohesive-potential-energy-cubic-crystal"'    // cr //  &
    '    "space-group" {'                          // cr //  &
    '        "source-value" "Fm-3m"'               // cr //  &
    '    }'                                        // cr //  &
    '    "cohesive-potential-energy" {'            // cr //  &
    '        "si-unit" "kg m^2 / s^2"'             // cr //  &
    '        "source-unit" "eV"'                   // cr //  &
    '        "si-value" 1.3859709e-20'             // cr //  &
    '        "source-value" 0.0865055077405508'    // cr //  &
    '    }'                                        // cr //  &
    '    "basis-atom-coordinates" {'               // cr //  &
    '        "source-value" ['                     // cr //  &
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
    '            ]'                                // cr //  &
    '        ]'                                    // cr //  &
    '    }'                                        // cr //  &
    '    "wyckoff-coordinates" {'                  // cr //  &
    '        "source-value" ['                     // cr //  &
    '            ['                                // cr //  &
    '                0'                            // cr //  &
    '                0'                            // cr //  &
    '                0'                            // cr //  &
    '            ]'                                // cr //  &
    '        ]'                                    // cr //  &
    '    }'                                        // cr //  &
    '    "species" {'                              // cr //  &
    '        "source-value" ['                     // cr //  &
    '            "Ar"'                             // cr //  &
    '            "Ar"'                             // cr //  &
    '            "Ar"'                             // cr //  &
    '            "Ar"'                             // cr //  &
    '        ]'                                    // cr //  &
    '    }'                                        // cr //  &
    '    "instance-id" 1'                          // cr //  &
'}'
end subroutine write_results

end module mod_process

!*******************************************************************************
!**
!**  PROGRAM runner
!**
!**  KIM compliant program to compute the energy of and forces and virial on an
!**  isolated cluster of Al particles
!**
!**  Works with the following NBC methods:
!**        NEIGH_RVEC_H
!**        NEIGH_PURE_H
!**        NEIGH_RVEC_F
!**        NEIGH_PURE_F
!**        MI_OPBC_H
!**        MI_OPBC_F
!**
!*******************************************************************************

!-------------------------------------------------------------------------------
!
! Main program
!
!-------------------------------------------------------------------------------
program runner
    use mod_global
    use KIM_API_F03
    use mod_coordinates
    use mod_neighborlist
    use mod_process, only : Cijkl,Cijkl_numerical,      &
                            Dijmkln,Dijmkln_numerical,  &
                            write_results,              &
                            setup_process_KIM_access
    implicit none

    integer(c_int), parameter :: nCellsPerSide  = 8
    integer(c_int), parameter :: ASpecies       = 1
    integer(c_int), parameter :: SizeOne        = 1

    integer(c_int)   :: N,dunit
    real(c_double)   :: LatticeSpacing,AtomVol
    character(len=3) :: CrystalStructure
    character(len=2) :: Species
    integer,external :: GetFreeUnit

    !
    ! neighbor list
    !
    type(neighObject_type), target :: neighObject
    real(c_double), parameter      :: cutpad = 0.75_cd ! in angstroms

    !
    ! KIM variables
    !
    character(len=KIM_KEY_STRING_LENGTH) :: testname,DependentTestname
    character(len=10000) testkimfile
    character(len=KIM_KEY_STRING_LENGTH) :: modelname
    character(len=KIM_KEY_STRING_LENGTH) :: NBC_Method
    ! 0- NEIGH_RVEC_H, 1- NEIGH_PURE_H, 2- NEIGH_RVEC_F, 3- NEIGH_PURE_F,
    ! 4- MI_OPBC_H,    5- MI_OPBC_F
    integer(c_int) nbc

    type(c_ptr)    :: pkim
    integer(c_int) :: ier, idum
    integer(c_int), pointer :: numberOfParticles;   type(c_ptr) :: pnParts
    integer(c_int), pointer :: numContrib;          type(c_ptr) :: pnumContrib
    integer(c_int), pointer :: numberOfSpecies;     type(c_ptr) :: pnOfSpecies
    integer(c_int), pointer :: particleSpecies(:);  type(c_ptr) :: pparticleSpecies
    real(c_double), pointer :: cutoff;              type(c_ptr) :: pcutoff
    real(c_double), pointer :: energy;              type(c_ptr) :: penergy
    real(c_double), pointer :: boxSideLengths(:);   type(c_ptr) :: pboxSideLengths
    integer :: numParticles_PrimitiveCell,ind

    Cijkl   = 0._cd
    Dijmkln = 0._cd

    ! Get KIM Model name to use
    read(*,*,err=100) testname
    read(*,*,err=100) modelname

    read(*,*,err=100) Species
    read(*,*,err=100) CrystalStructure
    read(*,*,err=100) DependentTestname
    read(*,*,err=100) LatticeSpacing


    if (CrystalStructure .eq. 'fcc') then
      N = 4*(nCellsPerSide)**3
      AtomVol = (LatticeSpacing**3)/4._cd
      numParticles_PrimitiveCell = 1
    elseif (CrystalStructure .eq. 'bcc') then
      N = 2*(nCellsPerSide)**3
      AtomVol = (LatticeSpacing**3)/2._cd
      numParticles_PrimitiveCell = 1
    elseif (CrystalStructure .eq. 'diamond') then
      N = 5*(nCellsPerSide)**3
      AtomVol = (LatticeSpacing**3)/4._cd
      numParticles_PrimitiveCell = 2
    else
      write(*,*) "Error: Crystal structure ", CrystalStructure," not supported."
      stop
    endif

    ! Initialize the KIM object
    testkimfile = '../../te/' // trim(testname) // '/descriptor.kim'
    ier = kim_api_file_init(pkim, testkimfile, modelname)
    if (ier.lt.KIM_STATUS_OK) then
        idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                 "kim_api_file_init", ier)
        stop
    endif

    ! Write to dependencies.edn the dependent testname
    dunit = GetFreeUnit()
    open(unit=dunit,file='output/dependencies.edn',status='replace')
    write(dunit,*) "[ " // "'" // trim(DependentTestname) // "'" // " ]"
    close(dunit)

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

    ! Allocate memory via the KIM system
    call kim_api_allocate(pkim, N, ASpecies, ier)
    if (ier.lt.KIM_STATUS_OK) then
        idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                 "kim_api_allocate", ier)
        stop
    endif

    ! Allocate and store pointers to neighbor list object and access function
    allocate(neighObject%neighborList(N+1, N))
    if (nbc.eq.0 .or. nbc.eq.2) then
        allocate(neighObject%RijList(DIM,N+1, N))
    endif
    !
    ier = kim_api_set_data(pkim, "neighObject", SizeOne, c_loc(neighObject))
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

    ! call model's init routine
    ier = kim_api_model_init(pkim)
    if (ier.lt.KIM_STATUS_OK) then
        idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                 "kim_api_model_init", ier)
        stop
    endif

    call setup_process_KIM_access(pkim)

    ! Unpack data from KIM object
    !
    call kim_api_getm_data(pkim, ier, &
       "numberOfParticles",           pnParts,          1,                                   &
       "numberContributingParticles", pnumContrib,      TRUEFALSE((nbc.eq.0).or.(nbc.eq.1).or.(nbc.eq.4)), &
       "numberOfSpecies",             pnOfSpecies,      1,                                   &
       "particleSpecies",             pparticleSpecies, 1,                                   &
       "coordinates",                 pcoor,            1,                                   &
       "cutoff",                      pcutoff,          1,                                   &
       "boxSideLengths",              pboxSideLengths,  TRUEFALSE((nbc.eq.4).or.(nbc.eq.5)), &
       "energy",                      penergy,          1)                              
    if (ier.lt.KIM_STATUS_OK) then
     idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                 "kim_api_getm_data", ier)
     stop
    endif
    call c_f_pointer(pnParts, numberOfParticles)
    if ((nbc.eq.0).or.(nbc.eq.1).or.(nbc.eq.4)) call c_f_pointer(pnumContrib,numContrib)
    call c_f_pointer(pnOfSpecies, numberOfSpecies)
    call c_f_pointer(pparticleSpecies, particleSpecies, [N])
    call c_f_pointer(pcoor, coords, [DIM,N])
    call c_f_pointer(pcutoff, cutoff)

    if ((nbc.eq.4).or.(nbc.eq.5)) call c_f_pointer(pboxSideLengths, &
                                                 boxSideLengths, [DIM])
    call c_f_pointer(penergy, energy)

    ! Set values
    numberOfParticles = N
    if ((nbc.eq.0).or.(nbc.eq.1).or.(nbc.eq.4)) numContrib = 1
    numberOfSpecies = ASpecies
    particleSpecies(:)    = kim_api_get_species_code(pkim, Species, ier)
    if (ier.lt.KIM_STATUS_OK) then
        idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                 "kim_api_get_species_code", ier)
        stop
    endif

    ! set up the cluster part positions
    if (CrystalStructure .eq. 'fcc') then
        call create_FCC_configuration(LatticeSpacing, nCellsPerSide, .false., coords, &
                                      MiddlePartId)
    elseif (CrystalStructure .eq. 'bcc') then
        call create_BCC_configuration(LatticeSpacing, nCellsPerSide, .false., coords, &
                                      MiddlePartId)
    elseif (CrystalStructure .eq. 'diamond') then
        call create_diamond_configuration(LatticeSpacing, nCellsPerSide,  coords )
    endif

    ! set boxSideLengths large enough to make the cluster isolated
    if (nbc.eq.4 .or. nbc.eq.5) boxSideLengths(:) = nCellsPerSide*LatticeSpacing



    ! The loop over the atoms of the unit cell starts here
    !do ind = 1,numParticles_PrimitiveCell
        ! Interchange coordinates of MiddlePartId with that of the index
        MiddlePartId = 1

        ! compute neighbor lists
        if (nbc.eq.0) then
            call NEIGH_RVEC_neighborlist(.true., N, coords, (cutoff+cutpad), &
                                              MiddlePartId, neighObject)
        elseif (nbc.eq.1) then
            call NEIGH_PURE_neighborlist(.true., N, coords, (cutoff+cutpad), &
                                              MiddlePartId, neighObject)
        elseif (nbc.eq.2) then
            call NEIGH_RVEC_neighborlist(.false., N, coords, (cutoff+cutpad), &
                                              MiddlePartId, neighObject)
        elseif (nbc.eq.3) then
            call NEIGH_PURE_neighborlist(.false., N, coords, (cutoff+cutpad), &
                                              MiddlePartId, neighObject)
        elseif (nbc.eq.4) then
            call MI_OPBC_neighborlist(.true., N, coords, (cutoff+cutpad), &
                                            boxSideLengths, MiddlePartId, neighObject)
        elseif (nbc.eq.5) then
            call MI_OPBC_neighborlist(.false., N, coords, (cutoff+cutpad), &
                                           boxSideLengths, MiddlePartId, neighObject)
        endif

        ! Call model compute
        ier = kim_api_model_compute(pkim)
        if (ier.lt.KIM_STATUS_OK) then
            idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                     "kim_api_model_compute", ier)
            stop
        endif
        ! print results to screen
        print '(80(''-''))'
        print '("This is Test             : ",A)', trim(testname)
        print '("Results for KIM Model    : ",A)', trim(modelname)
        print '("Crystal strucutre        : ",A)', CrystalStructure
        print '("Species                  : ",A)', species
        print '("Dependent Test           : ",A)', trim(DependentTestname)
        print '("Queried Lattice constant : ",F20.15)', LatticeSpacing
        print '("Using NBC                : ",A)', trim(NBC_Method)

        !-------------------------------------------------------------------------------
        ! Numerical checks
        !

        ! Switch off process_dedr and process_d2edr2
        call kim_api_set_compute(pkim,'process_dEdr',KIM_COMPUTE_FALSE,ier)
        call kim_api_set_compute(pkim,'process_d2Edr2',KIM_COMPUTE_FALSE,ier)
        if (ier.lt.KIM_STATUS_OK) then
            idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                     "kim_api_set_compute", ier)
            stop
        endif

        ! Perform the following non-uniform deformation which corresponts to zero strain
        ! and infinitesimal first strain gradient

        !call numerical_check_Dijmkln(pkim,nbc)
        call numerical_check_Cijkl(pkim,nbc)
    !enddo

    Cijkl             = Cijkl/AtomVol
    Dijmkln           = Dijmkln/AtomVol
    Cijkl_numerical   = Cijkl_numerical/AtomVol
    Dijmkln_numerical = Dijmkln_numerical/AtomVol

    write(*,*) Cijkl(1,1,1,1),"reached"
    stop

    call write_results

    ! Don't forget to free and/or deallocate
    deallocate(neighObject%neighborList)
    if (nbc.eq.0 .or. nbc.eq.2) then
        deallocate(neighObject%RijList)
    endif

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

    stop

    100 idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                 "Invalide query results", ier)
     stop

end program runner

!-------------------------------------------------------------------------------
!
! MI_OPBC_neighborlist : construct a half or full neighbor list using
!                                the particle coordinates in coords()
!
!-------------------------------------------------------------------------------
subroutine MI_OPBC_neighborlist(half, numberOfParticles, coords, rcut, &
                                        boxSideLengths, MiddlePartId, neighObject)
    use mod_global
    use KIM_API_F03
    use mod_neighborlist
    implicit none

    !-- Transferred variables
    logical,        intent(in)                  :: half
    integer(c_int), intent(in)                  :: numberOfParticles
    real(c_double), dimension(DIM,numberOfParticles), &
                  intent(in)                  :: coords
    real(c_double), intent(in)                  :: rcut
    real(c_double), dimension(DIM), intent(in)  :: boxSideLengths
    integer(c_int), intent(in)                  :: MiddlePartId
    type(neighObject_type),       intent(inout) :: neighObject

    !-- Local variables
    integer(c_int) i, j, a
    real(c_double) dx(DIM)
    real(c_double) r2
    real(c_double) rcut2

    rcut2 = rcut**2

    do i=1,numberOfParticles
        a = 1
        do j=1,numberOfParticles
            dx(:) = coords(:, j) - coords(:, i)
            where (abs(dx) > 0.5_cd*boxSideLengths)  ! apply PBC
               dx = dx - sign(boxSideLengths,dx)
            endwhere
            r2 = dot_product(dx, dx)
            if (r2.le.rcut2) then
               ! part j is a neighbor of part i
               if (half) then
                   if ( ((i.eq.MiddlePartId) .or. (j.eq.MiddlePartId)) .and. &
                        (i .lt. j) ) then
                      a = a+1
                      neighObject%neighborList(a,i) = j
                   endif
               else
                   if (i.eq.MiddlePartId .and. i.ne.j) then
                      a = a+1
                      neighObject%neighborList(a,i) = j
                   endif
               endif
            endif
        enddo
        ! part i has a-1 neighbors
        neighObject%neighborList(1,i) = a-1
    enddo

    return

end subroutine MI_OPBC_neighborlist

!-------------------------------------------------------------------------------
!
! NEIGH_PURE_neighborlist
!
!-------------------------------------------------------------------------------
subroutine NEIGH_PURE_neighborlist(half, numberOfParticles, coords, &
                                           cutoff, MiddlePartId, neighObject)
    use mod_global
    use KIM_API_F03
    use mod_neighborlist
    implicit none

    !-- Transferred variables
    logical,        intent(in)            :: half
    integer(c_int), intent(in)            :: numberOfParticles
    real(c_double), dimension(DIM,numberOfParticles), &
                  intent(in)            :: coords
    real(c_double), intent(in)            :: cutoff
    integer(c_int),         intent(in)    :: MiddlePartId
    type(neighObject_type), intent(inout) :: neighObject

    !-- Local variables
    integer(c_int) i, j, a
    real(c_double) dx(DIM)
    real(c_double) r2
    real(c_double) cutoff2

    cutoff2 = cutoff**2

    do i=1,numberOfParticles
        a = 1
        do j=1,numberOfParticles
            dx(:) = coords(:, j) - coords(:, i)
            r2 = dot_product(dx, dx)
            if (r2.le.cutoff2) then

               if (half) then
                   if ( ((i.eq.MiddlePartId) .or. (j.eq.MiddlePartId)) .and. &
                        (i .lt. j) ) then
                      a = a+1
                      neighObject%neighborList(a,i) = j
                   endif
               else
                   if (i.eq.MiddlePartId .and. i.ne.j) then
                      a = a+1
                      neighObject%neighborList(a,i) = j
                   endif
               endif

            endif
        enddo
        ! part i has a-1 neighbors
        neighObject%neighborList(1,i) = a-1
    enddo

    return

end subroutine NEIGH_PURE_neighborlist

!-------------------------------------------------------------------------------
!
! NEIGH_RVEC_neighborlist
!
!-------------------------------------------------------------------------------
subroutine NEIGH_RVEC_neighborlist(half, numberOfParticles, coords, &
                                           cutoff, MiddlePartId, neighObject)
    use mod_global
    use KIM_API_F03
    use mod_neighborlist
    implicit none

    !-- Transferred variables
    logical,        intent(in)            :: half
    integer(c_int), intent(in)            :: numberOfParticles
    real(c_double), dimension(DIM,numberOfParticles), &
                  intent(in)            :: coords
    real(c_double), intent(in)            :: cutoff
    integer(c_int),         intent(in)    :: MiddlePartId
    type(neighObject_type), intent(inout) :: neighObject

    !-- Local variables
    integer(c_int) i, j, a
    real(c_double) dx(DIM)
    real(c_double) r2
    real(c_double) cutoff2

    cutoff2 = cutoff**2

    do i=1,numberOfParticles
        a = 1
        do j=1,numberOfParticles
            dx(:) = coords(:, j) - coords(:, i)
            r2 = dot_product(dx, dx)
            if (r2.le.cutoff2) then
                if (half) then
                    if ( ((i.eq.MiddlePartId) .or. (j.eq.MiddlePartId)) .and. &
                        (i .lt. j) ) then
                        a = a+1
                        neighObject%neighborList(a,i) = j
                        neighObject%RijList(:,a-1,i) = dx
                    endif
                else
                    if (i.eq.MiddlePartId .and. i.ne.j) then
                      a = a+1
                      neighObject%neighborList(a,i) = j
                      neighObject%RijList(:,a-1,i) = dx
                    endif
                endif
            endif
         enddo
         ! part i has a-1 neighbors
         neighObject%neighborList(1,i) = a-1
    enddo

    return

end subroutine NEIGH_RVEC_neighborlist

!-------------------------------------------------------------------------------
!
! create_FCC_configuration subroutine
!
!  creates a cubic configuration of FCC particles with lattice spacing
!  `LatticeSpacing' and `nCellsPerSide' cells along each direction.
!
!  With periodic==.true. this will result in a total number of particles equal
!  to 4*(nCellsPerSide)**3 + 6*(nCellsPerSide)**2 + 3*(nCellsPerSide) + 1
!
!  With periodic==.false. this will result in a total number of particles equal
!  to 4*(nCellsPerSide)**3
!
!  Returns the Id of the particle situated in the middle of the configuration
!  (this particle will have the most neighbors.)
!
!-------------------------------------------------------------------------------
subroutine create_FCC_configuration(LatticeSpacing, nCellsPerSide, periodic, &
                                coords, MiddlePartId)
    use mod_global
    use KIM_API_F03
    implicit none

    !-- Transferred variables
    real(c_double), intent(in)  :: LatticeSpacing
    integer(c_int), intent(in)  :: nCellsPerSide
    logical,        intent(in)  :: periodic
    real(c_double), intent(out) :: coords(DIM,*)
    integer(c_int), intent(out) :: MiddlePartId
    !
    ! cluster setup variables
    !
    real(c_double) FCCshifts(DIM,4)
    real(c_double) latVec(DIM)
    integer(c_int) a, i, j, k, m

    ! Create a cubic FCC cluster
    !
    FCCshifts(1,1) = 0.0_cd
    FCCshifts(2,1) = 0.0_cd
    FCCshifts(3,1) = 0.0_cd
    FCCshifts(1,2) = 0.5_cd*LatticeSpacing
    FCCshifts(2,2) = 0.5_cd*LatticeSpacing
    FCCshifts(3,2) = 0.0_cd
    FCCshifts(1,3) = 0.5_cd*LatticeSpacing
    FCCshifts(2,3) = 0.0_cd
    FCCshifts(3,3) = 0.5_cd*LatticeSpacing
    FCCshifts(1,4) = 0.0_cd
    FCCshifts(2,4) = 0.5_cd*LatticeSpacing
    FCCshifts(3,4) = 0.5_cd*LatticeSpacing

    MiddlePartID = 1 ! Always put middle particle as #1
    a = 1            ! leave space for middle particle as particle #1
    do i=1,nCellsPerSide
     latVec(1) = (i-1)*LatticeSpacing
     do j=1,nCellsPerSide
        latVec(2) = (j-1)*LatticeSpacing
        do k=1,nCellsPerSide
           latVec(3) = (k-1)*LatticeSpacing
           do m=1,4
              a = a+1
              coords(:,a) = latVec + FCCshifts(:,m)
              if ((i.eq.nCellsPerside/2+1).and.(j.eq.nCellsPerSide/2+1) .and. &
                   (k.eq.nCellsPerSide/2+1) .and. (m.eq.1)) then
                 coords(:,1) = latVec + FCCshifts(:,m) ! put middle particle as #1
                 a = a - 1
              endif
           enddo
        enddo
        if (.not. periodic) then
            ! Add in the remaining three faces
            ! pos-x face
            latVec(1) = nCellsPerSide*LatticeSpacing
            latVec(2) = (i-1)*LatticeSpacing
            latVec(3) = (j-1)*LatticeSpacing
            a = a+1; coords(:,a) = latVec
            a = a+1; coords(:,a) = latVec + FCCshifts(:,4)
            ! pos-y face
            latVec(1) = (i-1)*LatticeSpacing
            latVec(2) = nCellsPerSide*LatticeSpacing
            latVec(3) = (j-1)*LatticeSpacing
            a = a+1; coords(:,a) = latVec
            a = a+1; coords(:,a) = latVec + FCCshifts(:,3)
            ! pos-z face
            latVec(1) = (i-1)*LatticeSpacing
            latVec(2) = (j-1)*LatticeSpacing
            latVec(3) = nCellsPerSide*LatticeSpacing
            a = a+1; coords(:,a) = latVec
            a = a+1; coords(:,a) = latVec + FCCshifts(:,2)
         endif
     enddo
     if (.not. periodic) then
         ! Add in the remaining three edges
         latVec(1) = (i-1)*LatticeSpacing
         latVec(2) = nCellsPerSide*LatticeSpacing
         latVec(3) = nCellsPerSide*LatticeSpacing
         a = a+1; coords(:,a) = latVec
         latVec(1) = nCellsPerSide*LatticeSpacing
         latVec(2) = (i-1)*LatticeSpacing
         latVec(3) = nCellsPerSide*LatticeSpacing
         a = a+1; coords(:,a) = latVec
         latVec(1) = nCellsPerSide*LatticeSpacing
         latVec(2) = nCellsPerSide*LatticeSpacing
         latVec(3) = (i-1)*LatticeSpacing
         a = a+1; coords(:,a) = latVec
      endif
    enddo
    if (.not. periodic) then
      ! Add in the remaining corner
      a = a+1; coords(:,a) = nCellsPerSide*LatticeSpacing
    endif

    return

end subroutine create_FCC_configuration

!-------------------------------------------------------------------------------
!
! create_BCC_configuration subroutine
!
!  creates a cubic configuration of BCC particles with lattice spacing
!  `LatticeSpacing' and `nCellsPerSide' cells along each direction.
!
!  With periodic==.true. this will result in a total number of particles equal
!  to 2*(nCellsPerSide)**3 + 3*(nCellsPerSide)**2 + 3*(nCellsPerSide) + 1
!
!  With periodic==.false. this will result in a total number of particles equal
!  to 2*(nCellsPerSide)**3
!
!  Returns the Id of the particle situated in the middle of the configuration
!  (this particle will have the most neighbors.)
!
!-------------------------------------------------------------------------------
subroutine create_BCC_configuration(LatticeSpacing, nCellsPerSide, periodic, &
                                coords, MiddlePartId)
    use mod_global
    use KIM_API_F03
    implicit none

    !-- Transferred variables
    real(c_double), intent(in)  :: LatticeSpacing
    integer(c_int), intent(in)  :: nCellsPerSide
    logical,        intent(in)  :: periodic
    real(c_double), intent(out) :: coords(DIM,*)
    integer(c_int), intent(out) :: MiddlePartId
    !
    ! cluster setup variables
    !
    real(c_double) FCCshifts(DIM,2)
    real(c_double) latVec(DIM)
    integer(c_int) a, i, j, k, m

    ! Create a cubic BCC cluster
    !
    FCCshifts(1,1) = 0.0_cd
    FCCshifts(2,1) = 0.0_cd
    FCCshifts(3,1) = 0.0_cd
    FCCshifts(1,2) = 0.5_cd*LatticeSpacing
    FCCshifts(2,2) = 0.5_cd*LatticeSpacing
    FCCshifts(3,2) = 0.5_cd*LatticeSpacing

    MiddlePartID = 1 ! Always put middle particle as #1
    a = 1            ! leave space for middle particle as particle #1
    do i=1,nCellsPerSide
     latVec(1) = (i-1)*LatticeSpacing
     do j=1,nCellsPerSide
        latVec(2) = (j-1)*LatticeSpacing
        do k=1,nCellsPerSide
           latVec(3) = (k-1)*LatticeSpacing
           do m=1,2
              a = a+1
              coords(:,a) = latVec + FCCshifts(:,m)
              if ((i.eq.nCellsPerside/2+1).and.(j.eq.nCellsPerSide/2+1) .and. &
                   (k.eq.nCellsPerSide/2+1) .and. (m.eq.1)) then
                 coords(:,1) = latVec + FCCshifts(:,m) ! put middle particle as #1
                 a = a - 1
              endif
           enddo
        enddo
        if (.not. periodic) then
            ! Add in the remaining three faces
            ! pos-x face
            latVec(1) = nCellsPerSide*LatticeSpacing
            latVec(2) = (i-1)*LatticeSpacing
            latVec(3) = (j-1)*LatticeSpacing
            a = a+1; coords(:,a) = latVec
            ! pos-y face
            latVec(1) = (i-1)*LatticeSpacing
            latVec(2) = nCellsPerSide*LatticeSpacing
            latVec(3) = (j-1)*LatticeSpacing
            a = a+1; coords(:,a) = latVec
            ! pos-z face
            latVec(1) = (i-1)*LatticeSpacing
            latVec(2) = (j-1)*LatticeSpacing
            latVec(3) = nCellsPerSide*LatticeSpacing
            a = a+1; coords(:,a) = latVec
         endif
     enddo
     if (.not. periodic) then
         ! Add in the remaining three edges
         latVec(1) = (i-1)*LatticeSpacing
         latVec(2) = nCellsPerSide*LatticeSpacing
         latVec(3) = nCellsPerSide*LatticeSpacing
         a = a+1; coords(:,a) = latVec
         latVec(1) = nCellsPerSide*LatticeSpacing
         latVec(2) = (i-1)*LatticeSpacing
         latVec(3) = nCellsPerSide*LatticeSpacing
         a = a+1; coords(:,a) = latVec
         latVec(1) = nCellsPerSide*LatticeSpacing
         latVec(2) = nCellsPerSide*LatticeSpacing
         latVec(3) = (i-1)*LatticeSpacing
         a = a+1; coords(:,a) = latVec
      endif
    enddo
    if (.not. periodic) then
      ! Add in the remaining corner
      a = a+1; coords(:,a) = nCellsPerSide*LatticeSpacing
    endif

    return

end subroutine create_BCC_configuration

subroutine create_diamond_configuration(LatticeSpacing, nCellsPerSide, &
                                        coords)
    use mod_global
    use KIM_API_F03
    implicit none

    !-- Transferred variables
    real(c_double), intent(in)  :: LatticeSpacing
    integer(c_int), intent(in)  :: nCellsPerSide
    real(c_double), intent(out) :: coords(DIM,*)
    !
    ! cluster setup variables
    !
    real(c_double) FCCshifts(DIM,5)
    real(c_double) latVec(DIM)
    integer(c_int) a, i, j, k, m

    ! Create a cubic FCC cluster
    !
    FCCshifts(1,1) = 0.0_cd
    FCCshifts(2,1) = 0.0_cd
    FCCshifts(3,1) = 0.0_cd
    FCCshifts(1,2) = 0.25_cd*LatticeSpacing
    FCCshifts(2,2) = 0.25_cd*LatticeSpacing
    FCCshifts(3,2) = 0.25_cd*LatticeSpacing
    FCCshifts(1,3) = 0.5_cd*LatticeSpacing
    FCCshifts(2,3) = 0.5_cd*LatticeSpacing
    FCCshifts(3,3) = 0.0_cd
    FCCshifts(1,4) = 0.5_cd*LatticeSpacing
    FCCshifts(2,4) = 0.0_cd
    FCCshifts(3,4) = 0.5_cd*LatticeSpacing
    FCCshifts(1,5) = 0.0_cd
    FCCshifts(2,5) = 0.5_cd*LatticeSpacing
    FCCshifts(3,5) = 0.5_cd*LatticeSpacing

    a = 2            ! leave space for middle particle as particle #1
    do i=1,nCellsPerSide
     latVec(1) = (i-1)*LatticeSpacing
     do j=1,nCellsPerSide
        latVec(2) = (j-1)*LatticeSpacing
        do k=1,nCellsPerSide
           latVec(3) = (k-1)*LatticeSpacing
           do m=1,5
              a = a+1
              coords(:,a) = latVec + FCCshifts(:,m)
              if ((i.eq.nCellsPerside/2+1).and.(j.eq.nCellsPerSide/2+1) .and. &
                   (k.eq.nCellsPerSide/2+1) ) then
                if (m .eq. 1) then
                     coords(:,1) = latVec + FCCshifts(:,m) ! put middle particle as #1
                elseif (m .eq. 2) then
                    coords(:,2) = latVec + FCCshifts(:,m)
                endif
                 a = a - 1
              endif
           enddo
        enddo
     enddo
 enddo
    return

end subroutine create_diamond_configuration

subroutine numerical_check_Cijkl(pkim,nbc)
    use mod_global
    use mod_coordinates
    use mod_process, only : Cijkl_numerical
    use mod_neighborlist, only : neighObject_type
    use KIM_API_F03
    implicit none
    ! Passed variables
    type(c_ptr),   intent(in) :: pkim
    integer(c_int),intent(in) :: nbc

    ! Local variables
    real(c_double), allocatable :: coords_old(:,:) 
    real(c_double)              :: energy_plus,energy_neg, &
                                   eps,gradU(DIM,DIM)
    integer(c_int) :: partIndex,xc,yc,p,q,r,s,t,u,ier,idum

    integer(c_int), pointer :: N;             type(c_ptr) :: pnParts
    real(c_double), pointer :: energy;        type(c_ptr) :: penergy
    real(c_double), pointer :: cutoff;        type(c_ptr) :: pcutoff
    type(neighObject_type), pointer :: neighObject;   type(c_ptr) :: pneighObject

    ! Unpack data from KIM object
    !
    call kim_api_getm_data(pkim, ier, &
       "numberOfParticles",     pnParts,          1,              &
       "neighObject",           pneighObject,     1,              &
       "cutoff",                pcutoff,          1,              &
       "energy",                penergy,          1)

    if (ier.lt.KIM_STATUS_OK) then
        idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                     "kim_api_getm_data", ier)
        stop
    endif

    call c_f_pointer(pnParts, N)
    call c_f_pointer(penergy, energy)
    call c_f_pointer(pcutoff, cutoff)
    call c_f_pointer(pneighObject, neighObject)

    allocate(coords_old(DIM,N))

    ! Initialize
    coords_old      = coords
    eps             = 1.e-4
    Cijkl_numerical = 0._cd

    do u=1,3
        energy_plus = 0._cd
        energy_neg  = 0._cd
        coords = coords_old

        select case(u)
        case(1)
            p=1; q=1; r=1; s=1
        case(2)
            p=1; q=2; r=1; s=2
        case(3)
            p=1; q=1; r=2; s=2
        end select

        do t=1,4
            coords = coords_old
            gradU  = 0._cd
            select case(t)
            case(1)
                if (p .eq. q) then
                    gradU(p,q) = eps
                else
                    gradU(p,q) = eps/2._cd
                    gradU(q,p) = eps/2._cd
                endif
                if (r .eq. s) then
                    gradU(r,s) = eps
                else
                    gradU(r,s) = eps/2._cd
                    gradU(s,r) = eps/2._cd
                endif
            case(2)
                if (p .eq. q) then
                    gradU(p,q) = eps
                else
                    gradU(p,q) = eps/2._cd
                    gradU(q,p) = eps/2._cd
                endif
                if (r .eq. s) then
                    gradU(r,s) = -eps
                else
                    gradU(r,s) = -eps/2._cd
                    gradU(s,r) = -eps/2._cd
                endif
                if (p .eq. r .and. q .eq. s) then
                    gradU(p,q) = 0._cd
                    gradU(q,p) = 0._cd
                endif
            case(3)
                if (p .eq. q) then
                    gradU(p,q) = -eps
                else
                    gradU(p,q) = -eps/2._cd
                    gradU(q,p) = -eps/2._Cd
                endif
                if (r .eq. s) then
                    gradU(r,s) = eps
                else
                    gradU(r,s) = eps/2._cd
                    gradU(s,r) = eps/2._cd
                endif
                if (p .eq. r .and. q .eq. s) then 
                    gradU(p,q) = 0._cd
                    gradU(q,p) = 0._cd
                endif
            case(4)
                if (p .eq. q) then
                    gradU(p,q) = -eps
                else
                    gradU(p,q) = -eps/2._cd
                    gradU(q,p) = -eps/2._cd
                endif
                if (r .eq. s) then
                    gradU(r,s) = -eps
                else
                    gradU(r,s) = -eps/2._cd
                    gradU(s,r) = -eps/2._cd
                endif
            endselect

            ! Deform
            do partIndex = 1,N
                do xc=1,3
                    do yc=1,3
                        coords(xc,partIndex) = coords(xc,partIndex) + gradU(xc,yc) * coords_old(yc,partIndex)
                    enddo
                enddo
            enddo

            if (nbc.eq.0) then
                call NEIGH_RVEC_neighborlist(.true., N, coords, cutoff, &
                                                  MiddlePartId, neighObject)
            elseif (nbc.eq.2) then
                call NEIGH_RVEC_neighborlist(.false., N, coords, cutoff, &
                                                  MiddlePartId, neighObject)
            endif
            
            ! Now compute energy
            ier = kim_api_model_compute(pkim)
            if (ier.lt.KIM_STATUS_OK) then
               idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                           "kim_api_model_compute", ier)
               stop
            endif

            if (t .eq. 1 .or. t .eq. 4) then
                energy_plus = energy_plus + energy
            elseif (t .eq. 2 .or. t .eq. 3) then
                energy_neg = energy_neg + energy
            endif
        enddo


        if ( p .eq. r .and. q .eq. s) then
            Cijkl_numerical(p,q,r,s) = (energy_plus - energy_neg)/(eps**2)
        else
            Cijkl_numerical(p,q,r,s) = (energy_plus - energy_neg)/(4._cd * eps**2)
        endif
    enddo

coords = coords_old
deallocate(coords_old)

end subroutine numerical_check_Cijkl

subroutine numerical_check_Dijmkln(pkim,nbc)
    use mod_global
    use KIM_API_F03
    use mod_coordinates
    use mod_neighborlist, only : neighObject_type
    use mod_process, only : Dijmkln_numerical
    implicit none
    
    ! Passed variables
    type(c_ptr),   intent(in)  :: pkim
    integer(c_int),intent(in) :: nbc

    ! Local variables
    real(c_double), allocatable            :: coords_old(:,:)
    integer(c_int)                         :: p,q,r,s,i,j,m,k,l,n,ier,idum
    real(c_double)                         :: energy_plus,energy_neg,eps
    real(c_double), dimension(DIM)         :: shift
    real(c_double), dimension(DIM,DIM,DIM) :: Cijk, Eijk

    integer(c_int), pointer :: nParts;       type(c_ptr) :: pnParts
    real(c_double), pointer :: energy;  type(c_ptr) :: penergy
    real(c_double), pointer :: cutoff;        type(c_ptr) :: pcutoff
    type(neighObject_type), pointer :: neighObject;   type(c_ptr) :: pneighObject


    ! Unpack data from KIM object
    !
    call kim_api_getm_data(pkim, ier, &
       "numberOfParticles",     pnParts,          1,   &
       "neighObject",           pneighObject,     1,   &
       "cutoff",                pcutoff,          1,   &
       "energy",                penergy,          1)

    if (ier.lt.KIM_STATUS_OK) then
        idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                     "kim_api_getm_data", ier)
        stop
    endif

    call c_f_pointer(pnParts, nParts)
    call c_f_pointer(penergy, energy)
    call c_f_pointer(pcutoff, cutoff)
    call c_f_pointer(pneighObject, neighObject)

    allocate(coords_old(DIM,nParts))
    coords_old      = coords
    eps             = 1.e-3
    Dijmkln_numerical = 0._cd

    shift(:) = coords_old(:,MiddlePartId)
    do p=1,nParts
        coords_old(:,p) = coords_old(:,p) - shift(:)
    enddo


    do q=1,11
        energy_plus = 0._cd
        energy_neg  = 0._cd

        select case(q)
        case(1)
            i = 1; j = 1; m = 1
            k = 1; l = 1; n = 1
        case(2)
            i = 1; j = 1; m = 1
            k = 1; l = 2; n = 2
        case(3)
            i = 1; j = 1; m = 1
            k = 2; l = 2; n = 1
        case(4)
            i = 1; j = 1; m = 2
            k = 1; l = 1; n = 2
        case(5)
            i = 1; j = 1; m = 2
            k = 1; l = 2; n = 1
        case(6)
            i = 1; j = 1; m = 2
            k = 2; l = 3; n = 3
        case(7)
            i = 1; j = 1; m = 2
            k = 3; l = 3; n = 2
        case(8)
            i = 1; j = 2; m = 1
            k = 1; l = 2; n = 1
        case(9)
            i = 1; j = 2; m = 1
            k = 2; l = 3; n = 3
        case(10)
            i = 1; j = 2; m = 3
            k = 1; l = 2; n = 3
        case(11)
            i = 1; j = 2; m = 3
            k = 1; l = 3; n = 2
        end select

        do p=1,4
            coords = coords_old
            Eijk = 0._cd
            Cijk = 0._cd

            select case(p)
                case(1)
                    if (i .eq. j) then
                        Eijk(i,j,m) = eps
                    else
                        Eijk(i,j,m) = eps/2._cd
                        Eijk(j,i,m) = eps/2._cd
                    endif
                    if (k .eq. l) then
                        Eijk(k,l,n) = eps
                    else
                        Eijk(k,l,n) = eps/2._cd
                        Eijk(l,k,n) = eps/2._cd
                    endif
                case(2)
                    if (i .eq. j) then
                        Eijk(i,j,m) = eps
                    else
                        Eijk(i,j,m) = eps/2._cd
                        Eijk(j,i,m) = eps/2._cd
                    endif
                    if (k .eq. l) then
                        Eijk(k,l,n) = -eps
                    else
                        Eijk(k,l,n) = -eps/2._cd
                        Eijk(l,k,n) = -eps/2._cd
                    endif
                    if ( i .eq. k .and. j .eq. l .and. m .eq. n) then
                        Eijk(i,j,m) = 0._cd
                        Eijk(j,i,m) = 0._cd
                    endif
                case(3)
                    if (i .eq. j) then
                        Eijk(i,j,m) = -eps
                    else
                        Eijk(i,j,m) = -eps/2._cd
                        Eijk(j,i,m) = -eps/2._cd
                    endif
                    if (k .eq. l) then
                        Eijk(k,l,n) = eps
                    else
                        Eijk(k,l,n) = eps/2._cd
                        Eijk(l,k,n) = eps/2._cd
                    endif
                    if ( i .eq. k .and. j .eq. l .and. m .eq. n) then
                        Eijk(i,j,m) = 0._cd
                        Eijk(j,i,m) = 0._cd
                    endif
                case(4)
                    if (i .eq. j) then
                        Eijk(i,j,m) = -eps
                    else
                        Eijk(i,j,m) = -eps/2._cd
                        Eijk(j,i,m) = -eps/2._cd
                    endif
                    if (k .eq. l) then
                        Eijk(k,l,n) = -eps
                    else
                        Eijk(k,l,n) = -eps/2._cd
                        Eijk(l,k,n) = -eps/2._cd
                    endif
            end select

            call Calculate_Cijk(Eijk,Cijk)
        
            do r=1,nParts
                do s=1,3
                   coords(s,r)         = coords(s,r) + 0.5d0*Cijk(s,1,1) * coords_old(1,r) * coords_old(1,r) +   &
                                                       0.5d0*Cijk(s,1,2) * coords_old(1,r) * coords_old(2,r) +   &
                                                       0.5d0*Cijk(s,1,3) * coords_old(1,r) * coords_old(3,r) +   &
                                                       0.5d0*Cijk(s,2,1) * coords_old(2,r) * coords_old(1,r) +   &
                                                       0.5d0*Cijk(s,2,2) * coords_old(2,r) * coords_old(2,r) +   &
                                                       0.5d0*Cijk(s,2,3) * coords_old(2,r) * coords_old(3,r) +   &
                                                       0.5d0*Cijk(s,3,1) * coords_old(3,r) * coords_old(1,r) +   &
                                                       0.5d0*Cijk(s,3,2) * coords_old(3,r) * coords_old(2,r) +   &
                                                       0.5d0*Cijk(s,3,3) * coords_old(3,r) * coords_old(3,r)
                enddo
            enddo

            if (nbc.eq.0) then
                call NEIGH_RVEC_neighborlist(.true., nParts, coords, cutoff, &
                                                  MiddlePartId, neighObject)
            elseif (nbc.eq.2) then
                call NEIGH_RVEC_neighborlist(.false., nParts, coords, cutoff, &
                                                  MiddlePartId, neighObject)
            endif
            
            ier = kim_api_model_compute(pkim)
            if (ier.lt.KIM_STATUS_OK) then
               idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                           "kim_api_model_compute", ier)
               stop
            endif

            if (p .eq. 1 .or. p .eq. 4) then
                energy_plus = energy_plus + energy
            elseif (p .eq. 2 .or. p .eq. 3) then
                energy_neg = energy_neg + energy
            endif

        enddo

        if ( i .eq. k .and. j .eq. l .and. m .eq. n) then
            Dijmkln_numerical(i,j,m,k,l,n) = ((energy_plus - energy_neg)/(eps**2))
        else
            Dijmkln_numerical(i,j,m,k,l,n) = ((energy_plus - energy_neg)/(4._cd * eps**2))
        endif
    enddo

    coords = coords_old

    deallocate(coords_old)


end subroutine numerical_check_Dijmkln

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

subroutine Calculate_Cijk(Eijk,Cijk)
    use mod_global
    implicit none

    ! Passed variables
    double precision,intent(in) :: Eijk(DIM,DIM,DIM)
    ! Local variables
    integer :: i,j,k
    double precision,intent(out) :: Cijk(DIM,DIM,DIM)

    Cijk = 0.d0
    do i=1,DIM
        do j=1,DIM
            do k=1,DIM
                Cijk(i,j,k) = Eijk(i,j,k) + Eijk(i,k,j) - Eijk(k,j,i) 
            enddo
        enddo
    enddo
    return
end subroutine Calculate_Cijk
