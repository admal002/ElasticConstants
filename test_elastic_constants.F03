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
! Copyright (c) 2013--2014, Regents of the University of Minnesota.
! All rights reserved.
!
! Contributors:
!    Ellad B. Tadmor
!    Ryan S. Elliott
!    Stephen M. Whalen
!

!
! Release: This file is part of the kim-api.git package.
!


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
! module mod coordinates:
!     Contains coordinates and the middle particle ID of the sysetem
!
module coordinates
    use mod_global
    implicit none
    real(c_double), pointer :: coords(:,:); type(c_ptr) :: pcoor
    integer(c_int)          :: MiddlePartId
end module coordinates

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

    public setup_neighborlist_KIM_access, get_neigh

    type neighObject_type
        integer(c_int), pointer :: neighborList(:,:)
        real(c_double), pointer :: RijList(:,:,:)
    end type neighObject_type

contains

!-------------------------------------------------------------------------------
!
! setup_neighborlist_KIM_access :
!
!    Store necessary pointers in KIM API object to access the neighbor list
!    data and methods.
!
!-------------------------------------------------------------------------------
subroutine setup_neighborlist_KIM_access(pkim, neighObject)
  implicit none

  !-- Transferred variables
  type(c_ptr),                    intent(in) :: pkim
  type(neighObject_type), target, intent(in) :: neighObject

  !-- Local variables
  integer(c_int), parameter :: SizeOne = 1
  integer(c_int)            :: ier, idum

  ! store location of neighObject variable
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

  return

end subroutine setup_neighborlist_KIM_access

!-------------------------------------------------------------------------------
!
! get_neigh neighbor list access function
!
! This function implements Locator and Iterator mode
!
!-------------------------------------------------------------------------------
function get_neigh(pkim,mode,request,part,numnei,pnei1part, &
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
  if (associated(neighObject%RijList)) then
    pRij = c_loc(neighObject%RijList(1,1,part))
  else
    pRij = c_null_ptr
  endif

  get_neigh = KIM_STATUS_OK
  return
end function get_neigh

end module mod_neighborlist

!-------------------------------------------------------------------------------
! module mod_process:
!     Module contains type and routines related to the processing of the
!     derivatives of the potential energy
!
module mod_process
use mod_global
implicit none

real(c_double) :: Cijkl(DIM,DIM,DIM,DIM),Dijmkln(DIM,DIM,DIM,DIM,DIM,DIM)
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
    use coordinates
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
    use coordinates
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
                               operatorT(i,j,m,k,l,n) = operatorT(i,j,m,k,l,n) +    &
                                                        tensor4(l,m,j,n)       +    &
                                                        tensor4(j,n,l,m)       +    &
                                                        tensor4(j,l,m,n)       +    &
                                                        tensor4(m,n,j,l)
                           endif

                           if (i .eq. l) then
                               operatorT(i,j,m,k,l,n) = operatorT(i,j,m,k,l,n) +    &
                                                        tensor4(j,k,m,n)       +    &
                                                        tensor4(j,n,k,m)       +    &
                                                        tensor4(k,m,j,n)       +    &
                                                        tensor4(m,n,j,k)
                           endif

                           if (i .eq. n) then
                               operatorT(i,j,m,k,l,n) = operatorT(i,j,m,k,l,n) -    &
                                                        tensor4(j,k,l,m)       -    &
                                                        tensor4(j,l,k,m)       -    &
                                                        tensor4(k,m,j,l)       -    &
                                                        tensor4(l,m,j,k)
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

end module mod_process

!*******************************************************************************
!**
!**  PROGRAM test_elastic_constants
!**
!**  KIM compliant program to compute Mindlin's First Strain Gradient elastic
!**  constants
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
program test_elastic_constants
    use mod_global
    use KIM_API_F03
    use coordinates
    use mod_neighborlist
    use mod_process, only : Cijkl,Dijmkln,setup_process_KIM_access
    implicit none

    !real(c_double), parameter :: FCCspacing     = 4.05_cd       ! For Al in angstroms
    !real(c_double), parameter :: FCCspacing     = 4.3265037_cd ! For Si in  angstroms
    real(c_double), parameter :: FCCspacing     = 5.29216_cd   ! For Ar in angstroms
    !real(c_double), parameter :: FCCspacing     = 3.6045828843981988    ! For custom lj

    integer(c_int), parameter :: nCellsPerSide  = 8
    integer(c_int), parameter :: ASpecies       = 1
    integer(c_int), parameter :: N              = 4*(nCellsPerSide)**3
    integer(c_int), parameter :: SizeOne        = 1

    !
    ! neighbor list
    !
    type(neighObject_type), target :: neighObject
    real(c_double), parameter      :: cutpad = 0.75_cd ! in angstroms

    !
    ! KIM variables
    !
    character(len=KIM_KEY_STRING_LENGTH) :: testname    = "test_elastic_constants"
    character(len=KIM_KEY_STRING_LENGTH) :: testkimfile = "descriptor.kim"
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
    real(c_double), pointer :: forces(:,:);         type(c_ptr) :: pforces
    real(c_double), pointer :: boxSideLengths(:);   type(c_ptr) :: pboxSideLengths
    character(len=2)        :: species

    Cijkl = 0._cd
    Dijmkln = 0._cd

    ! Get KIM Model name to use
    print '("Please enter a valid KIM model name: ")'

    !modelname='EAM_Dynamo_Zope_Mishin_Al__MO_664470114311_001'
    !species="Al"

    !modelname='Three_Body_Stillinger_Weber_Hauch_Brittle_Si__MO_119167353542_001'
    !species="Si"

    modelname='Pair_LJ_Smooth_Bernardes_Ar__MO_764178710049_000'
    species="Ar"

    !modelname='Pair_LJ_Custom'
    !species="Ar"

    ! Initialize the KIM object
    ier = kim_api_file_init(pkim, testkimfile, modelname)
    if (ier.lt.KIM_STATUS_OK) then
        idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                 "kim_api_file_init", ier)
        stop
    endif

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
       "energy",                      penergy,          1,                                   &
       "forces",                      pforces,          1)
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

    write(*,*) "Cutoff = ",cutoff

    if ((nbc.eq.4).or.(nbc.eq.5)) call c_f_pointer(pboxSideLengths, &
                                                 boxSideLengths, [DIM])
    call c_f_pointer(penergy, energy)
    call c_f_pointer(pforces, forces, [DIM,N])

    ! Set values
    numberOfParticles = N
    if ((nbc.eq.0).or.(nbc.eq.1).or.(nbc.eq.4)) numContrib = 1
    numberOfSpecies = ASpecies
    particleSpecies(:)    = kim_api_get_species_code(pkim, species, ier)
    if (ier.lt.KIM_STATUS_OK) then
        idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                 "kim_api_get_species_code", ier)
        stop
    endif

    ! set up the cluster part positions
    call create_FCC_configuration(FCCspacing, nCellsPerSide, .true., coords, &
                                MiddlePartId)
    ! set boxSideLengths large enough to make the system isolated
    if (nbc.eq.4 .or. nbc.eq.5) boxSideLengths(:) = nCellsPerSide*FCCspacing

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
    print '("This is Test          : ",A)', trim(testname)
    print '("Results for KIM Model : ",A)', trim(modelname)
    print '("Using NBC: ",A)', trim(NBC_Method)
    print '("Energy        = ",ES25.15)',energy

    Cijkl = Cijkl/((FCCspacing**3)/4._cd)
    Dijmkln = Dijmkln/((FCCspacing**3)/4._cd)

    call write_results(Cijkl,Dijmkln)

    !-------------------------------------------------------------------------------
    ! Numerical checks
    !

    ! Switch off process_dedr and process_d2edr2
    call kim_api_set_compute(pkim,'process_dEdr',KIM_COMPUTE_FALSE,ier)
    call kim_api_set_compute(pkim,'process_d2Edr2',KIM_COMPUTE_FALSE,ier)

    ! Perform the following non-uniform deformation which corresponts to zero strain
    ! and infinitesimal first strain gradient

    call nonuniform_deformation(pkim,Dijmkln,FCCspacing)
    call numerical_check(pkim,FCCSpacing)

    ! Don't forget to free and/or deallocate
    if (nbc.lt.6) deallocate(neighObject%neighborList)
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
end program test_elastic_constants

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
                   !if ( ((i.eq.MiddlePartId) .or. (j.eq.MiddlePartId)) .and. &
                   !     (i .lt. j) ) then
                   if ( j .gt. i ) then
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
!  `FCCspacing' and `nCellsPerSide' cells along each direction.
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
subroutine create_FCC_configuration(FCCspacing, nCellsPerSide, periodic, &
                                    coords, MiddlePartId)
    use mod_global
    use KIM_API_F03
    implicit none

    !-- Transferred variables
    real(c_double), intent(in)  :: FCCspacing
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
    FCCshifts(1,2) = 0.5_cd*FCCspacing
    FCCshifts(2,2) = 0.5_cd*FCCspacing
    FCCshifts(3,2) = 0.0_cd
    FCCshifts(1,3) = 0.5_cd*FCCspacing
    FCCshifts(2,3) = 0.0_cd
    FCCshifts(3,3) = 0.5_cd*FCCspacing
    FCCshifts(1,4) = 0.0_cd
    FCCshifts(2,4) = 0.5_cd*FCCspacing
    FCCshifts(3,4) = 0.5_cd*FCCspacing

    MiddlePartId = 1 ! Always put middle particle as #1
    a = 1            ! leave space for middle particle as particle #1
    do i=1,nCellsPerSide
        latVec(1) = (i-1)*FCCspacing
        do j=1,nCellsPerSide
            latVec(2) = (j-1)*FCCspacing
            do k=1,nCellsPerSide
                latVec(3) = (k-1)*FCCspacing
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
                latVec(1) = nCellsPerSide*FCCspacing
                latVec(2) = (i-1)*FCCspacing
                latVec(3) = (j-1)*FCCspacing
                a = a+1; coords(:,a) = latVec
                a = a+1; coords(:,a) = latVec + FCCshifts(:,4)
                ! pos-y face
                latVec(1) = (i-1)*FCCspacing
                latVec(2) = nCellsPerSide*FCCspacing
                latVec(3) = (j-1)*FCCspacing
                a = a+1; coords(:,a) = latVec
                a = a+1; coords(:,a) = latVec + FCCshifts(:,3)
                ! pos-z face
                latVec(1) = (i-1)*FCCspacing
                latVec(2) = (j-1)*FCCspacing
                latVec(3) = nCellsPerSide*FCCspacing
                a = a+1; coords(:,a) = latVec
                a = a+1; coords(:,a) = latVec + FCCshifts(:,2)
            endif
        enddo
        if (.not. periodic) then
            ! Add in the remaining three edges
            latVec(1) = (i-1)*FCCspacing
            latVec(2) = nCellsPerSide*FCCspacing
            latVec(3) = nCellsPerSide*FCCspacing
            a = a+1; coords(:,a) = latVec
            latVec(1) = nCellsPerSide*FCCspacing
            latVec(2) = (i-1)*FCCspacing
            latVec(3) = nCellsPerSide*FCCspacing
            a = a+1; coords(:,a) = latVec
            latVec(1) = nCellsPerSide*FCCspacing
            latVec(2) = nCellsPerSide*FCCspacing
            latVec(3) = (i-1)*FCCspacing
            a = a+1; coords(:,a) = latVec
        endif
    enddo
    if (.not. periodic) then
        ! Add in the remaining corner
        a = a+1; coords(:,a) = nCellsPerSide*FCCspacing
    endif

    return

end subroutine create_FCC_configuration

subroutine numerical_check(pkim,FCCSpacing)
    use mod_global
    use coordinates
    use KIM_API_F03
    implicit none
    ! Passed variables
    type(c_ptr),   intent(in) :: pkim
    real(c_double)            :: FCCSpacing

    ! Local variables
    real(c_double), allocatable :: coords_old(:,:) 
    real(c_double)              :: energy0,energy_plus,energy_neg, &
                                   eps,gradU(DIM,DIM),                 &
                                   Cijkl_numerical(DIM,DIM,DIM,DIM)
    integer(c_int) :: partIndex,xc,yc,p,q,r,s,t,u,ier,idum

    integer(c_int), pointer :: N;       type(c_ptr) :: pnParts
    real(c_double), pointer :: cutoff;  type(c_ptr) :: pcutoff
    real(c_double), pointer :: energy;  type(c_ptr) :: penergy

    ! Unpack data from KIM object
    !
    call kim_api_getm_data(pkim, ier, &
       "numberOfParticles",           pnParts,          1,                                   &
       "cutoff",                      pcutoff,          1,                                   &
       "energy",                      penergy,          1)

    if (ier.lt.KIM_STATUS_OK) then
        idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                     "kim_api_getm_data", ier)
        stop
    endif

    call c_f_pointer(pnParts, N)
    call c_f_pointer(pcutoff, cutoff)
    call c_f_pointer(penergy, energy)

    allocate(coords_old(DIM,N))

            ier = kim_api_model_compute(pkim)
            if (ier.lt.KIM_STATUS_OK) then
               idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                           "kim_api_model_compute", ier)
               stop
            endif

            energy0 = energy
    ! Initialize
    coords_old      = coords
    eps             = 1.e-4
    Cijkl_numerical = 0._cd

    do u=1,2
        if (u .eq. 1) then
            p=1
            q=1
        else
            p=1
            q=2
        endif

        do t=1,2
            coords = coords_old
            gradU  = 0._cd
            if (t .eq. 1) then
                gradU(p,q) = eps
            else
                gradU(p,q) = -eps
            endif

            ! Deform
            do partIndex = 1,N
                do xc=1,3
                    do yc=1,3
                        coords(xc,partIndex) = coords(xc,partIndex) + gradU(xc,yc) * coords_old(yc,partIndex)
                    enddo
                enddo
            enddo
            
            ier = kim_api_model_compute(pkim)
            if (ier.lt.KIM_STATUS_OK) then
               idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                           "kim_api_model_compute", ier)
               stop
            endif

            if (t .eq. 1) then
                energy_plus = energy
            else
                energy_neg = energy
            endif
        enddo

        ! We now have energy0, energy_plus and energy_neg
        Cijkl_numerical(p,q,p,q) = (energy_plus + energy_neg - 2._cd * energy0)/(eps**2)

        write(*,*) p,q,p,q, Cijkl_numerical(p,q,p,q)/((FCCspacing**3)/4._cd), &
                        (energy_plus-energy_neg)/((FCCspacing**3)/4._cd),N
    enddo


    p = 1
    q = 1
    r = 2
    s = 2
    do t=1,4
        coords = coords_old
        gradU  = 0._cd
        if (t .eq. 1) then
            gradU(p,q) = eps
            gradU(r,s) = eps
        elseif (t .eq. 2) then
            gradU(p,q) = eps
            gradU(r,s) = -eps
        elseif (t .eq. 3) then
            gradU(p,q) = -eps
            gradU(r,s) = eps
        else
            gradU(p,q) = -eps
            gradU(r,s) = -eps
        endif

        ! Deform
        do partIndex = 1,N
            do xc=1,3
                do yc=1,3
                    coords(xc,partIndex) = coords(xc,partIndex) + gradU(xc,yc) * coords_old(yc,partIndex)
                enddo
            enddo
        enddo
        
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

    ! We now have energy0, energy_plus and energy_neg
    Cijkl_numerical(p,q,r,s) = (energy_plus - energy_neg)/(4._cd * eps**2)
    write(*,*) p,q,r,s, Cijkl_numerical(p,q,r,s)/((FCCspacing**3)/4._cd)

    deallocate(coords_old)
end subroutine numerical_check

subroutine nonuniform_deformation(pkim,Dijmkln,FCCSpacing)
    use mod_global
    use KIM_API_F03
    use coordinates
    implicit none
    
    ! Passed variables
    type(c_ptr),   intent(in)  :: pkim
    real(c_double), intent(in) :: Dijmkln(DIM,DIM,DIM,DIM,DIM,DIM)
    real(c_double)             :: FCCSpacing

    ! Local variables
    real(c_double), allocatable            :: coords_old(:,:)
    integer(c_int)                         :: p,q,k,l,nsteps,ier,idum
    real(c_double)                         :: sgmax,EnergyDiff,energy0
    real(c_double), dimension(DIM,DIM,DIM) :: Cijk, Cijk_new, Eijk

    integer(c_int), pointer :: N;       type(c_ptr) :: pnParts
    real(c_double), pointer :: cutoff;  type(c_ptr) :: pcutoff
    real(c_double), pointer :: energy;  type(c_ptr) :: penergy

    ! Unpack data from KIM object
    !
    call kim_api_getm_data(pkim, ier, &
       "numberOfParticles",           pnParts,          1,                                   &
       "cutoff",                      pcutoff,          1,                                   &
       "energy",                      penergy,          1)

    if (ier.lt.KIM_STATUS_OK) then
        idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                     "kim_api_getm_data", ier)
        stop
    endif

    call c_f_pointer(pnParts, N)
    call c_f_pointer(pcutoff, cutoff)
    call c_f_pointer(penergy, energy)

    allocate(coords_old(DIM,N))
    coords_old      = coords
    energy0         = energy

    do p=2,N
        coords_old(:,p) = coords_old(:,p) - coords_old(:,1)
    enddo
    coords_old(:,MiddlePartId) = 0._cd

    nsteps = 100
    sgmax = 0.005d0
    do q = 1,nsteps !Deformation step
        do p=1,4
            Eijk = 0.d0
            Cijk = 0.d0
            select case(p)
                case(1)
                    Eijk(1,1,1) = 1.d0
                    EnergyDiff = 0.5d0*Dijmkln(1,1,1,1,1,1)*(q*sgmax/nsteps)**2 
                case(2)
                    Eijk(1,1,2) = 1.d0
                    EnergyDiff = 0.5d0*Dijmkln(1,1,2,1,1,2)*(q*sgmax/nsteps)**2 
                case(3)
                    Eijk(1,2,1) = 1.0d0
                    Eijk(2,1,1) = 1.0d0
                    EnergyDiff = 2.d0*Dijmkln(1,2,1,1,2,1)*(q*sgmax/nsteps)**2 
                case(4)
                    Eijk(1,2,3) = 1.0d0
                    Eijk(2,1,3) = 1.0d0
                    EnergyDiff = 2.d0*Dijmkln(1,2,3,1,2,3)*(q*sgmax/nsteps)**2 
            end select
            call Calculate_Cijk(Eijk,Cijk)


            Cijk_new = Cijk * q * sgmax / nsteps
            coords = coords_old
            do k=1,N
                do l=1,3
                   coords(l,k)         = coords(l,k) + 0.5d0*Cijk_new(l,1,1) * coords_old(1,k) * coords_old(1,k) +   &
                                                       0.5d0*Cijk_new(l,1,2) * coords_old(1,k) * coords_old(2,k) +   &
                                                       0.5d0*Cijk_new(l,1,3) * coords_old(1,k) * coords_old(3,k) +   &
                                                       0.5d0*Cijk_new(l,2,1) * coords_old(2,k) * coords_old(1,k) +   &
                                                       0.5d0*Cijk_new(l,2,2) * coords_old(2,k) * coords_old(2,k) +   &
                                                       0.5d0*Cijk_new(l,2,3) * coords_old(2,k) * coords_old(3,k) +   &
                                                       0.5d0*Cijk_new(l,3,1) * coords_old(3,k) * coords_old(1,k) +   &
                                                       0.5d0*Cijk_new(l,3,2) * coords_old(3,k) * coords_old(2,k) +   &
                                                       0.5d0*Cijk_new(l,3,3) * coords_old(3,k) * coords_old(3,k)
                enddo
            enddo

            ier = kim_api_model_compute(pkim)
            if (ier.lt.KIM_STATUS_OK) then
               idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                           "kim_api_model_compute", ier)
               stop
            endif

            write(93+p,*) q*sgmax/nsteps, (energy-energy0)* 4.d0/(FCCSpacing**3),EnergyDiff
        enddo
    enddo

    coords = coords_old

    deallocate(coords_old)

end subroutine nonuniform_deformation

subroutine write_results(Cijkl,Dijmkln)
    use mod_global
    implicit none
    ! Transferred variables
    real(c_double), intent(in) :: Cijkl(DIM,DIM,DIM,DIM), &
                                  Dijmkln(DIM,DIM,DIM,DIM,DIM,DIM)
    ! Local variables
    integer(c_int)             :: cunit,dunit
    integer,external           :: GetFreeUnit
    character(len=1)           :: cr

    cr = char(10)

    cunit = GetFreeUnit()
    open(unit=cunit,file='Cijkl.out',status='replace')
    dunit = GetFreeUnit()
    open(unit=dunit,file='Dijmkln.out',status='replace')
    
    write(cunit, *) "c11 = ",Cijkl(1,1,1,1)
    write(cunit, *) "c12 = ",Cijkl(1,1,2,2)
    write(cunit, *) "c44 = ",Cijkl(1,2,1,2)

    write(dunit,*) "D111111 = ",Dijmkln(1,1,1,1,1,1)
    write(dunit,*) "D111122 = ",Dijmkln(1,1,1,1,2,2)
    write(dunit,*) "D111221 = ",Dijmkln(1,1,1,2,2,1)
    write(dunit,*) "D112112 = ",Dijmkln(1,1,2,1,1,2)
    write(dunit,*) "D112121 = ",Dijmkln(1,1,2,1,2,1)
    write(dunit,*) "D112233 = ",Dijmkln(1,1,2,2,3,3)
    write(dunit,*) "D112332 = ",Dijmkln(1,1,2,3,3,2)
    write(dunit,*) "D121121 = ",Dijmkln(1,2,1,1,2,1)
    write(dunit,*) "D121233 = ",Dijmkln(1,2,1,2,3,3)
    write(dunit,*) "D123123 = ",Dijmkln(1,2,3,1,2,3)
    write(dunit,*) "D123132 = ",Dijmkln(1,2,3,1,3,2)

    close(cunit)
    close(dunit)

end subroutine write_results

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
