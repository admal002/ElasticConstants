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
!
! module mod_neighborlist :
!
!    Module contains type and routines related to neighbor list calculation
!
!-------------------------------------------------------------------------------

module mod_neighborlist

  use, intrinsic :: iso_c_binding
  use KIM_API_F03

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
integer(c_int) function get_neigh(pkim,mode,request,part,numnei,pnei1part, &
                                  pRij) bind(c)
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
  integer(c_int), parameter :: DIM = 3
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
use, intrinsic :: iso_c_binding
use KIM_API_F03
implicit none

integer(c_int),parameter :: DIM=3
real(c_double), allocatable :: prcoor(:,:)
real(c_double) :: Cnorm,lambda(DIM,DIM),Cijkl(DIM,DIM,DIM,DIM),Cijkl_global(DIM,DIM,DIM,DIM), Dijmkln(DIM,DIM,DIM,DIM,DIM,DIM)
contains
    subroutine setup_process_KIM_access(pkim,coords,DIM,N)
    implicit none
    ! Passed variables
    type(c_ptr), intent(in) :: pkim
    integer(c_int),intent(in) :: N,DIM
    real(c_double),intent(in) :: coords(DIM,N)

    !-- Local variables
    integer(c_int), parameter :: SizeOne = 1
    integer(c_int)            :: ier, idum,i


    allocate(prcoor(DIM,N))
    do i=1,N
        prcoor(:,i) = coords(:,i)
    enddo

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

    integer(c_int) function process_dedr(pkim,deidr,r,pdx,i,j)bind(c)
    implicit none

    ! Passed variables
    type(c_ptr),   intent(in) :: pkim,pdx
    real(c_double),intent(in) :: deidr,r
    integer(c_int),intent(in) :: i,j

    ! Local variables
    real(c_double), pointer :: dx(:)
    integer(c_int) :: a,b,c,d,e,f
    real(c_double) :: factor

    call c_f_pointer(pdx,dx,[3])

    factor = 0.5*deidr/r
    Cijkl_global = Cijkl_global - factor * tensor_prod(dx,dx,dx,dx)/(r**2)
    if (i .eq. 1 .or. j .eq. 1) then
        factor = 0.5*deidr/r
        Cijkl = Cijkl - factor * tensor_prod(dx,dx,dx,dx)/(r**2)
        Dijmkln = Dijmkln + 0.25*factor*tensor_prod_6(dx,dx,dx,dx,dx,dx)/(r**2)
       do f=1,DIM
           do e=1,DIM
               do d=1,DIM
                   do c=1,DIM
                       do b=1,DIM
                           do a=1,DIM
                         if (a .eq. d) then
                             Dijmkln(a,b,c,d,e,f) = Dijmkln(a,b,c,d,e,f) - 0.25*factor * dx(b)*dx(c)*dx(e)*dx(f)
                         endif
                         if (a .eq. e) then
                             Dijmkln(a,b,c,d,e,f) = Dijmkln(a,b,c,d,e,f) - 0.25*factor * dx(b)*dx(c)*dx(d)*dx(f)
                         endif
                         if (a .eq. f) then
                             Dijmkln(a,b,c,d,e,f) = Dijmkln(a,b,c,d,e,f) + 0.25*factor * dx(b)*dx(c)*dx(d)*dx(e)
                         endif
                         if (b .eq. d) then
                             Dijmkln(a,b,c,d,e,f) = Dijmkln(a,b,c,d,e,f) - 0.25*factor * dx(a)*dx(c)*dx(e)*dx(f)
                         endif
                         if (b .eq. e) then
                             Dijmkln(a,b,c,d,e,f) = Dijmkln(a,b,c,d,e,f) - 0.25*factor * dx(a)*dx(c)*dx(d)*dx(f)
                         endif
                         if (b .eq. f) then
                             Dijmkln(a,b,c,d,e,f) = Dijmkln(a,b,c,d,e,f) + 0.25*factor * dx(a)*dx(c)*dx(d)*dx(e)
                         endif
                         if (c .eq. d) then
                             Dijmkln(a,b,c,d,e,f) = Dijmkln(a,b,c,d,e,f) + 0.25*factor * dx(a)*dx(b)*dx(e)*dx(f)
                         endif
                         if (c .eq. e) then
                             Dijmkln(a,b,c,d,e,f) = Dijmkln(a,b,c,d,e,f) + 0.25*factor * dx(a)*dx(b)*dx(d)*dx(f)
                         endif
                         if (c .eq. f) then
                             Dijmkln(a,b,c,d,e,f) = Dijmkln(a,b,c,d,e,f) - 0.25*factor * dx(a)*dx(b)*dx(d)*dx(e)
                         endif

                           enddo
                       enddo
                   enddo
               enddo
           enddo
       enddo
    endif
        
    process_dedr = KIM_STATUS_OK

    end function process_dedr

    integer(c_int) function process_d2edr2(pkim,d2edr2,pr,pdx,pi,pj)bind(c)
    implicit none

    ! Passed variables
    type(c_ptr),   intent(in) :: pkim,pr,pdx
    real(c_double),intent(in) :: d2edr2
    type(c_ptr),   intent(in) :: pi,pj

    ! Local variables
    integer(c_int), pointer :: i(:);
    integer(c_int), pointer :: j(:);
    real(c_double), pointer :: dx(:,:)
    real(c_double), pointer :: r(:)
    real(c_double) :: R1(DIM),R2(DIM),R3(DIM),R4(DIM)

    !integer(c_int), pointer ::  i(:)
    !integer(c_int) j(2)
    call c_f_pointer(pi, i, [2])
    call c_f_pointer(pj, j, [2])
    call c_f_pointer(pdx,dx,[3,2])
    call c_f_pointer(pr,r,[2])

    R1 = dx(:,1)
    R2 = dx(:,2)
    if (j(1) .eq. j(2)) then
        ! this is a pair
        Cijkl_global = Cijkl_global + d2edr2*tensor_prod(R1,R1,R2,R2)/(r(1)*r(2))
    else
        ! this is 3-tuple
        Cijkl_global = Cijkl_global + 0.5*d2edr2*tensor_prod(R1,R1,R2,R2)/(r(1)*r(2))
    endif

    ! In order to understand the following code one has to understand how KIM API's process_d2edr2 functions. In order to obtain
    ! d2V/dr_{pq} dr_{rs}, the model may broadcast different values corresponding to 
    ! i[1]=p, j[1]=q, i[2]=r, j[2]=s
    ! i[1]=q, j[1]=q, i[2]=r, j[2]=s
    ! i[1]=p, j[1]=q, i[2]=s, j[2]=r
    ! i[1]=q, j[1]=q, i[2]=s, j[2]=r
    ! d2V/dr_{pq} dr_{rs} is equal to the sum of the four values obtained in the above broadcasts. Note that sum models may not
    ! broadcast four times. For example, the EAM model is written such that i[1] is always equal to i[2]. 
    ! Therefore, for instance, the value of d2V/dr_{pq} dr_{ps} is obtained in a single broadcast corresponding to i[1]=i[2]=p,j[1]=q,j[2]=s. On
    ! the other hand, the value of d2V/dr_{pq} dr_{pq} is obtained as a sum of two braodcasted values corresponding to
    ! i[1]=i[2]=p,j[1]=j[2]=q.

    if (i(1) .eq. 1 .or. i(2) .eq. 1 .or. j(1) .eq. 1 .or. j(2) .eq. 1) then
        if (i(1) .eq. 1) then
            !-------------------------------------------------------------------------------
            ! alpha = i(1)
            ! beta = j(1)
            ! gamma = i(2)
            ! delta = j(2)
            !
            Cijkl = Cijkl + 0.5d0*d2edr2*tensor_prod(R1,R1,R2,R2)/(r(1)*r(2))
            R3 = prcoor(:,i(1)) - prcoor(:,i(2))
            R4 = prcoor(:,i(1)) - prcoor(:,j(2)) 
            Dijmkln = &
                Dijmkln - (1.d0/8.d0) * d2edr2 * (tensor_prod_6(R1,R1,R3,R2,R2,R3)+tensor_prod_6(R1,R1,R4,R2,R2,R4))/(r(1)*r(2))
        elseif (j(1) .eq. 1) then
            !-------------------------------------------------------------------------------
            ! alpha = j(1)
            ! beta = i(1)
            ! gamma = j(2)
            ! delta = i(2)
            !
            Cijkl = Cijkl + 0.5d0*d2edr2*tensor_prod(R1,R1,R2,R2)/(r(1)*r(2))
            R3 = prcoor(:,j(1)) - prcoor(:,j(2))
            R4 = prcoor(:,j(1)) - prcoor(:,i(2))
            Dijmkln = &
                Dijmkln - (1.d0/8.d0) * d2edr2 * (tensor_prod_6(R1,R1,R3,R2,R2,R3)+tensor_prod_6(R1,R1,R4,R2,R2,R4))/(r(1)*r(2))
        endif
    endif

    process_d2edr2 = KIM_STATUS_OK
    end function process_d2edr2

    function tensor_prod(vec1,vec2,vec3,vec4)
    implicit none
    real(c_double),intent(in) :: vec1(DIM),vec2(DIM),vec3(DIM),vec4(DIM)
    real(c_double) :: tensor_prod(DIM,DIM,DIM,DIM)

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
!**        CLUSTER
!**
!*******************************************************************************

!-------------------------------------------------------------------------------
!
! Main program
!
!-------------------------------------------------------------------------------
program test_elastic_constants
  use, intrinsic :: iso_c_binding
  use KIM_API_F03
  use mod_neighborlist
  use mod_process, only : Cnorm,lambda,Cijkl,Cijkl_global,Dijmkln,setup_process_KIM_access,prcoor
  implicit none
  integer(c_int), parameter :: cd = c_double ! used for literal constants

  !real(c_double), parameter :: FCCspacing     = 4.32650374152 ! For Si in  angstroms
  !real(c_double), parameter :: FCCspacing     = 4.050d0 ! For Al in angstroms
  !real(c_double), parameter :: FCCspacing     = 5.29216d0 ! For Ar in angstroms
  !real(c_double), parameter :: FCCspacing     = 3.61500004445 ! For Cu in angstroms - mishin
  !real(c_double), parameter :: FCCspacing     = 3.63902985156d0! For Cu in angstroms - king
  real(c_double), parameter :: FCCspacing     = 3.61500009724d0! For Cu in angstroms - Ackland
  !real(c_double), parameter :: FCCspacing     = 3.61495904379d0! For Cu in angstroms - Zhou
  integer(c_int), parameter :: nCellsPerSide  = 7
  integer(c_int), parameter :: DIM            = 3
  integer(c_int), parameter :: ASpecies         = 1
  integer(c_int), parameter :: &
       N = 4*(nCellsPerSide)**3 + 6*(nCellsPerSide)**2 + 3*(nCellsPerSide) + 1
       !N = 4*(nCellsPerSide)**3
  integer(c_int), parameter :: SizeOne        = 1

  !
  ! neighbor list
  !
  type(neighObject_type), target :: neighObject
  real(c_double), parameter      :: cutpad = 0.750d0 ! in angstroms

  !
  ! KIM variables
  !
  character(len=KIM_KEY_STRING_LENGTH) :: testname = "test_elastic_constants"
  character(len=KIM_KEY_STRING_LENGTH) :: testkimfile = "descriptor.kim"
  character(len=KIM_KEY_STRING_LENGTH) :: modelname
  character(len=KIM_KEY_STRING_LENGTH) :: NBC_Method
  ! 0- NEIGH_RVEC_H, 1- NEIGH_PURE_H, 2- NEIGH_RVEC_F, 3- NEIGH_PURE_F,
  ! 4- MI_OPBC_H,    5- MI_OPBC_F,    6- CLUSTER
  integer(c_int) nbc

  type(c_ptr)    :: pkim
  integer(c_int) :: ier, idum
  integer(c_int) :: middleDum
  integer(c_int) :: I
  integer(c_int), pointer :: numberOfParticles;   type(c_ptr) :: pnParts
  integer(c_int), pointer :: numContrib;          type(c_ptr) :: pnumContrib
  integer(c_int), pointer :: numberOfSpecies;     type(c_ptr) :: pnOfSpecies
  integer(c_int), pointer :: particleSpecies(:);  type(c_ptr) :: pparticleSpecies
  real(c_double), pointer :: cutoff;              type(c_ptr) :: pcutoff
  real(c_double), pointer :: energy;              type(c_ptr) :: penergy
  real(c_double), pointer :: coords(:,:);         type(c_ptr) :: pcoor
  real(c_double), pointer :: forces(:,:);         type(c_ptr) :: pforces
  real(c_double), pointer :: boxSideLengths(:);   type(c_ptr) :: pboxSideLengths


  integer(c_int) :: p,q,r,s,t,u

  Cijkl = 0.0
  Cijkl_global = 0.0
  Dijmkln = 0.0

  ! Get KIM Model name to use
  print '("Please enter a valid KIM model name: ")'
  !read(*,*) modelname
  !modelname='Three_Body_Stillinger_Weber_Hauch_Brittle_Si__MO_119167353542_001'
  !modelname='EAM_Dynamo_Zope_Mishin_Al__MO_664470114311_001'
  !modelname='EAM_Dynamo_Mishin_Mehl_Cu__MO_346334655118_001'
  !modelname='EAM_Dynamo_Mendelev_King_Cu__MO_748636486270_001'
  modelname='EAM_Dynamo_Ackland_Tichy_Cu__MO_179025990738_001'
  !modelname='EAM_Dynamo_Zhou_Johnson_Cu__MO_127245782811_001'
  ! modelname='Pair_LJ_Smooth_Central_Model'

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
  elseif (index(NBC_Method,"CLUSTER").eq.1) then
     nbc = 6
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
  if (nbc.lt.6) allocate(neighObject%neighborList(N+1, N))
  if (nbc.eq.0 .or. nbc.eq.2) then
     allocate(neighObject%RijList(DIM,N+1, N))
  endif
  !
  if (nbc.ne.6) then
     ier = kim_api_set_data(pkim, "neighObject", SizeOne, c_loc(neighObject))
     if (ier.lt.KIM_STATUS_OK) then
        idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                    "kim_api_set_data", ier)
        stop
     endif
  endif

  if (nbc.ne.6) then
     ier = kim_api_set_method(pkim, "get_neigh", SizeOne, c_funloc(get_neigh))
     if (ier.lt.KIM_STATUS_OK) then
        idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                    "kim_api_set_method", ier)
        stop
     endif
  endif

  ! call model's init routine
  ier = kim_api_model_init(pkim)
  if (ier.lt.KIM_STATUS_OK) then
     idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                 "kim_api_model_init", ier)
     stop
  endif

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
  if ((nbc.eq.0).or.(nbc.eq.1).or.(nbc.eq.4)) call c_f_pointer(pnumContrib, &
                                                               numContrib)
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
  if ((nbc.eq.0).or.(nbc.eq.1).or.(nbc.eq.4)) numContrib = N
  numberOfSpecies = ASpecies
  particleSpecies(:)    = kim_api_get_species_code(pkim, "Cu", ier)
  if (ier.lt.KIM_STATUS_OK) then
     idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                 "kim_api_get_species_code", ier)
     stop
  endif

  ! set up the cluster part positions
  call create_FCC_configuration(FCCspacing, nCellsPerSide, .false., coords, &
                                middleDum)
  call setup_process_KIM_access(pkim,coords,DIM,N)

 ! set boxSideLengths large enough to make the cluster isolated
  if (nbc.eq.4 .or. nbc.eq.5) boxSideLengths(:) = 60000.0_cd

  ! compute neighbor lists
  if (nbc.eq.0) then
     call NEIGH_RVEC_cluster_neighborlist(.true., N, coords, (cutoff+cutpad), &
                                          neighObject)
  elseif (nbc.eq.1) then
     call NEIGH_PURE_cluster_neighborlist(.true., N, coords, (cutoff+cutpad), &
                                          neighObject)
  elseif (nbc.eq.2) then
     call NEIGH_RVEC_cluster_neighborlist(.false., N, coords, (cutoff+cutpad), &
                                          neighObject)
  elseif (nbc.eq.3) then
     call NEIGH_PURE_cluster_neighborlist(.false., N, coords, (cutoff+cutpad), &
                                          neighObject)
  elseif (nbc.eq.4) then
     call MI_OPBC_cluster_neighborlist(.true., N, coords, (cutoff+cutpad), &
                                       boxSideLengths, neighObject)
  elseif (nbc.eq.5) then
     call MI_OPBC_cluster_neighborlist(.false., N, coords, (cutoff+cutpad), &
                                       boxSideLengths, neighObject)
  else
     ! nothing to do for CLUSTER
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
  print '("Forces:")'
  print '("Part     ' // &
  'X                        ' // &
  'Y                        ' // &
  'Z                        ")'
  !print '(I2,"   ",3ES25.15)', (I,forces(:,I),I=1,N)
  print *
  print '("Energy        = ",ES25.15)', energy

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

  Cijkl = Cijkl/((FCCspacing**3)/4.0)
  Cijkl_global = Cijkl_global/((FCCspacing*nCellsPerSide)**3)
  Dijmkln = Dijmkln/((FCCspacing**3)/4.d0)

  Cnorm=0.d0
  lambda = 0.d0
  do u=1,DIM
      do t=1,DIM
          do s=1,DIM
            do r=1,DIM
                do q=1,DIM
                    do p=1,DIM
                        lambda(r,u) = lambda(r,u) + (Dijmkln(p,q,r,s,t,u)+Dijmkln(p,q,u,s,t,r))*Cijkl(p,q,s,t)/2.d0
                        if (u .eq. 1 .and. r .eq. 1) then
                            Cnorm = Cnorm + Cijkl(p,q,s,t)*Cijkl(p,q,s,t)
                            write(92,*) p,q,r,s,Cijkl(p,q,s,t)
                            write(93,*) p,q,r,s,Cijkl_global(p,q,s,t) 
                        endif
                        write(94,*) p,q,r,s,t,u,Dijmkln(p,q,r,s,t,u)
                    enddo
                enddo
            enddo
          enddo
      enddo
  enddo

  do r=1,DIM
      do u=1,DIM
          write(95,*) r,u,lambda(r,u)/Cnorm
      enddo
  enddo
  deallocate(prcoor)
  stop
end program test_elastic_constants

!-------------------------------------------------------------------------------
!
! MI_OPBC_cluster_neighborlist : construct a half or full neighbor list using
!                                the particle coordinates in coords()
!
!-------------------------------------------------------------------------------
subroutine MI_OPBC_cluster_neighborlist(half, numberOfParticles, coords, rcut, &
                                        boxSideLengths, neighObject)
  use, intrinsic :: iso_c_binding
  use KIM_API_F03
  use mod_neighborlist
  implicit none
  integer(c_int), parameter :: cd = c_double ! used for literal constants

  !-- Transferred variables
  logical,        intent(in)                  :: half
  integer(c_int), intent(in)                  :: numberOfParticles
  real(c_double), dimension(3,numberOfParticles), &
                  intent(in)                  :: coords
  real(c_double), intent(in)                  :: rcut
  real(c_double), dimension(3), intent(in)    :: boxSideLengths
  type(neighObject_type),       intent(inout) :: neighObject

  !-- Local variables
  integer(c_int) i, j, a
  real(c_double) dx(3)
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
           if (i.ne.j) then
              if ( (j .gt. i) .or. ((.not. half) .AND. (i.ne.j)) ) then
                  ! part j is a neighbor of part i
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

end subroutine MI_OPBC_cluster_neighborlist

!-------------------------------------------------------------------------------
!
! NEIGH_PURE_cluster_neighborlist
!
!-------------------------------------------------------------------------------
subroutine NEIGH_PURE_cluster_neighborlist(half, numberOfParticles, coords, &
                                           cutoff, neighObject)
  use, intrinsic :: iso_c_binding
  use KIM_API_F03
  use mod_neighborlist
  implicit none

  !-- Transferred variables
  logical,        intent(in)            :: half
  integer(c_int), intent(in)            :: numberOfParticles
  real(c_double), dimension(3,numberOfParticles), &
                  intent(in)            :: coords
  real(c_double), intent(in)            :: cutoff
  type(neighObject_type), intent(inout) :: neighObject

  !-- Local variables
  integer(c_int) i, j, a
  real(c_double) dx(3)
  real(c_double) r2
  real(c_double) cutoff2

  cutoff2 = cutoff**2

  do i=1,numberOfParticles
     a = 1
     do j=1,numberOfParticles
        dx(:) = coords(:, j) - coords(:, i)
        r2 = dot_product(dx, dx)
        if (r2.le.cutoff2) then
           ! part j is a neighbor of part i
           if ( (j .gt. i) .OR. ((.not. half) .AND. (i.ne.j)) ) then
               a = a+1
               neighObject%neighborList(a,i) = j
           endif
        endif
     enddo
     ! part i has a-1 neighbors
     neighObject%neighborList(1,i) = a-1
  enddo

  return

end subroutine NEIGH_PURE_cluster_neighborlist

!-------------------------------------------------------------------------------
!
! NEIGH_RVEC_cluster_neighborlist
!
!-------------------------------------------------------------------------------
subroutine NEIGH_RVEC_cluster_neighborlist(half, numberOfParticles, coords, &
                                           cutoff, neighObject)
  use, intrinsic :: iso_c_binding
  use KIM_API_F03
  use mod_neighborlist
  implicit none

  !-- Transferred variables
  logical,        intent(in)            :: half
  integer(c_int), intent(in)            :: numberOfParticles
  real(c_double), dimension(3,numberOfParticles), &
                  intent(in)            :: coords
  real(c_double), intent(in)            :: cutoff
  type(neighObject_type), intent(inout) :: neighObject

  !-- Local variables
  integer(c_int) i, j, a
  real(c_double) dx(3)
  real(c_double) r2
  real(c_double) cutoff2

  cutoff2 = cutoff**2

  do i=1,numberOfParticles
     a = 1
     do j=1,numberOfParticles
        dx(:) = coords(:, j) - coords(:, i)
        r2 = dot_product(dx, dx)
        if (r2.le.cutoff2) then
           if ((half .and. i.lt.j) .or. (.not.half .and. i.ne.j)) then
              ! part j is a neighbor of part i
              a = a+1
              neighObject%neighborList(a,i) = j
              neighObject%RijList(:,a-1,i) = dx
           endif
        endif
     enddo
     ! part i has a-1 neighbors
     neighObject%neighborList(1,i) = a-1
  enddo

  return

end subroutine NEIGH_RVEC_cluster_neighborlist

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
  use, intrinsic :: iso_c_binding
  use KIM_API_F03
  implicit none
  integer(c_int), parameter :: cd = c_double ! used for literal constants

  !-- Transferred variables
  real(c_double), intent(in)  :: FCCspacing
  integer(c_int), intent(in)  :: nCellsPerSide
  logical,        intent(in)  :: periodic
  real(c_double), intent(out) :: coords(3,*)
  integer(c_int), intent(out) :: MiddlePartId
  !
  ! cluster setup variables
  !
  real(c_double) FCCshifts(3,4)
  real(c_double) latVec(3)
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

  MiddlePartID = 1 ! Always put middle particle as #1
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
