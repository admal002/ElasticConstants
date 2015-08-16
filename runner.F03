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
!    Nikhil Chandra Admal
!

!
! Release: This file is part of the kim-api.git package.
!
!-------------------------------------------------------------------------------
! Notes:
! 0. To do: Include the name of the test driver in the subroutine
!           WriteDescriptorFile
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
real(c_double) :: Cijkl(DIM,DIM,DIM,DIM),Cijkl_global(DIM,DIM,DIM,DIM), Dijmkln(DIM,DIM,DIM,DIM,DIM,DIM)
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
!**        CLUSTER
!**
!*******************************************************************************

!-------------------------------------------------------------------------------
!
! Main program
!
!-------------------------------------------------------------------------------
program runner
  use, intrinsic :: iso_c_binding
  use KIM_API_F03
  use mod_neighborlist
  use mod_process, only : Cijkl,Cijkl_global,Dijmkln,setup_process_KIM_access,prcoor
  implicit none
  integer(c_int), parameter :: cd = c_double ! used for literal constants

  integer(c_int), parameter :: nCellsPerSide  = 7
  integer(c_int), parameter :: DIM            = 3
  integer(c_int), parameter :: ASpecies       = 1
  integer(c_int), parameter :: SizeOne        = 1

  integer(c_int)   :: N
  real(c_double)   :: LatticeSpacing,AtomVol
  character(len=3) :: CrystalStructure
  character(len=2) :: Species
  character(len=6) :: QuerySpecies

  !
  ! neighbor list
  !
  type(neighObject_type), target :: neighObject
  real(c_double), parameter      :: cutpad = 0.750d0 ! in angstroms

  !
  ! KIM variables
  !
  character(len=KIM_KEY_STRING_LENGTH) :: testname
  character(len=10000) testkimfile
  character(len=KIM_KEY_STRING_LENGTH) :: modelname
  character(len=KIM_KEY_STRING_LENGTH) :: NBC_Method
  ! 0- NEIGH_RVEC_H, 1- NEIGH_PURE_H, 2- NEIGH_RVEC_F, 3- NEIGH_PURE_F,
  ! 4- MI_OPBC_H,    5- MI_OPBC_F,    6- CLUSTER
  integer(c_int) nbc

  type(c_ptr)    :: pkim
  integer(c_int) :: ier, idum
  integer(c_int) :: middleDum
  integer(c_int), pointer :: numberOfParticles;   type(c_ptr) :: pnParts
  integer(c_int), pointer :: numContrib;          type(c_ptr) :: pnumContrib
  integer(c_int), pointer :: numberOfSpecies;     type(c_ptr) :: pnOfSpecies
  integer(c_int), pointer :: particleSpecies(:);  type(c_ptr) :: pparticleSpecies
  real(c_double), pointer :: cutoff;              type(c_ptr) :: pcutoff
  real(c_double), pointer :: energy;              type(c_ptr) :: penergy
  real(c_double), pointer :: coords(:,:);         type(c_ptr) :: pcoor
  real(c_double), pointer :: boxSideLengths(:);   type(c_ptr) :: pboxSideLengths

  Cijkl = 0.0
  Cijkl_global = 0.0
  Dijmkln = 0.0

  ! Get KIM Model name to use
  read(*,*) testname
  read(*,*) modelname

  read(*,*) CrystalStructure
  read(*,*) LatticeSpacing
  read(*,*) QuerySpecies

  Species = QuerySpecies(3:4)

  if (CrystalStructure .eq. 'fcc') then
      N = 4*(nCellsPerSide)**3 + 6*(nCellsPerSide)**2 + 3*(nCellsPerSide) + 1
      AtomVol = (LatticeSpacing**3)/4._cd
  elseif (CrystalStructure .eq. 'bcc') then
      N = 2*(nCellsPerSide)**3 + 3*(nCellsPerSide)**2 + 3*(nCellsPerSide) + 1
      AtomVol = (LatticeSpacing**3)/2._cd
  else
      write(*,*) "Error: Crystal structure ", CrystalStructure," not supported."
      stop
  endif

  call WriteDescriptorFile(testkimfile,Species)
  ! Initialize the KIM object
  ier = kim_api_string_init(pkim, trim(testkimfile), modelname)
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
       "energy",                      penergy,          1)                              
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

  if ((nbc.eq.4).or.(nbc.eq.5)) call c_f_pointer(pboxSideLengths, &
                                                 boxSideLengths, [DIM])
  call c_f_pointer(penergy, energy)

  ! Set values
  numberOfParticles = N
  if ((nbc.eq.0).or.(nbc.eq.1).or.(nbc.eq.4)) numContrib = N
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
                                    middleDum)
  elseif (CrystalStructure .eq. 'bcc') then
      call create_BCC_configuration(LatticeSpacing, nCellsPerSide, .false., coords, &
                                    middleDum)
  endif
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
  print '("This is Test             : ",A)', trim(testname)
  print '("Results for KIM Model    : ",A)', trim(modelname)
  print '("Crystal strucutre        : ",A)', CrystalStructure
  print '("Species                  : ",A)', species
  print '("Queried Lattice constant : ",F20.15)', LatticeSpacing
  print '("Using NBC                : ",A)', trim(NBC_Method)

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

  Cijkl = Cijkl/AtomVol
  Cijkl_global = Cijkl_global/((LatticeSpacing*nCellsPerSide)**3)
  Dijmkln = Dijmkln/AtomVol

  call write_results(Cijkl,Cijkl_global,Dijmkln)
  deallocate(prcoor)
  stop
end program runner

subroutine write_results(Cijkl,Cijkl_global,Dijmkln)
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), parameter  :: DIM = 3
    ! Transferred variables
    real(c_double), intent(in) :: Cijkl(DIM,DIM,DIM,DIM),Cijkl_global(DIM,DIM,DIM,DIM), Dijmkln(DIM,DIM,DIM,DIM,DIM,DIM)
    ! Local variables
    integer(c_int)             :: cunit,cgunit,dunit,runit
    integer,external           :: GetFreeUnit
    character(len=1)           :: cr

    cr = char(10)

    cunit = GetFreeUnit()
    open(unit=cunit,file='output/Cijkl.out',status='replace')
    cgunit = GetFreeUnit()
    open(unit=cgunit,file='output/Cijkl_global.out',status='replace')
    dunit = GetFreeUnit()
    open(unit=dunit,file='output/Dijmkln.out',status='replace')
    
    write(cunit,*) "c11 = ",Cijkl(1,1,1,1)
    write(cgunit,*) "c11 = ",Cijkl_global(1,1,1,1)
    write(cunit,*) "c12 = ",Cijkl(1,1,2,2)
    write(cgunit,*) "c12 = ",Cijkl_global(1,1,2,2)
    write(cunit,*) "c44 = ",Cijkl(1,2,1,2)
    write(cgunit,*) "c44 = ",Cijkl_global(1,2,1,2)

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
    close(cgunit)
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

subroutine WriteDescriptorFile(descriptorfile,species)
use, intrinsic :: iso_c_binding
implicit none

!-- Transferred variables
character(len=10000), intent(out)  :: descriptorfile
character(len=2),     intent(in)   :: species

!-- Local variables
integer(c_int) :: i
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
   '# Author: Automatically generated by the test driver'                             // cr // &
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
   'CLUSTER                 flag'                                                     // cr // &
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
                                                                                         cr // &
   'process_dEdr                  method       none       []'                         // cr // &
                                                                                         cr // &
   'process_d2Edr2                method       none       []'                         // cr // &
                                                                                         cr // &
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
   divider                                                                            // cr // char(0)

return

end subroutine WriteDescriptorFile

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
  use, intrinsic :: iso_c_binding
  use KIM_API_F03
  implicit none
  integer(c_int), parameter :: cd = c_double ! used for literal constants

  !-- Transferred variables
  real(c_double), intent(in)  :: LatticeSpacing
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
  use, intrinsic :: iso_c_binding
  use KIM_API_F03
  implicit none
  integer(c_int), parameter :: cd = c_double ! used for literal constants

  !-- Transferred variables
  real(c_double), intent(in)  :: LatticeSpacing
  integer(c_int), intent(in)  :: nCellsPerSide
  logical,        intent(in)  :: periodic
  real(c_double), intent(out) :: coords(3,*)
  integer(c_int), intent(out) :: MiddlePartId
  !
  ! cluster setup variables
  !
  real(c_double) FCCshifts(3,2)
  real(c_double) latVec(3)
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
