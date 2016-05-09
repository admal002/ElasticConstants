#include "KIM_API_status.h"
#define THIS_FILE_NAME __FILE__
#define TRUEFALSE(TRUTH) merge(1,0,(TRUTH))

module mod_atomistic
    use mod_global
    use mod_fsge, only : Fsge

    type NeighObject_type
        !integer(c_int), pointer :: neighborList(:,:)
        !real(c_double), pointer :: RijList(:,:,:)
        integer(c_int),allocatable :: neighborList(:,:)
        real(c_double),allocatable :: RijList(:,:,:)
    end type NeighObject_type

    type ::  Atomistic
        real(c_double),allocatable   :: coords(:,:)
        integer(c_int)               :: numParticles,numContributingParticles,nbc,numSpecies,&
                                        processdEdrFlag,processd2Edr2Flag
        type(neighObject_type)       :: neighObject
        integer(c_int),allocatable   :: particleSpecies(:)
        real(c_double)               :: boxSideLengths(DIM),cutoff,energy
        type(Fsge)                   :: fsgeType,fsgeTypeNumerical
    end type Atomistic

    type(Atomistic), pointer :: pA
    contains

        subroutine constructAtomistic(aS,numParticles,nbc)
            use mod_fsge, only : initializeFsge
            implicit none

            integer(c_int), intent(in)    :: numParticles,nbc
            type(Atomistic),intent(inout) :: aS

            aS%numParticles = numParticles

            ! Allocate memory to atomistic system
            allocate(aS%coords(DIM,numParticles))
            allocate(aS%particleSpecies(numParticles))

            allocate(aS%neighObject%neighborList(numParticles+1, numParticles))
            if (nbc.eq.0 .or. nbc.eq.2) then
                allocate(aS%neighObject%RijList(DIM,numParticles+1, numParticles))
            endif
            call initializeFsge(aS%fsgeType)
            call initializeFsge(aS%fsgeTypeNumerical)

        end subroutine constructAtomistic

        function process_dedr(pkim,deidr,r,pdx,i,j)bind(c)
        use KIM_API_F03
        use mod_matrix
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

        R3 = (pA%coords(:,i) + pA%coords(:,j))/2._cd

        pA%fsgeType%Tijm = pA%fsgeType%Tijm + deidr * tensor_prod(dx,dx,R3)/r
        pA%fsgeType%Eijkln = pA%fsgeType%Eijkln - deidr*tensor_prod(dx,dx,dx,dx,R3)/(r**3)
        pA%fsgeType%Cijkl = pA%fsgeType%Cijkl - deidr * tensor_prod(dx,dx,dx,dx)/(r**3)
        pA%fsgeType%Dijmkln = pA%fsgeType%Dijmkln - deidr * tensor_prod(dx,dx,R3,dx,dx,R3)/(r**3)
        pA%fsgeType%Dijmkln = pA%fsgeType%Dijmkln + deidr * tensorG_dot_tensorG(tensorG(dx,R3)) / r

        process_dedr = KIM_STATUS_OK
        return
        end function process_dedr

        function process_d2edr2(pkim,d2edr2,pr,pdx,pi,pj)bind(c)
        use KIM_API_F03
        use mod_matrix
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
        R3 = (pA%coords(:,i(1)) + pA%coords(:,j(1)))/2._cd
        R4 = (pA%coords(:,i(2)) + pA%coords(:,j(2)))/2._cd

        pA%fsgeType%Eijkln = pA%fsgeType%Eijkln + d2edr2*tensor_prod(R1,R1,R2,R2,R4)/(r(1)*r(2))
        pA%fsgeType%Cijkl  = pA%fsgeType%Cijkl + d2edr2*tensor_prod(R1,R1,R2,R2)/(r(1)*r(2))
        pA%fsgeType%Dijmkln= pA%fsgeType%Dijmkln + d2edr2*tensor_prod(R1,R1,R3,R2,R2,R4)/(r(1)*r(2))


        process_d2edr2 = KIM_STATUS_OK
        return
        end function process_d2edr2

        function tensorG(dx,R3)
        implicit none

        real(c_double) ::tensorG(DIM,DIM,DIM,DIM)
        real(c_double), intent(in) :: dx(DIM), R3(DIM)
        integer :: p,i,j,m

        tensorG = 0._cd

        do m=1,DIM
            do j=1,DIM
                do i=1,DIM
                    do p=1,DIM
                        if (p .eq. i) tensorG(p,i,j,m) = tensorG(p,i,j,m) + dx(j)*R3(m) + dx(m)*R3(j)
                        if (p .eq. j) tensorG(p,i,j,m) = tensorG(p,i,j,m) + dx(i)*R3(m) + dx(m)*R3(i)
                        if (p .eq. m) tensorG(p,i,j,m) = tensorG(p,i,j,m) - dx(i)*R3(j) - dx(j)*R3(i)
                    enddo
                enddo
            enddo
        enddo
        tensorG = 0.5_cd * tensorG
        end function tensorG

        function tensorG_dot_tensorG(tensorG)
        implicit none

        real(c_double) ::tensorG_dot_tensorG(DIM,DIM,DIM,DIM,DIM,DIM)
        real(c_double), intent(in) :: tensorG(DIM,DIM,DIM,DIM)
        integer :: i,j,m,k,l,n,p

        tensorG_dot_tensorG = 0._cd
        do i = 1,DIM
            do j = 1,DIM
                do m = 1,DIM
                    do k = 1,DIM
                        do l = 1,DIM
                            do n = 1,DIM
                                do p = 1,DIM
                                    tensorG_dot_tensorG(i,j,m,k,l,n) = tensorG_dot_tensorG(i,j,m,k,l,n) + &
                                                                        tensorG(p,i,j,m) * tensorG(p,k,l,n)
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        end function tensorG_dot_tensorG

        subroutine constructNeighList(aS)
            implicit none
            type(Atomistic) :: aS

            select case(aS%nbc)
            case(0)
                call NEIGH_RVEC_neighborlist(.true.,aS)
            case(1)
                call NEIGH_PURE_neighborlist(.true.,aS)
            case(2)
                call NEIGH_RVEC_neighborlist(.false.,aS)
            case(3)
                call NEIGH_PURE_neighborlist(.false.,aS)
            case(4)
                call MI_OPBC_neighborlist(.true.,aS)
            case(5)
                call MI_OPBC_neighborlist(.false.,aS)
            end select

        end subroutine constructNeighList

    subroutine NEIGH_PURE_neighborlist(half,aS)
        use mod_global
        implicit none

        !-- Transferred variables
        type(Atomistic)             :: aS
        logical,        intent(in)  :: half

        !-- Local variables
        integer(c_int) i, j, a
        real(c_double) dx(DIM)
        real(c_double) r2
        real(c_double) cutoff2

        cutoff2 = aS%cutoff**2

        do i=1,aS%numParticles
            a = 1
            do j=1,aS%numParticles
                dx(:) = aS%coords(:, j) - aS%coords(:, i)
                r2 = dot_product(dx, dx)
                if (r2.le.cutoff2) then

                   if (half) then
                       if ( ((i.le. aS%numContributingParticles) .or. &
                             (j.le. aS%numContributingParticles)) .and. &
                            (i .lt. j) ) then
                          a = a+1
                          aS%neighObject%neighborList(a,i) = j
                       endif
                   else
                       if (i .le. aS%numContributingParticles .and. i.ne.j) then
                          a = a+1
                          aS%neighObject%neighborList(a,i) = j
                       endif
                   endif

                endif
            enddo
            ! part i has a-1 neighbors
            aS%neighObject%neighborList(1,i) = a-1
        enddo

        return

    end subroutine NEIGH_PURE_neighborlist

    subroutine MI_OPBC_neighborlist(half,aS)
        use mod_global
        implicit none

        !-- Transferred variables
        logical,        intent(in)            :: half
        type(Atomistic)                       :: aS

        !-- Local variables
        integer(c_int) i, j, a
        real(c_double) dx(DIM)
        real(c_double) r2
        real(c_double) rcut2

        rcut2 = aS%cutoff**2

        do i=1,aS%numParticles
            a = 1
            do j=1,aS%numParticles
                dx(:) = aS%coords(:, j) - aS%coords(:, i)
                where (abs(dx) > 0.5_cd*aS%boxSideLengths)  ! apply PBC
                   dx = dx - sign(aS%boxSideLengths,dx)
                endwhere
                r2 = dot_product(dx, dx)
                if (r2.le.rcut2) then
                   ! part j is a neighbor of part i
                   if (half) then
                       if ( ((i .le. aS%numContributingParticles) .or. &
                             (j .le. aS%numContributingParticles)) .and. &
                            (i .lt. j) ) then
                          a = a+1
                          aS%neighObject%neighborList(a,i) = j
                       endif
                   else
                       if (i .le. aS%numContributingParticles .and. i.ne.j) then
                          a = a+1
                          aS%neighObject%neighborList(a,i) = j
                       endif
                   endif
                endif
            enddo
            ! part i has a-1 neighbors
            aS%neighObject%neighborList(1,i) = a-1
        enddo

        return

    end subroutine MI_OPBC_neighborlist

    subroutine NEIGH_RVEC_neighborlist(half,aS)
        use mod_global
        implicit none

        !-- Transferred variables
        logical,        intent(in)            :: half
        type(Atomistic)                       :: aS

        !-- Local variables
        integer(c_int) i, j, a
        real(c_double) dx(DIM)
        real(c_double) r2
        real(c_double) cutoff2

        cutoff2 = aS%cutoff**2

        do i=1,aS%numParticles
            a = 1
            do j=1,aS%numParticles
                dx(:) = aS%coords(:, j) - aS%coords(:, i)
                r2 = dot_product(dx, dx)
                if (r2.le.cutoff2) then
                    if (half) then
                        if ( ((i .le. aS%numContributingParticles) .or.  &
                              (j .le. aS%numContributingParticles)) .and. &
                              (i .lt. j) ) then
                            a = a+1
                            aS%neighObject%neighborList(a,i) = j
                            aS%neighObject%RijList(:,a-1,i) = dx
                        endif
                    else
                        if (i .le. aS%numContributingParticles .and. i.ne.j) then
                          a = a+1
                          aS%neighObject%neighborList(a,i) = j
                          aS%neighObject%RijList(:,a-1,i) = dx
                        endif
                    endif
                endif
             enddo
             ! part i has a-1 neighbors
             aS%neighObject%neighborList(1,i) = a-1
        enddo

        return

    end subroutine NEIGH_RVEC_neighborlist

        subroutine destroyAtomistic(aS)
            implicit none
            type(Atomistic) :: aS

            if (allocated(aS%coords)) deallocate(aS%coords)
            if (allocated(aS%particleSpecies)) deallocate(aS%particleSpecies)
            if (allocated(aS%neighObject%neighborList)) deallocate(aS%neighObject%neighborList)
            if (allocated(aS%neighObject%RijList)) deallocate(aS%neighObject%RijList)
            if (allocated(aS%fsgeType%IndComponentsC)) deallocate(aS%fsgeType%IndComponentsC)
            if (allocated(aS%fsgeType%IndComponentsD)) deallocate(aS%fsgeType%IndComponentsD)
            if (allocated(aS%fsgeTypeNumerical%IndComponentsC)) deallocate(aS%fsgeTypeNumerical%IndComponentsC)
            if (allocated(aS%fsgeTypeNumerical%IndComponentsD)) deallocate(aS%fsgeTypeNumerical%IndComponentsD)
        end subroutine destroyAtomistic

end module mod_atomistic



