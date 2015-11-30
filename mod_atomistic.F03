#include "KIM_API_status.h"
#define THIS_FILE_NAME __FILE__
#define TRUEFALSE(TRUTH) merge(1,0,(TRUTH))

module mod_atomistic
    use mod_global

    type NeighObject_type
        integer(c_int), pointer :: neighborList(:,:)
        real(c_double), pointer :: RijList(:,:,:)
    end type NeighObject_type

    type ::  Atomistic
        real(c_double),allocatable   :: coords(:,:)
        integer(c_int)               :: numParticles,numContributingParticles,nbc,numSpecies
        type(neighObject_type)       :: neighObject
        integer(c_int),allocatable   :: particleSpecies(:)
        real(c_double)               :: boxSideLengths,cutoff,energy
    end type Atomistic

    contains

        subroutine constructAtomistic(pAS,numParticles,nbc)
            implicit none

            integer(c_int), intent(in) :: numParticles,nbc
            type(Atomistic)            :: pAS

            pAS%numParticles = numParticles

            ! Allocate memory to atomistic system
            allocate(pAS%coords(DIM,numParticles))
            allocate(pAS%particleSpecies(numParticles))

            allocate(pAS%neighObject%neighborList(numParticles+1, numParticles))
            if (nbc.eq.0 .or. nbc.eq.2) then
                allocate(pAS%neighObject%RijList(DIM,numParticles+1, numParticles))
            endif

        end subroutine constructAtomistic

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
                       if ( ((i.eq.1) .or. (j.eq.1)) .and. &
                            (i .lt. j) ) then
                          a = a+1
                          aS%neighObject%neighborList(a,i) = j
                       endif
                   else
                       if (i.eq.1 .and. i.ne.j) then
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
                       if ( ((i.eq. 1) .or. (j.eq. 1)) .and. &
                            (i .lt. j) ) then
                          a = a+1
                          aS%neighObject%neighborList(a,i) = j
                       endif
                   else
                       if (i.eq. 1 .and. i.ne.j) then
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
                        if ( ((i.eq. 1) .or. (j.eq. 1)) .and. &
                            (i .lt. j) ) then
                            a = a+1
                            aS%neighObject%neighborList(a,i) = j
                            aS%neighObject%RijList(:,a-1,i) = dx
                        endif
                    else
                        if (i.eq. 1 .and. i.ne.j) then
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

            deallocate(aS%coords)
            deallocate(aS%particleSpecies)
            deallocate(aS%neighObject%neighborList)
            deallocate(aS%neighObject%RijList)
        end subroutine destroyAtomistic
end module mod_atomistic



