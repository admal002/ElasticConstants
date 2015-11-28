module mod_fsge
    use mod_global
    implicit none

    type Fsge
        real(c_double) :: Cijkl(DIM,DIM,DIM,DIM),           &
                          Dijmkln(DIM,DIM,DIM,DIM,DIM,DIM), &
                          Tijm(DIM,DIM,DIM),                &
                          Eijkln(DIM,DIM,DIM,DIM,DIM)
        integer(c_int), allocatable :: IndComponents(:,:)
    end type fsge

    contains
        subroutine contructCubicFsgeMaterial(fsgeType)
            implicit none

            type(Fsge) :: fsgeType

            
        end subroutine constructCubicFsgeMaterial

        subroutine voigt

        end subroutine voigt

        subroutine writeFsgeMaterial
        end subroutine writeFsgeMaterial

        subroutine destroyFsgeMaterial(fsgeType)
            implicit none
            type(Fsge) :: fsgeType

            deallocate(fsgeType%IndComponents)
        end subroutine destroyFsgeMaterial
end module mod_fsge
