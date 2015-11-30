module mod_crystal
    use mod_global
    use mod_lattice, only : Lattice
    implicit none
    type Crystal
        type(Lattice)    :: LatticeType
        integer          :: numSpecies
        character(len=2),allocatable :: basisAtomsSpecies(:)
        integer(c_int),allocatable   :: basisAtomsSpeciesCode(:)
        real(c_double),allocatable   :: basisAtomsCoords(:,:)
    end type Crystal

    contains

        subroutine allocateCrystal(crystalType,numBasisAtoms)
            implicit none

            type(Crystal),intent(inout) :: crystalType
            integer, intent(in)         :: numBasisAtoms

            allocate(crystalType%basisAtomsSpecies(numBasisAtoms))
            allocate(crystalType%basisAtomsSpeciesCode(numBasisAtoms))
            allocate(crystalType%basisAtomsCoords(DIM,numBasisAtoms))

        end subroutine allocateCrystal

        subroutine fillCrystal(crystalType,numSpecies,basisAtomsCoords,basisAtomsSpecies)
            implicit none

            type(Crystal),intent(inout) :: crystalType
            integer, intent(in)         :: numSpecies
            real(c_double),intent(inout):: basisAtomsCoords(:,:)
            character(len=*)            :: basisAtomsSpecies(:)

            crystalType%numSpecies    = numSpecies
            crystalType%basisAtomsSpecies = basisAtomsSpecies
            crystalType%basisAtomsCoords = basisAtomsCoords
        end subroutine fillCrystal

        subroutine destroyCrystal(crystalType)
            implicit none

            type(Crystal),intent(inout) :: crystalType

            deallocate(crystalType%basisAtomsSpecies)
            deallocate(crystalType%basisAtomsCoords)
        end subroutine destroyCrystal

        subroutine constructGlobalCrystal(numCellsPerSide,cT,coords,speciesCode)
            implicit none
            integer, intent(in)         :: numCellsPerSide
            type(Crystal), intent(in)   :: cT
            real(c_double), intent(out) :: coords(DIM,*)
            integer(c_int)              :: speciesCode(*)

            ! Local variables
            integer(c_int) :: i,j,k,m,a
            real(c_double) :: latVec(DIM)

            ! make sure that crystalType is ready before proceeding
            ! Nikhil
            a = 1
            do i=1,numCellsPerSide
                do j=1,numCellsPerSide
                    do k=1,numCellsPerSide
                        latVec = (i-1) * cT%latticeType%latVec(:,1) + &
                                 (j-1) * cT%latticeType%latVec(:,2) + &
                                 (k-1) * cT%latticeType%latVec(:,3) 
                        do m=1,size(cT%basisAtomsCoords,2)
                            a = a+1
                            coords(:,a) = latVec + cT%basisAtomsCoords(:,m)
                            speciesCode(a) = cT%basisAtomsSpeciesCode(m)
                            if ((i.eq.numCellsPerside/2+1).and.(j.eq.numCellsPerSide/2+1) .and. &
                                (k.eq.numCellsPerSide/2+1) ) then
                                coords(:,m) = latVec + cT%basisAtomsCoords(:,m) ! put middle particle as #1
                                speciesCode(m) = cT%basisAtomsSpeciesCode(m)
                                a = a - 1
                            endif
                       enddo
                    enddo
                enddo
            enddo
        end subroutine constructGlobalCrystal

end module mod_crystal

