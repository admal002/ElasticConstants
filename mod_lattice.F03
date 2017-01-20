module mod_lattice
    use mod_global
    implicit none
    type Lattice
        real(c_double), dimension(DIM,DIM) :: latVec
        real(c_double) :: latVol
    end type Lattice

    contains

        function constructFCCLattice(latticeSpacing)
            implicit none
            ! Passed variables
            type(Lattice)  :: constructFCCLattice
            real(c_double), intent(in) :: latticeSpacing

            ! Local variables
            real(c_double), dimension(DIM,DIM) :: latVec
            
            latVec(:,1) = latticeSpacing*(/0.0_cd, 0.5_cd, 0.5_cd/)
            latVec(:,2) = latticeSpacing*(/0.5_cd, 0.0_cd, 0.5_cd/)
            latVec(:,3) = latticeSpacing*(/0.5_cd, 0.5_cd, 0.0_cd/)

            constructFCCLattice%latVec = latVec
            constructFCCLattice%latVol = (latticeSpacing**3)/4._cd

        end function constructFCCLattice

        function constructBCCLattice(latticeSpacing)
            implicit none
            type(Lattice)  :: constructBCCLattice
            real(c_double), intent(in) :: latticeSpacing

            ! Local variables
            real(c_double), dimension(DIM,DIM) :: latVec
            
            latVec(:,1) = latticeSpacing*(/-0.5_cd, 0.5_cd, 0.5_cd/)
            latVec(:,2) = latticeSpacing*(/0.5_cd, -0.5_cd, 0.5_cd/)
            latVec(:,3) = latticeSpacing*(/0.5_cd, 0.5_cd, -0.5_cd/)

            constructBCCLattice%latVec = latVec
            constructBCCLattice%latVol = (latticeSpacing**3)/2._cd

        end function constructBCCLattice

        function constructHCPLattice(latticeSpacing,latticeSpacing_z)
            implicit none
            type(Lattice)  :: constructHCPLattice
            real(c_double), intent(in) :: latticeSpacing,latticeSpacing_z

            ! Local variables
            real(c_double), dimension(DIM,DIM) :: latVec
            
            ! A1 passes through atom
            latVec(:,1) = latticeSpacing*(/1.0_cd, 0.0_cd, 0.0_cd/)
            latVec(:,2) = latticeSpacing*(/-0.5_cd, 0.5_cd*sqrt(3._cd), 0.0_cd/)

            latVec(:,3) = latticeSpacing_z*(/0.0_cd, 0.0_cd, 1.0_cd/)

            constructHCPLattice%latVec = latVec

            ! Note that this is actually a hexagonal lattice
            ! Therefore, there is a factor of 2 instead of 4 in the lattice
            ! volume
            constructHCPLattice%latVol = sqrt(3._cd) * (latticeSpacing**2) * latticeSpacing_z/2._cd

        end function constructHCPLattice

        function constructOrthLattice(latticeSpacing_x,latticeSpacing_y,latticeSpacing_z)
            implicit none
            type(Lattice)  :: constructOrthLattice
            real(c_double), intent(in) :: latticeSpacing_x,latticeSpacing_y,latticeSpacing_z

            ! Local variables
            real(c_double), dimension(DIM,DIM) :: latVec
            
            ! A2 passes through atom
            latVec(:,1) = latticeSpacing_x*(/1._cd, 0._cd, 0._cd/)
            latVec(:,2) = latticeSpacing_y*(/0._cd, 1._cd, 0._cd/)
            latVec(:,3) = latticeSpacing_z*(/0._cd, 0._cd, 1._cd/)

            constructOrthLattice%latVec = latVec

            constructOrthLattice%latVol = latticeSpacing_x*latticeSpacing_y*latticeSpacing_z

        end function constructOrthLattice
end module mod_lattice


