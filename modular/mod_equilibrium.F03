#include "KIM_API_status.h"
#define THIS_FILE_NAME __FILE__
module mod_equilibrium
    use mod_global
    use mod_crystal, only : Crystal
    use mod_atomistic
    use mod_lattice, only : constructFCCLattice, constructBCCLattice, constructHCPLattice
    use mod_kim

    contains

        function equilibriumCrystal(cT,crystalStructure,testname,modelname)
            implicit none

            ! Passed variables
            character(len=KIM_KEY_STRING_LENGTH) :: testname,modelname
            character(len=*) :: crystalStructure
            type(Crystal) :: cT
            type(Crystal) :: equilibriumCrystal

            ! Local variables
            type(c_ptr)    :: pkim
            type(Atomistic),pointer :: atomisticSystem
            integer(c_int) :: nbc,numBasisAtoms,numParticles
            integer(c_int), parameter :: numCellsPerSide  = 8
            real(c_double) :: latticeSpacing, latticeSpacing_z


            ! Perform the hand-shake and determine the boundary conditions
            call setup_kim(pkim,testname,modelname,nbc)

            numBasisAtoms = size(cT%basisAtomsCoords(DIM,:))
            numParticles = numBasisAtoms * (numCellsPerSide**(DIM))

            ! Allocate memory to the kim object
            call allocate_kim(pkim,numParticles,cT%numSpecies)

            ! Initialize KIM
            call initialize_kim(pkim)


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Now prepare the atomistic system

            atomisticSystem = constructAtomistic(numParticles,nbc)

            ! Broadcast the atomistic system to the kim api
            call broadcastAtomistic_kim(pkim,atomisticSystem)


            ! Fill in the coordinates using the relaxed crystal
            call constructGlobalCrystal(numCellsPerSide,cT,atomisticSystem%coords,atomisticSystem%particleSpecies)
            atomisticSystem%nbc = nbc
            atomisticSystem%numContributingParticles = 1
            atomisticSystem%boxSideLengths  = 1000._cd
            atomisticSystem%numSpecies = cT%numSpecies
            atomisticSystem%cutoff = getCutoff(pkim)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


            ! start loop here
                call constructGlobalCrystal(numCellsPerSide,                &
                                            cT,                             &
                                            atomisticSystem%coords,         &
                                            atomisticSystem%particleSpecies)
                call constructNeighList(atomisticSystem)

                ! Call model compute
                ier = kim_api_model_compute(pkim)
                if (ier.lt.KIM_STATUS_OK) then
                    idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                             "kim_api_model_compute", ier)
                    stop
                endif

                ! If converged, exit the loop
                ! Update the crystal
                select case(trim(crystalStructure))
                ! update the lattice spacing based on crystalStructure
                case('fcc')
                    cT%LatticeType = constructFCCLattice(latticeSpacing)
                    cT%basisAtomsCoords(:,1) = (/0._cd,0._cd,0._cd/)
                case('bcc')
                    cT%LatticeType = constructBCCLattice(latticeSpacing)
                    cT%basisAtomsCoords(:,1) = (/0._cd,0._cd,0._cd/)
                case('hcp')
                    cT%LatticeType = constructHCPLattice(latticeSpacing,latticeSpacing_z)
                    cT%basisAtomsCoords(:,1) = (/0._cd,0._cd,0._cd/)
                    cT%basisAtomsCoords(:,2) = (/latticeSpacing*(2._cd/3._cd),  &
                                                latticeSpacing*(1._cd/3._cd),   &
                                                latticeSpacing_z*0.5_cd/) 
                case('diamond')
                    cT%LatticeType = constructFCCLattice(latticeSpacing)
                    cT%basisAtomsCoords(:,1) = (/0._cd,0._cd,0._cd/)
                    cT%basisAtomsCoords(:,2) = (/0.25_cd,0.25_cd,0.25_cd/)*latticeSpacing
                end select 
            !loop

            equilibriumCrystal = cT

            call destroy_kim(pkim)
            call destroyAtomistic(atomisticSystem)

        end function equilibriumCrystal

end module mod_equilibrium
