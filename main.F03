#include "KIM_API_status.h"
#define THIS_FILE_NAME __FILE__
program main
    use mod_global
    use mod_crystal
    use mod_lattice,     only : constructFCCLattice,constructBCCLattice,constructHCPLattice
    use mod_equilibrium, only : equilibriumCrystal
    use mod_atomistic,   only : Atomistic, constructAtomistic, destroyAtomistic
    use mod_kim,         only : setup_kim, allocate_kim, initialize_kim, getCutoff, destroy_kim
    use KIM_API_F03
    implicit none

    integer(c_int), parameter :: numCellsPerSide  = 8
    character(len=KIM_KEY_STRING_LENGTH) :: testname,modelname

    type(c_ptr)      :: pkim
    integer(c_int)   :: idum,ier

    ! Crystal related variables
    type(Crystal)    :: crystalType, crystalTypeRelaxed
    integer(c_int)   :: numBasisAtoms,numSpecies
    character(len=20):: crystalStructure
    character(len=2) :: species

    ! Atomistic system related variables
    type(Atomistic),pointer     :: atomisticSystem
    real(c_double)              :: rcutoff
    real(c_double), allocatable :: coords(:,:)
    integer(c_int)              :: numParticles,nbc


    read(*,*,err=100) testname
    read(*,*,err=100) modelname
    read(*,*,err=100) species
    read(*,*,err=100) crystalStructure

    ! Perform the hand-shake and determine the boundary conditions
    call setup_kim(pkim,testname,modelname,nbc)

    select case(trim(crystalStructure))
    case('fcc')
        numBasisAtoms = 1
    case('bcc')
        numBasisAtoms = 1
    case('hcp')
        numBasisAtoms = 2
    case('diamond')
        numBasisAtoms = 2
    end select
    call allocateCrystal(crystalType,numBasisAtoms)

    numParticles = numBasisAtoms * (numCellsPerSide**(DIM))
    numSpecies = 1

    ! Allocate memory to the kim object
    call allocate_kim(pkim,numParticles,numSpecies)

    ! Initialize KIM
    call initialize_kim(pkim)

    rcutoff = getCutoff(pkim)

    ! Now define the initial crystal using a lattice spacing of rcutoff/2
    crystalType%numSpecies = numSpecies
    crystalType%basisAtomsSpecies(:) = species
    select case(trim(crystalStructure))
    case('fcc')
        crystalType%LatticeType = constructFCCLattice(rcutoff/2._cd)
        crystalType%basisAtomsCoords(:,1) = (/0._cd,0._cd,0._cd/)
    case('bcc')
        crystalType%LatticeType = constructBCCLattice(rcutoff/2._cd)
        crystalType%basisAtomsCoords(:,1) = (/0._cd,0._cd,0._cd/)
    case('hcp')
        crystalType%LatticeType = constructHCPLattice(rcutoff/2._cd,rcutoff/2._cd)
        crystalType%basisAtomsCoords(:,1) = (/0._cd,0._cd,0._cd/)
        crystalType%basisAtomsCoords(:,2) = (/2._cd/3._cd,1._cd/3._cd,0.5_cd/)*rcutoff/2._cd
    case('diamond')
        crystalType%LatticeType = constructFCCLattice(rcutoff/2._cd)
        crystalType%basisAtomsCoords(:,1) = (/0._cd,0._cd,0._cd/)
        crystalType%basisAtomsCoords(:,2) = (/0.25_cd,0.25_cd,0.25_cd/)*rcutoff/2._cd
    end select 

    ! Relax the crystal
    crystalTypeRelaxed = equilibriumCrystal(crystalType,crystalStructure,testname,modelname)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now prepare the atomistic system using the relaxed crystal

    ! allocate memory
    atomisticSystem = constructAtomistic(numParticles,nbc)

    ! Broadcast the atomistic system to the kim api
    call broadcastAtomistic_kim(pkim,atomisticSystem)


    ! Fill in the coordinates using the relaxed crystal
    call constructGlobalCrystal(numCellsPerSide,crystalTypeRelaxed,atomisticSystem%coords,atomisticSystem%particleSpecies)
    atomisticSystem%nbc = nbc
    atomisticSystem%numContributingParticles = 1
    atomisticSystem%boxSideLengths  = 1000._cd
    atomisticSystem%numSpecies = numSpecies
    atomisticSystem%cutoff = rcutoff
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Write results
    call write_results_kim

    ! print results to screen
    print '(80(''-''))'
    print '("This is Test             : ",A)', trim(testname)
    print '("Results for KIM Model    : ",A)', trim(modelname)
    print '("Crystal strucutre        : ",A)', trim(crystalStructure)
    print '("Species                  : ",A)', species

    call destroyAtomistic(atomisticSystem)
    call destroy_kim(pkim)


    100 stop

end program main

