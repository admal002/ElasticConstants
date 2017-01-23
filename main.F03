#include "KIM_API_status.h"
#define THIS_FILE_NAME __FILE__
program main
    use KIM_API_F03
    use mod_global
    use mod_crystal
    use mod_lattice
    use mod_equilibrium, only : equilibriumCrystal
    use mod_atomistic,   only : Atomistic, constructAtomistic, destroyAtomistic,constructNeighList,pA
    use mod_kim,         only : setup_kim, allocate_kim, initialize_kim, getCutoff_kim, &
                                destroy_kim, write_results_kim, broadcastAtomistic_kim,&
                                getSpeciesCode_kim,compute_kim,setCompute_kim,WriteDescriptorFile
    use mod_fsge,        only : allocateIndComponents,writeFsgeMaterial
    use KIM_API_F03
    implicit none

    integer(c_int), parameter :: numCellsPerSide  = 8
    character(len=KIM_KEY_STRING_LENGTH) :: testname,modelname

    type(c_ptr)      :: pkim

    ! Crystal related variables
    type(Crystal)    :: crystalType, crystalTypeRelaxed
    integer(c_int)   :: numBasisAtoms,numSpecies,i
    character(len=20):: crystalStructure
    character(len=2) :: species

    ! Atomistic system related variables
    type(Atomistic),target      :: atomisticSystem
    type(Atomistic),pointer     :: patomisticSystem
    real(c_double)              :: rcutoff,coord(DIM),latVol
    integer(c_int)              :: numParticles,nbc,speciesCode


    read(*,*,err=100) testname
    read(*,*,err=100) modelname
    read(*,*,err=100) species
    read(*,*,err=100) crystalStructure
    read(*,*,err=100) analytical

    print '(80(''-''))'
    print '("This is Test                : ",A)', trim(testname)
    print '("Results for KIM Model       : ",A)', trim(modelname)
    print '("Crystal strucutre           : ",A)', trim(crystalStructure)
    print '("Species                     : ",A)', species

    ! Initialize the KIM object
    call WriteDescriptorFile(pkim,species,modelname)

    ! Perform the hand-shake and determine the boundary conditions and process
    ! flags
    call setup_kim(pkim,nbc)

    select case(trim(crystalStructure))
    case('fcc')
        numBasisAtoms = 1
    case('bcc')
        numBasisAtoms = 1
    case('hex')
        numBasisAtoms = 1
    case('hcp')
        numBasisAtoms = 2
    case('diamond')
        numBasisAtoms = 2
    case('orth')
        numBasisAtoms = 2
    end select
    call allocateCrystal(crystalType,numBasisAtoms)

    numParticles = numBasisAtoms * (numCellsPerSide**(DIM))
    numSpecies = 1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate memory and initialize the kim api object for the sole purpose of
    ! getting the cutoff and species codes.

    call allocate_kim(pkim,numParticles,numSpecies)
    call initialize_kim(pkim,rcutoff)
    rcutoff = getCutoff_kim(pkim)
    speciesCode = getSpeciesCode_kim(pkim,species)
    call destroy_kim(pkim)

    !!!!!!!!!!!! We are done with the pkim object for now !!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Prepare a crystal of intial size rcutoff/2, and find the relaxed size

    crystalType%numSpecies = numSpecies
    crystalType%basisAtomsSpecies(:) = species
    crystalType%basisAtomsSpeciesCode(:) = speciesCode
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
        crystalType%basisAtomsCoords(:,2) =                     &
            (2._cd/3._cd)*crystalType%LatticeType%latVec(:,1) + &
            (1._cd/3._cd)*crystalType%LatticeType%latVec(:,2) + &
            (1._cd/2._cd)*crystalType%LatticeType%latVec(:,3)
    case('hex')
        crystalType%LatticeType = constructHCPLattice(rcutoff/3._cd,rcutoff/3._cd)
        crystalType%basisAtomsCoords(:,1) = (/0._cd,0._cd,0._cd/)
    case('diamond')
        crystalType%LatticeType = constructFCCLattice(rcutoff/2._cd)
        crystalType%basisAtomsCoords(:,1) = (/0._cd,0._cd,0._cd/)
        crystalType%basisAtomsCoords(:,2) = (/0.25_cd,0.25_cd,0.25_cd/)*rcutoff/2._cd
    case('orth')
        crystalType%LatticeType = constructOrthLattice(rcutoff*0.4_cd,rcutoff*0.6_cd,rcutoff*0.4_cd)
        ! B1   =  - y*A1 + y*A2 + 0.25*A3     =  + y*b*e2 + 0.25*c e3 
        ! B2   =  + y*A1 - y*A2 - 0.25*A3     =  - y*b*e2 - 0.25*c e3 
        ! Taking y=0.1
        crystalType%basisAtomsCoords(:,1) = (/0._cd, 0.1_cd*0.6_cd, 0.25_cd*0.4_cd/) * rcutoff
        crystalType%basisAtomsCoords(:,2) = (/0._cd,-0.1_cd*0.6_cd,-0.25_cd*0.4_cd/) * rcutoff
    end select 

    crystalTypeRelaxed = equilibriumCrystal(crystalType,crystalStructure,species,modelname)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now prepare the atomistic system using the relaxed crystal and broadcast
    ! it

    ! allocate memory
    call constructAtomistic(atomisticSystem,numParticles,nbc)
    patomisticSystem => atomisticSystem

    call WriteDescriptorFile(pkim,species,modelname)
    call setup_kim(pkim,nbc)
    call initialize_kim(pkim,patomisticSystem%cutoff)

    ! Fill in the coordinates using the relaxed crystal
    call constructglobalcrystal(numcellsperside,crystalTypeRelaxed,atomisticsystem%coords,atomisticsystem%particlespecies)
    atomisticSystem%nbc = nbc
    atomisticSystem%numContributingParticles = 1
    atomisticSystem%boxSideLengths  = 1000._cd
    atomisticSystem%numSpecies = numSpecies
    if (analytical .eq. 1) then 
        atomisticSystem%processdEdrFlag   = 1
        atomisticSystem%processd2Edr2Flag = 1
    else
        atomisticSystem%processdEdrFlag   = 0
        atomisticSystem%processd2Edr2Flag = 0
    endif

    ! Broadcast the atomistic system
    call broadcastAtomistic_kim(pkim,patomisticSystem)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! We are now ready to compute

    call allocateIndComponents(atomisticSystem%fsgeType,crystalStructure)
    call allocateIndComponents(atomisticSystem%fsgeTypeNumerical,crystalStructure)

    do i=1,numBasisAtoms
        coord = atomisticSystem%coords(:,1)
        atomisticSystem%coords(:,1) = atomisticSystem%coords(:,i)
        atomisticSystem%coords(:,i) = coord

        call constructNeighList(atomisticSystem)

        if (analytical .eq. 1) then
            pA => atomisticSystem
            call compute_kim(pkim)
            nullify(pA)
        else 
            call numerical_check_Cijkl(pkim,atomisticSystem)
            call numerical_check_Dijmkln(pkim,atomisticSystem)
        endif

        coord = atomisticSystem%coords(:,1)
        atomisticSystem%coords(:,1) = atomisticSystem%coords(:,i)
        atomisticSystem%coords(:,i) = coord
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Write results

    latVol = crystalTypeRelaxed%LatticeType%latVol
    atomisticSystem%fsgeType%Cijkl   = atomisticSystem%fsgeType%Cijkl/latVol
    atomisticSystem%fsgeType%Dijmkln = atomisticSystem%fsgeType%Dijmkln/latVol
    atomisticSystem%fsgeType%Tijm    = atomisticSystem%fsgeType%Tijm/latVol
    atomisticSystem%fsgeType%Eijkln  = atomisticSystem%fsgeType%Eijkln/latVol

    atomisticSystem%fsgeTypeNumerical%Cijkl   = atomisticSystem%fsgeTypeNumerical%Cijkl/latVol
    atomisticSystem%fsgeTypeNumerical%Dijmkln = atomisticSystem%fsgeTypeNumerical%Dijmkln/latVol
    atomisticSystem%fsgeTypeNumerical%Tijm    = atomisticSystem%fsgeTypeNumerical%Tijm/latVol
    atomisticSystem%fsgeTypeNumerical%Eijkln  = atomisticSystem%fsgeTypeNumerical%Eijkln/latVol

    if (analytical .eq. 1) then
        call writeFsgeMaterial(atomisticSystem%fsgeType)
    else
        call writeFsgeMaterial(atomisticSystem%fsgeTypeNumerical)
    endif

    if (analytical .eq. 1) then
        call write_results_kim(species,crystalStructure,crystalTypeRelaxed,atomisticSystem%fsgeType)
    else
        call write_results_kim(species,crystalStructure,crystalTypeRelaxed,atomisticSystem%fsgeTypeNumerical)
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! We are done. Deallocate memory

    call destroyCrystal(crystalType)
    call destroyCrystal(crystalTypeRelaxed)
    call destroyAtomistic(atomisticSystem)
    call destroy_kim(pkim)
    nullify(patomisticSystem)


    100 stop

end program main

subroutine numerical_check_Cijkl(pkim,aS)
    use mod_global
    use mod_kim, only : setCompute_kim,compute_kim
    use mod_fsge, only : voigtStrain
    use mod_atomistic, only : Atomistic, constructNeighList
    implicit none

    ! Passed variables
    type(c_ptr),   intent(in)     :: pkim
    type(Atomistic),intent(inout) :: aS

    ! Local variables
    real(c_double), allocatable :: coords_old(:,:) 
    real(c_double)              :: energy_plus,energy_neg, &
                                   eps,gradU(DIM,DIM)
    integer(c_int) :: partIndex,xc,yc,p,q,r,s,t,u,indx(2),nParts

    ! Set the process flags in the kim api to 0, as we do not need them
    if (aS%processdEdrFlag .eq. 1)   call setCompute_kim(pkim,'process_dEdr',0)
    if (aS%processd2Edr2Flag .eq. 1) call setCompute_kim(pkim,'process_d2Edr2',0)

    nParts = aS%numParticles
    allocate(coords_old(DIM,nParts))

    ! Initialize
    coords_old      = aS%coords
    eps             = 5.e-3


    do u=1,size(aS%fsgeTypeNumerical%IndComponentsC,1)

        energy_plus = 0._cd
        energy_neg  = 0._cd

        indx = voigtStrain(aS%fsgeTypeNumerical%IndComponentsC(u,1))
        p = indx(1); q = indx(2)
        indx = voigtStrain(aS%fsgeTypeNumerical%IndComponentsC(u,2))
        r = indx(1); s = indx(2)

        do t=1,4
            aS%coords = coords_old
            gradU  = 0._cd
            select case(t)
            case(1)
                if (p .eq. q) then
                    gradU(p,q) = eps
                else
                    gradU(p,q) = eps/2._cd
                    gradU(q,p) = eps/2._cd
                endif
                if (r .eq. s) then
                    gradU(r,s) = eps
                else
                    gradU(r,s) = eps/2._cd
                    gradU(s,r) = eps/2._cd
                endif
            case(2)
                if (p .eq. q) then
                    gradU(p,q) = eps
                else
                    gradU(p,q) = eps/2._cd
                    gradU(q,p) = eps/2._cd
                endif
                if (r .eq. s) then
                    gradU(r,s) = -eps
                else
                    gradU(r,s) = -eps/2._cd
                    gradU(s,r) = -eps/2._cd
                endif
                if (p .eq. r .and. q .eq. s) then
                    gradU(p,q) = 0._cd
                    gradU(q,p) = 0._cd
                endif
            case(3)
                if (p .eq. q) then
                    gradU(p,q) = -eps
                else
                    gradU(p,q) = -eps/2._cd
                    gradU(q,p) = -eps/2._Cd
                endif
                if (r .eq. s) then
                    gradU(r,s) = eps
                else
                    gradU(r,s) = eps/2._cd
                    gradU(s,r) = eps/2._cd
                endif
                if (p .eq. r .and. q .eq. s) then 
                    gradU(p,q) = 0._cd
                    gradU(q,p) = 0._cd
                endif
            case(4)
                if (p .eq. q) then
                    gradU(p,q) = -eps
                else
                    gradU(p,q) = -eps/2._cd
                    gradU(q,p) = -eps/2._cd
                endif
                if (r .eq. s) then
                    gradU(r,s) = -eps
                else
                    gradU(r,s) = -eps/2._cd
                    gradU(s,r) = -eps/2._cd
                endif
            endselect

            ! Deform
            do partIndex = 1,nParts
                do xc=1,3
                    do yc=1,3
                        aS%coords(xc,partIndex) = aS%coords(xc,partIndex) + gradU(xc,yc) * coords_old(yc,partIndex)
                    enddo
                enddo
            enddo

            call constructNeighList(aS)
            call compute_kim(pkim)

            if (t .eq. 1 .or. t .eq. 4) then
                energy_plus = energy_plus + aS%energy
            elseif (t .eq. 2 .or. t .eq. 3) then
                energy_neg = energy_neg + aS%energy
            endif
        enddo


        if ( p .eq. r .and. q .eq. s) then
            aS%fsgeTypeNumerical%Cijkl(p,q,r,s) = aS%fsgeTypeNumerical%Cijkl(p,q,r,s) + (energy_plus - energy_neg)/(eps**2)
        else
            aS%fsgeTypeNumerical%Cijkl(p,q,r,s) = aS%fsgeTypeNumerical%Cijkl(p,q,r,s) + (energy_plus - energy_neg)/(4._cd * eps**2)
        endif
    enddo

aS%coords = coords_old
deallocate(coords_old)

end subroutine numerical_check_Cijkl

subroutine numerical_check_Dijmkln(pkim,aS)
    use mod_global
    use mod_kim, only : setCompute_kim,compute_kim
    use mod_fsge, only : voigtStrainGradient
    use mod_atomistic, only : Atomistic, constructNeighList
    implicit none
    
    ! Passed variables
    type(c_ptr),   intent(in)     :: pkim
    type(Atomistic),intent(inout) :: aS

    ! Local variables
    real(c_double), allocatable            :: coords_original(:,:)
    integer(c_int)                         :: p,q,r,s,i,j,m,k,l,n,indx(3),nParts
    real(c_double)                         :: energy_plus,energy_neg,eps
    real(c_double), dimension(DIM,DIM,DIM) :: Cijk, Eijk

    ! Set the process flags in the kim api to 0, as we do not need them
    if (aS%processdEdrFlag .eq. 1)   call setCompute_kim(pkim,'process_dEdr',0)
    if (aS%processd2Edr2Flag .eq. 1) call setCompute_kim(pkim,'process_d2Edr2',0)

    nParts = aS%numParticles
    allocate(coords_original(DIM,nParts))

    coords_original = aS%coords
    eps             = 5.e-3

    do q=1,size(aS%fsgeTypeNumerical%IndComponentsD,1)
        energy_plus = 0._cd
        energy_neg  = 0._cd

        indx = voigtStrainGradient(aS%fsgeTypeNumerical%IndComponentsD(q,1))
        i = indx(1); j = indx(2); m = indx(3)
        indx = voigtStrainGradient(aS%fsgeTypeNumerical%IndComponentsD(q,2))
        k = indx(1); l = indx(2); n = indx(3)

        do p=1,4
            aS%coords = coords_original
            Eijk = 0._cd
            Cijk = 0._cd

            select case(p)
                case(1)
                    if (i .eq. j) then
                        Eijk(i,j,m) = eps
                    else
                        Eijk(i,j,m) = eps/2._cd
                        Eijk(j,i,m) = eps/2._cd
                    endif
                    if (k .eq. l) then
                        Eijk(k,l,n) = eps
                    else
                        Eijk(k,l,n) = eps/2._cd
                        Eijk(l,k,n) = eps/2._cd
                    endif
                case(2)
                    if (i .eq. j) then
                        Eijk(i,j,m) = eps
                    else
                        Eijk(i,j,m) = eps/2._cd
                        Eijk(j,i,m) = eps/2._cd
                    endif
                    if (k .eq. l) then
                        Eijk(k,l,n) = -eps
                    else
                        Eijk(k,l,n) = -eps/2._cd
                        Eijk(l,k,n) = -eps/2._cd
                    endif
                    if ( i .eq. k .and. j .eq. l .and. m .eq. n) then
                        Eijk(i,j,m) = 0._cd
                        Eijk(j,i,m) = 0._cd
                    endif
                case(3)
                    if (i .eq. j) then
                        Eijk(i,j,m) = -eps
                    else
                        Eijk(i,j,m) = -eps/2._cd
                        Eijk(j,i,m) = -eps/2._cd
                    endif
                    if (k .eq. l) then
                        Eijk(k,l,n) = eps
                    else
                        Eijk(k,l,n) = eps/2._cd
                        Eijk(l,k,n) = eps/2._cd
                    endif
                    if ( i .eq. k .and. j .eq. l .and. m .eq. n) then
                        Eijk(i,j,m) = 0._cd
                        Eijk(j,i,m) = 0._cd
                    endif
                case(4)
                    if (i .eq. j) then
                        Eijk(i,j,m) = -eps
                    else
                        Eijk(i,j,m) = -eps/2._cd
                        Eijk(j,i,m) = -eps/2._cd
                    endif
                    if (k .eq. l) then
                        Eijk(k,l,n) = -eps
                    else
                        Eijk(k,l,n) = -eps/2._cd
                        Eijk(l,k,n) = -eps/2._cd
                    endif
            end select

            call Calculate_Cijk(Eijk,Cijk)
        
            do r=1,nParts
                do s=1,DIM
                   aS%coords(s,r)   = aS%coords(s,r) + 0.5_cd*Cijk(s,1,1) * coords_original(1,r) * coords_original(1,r) +   &
                                                       0.5_cd*Cijk(s,1,2) * coords_original(1,r) * coords_original(2,r) +   &
                                                       0.5_cd*Cijk(s,1,3) * coords_original(1,r) * coords_original(3,r) +   &
                                                       0.5_cd*Cijk(s,2,1) * coords_original(2,r) * coords_original(1,r) +   &
                                                       0.5_cd*Cijk(s,2,2) * coords_original(2,r) * coords_original(2,r) +   &
                                                       0.5_cd*Cijk(s,2,3) * coords_original(2,r) * coords_original(3,r) +   &
                                                       0.5_cd*Cijk(s,3,1) * coords_original(3,r) * coords_original(1,r) +   &
                                                       0.5_cd*Cijk(s,3,2) * coords_original(3,r) * coords_original(2,r) +   &
                                                       0.5_cd*Cijk(s,3,3) * coords_original(3,r) * coords_original(3,r)
                enddo
            enddo

            call constructNeighList(aS)
            call compute_kim(pkim)

            if (p .eq. 1 .or. p .eq. 4) then
                energy_plus = energy_plus + aS%energy
            elseif (p .eq. 2 .or. p .eq. 3) then
                energy_neg = energy_neg + aS%energy
            endif

        enddo

        if ( i .eq. k .and. j .eq. l .and. m .eq. n) then
            aS%fsgeTypeNumerical%Dijmkln(i,j,m,k,l,n) = &
                aS%fsgeTypeNumerical%Dijmkln(i,j,m,k,l,n) + ((energy_plus - energy_neg)/(eps**2))
        else
            aS%fsgeTypeNumerical%Dijmkln(i,j,m,k,l,n) = &
                aS%fsgeTypeNumerical%Dijmkln(i,j,m,k,l,n) + ((energy_plus - energy_neg)/(4._cd * eps**2))
        endif
    enddo

    aS%coords = coords_original
    deallocate(coords_original)

end subroutine numerical_check_Dijmkln

subroutine numerical_check_Tijm(pkim,aS)
    use mod_global
    use mod_kim, only : setCompute_kim,compute_kim
    use mod_fsge, only : voigtStrainGradient
    use mod_atomistic, only : Atomistic, constructNeighList
    implicit none
    
    ! Passed variables
    type(c_ptr),   intent(in)     :: pkim
    type(Atomistic),intent(inout) :: aS

    ! Local variables
    real(c_double), allocatable            :: coords_original(:,:)
    integer(c_int)                         :: p,q,r,s,i,j,m,indx(3),nParts
    real(c_double)                         :: energy_plus,energy_neg,eps
    real(c_double), dimension(DIM,DIM,DIM) :: Cijk, Eijk

    ! Set the process flags in the kim api to 0, as we do not need them
    if (aS%processdEdrFlag .eq. 1)   call setCompute_kim(pkim,'process_dEdr',0)
    if (aS%processd2Edr2Flag .eq. 1) call setCompute_kim(pkim,'process_d2Edr2',0)

    nParts = aS%numParticles

    allocate(coords_original(DIM,nParts))

    coords_original = aS%coords
    eps             = 5.e-3

    do q=1,18
        energy_plus = 0._cd
        energy_neg  = 0._cd

        indx = voigtStrainGradient(q)
        i = indx(1); j = indx(2); m = indx(3)

        do p=1,2
            aS%coords = coords_original
            Eijk = 0._cd
            Cijk = 0._cd

            select case(p)
                case(1)
                    if (i .eq. j) then
                        Eijk(i,j,m) = eps
                    else
                        Eijk(i,j,m) = eps/2._cd
                        Eijk(j,i,m) = eps/2._cd
                    endif
                case(2)
                    if (i .eq. j) then
                        Eijk(i,j,m) = -eps
                    else
                        Eijk(i,j,m) = -eps/2._cd
                        Eijk(j,i,m) = -eps/2._cd
                    endif
            end select

            call Calculate_Cijk(Eijk,Cijk)
        
            do r=1,nParts
                do s=1,3
                   aS%coords(s,r)   = aS%coords(s,r) + 0.5d0*Cijk(s,1,1) * coords_original(1,r) * coords_original(1,r) +   &
                                                       0.5d0*Cijk(s,1,2) * coords_original(1,r) * coords_original(2,r) +   &
                                                       0.5d0*Cijk(s,1,3) * coords_original(1,r) * coords_original(3,r) +   &
                                                       0.5d0*Cijk(s,2,1) * coords_original(2,r) * coords_original(1,r) +   &
                                                       0.5d0*Cijk(s,2,2) * coords_original(2,r) * coords_original(2,r) +   &
                                                       0.5d0*Cijk(s,2,3) * coords_original(2,r) * coords_original(3,r) +   &
                                                       0.5d0*Cijk(s,3,1) * coords_original(3,r) * coords_original(1,r) +   &
                                                       0.5d0*Cijk(s,3,2) * coords_original(3,r) * coords_original(2,r) +   &
                                                       0.5d0*Cijk(s,3,3) * coords_original(3,r) * coords_original(3,r)
                enddo
            enddo

            call constructNeighList(aS)
            call compute_kim(pkim)

            if (p .eq. 1) then
                energy_plus = energy_plus + aS%energy
            elseif (p .eq. 2) then
                energy_neg = energy_neg + aS%energy
            endif

        enddo

        aS%fsgeTypeNumerical%Tijm(i,j,m) = aS%fsgeTypeNumerical%Tijm(i,j,m) + (energy_plus - energy_neg)/(2._cd*eps)
    enddo

    aS%coords = coords_original

    deallocate(coords_original)


end subroutine numerical_check_Tijm

subroutine Calculate_Cijk(Eijk,Cijk)
    use mod_global
    implicit none

    ! Passed variables
    real(c_double),intent(in) :: Eijk(DIM,DIM,DIM)
    ! Local variables
    integer :: i,j,k
    real(c_double),intent(out) :: Cijk(DIM,DIM,DIM)

    Cijk = 0._cd
    do i=1,DIM
        do j=1,DIM
            do k=1,DIM
                Cijk(i,j,k) = Eijk(i,j,k) + Eijk(i,k,j) - Eijk(k,j,i) 
            enddo
        enddo
    enddo
    return
end subroutine Calculate_Cijk

