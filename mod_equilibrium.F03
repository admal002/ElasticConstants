#include "KIM_API_status.h"
#define THIS_FILE_NAME __FILE__
module mod_equilibrium
    use mod_global
    use mod_crystal, only : Crystal,constructGlobalCrystal,scaleCrystal,allocateCrystal,destroyCrystal
    use mod_atomistic
    use mod_lattice, only : constructFCCLattice, constructOrthLattice, constructBCCLattice, constructHCPLattice
    use mod_kim

    integer(c_int), parameter :: numCellsPerSide  = 8

    type(c_ptr)   :: pkim

    type(Crystal)    :: crystalType
    character(len=20) :: crystalStructure

    type(Atomistic),target  :: atomisticSystem
    type(Atomistic),pointer :: patomisticSystem

    contains

        function equilibriumCrystal(cT,cS,species,modelname)
            implicit none

            ! Passed variables
            character(len=KIM_KEY_STRING_LENGTH) :: modelname
            character(len=2) :: species
            character(len=*) :: cS
            type(Crystal) :: cT
            type(Crystal) :: equilibriumCrystal

            ! Local variables
            integer(c_int) :: nbc,numBasisAtoms,numParticles
            integer(c_int) :: np,mp,ndim
            real(c_double) :: ftol
            real(c_double),allocatable :: p(:,:),p0(:),y(:)
            integer(c_int) :: iter,dimn
            
            select case(trim(cS))
            ! update the lattice spacing based on crystalStructure
            case('fcc' , 'bcc' , 'diamond')
                ndim=1
            case('hcp', 'hex')
                ndim=2
            case('orth')
                ndim=4
            case default
                print '("mod_equilibrium: Crystal structure not supported")'
            end select
            np = ndim
            mp = ndim+1
            ftol = epsilon(1._cd)
            allocate(p(mp,np))
            allocate(p0(np))
            allocate(y(mp))
            
            crystalStructure = cS

            ! Perform the hand-shake and determine the boundary conditions
            call WriteDescriptorFile(pkim,species,modelname)
            call setup_kim(pkim,nbc)

            numBasisAtoms = size(cT%basisAtomsCoords,2)
            numParticles = numBasisAtoms * (numCellsPerSide**(DIM))

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Now prepare the atomistic system

            call constructAtomistic(atomisticSystem,numParticles,nbc)
            patomisticSystem => atomisticSystem
            call initialize_kim(pkim,patomisticSystem%cutoff)

            ! Fill in the coordinates
            call constructGlobalCrystal(numCellsPerSide,cT,atomisticSystem%coords,  &
                                        atomisticSystem%particleSpecies)
            atomisticSystem%nbc = nbc
            atomisticSystem%numContributingParticles = size(cT%basisAtomsCoords,2)
            atomisticSystem%boxSideLengths  = 1000._cd
            atomisticSystem%numSpecies = cT%numSpecies
            if (analytical .eq. 1) then 
                atomisticSystem%processdEdrFlag   = 1
                atomisticSystem%processd2Edr2Flag = 1
            else
                atomisticSystem%processdEdrFlag   = 0
                atomisticSystem%processd2Edr2Flag = 0
            endif

            ! Broadcast the atomistic system
            call broadcastAtomistic_kim(pkim,patomisticSystem)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


            call allocateCrystal(crystalType,numBasisAtoms)
            crystalType = cT

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
            ! Start relaxation using the simplex method

            write(*,*) 'Relaxing the crystal to its equilibrium lattice constant'
            write(*,*) '--------------------------------------------------------'
            if(trim(cS) .eq. 'fcc' .or. trim(cS) .eq. 'bcc' .or. trim(cS) .eq.  'diamond') then
                write(*,'(A20,2x,A20)') 'Lattice spacing','energy'
            elseif (trim(cS) .eq. 'hcp' .or. trim(cS) .eq. 'hex') then
                write(*,'(A20,2x,A20,2x,A10)') 'Lattice spacing a','Lattice spacing c','energy'
            elseif (trim(cS) .eq. 'orth') then
                write(*,'(A20,2x,A20,2x,A20,2x,A20,2x,A10)') 'Lattice spacing a', &
                                               'Lattice spacing b', &
                                               'Lattice spacing c', &
                                               'Lattice shift y','energy'
            endif
            ! Set the process flags in the kim api to 0, as we do not need them
            if (atomisticSystem%processdEdrFlag .eq. 1)   call setCompute_kim(pkim,'process_dEdr',0)
            if (atomisticSystem%processd2Edr2Flag .eq. 1) call setCompute_kim(pkim,'process_d2Edr2',0)

            ! Calculate energy for initial crystal
            call constructNeighList(atomisticSystem)
            p(1,:) = atomisticSystem%cutoff/2._cd
            if(trim(cS) .eq. 'orth') then 
                p(1,1:3) = (/0.4_cd,0.6_cd,0.4_cd/)*atomisticSystem%cutoff
                p(1,4) = 0.1_cd
            endif
            p0 = p(1,:)
            pA => atomisticSystem
            call compute_kim(pkim)
            nullify(pA)
            y(1) = atomisticSystem%energy

            do dimn=1,ndim
                select case(trim(cS))
                case('fcc'  , 'bcc', 'diamond')
                    crystalType = scaleCrystal(cT,2._cd,2._cd,2._cd)
                    p(dimn+1,:) = 2._cd*p(1,:)
                case('hcp','hex')
                    select case(dimn)
                    case (1)
                        crystalType = scaleCrystal(cT,2._cd,2._cd,1._cd)
                        p(dimn+1,1) = 2._cd*p(1,1)
                        p(dimn+1,2) = p(1,2)
                    case(2)
                        crystalType = scaleCrystal(cT,1._cd,1._cd,2._cd)
                        p(dimn+1,1) = p(1,1)
                        p(dimn+1,2) = 2._cd*p(1,2)
                    end select
                case('orth')
                    select case(dimn)
                    case(1)
                        crystalType = scaleCrystal(cT,2._cd,1._cd,1._cd)
                        p(dimn+1,1) = 2._cd*p(1,1)
                        p(dimn+1,2) = p(1,2)
                        p(dimn+1,3) = p(1,3)
                        p(dimn+1,4) = p(1,4)
                    case(2)
                        crystalType = scaleCrystal(cT,1._cd,2._cd,1._cd)
                        p(dimn+1,1) = p(1,1)
                        p(dimn+1,2) = 2._cd*p(1,2)
                        p(dimn+1,3) = p(1,3)
                        p(dimn+1,4) = p(1,4)
                    case(3)
                        crystalType = scaleCrystal(cT,1._cd,1._cd,2._cd)
                        p(dimn+1,1) = p(1,1)
                        p(dimn+1,2) = p(1,2)
                        p(dimn+1,3) = 2._cd*p(1,3)
                        p(dimn+1,4) = p(1,4)
                    case(4)
                        crystalType = scaleCrystal(cT,1._cd,1._cd,1._cd)
                        crystalType%basisAtomsCoords(2,:) = 2._cd*crystalType%basisAtomsCoords(2,:)
                        p(dimn+1,1) = p(1,1)
                        p(dimn+1,2) = p(1,2)
                        p(dimn+1,3) = p(1,3)
                        p(dimn+1,4) = 2._cd*p(1,4)

                    end select
                end select

                call constructGlobalCrystal(numCellsPerSide,crystalType,atomisticSystem%coords,atomisticSystem%particleSpecies)
                call constructNeighList(atomisticSystem)
                pA => atomisticSystem
                call compute_kim(pkim)
                nullify(pA)
                y(dimn+1) = atomisticSystem%energy

            enddo

            call amoeba(p,y,mp,np,ndim,ftol,iter)
            print '("Equilibrium lattice spacing : ",3f15.10)', p(1,:)

            call allocateCrystal(equilibriumCrystal,numBasisAtoms)
            equilibriumCrystal = cT

            select case(trim(cS))
            case('fcc' , 'bcc' , 'diamond')
                equilibriumCrystal = scaleCrystal(cT,p(1,1)/p0(1),p(1,1)/p0(1),p(1,1)/p0(1))
            case('hcp','hex')
                equilibriumCrystal = scaleCrystal(cT,p(1,1)/p0(1),p(1,1)/p0(1),p(1,2)/p0(2))
            case('orth')
                equilibriumCrystal = scaleCrystal(cT,p(1,1)/p0(1),p(1,2)/p0(2),p(1,3)/p0(3))
                equilibriumCrystal%basisAtomsCoords(2,:) = equilibriumCrystal%basisAtomsCoords(2,:)*p(1,4)/p0(4)
            end select

            deallocate(y)
            deallocate(p)
            deallocate(p0)
            call destroy_kim(pkim)
            call destroyAtomistic(atomisticSystem)
            call destroyCrystal(crystalType)

        end function equilibriumCrystal

        function funk(latticeSpacing)
            implicit none
            real(c_double), intent(in) :: latticeSpacing(*)
            real(c_double)             :: funk

            real(c_double), parameter :: latticeSpacingMin = 1e-3,  &
                                         energyMax = 1e5

            if (trim(crystalStructure) .eq. 'fcc' .or.  &
                trim(crystalStructure) .eq. 'bcc' .or.  &
                trim(crystalStructure) .eq. 'diamond') then
                if (latticeSpacing(1) .lt. latticeSpacingMin) then
                    funk = energyMax;
                    return
                endif
            elseif (trim(crystalStructure) .eq. 'hcp' .or. trim(crystalStructure) .eq. 'hex') then
                if (min(latticeSpacing(1),latticeSpacing(2)) .lt. latticeSpacingMin) then
                    funk = energyMax;
                    return
                endif
            elseif (trim(crystalStructure) .eq. 'orth') then
                if (min(latticeSpacing(1),latticeSpacing(2),latticeSpacing(3),latticeSpacing(4)) .lt. latticeSpacingMin) then
                    funk = energyMax;
                    return
                endif
            endif

            select case(trim(crystalStructure))
            ! update the lattice spacing based on crystalStructure
            case('fcc')
                crystalType%LatticeType = constructFCCLattice(latticeSpacing(1))
                crystalType%basisAtomsCoords(:,1) = (/0._cd,0._cd,0._cd/)
            case('bcc')
                crystalType%LatticeType = constructBCCLattice(latticeSpacing(1))
                crystalType%basisAtomsCoords(:,1) = (/0._cd,0._cd,0._cd/)
            case('hcp')
                crystalType%LatticeType = constructHCPLattice(latticeSpacing(1),latticeSpacing(2))
                crystalType%basisAtomsCoords(:,1) = (/0._cd,0._cd,0._cd/)
                crystalType%basisAtomsCoords(:,2) = (/latticeSpacing(1)*0.5_cd,  &
                                            latticeSpacing(1)*(0.5_cd/sqrt(3._cd)),   &
                                            latticeSpacing(2)*0.5_cd/) 
            case('hex')
                crystalType%LatticeType = constructHCPLattice(latticeSpacing(1),latticeSpacing(2))
                crystalType%basisAtomsCoords(:,1) = (/0._cd,0._cd,0._cd/)
            case('diamond')
                crystalType%LatticeType = constructFCCLattice(latticeSpacing(1))
                crystalType%basisAtomsCoords(:,1) = (/0._cd,0._cd,0._cd/)
                crystalType%basisAtomsCoords(:,2) = (/0.25_cd,0.25_cd,0.25_cd/)*latticeSpacing(1)
            case('orth')
                crystalType%LatticeType = constructOrthLattice(latticeSpacing(1),   &
                                                               latticeSpacing(2),   &
                                                               latticeSpacing(3))
                ! B1   =  - y*A1 + y*A2 + 0.25*A3     =  + y*b*e2 + 0.25*c e3 
                ! B2   =  + y*A1 - y*A2 - 0.25*A3     =  - y*b*e2 - 0.25*c e3 
                crystalType%basisAtomsCoords(:,1) =  (/0._cd,latticeSpacing(4)*latticeSpacing(2),0.25_cd*latticeSpacing(3)/)
                crystalType%basisAtomsCoords(:,2) = -(/0._cd,latticeSpacing(4)*latticeSpacing(2),0.25_cd*latticeSpacing(3)/)
            end select 

            call constructGlobalCrystal(numCellsPerSide,                &
                                        crystalType,                    &
                                        atomisticSystem%coords,         &
                                        atomisticSystem%particleSpecies)
            call constructNeighList(atomisticSystem)

            ! This is actually a fix because I am unable to set the process
            ! functions as type-bound procedures to the Atomistic type
            ! atomisticSystem. Instead, I use the pointer pA to point to
            ! atomisticSystem to tell the "model compute" to build the elastic
            ! tensors of the atomisticSystem.

            pA => atomisticSystem
            call compute_kim(pkim)
            nullify(pA)

            if(trim(crystalStructure) .eq. 'fcc' .or.       &
               trim(crystalStructure) .eq. 'bcc' .or.       &
               trim(crystalStructure) .eq.  'diamond') then
                write(*,'(f20.5,2x,f20.5,2x,f20.10)') latticeSpacing(1:1),atomisticSystem%energy
            elseif (trim(crystalStructure) .eq. 'hcp' .or.  &
                    trim(crystalStructure) .eq. 'hex') then
                write(*,'(f20.5,2x,f20.5,2x,f20.10)') latticeSpacing(1:2),atomisticSystem%energy
            elseif (trim(crystalStructure) .eq. 'orth') then
                write(*,'(f20.5,2x,f20.5,2x,f20.10,2x,f20.10,2x,f20.10)') latticeSpacing(1:4),atomisticSystem%energy
            endif
            funk = atomisticSystem%energy
        end function funk

        SUBROUTINE amoeba(p,y,mp,np,ndim,ftol,iter) 
            implicit none
            INTEGER(c_int) iter,mp,ndim,np,NMAX,ITMAX
            REAL(c_double) ftol,p(mp,np),y(mp),TINY
            PARAMETER (NMAX=20,ITMAX=1000,TINY=1.e-10) 
            INTEGER(c_int) i,ihi,ilo,inhi,j,m,n
            REAL(c_double) rtol,sum,swap,ysave,ytry,psum(NMAX)

            iter=0
            1 do n=1,ndim 
                sum=0.
                do m=1,ndim+1
                    sum=sum+p(m,n) 
                enddo
                psum(n)=sum 
            enddo
            2 ilo=1
            if (y(1).gt.y(2)) then
                ihi=1
                inhi=2 
            else
                ihi=2
                inhi=1 
            endif
            do i=1,ndim+1 
                if(y(i).le.y(ilo)) ilo=i 
                if(y(i).gt.y(ihi)) then
                    inhi=ihi
                    ihi=i
                else if(y(i).gt.y(inhi)) then
                    if(i.ne.ihi) inhi=i
                endif 
            enddo
            rtol=2.*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
            if (rtol.lt.ftol) then 
                swap=y(1)
                y(1)=y(ilo) 
                y(ilo)=swap 
                do n=1,ndim
                    swap=p(1,n) 
                    p(1,n)=p(ilo,n) 
                    p(ilo,n)=swap
                enddo
                return 
            endif
            if (iter.ge.ITMAX) then
                print '("ITMAX exceeded in amoeba")'
                print '("Check convergence. If there seems to be convergence, then increase ITMAX")'
                stop
            endif
            iter=iter+2
            ytry=amotry(p,y,psum,mp,np,ndim,ihi,-1._cd) 
            if (ytry.le.y(ilo)) then
                ytry=amotry(p,y,psum,mp,np,ndim,ihi,2._cd)
            else if (ytry.ge.y(inhi)) then
                ysave=y(ihi)
                ytry=amotry(p,y,psum,mp,np,ndim,ihi,0.5_cd)
                if (ytry.ge.ysave) then 
                    do i=1,ndim+1
                        if(i.ne.ilo) then
                            do j=1,ndim 
                                psum(j)=0.5*(p(i,j)+p(ilo,j))
                                p(i,j)=psum(j)
                            enddo
                            y(i)=funk(psum) 
                        endif
                    enddo 
                    iter=iter+ndim
                    goto 1
                endif 
            else
                iter=iter-1 
            endif
            goto 2 
        END subroutine amoeba

        FUNCTION amotry(p,y,psum,mp,np,ndim,ihi,fac) 
            implicit none
            INTEGER(c_int) ihi,mp,ndim,np,NMAX
            REAL(c_double) amotry,fac,p(mp,np),psum(np),y(mp)
            PARAMETER (NMAX=20)
            INTEGER(c_int) j
            REAL(c_double) fac1,fac2,ytry,ptry(NMAX) 
            fac1=(1.-fac)/ndim 
            fac2=fac1-fac
            do j=1,ndim
                ptry(j)=psum(j)*fac1-p(ihi,j)*fac2 
            enddo
            ytry=funk(ptry)
            if (ytry.lt.y(ihi)) then
                y(ihi)=ytry 
                do j=1,ndim
                    psum(j)=psum(j)-p(ihi,j)+ptry(j)
                    p(ihi,j)=ptry(j) 
                enddo
            endif 
            amotry=ytry 
            return
        END function amotry

end module mod_equilibrium
