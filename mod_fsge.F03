module mod_Fsge
    use mod_global
    implicit none

    type Fsge
        real(c_double) :: Cijkl(DIM,DIM,DIM,DIM),           &
                          Dijmkln(DIM,DIM,DIM,DIM,DIM,DIM), &
                          Tijm(DIM,DIM,DIM),                &
                          Eijkln(DIM,DIM,DIM,DIM,DIM)
        integer(c_int), allocatable :: IndComponentsC(:,:), &
                                       IndComponentsD(:,:)
    end type Fsge

    contains
        subroutine initializeFsge(fsgeType)
            implicit none
            type(Fsge),intent(out) :: fsgeType

            fsgeType%Cijkl   = 0._cd
            fsgeType%Dijmkln = 0._cd
            fsgeType%Tijm    = 0._cd
            fsgeType%Eijkln  = 0._cd

        end subroutine initializeFsge

        subroutine allocateIndComponents(fsgeType,crystalStructure)
            implicit none
            type(Fsge),intent(inout) :: fsgeType
            character(len=*),intent(in) :: crystalStructure

            integer :: i,j,m

            if (trim(crystalStructure) .eq. 'fcc' .or. &
                trim(crystalStructure) .eq. 'bcc') then
                allocate(fsgeType%IndComponentsC(3,2))
                fsgeType%IndComponentsC(1,:) = (/1,1/)
                fsgeType%IndComponentsC(2,:) = (/1,2/)
                fsgeType%IndComponentsC(3,:) = (/4,4/)
                allocate(fsgeType%IndComponentsD(11,2))
                fsgeType%IndComponentsD(1,:)  = (/1,1/); fsgeType%IndComponentsD(2,:)  = (/1,2/)
                fsgeType%IndComponentsD(3,:)  = (/1,3/); fsgeType%IndComponentsD(4,:)  = (/2,2/)
                fsgeType%IndComponentsD(5,:)  = (/2,3/); fsgeType%IndComponentsD(6,:)  = (/2,4/)
                fsgeType%IndComponentsD(7,:)  = (/2,5/); fsgeType%IndComponentsD(8,:)  = (/3,3/)
                fsgeType%IndComponentsD(9,:)  = (/3,5/); fsgeType%IndComponentsD(10,:) = (/16,16/)
                fsgeType%IndComponentsD(11,:) = (/16,17/)
            elseif (trim(crystalStructure) .eq. 'hex') then
                allocate(fsgeType%IndComponentsC(5,2))
                fsgeType%IndComponentsC(1,:) = (/1,1/)
                fsgeType%IndComponentsC(2,:) = (/1,2/)
                fsgeType%IndComponentsC(3,:) = (/1,3/)
                fsgeType%IndComponentsC(4,:) = (/3,3/)
                fsgeType%IndComponentsC(5,:) = (/5,5/)
                allocate(fsgeType%IndComponentsD(22,2))
                fsgeType%IndComponentsD(1,:)  = (/1,1/); fsgeType%IndComponentsD(2,:)  = (/6,6/)
                fsgeType%IndComponentsD(3,:)  = (/6,7/); fsgeType%IndComponentsD(4,:)  = (/6,8/)
                fsgeType%IndComponentsD(5,:)  = (/6,9/); fsgeType%IndComponentsD(6,:)  = (/6,10/)
                fsgeType%IndComponentsD(7,:)  = (/7,7/); fsgeType%IndComponentsD(8,:)  = (/8,9/)
                fsgeType%IndComponentsD(9,:)  = (/8,10/); fsgeType%IndComponentsD(10,:) = (/9,9/)
                fsgeType%IndComponentsD(11,:) = (/9,10/); fsgeType%IndComponentsD(12,:) = (/10,10/)
                fsgeType%IndComponentsD(13,:) = (/11,11/); fsgeType%IndComponentsD(14,:) = (/11,12/)
                fsgeType%IndComponentsD(15,:) = (/11,13/); fsgeType%IndComponentsD(16,:) = (/12,12/)
                fsgeType%IndComponentsD(17,:) = (/12,13/); fsgeType%IndComponentsD(18,:) = (/13,13/)
                fsgeType%IndComponentsD(19,:) = (/16,16/); fsgeType%IndComponentsD(20,:) = (/16,17/)
                fsgeType%IndComponentsD(21,:) = (/17,17/); fsgeType%IndComponentsD(22,:) = (/17,18/)
            else
                allocate(fsgeType%IndComponentsC(21,2))
                m=1
                do i=1,6
                    do j=i,6
                        fsgeType%IndComponentsC(m,:) = (/i,j/)
                        m=m+1
                    enddo
                enddo

                allocate(fsgeType%IndComponentsD(171,2))
                m=1
                do i=1,18
                    do j=i,18
                        fsgeType%IndComponentsD(m,:) = (/i,j/)
                        m=m+1
                    enddo
                enddo
            endif
        end subroutine allocateIndComponents

        function voigtStrain(i)
        implicit none

        integer(c_int)             :: voigtStrain(2)
        integer(c_int), intent(in) :: i

        select case (i)
        case (1)
            voigtStrain = (/1,1/)
        case (2)
            voigtStrain = (/2,2/)
        case (3)
            voigtStrain = (/3,3/)
        case (4)
            voigtStrain = (/2,3/)
        case (5)
            voigtStrain = (/1,3/)
        case (6)
            voigtStrain = (/1,2/)
        end select

        end function voigtStrain

        function voigtStrainGradient(i)
        implicit none

        integer(c_int)             :: voigtStrainGradient(3)
        integer(c_int), intent(in) :: i

        select case (i)
        case (1)
            voigtStrainGradient = (/1,1,1/)
        case (2)
            voigtStrainGradient = (/2,2,1/)
        case (3)
            voigtStrainGradient = (/1,2,2/)
        case (4)
            voigtStrainGradient = (/3,3,1/)
        case (5)
            voigtStrainGradient = (/1,3,3/)
        case (6)
            voigtStrainGradient = (/2,2,2/)
        case (7)
            voigtStrainGradient = (/1,1,2/)
        case (8)
            voigtStrainGradient = (/1,2,1/)
        case (9)
            voigtStrainGradient = (/3,3,2/)
        case (10)
            voigtStrainGradient = (/2,3,3/)
        case (11)
            voigtStrainGradient = (/3,3,3/)
        case (12)
            voigtStrainGradient = (/1,1,3/)
        case (13)
            voigtStrainGradient = (/1,3,1/)
        case (14)
            voigtStrainGradient = (/2,2,3/)
        case (15)
            voigtStrainGradient = (/2,3,2/)
        case (16)
            voigtStrainGradient = (/1,2,3/)
        case (17)
            voigtStrainGradient = (/1,3,2/)
        case (18)
            voigtStrainGradient = (/2,3,1/)
        end select
        end function voigtStrainGradient

        subroutine writeFsgeMaterial(fsgeType)
            implicit none
            ! Passed variables
            type(Fsge), intent(in) :: fsgeType
            ! Local variables
            integer(c_int)   :: cpntsC(4), cpntsD(6),cpntsT(3),cpntsE(5)
            integer(c_int)   :: cunit,dunit,tunit,eunit,i,j,m,k,l,n
            integer,external :: GetFreeUnit
            character(len=1) :: cr

            cr = char(10)

            cunit = GetFreeUnit()
            open(unit=cunit,file='output/Cijkl.out',status='replace',position='append')
            do i=1,size(fsgeType%IndComponentsC,1)
                cpntsC(1:2) = voigtStrain(fsgeType%IndComponentsC(i,1))
                cpntsC(3:4) = voigtStrain(fsgeType%IndComponentsC(i,2))
                write(cunit,*) cpntsC(1),cpntsC(2),cpntsC(3),cpntsC(4),fsgeType%Cijkl(cpntsC(1),cpntsC(2),cpntsC(3),cpntsC(4))
            enddo
            close(cunit)

            cunit = GetFreeUnit()
            open(unit=cunit,file='output/Cijkl_ext.out',status='replace',position='append')
            do i=1,DIM
                do j=1,DIM
                    do k=1,DIM
                        do l=1,DIM
                            write(cunit,*) i,j,k,l,fsgeType%Cijkl(i,j,k,l)
                        enddo
                    enddo
                enddo
            enddo
            close(cunit)

            dunit = GetFreeUnit()
            open(unit=dunit,file='output/Dijmkln.out',status='replace',position='append')
            do i=1,size(fsgeType%IndComponentsD,1)
                cpntsD(1:3) = voigtStrainGradient(fsgeType%IndComponentsD(i,1))
                cpntsD(4:6) = voigtStrainGradient(fsgeType%IndComponentsD(i,2))
                write(dunit,*) cpntsD(1),cpntsD(2),cpntsD(3),cpntsD(4),cpntsD(5),cpntsD(6), &
                        fsgeType%Dijmkln(cpntsD(1),cpntsD(2),cpntsD(3),cpntsD(4),cpntsD(5),cpntsD(6))
            enddo
            close(dunit)

            dunit = GetFreeUnit()
            open(unit=dunit,file='output/Dijmkln_ext.out',status='replace',position='append')
            do i=1,18
                cpntsD(1:3) = voigtStrainGradient(i)
                do j=i,18
                    cpntsD(4:6) = voigtStrainGradient(j)
                    write(dunit,*) fsgeType%Dijmkln(cpntsD(1),cpntsD(2),cpntsD(3),cpntsD(4),cpntsD(5),cpntsD(6))
                enddo
            enddo
            !do i=1,DIM
            !    do j=1,DIM
            !        do m=1,DIM
            !            do k=1,DIM
            !                do l=1,DIM
            !                    do n=1,DIM
            !                        write(dunit,*) i,j,m,k,l,n,fsgeType%Dijmkln(i,j,m,k,l,n)
            !                    enddo
            !                enddo
            !            enddo
            !        enddo
            !    enddo
            !enddo
            close(dunit)

            tunit = GetFreeUnit()
            open(unit=tunit,file='output/Tijm.out',status='replace',position='append')
            do i=1,18
                cpntsT = voigtStrainGradient(i)
                write(tunit,*) fsgeType%Tijm(cpntsT(1),cpntsT(2),cpntsT(3))
            enddo
            close(tunit)

            eunit = GetFreeUnit()
            open(unit=eunit,file='output/Eijkln.out',status='replace',position='append')
            do i=1,6
                cpntsE(1:2) = voigtStrain(i)
                do j=1,18
                    cpntsE(3:5) = voigtStrainGradient(j)
                    write(eunit,*) fsgeType%Eijkln(cpntsE(1),cpntsE(2),cpntsE(3),cpntsE(4),cpntsE(5))
                enddo
            enddo
            close(eunit)

        end subroutine writeFsgeMaterial

        subroutine destroyFsgeMaterial(fsgeType)
            implicit none
            type(Fsge) :: fsgeType

            if (allocated(fsgeType%IndComponentsC)) deallocate(fsgeType%IndComponentsC)
            if (allocated(fsgeType%IndComponentsD)) deallocate(fsgeType%IndComponentsD)
        end subroutine destroyFsgeMaterial

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

end module mod_fsge
