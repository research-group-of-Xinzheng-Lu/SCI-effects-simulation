!Main subroutines for MDOF input and matrix generation. by Y TIAN
subroutine MDOFInputd(mdofid)
	use MDOF_MyData
	implicit none
    Character(len=200) tempchar, StrucTypeName
    integer mdofid
    integer I,J,tempint
    integer MDOFmotemp1,MDOFmotemp2
    real*8 tempreal

    MDOFmotemp1=1+(mdofid+1)*10
    MDOFmotemp2=2+(mdofid+1)*10

    open(MDOFmotemp1,file=MDOFfinfo) 
	    read(MDOFmotemp1,*)  	
    do I=1, NumofBlgs 		
        read(MDOFmotemp1,"(3I10,2E15.7)") Blg(I)%ID, &
            & Blg(I)%StrucType, Blg(I)%NDOF, Blg(I)%Floor_h, Blg(I)%Area		
        Blg(I)%Height=Blg(I)%NDOF*Blg(I)%Floor_h        
            Blgdof(I)=Blg(I)%NDOF
            allocate(Blg(I)%props(Blg(I)%NDOF,10,1))
            allocate(Blg(I)%dval(Blgdof(I),4))
            allocate(Blg(I)%BlgM(Blgdof(I),Blgdof(I)))
            allocate(Blg(I)%BlgK(Blgdof(I),Blgdof(I)))
            allocate(Blg(I)%BlgC(Blgdof(I),Blgdof(I)))
            allocate(Blg(I)%BlgM_Inv(Blgdof(I),Blgdof(I)))
            allocate(Blg(I)%MDOFEt(Blgdof(I),2)) 
            allocate(Blg(I)%MDOFstatev(Blgdof(I),11,2))
            allocate(Blg(I)%MDOFspd(Blgdof(I),2))
            allocate(Blg(I)%MDOFyield(Blgdof(I),2))
            allocate(Blg(I)%MDOFIDR(Blgdof(I),2))
            allocate(Blg(I)%dMDOFIDR(Blgdof(I),2))
            allocate(Blg(I)%MDOFItF(Blgdof(I),2))
            allocate(Blg(I)%tempMDOFU0(Blgdof(I),2))
            allocate(Blg(I)%tempMDOFU1(Blgdof(I),2))
            allocate(Blg(I)%tempMDOFU(Blgdof(I),2))
            allocate(Blg(I)%BlgMaxIDR(Blgdof(I),2))
            allocate(Blg(I)%MDOFIDeath(Blgdof(I),2))
            allocate(Blg(I)%MDOFstate(Blg(I)%NDOF,2))
            allocate(Blg(I)%tempMDOFA1(Blg(I)%NDOF,2))
            Blg(I)%BlgM=0; Blg(I)%BlgK=0; Blg(I)%BlgC=0; Blg(I)%BlgM_Inv=0; 
            Blg(I)%MDOFEt=0; Blg(I)%MDOFstatev=0; Blg(I)%MDOFspd=0; Blg(I)%MDOFyield=0; 
            Blg(I)%MDOFIDR=0; Blg(I)%dMDOFIDR=0; Blg(I)%MDOFItF=0;
            Blg(I)%tempMDOFU0=0; Blg(I)%tempMDOFU1=0; Blg(I)%tempMDOFU=0; 
            Blg(I)%BlgMaxIDR=0; Blg(I)%MDOFIDeath=0; Blg(I)%MDOFstate=0
        do J=1, Blgdof(I)
		    read(MDOFmotemp1,"(10E15.7)")Blg(I)%props(J,1:10,1)
        end do
    end do
    do I=1, NumofBlgs
        read(MDOFmotemp1,*) 
        do J=1, Blgdof(I)
		    read(MDOFmotemp1,"(4E15.7)")Blg(I)%dval(J,1:4)
        end do 
    end do    
    read(MDOFmotemp1,*) 
    do I=1, NumofBlgs
		read(MDOFmotemp1,"(I10,5E15.7)")tempint,&
            & Blg(I)%T1, Blg(I)%T2,Blg(I)%Calpha,Blg(I)%Cbeta, Blg(I)%TN
    end do

	close(MDOFmotemp1)
    return
end subroutine MDOFInputd
!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!               
subroutine MDOFMatrx(dtsite)
	use MDOF_MyData
	implicit none
    integer I,J, tempint
	real*8 tempreal,damping_temp
    real*8,allocatable:: tempprops(:,:,:)
    real*8,intent(in):: dtsite
    do I=1, NumofBlgs
        Blg(I)%dt=Blg(I)%TN/(dtan(1.d0)*4.d0)
        Blg(I)%ndt=ceiling(dtsite/Blg(I)%dt)
        Blg(I)%dt=dtsite/Blg(I)%ndt
        Blg(I)%dt2=Blg(I)%dt*Blg(I)%dt
            !!!!!!!!!!!MassMatrix!!!!!!!!!!!!
            do J=1,Blg(I)%NDOF
                Blg(I)%BlgM(J,J)=MasspArea*Blg(I)%Area
            end do
            !!!!!!!!!!!MassMatrix!!!!!!!!!!!!
            !!!!!!!!!!!StiffMatrix!!!!!!!!!!!!
            do J=1,Blg(I)%NDOF
                Blg(I)%BlgK(J,J)=Blg(I)%BlgK(J,J)+Blg(I)%props(J,1,1) 
                if(J>1) then
                    Blg(I)%BlgK(J-1,J-1)=Blg(I)%BlgK(J-1,J-1)+Blg(I)%props(J,1,1)
                    Blg(I)%BlgK(J,J-1)=Blg(I)%BlgK(J,J-1)-Blg(I)%props(J,1,1)
                    Blg(I)%BlgK(J-1,J)=Blg(I)%BlgK(J-1,J)-Blg(I)%props(J,1,1)
                end if
            end do
            !!!!!!!!!!!StiffMatrix!!!!!!!!!!!!
            !!!!!!!!!!!DampingMatrix!!!!!!!!!!!!
            Blg(I)%BlgC(1:Blg(I)%NDOF,1:Blg(I)%NDOF)= &
            & Blg(I)%Calpha*Blg(I)%BlgM(1:Blg(I)%NDOF,1:Blg(I)%NDOF)+Blg(I)%Cbeta*Blg(I)%BlgK(1:Blg(I)%NDOF,1:Blg(I)%NDOF)
            !!!!!!!!!!!DampingMatrix!!!!!!!!!!!!
        deallocate(Blg(I)%BlgK)
	enddo
	return
end subroutine MDOFMatrx

!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MDOFDamageState(BlgMaxIDR,BlgNDOF,dof,BlgStrucType,damagetemp,BlgDValtmp)
	use MDOF_MyData
	integer*4:: BlgStrucType,tempii,BlgNDOF,dof
	real*8:: BlgMaxIDR(dof),BlgDValtmp(dof,4)
    integer*4:: damagetemp(BlgNDOF)
	do tempii=1,BlgNDOF
		if(BlgMaxIDR(tempii).LE.BlgDValtmp(tempii,1))then
			damagetemp(tempii)=0
		elseif((BlgMaxIDR(tempii).gt.BlgDValtmp(tempii,1)).and.(BlgMaxIDR(tempii).LE.BlgDValtmp(tempii,2)))then
			damagetemp(tempii)=1
		elseif((BlgMaxIDR(tempii).gt.BlgDValtmp(tempii,2)).and.(BlgMaxIDR(tempii).LE.BlgDValtmp(tempii,3)))then
			damagetemp(tempii)=2
		elseif((BlgMaxIDR(tempii).gt.BlgDValtmp(tempii,3)).and.(BlgMaxIDR(tempii).LE.BlgDValtmp(tempii,4)))then
			damagetemp(tempii)=3
		elseif(BlgMaxIDR(tempii).GT.BlgDValtmp(tempii,4))then
			damagetemp(tempii)=4
		endif
	enddo
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GenBFile
    use MDOF_MyData
    use speed_timeloop
	implicit none
    logical :: BRornot
    
    MDOFmonitor=701+mpi_id
    inquire(file="BRCtrl.txt",exist=BRornot)
    if(BRornot) then
    open(MDOFmonitor,file="BRCtrl.txt")
        read(MDOFmonitor,*)BRoutput(1:6) 
    close(MDOFmonitor)
    else
    BRoutput(1:4)=1
    BRoutput(5:6)=0
    end if

    Maxdof=0
    write(MDOFstr,*) mpi_id+1
    write(MDOFfinfo,*) trim(adjustl(mpi_file)),'/BuildingInfo',trim(adjustl(MDOFstr)),'.txt'
    MDOFfinfo=trim(adjustl(MDOFfinfo))
    write(MDOFftest,*) trim(adjustl(mpi_file)),'/mpitest',trim(adjustl(MDOFstr)),'.txt'
    MDOFftest=trim(adjustl(MDOFftest))
    write(MDOFftest0,*) trim(adjustl(mpi_file)),'/mpitest0.txt'
    MDOFftest0=trim(adjustl(MDOFftest0))

    if (BRoutput(1).eq.1) then
    write(MDOFfdisp1,*) trim(adjustl(monitor_file)),'/TDisX-',trim(adjustl(MDOFstr)),'.txt'
    MDOFfdisp1=trim(adjustl(MDOFfdisp1))
    write(MDOFfdisp2,*) trim(adjustl(monitor_file)),'/TDisY-',trim(adjustl(MDOFstr)),'.txt'
    MDOFfdisp2=trim(adjustl(MDOFfdisp2))
    open(MDOFmonitor,file=MDOFfdisp1,status='replace')
    close(MDOFmonitor)
    open(MDOFmonitor,file=MDOFfdisp2,status='replace')
    close(MDOFmonitor)
    end if
    if (BRoutput(1).eq.2) then
    write(MDOFfdisp1,*) trim(adjustl(monitor_file)),'/SDisX-',trim(adjustl(MDOFstr)),'.txt'
    MDOFfdisp1=trim(adjustl(MDOFfdisp1))
    write(MDOFfdisp2,*) trim(adjustl(monitor_file)),'/SDisY-',trim(adjustl(MDOFstr)),'.txt'
    MDOFfdisp2=trim(adjustl(MDOFfdisp2))
    open(MDOFmonitor,file=MDOFfdisp1,status='replace')
    close(MDOFmonitor)
    open(MDOFmonitor,file=MDOFfdisp2,status='replace')
    close(MDOFmonitor)
    end if

    if (BRoutput(2).eq.1) then
    write(MDOFfacc1,*) trim(adjustl(monitor_file)),'/TAccX-',trim(adjustl(MDOFstr)),'.txt'
    MDOFfacc1=trim(adjustl(MDOFfacc1))
    write(MDOFfacc2,*) trim(adjustl(monitor_file)),'/TAccY-',trim(adjustl(MDOFstr)),'.txt'
    MDOFfacc2=trim(adjustl(MDOFfacc2))
    open(MDOFmonitor,file=MDOFfacc1,status='replace')
    close(MDOFmonitor)
    open(MDOFmonitor,file=MDOFfacc2,status='replace')
    close(MDOFmonitor)
    end if
    if (BRoutput(2).eq.2) then
    write(MDOFfacc1,*) trim(adjustl(monitor_file)),'/SAccX-',trim(adjustl(MDOFstr)),'.txt'
    MDOFfacc1=trim(adjustl(MDOFfacc1))
    write(MDOFfacc2,*) trim(adjustl(monitor_file)),'/SAccY-',trim(adjustl(MDOFstr)),'.txt'
    MDOFfacc2=trim(adjustl(MDOFfacc2))
    open(MDOFmonitor,file=MDOFfacc1,status='replace')
    close(MDOFmonitor)
    open(MDOFmonitor,file=MDOFfacc2,status='replace')
    close(MDOFmonitor)
    end if
    
    if (BRoutput(3).eq.1) then
    write(MDOFfdamag1,*) trim(adjustl(monitor_file)),'/TDmgX-',trim(adjustl(MDOFstr)),'.txt'
    MDOFfdamag1=trim(adjustl(MDOFfdamag1))
    write(MDOFfdamag2,*) trim(adjustl(monitor_file)),'/TDmgY-',trim(adjustl(MDOFstr)),'.txt'
    MDOFfdamag2=trim(adjustl(MDOFfdamag2))
    open(MDOFmonitor,file=MDOFfdamag1,status='replace')
    close(MDOFmonitor)
    open(MDOFmonitor,file=MDOFfdamag2,status='replace')
    close(MDOFmonitor)
    end if
    if (BRoutput(3).eq.2) then
    write(MDOFfdamag1,*) trim(adjustl(monitor_file)),'/SDmgX-',trim(adjustl(MDOFstr)),'.txt'
    MDOFfdamag1=trim(adjustl(MDOFfdamag1))
    write(MDOFfdamag2,*) trim(adjustl(monitor_file)),'/SDmgY-',trim(adjustl(MDOFstr)),'.txt'
    MDOFfdamag2=trim(adjustl(MDOFfdamag2))
    open(MDOFmonitor,file=MDOFfdamag1,status='replace')
    close(MDOFmonitor)
    open(MDOFmonitor,file=MDOFfdamag2,status='replace')
    close(MDOFmonitor)
    end if

    if (BRoutput(4).eq.1) then
    write(MDOFfaccip1,*) trim(adjustl(monitor_file)),'/AcIpX-',trim(adjustl(MDOFstr)),'.txt'
    MDOFfaccip1=trim(adjustl(MDOFfaccip1))
    write(MDOFfaccip2,*) trim(adjustl(monitor_file)),'/AcIpY-',trim(adjustl(MDOFstr)),'.txt'
    MDOFfaccip2=trim(adjustl(MDOFfaccip2))
    open(MDOFmonitor,file=MDOFfaccip1,status='replace')
    close(MDOFmonitor)
    open(MDOFmonitor,file=MDOFfaccip2,status='replace')
    close(MDOFmonitor)
    end if

    if (BRoutput(5).eq.1) then
    write(MDOFfreat1,*) trim(adjustl(monitor_file)),'/ReactionX-',trim(adjustl(MDOFstr)),'.txt'
    MDOFfreat1=trim(adjustl(MDOFfreat1))
    write(MDOFfreat2,*) trim(adjustl(monitor_file)),'/ReactionY-',trim(adjustl(MDOFstr)),'.txt'
    MDOFfreat2=trim(adjustl(MDOFfreat2))
    open(MDOFmonitor,file=MDOFfreat1,status='replace')
    close(MDOFmonitor)
    open(MDOFmonitor,file=MDOFfreat2,status='replace')
    close(MDOFmonitor)
    end if
    
    if (BRoutput(6).eq.1) then
    write(MDOFfreat4,*) trim(adjustl(monitor_file)),'/MomentX-',trim(adjustl(MDOFstr)),'.txt'
    MDOFfreat4=trim(adjustl(MDOFfreat4))
    write(MDOFfreat5,*) trim(adjustl(monitor_file)),'/MomentY-',trim(adjustl(MDOFstr)),'.txt'
    MDOFfreat5=trim(adjustl(MDOFfreat5))
    open(MDOFmonitor,file=MDOFfreat4,status='replace')
    close(MDOFmonitor)
    open(MDOFmonitor,file=MDOFfreat5,status='replace')
    close(MDOFmonitor)
    end if
    
    
    open(MDOFmonitor,file='Config.txt')
        read(MDOFmonitor,*)
        read(MDOFmonitor,*) configtmp,MasspArea
	    read(MDOFmonitor,*) kclose
    close(MDOFmonitor)

    open(MDOFmonitor,file=MDOFftest)
        read(MDOFmonitor,*) MDOFbf,MDOFef
    close(MDOFmonitor)
        
    if(kclose.ne.0) kclose=1.d0
    open(MDOFmonitor,file=MDOFfinfo)
        read(MDOFmonitor,*) NumofBlgs,Maxdof
    close(MDOFmonitor)
    if(NumofBlgs.gt.0) then
        allocate(Blg(NumofBlgs))
        allocate(Blgdof(NumofBlgs))
        allocate(MDOFacl(NumofBlgs,3))
        allocate(NewOrNot(NumofBlgs))
        Blgdof=0
        NewOrNot=0
        call MDOFInputd(mpi_id)
        call MDOFMatrx(deltat)
    end if
    open(666,file='BuildingInfo.txt')
        read(666,*) MDOFMnum
    close(666)
    allocate(MDOFinput(3*MDOFMnum))
    allocate(MDOFinputfk(3*MDOFMnum))
    allocate(MDOFinputbuffer(mpi_np*MDOFMnum*3))
    allocate(MDOFinputfkbuffer(mpi_np*MDOFMnum*3))
    allocate(MDOFforceinput(6*MDOFMnum)) !(Fx,Fy,Fz,Mx,My,Mz)
    MDOFforceinput=0
    allocate(MDOFforceinputbuffer(6*MDOFMnum*mpi_np))   
    allocate(MDOFftranfg(mpi_np*2))
    open(666,file=MDOFftest0)
    do countyi=1,mpi_np
        read(666,*) MDOFftranfg(countyi),MDOFftranfg(countyi+mpi_np)
    end do
    close(666)

    return
end subroutine GenBFile
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GenBResultFile(itstmp,tt1tmp)
    use MDOF_MyData
    use speed_timeloop
	implicit none
    integer*4 :: tempj
    integer*4 :: itstmp
    real*8 :: tt1tmp

    MDOFmonitor1=10*(mpi_id+1)+7

    if (BRoutput(1).ge.1) then
        open(MDOFmonitor1,file=MDOFfdisp1,position='append')
        write(MDOFmonitor1,"(E16.7)",advance='NO') tt1tmp
        do MDOFBlgi=1,NumofBlgs 
            if (BRoutput(1).eq.1) then
                write(MDOFmonitor1,"(E16.7)",advance='NO') Blg(MDOFBlgi)%tempMDOFU1(Blg(MDOFBlgi)%NDOF,1)
            elseif (BRoutput(1).eq.2) then
                do tempj=1,Blg(MDOFBlgi)%NDOF
                    write(MDOFmonitor1,"(E16.7)",advance='NO') Blg(MDOFBlgi)%tempMDOFU1(tempj,1)
                end do
            end if
        end do
        write(MDOFmonitor1,"(A1)") " "
        close(MDOFmonitor1)
        !!!!!
        open(MDOFmonitor1,file=MDOFfdisp2,position='append')
        write(MDOFmonitor1,"(E16.7)",advance='NO') tt1tmp
        do MDOFBlgi=1,NumofBlgs 
            if (BRoutput(1).eq.1) then
                write(MDOFmonitor1,"(E16.7)",advance='NO') Blg(MDOFBlgi)%tempMDOFU1(Blg(MDOFBlgi)%NDOF,2)
            elseif (BRoutput(1).eq.2) then
                do tempj=1,Blg(MDOFBlgi)%NDOF
                    write(MDOFmonitor1,"(E16.7)",advance='NO') Blg(MDOFBlgi)%tempMDOFU1(tempj,2)
                end do
            end if
        end do
        write(MDOFmonitor1,"(A1)") " "
        close(MDOFmonitor1)
    end if

    if (BRoutput(2).ge.1) then
        open(MDOFmonitor1,file=MDOFfacc1,position='append')
        write(MDOFmonitor1,"(E16.7)",advance='NO') tt1tmp
        do MDOFBlgi=1,NumofBlgs 
            if (BRoutput(2).eq.1) then
                write(MDOFmonitor1,"(E16.7)",advance='NO') Blg(MDOFBlgi)%tempMDOFA1(Blg(MDOFBlgi)%NDOF,1)
            elseif (BRoutput(2).eq.2) then
                do tempj=1,Blg(MDOFBlgi)%NDOF
                    write(MDOFmonitor1,"(E16.7)",advance='NO') Blg(MDOFBlgi)%tempMDOFA1(tempj,1)
                end do
            end if
        end do
        write(MDOFmonitor1,"(A1)") " "
        close(MDOFmonitor1)
        !!!!!
        open(MDOFmonitor1,file=MDOFfacc2,position='append')
        write(MDOFmonitor1,"(E16.7)",advance='NO') tt1tmp
        do MDOFBlgi=1,NumofBlgs 
            if (BRoutput(2).eq.1) then
                write(MDOFmonitor1,"(E16.7)",advance='NO') Blg(MDOFBlgi)%tempMDOFA1(Blg(MDOFBlgi)%NDOF,2)
            elseif (BRoutput(2).eq.2) then
                do tempj=1,Blg(MDOFBlgi)%NDOF
                    write(MDOFmonitor1,"(E16.7)",advance='NO') Blg(MDOFBlgi)%tempMDOFA1(tempj,2)
                end do
            end if
        end do
        write(MDOFmonitor1,"(A1)") " "
        close(MDOFmonitor1)
    end if

    if ((BRoutput(3).eq.1).and.(itstmp.eq.nts)) then
        open(MDOFmonitor1,file=MDOFfdamag1,position='append')
        do MDOFBlgi=1,NumofBlgs 
            call MDOFDamageState(Blg(MDOFBlgi)%BlgMaxIDR(1:Blgdof(MDOFBlgi),1), &
                Blg(MDOFBlgi)%NDOF,Blgdof(MDOFBlgi),Blg(MDOFBlgi)%StrucType, &
                Blg(MDOFBlgi)%MDOFstate(1:Blg(MDOFBlgi)%NDOF,1),Blg(MDOFBlgi)%dval(1:Blgdof(MDOFBlgi),1:4))
            write(MDOFmonitor1,*) Blg(MDOFBlgi)%ID,Blg(MDOFBlgi)%MDOFstate(1:Blg(MDOFBlgi)%NDOF,1),Blg(MDOFBlgi)%BlgMaxIDR(1:Blgdof(MDOFBlgi),1)
        end do
        close(MDOFmonitor1)
        open(MDOFmonitor1,file=MDOFfdamag2,position='append')
        do MDOFBlgi=1,NumofBlgs 
            call MDOFDamageState(Blg(MDOFBlgi)%BlgMaxIDR(1:Blgdof(MDOFBlgi),2), &
                Blg(MDOFBlgi)%NDOF,Blgdof(MDOFBlgi),Blg(MDOFBlgi)%StrucType, &
                Blg(MDOFBlgi)%MDOFstate(1:Blg(MDOFBlgi)%NDOF,2),Blg(MDOFBlgi)%dval(1:Blgdof(MDOFBlgi),1:4))
            write(MDOFmonitor1,*) Blg(MDOFBlgi)%ID,Blg(MDOFBlgi)%MDOFstate(1:Blg(MDOFBlgi)%NDOF,2),Blg(MDOFBlgi)%BlgMaxIDR(1:Blgdof(MDOFBlgi),2)
        end do
        close(MDOFmonitor1)
    else if (BRoutput(3).eq.2)  then
        open(MDOFmonitor1,file=MDOFfdamag1,position='append')
        write(MDOFmonitor1,"(E16.7)") tt1tmp
        do MDOFBlgi=1,NumofBlgs 
            call MDOFDamageState(Blg(MDOFBlgi)%BlgMaxIDR(1:Blgdof(MDOFBlgi),1), &
                Blg(MDOFBlgi)%NDOF,Blgdof(MDOFBlgi),Blg(MDOFBlgi)%StrucType, &
                Blg(MDOFBlgi)%MDOFstate(1:Blg(MDOFBlgi)%NDOF,1),Blg(MDOFBlgi)%dval(1:Blgdof(MDOFBlgi),1:4))
            write(MDOFmonitor1,*) Blg(MDOFBlgi)%ID,Blg(MDOFBlgi)%MDOFstate(1:Blg(MDOFBlgi)%NDOF,1),Blg(MDOFBlgi)%BlgMaxIDR(1:Blgdof(MDOFBlgi),1)
        end do
        close(MDOFmonitor1)
        open(MDOFmonitor1,file=MDOFfdamag2,position='append')
        write(MDOFmonitor1,"(E16.7)") tt1tmp
        do MDOFBlgi=1,NumofBlgs 
            call MDOFDamageState(Blg(MDOFBlgi)%BlgMaxIDR(1:Blgdof(MDOFBlgi),2), &
                Blg(MDOFBlgi)%NDOF,Blgdof(MDOFBlgi),Blg(MDOFBlgi)%StrucType, &
                Blg(MDOFBlgi)%MDOFstate(1:Blg(MDOFBlgi)%NDOF,2),Blg(MDOFBlgi)%dval(1:Blgdof(MDOFBlgi),1:4))
            write(MDOFmonitor1,*) Blg(MDOFBlgi)%ID,Blg(MDOFBlgi)%MDOFstate(1:Blg(MDOFBlgi)%NDOF,2),Blg(MDOFBlgi)%BlgMaxIDR(1:Blgdof(MDOFBlgi),2)
        end do
        close(MDOFmonitor1)
    end if

    if (BRoutput(4).eq.1) then
        open(MDOFmonitor1,file=MDOFfaccip1,position='append')
        write(MDOFmonitor1,"(E16.7)",advance='NO') tt1tmp
        do MDOFBlgi=1,NumofBlgs 
            write(MDOFmonitor1,"(E16.7)",advance='NO') -MDOFacl(MDOFBlgi,1)
        end do
        write(MDOFmonitor1,"(A1)") " "
        close(MDOFmonitor1)
        open(MDOFmonitor1,file=MDOFfaccip2,position='append')
        write(MDOFmonitor1,"(E16.7)",advance='NO') tt1tmp
        do MDOFBlgi=1,NumofBlgs 
            write(MDOFmonitor1,"(E16.7)",advance='NO') -MDOFacl(MDOFBlgi,2)
        end do
        write(MDOFmonitor1,"(A1)") " "
        close(MDOFmonitor1)
    end if

    if (BRoutput(5).eq.1) then
        open(MDOFmonitor1,file=MDOFfreat1,position='append')
        write(MDOFmonitor1,"(E16.7)",advance='NO') tt1tmp
        do MDOFBlgi=1,NumofBlgs 
            write(MDOFmonitor1,"(E16.7)",advance='NO') MDOFforceinput(1+6*(MDOFftranfg(mpi_id+1)+MDOFBlgi-2))
        end do
        write(MDOFmonitor1,"(A1)") " "
        close(MDOFmonitor1)
        open(MDOFmonitor1,file=MDOFfreat2,position='append')
        write(MDOFmonitor1,"(E16.7)",advance='NO') tt1tmp
        do MDOFBlgi=1,NumofBlgs 
            write(MDOFmonitor1,"(E16.7)",advance='NO') MDOFforceinput(2+6*(MDOFftranfg(mpi_id+1)+MDOFBlgi-2))
        end do
        write(MDOFmonitor1,"(A1)") " "
        close(MDOFmonitor1)
    end if

    if (BRoutput(6).eq.1) then
        open(MDOFmonitor1,file=MDOFfreat4,position='append')
        write(MDOFmonitor1,"(E16.7)",advance='NO') tt1tmp
        do MDOFBlgi=1,NumofBlgs 
            write(MDOFmonitor1,"(E16.7)",advance='NO') MDOFforceinput(4+6*(MDOFftranfg(mpi_id+1)+MDOFBlgi-2))
        end do
        write(MDOFmonitor1,"(A1)") " "
        close(MDOFmonitor1)
        open(MDOFmonitor1,file=MDOFfreat5,position='append')
        write(MDOFmonitor1,"(E16.7)",advance='NO') tt1tmp
        do MDOFBlgi=1,NumofBlgs 
            write(MDOFmonitor1,"(E16.7)",advance='NO') MDOFforceinput(5+6*(MDOFftranfg(mpi_id+1)+MDOFBlgi-2))
        end do
        write(MDOFmonitor1,"(A1)") " "
        close(MDOFmonitor1)
    end if

    return
end subroutine GenBResultFile
