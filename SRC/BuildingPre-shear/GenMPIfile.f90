subroutine GenMPIfile

use GenModelPara

implicit none
integer*4,allocatable :: flg(:),maxdof(:),flgdof(:),readline(:,:)
integer*4 :: noteflag,i,j,k,begin
real*8 :: partflag,partflag2
character(len=100) :: filestr,filestr2,filestr3,filestr4
allocate(maxdof(nbld))
maxdof=0
do i=1,nbld
    if(i.gt.1) then
        Building(i)%stpflag=Building(i-1)%stpflag+Building(i)%NDOF*Building(i)%NDOF
        maxdof(i)=max(maxdof(i-1),Building(i)%NDOF)
    elseif(i.eq.1) then
        Building(i)%stpflag=Building(i)%NDOF*Building(i)%NDOF
        maxdof(i)=Building(i)%NDOF
    end if
end do

allocate(flg(np))
allocate(flgdof(np))
allocate(readline(np,3))
partflag=real(Building(nbld)%stpflag)/np  !average computation load
noteflag=1
    flg=0
    flgdof=0
    readline=0
    do j=1,np
        partflag2=partflag*j
        if(noteflag.le.nbld) then
        do i=noteflag,nbld
            if(partflag2.lt.(Building(i)%stpflag)) then
                if(i.gt.1) then
                flg(j)=i-1
               ! flgdof(j)=maxdof(i-1)
                noteflag=i
                elseif(i.eq.1) then
                flg(j)=0
                !flgdof(j)=0
                noteflag=i
                end if
                exit
            elseif(partflag2.eq.(Building(i)%stpflag)) then
                flg(j)=i
                !flgdof(j)=maxdof(i)
                noteflag=i+1
                exit
            end if
        end do
        end if
    end do
    !!prevent calculating error
    flg(np)=nbld
    !flgdof(np)=maxdof(nbld)
    open(11,file="./FILES_MPI/mpitest0.txt",status="replace")
        do i=1,np
        write(filestr,*) i
        write(filestr4,*) './FILES_MPI/BuildingInfo',trim(adjustl(filestr)),'.txt'
        write(filestr,*) './FILES_MPI/mpitest',trim(adjustl(filestr)),'.txt'

        open(22,file=filestr,status='replace')
        open(55,file=filestr4,status='replace')
        if(i.eq.1) then
            write(11,"(2I<digitout>)")1,flg(i)  !Beginning Id, number of buildings in i th core
            write(22,"(2I<digitout>)")1,flg(i)  !Beginning Id, end Id
            if(flg(i).ge.1) then
            flgdof(i)=maxval(Building(1:flg(i))%NDOF)
            else
            flgdof(i)=0
            end if
            write(55,"(2I<digitout>)")flg(i),flgdof(i)
        elseif(i.gt.1) then
            write(11,"(2I<digitout>)")flg(i-1)+1,flg(i)-flg(i-1)  !The beginning Building Id for i th core is flg(i-1)+1; there are flg(i)-flg(i-1) buildings
            write(22,"(2I<digitout>)")flg(i-1)+1,flg(i)  !Beginning Id, end Id
            if(flg(i).ge.(flg(i-1)+1))then
            flgdof(i)=maxval(Building(flg(i-1)+1:flg(i))%NDOF)
            else
            flgdof(i)=0
            end if
            write(55,"(2I<digitout>)")flg(i)-flg(i-1),flgdof(i)
        end if
        close(22)
        close(33)
        close(44)
        close(55)
        end do
    close(11)
    if(np.gt.1)then
    k=1
        if(flg(k).gt.0) then
        write(filestr,*) k
        write(filestr2,*) './FILES_MPI/BuildingInfo',trim(adjustl(filestr)),'.txt'
        open(11,file=filestr2,position='append')
        do i=1,flg(k)
            write(11, "(3I<digitout>,2E<digitoute1>.<digitoute2>)") Building(I)%ID, Building(I)%StrucType, &
            & Building(I)%NDOF,Building(I)%Floor_h, Building(I)%Area
            do J=1, Building(I)%NDOF
		        write(11, "(10E<digitoute1>.<digitoute2>)")Building(I)%props(J,1:10)
            enddo
        end do
        do i=1,flg(k)
            write(11, "(3I<digitout>,2E<digitoute1>.<digitoute2>)") Building(I)%ID, Building(I)%StrucType, &
            & Building(I)%NDOF,Building(I)%Floor_h, Building(I)%Area
            do J=1, Building(I)%NDOF
		        write(11, "(4E<digitoute1>.<digitoute2>)")Building(I)%DamageCriteria(J,1:4)
            enddo
        end do
        write(11, "(A6)") "Period"
        do i=1,flg(k)
            write(11,"(I<digitout>,5E<digitoute1>.<digitoute2>)") Building(I)%ID, Building(I)%T1, Building(I)%T2,&
            & Building(I)%Calpha,Building(I)%Cbeta,Building(I)%TN
        end do
        close(11)
        end if
    do k=2,np
        if((flg(k)-flg(k-1)).gt.0) then
        write(filestr,*) k
        write(filestr2,*) './FILES_MPI/BuildingInfo',trim(adjustl(filestr)),'.txt'
        open(11,file=filestr2,position='append')
        do i=flg(k-1)+1,flg(k)
            write(11, "(3I<digitout>,2E<digitoute1>.<digitoute2>)") Building(I)%ID, Building(I)%StrucType, &
            & Building(I)%NDOF,Building(I)%Floor_h, Building(I)%Area
            do J=1, Building(I)%NDOF
		        write(11, "(10E<digitoute1>.<digitoute2>)")Building(I)%props(J,1:10)
            enddo
        end do
        do i=flg(k-1)+1,flg(k)
            write(11, "(3I<digitout>,2E<digitoute1>.<digitoute2>)") Building(I)%ID, Building(I)%StrucType, &
            & Building(I)%NDOF,Building(I)%Floor_h, Building(I)%Area
            do J=1, Building(I)%NDOF
		        write(11, "(4E<digitoute1>.<digitoute2>)")Building(I)%DamageCriteria(J,1:4)
            enddo
        end do
        write(11, "(A6)") "Period"
        do i=flg(k-1)+1,flg(k)
            write(11,"(I<digitout>,5E<digitoute1>.<digitoute2>)") Building(I)%ID, Building(I)%T1, Building(I)%T2,&
            & Building(I)%Calpha,Building(I)%Cbeta,Building(I)%TN
        end do
        close(11)
        end if
    end do 
    else
        open(11,file='./FILES_MPI/BuildingInfo1.txt',position='append')
        do i=1,nbld
            write(11, "(3I<digitout>,2E<digitoute1>.<digitoute2>)") Building(I)%ID, Building(I)%StrucType, &
            & Building(I)%NDOF,Building(I)%Floor_h, Building(I)%Area
            do J=1, Building(I)%NDOF
		        write(11, "(10E<digitoute1>.<digitoute2>)")Building(I)%props(J,1:10)
            enddo
        end do
        do i=1,nbld
            write(11, "(3I<digitout>,2E<digitoute1>.<digitoute2>)") Building(I)%ID, Building(I)%StrucType, &
            & Building(I)%NDOF,Building(I)%Floor_h, Building(I)%Area
            do J=1, Building(I)%NDOF
		        write(11, "(4E<digitoute1>.<digitoute2>)")Building(I)%DamageCriteria(J,1:4)
            enddo
        end do
        write(11, "(A6)") "Period"
        do i=1,nbld
            write(11,"(I<digitout>,5E<digitoute1>.<digitoute2>)") Building(I)%ID, Building(I)%T1, &
            & Building(I)%T2,Building(I)%Calpha,Building(I)%Cbeta,Building(I)%TN
        end do
        close(11)
    end if


return
end subroutine GenMPIfile