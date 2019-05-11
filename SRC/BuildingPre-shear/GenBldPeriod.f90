subroutine GenBldPeriod

use GenModelPara

implicit none

real*8,allocatable:: BlgM(:,:,:),BlgK(:,:,:),Cdamp(:)
integer Maxdof,NumofBlgs,I
real*8 Comiga1,Comiga2

Maxdof=0
Cdamp=0

open(11,file="BuildingInfo.txt",position="rewind")
    read(11,*) nbld
    call judgelimit(nbld,0.d0,0.d0,limittag)
    if(limittag.eq.1) then
        close(11) 
        pause 
        stop
    end if
close(11)

allocate(Building(nbld))
allocate(Cdamp(nbld))
call MDOFInput(nbld,Building,Maxdof)
allocate(BlgM(nbld,Maxdof,Maxdof))
allocate(BlgK(nbld,Maxdof,Maxdof))
BlgM=0
BlgK=0
write(*,"(A13)") "Preparing ..."
call MDOFMatrx(nbld,Building,Maxdof,BlgM,BlgK,Cdamp)
deallocate(BlgM)
deallocate(BlgK)
open(11,file='BuildingPeriod.txt',status='replace')
write(11,"(2I<digitout>)")nbld, Maxdof
do I=1,nbld
    Comiga1=2.*3.1415926535/Building(I)%T1
    Comiga2=2.*3.1415926535/Building(I)%T2; ! 计算圆频率
    Building(I)%Calpha=2.*Cdamp(I)*Comiga1*Comiga2/(Comiga1+Comiga2)
    Building(I)%Cbeta=2.*Cdamp(I)/(Comiga1+Comiga2); ! 计算瑞雷阻尼系数
    write(11,"(I<digitout>,5E<digitoute1>.<digitoute2>)") Building(I)%ID, Building(I)%T1,&
    & Building(I)%T2,Building(I)%Calpha,Building(I)%Cbeta,Building(I)%TN
end do
close(11)
call GenMPIfile


return


end subroutine GenBldPeriod