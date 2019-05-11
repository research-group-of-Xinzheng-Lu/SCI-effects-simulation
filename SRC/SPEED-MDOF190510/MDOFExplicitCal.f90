!Computation for MDOF shear model. By Group of XZ LU
!NumofBlgs: Total number of buildings in the ith core
!BlgID: ID
!MDOFacltmp: Acceleration input
!Blg: All building informations are contained here.
!tempMDOFU0,tempMDOFU1,tempMDOFU: Displacement at n-1,n,n+1

subroutine ExplicitCalShear (BlgID,dof,MDOFacltmp,direction, &
                        NewOrNottmp)
use MDOF_MyData
implicit none
integer I,J,BlgID,dof
real*8 E(dof),dE(dof),S(dof)
real*8 MDOFPeff0(dof)
integer NewOrNottmp,direction
real*8 MDOFacltmp,MDOFPeff(dof)

do I=1,dof
    MDOFPeff(I)=MDOFacltmp*Blg(BlgID)%BlgM(I,I)
end do

    do I=1, dof-1; 
	    MDOFPeff0(I)=MDOFPeff(I)-Blg(BlgID)%MDOFItF(I,direction)+Blg(BlgID)%MDOFItF(I+1,direction);
	end do
    MDOFPeff0(dof)=MDOFPeff(dof)-Blg(BlgID)%MDOFItF(dof,direction);
    call Central_Difference(Blg(BlgID)%BlgM(1:dof,1:dof),Blg(BlgID)%BlgC(1:dof,1:dof), &
        Blg(BlgID)%BlgM_Inv(1:dof,1:dof), Blg(BlgID)%dt, &
        Blg(BlgID)%tempMDOFU(1:dof,direction), Blg(BlgID)%tempMDOFU1(1:dof,direction), &
        Blg(BlgID)%tempMDOFU0(1:dof,direction), MDOFPeff0, dof, NewOrNottmp); 
    Blg(BlgID)%dMDOFIDR(1,direction)=Blg(BlgID)%tempMDOFU(1,direction)-Blg(BlgID)%MDOFIDR(1,direction)
    do J=2, dof
        Blg(BlgID)%dMDOFIDR(J,direction)=(Blg(BlgID)%tempMDOFU(J,direction)-Blg(BlgID)%tempMDOFU(J-1,direction))-Blg(BlgID)%MDOFIDR(J,direction)
    end do
	Blg(BlgID)%MDOFIDR(1,direction)=Blg(BlgID)%tempMDOFU(1,direction); 
	do J=2, dof
        Blg(BlgID)%MDOFIDR(J,direction)=Blg(BlgID)%tempMDOFU(J,direction)-Blg(BlgID)%tempMDOFU(J-1,direction)
       ! write(*,*) MDOFIDR(J)
    end do
    do J=1,	dof
        Blg(BlgID)%BlgMaxIDR(J,direction)=max(Blg(BlgID)%BlgMaxIDR(J,direction),abs(Blg(BlgID)%MDOFIDR(J,direction)/Blg(BlgID)%Floor_h))
    end do
    Blg(BlgID)%tempMDOFA1(1:dof,direction)=Blg(BlgID)%tempMDOFU(1:dof,direction)+Blg(BlgID)%tempMDOFU0(1:dof,direction)-2.d0*Blg(BlgID)%tempMDOFU1(1:dof,direction)
    Blg(BlgID)%tempMDOFA1(1:dof,direction)=Blg(BlgID)%tempMDOFA1(1:dof,direction)/Blg(BlgID)%dt2-MDOFacltmp
    Blg(BlgID)%tempMDOFU0(1:dof,direction)= Blg(BlgID)%tempMDOFU1(1:dof,direction)
    Blg(BlgID)%tempMDOFU1(1:dof,direction)= Blg(BlgID)%tempMDOFU(1:dof,direction)
    Blg(BlgID)%tempMDOFU(1:dof,direction)=0 

    do I=1,dof
        E(I)=Blg(BlgID)%MDOFIDR(I,direction)
        dE(I)=Blg(BlgID)%dMDOFIDR(I,direction)
        S(I)=Blg(BlgID)%MDOFItF(I,direction)
        if((abs(Blg(BlgID)%MDOFIDR(I,direction))/Blg(BlgID)%Floor_h).GT.Blg(BlgID)%dval(I,4))then
	        Blg(BlgID)%MDOFIDeath(I,direction)=1
	    endif
        call ksteel02(Blg(BlgID)%props(I,1:10,1),s(I),e(I),de(I),Blg(BlgID)%MDOFEt(I,direction), &
             Blg(BlgID)%MDOFstatev(I,1:11,direction),Blg(BlgID)%MDOFspd(I,direction), &
             Blg(BlgID)%MDOFyield(I,direction), &
             Blg(BlgID)%MDOFIDeath(I,direction),Blg(BlgID)%BlgM(1:dof,1:dof),dof) 
		Blg(BlgID)%MDOFItF(I,direction)=S(I) 
    end do




end subroutine ExplicitCalShear