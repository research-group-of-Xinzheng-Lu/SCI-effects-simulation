subroutine MDOFMatrx(NumofBlgs,Blg,Maxdof,BlgM,BlgK,Cdamp)
	use GenModelPara
	implicit none
    integer I,J, Maxdof,NumofBlgs
    real*8,allocatable:: tempprops(:,:,:)
    real*8 BlgM(NumofBlgs,Maxdof,Maxdof),BlgK(NumofBlgs,Maxdof,Maxdof)
	real*8,allocatable:: BlgT(:),BlgVec(:,:)
	real*8 Cdamp(NumofBlgs)
    type(BlgInfo):: Blg(NumofBlgs)
    
    do I=1,NumofBlgs
    if(damping_bld.gt.0) then
        Cdamp(I)=damping_bld
    else
       Cdamp(I)=0.05
    end if
    end do

    
    do I=1, NumofBlgs
            !!!!!!!!!!!MassMatrix!!!!!!!!!!!!
            do J=1,Blg(I)%NDOF
                BlgM(I,J,J)=MasspArea*Blg(I)%Area
            end do
            !!!!!!!!!!!MassMatrix!!!!!!!!!!!!
            !!!!!!!!!!!StiffnessMatrix!!!!!!!!!!!!
            do J=1,Blg(I)%NDOF
                BlgK(I,J,J)=BlgK(I,J,J)+Blg(I)%props(J,1) 
                if(J>1) then
                    BlgK(I,J-1,J-1)=BlgK(I,J-1,J-1)+Blg(I)%props(J,1)
                    BlgK(I,J,J-1)=BlgK(I,J,J-1)-Blg(I)%props(J,1)
                    BlgK(I,J-1,J)=BlgK(I,J-1,J)-Blg(I)%props(J,1)
                end if
            end do
            !!!!!!!!!!!StiffnessMatrix!!!!!!!!!!!!
            allocate(BlgT(Blg(I)%NDOF))
            allocate(BlgVec(Blg(I)%NDOF,Blg(I)%NDOF))
            BlgT=0
            BlgVec=0
            call ModalAnalysis (BlgK(I,1:Blg(I)%NDOF,1:Blg(I)%NDOF), BlgM(I,1:Blg(I)%NDOF,1:Blg(I)%NDOF), Blg(I)%NDOF, BlgT, BlgVec)	!Ä£Ì¬·ÖÎö
            Blg(I)%T1=BlgT(1);
            Blg(I)%TN=BlgT(Blg(I)%NDOF);
	            if(Blg(I)%NDOF.gt.1)then
		            Blg(I)%T2=BlgT(2)
	            else
		            Blg(I)%T2=Blg(I)%T1/3.0  
	            endif 
            deallocate(BlgT)
            deallocate(BlgVec)    
	enddo
	return
end subroutine MDOFMatrx