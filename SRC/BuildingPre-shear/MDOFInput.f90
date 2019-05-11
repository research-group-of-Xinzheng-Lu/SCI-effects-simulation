subroutine MDOFInput(NumofBlgs,Blg,Maxdof)
	use GenModelPara
	implicit none
    Character(len=200) tempchar, StrucTypeName
    integer Maxdof,NumofBlgs
    integer I,J,tempint
    real*8 tempreal
	type(BlgInfo):: Blg(NumofBlgs)

    Maxdof=0

    open(11,file='BuildingInfo.txt',position="rewind") 
	    read(11,*)  
        do I=1, NumofBlgs
 		    read(11, *) Blg(I)%ID, Blg(I)%StrucType, Blg(I)%NDOF,Blg(I)%Floor_h, Blg(I)%Area !, Blg(I)%position(1:2)
            Blg(I)%position(3)=0.d0
		    Blg(I)%Height=Blg(I)%NDOF*Blg(I)%Floor_h		
            do J=1, Blg(I)%NDOF
		        read(11, *)Blg(I)%props(J,1:10)
            enddo
            Maxdof=max(Maxdof,Blg(I)%NDOF)
        end do
        do I=1, NumofBlgs
 		    read(11, *) 		
            do J=1, Blg(I)%NDOF
		        read(11,*)Blg(I)%DamageCriteria(J,1:4)
            enddo
        end do
	close(11)
    write(*,"(A25)") "[Info] BuildingInfo read."
    return
end subroutine MDOFInput