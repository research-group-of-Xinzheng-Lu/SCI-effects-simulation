!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module GenModelPara
    implicit none
    !Codelimit
    integer*4, parameter :: codelimit = 0      !0: No limit; non0: Building number and Site size Limit
    integer*4, parameter :: digitout = 10
    integer*4, parameter :: digitoute1 = 15
    integer*4, parameter :: digitoute2 = 7
    integer*4 :: limittag 
    integer*4 :: nbldlimit
    !Building
    integer*4 :: nbld                           !total number of buildings
    integer*4 :: shapeornot,nposfunc
    real*8 :: MasspArea,damping_bld
    type BlgInfo
        integer*4 :: ID	
        integer*4 :: NDOF	!num of story
        integer*4 :: stpflag
        integer*4 :: StrucType
        integer*4 :: posfunc
        real*8 :: position(3)
        real*8 :: height	!total height
        real*8 :: Floor_h	!story height
        real*8 :: T1, T2, TN
        real*8 :: Area	
        real*8 :: props(200,10)			
        real*8 :: DamageCriteria(200,4)
        real*8 :: damping
        real*8 :: Calpha,Cbeta
    end type
    type(BlgInfo),allocatable:: Building(:)
    !Analysis
    integer*4 :: np
    logical :: exist_f,exist_c
end module GenModelPara

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine LimitModule
    use GenModelPara
    implicit none
    if (codelimit.eq.0) then
        nbldlimit=0
    else
        nbldlimit=50
    end if
end subroutine LimitModule
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ModalAnalysis(K,M,NDOF,T,vec)
	use mkl95_lapack
	implicit none 
	integer NDOF, I,tempi,tempj
	real*8 K(1,NDOF,NDOF), M(1,NDOF,NDOF), M_temp(NDOF,NDOF)
	real*8 T(NDOF), vec(NDOF, NDOF)
	intent(in) K, M, NDOF
	intent(out) T, vec
	do tempi=1,NDOF
        do tempj=1,NDOF
            Vec(tempi,tempj)=K(1,tempi,tempj)
	        M_Temp(tempi,tempj)=M(1,tempi,tempj)
        end do
    end do
 	call sygv(vec, M_temp, T, 1, "V", "U")
    do I=1, NDOF
		T(I)=2.0*3.1415926/sqrt(T(I))
		vec(1:NDOF,I)=vec(1:NDOF,I)/Vec(NDOF,I);
	end do
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine outputbug(bugid)
    use GenModelPara
    implicit none
    integer*4 :: bugid
    selectcase(bugid)
    case(0)
        write(*,"(A12)") "[Error] Bye."
    case(1)
        write(*,"(A36)") "[Error] BuildingInfo.txt is missing."
    case(2)
        write(*,"(A64)") "[Warning] Config.txt is missing. Default values will be adopted."
    case(3)
        write(*,"(A102)") "[Error] We are sorry that you do not have the authority to calculate more than",nbldlimit," buildings."
    end select
    return
end subroutine outputbug
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine checkfile
    use GenModelPara
    implicit none
    inquire(file="BuildingInfo.txt",exist=exist_f)
    if (exist_f) then
        continue
    else
        call outputbug(1)
        pause 
        stop
    end if
    inquire(file="Config.txt",exist=exist_c)
    if (exist_c) then
        call readconfig
        continue
    else
        call outputbug(2)
        call defaultconfig
        continue
    end if
    return
end subroutine checkfile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine judgelimit(a,b,c,tag)
    use GenModelPara
    implicit none
    integer*4::a,tag
    real*8::b,c
    tag=0
    if(nbldlimit.gt.0) then
        if(a.gt.nbldlimit) then
            call outputbug(3)
            tag=1
        end if
    end if
end subroutine judgelimit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine defaultconfig
    use GenModelPara
    implicit none
    damping_bld=-1
    MasspArea=1000.d0
    np=4
end subroutine defaultconfig
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readconfig
    use GenModelPara
    implicit none
    integer*4 :: i
    open(11,file="Config.txt",position="rewind")
        read(11,*) np
        read(11,*) damping_bld,MasspArea
    close(11)
    return
end subroutine readconfig