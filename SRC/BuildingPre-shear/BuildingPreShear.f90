!  BuildingPre.f90 
!
!  FUNCTIONS:
!  BuildingPre - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: BuildingPre
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program BuildingPreShear

    use GenModelPara
    
    implicit none
    call LimitModule   
    call checkfile    
    call GenFile
    call GenBldPeriod

    pause
    stop

    end program BuildingPreShear

