!    Copyright (C) 2012 The SPEED FOUNDATION
!    Author: Ilario Mazzieri
!
!    This file is part of SPEED.
!
!    SPEED is free software; you can redistribute it and/or modify it
!    under the terms of the GNU Affero General Public License as
!    published by the Free Software Foundation, either version 3 of the
!    License, or (at your option) any later version.
!
!    SPEED is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Affero General Public License for more details.
!
!    You should have received a copy of the GNU Affero General Public License
!    along with SPEED.  If not, see <http://www.gnu.org/licenses/>.


!> @brief Makes Butcher array for Runge-Kutta scheme 
!! @author Ilario Mazzieri
!> @note Compute low storage Butcher array for Runge-Kutta(p,s) time integration scheme 
!  where p is the order of the scheme and s is the number of stages
!  Schemes Implemented: RK(2,2), RK(3,3), RK(4,4)  
!> @date September, 2013 
!> @version 1.0
!> @param[in] order  p
!> @param[in] stages  s
!> @param[out] A_low  low storage Butcher array A  
!> @param[out] B_low  low storage Butcher array B
!> @param[out] c  Butcher array c  



     subroutine MAKE_BUTCHERARRAY(order, stages, A_low, B_low, c)
          
     implicit none
     
     integer*4:: order, stages
     real*8, dimension(stages,stages) :: A
     real*8, dimension(stages) :: b, c, A_low, B_low

     A = 0; b = 0; c = 0; A_low = 0; B_low = 0;
     
     select case(order)   

        case(2)               
           ! RK(2,2)
           select case(stages)
              case(2)
                 c(2) = 1
                 b(1) = 0.5; b(2)=0.5
                 A(2,1) = 1
                 
                 call MAKE_LOW_STORAGE_COEFFICIENTS(stages,A,b,c,A_low,B_low)
                 
              case default   
                  write(*,'(A)') 'RK scheme not implemented'
           end select               
        
        case(3)
           ! RK(3,3)
           select case(stages)
              case(3)
                 c(2) = 0.5; c(3) = 1
                 b(1) = 1.d0/6; b(2) = 2.d0/3; b(3) = 1.d0/6
                 A(2,1) = 0.5; 
                 A(3,1) = -1; A(3,2) = 2
                 
                 call MAKE_LOW_STORAGE_COEFFICIENTS(stages,A,b,c,A_low,B_low)
                                  

              case default   
                  write(*,'(A)') 'RK scheme not implemented'
           end select                       
        
        case(4)

           select case(stages)
              ! RK(4,4)
              case(4)
                 c(2) = 0.5; c(3) = 0.5; c(4) = 1
                 b(1) = 1.d0/6; b(2) = 1.d0/3; b(3) = 1.d0/3; b(4) = 1.d0/6
                 A(2,1) = 0.5; 
                 A(3,2) = 0.5;
                 A(4,3) = 1;
                 
                 call MAKE_LOW_STORAGE_COEFFICIENTS(stages,A,b,c,A_low,B_low)
               
                 
              case default   
                  write(*,'(A)') 'RK scheme not implemented'
                  
           end select                
       
       
      
       
       
        case default
           write(*,'(A)') 'RK scheme not implemented'
           
      end select    
           
     
     end subroutine MAKE_BUTCHERARRAY     
     
!**************************************************************************************************


!> @brief Computes Low-storage memory coefficients for Runge-Kutta scheme 
!! @author Ilario Mazzieri
!> @note Compute low storage Butcher array for Runge-Kutta(p,s) time integration scheme 
!!  where p is the order of the scheme and s is the number of stages  
!> @date September, 2013 
!> @version 1.0
!> @param[in] stages  s
!> @param[in] A Butcher array
!> @param[in] b Butcher array
!> @param[in] c Butcher array
!> @param[out] A_low  low storage Butcher array A  
!> @param[out] B_low  low storage Butcher array B


     
     subroutine MAKE_LOW_STORAGE_COEFFICIENTS(stages,A,b,c,A_low,B_low)
     
     implicit none
     
     integer*4:: stages
     real*8, dimension(stages,stages) :: A
     real*8, dimension(stages) :: b, c, A_low, B_low
     integer*4 :: i
     

     ! Compute Low storage coefficients according to the
     ! paper Toulorge and Desmet: "Optimal Runge-Kutta schemes for discontinuous
     ! Galerkin space discretizations applied to wave propagation problems"
     
     A_low(1) = 0;
     
     do i = 1, stages -1
        B_low(i) = a(i+1,i)
     enddo
     B_low(stages) = b(stages)
    
     do i = 2, stages 
        if (b(i) .ne. 0.d0) then
           A_low(i) = (b(i-1) - B_low(i-1))/b(i);
        else
           A_low(i) = (A(i+1,i-1) - c(i))/B_low(i);
        endif
     enddo   
        
     
     end subroutine MAKE_LOW_STORAGE_COEFFICIENTS
     
     
     
