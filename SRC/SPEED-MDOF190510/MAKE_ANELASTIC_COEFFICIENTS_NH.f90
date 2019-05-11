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

!> @brief Compute anelastic coefficients for Standard Linear Solid model, when the not-honoring strategy is applied
!! USE THE LIBRARY QR_SOLVE!
!! @author Ilario Mazzieri
!> @date November, 2014 
!> @version 1.0
!> @param[in] nmat number of materials 
!> @param[in] N_SLS number of subdivision in the standard linear solid (dafult N_SLS=3)
!!            possible values are N_SLS=3,4,5 (see MODULES.f90)
!> @param[in] lambda_el Lame 1st modulus given node by node
!> @param[in] mu_el Lame 2nd modulus given node by node
!> @param[in] QS_nh quality factors for s-waves
!> @param[in] QP_nh quality factors for p-waves
!> @param[in] f_val value for the reference frequency
!> @param[in] mpi_id id of the mpi-process 
!> @param[out] Y_lambda_nh anelastic coefficient (Lame 1st modulus)
!> @param[out] Y_mu_nh anelastic coefficient (Lame 2nd modulus)
!> @param[out] frequency_range frequency range where quality factors are assumed to be constant

      subroutine MAKE_ANELASTIC_COEFFICIENTS_NH(nn, N_SLS, lambda_el, mu_el,&
                                            rho_el, QS_nh, QP_nh, f_val, &
                                            Y_lambda_nh, Y_mu_nh, frequency_range, mpi_id)
      
      use speed_exit_codes                  

      implicit none
 
      integer*4 :: nn, mpi_id, N_SLS, i, j, k, &
                   PIVOT_QP(N_SLS),PIVOT_QS(N_SLS)
                   
      real*8 :: f_val, fmin, fmax, esp1, esp2, deltax, &
                vp2, vs2, pi
      
      real*8 :: QS_nh, QP_nh, &
                Y_lambda_nh(1,N_SLS),&
                Y_mu_nh(1,N_SLS), &
                frequency_range(N_SLS),&
                frequency_range_sampling(2*N_SLS-1),&
                YS(N_SLS),YP(N_SLS), &
                A_QP(2*N_SLS-1,N_SLS),RHS_QP(2*N_SLS-1),A_QS(2*N_SLS-1,N_SLS),RHS_QS(2*N_SLS-1),&
                rho, rho_el(nn,nn,nn), &
                lambda, lambda_el(nn,nn,nn), &
                mu, mu_el(nn,nn,nn), &
                rho_all, mu_all, lambda_all 
      
      rho_all = 0.d0; lambda_all = 0.d0; mu_all = 0.d0
      
      do i =1, nn
         do j =1, nn
            do k =1, nn
               mu_all = mu_all + mu_el(i,j,k)
               lambda_all = lambda_all + lambda_el(i,j,k)
               rho_all = rho_all + rho_el(i,j,k)
            enddo
         enddo
      enddo      
                         
      mu = mu_all/(nn**3)
      lambda = lambda_all/(nn**3)                
      rho = rho_all/(nn**3) 
 
      pi = 4.d0*datan(1.d0);
      
      if (f_val .le. 10) then
         fmin = 0.05d0; !fmax = 5.d0;
      elseif(f_val .le. 100) then
         fmin = 5.d0; !fmax = 100.d0;   
      elseif(f_val .le. 1000) then
         fmin = 50.d0; !fmax = 1000.d0;   
      elseif(f_val .le. 10000) then
         fmin = 500.d0; !fmax = 10000.d0;   
      elseif(f_val .le. 100000) then
         fmin = 5000.d0; !fmax = 100000.d0;
      endif
      
      fmax=100*fmin;
      
      !if(mpi_id .eq. 0) then
      !   write(*,*) 'Frequency range where Quality Factor is assumed to be constant is:'
      !   write(*,*) 'FMIN =' ,fmin, '  FMAX =', fmax;
      !endif
            
     esp1 = LOG10(fmin); esp2 = LOG10(fmax);
     deltax = (esp2-esp1)/(N_SLS-1);
     

     if (N_SLS .eq. 3) then

        frequency_range(1) = 2.*pi*fmin;
        frequency_range(2) = 2.*pi*10*fmin;
        frequency_range(3) = 2.*pi*fmax;
        frequency_range_sampling(1) = frequency_range(1);
        frequency_range_sampling(2) = 0.5*(frequency_range(1)+frequency_range(2))
        frequency_range_sampling(3) = frequency_range(2);
        frequency_range_sampling(4) = 0.5*(frequency_range(2)+frequency_range(3))
        frequency_range_sampling(5) = frequency_range(3);
                

     elseif (N_SLS .eq. 4) then

        frequency_range(1) = 2.*pi*fmin;
        frequency_range(2) = 2.*pi*10**(esp1+deltax)
        frequency_range(3) = 2.*pi*10**(esp1+2*deltax)
        frequency_range(4) = 2.*pi*fmax;
        frequency_range_sampling(1) = frequency_range(1);
        frequency_range_sampling(2) = 0.5*(frequency_range(1)+frequency_range(2))
        frequency_range_sampling(3) = frequency_range(2);
        frequency_range_sampling(4) = 0.5*(frequency_range(2)+frequency_range(3))
        frequency_range_sampling(5) = frequency_range(3);
        frequency_range_sampling(6) = 0.5*(frequency_range(3)+frequency_range(4))
        frequency_range_sampling(7) = frequency_range(4);


      elseif (N_SLS .eq. 5) then

        frequency_range(1) = 2.*pi*fmin;
        frequency_range(2) = 2.*pi*10**(esp1+deltax)
        frequency_range(3) = 2.*pi*10**(esp1+2*deltax)
        frequency_range(4) = 2.*pi*10**(esp1+3*deltax)
        frequency_range(5) = 2.*pi*fmax;
        frequency_range_sampling(1) = frequency_range(1);
        frequency_range_sampling(2) = 0.5*(frequency_range(1)+frequency_range(2))
        frequency_range_sampling(3) = frequency_range(2);
        frequency_range_sampling(4) = 0.5*(frequency_range(2)+frequency_range(3))
        frequency_range_sampling(5) = frequency_range(3);
        frequency_range_sampling(6) = 0.5*(frequency_range(3)+frequency_range(4))
        frequency_range_sampling(7) = frequency_range(4);
        frequency_range_sampling(6) = 0.5*(frequency_range(4)+frequency_range(5))
        frequency_range_sampling(7) = frequency_range(5);


      endif           


 
 

         vp2 = (lambda + 2.d0*mu)/rho;  vs2 = mu/rho;
  
         Y_mu_nh(1,:) = 0.d0; Y_lambda_nh(1,:) = 0.d0;
         if (QP_nh .ne. 0.d0 .and. QS_nh .ne. 0.d0) then
 
           do i = 1, 2*N_SLS-1
              do j = 1, N_SLS

                A_QP(i,j) = (frequency_range(j)*frequency_range_sampling(i) + frequency_range(j)**2*(1.d0/QP_nh)) &
                       / (frequency_range(j)**2 + frequency_range_sampling(i)**2);

                A_QS(i,j) = (frequency_range(j)*frequency_range_sampling(i) + frequency_range(j)**2*(1.d0/QS_nh)) &
                      / (frequency_range(j)**2 + frequency_range_sampling(i)**2);
                      
              enddo
            
              RHS_QP(i) = 1.d0/QP_nh;
              RHS_QS(i) = 1.d0/QS_nh;
            
           enddo
         
           !call FACTORIZE_MATRIX(A_QP,N_SLS,PIVOT_QP)
           !call FACTORIZE_MATRIX(A_QS,N_SLS,PIVOT_QS)
         
           !call DIRECT_LU_SOLVER(A_QP,RHS_QP,N_SLS,YP,PIVOT_QP)
           !call DIRECT_LU_SOLVER(A_QS,RHS_QS,N_SLS,YS,PIVOT_QS)
                  
           call QR_SOLVE(2*N_SLS-1,N_SLS, A_QP, RHS_QP, YP) 
           call QR_SOLVE(2*N_SLS-1,N_SLS, A_QS, RHS_QS, YS) 
           
                             
           Y_mu_nh(1,:) = YS;
           Y_lambda_nh(1,:) = (vp2*YP - 2.d0*vs2*YS)/(vp2-2.d0*vs2);
      
         endif
     

!         do i = 1, N_SLS
!            if (Y_lambda_nh(1,i) .lt. 0.d0 ) then 
!              write(*,*) 'ERROR! ANELASTIC COEFFICIENTS HAVE TO BE POSITIVE!'
!              write(*,*) 'Y_LAMBDA = ', Y_lambda_nh(1,i)
!              write(*,*) 'SETTING Y_LAMBDA = 0'
!                Y_lambda_nh(1,i) = 0.d0
!            endif
!            if (Y_mu_nh(1,i) .lt. 0.d0) then
!              write(*,*) 'ERROR! ANELASTIC COEFFICIENTS HAVE TO BE POSITIVE!'
!              write(*,*) 'Y_MU = ', Y_mu_nh(1,i)
!              write(*,*) 'SETTING Y_MU = 0'
!                Y_mu_nh(1,i) = 0.d0           
!              write(*,*) '----------------------------------------------------'              
!            endif
!         enddo

         

       return
      
       end subroutine MAKE_ANELASTIC_COEFFICIENTS_NH
 
 

