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

!> @brief Make damping matrices for non-linear elastic materials.
!! @author Ilario Mazzieri
!> @date November, 2014 
!> @version 1.2
!> @param[in] nn nuber of 1D GLL nodes
!> @param[in] ct 1D GLL nodes
!> @param[in] ww 1D GLL weights
!> @param[in] dd spectral derivative matrix
!> @param[in] rho material density
!> @param[in] gamma0 damping coefficient
!> @param[in] alfa11 costant values for the bilinear map
!> @param[in] alfa12 costant values for the bilinear map
!> @param[in] alfa13 costant values for the bilinear map
!> @param[in] alfa21 costant values for the bilinear map
!> @param[in] alfa22 costant values for the bilinear map
!> @param[in] alfa23 costant values for the bilinear map
!> @param[in] alfa31 costant values for the bilinear map
!> @param[in] alfa32 costant values for the bilinear map
!> @param[in] alfa33 costant values for the bilinear map
!> @param[in] beta11 costant values for the bilinear map 
!> @param[in] beta12 costant values for the bilinear map
!> @param[in] beta13 costant values for the bilinear map
!> @param[in] beta21 costant values for the bilinear map
!> @param[in] beta22 costant values for the bilinear map
!> @param[in] beta23 costant values for the bilinear map
!> @param[in] beta31 costant values for the bilinear map
!> @param[in] beta32 costant values for the bilinear map
!> @param[in] beta33 costant values for the bilinear map
!> @param[in] gamma1 costant values for the bilinear map
!> @param[in] gamma2 costant values for the bilinear map
!> @param[in] gamma3 costant values for the bilinear map
!> @param[in] delta1 costant values for the bilinear map
!> @param[in] delta2 costant values for the bilinear map
!> @param[in] delta3 costant values for the bilinear map
!> @param[in] cs_nnz length of cs
!> @param[in] cs connectivity vector
!> @param[in] vcase  depth for non linear block 
!> @param[in] R_el  reduction factor associated to strain invariant 
!> @param[in] fpeak  peak frequency
!> @param[in] nf number of functions 
!> @param[in] func_type function type
!> @param[in] func_indx  indices for the data 
!> @param[in] nfdata number of data for each function
!> @param[in] func_data data for the calculation (depending on type) 
!> @param[in] t_stress  time
!> @param[in] tag_func  function label  
!> @param[in] yon  1 if the node is non-linear elastic, 0 otherwise
!> @param[in] tcase  not-honoring test case
!> @param[in] nnod_loc  number of local nodes
!> @param[in] vs_tria  vs values for each node
!> @param[out] mc_el non-linear damping matrix
!> @param[out] mck_el non-linear damping matrix 
!> @param[out] gamma_new gamma value for debug mode

      subroutine MAKE_DAMPING_MATRIX_NLE(nn,ct,ww,dd,rho,gamma0,&                                        
                           alfa11,alfa12,alfa13,alfa21,alfa22,alfa23,&        
                           alfa31,alfa32,alfa33,beta11,beta12,beta13,&        
                           beta21,beta22,beta23,beta31,beta32,beta33,&        
                           gamma1,gamma2,gamma3,delta1,delta2,delta3,&        
                           mc_el,mck_el,&
                           cs_nnz,cs,ielem,&                                
                           vcase,R_el,fpeak, &
                           func_type,func_indx,func_data,nfdata,nf,t_stress,tag_func,yon,gamma_new,&
                           tcase, nnod_loc, vs_tria)

      
      implicit none 
      
      integer*4 :: nn 
      integer*4 :: i,j,k 
      integer*4 :: nnt,tagmat                        
      integer*4 :: cs_nnz                                
      integer*4 :: vcase,nfdata                        
      integer*4 :: nf,fn,ielem,tcase,nnod_loc,in,is                                                                       

      integer*4, dimension(nf) :: func_type                                
      integer*4, dimension(nf +1) :: func_indx                        
      integer*4, dimension(0:cs_nnz) :: cs        
      integer*4, dimension(nf) :: tag_func                                
 
      integer*4, dimension(nn,nn,nn) :: yon

      real*8 :: alfa11,alfa12,alfa13,alfa21,alfa22,alfa23,alfa31,alfa32,alfa33 
      real*8 :: beta11,beta12,beta13,beta21,beta22,beta23,beta31,beta32,beta33 
      real*8 :: gamma1,gamma2,gamma3,delta1,delta2,delta3 
      real*8 :: dxdx,dxdy,dxdz,dydx,dydy,dydz,dzdx,dzdy,dzdz,det_j 
      real*8 :: DEPTH
      real*8 :: lambda, mu
      real*8 :: t_stress                                                                
      real*8 :: get_func_value        
      real*8 :: csi,fpeak,Q,PI,gamma

      real*8, dimension(nfdata) :: func_data                                        
      real*8, dimension(nn) :: ct,ww 
      real*8, dimension(nnod_loc) :: vs_tria 
            
      real*8, dimension(nn,nn) :: dd 
      
      real*8, dimension(nn,nn,nn) :: rho,gamma0,gamma_new 
      real*8, dimension(nn,nn,nn) :: mc_el,mck_el 
      real*8, dimension(nn,nn,nn) :: R_el

      
      gamma_new = gamma0;
      
!     ELEMENT NODAL MASS CALCULATION 
      
      do k = 1,nn 
         do j = 1,nn 
            do i = 1,nn

                if (yon(i,j,k).eq.1) then

                   dxdx = alfa11 +beta12*ct(k) +beta13*ct(j) & 
                                        + gamma1*ct(j)*ct(k) 
                   dydx = alfa21 +beta22*ct(k) +beta23*ct(j) & 
                                        + gamma2*ct(j)*ct(k) 
                   dzdx = alfa31 +beta32*ct(k) +beta33*ct(j) & 
                                        + gamma3*ct(j)*ct(k) 
                   dxdy = alfa12 +beta11*ct(k) +beta13*ct(i) & 
                                       + gamma1*ct(k)*ct(i) 
                   dydy = alfa22 +beta21*ct(k) +beta23*ct(i) & 
                                        + gamma2*ct(k)*ct(i) 
                   dzdy = alfa32 +beta31*ct(k) +beta33*ct(i) & 
                                        + gamma3*ct(k)*ct(i) 
                   dxdz = alfa13 +beta11*ct(j) +beta12*ct(i) & 
                                        + gamma1*ct(i)*ct(j) 
                   dydz = alfa23 +beta21*ct(j) +beta22*ct(i) & 
                                        + gamma2*ct(i)*ct(j) 
                   dzdz = alfa33 +beta31*ct(j) +beta32*ct(i) & 
                                        + gamma3*ct(i)*ct(j) 
                   det_j = dxdz * (dydx*dzdy - dzdx*dydy) & 
                                        - dydz * (dxdx*dzdy - dzdx*dxdy) & 
                                        + dzdz * (dxdx*dydy - dydx*dxdy) 

                  is = nn*nn*(k -1) +nn*(j -1) +i        
                  in = cs(cs(ielem -1) +is)
                  !-------------------------------------------------------------------
                  !
                  if(tcase .ne. 16) then
                    do fn = 1,nf                                                                
                       if (vcase.eq.tag_func(fn)) then    
                          if (func_type(fn) .eq. 61) then                                                 
                             csi = get_func_value(nf,func_type,func_indx,func_data, nfdata,&  
                                                 fn,R_el(i,j,k)) 
                                                                  
                                                                                                       
                              Q = 1.d0/(2.d0*csi)
                              PI = 4.0d0 * datan(1.0d0)
                              gamma = PI * fpeak / Q
                                                
 
                              if (gamma0(i,j,k) .ge. gamma) gamma = gamma0(i,j,k)
                              gamma_new(i,j,k) = gamma                
                                                                
                              mc_el(i,j,k) = 2.d0 * gamma * rho(i,j,k) &
                                                                * det_j * ww(i) * ww(j) * ww(k)      
                              mck_el(i,j,k) = (gamma**2) * rho(i,j,k) & 
                                                                * det_j * ww(i) * ww(j) * ww(k)    
                                 
                           endif                                                                        
                       endif                                                                                
                     enddo
                   elseif (tcase .eq. 16) then
                     
                     do fn = 1,nf                        
                         if (vcase .eq. tag_func(fn)) then                                                            
                             if (func_type(fn) .eq. 61 .and. vs_tria(in) .lt. 325.d0) then     

                               csi = get_func_value(nf,func_type,func_indx,func_data, nfdata,&  
                                                 fn,R_el(i,j,k)) 
                                                                  
                                                                                                       
                                Q = 1.d0/(2.d0*csi)
                                PI = 4.0d0 * datan(1.0d0)
                                gamma = PI * fpeak / Q
                                                
 
                                if (gamma0(i,j,k) .ge. gamma) gamma = gamma0(i,j,k)
                                gamma_new(i,j,k) = gamma                
                                                                
                                mc_el(i,j,k) = 2.d0 * gamma * rho(i,j,k) &
                                                                  * det_j * ww(i) * ww(j) * ww(k)      
                                mck_el(i,j,k) = (gamma**2) * rho(i,j,k) & 
                                                                  * det_j * ww(i) * ww(j) * ww(k)    


                             elseif (func_type(fn) .eq. 63 .and. vs_tria(in) .lt. 450.d0) then     
                               csi = get_func_value(nf,func_type,func_indx,func_data, nfdata,&  
                                                 fn,R_el(i,j,k)) 
                                                                  
                                                                                                       
                                Q = 1.d0/(2.d0*csi)
                                PI = 4.0d0 * datan(1.0d0)
                                gamma = PI * fpeak / Q
                                                
 
                                if (gamma0(i,j,k) .ge. gamma) gamma = gamma0(i,j,k)
                                gamma_new(i,j,k) = gamma                
                                                                
                                mc_el(i,j,k) = 2.d0 * gamma * rho(i,j,k) &
                                                                * det_j * ww(i) * ww(j) * ww(k)      
                                mck_el(i,j,k) = (gamma**2) * rho(i,j,k) & 
                                                                * det_j * ww(i) * ww(j) * ww(k)    
                               
                                                                                                                                    
                              endif     
                          endif                                                                     
                        enddo
                     
                   endif!tcase = 16
           
                        !
                        !-------------------------------------------------------------------

                endif !if (yon(i,j,k).eq.1) then

            enddo 
         enddo 
      enddo 

      
      return 
      
      end subroutine MAKE_DAMPING_MATRIX_NLE 

