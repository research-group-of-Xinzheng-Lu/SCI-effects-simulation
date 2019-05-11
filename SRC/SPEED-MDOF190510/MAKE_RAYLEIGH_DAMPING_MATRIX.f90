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

!> @brief Compute the Rayleigh damping matrix C = A0 M + A1 K,
!! M = mass matrix, K= stiffness matrix
!! @author Ilario Mazzieri
!> @date November, 2014 
!> @version 1.0
!> @param[in] nn polynomial degree + 1
!> @param[in] ct LGL nodes
!> @param[in] ww LGL weights
!> @param[in] dd matrix of spectral derivates 
!> @param[in] rho mass density
!> @param[in] A11,...,DD3 coefficients for the bilinear map
!> @param[out] mc_el Rayleigh damping matrix

      subroutine MAKE_RAYLEIGH_DAMPING_MATRIX(nn,ct,ww,dd,rho,&                                        
                           A11,A12,A13,A21,A22,A23,&        
                           A31,A32,A33,B11,B12,B13,&        
                           B21,B22,B23,B31,B32,B33,&        
                           GG1,GG2,GG3,DD1,DD2,DD3,&        
                           mc_el)
      
 
      implicit none 
      
      integer*4 :: nn 
      integer*4 :: i,j,k 

      real*8 :: A11,A12,A13,A21,A22,A23,A31,A32,A33 
      real*8 :: B11,B12,B13,B21,B22,B23,B31,B32,B33 
      real*8 :: GG1,GG2,GG3,DD1,DD2,DD3 
      real*8 :: dxdx,dxdy,dxdz,dydx,dydy,dydz,dzdx,dzdy,dzdz,det_j 

      real*8, dimension(nn) :: ct,ww 

      real*8, dimension(nn,nn) :: dd 

      real*8, dimension(nn,nn,nn) :: rho,gamma
      real*8, dimension(nn,nn,nn) :: mc_el,mck_el 

     
      
      do k = 1,nn 
         do j = 1,nn 
            do i = 1,nn 
               dxdx = A11 +B12*ct(k) +B13*ct(j) & 
                    + GG1*ct(j)*ct(k) 
               dydx = A21 +B22*ct(k) +B23*ct(j) & 
                    + GG2*ct(j)*ct(k) 
               dzdx = A31 +B32*ct(k) +B33*ct(j) & 
                    + GG3*ct(j)*ct(k) 
               
               dxdy = A12 +B11*ct(k) +B13*ct(i) & 
                    + GG1*ct(k)*ct(i) 
               dydy = A22 +B21*ct(k) +B23*ct(i) & 
                    + GG2*ct(k)*ct(i) 
               dzdy = A32 +B31*ct(k) +B33*ct(i) & 
                    + GG3*ct(k)*ct(i) 
               
               dxdz = A13 +B11*ct(j) +B12*ct(i) & 
                    + GG1*ct(i)*ct(j) 
               dydz = A23 +B21*ct(j) +B22*ct(i) & 
                    + GG2*ct(i)*ct(j) 
               dzdz = A33 +B31*ct(j) +B32*ct(i) & 
                    + GG3*ct(i)*ct(j) 
               
               det_j = dxdz * (dydx*dzdy - dzdx*dydy) & 
                     - dydz * (dxdx*dzdy - dzdx*dxdy) & 
                     + dzdz * (dxdx*dydy - dydx*dxdy) 


             
                    mc_el(i,j,k) = rho(i,j,k) &                                                
                                 * det_j * ww(i) * ww(j) * ww(k) 
                                   
                   

            enddo 
         enddo 
      enddo 
      
      return 
      
      end subroutine MAKE_RAYLEIGH_DAMPING_MATRIX

