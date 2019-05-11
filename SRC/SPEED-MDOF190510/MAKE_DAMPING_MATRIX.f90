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

!> @brief Make damping matrices.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nn number of 1-D GLL nodes
!> @param[in] xq GLL nodes
!> @param[in] wq GLL weights  
!> @param[in] dd spectral derivatives matrix
!> @param[in] rho mass density
!> @param[in] gamma damping coefficient
!> @param[in] A11 costant values for the bilinear map
!> @param[in] A12 costant values for the bilinear map
!> @param[in] A13 costant values for the bilinear map
!> @param[in] A21 costant values for the bilinear map
!> @param[in] A22 costant values for the bilinear map
!> @param[in] A23 costant values for the bilinear map
!> @param[in] A31 costant values for the bilinear map
!> @param[in] A32 costant values for the bilinear map
!> @param[in] A33 costant values for the bilinear map
!> @param[in] B11 costant values for the bilinear map 
!> @param[in] B12 costant values for the bilinear map
!> @param[in] B13 costant values for the bilinear map
!> @param[in] B21 costant values for the bilinear map
!> @param[in] B22 costant values for the bilinear map
!> @param[in] B23 costant values for the bilinear map
!> @param[in] B31 costant values for the bilinear map
!> @param[in] B32 costant values for the bilinear map
!> @param[in] B33 costant values for the bilinear map
!> @param[in] GG1 costant values for the bilinear map
!> @param[in] GG2 costant values for the bilinear map
!> @param[in] GG3 costant values for the bilinear map
!> @param[in] DD1 costant values for the bilinear map
!> @param[in] DD2 costant values for the bilinear map
!> @param[in] DD3 costant values for the bilinear map
!> @param[out] mc_el viscous forces proportional to the velocity field
!> @param[out] mck_el viscous forces proportional to the displacement field


      subroutine MAKE_DAMPING_MATRIX(nn,ct,ww,dd,rho,gamma,&                                        
                           A11,A12,A13,A21,A22,A23,&        
                           A31,A32,A33,B11,B12,B13,&        
                           B21,B22,B23,B31,B32,B33,&        
                           GG1,GG2,GG3,DD1,DD2,DD3,&        
                           mc_el,mck_el)
      
 
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

                    mc_el(i,j,k) = 2* gamma(i,j,k) * rho(i,j,k) &                                               
                                 * det_j * ww(i) * ww(j) * ww(k) 
                                   
                   mck_el(i,j,k) = (gamma(i,j,k)**2) * rho(i,j,k) & 
                                * det_j * ww(i) * ww(j) * ww(k) 


            enddo 
         enddo 
      enddo 
      
      return 
      
      end subroutine MAKE_DAMPING_MATRIX 

