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

!> @brief Makes mass matrix
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nn number of 1-D GLL nodes
!> @param[in] xq GLL nodes
!> @param[in] wq GLL weights  
!> @param[in] dd spectral derivatives matrix
!> @param[in] rho mass density
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
!> @param[out] mass  diagonal mass matrix 

      subroutine MAKE_MASS_MATRIX(nn,xq,wq,dd,rho,&
                           AA11,AA12,AA13,AA21,AA22,AA23,&
                           AA31,AA32,AA33,BB11,BB12,BB13,&
                           BB21,BB22,BB23,BB31,BB32,BB33,&
                           GG1,GG2,GG3,DD1,DD2,DD3,&
                           mass)
      
      
      implicit none
      
      integer*4 :: nn
      integer*4 :: i,j,k

      real*8 :: AA11,AA12,AA13,AA21,AA22,AA23,AA31,AA32,AA33
      real*8 :: BB11,BB12,BB13,BB21,BB22,BB23,BB31,BB32,BB33
      real*8 :: GG1,GG2,GG3,DD1,DD2,DD3
      real*8 :: dxdx,dxdy,dxdz,dydx,dydy,dydz,dzdx,dzdy,dzdz,det_j

      real*8, dimension(nn) :: xq, wq

      real*8, dimension(nn,nn) :: dd

      real*8, dimension(nn,nn,nn) :: rho 
      real*8, dimension(nn,nn,nn) :: mass
      
      
!     ELEMENT NODAL MASS CALCULATION
      
      do k = 1,nn
         do j = 1,nn
            do i = 1,nn
               dxdx = AA11 + BB12*xq(k) + BB13*xq(j) + GG1*xq(j)*xq(k)
               dydx = AA21 + BB22*xq(k) + BB23*xq(j) + GG2*xq(j)*xq(k)
               dzdx = AA31 + BB32*xq(k) + BB33*xq(j) + GG3*xq(j)*xq(k)
               
               dxdy = AA12 + BB11*xq(k) + BB13*xq(i) + GG1*xq(k)*xq(i)
               dydy = AA22 + BB21*xq(k) + BB23*xq(i) + GG2*xq(k)*xq(i)
               dzdy = AA32 + BB31*xq(k) + BB33*xq(i) + GG3*xq(k)*xq(i)
               
               dxdz = AA13 + BB11*xq(j) + BB12*xq(i) + GG1*xq(i)*xq(j)
               dydz = AA23 + BB21*xq(j) + BB22*xq(i) + GG2*xq(i)*xq(j)
               dzdz = AA33 + BB31*xq(j) + BB32*xq(i) + GG3*xq(i)*xq(j)
               
               det_j = dxdz * (dydx*dzdy - dzdx*dydy) &
                     - dydz * (dxdx*dzdy - dzdx*dxdy) &
                     + dzdz * (dxdx*dydy - dydx*dxdy)
               
               mass(i,j,k) = rho(i,j,k) * det_j * wq(i)*wq(j)*wq(k)
            enddo
         enddo
      enddo
      
      return
      
      end subroutine MAKE_MASS_MATRIX

