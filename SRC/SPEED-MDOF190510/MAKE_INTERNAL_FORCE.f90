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
!    along with SPEED.  If not, see <http://wqw.gnu.org/licenses/>.

!> @brief Makes internal forces.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nn number of 1-D GLL nodes
!> @param[in] xq GLL nodes
!> @param[in] wq GLL weights
!> @param[in] dd spectral derivative matrix
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
!> @param[in] sxx nodal values for the stress tensor
!> @param[in] syy nodal values for the stress tensor
!> @param[in] szz nodal values for the stress tensor
!> @param[in] syz nodal values for the stress tensor
!> @param[in] szx nodal values for the stress tensor
!> @param[in] sxy nodal values for the stress tensor
!> @param[out] fx x-componnent for internal forces
!> @param[out] fy y-componnent for internal forces
!> @param[out] fz z-componnent for internal forces

      subroutine MAKE_INTERNAL_FORCE(nn,xq,wq,dd,&
                               AA11,AA12,AA13,AA21,AA22,AA23,&
                               AA31,AA32,AA33,BB11,BB12,BB13,&
                               BB21,BB22,BB23,BB31,BB32,BB33,&
                               GG1,GG2,GG3,DD1,DD2,DD3,&
                               sxx,syy,szz,syz,szx,sxy,fx,fy,fz)
      
      implicit none
      
      integer*4 :: nn
      integer*4 :: i,j,k,p,q,r

      real*8 :: AA11,AA12,AA13,AA21,AA22,AA23,AA31,AA32,AA33
      real*8 :: BB11,BB12,BB13,BB21,BB22,BB23,BB31,BB32,BB33
      real*8 :: GG1,GG2,GG3,DD1,DD2,DD3
      real*8 :: dxdx,dxdy,dxdz,dydx,dydy,dydz,dzdx,dzdy,dzdz,det_j
      real*8 :: duxdx,duxdy,duxdz,duydx,duydy,duydz,duzdx,duzdy,duzdz
      real*8 :: t1ux,t2ux,t3ux,t1uy,t2uy,t3uy,t1uz,t2uz,t3uz

      real*8, dimension(nn) :: xq,wq

      real*8, dimension(nn,nn) :: dd

      real*8, dimension(nn,nn,nn) :: sxx,syy,szz,syz,szx,sxy
      real*8, dimension(nn,nn,nn) :: fx,fy,fz
      
!     FORCE CALCULATION
      
      do r = 1,nn
         do q = 1,nn
            do p = 1,nn
               t1ux = 0.d0;   t1uy = 0.d0;   t1uz = 0.d0
               t2ux = 0.d0;   t2uy = 0.d0;   t2uz = 0.d0
               t3ux = 0.d0;   t3uy = 0.d0;   t3uz = 0.d0
               
               do i = 1,nn
                  dxdy = AA12 + BB11*xq(r) + BB13*xq(i) + GG1*xq(r)*xq(i)
                  dydy = AA22 + BB21*xq(r) + BB23*xq(i) + GG2*xq(r)*xq(i)
                  dzdy = AA32 + BB31*xq(r) + BB33*xq(i) + GG3*xq(r)*xq(i)
                  
                  dxdz = AA13 + BB11*xq(q) + BB12*xq(i) + GG1*xq(i)*xq(q)
                  dydz = AA23 + BB21*xq(q) + BB22*xq(i) + GG2*xq(i)*xq(q)
                  dzdz = AA33 + BB31*xq(q) + BB32*xq(i) + GG3*xq(i)*xq(q)
                  
                  t1ux = t1ux + wq(i)*wq(q)*wq(r) * dd(i,p) * &
                         ((dydy*dzdz - dydz*dzdy) * sxx(i,q,r) + &
                          (dzdy*dxdz - dzdz*dxdy) * sxy(i,q,r) + &
                          (dxdy*dydz - dxdz*dydy) * szx(i,q,r))
                  t1uy = t1uy + wq(i)*wq(q)*wq(r) * dd(i,p) * &
                         ((dydy*dzdz - dydz*dzdy) * sxy(i,q,r) + &
                          (dzdy*dxdz - dzdz*dxdy) * syy(i,q,r) + &
                          (dxdy*dydz - dxdz*dydy) * syz(i,q,r))
                  t1uz = t1uz + wq(i)*wq(q)*wq(r) * dd(i,p) * &
                         ((dydy*dzdz - dydz*dzdy) * szx(i,q,r) + &
                          (dzdy*dxdz - dzdz*dxdy) * syz(i,q,r) + &
                          (dxdy*dydz - dxdz*dydy) * szz(i,q,r))
               enddo
               
               do j = 1,nn
                  dxdx = AA11 + BB12*xq(r) + BB13*xq(j) + GG1*xq(j)*xq(r)
                  dydx = AA21 + BB22*xq(r) + BB23*xq(j) + GG2*xq(j)*xq(r)
                  dzdx = AA31 + BB32*xq(r) + BB33*xq(j) + GG3*xq(j)*xq(r)
                  
                  dxdz = AA13 + BB11*xq(j) + BB12*xq(p) + GG1*xq(p)*xq(j)
                  dydz = AA23 + BB21*xq(j) + BB22*xq(p) + GG2*xq(p)*xq(j)
                  dzdz = AA33 + BB31*xq(j) + BB32*xq(p) + GG3*xq(p)*xq(j)
                  
                  t2ux = t2ux + wq(p)*wq(j)*wq(r) * dd(j,q) * &
                         ((dydz*dzdx - dydx*dzdz) * sxx(p,j,r) + &
                          (dzdz*dxdx - dzdx*dxdz) * sxy(p,j,r) + &
                          (dxdz*dydx - dxdx*dydz) * szx(p,j,r))
                  t2uy = t2uy + wq(p)*wq(j)*wq(r) * dd(j,q) * &
                         ((dydz*dzdx - dydx*dzdz) * sxy(p,j,r) + &
                          (dzdz*dxdx - dzdx*dxdz) * syy(p,j,r) + &
                          (dxdz*dydx - dxdx*dydz) * syz(p,j,r))
                  t2uz = t2uz + wq(p)*wq(j)*wq(r) * dd(j,q) * &
                         ((dydz*dzdx - dydx*dzdz) * szx(p,j,r) + &
                          (dzdz*dxdx - dzdx*dxdz) * syz(p,j,r) + &
                          (dxdz*dydx - dxdx*dydz) * szz(p,j,r))
               enddo
               
               do k = 1,nn
                  dxdx = AA11 + BB12*xq(k) + BB13*xq(q) + GG1*xq(q)*xq(k)
                  dydx = AA21 + BB22*xq(k) + BB23*xq(q) + GG2*xq(q)*xq(k)
                  dzdx = AA31 + BB32*xq(k) + BB33*xq(q) + GG3*xq(q)*xq(k)
                  
                  dxdy = AA12 + BB11*xq(k) + BB13*xq(p) + GG1*xq(k)*xq(p)
                  dydy = AA22 + BB21*xq(k) + BB23*xq(p) + GG2*xq(k)*xq(p)
                  dzdy = AA32 + BB31*xq(k) + BB33*xq(p) + GG3*xq(k)*xq(p)
                  
                  t3ux = t3ux + wq(p)*wq(q)*wq(k) * dd(k,r) * &
                         ((dydx*dzdy - dydy*dzdx) * sxx(p,q,k) + &
                          (dzdx*dxdy - dzdy*dxdx) * sxy(p,q,k) + &
                          (dxdx*dydy - dxdy*dydx) * szx(p,q,k))
                  t3uy = t3uy + wq(p)*wq(q)*wq(k) * dd(k,r) * &
                         ((dydx*dzdy - dydy*dzdx) * sxy(p,q,k) + &
                          (dzdx*dxdy - dzdy*dxdx) * syy(p,q,k) + &
                          (dxdx*dydy - dxdy*dydx) * syz(p,q,k))
                  t3uz = t3uz + wq(p)*wq(q)*wq(k) * dd(k,r) * &
                         ((dydx*dzdy - dydy*dzdx) * szx(p,q,k) + &
                          (dzdx*dxdy - dzdy*dxdx) * syz(p,q,k) + &
                          (dxdx*dydy - dxdy*dydx) * szz(p,q,k))
               enddo
               
               fx(p,q,r) = t1ux + t2ux + t3ux
               fy(p,q,r) = t1uy + t2uy + t3uy
               fz(p,q,r) = t1uz + t2uz + t3uz
            enddo
         enddo
      enddo
      
      
      return
      
      end subroutine MAKE_INTERNAL_FORCE

