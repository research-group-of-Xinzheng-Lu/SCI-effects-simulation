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

!> @brief Findes node coordinates in the reference element.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nb_el number of elements
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
!> @param[in] x_true x-coordinate of the node 
!> @param[in] y_true y-coordinate of the node 
!> @param[in] z_true z-coordinate of the node 
!> @param[in] highest  elevation of the elemnt
!> @param[in] elevation  maximal elevation 
!> @param[out] el  element index
!> @param[out] xref  x-coordinate in the refernce element
!> @param[out] yref  y-coordinate in the refernce element
!> @param[out] zref  z-coordinate in the refernce element

      subroutine GET_PNT_POS_PGM(nb_el,&
                                 AA11,AA12,AA13,AA21,AA22,AA23,&
                                 AA31,AA32,AA33,BB11,BB12,BB13,&
                                 BB21,BB22,BB23,BB31,BB32,BB33,&
                                 GG1,GG2,GG3,DD1,DD2,DD3,&
                                 x_true,y_true,z_true,el,xref,yref,zref,&
                                 highest,elevation)
      
      implicit none
      
      integer*4 :: nb_el,el
      integer*4 :: inr,ie

      real*8 :: x_true,y_true,z_true,xref,yref,zref
      real*8 :: csi,eta,zeta,delta_csi,delta_eta,delta_zeta,dist,dist_min
      real*8 :: termX,termY,termZ,det_j
      real*8 :: dxdx,dxdy,dxdz
      real*8 :: dydx,dydy,dydz
      real*8 :: dzdx,dzdy,dzdz
      real*8 :: a11,a12,a13
      real*8 :: a21,a22,a23
      real*8 :: a31,a32,a33
      real*8 :: elevation

      real*8, dimension(nb_el) :: AA11,AA12,AA13
      real*8, dimension(nb_el) :: AA21,AA22,AA23
      real*8, dimension(nb_el) :: AA31,AA32,AA33
      real*8, dimension(nb_el) :: BB11,BB12,BB13
      real*8, dimension(nb_el) :: BB21,BB22,BB23
      real*8, dimension(nb_el) :: BB31,BB32,BB33
      real*8, dimension(nb_el) :: GG1,GG2,GG3
      real*8, dimension(nb_el) :: DD1,DD2,DD3
      real*8, dimension(nb_el) :: highest
      
      dist_min = 1.0E20
      xref = 0.0;      yref = 0.0;      zref = 0.0
      el = 1
      
      do ie = 1,nb_el
         if (highest(ie) .gt. elevation) then
            csi = 0.0;             eta = 0.0;     zeta = 0.0
            
            do inr = 1, 100
                dxdx = AA11(ie) + BB12(ie)*zeta + BB13(ie)*eta + GG1(ie)*eta*zeta
                dydx = AA21(ie) + BB22(ie)*zeta + BB23(ie)*eta + GG2(ie)*eta*zeta
                dzdx = AA31(ie) + BB32(ie)*zeta + BB33(ie)*eta + GG3(ie)*eta*zeta
          
                dxdy = AA12(ie) + BB11(ie)*zeta + BB13(ie)*csi + GG1(ie)*zeta*csi
                dydy = AA22(ie) + BB21(ie)*zeta + BB23(ie)*csi + GG2(ie)*zeta*csi
                dzdy = AA32(ie) + BB31(ie)*zeta + BB33(ie)*csi + GG3(ie)*zeta*csi
          
                dxdz = AA13(ie) + BB11(ie)*eta + BB12(ie)*csi  + GG1(ie)*csi*eta
                dydz = AA23(ie) + BB21(ie)*eta + BB22(ie)*csi  + GG2(ie)*csi*eta
                dzdz = AA33(ie) + BB31(ie)*eta + BB32(ie)*csi  + GG3(ie)*csi*eta
          
                det_j = dxdz * (dydx*dzdy - dzdx*dydy) - dydz * (dxdx*dzdy - dzdx*dxdy) &
                          + dzdz * (dxdx*dydy - dydx*dxdy)
          
          
                a11 = +dydy*dzdz -dydz*dzdy
                a12 = -dydx*dzdz +dydz*dzdx
                a13 = +dydx*dzdy -dydy*dzdx
          
                a21 = -dxdy*dzdz +dxdz*dzdy
                a22 = +dxdx*dzdz -dxdz*dzdx
                a23 = -dxdx*dzdy +dxdy*dzdx
          
                a31 = +dxdy*dydz -dxdz*dydy
                a32 = -dxdx*dydz +dxdz*dydx
                a33 = +dxdx*dydy -dxdy*dydx
          
                termX = AA11(ie)*csi + AA12(ie)*eta &
                   + AA13(ie)*zeta + BB11(ie)*eta*zeta &
                   + BB12(ie)*csi*zeta + BB13(ie)*csi*eta &
                   + GG1(ie)*csi*eta*zeta + DD1(ie) -x_true
          
                termY = AA21(ie)*csi + AA22(ie)*eta &
                   + AA23(ie)*zeta + BB21(ie)*eta*zeta &
                   + BB22(ie)*csi*zeta + BB23(ie)*csi*eta &
                   + GG2(ie)*csi*eta*zeta + DD2(ie) -y_true
            
                termZ = AA31(ie)*csi + AA32(ie)*eta &
                   + AA33(ie)*zeta + BB31(ie)*eta*zeta &
                   + BB32(ie)*csi*zeta + BB33(ie)*csi*eta &
                   + GG3(ie)*csi*eta*zeta + DD3(ie) -z_true
            
                    delta_csi  = -(a11*termX +a21*termY +a31*termZ) / det_j
                delta_eta  = -(a12*termX +a22*termY +a32*termZ) / det_j
                delta_zeta = -(a13*termX +a23*termY +a33*termZ) / det_j
            
                csi = csi + delta_csi
                eta = eta + delta_eta
                zeta = zeta + delta_zeta
                    
                    if (csi.lt.-1.0) csi = -1.0
                if (csi.gt. 1.0) csi =  1.0
                if (eta.lt.-1.0) eta = -1.0
                if (eta.gt. 1.0) eta =  1.0
                if (zeta.lt.-1.0) zeta = -1.0
                if (zeta.gt. 1.0) zeta =  1.0
            
            enddo
         
                 dist = termX*termX +termY*termY +termZ*termZ
         
                 if (dist.lt.dist_min) then
                    el = ie
                    xref = csi
                    yref = eta
                    zref = zeta
                    dist_min = dist
                 endif

        endif !if (highest(ie).gt.elevation) then
         
      enddo !do ie = 1,nb_el
      
      return
      
      end subroutine GET_PNT_POS_PGM

