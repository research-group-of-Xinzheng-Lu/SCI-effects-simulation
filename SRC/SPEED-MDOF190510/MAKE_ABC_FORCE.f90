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

!> @brief Computes ABC's. 
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nn number of 1D GLL nodes
!> @param[in] xq GLL nodes
!> @param[in] wq GLL weights
!> @param[in] dd spectral derivative matrix
!> @param[in] rho
!> @param[in] lambda
!> @param[in] mu material properties given node by node
!> @param[in] AA11  costant value for the bilinear map
!> @param[in] AA12  costant value for the bilinear map
!> @param[in] AA13  costant value for the bilinear map
!> @param[in] AA21  costant value for the bilinear map
!> @param[in] AA22  costant value for the bilinear map
!> @param[in] AA23  costant value for the bilinear map
!> @param[in] AA31  costant value for the bilinear map
!> @param[in] AA32  costant value for the bilinear map
!> @param[in] AA33  costant value for the bilinear map
!> @param[in] BB11  costant value for the bilinear map
!> @param[in] BB12  costant value for the bilinear map
!> @param[in] BB13  costant value for the bilinear map
!> @param[in] BB21  costant value for the bilinear map
!> @param[in] BB22  costant value for the bilinear map
!> @param[in] BB23  costant value for the bilinear map
!> @param[in] BB31  costant value for the bilinear map
!> @param[in] BB32  costant value for the bilinear map
!> @param[in] BB33  costant value for the bilinear map
!> @param[in] GG1  costant value for the bilinear map
!> @param[in] GG2  costant value for the bilinear map
!> @param[in] GG3  costant value for the bilinear map
!> @param[in] DD1  costant value for the bilinear map
!> @param[in] DD2  costant value for the bilinear map
!> @param[in] DD3  costant value for the bilinear map
!> @param[in] ia index for selecting the element absorbing face
!> @param[in] ib index for selecting the element absorbing face
!> @param[in] ja index for selecting the element absorbing face
!> @param[in] jb index for selecting the element absorbing face
!> @param[in] ka index for selecting the element absorbing face
!> @param[in] kb index for selecting the element absorbing face
!> @param[in] ux x-displacement
!> @param[in] uy y-displacement
!> @param[in] uz z-displacement
!> @param[in] vx x-velocity
!> @param[in] vy y-velocity
!> @param[in] vz z-velocity
!> @param[out] fx x-absorbing forces
!> @param[out] fy y-absorbing forces
!> @param[out] fz z-absorbing forces


      subroutine MAKE_ABC_FORCE(nn,xq,wq,dd,rho,lambda,mu,&
                                AA11,AA12,AA13,AA21,AA22,AA23,&
                                AA31,AA32,AA33,BB11,BB12,BB13,&
                                BB21,BB22,BB23,BB31,BB32,BB33,&
                                GG1,GG2,GG3,DD1,DD2,DD3,&
                                ia,ib,ja,jb,ka,kb,ux,uy,uz,vx,vy,vz,fx,fy,fz)
      
!
!************************************************************************************************** 
      
      implicit none
      
      integer*4 :: nn
      integer*4 :: ia,ib,ja,jb,ka,kb
      integer*4 :: i,j,k,p,q,r
      
      real*8 :: AA11,AA12,AA13,AA21,AA22,AA23,AA31,AA32,AA33
      real*8 :: BB11,BB12,BB13,BB21,BB22,BB23,BB31,BB32,BB33
      real*8 :: GG1,GG2,GG3,DD1,DD2,DD3
      real*8 :: dxdx,dxdy,dxdz,dydx,dydy,dydz,dzdx,dzdy,dzdz,det_j
      real*8 :: duxdx,duxdy,duxdz,duydx,duydy,duydz,duzdx,duzdy,duzdz
      real*8 :: dut1dt1,dut2dt2,dut3dt1,dut3dt2
      real*8 :: dxdu,dxdv,dydu,dydv,dzdu,dzdv
      real*8 :: s1,s2,s3,vt1,vt2,vt3,r8t
      real*8 :: t1x,t1y,t1z,t2x,t2y,t2z,t3x,t3y,t3z
      real*8 :: t1ux,t2ux,t3ux,t1uy,t2uy,t3uy,t1uz,t2uz,t3uz
           
      real*8, dimension(nn) :: xq,wq

      real*8, dimension(nn,nn) :: dd

      real*8, dimension(nn,nn,nn) :: rho,lambda,mu 
      real*8, dimension(nn,nn,nn) :: ux,uy,uz
      real*8, dimension(nn,nn,nn) :: vx,vy,vz
      real*8, dimension(nn,nn,nn) :: fx,fy,fz
      real*8, dimension(nn,nn,nn) :: alfa,beta 
          
      do r = 1,nn
         do q = 1,nn
            do p = 1,nn
               alfa(p,q,r) = dsqrt((lambda(p,q,r) + 2.0d0*mu(p,q,r))/rho(p,q,r))
               beta(p,q,r) = dsqrt(mu(p,q,r)/rho(p,q,r))
            enddo
         enddo
      enddo
      
      do r = 1,nn
         do q = 1,nn
            do p = 1,nn
               fx(p,q,r) = 0.0d0;  fy(p,q,r) = 0.0d0;   fz(p,q,r) = 0.0d0
            enddo
         enddo
      enddo
      
      do r = ka,kb
         do q = ja,jb
            do p = ia,ib
               t1ux = 0.d0; t1uy = 0.d0; t1uz = 0.d0
               t2ux = 0.d0; t2uy = 0.d0; t2uz = 0.d0
               t3ux = 0.d0; t3uy = 0.d0; t3uz = 0.d0
               
               dxdx = AA11 + BB12*xq(r) + BB13*xq(q) + GG1*xq(q)*xq(r)
               dydx = AA21 + BB22*xq(r) + BB23*xq(q) + GG2*xq(q)*xq(r)
               dzdx = AA31 + BB32*xq(r) + BB33*xq(q) + GG3*xq(q)*xq(r)
               
               dxdy = AA12 + BB11*xq(r) + BB13*xq(p) + GG1*xq(r)*xq(p)
               dydy = AA22 + BB21*xq(r) + BB23*xq(p) + GG2*xq(r)*xq(p)
               dzdy = AA32 + BB31*xq(r) + BB33*xq(p) + GG3*xq(r)*xq(p)
               
               dxdz = AA13 + BB11*xq(q) + BB12*xq(p) + GG1*xq(p)*xq(q)
               dydz = AA23 + BB21*xq(q) + BB22*xq(p) + GG2*xq(p)*xq(q)
               dzdz = AA33 + BB31*xq(q) + BB32*xq(p) + GG3*xq(p)*xq(q)
               
               det_j = dxdz * (dydx*dzdy - dzdx*dydy) &
                     - dydz * (dxdx*dzdy - dzdx*dxdy) &
                     + dzdz * (dxdx*dydy - dydx*dxdy)
               
               
               do i = 1,nn
                  t1ux = t1ux + ux(i,q,r) * dd(p,i)
                  t1uy = t1uy + uy(i,q,r) * dd(p,i)
                  t1uz = t1uz + uz(i,q,r) * dd(p,i)
               enddo
               
               do j = 1,nn
                  t2ux = t2ux + ux(p,j,r) * dd(q,j)
                  t2uy = t2uy + uy(p,j,r) * dd(q,j)
                  t2uz = t2uz + uz(p,j,r) * dd(q,j)
               enddo
               
               do k = 1,nn
                  t3ux = t3ux + ux(p,q,k) * dd(r,k)
                  t3uy = t3uy + uy(p,q,k) * dd(r,k)
                  t3uz = t3uz + uz(p,q,k) * dd(r,k)
               enddo
               
               
               duxdx = 1.0d0 / det_j *(&
                    ((dydy*dzdz - dydz*dzdy) * t1ux) + &
                    ((dydz*dzdx - dydx*dzdz) * t2ux) + &
                    ((dydx*dzdy - dydy*dzdx) * t3ux))
               
               duydx = 1.0d0 / det_j *(&
                    ((dydy*dzdz - dydz*dzdy) * t1uy) + &
                    ((dydz*dzdx - dydx*dzdz) * t2uy) + &
                    ((dydx*dzdy - dydy*dzdx) * t3uy))
               
               duzdx = 1.0d0 / det_j *(&
                    ((dydy*dzdz - dydz*dzdy) * t1uz) + &
                    ((dydz*dzdx - dydx*dzdz) * t2uz) + &
                    ((dydx*dzdy - dydy*dzdx) * t3uz))
               
               duxdy = 1.0d0 / det_j *(&
                    ((dzdy*dxdz - dzdz*dxdy) * t1ux) + &
                    ((dzdz*dxdx - dzdx*dxdz) * t2ux) + &
                    ((dzdx*dxdy - dzdy*dxdx) * t3ux))
               
               duydy = 1.0d0 / det_j *(&
                    ((dzdy*dxdz - dzdz*dxdy) * t1uy) + &
                    ((dzdz*dxdx - dzdx*dxdz) * t2uy) + &
                    ((dzdx*dxdy - dzdy*dxdx) * t3uy))
               
               duzdy = 1.0d0 / det_j *(&
                    ((dzdy*dxdz - dzdz*dxdy) * t1uz) + &
                    ((dzdz*dxdx - dzdx*dxdz) * t2uz) + &
                    ((dzdx*dxdy - dzdy*dxdx) * t3uz))
               
               duxdz = 1.0d0 / det_j *(&
                    ((dxdy*dydz - dxdz*dydy) * t1ux) + &
                    ((dxdz*dydx - dxdx*dydz) * t2ux) + &
                    ((dxdx*dydy - dxdy*dydx) * t3ux))
               
               duydz = 1.0d0 / det_j *(&
                    ((dxdy*dydz - dxdz*dydy) * t1uy) + &
                    ((dxdz*dydx - dxdx*dydz) * t2uy) + &
                    ((dxdx*dydy - dxdy*dydx) * t3uy))
               
               duzdz = 1.0d0 / det_j *(&
                    ((dxdy*dydz - dxdz*dydy) * t1uz) + &
                    ((dxdz*dydx - dxdx*dydz) * t2uz) + &
                    ((dxdx*dydy - dxdy*dydx) * t3uz))
               
               if ((ia.eq.ib).and.(ia.eq.1)) then
                  dxdu = dxdy
                  dydu = dydy
                  dzdu = dzdy
                  
                  dxdv = -dxdz
                  dydv = -dydz
                  dzdv = -dzdz
               endif
               
               if ((ja.eq.jb).and.(ja.eq.1)) then
                  dxdu = dxdz
                  dydu = dydz
                  dzdu = dzdz
                  
                  dxdv = -dxdx
                  dydv = -dydx
                  dzdv = -dzdx
               endif
               
               if ((ka.eq.kb).and.(ka.eq.1)) then
                  dxdu = dxdx
                  dydu = dydx
                  dzdu = dzdx
                  
                  dxdv = -dxdy
                  dydv = -dydy
                  dzdv = -dzdy
               endif
               
               if ((ia.eq.ib).and.(ia.eq.nn)) then
                  dxdu = dxdy
                  dydu = dydy
                  dzdu = dzdy
                  
                  dxdv = dxdz
                  dydv = dydz
                  dzdv = dzdz
               endif
               
               if ((ja.eq.jb).and.(ja.eq.nn)) then
                  dxdu = dxdz
                  dydu = dydz
                  dzdu = dzdz
                  
                  dxdv = dxdx
                  dydv = dydx
                  dzdv = dzdx
               endif
               
               if ((ka.eq.kb).and.(ka.eq.nn)) then
                  dxdu = dxdx
                  dydu = dydx
                  dzdu = dzdx
                  
                  dxdv = dxdy
                  dydv = dydy
                  dzdv = dzdy
               endif
               
               t1x = dxdu
               t1y = dydu
               t1z = dzdu
               
               t2x = dxdv
               t2y = dydv
               t2z = dzdv
               
               t3x = t1y*t2z - t1z*t2y
               t3y = t1z*t2x - t1x*t2z
               t3z = t1x*t2y - t1y*t2x
               
               r8t = dsqrt(t1x*t1x +t1y*t1y +t1z*t1z)
               t1x = t1x/r8t
               t1y = t1y/r8t
               t1z = t1z/r8t
               
               r8t = dsqrt(t3x*t3x +t3y*t3y +t3z*t3z)
               t3x = t3x/r8t
               t3y = t3y/r8t
               t3z = t3z/r8t
               
               t2x = t3y*t1z - t3z*t1y
               t2y = t3z*t1x - t3x*t1z
               t2z = t3x*t1y - t3y*t1x
               
               
               dut1dt1 = t1x*t1x*duxdx + t1x*t1y*duydx + t1x*t1z*duzdx &
                       + t1y*t1x*duxdy + t1y*t1y*duydy + t1y*t1z*duzdy &
                       + t1z*t1x*duxdz + t1z*t1y*duydz + t1z*t1z*duzdz
               
               dut2dt2 = t2x*t2x*duxdx + t2x*t2y*duydx + t2x*t2z*duzdx &
                       + t2y*t2x*duxdy + t2y*t2y*duydy + t2y*t2z*duzdy &
                       + t2z*t2x*duxdz + t2z*t2y*duydz + t2z*t2z*duzdz
               
               dut3dt1 = t1x*t3x*duxdx + t1x*t3y*duydx + t1x*t3z*duzdx &
                       + t1y*t3x*duxdy + t1y*t3y*duydy + t1y*t3z*duzdy &
                       + t1z*t3x*duxdz + t1z*t3y*duydz + t1z*t3z*duzdz
               
               dut3dt2 = t2x*t3x*duxdx + t2x*t3y*duydx + t2x*t3z*duzdx &
                       + t2y*t3x*duxdy + t2y*t3y*duydy + t2y*t3z*duzdy &
                       + t2z*t3x*duxdz + t2z*t3y*duydz + t2z*t3z*duzdz
               
               vt1 = t1x*vx(p,q,r) + t1y*vy(p,q,r) + t1z*vz(p,q,r)
               vt2 = t2x*vx(p,q,r) + t2y*vy(p,q,r) + t2z*vz(p,q,r)
               vt3 = t3x*vx(p,q,r) + t3y*vy(p,q,r) + t3z*vz(p,q,r)
               
               s1 = (mu(p,q,r)*(2.0d0*beta(p,q,r) -alfa(p,q,r))) &
                                /beta(p,q,r) * dut3dt1 - mu(p,q,r)/beta(p,q,r) * vt1
               s2 = (mu(p,q,r)*(2.0d0*beta(p,q,r) -alfa(p,q,r))) &
                                /beta(p,q,r) * dut3dt2 - mu(p,q,r)/beta(p,q,r) * vt2
               s3 = (lambda(p,q,r)*beta(p,q,r) +2.0d0*mu(p,q,r)* &
                                (beta(p,q,r) -alfa(p,q,r)))/alfa(p,q,r) &
                   * (dut1dt1 +dut2dt2) - (lambda(p,q,r) +2.0d0*mu(p,q,r))/alfa(p,q,r) * vt3
               
               det_j = dsqrt(((dxdu*dxdu +dydu*dydu +dzdu*dzdu) &
                     * (dxdv*dxdv +dydv*dydv +dzdv*dzdv)) &
                     -((dxdu*dxdv +dydu*dydv +dzdu*dzdv) &
                     * (dxdu*dxdv +dydu*dydv +dzdu*dzdv)))
                             
               fx(p,q,r) = fx(p,q,r) - wq(p)*wq(q)*wq(r)*det_j/wq(1) &
                         * (s1*t1x + s2*t2x + s3*t3x)
               fy(p,q,r) = fy(p,q,r) - wq(p)*wq(q)*wq(r)*det_j/wq(1) &
                         * (s1*t1y + s2*t2y + s3*t3y)
               fz(p,q,r) = fz(p,q,r) - wq(p)*wq(q)*wq(r)*det_j/wq(1) &
                         * (s1*t1z + s2*t2z + s3*t3z)
               
            enddo
         enddo
      enddo
      
      return
      
      end subroutine MAKE_ABC_FORCE

