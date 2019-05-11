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

!> @brief Computes nodal values for ratiotional, strain and stress tensor.
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] nn number of 1-D Legendre nodes
!> @param[in] ct GLL nodes
!> @param[in] ww GLL weights
!> @param[in] dd matrix of spectral derivatives
!> @param[in] alfa11 costant for the bilinear map
!> @param[in] alfa12 costant for the bilinear map
!> @param[in] alfa13 costant for the bilinear map
!> @param[in] alfa21 costant for the bilinear map
!> @param[in] alfa22 costant for the bilinear map
!> @param[in] alfa23 costant for the bilinear map
!> @param[in] alfa31 costant for the bilinear map
!> @param[in] alfa32 costant for the bilinear map
!> @param[in] alfa33 costant for the bilinear map
!> @param[in] beta11 costant for the bilinear map
!> @param[in] beta12 costant for the bilinear map
!> @param[in] beta13 costant for the bilinear map
!> @param[in] beta21 costant for the bilinear map
!> @param[in] beta22 costant for the bilinear map
!> @param[in] beta23 costant for the bilinear map
!> @param[in] beta31 costant for the bilinear map
!> @param[in] beta32 costant for the bilinear map
!> @param[in] beta33 costant for the bilinear map
!> @param[in] gamma1  costant for the bilinear map
!> @param[in] gamma2  costant for the bilinear map
!> @param[in] gamma3  costant for the bilinear map
!> @param[in] delta1  costant for the bilinear map
!> @param[in] delta2  costant for the bilinear map
!> @param[in] delta3  costant for the bilinear map
!> @param[in] ux x-displacement
!> @param[in] uy y-displacement
!> @param[in] uy z-displacement
!> @param[in] div duxdx + duydy + duzdz
!> @param[in] lambda nodal values of Lame coefficient lambda  
!> @param[in] mu nodal values of Lame coefficient mu
!> @param[out] rotx  nodal values of rotational tensor x-axis
!> @param[out] roty  nodal values of rotational tensor y-axis 
!> @param[out] rotz  nodal values of rotational tensor z-axis
!> @param[out] strainxx nodal values of strain tensor 
!> @param[out] strainyy nodal values of strain tensor
!> @param[out] strainzz nodal values of strain tensor
!> @param[out] strainxy nodal values of strain tensor
!> @param[out] strainyz nodal values of strain tensor
!> @param[out] strainzx nodal values of strain tensor
!> @param[out] stressxx nodal values of stress tensor
!> @param[out] stressyy nodal values of stress tensor
!> @param[out] stresszz nodal values of stress tensor
!> @param[out] stressxy nodal values of stress tensor
!> @param[out] stressyz nodal values of stress tensor
!> @param[out] stresszx nodal values of stress tensor

      subroutine MAKE_STRAIN_ROT_STRESS_TENSOR(nn,ct,ww,dd,&
                                        alfa11,alfa12,alfa13,alfa21,alfa22,alfa23,&
                                        alfa31,alfa32,alfa33,beta11,beta12,beta13,&
                                        beta21,beta22,beta23,beta31,beta32,beta33,&
                                        gamma1,gamma2,gamma3,delta1,delta2,delta3,&
                                        ux,uy,uz,div,&
                                        lambda,mu,&
                                        rotx,roty,rotz,&
                                        strainxx,strainyy,strainzz,&
                                        strainxy,strainyz,strainzx,&
                                        stressxx,stressyy,stresszz,&
                                        stressxy,stressyz,stresszx)
      

    
      implicit none
      
      integer*4 :: nn
      integer*4 :: i,j,k,p,q,r,L

      real*8 :: alfa11,alfa12,alfa13,alfa21,alfa22,alfa23,alfa31,alfa32,alfa33
      real*8 :: beta11,beta12,beta13,beta21,beta22,beta23,beta31,beta32,beta33
      real*8 :: gamma1,gamma2,gamma3,delta1,delta2,delta3,dphi,phi
      real*8 :: dxdx,dxdy,dxdz,dydx,dydy,dydz,dzdx,dzdy,dzdz,det_j
      real*8 :: duxdx,duxdy,duxdz,duydx,duydy,duydz,duzdx,duzdy,duzdz
      real*8 :: t1ux,t2ux,t3ux,t1uy,t2uy,t3uy,t1uz,t2uz,t3uz, dpdcsi, dpdeta, dpdzeta
      
      real*8, dimension(nn) :: ct,ww

      real*8, dimension(nn,nn) :: dd

      real*8, dimension(nn,nn,nn) :: ux,uy,uz
      real*8, dimension(nn,nn,nn) :: div,rotx,roty,rotz
      real*8, dimension(nn,nn,nn) :: strainxx,strainyy,strainzz
      real*8, dimension(nn,nn,nn) :: strainxy,strainyz,strainzx
      real*8, dimension(nn,nn,nn) :: stressxx,stressyy,stresszz
      real*8, dimension(nn,nn,nn) :: stressxy,stressyz,stresszx
      real*8, dimension(nn,nn,nn) :: lambda,mu
      
      real*8, dimension(nn**3) :: derx, dery, derz, ux_el, uy_el, uz_el 
      
!     STRESS CALCULATION
      

          
      do r = 1,nn
         do q = 1,nn
            do p = 1,nn
               !M = p + (q-1)*nn + (r-1)*nn*nn

                do k = 1, nn
                  do j = 1, nn
                     do i = 1, nn
                        L = i + (j-1)*nn + (k-1)*nn*nn
                    
                      ux_el(L) = ux(i,j,k)
                      uy_el(L) = uy(i,j,k)
                      uz_el(L) = uz(i,j,k)

               
                       dxdx = alfa11 + beta12*ct(r) + beta13*ct(q) &
                            + gamma1*ct(q)*ct(r)
                       dydx = alfa21 + beta22*ct(r) + beta23*ct(q) &
                            + gamma2*ct(q)*ct(r)
                       dzdx = alfa31 + beta32*ct(r) + beta33*ct(q) &
                            + gamma3*ct(q)*ct(r)
                       
                       dxdy = alfa12 + beta11*ct(r) + beta13*ct(p) &
                            + gamma1*ct(r)*ct(p)
                       dydy = alfa22 + beta21*ct(r) + beta23*ct(p) &
                            + gamma2*ct(r)*ct(p)
                       dzdy = alfa32 + beta31*ct(r) + beta33*ct(p) &
                            + gamma3*ct(r)*ct(p)
                       
                       dxdz = alfa13 + beta11*ct(q) + beta12*ct(p) &
                            + gamma1*ct(p)*ct(q)
                       dydz = alfa23 + beta21*ct(q) + beta22*ct(p) &
                            + gamma2*ct(p)*ct(q)
                       dzdz = alfa33 + beta31*ct(q) + beta32*ct(p) &
                            + gamma3*ct(p)*ct(q)
                       
                       det_j = dxdz * (dydx*dzdy - dzdx*dydy) &
                             - dydz * (dxdx*dzdy - dzdx*dxdy) &
                             + dzdz * (dxdx*dydy - dydx*dxdy)
                       
               
                      dpdcsi   = dphi(nn-1, ct(p), ct(i)) * phi(nn-1, ct(q), ct(j)) * phi(nn-1, ct(r), ct(k))               
                          dpdeta   = phi(nn-1, ct(p), ct(i)) * dphi(nn-1, ct(q), ct(j)) * phi(nn-1, ct(r), ct(k))
                    dpdzeta = phi(nn-1, ct(p), ct(i)) * phi(nn-1, ct(q), ct(j)) * dphi(nn-1, ct(r), ct(k))


                         derx(L) = (dydx*(dzdy*dpdzeta - dzdz*dpdeta) - dydy*(dzdx*dpdzeta - dzdz*dpdcsi) + &
                             dydz*(dzdx*dpdeta - dzdy*dpdcsi))/det_j
 
                    dery(L) = (dzdx*(dxdy*dpdzeta - dxdz*dpdeta) - dzdy*(dxdx*dpdzeta - dxdz*dpdcsi) + &
                             dzdz*(dxdx*dpdeta - dxdy*dpdcsi))/det_j

                    derz(L) = (dxdx*(dydy*dpdzeta - dydz*dpdeta) - dxdy*(dydx*dpdzeta - dydz*dpdcsi) +  &
                             dxdz*(dydx*dpdeta - dydy*dpdcsi))/det_j
                    

                    enddo
                  enddo  
                enddo      
                          
                          
                          
               strainxx(p,q,r) = dot_product(derx,ux_el)
               strainyy(p,q,r) = dot_product(dery,uy_el)
               strainyy(p,q,r) = dot_product(derz,uz_el)
               strainxy(p,q,r) = 0.5d0*(dot_product(derx,uy_el) + dot_product(dery,ux_el))
               strainyz(p,q,r) = 0.5d0*(dot_product(derz,uy_el) + dot_product(dery,uz_el))
               strainzx(p,q,r) = 0.5d0*(dot_product(derx,uz_el) + dot_product(derz,ux_el))

               rotx(p,q,r) = 0.5d0*(dot_product(dery,uz_el) - dot_product(derz,uy_el))
               roty(p,q,r) = 0.5d0*(dot_product(derz,ux_el) - dot_product(derx,uz_el))
               rotz(p,q,r) = 0.5d0*(dot_product(derx,uy_el) - dot_product(dery,ux_el))

               stressxx(p,q,r) = lambda(p,q,r) * (dot_product(derx,ux_el) + dot_product(dery,uy_el) &
                                + dot_product(derz,uz_el)) + 2.0d0*mu(p,q,r) * dot_product(derx,ux_el)
               stressyy(p,q,r) = lambda(p,q,r) * (dot_product(derx,ux_el) + dot_product(dery,uy_el) &
                                + dot_product(derz,uz_el)) + 2.0d0*mu(p,q,r) * dot_product(dery,uy_el)
               stresszz(p,q,r) = lambda(p,q,r) * (dot_product(derx,ux_el) + dot_product(dery,uy_el) &
                                + dot_product(derz,uz_el)) + 2.0d0*mu(p,q,r) * dot_product(derz,uz_el)
               
               stressyz(p,q,r) = mu(p,q,r) * (dot_product(derz,uy_el) + dot_product(dery,uz_el))
               stresszx(p,q,r) = mu(p,q,r) * (dot_product(derx,uz_el) + dot_product(derz,ux_el))
               stressxy(p,q,r) = mu(p,q,r) * (dot_product(derx,uy_el) + dot_product(dery,ux_el))
                            
                              
            enddo
         enddo
      enddo
               
      
      return
      
      end subroutine MAKE_STRAIN_ROT_STRESS_TENSOR
