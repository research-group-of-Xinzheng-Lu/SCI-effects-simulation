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

!> @brief Computes the stress tensor.
!! @author Ilario Mazzieri
!> @date November, 2014
!> @version 1.0
!> @param[in] nn number of 1-D Legendre nodes
!> @param[in] lambda nodal values of Lame coefficient lambda 
!> @param[in] mu nodal values of Lame coefficient mu
!> @param[in] duxdx nodal values for spatial derivatives of the displacement 
!> @param[in] duydx nodal values for spatial derivatives of the displacement
!> @param[in] duzdx nodal values for spatial derivatives of the displacement
!> @param[in] duxdy nodal values for spatial derivatives of the displacement
!> @param[in] duydy nodal values for spatial derivatives of the displacement
!> @param[in] duzdy nodal values for spatial derivatives of the displacement
!> @param[in] duxdz nodal values for spatial derivatives of the displacement
!> @param[in] duydz nodal values for spatial derivatives of the displacement      
!> @param[in] duzdz nodal values for spatial derivatives of the displacement
!> @param[in] strain_visc_xx nodal values for viscoelastic stress tensor
!> @param[in] strain_visc_xy nodal values for viscoelastic stress tensor
!> @param[in] strain_visc_xz nodal values for viscoelastic stress tensor
!> @param[in] strain_visc_yy nodal values for viscoelastic stress tensor
!> @param[in] strain_visc_yz nodal values for viscoelastic stress tensor
!> @param[in] strain_visc_zz nodal values for viscoelastic stress tensor
!> @param[in] Ylambda anelastic coefficients (1st Lame modulus)
!> @param[in] Ymu anelastic coefficients  (2nd Lame modulus)
!> @param[in] ie  local element id
!> @param[out] sxx nodal values for the stress tensor
!> @param[out] syy nodal values for the stress tensor
!> @param[out] szz nodal values for the stress tensor
!> @param[out] syz nodal values for the stress tensor
!> @param[out] szx nodal values for the stress tensor
!> @param[out] sxy nodal values for the stress tensor


      subroutine MAKE_STRESS_TENSOR_DAMPED(nn,lambda,mu,&
                               duxdx,duydx,duzdx,&
                               duxdy,duydy,duzdy,&
                               duxdz,duydz,duzdz,&
                               strain_visc_xx, strain_visc_xy, strain_visc_xz, &
                               strain_visc_yy, strain_visc_yz, strain_visc_zz, &
                               Ylambda,Ymu, N_SLS, & 
                               sxx,syy,szz,&
                               syz,szx,sxy,ie)
      
      
      implicit none
      
      integer*4 :: nn,N_SLS
      integer*4 :: p,q,r,ie
      
      real*8, dimension(N_SLS) :: Ylambda, Ymu, temp
      real*8, dimension(nn,nn,nn) :: lambda,mu
      real*8, dimension(nn,nn,nn) :: sxx,syy,szz,syz,szx,sxy
      real*8, dimension(nn,nn,nn) :: duxdx,duxdy,duxdz,duydx,duydy,duydz,duzdx,duzdy,duzdz
      real*8, dimension(nn,nn,nn,N_SLS) :: strain_visc_xx, strain_visc_xy, strain_visc_xz, &
                                           strain_visc_yy, strain_visc_yz, strain_visc_zz

      
      do r = 1,nn
         do q = 1,nn
            do p = 1,nn
            
               
               temp = strain_visc_xx(p,q,r,:) + strain_visc_yy(p,q,r,:) + strain_visc_zz(p,q,r,:) 
                   
               sxx(p,q,r) = lambda(p,q,r) * (duxdx(p,q,r) +duydy(p,q,r) +duzdz(p,q,r) &
                                        - DOT_PRODUCT(Ylambda,temp)) &
                                        +2.0d0*mu(p,q,r) * (duxdx(p,q,r) - DOT_PRODUCT(Ymu,strain_visc_xx(p,q,r,:)))
                                        
               syy(p,q,r) = lambda(p,q,r) * (duxdx(p,q,r) +duydy(p,q,r) +duzdz(p,q,r) &
                                        - DOT_PRODUCT(Ylambda,temp)) &
                                        +2.0d0*mu(p,q,r) * (duydy(p,q,r) - DOT_PRODUCT(Ymu,strain_visc_yy(p,q,r,:)))
                                        
               szz(p,q,r) = lambda(p,q,r) * (duxdx(p,q,r) +duydy(p,q,r) +duzdz(p,q,r) &
                                        - DOT_PRODUCT(Ylambda,temp)) &
                                        +2.0d0*mu(p,q,r) * (duzdz(p,q,r) - DOT_PRODUCT(Ymu,strain_visc_zz(p,q,r,:)))
               
               syz(p,q,r) = mu(p,q,r) * (duydz(p,q,r) + duzdy(p,q,r) -  DOT_PRODUCT(Ymu,strain_visc_yz(p,q,r,:)))
               szx(p,q,r) = mu(p,q,r) * (duzdx(p,q,r) + duxdz(p,q,r) -  DOT_PRODUCT(Ymu,strain_visc_xz(p,q,r,:)))
               sxy(p,q,r) = mu(p,q,r) * (duxdy(p,q,r) + duydx(p,q,r) -  DOT_PRODUCT(Ymu,strain_visc_xy(p,q,r,:)))
               
              ! write(*,*) p,q,r
              ! write(*,*) temp
              ! write(*,*) DOT_PRODUCT(Ylambda,temp)
              ! write(*,*) DOT_PRODUCT(Ymu,strain_visc_xx(p,q,r,:))
              ! write(*,*) DOT_PRODUCT(Ymu,strain_visc_yy(p,q,r,:))
              ! write(*,*) DOT_PRODUCT(Ymu,strain_visc_zz(p,q,r,:))
              ! write(*,*) DOT_PRODUCT(Ymu,strain_visc_yz(p,q,r,:))
              ! write(*,*) DOT_PRODUCT(Ymu,strain_visc_xz(p,q,r,:))
              ! write(*,*) DOT_PRODUCT(Ymu,strain_visc_xy(p,q,r,:))
              ! write(*,*) '--------------------------------------'
            enddo
         enddo
      enddo
      
      return
      
      end subroutine MAKE_STRESS_TENSOR_DAMPED

