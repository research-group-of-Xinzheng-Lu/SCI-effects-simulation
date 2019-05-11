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

!> @brief Evaluates Lengendre basis function on a specific point. 
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] ng degree of the polynomial
!> @param[in] csii i-th Legendre function 
!> @param[in] csij node where the Legendre function is evaluated
!> @param[out] PHI  PHI^{ng}_i(xi_j) i-th (1-D) basis function of degree ng evaluate at 
!!                          node csij 


      real*8 function PHI(ng,csii, csij)
      
      implicit real*8(a-h,o-z)
      integer*4 :: ng

      call GET_LEGENDRE_VALUE_2DER(p2der2,p2,p2der,p1,p1der,ng,csii)          
      call GET_LEGENDRE_VALUE_2DER(q2der2,q2,q2der,q1,q1der,ng,csij)   

      if(abs(csii - csij).le. 1.e-8) then
          PHI = 1.d0  
      else
          PHI =  (1.d0-csij**2.d0)*q2der /(ng*(ng+1.d0)*p2*(csii-csij))
          if(abs(PHI) .lt. 1.e-10) then
             PHI = 0.d0
          endif
      endif

      return
      end function PHI
