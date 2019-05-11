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

!> @brief Computes derivative of Legendre basis function.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] ng polynomial degree
!> @param[in] csii GLL node
!> @param[in] csij GLL node identifying the basis function
!> @param[out] DPHI derivative of the j-th (csij) basis function of degree ng at node csii 
!!              DPHI(ng,csii,csij) = DPHI_csij/dcsi (csii) 

      function DPHI(ng,csii,csij)
                                                                    
     implicit real*8 (a-h,o-z)                            

     integer*4 :: ng 
     
     call GET_LEGENDRE_VALUE_2DER(p2der2,p2,p2der,p1,p1der,ng,csii)     
     call GET_LEGENDRE_VALUE_2DER(q2der2,q2,q2der,q1,q1der,ng,csij)   

     if(abs(csij + 1.d0).le. 1.e-8 .and. abs(csii + 1.d0).le. 1.e-8) then
         DPHI = -(ng+1.d0)*ng/4.d0 
     elseif(abs(csij - 1.d0).le. 1.e-8 .and. abs(csii - 1.d0).le. 1.e-8) then
         DPHI = (ng+1.d0)*ng/4.d0
     elseif(abs(csii - csij).le. 1.e-8) then
         DPHI = 0.d0   
     else
         anum = (csij-csii)*ng*(ng+1.d0)*q2 + (1.d0-csij**2.d0)*q2der
         aden = p2*ng*(ng+1.d0)*(csii - csij)**2.d0
         DPHI =  anum/aden
     endif
           
     if(abs(DPHI) .lt. 1.e-10) then
         DPHI = 0.d0
     endif

     return    
     end function DPHI

