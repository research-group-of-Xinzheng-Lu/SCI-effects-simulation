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

!> @brief Computes Legendre polynomial and its derivative on a given point x.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] n  degree of the polynomial
!> @param[in] x point
!> @param[out] p2  Legendre polynomial of degree n evaluated in x
!> @param[out] p2der  derivative of p2 evaluated in x
!> @param[out] p1  Legendre polynomial of degree n-1 evaluated in x
!> @param[out] p1der  derivative of p1 evaluated in x


      subroutine GET_LEGENDRE_VALUE(p2,p2der,p1,p1der,n,x)
                                                                    
      implicit real*8 (a-h,o-z)                            
      
      p2 = 1.d0                                                    
      p2der = 0.d0                                                    
      if (n .eq. 0) return                                          
      p1 = p2
      p2 = x
      p1der = p2der
      p2der = 1.d0
      if (n .eq. 1) return                               
      do k = 1, n-1
        p0 = p1
        p1 = p2
        p0der = p1der
        p1der = p2der
        dk = dfloat(k)                                         
        a1 = (2.d0*dk+1.d0) / (dk+1.d0)
        a2 = -dk / (dk+1.d0)
        p2 = a1*p1*x + a2*p0
        p2der = a1*p1 + a1*p1der*x + a2*p0der
      enddo
      
      return                    
      end subroutine GET_LEGENDRE_VALUE
      
