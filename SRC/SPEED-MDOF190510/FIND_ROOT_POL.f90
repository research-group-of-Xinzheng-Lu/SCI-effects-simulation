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

!> @brief Find a root of a polynomial of degree n.
!! @author Ilario Mazzieri
!> @date September, 2013 - Creation
!> @version 1.0
!> @param[in] x1 left bound of the interval [x1,x2]
!> @param[in] x2 right bound of the interval [x1,x2]
!> @param[in] xacc  acceleration parameter
!> @param[in] n  polynomial degree
!> @param[out] FIND_ROOT_POL  zero of the polynomial function of order n 


      real*8 function FIND_ROOT_POL(x1,x2,xacc,n)

      use speed_exit_codes
      
      implicit real*8(a-h,o-z)
      
      parameter (jmax=1000)
      
      call GET_LEGENDRE_VALUE(p2,fmid,p1,p1der,n,x2)
      call GET_LEGENDRE_VALUE(p2,f,p1,p1der,n,x1)

      if (f*fmid .ge. 0.d0) then
        write(*,*) 'root non-bracketed !'
        call EXIT(EXIT_ROOT)
      endif
      
      if (f .lt. 0.d0) then
        FIND_ROOT_POL = x1
        dx = x2 - x1
      else
        FIND_ROOT_POL = x2
        dx = x1 - x2
      endif
      
      do j = 1,jmax
        dx = 0.5d0*dx
        xmid = FIND_ROOT_POL + dx
        call GET_LEGENDRE_VALUE(p2,fmid,p1,p1der,n,xmid)
        if (fmid .le. 0.d0) FIND_ROOT_POL = xmid
        
        if (dabs(dx) .lt. xacc .and. dabs(fmid) .le. xacc) then
          FIND_ROOT_POL = xmid
          return
        endif
        
      enddo
      
      end function FIND_ROOT_POL
      

