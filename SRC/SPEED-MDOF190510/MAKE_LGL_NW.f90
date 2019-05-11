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

!> @brief Makes Gauss-Legendre-Lobatto nodes, weigths and spectral derivatives.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nb_pnt  polynomial  degree
!> @param[out] xq  LGL nodes
!> @param[out] wq  LGL weights
!> @param[out] dd  matrix for  spectral derivatives

     subroutine MAKE_LGL_NW(nb_pnt, xq, wq, dd)
      
     implicit real*8 (a-h,o-z)
      
      parameter (nstep = 1000, acc = 1.d-15)
      
      
      integer*4, intent(inout) :: nb_pnt
      real*8, dimension(nb_pnt) :: xq, wq      
      real*8, dimension(nb_pnt) :: dd(nb_pnt,nb_pnt)
              
      n = nb_pnt-1
      xq(1) = -1.d0
      xq(nb_pnt) = 1.d0
      n2=idint(0.5d0*nb_pnt)
      dx=2.d0/dfloat(nstep)

      x = -1.d0
      iroot = 1

      call GET_LEGENDRE_VALUE(p2,a1,p1,p1der,n,x)
      
      do while (iroot .lt. n2)
        x = x + dx
        call GET_LEGENDRE_VALUE(p2,a2,p1,p1der,n,x)
        
        if (dabs(a2) .le. acc) then
          iroot = iroot + 1
          xq(iroot) = a2
        else
          aa = a1 * a2
          if (aa .lt. 0.d0) then
            iroot = iroot + 1
            xq(iroot) = FIND_ROOT_POL(x-dx,x,acc,n)
          endif
        endif
        
        a1 = a2
      enddo
      
      n_filt = 2*n2
      if (n_filt .ne. nb_pnt) then         ! nb_pnt odd
        xq(n2+1) = 0.d0
        do i = 1, n2-1
          xq(nb_pnt-i) = -xq(i+1)
        enddo
      else                           ! nb_pnt even
        do i = 1, n2-1
          xq(nb_pnt-i) = -xq(i+1)
        enddo
      endif
      
      
      xn = dfloat(n)
      acost = 2.d0/(xn*(xn+1.d0))
      
      do j = 1,nb_pnt
         call GET_LEGENDRE_VALUE(p2,p2der,p1,p1der,n,xq(j))
         den = p2*p2
         wq(j) = acost/den
         
         do i = 1, nb_pnt
            if (i .ne. j) then
               call GET_LEGENDRE_VALUE(pnum,p2der,p1,p1der,n,xq(i))
               den = p2 * (xq(i)-xq(j))
               dd(i,j) = pnum / den
            endif
         enddo
         
      enddo
      
      do j = 2, nb_pnt-1
         dd(j,j) = 0.d0
      enddo
      
      dd(1,1) = -0.25d0 * xn * (xn+1.d0)
      dd(nb_pnt,nb_pnt) = 0.25d0 * xn * (xn+1.d0)
     
      
      return
      
      end subroutine MAKE_LGL_NW
