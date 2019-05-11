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

!> @brief Makes Gauss-Legendre-Lobatto nodes.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] np number of points
!> @param[out] ct  LGL nodes

      subroutine MAKE_LGL_NODES(np,ct)
      
      implicit real*8 (a-h,o-z)

      integer*4 :: np
      
      real*8, dimension(np) :: ct
          
      parameter (nstep=1000,acc=1.d-15)
      
      n=np-1
      ct(1)=-1.d0
      ct(np)=1.d0
      n2=idint(0.5d0*np)
      dx=2.d0/dfloat(nstep)
      x=-1.d0
      iroot=1
      call GET_LEGENDRE_VALUE (p2,a1,p1,p1der,n,x)
      do while (iroot.lt.n2)
        x=x+dx
        call GET_LEGENDRE_VALUE (p2,a2,p1,p1der,n,x)
        if (dabs(a2).le.acc) then
          iroot=iroot+1
          ct(iroot)=a2
        else
          aa=a1*a2
          if (aa.lt.0.d0) then
            iroot=iroot+1
            ct(iroot)=FIND_ROOT_POL(x-dx,x,acc,n)
          endif
        endif
        a1=a2
      enddo
      
      n_filt=2*n2
      if (n_filt.ne.np) then         ! np odd
        ct(n2+1)=0.d0
        do i=1,n2-1
          ct(np-i)=-ct(i+1)
        enddo
      else                           ! np even
        do i=1,n2-1
          ct(np-i)=-ct(i+1)
        enddo
      endif
      

      
      return
      
      end subroutine MAKE_LGL_NODES
