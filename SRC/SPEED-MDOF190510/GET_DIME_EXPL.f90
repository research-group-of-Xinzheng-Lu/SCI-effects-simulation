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

!> @brief Computes local number of nodes where explosive source is imposed.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] Xipo x-coordinate of Hypocenter
!> @param[in] Yipo y-coordinate of Hypocenter
!> @param[in] Zipo z-coordinate of Hypocenter
!> @param[in] X1 vertex coordinate (x) of the triangular fault
!> @param[in] Y1 vertex coordinate (y) of the triangular fault
!> @param[in] Z1 vertex coordinate (z) of the triangular fault
!> @param[in] X2 vertex coordinate (x) of the triangular fault
!> @param[in] Y2 vertex coordinate (y) of the triangular fault
!> @param[in] Z2 vertex coordinate (z) of the triangular fault
!> @param[in] X3 vertex coordinate (x) of the triangular fault
!> @param[in] Y3 vertex coordinate (y) of the triangular fault
!> @param[in] Z3 vertex coordinate (z) of the triangular fault
!> @param[in] nnod number of local nodes
!> @param[in] xs x-coordinate of spectral nodes
!> @param[in] ys y-coordinate of spectral nodes
!> @param[in] zs z-coordinate of spectral nodes
!> @param[in] nnloc (obsolete)
!> @param[out] node_expl  local number of explosive nodes

      subroutine GET_DIME_EXPL(Xipo,Yipo,Zipo, &
                                 X1,Y1,Z1, &
                                 X2,Y2,Z2, &
                                 X3,Y3,Z3, &
                                 nnod,xs,ys,zs, &
                                 node_expl,nnloc)
         
      implicit none
      
      integer*4 :: nnod,node_expl,nnloc
      integer*4 :: isn

      real*8 :: Xipo,Yipo,Zipo,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,ux,uy,uz,vx,vy,vz,tol
      real*8 :: Xmax,Xmin,Ymax,Ymin,Zmin,Zmax,a,b,c, den
      real*8, dimension(nnloc) :: xs,ys,zs

      node_expl = 0
      tol = 5.0d0
      
      ux=(X1-X2)/sqrt((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
      uy=(Y1-Y2)/sqrt((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
      uz=(Z1-Z2)/sqrt((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
      vx=(X3-X2)/sqrt((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2)
      vy=(Y3-Y2)/sqrt((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2)
      vz=(Z3-Z2)/sqrt((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2)
      
      a = uy * vz - uz * vy 
      b = uz * vx - ux * vz
      c = ux * vy - uy * vx
      den = sqrt(a**2 + b**2 + c**2)

      !Check on the relative position of fault end points
      !In this preliminary part we check if P1 point is the lower-left fault point or not
      !If this check would fail the two points are swapped
  
      Xmax=max(X1,X2,X3)+tol
      Xmin=min(X1,X2,X3)-tol
      Ymax=max(Y1,Y2,Y3)+tol
      Ymin=min(Y1,Y2,Y3)-tol
      Zmax=max(Z1,Z2,Z3)+tol
      Zmin=min(Z1,Z2,Z3)-tol
  
  
      do isn = 1,nnod
         if ( dabs( a*(xs(isn)-Xipo) + b*(ys(isn)-Yipo) + c*(zs(isn)-Zipo)/den ).le.tol) then
            if ((xs(isn).ge.Xmin).and.(xs(isn).le.Xmax)) then
                if ((ys(isn).ge.Ymin).and.(ys(isn).le.Ymax)) then
                    if ((zs(isn).ge.Zmin).and.(zs(isn).le.Zmax)) then
                        node_expl = node_expl+1
                    endif
                endif
            endif
         endif         
      enddo

      return
      end subroutine GET_DIME_EXPL

