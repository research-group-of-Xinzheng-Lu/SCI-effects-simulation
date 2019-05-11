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


!> @brief Computes the nearest node with respet to (xt,yt,zt). 
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] n number of local nodes
!> @param[in] xs vertex x-coordinate of local nodes
!> @param[in] ys vertex y-coordinate of local nodes
!> @param[in] zs vertex z-coordinate of local nodes
!> @param[in] xt x-coor current vertex 
!> @param[in] yt y-coor current vertex 
!> @param[in] zt z-coor current vertex 
!> @param[out] nt index of the  nearest node from (xt,yt,zt) 
!> @param[out] dist_min minimal distance from (xt,yt,zt)

       subroutine GET_NEAREST_NODE(n,xs,ys,zs,xt,yt,zt,nt,dist_min)

      implicit none
      
      integer*4 :: n,nt, i
      integer*4, dimension(1) :: pos

      !real*8 :: dx,dy,dz,d2,d2min
      real*8 :: xt,yt,zt
      real*8 :: dist_min
      
      real*8, dimension(n) :: xs,ys,zs
      real*8, dimension(n) :: dist_all
      
      !d2min = 1.0d30
      
      !do i = 1,n
      !   dx = xs(i) - xt
      !   dy = ys(i) - yt
      !   dz = zs(i) - zt
      !   d2 = dx*dx + dy*dy + dz*dz
      !   if (d2.lt.d2min) then
      !      d2min = d2
      !      nt = i
      !      
      !   endif
      !enddo

      dist_all = (xs-xt)*(xs-xt) + (ys-yt)*(ys-yt) + (zs-zt)*(zs-zt)
      dist_min = minval(dist_all)          
      pos = minloc(dist_all); nt = pos(1)
      
           
      return
      end subroutine GET_NEAREST_NODE

