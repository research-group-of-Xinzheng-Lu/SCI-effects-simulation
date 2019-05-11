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

!> @brief Computes the nearest node with respet to (xt,yt,zt) 
!! starting from a given elevation.  
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] nloc number of local nodes
!> @param[in] xs vertex x-coordinate of local nodes
!> @param[in] ys vertex y-coordinate of local nodes
!> @param[in] zs vertex z-coordinate of local nodes
!> @param[in] xt x-coor current vertex 
!> @param[in] yt y-coor current vertex 
!> @param[in] zt z-coor current vertex 
!> @param[in] elevation z-value from wich we start the search
!> @param[out] nt index of the  nearest node from (xt,yt,zt) 
!> @param[out] dist_min minimal distance from (xt,yt,zt)
 
      subroutine GET_NEAREST_NODE_PGM(nloc,xs,ys,zs,xt,yt,zt,nt,dist_min,elevation)
      
      implicit none
      
      integer*4 :: i, nloc, nt

      real*8 :: xt,yt,zt
      real*8 :: elevation
      real*8 :: dist,dist_min

      real*8, dimension(nloc) :: xs, ys, zs
      
      dist_min = 1.0d30
      nt = 1

      do i = 1,nloc

         if (zs(i) .gt. elevation) then
            
            dist = (xs(i) - xt)*(xs(i) - xt) + (ys(i) - yt)*(ys(i) - yt) + (zs(i) - zt)*(zs(i) - zt) 

            if (dist .lt. dist_min) then
               dist_min = dist; nt = i
            endif

        endif

      enddo
      
      return
      end subroutine GET_NEAREST_NODE_PGM

