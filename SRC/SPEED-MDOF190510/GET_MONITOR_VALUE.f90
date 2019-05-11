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

!> @brief Computes the  mean value. 
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] nb_nod number of 1D GLL nodes
!> @param[in] xq coor GLL node
!> @param[in] val value in each GLL node
!> @param[in] xref x-coordinate of the monitor
!> @param[in] yref y-coordinate of the monitor
!> @param[in] zref z-coordinate of the monitor
!> @param[out] res mean value (interpolation)


      subroutine GET_MONITOR_VALUE(nb_nod, xq, val, xref, yref, zref, res)
            
      implicit none
      
      integer*4 :: nb_nod
      integer*4 :: i,j,k,h
      
      real*8 :: xref,yref,zref,res
      real*8 :: fx,fy,fz

      real*8, dimension(nb_nod) :: xq

      real*8, dimension(nb_nod,nb_nod,nb_nod) :: val

      res = 0.0d0
      
      do k = 1,nb_nod
         fz = 1.0d0
         do h = 1,nb_nod
            if (h.ne.k) fz = fz * (zref - xq(h)) / (xq(k) - xq(h))
         enddo
         
         do j = 1,nb_nod
            fy = 1.0d0
            do h = 1,nb_nod
               if (h.ne.j) fy = fy * (yref - xq(h)) / (xq(j) - xq(h))
            enddo
            
            do i = 1,nb_nod
               fx = 1.0d0
               do h = 1,nb_nod
                  if (h.ne.i) fx = fx * (xref - xq(h)) / (xq(i) - xq(h))
               enddo
               
               res = res + val(i,j,k) * fx * fy * fz
            enddo
         enddo
      enddo
      
      end subroutine GET_MONITOR_VALUE

