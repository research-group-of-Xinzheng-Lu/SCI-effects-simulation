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

!> @brief Assignes labels for multi not-honoring technique.
!> @warning Not yet validated. Do not use it.
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] nn_s number of nodes 
!> @param[in] zz_elevation  depth of the nodes with respect to the surface
!> @param[in] j = tag
!> @param[in] xx_s  x-coordinate of spectral nodes
!> @param[out] sub_tag_all labels for multi not-honoring

      subroutine MAKE_SUBTAG_ALLUVIAL(nn_s, zz_elevation, j, sub_tag_all, xx_s)
            
      implicit none
      
      integer*4 :: nn_s,isn,j

      integer*4, dimension(nn_s) :: sub_tag_all
      
      real*8, dimension(nn_s) :: zz_elevation
      real*8, dimension(nn_s) :: xx_s
      

      do isn = 1, nn_s
         if ((zz_elevation(isn) .ge. 0.0d0) .and. (sub_tag_all(isn) .eq. 4)) then

                !if ((xx_s(isn) .ge. 685956).and.(xx_s(isn) .le. 701105)) then
                        sub_tag_all(isn) = j
                !else
                !        sub_tag_all(isn) = j+4
                !endif
        endif
     enddo
      
      return
      
      end subroutine MAKE_SUBTAG_ALLUVIAL

