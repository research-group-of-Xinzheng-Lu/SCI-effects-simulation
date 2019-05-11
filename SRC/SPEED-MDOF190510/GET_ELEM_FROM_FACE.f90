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

!> @brief Find element index from verteces number.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nb_nz_el number of element
!> @param[in] list_el vector containing spectral node repetition
!> @param[in] v1 node index of first vertex of the element (quad) 
!> @param[in] v2 node index of second vertex of the element (quad)
!> @param[in] v3 node index of third vertex of the element (quad)
!> @param[out] ie element index

      subroutine GET_ELEM_FROM_FACE(nb_nz_el,list_el,v1,v2,v3,ie)
      
      
      implicit none
      
      integer*4 :: nb_nz_el,v1,v2,v3,ie
      integer*4 :: i,j,k

      integer*4, dimension(0:nb_nz_el) :: list_el
      
      ie = 0
      
      do i = list_el(v1 -1),list_el(v1) -1
         do j = list_el(v2 -1),list_el(v2) -1
            do k = list_el(v3 -1),list_el(v3) -1
            
               if ((list_el(i).eq.list_el(j)).and.(list_el(j).eq.list_el(k))) ie = list_el(i)
               
            enddo
         enddo
      enddo
      
      return
      
      end subroutine GET_ELEM_FROM_FACE

