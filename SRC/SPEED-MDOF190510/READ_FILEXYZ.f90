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

!> @brief Reads mesh of triangular elements (topography or alluvial basin)
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] filec file name for the mesh
!> @param[in] num_nodes  number of nodes
!> @param[in] num_tria  number of triangular fault 
!> @param[out] x x-coordinate for the node
!> @param[out] y y-coordinate for the node
!> @param[out] z z-coordinate for the node
!> @param[out] node1 1-vertex of the triangular fault
!> @param[out] node2 2-vertex of the triangular fault
!> @param[out] node3 3-vertex of the triangular fault
!> @param[out] max_length  max diameter of the triangles

      subroutine READ_FILEXYZ(filec,num_nodes,num_tria, &
                                x,y,z,&
                                node1,node2,node3,&
                                max_length)

      implicit none

      integer*4 :: num_nodes,num_tria
      integer*4 :: i
      integer*4 :: trash
      real*8, dimension(num_nodes) :: x,y,z
      integer*4, dimension(num_tria) :: node1,node2,node3
      character*70 :: filec
      integer*4 :: ileft,iright
      integer*4 :: status
      character*100000 :: input_line
      real*8 :: length_edge1,length_edge2,length_edge3
      real*8 :: max_edge_length,max_length

      max_length = 0.0d0

      open(20,file=filec)
      read(20,'(A)',IOSTAT = status) input_line
      ileft = 1
      iright = len(input_line)
      read(input_line(ileft:iright),*) num_nodes,num_tria 

      do i = 1,num_nodes
        read(20,'(A)',IOSTAT = status) input_line
        if (status.ne.0) exit
                  ileft = 1
                iright = len(input_line)
                read(input_line(ileft:iright),*)trash,x(i),y(i),z(i) 
      enddo

      do i = 1,num_tria
        read(20,'(A)',IOSTAT = status) input_line
        if (status.ne.0) exit
                ileft = 1
                iright = len(input_line)
                read(input_line(ileft:iright),*)trash,node1(i),node2(i),node3(i) 
      enddo 
                    
      close(20)

     ! Compute max tria edge length

      do i = 1,num_tria

        length_edge1 = ( x(node1(i)) - x(node2(i)) )**2 + ( y(node1(i)) - y(node2(i)) )**2
        length_edge2 = ( x(node2(i)) - x(node3(i)) )**2 + ( y(node2(i)) - y(node3(i)) )**2
        length_edge3 = ( x(node3(i)) - x(node1(i)) )**2 + ( y(node3(i)) - y(node1(i)) )**2

        max_edge_length = max(length_edge1,length_edge2,length_edge3)

        if (max_length.lt.max_edge_length) then
                 max_length = max_edge_length
        endif
                
      enddo
       
        max_length = sqrt(max_length)

      return

      end subroutine READ_FILEXYZ

