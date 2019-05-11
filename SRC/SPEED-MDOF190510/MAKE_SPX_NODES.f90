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

!> @brief Makes a pointer for spectral connectivity vector.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nb_node  GLL nodes
!> @param[in] spx_con spectral connectivity vector
!> @param[in] node_wgt node_wgt(i) contains the multiplicity of the GLL node i
!!                  (e.g. if the GLL node i is shared by 4 elements node_wgt(i) = 4  
!> @param[in] nnz  number of GLL nodes + number of GLL nodes including repetitions + 1
!> @param[out] node_pointer as explained in MAKE_GRID_NODES.f90

      subroutine MAKE_SPX_NODES(nb_node,spx_con_nnz,spx_con,node_wgt,nnz,node_pointer)
      
  
      
      implicit none
      
      integer*4 :: nb_node,spx_con_nnz,nnz
      integer*4 :: i,j,ie,nb_elem,nn3

      integer*4, dimension(:), allocatable :: i4count
      integer*4, dimension(0:spx_con_nnz) :: spx_con
      integer*4, dimension(nb_node) :: node_wgt
      integer*4, dimension(0:nnz) :: node_pointer
      
      allocate(i4count(nb_node))
      
      node_pointer(0) = nb_node + 1
      do j = 1, nb_node
         node_pointer(j) = node_pointer(j -1) + node_wgt(j)
      enddo
      
      do j = 1, nb_node
         i4count(j) = node_pointer(j -1)
      enddo
      
      nb_elem = spx_con(0) - 1
      do ie = 1, nb_elem
      
         nn3 = spx_con(ie) - spx_con(ie -1) -1
         do i = 1, nn3
      
            j = spx_con(spx_con(ie -1) +i)       
            node_pointer(i4count(j)) = ie
            i4count(j) = i4count(j) + 1
         enddo
      enddo
      
      deallocate(i4count)
      
      return
      
      end subroutine MAKE_SPX_NODES

