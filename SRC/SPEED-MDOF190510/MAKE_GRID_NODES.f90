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

!> @brief Makes pointer for connectivity of mesh nodes. 
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nb_node  grid nodes  
!> @param[in] nb_elem  number of elements (hexes)
!> @param[in]  con_matrix   con_matrix(i,1) -> element 
!!                     con_matrix(i,2),...,con_matrix(i,9) -> grid nodes of element con_matrix(i,1) 
!> @param[in] node_wgt  node_wgt(i) contains the multiplicity of the grid node i
!!                  (e.g. if the node node i is shared by 4 elements node_wgt(i) = 4  
!> @param[in] nnz  number of grid nodes + number of grid nodes including repetitions + 1
!!
!> @param[out] node_pointer  pointer for mesh nodes connectivity 
!  
!    E.G. two elements grid   
!    1------2------3
!    |      |      |       nnz = 15
!    | el.1 | el.2 |       node_wgt(1) = node_wgt(3) = node_wgt(4) = node_wgt(6) = 1
!    4------5------6       node_wgt(2) = node_wgt(5) = 2
!      
!    BEFORE
!    node_pointer(0) = 7                      i4count(1) = 7     
!    node_pointer(1) = 8                      i4count(2) = 8 
!    node_pointer(2) = 10                     i4count(3) = 10
!    node_pointer(3) = 11                     i4count(4) = 11
!    node_pointer(4) = 12                     i4count(5) = 12
!    node_pointer(5) = 14                     i4count(6) = 14
!    node_pointer(6) = 15                     i4count(7) = 15
!    node_pointer(7) = ... = node_pointer(15) = 0
!
!
!    AFTER
!    node_pointer(0),...,node_pointer(6) unchanged
!    node_pointer(7) = 1                      i4count(1) = 8 
!    node 1 el. 1    
!    node_pointer(8) = 1                      i4count(2) = 10 
!    node_pointer(9) = 2                      i4count(3) = 11
!    node 2 el. 1 & 2
!    node_pointer(10) = 2                     i4count(4) = 12
!    node 3 el. 2
!    node_pointer(11) = 1                     i4count(5) = 14
!    node 4 el. 1    
!    node_pointer(12) = 1                     i4count(6) = 15
!    node_pointer(13) = 2                     i4count(7) = 16
!    node 5 el. 1 & 2
!    node_pointer(14) = 2
!    node 6 el. 2        
!    node_pointer(15) = 2        
!    node 7 el. 2
!
!**************************************************************************************************    


      subroutine MAKE_GRID_NODES(nb_node,nb_elem,con_matrix,node_wgt,nnz,node_pointer)
      
      
      implicit none
      
      integer*4 :: nb_node,nb_elem,nnz
      integer*4 :: i,j,ie
      
      integer*4, dimension(:), allocatable :: i4count
      integer*4, dimension(nb_node) :: node_wgt
      integer*4, dimension(0:nnz) :: node_pointer
      
      integer*4, dimension(nb_elem,9) :: con_matrix
      

      allocate(i4count(nb_node))
      
      node_pointer(0) = nb_node +1
      do j = 1,nb_node
         node_pointer(j) = node_pointer(j-1) + node_wgt(j)
      enddo
      
      do j = 1,nb_node
         i4count(j) = node_pointer(j-1)
      enddo
      
      
      do ie = 1,nb_elem
         do i = 1,8
            j = con_matrix(ie,i +1)
            node_pointer(i4count(j)) = ie      
            i4count(j) = i4count(j) +1
         enddo
      enddo
      
      deallocate(i4count)
      
      return
      
      end subroutine MAKE_GRID_NODES

