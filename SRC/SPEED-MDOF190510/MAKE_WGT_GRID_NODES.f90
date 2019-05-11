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

!> @brief Computes multeplicity for mesh nodes. 
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nb_node  grid nodes
!> @param[in] nb_elem  number of elements (hexes)
!> @param[in] con_matrix con_matrix(i,1) -> element,  
!!                     con_matrix(i,2),...,con_matrix(i,9) -> grid nodes of element con_matrix(i,1) 
!> @param[out] node_wgt  node_wgt(i) matrix containing the multiplicity of the i-th grid node 
!!                  (e.g. if the node i is shared by 4 elements node_wgt(i) = 4  
!> @param[out] nb_nz number of grid nodes + number of grid nodes including repetitions + 1

 
      subroutine MAKE_WGT_GRID_NODES(nb_node,nb_elem,con_matrix,node_wgt,nb_nz)
      
      implicit none
      
      integer*4 :: nb_node,nb_elem,nb_nz
      integer*4 :: i,j,ie
      
      integer*4, dimension(nb_node) :: node_wgt
      integer*4, dimension(nb_elem,9) :: con_matrix

      
      nb_nz = nb_node + 1
      node_wgt = 0
      
      do ie = 1,nb_elem
         do i = 1,8
            j = con_matrix(ie,i+1)         
            node_wgt(j) = node_wgt(j) + 1
            nb_nz = nb_nz + 1
         enddo
      enddo
      
      return
      
      end subroutine MAKE_WGT_GRID_NODES

