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

!> @brief Computes multeplicity for spectral nodes. 
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nb_node  GLL nodes
!> @param[in] spx_con spectral connectivity vector
!> @param[out] node_wgt node_wgt(i) contains the multiplicity of the i-th GLL node
!!                  (e.g. if the node i is shared by 4 elements node_wgt(i) = 4  
!> @param[out] nb_nz  number of GLL nodes + number of GLL nodes including repetitions + 1

      subroutine MAKE_WGT_SPX_NODES(nb_node,spx_con_nnz,spx_con,node_wgt,nb_nz)
      
      
      implicit none
      
      integer*4 :: nb_node,spx_con_nnz,nb_nz
      integer*4 :: i,j,ie,nb_elem,nn3
      
      integer*4, dimension(0:spx_con_nnz) :: spx_con
      integer*4, dimension(nb_node) :: node_wgt
      
      nb_nz = nb_node +1
      node_wgt = 0
      
      nb_elem = spx_con(0) -1
      do ie = 1, nb_elem
         nn3 = spx_con(ie) - spx_con(ie -1) -1
         do i = 1,nn3
            j = spx_con(spx_con(ie -1) +i)
            
            node_wgt(j) = node_wgt(j) +1
            nb_nz = nb_nz +1
         enddo
      enddo
      
      return
      
      end subroutine MAKE_WGT_SPX_NODES

