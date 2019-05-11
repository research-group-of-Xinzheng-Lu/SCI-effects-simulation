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

!> @brief Computes total number of boundary nodes 
!! and second derivative on a given point x.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nb_nodes number of nodes
!> @param[in] nb_nz length of con_spc
!> @param[in] con_spc  spectral connectivity vector
!> @param[in] nb_load  load number for boundary condition
!> @param[in] lab_bc  label for boundary condition 
!> @param[in] local_n_num = local node numbering (dummy)
!> @param[in] node_indicator  node_indicator(i) = 1 if the node belongs to the boundary
!> @param[out] nb_nodes_face_total number of nodes lying on the boundary

      subroutine GET_NODE_FROM_FACE(nb_nodes, nb_nz, con_spc, nb_load, lab_bc,&
                                    nb_nodes_face_total, node_indicator, local_n_num)
      
      implicit none
      
      integer*4 :: nb_nodes,nb_nz,nb_load
      integer*4 :: nb_nodes_face_total
      integer*4 :: i,j,ie, nb_el, nb_nodes_face, cek_el
      
      integer*4, dimension(0:nb_nz) :: con_spc
      integer*4, dimension(nb_load) :: lab_bc
      integer*4, dimension(nb_nodes) :: node_indicator, local_n_num

      
      node_indicator = 0
      nb_nodes_face_total = 0
      nb_el = con_spc(0) - 1
      
      if (nb_nz .gt. 0) then
         
         do ie = 1,nb_el
            nb_nodes_face = con_spc(ie) - con_spc(ie -1) -1
            
            cek_el = 0
            do j = 1,nb_load
               if (con_spc(con_spc(ie -1)) .eq. lab_bc(j)) cek_el = 1
            enddo
            
            if (cek_el .ne. 0) then          
               do i = 1, nb_nodes_face                      
                  node_indicator(con_spc(con_spc(ie -1) +i)) = 1
               enddo
            endif
         enddo
         
         do i = 1,nb_nodes
            if (node_indicator(i) .ne. 0) nb_nodes_face_total = nb_nodes_face_total + 1
         enddo
      endif
      
      return
      
      end subroutine GET_NODE_FROM_FACE

