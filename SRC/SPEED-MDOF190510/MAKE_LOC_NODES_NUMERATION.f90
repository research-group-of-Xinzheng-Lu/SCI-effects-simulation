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

!> @brief Makes local numeration vector.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nnloc number of local nodes
!> @param[in,out] local_n_num  local node numbering
!> @param[in] nrcv number of nodes to receive
!> @param[in] node_rcv list of nodes to receive


      subroutine  MAKE_LOC_NODES_NUMERATION(nnloc, local_n_num, nrcv, node_rcv) 

      implicit none

      integer*4 :: nnloc, nrcv
      integer*4 :: i,j

      integer*4, dimension(nnloc) :: copy_n_num
      integer*4, dimension(nrcv) :: node_rcv
      integer*4, dimension(nnloc), intent(inout) :: local_n_num
      
      do i = 1, nrcv
        do j = 1, nnloc
            if(node_rcv(i) .eq. local_n_num(j)) local_n_num(j) = 0
        enddo
      enddo
      
      i = 1
      do j = 1, nnloc
         if(local_n_num(j) .ne. 0) then 
                  copy_n_num(i) = local_n_num(j) 
                  i = i + 1
         endif
      enddo
      
      copy_n_num(nnloc - nrcv + 1 : nnloc) = node_rcv
      local_n_num = copy_n_num      
      
      end subroutine MAKE_LOC_NODES_NUMERATION
