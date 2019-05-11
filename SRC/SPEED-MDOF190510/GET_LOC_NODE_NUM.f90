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


!> @brief Computes number of local nodes.
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] nnz_loc length of cs_loc
!> @param[in] cs_loc local connectivity vector
!> @param[in] id  number of mpi process
!> @param[out] nn_loc  number of local unknowns

      subroutine GET_LOC_NODE_NUM(nnz_loc, cs_loc, nn_loc, id)

      use qsort_c_module
      
      implicit none

      integer*4 :: nnz_loc
      integer*4 :: ie, ne, i, j, ind_nd, id

      integer*4, intent(out) :: nn_loc

      integer*4, dimension(0:nnz_loc), intent(in) :: cs_loc
      integer*4, dimension(0:nnz_loc) :: cs_loc_copy

            
      cs_loc_copy = cs_loc
      ne = cs_loc_copy(0) - 1
           
      
      do ie = 1,ne
         cs_loc_copy(cs_loc_copy(ie-1)) = 0       
      enddo

      !QUICKSORT FOR REORDING
      !------------------------------------------------------------
      nn_loc = 0
      do ie = 0, cs_loc_copy(0)
         cs_loc_copy(ie) = 0
      enddo         
      
      call QsortC(cs_loc_copy)
      
      
      do i = 1, nnz_loc
         if (cs_loc_copy(i) .ne. 0) then
               if(cs_loc_copy(i) .ne. cs_loc_copy(i-1)) then
                  nn_loc = nn_loc + 1
               endif                
         endif           
      enddo
            
      !-------------------------------------------------------------
  
      !-------------------------------------------------------------- 
!      do i = cs_loc_copy(0) + 1, nnz_loc
!         ind_nd = cs_loc_copy(i)
!      
!         do j = i + 1, nnz_loc
!            if (ind_nd .eq. cs_loc_copy(j)) cs_loc_copy(j) = 0          
!         enddo
!      enddo
!      
!      do ie = 0, cs_loc_copy(0)
!         cs_loc_copy(ie) = 0
!      enddo         
!      
!      
!      nn_loc = 0
!      do i = 0, nnz_loc - 1
!         if(cs_loc_copy(i) .ne. 0) nn_loc = nn_loc + 1
!      enddo   
      !-------------------------------------------------------------
      
            
      end subroutine GET_LOC_NODE_NUM
