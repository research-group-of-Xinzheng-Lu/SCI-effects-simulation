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

!> @brief Local numeration for DG elements.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nnloc number of local elements
!> @param[in] loc_n_num local node numeration
!> @param[in] neloc number of local elements
!> @param[in] loc_el_num local element numeration
!> @param[in,out] cs_nnz length of cs_loc 
!> @param[in,out] cs_loc DG spectral connectivity vector (IN -> global num, OUT -> local num.)


      subroutine MAKE_DG_LOC_NUMERATION(nnloc, loc_n_num, &
                                 cs_nnz, cs_loc, neloc, loc_el_num, mpi_id)

      implicit none
      
      integer*4 :: ne_loc, ieloc, neloc       
      integer*4 :: ie, j, iglo, iloc
      integer*4 :: cs_nnz, nnloc, mpi_id
      
      integer*4, dimension(0:cs_nnz), intent(inout) :: cs_loc
      integer*4, dimension(nnloc), intent(inout) :: loc_n_num
      integer*4, dimension(neloc), intent(inout) :: loc_el_num
      
      
      ne_loc = cs_loc(0) - 1

      
      do ie = 1, ne_loc

         call GET_INDLOC_FROM_INDGLO(loc_el_num, neloc, cs_loc(cs_loc(ie-1)), ieloc)

         if(ieloc .ne. 0)  then
         
            cs_loc(cs_loc(ie-1)) = 0   
            do j = cs_loc(ie -1) + 1, cs_loc(ie) - 1
               cs_loc(j) = 0  
            enddo
         
         else   
               
            do j = cs_loc(ie -1) + 1, cs_loc(ie) - 1
               iglo = cs_loc(j)
               call GET_INDLOC_FROM_INDGLO(loc_n_num, nnloc, iglo, iloc)       
               cs_loc(j) = iloc
            enddo
         
         endif      
         
      enddo      
      
     
      
      end subroutine MAKE_DG_LOC_NUMERATION
