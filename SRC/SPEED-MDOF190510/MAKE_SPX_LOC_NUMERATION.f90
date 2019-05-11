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


!> @brief Makes local numeration of spectral nodes.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nnloc number of local nodes
!> @param[in] loc_n_num local node numeration
!> @param[in] cs_nnz length of cs_loc
!> @param[in] cs_nnz_bc length of cs_bc_loc
!> @param[in,out] cs_loc  spectral connectivity vector (IN -> global num, OUT -> local num.)
!> @param[in,out] cs_bc_loc  spectral boundary connectivity vector 
!!                                 (IN -> global num, OUT -> local num.)


      subroutine MAKE_SPX_LOC_NUMERATION(nnloc, loc_n_num, &
                                 cs_nnz, cs_loc, &
                                 cs_nnz_bc, cs_bc_loc)


      implicit none
      
      integer*4 :: cs_nnz, cs_nnz_bc, nnloc
      integer*4, dimension(0:cs_nnz), intent(inout) :: cs_loc
      integer*4, dimension(0:cs_nnz_bc), intent(inout) :: cs_bc_loc
      integer*4, dimension(nnloc), intent(inout) :: loc_n_num

      integer*4 :: ne_loc, ne_loc_bc       
      integer*4 :: ie, j, iglo, iloc
      
      
      ne_loc = cs_loc(0) - 1
      
      do ie = 1, ne_loc
         do j = cs_loc(ie -1) + 1, cs_loc(ie) - 1
            
            iglo = cs_loc(j)
            call GET_INDLOC_FROM_INDGLO(loc_n_num, nnloc, iglo, iloc)       
            cs_loc(j) = iloc
            
            if(iloc .eq. 0) write(*,*) '1 Error in MAKE_SPX_LOC_NUMERATION'
         enddo
      enddo      
      
      ne_loc_bc = cs_bc_loc(0) - 1
      
      do ie = 1, ne_loc_bc
         do j = cs_bc_loc(ie -1) + 1, cs_bc_loc(ie) - 1
            
            iglo = cs_bc_loc(j)
            call GET_INDLOC_FROM_INDGLO(loc_n_num, nnloc, iglo, iloc)       
            cs_bc_loc(j) = iloc
            
            if(iloc .eq. 0) write(*,*) '2 Error in MAKE_SPX_LOC_NUMERATION'
         enddo
      enddo      
      
      
      
      end subroutine MAKE_SPX_LOC_NUMERATION
      
