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

!> @brief Makes local connectivity vector.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nel_loc  number of local elements
!> @param[in] nel  number of total elements
!> @param[in] el_dom  metis output (el_dom(i) = k --> element i belongs to mpi process k)
!> @param[in] tag_mat material label
!> @param[in] sdeg_mat  polynomial degree vector
!> @param[in] mpi_id  mpi process id
!> @param[out] cs_loc local spectral connectivity vector

   subroutine MAKE_SPX_CON_LOC(nel_loc, nel, el_dom, &
                                       nnz,  cs, nnz_loc, cs_loc, &
                                       nm, tag_mat, sdeg_mat, mpi_id)
                                  

    implicit none

    integer*4 :: j,k,ie,nm, nn
    integer*4 :: nel_loc, nel, nnz, mpi_id
    integer*4 :: nnz_loc

    integer*4, dimension(0:nnz) :: cs
    integer*4, dimension(nel) :: el_dom
    integer*4, dimension(nm) :: tag_mat, sdeg_mat
    integer*4, dimension(0:nnz_loc) :: cs_loc
    

    cs_loc = 0
    cs_loc(0) = nel_loc+1
    k = 0           
    do ie = 1, nel
         do j = 1, nm
            if (( tag_mat(j).eq. cs(cs(ie-1)) ) .and. (el_dom(ie) .eq. mpi_id)) then 
                nn = sdeg_mat(j) +1
                k = k + 1
                cs_loc(k) = cs_loc(k-1) + nn*nn*nn +1
                cs_loc(cs_loc(k-1)) = tag_mat(j)
                cs_loc(cs_loc(k-1) + 1 : cs_loc(k)-1) = cs(cs(ie-1)+1 : cs(ie)-1)   
            endif          
         enddo
      enddo
                                  
                                  
    end subroutine MAKE_SPX_CON_LOC                                  

