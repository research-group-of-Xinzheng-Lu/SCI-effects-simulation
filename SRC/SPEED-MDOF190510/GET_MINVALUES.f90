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

!> @brief Computes positions of minimum values. 
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] n_glo number of global indices
!> @param[in] i_glo vector conteining global indices
!> @param[in] v_glo vector conteining global values
!> @param[in] np number of MPI processors
!> @param[in] n_loc number of local indices
!> @param[out] i_loc vector containing local indices corresponding to global minimal values


     subroutine  GET_MINVALUES(i_glo, v_glo, n_glo, i_loc, n_loc, np)
     
     implicit none
     
     integer*4 :: n_glo, n_loc, np
     integer*4 :: i,j

     integer*4, dimension(1) :: pos
     integer*4, dimension(np) :: ind     
     integer*4, dimension(n_glo) :: i_glo
     integer*4, dimension(n_loc) :: i_loc
          
     real*8, dimension(n_glo) :: v_glo
     real*8, dimension(n_loc) :: v_loc
     real*8, dimension(np) :: val
     

     do i = 1, n_loc
        do j = 1, np
          
          ind(j) = i_glo(n_loc*(j-1) + i)
          val(j) = v_glo(n_loc*(j-1) + i)
        enddo  
        
       pos = minloc(val)
       i_loc(i) = ind(pos(1))
            
     enddo         

     end subroutine  GET_MINVALUES
