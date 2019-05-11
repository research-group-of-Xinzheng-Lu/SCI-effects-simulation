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

!> @brief Checks if a an element is present on a vector and give its position.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] n dimension of the vector v
!> @param[in] v vector
!> @param[in] index to find
!> @param[out] tt 1 if ind is in v, 0 otherwise 
!> @param[out] pos  position of ind in v

    subroutine CHECK_MPI(n,v,ind,tt,pos)
    
    implicit none
    
    integer*4 :: i, ne
    integer*4, intent(in) ::n,ind
    integer*4, dimension(0:n), intent(in)  :: v
    integer*4, intent(out) :: tt, pos
    

    tt = 0
    pos = - 1
    i = 1
    ne = v(0) - 1
    
      
    do while (i.le.ne)
       if (v(v(i-1)) .eq. ind) then
          tt = 1
          pos = i
          return
       endif   
       i=i+1
       
    enddo    
    
    end subroutine CHECK_MPI
