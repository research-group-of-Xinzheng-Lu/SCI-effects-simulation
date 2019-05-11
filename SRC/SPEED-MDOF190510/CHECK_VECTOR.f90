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

!> @brief Checks if an element is present in a vector.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] n dimension of v
!> @param[in] v vector
!> @param[in] ind index to find
!> @param[out] tt 1 if ind is in v, 0 otherwise 

    subroutine CHECK_VECTOR(n,v,ind,tt)
    
    implicit none
    
    integer*4 :: i
    
    integer*4, intent(in) ::n,ind
    integer*4, dimension(n), intent(in)  :: v
    integer*4, intent(out) :: tt

    tt = 0
    i = 1
    
    do while (i.le.n)
       if (v(i) .eq. ind) then
          tt = 1
          return
       endif   
       i=i+1
    enddo    
    
    end subroutine CHECK_VECTOR
