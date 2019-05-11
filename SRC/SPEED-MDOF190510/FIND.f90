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

!> @brief Searches for element in a vector
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] n dimension if vec
!> @param[in] vec vector
!> @param[in] tar target to find
!> @param[out] FIND  1 if tar is in vec, 0 otherwise 

      function FIND(n,vec,tar)
                                                                     
      implicit none

      integer*4 :: n, i, tar, FIND
      integer*4, dimension(n) :: vec

      do i = 1, n
        if(vec(i) .eq. tar) then
          FIND = 1
          return
        endif
      enddo

      FIND = 0 
      return   
       
      end function FIND

