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

!> @brief Find element position on a vector.
!! @author Ilario Mazzieri
!> @date September, 2013 - Creation
!> @version 1.0
!> @param[in] n
!> @param[in] vec vec vector where seraching tar  
!> @param[in] tar target
!> @param[ount] ind position of tar in vec, 0 otherwise

      subroutine FIND_POS(n,vec,tar,ind)
                                                                     
      implicit none

      integer*4 :: n, i, tar, ind
      integer*4, dimension(n) :: vec

      ind = 0  
      do i = 1, n
        if(vec(i) .eq. tar) then
          ind = i
          return
        endif
      enddo


      return   
       
      end subroutine FIND_POS

