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

!> @brief Returns local id from global id.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] local_el local element numeration local_el(i) = global numeration of the 
!!                           i-th local element
!> @param[in] nel_loc  number of local elements
!> @param[in] ie global id for the element  
!> @param[out] ic  local number for the element, 0 if the element is not present

     subroutine GET_INDLOC_FROM_INDGLO(local_el, nel_loc, ie, ic)
     

     implicit none

     integer*4 :: nel_loc, ie, ic, i
     integer*4, dimension(nel_loc) :: local_el
     
     do i = 1, nel_loc
        if(ie .eq. local_el(i))  then
           ic = i
           return
        endif   
       
     enddo

     ic = 0     
     return
    
     end subroutine GET_INDLOC_FROM_INDGLO
