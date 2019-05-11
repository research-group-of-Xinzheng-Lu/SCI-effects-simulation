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

!> @brief Find the face (quad) of a DG element (hex) 
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] faces faces(1,i) = block id of the element
!!                             faces(2,i) = element id
!!                             faces(3,i) = face id    
!!
!> @param[in] nel_dg_glo number of dg elements
!> @param[in] ik element id to find
!> @param[in] ind  face id to find
!> @param[out]  face_found  0 if not found, number of the face (from 1 to 6) if found


    subroutine GET_FACE_DG(faces, nel_dg_glo, ik, ind, face_found)


    implicit none
    
    integer*4 :: nel_dg_glo, ind, ik, face_found, i

    integer*4, dimension(3,nel_dg_glo) :: faces
    
    face_found = 0 
    do i = 1, nel_dg_glo
    
       if(faces(2,i) .eq. ik .and. faces(3,i) .eq. ind ) then
          face_found = i
          return
       endif
    enddo     
    
        
    
    end subroutine GET_FACE_DG
