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

!> @brief Verifies if 2 different faces are opposite or not.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nor_x  x-comp of the normal vector for the current element 
!> @param[in] nor_y  y-comp of the normal vector for the current element 
!> @param[in] nor_z  z-comp of the normal vector for the current element 
!> @param[in]  material  block id for the neigh. element
!> @param[in] element  element id for the neigh. element
!> @param[in] face  face id for the neigh. element
!> @param[in] nel_glo  number of global element dg
!> @param[in] nomalxyz  nomals for dg elements
!> @param[in] mat_el_fac id for dg elements
!> @param[out] yon  1 if the normals are opposite 0 otherwise


      subroutine CHECK_NORMAL(nor_x, nor_y, nor_z, &
                              material, element,  face, &
                              nel_glo, normalxyz, mat_el_fac, yon)
      
      implicit none
      
      integer*4 :: i
      integer*4 :: material, element, face, yon, nel_glo
      integer*4, dimension(nel_glo,3) :: mat_el_fac
      
      real*8 :: nor_x, nor_y, nor_z, res                        
      real*8, dimension(nel_glo,3) :: normalxyz
                                 
      yon = 0              
      
      do i = 1, nel_glo
      
         if(material .eq. mat_el_fac(i,1) .and. &
            element .eq. mat_el_fac(i,2)  .and. &
            face .eq. mat_el_fac(i,3) ) then
      
            res = dsqrt(dabs(nor_x + normalxyz(i,1))**2 + dabs(nor_y + normalxyz(i,2))**2 + dabs(nor_z + normalxyz(i,3))**2)
            
            if(res .le. 5.d-2) yon = 1
            
         endif
      
      enddo                     
                    
                    
                    
      end subroutine CHECK_NORMAL             
