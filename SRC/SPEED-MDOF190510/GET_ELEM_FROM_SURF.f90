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


!> @brief Find element index from verteces number.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nnz_bc length of cs_bc
!> @param[in] cs_bc spectral connectivity vector
!> @param[in] n1  node index of 1 vertex of the element (quad)
!> @param[in] n2  node index of 1 vertex of the element (quad)
!> @param[in] n3  node index of 1 vertex of the element (quad)
!> @param[in] n4  node index of 1 vertex of the element (quad)
!> @param[out] ie_surf index of the element found

     subroutine GET_ELEM_FROM_SURF(nnz_bc, cs_bc, n1, n2, n3, n4, ie_surf)
     

     implicit none

     integer*4 :: i, ne, nel, find
     integer*4 :: nnz_bc, n1, n2, n3, n4, ie_surf
     
     integer*4, dimension(:), allocatable :: vet
     integer*4, dimension(0:nnz_bc) :: cs_bc
     

     
     ne = cs_bc(0) - 1 
     
     do i = 1, ne
        nel = cs_bc(i) - cs_bc(i-1)  - 1 
        allocate(vet(nel))
        
        vet = cs_bc(cs_bc(i-1) + 1 : cs_bc(i) - 1 )
        
        if (find(nel, vet, n1) .eq. 1 .and. &
            find(nel, vet, n2) .eq. 1 .and. &
            find(nel, vet, n3) .eq. 1 .and. &
            find(nel, vet, n4) .eq. 1 ) then
           
            ie_surf = i
            return
        endif    
        
        deallocate(vet)   
     enddo
     
     ie_surf = 0
     return
     
     
     end subroutine GET_ELEM_FROM_SURF
