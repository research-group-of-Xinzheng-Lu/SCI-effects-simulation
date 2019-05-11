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

!> @brief Findes index of a boundary element. 
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nfac number of boundary elements
!> @param[in] cn_bc connectivity for boundary elements
!> @param[in] v1  index of 1-vertex
!> @param[in] v2  index of 2-vertex
!> @param[in] v3  index of 3-vertex
!> @param[in] v4  index of 4-vertex  
!> @param[in] tag_dg label for dg elements 
!> @param[in] nl_dg  number of dg labels
!> @param[out] ind  index of the element found

      subroutine GET_TAG_BC(cn_bc, nfac, &
                             v1, v2, v3, v4, &
                             tag_dg, nl_dg, ind)

      implicit none
      
      integer*4 :: nfac, v1, v2, v3, v4, nl_dg, ind
      integer*4 :: i, j
      
      integer*4 :: tag_dg(nl_dg)
      
      integer*4, dimension(nfac,5) :: cn_bc


      ind = 0
      
      do i = 1, nfac
      
         !ordering 1 2 3 4
         if(cn_bc(i,2) .eq. v1 .and. & 
            cn_bc(i,3) .eq. v2 .and. &
            cn_bc(i,4) .eq. v3 .and. &
            cn_bc(i,5) .eq. v4) then
            
               do j = 1, nl_dg
                  if(cn_bc(i,1) .eq. tag_dg(j)) then
                     ind = j
                     return
                  endif 
               enddo
            
         endif   
         
         !ordering 2 3 4 1         
         if(cn_bc(i,2) .eq. v2 .and. & 
            cn_bc(i,3) .eq. v3 .and. &
            cn_bc(i,4) .eq. v4 .and. &
            cn_bc(i,5) .eq. v1) then
            
               do j = 1, nl_dg
                  if(cn_bc(i,1) .eq. tag_dg(j)) then
                     ind = j
                     return
                  endif 
               enddo
            
         endif 
         
         !ordering 3 4 1 2
         if(cn_bc(i,2) .eq. v3 .and. & 
            cn_bc(i,3) .eq. v4 .and. &
            cn_bc(i,4) .eq. v1 .and. &
            cn_bc(i,5) .eq. v2) then
            
               do j = 1, nl_dg
                  if(cn_bc(i,1) .eq. tag_dg(j)) then
                     ind = j
                     return
                  endif 
               enddo
            
         endif 
         
         !ordering 4 1 2 3
         if(cn_bc(i,2) .eq. v4 .and. & 
            cn_bc(i,3) .eq. v1 .and. &
            cn_bc(i,4) .eq. v2 .and. &
            cn_bc(i,5) .eq. v3) then
            
               do j = 1, nl_dg
                  if(cn_bc(i,1) .eq. tag_dg(j)) then
                     ind = j
                     return
                  endif 
               enddo
            
         endif 
         
         !ordering 1 4 3 2
         if(cn_bc(i,2) .eq. v1 .and. & 
            cn_bc(i,3) .eq. v4 .and. &
            cn_bc(i,4) .eq. v3 .and. &
            cn_bc(i,5) .eq. v2) then
            
               do j = 1, nl_dg
                  if(cn_bc(i,1) .eq. tag_dg(j)) then
                     ind = j
                     return
                  endif 
               enddo
            
         endif 

         !ordering 2 1 4 3
         if(cn_bc(i,2) .eq. v2 .and. & 
            cn_bc(i,3) .eq. v1 .and. &
            cn_bc(i,4) .eq. v4 .and. &
            cn_bc(i,5) .eq. v3) then
            
               do j = 1, nl_dg
                  if(cn_bc(i,1) .eq. tag_dg(j)) then
                     ind = j
                     return
                  endif 
               enddo
            
         endif 

         !ordering 3 2 1 4
         if(cn_bc(i,2) .eq. v3 .and. & 
            cn_bc(i,3) .eq. v2 .and. &
            cn_bc(i,4) .eq. v1 .and. &
            cn_bc(i,5) .eq. v4) then
            
               do j = 1, nl_dg
                  if(cn_bc(i,1) .eq. tag_dg(j)) then
                     ind = j
                     return
                  endif 
               enddo
            
         endif 

         !ordering 4 3 2 1
         if(cn_bc(i,2) .eq. v4 .and. & 
            cn_bc(i,3) .eq. v3 .and. &
            cn_bc(i,4) .eq. v2 .and. &
            cn_bc(i,5) .eq. v1) then
            
               do j = 1, nl_dg
                  if(cn_bc(i,1) .eq. tag_dg(j)) then
                     ind = j
                     return
                  endif 
               enddo
            
         endif 
         
         
            
      enddo


      end subroutine GET_TAG_BC
