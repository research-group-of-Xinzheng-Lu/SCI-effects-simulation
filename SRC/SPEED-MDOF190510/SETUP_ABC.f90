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

!> @brief Setup for ABC faces. Stores data on array ielem_abc. 
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] nm number of materials
!> @param[in] sd polynomial degree vector
!> @param[in] tag_mat label for materials
!> @param[in] cs_nnz_loc length of cs_loc
!> @param[in] cs_loc  local spectral connectivity vector
!> @param[in] nn_loc number of local elements
!> @param[in] local_n_num local node numeration
!> @param[in] ne_loc number of local elements
!> @param[in] local_el_num local element numeration
!> @param[in] nel_abc  number of elements where ABC conditions are applied
!> @param[in] i4count  vector identifying nodes where ABC conditions are applied
!> @param[in] ielem_abc matrix where abc elements are stored: 
!!      ielem_abc(nel_abc,1) = global id, 
!!      ielem_abc(nel_abc,2),...,ielem_abc(nel_abc,7)  indices for selecting the ABC face 


      subroutine SETUP_ABC(nm, sd, tag_mat, cs_nnz_loc, cs_loc, &
                           nn_loc, local_n_num, ne_loc, local_el_num, &
                           nel_abc, ielem_abc, i4count)

     implicit none
     
     integer*4 :: nm, cs_nnz_loc, nn_loc, ne_loc, nel_abc
     integer*4 :: im, nn, ie, ned
     integer*4 :: ne1, ne2, ne3, ne4, ic1, ic2, ic3, ic4
     
     integer*4, dimension(nm) :: tag_mat, sd
     integer*4, dimension(0:cs_nnz_loc) :: cs_loc
     integer*4, dimension(nn_loc) :: local_n_num, i4count
     integer*4, dimension(ne_loc) :: local_el_num
     integer*4, dimension(nel_abc,7):: ielem_abc

!      write(*,*) 'dentrodopo', nel_abc

      ned = cs_loc(0) - 1

      nel_abc = 0      
      
      do im = 1, nm

         nn = sd(im) +1         

         do ie = 1, ned
      
            if (cs_loc(cs_loc(ie -1) +0) .eq. tag_mat(im)) then
            
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(1 -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(1 -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(1 -1) +1)              
               
                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then
                  nel_abc = nel_abc +1
                  ielem_abc(nel_abc,1) = local_el_num(ie)
                  ielem_abc(nel_abc,2) = 1
                  ielem_abc(nel_abc,3) = 1
                  ielem_abc(nel_abc,4) = 1
                  ielem_abc(nel_abc,5) = nn
                  ielem_abc(nel_abc,6) = 1
                  ielem_abc(nel_abc,7) = nn
               endif               
               
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(1 -1) +1)
                                           
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then
                  nel_abc = nel_abc +1
                  ielem_abc(nel_abc,1) = local_el_num(ie)
                  ielem_abc(nel_abc,2) = 1
                  ielem_abc(nel_abc,3) = nn
                  ielem_abc(nel_abc,4) = 1
                  ielem_abc(nel_abc,5) = 1
                  ielem_abc(nel_abc,6) = 1
                  ielem_abc(nel_abc,7) = nn
               endif               

               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(1 -1) +1)

                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then
                  nel_abc = nel_abc +1
                  ielem_abc(nel_abc,1) = local_el_num(ie)
                  ielem_abc(nel_abc,2) = 1
                  ielem_abc(nel_abc,3) = nn
                  ielem_abc(nel_abc,4) = 1
                  ielem_abc(nel_abc,5) = nn
                  ielem_abc(nel_abc,6) = 1
                  ielem_abc(nel_abc,7) = 1
               endif               

               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(nn -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(nn -1) +1)
                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then
               
                  nel_abc = nel_abc +1
                  ielem_abc(nel_abc,1) = local_el_num(ie)
                  ielem_abc(nel_abc,2) = nn
                  ielem_abc(nel_abc,3) = nn
                  ielem_abc(nel_abc,4) = 1
                  ielem_abc(nel_abc,5) = nn
                  ielem_abc(nel_abc,6) = 1
                  ielem_abc(nel_abc,7) = nn
               endif               

               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(1 -1) +1)

               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then
                  nel_abc = nel_abc +1
                  ielem_abc(nel_abc,1) = local_el_num(ie)
                  ielem_abc(nel_abc,2) = 1
                  ielem_abc(nel_abc,3) = nn
                  ielem_abc(nel_abc,4) = nn
                  ielem_abc(nel_abc,5) = nn
                  ielem_abc(nel_abc,6) = 1
                  ielem_abc(nel_abc,7) = nn
               endif               

               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(1 -1) +1)

                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then
                   nel_abc = nel_abc +1
                  ielem_abc(nel_abc,1) = local_el_num(ie)
                  ielem_abc(nel_abc,2) = 1
                  ielem_abc(nel_abc,3) = nn
                  ielem_abc(nel_abc,4) = 1
                  ielem_abc(nel_abc,5) = nn
                  ielem_abc(nel_abc,6) = nn
                  ielem_abc(nel_abc,7) = nn
               endif               
            endif
         enddo
      enddo



      return

      end subroutine SETUP_ABC
