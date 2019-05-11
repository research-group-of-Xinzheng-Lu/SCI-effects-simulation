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

!> @brief Computes the highest node of the mesh (z-dir)
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] nn_loc  number of local nodes   
!> @param[in] ne_loc  number of local element
!> @param[in] zz_loc  z-coordinate for the nodes
!> @param[in] local_n_num vector for local node numbering 
!> @param[in] nnz_loc length of cs_loc
!> @param[in] cs_loc local spectral connectivity
!> @param[in] nm number of materials
!> @param[in] sd polynomial degree vector 
!> @param[in] tm labels for material vector 


      subroutine GET_HIGHEST_NODE(nn_loc, ne_loc, zz_loc, loc_n_num, nnz_loc, cs_loc, &
                                  nm, tm, sd, highest)


      implicit none
      
      integer*4 :: nn_loc, nnz_loc, ne_loc
      integer*4 :: ie, nn, nm, im
      integer*4 :: n1,n2,n3,n4,n5,n6,n7,n8 
      integer*4 :: ic1,ic2,ic3,ic4,ic5,ic6,ic7,ic8

      integer*4, dimension(nn_loc) :: loc_n_num
      integer*4, dimension(nm) :: tm
      integer*4, dimension(nm) :: sd

      integer*4, dimension(0:nnz_loc) :: cs_loc

      real*8 :: zz1,zz2,zz3,zz4,zz5,zz6,zz7,zz8

      real*8, dimension(:), allocatable :: ct,ww
      real*8, dimension(nn_loc) :: zz_loc
      real*8, dimension(ne_loc) :: highest 
      
      real*8, dimension(:,:), allocatable :: dd
      
      
      nn = 2
      allocate(ct(nn),ww(nn),dd(nn,nn))
      call MAKE_LGL_NW(nn,ct,ww,dd)

      do im = 1,nm
         if ((sd(im) +1).ne.nn) then
            deallocate(ct,ww,dd)
            nn = sd(im) +1
            allocate(ct(nn),ww(nn),dd(nn,nn))
            call MAKE_LGL_NW(nn,ct,ww,dd)
         endif
         
         do ie = 1,ne_loc
            if (cs_loc(cs_loc(ie -1) +0).eq.tm(im)) then

                n1 = nn*nn*(1 -1) +nn*(1 -1) +1
                n2 = nn*nn*(1 -1) +nn*(1 -1) +nn
                n3 = nn*nn*(1 -1) +nn*(nn -1) +nn
                n4 = nn*nn*(1 -1) +nn*(nn -1) +1
                n5 = nn*nn*(nn -1) +nn*(1 -1) +1
                n6 = nn*nn*(nn -1) +nn*(1 -1) +nn
                n7 = nn*nn*(nn -1) +nn*(nn -1) +nn
                n8 = nn*nn*(nn -1) +nn*(nn -1) +1
               
                ic1 = cs_loc(cs_loc(ie -1) +n1)
                ic2 = cs_loc(cs_loc(ie -1) +n2)
                ic3 = cs_loc(cs_loc(ie -1) +n3)
                ic4 = cs_loc(cs_loc(ie -1) +n4)
                ic5 = cs_loc(cs_loc(ie -1) +n5)
                ic6 = cs_loc(cs_loc(ie -1) +n6)
                ic7 = cs_loc(cs_loc(ie -1) +n7)
                ic8 = cs_loc(cs_loc(ie -1) +n8)
                    
                zz1 = zz_loc(ic1)
                zz2 = zz_loc(ic2)
                zz3 = zz_loc(ic3)
                zz4 = zz_loc(ic4)
                zz5 = zz_loc(ic5)
                zz6 = zz_loc(ic6)
                zz7 = zz_loc(ic7)
                zz8 = zz_loc(ic8)

                highest(ie) = max(zz1,zz2,zz3,zz4,zz5,zz6,zz7,zz8)
                
             endif  
         enddo 
      enddo     


      deallocate(ct,ww,dd)
      
      return
      
      end subroutine GET_HIGHEST_NODE
