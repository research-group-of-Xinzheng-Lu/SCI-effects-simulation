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


!> @brief Fills array check_ns for explosive force.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] cs_nnz length of cs  
!> @param[in] cs spectral connectivity vector
!> @param[in] nm number of materials
!> @param[in] tm label for materias 
!> @param[in] sd polynomial degree vector
!> @param[in] nl_expl  number of explosive nodes
!> @param[in] num_ne  number of explosive nodes
!> @param[in] max_num_ne  max number of explosive nodes
!> @param[in] sour_ne explosive source
!> @param[in] dist_sour_ne  distance of the node from the explosive source
!> @param[in] length_cne  length check explosive nodes (useless)
!> @param[in] fun_expl  function associeted with the explosive load
!> @param[in] nf  number of functions
!> @param[in] tag_func label for explosive functions
!> @param[in] val_expl values for the explosive loads
!> @param[in] nn_loc number of local nodes
!> @param[in] local_n_num  local numeration vector
!> @param[out] check_ne  check_ne(i,1)  explosive source node, 
!!      check_ne(i,2)  explosive function for the node i,
!!      check_ne(i,3)  i for explosive load,
!!      check_ne(i,4)  number of the element       
!> @param[out]  check_dist_ne check_dist_ne(i,1) = dist_sour_ne / val_expl

      subroutine CHECK_EXPL(cs_nnz,cs,&
                          nm,tm,sd,&
                          nl_expl,&
                          num_ne,max_num_ne,&
                          sour_ne,dist_sour_ne,&
                          check_ne,check_dist_ne,&
                          length_cne,&
                          fun_expl,nf,tag_func,val_expl,&
                          nn_loc, local_n_num)

      implicit none

      integer*4 :: cs_nnz,nm,ne,nl_expl,nf, nn_loc
      integer*4 :: max_num_ne,im,ie,iexpl,nn,is,in,ip,i,j,k,h,fn
      integer*4 :: length_cne
      
      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(nm) :: tm,sd
      integer*4, dimension(nl_expl) :: num_ne
      integer*4, dimension(nl_expl) :: fun_expl
      integer*4, dimension(nn_loc) :: local_n_num
      integer*4, dimension(nf) :: tag_func

      integer*4, dimension(length_cne,4) :: check_ne
      integer*4, dimension(max_num_ne,nl_expl) :: sour_ne
      
      real*8 :: vel_prop 

      real*8, dimension(:), allocatable :: ct,ww

      real*8, dimension(:,:), allocatable :: dd
      real*8, dimension(nl_expl,20) :: val_expl
      real*8, dimension(max_num_ne,nl_expl) :: dist_sour_ne
      real*8, dimension(length_cne,1) :: check_dist_ne


      
      h = 0
      nn = 2
      ne = cs(0) -1

      allocate(ct(nn),ww(nn),dd(nn,nn))
      call MAKE_LGL_NW(nn,ct,ww,dd)
      
      do im = 1,nm
         if ((sd(im) +1).ne.nn) then
            deallocate(ct,ww,dd)
            nn = sd(im) +1
            allocate(ct(nn),ww(nn),dd(nn,nn))
            call MAKE_LGL_NW(nn,ct,ww,dd)
         endif
         
         do ie = 1,ne
            if (cs(cs(ie -1) +0).eq.tm(im)) then
                do k = 1,nn
                   do j = 1,nn
                      do i = 1,nn
                         is = nn*nn*(k -1) +nn*(j -1) +i
                         in = cs(cs(ie -1) +is)

                         do iexpl = 1,nl_expl
                             do ip = 1,num_ne(iexpl)

                               if (local_n_num(in) .eq. sour_ne(ip,iexpl)) then

                                  h = h + 1
                                  check_ne(h,1) = sour_ne(ip,iexpl)
                                  check_ne(h,2) = fun_expl(iexpl)
                                  check_ne(h,3) = iexpl
                                  check_ne(h,4) = ie
                                  check_dist_ne(h,1) = dist_sour_ne(ip,iexpl) / val_expl(iexpl,19)
                                     
                               endif
                             enddo  !ip
                          enddo  !iexpl
                       enddo  !i
                    enddo  !j
                  enddo  !k
             

             endif  !if (cs....)

            enddo  !ie = 1,ne

          enddo !im = 1,nm 

      return

      end subroutine CHECK_EXPL

