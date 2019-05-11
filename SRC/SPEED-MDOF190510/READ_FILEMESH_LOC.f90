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

!> @brief Reads the local mesh.
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] file_grid file name
!> @param[in] nn_loc number of local nodes
!> @param[in] loc_n_num local node numbering
!> @param[out] xx_loc x-coordinate of local nodes
!> @param[out] yy_loc y-coordinate of local nodes
!> @param[out] zz_loc z-coordinate of local nodes
 
      subroutine READ_FILEMESH_LOC(file_grid, nn_loc, xx_loc, yy_loc, zz_loc, &
                            loc_n_num)

      implicit none
      
      character*70 :: file_grid
      character*110 :: input_line

      integer*4 :: nn_loc
      integer*4 :: inode
      integer*4 :: nnode,nelem,i,ic,status

      integer*4, dimension(nn_loc) :: loc_n_num

      real*8 :: xx,yy,zz      
      real*8, dimension(nn_loc), intent(inout) :: xx_loc(nn_loc), yy_loc(nn_loc), zz_loc(nn_loc)
      
      inode = 0
      status = 0 
      
      open(23,file=file_grid)
      
      do 
         read(23,'(A)') input_line
         if (input_line(1:1) .ne. '#') exit
      enddo
      
      read(input_line,*) nnode, nelem
      
      do i = 1,nnode
        read(23,*) inode, xx, yy, zz
        if (inode.ne.i) then
          status = 1
        endif
        
        call GET_INDLOC_FROM_INDGLO(loc_n_num, nn_loc, inode, ic)
        
        
        if (ic .ne. 0) then
          xx_loc(ic) = xx
          yy_loc(ic) = yy
          zz_loc(ic) = zz
        endif
        
      enddo
      
      
      close(23)
      
      return
      
      end subroutine READ_FILEMESH_LOC

