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

!> @brief Reads LS.input file.
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] filec file name
!> @param[in] num_nodes  number of monitored nodes 
!> @param[out] x x-coord. of the monitor 
!> @param[out] y y-coord. of the monitor 
!> @param[out] z z-coord. of the monitor 

      subroutine READ_FILEPG(filec,num_nodes, x,y,z)

      implicit none

      character*70 :: filec
      character*100000 :: input_line

      integer*4 :: i,trash, num_nodes
      integer*4 :: ileft,iright, status
      
      real*8 :: length_edge1,length_edge2,length_edge3

      real*8, dimension(num_nodes) :: x,y,z


      open(20,file=filec)
      read(20,'(A)',IOSTAT = status) input_line
      ileft = 1
      iright = len(input_line)
      read(input_line(ileft:iright),*)num_nodes 

      do i = 1,num_nodes
        read(20,'(A)',IOSTAT = status) input_line
        if (status.ne.0) exit
        ileft = 1
        iright = len(input_line)
        read(input_line(ileft:iright),*)trash,x(i),y(i),z(i) 
      enddo
                    
      close(20)

      return

      end subroutine READ_FILEPG

