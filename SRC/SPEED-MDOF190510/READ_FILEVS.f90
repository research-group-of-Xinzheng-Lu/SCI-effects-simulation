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

!> @brief Reads vs30 and thickness from VS_RS.out 
!! @author Ilario Mazzieri
!> @date September, 2015
!> @version 1.0
!> @param[in] filec file name for the mesh
!> @param[in] num_tria  number of triangular elements
!> @param[out] velocity vs30 
!> @param[out] thickness thickness of soft sediments

      subroutine READ_FILEVS(filec,num_tria,velocity,sediments)

      implicit none

      integer*4 :: num_tria
      integer*4 :: i
      integer*4 :: trash

      character*70 :: filec
      real*8, dimension(num_tria) :: velocity, sediments

      open(20,file=filec)

      do i = 1,num_tria
        read(20,*) trash, velocity(i), sediments(i)
      enddo
                    
      close(20)

      return

      end subroutine READ_FILEVS

