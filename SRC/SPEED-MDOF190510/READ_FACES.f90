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

!> @brief Reads FACES.input.
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] filename  file name
!> @param[in] n number of DG faces 
!> @param[in] mpi_id  mpi processor identity
!> @param[in] mpi_com  mpi communicator
!> @param[out] fac material, element, face for a DG face
!> @param[out] nodes constants for the bilienar map from (-1,1)^3 to the current element

     subroutine READ_FACES(filename, fac, nodes, n, mpi_id, mpi_com)

     implicit none

     character*70 :: filename

     integer*4 :: j, n, mpi_id, mpi_com, error, unitname

     integer*4, dimension(3,n), intent(inout) :: fac

     real*8, dimension(25,n), intent(inout) :: nodes
    
     unitname = 400 + mpi_id
     open(unitname, file=filename)
     do j = 1,n
         read(unitname,*) &
                         fac(1,j), fac(2,j), fac(3,j), &
                         nodes(1,j), nodes(2,j), nodes(3,j), & 
                         nodes(4,j), nodes(5,j), nodes(6,j), & 
                         nodes(7,j), nodes(8,j), nodes(9,j), & 
                         nodes(10,j), nodes(11,j), nodes(12,j), & 
                         nodes(13,j), nodes(14,j), nodes(15,j), & 
                         nodes(16,j), nodes(17,j), nodes(18,j), & 
                         nodes(19,j), nodes(20,j), nodes(21,j), & 
                         nodes(22,j), nodes(23,j), nodes(24,j), & 
                         nodes(25,j) 

      enddo
      
      close(unitname)


      return

      end subroutine READ_FACES
