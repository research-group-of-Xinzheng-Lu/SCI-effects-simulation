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

!> @brief Inizialization for parallel computation.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] comm mpi common world
!> @param[in] np number of mpi processes
!> @param[in] id mpi processor index
!> @param[out] mpi error

  subroutine INITIALIZATION(comm,np,id,ierr)
     
     implicit none
     include 'SPEED.MPI'
     
     integer*4 :: comm,np,id,ierr
     
     call MPI_INIT(ierr)
     call MPI_COMM_RANK(comm,id,ierr)
     call MPI_COMM_SIZE(comm,np,ierr)
     
     if (ierr.ne.0) then
        write(*,*)'MPI Initialization error - proc : ',id
     endif
     
     speed_tag = speed_tag_min
     
     return
     
     end subroutine INITIALIZATION
