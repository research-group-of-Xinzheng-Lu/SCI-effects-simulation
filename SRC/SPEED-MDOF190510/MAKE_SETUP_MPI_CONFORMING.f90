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

!> @brief Defines buffers for MPI communication.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0

     subroutine MAKE_SETUP_MPI_CONFORMING()
     
      use max_var
      use speed_par
      

      implicit none
 
      include 'SPEED.MPI'      


      nsend = 0
      nrecv = 0
      
      allocate(proc_send(mpi_np))
      allocate(proc_recv(mpi_np))

      !if(mpi_id .eq.0) then
      !   start1=MPI_WTIME()
      !endif
      
      call SETUP_MPI(nnod,node_domain,con_nnz_loc,con_spx_loc, &
                      nsend,node_send,nrecv,node_recv, &
                      mpi_np,proc_send,proc_recv,mpi_id,mpi_file,mpi_comm)    

                
      allocate(node_send(nsend))
      allocate(node_recv(nrecv))
      
      call SETUP_MPI(nnod,node_domain,con_nnz_loc,con_spx_loc, &
                      nsend,node_send,nrecv,node_recv, &
                      mpi_np,proc_send,proc_recv,mpi_id,mpi_file,mpi_comm)   
      
      deallocate(node_domain)
      call MPI_BARRIER(mpi_comm, mpi_ierr)
     
     ! if(mpi_id .eq.0) then
     !    start2=MPI_WTIME()
     !    write (*,*) "Time in routine setup_mpi ",start2-start1,"sec"
     ! endif


     end subroutine MAKE_SETUP_MPI_CONFORMING
