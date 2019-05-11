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
     
!> @brief Routine used to setup the communication buffer structure.
!> @note If run with nsend and nrecv equal to 
!!     zero gives back the correct values of these to allocate node_send and node_recv, otherwise 
!!     fills node_send and node_recv with the indices of the nodes to be sent to and received from 
!!     the other processors, and proc_send and proc_recv with the number of nodes to be sent and
!!     received.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0      
!> @param[in] nnode  number of nodes
!> @param[in] node_proc  vector containing the processor to which the node is assigned
!> @param[in] ncs  length of the connectivity vector
!> @param[in] cs  connectivity vector
!> @param[in] id  processor identifyer
!> @param[in] filempi folder where to find file *.mpi
!> @param[in] mpi_comm mpi common world
!> @param[out] nsend number of nodes to be sent to other processors
!> @param[out] node_send   ordered vector of nodes to be sent
!> @param[out] nrecv       number of nodes to be received from other procs
!> @param[out] node_recv   ordered vector of nodes to be received
!> @param[out] nproc       number of processors
!> @param[out] proc_send   vector containing the number of the nodes to be sent to each proc.
!> @param[out] proc_recv   vector containing the number of the nodes to be received from each proc.
      
     subroutine SETUP_MPI(nnode,node_proc,ncs,cs,&
                           nsend,node_send,nrecv,node_recv,&
                           nproc,proc_send,proc_recv,id,filempi,mpi_comm)  
     
     
     implicit none
     include 'SPEED.MPI'
     
     character*70 :: filename, filempi, filename_new

     integer*4 :: nnode,ncs,nsend,nrecv,nproc,id,mpi_ierr  
     integer*4 :: i,j,k,ie,ip,ic,unitfile,ncs_mpi,nelem_mpi,nelem_loc, mpi_comm

     integer*4, dimension(:), allocatable :: cs_mpi 
     integer*4, dimension(:), allocatable :: i4vec
     integer*4, dimension(nnode) :: node_proc
     integer*4, dimension(0:ncs) :: cs
     integer*4, dimension(nsend) :: node_send
     integer*4, dimension(nrecv) :: node_recv
     integer*4, dimension(nproc) :: proc_send, proc_recv
     
     !real*8 :: start1,start2
     allocate(i4vec(ncs))
     nelem_loc = cs(0) - 1

     
     do ip = 1, nproc
        if ((ip -1).eq.id) then
           proc_recv(ip) = 0
        else

           i4vec = 0
           ic = 0
           do ie = 1,nelem_loc
                 do i = cs(ie -1) +1, cs(ie) -1
                    if (node_proc(cs(i)).eq.(ip -1)) then
                       ic = ic +1
                       i4vec(ic) = cs(i)
                    endif
                 enddo
           enddo
           
           !ordering nodes
           do i = 1, ic
              do j = i, ic
                 if (i4vec(j).lt.i4vec(i)) then
                    k = i4vec(i)
                    i4vec(i) = i4vec(j)
                    i4vec(j) = k
                 endif
              enddo
           enddo
          
           !counting nodes to receive
           j = 1
           do i = 2, ic
              if (i4vec(i).ne.i4vec(j)) then
                 j = j +1
                 i4vec(j) = i4vec(i)
              endif
           enddo
           
           if (ic.eq.0) j = 0
           proc_recv(ip) = j
           
           !fill node_recv vector
           if (nrecv.ne.0) then
              ic = 0
              do i = 1, ip -1
                 ic = ic +proc_recv(i)
              enddo
              
              do i = 1, proc_recv(ip)
                 node_recv(ic +i) = i4vec(i)
              enddo
           endif
         
        endif
     enddo
     
     if (nrecv.eq.0) then
        do i = 1, nproc
           nrecv = nrecv +proc_recv(i)
        enddo
     endif
 
     deallocate(i4vec) 
 
     do ip = 1, nproc
     
      ! if(id .eq. 0) then

     !   if(id .eq. ip-1) start1=MPI_WTIME()
      !    filename = 'cons000000.mpi'
      !    unitfile = 40                                 
      !    if (ip-1 .lt. 10) then                                        
      !        write(filename(10:10),'(i1)') ip-1                
      !    else if (ip-1 .lt. 100) then                                
      !        write(filename(9:10),'(i2)') ip-1
      !    else if (ip-1 .lt. 1000) then                                
      !        write(filename(8:10),'(i3)') ip-1
      !    else if (ip-1 .lt. 10000) then                                
      !        write(filename(7:10),'(i4)') ip-1
      !    else if (ip-1 .lt. 100000) then                                
      !        write(filename(6:10),'(i5)') ip-1
      !    else if (ip-1 .lt. 1000000) then                                
      !        write(filename(5:10),'(i6)') ip-1
      !    endif                                                              

      !    if(len_trim(filempi) .ne. 70) then                                                                                  
      !       filename_new = filempi(1:len_trim(filempi)) // '/' // filename
      !    else
      !       filename_new = filename    
      !    endif

      !    open(unitfile,file=filename_new)        
      !    read(unitfile,*) ncs_mpi
      !    allocate(cs_mpi(0:ncs_mpi))

                      
      !    do i = 0, ncs_mpi
      !               read(unitfile,*) cs_mpi(i)
      !    enddo
      !    close(unitfile)                
       
      ! endif
       
       !call MPI_SCATTER(ncs_mpi, 1, SPEED_INTEGER, ncs_mpi, 1, SPEED_INTEGER, 0, mpi_comm)
       !call MPI_BCAST(ncs_mpi,1,SPEED_INTEGER,0,MPI_COMM_WORLD,mpi_ierr) 
       if(id .eq. ip-1) ncs_mpi = ncs; 
       call MPI_BCAST(ncs_mpi,1,SPEED_INTEGER,ip-1,MPI_COMM_WORLD,mpi_ierr) 
       
       !if(id .ne. 0)  allocate(cs_mpi(0:ncs_mpi))
       allocate(cs_mpi(0:ncs_mpi))
       if(id .eq. ip-1) cs_mpi = cs;
       
       call MPI_BCAST(cs_mpi,ncs_mpi+1,SPEED_INTEGER,ip-1,MPI_COMM_WORLD,mpi_ierr)

       !call MPI_BCAST(cs_mpi,ncs_mpi+1,SPEED_INTEGER,0,MPI_COMM_WORLD,mpi_ierr)
       !call MPI_SCATTER(cs_mpi,ncs_mpi+1, SPEED_INTEGER, cs_mpi, ncs_mpi+1, SPEED_INTEGER, 0, mpi_comm)
            
     !  if(id .eq. ip-1) then
     !   start2=MPI_WTIME()
     !   write (*,*) "Time to send all .mpi files to all proc",start2-start1,"sec"
     !  endif
       
       
       allocate(i4vec(ncs_mpi))
       nelem_mpi = cs_mpi(0)-1
       
       
             
       if ((ip -1).eq.id) then
           proc_send(ip) = 0
       else
           i4vec = 0

           ic = 0
           do ie = 1,nelem_mpi
                 do i = cs_mpi(ie -1) +1, cs_mpi(ie) -1
                    if (node_proc(cs_mpi(i)).eq.id) then   
                       ic = ic +1
                       i4vec(ic) = cs_mpi(i)
                    endif
                 enddo
           enddo
           
           do i = 1, ic
              do j = i, ic
                 if (i4vec(j).lt.i4vec(i)) then
                    k = i4vec(i)
                    i4vec(i) = i4vec(j)
                    i4vec(j) = k
                 endif
              enddo
           enddo
           
           j = 1
           do i = 2, ic
              if (i4vec(i).ne.i4vec(j)) then
                 j = j +1
                 i4vec(j) = i4vec(i)
              endif
           enddo
           
           if (ic .eq. 0) j = 0
           proc_send(ip) = j
           
           if (nsend.ne.0) then
              ic = 0
              do i = 1, ip -1
                 ic = ic +proc_send(i)
              enddo
              
              do i = 1, proc_send(ip)
                 node_send(ic +i) = i4vec(i)
              enddo
           endif
       endif

       deallocate(i4vec, cs_mpi)

     enddo
     
     if (nsend.eq.0) then
        do i = 1, nproc
           nsend = nsend +proc_send(i)
        enddo
     endif
     
     return
     
     end subroutine SETUP_MPI
