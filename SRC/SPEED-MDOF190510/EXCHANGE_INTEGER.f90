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

!> @brief Exchanges integers between MPI processes.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nsend number of values to be sent
!> @param[in] buff_send  buffer containing the values to be sent
!> @param[in] nrecv  number of values to be received
!> @param[in] buff_recv  buffer for the values to be received
!> @param[in] nproc  number of processors
!> @param[in] proc_send  vector containing the number of values to be sent to each proc.
!> @param[in] proc_recv  vector containing the number of values to be received from each proc.
!> @param[in] comm  MPI communicator
!> @param[in] status  MPI status
!> @param[in] ierr  MPI error tag
!> @param[in] myid MPI id
     
     subroutine EXCHANGE_INTEGER(nsend,buff_send,nrecv,buff_recv,&
                                     nproc,proc_send,proc_recv,&
                                     comm,status,ierr,myid)
    
     implicit none
     include 'SPEED.MPI'
     
     integer*4 :: nsend,nrecv,nproc,comm,ierr,tag,myid
     integer*4 :: ip,ir,is
     
     integer*4, dimension(nproc) :: proc_send,proc_recv
     integer*4, dimension(SPEED_STATUS_SIZE) :: status
     integer*4, dimension(nproc) :: request,requests

!     real*8, dimension(*) :: buff_send,buff_recv
     integer*4, dimension(nsend) :: buff_send
     integer*4, dimension(nrecv) :: buff_recv
    
     tag=2120
     
     
     
     is = 1
     do ip = 1, nproc
         if (proc_send(ip).gt.0) then
            call MPI_ISEND(buff_send(is:(is +proc_send(ip) -1)),&
                 proc_send(ip),SPEED_INTEGER,(ip -1),&
                tag,comm,requests(ip),ierr)
         endif
         is = is +proc_send(ip)
      enddo

      ir = 1
      do ip = 1, nproc
         if (proc_recv(ip).gt.0) then
            call MPI_IRECV(buff_recv(ir:(ir +proc_recv(ip) -1)),&
                 proc_recv(ip),SPEED_INTEGER,(ip -1),&
                 MPI_ANY_TAG,comm,request(ip),ierr)
         endif
         ir = ir +proc_recv(ip)
      enddo
    
     do ip = 1, nproc 
         if (proc_send(ip).gt.0) then 
            call MPI_WAIT(requests(ip),status,ierr)
         endif
     enddo 
     do ip = 1, nproc
        if (proc_recv(ip).gt.0) then
           call MPI_WAIT(request(ip),status,ierr)
        endif
     enddo
     
     speed_tag = speed_tag +1
     if (speed_tag.gt.speed_tag_max) speed_tag = speed_tag_min
     
     return
     
     end subroutine EXCHANGE_INTEGER

