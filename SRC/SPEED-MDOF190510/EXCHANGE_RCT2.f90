!> @brief Exchanges acc between MPI processes. by Y TIAN
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
     
subroutine EXCHANGE_RCT2(nsend,buff_send,nrecv,buff_recv,&
                                nproc,proc_send,proc_recv,&
                                comm,status,ierr,myid,MDOFftran)
     
implicit none
include 'SPEED.MPI'
     
integer*4 :: nsend,nrecv,nproc,comm,ierr,tag,myid  
integer*4 :: ip,ir,is
integer*8,dimension(2*nproc):: MDOFftran
integer*4 :: proc_send,proc_recv
integer*4, dimension(SPEED_STATUS_SIZE) :: status
integer*4, dimension(nproc) :: request,requests
real*8, dimension(nsend) :: buff_send
real*8, dimension(nrecv) :: buff_recv
    
tag=2120
     
do ip=1,nproc-1
if(myid.eq.0) then
if(MDOFftran(nproc+ip+1).gt.0) then
    call MPI_ISEND(buff_send(1:proc_send),&
            proc_send,SPEED_DOUBLE,ip,&
        tag,comm,requests(ip),ierr)
    call MPI_WAIT(requests(ip),status,ierr)
end if
elseif(myid.eq.ip) then
    call MPI_IRECV(buff_recv(1:proc_recv),&
            proc_recv,SPEED_DOUBLE,0,&
            MPI_ANY_TAG,comm,request(ip),ierr)
    call MPI_WAIT(request(ip),status,ierr)
end if
end do

return
     
end subroutine EXCHANGE_RCT2

