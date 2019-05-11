!> @brief Exchanges acc between MPI processes. by Y TIAN
     
subroutine EXCHANGE_ACC(nsend,buff_send,nrecv,buff_recv,&
                                nproc,proc_send,proc_recv,&
                                comm,status,ierr,myid)
     
implicit none
include 'SPEED.MPI'
     
integer*4 :: nsend,nrecv,nproc,comm,ierr,tag,myid  
integer*4 :: ip,ir,is
     
integer*4 :: proc_send,proc_recv
integer*4, dimension(SPEED_STATUS_SIZE) :: status
integer*4, dimension(nproc) :: request,requests
real*8, dimension(nsend) :: buff_send
real*8, dimension(nproc*nrecv) :: buff_recv
    
tag=2120
     
do ip=1,nproc-1
if(myid.eq.ip) then
    call MPI_ISEND(buff_send(1:proc_send),&
            proc_send,SPEED_DOUBLE,0,&
        tag,comm,requests(ip),ierr)
    call MPI_WAIT(requests(ip),status,ierr)
elseif(myid.eq.0) then
    call MPI_IRECV(buff_recv((1+proc_recv*(ip-1)):proc_recv*ip),&
            proc_recv,SPEED_DOUBLE,ip,&
            MPI_ANY_TAG,comm,request(ip),ierr)
    call MPI_WAIT(request(ip),status,ierr)
end if
end do
    


     
return
     
     end subroutine EXCHANGE_ACC

