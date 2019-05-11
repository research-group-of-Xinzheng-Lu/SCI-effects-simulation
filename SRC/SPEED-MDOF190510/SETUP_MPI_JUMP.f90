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
!> @note If run with nsend_jump and 
!!     nrecv_jump equal to zero gives back the correct values of these to allocate node_send_jump 
!!     and node_recv_jump, otherwise fills node_send_jump and node_recv_jump with the indices of the
!!     nodes to be sent to and received from the other processors, and proc_send_jump and 
!!     proc_recv_jump with the number of nodes to be sent and received
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] ne_dg_loc  number of local dg elements
!> @param[in] el_new  structure containing dg elements
!> @param[in] nnode  number of local nodes
!> @param[in] node_proc  vector containing the processor to which the node is assigned
!> @param[in] nel_loc   number of local elements
!> @param[in] loc_el_num   local element numeration
!> @param[in] nn_loc  number of local nodes
!> @param[in] loc_n_num  local node numeration
!> @param[in] ncs  length of the connectivity vector
!> @param[in] cs  local connectivity vector
!> @param[in] nproc  number of processors
!> @param[in] id  processor identifyer
!> @param[in] mpi_file folder where *.mpi files are stores
!> @param[out] nsend_jump  number of nodes to be sent to other processors
!> @param[out] node_send_jump  ordered vector of nodes to be sent
!> @param[out] nrecv_jump  number of nodes to be received from other procs
!> @param[out] node_recv_jump  ordered vector of nodes to be received
!> @param[out] proc_send_jump  vector containing the number of the nodes to be sent to each proc.
!> @param[out] proc_recv_jump  vector containing the number of the nodes to be received from each proc.

   
     subroutine SETUP_MPI_JUMP(ne_dg_loc, el_new, nnode, node_proc, nel_loc, loc_el_num, &
                           nnloc, loc_n_num, ncs,cs,&
                           nsend_jump, node_send_jump, nrecv_jump, node_recv_jump,&
                           nproc, proc_send_jump, proc_recv_jump, id, mpi_file)
     
     use max_var
     use DGJUMP
     
     implicit none

     type(el4loop), dimension(ne_dg_loc), intent(in) :: el_new
     
     character*70 :: filename, mpi_file, filename_new

     integer*4 :: nel_loc, ncs, nsend_jump, nrecv_jump, nproc, id, ne_dg_loc, nnode, nnloc
     integer*4 :: i,j,k,ie,ip,ic,ieloc
     integer*4 :: it,ir,im, tt, ttt
     integer*4 :: unitfile, ncs_mpi, nelem_mpi, pos     

     integer*4, dimension(:), allocatable :: i4vec
     integer*4, dimension(:), allocatable :: cs_mpi
     integer*4, dimension(nel_loc) :: loc_el_num
     integer*4, dimension(nnloc) :: loc_n_num(nnloc)
     integer*4, dimension(nnode) :: node_proc
     integer*4, dimension(0:ncs) :: cs
     integer*4, dimension(nsend_jump) :: node_send_jump
     integer*4, dimension(nrecv_jump) :: node_recv_jump
     integer*4, dimension(nproc) :: proc_send_jump, proc_recv_jump
     
     do ip = 1, nproc
     
        if ((ip -1).eq.id) then
           proc_recv_jump(ip) = 0
        else
           
           filename = 'DGCS000000.mpi'
           unitfile = 40                                 
           if (ip-1 .lt. 10) then                                        
               write(filename(10:10),'(i1)') ip-1                
           else if (ip-1 .lt. 100) then                                
               write(filename(9:10),'(i2)') ip-1
           else if (ip-1 .lt. 1000) then                                
               write(filename(8:10),'(i3)') ip-1
           else if (ip-1 .lt. 10000) then                                
               write(filename(7:10),'(i4)') ip-1
           else if (ip-1 .lt. 100000) then                                
               write(filename(6:10),'(i5)') ip-1
           else if (ip-1 .lt. 1000000) then                                
               write(filename(5:10),'(i6)') ip-1
           endif
           
           if(len_trim(mpi_file) .ne. 70) then                                                                                  
              filename_new = mpi_file(1:len_trim(mpi_file)) // '/' // filename
           else
              filename_new = filename
           endif
          
           open(unitfile,file=filename_new)        
           read(unitfile,*) ncs_mpi
           allocate(cs_mpi(0:ncs_mpi))
                      
           do i = 0, ncs_mpi
               read(unitfile,*) cs_mpi(i)
           enddo
           close(unitfile)                
       

           allocate(i4vec(ncs_mpi))
           nelem_mpi = cs_mpi(0)-1

           i4vec = 0           
           ic = 0
           
           
           do im = 1, ne_dg_loc
           
              do it = 1, el_new(im)%num_of_ne
    
                 ie = el_new(im)%el_conf(it,1)
                 
                 call CHECK_MPI(ncs_mpi, cs_mpi, ie, tt, pos) 
                 
                 if(tt .eq. 1) then
                 
                    do i = cs_mpi(pos -1) +1, cs_mpi(pos) -1
                     
                      ! if (node_proc(cs_mpi(i)) .eq. (ip -1)) then
                           ttt = 0
                           if(ic .ge. 1) call CHECK_VECTOR(ic, i4vec(1:ic), cs_mpi(i), ttt)
                           if(ttt.eq.0) then
                             ic = ic +1
                             i4vec(ic) = cs_mpi(i)
                           endif  
                      ! endif
                    enddo
                 
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
           proc_recv_jump(ip) = j
           
           if (nrecv_jump .ne. 0) then
              ic = 0
              do i = 1, ip -1
                 ic = ic + proc_recv_jump(i)
              enddo
              
              do i = 1, proc_recv_jump(ip)
                 node_recv_jump(ic +i) = i4vec(i)
              enddo
           endif
         
           deallocate(cs_mpi, i4vec) 
        endif
             
     enddo
     
     if (nrecv_jump.eq.0) then
        do i = 1, nproc
           nrecv_jump = nrecv_jump + proc_recv_jump(i)
        enddo
     endif
     
!-----------------------nodes to send------------------------------------
     
     allocate(i4vec(ncs))
     
     do ip = 1, nproc
        if ((ip -1).eq.id) then
           proc_send_jump(ip) = 0
        else

           filename = 'DGCS000000.mpi'
           unitfile = 40                                 
           if (ip-1 .lt. 10) then                                        
               write(filename(10:10),'(i1)') ip-1                
           else if (ip-1 .lt. 100) then                                
               write(filename(9:10),'(i2)') ip-1
           else if (ip-1 .lt. 1000) then                                
               write(filename(8:10),'(i3)') ip-1
           else if (ip-1 .lt. 10000) then                                
               write(filename(7:10),'(i4)') ip-1
           else if (ip-1 .lt. 100000) then                                
               write(filename(6:10),'(i5)') ip-1
           else if (ip-1 .lt. 1000000) then                                
               write(filename(5:10),'(i6)') ip-1
           endif                                                              

           if(len_trim(mpi_file) .ne. 70) then                                                                                  
             filename_new = mpi_file(1:len_trim(mpi_file)) // '/' // filename
           else
             filename_new = filename
           endif
          
           open(unitfile,file=filename_new)        
           read(unitfile,*) ncs_mpi
           allocate(cs_mpi(0:ncs_mpi))
                      
           do i = 0, ncs_mpi
               read(unitfile,*) cs_mpi(i)
           enddo
           close(unitfile)                


           i4vec = 0
           ic = 0
           
           do im = 1, ne_dg_loc
              ir = el_new(im)%ind
              
              do it = 1, el_new(im)%num_of_ne    
                 ie = el_new(im)%el_conf(it,1)
                 
                 call CHECK_MPI(ncs_mpi, cs_mpi, ie, tt, pos) 
                 if(tt .eq. 1) then
                 
                    call GET_INDLOC_FROM_INDGLO(loc_el_num, nel_loc, ir, ieloc)
                              
                    do i = cs(ieloc -1) +1, cs(ieloc) -1

                      ! if (node_proc(loc_n_num(cs(i))) .eq. id) then
                           ttt = 0
                           if(ic .ge. 1)  call CHECK_VECTOR(ic, i4vec(1:ic), loc_n_num(cs(i)), ttt)
                           
                           if(ttt.eq.0) then
                                ic = ic +1
                                i4vec(ic) = loc_n_num(cs(i))
                           endif    
                      ! endif
                    enddo
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
                      
           
           if (ic.eq.0) j = 0
           proc_send_jump(ip) = j
           
           if (nsend_jump.ne.0) then
              ic = 0
              do i = 1, ip -1
                 ic = ic +proc_send_jump(i)
              enddo
              
              do i = 1, proc_send_jump(ip)
                 node_send_jump(ic +i) = i4vec(i)
              enddo
           endif
           
           deallocate(cs_mpi)
           
        endif
     enddo
     
     if (nsend_jump.eq.0) then
        do i = 1, nproc
           nsend_jump = nsend_jump +proc_send_jump(i)
        enddo
     endif
     

     deallocate(i4vec)
     
     return
     
     end subroutine SETUP_MPI_JUMP
