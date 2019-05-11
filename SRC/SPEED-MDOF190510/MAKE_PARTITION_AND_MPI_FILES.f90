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

!> @brief Makes Partitioning and writes files *.mpi
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0

!> @note A new partition is made if elementdomain.mpi is absent in the
!! workig directory or in the folder given in SPEED.input (MPIFILE). 

      subroutine MAKE_PARTITION_AND_MPI_FILES()
      
      use max_var
      use speed_par
      

      implicit none
 
      include 'SPEED.MPI'      
            
!***************************************************************************************************************
!     Partitioning
                                       
      if (mpi_id.eq.0)  write(*,'(A)') '---------------------Partitioning----------------------'
      

      
      if (mpi_id.eq.0) then 
              if(len_trim(mpi_file) .ne. 70) then                                                                
                 file_part = mpi_file(1:len_trim(mpi_file)) // '/elemdomain.mpi'
              else 
                 file_part = 'elemdomain.mpi'     
              endif      
      
              inquire(file=file_part,exist=filefound)
      
!              write(*,*) filefound, file_part, mpi_id
!        read(*,*)

              if(mpi_np .gt. 1) then 
         
           if(filefound .eqv. .TRUE.) then
      
              if (mpi_id.eq.0)  write(*,'(A)') 'Reading existing partitioning...'
                       unit_part = 400                                 
                 open(unit_part,file=file_part)        
                 read(unit_part,*) trash
                 do i = 1, nelem
                       read(unit_part,*) trash, elem_domain(i)
                 enddo
                 close(unit_part)
                 if (mpi_id.eq.0)  write(*,'(A)') 'Read.'
           
              else
      
                 write(*,'(A)') 'Making new partitioning...'
                 call MESH_PARTITIONING(mpi_file,nelem,nnod_macro,mpi_np,con,0)      
!                call MPI_BARRIER(mpi_comm, mpi_ierr)
                 !sustitute this part adding output elem_domain in mesh partitioning
                 unit_part = 400                                 
                 open(unit_part, file = file_part)        
                 read(unit_part,*) trash
                 do i = 1, nelem
                    read(unit_part,*) trash, elem_domain(i)
                 enddo
                 close(unit_part)
                 write(*,'(A)') 'Made.'
      
              endif
         endif         

      endif!if (mpi_id .eq. 0)
      
      call MPI_BARRIER(mpi_comm, mpi_ierr)
      call MPI_BCAST(elem_domain,nelem,SPEED_INTEGER,0,MPI_COMM_WORLD,mpi_ierr)
      
!************************************************************************************************************** 
      
      allocate(elem_index(nelem))


      nelem_dom = 0
      do ie = 1,nelem
         elem_index(ie) = 0
         if (elem_domain(ie).eq.mpi_id) then
            nelem_dom = nelem_dom +1
            elem_index(ie) = nelem_dom
         endif
      enddo
      
 
      
!***************************************************************************************************************
!     SPECTRAL CONNECTIVITY --- ONLY FOR MASTER PROCESS ---
      
      if (mpi_id.eq.0) then
         write(*,'(A)') 
         write(*,'(A)')'--------------Making Spectral connectivity-------------'
      endif
      
      call MPI_BARRIER(mpi_comm, mpi_ierr)
            
      if(mpi_id .eq. 0) then 
        
        file_mpi = 'nodedomain.mpi';
        
        if(len_trim(mpi_file) .ne. 70) then                                                                                  
           file_mpi_new = mpi_file(1:len_trim(mpi_file)) // '/' // file_mpi
        else 
           file_mpi_new = file_mpi     
        endif
      
        inquire(file=file_mpi_new, exist=filefound) 

        if(filefound .eqv. .FALSE.) then
      
                allocate(node_weight(nnod_macro))
      
                 call MAKE_WGT_GRID_NODES(nnod_macro,nelem,con,node_weight,nnz_node_weight)
          
                allocate(node_pointer(0:nnz_node_weight))
      
                call MAKE_GRID_NODES(nnod_macro,nelem,con,node_weight,nnz_node_weight,node_pointer)
      
                deallocate(node_weight)
    
                con_nnz = nelem +1
                do ie = 1,nelem
                     do j = 1,nmat
                        if (tag_mat(j).eq.con(ie,1)) nn = sdeg_mat(j) +1
                     enddo
                     con_nnz = con_nnz + nn*nn*nn +1
                enddo
      
                allocate(con_spx(0:con_nnz))
      
                call MAKE_SPX_CON(nelem,con,nmat,tag_mat,sdeg_mat,&
                                      nnz_node_weight,node_pointer,con_nnz,con_spx,nnod)


                deallocate(node_pointer);  allocate(node_weight(nnod))
      
                call MAKE_WGT_SPX_NODES(nnod,con_nnz,con_spx,node_weight,nnz_node_weight)
      
                allocate(node_pointer(0:nnz_node_weight))

                call MAKE_SPX_NODES(nnod,con_nnz,con_spx,node_weight,nnz_node_weight,node_pointer)
      
                deallocate(node_weight)
  
          write(*,'(A)') 'Made.'
        
!     SPECTRAL CONNECTIVITY END --- ONLY FOR MASTER PROCESS ---
!***************************************************************************************************************
!***************************************************************************************************************
!     NODE DOMAIN FOR MPI

          allocate(node_domain(nnod))
      
          do ie = 1, nelem
             do i = con_spx(ie -1) +1, con_spx(ie) -1
                node_domain(con_spx(i)) = elem_domain(ie)
             enddo
          enddo

          file_mpi = 'nodedomain.mpi'; unit_mpi = 40
        
          if(len_trim(mpi_file) .ne. 70) then                                                                                  
             file_mpi_new = mpi_file(1:len_trim(mpi_file)) // '/' // file_mpi
          else 
             file_mpi_new = file_mpi     
          endif

          open(unit_mpi,file=file_mpi_new)        
          write(unit_mpi,*) nnod
          do i = 1, nnod
             write(unit_mpi,*) node_domain(i)
          enddo
          close(unit_mpi)
          
          deallocate(node_domain)                  

!     NODE DOMAIN FOR MPI
!***************************************************************************************************************
!***************************************************************************************************************
!     SPECTRAL LOCAL CONNECTIVITY  
          
          write(*,'(A)')                                
          write(*,'(A)') '---------------Making local connectivity---------------'
           
          do ip = 0, mpi_np-1 
         
              nelem_loc = 0 
              do in = 1,nelem
                 if(elem_domain(in) .eq. ip) nelem_loc = nelem_loc +1
                    enddo        
      
              con_nnz_loc = nelem_loc +1
              do ie = 1, nelem
                 do j = 1, nmat         
                    if ((tag_mat(j).eq. con_spx(con_spx(ie-1))) .and. (elem_domain(ie) .eq. ip)) then
                       nn = sdeg_mat(j) +1
                       con_nnz_loc = con_nnz_loc + nn*nn*nn +1
                    endif
                 enddo
              enddo

              allocate(con_spx_loc(0:con_nnz_loc))
    
              call MAKE_SPX_CON_LOC(nelem_loc, nelem, elem_domain, &
                        con_nnz,  con_spx, con_nnz_loc, con_spx_loc, &
                        nmat, tag_mat, sdeg_mat, ip)



         
              file_mpi = 'cons000000.mpi'
              unit_mpi = 40                                 
              if (ip.lt. 10) then                                        
                  write(file_mpi(10:10),'(i1)') ip                
              else if (ip .lt. 100) then                                
                  write(file_mpi(9:10),'(i2)') ip                
              else if (ip .lt. 1000) then                                
                  write(file_mpi(8:10),'(i3)') ip                
              else if (ip .lt. 10000) then                                
                  write(file_mpi(7:10),'(i4)') ip                
              else if (ip .lt. 100000) then        
                  write(file_mpi(6:10),'(i5)') ip                
              else if (ip .lt. 1000000) then                                
                  write(file_mpi(5:10),'(i6)') ip                
              endif
                
              if(len_trim(mpi_file) .ne. 70) then                                         
                 file_mpi_new = mpi_file(1:len_trim(mpi_file)) // '/' // file_mpi
              else 
                 file_mpi_new = file_mpi     
              endif
             
                open(unit_mpi,file=file_mpi_new)
              write(unit_mpi,*) con_nnz_loc                
              do i = 0, con_nnz_loc
                       write(unit_mpi,*) con_spx_loc(i)
              enddo
              close(unit_mpi)        

              ic = 0
              do in = 1,nelem
                 if (elem_domain(in) .eq. ip) ic = ic +1
              enddo
              nelem_loc = ic
     
              allocate(local_el_num(nelem_loc))
      
              ic = 0
              do in = 1,nelem
                 if (elem_domain(in) .eq. ip) then
                    ic = ic +1
                    local_el_num(ic) = in
                 endif
              enddo

                   con_nnz_bc_loc = 0
                   nface_loc = 0
    
              ic = 0 
              if (nface.gt.0) then

                do i = 1,nface

                  call GET_ELEM_FROM_FACE(nnz_node_weight,node_pointer,con_bc(i,2),con_bc(i,3),con_bc(i,5),ie)
                  call GET_INDLOC_FROM_INDGLO(local_el_num, nelem_loc, ie, ic)

                  if (ic .ne. 0) then !I found global index in local_el_num

                    do j = 1,nmat
            
                      if (tag_mat(j) .eq. con_spx_loc(con_spx_loc(ic-1))) then
                          nn = sdeg_mat(j) +1  
                          con_nnz_bc_loc = con_nnz_bc_loc +nn*nn +1
                          nface_loc = nface_loc + 1                                             
                      endif
                   
                    enddo
                  endif    
                
                enddo
         
                con_nnz_bc_loc = con_nnz_bc_loc + nface_loc + 1 

                allocate(con_spx_bc_loc(0:con_nnz_bc_loc))
                con_spx_bc_loc(0) = nface_loc + 1
         
                call MAKE_SPX_CON_LOC_BOUND(con_nnz_loc,con_spx_loc, nface,con_bc, &
                                         nmat,tag_mat,sdeg_mat,nnz_node_weight,node_pointer, &
                                         con_nnz_bc_loc,con_spx_bc_loc, local_el_num, nelem_loc, &
                                         nface_loc)

              endif


              file_mpi = 'conb000000.mpi'
              unit_mpi = 40                                 
              if (ip.lt. 10) then                                        
                  write(file_mpi(10:10),'(i1)') ip                
              else if (ip .lt. 100) then                                
                  write(file_mpi(9:10),'(i2)') ip                
              else if (ip .lt. 1000) then                                
                  write(file_mpi(8:10),'(i3)') ip                
              else if (ip .lt. 10000) then                                
                  write(file_mpi(7:10),'(i4)') ip                
              else if (ip .lt. 100000) then        
                  write(file_mpi(6:10),'(i5)') ip                
              else if (ip .lt. 1000000) then                                
                  write(file_mpi(5:10),'(i6)') ip                
              endif

              if(len_trim(mpi_file) .ne. 70) then             
                 file_mpi_new = mpi_file(1:len_trim(mpi_file)) // '/' // file_mpi
              else 
                 file_mpi_new = file_mpi     
              endif

              open(unit_mpi,file=file_mpi_new)                        
              write(unit_mpi,*) con_nnz_bc_loc                
              do i = 0, con_nnz_bc_loc
                       write(unit_mpi,*) con_spx_bc_loc(i)
              enddo
              close(unit_mpi)        


           deallocate(con_spx_loc, con_spx_bc_loc, local_el_num)

          enddo

          deallocate(con_spx,node_pointer)
        
        endif

      endif 

      if(mpi_id.eq.0) write(*,'(A)') 'Made.'
      call MPI_BARRIER(mpi_comm, mpi_ierr)
      
!***************************************************************************************************************
!     LOCAL NUMERATION FOR ELEMENTS  
   
      ic = 0
      do in = 1,nelem
       
         if (elem_domain(in) .eq. mpi_id) then
            ic = ic +1
         endif
      enddo
      nelem_loc = ic
     
      allocate(local_el_num(nelem_loc))
      
      ic = 0
      do in = 1,nelem
       
         if (elem_domain(in) .eq. mpi_id) then
            ic = ic +1
            local_el_num(ic) = in
         endif
      enddo

      call MPI_BARRIER(mpi_comm, mpi_ierr)
      
!     LOCAL NUMERATION FOR ELEMENTS  END
!***************************************************************************************************************

       file_mpi = 'cons000000.mpi';   unit_mpi = 40                                 

       if (mpi_id .lt. 10) then                                        
          write(file_mpi(10:10),'(i1)') mpi_id                
       else if (mpi_id .lt. 100) then                                
          write(file_mpi(9:10),'(i2)') mpi_id                
       else if (mpi_id .lt. 1000) then                                
          write(file_mpi(8:10),'(i3)') mpi_id                
       else if (mpi_id .lt. 10000) then                                
          write(file_mpi(7:10),'(i4)') mpi_id                
       else if (mpi_id .lt. 100000) then        
          write(file_mpi(6:10),'(i5)') mpi_id                
       else if (mpi_id .lt. 1000000) then                                
          write(file_mpi(5:10),'(i6)') mpi_id                
       endif

       if(len_trim(mpi_file) .ne. 70) then                                                                                  
          file_mpi_new = mpi_file(1:len_trim(mpi_file)) // '/' // file_mpi
       else 
          file_mpi_new = file_mpi     
       endif


       open(unit_mpi,file=file_mpi_new)        
       read(unit_mpi,*) con_nnz_loc
       allocate(con_spx_loc(0:con_nnz_loc))

                      
       do i = 0, con_nnz_loc
            read(unit_mpi,*) con_spx_loc(i)
       enddo
       close(unit_mpi)        


       file_mpi = 'conb000000.mpi';  unit_mpi = 40                                 

       if (mpi_id .lt. 10) then                                        
          write(file_mpi(10:10),'(i1)') mpi_id                
       else if (mpi_id .lt. 100) then                                
          write(file_mpi(9:10),'(i2)') mpi_id                
       else if (mpi_id .lt. 1000) then                                
          write(file_mpi(8:10),'(i3)') mpi_id                
       else if (mpi_id .lt. 10000) then                                
          write(file_mpi(7:10),'(i4)') mpi_id                
       else if (mpi_id .lt. 100000) then        
          write(file_mpi(6:10),'(i5)') mpi_id                
       else if (mpi_id .lt. 1000000) then                                
          write(file_mpi(5:10),'(i6)') mpi_id                
       endif

       if(len_trim(mpi_file) .ne. 70) then                                                                                  
          file_mpi_new = mpi_file(1:len_trim(mpi_file)) // '/' // file_mpi
       else 
          file_mpi_new = file_mpi     
       endif
    
       open(unit_mpi,file=file_mpi_new)        
       read(unit_mpi,*) con_nnz_bc_loc
       allocate(con_spx_bc_loc(0:con_nnz_bc_loc))

                      
      do i = 0, con_nnz_bc_loc
            read(unit_mpi,*) con_spx_bc_loc(i)
      enddo
      close(unit_mpi)        


!     SPECTRAL LOCAL CONNECTIVITY  
!***************************************************************************************************************
!***************************************************************************************************************
!     COUNTING LOCAL NODES & LOCAL NODE NUMBERING
      if(mpi_id .eq.0) write(*,'(A)')
      if(mpi_id .eq.0) write(*,'(A)') '------------------Counting local nodes-----------------'

      call GET_LOC_NODE_NUM(con_nnz_loc, con_spx_loc, nnod_loc, mpi_id)

      if(mpi_id .eq.0) write(*,'(A)') 'Made.'

      allocate(local_node_num(nnod_loc))
      local_node_num = 0
      
      
      if(mpi_id .eq.0) write(*,'(A)')
      if(mpi_id .eq.0) write(*,'(A)') '--------------Make local nodes numbering---------------'
      
      
      call MAKE_REN_LOC_NODE(con_nnz_loc, con_spx_loc, nnod_loc, mpi_id, &
                             local_node_num)
      
!      ic = 0
!!      do im = 1,nmat
!!         nn = sdeg_mat(im) +1
!         do ie = 1,nelem_loc
!            im = con_spx_loc(con_spx_loc(ie -1) +0); nn = sdeg_mat(im) +1
!
!!            if (con_spx_loc(con_spx_loc(ie -1) +0) .eq. tag_mat(im)) then
!                do i = 1, nn*nn*nn 
!                   if (find(nnod_loc, local_node_num, con_spx_loc(con_spx_loc(ie -1) + i)) .ne. 1) then
!                      ic = ic + 1
!                      local_node_num(ic) = con_spx_loc(con_spx_loc(ie -1) + i)
!                   endif
!                 enddo
!!             endif         
!          enddo
!!       enddo         
               
                   
     if(mpi_id .eq.0) write(*,'(A)') 'Made.'
     if(mpi_id .eq.0) write(*,'(A)')

     file_mpi = 'nodedomain.mpi'; unit_mpi = 40                                 
     if(len_trim(mpi_file) .ne. 70) then                                                                                  
        file_mpi_new = mpi_file(1:len_trim(mpi_file)) // '/' // file_mpi
     else 
        file_mpi_new = file_mpi     
     endif

     if(mpi_id .eq.0) then
        !start1=MPI_WTIME() 
       open(unit_mpi,file=file_mpi_new)        
       read(unit_mpi,*) nnod
     endif
     call MPI_BCAST(nnod,1,SPEED_INTEGER,0,MPI_COMM_WORLD,mpi_ierr)

     allocate(node_domain(nnod))
     if(mpi_id .eq.0) then
       do i = 1, nnod
            read(unit_mpi,*) node_domain(i)
       enddo
       close(unit_mpi)           
       !write (*,*) "File nodedomain.mpi read, broadcast to slaves" 
     endif
     call MPI_BCAST(node_domain,nnod,SPEED_INTEGER,0,MPI_COMM_WORLD,mpi_ierr)

       !write (*,*) "Vector received by process",mpi_id

       
      ic = 0
      do in = 1, nnod
         if (node_domain(in).eq. mpi_id) ic = ic +1
      enddo   
      nnode_dom = ic
           
  !    write(*,'(A,I3,I8)')'Nodes active on proc : ', mpi_id, ic      
      write(*,'(A,I6,I15)') 'Total nodes on Proc. :  ', mpi_id, nnode_dom
      
 
      call MPI_BARRIER(mpi_comm, mpi_ierr)

      !if(mpi_id .eq.0) then
         !start2=MPI_WTIME()
         !write (*,*) "Time to send nodedomain.mpi to all &
         !& proc",start2-start1,"sec"
      !endif

!     COUNTING LOCAL NODES & LOCAL NODE NUMBERING
!***********************************************************************************

      end subroutine MAKE_PARTITION_AND_MPI_FILES
