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


!> @brief Setup for DG interface conditions 
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
      subroutine MAKE_DG_INTERFACE_CONDITIONS()
        
              
      use max_var
      use str_mesh 
      use str_mesh_scratch                   
      use DGJUMP

      use speed_par
      use speed_par_dg

      implicit none
 
      include 'SPEED.MPI' 
      

      
      
      nelem_dg = 0     
      nnode_dg = 0
      allocate(i4count(nnod_loc))
      i4count = 0

    
      call GET_NODE_FROM_FACE(nnod_loc, con_nnz_bc_loc, con_spx_bc_loc, nload_dg_el,tag_dg_el,&
                          nnode_dg, i4count, local_node_num)


      ! output: - number of dg local elements 
      !         - number of dg global elements 

      call GET_DIME_DG(nmat, sdeg_mat, tag_mat, con_nnz_loc, con_spx_loc, &
                     nnod_loc, local_node_num, nelem_loc, local_el_num, &
                     xx_spx_loc, yy_spx_loc, zz_spx_loc, &
                     nelem_dg, nelem_dg_glo, i4count, mpi_id, mpi_comm)

      if (mpi_id.eq.0 .and. nelem_dg_glo .gt. 0) write(*,'(A)') '-------------Setup DG interface conditions-------------'  

      allocate(faces(3,nelem_dg), area_nodes(25,nelem_dg))  
      faces = 0; area_nodes = 0.d0
      
      if(nelem_dg_glo .gt. 0) then

              file_face = 'FACES.input'

              call SETUP_DG(nmat, sdeg_mat, tag_mat, con_nnz_loc, con_spx_loc, &
                     nnod_loc, local_node_num, nelem_loc, local_el_num, &
                     xx_spx_loc, yy_spx_loc, zz_spx_loc, &
                     nelem_dg, nelem_dg_glo, i4count, mpi_id, mpi_comm,&
                     alfa11, alfa12, alfa13, &
                     alfa21, alfa22, alfa23, &
                     alfa31, alfa32, alfa33, &
                     beta11, beta12, beta13, &
                     beta21, beta22, beta23, &
                     beta31, beta32, beta33, &
                     gamma1, gamma2, gamma3, &
                     delta1, delta2, delta3, &                     
                     faces, area_nodes)


              write(*,'(A,I10,A4,I10)') 'DG Faces on Proc.', mpi_id, ' are', nelem_dg
              allocate(count_faces(mpi_np))
              call MPI_BARRIER(mpi_comm, mpi_ierr)                  
              call MPI_ALLGATHER(nelem_dg, 1, SPEED_INTEGER, count_faces, 1, SPEED_INTEGER, mpi_comm, mpi_ierr)
      
      

              inquire(file=file_face, exist=filefound)
              if(filefound .eqv. .FALSE.) then
                 
                 if (mpi_id .eq. 0) write(*,'(A)') 'Writing FACES.input'
                 call WRITE_FILE_FACES(mpi_file, file_face, faces, area_nodes, nelem_dg, &
                                             nelem_dg_glo, mpi_id, mpi_comm, count_faces, mpi_np)
              endif    
      
              deallocate(faces, area_nodes, count_faces)  
              allocate(faces(3,nelem_dg_glo), area_nodes(25,nelem_dg_glo))
              faces = 0 
              area_nodes = 0.d0
                  
              file_face = 'FACES.input'         
              call READ_FACES(file_face, faces, area_nodes, nelem_dg_glo, mpi_id, mpi_comm)
                   
              head_file = 'DGFS.input'
              inquire(file=head_file, exist=filefound)
              if(filefound .eqv. .FALSE.) then     
              
                  if (mpi_id .eq. 0) write(*,'(A)') 'Writing DGFS.input'
                   
                  ! output: - write DG_INTERFACE.out file
        
                  allocate(dg_els(nelem_dg))       
                  allocate(scratch_dg_els(nelem_dg)) 
        
                  call SETUP_DG_ELEM(nmat, sdeg_mat, tag_mat, con_nnz_loc, con_spx_loc, &
                                   nnod_loc, local_node_num, nelem_loc, local_el_num, &
                                   xx_spx_loc,yy_spx_loc,zz_spx_loc,&
                                   nelem_dg, nelem_dg_glo, &
                                   i4count, mpi_id, mpi_comm, mpi_np, &
                                   alfa11, alfa12, alfa13, &
                                   alfa21, alfa22, alfa23, &
                                   alfa31, alfa32, alfa33, &
                                   beta11, beta12, beta13, &
                                   beta21, beta22, beta23, &
                                   beta31, beta32, beta33, &
                                   gamma1, gamma2, gamma3, &
                                   delta1, delta2, delta3, &
                                   dg_els, scratch_dg_els, &
                                   tag_dg_el, tag_dg_yn, nload_dg_el, &
                                   con_bc, nface, mpi_file)

        
                   call WRITE_FILE_DGFS(nmat, sdeg_mat, tag_mat, con_nnz_loc, con_spx_loc, &
                                   nnod_loc, local_node_num, nelem_loc, local_el_num, &
                                   xx_spx_loc,yy_spx_loc,zz_spx_loc,&
                                   nelem_dg, nelem_dg_glo, &
                                  mpi_id, mpi_comm, mpi_np, &
                                   alfa11, alfa12, alfa13, &
                                   alfa21, alfa22, alfa23, &
                                   alfa31, alfa32, alfa33, &
                                   beta11, beta12, beta13, &
                                   beta21, beta22, beta23, &
                                   beta31, beta32, beta33, &
                                   gamma1, gamma2, gamma3, &
                                   delta1, delta2, delta3, &
                                   faces, area_nodes, dg_els, scratch_dg_els, &
                                   head_file, mpi_file)
        
                    deallocate(dg_els, scratch_dg_els) 
        
        
               endif
        
           endif

           
          if(nelem_dg_glo .gt. 0) then
 

           allocate(el_new(nelem_dg))  
           
                                            
           if (mpi_id .eq. 0) write(*,'(A)') 
           if (mpi_id .eq. 0) write(*,'(A)') '-----------------Making DG interfaces------------------' 

           

           call MAKE_DG_INTERFACE(nmat, sdeg_mat, tag_mat, prop_mat, con_nnz_loc, con_spx_loc, &
                           nnod_loc, local_node_num, nelem_loc, local_el_num, &
                           xx_spx_loc,yy_spx_loc,zz_spx_loc,&
                           nelem_dg, nelem_dg_glo, &
                           i4count, mpi_id, mpi_comm, mpi_np,&
                           alfa11, alfa12, alfa13, &
                           alfa21, alfa22, alfa23, &
                           alfa31, alfa32, alfa33, &
                           beta11, beta12, beta13, &
                           beta21, beta22, beta23, &
                           beta31, beta32, beta33, &
                           gamma1, gamma2, gamma3, &
                           delta1, delta2, delta3, dg_c, pen_c, & 
                           faces, area_nodes, el_new, head_file,testmode)
                                      
               ! 0 - file DGCSXXXXXX.mpi contains spectral connectivity for interface elements

          write(*,'(A,I6,A)') 'Proc:', mpi_id, ' Made.'
           
          call MPI_BARRIER(mpi_comm, mpi_ierr) 
          if (mpi_id .eq. 0) write(*,'(A)') 'Made All.' 
          
!***********************************************************************************
!     START SETUP 4 MPI COMUNICATION (ONLY BETWEEN NON-CONFORMING REGIONS)
!***********************************************************************************

          if (mpi_id.eq.0) write(*,'(A)')
          if (mpi_id.eq.0) write(*,'(A)')'----------Setting mpi parameters for DG jumps----------'

                 call SETUP_DG4MPI(el_new, nelem_dg, nelem_dg_glo, &
                               mpi_id, mpi_comm, mpi_np, con_nnz_loc, con_spx_loc, &
                               nelem_loc, local_el_num, &
                               local_node_num, nnod_loc, con_nnz_dg, total_els, mpi_file)

          file_mpi = 'nodedomain.mpi';   unit_mpi = 40                                 

          if(len_trim(mpi_file) .ne. 70) then                                                                                  
            file_mpi_new = mpi_file(1:len_trim(mpi_file)) // '/' // file_mpi
          else
            file_mpi_new = file_mpi
          endif

          open(unit_mpi,file=file_mpi_new)        
          read(unit_mpi,*) nnod

          allocate(node_domain(nnod))
          do i = 1, nnod
            read(unit_mpi,*) node_domain(i)
          enddo
          close(unit_mpi)            

          nsend_jump = 0
          nrecv_jump = 0
     
          allocate(proc_send_jump(mpi_np))
          allocate(proc_recv_jump(mpi_np)) 


          call SETUP_MPI_JUMP(nelem_dg, el_new, nnod, node_domain, nelem_loc, local_el_num, &
                           nnod_loc, local_node_num, con_nnz_loc, con_spx_loc,&
                           nsend_jump, node_send_jump, nrecv_jump, node_recv_jump,&
                           mpi_np, proc_send_jump, proc_recv_jump, mpi_id, mpi_file)

       
         allocate(node_send_jump(nsend_jump))
         allocate(node_recv_jump(nrecv_jump))


         call SETUP_MPI_JUMP(nelem_dg, el_new, nnod, node_domain, nelem_loc, local_el_num, &
                           nnod_loc, local_node_num, con_nnz_loc, con_spx_loc,&
                           nsend_jump, node_send_jump, nrecv_jump, node_recv_jump,&
                           mpi_np, proc_send_jump, proc_recv_jump, mpi_id, mpi_file)



         deallocate(node_domain)

         allocate(con_spx_dg(0:con_nnz_dg))
         con_spx_dg = 0
         call MAKE_SPX_CON_DG(con_spx_dg, con_nnz_dg, total_els, mpi_np, mpi_file)
       
         torf = 0    
         if(torf .eq. 1) then
             write(*,'(A)') 'ERROR: Repetitions on vector node_recv_jump' 
         else
            allocate(local_node_num_dg(nrecv_jump))
            do i = 1, nrecv_jump
               local_node_num_dg(i) = node_recv_jump(i)
               node_recv_jump(i) = i
            enddo    

             call MAKE_DG_LOC_NUMERATION(nrecv_jump, local_node_num_dg, &
                            con_nnz_dg, con_spx_dg, nelem_loc, local_el_num, mpi_id)            
            
         endif
     
      endif



      deallocate(faces, area_nodes)


      deallocate(i4count)
      if (mpi_id.eq.0 .and. nelem_dg_glo .gt. 0) write(*,'(A)') 'Made.'                  
      

      end subroutine MAKE_DG_INTERFACE_CONDITIONS

