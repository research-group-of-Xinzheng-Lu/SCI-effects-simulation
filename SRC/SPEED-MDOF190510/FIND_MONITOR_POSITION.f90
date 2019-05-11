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

!> @brief Find monitor position an writes MLST.input or MPGM. input files.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0


      subroutine FIND_MONITOR_POSITION()
      
      use max_var
      use speed_par
      

      implicit none
 
      include 'SPEED.MPI' 
      
      if (num_pgm.eq.1) then
         if (nelem_loc .gt. 0) then
             allocate (highest_pgm_loc(nelem_loc))
         endif
         call GET_HIGHEST_NODE(nnod_loc, nelem_loc, zz_spx_loc,&
                               local_node_num, con_nnz_loc, con_spx_loc, &
                               nmat, tag_mat, sdeg_mat, highest_pgm_loc)
      endif

      if (num_lst.eq.1) then
         if (nelem_loc .gt. 0) then
            allocate (highest_lst_loc(nelem_loc))
        endif
        call GET_HIGHEST_NODE(nnod_loc, nelem_loc, zz_spx_loc,&
                             local_node_num, con_nnz_loc, con_spx_loc, &
                             nmat, tag_mat, sdeg_mat, highest_lst_loc)
      endif


      if (mpi_id.eq.0) then
         write(*,'(A)') 
         write(*,'(A)') '--------------------Monitored points------------------'
      endif

!****************************************************************************************************
!                             MONITORED NODES - Peak Ground Map	-  
!****************************************************************************************************

      nmonitors_pgm = 0
      if (num_pgm.eq.1) then
                 !ATTENTION THIS PART SHUOLD BE MODIFIED ACCORDING TO THE CORRESPONDING LS.input PART
                 !IN ORDER TO AVOID MEMORY OVERFLOW  
		 file_PG = 'PG.input'
		 call READ_DIME_FILEPG(file_PG,nmonitors_pgm)     

		 allocate(n_monitor_pgm(nmonitors_pgm),el_monitor_pgm(nmonitors_pgm),dist_monitor_pgm(nmonitors_pgm))
		 allocate(x_monitor_pgm(nmonitors_pgm),y_monitor_pgm(nmonitors_pgm),z_monitor_pgm(nmonitors_pgm))
		 allocate(xr_monitor_pgm(nmonitors_pgm),yr_monitor_pgm(nmonitors_pgm),zr_monitor_pgm(nmonitors_pgm))



		 call READ_FILEPG(file_PG,nmonitors_pgm,x_monitor_pgm,y_monitor_pgm,z_monitor_pgm)		

	 	 if (file_mon_pgm.eq.0) then ! NO input file with the position of LST monitors

			do i = 1,nmonitors_pgm
		           call GET_NEAREST_NODE_PGM(nnod_loc, xx_spx_loc, yy_spx_loc, zz_spx_loc,&
					          x_monitor_pgm(i), y_monitor_pgm(i), z_monitor_pgm(i),&
						  n_monitor_pgm(i), dist_monitor_pgm(i), depth_search_mon_pgm)	
              
           
			   call GET_PNT_POS_PGM(nelem_loc,&
					 alfa11,alfa12,alfa13,alfa21,alfa22,alfa23,&
					 alfa31,alfa32,alfa33,beta11,beta12,beta13,&
	 				 beta21,beta22,beta23,beta31,beta32,beta33,&
					 gamma1,gamma2,gamma3,delta1,delta2,delta3,&
					 x_monitor_pgm(i),y_monitor_pgm(i),z_monitor_pgm(i),&
					 el_monitor_pgm(i),xr_monitor_pgm(i),yr_monitor_pgm(i),zr_monitor_pgm(i),&
					 highest_pgm_loc, depth_search_mon_pgm)  

               
                          xr_monitor_pgm(i) = xx_spx_loc(n_monitor_pgm(i))
                          yr_monitor_pgm(i) = yy_spx_loc(n_monitor_pgm(i))
                          zr_monitor_pgm(i) = zz_spx_loc(n_monitor_pgm(i))
                          n_monitor_pgm(i) = local_node_num(n_monitor_pgm(i))
                          el_monitor_pgm(i) = local_el_num(el_monitor_pgm(i))

			enddo
			
                allocate(dist_glo(nmonitors_pgm*mpi_np),n_glo(nmonitors_pgm*mpi_np), el_glo(nmonitors_pgm*mpi_np))
                allocate(xr_glo(nmonitors_pgm*mpi_np), yr_glo(nmonitors_pgm*mpi_np), zr_glo(nmonitors_pgm*mpi_np)) 

                call MPI_BARRIER(mpi_comm, mpi_ierr)
                        
                call MPI_ALLGATHER(dist_monitor_pgm, nmonitors_pgm, SPEED_DOUBLE, dist_glo, nmonitors_pgm, &
                                   SPEED_DOUBLE, mpi_comm, mpi_ierr)
                call MPI_ALLGATHER(xr_monitor_pgm, nmonitors_pgm, SPEED_DOUBLE, xr_glo, nmonitors_pgm, &
                                   SPEED_DOUBLE, mpi_comm, mpi_ierr)
                call MPI_ALLGATHER(yr_monitor_pgm, nmonitors_pgm, SPEED_DOUBLE, yr_glo, nmonitors_pgm, &
                                   SPEED_DOUBLE, mpi_comm, mpi_ierr)
                call MPI_ALLGATHER(zr_monitor_pgm, nmonitors_pgm, SPEED_DOUBLE, zr_glo, nmonitors_pgm, &
                                   SPEED_DOUBLE, mpi_comm, mpi_ierr)

                call MPI_ALLGATHER(n_monitor_pgm, nmonitors_pgm, SPEED_INTEGER, n_glo, nmonitors_pgm, &
                                   SPEED_INTEGER, mpi_comm, mpi_ierr)
                call MPI_ALLGATHER(el_monitor_pgm, nmonitors_pgm, SPEED_INTEGER, el_glo, nmonitors_pgm, &
                                   SPEED_INTEGER, mpi_comm, mpi_ierr)
                
                call GET_MINVALUES(n_glo, dist_glo, nmonitors_pgm*mpi_np, n_monitor_pgm, nmonitors_pgm, mpi_np)
        

                j = 1
                do while(j .le. nmonitors_pgm)
                   call GET_INDLOC_FROM_INDGLO(n_glo, nmonitors_pgm*mpi_np, n_monitor_pgm(j), ic)         
                   xr_monitor_pgm(j) = xr_glo(ic)
                   yr_monitor_pgm(j) = yr_glo(ic)
                   zr_monitor_pgm(j) = zr_glo(ic)
                   el_monitor_pgm(j) = el_glo(ic)                                   
                   j=j+1
                enddo

                deallocate(dist_glo, n_glo, el_glo, xr_glo, yr_glo, zr_glo, dist_monitor_pgm)

		file_MPGM = 'MPGM.input'
		call WRITE_FILE_MPGM(file_MPGM, nmonitors_pgm, n_monitor_pgm, el_monitor_pgm,&
		                              xr_monitor_pgm, yr_monitor_pgm, zr_monitor_pgm)

	 else ! YES, it exists an input file with the position of PGM monitors

		file_MPGM = 'MPGM.input'
		call READ_FILE_MPGM(file_MPGM, nmonitors_pgm, n_monitor_pgm, el_monitor_pgm,&
		                              xr_monitor_pgm, yr_monitor_pgm, zr_monitor_pgm)
 
	 endif

		 deallocate(highest_pgm_loc)


	  if (nmonitors_pgm .gt. 0) then
	       allocate(max_u(nmonitors_pgm,9))
	       allocate(max_v(nmonitors_pgm,9))
	       allocate(max_a(nmonitors_pgm,9))
	       allocate(max_o(nmonitors_pgm,3))                                                              
	       max_u = 0.0d0 
               max_v = 0.0d0
	       max_a = 0.0d0
	       max_o = 0.0d0 
	   endif

      endif

      
      if (mpi_id.eq.0) then
         write(*,'(A)')
         write(*,'(A,I10)') 'Monitored points Peak Ground Map: ',nmonitors_pgm
         write(*,'(A,I2)') 'File MPGM : ',file_mon_pgm
      endif



!*****************************************************************************************
!                      MONITORED NODES - Monitor LiST -
!*****************************************************************************************


      nmonitors_lst = 0
      if (num_lst .eq. 1) then

         
	 file_LS = 'LS.input'
	 call READ_DIME_FILEPG(file_LS,nmonitors_lst)     

	 allocate(n_monitor_lst(nmonitors_lst),el_monitor_lst(nmonitors_lst),dist_monitor_lst(nmonitors_lst))
	 allocate(x_monitor_lst(nmonitors_lst),y_monitor_lst(nmonitors_lst),z_monitor_lst(nmonitors_lst))
	 allocate(xr_monitor_lst(nmonitors_lst),yr_monitor_lst(nmonitors_lst),zr_monitor_lst(nmonitors_lst))
	 

	 call READ_FILEPG(file_LS,nmonitors_lst,x_monitor_lst,y_monitor_lst,z_monitor_lst)		
          
	 if (file_mon_lst.eq.0) then ! NO input file with the position of LST monitors
                allocate(x_monitor_real(nmonitors_lst), y_monitor_real(nmonitors_lst), z_monitor_real(nmonitors_lst))

		do i = 1,nmonitors_lst
			call GET_NEAREST_NODE_PGM(nnod_loc, xx_spx_loc, yy_spx_loc, zz_spx_loc,&
					          x_monitor_lst(i), y_monitor_lst(i), z_monitor_lst(i),&
						  n_monitor_lst(i), dist_monitor_lst(i), depth_search_mon_lst)	
              
			call GET_PNT_POS_PGM(nelem_loc,&
					 alfa11,alfa12,alfa13,alfa21,alfa22,alfa23,&
					 alfa31,alfa32,alfa33,beta11,beta12,beta13,&
	 				 beta21,beta22,beta23,beta31,beta32,beta33,&
					 gamma1,gamma2,gamma3,delta1,delta2,delta3,&
					 x_monitor_lst(i),y_monitor_lst(i),z_monitor_lst(i),&
					 el_monitor_lst(i),xr_monitor_lst(i),yr_monitor_lst(i),zr_monitor_lst(i),&
					 highest_lst_loc, depth_search_mon_lst)  

                
                          x_monitor_real(i) = xx_spx_loc(n_monitor_lst(i))
                          y_monitor_real(i) = yy_spx_loc(n_monitor_lst(i))
                          z_monitor_real(i) = zz_spx_loc(n_monitor_lst(i))
                          n_monitor_lst(i) = local_node_num(n_monitor_lst(i))
                          el_monitor_lst(i) = local_el_num(el_monitor_lst(i))
                          

		enddo



                
                j = 1
                do while(j .le. nmonitors_lst)
                
                	allocate(dist_glo(mpi_np),n_glo(mpi_np), el_glo(mpi_np))
                        allocate(xr_glo(mpi_np), yr_glo(mpi_np), zr_glo(mpi_np)) 
                        allocate(x_glo_real(mpi_np), y_glo_real(mpi_np),z_glo_real(mpi_np)) 

                        call MPI_BARRIER(mpi_comm, mpi_ierr)
                        
                call MPI_ALLGATHER(dist_monitor_lst(j), 1, SPEED_DOUBLE, dist_glo, 1, SPEED_DOUBLE, mpi_comm, mpi_ierr)
                call MPI_ALLGATHER(xr_monitor_lst(j), 1, SPEED_DOUBLE, xr_glo, 1, SPEED_DOUBLE, mpi_comm, mpi_ierr)
                call MPI_ALLGATHER(yr_monitor_lst(j), 1, SPEED_DOUBLE, yr_glo, 1, SPEED_DOUBLE, mpi_comm, mpi_ierr)
                call MPI_ALLGATHER(zr_monitor_lst(j), 1, SPEED_DOUBLE, zr_glo, 1, SPEED_DOUBLE, mpi_comm, mpi_ierr)
                call MPI_ALLGATHER(n_monitor_lst(j), 1, SPEED_INTEGER, n_glo, 1, SPEED_INTEGER, mpi_comm, mpi_ierr)
                call MPI_ALLGATHER(el_monitor_lst(j), 1, SPEED_INTEGER, el_glo, 1, SPEED_INTEGER, mpi_comm, mpi_ierr)
                call MPI_ALLGATHER(x_monitor_real(j), 1, SPEED_DOUBLE, x_glo_real, 1, SPEED_DOUBLE, mpi_comm, mpi_ierr)
                call MPI_ALLGATHER(y_monitor_real(j), 1, SPEED_DOUBLE, y_glo_real, 1, SPEED_DOUBLE, mpi_comm, mpi_ierr)
                call MPI_ALLGATHER(z_monitor_real(j), 1, SPEED_DOUBLE, z_glo_real, 1, SPEED_DOUBLE, mpi_comm, mpi_ierr)
                                   
                                   
                
                call GET_MINVALUES(n_glo, dist_glo, mpi_np, n_monitor_lst(j), 1, mpi_np)
        
                call GET_INDLOC_FROM_INDGLO(n_glo, mpi_np, n_monitor_lst(j), ic)              

                   x_monitor_real(j) = x_glo_real(ic)
                   y_monitor_real(j) = y_glo_real(ic)
                   z_monitor_real(j) = z_glo_real(ic)
                   xr_monitor_lst(j) = xr_glo(ic)
                   yr_monitor_lst(j) = yr_glo(ic)
                   zr_monitor_lst(j) = zr_glo(ic)
                   el_monitor_lst(j) = el_glo(ic)                                   
                   j=j+1

                deallocate(dist_glo, n_glo, el_glo, xr_glo, yr_glo, zr_glo, &
                           x_glo_real, y_glo_real, z_glo_real)
                enddo
                deallocate(dist_monitor_lst)

                if(mpi_id.eq. 0) then 
			file_MLST = 'MLST.input'
			call WRITE_FILE_MPGM(file_MLST, nmonitors_lst, n_monitor_lst, el_monitor_lst, &
                                              xr_monitor_lst, yr_monitor_lst, zr_monitor_lst)
			filename  = 'MLST.position'
			call WRITE_FILE_MPGM(filename, nmonitors_lst, n_monitor_lst, el_monitor_lst, &
                                              x_monitor_real, y_monitor_real, z_monitor_real)
                                              
                                              
                endif                                              
                deallocate(x_monitor_real, y_monitor_real, z_monitor_real)

	 else ! YES, it exists an input file with the position of LST monitors

		file_MLST = 'MLST.input'
		call READ_FILE_MPGM(file_MLST, nmonitors_lst, n_monitor_lst, el_monitor_lst, &
					      xr_monitor_lst, yr_monitor_lst, zr_monitor_lst)
 
	 endif


	 deallocate(highest_lst_loc)	 

      else
         write(*,'(A)') 'MLST key not found!'
      
      endif
      
      
      if (mpi_id.eq.0) then
         write(*,'(A)')
         write(*,'(A,I10)') 'Monitored points LST: ',nmonitors_lst
         write(*,'(A,I2)') 'File MLST : ',file_mon_lst
      endif

      call MPI_BARRIER(mpi_comm, mpi_ierr)

!----------------------------------------------------------------------------------
!	WRITING MONITOR.INFO  
!----------------------------------------------------------------------------------
      
      if(mpi_id .eq. 0) then
      
        allocate(monit_files(mpi_np)) 
        monit_files = 0
             
        do i = 1, nmonitors_lst
           monit_files(elem_domain(el_monitor_lst(i))+1) = monit_files(elem_domain(el_monitor_lst(i))+1) + 1
        enddo
       
        if(len_trim(monitor_file) .ne. 70) then              
           monitor_file_new = monitor_file(1:len_trim(monitor_file)) // '/MONITOR.INFO'
        else
           monitor_file_new = 'MONITOR.INFO'  
        endif
     
     
        open(unit=50,file=monitor_file_new)
        write(50,*) tstop
        write(50,*) deltat
        write(50,*) ndt_mon_lst
        write(50,*) mpi_np
        do i = 1, mpi_np
           write(50,*) monit_files(i)
        enddo 
        close(50)
      
        deallocate(monit_files)
      endif      



      
      end subroutine FIND_MONITOR_POSITION
