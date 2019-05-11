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

!> @brief Writes file DGFS.input, containing infos for computing integrals on DG interfaces.
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] nm  number of materials
!> @param[in] sd polynomial degree vector 
!> @param[in] tag_mat label for materials
!> @param[in] cs_nnz_loc length of cs_loc
!> @param[in] cs_loc local spectral connectivity
!> @param[in] nn_loc number of local nodes
!> @param[in] local_n_num local node numeration
!> @param[in] number of local elements
!> @param[in] local_el_num local element numeration
!> @param[in] xs x-coordinate of GLL nodes
!> @param[in] ys y-coordinate of GLL nodes
!> @param[in] zs z-coordinate of GLL nodes
!> @param[in] nel_dg_loc  number of local DG elements
!> @param[in] nel_dg_glo  number of global DG elements
!> @param[in] i4count  vector containing info for DG interface (1 in the i-th position if the i-th 
!!                node lies on a DG interface)
!> @param[in] mpi_id  mpi processor identity
!> @param[in]  mpi_comm  mpi communicator
!> @param[in] mpi_np  mpi processor number
!> @param[in] alfa11 bilinear mapping from (-1,1)^3 to the current el
!> @param[in] alfa12 bilinear mapping from (-1,1)^3 to the current el
!> @param[in] alfa13 bilinear mapping from (-1,1)^3 to the current el
!> @param[in] alfa21 bilinear mapping from (-1,1)^3 to the current el
!> @param[in] alfa22 bilinear mapping from (-1,1)^3 to the current el
!> @param[in] alfa23 bilinear mapping from (-1,1)^3 to the current el
!> @param[in] alfa31 bilinear mapping from (-1,1)^3 to the current el
!> @param[in] alfa32 bilinear mapping from (-1,1)^3 to the current el
!> @param[in] alfa33 bilinear mapping from (-1,1)^3 to the current el
!> @param[in] beta11 bilinear mapping from (-1,1)^3 to the current el
!> @param[in] beta12 bilinear mapping from (-1,1)^3 to the current el
!> @param[in] beta13 bilinear mapping from (-1,1)^3 to the current el
!> @param[in] beta21 bilinear mapping from (-1,1)^3 to the current el
!> @param[in] beta22 bilinear mapping from (-1,1)^3 to the current el
!> @param[in] beta23 bilinear mapping from (-1,1)^3 to the current el
!> @param[in] beta31 bilinear mapping from (-1,1)^3 to the current el
!> @param[in] beta32 bilinear mapping from (-1,1)^3 to the current el
!> @param[in] beta33 bilinear mapping from (-1,1)^3 to the current el
!> @param[in] gamma1 bilinear mapping from (-1,1)^3 to the current el
!> @param[in] gamma2 bilinear mapping from (-1,1)^3 to the current el
!> @param[in] gamma3 bilinear mapping from (-1,1)^3 to the current el
!> @param[in] delta1 bilinear mapping from (-1,1)^3 to the current el
!> @param[in] delta2 bilinear mapping from (-1,1)^3 to the current el
!> @param[in] delta3 bilinear mapping from (-1,1)^3 to the current el
!> @param[in] faces  data structure containing info about DG faces (material, element, face)
!> @param[in] area_nodes  data structure containing info about DG faces (area, constants for the 
!!                   bilinear mapping)
!> @param[in] dg_els  data structure for dg interface elements --> see module.f90
!> @param[in] scratch_dg_els  temporary data structure for dg interface elements --> see module.f90!
!> @param[in] mpi_file directory where reading *.mpi files
!> @param[out] filename (DGFS.input)  file containing info for DG interface integrals
 
 
     subroutine WRITE_FILE_DGFS(nm, sd, tag_mat, cs_nnz_loc, cs_loc, &
                           nn_loc, local_n_num, ne_loc, local_el_num, &
                           xs,ys,zs,&
                           nel_dg_loc, nel_dg_glo, &
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
                           filename,mpi_file)

     use max_var
     use str_mesh 
     use str_mesh_scratch                   

     implicit none

                    
     include 'SPEED.MPI'

     
     type(ELEMENT), dimension(nel_dg_loc), intent(inout) :: dg_els
     type(scratch_ELEMENT), dimension(nel_dg_loc), intent(inout) :: scratch_dg_els
 
     character*70 :: filename, filempi, filenorm, mpi_file, filempi_new
     
     integer*4 :: nm, cs_nnz_loc, nn_loc, ne_loc, nel_dg_loc, nel_dg_glo
     integer*4 :: mpi_comm, mpi_id, mpi_np
     integer*4 :: im, nn, ie, ned, yon
     integer*4 :: ne1, ne2, ne3, ne4, ic1, ic2, ic3, ic4
     integer*4 :: ne5, ne6, ne7, ne8, ic5, ic6, ic7, ic8
     integer*4 :: el_conf, face_conf, face_found, unitmpi, unitname, unitnorm
     integer*4 :: ip, k, j, i, it, ih, ik, tt, indice, node_not_ass, ic, dim1, dim2, jstart
     integer*4 :: error, ncol
     integer*4 :: mpierror, nofel

     integer*4, dimension(2) :: dims, dimsfi
     integer*4, dimension(nm) :: tag_mat, sd
     integer*4, dimension(0:cs_nnz_loc) :: cs_loc
     integer*4, dimension(nn_loc) :: local_n_num
     integer*4, dimension(ne_loc) :: local_el_num
     integer*4, dimension(mpi_np) :: count_dg
     
     integer*4, dimension(:,:), allocatable :: DG_CON_ALL 
     integer*4, dimension(3,nel_dg_glo) :: faces
     integer*4, dimension(nel_dg_glo,3) :: mat_el_fac

     real*8 :: normal_x, normal_y, normal_z
     real*8 :: c_alfa11, c_alfa12, c_alfa13, c_alfa21, c_alfa22, c_alfa23, c_alfa31, c_alfa32, c_alfa33
     real*8 :: c_beta11, c_beta12, c_beta13, c_beta21, c_beta22, c_beta23, c_beta31, c_beta32, c_beta33
     real*8 :: c_gamma1, c_gamma2, c_gamma3, c_delta1, c_delta2, c_delta3
     real*8 :: xnod, ynod, znod, csi, eta, zeta, xnod1, ynod1, znod1
     real*8 :: val1, val2, val3, val4, val5, val6, valmin
     real*8 :: coef_a, coef_b, coef_t, det_trasf

     real*8, dimension(nn_loc) :: xs,ys,zs
     
     real*8, dimension(:), allocatable :: ctgl,wwgl
     real*8, dimension(:), allocatable :: ct, ww
     real*8, dimension(ne_loc) :: alfa11,alfa12,alfa13
     real*8, dimension(ne_loc) :: alfa21,alfa22,alfa23
     real*8, dimension(ne_loc) :: alfa31,alfa32,alfa33
     real*8, dimension(ne_loc) :: beta11,beta12,beta13
     real*8, dimension(ne_loc) :: beta21,beta22,beta23
     real*8, dimension(ne_loc) :: beta31,beta32,beta33
     real*8, dimension(ne_loc) :: gamma1,gamma2,gamma3
     real*8, dimension(ne_loc) :: delta1,delta2,delta3

     real*8, dimension(:,:), allocatable :: dd
     real*8, dimension(:,:), allocatable :: DG_PW_ALL
     real*8, dimension(25,nel_dg_glo) :: area_nodes
     real*8, dimension(3,3) :: rot_mat
     real*8, dimension(nel_dg_glo,3) :: normalxyz
     
     filempi = 'DGFS000000.mpi'
     unitmpi = 50
     
     if (mpi_id .lt. 10) then                                        
        write(filempi(10:10),'(i1)') mpi_id                
     else if (mpi_id .lt. 100) then                                
        write(filempi(9:10),'(i2)') mpi_id                
     else if (mpi_id .lt. 1000) then                                
        write(filempi(8:10),'(i3)') mpi_id                
     else if (mpi_id .lt. 10000) then                                
        write(filempi(7:10),'(i4)') mpi_id                
     else if (mpi_id .lt. 100000) then        
        write(filempi(6:10),'(i5)') mpi_id                
     else if (mpi_id .lt. 1000000) then                                
        write(filempi(5:10),'(i6)') mpi_id                
     endif

     if(len_trim(mpi_file) .ne. 70) then                                                                                  
         filempi_new = mpi_file(1:len_trim(mpi_file)) // '/' // filempi
     else 
         filempi_new = filempi     
     endif
                                      
     open(unitmpi,file=filempi_new)                        


     filenorm = 'NORMALL.input';  unitnorm = 40

     open(unitnorm,file=filenorm)                        
     read(unitnorm,*) nofel
     do i = 1, nofel
       read(unitnorm,*) mat_el_fac(i,1),mat_el_fac(i,2),mat_el_fac(i,3) ,&
                       normalxyz(i,1), normalxyz(i,2), normalxyz(i,3)
     enddo
     
     close(unitnorm)
     
                        
       
       
     ic = 0
             
     ned = cs_loc(0) - 1
     node_not_ass = 0

     do it = 1, nel_dg_loc   
                  
        if(dg_els(it)%proj_yn .eq. 1) then 


           
           do i = 1, (dg_els(it)%quad_rule)**2

              tt = 0
              ih = 1
                             
              do while (tt.eq.0 .and. ih.le. nel_dg_glo)   
              
                 if( faces(1,ih) .ne. dg_els(it)%mat) then
                              
                    !ik = faces(2,ih)
                    !il = faces(3,ih)
                      
                    !CHECK THE NORMAL !!!
                    call CHECK_NORMAL(dg_els(it)%nx, dg_els(it)%ny, dg_els(it)%nz, &
                                       faces(1,ih), faces(2,ih), faces(3,ih), &
                                       nel_dg_glo, normalxyz, mat_el_fac, yon)
                    
                    
                    if(yon .eq. 1) then 
   
                       xnod =  scratch_dg_els(it)%x_nq(i)
                       ynod =  scratch_dg_els(it)%y_nq(i)
                       znod =  scratch_dg_els(it)%z_nq(i)

                       call MAKE_BILINEAR_MAP(area_nodes(2:25,ih), &
                                             c_alfa11, c_alfa12, c_alfa13, &
                                             c_alfa21, c_alfa22, c_alfa23, &
                                             c_alfa31, c_alfa32, c_alfa33, &  
                                             c_beta11, c_beta12, c_beta13, & 
                                             c_beta21, c_beta22, c_beta23, & 
                                             c_beta31, c_beta32, c_beta33, &
                                             c_gamma1, c_gamma2, c_gamma3, &
                                             c_delta1, c_delta2, c_delta3)
 
                       call NEWTON_RAPSON(xnod, ynod, znod, &
                                           c_alfa11, c_alfa12, c_alfa13, &
                                           c_alfa21, c_alfa22, c_alfa23, &
                                           c_alfa31, c_alfa32, c_alfa33, &  
                                           c_beta11, c_beta12, c_beta13, & 
                                           c_beta21, c_beta22, c_beta23, & 
                                           c_beta31, c_beta32, c_beta33, &
                                           c_gamma1, c_gamma2, c_gamma3, &
                                           c_delta1, c_delta2, c_delta3, & 
                                           tt, csi, eta, zeta, nofinr, mpi_id,&
                                           dg_els(it)%ind_el,faces(2,ih), 1.d-6, 1.01d0, 1)
                                           
                                           
                        if(tt == 1) then
                       
                           ic = ic + 1
                            
                           write(unitmpi,"(1I2,1X,1I12,1X,1I2,1X,1I2,1X,1I12,1X,1I2,3(1X,ES25.16),3(1X,ES25.16))") &
                           dg_els(it)%mat, dg_els(it)%ind_el, dg_els(it)%face_el, &
                           faces(1,ih), faces(2,ih), faces(3,ih), &
                           xnod, ynod, znod, &
                           dg_els(it)%wx_pl(i),dg_els(it)%wy_pl(i),dg_els(it)%wz_pl(i) 
                           
                        endif ! if (tt == 1)                  
                                           
 
 
                    else
                    
                       tt = 0 
                    
                  
                    endif ! (yon .eq. 1)      
                   
                 endif ! (faces(ih,1) .ne. dg_els(it)%mat)         
              
                 ih = ih + 1
                  
              enddo  ! while tt == 0
    
              !  CHECK ON NODE NOT ASSIGNED
              if(tt .eq. 0 .and. ih .gt. nel_dg_glo) then
                  node_not_ass = node_not_ass + 1
                  write(*,*) mpi_id ,'NODE', i,' NOT ASSIGNED', ' el', dg_els(it)%ind_el
              endif
                 

           enddo ! i         
       
        endif !if proj_yn == 1
       
     enddo ! it 

     close(unitmpi)        

     
     dims(1) = 6
     dims(2) = ic
     
     if(ic .ne. 0) then 
        write(*,'(A,I5,A,I6)') 'Proc. :', mpi_id, ' not assigned nodes : ', node_not_ass
        write(*,'(A,I5,A,I6)') 'Proc. :', mpi_id, ' assigned nodes     : ', ic        
     else
        write(*,'(A,I5)') 'No projection for Proc. :', mpi_id             
     endif                                                              

!     write(unitmpi,*) 6, ic                
!     do i = 1, ic
!         write(unitmpi,"(1I2,1X,1I12,1X,1I2,1X,1I2,1X,1I12,1X,1I2,3(1X,ES16.9),3(1X,ES16.9))") &
!                 DG_CON(1,i), DG_CON(2,i), DG_CON(3,i), DG_CON(4,i), DG_CON(5,i), DG_CON(6,i), & 
!                 DG_PW(1,i), DG_PW(2,i), DG_PW(3,i), DG_PW(4,i), DG_PW(5,i), DG_PW(6,i) 
!         write(unit_mpi,"(1I2,1X,1I12,1X,1I2,1X,1I2,1X,1I12,1X,1I2,3(1X,ES25.16),3(1X,ES25.16))") &
!          write(*,*)              
!     enddo



     call MPI_BARRIER(mpi_comm, mpierror)                  
     call MPI_ALLGATHER(dims(2), 1, SPEED_INTEGER, count_dg, 1, SPEED_INTEGER, mpi_comm, mpierror)

     dimsfi(1) = 6
     dimsfi(2) = sum(count_dg)     
     
 !    write(*,*) 'qui in fondo'


    if(mpi_id .eq. 0) then 
    
       allocate(DG_CON_ALL(dimsfi(1),dimsfi(2)), DG_PW_ALL(dimsfi(1),dimsfi(2)))

       do i = 1, mpi_np
       
          
          if (count_dg(i) .ne. 0) then
          
              filempi = 'DGFS000000.mpi'
              unitmpi = 40                                 
              if (i-1 .lt. 10) then                                        
                  write(filempi(10:10),'(i1)') i-1                
              else if (i-1 .lt. 100) then                                
                  write(filempi(9:10),'(i2)') i-1
              else if (i-1 .lt. 1000) then                                
                  write(filempi(8:10),'(i3)') i-1
              else if (i-1 .lt. 10000) then                                
                  write(filempi(7:10),'(i4)') i-1
              else if (i-1 .lt. 100000) then                                
                  write(filempi(6:10),'(i5)') i-1
              else if (i-1 .lt. 1000000) then                                
                  write(filempi(5:10),'(i6)') i-1
              endif                                                              

              if(len_trim(mpi_file) .ne. 70) then                                                                                  
                 filempi_new = mpi_file(1:len_trim(mpi_file)) // '/' // filempi
              else 
                 filempi_new = filempi     
              endif
              
              open(unitmpi,file=filempi_new)        
              dim1 = 6
              dim2 = count_dg(i)
              
              if(i.eq. 1) then
                 jstart = 0
              else
                 jstart = sum(count_dg(1:i-1))
              endif
                  
              do j = 1, dim2
                read(unitmpi,*)  DG_CON_ALL(1,j+jstart), DG_CON_ALL(2,j+jstart), DG_CON_ALL(3,j+jstart), &
                                  DG_CON_ALL(4,j+jstart), DG_CON_ALL(5,j+jstart), DG_CON_ALl(6,j+jstart), & 
                                  DG_PW_ALL(1,j+jstart), DG_PW_ALL(2,j+jstart), DG_PW_ALL(3,j+jstart), &
                                  DG_PW_ALL(4,j+jstart), DG_PW_ALL(5,j+jstart), DG_PW_ALL(6,j+jstart)              

             !   write(*,*)  DG_CON_ALL(1,j+jstart), DG_CON_ALL(2,j+jstart), DG_CON_ALL(3,j+jstart), &
             !                     DG_CON_ALL(4,j+jstart), DG_CON_ALL(5,j+jstart), DG_CON_ALl(6,j+jstart), & 
             !                     DG_PW_ALL(1,j+jstart), DG_PW_ALL(2,j+jstart), DG_PW_ALL(3,j+jstart), &
             !                     DG_PW_ALL(4,j+jstart), DG_PW_ALL(5,j+jstart), DG_PW_ALL(6,j+jstart)              



              enddo
          
              close(unitmpi)
          
          endif
      
       enddo

      unitname = 500 
      open(unitname,file=filename)
      write(unitname,"(1I15)") dimsfi(2)

      do j = 1, dimsfi(2)

        write(unitname,"(1I2,1X,1I12,1X,1I2,1X,1I2,1X,1I12,1X,1I2,3(1X,ES25.16),3(1X,ES25.16))") &
                                  DG_CON_ALL(1,j), DG_CON_ALL(2,j), DG_CON_ALL(3,j), &
                                  DG_CON_ALL(4,j), DG_CON_ALL(5,j), DG_CON_ALl(6,j), & 
                                  DG_PW_ALL(1,j), DG_PW_ALL(2,j), DG_PW_ALL(3,j), &
                                  DG_PW_ALL(4,j), DG_PW_ALL(5,j), DG_PW_ALL(6,j)     



      enddo

      close(unitname)

         
       deallocate(DG_CON_ALL, DG_PW_ALL)  

    endif

    call MPI_BARRIER(mpi_comm, error) 

        
    return
    
    
     end subroutine WRITE_FILE_DGFS
