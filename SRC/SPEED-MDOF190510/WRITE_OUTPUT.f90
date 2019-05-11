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

!> @brief Writes output results.
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] nmonitlst  nmber of monitors
!> @param[in] mpi_id  id mpi process
!> @param[in] elem_mlst  list of elements for monitors
!> @param[in] ne_loc number of local elements
!> @param[in] local_el_num  local element numeration
!> @param[in] cs_nnz_loc length of cs_loc
!> @param[in] cs_loc  local connectivity vector
!> @param[in] nmat number of materials
!> @param[in] sdeg_mat  polynomial degree vector 
!> @param[in] nnod_loc number of local nodes
!> @param[in] u2 displacement variables time t+1
!> @param[in] u1 displacement variables time t
!> @param[in] u0 displacement variables time t-1
!> @param[in] v1 velocity variables time t
!> @param[in] stress stress tensor
!> @param[in] strain strain tensor
!> @param[in] omega rotational tensor
!> @param[in] xr_mlst monitors x-coordinate
!> @param[in] yr_mlst monitors y-coordinate
!> @param[in] zr_mlst  monitors coordinate
!> @param[in] tt1 current time
!> @param[in] dt time step
!> @param[in] dt2 dt^2
!> @param[in] option_out_var  option for output
!> @param[in] count_mon  number of current monitors
!> @param[in] monitor_file directory where saving files
!> @param[in] dbg flag for debugging mode 
!> @param[in] mu second elastic coefficient to print
!> @param[in] gamma damping coefficient to print
!> @param[in] b_instabilitycontrol flag to check numerical instability
!> @param[in] instability_maxval  check numerical instability
!> @param[out] b_instability_abort if true, simulation is unstable
!> @param[out] --- Files MONITORXXXXXX.D (disp),MONITORXXXXXX.V (vel),MONITORXXXXXX.A (acc)
!            MONITORXXXXXX.S (stress),MONITORXXXXXX.E (strain),MONITORXXXXXX.O (rot.) 


     subroutine WRITE_OUTPUT(nmonitlst, mpi_id, elem_mlst, local_el_num, ne_loc, &
                              cs_loc, cs_nnz_loc, sdeg_mat, nmat, &
                              u2, u1, u0, v1, nnod_loc, &
                              stress, strain, omega, &
                              xr_mlst, yr_mlst, zr_mlst, & 
                              tt1, dt, dt2, &
                              option_out_var, count_mon, monitor_file, &
                              dbg, mu, gamma, b_instabilitycontrol, instability_maxval, b_instability_abort)
                              
     
     implicit none
     
     character*70 :: monitor_file
     character*70 :: monitor_file1,monitor_file2,monitor_file3,monitor_file4,monitor_file5,monitor_file6
     character*70 :: monitor_file7,monitor_file8
     character*14 :: file_monitor1,file_monitor2,file_monitor3,file_monitor4,file_monitor5,file_monitor6
     character*14 :: file_monitor7,file_monitor8
     integer*4 :: tycount !!!190510
     integer*4 :: unit_monitor1,unit_monitor2,unit_monitor3,unit_monitor4,unit_monitor5,unit_monitor6
     integer*4 :: unit_monitor7,unit_monitor8
     integer*4 :: imon,ielem,ie,im,nn,k,j,i,is,in,iaz,count_mon,ishift,jshift,kshift,dbg
     integer*4 :: nmonitlst, mpi_id, ne_loc, nmat, nnod_loc, cs_nnz_loc
     integer*4, dimension(0:cs_nnz_loc) :: cs_loc
     integer*4, dimension(nmonitlst) :: elem_mlst
     integer*4, dimension(nmat) :: sdeg_mat
     integer*4, dimension(ne_loc) :: local_el_num
     integer*4, dimension(6) :: option_out_var
     integer*4 :: mpi_ierr
     
     real*8 :: tt1, dt, dt2
     real*8 :: uxm,uym,uzm
     real*8 :: vxm,vym,vzm
     real*8 :: axm,aym,azm
     real*8 :: variable1m,variable2m,variable3m,variable4m,variable5m,variable6m

     real*8, dimension(:), allocatable :: variables_1,variables_2,variables_3
     real*8, dimension(:), allocatable :: variables_4,variables_5,variables_6
     real*8, dimension(:), allocatable :: variables_7,variables_8
          
     real*8, dimension(:), allocatable :: ct, ww
     real*8, dimension(:,:), allocatable :: dd
     real*8, dimension(:,:,:), allocatable :: ux_el, uy_el, uz_el
     real*8, dimension(:,:,:), allocatable :: variable1_el,variable2_el,variable3_el
     real*8, dimension(:,:,:), allocatable :: variable4_el,variable5_el,variable6_el
          
     real*8, dimension(3*nnod_loc) :: u1, u0, u2, v1, omega
     real*8, dimension(6*nnod_loc) :: stress, strain
     real*8, dimension(nnod_loc) :: mu, gamma
     real*8, dimension(nmonitlst) ::xr_mlst, yr_mlst, zr_mlst

     ! Instability control
     real*8, intent(in)  :: instability_maxval
     logical, intent(in) :: b_instabilitycontrol
     logical, intent(out) :: b_instability_abort

     
     unit_monitor1 = 10*(mpi_id+1) + 1
     unit_monitor2 = 10*(mpi_id+1) + 2
     unit_monitor3 = 10*(mpi_id+1) + 3
     unit_monitor4 = 10*(mpi_id+1) + 4
     unit_monitor5 = 10*(mpi_id+1) + 5
     unit_monitor6 = 10*(mpi_id+1) + 6
     
     if (dbg .eq. 1) then
         unit_monitor7 = 100*(mpi_id+1) + 1
         unit_monitor8 = 100*(mpi_id+1) + 2
     endif    
          
     if (option_out_var(1) .eq. 1) file_monitor1 = 'MONITOR00000.D'
     if (option_out_var(2) .eq. 1) file_monitor2 = 'MONITOR00000.V'
     if (option_out_var(3) .eq. 1) file_monitor3 = 'MONITOR00000.A'
     if (option_out_var(4) .eq. 1) file_monitor4 = 'MONITOR00000.S'
     if (option_out_var(5) .eq. 1) file_monitor5 = 'MONITOR00000.E'
     if (option_out_var(6) .eq. 1) file_monitor6 = 'MONITOR00000.O'
     if (dbg .eq. 1) then
         file_monitor7 = 'MONITOR00000.M'
         file_monitor8 = 'MONITOR00000.G'     
      endif    
          
            
     if (mpi_id .lt. 10) then                                        
           write(file_monitor1(12:12),'(i1)') mpi_id;
           write(file_monitor2(12:12),'(i1)') mpi_id;
           write(file_monitor3(12:12),'(i1)') mpi_id;
           write(file_monitor4(12:12),'(i1)') mpi_id;
           write(file_monitor5(12:12),'(i1)') mpi_id;
           write(file_monitor6(12:12),'(i1)') mpi_id;
           write(file_monitor7(12:12),'(i1)') mpi_id;
           write(file_monitor8(12:12),'(i1)') mpi_id;
                           
     else if (mpi_id .lt. 100) then                                
           write(file_monitor1(11:12),'(i2)') mpi_id;                       
           write(file_monitor2(11:12),'(i2)') mpi_id;                       
           write(file_monitor3(11:12),'(i2)') mpi_id;                       
           write(file_monitor4(11:12),'(i2)') mpi_id;                       
           write(file_monitor5(11:12),'(i2)') mpi_id;                       
           write(file_monitor6(11:12),'(i2)') mpi_id;                       
           write(file_monitor7(11:12),'(i1)') mpi_id;
           write(file_monitor8(11:12),'(i1)') mpi_id;
           
     else if (mpi_id .lt. 1000) then                                
           write(file_monitor1(10:12),'(i3)') mpi_id;                
           write(file_monitor2(10:12),'(i3)') mpi_id;                
           write(file_monitor3(10:12),'(i3)') mpi_id;                
           write(file_monitor4(10:12),'(i3)') mpi_id;                
           write(file_monitor5(10:12),'(i3)') mpi_id;                
           write(file_monitor6(10:12),'(i3)') mpi_id;                
           write(file_monitor7(10:12),'(i3)') mpi_id;                
           write(file_monitor8(10:12),'(i3)') mpi_id;                

     else if (mpi_id .lt. 10000) then                                
           write(file_monitor1(9:12),'(i4)') mpi_id;                
           write(file_monitor2(9:12),'(i4)') mpi_id;                
           write(file_monitor3(9:12),'(i4)') mpi_id;                
           write(file_monitor4(9:12),'(i4)') mpi_id;                
           write(file_monitor5(9:12),'(i4)') mpi_id;                
           write(file_monitor6(9:12),'(i4)') mpi_id;                
           write(file_monitor7(9:12),'(i4)') mpi_id;                
           write(file_monitor8(9:12),'(i4)') mpi_id;                

     else                                                        
           write(file_monitor1(8:12),'(i5)') mpi_id                
           write(file_monitor2(8:12),'(i5)') mpi_id                
           write(file_monitor3(8:12),'(i5)') mpi_id                
           write(file_monitor4(8:12),'(i5)') mpi_id                
           write(file_monitor5(8:12),'(i5)') mpi_id                
           write(file_monitor6(8:12),'(i5)') mpi_id
           write(file_monitor7(8:12),'(i5)') mpi_id
           write(file_monitor8(8:12),'(i5)') mpi_id
                                      
     endif                                                
        
     if (option_out_var(1) .eq. 1) then
           if(len_trim(monitor_file) .ne. 70) then
             monitor_file1 = monitor_file(1:len_trim(monitor_file)) // '/' // file_monitor1
           else
             monitor_file1 = file_monitor1
           endif     
           allocate(variables_1(0:3*count_mon))
           open(unit_monitor1,file=monitor_file1,position='APPEND')  
           variables_1(0) = tt1
     endif     
     if (option_out_var(2) .eq. 1) then
           if(len_trim(monitor_file) .ne. 70) then
             monitor_file2 = monitor_file(1:len_trim(monitor_file)) // '/' // file_monitor2
           else
             monitor_file2 = file_monitor2
           endif     
           allocate(variables_2(0:3*count_mon))
           open(unit_monitor2,file=monitor_file2,position='APPEND')  
           variables_2(0) = tt1
     endif
     if (option_out_var(3) .eq. 1) then
           if(len_trim(monitor_file) .ne. 70) then
             monitor_file3 = monitor_file(1:len_trim(monitor_file)) // '/' // file_monitor3
           else
             monitor_file3 = file_monitor3
           endif     
           allocate(variables_3(0:3*count_mon))
           open(unit_monitor3,file=monitor_file3,position='APPEND')  
           variables_3(0) = tt1
     endif      
     if (option_out_var(4) .eq. 1) then
           if(len_trim(monitor_file) .ne. 70) then
             monitor_file4 = monitor_file(1:len_trim(monitor_file)) // '/' // file_monitor4
           else
             monitor_file4 = file_monitor4
           endif     
           allocate(variables_4(0:6*count_mon))
           open(unit_monitor4,file=monitor_file4,position='APPEND')  
           variables_4(0) = tt1
     endif      
     if (option_out_var(5) .eq. 1) then
           if(len_trim(monitor_file) .ne. 70) then                                  
             monitor_file5 = monitor_file(1:len_trim(monitor_file)) // '/' // file_monitor5
           else
             monitor_file5 = file_monitor5
           endif     
           allocate(variables_5(0:6*count_mon))
           open(unit_monitor5,file=monitor_file5,position='APPEND')  
           variables_5(0) = tt1
     endif      
     if (option_out_var(6) .eq. 1) then
           if(len_trim(monitor_file) .ne. 70) then                                                                   
             monitor_file6 = monitor_file(1:len_trim(monitor_file)) // '/' // file_monitor6
           else
             monitor_file6 = file_monitor6
           endif     
           allocate(variables_6(0:3*count_mon))
           open(unit_monitor6,file=monitor_file6,position='APPEND')  
           variables_6(0) = tt1
     endif
     
     if (dbg .eq. 1) then
           if(len_trim(monitor_file) .ne. 70) then                                                                   
             monitor_file7 = monitor_file(1:len_trim(monitor_file)) // '/' // file_monitor7
             monitor_file8 = monitor_file(1:len_trim(monitor_file)) // '/' // file_monitor8
           else
             monitor_file7 = file_monitor7
             monitor_file8 = file_monitor8
           endif     
           allocate(variables_7(0:count_mon),variables_8(0:count_mon))
           open(unit_monitor7,file=monitor_file7,position='APPEND')  
           open(unit_monitor8,file=monitor_file8,position='APPEND')  
           variables_7(0) = tt1;    variables_8(0) = tt1
           
          
     endif
    


     ishift = 0
     jshift = 0
     kshift = 0
     
     do imon = 1, nmonitlst
        
        ielem = elem_mlst(imon)     !ie = index hexahedra containing monitor
        call GET_INDLOC_FROM_INDGLO(local_el_num, ne_loc, ielem, ie)                                             
        if (ie .ne. 0) then
            
            im = cs_loc(cs_loc(ie -1) +0)
            nn = sdeg_mat(im) +1

            allocate(ct(nn),ww(nn),dd(nn,nn))
            allocate(ux_el(nn,nn,nn),uy_el(nn,nn,nn),uz_el(nn,nn,nn))
            allocate(variable1_el(nn,nn,nn),variable2_el(nn,nn,nn),variable3_el(nn,nn,nn))        
            allocate(variable4_el(nn,nn,nn),variable5_el(nn,nn,nn),variable6_el(nn,nn,nn))        
            call MAKE_LGL_NW(nn,ct,ww,dd)

            do k = 1,nn
               do j = 1,nn
                  do i = 1,nn
                     is = nn*nn*(k -1) +nn*(j -1) +i
                     in = cs_loc(cs_loc(ie -1) + is)

                     iaz = 3*(in -1) +1
                     ux_el(i,j,k) = u1(iaz)
                     iaz = 3*(in -1) +2
                     uy_el(i,j,k) = u1(iaz)
                     iaz = 3*(in -1) +3
                     uz_el(i,j,k) = u1(iaz)
                  enddo
               enddo
            enddo

            if (option_out_var(1).eq.1) then  
            !-------------------------------------------------------------
            ! DISPLACEMENT
             !-------------------------------------------------------------
               
                    call GET_MONITOR_VALUE(nn,ct,ux_el,&
                                xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),uxm)
              
                    call GET_MONITOR_VALUE(nn,ct,uy_el,&
                                xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),uym)
               
                    call GET_MONITOR_VALUE(nn,ct,uz_el,&
                                xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),uzm)

                    if (dabs(uxm).lt.1.0e-30) uxm = 0.0e+00
                    if (dabs(uym).lt.1.0e-30) uym = 0.0e+00
                    if (dabs(uzm).lt.1.0e-30) uzm = 0.0e+00


                variables_1(ishift+1) = uxm
                variables_1(ishift+2) = uym
                variables_1(ishift+3) = uzm
            !-------------------------------------------------------------
             endif
                
            if (option_out_var(2).eq.1) then  
            !-------------------------------------------------------------
            ! VELOCITY
            !-------------------------------------------------------------
                   do k = 1,nn
                   do j = 1,nn
                      do i = 1,nn
                         is = nn*nn*(k -1) +nn*(j -1) +i
                         in = cs_loc(cs_loc(ie -1) + is)
                         
                         iaz = 3*(in -1) +1
                         ux_el(i,j,k) = (u2(iaz) - u0(iaz)) / (2.0d0*dt)
                         iaz = 3*(in -1) +2
                         uy_el(i,j,k) = (u2(iaz) - u0(iaz)) / (2.0d0*dt)
                         iaz = 3*(in -1) +3
                         uz_el(i,j,k) = (u2(iaz) - u0(iaz)) / (2.0d0*dt)
                       enddo
                  enddo
                enddo
               
                call GET_MONITOR_VALUE(nn,ct,ux_el,&
                                xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),vxm)
               
                call GET_MONITOR_VALUE(nn,ct,uy_el,&
                                xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),vym)
               
                call GET_MONITOR_VALUE(nn,ct,uz_el,&
                                xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),vzm)

                if (dabs(vxm).lt.1.0e-30) vxm = 0.0e+00
                if (dabs(vym).lt.1.0e-30) vym = 0.0e+00
                if (dabs(vzm).lt.1.0e-30) vzm = 0.0e+00


                variables_2(ishift+1) = vxm
                variables_2(ishift+2) = vym
                variables_2(ishift+3) = vzm
            !-------------------------------------------------------------
            endif        
                        
            if (option_out_var(3).eq.1) then  
            !-------------------------------------------------------------
            ! ACCELERATION
            !-------------------------------------------------------------
                    do k = 1,nn
                      do j = 1,nn
                         do i = 1,nn
                            is = nn*nn*(k -1) +nn*(j -1) +i
                            in = cs_loc(cs_loc(ie -1) + is)

                            iaz = 3*(in -1) +1
                            ux_el(i,j,k) = (u2(iaz) -2.0*u1(iaz) +u0(iaz)) / dt2
                            iaz = 3*(in -1) +2
                            uy_el(i,j,k) = (u2(iaz) -2.0*u1(iaz) +u0(iaz)) / dt2
                            iaz = 3*(in -1) +3
                            uz_el(i,j,k) = (u2(iaz) -2.0*u1(iaz) +u0(iaz)) / dt2
                           enddo
                      enddo
                   enddo
               
                   call GET_MONITOR_VALUE(nn,ct,ux_el,&
                                xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),axm)
             
                   call GET_MONITOR_VALUE(nn,ct,uy_el,&
                                xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),aym)
            
                   call GET_MONITOR_VALUE(nn,ct,uz_el,&
                                xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),azm)
           
                   if (dabs(axm).lt.1.0e-30) axm = 0.0e+00
                   if (dabs(aym).lt.1.0e-30) aym = 0.0e+00
                   if (dabs(azm).lt.1.0e-30) azm = 0.0e+00

                                                
                variables_3(ishift+1) = axm
                variables_3(ishift+2) = aym
                variables_3(ishift+3) = azm
            !-------------------------------------------------------------
              endif

            if (option_out_var(4).eq.1) then  
            !-------------------------------------------------------------
            ! STRESS
            !-------------------------------------------------------------
                do k = 1,nn
                   do j = 1,nn
                      do i = 1,nn
                         is = nn*nn*(k -1) +nn*(j -1) +i
                         in = cs_loc(cs_loc(ie -1) + is)

                         iaz = 6*(in -1) +1
                         variable1_el(i,j,k) = stress(iaz)
                         iaz = 6*(in -1) +2
                         variable2_el(i,j,k) = stress(iaz)
                         iaz = 6*(in -1) +3
                         variable3_el(i,j,k) = stress(iaz)
                         iaz = 6*(in -1) +4
                         variable4_el(i,j,k) = stress(iaz)
                         iaz = 6*(in -1) +5
                         variable5_el(i,j,k) = stress(iaz)
                         iaz = 6*(in -1) +6
                         variable6_el(i,j,k) = stress(iaz)
                      enddo
                   enddo
                 enddo
              
                 call GET_MONITOR_VALUE(nn,ct,variable1_el,&
                                        xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),variable1m)
                 call GET_MONITOR_VALUE(nn,ct,variable2_el,&
                                        xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),variable2m)
                 call GET_MONITOR_VALUE(nn,ct,variable3_el,&
                                        xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),variable3m)
                 call GET_MONITOR_VALUE(nn,ct,variable4_el,&
                                        xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),variable4m)
                 call GET_MONITOR_VALUE(nn,ct,variable5_el,&
                                        xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),variable5m)
                 call GET_MONITOR_VALUE(nn,ct,variable6_el,&
                                        xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),variable6m)

                 if (dabs(variable1m).lt.1.0e-30) variable1m = 0.0e+00
                 if (dabs(variable2m).lt.1.0e-30) variable2m = 0.0e+00
                 if (dabs(variable3m).lt.1.0e-30) variable3m = 0.0e+00
                 if (dabs(variable4m).lt.1.0e-30) variable4m = 0.0e+00
                 if (dabs(variable5m).lt.1.0e-30) variable5m = 0.0e+00
                 if (dabs(variable6m).lt.1.0e-30) variable6m = 0.0e+00

                  ! Instability control
                  ! Check if one of the monitored values is out of bounds
                  if (b_instabilitycontrol) then
                    if (dabs(variable1m) .gt. instability_maxval .or. dabs(variable2m) .gt. instability_maxval .or. &
                        dabs(variable3m) .gt. instability_maxval .or. dabs(variable4m) .gt. instability_maxval .or. &
                        dabs(variable5m) .gt. instability_maxval .or. dabs(variable6m) .gt. instability_maxval) then
                          write(*,*)             '+======================================================================================+'
                          write(*,'(A,E12.4,A)') ' | Instability control: monitored stress out of bounds (maxval=', instability_maxval, '), aborting! |'
                          write(*,*)             '+======================================================================================+'
                          b_instability_abort = .TRUE.
                          return
                    endif
                  endif
                
                variables_4(jshift+1) = variable1m
                variables_4(jshift+2) = variable2m
                variables_4(jshift+3) = variable3m
                variables_4(jshift+4) = variable4m
                variables_4(jshift+5) = variable5m
                variables_4(jshift+6) = variable6m

            !-------------------------------------------------------------
            endif

            if (option_out_var(5).eq.1) then  
            !-------------------------------------------------------------
            ! STRAIN
            !-------------------------------------------------------------
                do k = 1,nn
                   do j = 1,nn
                      do i = 1,nn
                         is = nn*nn*(k -1) +nn*(j -1) +i
                         in = cs_loc(cs_loc(ie -1) + is)

                         iaz = 6*(in -1) +1
                         variable1_el(i,j,k) = strain(iaz)
                         iaz = 6*(in -1) +2
                         variable2_el(i,j,k) = strain(iaz)
                         iaz = 6*(in -1) +3
                         variable3_el(i,j,k) = strain(iaz)
                         iaz = 6*(in -1) +4
                         variable4_el(i,j,k) = strain(iaz)
                         iaz = 6*(in -1) +5
                         variable5_el(i,j,k) = strain(iaz)
                         iaz = 6*(in -1) +6
                         variable6_el(i,j,k) = strain(iaz)
                      enddo
                   enddo
                enddo
               
                call GET_MONITOR_VALUE(nn,ct,variable1_el,&
                                xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),variable1m)
               
                       call GET_MONITOR_VALUE(nn,ct,variable2_el,&
                                xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),variable2m)
               
                call GET_MONITOR_VALUE(nn,ct,variable3_el,&
                                xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),variable3m)

                call GET_MONITOR_VALUE(nn,ct,variable4_el,&
                                xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),variable4m)

                call GET_MONITOR_VALUE(nn,ct,variable5_el,&
                                xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),variable5m)

                call GET_MONITOR_VALUE(nn,ct,variable6_el,&
                                xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),variable6m)

                if (dabs(variable1m).lt.1.0e-30) variable1m = 0.0e+00
                if (dabs(variable2m).lt.1.0e-30) variable2m = 0.0e+00
                if (dabs(variable3m).lt.1.0e-30) variable3m = 0.0e+00
                if (dabs(variable4m).lt.1.0e-30) variable4m = 0.0e+00
                if (dabs(variable5m).lt.1.0e-30) variable5m = 0.0e+00
                if (dabs(variable6m).lt.1.0e-30) variable6m = 0.0e+00


                variables_5(jshift+1) = variable1m
                variables_5(jshift+2) = variable2m
                variables_5(jshift+3) = variable3m
                variables_5(jshift+4) = variable4m
                variables_5(jshift+5) = variable5m
                variables_5(jshift+6) = variable6m

            !-------------------------------------------------------------
            endif

            if (option_out_var(6).eq.1) then  
            !-------------------------------------------------------------
            ! OMEGA
            !-------------------------------------------------------------
                do k = 1,nn
                   do j = 1,nn
                      do i = 1,nn
                         is = nn*nn*(k -1) +nn*(j -1) +i
                         in = cs_loc(cs_loc(ie -1) + is)

                          iaz = 3*(in -1) +1
                         variable1_el(i,j,k) = omega(iaz)
                         iaz = 3*(in -1) +2
                         variable2_el(i,j,k) = omega(iaz)
                         iaz = 3*(in -1) +3
                         variable3_el(i,j,k) = omega(iaz)
                      enddo
                   enddo
                 enddo
               
                 call GET_MONITOR_VALUE(nn,ct,variable1_el,&
                                xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),variable1m)
               
                 call GET_MONITOR_VALUE(nn,ct,variable2_el,&
                                xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),variable2m)
               
                 call GET_MONITOR_VALUE(nn,ct,variable3_el,&
                                xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),variable3m)

                 if (dabs(variable1m).lt.1.0e-30) variable1m = 0.0e+00
                 if (dabs(variable2m).lt.1.0e-30) variable2m = 0.0e+00
                 if (dabs(variable3m).lt.1.0e-30) variable3m = 0.0e+00

                variables_6(ishift+1) = variable1m
                variables_6(ishift+2) = variable2m
                variables_6(ishift+3) = variable3m

            endif
            
            if (dbg.eq.1) then  
            !-------------------------------------------------------------
            ! MU, GAMMA NLE
            !-------------------------------------------------------------
                do k = 1,nn
                   do j = 1,nn
                      do i = 1,nn
                         is = nn*nn*(k -1) +nn*(j -1) +i
                         in = cs_loc(cs_loc(ie -1) + is)

                         
                         variable1_el(i,j,k) = mu(in)
                         variable2_el(i,j,k) = gamma(in)
                      enddo
                   enddo
                 enddo
               
                 call GET_MONITOR_VALUE(nn,ct,variable1_el,&
                                xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),variable1m)
               
                 call GET_MONITOR_VALUE(nn,ct,variable2_el,&
                                xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),variable2m)
               

                 if (dabs(variable1m).lt.1.0e-30) variable1m = 0.0e+00
                 if (dabs(variable2m).lt.1.0e-30) variable2m = 0.0e+00


                variables_7(kshift+1) = variable1m
                variables_8(kshift+1) = variable2m

            endif

            deallocate(ct,ww,dd)
            deallocate(ux_el,uy_el,uz_el)
            deallocate(variable1_el,variable2_el,variable3_el)        
            deallocate(variable4_el,variable5_el,variable6_el)        
            ishift = ishift + 3
            jshift = jshift + 6
            kshift = kshift + 1

        endif
            
        
     enddo
         
      
         
     
     if (option_out_var(1) .eq. 1) then
        do tycount=0,3*count_mon
            write(unit_monitor1,"(E16.7,A1)",advance='NO') variables_1(tycount)," "
        end do        
            write(unit_monitor1,"(A1)") " "
        close(unit_monitor1)
        deallocate(variables_1)
     endif   
     if (option_out_var(2) .eq. 1) then 
        do tycount=0,3*count_mon
            write(unit_monitor2,"(E16.7,A1)",advance='NO') variables_2(tycount)," "
        end do
            write(unit_monitor2,"(A1)") " "
        close(unit_monitor2)
        deallocate(variables_2)
     endif   
     if (option_out_var(3) .eq. 1) then
        do tycount=0,3*count_mon
            write(unit_monitor3,"(E16.7,A1)",advance='NO') variables_3(tycount)," "
        end do
            write(unit_monitor3,"(A1)") " "
        close(unit_monitor3)
        deallocate(variables_3)
     endif   
     if (option_out_var(4) .eq. 1) then
        do tycount=0,6*count_mon
            write(unit_monitor4,"(E16.7,A1)",advance='NO') variables_4(tycount)," "
        end do
            write(unit_monitor4,"(A1)") " "
        close(unit_monitor4)
        deallocate(variables_4)
     endif   
     if (option_out_var(5) .eq. 1) then
        do tycount=0,6*count_mon
            write(unit_monitor5,"(E16.7,A1)",advance='NO') variables_5(tycount)," "
        end do
            write(unit_monitor5,"(A1)") " "
        close(unit_monitor5)
        deallocate(variables_5)
     endif   
     if (option_out_var(6) .eq. 1) then
        do tycount=0,3*count_mon
            write(unit_monitor6,"(E16.7,A1)",advance='NO') variables_6(tycount)," "
        end do
            write(unit_monitor5,"(A1)") " "
        close(unit_monitor6)
        deallocate(variables_6)
     endif   
     
    
    
     if (dbg .eq. 1) then
        do tycount=0,count_mon
            write(unit_monitor7,"(E16.7,A1)",advance='NO') variables_7(tycount)," "
        end do
            write(unit_monitor5,"(A1)") " "
        close(unit_monitor7)
        deallocate(variables_7)
        do tycount=0,count_mon
            write(unit_monitor8,"(E16.7,A1)",advance='NO') variables_8(tycount)," "
        end do
            write(unit_monitor5,"(A1)") " "
        close(unit_monitor8)
        deallocate(variables_8)
     endif          

     
        
     end subroutine WRITE_OUTPUT
