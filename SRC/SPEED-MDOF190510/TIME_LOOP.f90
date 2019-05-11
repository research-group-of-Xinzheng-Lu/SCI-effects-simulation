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

!> @brief Performs time integration with the leap-frog or
!! the Runge-Kutta scheme.
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in,out] el_new structure for DG interfaces
!> @note Use Hypbrid programming MPI-OpenMP. 

      subroutine TIME_LOOP(el_new)                      
                                                     
      use max_var
      use DGJUMP
      use speed_timeloop
      !!!ty!!!
      use MDOF_MyData
      !!!ty!!!
      implicit none

      include 'SPEED.MPI'

      type(el4loop), dimension(nelem_dg), intent(in) :: el_new 
      
! GENERAL
      character*70 :: file_monitor, file_outU1, file_outU0, file_outV1, filePG
      character*70 :: file_fe0, file_fe1
      integer*4 :: unit_fe0, unit_fe1

      integer*4 :: mm, iface, ifacene, iene, imne, ie_curr, iene_curr, trofa, posit, jshift, &
                   ic, ih, ik, il, it, ipos, ineg, &
                   isnap, its, fn, nn, ip, im, imon, iaz, &
                   i, j, k, is, in, id, istage, itime, &
                   ie, ielem, count_monitor, find_tag
      ! Elapsed time print-out iteration divisor
      integer*4 :: it_divisor = 1000

      integer*4, dimension(:), allocatable :: node_counter


      real*8 :: start_time, final_time, time_in_seconds, &
                start, finish, &
                tt0, tt1, tt2, deltat2, tt_int, &
                GET_FUNC_VALUE

      real*8 :: tol = 1.d-1

 
      real*8, dimension(:), allocatable     :: jumpx, jumpy, jumpz
      real*8, dimension(:), allocatable     :: func_value

! ELEMENT BY ELEMENT
      real*8                                :: rho, lambda, mu, tref

      real*8, dimension(:), allocatable     :: ct,ww
      real*8, dimension(:,:), allocatable   :: dd
      real*8, dimension(:,:,:), allocatable :: ux_el, uy_el, uz_el, &
                                               duxdx_el, duydx_el, duzdx_el, &
                                               duxdy_el, duydy_el, duzdy_el, &
                                               duxdz_el, duydz_el, duzdz_el, &
                                               sxx_el, syy_el, szz_el, &
                                               syz_el, szx_el, sxy_el, &
                                               fx_el, fy_el, fz_el, &
                                               mv_el, &
                                               vx_el, vy_el, vz_el, &
                                               div_el, rotx_el, roty_el, rotz_el, &
                                               strainxx_el, strainyy_el, strainzz_el, &
                                               strainxy_el, strainyz_el, strainzx_el, &
                                               stressxx_el, stressyy_el, stresszz_el, &
                                               stressxy_el, stressyz_el, stresszx_el, &
                                               rho_el, lambda_el, mu_el, gamma_el             

! MPI
     integer*4, dimension(mpi_np)         :: send_length, recv_length, small_send_length, small_recv_length
     integer*4, dimension(mpi_np)         :: send_length_jump, recv_length_jump, double_send_length, double_recv_length

     integer*4, dimension(:), allocatable :: small_send_buffer,small_recv_buffer

     real*8, dimension(:), allocatable    :: u_1, u0, u1, u2, v1, v2, jump, jump_minus, &
                                          u_p, u_m, jump_all, &
                                          fk, fe, mv, fkold,&
                                          send_buffer, recv_buffer, &
                                          double_send_buffer, double_recv_buffer, &
                                          send_buffer_jump, recv_buffer_jump

! DAMPING
      real*8                                  :: gamma, QS_nh, QP_nh
      real*8, dimension(1,N_SLS)              :: Y_lambda_nh, Y_mu_nh
      real*8, dimension(:), allocatable       :: mc, mck, damp_term,mckty !mckty for stifness damping
      real*8, dimension(:,:,:), allocatable   :: mc_el, mck_el
      real*8, dimension(:,:), allocatable     :: strain_visc
      real*8, dimension(:,:,:,:), allocatable :: strain_visc_xx,strain_visc_xy,&
                                                 strain_visc_xz,strain_visc_yy,&
                                                 strain_visc_yz,strain_visc_zz


! PEAK GROUND MAP
      real*8, dimension (:,:), allocatable     :: max_u,max_v,max_a, max_o
      
! NON LINEAR ELASTIC MATERIAL
      integer*4                                :: im_nle, which_mat_nle, y_or_n_node
      integer*4, dimension(:), allocatable     :: con_spx_loc_nle
      integer*4, dimension(:), allocatable     :: node_nle_4_mpi     
      integer*4, dimension(:,:,:), allocatable :: yes_or_not_nle_el
      real*8                                   :: Depth_nle
      real*8, dimension(:,:,:), allocatable    :: R_el
      real*8, dimension(:), allocatable        :: Depth_nle_el

! OUTPUT
      real*8, dimension(:), allocatable        :: strain, omega, stress

! RUNGE-KUTTA
      real*8, dimension(:), allocatable        :: A, b,c, u_int, v_int
      real*8, dimension(:,:), allocatable      :: Kappau, Kappav

! DEBUG
      real*8, dimension(:,:,:), allocatable    :: gamma_new
      real*8, dimension(:), allocatable        :: mu_nle, gamma_nle

! INSTABILITY CONTROL
       ! Set to TRUE when instability control is enabled and instability is detected.
      logical                                  :: b_instability_abort
      logical, dimension(mpi_np)               :: b_instability_abort_all
       
       b_instability_abort = .FALSE.

!---------------------------------------------------------------------------

      if(mpi_id .eq. 0) start = MPI_WTIME()

      isnap = 1      
      allocate(func_value(nfunc))            
     
!---------------------------------------------------------------------------
! Monitor LiST                                                                                
!---------------------------------------------------------------------------

      if (nmonitors_lst.ge.1) then                
           
         call MAKE_MONITOR_FILES(nmonitors_lst,el_monitor_lst,local_el_num, nelem_loc,&
                                 count_monitor, file_monitor,monitor_file, mpi_id, &
                                 opt_out_var)

      endif 

       if(nmonitors_pgm .ge. 1) then
          allocate(max_u(nmonitors_pgm,9), max_v(nmonitors_pgm,9), max_a(nmonitors_pgm,9), max_o(nmonitors_pgm,3))
          max_u = 0.d0; max_v = 0.d0; max_a = 0.d0; max_o = 0.d0
       endif 

!---------------------------------------------------------------------------
!     MPI INITIALIZATION  ---  NEW NUMERATION recv & send
!---------------------------------------------------------------------------
      
      allocate(send_buffer(3*nsend));      send_buffer= 0.0d0
      allocate(recv_buffer(3*nrecv));      recv_buffer = 0.0d0
      
      do i = 1,mpi_np
         send_length(i) = 3*proc_send(i);  recv_length(i) = 3*proc_recv(i)
         if(opt_out_var(4) .eq. 1 .or. opt_out_var(5) .eq. 1 .or. damping_type .eq. 2) then  
            double_send_length(i) = 6*proc_send(i); double_recv_length(i) = 6*proc_recv(i)
         endif
      enddo

      do i = 1,nrecv
         in = node_recv(i)
         call GET_INDLOC_FROM_INDGLO(local_node_num, nnod_loc, in, id)
         node_recv(i) = id
      enddo
      
      do i = 1,nsend
         in = node_send(i)
         call GET_INDLOC_FROM_INDGLO(local_node_num, nnod_loc, in, id)
         node_send(i) = id
      enddo   
      
!---------------------------------------------------------------------------
!     MPI INITIALIZATION  ---  NON CONFORMING DG 
!---------------------------------------------------------------------------

    if (nelem_dg_glo .gt. 0) then      

       allocate(send_buffer_jump(3*nsend_jump), recv_buffer_jump(3*nrecv_jump))
       send_buffer_jump = 0.0d0;       recv_buffer_jump = 0.0d0
      
      do i = 1,mpi_np
         send_length_jump(i) = 3*proc_send_jump(i)
         recv_length_jump(i) = 3*proc_recv_jump(i)
      enddo      

      allocate(jump_minus(3*nrecv_jump))      

      do i = 1,nsend_jump
         in = node_send_jump(i)
         call GET_INDLOC_FROM_INDGLO(local_node_num, nnod_loc, in, id)
         node_send_jump(i) = id
      enddo        

    endif  
      
!---------------------------------------------------------------------------
!     UNKNOWNS INITIALIZATION
!---------------------------------------------------------------------------
     
      allocate(fe(3*nnod_loc),fk(3*nnod_loc),mv(3*nnod_loc)); fe = 0.0d0; fk = 0.0d0; mv = 0.0d0
      allocate(fkold(3*nnod_loc));fkold=0.0d0
      if(rk_scheme .eq. 'RUNGEKUTTA') then
        allocate(u1(3*nnod_loc),u2(3*nnod_loc),v1(3*nnod_loc),v2(3*nnod_loc),jump(3*nnod_loc))
        u1 = 0.0d0; u2 = 0.0d0; v1 = 0.0d0; jump = 0.0d0; 
        if(testmode .eq. 1) call GET_INITIAL_CONDITIONS(nnod_loc, u2, u1, v1, nelem_loc, con_spx_loc, con_nnz_loc,&
                                                        nmat,prop_mat, sdeg_mat, local_node_num, &
                                                        xx_spx_loc,yy_spx_loc,zz_spx_loc, mpi_id,deltat)
        u2 = 0.d0                                                
      else 
        allocate(u0(3*nnod_loc),u1(3*nnod_loc),u2(3*nnod_loc),v1(3*nnod_loc),jump(3*nnod_loc))
        u0 = 0.0d0; u1 = 0.0d0; u2 = 0.0d0; v1 = 0.0d0; jump = 0.0d0;
        if(testmode .eq. 1) call GET_INITIAL_CONDITIONS(nnod_loc, u0, u1, v1, nelem_loc, con_spx_loc, con_nnz_loc,&
                                                        nmat, prop_mat, sdeg_mat, local_node_num, &
                                                        xx_spx_loc,yy_spx_loc,zz_spx_loc, mpi_id,deltat) 
 
      
      endif

!---------------------------------------------------------------------------
!     UNKNOWNS INITIALIZATION FOR DAMPING
!---------------------------------------------------------------------------

      if (make_damping_yes_or_not.eq. 1 .and. damping_type .eq. 1) then                        
              allocate(mc(3*nnod_loc), mck(3*nnod_loc),mckty(3*nnod_loc)); mc = 0.d0; mck = 0.d0; mckty=0.d0  
      elseif(damping_type .eq. 2) then  
             allocate(strain_visc(6*nnod_loc,N_SLS)); strain_visc = 0.d0          
      elseif(damping_type .eq. 3) then  
              allocate(mc(3*nnod_loc),damp_term((3*nnod_loc))); mc = 0.d0; damp_term=0.d0;
      endif                                                        

!---------------------------------------------------------------------------
!     UNKNOWNS INITIALIZATION FOR OUTOUT
!---------------------------------------------------------------------------

      if(opt_out_var(4) .eq. 1) then
         allocate(stress(6*nnod_loc)); stress = 0.d0
      endif   
      if(opt_out_var(5) .eq. 1 .or. damping_type .eq. 2) then
         allocate(strain(6*nnod_loc)); strain = 0.d0
      endif
      if(opt_out_var(6) .eq. 1) then
         allocate(omega(3*nnod_loc)); omega = 0.d0
      endif    

!---------------------------------------------------------------------------
!     UNKNOWNS INITIALIZATION FOR DEBUG
!---------------------------------------------------------------------------

      if(debug .eq. 1) then
         allocate(mu_nle(nnod_loc), gamma_nle(nnod_loc))
         mu_nle = 0.d0; gamma_nle = 0.d0
      endif 

!---------------------------------------------------------------------------
!     NODE COUNTER INITIALIZATION FOR OUTPUT
!---------------------------------------------------------------------------

      
      if(opt_out_var(4) .eq. 1 .or. opt_out_var(5) .eq. 1 .or. opt_out_var(6) .eq. 1 &
         .or. damping_type .eq. 2) then

        allocate(node_counter(nnod_loc)) ;  node_counter = 0      
         
        do ie = 1,nelem_loc
           im = con_spx_loc(con_spx_loc(ie -1)); nn = sdeg_mat(im) +1
           do k = 1,nn
              do j = 1,nn
                 do i = 1,nn
                    is = nn*nn*(k -1) +nn*(j -1) +i; in = con_spx_loc(con_spx_loc(ie -1) + is)
                    node_counter(in) = node_counter(in) + 1
                 enddo
              enddo
           enddo
        enddo
      
        allocate(small_recv_buffer(nrecv), small_send_buffer(nsend))
        small_recv_buffer = 0; small_send_buffer = 0

        do i = 1,mpi_np
           small_send_length(i) = proc_send(i); small_recv_length(i) = proc_recv(i)
        enddo
 
        do i = 1,nsend
           small_send_buffer(i) = node_counter(node_send(i))
        enddo

        call EXCHANGE_INTEGER(nsend,small_send_buffer,nrecv,small_recv_buffer,&
                                   mpi_np,small_send_length,small_recv_length,&
                                   mpi_comm,mpi_stat,mpi_ierr,mpi_id)
         
        do i = 1,nrecv
           node_counter(node_recv(i)) = node_counter(node_recv(i)) + small_recv_buffer(i)
        enddo

        deallocate(small_recv_buffer, small_send_buffer)
              
        do i = 1, nnod_loc
           if(node_counter(i) .eq. 0) write(*,'(A)') 'Error in node counter!!!'
        enddo
      endif

     !  if (mpi_id .eq. 0) write(*,*) node_counter
!---------------------------------------------------------------------------
!     STOP & GO <--->  READING BACKUP FILES (If present)
!---------------------------------------------------------------------------

     if(len_trim(bkp_file) .ne. 70) then                                                                                  
        !file_outU0 = bkp_file(1:len_trim(bkp_file)) // '/U0_'
        file_outU1 = bkp_file(1:len_trim(bkp_file)) // '/U1_'
        file_outV1 = bkp_file(1:len_trim(bkp_file)) // '/V1_'   
     else
        file_outU0 = 'U0_'; file_outU1 = 'U1_'; file_outV1 = 'V1_' 
     endif
        
      if (tstart .gt. 0.0d0) then
         !call READ_FILEOUT(file_outU0,initial_snap,mpi_id,3*nnod_loc,u0)
         call READ_FILEOUT(file_outU1,initial_snap,mpi_id,3*nnod_loc,u1)
         call READ_FILEOUT(file_outV1,initial_snap,mpi_id,3*nnod_loc,v1)
      endif

!---------------------------------------------------------------------------
!     MAKING MASS MATRIX
!---------------------------------------------------------------------------

      if (nelem_loc.gt.0) then
      
!$OMP PARALLEL &
!$OMP PRIVATE(ie,im,nn,ct,ww,dd,mv_el,rho_el,lambda_el,mu_el,gamma_el,k,j,i,is,in,iaz)

!$OMP DO    
         do ie = 1, nelem_loc

            im = con_spx_loc(con_spx_loc(ie - 1));  nn = sdeg_mat(im) + 1
            
            allocate(ct(nn), ww(nn), dd(nn,nn), mv_el(nn,nn,nn))
            allocate(rho_el(nn,nn,nn),lambda_el(nn,nn,nn),mu_el(nn,nn,nn),gamma_el(nn,nn,nn))
            call MAKE_LGL_NW(nn,ct,ww,dd)              
      
            rho_el = prop_mat(im,1); lambda_el = prop_mat(im,2)
            mu_el = prop_mat(im,3);  gamma_el = prop_mat(im,4)

            ! < STANDARD CASE:                                 >
            ! < if you DON'T go THROUGH the following loop the >
            ! < mechanical properties are given as usual       >
            ! < BLOCK BY BLOCK                                 >
                          
              if (n_case.gt.0) then                                                        
                if (val_case.eq.tag_mat(im)) then                                        
                        call MAKE_ELTENSOR_FOR_CASES(tag_case, val_case,&                        
                                            nn, rho_el, lambda_el, mu_el, gamma_el,&        
                                            nnod_loc, zs_elev, zs_all, vs_tria, thick, &        
                                            con_nnz_loc, con_spx_loc, ie,&        
                                            sub_tag_all, zz_spx_loc, mpi_id, local_node_num, &
                                            damping_type, QS_nh, QP_nh, &
                                            xx_spx_loc, yy_spx_loc,0)

                endif        
             endif                
                
            call MAKE_MASS_MATRIX(nn,ct,ww,dd,rho_el,&
                     alfa11(ie),alfa12(ie),alfa13(ie),&
                     alfa21(ie),alfa22(ie),alfa23(ie),&
                     alfa31(ie),alfa32(ie),alfa33(ie),&
                     beta11(ie),beta12(ie),beta13(ie),&
                     beta21(ie),beta22(ie),beta23(ie),&
                     beta31(ie),beta32(ie),beta33(ie),&
                     gamma1(ie),gamma2(ie),gamma3(ie),&
                     delta1(ie),delta2(ie),delta3(ie),&
                     mv_el)             
                     
                      
                       
            do k = 1,nn
               do j = 1,nn
                  do i = 1,nn
                     is = nn*nn*(k -1) +nn*(j -1) +i
                     in = con_spx_loc(con_spx_loc(ie -1) + is) 
         
                     iaz = 3*(in -1) +1;  mv(iaz) = mv(iaz) + mv_el(i,j,k)                              
                     iaz = 3*(in -1) +2;  mv(iaz) = mv(iaz) + mv_el(i,j,k)                              
                     iaz = 3*(in -1) +3;  mv(iaz) = mv(iaz) + mv_el(i,j,k)

                  enddo
               enddo
            enddo                  

            deallocate(ct, ww, dd, mv_el, rho_el, lambda_el, mu_el, gamma_el)

         enddo    ! loop on elements
!$OMP END DO
!$OMP END PARALLEL
      endif


      do i = 1,nrecv
         in = node_recv(i)
         recv_buffer(3*(i -1) +1) = mv(3*(in -1) +1)
         recv_buffer(3*(i -1) +2) = mv(3*(in -1) +2)
         recv_buffer(3*(i -1) +3) = mv(3*(in -1) +3)
      enddo

      call EXCHANGE_DOUBLE(3*nrecv,recv_buffer,3*nsend,send_buffer,&
                                mpi_np,recv_length,send_length,&
                                mpi_comm,mpi_stat,mpi_ierr,mpi_id)
        
      do i = 1,nsend
         in = node_send(i)
         mv(3*(in -1) +1) = mv(3*(in -1) +1) +send_buffer(3*(i -1) +1)
         mv(3*(in -1) +2) = mv(3*(in -1) +2) +send_buffer(3*(i -1) +2)
         mv(3*(in -1) +3) = mv(3*(in -1) +3) +send_buffer(3*(i -1) +3)
      enddo
      
!---------------------------------------------------------------------------
!     MAKING DAMPING MATRIX
!---------------------------------------------------------------------------

      if (make_damping_yes_or_not .eq. 1 .and. damping_type .ne. 2) then
              if (nelem_loc .gt. 0) then

!$OMP PARALLEL &
!$OMP PRIVATE(ie,im,nn,ct,ww,dd,rho_el,lambda_el,mu_el,gamma_el,mc_el,mck_el,i,j,k,is,in,iaz)
 
!$OMP DO        
              do ie = 1,nelem_loc                                        
                 
                 im = con_spx_loc(con_spx_loc(ie -1));  nn = sdeg_mat(im) +1
                 
                 allocate(ct(nn),ww(nn),dd(nn,nn))
                 allocate(rho_el(nn,nn,nn),lambda_el(nn,nn,nn),mu_el(nn,nn,nn),gamma_el(nn,nn,nn)) 
                 call MAKE_LGL_NW(nn,ct,ww,dd) 
                      
                 rho_el = prop_mat(im,1); lambda_el = prop_mat(im,2)
                 mu_el = prop_mat(im,3);  gamma_el = prop_mat(im,4)
                                     
                 if (damping_type .eq. 1) allocate(mc_el(nn,nn,nn),mck_el(nn,nn,nn))                                        
                 if (damping_type .eq. 3) allocate(mc_el(nn,nn,nn))                                        
                 
                 
                 ! < STANDARD CASE:                                 >
                 ! < if you DON'T go THROUGH the following loop the >
                 ! < mechanical properties are given as usual       >
                 ! < BLOCK BY BLOCK                                 >
                              
                 if (n_case.gt.0) then                                                        
                        if (val_case.eq.tag_mat(im)) then                                        
                                call MAKE_ELTENSOR_FOR_CASES(tag_case, val_case,&                        
                                        nn, rho_el, lambda_el, mu_el, gamma_el,&        
                                        nnod_loc, zs_elev, zs_all, vs_tria, thick, &        
                                        con_nnz_loc, con_spx_loc, ie,&        
                                        sub_tag_all, zz_spx_loc, mpi_id, local_node_num, &
                                        damping_type, QS_nh, QP_nh, &
                                        xx_spx_loc, yy_spx_loc,0)
                        endif        
                 endif                
                 
                 if (damping_type .eq. 1)  then                  

                    call MAKE_DAMPING_MATRIX(nn,ct,ww,dd,&                        
                                          rho_el,gamma_el,&                
                                          alfa11(ie),alfa12(ie),alfa13(ie),&        
                                          alfa21(ie),alfa22(ie),alfa23(ie),&        
                                          alfa31(ie),alfa32(ie),alfa33(ie),&        
                                          beta11(ie),beta12(ie),beta13(ie),&        
                                          beta21(ie),beta22(ie),beta23(ie),&        
                                          beta31(ie),beta32(ie),beta33(ie),&        
                                          gamma1(ie),gamma2(ie),gamma3(ie),&        
                                          delta1(ie),delta2(ie),delta3(ie),&        
                                          mc_el,mck_el)                  
                    do k = 1,nn                                                        
                       do j = 1,nn                                                
                          do i = 1,nn                                                
                             is = nn*nn*(k -1) +nn*(j -1) +i                        
                             in = con_spx_loc(con_spx_loc(ie -1) + is)    
                                   
                             iaz = 3*(in -1) +1
                             mc(iaz) = mc(iaz) +mc_el(i,j,k)                        
                             mck(iaz) = mck(iaz) +mck_el(i,j,k)                                              
                                              
                             iaz = 3*(in -1) +2                        
                             mc(iaz) = mc(iaz) +mc_el(i,j,k)                
                             mck(iaz) = mck(iaz) +mck_el(i,j,k)                                      
                                              
                             iaz = 3*(in -1) +3                        
                             mc(iaz) = mc(iaz) +mc_el(i,j,k)                
                             mck(iaz) = mck(iaz) +mck_el(i,j,k)        

                          enddo                                                
                       enddo                                                
                    enddo                                                        

                 elseif(damping_type .eq. 3) then 
                 
                    call MAKE_RAYLEIGH_DAMPING_MATRIX(nn,ct,ww,dd,&                        
                                          rho_el,&                
                                          alfa11(ie),alfa12(ie),alfa13(ie),&        
                                          alfa21(ie),alfa22(ie),alfa23(ie),&        
                                          alfa31(ie),alfa32(ie),alfa33(ie),&        
                                          beta11(ie),beta12(ie),beta13(ie),&        
                                          beta21(ie),beta22(ie),beta23(ie),&        
                                          beta31(ie),beta32(ie),beta33(ie),&        
                                          gamma1(ie),gamma2(ie),gamma3(ie),&        
                                          delta1(ie),delta2(ie),delta3(ie),&        
                                          mc_el)  
                
                    do k = 1,nn                                                        
                       do j = 1,nn                                                
                          do i = 1,nn                                                
                             is = nn*nn*(k -1) +nn*(j -1) +i                        
                             in = con_spx_loc(con_spx_loc(ie -1) + is)    
                                   
                             iaz = 3*(in -1) +1; mc(iaz) = mc(iaz) + A0_ray(im)*mc_el(i,j,k);
                             iaz = 3*(in -1) +2; mc(iaz) = mc(iaz) + A0_ray(im)*mc_el(i,j,k);                
                             iaz = 3*(in -1) +3; mc(iaz) = mc(iaz) + A0_ray(im)*mc_el(i,j,k);                

                          enddo                                                
                       enddo                                                
                    enddo                                                        

                endif                              
                                          
                 deallocate(ct,ww,dd,rho_el,lambda_el,mu_el,gamma_el,mc_el)
                 if (damping_type .eq. 1 ) deallocate(mck_el)

              enddo   !loop on lements      
!$OMP END DO
!$OMP END PARALLEL
                                
      endif                        

      if (damping_type .eq. 1 .or. damping_type .eq. 3) then
         do i = 1,nrecv
            in = node_recv(i)
            recv_buffer(3*(i -1) +1) = mc(3*(in -1) +1)
            recv_buffer(3*(i -1) +2) = mc(3*(in -1) +2)
            recv_buffer(3*(i -1) +3) = mc(3*(in -1) +3)
         enddo
      
         call EXCHANGE_DOUBLE(3*nrecv,recv_buffer,3*nsend,send_buffer,&
                                mpi_np,recv_length,send_length,&
                                mpi_comm,mpi_stat,mpi_ierr,mpi_id)
      
         do i = 1,nsend
            in = node_send(i)
            mc(3*(in -1) +1) = mc(3*(in -1) +1) +send_buffer(3*(i -1) +1)
            mc(3*(in -1) +2) = mc(3*(in -1) +2) +send_buffer(3*(i -1) +2)
            mc(3*(in -1) +3) = mc(3*(in -1) +3) +send_buffer(3*(i -1) +3)
         enddo
      
      endif
      
      if (damping_type .eq. 1) then
      
         do i = 1,nrecv
            in = node_recv(i)
            recv_buffer(3*(i -1) +1) = mck(3*(in -1) +1)
            recv_buffer(3*(i -1) +2) = mck(3*(in -1) +2)
            recv_buffer(3*(i -1) +3) = mck(3*(in -1) +3)
         enddo
      
         call EXCHANGE_DOUBLE(3*nrecv,recv_buffer,3*nsend,send_buffer,&
                                mpi_np,recv_length,send_length,&
                                mpi_comm,mpi_stat,mpi_ierr,mpi_id)
      
         do i = 1,nsend
            in = node_send(i)
            mck(3*(in -1) +1) = mck(3*(in -1) +1) +send_buffer(3*(i -1) +1)
            mck(3*(in -1) +2) = mck(3*(in -1) +2) +send_buffer(3*(i -1) +2)
            mck(3*(in -1) +3) = mck(3*(in -1) +3) +send_buffer(3*(i -1) +3)
         enddo
         
      endif
!---------------------------------------------------------------------------
   endif    

!---------------------------------------------------------------------------
!     SET VARAIBLES FOR LOW STORAGE RUNGE-KUTTA SCHEME

      if (rk_scheme .eq. 'RUNGEKUTTA') then 

        allocate(A(rk_stages), B(rk_stages), c(rk_stages)) 
        A=0.d0; B=0.d0; c=0.d0
        call MAKE_BUTCHERARRAY(rk_order, rk_stages, A,B,c)

      else  !LEAP-FROG
        allocate(A(1),b(1),c(1))
        rk_stages=1; rk_order=0;
        A=0.d0 ; b=0.d0; c=0.d0;
      endif
     
      allocate(u_int(3*nnod_loc), v_int(3*nnod_loc)); u_int = 0.d0; v_int = 0.d0;          

!---------------------------------------------------------------------------
!     SET VARAIBLES FOR NON-LINEAR ELASTIC FORCES
!     NEW NON LINEAR IMPLEMETATION
      if (nmat_nle .gt. 0)  then
          allocate(con_spx_loc_nle(0:con_nnz_loc));con_spx_loc_nle = con_spx_loc
          allocate(Depth_nle_el(nelem_loc)); Depth_nle_el = 0.d0;
          allocate(node_nle_4_mpi(1:nnod_loc)); node_nle_4_mpi = 0;
                                                  
        do ie = 1,nelem_loc
           im = con_spx_loc(con_spx_loc(ie -1)); nn = sdeg_mat(im) +1

           
            which_mat_nle = 0                                                                        
            do im_nle = 1,nmat_nle                                                                        
               if (tag_mat(im).eq.tag_mat_nle(im_nle)) then                                        
                 which_mat_nle = im_nle                                                        
                 Depth_nle = val_mat_nle(im_nle,1)                                        
               endif                                                                                
            enddo                               
            con_spx_loc_nle(con_spx_loc_nle(ie -1)) = which_mat_nle
            Depth_nle_el(ie) = Depth_nle;
            if (which_mat_nle .eq. 0) con_spx_loc_nle(con_spx_loc_nle(ie-1)+1:con_spx_loc_nle(ie)-1) = 0;
            if (which_mat_nle .ne. 0) then                                                        
                do k = 1,nn                                                                
                   do j = 1,nn                                                        
                       do i = 1,nn                                                
                          is = nn*nn*(k -1) +nn*(j -1) +i                        
                          in = con_spx_loc(con_spx_loc(ie -1) + is)                                
                          con_spx_loc_nle(con_spx_loc_nle(ie -1) + is) = 1
                          node_nle_4_mpi(in) = 1;
                          if (n_case .gt. 0 .and. tag_case .ne. 16 .and. tag_case .ne. 21) then        
                                if (Depth_nle_el(ie) .lt. zs_elev(in) .or. zs_elev(in) .lt. 0.d0) then
                                             con_spx_loc_nle(con_spx_loc_nle(ie -1) + is) = 0; node_nle_4_mpi(in) = 0
                                endif             
                                if ((tag_case .gt. 0).and.(zs_all(in) .lt. 0.0d0)) then
                                             con_spx_loc_nle(con_spx_loc_nle(ie -1) + is) = 0; node_nle_4_mpi(in) = 0
                                endif
                              
                          elseif(n_case .gt. 0 .and. tag_case .eq. 16) then
                                if (Depth_nle_el(ie) .lt. zs_elev(in) .or. zs_elev(in) .lt. 0.d0) then
                                             con_spx_loc_nle(con_spx_loc_nle(ie -1) + is) = 0; node_nle_4_mpi(in) = 0
                                endif
                                if (vs_tria(in) .ge. 450.d0)  then
                                             con_spx_loc_nle(con_spx_loc_nle(ie -1) + is) = 0; node_nle_4_mpi(in) = 0
                                endif             
                         
                          elseif(n_case .gt. 0 .and. tag_case .eq. 21) then
                                if (Depth_nle_el(ie) .lt. zs_elev(in) .or. zs_elev(in) .lt. 0.d0) then
                                             con_spx_loc_nle(con_spx_loc_nle(ie -1) + is) = 0; node_nle_4_mpi(in) = 0
                                endif
                                if (vs_tria(in) .ge. 600.d0 .or. (zs_all(in) .lt. 0.0d0)) then
                                             con_spx_loc_nle(con_spx_loc_nle(ie -1) + is) = 0; node_nle_4_mpi(in) = 0
                                endif             
                          endif
                       enddo                                                                                
                    enddo                                                                                        
                enddo                                                                                                
          endif             
        enddo
      endif

!      write(*,*) 'nle', con_spx_loc_nle
        !!!!!!!!!!!!!!!!!!ty!!!!!!!!!!!!!!!
    !    write(mdoftest2,*) mpi_id
    !    write(mdoftest,*) 'tytest',trim(adjustl(mdoftest2)),'.txt'
    !    mdoftest=TRIM(ADJUSTL(mdoftest))
        !!!!!!!!!!!!!!!!!ty!!!!!!!!!!!!!!!!!
        
!---------------------------------------------------------------------------
!     Make MCKTY
!---------------------------------------------------------------------------
            do id = 1, 3*nnod_loc 
               mckty(id)=mc(id)/mv(id)/(64.d0*datan(1.d0)*datan(1.d0)*(fmax/3.d0)*(fmax/3.d0))  
               !This is not accurate when different material share the same node
               !because we should use stiffness as the weight for stiffness matrix instead of mass
            enddo

!---------------------------------------------------------------------------
!     BEGINNING OF THE TIME LOOP
!---------------------------------------------------------------------------
      !!!ty
      allocate(utyty1(3*count_monitor))
      allocate(utyty2(3*count_monitor))
      allocate(utyty3(3*count_monitor))
      allocate(fty1(3*count_monitor),fty2(3*count_monitor))
      utyty1=0
      utyty2=0
      utyty3=0
      fty1=0
      fty2=0
      !!!ty


      itime = 1 !CONVERGENCE TEST 
          

      if(mpi_id .eq. 0)  write(*,'(A)') '--------------Starting the loop------------------------'
      
      tt0 = tstart - dfloat(1)*deltat; tt1 = tstart; tt2 = tstart + dfloat(1) * deltat
      deltat = tt1 - tt0;  deltat2 = deltat*deltat;
      its = 0;          
      !!!!!!!!!!!!!ty!!!!!!!!!!!!!!!!!!!!!!!!!
      call GenBFile
      !!!!!!!!!!!!!ty!!!!!!!!!!!!!!!!!!!!!!!!!
      do its = 0,nts
          
         if (opt_out_var(4) .eq. 1) stress = 0.d0
         if (opt_out_var(5) .eq. 1 .or. damping_type .eq. 2) strain = 0.d0
         if (opt_out_var(6) .eq. 1) omega = 0.d0

         if (mpi_id.eq.0) write(*,'(A,E14.6)')'TIME = ',tt1 
         if (its.le. 200) then
            start_time = MPI_WTIME()
         elseif (mod(its, it_divisor) .eq. 0) then
            ! Print elapsed time since 200-th iteration ivery it_divisor iterations
            if (mpi_id .eq. 0) then
              final_time = MPI_WTIME()
              write(*,'(A5,I9,A16,F12.3)') '[Iter ', its, '] Elapsed time: ', final_time - start_time
            endif
         endif

         u_int = 0.d0; v_int = 0.d0;
         !!!!!!!!!!!!!!!!!!ty!!!!!!!!!!!!!!!!!

         !!!!!!!!!!!!!!!!!!ty!!!!!!!!!!!!!!!!!    
!---------------------------------------------------------------------------
!     LOOP FOR RK SCHEMES
!---------------------------------------------------------------------------
        ! write(*,*) mpi_id, nnod_loc
         do istage = 1, rk_stages    

            tt_int = tt1 + c(istage)*deltat
                              
            fk = 0.0d0; jump = 0.d0; fe = 0.0d0; 
            if (damping_type .eq. 3) damp_term = 0.d0;

            
            do id = 1, nnod_loc
            
               do fn = 1,nfunc 
                  if(nload_traX_el .gt. 0) then !TRAVELING POINT LOAD X
                     call FIND_POS(nload_traX_el,fun_traX_el,tag_func(fn),find_tag)
                     do i = 1, nb_tra_load
                       if(local_node_num(id) .eq. node_tra(i) .and. find_tag .ne. 0) then
                           iaz = 3*(id -1) +1; fe(iaz) = fe(iaz) + Fel(fn,(3*(id -1) +1),1) * & 
                                               GET_FUNC_VALUE(nfunc,func_type,func_indx,func_data,&
                                                              nfunc_data,fn,tt_int,dist_tra(i),val_traX_el(find_tag,1), &
                                                              MDOFMnum,MDOFforceinput,Fel(fn,(3*(id -1) +1),2), &
                                                              Fel(fn,(3*(id -1) +1),3),its)
                       endif
                     enddo
                  elseif(nload_traY_el .gt. 0) then  !TRAVELING POINT LOAD X    
                     call FIND_POS(nload_traY_el,fun_traY_el,tag_func(fn),find_tag)
                     do i = 1, nb_tra_load
                       if(local_node_num(id) .eq. node_tra(i) .and. find_tag .ne. 0) then
                          iaz = 3*(id -1) +2; fe(iaz) = fe(iaz) + Fel(fn,(3*(id -1) +2),1) * & 
                                               GET_FUNC_VALUE(nfunc,func_type,func_indx,func_data,&
                                                              nfunc_data,fn,tt_int,dist_tra(i),val_traY_el(find_tag,1), &
                                                              MDOFMnum,MDOFforceinput,Fel(fn,(3*(id -1) +2),2), &
                                                              Fel(fn,(3*(id -1) +2),3),its)
                       endif
                     enddo
                  elseif(nload_traZ_el .gt. 0) then  !TRAVELING POINT LOAD Z 
                     call FIND_POS(nload_traZ_el,fun_traZ_el,tag_func(fn),find_tag)
                     do i = 1, nb_tra_load
                       if(local_node_num(id) .eq. node_tra(i) .and. find_tag .ne. 0) then
                          iaz = 3*(id -1) +3; fe(iaz) = fe(iaz) + Fel(fn,(3*(id -1) +3),1) * & 
                                               GET_FUNC_VALUE(nfunc,func_type,func_indx,func_data,&
                                                              nfunc_data,fn,tt_int,dist_tra(i),val_traZ_el(find_tag,1), &
                                                              MDOFMnum,MDOFforceinput,Fel(fn,(3*(id -1) +3),2), &
                                                              Fel(fn,(3*(id -1) +3),3),its)
                       endif
                     enddo
                       

                  else  !EXTERNAL FORCES   
                !  write(*,*) "ty",id,Fel(fn,(3*(id -1) +1),1:3)
                    iaz = 3*(id -1) +1; fe(iaz) = fe(iaz) + Fel(fn,(3*(id -1) +1),1) * & 
                                 GET_FUNC_VALUE(nfunc,func_type,func_indx,func_data,nfunc_data,fn,tt_int,0,0, &
                                                              MDOFMnum,MDOFforceinput,Fel(fn,(3*(id -1) +1),2), &
                                                              Fel(fn,(3*(id -1) +1),3),its)
                    iaz = 3*(id -1) +2; fe(iaz) = fe(iaz) + Fel(fn,(3*(id -1) +2),1) * & 
                                 GET_FUNC_VALUE(nfunc,func_type,func_indx,func_data,nfunc_data,fn,tt_int,0,0, &
                                                              MDOFMnum,MDOFforceinput,Fel(fn,(3*(id -1) +2),2), &
                                                              Fel(fn,(3*(id -1) +2),3),its)
                    iaz = 3*(id -1) +3; fe(iaz) = fe(iaz) + Fel(fn,(3*(id -1) +3),1) * & 
                                 GET_FUNC_VALUE(nfunc,func_type,func_indx,func_data,nfunc_data,fn,tt_int,0,0, &
                                                              MDOFMnum,MDOFforceinput,Fel(fn,(3*(id -1) +3),2), &
                                                              Fel(fn,(3*(id -1) +3),3),its)
                  endif
               enddo
            enddo

            do i = 1,nrecv
               in = node_recv(i)
               recv_buffer(3*(i -1) +1) = fe(3*(in -1) +1)
               recv_buffer(3*(i -1) +2) = fe(3*(in -1) +2)
               recv_buffer(3*(i -1) +3) = fe(3*(in -1) +3)
            enddo

            call EXCHANGE_DOUBLE(3*nrecv,recv_buffer,3*nsend,send_buffer,&
                                mpi_np,recv_length,send_length,&
                                mpi_comm,mpi_stat,mpi_ierr,mpi_id)
        
            do i = 1,nsend
               in = node_send(i)
               !!!!!!!!!!???如果一个节点同时属于不同的core，那么它上面的节点力会被叠加（内力可以，但是这是外力），这样应该是错误的？
               fe(3*(in -1) +1) =fe(3*(in -1) +1) +send_buffer(3*(i -1) +1)
               fe(3*(in -1) +2) =fe(3*(in -1) +2) +send_buffer(3*(i -1) +2)
               fe(3*(in -1) +3) =fe(3*(in -1) +3) +send_buffer(3*(i -1) +3)
            enddo
            
          
            
                          
!---------------------------------------------------------------------------
!     EXCHANGE U1 & V1 VALUES
!---------------------------------------------------------------------------
            do i = 1,nsend
               in = node_send(i)
               send_buffer(3*(i -1) +1) = u1(3*(in -1) +1)
               send_buffer(3*(i -1) +2) = u1(3*(in -1) +2)
               send_buffer(3*(i -1) +3) = u1(3*(in -1) +3)
            enddo
                
            call EXCHANGE_DOUBLE(3*nsend,send_buffer,3*nrecv,recv_buffer,&
                                 mpi_np,send_length,recv_length,&
                                 mpi_comm,mpi_stat,mpi_ierr,mpi_id)
                  
            do i = 1,nrecv
               in = node_recv(i)
               u1(3*(in -1) +1) = recv_buffer(3*(i -1) +1) 
               u1(3*(in -1) +2) = recv_buffer(3*(i -1) +2)
               u1(3*(in -1) +3) = recv_buffer(3*(i -1) +3)
            enddo
         
            do i = 1,nsend
               in = node_send(i)
               send_buffer(3*(i -1) +1) = v1(3*(in -1) +1)
               send_buffer(3*(i -1) +2) = v1(3*(in -1) +2)
               send_buffer(3*(i -1) +3) = v1(3*(in -1) +3)
            enddo
                
            call EXCHANGE_DOUBLE(3*nsend,send_buffer,3*nrecv,recv_buffer,&
                                 mpi_np,send_length,recv_length,&
                                 mpi_comm,mpi_stat,mpi_ierr,mpi_id)
                  
            do i = 1,nrecv
               in = node_recv(i)
               v1(3*(in -1) +1) = recv_buffer(3*(i -1) +1)
               v1(3*(in -1) +2) = recv_buffer(3*(i -1) +2)
               v1(3*(in -1) +3) = recv_buffer(3*(i -1) +3)
            enddo
            
            

!---------------------------------------------------------------------------
!     MAKE INTERNAL FORCES
!---------------------------------------------------------------------------

         

         if (nelem_loc.gt.0) then

 !     write(*,*) 'nmat_nle', nmat_nle

!$OMP PARALLEL &
!$OMP PRIVATE(ie,im,nn,ct,ww,dd,ux_el,uy_el,uz_el) &
!$OMP PRIVATE(duxdx_el,duydx_el,duzdx_el,duxdy_el,duydy_el,duzdy_el,duxdz_el,duydz_el,duzdz_el) &
!$OMP PRIVATE(sxx_el,syy_el,szz_el,syz_el,szx_el,sxy_el,fx_el,fy_el,fz_el,rho_el,lambda_el) &
!$OMP PRIVATE(mu_el,gamma_el,mc_el,mck_el,R_el,yes_or_not_nle_el) &
!$OMP PRIVATE(im_nle,which_mat_nle,Depth_nle,gamma_new) &
!$OMP PRIVATE(k,j,i,is,in,iaz,strain_visc_xx,strain_visc_xy) &
!$OMP PRIVATE(strain_visc_xz,strain_visc_yy,strain_visc_yz,strain_visc_zz,QS_nh,QP_nh) &
!$OMP PRIVATE(Y_lambda_nh,Y_mu_nh)
 
!$OMP DO

            do ie = 1,nelem_loc
                  
                  im = con_spx_loc(con_spx_loc(ie -1)); nn = sdeg_mat(im) +1

                  allocate(ct(nn),ww(nn),dd(nn,nn))
                  allocate(ux_el(nn,nn,nn),uy_el(nn,nn,nn),uz_el(nn,nn,nn))
                  allocate(duxdx_el(nn,nn,nn),duydx_el(nn,nn,nn),duzdx_el(nn,nn,nn)) 
                  allocate(duxdy_el(nn,nn,nn),duydy_el(nn,nn,nn),duzdy_el(nn,nn,nn)) 
                  allocate(duxdz_el(nn,nn,nn),duydz_el(nn,nn,nn),duzdz_el(nn,nn,nn)) 
                  allocate(sxx_el(nn,nn,nn),syy_el(nn,nn,nn),szz_el(nn,nn,nn))
                  allocate(syz_el(nn,nn,nn),szx_el(nn,nn,nn),sxy_el(nn,nn,nn))
                  allocate(fx_el(nn,nn,nn),fy_el(nn,nn,nn),fz_el(nn,nn,nn))
                  allocate(rho_el(nn,nn,nn),lambda_el(nn,nn,nn),mu_el(nn,nn,nn),gamma_el(nn,nn,nn)) 
                  call MAKE_LGL_NW(nn,ct,ww,dd)

                  if (make_damping_yes_or_not .eq. 1 .and. damping_type .eq. 1) allocate(mc_el(nn,nn,nn),mck_el(nn,nn,nn))
                  
                  if (damping_type .eq. 2) & 
                    allocate(strain_visc_xx(nn,nn,nn,N_SLS),strain_visc_xy(nn,nn,nn,N_SLS),&
                             strain_visc_xz(nn,nn,nn,N_SLS),strain_visc_yy(nn,nn,nn,N_SLS),&
                             strain_visc_yz(nn,nn,nn,N_SLS),strain_visc_zz(nn,nn,nn,N_SLS))
                             
                  if (nmat_nle .gt. 0) then
                          allocate(R_el(nn,nn,nn), yes_or_not_nle_el(nn,nn,nn),gamma_new(nn,nn,nn))        
                          R_el = 0; yes_or_not_nle_el = 0; gamma_new = 0;
                          
                          !write(*,*) 'con', 'ie', con_spx_loc_nle(con_spx_loc_nle(ie -1)) 
                          !write(*,*) 'con', 'ie', con_spx_loc_nle(con_spx_loc_nle(ie)) 
                          !read(*,*)
                          if (con_spx_loc_nle(con_spx_loc_nle(ie -1)) .ne. 0) then                                                        
                              do k = 1,nn                                                                
                                 do j = 1,nn                                                        
                                    do i = 1,nn                                                
                                       is = nn*nn*(k -1) +nn*(j -1) +i                
                                       !write(*,*) con_spx_loc_nle(con_spx_loc_nle(ie -1) + is) 
                                       !read(*,*)        
                                       yes_or_not_nle_el(i,j,k) = con_spx_loc_nle(con_spx_loc_nle(ie -1) + is)   
                                       in = con_spx_loc(con_spx_loc(ie -1) + is)

                                       if (yes_or_not_nle_el(i,j,k) .eq. 1) then
                                           iaz = 3*(in -1) +1; mc(iaz) = 0.0d0; mck(iaz) = 0.0d0
                                           iaz = 3*(in -1) +2; mc(iaz) = 0.0d0; mck(iaz) = 0.0d0            
                                           iaz = 3*(in -1) +3; mc(iaz) = 0.0d0; mck(iaz) = 0.0d0 
                                       endif  
                               
                                       if (debug .eq. 1) then
                                          mu_nle(in) = 0.0d0; gamma_nle(in) = 0.0d0
                                       endif                                     
                                    enddo                                                                                
                                 enddo                                                                                        
                             enddo  
                          endif
                   endif       
                  ! write(*,*) nmat_nle,'ie',ie, yes_or_not_nle_el
                  ! read(*,*)
                   
                                   
                   !|
                   !+----------------------------------------------------------------------------------------
                   rho_el = prop_mat(im,1); lambda_el = prop_mat(im,2)                                        
                   mu_el = prop_mat(im,3);  gamma_el = prop_mat(im,4)                                        
                     
                   do k = 1,nn
                      do j = 1,nn
                         do i = 1,nn
                            is = nn*nn*(k -1) +nn*(j -1) +i
                            in = con_spx_loc(con_spx_loc(ie -1) + is)
    
                            if (damping_type .eq. 3) then 
    
                              iaz = 3*(in -1) +1; ux_el(i,j,k) = u1(iaz) + A1_ray(im)*v1(iaz)
                              iaz = 3*(in -1) +2; uy_el(i,j,k) = u1(iaz) + A1_ray(im)*v1(iaz)
                              iaz = 3*(in -1) +3; uz_el(i,j,k) = u1(iaz) + A1_ray(im)*v1(iaz)
                            else
                              iaz = 3*(in -1) +1; ux_el(i,j,k) = u1(iaz)
                              iaz = 3*(in -1) +2; uy_el(i,j,k) = u1(iaz)
                              iaz = 3*(in -1) +3; uz_el(i,j,k) = u1(iaz)
                            endif
                            
                            if(damping_type .eq. 2) then 
                               strain_visc_xx(i,j,k,:) = strain_visc(6*(in -1) +1,:)
                               strain_visc_yy(i,j,k,:) = strain_visc(6*(in -1) +2,:)
                               strain_visc_zz(i,j,k,:) = strain_visc(6*(in -1) +3,:)
                               strain_visc_xy(i,j,k,:) = strain_visc(6*(in -1) +4,:)
                               strain_visc_yz(i,j,k,:) = strain_visc(6*(in -1) +5,:)
                               strain_visc_xz(i,j,k,:) = strain_visc(6*(in -1) +6,:)
                            endif   
                          enddo
                       enddo
                    enddo

                    call MAKE_STRAIN_TENSOR(nn,ct,ww,dd,&                                
                                               alfa11(ie),alfa12(ie),alfa13(ie),&        
                                               alfa21(ie),alfa22(ie),alfa23(ie),&        
                                               alfa31(ie),alfa32(ie),alfa33(ie),&        
                                               beta11(ie),beta12(ie),beta13(ie),&        
                                               beta21(ie),beta22(ie),beta23(ie),&        
                                               beta31(ie),beta32(ie),beta33(ie),&        
                                               gamma1(ie),gamma2(ie),gamma3(ie),&        
                                               delta1(ie),delta2(ie),delta3(ie),&        
                                               ux_el,uy_el,uz_el,&                                
                                               duxdx_el,duydx_el,duzdx_el,&                        
                                               duxdy_el,duydy_el,duzdy_el,&                        
                                               duxdz_el,duydz_el,duzdz_el)
                                                                                                    
                    ! < STANDARD CASE:                                 >
                    ! < if you DON'T go THROUGH the following loop the >
                    ! < mechanical properties are given as usual       >
                    ! < BLOCK BY BLOCK                                 >
                                        
                    if (n_case .gt. 0) then                                                        
                        if (val_case .eq. tag_mat(im)) then                                        
                                call MAKE_ELTENSOR_FOR_CASES(tag_case, val_case,&                        
                                                nn, rho_el, lambda_el, mu_el, gamma_el,&        
                                                nnod_loc, zs_elev, zs_all, vs_tria, thick, &        
                                                con_nnz_loc, con_spx_loc, ie,&        
                                                sub_tag_all, zz_spx_loc, mpi_id, local_node_num,&
                                                damping_type, QS_nh, QP_nh,&
                                                xx_spx_loc, yy_spx_loc,0)
                                                
                         endif        
                    endif                

                    !+----------------------------------------------------------------------------------------
                    !|
                    ! NON LINEAR ELASTIC IMPLEMENTATION:
                    ! MODIFYING lambda_el, mu_el, gamma_el
        
!                   if (which_mat_nle.ne.0 .and. damping_type .eq. 1) then    
                    if (nmat_nle .gt. 0 .and. damping_type .eq. 1) then 
                          
                        which_mat_nle = con_spx_loc_nle(con_spx_loc_nle(ie-1))       
                       ! write(*,*) which_mat_nle,ie
                       ! read(*,*)
                           
                        if (which_mat_nle .ne. 0) then                                                
                                                                                                                                              
                        ! Invariants of the strain tensor sigma_ij
                        ! & main strain computations
                   
                           call MAKE_INVARIANTS_AND_MAIN_STRAIN(nn,&                        
                                                 duxdx_el,duydx_el,duzdx_el,&                                
                                                 duxdy_el,duydy_el,duzdy_el,&                                
                                                 duxdz_el,duydz_el,duzdz_el,&                                
                                                 R_el,yes_or_not_nle_el)                                
                                                                                                  
                           call MAKE_ELTENSOR_FOR_CASES_NLE(prop_mat_nle(which_mat_nle,1),R_el, &  
                                                 nn,rho_el,lambda_el,mu_el,gamma_el,&
                                                 con_nnz_loc,con_spx_loc,ie,&
                                                 func_type,func_indx,func_data,nfunc_data,&
                                                 nfunc,tt_int,tag_func,yes_or_not_nle_el,tag_case,&
                                                 nnod_loc,vs_tria)
                                                 
                                                                                                                                          
                           call MAKE_DAMPING_MATRIX_NLE(nn,ct,ww,dd,&                        
                                                   rho_el,gamma_el,&                                        
                                                   alfa11(ie),alfa12(ie),alfa13(ie),&                
                                                   alfa21(ie),alfa22(ie),alfa23(ie),&                
                                                   alfa31(ie),alfa32(ie),alfa33(ie),&                
                                                   beta11(ie),beta12(ie),beta13(ie),&                
                                                   beta21(ie),beta22(ie),beta23(ie),&                
                                                   beta31(ie),beta32(ie),beta33(ie),&                
                                                   gamma1(ie),gamma2(ie),gamma3(ie),&                
                                                   delta1(ie),delta2(ie),delta3(ie),&                
                                                   mc_el,mck_el,&                                        
                                                   con_nnz_loc,con_spx_loc,ie,&        
                                                   prop_mat_nle(which_mat_nle,1),R_el,fpeak, &
                                                   func_type,func_indx,func_data,nfunc_data,&
                                                   nfunc,tt_int,tag_func,yes_or_not_nle_el,gamma_new, &
                                                   tag_case,nnod_loc,vs_tria)
                                 
                                 
                            do k = 1,nn                                                
                               do j = 1,nn                                        
                                  do i = 1,nn                                
                                                                   
                                     if (yes_or_not_nle_el(i,j,k) .eq. 1) then
                                                           
                                        is = nn*nn*(k -1) +nn*(j -1) +i                
                                        in = con_spx_loc(con_spx_loc(ie -1) + is)                        
                                
                                        iaz = 3*(in -1) +1                        
                                        mc(iaz) = mc(iaz) +mc_el(i,j,k)        
                                        mck(iaz) = mck(iaz) +mck_el(i,j,k)  
                                        if (debug .eq. 1) mu_nle(in) = mu_el(i,j,k)
                                        if (debug .eq. 1) gamma_nle(in) = gamma_new(i,j,k)
                                                 
                                        iaz = 3*(in -1) +2                        
                                        mc(iaz) = mc(iaz) +mc_el(i,j,k)        
                                        mck(iaz) = mck(iaz) +mck_el(i,j,k)   
                                        iaz = 3*(in -1) +3                        
                                        mc(iaz) = mc(iaz) + mc_el(i,j,k)        
                                        mck(iaz) = mck(iaz) + mck_el(i,j,k)        
                                                                
                                     endif !if (yes_or_not_nle_el(i,j,k).eq.1) then
                                                                
                                  enddo                                                        
                              enddo                                                                
                            enddo                                                                         
                        endif !which_mat_nle .ne. 0
                    endif  !nmat_nle .ne. 0                                                                                  
                    !|
                    !+----------------------------------------------------------------------------------------
                                                                         
                    !+----------------------------------------------------------------------------------------
                    !|
                    ! HERE VISCOPLASTIC IMPLEMENTATION:
                    ! MODIFYING lambda_el and mu_el
                    !|
                    !+----------------------------------------------------------------------------------------

                    if (damping_type .eq. 2) then 


                       if (n_case.gt.0) then                                                        
                           if (val_case.eq.tag_mat(im)) then 
                              call MAKE_ANELASTIC_COEFFICIENTS_NH(nn, N_SLS, lambda_el, mu_el,&
                                            rho_el, QS_nh, QP_nh, fmax, &
                                            Y_lambda_nh, Y_mu_nh, frequency_range, mpi_id)
                           endif
                       else
                           Y_lambda_nh(1,:) = Y_lambda(im,:);
                           Y_mu_nh(1,:) = Y_mu(im,:);
                       
                       endif
                           
                      
                       call MAKE_STRESS_TENSOR_DAMPED(nn,lambda_el,mu_el,&                                        
                                   duxdx_el,duydx_el,duzdx_el,&                                
                                   duxdy_el,duydy_el,duzdy_el,&                                
                                   duxdz_el,duydz_el,duzdz_el,&
                                   strain_visc_xx, strain_visc_xy, strain_visc_xz, &
                                   strain_visc_yy, strain_visc_yz, strain_visc_zz, &
                                   Y_lambda_nh(1,:),Y_mu_nh(1,:), N_SLS, &                                
                                   sxx_el,syy_el,szz_el,&                                
                                   syz_el,szx_el,sxy_el,ie)        
                                   
              
                                   
                        do k = 1,nn
                           do j = 1,nn
                              do i = 1,nn
                                 is = nn*nn*(k -1) +nn*(j -1) +i
                                 in = con_spx_loc(con_spx_loc(ie -1) + is)
                                
                                 iaz = 6*(in -1) +1; strain(iaz) = strain(iaz) + duxdx_el(i,j,k)/node_counter(in) !xx
                                !if (local_el_num(ie) .eq. 453) write(*,*) mpi_id, &
                                !                              local_node_num(in), duxdx_el(i,j,k)
                                 iaz = 6*(in -1) +2; strain(iaz) = strain(iaz) + duydy_el(i,j,k)/node_counter(in) !yy
                                !if (local_el_num(ie) .eq. 453) write(*,*) mpi_id, &
                                !                              local_node_num(in), duydy_el(i,j,k)
                                 iaz = 6*(in -1) +3; strain(iaz) = strain(iaz) + duzdz_el(i,j,k)/node_counter(in) !zz
                                !if (local_el_num(ie) .eq. 453) write(*,*) mpi_id, &
                                !                              local_node_num(in), duzdz_el(i,j,k)
                                 iaz = 6*(in -1) +4; strain(iaz) = strain(iaz) &
                                                                  + 0.5*(duxdy_el(i,j,k) + duydx_el(i,j,k))/node_counter(in) !xy
                                !if (local_el_num(ie) .eq. 453) write(*,*) mpi_id, &
                                !                              local_node_num(in), duxdy_el(i,j,k)
                                 iaz = 6*(in -1) +5; strain(iaz) = strain(iaz)  & 
                                                                  + 0.5*(duydz_el(i,j,k) + duzdy_el(i,j,k))/node_counter(in) !yz
                                !if (local_el_num(ie) .eq. 453) write(*,*) mpi_id, &
                                !                              local_node_num(in), duxdz_el(i,j,k)
                                 iaz = 6*(in -1) +6; strain(iaz) = strain(iaz)  &
                                                                 + 0.5*(duxdz_el(i,j,k) + duzdx_el(i,j,k))/node_counter(in) !zx
                                !if (local_el_num(ie) .eq. 453) write(*,*) mpi_id, &
                                !                              local_node_num(in), duxdz_el(i,j,k)                         
                             enddo
                          enddo
                       enddo
                                   
                                                           
                   
                    else 
                        call MAKE_STRESS_TENSOR(nn,lambda_el,mu_el,&                                        
                                   duxdx_el,duydx_el,duzdx_el,&                                
                                   duxdy_el,duydy_el,duzdy_el,&                                
                                   duxdz_el,duydz_el,duzdz_el,&                                
                                   sxx_el,syy_el,szz_el,&                                
                                   syz_el,szx_el,sxy_el)                                
                    endif
                 
                                  
                 
                    if (nload_sism_el.gt.0) then                                                
                        
                        if(length_check_node_sism .gt. 0) then

                                call MAKE_SEISMIC_MOMENT(nn,ct,ww,dd,&                        
                                                sxx_el,syy_el,szz_el,syz_el,szx_el,sxy_el,&        
                                                nload_sism_el,&                                        
                                                check_node_sism,check_dist_node_sism,&                         
                                                length_check_node_sism,ie,factor_seismic_moment,&                        
                                                tau_seismic_moment,&                                         
                                                func_type,func_indx,func_data,nfunc_data,nfunc,tt_int,&         
                                                con_nnz_loc,con_spx_loc,tag_func, &
                                                nelem_loc, local_el_num, &
                                                nnod_loc, local_node_num)                                  
                        endif
                   endif                                                                
                  
                   if (nload_expl_el.gt.0) then                                                
                        
                        if(length_check_node_expl .gt. 0) then

                                call MAKE_EXPL_SOURCE(nn,ct,ww,dd,&                
                                                sxx_el,syy_el,szz_el,syz_el,szx_el,sxy_el,&        
                                                nload_expl_el,&                                        
                                                check_node_expl,check_dist_node_expl,&                        
                                                length_check_node_expl,ie,factor_explosive_source,&                        
                                                func_type,func_indx,func_data,nfunc_data,nfunc,tt_int,&         
                                                con_nnz_loc,con_spx_loc,tag_func, &
                                                nelem_loc, local_el_num, &
                                                nnod_loc, local_node_num)                                  
                        endif
                   endif                                                                

                 if(opt_out_var(4) .eq. 1) then 
                   do k = 1,nn
                      do j = 1,nn
                         do i = 1,nn
                           is = nn*nn*(k -1) +nn*(j -1) +i
                           in = con_spx_loc(con_spx_loc(ie -1) + is)                         
                           iaz = 6*(in -1) +1; stress(iaz) = stress(iaz) + sxx_el(i,j,k)/node_counter(in) 
                           iaz = 6*(in -1) +2; stress(iaz) = stress(iaz) + syy_el(i,j,k)/node_counter(in) 
                           iaz = 6*(in -1) +3; stress(iaz) = stress(iaz) + szz_el(i,j,k)/node_counter(in) 
                           iaz = 6*(in -1) +4; stress(iaz) = stress(iaz) + sxy_el(i,j,k)/node_counter(in) 
                           iaz = 6*(in -1) +5; stress(iaz) = stress(iaz) + syz_el(i,j,k)/node_counter(in) 
                           iaz = 6*(in -1) +6; stress(iaz) = stress(iaz) + szx_el(i,j,k)/node_counter(in) 
                         enddo
                      enddo
                   enddo             
                 endif
                 if(opt_out_var(5) .eq. 1 .and. damping_type .eq. 1) then
                       do k = 1,nn
                          do j = 1,nn
                             do i = 1,nn
                                is = nn*nn*(k -1) +nn*(j -1) +i
                                in = con_spx_loc(con_spx_loc(ie -1) + is)                              
                                iaz = 6*(in -1) +1; strain(iaz) = strain(iaz) + duxdx_el(i,j,k) /node_counter(in)
                                iaz = 6*(in -1) +2; strain(iaz) = strain(iaz) + duydy_el(i,j,k) /node_counter(in)
                                iaz = 6*(in -1) +3; strain(iaz) = strain(iaz) + duzdz_el(i,j,k) /node_counter(in)
                                iaz = 6*(in -1) +4; strain(iaz) = strain(iaz) + 0.5*(duxdy_el(i,j,k) + duydx_el(i,j,k)) &
                                                                                    /node_counter(in)
                                iaz = 6*(in -1) +5; strain(iaz) = strain(iaz) + 0.5*(duydz_el(i,j,k) + duzdy_el(i,j,k)) &
                                                                                    /node_counter(in)
                                iaz = 6*(in -1) +6; strain(iaz) = strain(iaz) + 0.5*(duxdz_el(i,j,k) + duzdx_el(i,j,k)) &
                                                                                    /node_counter(in)
                            enddo
                         enddo
                      enddo
                 endif           

                 if(opt_out_var(6) .eq. 1)  then
                     do k = 1,nn
                        do j = 1,nn
                           do i = 1,nn
                              is = nn*nn*(k -1) +nn*(j -1) +i
                              in = con_spx_loc(con_spx_loc(ie -1) + is)                              
                              iaz = 3*(in -1) +1; omega(iaz) = omega(iaz) + 0.5*(duxdy_el(i,j,k) - duydx_el(i,j,k)) &
                                                                                    /node_counter(in)
                              iaz = 3*(in -1) +2; omega(iaz) = omega(iaz) + 0.5*(duydz_el(i,j,k) - duzdy_el(i,j,k)) &
                                                                                    /node_counter(in)
                              iaz = 3*(in -1) +3; omega(iaz) = omega(iaz) + 0.5*(duxdz_el(i,j,k) - duzdx_el(i,j,k)) &
                                                                                    /node_counter(in)
                           enddo
                        enddo
                     enddo
                 endif



                   call MAKE_INTERNAL_FORCE(nn,ct,ww,dd,&
                               alfa11(ie),alfa12(ie),alfa13(ie),&
                               alfa21(ie),alfa22(ie),alfa23(ie),&
                               alfa31(ie),alfa32(ie),alfa33(ie),&
                               beta11(ie),beta12(ie),beta13(ie),&
                               beta21(ie),beta22(ie),beta23(ie),&
                               beta31(ie),beta32(ie),beta33(ie),&
                               gamma1(ie),gamma2(ie),gamma3(ie),&
                               delta1(ie),delta2(ie),delta3(ie),&
                               sxx_el,syy_el,szz_el,syz_el,szx_el,sxy_el,&
                               fx_el,fy_el,fz_el)
                               
                 

                    do k = 1,nn
                       do j = 1,nn
                          do i = 1,nn
                             is = nn*nn*(k -1) +nn*(j -1) +i
                             in = con_spx_loc(con_spx_loc(ie -1) + is)
                  
                             iaz = 3*(in -1) +1; fk(iaz) = fk(iaz) +fx_el(i,j,k); 
                             iaz = 3*(in -1) +2; fk(iaz) = fk(iaz) +fy_el(i,j,k)
                             iaz = 3*(in -1) +3; fk(iaz) = fk(iaz) +fz_el(i,j,k)

                             if (damping_type .eq. 3) then 
                                iaz = 3*(in -1) +1; damp_term(iaz) = damp_term(iaz) + ux_el(i,j,k) 
                                iaz = 3*(in -1) +2; damp_term(iaz) = damp_term(iaz) + uy_el(i,j,k) 
                                iaz = 3*(in -1) +3; damp_term(iaz) = damp_term(iaz) + uz_el(i,j,k) 
                             endif
                             
                          enddo
                       enddo
                    enddo
                    
                                
                    deallocate(ct,ww,dd,ux_el,uy_el,uz_el)
                    deallocate(duxdx_el,duydx_el,duzdx_el,duxdy_el,duydy_el,duzdy_el,duxdz_el,duydz_el,duzdz_el) 
                    deallocate(sxx_el,syy_el,szz_el,syz_el,szx_el,sxy_el,fx_el,fy_el,fz_el)
                    deallocate(rho_el,lambda_el,mu_el,gamma_el) 
                    if(damping_type .eq. 2) deallocate(strain_visc_xx,strain_visc_xy,strain_visc_xz,&
                                                       strain_visc_yy,strain_visc_yz,strain_visc_zz)
                        
                    if (make_damping_yes_or_not .eq. 1 .and. damping_type .eq. 1)  deallocate(mc_el, mck_el)  
                    if (nmat_nle .gt. 0)  deallocate(R_el,yes_or_not_nle_el,gamma_new)        

            enddo   !loop on elements
!$OMP END DO
!$OMP END PARALLEL

         endif




!---------------------------------------------------------------------------
!     MAKE INTERNAL FORCES FOR ABC
!---------------------------------------------------------------------------
          
         if (nelem_abc.gt.0) then

!$OMP PARALLEL &
!$OMP PRIVATE(ie,ielem, ie_curr, im, nn, ct,ww,dd, ux_el, uy_el, uz_el, fx_el, fy_el, fz_el) &
!$OMP PRIVATE(vx_el, vy_el, vz_el, rho_el, lambda_el, mu_el, gamma_el) &
!$OMP PRIVATE( k, j, i, is, in, iaz, QS_nh, QP_nh)
 
!$OMP DO

             do ie = 1,nelem_abc
               
               ielem = ielem_abc(ie,1)
               call GET_INDLOC_FROM_INDGLO(local_el_num, nelem_loc, ielem, ie_curr)
                  
               im = con_spx_loc(con_spx_loc(ie_curr -1)); nn = sdeg_mat(im) +1
                      
               allocate(ct(nn),ww(nn),dd(nn,nn))
               allocate(ux_el(nn,nn,nn),uy_el(nn,nn,nn),uz_el(nn,nn,nn))
               allocate(fx_el(nn,nn,nn),fy_el(nn,nn,nn),fz_el(nn,nn,nn))
               allocate(vx_el(nn,nn,nn),vy_el(nn,nn,nn),vz_el(nn,nn,nn))
               allocate(rho_el(nn,nn,nn),lambda_el(nn,nn,nn),mu_el(nn,nn,nn),gamma_el(nn,nn,nn)) 
               
               call MAKE_LGL_NW(nn,ct,ww,dd)

               rho_el = prop_mat(im,1); lambda_el = prop_mat(im,2)                                        
               mu_el = prop_mat(im,3);  gamma_el = prop_mat(im,4)                                        
             
               do k = 1,nn
                  do j = 1,nn
                     do i = 1,nn
                        is = nn*nn*(k -1) +nn*(j -1) +i
                        in = con_spx_loc(con_spx_loc(ie_curr -1) + is)
                              
                        iaz = 3*(in -1) +1
                        ux_el(i,j,k) = u1(iaz) !u1(iaz)
                        vx_el(i,j,k) = v1(iaz) !(3.0d0*u1(iaz) - 4.0d0*u0(iaz) + 1.0d0*u_1(iaz)) / (2.0d0*deltat) 
                        iaz = 3*(in -1) +2
                        uy_el(i,j,k) = u1(iaz) !u1(iaz)
                        vy_el(i,j,k) = v1(iaz) !(3.0d0*u1(iaz) - 4.0d0*u0(iaz) + 1.0d0*u_1(iaz)) / (2.0d0*deltat) 
                        iaz = 3*(in -1) +3
                        uz_el(i,j,k) = u1(iaz) !u1(iaz)
                        vz_el(i,j,k) = v1(iaz) !(3.0d0*u1(iaz) - 4.0d0*u0(iaz) + 1.0d0*u_1(iaz)) / (2.0d0*deltat) 

                     enddo
                  enddo
               enddo
                                         
               ! < STANDARD CASE:                                 >
               ! < if you DON'T go THROUGH the following loop the >
               ! < mechanical properties are given as usual       >
               ! < BLOCK BY BLOCK                                 >
                                         
               if (n_case.gt.0) then                                                        
                  if (val_case.eq.tag_mat(im)) then                                        
                        call MAKE_ELTENSOR_FOR_CASES(tag_case, val_case,&                        
                                                nn, rho_el, lambda_el, mu_el, gamma_el,&        
                                                nnod_loc, zs_elev, zs_all, vs_tria, thick,  &        
                                                con_nnz_loc, con_spx_loc, ie,&        
                                                sub_tag_all, zz_spx_loc, mpi_id, local_node_num, &
                                                damping_type, QS_nh, QP_nh, &
                                                xx_spx_loc, yy_spx_loc,0)

                  endif        
                endif                

                     
                call MAKE_ABC_FORCE(nn,ct,ww,dd,&
                                   rho_el,lambda_el,mu_el,&                         
                                  alfa11(ie_curr),alfa12(ie_curr),alfa13(ie_curr),&
                                  alfa21(ie_curr),alfa22(ie_curr),alfa23(ie_curr),&
                                  alfa31(ie_curr),alfa32(ie_curr),alfa33(ie_curr),&
                                  beta11(ie_curr),beta12(ie_curr),beta13(ie_curr),&
                                  beta21(ie_curr),beta22(ie_curr),beta23(ie_curr),&
                                  beta31(ie_curr),beta32(ie_curr),beta33(ie_curr),&
                                  gamma1(ie_curr),gamma2(ie_curr),gamma3(ie_curr),&
                                  delta1(ie_curr),delta2(ie_curr),delta3(ie_curr),&
                                  ielem_abc(ie,2),ielem_abc(ie,3),&
                                  ielem_abc(ie,4),ielem_abc(ie,5),&
                                  ielem_abc(ie,6),ielem_abc(ie,7),&
                                  ux_el,uy_el,uz_el,vx_el,vy_el,vz_el,&
                                  fx_el,fy_el,fz_el)
                     
                do k = 1,nn
                   do j = 1,nn
                      do i = 1,nn
                         is = nn*nn*(k -1) +nn*(j -1) +i
                         in = con_spx_loc(con_spx_loc(ie_curr -1) + is)

                         iaz = 3*(in -1) +1; fk(iaz) = fk(iaz) +fx_el(i,j,k)
                         iaz = 3*(in -1) +2; fk(iaz) = fk(iaz) +fy_el(i,j,k)
                         iaz = 3*(in -1) +3; fk(iaz) = fk(iaz) +fz_el(i,j,k)
                      enddo
                   enddo
                 enddo
               
                 deallocate(ct,ww,dd,ux_el,uy_el,uz_el)               
                 deallocate(rho_el,lambda_el,mu_el,gamma_el) 
                 deallocate(vx_el,vy_el,vz_el)
                 deallocate(fx_el,fy_el,fz_el)
                 
             enddo   !loop on materials
!$OMP END DO
!$OMP END PARALLEL

         endif

!---------------------------------------------------------------------------
!     MAKE DG JUMPS
!---------------------------------------------------------------------------

      if(nelem_dg_glo .gt. 0) then
         
          jump_minus = 0.d0

          do i = 1,nsend_jump
             in = node_send_jump(i)
                         
             if (testmode .eq. 1) then 
               send_buffer_jump(3*(i -1) +1) = u1(3*(in -1) +1) !+ v1(3*(in -1) +1) 
               send_buffer_jump(3*(i -1) +2) = u1(3*(in -1) +2) !+ v1(3*(in -1) +2) 
               send_buffer_jump(3*(i -1) +3) = u1(3*(in -1) +3) !+ v1(3*(in -1) +3)          
             elseif (damping_type .eq. 3) then
               send_buffer_jump(3*(i -1) +1) = damp_term(3*(in -1) +1)           
               send_buffer_jump(3*(i -1) +2) = damp_term(3*(in -1) +2)
               send_buffer_jump(3*(i -1) +3) = damp_term(3*(in -1) +3)
             else
               send_buffer_jump(3*(i -1) +1) = u1(3*(in -1) +1)           
               send_buffer_jump(3*(i -1) +2) = u1(3*(in -1) +2)
               send_buffer_jump(3*(i -1) +3) = u1(3*(in -1) +3)
             endif
             
          enddo
      
         call EXCHANGE_DOUBLE(3*nsend_jump, send_buffer_jump, 3*nrecv_jump, recv_buffer_jump,&
                                mpi_np, send_length_jump, recv_length_jump,&
                                mpi_comm, mpi_stat, mpi_ierr,mpi_id)
     
          do i = 1,nrecv_jump
             in = node_recv_jump(i)
             jump_minus(3*(i -1) +1) = recv_buffer_jump(3*(in -1) +1)           
             jump_minus(3*(i -1) +2) = recv_buffer_jump(3*(in -1) +2)
             jump_minus(3*(i -1) +3) = recv_buffer_jump(3*(in -1) +3)
          enddo

      
!$OMP PARALLEL &
!$OMP PRIVATE(ie,ielem, nn, ie_curr, u_p, jumpx, jumpy, jumpz, i, j, k, is, in, iaz) &
!$OMP PRIVATE(u_m, jshift, ic, iene, imne, mm, trofa, posit, id, iene_curr, im) &
!$OMP PRIVATE(jump_all)
 
!$OMP DO      

            do ie = 1, nelem_dg

               ielem = el_new(ie)%ind; nn = el_new(ie)%deg
               im = el_new(ie)%mate;
               call GET_INDLOC_FROM_INDGLO(local_el_num, nelem_loc, ielem, ie_curr)

               allocate(u_p(3*nn**3),jumpx(nn**3),jumpy(nn**3),jumpz(nn**3))

               u_p = 0.d0; jumpx = 0.d0; jumpy = 0.d0; jumpz = 0.d0

               do k = 1,nn
                  do j = 1,nn
                     do i = 1,nn
                        is = nn*nn*(k -1) +nn*(j -1) +i
                        in = con_spx_loc(con_spx_loc(ie_curr -1) + is)
                          
                        if (testmode .eq. 1) then 
                           iaz = 3*(in -1) +1; u_p(is) = u1(iaz) !+ v1(iaz)
                           iaz = 3*(in -1) +2; u_p(nn**3+is) = u1(iaz) !+ v1(iaz)
                           iaz = 3*(in -1) +3; u_p(2*nn**3 + is) = u1(iaz) !+ v1(iaz)
                        elseif(damping_type .eq. 3) then 
                           iaz = 3*(in -1) +1; u_p(is) = u1(iaz) + A1_ray(im)*v1(iaz)
                           iaz = 3*(in -1) +2; u_p(nn**3+is) = u1(iaz) + A1_ray(im)*v1(iaz)
                           iaz = 3*(in -1) +3; u_p(2*nn**3 + is) = u1(iaz) + A1_ray(im)*v1(iaz)
                        else
                           iaz = 3*(in -1) +1; u_p(is) = u1(iaz) 
                           iaz = 3*(in -1) +2; u_p(nn**3+is) = u1(iaz) 
                           iaz = 3*(in -1) +3; u_p(2*nn**3 + is) = u1(iaz)
                        endif

                     enddo
                  enddo
                enddo


               allocate(u_m(el_new(ie)%nnz_col)); u_m = 0.d0; jshift = 0 
                
               do ic = 1, el_new(ie)%num_of_ne 
                  
                   iene = el_new(ie)%el_conf(ic,1)       
                   imne = el_new(ie)%el_conf(ic,0)
                       
                   mm = 2
                   do i = 1, nmat 
                      if(tag_mat(i) .eq. imne ) then
                         mm = sdeg_mat(i) +1
                      endif
                   enddo   

                   trofa  = 0 
                   if(mpi_np .gt. 1) call CHECK_MPI(con_nnz_dg, con_spx_dg, iene, trofa, posit)

                   if(trofa .eq. 1) then
                     do k = 1,mm
                        do j = 1,mm
                           do i = 1,mm
                              is =  mm*mm*(k -1) +mm*(j -1) +i
                              id = con_spx_dg(con_spx_dg(posit-1)+is)
                                                                                
                              iaz = 3*(id -1) +1; u_m(is + jshift) = jump_minus(iaz)
                              iaz = 3*(id -1) +2; u_m(mm**3 + is + jshift) = jump_minus(iaz)
                              iaz = 3*(id -1) +3; u_m(2*(mm**3) + is + jshift) = jump_minus(iaz)  
                            enddo
                         enddo   
                      enddo
                   else                                     

                     call GET_INDLOC_FROM_INDGLO(local_el_num, nelem_loc, iene, iene_curr)

                     do k = 1,mm
                        do j = 1,mm
                           do i = 1,mm
                              is =  mm*mm*(k -1) +mm*(j -1) +i
                              id = con_spx_loc(con_spx_loc(iene_curr-1)+is)
     
                              if (testmode .eq. 1) then
                                    iaz = 3*(id -1) +1; u_m(is + jshift) = u1(iaz) !+ v1(iaz) 
                                    iaz = 3*(id -1) +2; u_m(mm**3 + is + jshift) = u1(iaz) !+ v1(iaz)  
                                    iaz = 3*(id -1) +3; u_m(2*(mm**3) + is + jshift) = u1(iaz) !+ v1(iaz)
                              elseif (damping_type .eq. 3) then
                                    iaz = 3*(id -1) +1; u_m(is + jshift) = u1(iaz) + A1_ray(imne)*v1(iaz) 
                                    iaz = 3*(id -1) +2; u_m(mm**3 + is + jshift) = u1(iaz) + A1_ray(imne)*v1(iaz)  
                                    iaz = 3*(id -1) +3; u_m(2*(mm**3) + is + jshift) = u1(iaz) + A1_ray(imne)*v1(iaz)
                              else
                                    iaz = 3*(id -1) +1; u_m(is + jshift) = u1(iaz)
                                    iaz = 3*(id -1) +2; u_m(mm**3 + is + jshift) = u1(iaz)
                                    iaz = 3*(id -1) +3; u_m(2*(mm**3) + is + jshift) = u1(iaz)
                              endif  
                           enddo
                        enddo
                      enddo
                    endif

                   jshift = jshift + 3*mm**3 

                enddo

                allocate(jump_all(3*nn**3))
                jump_all = 0.d0


                call MATMUL_SPARSE(el_new(ie)%matMin, el_new(ie)%nnz_minus, el_new(ie)%JMin, &
                                  el_new(ie)%IMin, jump_all, 3*nn**3,  &
                                  u_m, el_new(ie)%nnz_col,1) 
                
                jumpx = jumpx + jump_all(1:nn**3)
                jumpy = jumpy + jump_all(nn**3+1 : 2*nn**3)
                jumpz = jumpz + jump_all(2*nn**3+1 : 3*nn**3)

                deallocate(u_m, jump_all)
                  
                allocate(jump_all(3*nn**3)); jump_all = 0.d0   

                call MATMUL_SPARSE(el_new(ie)%matPlus, el_new(ie)%nnz_plus,el_new(ie)%JPlus, &
                                  el_new(ie)%IPlus, jump_all, 3*nn**3, &
                                  u_p, 3*nn**3,1) 
                
                
                jumpx = jumpx + jump_all(1:nn**3)
                jumpy = jumpy + jump_all(nn**3+1 : 2*nn**3)
                jumpz = jumpz + jump_all(2*nn**3+1 : 3*nn**3)

                 do k = 1,nn
                    do j = 1,nn
                       do i = 1,nn
                          is = nn*nn*(k -1) +nn*(j -1) +i
                          in = con_spx_loc(con_spx_loc(ie_curr -1) + is)
                           
                          iaz = 3*(in -1) + 1; jump(iaz) = jump(iaz) + jumpx(is)    
                          iaz = 3*(in -1) + 2; jump(iaz) = jump(iaz) + jumpy(is)
                          iaz = 3*(in -1) + 3; jump(iaz) = jump(iaz) + jumpz(is)

                       enddo
                    enddo
                 enddo
                    
                 deallocate(jumpx,jumpy,jumpz)
                 deallocate(u_p, jump_all)

              enddo

!$OMP END DO
!$OMP END PARALLEL

      endif   !if (nelem_dg_glo .gt. 0)


      
!---------------------------------------------------------------------------
!    EXCHANGE DAMPING MATRIX FOR NONLINEAR ELASTICITY
!---------------------------------------------------------------------------
                                                     

     if (nmat_nle.gt.0) then                

          do i = 1,nrecv
                 in = node_recv(i)
                 recv_buffer(3*(i -1) +1) = mc(3*(in -1) +1)
                 recv_buffer(3*(i -1) +2) = mc(3*(in -1) +2)
                 recv_buffer(3*(i -1) +3) = mc(3*(in -1) +3)
          enddo
      
                call EXCHANGE_DOUBLE(3*nrecv,recv_buffer,3*nsend,send_buffer,&
                                        mpi_np,recv_length,send_length,&
                                        mpi_comm,mpi_stat,mpi_ierr,mpi_id)
      
          do i = 1,nsend
                 in = node_send(i)
!                 y_or_n_node = 1
!                 
!                 if (n_case.gt.0 .and. tag_case .ne. 16) then        
!                    if (Depth_nle .lt. zs_elev(in) .or. zs_elev(in) .lt. 0) y_or_n_node = 0
!                    if ((tag_case .gt. 0).and.(zs_all(in) .lt. 0.0d0)) y_or_n_node = 0
                              
!                 elseif(n_case .gt. 0 .and. tag_case .eq. 16) then
!                    if (Depth_nle .lt. zs_elev(in) .or. zs_elev(in) .lt. 0) y_or_n_node = 0
!                    if (vs_tria(in) .ge. 450.d0) y_or_n_node = 0
                              
!                 endif

                                               
!                 if(y_or_n_node .eq. 1) then
                 if(node_nle_4_mpi(in) .eq. 1) then
                    mc(3*(in -1) +1) = mc(3*(in -1) +1) +send_buffer(3*(i -1) +1)
                    mc(3*(in -1) +2) = mc(3*(in -1) +2) +send_buffer(3*(i -1) +2)
                    mc(3*(in -1) +3) = mc(3*(in -1) +3) +send_buffer(3*(i -1) +3)
                 endif   
          enddo
      
      
          do i = 1,nrecv
                  in = node_recv(i)
                 recv_buffer(3*(i -1) +1) = mck(3*(in -1) +1)
                 recv_buffer(3*(i -1) +2) = mck(3*(in -1) +2)
                 recv_buffer(3*(i -1) +3) = mck(3*(in -1) +3)
          enddo
      
          call EXCHANGE_DOUBLE(3*nrecv,recv_buffer,3*nsend,send_buffer,&
                                        mpi_np,recv_length,send_length,&
                                        mpi_comm,mpi_stat,mpi_ierr,mpi_id)
    
          do i = 1,nsend
                 in = node_send(i)
                 if(node_nle_4_mpi(in) .eq. 1) then
                    mck(3*(in -1) +1) = mck(3*(in -1) +1) +send_buffer(3*(i -1) +1)
                    mck(3*(in -1) +2) = mck(3*(in -1) +2) +send_buffer(3*(i -1) +2)
                    mck(3*(in -1) +3) = mck(3*(in -1) +3) +send_buffer(3*(i -1) +3)
                 endif   
          enddo

     endif !if (nmat_nle.gt.0) then        
        
   
         do i = 1,nrecv
            in = node_recv(i)
            recv_buffer(3*(i -1) +1) = fk(3*(in -1) +1)
            recv_buffer(3*(i -1) +2) = fk(3*(in -1) +2)
            recv_buffer(3*(i -1) +3) = fk(3*(in -1) +3)
         enddo
         
         call EXCHANGE_DOUBLE(3*nrecv,recv_buffer,3*nsend,send_buffer,&
                                   mpi_np,recv_length,send_length,&
                                   mpi_comm,mpi_stat,mpi_ierr,mpi_id)
         
         do i = 1,nsend
            in = node_send(i)
            fk(3*(in -1) +1) = fk(3*(in -1) +1) +send_buffer(3*(i -1) +1)
            fk(3*(in -1) +2) = fk(3*(in -1) +2) +send_buffer(3*(i -1) +2)
            fk(3*(in -1) +3) = fk(3*(in -1) +3) +send_buffer(3*(i -1) +3)
         enddo


!---------------------------------------------------------------------------
!     EXCHANGE JUMPS
!---------------------------------------------------------------------------

     if (nelem_dg_glo.gt.0) then    
   
        do i = 1,nrecv
           in = node_recv(i)
           recv_buffer(3*(i -1) +1) = jump(3*(in -1) +1)
           recv_buffer(3*(i -1) +2) = jump(3*(in -1) +2)
           recv_buffer(3*(i -1) +3) = jump(3*(in -1) +3)
        enddo
      
        call EXCHANGE_DOUBLE(3*nrecv,recv_buffer,3*nsend,send_buffer,&
                                   mpi_np,recv_length,send_length,&
                                   mpi_comm,mpi_stat,mpi_ierr,mpi_id)
      
        do i = 1,nsend
           in = node_send(i)
           jump(3*(in -1) +1) = jump(3*(in -1) +1) +send_buffer(3*(i -1) +1)
           jump(3*(in -1) +2) = jump(3*(in -1) +2) +send_buffer(3*(i -1) +2)
           jump(3*(in -1) +3) = jump(3*(in -1) +3) +send_buffer(3*(i -1) +3)
        enddo
      
      endif    
  
           
        if (rk_scheme .eq. 'RUNGEKUTTA') then         
           !Low storage Kappau = dq for u
           if (istage .eq. 1) then
              u_int = deltat*v1
              if(make_damping_yes_or_not.eq. 1 .and. damping_type .eq. 1) &
                                               v_int = -deltat*(fk + jump + mck*u1 + mc*v1 - fe)/mv
              if(make_damping_yes_or_not.eq. 1 .and. damping_type .eq. 3) &
                                               v_int = -deltat*(fk + jump + mc*v1 - fe)/mv
              if(make_damping_yes_or_not.eq. 0 .or. damping_type .eq. 2) v_int = -deltat*(fk + jump - fe)/mv

           else     
              u_int = A(istage)*u_int + deltat*v1
              if(make_damping_yes_or_not.eq. 1 .and. damping_type .eq. 1) &
                                               v_int = A(istage)*v_int - deltat*(fk + jump + mck*u1 + mc*v1 - fe)/mv
              if(make_damping_yes_or_not.eq. 1 .and. damping_type .eq. 3) &
                                               v_int = A(istage)*v_int - deltat*(fk + jump + mc*v1 - fe)/mv
              if(make_damping_yes_or_not.eq. 0 .or. damping_type .eq. 2) v_int = A(istage)*v_int - deltat*(fk + jump - fe)/mv

           endif
           u1 = u1 + B(istage)*u_int
           v1 = v1 + B(istage)*v_int
           
        endif

     enddo   !istage


      if(rk_scheme .eq. 'RUNGEKUTTA') then
         u2 = u1;  v2 = v1

      else
      

      
!---------------------------------------------------------------------------
!     TIME STEP FOR LEAP FROG
!---------------------------------------------------------------------------
!         if(its .eq. 0 .and. testmode .ne. 1) then

!$OMP PARALLEL DO PRIVATE(iaz,id)

                         
            do id = 1, nnod_loc 
           
               iaz = 3*(id -1) +1
               if (make_damping_yes_or_not .eq. 0 .or. damping_type .eq. 2) then
                   u2(iaz) = 2.0d0 * u1(iaz) - u0(iaz) &
                     + (fe(iaz) - fk(iaz) - jump(iaz)) / mv(iaz) *deltat2
               else
                   u2(iaz) = ( fe(iaz) - fk(iaz) - jump(iaz) - mck(iaz) * u1(iaz) & 
                       -mckty(iaz)/deltat*(fk(iaz)-fkold(iaz))*kclose/2 &                      
                       + (mc(iaz)*(2-kclose)/(4*deltat)) * u0(iaz) + (mv(iaz)/deltat2) &                         
                       * ( 2.0d0 * u1(iaz) - u0(iaz) ) ) / (mv(iaz)/deltat2 + mc(iaz)*(2-kclose)/(4*deltat) )
                       
               endif
               if(abs(u2(iaz)).lt. 1.0e-30) u2(iaz) = 0.d0  
               
               iaz = 3*(id -1) +2
               if (make_damping_yes_or_not.eq. 0 .or. damping_type .eq. 2) then
                  u2(iaz) = 2.0d0 * u1(iaz) - u0(iaz) &
                     + (fe(iaz) - fk(iaz) - jump(iaz)) / mv(iaz) *deltat2
               else
                  u2(iaz) = ( fe(iaz) - fk(iaz) - jump(iaz) - mck(iaz) * u1(iaz) &  
                       -mckty(iaz)/deltat*(fk(iaz)-fkold(iaz))*kclose/2 &                              
                       + (mc(iaz)*(2-kclose)/(4*deltat)) * u0(iaz) + (mv(iaz)/deltat2) &                         
                       * ( 2.0d0 * u1(iaz) - u0(iaz) ) ) &                                        
                       / (mv(iaz)/deltat2 + mc(iaz)*(2-kclose)/(4*deltat) )                                        
               endif
               if(abs(u2(iaz)).lt. 1.0e-30) u2(iaz) = 0.d0           

               iaz = 3*(id -1) +3
               if (make_damping_yes_or_not.eq. 0 .or. damping_type .eq. 2) then 
                  u2(iaz) = 2.0d0 * u1(iaz) - u0(iaz) &
                     + (fe(iaz) - fk(iaz) - jump(iaz)) / mv(iaz) *deltat2
               else
                  u2(iaz) = ( fe(iaz) - fk(iaz) - jump(iaz) - mck(iaz) * u1(iaz) & 
                         -mckty(iaz)/deltat*(fk(iaz)-fkold(iaz))*kclose/2 &                               
                         + (mc(iaz)*(2-kclose)/(4*deltat)) * u0(iaz) + (mv(iaz)/deltat2) &                         
                         * ( 2.0d0 * u1(iaz) - u0(iaz) ) )  / (mv(iaz)/deltat2 + mc(iaz)*(2-kclose)/(4*deltat) )                       
               endif
               if(abs(u2(iaz)).lt. 1.0e-30) u2(iaz) = 0.d0 


            enddo

            fkold=fk

!$END OMP PARALLEL DO
!         endif
                                                                                          
      endif                                                                                        

      do fn = 1,nfunc
            func_value(fn) = GET_FUNC_VALUE(nfunc,func_type,func_indx,func_data,nfunc_data,fn,tt2,0,0, &
                                                              MDOFMnum,MDOFforceinput,0,0,its)
      enddo
         
!---------------------------------------------------------------------------
!     SET DIRICHLET BOUNDARY CONDITIONS
!---------------------------------------------------------------------------

      if (nnode_dirX.gt.0) then
         do i = 1,nnode_dirX
            in = inode_dirX(i)
            iaz = 3*(in -1) +1; 
            u2(iaz) = 0.0d0; !180807 ty
            if(rk_scheme .eq. 'RUNGEKUTTA') v2(iaz) = 0.d0
            do fn = 1,nfunc
                if((Fel(fn,(3*(in -1) +1),2).eq.0.d0).and.(Fel(fn,(3*(in -1) +1),3).eq.-1.d0) ) then
                u2(iaz) = u2(iaz) + Fel(fn,(3*(in -1) +1),1) * func_value(fn)
                end if
            enddo
         enddo
      endif
      
      if (nnode_dirY.gt.0) then
         do i = 1,nnode_dirY
            in = inode_dirY(i)
            iaz = 3*(in -1) +2;
            u2(iaz) = 0.0d0; !180807 ty
            if(rk_scheme .eq. 'RUNGEKUTTA') v2(iaz) = 0.d0
            do fn = 1,nfunc
                if((Fel(fn,(3*(in -1) +2),2).eq.0.d0).and.(Fel(fn,(3*(in -1) +2),3).eq.-1.d0) ) then
                u2(iaz) = u2(iaz) + Fel(fn,(3*(in -1) +2),1) * func_value(fn)
                end if
            enddo
         enddo
      endif
      
      if (nnode_dirZ.gt.0) then
         do i = 1,nnode_dirZ
            in = inode_dirZ(i)
            iaz = 3*(in -1) +3; 
            u2(iaz) = 0.0d0; !180807 ty
            if(rk_scheme .eq. 'RUNGEKUTTA') v2(iaz) = 0.d0
            do fn = 1,nfunc
               if((Fel(fn,(3*(in -1) +3),2).eq.0.d0).and.(Fel(fn,(3*(in -1) +3),3).eq.-1.d0) ) then
               u2(iaz) = u2(iaz) + Fel(fn,(3*(in -1) +3),1) * func_value(fn)
               end if
            enddo
         enddo
      endif
 

!---------------------------------------------------------------------------
!     EXCHANGE STRESS-STRAIN-ROTATION FOR OUTPUT
!---------------------------------------------------------------------------

!! OMP DO COULD BE INSERTED HERE

         if(opt_out_var(4) .eq. 1) then
           allocate(double_send_buffer(6*nsend)); double_send_buffer = 0.0d0
           allocate(double_recv_buffer(6*nrecv)); double_recv_buffer = 0.0d0
         
            do i = 1,nsend
               in = node_send(i)
               double_send_buffer(6*(i -1) +1) = stress(6*(in -1) +1)
               double_send_buffer(6*(i -1) +2) = stress(6*(in -1) +2)
               double_send_buffer(6*(i -1) +3) = stress(6*(in -1) +3)
               double_send_buffer(6*(i -1) +4) = stress(6*(in -1) +4)
               double_send_buffer(6*(i -1) +5) = stress(6*(in -1) +5)
               double_send_buffer(6*(i -1) +6) = stress(6*(in -1) +6)
            enddo
         
        call EXCHANGE_DOUBLE(6*nsend,double_send_buffer,6*nrecv,double_recv_buffer,&
                                   mpi_np,double_send_length,double_recv_length,&
                                   mpi_comm,mpi_stat,mpi_ierr,mpi_id)
         
            do i = 1,nrecv
               in = node_recv(i)
               stress(6*(in -1) +1) = double_recv_buffer(6*(i -1) +1)
               stress(6*(in -1) +2) = double_recv_buffer(6*(i -1) +2)
               stress(6*(in -1) +3) = double_recv_buffer(6*(i -1) +3)
               stress(6*(in -1) +4) = double_recv_buffer(6*(i -1) +4)
               stress(6*(in -1) +5) = double_recv_buffer(6*(i -1) +5)
               stress(6*(in -1) +6) = double_recv_buffer(6*(i -1) +6)
            enddo
            
            deallocate(double_send_buffer,double_recv_buffer)
        
         endif


         if(opt_out_var(5) .eq. 1 .or. damping_type .eq. 2) then
         
           allocate(double_send_buffer(6*nsend)); double_send_buffer = 0.0d0
           allocate(double_recv_buffer(6*nrecv)); double_recv_buffer = 0.0d0

            do i = 1,nsend
               in = node_send(i)
               double_send_buffer(6*(i -1) +1) = strain(6*(in -1) +1)
               double_send_buffer(6*(i -1) +2) = strain(6*(in -1) +2)
               double_send_buffer(6*(i -1) +3) = strain(6*(in -1) +3)
               double_send_buffer(6*(i -1) +4) = strain(6*(in -1) +4)
               double_send_buffer(6*(i -1) +5) = strain(6*(in -1) +5)
               double_send_buffer(6*(i -1) +6) = strain(6*(in -1) +6)
            enddo
         
        call EXCHANGE_DOUBLE(6*nsend,double_send_buffer,6*nrecv,double_recv_buffer,&
                                   mpi_np,double_send_length,double_recv_length,&
                                   mpi_comm,mpi_stat,mpi_ierr,mpi_id)
         
            do i = 1,nrecv
               in = node_recv(i)
               strain(6*(in -1) +1) = double_recv_buffer(6*(i -1) +1)
               strain(6*(in -1) +2) = double_recv_buffer(6*(i -1) +2)
               strain(6*(in -1) +3) = double_recv_buffer(6*(i -1) +3)
               strain(6*(in -1) +4) = double_recv_buffer(6*(i -1) +4)
               strain(6*(in -1) +5) = double_recv_buffer(6*(i -1) +5)
               strain(6*(in -1) +6) = double_recv_buffer(6*(i -1) +6)
            enddo
            

            deallocate(double_send_buffer,double_recv_buffer)
        
         endif

         if(opt_out_var(6) .eq. 1) then
         
            do i = 1,nsend
               in = node_send(i)
               send_buffer(3*(i -1) +1) = omega(3*(in -1) +1)
               send_buffer(3*(i -1) +2) = omega(3*(in -1) +2)
               send_buffer(3*(i -1) +3) = omega(3*(in -1) +3)
            enddo
         
         call EXCHANGE_DOUBLE(3*nsend,send_buffer,3*nrecv,recv_buffer,&
                                   mpi_np,send_length,recv_length,&
                                   mpi_comm,mpi_stat,mpi_ierr,mpi_id)
         
            do i = 1,nrecv
               in = node_recv(i)
               omega(3*(in -1) +1) = recv_buffer(3*(i -1) +1)
               omega(3*(in -1) +2) = recv_buffer(3*(i -1) +2)
               omega(3*(in -1) +3) = recv_buffer(3*(i -1) +3)
               
            enddo
         endif


!---------------------------------------------------------------------------
!     UPDATE DAMPING TENSOR 
!---------------------------------------------------------------------------

      if (damping_type .eq. 2) then
          
          call UPDATE_DAMPING_TENSOR(nnod_loc, N_SLS, frequency_range, &
                                     strain_visc, strain, deltat)
      endif
      
!---------------------------------------------------------------------------
!     WRITE MONITORS IN OUTPUT 
!---------------------------------------------------------------------------

      if (nmonitors_lst .ge. 1) then
                        
        if (mod(its,ndt_mon_lst) .eq. 0  .and. count_monitor .gt. 0 ) then 

            call WRITE_OUTPUT(nmonitors_lst, mpi_id, el_monitor_lst, local_el_num, nelem_loc, &
                              con_spx_loc, con_nnz_loc, sdeg_mat, nmat, &
                              u2, u1, u0, v1, nnod_loc, &
                              stress, strain, omega, &
                              xr_monitor_lst, yr_monitor_lst, zr_monitor_lst, & 
                              tt1, deltat, deltat2, &
                              opt_out_var, count_monitor, monitor_file,debug, &
                              mu_nle, gamma_nle, b_instabilitycontrol, instability_maxval, b_instability_abort)            
         endif         
                        
      endif


      !!!!!!!!!TY!!!!!!!!!!!!!!

      MDOFinput=0
      MDOFinputbuffer=0
      MDOFinputfk=0
      MDOFinputfkbuffer=0

      call WRITE_OUTPUT_MDOF(nmonitors_lst, mpi_id, el_monitor_lst, local_el_num, nelem_loc, &
                              con_spx_loc, con_nnz_loc, sdeg_mat, nmat, &
                              u2, u1, u0, v1, nnod_loc, &
                           !   stress, strain, omega, &
                              xr_monitor_lst, yr_monitor_lst, zr_monitor_lst, & 
                              tt1, deltat, deltat2, &
                              count_monitor,  &
                              MDOFinput,mpi_np, &
                              utyty1,utyty2,utyty3,MDOFMnum,fk,fty1,fty2,MDOFinputfk)  

      if(mpi_np.gt.1) then
        call EXCHANGE_ACC(MDOFMnum*3,MDOFinput,MDOFMnum*3,MDOFinputbuffer,&        ! non 0 => 0
                                   mpi_np,MDOFMnum*3,MDOFMnum*3,&
                                   mpi_comm,mpi_stat,mpi_ierr,mpi_id)
      end if
      
      if (mpi_id .eq. 0) then
        do MDOFexi=1,MDOFMnum*3 
            do MDOFexj=1,mpi_np  
                MDOFinput(MDOFexi)=MDOFinput(MDOFexi)+MDOFinputbuffer(MDOFMnum*3*(MDOFexj-1)+MDOFexi)
            end do
        end do

      end if
      MDOFinputbuffer=0  
      if(mpi_np.gt.1) then
        call EXCHANGE_ACC2(MDOFMnum*3,MDOFinput,MDOFMnum*3,MDOFinputbuffer,&        ! 0 => non 0
                                   mpi_np,MDOFMnum*3,MDOFMnum*3,&
                                   mpi_comm,mpi_stat,mpi_ierr,mpi_id)
      end if
      if (mpi_id .ne. 0) then
         MDOFinput(1:MDOFMnum*3)=MDOFinputbuffer(1:MDOFMnum*3)
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(mpi_np.gt.1) then
        call EXCHANGE_ACC(MDOFMnum*3,MDOFinputfk,MDOFMnum*3,MDOFinputfkbuffer,&        
                                   mpi_np,MDOFMnum*3,MDOFMnum*3,&
                                   mpi_comm,mpi_stat,mpi_ierr,mpi_id)
      end if
      
      if (mpi_id .eq. 0) then
        do MDOFexi=1,MDOFMnum*3   
            do MDOFexj=1,mpi_np  
                MDOFinputfk(MDOFexi)=MDOFinputfk(MDOFexi)+MDOFinputfkbuffer(MDOFMnum*3*(MDOFexj-1)+MDOFexi)
            end do
        end do
      end if
      MDOFinputfkbuffer=0  
      if(mpi_np.gt.1) then
        call EXCHANGE_ACC2(MDOFMnum*3,MDOFinputfk,MDOFMnum*3,MDOFinputfkbuffer,&       
                                   mpi_np,MDOFMnum*3,MDOFMnum*3,&
                                   mpi_comm,mpi_stat,mpi_ierr,mpi_id)
      end if
      if (mpi_id .ne. 0) then
         MDOFinputfk(1:MDOFMnum*3)=MDOFinputfkbuffer(1:MDOFMnum*3)
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!TY
      MDOFforceinput=0  
      if(NumofBlgs .gt. 0) then 
      do MDOFBlgi=1,NumofBlgs !MDOFBlgi
        MDOFacl(MDOFBlgi,1:3)=(-1.d0)*MDOFinput((3*(MDOFbf+MDOFBlgi-1)-2):(3*(MDOFbf+MDOFBlgi-1)))
            do MDOFBlgti=1,Blg(MDOFBlgi)%ndt
            !!!!X
            call ExplicitCalShear (MDOFBlgi,Blgdof(MDOFBlgi),MDOFacl(MDOFBlgi,1), &
                        1,NewOrNot(MDOFBlgi))
            !!!!Y
            call ExplicitCalShear (MDOFBlgi,Blgdof(MDOFBlgi),MDOFacl(MDOFBlgi,2), &
                        2,NewOrNot(MDOFBlgi))
            end do            
            MDOFforceinput(1+6*(MDOFftranfg(mpi_id+1)+MDOFBlgi-2))=Blg(MDOFBlgi)%MDOFItF(1,1)  
            MDOFforceinput(2+6*(MDOFftranfg(mpi_id+1)+MDOFBlgi-2))=Blg(MDOFBlgi)%MDOFItF(1,2) 
            MDOFforceinput(4+6*(MDOFftranfg(mpi_id+1)+MDOFBlgi-2))=MDOFforceinput(2+6*(MDOFftranfg(mpi_id+1)+MDOFBlgi-2))*Blg(MDOFBlgi)%Floor_h/2.0d0   !Mx, calculated by Fy
            MDOFforceinput(5+6*(MDOFftranfg(mpi_id+1)+MDOFBlgi-2))=MDOFforceinput(1+6*(MDOFftranfg(mpi_id+1)+MDOFBlgi-2))*Blg(MDOFBlgi)%Floor_h/2.0d0  !My, calculated by Fx  
      end do  !MDOFBlgi
      
      if(mod(its,ndt_mon_lst) .eq. 0) then
      call GenBResultFile(its,tt1)      
      end if
    end if   !NumofBlgs.gt.0
      
     
      MDOFforceinputbuffer=0
      if(mpi_np.gt.1) then
        call EXCHANGE_ACC(MDOFMnum*6,MDOFforceinput,MDOFMnum*6,MDOFforceinputbuffer,&        !non 0 => 0
                                   mpi_np,MDOFMnum*6,MDOFMnum*6,&
                                   mpi_comm,mpi_stat,mpi_ierr,mpi_id)
      end if
      if (mpi_id .eq. 0) then
        do MDOFexi=1,MDOFMnum*6 
            do MDOFexj=1,mpi_np  
                MDOFforceinput(MDOFexi)=MDOFforceinput(MDOFexi)+MDOFforceinputbuffer(MDOFMnum*6*(MDOFexj-1)+MDOFexi)
            end do
        end do
      end if
      MDOFforceinputbuffer=0  
      if(mpi_np.gt.1) then
        call EXCHANGE_ACC2(MDOFMnum*6,MDOFforceinput,MDOFMnum*6,MDOFforceinputbuffer,&        !0 => non 0
                                   mpi_np,MDOFMnum*6,MDOFMnum*6,&
                                   mpi_comm,mpi_stat,mpi_ierr,mpi_id)
      end if
      if (mpi_id .ne. 0) then
                 MDOFforceinput(1:MDOFMnum*6)=MDOFforceinputbuffer(1:MDOFMnum*6)
      end if
  

      if (b_instabilitycontrol) then
        if (b_instability_abort) then
          write(*,*) 'Worker ', mpi_id, ': instability detected. Signalling to other workers.'
        endif

        ! Send & read all workers statuses
        call MPI_ALLGATHER(b_instability_abort, 1, MPI_LOGICAL, b_instability_abort_all, 1, MPI_LOGICAL, mpi_comm, mpi_ierr)

        ! Abort if one worker fails
        if (any(b_instability_abort_all)) then
          if (mpi_id .eq. 0) then
            write(*,*) 'Instability detected, aborting.'           
          endif

          final_time = MPI_WTIME()
          call MPI_FINALIZE(mpi_ierr)

          if (mpi_id .eq. 0) then
            write(*,*) 'Jobs aborted. '
            write(*,*) 'Elapsed wall clock time: ', final_time - start_time, ' seconds'
            write(*,*) 'Instability detected at iteration = ', its, ', TIME =', tt1
          endif
          call EXIT(EXIT_INSTAB)

        endif

      endif
      
      
!---------------------------------------------------------------------------
!     WRITE BACKUP FILES
!---------------------------------------------------------------------------

    if ((nsnaps .gt. 0) .and. (isnap .le. nsnaps)) then
            
        if (its .ge. itersnap(isnap)) then

           if (mpi_id.eq.0) write(*,*) 'Writing Backup files at TIME : ',tt2

               call WRITE_FILEOUT_GRID(file_outU1,isnap,mpi_id,3*nnod_loc,u1,& 
                                        xx_spx_loc,yy_spx_loc,zz_spx_loc,&
                                        local_node_num)
               call WRITE_FILEOUT_GRID(file_outV1,isnap,mpi_id,3*nnod_loc,v1,& 
                                        xx_spx_loc,yy_spx_loc,zz_spx_loc,&
                                        local_node_num)
               
               isnap = isnap + 1
               
        endif
    endif

!---------------------------------------------------------------------------
!     UPDATE VARIABLES
!---------------------------------------------------------------------------
       
      if (rk_scheme .eq. 'RUNGEKUTTA') then
         v1 = v2
         u1 = u2
      else
         v1 = (u2-u0)/(2.0d0*deltat) 
         u0 = u1
         u1 = u2
      endif
          
      if (its .le. 200) then 
          call MPI_BARRIER(mpi_comm, mpi_ierr)
          final_time = MPI_WTIME()
          if (mpi_id .eq. 0) write(*,*) 'Time for a single step', final_time - start_time
      endif

         
!---------------------------------------------------------------------------
!     COMPUTE ERROR
!---------------------------------------------------------------------------
   if (itime .le. ntime_err) then 
         
      if (testmode .eq. 1 .and. abs(time_error(itime)-tt1) .le. deltat) then 
    
         call COMPUTE_ENERGY_ERROR(nnod_loc, u1, v1, tt1, nelem_loc, con_spx_loc, con_nnz_loc,&
                              nmat, prop_mat, sdeg_mat, tag_mat, &
                              local_node_num, local_el_num, nelem_dg,&
                              alfa11,alfa12,alfa13,&
                              alfa21,alfa22,alfa23,&
                              alfa31,alfa32,alfa33,&
                              beta11,beta12,beta13,&
                              beta21,beta22,beta23,&
                              beta31,beta32,beta33,&
                              gamma1,gamma2,gamma3,&
                              delta1,delta2,delta3,&
                              xx_spx_loc, yy_spx_loc, zz_spx_loc, mpi_id, el_new, &
                              nelem_dg_glo, &
                              nsend_jump,node_send_jump, &
                              nrecv_jump,node_recv_jump, &
                              send_length_jump, recv_length_jump, &  
                              mpi_np, mpi_comm, &
                              con_nnz_dg, con_spx_dg) 
      
        itime = itime + 1

     endif
  endif     
  
     tt0 = tt0 + deltat; tt1 = tt1 + deltat; tt2 = tt2 + deltat;
  enddo   !loop on its
      
!---------------------------------------------------------------------------
!     END TIME LOOP
!---------------------------------------------------------------------------
    
      if(mpi_id .eq. 0) then 
         finish = MPI_WTIME();  time_in_seconds = finish - start                
      endif

      if (mpi_id.eq.0) then
         write(*,'(A)')
         write(*,'(A,F20.3,A)')'Time loop completed in ',time_in_seconds,' s'
         write(*,'(A)')'-------------------------------------------------------'
         write(*,'(A)')
      endif
      
!---------------------------------------------------------------------------
!     WRITING PEAK GROUD MAP
!---------------------------------------------------------------------------

      if (nmonitors_pgm.ge.1) then        
               
        if(mpi_id.eq.0) then
           write(*,'(A)')
           write(*,'(A)')'--------------Writing Peak Ground Maps-----------------'
           write(*,'(A)')
        endif 
                                
        filePG='PGD'
        call WRITE_FILEOUT_PG(monitor_file,filePG,mpi_id,nmonitors_pgm,max_u,9)
        filePG='PGV'
        call WRITE_FILEOUT_PG(monitor_file,filePG,mpi_id,nmonitors_pgm,max_v,9)
        filePG='PGA'
        call WRITE_FILEOUT_PG(monitor_file,filePG,mpi_id,nmonitors_pgm,max_a,9)                        
        filePG='PGO'
        call WRITE_FILEOUT_PG(monitor_file,filePG,mpi_id,nmonitors_pgm,max_o,3)          
      endif                                                                                

      
!---------------------------------------------------------------------------
!     DEALLOCATION 
!---------------------------------------------------------------------------
        
    !!!!!!!!!!!!!!!!!!!!!!!!!!ty!!!!!!!!!!!!!!!!!!!!!!!!    
    if(NumofBlgs.gt.0) then
    deallocate(Blg)
    deallocate(Blgdof)
    deallocate(MDOFacl)
    deallocate(NewOrNot)
    end if
    deallocate(MDOFinput) 
    deallocate(MDOFinputbuffer)
    deallocate(utyty1)
    deallocate(utyty2)
    deallocate(utyty3)
    !!!!!!!!!!!!!!!!!!!!!!!!!!ty!!!!!!!!!!!!!!!!!!!!!!!!  
      
      if (debug .eq. 1) deallocate(mu_nle,gamma_nle)
      if (damping_type .eq. 2) deallocate(strain_visc)
      if (damping_type .eq. 3) deallocate(damp_term)
      
      
      deallocate(func_value, u1, u2, v1, fe, fk, fkold, mv, jump,send_buffer,recv_buffer,v_int,u_int)
      if (rk_scheme .ne. 'RUNGEKUTTA') deallocate(u0)

      if (nelem_dg_glo .gt. 0) deallocate(send_buffer_jump, recv_buffer_jump, jump_minus)
      if (opt_out_var(4) .eq. 1 .or. opt_out_var(5) .eq. 1 .or. opt_out_var(6) .eq. 1) then 
         deallocate(node_counter)
      endif   
      if (rk_scheme .eq. 'RUNGEKUTTA') deallocate(v2, A,b,c)
      if (nmonitors_pgm .ge. 1) deallocate(max_u, max_v, max_a, max_o)
      if (nmat_nle.gt.0) deallocate(con_spx_loc_nle,Depth_nle_el,node_nle_4_mpi)
      

      !ADAPTIVE TIME STEP 
      !if (time_adpt .eq. 1) deallocate(fe_1,fe_2,fe_3,fext,fext_1,fext_2,fext_3,fk_1,fk_2,eta_1,eta,dg,d2g)


     end subroutine TIME_LOOP
