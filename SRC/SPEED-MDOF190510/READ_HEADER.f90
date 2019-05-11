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

!> @brief Inizialization for parallel computation.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] headerfile file name for the header
!!
!> @param[out] gridfile  file name containing grid description
!> @param[out] matefile  file name containing material description, boundary confitions, 
!!                 forcing terms
!> @param[out] time_step  time integration step 
!> @param[out] start_time  start time
!> @param[out] stop_time  stop time
!> @param[out] option_out_var  options for the output of the variables  d v a e s w
!!                       1 - on, 0 - off                      ex: 1 0 0 0 0 0 (only displ.)                                
!> @param[out] option_out_form  options for the output format
!> @param[out] option_out_data  options for the ouptut data
!> @param[out] nsnapshots  number of snapshots of the solution
!> @param[out] t_snapshot  time of the snapshot
!> @param[out] ndt_monitor  number of time steps between two savings
!> @param[out] deltat_fixed   fixed time step or not ('YES' or 'NOT') 
!> @param[out] depth_search_mon_pgm  depth for searching monitors for peak ground maps
!> @param[out] ndt_mon_pgm  number of time steps between two savings for peak ground maps
!> @param[out] n_pgm  number of peak ground map
!> @param[out] rotation_angle_mon_pgm
!> @param[out] monfile_pgm  1 of the file PGDM already exists 0 otherwise
!> @param[out] depth_search_mon_lst  depth for searching monitors 
!> @param[out] n_lst  number of monitor list
!> @param[out] monfile_lst  1 of the file MLST already exists 0 otherwise
!> @param[out] dg_c  select DG method (-1 SIPG, 0-IIPG, 1-NIPG)
!> @param[out] pen_c  penality parameter for DG discretization
!> @param[out] scheme  time integration scheme - if absent LEAP_FROG scheme is adopdted, if keyword
!!               RUNGEKUTTA is specified, Runge-Kutta scehme is employed
!> @param[out] order  order of accuracy for Runge-Kutta scheme
!> @param[out] stages  stages for Runge-Kutta scheme
!> @param[out] testmode  a special test case is addressed
!> @param[out] time_err(ntime_err)  time instants for which the norms of the error are computed
!!                            (active if testmode keyword is present)
!> @param[out] damping_type  choose damping scheme 1-K&K, 2-SLS, 3-Rayleigh
!> @param[out] b_failoncoeffs  if damping=2, quit if anelastic coefficients are negative.
!> @param[out] b_setuponly  if true, quit before the time loop
!> @param[out] b_failCFL  if true, fails if CFL does not hold
!> @param[out] b_instabilitycontrol  if true, quit whenever the simulation diverges
!> @param[out] instability_maxval  max abs. value for monitored Stress

      subroutine READ_HEADER(headerfile,gridfile,matefile,mpifile,monfile,bkpfile,&   
                             time_step,start_time,stop_time,&
                             option_out_var,&          
                             option_out_data,option_out_form,&
                             nsnapshots,t_snapshot,&
                             ndt_monitor,&                
                             deltat_fixed,&                
                             depth_search_mon_pgm,ndt_mon_pgm,n_pgm,&        
                             rotation_angle_mon_pgm,monfile_pgm,&        
                             depth_search_mon_lst,n_lst, &
                             monfile_lst,dg_c,pen_c, scheme, order, stages, testmode, &
                             ntime_err, time_err, damping_type, &
                             b_failoncoeffs, b_setuponly, b_failCFL, b_instabilitycontrol, instability_maxval)

      use speed_exit_codes

      implicit none
      
      character*70 :: headerfile,gridfile,matefile,mpifile,monfile,bkpfile 
      character*70 :: inline
      character*10 :: scheme        
      character*8 :: keyword
      character*3 :: deltat_fixed        

      integer*4 :: nsnapshots, testmode
      integer*4 :: option_out_data,option_out_form
      integer*4 :: status
      integer*4 :: ileft,iright, order, stages
      integer*4 :: arglen      
      integer*4 :: i,j,im,is
      integer*4 :: ndt_monitor                
      integer*4 :: ndt_mon_pgm,n_pgm                                
      integer*4 :: monfile_pgm                                        
      integer*4 :: n_lst, ntime_err, itime
      integer*4 :: monfile_lst, damping_type
      integer*4 :: file_row = 0

      integer*4, dimension (6) :: option_out_var 

      real*8 :: val,dg_c, pen_c
      real*8 :: time_step,start_time,stop_time
      real*8 :: depth_search_mon_pgm                                
      real*8 :: rotation_angle_mon_pgm                                
      real*8 :: depth_search_mon_lst

      real*8, dimension(nsnapshots) :: t_snapshot
      real*8, dimension(ntime_err) :: time_err

      ! If .TRUE., fails when a mandatory argument is missing
      logical :: b_fail_unset_args = .FALSE.

      ! New settings
      ! Instability
      real*8 :: instability_maxval
      logical :: b_instabilitycontrol
      logical :: b_failoncoeffs, b_setuponly, b_failCFL


      im = 0;       is = 0;   itime = 0;
      n_pgm = 0;    monfile_pgm = 0                 
      n_lst = 0;    monfile_lst = 0                 

      time_err = 0      


      open(40,file=headerfile)
      
      do    
         read(40,'(A)',IOSTAT = status) inline
         file_row = file_row + 1
         
         if (status.ne.0) exit
         
         ! Skip comments
         if (inline(1:1) .eq. ' ') then
            cycle
         endif

         !!!! Parse ketword arguments
         ileft = 1
         iright = len_trim(inline)
         
         ! Compute index to first non-keyword argument
         do i = 1,iright
            if (inline(i:i) .eq.' ') exit
         enddo
         ileft = i + 1
         keyword = inline(1:(ileft-2))

         ! Now ileft points after the first blank
         arglen = len_trim(inline(ileft:iright))

!         write(*,*) 'keyword=', keyword, '  inline= ', inline
!         write(*,*) 'ileft= ', ileft, '   iright= ', iright, '   arglen=', arglen
!         write(*,*) 'inline(ileft:iright)=', inline(ileft:iright)
!         write(*,*)

         ! Remove this for keywords without arguments!
         if (arglen .eq. 0) then
            write(*,'(A,I3,A,A,A)') 'FATAL in SPEED.input, row', file_row, ': no argument given for command "', trim(keyword), '"'
            call EXIT(EXIT_SYNTAX_ERROR)
         endif

!         write(*,*) 'Comparing ', inline(1:(ileft-2)), ', ileft=', ileft

         select case (keyword)

           case('GRIDFILE')
            read(inline(ileft:iright),*) gridfile

           case('DGMETHOD') 
            read(inline(ileft:iright),*) dg_c

           case('PENALIZC') 
            read(inline(ileft:iright),*) pen_c
         
           case('TIMESTEP')
            read(inline(ileft:iright),*) time_step
         
           case('STARTIME')
            read(inline(ileft:iright),*) start_time
         
           case('STOPTIME')
            read(inline(ileft:iright),*) stop_time
         
           case('SNAPSHOT')
            is = is +1
            read(inline(ileft:iright),*) t_snapshot(is)

           case('TIMEFIXE')                                
            read(inline(ileft:iright),*) time_step                                        
            deltat_fixed = 'yes'                                                        

           case('TMONITOR')                                        
            read(inline(ileft:iright),*) ndt_monitor                
            
           case('TIMESCHM') 
            read(inline(ileft:iright),*) scheme, order, stages
             
           case('TESTMODE') 
            testmode = 1    

           case('MATFILE')
            read(inline(ileft:iright),*) matefile
           
           case('MPIFILE')
            read(inline(ileft:iright),*) mpifile   
           
           case('MONFILE')
            read(inline(ileft:iright),*) monfile   
           
           case('BKPFILE')
            read(inline(ileft:iright),*) bkpfile   

           case('OPTIOUT')
            read(inline(ileft:iright),*) option_out_var(1),option_out_var(2),&  
                        option_out_var(3), option_out_var(4),option_out_var(5),option_out_var(6),&  
                        option_out_form, option_out_data
                        
           case('TIMEERR') 
            itime = itime + 1                
            read(inline(ileft:iright),*) time_err(itime)  

           case('DAMPING')
            read(inline(ileft:iright),*) damping_type   
    

           case('PGDM')                                        
            n_pgm = 1                                                        
            read(inline(ileft:iright),*) depth_search_mon_pgm, &
                        ndt_mon_pgm,rotation_angle_mon_pgm,monfile_pgm                        
           
           case('MLST')                        
            n_lst = 1                                
            read(inline(ileft:iright),*) depth_search_mon_lst,monfile_lst
    
           !! New flags

           case('FAILCFL')
            ! If specified as "FAILCFL T", quit if CFL condition does not hold
            read(inline(ileft:(ileft+arglen)),'(L)') b_failCFL

           case('FAILINST')
            ! If specified as "FAILINST T", enable instability control
            read(inline(ileft:iright),*) b_instabilitycontrol

           case('SETUPONL')
            ! If specified as "SETUPONL T", quit before starting the time loop
            read(inline(ileft:iright),*) b_setuponly

           case('FAILCOEF')
            ! If specified as "FAILCOEF T", quit if any of the computed
            ! anelastic coefficients is negative [damping 2]
            read(inline(ileft:iright),*) b_failoncoeffs

           case('INSTVAL')
            ! If specified, overwrite preset instability threshold
            read(inline(ileft:iright),*) instability_maxval
            !b_instabilitycontrol = .true.

           ! Ignored keywords (for compatibility)
           case('r_f')
            write(*,'(A,I3,A,A)') 'In SPEED.input, row', file_row, ': ignored field ', keyword

           ! Fail if keyword is not recognised
           case default
            write(*,'(A,I3,A,A)') 'FATAL in SPEED.input, row', file_row, ': unknown keyword ', keyword
            call EXIT(EXIT_SYNTAX_ERROR)

         end select  

      enddo
      
      
! Snapshots reordering
      
      if (nsnapshots.gt.1) then
         do i = 1,nsnapshots-1
            do j = i+1, nsnapshots
               if (t_snapshot(i).gt.t_snapshot(j)) then
                  val = t_snapshot(i)
                  t_snapshot(i) = t_snapshot(j)
                  t_snapshot(j) = val
               endif
            enddo
         enddo
      endif
      
      close(40)
      
      return
      end subroutine READ_HEADER

