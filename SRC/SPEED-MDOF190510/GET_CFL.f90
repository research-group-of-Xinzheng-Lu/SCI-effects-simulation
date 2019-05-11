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

!> @brief Computes CFL condition when mechanical properties are given
!! block by block.
!! @author Ilario Mazzieri
!> @date September, 2013 - Creation
!> @version 1.0
!> @param[in,out]  time_step  deltat used for the time integration
!> @param[in] nn_loc number of local nodes
!> @param[in] loc_n_num local node numeration
!> @param[in] nm number of materials
!> @param[in] sdeg   polynomial degree vector 
!> @param[in] pm  material property vector (1-rho, 2-lambda, 3-mu, 4-zeta)
!> @param[in] xx_loc x-coordinate of local nodes
!> @param[in] yy_loc y-coordinate of local nodes
!> @param[in] zz_loc z-coordinate of local nodes
!> @param[in] cs_nnz_loc length of cs_loc
!> @param[in] cs_loc local connectivity vector
!> @param[in] time_step_cfl  maximum deltat available in order to avoid instability
!> @param[in] fmax  maximum frequency of the wave
!> @param[in] deltat_fixed  'yes' deltat fixed even if is extremely small for the time integration
!                     'not' deltat will be changed if is extremely small for the time integration
!> @param[in] mpi_comm  mpi common world
!> @param[in] mpi_np  total mpi processor number
!> @param[in] mpi_id  number of mpi process 
!> @param[in] b_failCFL  fail if CFL does not hold
!> @note Compute the CFL condition, i.e., deltat <= GET_CFL = h_min/(vp_max*N^2) 
!!  and print on the screen what is the percent of GET_CFL you are using for the time integration

       subroutine GET_CFL(time_step, nn_loc, loc_n_num, nm, tm, pm, sdeg, &
                             xx_loc, yy_loc, zz_loc, cs_nnz_loc, cs_loc, & 
                             time_step_cfl, fmax, deltat_fixed, &
                             mpi_comm, mpi_np, mpi_id, b_failCFL)   

        use speed_exit_codes

        implicit none
         
        include 'SPEED.MPI'

        character*3 :: deltat_fixed

        integer*4 :: nm, im, nn_loc, cs_nnz_loc, mpi_comm, mpi_np, mpi_err, mpi_id
        integer*4 :: nel_loc, ic1, ic2, n1, n2, istart, ifin
        integer*4 :: ie, i, j, nn, mcode, smcode, sdeg_deltat_cfl, sdeg_npoints
        
        integer*4, dimension(1) :: pos
        integer*4, dimension(nm) :: tm, sdeg
        integer*4, dimension(nn_loc) :: loc_n_num

        integer*4, dimension(0:cs_nnz_loc) :: cs_loc

        real*8 :: time_step, time_step_cfl, fmax        
        real*8 :: length_min,length,vs_length_min,vs_length,length_vp_min,length_vp
        real*8 :: percent_deltat
        real*8 :: num_of_points,vs_npoints,vp_deltat_cfl
        real*8 :: rho,lambda,mu,vs,vp
        real*8 :: x1_length_vp_min,y1_length_vp_min,z1_length_vp_min
        real*8 :: x2_length_vp_min,y2_length_vp_min,z2_length_vp_min

        real*8, dimension(:), allocatable :: ct,ww
        real*8, dimension(:), allocatable :: time_step_glo
        real*8, dimension(nn_loc) :: xx_loc, yy_loc, zz_loc
              
        real*8, dimension(:,:), allocatable :: dd
        real*8, dimension(nm,4) :: pm

        logical :: b_failCFL, b_CFL_failure = .FALSE.

        nel_loc = cs_loc(0) - 1

        length_vp_min = 1.d10
        vs_length_min = 1.d10
       

    
          
         do ie = 1, nel_loc

               im = cs_loc(cs_loc(ie -1) +0) 
               nn = sdeg(im) +1

               smcode=sdeg(im)        
               rho = pm(im,1)   
               lambda = pm(im,2)
               mu = pm(im,3)
               vs = dsqrt(mu/rho)
               vp = dsqrt((lambda+2*mu)/rho)

               n1 = nn*nn*(1 -1) +nn*(1 -1) +1
               n2 = nn*nn*(nn -1) +nn*(nn -1) +nn
               !ic1 = cs_loc(cs_loc(ie -1) +n1)
               !ic2 = cs_loc(cs_loc(ie -1) +n2)

               istart = cs_loc(ie -1) +n1
               ifin = cs_loc(ie -1) +n2

 
!               write(*,*) ic1, ic2
!               write(*,*) xx_loc(ic1),yy_loc(ic1),zz_loc(ic1)
!               write(*,*) xx_loc(ic2),yy_loc(ic2),zz_loc(ic2)
!               read(*,*) 
!               write(*,*) cs_loc(ie -1) +n2, cs_loc(cs_loc(ie -1) +n2), cs_loc
!               read(*,*)
               
               do ic1 = istart, ifin 
                  do ic2 = ic1 + 1, ifin               
               
                  if(cs_loc(ic1) .ne. 0 .and. cs_loc(ic2) .ne. 0) then
               
                     length=dsqrt((xx_loc(cs_loc(ic1))-xx_loc(cs_loc(ic2)))**2 + &
                                  (yy_loc(cs_loc(ic1))-yy_loc(cs_loc(ic2)))**2 + &
                                  (zz_loc(cs_loc(ic1))-zz_loc(cs_loc(ic2)))**2)
                        
                     length_vp = length/vp
                     !write(*,*) cs_loc(ic1),cs_loc(ic2), length, vp
                     

                     if (length_vp_min .gt. length_vp) then
                        length_vp_min = length_vp
                        x1_length_vp_min = xx_loc(cs_loc(ic1))
                        y1_length_vp_min = yy_loc(cs_loc(ic1))
                        z1_length_vp_min = zz_loc(cs_loc(ic1))
                        x2_length_vp_min = xx_loc(cs_loc(ic2))
                        y2_length_vp_min = yy_loc(cs_loc(ic2))
                        z2_length_vp_min = zz_loc(cs_loc(ic2))
                        vp_deltat_cfl = vp
                        sdeg_deltat_cfl = smcode
                        length_min = length
                     endif


                     vs_length = vs/length

                     if (vs_length_min.gt.vs_length) then
                        vs_length_min = vs_length
                        vs_npoints = vs
                        sdeg_npoints = smcode
                     endif
      

                                              
                  endif
 
               enddo
             enddo   
             
             !read(*,*)     
         enddo
        
        nn = sdeg_deltat_cfl + 1

        !write(*,*) length_vp_min,length_min,vp
        !read(*,*)

        time_step_cfl = length_vp_min  !/((nn-1)**2)
        
        allocate(time_step_glo(mpi_np))
        
        call MPI_BARRIER(mpi_comm, mpi_err)
        call MPI_ALLGATHER(time_step_cfl, 1, SPEED_DOUBLE, time_step_glo, 1, SPEED_DOUBLE, mpi_comm, mpi_err)

      
        pos = minloc(time_step_glo)
        
        
        deallocate(time_step_glo)

        if((mpi_id + 1) .eq. pos(1)) then
        
           write(*,'(A)')' '
           write(*,'(A)') '--------Stability for time advancing algorithm --------'
           write(*,'(A,E12.4)') 'Min. el. length :', length_vp_min*vp_deltat_cfl
           write(*,'(A,E12.4)') 'Min. Vp         :', vp_deltat_cfl
           write(*,'(A,E12.4,E12.4,E12.4)')'Min. node 1     : ',x1_length_vp_min,y1_length_vp_min,z1_length_vp_min
           write(*,'(A,E12.4,E12.4,E12.4)')'Min. node 2     : ',x2_length_vp_min,y2_length_vp_min,z2_length_vp_min
           write(*,'(A)') '-------------------------------------------------------'
           if (time_step.le.time_step_cfl) then
                percent_deltat=time_step/time_step_cfl*100
                if (percent_deltat.le.1) then
                        write(*,'(A,E12.4)')'WARNING!!! deltat CFL =',&
                        time_step_cfl
                        write(*,'(A,F6.2,A)')'You are using ',percent_deltat,&
                                '% of critical time step.'
                        write(*,'(A)')'This time step is excedingly lower the deltat CFL'
                        if (deltat_fixed.eq.'not') then  
                                write(*,'(A)')'deltat chosen will be substituted with 1% of deltat CFL'
                                time_step=time_step_cfl*0.01
                                write(*,'(A,E12.4)')'deltat chosen :',time_step
                        endif
                elseif (percent_deltat.le.25) then
                        write(*,'(A,E12.4)')'OK!!! deltat CFL =',&
                        time_step_cfl
                        write(*,'(A,F6.2,A)')'You are using ',percent_deltat,&
                                '% of critical time step.'
                else
                        write(*,'(A,E12.4)')'WARNING!!! deltat CFL =',&
                        time_step_cfl
                        write(*,'(A,F6.2,A)')'You are using ',percent_deltat,&
                                '% of critical time step.'
                        write(*,'(A)')'This could be not enough for the correct working of ABC'
                        if (deltat_fixed.eq.'not') then
                                write(*,'(A)')'deltat chosen will be substituted with 25% of deltat CFL'
                                time_step=time_step_cfl*0.25
                                write(*,'(A,E12.4)')'deltat chosen :',time_step
                        endif
                endif
           elseif (time_step.gt.time_step_cfl) then
                write(*,'(2X,A,E12.4)')'ERROR!!! deltat CFL = ',&
                time_step_cfl,' < deltat = ',time_step
                write(*,'(A)')'The time advancing scheme will be unstable!'
                if (deltat_fixed.eq.'not') then
                        write(*,'(A)')'deltat chosen will be substituted with 25% of deltat CFL'
                        time_step=time_step_cfl*0.25
                        write(*,'(A,E12.4)')'deltat chosen :',time_step
                else
                    b_CFL_failure = .TRUE.
                endif
           endif
           write(*,'(A)')'-------------------------------------------------------'
           write(*,'(A)')' '
        
           if (b_failCFL .and.b_CFL_failure) then
              write(*,*) 'CFL does NOT hold, asked to quit. Program finished.'
              call MPI_FINALIZE(mpi_err)
              call EXIT(EXIT_CFL)
           endif

        !Writing of the results for the maximum number of points per wave length
        
          allocate(ct(nn),ww(nn),dd(nn,nn))
          call MAKE_LGL_NW(nn,ct,ww,dd)

          num_of_points = vs_length_min*(ct(1)-ct(nn))/(ct(int(nn/2))-ct(int(nn/2)+1))/fmax
                       
          write(*,'(A)')'-----------Number of points per wave length------------'
          write(*,'(A,E12.4)')'Max. el. length       :',1/vs_length_min*vs_npoints
          write(*,'(A,E12.4)')'Max. Vs               :',vs_npoints
          write(*,'(A,E12.4)')'Points per wavelength :',num_of_points
          write(*,'(A)') ' '

          deallocate(ct,ww,dd)
        endif   
        
        call MPI_BARRIER(mpi_comm, mpi_err)
        
        return
               
        end subroutine GET_CFL

