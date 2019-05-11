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


!> @brief Set initial conditions.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nnod_loc  number of local nodes
!> @param[in] ne_loc  number of local elements
!> @param[in] length of cs_loc
!> @param[in] cs_loc local connectivity vector
!> @param[in] nm number of materials
!> @param[in] prop_mat material properties (rho,lambda,mu,gamma)
!> @param[in] sdeg_mat polinomial degree vector
!> @param[in] local_n_num
!> @param[in] xs_loc x-coordinate of spectral nodes
!> @param[in] ys_loc y-coordinate of spectral nodes
!> @param[in] zs_loc z-coordinate of spectral nodes
!> @param[in] mpi_id  MPI process id
!> @param[in] dt time step 
!> @param[out] u0 displacement at time -dt
!> @param[out] u1 displacement at time 0
!> @param[out] v1 velocity at time 0  
!  
!************************************************************************************************** 


     subroutine GET_INITIAL_CONDITIONS(nnod_loc, u0, u1, v1, ne_loc, cs_loc, cs_nnz_loc,&
                                       nm,prop_mat, sdeg_mat, loc_n_num, &
                                       xs_loc,ys_loc,zs_loc,&
                                       mpi_id,dt)
     
     
     implicit none
     
     character*100000 :: input_line
     character*4 :: keyword

     integer*4 :: nnod_loc, mpi_id, cs_nnz_loc
     integer*4 :: nm,  ne_loc, status
     integer*4 :: ie, im, nn, fn,is,in,iaz,i,j,k,ileft,iright,nval
     integer*4, dimension(0:cs_nnz_loc) :: cs_loc

     integer*4, dimension(nnod_loc) :: loc_n_num
     integer*4, dimension(nm) :: sdeg_mat
     
     
     real*8 :: time, term1,term2,term3, dt, pi
     real*8 :: a1,a2,a3,b1,b2,b3,c1,c2,c3,w1, rho,lambda,mu, cost1,cost2,cost3
     real*8 :: dxdx,dxdy,dxdz,dydx,dydy,dydz,dzdx,dzdy,dzdz,det_j, u1_ex,u2_ex,u3_ex

     real*8, dimension(nm,4) :: prop_mat
     real*8, dimension(ne_loc) :: alfa11,alfa12,alfa13,alfa21,alfa22,alfa23,alfa31,alfa32,alfa33
     real*8, dimension(ne_loc) :: beta11,beta12,beta13,beta21,beta22,beta23,beta31,beta32,beta33
     real*8, dimension(ne_loc) :: gamma1,gamma2,gamma3,delta1,delta2,delta3
     real*8, dimension(3*nnod_loc) :: u0, u1, v1
     real*8, dimension(nnod_loc) :: xs_loc, ys_loc, zs_loc
     real*8, dimension(:), allocatable :: val_data
 
     
     
     pi = 4.d0*datan(1.d0)

!     open(file='INITCOND.input',unit=40+mpi_id)
!     
!     do 
!         read(40+mpi_id,'(A)',IOSTAT = status) input_line
!         
!         if (status.ne.0) exit
!         
!         keyword = input_line(1:4)
!         
!         ileft = 0
!         iright = len(input_line)
!         do i = 1,iright
!            if (input_line(i:i).eq.' ') exit
!         enddo
!         ileft = i
!         
!         if (keyword .eq. 'DATA') then
!            read(input_line(ileft:iright),*) nval
!         endif   
!      enddo
!
!     close(40+mpi_id)
!     allocate(val_data(nval))
!
!     open(file='INITCOND.input',unit=40+mpi_id)
!     
!     do 
!         read(40+mpi_id,'(A)',IOSTAT = status) input_line
!         
!         if (status.ne.0) exit
!         
!         keyword = input_line(1:4)
!         
!         ileft = 0
!         iright = len(input_line)
!         do i = 1,iright
!            if (input_line(i:i).eq.' ') exit
!         enddo
!         ileft = i
!         
!         if (keyword .eq. 'DATA') then
!            read(input_line(ileft:iright),*) nval, (val_data(j), j = 1, nval)
!         endif   
!      enddo
!
!     close(40+mpi_id)


     do ie = 1,ne_loc
 
                im = cs_loc(cs_loc(ie -1)) 
                nn = sdeg_mat(im) +1
!                a1 = val_data(1); b1 = val_data(2); c1 = val_data(3)
!                a2 = val_data(4); b2 = val_data(5); c2 = val_data(6)  
!                a3 = val_data(7); b3 = val_data(8); c3 = val_data(9)  
!                w1 = val_data(10);
                                            
                do k = 1,nn
                   do j = 1,nn
                      do i = 1,nn
                                 
                         is = nn*nn*(k -1) +nn*(j -1) +i
                         in = cs_loc(cs_loc(ie -1) +is)
                         iaz = 3*(in -1) + 1
             
             !PAPER WITH BLANCA INITIAL CONDITIONS
            u0(iaz+0) = + dsin(3.d0*pi*dt)*(dsin(pi*xs_loc(in)))**2.d0*dsin(2.d0*pi*ys_loc(in))*dsin(2.d0*pi*zs_loc(in));
            u0(iaz+1) = - dsin(3.d0*pi*dt)*(dsin(pi*ys_loc(in)))**2.d0*dsin(2.d0*pi*xs_loc(in))*dsin(2.d0*pi*zs_loc(in));
            u0(iaz+2) = - dsin(3.d0*pi*dt)*(dsin(pi*zs_loc(in)))**2.d0*dsin(2.d0*pi*xs_loc(in))*dsin(2.d0*pi*ys_loc(in));

            u1(iaz+0) = - 0.d0*(dsin(pi*xs_loc(in)))**2.d0*dsin(2.d0*pi*ys_loc(in))*dsin(2.d0*pi*zs_loc(in));
            u1(iaz+1) = + 0.d0*(dsin(pi*ys_loc(in)))**2.d0*dsin(2.d0*pi*xs_loc(in))*dsin(2.d0*pi*zs_loc(in));
            u1(iaz+2) = + 0.d0*(dsin(pi*zs_loc(in)))**2.d0*dsin(2.d0*pi*xs_loc(in))*dsin(2.d0*pi*ys_loc(in));

            v1(iaz+0) = - 3.d0*pi*(dsin(pi*xs_loc(in)))**2.d0*dsin(2.d0*pi*ys_loc(in))*dsin(2.d0*pi*zs_loc(in));
            v1(iaz+1) = + 3.d0*pi*(dsin(pi*ys_loc(in)))**2.d0*dsin(2.d0*pi*xs_loc(in))*dsin(2.d0*pi*zs_loc(in));
            v1(iaz+2) = + 3.d0*pi*(dsin(pi*zs_loc(in)))**2.d0*dsin(2.d0*pi*xs_loc(in))*dsin(2.d0*pi*ys_loc(in));


             !STATIONARY POLYNOMIAL SOLUTION INITIAL CONDITIONS
             !u0(iaz+0) = -dsin(3.d0*pi*dt)*((xs_loc(in)*(xs_loc(in)-1.d0)) &
             !             *(ys_loc(in)*(ys_loc(in)-1.d0)) &
             !             *(zs_loc(in)*(zs_loc(in)-1.d0)));
             !u0(iaz+1) = 0.d0;
             !u0(iaz+2) = 0.d0;
             !
             !u1(iaz+0) = dsin(3.d0*pi*0.d0)*((xs_loc(in)*(xs_loc(in)-1.d0)) &
             !             *(ys_loc(in)*(ys_loc(in)-1.d0)) &
             !             *(zs_loc(in)*(zs_loc(in)-1.d0)));
             !u1(iaz+1) = + 0.d0
             !u1(iaz+2) = + 0.d0
             !
             !v1(iaz+0) = 3.d0*pi*((xs_loc(in)*(xs_loc(in)-1.d0)) &
             !             *(ys_loc(in)*(ys_loc(in)-1.d0)) &
             !             *(zs_loc(in)*(zs_loc(in)-1.d0)));
             !
             !v1(iaz+1) = + 0.d0
             !v1(iaz+2) = + 0.d0


            
                      enddo            
                   enddo
                enddo
     
     enddo     
     
     
     ! deallocate(val_data)
     
     end subroutine GET_INITIAL_CONDITIONS
