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

!> @brief Computes norms of error in test mode case.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nnod_loc number of local nodes
!> @param[in] u1 current displacement
!> @param[in] v1 current velocity
!> @param[in] time discrete instant time 
!> @param[in] ne_loc number of local elements
!> @param[in] cs_loc spectral local connectivity
!> @param[in] cs_nnz_loc length of cs_loc
!> @param[in] nm number of materials
!> @param[in] prop_mat material properties (rho, lambda, mu, gamma)
!> @param[in] sdeg_mat polynomial degree vector
!> @param[in] tag_mat label for material properties
!> @param[in] loc_n_num  local node numbering (local to global)
!> @param[in] local_el_num local element numbering (local to global) 
!> @param[in] nelem_dg_local number of DG elements 
!> @param[in] alfa11 costant value for the bilinear map from (-1,1)^3 to the current element
!> @param[in] alfa12 costant value for the bilinear map from (-1,1)^3 to the current element
!> @param[in] alfa13 costant value for the bilinear map from (-1,1)^3 to the current element
!> @param[in] alfa21 costant value for the bilinear map from (-1,1)^3 to the current element
!> @param[in] alfa22 costant value for the bilinear map from (-1,1)^3 to the current element
!> @param[in] alfa23 costant value for the bilinear map from (-1,1)^3 to the current element
!> @param[in] alfa31 costant value for the bilinear map from (-1,1)^3 to the current element
!> @param[in] alfa32 costant value for the bilinear map from (-1,1)^3 to the current element
!> @param[in] alfa33 costant value for the bilinear map from (-1,1)^3 to the current element 
!> @param[in] beta11 costant value for the bilinear map from (-1,1)^3 to the current element
!> @param[in] beta12 costant value for the bilinear map from (-1,1)^3 to the current element
!> @param[in] beta13 costant value for the bilinear map from (-1,1)^3 to the current element
!> @param[in] beta21 costant value for the bilinear map from (-1,1)^3 to the current element
!> @param[in] beta22 costant value for the bilinear map from (-1,1)^3 to the current element
!> @param[in] beta23 costant value for the bilinear map from (-1,1)^3 to the current element
!> @param[in] beta31 costant value for the bilinear map from (-1,1)^3 to the current element
!> @param[in] beta32 costant value for the bilinear map from (-1,1)^3 to the current element
!> @param[in] beta33 costant value for the bilinear map from (-1,1)^3 to the current element
!> @param[in] gamma1 costant value for the bilinear map from (-1,1)^3 to the current element
!> @param[in] gamma2 costant value for the bilinear map from (-1,1)^3 to the current element
!> @param[in] gamma3 costant value for the bilinear map from (-1,1)^3 to the current element
!> @param[in] delta1 costant value for the bilinear map from (-1,1)^3 to the current element
!> @param[in] delta2 costant value for the bilinear map from (-1,1)^3 to the current element
!> @param[in] delta3 costant value for the bilinear map from (-1,1)^3 to the current element
!> @param[in] xs_loc x-coordinate spectral node 
!> @param[in] ys_loc y-coordinate spectral node
!> @param[in] zs_loc z-coordinate spectral node
!> @param[in] mpi_id MPI id for process
!> @param[in] el_new struct for DG elements
!> @param[in] nelem_dg_global number of total DG elemnents
!> @param[in] nsd_jump number of node to send for jumps
!> @param[in] node_sd_jump list of node to send for jums
!> @param[in] nrv_jump number of node to receive for jumps
!> @param[in] node_rv_jump list of node to receive for jums
!> @param[in] send_length_jump number of node to send for jumps to each MPI process
!> @param[in] recv_length_jump number of node to receive for jumps to each MPI process
!> @param[in] mpi_np  number of MPI process
!> @param[in] mpi_comm MPI common world
!> @param[in] cs_nnz_dg length of cs_dg 
!> @param[in] cs_dg spectral connectivity for DG elements 
!> @param[out] ... Write EN.ERR file containing the computed error norms. 

!> @warning Only for serial computation. Need to be tested for parallel computation

     subroutine COMPUTE_ENERGY_ERROR(nnod_loc, u1, v1, time, ne_loc, cs_loc, cs_nnz_loc,&
                              nm,prop_mat, sdeg_mat,tag_mat, &
                              loc_n_num, local_el_num, nelem_dg_local, &
                              alfa11,alfa12,alfa13,&
                              alfa21,alfa22,alfa23,&
                              alfa31,alfa32,alfa33,& 
                              beta11,beta12,beta13,&
                              beta21,beta22,beta23,&
                              beta31,beta32,beta33,&
                              gamma1,gamma2,gamma3,&
                              delta1,delta2,delta3,&
                              xs_loc,ys_loc,zs_loc,&
                              mpi_id,el_new, & 
                              nelem_dg_global, &
                              nsd_jump,node_sd_jump, &
                              nrv_jump,node_rv_jump, &
                              send_length_jump, recv_length_jump, &  
                              mpi_np, mpi_comm, &
                              cs_nnz_dg, cs_dg)
 

     use max_var
     use DGJUMP
     use speed_exit_codes
     
     implicit none

     include 'SPEED.MPI'
      
     type(el4loop), dimension(ne_loc), intent(in) :: el_new 
     
     character*100000 :: input_line
     character*4 :: keyword
     character*70 :: file_mpi
     
     integer*4 :: nnod_loc, nfunc, mpi_id, cs_nnz_loc, nelem_dg_global, nsd_jump, nrv_jump 
     integer*4 :: nm, unit_mpi, ne_loc, status, ic, id, ie_curr, jb,kb,ib,isb
     integer*4 :: mpi_np, mpi_comm, mpi_ierr, trofa, posit, cs_nnz_dg
     integer*4 :: ielem, iene, iene_curr, imne, jshift,mm, nelem_dg_local
     integer*4 :: ie, im, nn, is,in,iaz,i,j,k,ileft,iright,nval
     integer*4, dimension(0:cs_nnz_loc) :: cs_loc
     integer*4, dimension(0:cs_nnz_dg) ::  cs_dg
     integer*4, dimension(nnod_loc) :: loc_n_num
     integer*4, dimension(ne_loc) :: local_el_num
     integer*4, dimension(nsd_jump):: node_sd_jump                     
     integer*4, dimension(nrv_jump):: node_rv_jump
     integer*4, dimension(nm) :: sdeg_mat
     integer*4, dimension(nm) :: tag_mat
     integer*4, dimension(SPEED_STATUS_SIZE) :: mpi_stat
     integer*4, dimension(mpi_np) :: send_length_jump, recv_length_jump     
               
     real*8 :: time, term1,term2,term3, L2_err, EN_err, term_ener,term_l2 
     real*8 :: Linf_err, EN_err_tot, L2_err_tot, Linf_err_tot, xp, yp, zp, phi

     real*8 :: L2_err_vel, H1_err, term_H1, H1_err_tot, term_l2_vel
     real*8 :: Linf_err_vel, L2_err_vel_tot, Linf_err_vel_tot
     real*8 :: e11_ex, e12_ex, e13_ex, e22_ex, e23_ex, e33_ex
     real*8 :: sxx_ex, syy_ex, szz_ex, sxy_ex, szx_ex, syz_ex


     real*8 :: a1,a2,a3,b1,b2,b3,c1,c2,c3,w1, rho,lambda,mu, cost1,cost2,cost3,pi,x,y,z, u1_ctm, u2_ctm, u3_ctm
     real*8 :: dxdx,dxdy,dxdz,dydx,dydy,dydz,dzdx,dzdy,dzdz,det_j, u1_ex,u2_ex,u3_ex, v1_ex, v2_ex, v3_ex

     real*8, dimension(4) :: Linf_array, Linf_array_vel
     real*8, dimension(nm,4) :: prop_mat
     real*8, dimension(ne_loc) :: alfa11,alfa12,alfa13,alfa21,alfa22,alfa23,alfa31,alfa32,alfa33
     real*8, dimension(ne_loc) :: beta11,beta12,beta13,beta21,beta22,beta23,beta31,beta32,beta33
     real*8, dimension(ne_loc) :: gamma1,gamma2,gamma3,delta1,delta2,delta3
     real*8, dimension(3*nnod_loc) :: u1, v1
     real*8, dimension(nnod_loc) :: xs_loc, ys_loc, zs_loc
     real*8, dimension(3*nsd_jump) :: send_buffer_error
     real*8, dimension(3*nrv_jump) :: recv_buffer_error, jump_minus_error 
     
     real*8, dimension(:), allocatable :: ct, ww, ctm, wwm, u_p, u_m, jump_all
     real*8, dimension(:,:), allocatable :: dd, ddm     
     real*8, dimension(:), allocatable :: val_data, basis_functions
     real*8, dimension(:,:,:), allocatable :: ux_el,uy_el,uz_el
     real*8, dimension(:,:,:), allocatable :: duxdx_el,duydx_el,duzdx_el 
     real*8, dimension(:,:,:), allocatable :: duxdy_el,duydy_el,duzdy_el 
     real*8, dimension(:,:,:), allocatable :: duxdz_el,duydz_el,duzdz_el 
     real*8, dimension(:,:,:), allocatable :: sxx_el,syy_el,szz_el
     real*8, dimension(:,:,:), allocatable :: syz_el,szx_el,sxy_el
     real*8, dimension(:,:,:), allocatable :: lambda_el,mu_el 

     pi = 4.d0*datan(1.d0)
     


     L2_err = 0.d0; EN_err = 0.d0; Linf_err = 0.d0;
     L2_err_tot = 0.d0; EN_err_tot = 0.d0; Linf_err_tot = 0.d0;
     
     L2_err_vel = 0.d0;  H1_err = 0.d0;  Linf_err_vel = 0.d0;
     L2_err_vel_tot = 0.d0;  H1_err_tot = 0.d0; Linf_err_vel_tot = 0.d0;

     do ie = 1,ne_loc
 
         im = cs_loc(cs_loc(ie -1)) 
         nn = sdeg_mat(im) +1;
         mm = nn ! 8 ! gll nodes for computing norms
         allocate(ct(nn),ww(nn),dd(nn,nn));  allocate(ctm(mm),wwm(mm),ddm(mm,mm));
         
         allocate(ux_el(nn,nn,nn),uy_el(nn,nn,nn),uz_el(nn,nn,nn))
     	 allocate(duxdx_el(nn,nn,nn),duydx_el(nn,nn,nn),duzdx_el(nn,nn,nn)) 
    	 allocate(duxdy_el(nn,nn,nn),duydy_el(nn,nn,nn),duzdy_el(nn,nn,nn)) 
         allocate(duxdz_el(nn,nn,nn),duydz_el(nn,nn,nn),duzdz_el(nn,nn,nn)) 
     	 allocate(sxx_el(nn,nn,nn),syy_el(nn,nn,nn),szz_el(nn,nn,nn))
      	 allocate(syz_el(nn,nn,nn),szx_el(nn,nn,nn),sxy_el(nn,nn,nn))
      	 allocate(lambda_el(nn,nn,nn),mu_el(nn,nn,nn)) 

         call MAKE_LGL_NW(nn,ct,ww,dd)
         call MAKE_LGL_NW(mm,ctm,wwm,ddm)
 
                      
         rho = prop_mat(im,1); lambda = prop_mat(im,2); mu = prop_mat(im,3)

 
                
         do k = 1,nn
            do j = 1,nn
               do i = 1,nn
                  is = nn*nn*(k -1) +nn*(j -1) +i
                  in = cs_loc(cs_loc(ie -1) + is)

     
                  iaz = 3*(in -1) +1; ux_el(i,j,k) = u1(iaz) 
                  iaz = 3*(in -1) +2; uy_el(i,j,k) = u1(iaz)  
                  iaz = 3*(in -1) +3; uz_el(i,j,k) = u1(iaz) 

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
 		     	  
        

        lambda_el = lambda; mu_el = mu;
        
        
        call MAKE_STRESS_TENSOR(nn,lambda_el,mu_el,&					
                                duxdx_el,duydx_el,duzdx_el,&				
	                        duxdy_el,duydy_el,duzdy_el,&				
				duxdz_el,duydz_el,duzdz_el,&				
                                sxx_el,syy_el,szz_el,&				
				syz_el,szx_el,sxy_el)


         do k = 1,mm
            do j = 1,mm
               do i = 1,mm
                      dxdx = alfa11(ie) +beta12(ie)*ctm(k) &
                             + beta13(ie)*ctm(j) +gamma1(ie)*ctm(j)*ctm(k)
                      dydx = alfa21(ie) +beta22(ie)*ctm(k) &
                             + beta23(ie)*ctm(j) +gamma2(ie)*ctm(j)*ctm(k)
                      dzdx = alfa31(ie) +beta32(ie)*ctm(k) &
                             + beta33(ie)*ctm(j) +gamma3(ie)*ctm(j)*ctm(k)
                              
                      dxdy = alfa12(ie) +beta11(ie)*ctm(k) &
                             + beta13(ie)*ctm(i) +gamma1(ie)*ctm(k)*ctm(i)
                      dydy = alfa22(ie) +beta21(ie)*ctm(k) &
                             + beta23(ie)*ctm(i) +gamma2(ie)*ctm(k)*ctm(i)
                      dzdy = alfa32(ie) +beta31(ie)*ctm(k) &
                             + beta33(ie)*ctm(i) +gamma3(ie)*ctm(k)*ctm(i)
                                 
                      dxdz = alfa13(ie) +beta11(ie)*ctm(j) &
                           + beta12(ie)*ctm(i) +gamma1(ie)*ctm(i)*ctm(j)
                      dydz = alfa23(ie) +beta21(ie)*ctm(j) &
                           + beta22(ie)*ctm(i) +gamma2(ie)*ctm(i)*ctm(j)
                      dzdz = alfa33(ie) +beta31(ie)*ctm(j) &
                           + beta32(ie)*ctm(i) +gamma3(ie)*ctm(i)*ctm(j)
                                 
                      det_j = dxdz * (dydx*dzdy - dzdx*dydy) &
                            - dydz * (dxdx*dzdy - dzdx*dxdy) &
                            + dzdz * (dxdx*dydy - dydx*dxdy)
                                 
                      !Local to global node transformation
                      xp = alfa11(ie)*ctm(i) + alfa12(ie)*ctm(j) &
                           + alfa13(ie)*ctm(k) + beta11(ie)*ctm(j)*ctm(k) &
                           + beta12(ie)*ctm(i)*ctm(k) + beta13(ie)*ctm(i)*ctm(j) &
                           + gamma1(ie)*ctm(i)*ctm(j)*ctm(k) + delta1(ie)
                        
                      yp = alfa21(ie)*ctm(i) + alfa22(ie)*ctm(j) &
                           + alfa23(ie)*ctm(k) + beta21(ie)*ctm(j)*ctm(k) &
                           + beta22(ie)*ctm(i)*ctm(k) + beta23(ie)*ctm(i)*ctm(j) &
                           + gamma2(ie)*ctm(i)*ctm(j)*ctm(k) + delta2(ie)
                        
                      zp = alfa31(ie)*ctm(i) + alfa32(ie)*ctm(j) &
                           + alfa33(ie)*ctm(k) + beta31(ie)*ctm(j)*ctm(k) &
                           + beta32(ie)*ctm(i)*ctm(k) + beta33(ie)*ctm(i)*ctm(j) &
                           + gamma3(ie)*ctm(i)*ctm(j)*ctm(k) + delta3(ie)
                        
                      !STATIONARY POLYNOMIAL SOLUTIONS  
                      !u1_ex = dsin(3.d0*pi*time)*((xp*(xp-1.d0))*(yp*(yp-1.d0))*(zp*(zp-1.d0)));
                      !u2_ex = 0.d0;
                      !u3_ex = 0.d0;


                      is = nn*nn*(k -1) +nn*(j -1) +i
                      in = cs_loc(cs_loc(ie -1) +is)

      !POPLYNOMIAL SOLUTION
      !u1_ex = dsin(3.d0*pi*time)*((xs_loc(in)*(xs_loc(in)-1.d0))**2.d0 &
      !                                      *(ys_loc(in)*(ys_loc(in)-1.d0))**2.d0 &
      !                                      *(zs_loc(in)*(zs_loc(in)-1.d0))**2.d0);
      !u2_ex = dsin(3.d0*pi*time)*((xs_loc(in)*(xs_loc(in)-1.d0))**2.d0 &
      !                                      *(ys_loc(in)*(ys_loc(in)-1.d0))**2.d0 &
      !                                      *(zs_loc(in)*(zs_loc(in)-1.d0))**2.d0);
      !u3_ex = dsin(3.d0*pi*time)*((xs_loc(in)*(xs_loc(in)-1.d0))**2.d0 &
      !                                      *(ys_loc(in)*(ys_loc(in)-1.d0))**2.d0 &
      !                                      *(zs_loc(in)*(zs_loc(in)-1.d0))**2.d0);
      !
      !v1_ex = + 3.d0*pi*dcos(3.d0*pi*time)*((xs_loc(in)*(xs_loc(in)-1.d0))**2.d0 &
      !                                      *(ys_loc(in)*(ys_loc(in)-1.d0))**2.d0 &
      !                                      *(zs_loc(in)*(zs_loc(in)-1.d0))**2.d0);
      !v2_ex = + 3.d0*pi*dcos(3.d0*pi*time)*((xs_loc(in)*(xs_loc(in)-1.d0))**2.d0 &
      !                                      *(ys_loc(in)*(ys_loc(in)-1.d0))**2.d0 &
      !                                      *(zs_loc(in)*(zs_loc(in)-1.d0))**2.d0);
      !v3_ex = + 3.d0*pi*dcos(3.d0*pi*time)*((xs_loc(in)*(xs_loc(in)-1.d0))**2.d0 &
      !                                      *(ys_loc(in)*(ys_loc(in)-1.d0))**2.d0 &
      !                                      *(zs_loc(in)*(zs_loc(in)-1.d0))**2.d0);




      !PAPER WITH BLANCA
      u1_ex = -dsin(3.d0*pi*time)*(dsin(pi*xp))**2.d0*dsin(2.d0*pi*yp)*dsin(2.d0*pi*zp);
      u2_ex = +dsin(3.d0*pi*time)*(dsin(pi*yp))**2.d0*dsin(2.d0*pi*xp)*dsin(2.d0*pi*zp);
      u3_ex = +dsin(3.d0*pi*time)*(dsin(pi*zp))**2.d0*dsin(2.d0*pi*xp)*dsin(2.d0*pi*yp);

      v1_ex = - 3.d0*pi*dcos(3.d0*pi*time)*(dsin(pi*xs_loc(in)))**2.d0*dsin(2.d0*pi*ys_loc(in))*dsin(2.d0*pi*zs_loc(in));
      v2_ex = + 3.d0*pi*dcos(3.d0*pi*time)*(dsin(pi*ys_loc(in)))**2.d0*dsin(2.d0*pi*xs_loc(in))*dsin(2.d0*pi*zs_loc(in));
      v3_ex = + 3.d0*pi*dcos(3.d0*pi*time)*(dsin(pi*zs_loc(in)))**2.d0*dsin(2.d0*pi*xs_loc(in))*dsin(2.d0*pi*ys_loc(in));

      x = xs_loc(in); y = ys_loc(in); z = zs_loc(in);         
      e11_ex = -2.d0*pi*dcos(pi*x)*dsin(3.d0*pi*time)*dsin(pi*x)*dsin(2.d0*pi*y)*dsin(2.d0*pi*z)
      e22_ex =  2.d0*pi*dcos(pi*y)*dsin(3.d0*pi*time)*dsin(2.d0*pi*x)*dsin(pi*y)*dsin(2.d0*pi*z)

      e33_ex =  2.d0*pi*dcos(pi*z)*dsin(3.d0*pi*time)*dsin(2.d0*pi*x)*dsin(2.d0*pi*y)*dsin(pi*z)

      e12_ex =  pi*dcos(2.d0*pi*x)*dsin(3.d0*pi*time)*dsin(pi*y)**2.d0*dsin(2.d0*pi*z) &
                              - pi*dcos(2.d0*pi*y)*dsin(3.d0*pi*time)*dsin(pi*x)**2.d0*dsin(2.d0*pi*z)

      e13_ex = pi*dcos(2.d0*pi*x)*dsin(3.d0*pi*time)*dsin(2.d0*pi*y)*dsin(pi*z)**2.d0 &
                             - pi*dcos(2.d0*pi*z)*dsin(3.d0*pi*time)*dsin(pi*x)**2.d0*dsin(2.d0*pi*y)

      e23_ex = pi*dcos(2.d0*pi*z)*dsin(3.d0*pi*time)*dsin(2.d0*pi*x)*dsin(pi*y)**2.d0 &
                             + pi*dcos(2.d0*pi*y)*dsin(3.d0*pi*time)*dsin(2.d0*pi*x)*dsin(pi*z)**2.d0

                      
      sxx_ex = lambda*(2.d0*pi*dcos(pi*y)*dsin(3.d0*pi*time)*dsin(2.d0*pi*x)*dsin(pi*y)*dsin(2.d0*pi*z) &
                               - 2.d0*pi*dcos(pi*x)*dsin(3.d0*pi*time)*dsin(pi*x)*dsin(2.d0*pi*y)*dsin(2.d0*pi*z) &
                               + 2.d0*pi*dcos(pi*z)*dsin(3.d0*pi*time)*dsin(2.d0*pi*x)*dsin(2.d0*pi*y)*dsin(pi*z)) &
                               - 4.d0*pi*mu*dcos(pi*x)*dsin(3.d0*pi*time)*dsin(pi*x)*dsin(2.d0*pi*y)*dsin(2.d0*pi*z)

      syy_ex = lambda*(2.d0*pi*dcos(pi*y)*dsin(3.d0*pi*time)*dsin(2.d0*pi*x)*dsin(pi*y)*dsin(2.d0*pi*z) &
                               - 2.d0*pi*dcos(pi*x)*dsin(3.d0*pi*time)*dsin(pi*x)*dsin(2.d0*pi*y)*dsin(2.d0*pi*z) &
                               + 2.d0*pi*dcos(pi*z)*dsin(3.d0*pi*time)*dsin(2.d0*pi*x)*dsin(2.d0*pi*y)*dsin(pi*z)) &
                               + 4.d0*pi*mu*dcos(pi*y)*dsin(3.d0*pi*time)*dsin(2.d0*pi*x)*dsin(pi*y)*dsin(2.d0*pi*z)

      szz_ex = lambda*(2*pi*dcos(pi*y)*dsin(3.d0*pi*time)*dsin(2.d0*pi*x)*dsin(pi*y)*dsin(2.d0*pi*z) &
                               - 2.d0*pi*dcos(pi*x)*dsin(3.d0*pi*time)*dsin(pi*x)*dsin(2.d0*pi*y)*dsin(2.d0*pi*z) &
                               + 2.d0*pi*dcos(pi*z)*dsin(3.d0*pi*time)*dsin(2.d0*pi*x)*dsin(2.d0*pi*y)*dsin(pi*z)) &
                               + 4.d0*pi*mu*dcos(pi*z)*dsin(3.d0*pi*time)*dsin(2.d0*pi*x)*dsin(2.d0*pi*y)*dsin(pi*z)
 
      sxy_ex = mu*(2.d0*pi*dcos(2.d0*pi*x)*dsin(3.d0*pi*time)*dsin(pi*y)**2.d0*dsin(2.d0*pi*z) &
                               - 2.d0*pi*dcos(2.d0*pi*y)*dsin(3.d0*pi*time)*dsin(pi*x)**2.d0*dsin(2.d0*pi*z))

      szx_ex = mu*(2.d0*pi*dcos(2.d0*pi*x)*dsin(3.d0*pi*time)*dsin(2.d0*pi*y)*dsin(pi*z)**2.d0 &
                               - 2.d0*pi*dcos(2.d0*pi*z)*dsin(3.d0*pi*time)*dsin(pi*x)**2.d0*dsin(2.d0*pi*y))

      syz_ex = mu*(2.d0*pi*dcos(2.d0*pi*z)*dsin(3.d0*pi*time)*dsin(2.d0*pi*x)*dsin(pi*y)**2.d0 &
                               + 2.d0*pi*dcos(2.d0*pi*y)*dsin(3.d0*pi*time)*dsin(2.d0*pi*x)*dsin(pi*z)**2.d0)


                      u1_ctm = 0.d0; u2_ctm = 0.d0; u3_ctm = 0.d0
                      
                      allocate(basis_functions(mm**3)) ! basis function evaluation
                      do kb = 1, nn
                         do jb = 1, nn
                            do ib =1, nn
                               isb = nn*nn*(kb -1) + nn*(jb -1) +ib
                               in = cs_loc(cs_loc(ie -1) +isb)
                               iaz = 3*(in -1) + 1
                               
                               basis_functions(isb) = phi(nn-1, ct(ib), ctm(i)) * &
                                                      phi(nn-1, ct(jb), ctm(j)) * & 
                                                      phi(nn-1, ct(kb), ctm(k))

                               u1_ctm =  u1_ctm +  basis_functions(isb)*u1(iaz);
                               u2_ctm =  u2_ctm +  basis_functions(isb)*u1(iaz+1);
                               u3_ctm =  u3_ctm +  basis_functions(isb)*u1(iaz+2);
                                                                                    
                            enddo
                         enddo
                      enddo                            
                      deallocate(basis_functions)


                      
                    

                      !L2 NORM OLD
                      is = nn*nn*(k -1) +nn*(j -1) +i
                      in = cs_loc(cs_loc(ie -1) +is)

                      iaz = 3*(in -1) + 1

                      term_l2_vel =  det_j * ww(i) * ww(j) * ww(k) * &
                        ( (v1_ex - v1(iaz))**2.d0 + (v2_ex - v1(iaz+1))**2.d0 +(v3_ex - v1(iaz+2))**2.d0)   

                      
                      L2_err_vel = L2_err_vel + term_l2_vel 
                      
                      

                      !LINF NORM - DISPL & VEL
                      Linf_array = (/ Linf_err,abs(u1_ex - u1_ctm),abs(u2_ex - u2_ctm),abs(u3_ex - u3_ctm) /) 
                      Linf_err = maxval(Linf_array)


                      !L2 NORM NEW
                      term_l2 =  det_j * wwm(i) * wwm(j) * wwm(k) * &
                        ( (u1_ex - u1_ctm)**2.d0 + (u2_ex - u2_ctm)**2.d0 +(u3_ex - u3_ctm)**2.d0)                      
                      L2_err = L2_err + term_l2 



                      !ENERGY NORM
                      iaz = 3*(in -1) + 1
                      term_ener =  det_j * ww(i) * ww(j) * ww(k) * ( &
                                   rho * (v1_ex - v1(iaz+0))**2.d0 &
                                 + rho * (v2_ex - v1(iaz+1))**2.d0 & 
                                 + rho * (v3_ex - v1(iaz+2))**2.d0 &
                                 + (sxx_el(i,j,k) - sxx_ex) * (duxdx_el(i,j,k) - e11_ex) &
                                 + (syy_el(i,j,k) - syy_ex) * (duydy_el(i,j,k) - e22_ex) &
                                 + (szz_el(i,j,k) - szz_ex) * (duzdz_el(i,j,k) - e33_ex) &
                                 + 2.d0*(sxy_el(i,j,k) - sxy_ex) * (0.5d0*duxdy_el(i,j,k) + 0.5d0*duydx_el(i,j,k) - e12_ex) &
                                 + 2.d0*(szx_el(i,j,k) - szx_ex) * (0.5d0*duxdz_el(i,j,k) + 0.5d0*duzdx_el(i,j,k) - e13_ex) & 
                                 + 2.d0*(syz_el(i,j,k) - syz_ex) * (0.5d0*duydz_el(i,j,k) + 0.5d0*duzdy_el(i,j,k) - e23_ex) &
                                )
                                  
                                  
                      EN_err = EN_err + term_ener 
                      
                      term_H1 =  det_j * ww(i) * ww(j) * ww(k) * ( &
                                 + (sxx_el(i,j,k) - sxx_ex) * (duxdx_el(i,j,k) - e11_ex) &
                                 + (syy_el(i,j,k) - syy_ex) * (duydy_el(i,j,k) - e22_ex) &
                                 + (szz_el(i,j,k) - szz_ex) * (duzdz_el(i,j,k) - e33_ex) &
                                 + 2.d0*(sxy_el(i,j,k) - sxy_ex) * (0.5d0*duxdy_el(i,j,k) + 0.5d0*duydx_el(i,j,k) - e12_ex) &
                                 + 2.d0*(szx_el(i,j,k) - szx_ex) * (0.5d0*duxdz_el(i,j,k) + 0.5d0*duzdx_el(i,j,k) - e13_ex) & 
                                 + 2.d0*(syz_el(i,j,k) - syz_ex) * (0.5d0*duydz_el(i,j,k) + 0.5d0*duzdy_el(i,j,k) - e23_ex) &
                                )

                      H1_err = H1_err + term_H1
            
                      enddo            
                   enddo
                enddo
               

         deallocate(ct,ww,dd, ctm, wwm, ddm)
         deallocate(ux_el,uy_el,uz_el)
     	 deallocate(duxdx_el,duydx_el,duzdx_el,duxdy_el,duydy_el,duzdy_el,duxdz_el,duydz_el,duzdz_el) 
     	 deallocate(sxx_el,syy_el,szz_el,syz_el,szx_el,sxy_el)
         deallocate(lambda_el,mu_el) 

     
     enddo     
     
    
     
!----------------------------------------------------------------------------------------------------------
!    NOW COMPUTE NON CONFORMING TERM S_e int_e [u-uh][u-uh] 
!----------------------------------------------------------------------------------------------------------     
     if(nelem_dg_global .gt. 0) then

          jump_minus_error = 0.d0

          do i = 1,nsd_jump
             in = node_sd_jump(i)
                  !PAPER WITH BLANCA 
                  u1_ex = -dsin(3.d0*pi*time)*(dsin(pi*xs_loc(in)))**2.d0*dsin(2.d0*pi*ys_loc(in))*dsin(2.d0*pi*zs_loc(in));
                  u2_ex = +dsin(3.d0*pi*time)*(dsin(pi*ys_loc(in)))**2.d0*dsin(2.d0*pi*xs_loc(in))*dsin(2.d0*pi*zs_loc(in));
                  u3_ex = +dsin(3.d0*pi*time)*(dsin(pi*zs_loc(in)))**2.d0*dsin(2.d0*pi*xs_loc(in))*dsin(2.d0*pi*ys_loc(in));

 
                  send_buffer_error(3*(i -1) +1) = abs(u1(3*(in -1) +1) - u1_ex)            
                  send_buffer_error(3*(i -1) +2) = abs(u1(3*(in -1) +2) - u2_ex)
                  send_buffer_error(3*(i -1) +3) = abs(u1(3*(in -1) +3) - u3_ex)
             
          enddo
      
         call EXCHANGE_DOUBLE(3*nsd_jump, send_buffer_error, 3*nrv_jump, recv_buffer_error,&
                                mpi_np, send_length_jump, recv_length_jump,&
                                mpi_comm, mpi_stat, mpi_ierr,mpi_id)
     
          do i = 1,nrv_jump
             in = node_rv_jump(i)
             jump_minus_error(3*(i -1) +1) = recv_buffer_error(3*(in -1) +1)           
             jump_minus_error(3*(i -1) +2) = recv_buffer_error(3*(in -1) +2)
             jump_minus_error(3*(i -1) +3) = recv_buffer_error(3*(in -1) +3)
          enddo

      endif


     
      do ie = 1, nelem_dg_local

        ielem = el_new(ie)%ind; nn = el_new(ie)%deg
        call GET_INDLOC_FROM_INDGLO(local_el_num, ne_loc, ielem, ie_curr)

        allocate(u_p(3*nn**3)); u_p = 0.d0; 

               do k = 1,nn
                  do j = 1,nn
                     do i = 1,nn
                        is = nn*nn*(k -1) +nn*(j -1) +i
                        in = cs_loc(cs_loc(ie_curr -1) + is)
                  
                  !PAPER WITH BLANCA        
                  u1_ex = -dsin(3.d0*pi*time)*(dsin(pi*xs_loc(in)))**2.d0*dsin(2.d0*pi*ys_loc(in))*dsin(2.d0*pi*zs_loc(in));
                  u2_ex = +dsin(3.d0*pi*time)*(dsin(pi*ys_loc(in)))**2.d0*dsin(2.d0*pi*xs_loc(in))*dsin(2.d0*pi*zs_loc(in));
                  u3_ex = +dsin(3.d0*pi*time)*(dsin(pi*zs_loc(in)))**2.d0*dsin(2.d0*pi*xs_loc(in))*dsin(2.d0*pi*ys_loc(in));

            
                        iaz = 3*(in -1) +1; u_p(is) = abs(u1(iaz) - u1_ex)
                        iaz = 3*(in -1) +2; u_p(nn**3+is) = abs(u1(iaz) - u2_ex)
                        iaz = 3*(in -1) +3; u_p(2*nn**3 + is) = abs(u1(iaz) - u3_ex)

                        
                     enddo
                  enddo
                enddo


               allocate(u_m(el_new(ie)%nnz_col_only_uv)); u_m = 0.d0; jshift = 0 
                
               do ic = 1, el_new(ie)%num_of_ne 
                  
                   iene = el_new(ie)%el_conf(ic,1)       
                   imne = el_new(ie)%el_conf(ic,0)
                       
                   mm = 2
                   do i = 1, nm 
                      if(tag_mat(i) .eq. imne ) then
                         mm = sdeg_mat(i) +1
                      endif
                   enddo   

                   trofa  = 0 
                   if(mpi_np .gt. 1) call CHECK_MPI(cs_nnz_dg, cs_dg, iene, trofa, posit)

                   if(trofa .eq. 1) then
                     do k = 1,mm
                        do j = 1,mm
                           do i = 1,mm
                              is =  mm*mm*(k -1) +mm*(j -1) +i
                              id = cs_dg(cs_dg(posit-1)+is)
                              iaz = 3*(id -1) +1; u_m(is + jshift) = jump_minus_error(iaz)
                              iaz = 3*(id -1) +2; u_m(mm**3 + is + jshift) = jump_minus_error(iaz)
                              iaz = 3*(id -1) +3; u_m(2*(mm**3) + is + jshift) = jump_minus_error(iaz)  
                            enddo
                         enddo   
                      enddo
                   else                                     

                     call GET_INDLOC_FROM_INDGLO(local_el_num, ne_loc, iene, iene_curr)

                     do k = 1,mm
                        do j = 1,mm
                           do i = 1,mm
                              is =  mm*mm*(k -1) +mm*(j -1) +i
                              id = cs_loc(cs_loc(iene_curr-1)+is)

                  !PAPER WITH BLANCA  
                  u1_ex = -dsin(3.d0*pi*time)*(dsin(pi*xs_loc(id)))**2.d0*dsin(2.d0*pi*ys_loc(id))*dsin(2.d0*pi*zs_loc(id));
                  u2_ex = +dsin(3.d0*pi*time)*(dsin(pi*ys_loc(id)))**2.d0*dsin(2.d0*pi*xs_loc(id))*dsin(2.d0*pi*zs_loc(id));
                  u3_ex = +dsin(3.d0*pi*time)*(dsin(pi*zs_loc(id)))**2.d0*dsin(2.d0*pi*xs_loc(id))*dsin(2.d0*pi*ys_loc(id));


                              iaz = 3*(id -1) +1; u_m(is + jshift) = abs(u1(iaz) - u1_ex)
                              iaz = 3*(id -1) +2; u_m(mm**3 + is + jshift) = abs(u1(iaz) - u2_ex)
                              iaz = 3*(id -1) +3; u_m(2*(mm**3) + is + jshift) = abs(u1(iaz) - u3_ex)  

 
                           enddo
                        enddo
                      enddo
                    endif

                      !  write(*,*) '--------------- fine +- ------------------------------'
                      !  read(*,*)


                   jshift = jshift + 3*mm**3 

                enddo

                allocate(jump_all(3*nn**3))
                jump_all = 0.d0


                call MATMUL_SPARSE(el_new(ie)%matMin_only_uv, el_new(ie)%nnz_minus_only_uv, el_new(ie)%JMin_only_uv, &
                                  el_new(ie)%IMin_only_uv, jump_all, 3*nn**3,  &
                                  u_m, el_new(ie)%nnz_col_only_uv,1) 
                
                    
                EN_err = EN_err + dot_product(jump_all,u_p)
                H1_err = H1_err + dot_product(jump_all,u_p)
               

                deallocate(u_m, jump_all)
                allocate(jump_all(3*nn**3)); jump_all = 0.d0   

                call MATMUL_SPARSE(el_new(ie)%matPlus_only_uv, el_new(ie)%nnz_plus_only_uv,el_new(ie)%JPlus_only_uv, &
                                  el_new(ie)%IPlus_only_uv, jump_all, 3*nn**3, &
                                  u_p, 3*nn**3,1) 
                
               
                EN_err = EN_err + dot_product(jump_all,u_p)
                H1_err = H1_err + dot_product(jump_all,u_p)


                 deallocate(u_p, jump_all)

      enddo
      
      
!----------------------------------------------------------------------------------------------------------
!    NOW COMPUTE NON CONFORMING TERM S_e int_e [v-vh][v-vh] v=du/dt
!----------------------------------------------------------------------------------------------------------     
     if(nelem_dg_global .gt. 0) then

          jump_minus_error = 0.d0

          do i = 1,nsd_jump
             in = node_sd_jump(i)
                  !PAPER WITH BLANCA 
      v1_ex = - 3.d0*pi*dcos(3.d0*pi*time)*(dsin(pi*xs_loc(in)))**2.d0*dsin(2.d0*pi*ys_loc(in))*dsin(2.d0*pi*zs_loc(in));
      v2_ex = + 3.d0*pi*dcos(3.d0*pi*time)*(dsin(pi*ys_loc(in)))**2.d0*dsin(2.d0*pi*xs_loc(in))*dsin(2.d0*pi*zs_loc(in));
      v3_ex = + 3.d0*pi*dcos(3.d0*pi*time)*(dsin(pi*zs_loc(in)))**2.d0*dsin(2.d0*pi*xs_loc(in))*dsin(2.d0*pi*ys_loc(in));
                   
 
                  send_buffer_error(3*(i -1) +1) = abs(v1(3*(in -1) +1) - v1_ex)            
                  send_buffer_error(3*(i -1) +2) = abs(v1(3*(in -1) +2) - v2_ex)
                  send_buffer_error(3*(i -1) +3) = abs(v1(3*(in -1) +3) - v3_ex)
             
          enddo
      
         call EXCHANGE_DOUBLE(3*nsd_jump, send_buffer_error, 3*nrv_jump, recv_buffer_error,&
                                mpi_np, send_length_jump, recv_length_jump,&
                                mpi_comm, mpi_stat, mpi_ierr,mpi_id)
     
          do i = 1,nrv_jump
             in = node_rv_jump(i)
             jump_minus_error(3*(i -1) +1) = recv_buffer_error(3*(in -1) +1)           
             jump_minus_error(3*(i -1) +2) = recv_buffer_error(3*(in -1) +2)
             jump_minus_error(3*(i -1) +3) = recv_buffer_error(3*(in -1) +3)
          enddo

      endif

     
      do ie = 1, nelem_dg_local

        ielem = el_new(ie)%ind; nn = el_new(ie)%deg
        call GET_INDLOC_FROM_INDGLO(local_el_num, ne_loc, ielem, ie_curr)

        allocate(u_p(3*nn**3)); u_p = 0.d0; 

               do k = 1,nn
                  do j = 1,nn
                     do i = 1,nn
                        is = nn*nn*(k -1) +nn*(j -1) +i
                        in = cs_loc(cs_loc(ie_curr -1) + is)
                  
                  !PAPER WITH BLANCA        
      v1_ex = - 3.d0*pi*dcos(3.d0*pi*time)*(dsin(pi*xs_loc(in)))**2.d0*dsin(2.d0*pi*ys_loc(in))*dsin(2.d0*pi*zs_loc(in));
      v2_ex = + 3.d0*pi*dcos(3.d0*pi*time)*(dsin(pi*ys_loc(in)))**2.d0*dsin(2.d0*pi*xs_loc(in))*dsin(2.d0*pi*zs_loc(in));
      v3_ex = + 3.d0*pi*dcos(3.d0*pi*time)*(dsin(pi*zs_loc(in)))**2.d0*dsin(2.d0*pi*xs_loc(in))*dsin(2.d0*pi*ys_loc(in));

            
                        iaz = 3*(in -1) +1; u_p(is) = abs(v1(iaz) - v1_ex)
                        iaz = 3*(in -1) +2; u_p(nn**3+is) = abs(v1(iaz) - v2_ex)
                        iaz = 3*(in -1) +3; u_p(2*nn**3 + is) = abs(v1(iaz) - v3_ex)

                        
                     enddo
                  enddo
                enddo


               allocate(u_m(el_new(ie)%nnz_col_only_uv)); u_m = 0.d0; jshift = 0 
                
               do ic = 1, el_new(ie)%num_of_ne 
                  
                   iene = el_new(ie)%el_conf(ic,1)       
                   imne = el_new(ie)%el_conf(ic,0)
                       
                   mm = 2
                   do i = 1, nm 
                      if(tag_mat(i) .eq. imne ) then
                         mm = sdeg_mat(i) +1
                      endif
                   enddo   

                   trofa  = 0 
                   if(mpi_np .gt. 1) call CHECK_MPI(cs_nnz_dg, cs_dg, iene, trofa, posit)

                   if(trofa .eq. 1) then
                     do k = 1,mm
                        do j = 1,mm
                           do i = 1,mm
                              is =  mm*mm*(k -1) +mm*(j -1) +i
                              id = cs_dg(cs_dg(posit-1)+is)
                              iaz = 3*(id -1) +1; u_m(is + jshift) = jump_minus_error(iaz)
                              iaz = 3*(id -1) +2; u_m(mm**3 + is + jshift) = jump_minus_error(iaz)
                              iaz = 3*(id -1) +3; u_m(2*(mm**3) + is + jshift) = jump_minus_error(iaz)  
                            enddo
                         enddo   
                      enddo
                   else                                     

                     call GET_INDLOC_FROM_INDGLO(local_el_num, ne_loc, iene, iene_curr)

                     do k = 1,mm
                        do j = 1,mm
                           do i = 1,mm
                              is =  mm*mm*(k -1) +mm*(j -1) +i
                              id = cs_loc(cs_loc(iene_curr-1)+is)

                  !PAPER WITH BLANCA  
      v1_ex = - 3.d0*pi*dcos(3.d0*pi*time)*(dsin(pi*xs_loc(id)))**2.d0*dsin(2.d0*pi*ys_loc(id))*dsin(2.d0*pi*zs_loc(id));
      v2_ex = + 3.d0*pi*dcos(3.d0*pi*time)*(dsin(pi*ys_loc(id)))**2.d0*dsin(2.d0*pi*xs_loc(id))*dsin(2.d0*pi*zs_loc(id));
      v3_ex = + 3.d0*pi*dcos(3.d0*pi*time)*(dsin(pi*zs_loc(id)))**2.d0*dsin(2.d0*pi*xs_loc(id))*dsin(2.d0*pi*ys_loc(id));



                              iaz = 3*(id -1) +1; u_m(is + jshift) = abs(v1(iaz) - v1_ex)
                              iaz = 3*(id -1) +2; u_m(mm**3 + is + jshift) = abs(v1(iaz) - v2_ex)
                              iaz = 3*(id -1) +3; u_m(2*(mm**3) + is + jshift) = abs(v1(iaz) - v3_ex)  

 
                           enddo
                        enddo
                      enddo
                    endif

                      !  write(*,*) '--------------- fine +- ------------------------------'
                      !  read(*,*)


                   jshift = jshift + 3*mm**3 

                enddo

                allocate(jump_all(3*nn**3))
                jump_all = 0.d0


                call MATMUL_SPARSE(el_new(ie)%matMin_only_uv, el_new(ie)%nnz_minus_only_uv, el_new(ie)%JMin_only_uv, &
                                  el_new(ie)%IMin_only_uv, jump_all, 3*nn**3,  &
                                  u_m, el_new(ie)%nnz_col_only_uv,1) 
                
                    
                EN_err = EN_err + dot_product(jump_all,u_p)
                !write(*,*) dot_product(jump_all,u_p)                  
                !H1_err = H1_err + dot_product(jump_all,u_p)
               

                deallocate(u_m, jump_all)
                allocate(jump_all(3*nn**3)); jump_all = 0.d0   

                call MATMUL_SPARSE(el_new(ie)%matPlus_only_uv, el_new(ie)%nnz_plus_only_uv,el_new(ie)%JPlus_only_uv, &
                                  el_new(ie)%IPlus_only_uv, jump_all, 3*nn**3, &
                                  u_p, 3*nn**3,1) 
                
               
                EN_err = EN_err + dot_product(jump_all,u_p)
                !write(*,*) dot_product(jump_all,u_p)
                !read(*,*)
                !H1_err = H1_err + dot_product(jump_all,u_p)


                 deallocate(u_p, jump_all)

      enddo


      

     Linf_err_tot = dabs(Linf_err)
     !Linf_err_vel_tot = dabs(Linf_err_vel)

     L2_err_tot = dabs(L2_err)
     L2_err_vel_tot = dabs(L2_err_vel)
     H1_err_tot = dabs(H1_err)
     EN_err_tot = dabs(En_err)

     if (EN_err_tot .ge. 100.d0) call EXIT(EXIT_ENERGY_ERROR)

     if(mpi_id .eq. 0) then  
        file_mpi = 'EN.ERR'; unit_mpi = 40;
        open(unit_mpi,file=file_mpi, position='APPEND')
        write(unit_mpi,*)  sqrt(H1_err_tot), sqrt(EN_err_tot)!, sqrt(L2_err_tot),sqrt(L2_err_vel_tot)
                            
        close(unit_mpi)
     endif


     
     
     
     
     
     
     
     end subroutine COMPUTE_ENERGY_ERROR
