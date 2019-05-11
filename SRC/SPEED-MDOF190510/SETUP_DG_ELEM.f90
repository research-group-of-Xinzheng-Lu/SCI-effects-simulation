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

!> @brief Setup for DG data structure.
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] nm  number of materials
!> @param[in] sd spectral degree vector 
!> @param[in] tag_mat label for materials
!> @param[in] cs_nnz_loc length of cs_loc
!> @param[in] cs_loc local spectral connectivity
!> @param[in] nn_loc number of local nodes
!> @param[in] local_n_num local nodes numeration
!> @param[in] ne_loc number of local elements
!> @param[in] local_el_num local elements numeration
!> @param[in] xs x-coordinate of GLL nodes
!> @param[in] ys y-coordinate of GLL nodes
!> @param[in] zs z-coordinate of GLL nodes
!> @param[in] nel_dg_loc number of local DG elements
!> @param[in] nel_dg_glo number of global DG elements
!> @param[in] i4count  vector containing info for DG interface (1 in the i-th position if the i-th 
!!                node lies on a DG interface)
!> @param[in] mpi_id mpi processor identity
!> @param[in] mpi_comm mpi communicator
!> @param[in] alfa11 bilinear mapping constant
!> @param[in] alfa12 bilinear mapping constant
!> @param[in] alfa13 bilinear mapping constant
!> @param[in] alfa21 bilinear mapping constant
!> @param[in] alfa22 bilinear mapping constant
!> @param[in] alfa23 bilinear mapping constant
!> @param[in] alfa31 bilinear mapping constant
!> @param[in] alfa32 bilinear mapping constant
!> @param[in] alfa33 bilinear mapping constant
!> @param[in] beta11 bilinear mapping constant
!> @param[in] beta12 bilinear mapping constant
!> @param[in] beta13 bilinear mapping constant
!> @param[in] beta21 bilinear mapping constant
!> @param[in] beta22 bilinear mapping constant
!> @param[in] beta23 bilinear mapping constant
!> @param[in] beta31 bilinear mapping constant
!> @param[in] beta32 bilinear mapping constant
!> @param[in] beta33 bilinear mapping constant
!> @param[in] gamma1 bilinear mapping constant
!> @param[in] gamma2 bilinear mapping constant
!> @param[in] gamma3 bilinear mapping constant
!> @param[in] delta1 bilinear mapping constant
!> @param[in] delta2 bilinear mapping constant
!> @param[in] delta3 bilinear mapping constant
!> @param[in] tag_dg_el label for DG interfaces
!> @param[in] tag_dg_yn label for projecting nodes from a DG surface
!> @param[in] nload_dg number of DG interfaces
!> @param[in] con_bc connectivity matrix for boundary faces
!> @param[in] nface number of boundary faces
!> @param[in] mpi_file folder where *.mpi file are stored
!> @param[out] dg_els  data structure for dg interface elements 
!> @param[out] scratch_dg_els  temporary data structure for dg interface elements 

     subroutine SETUP_DG_ELEM(nm, sd, tag_mat, cs_nnz_loc, cs_loc, &
                           nn_loc, local_n_num, ne_loc, local_el_num, &
                           xs,ys,zs,&
                           nel_dg_loc, nel_dg_glo, &
                           i4count, mpi_id, mpi_comm, mpi_np,&
                           alfa11, alfa12, alfa13, &
                           alfa21, alfa22, alfa23, &
                           alfa31, alfa32, alfa33, &
                           beta11, beta12, beta13, &
                           beta21, beta22, beta23, &
                           beta31, beta32, beta33, &
                           gamma1, gamma2, gamma3, &
                           delta1, delta2, delta3, &
                           dg_els, scratch_dg_els, &
                           tag_dg_el, tag_dg_yn, nload_dg, &
                           con_bc, nface,mpi_file)
                     

     use max_var
     use str_mesh 
     use str_mesh_scratch                   


     implicit none
                   
     include 'SPEED.MPI'
     
     type(ELEMENT), dimension(nel_dg_loc), intent(inout) :: dg_els
     type(scratch_ELEMENT), dimension(nel_dg_loc), intent(inout) :: scratch_dg_els
  
     character*70 :: filempi, filename, cmd, mpi_file, filempi_new
     
     integer*4 :: nm, cs_nnz_loc, nn_loc, ne_loc, nel_dg_loc, nel_dg_glo
     integer*4 :: mpi_comm, mpi_id, mpi_ierr, mpi_np, nload_dg, tag_ind, nface
     integer*4 :: im, nn, ie, ned, unitmpi, unitname, nofel, nel_dg_proc
     integer*4 :: ne1, ne2, ne3, ne4, ic1, ic2, ic3, ic4
     integer*4 :: ne5, ne6, ne7, ne8, ic5, ic6, ic7, ic8
     integer*4 :: el_conf, face_conf, face_found, imate, iele, iface
     integer*4 :: ip, k, j, i, it, ih, ik, tt, indice, node_not_ass

     integer*4, dimension(nm) :: tag_mat, sd
     integer*4, dimension(0:cs_nnz_loc) :: cs_loc
     integer*4, dimension(nn_loc) :: local_n_num, i4count
     integer*4, dimension(ne_loc) :: local_el_num
     integer*4, dimension(nload_dg) :: tag_dg_el, tag_dg_yn

     integer*4, dimension(:,:), allocatable :: ielem_dg 
     integer*4, dimension(:,:), allocatable :: mat_el_face
     integer*4, dimension(nface,5) :: con_bc

     real*8 :: normal_x, normal_y, normal_z
     real*8 :: c_alfa11, c_alfa12, c_alfa13, c_alfa21, c_alfa22, c_alfa23, c_alfa31, c_alfa32, c_alfa33
     real*8 :: c_beta11, c_beta12, c_beta13, c_beta21, c_beta22, c_beta23, c_beta31, c_beta32, c_beta33
     real*8 :: c_gamma1, c_gamma2, c_gamma3, c_delta1, c_delta2, c_delta3
     real*8 :: xnod, ynod, znod, csi, eta, zeta, xnod1, ynod1, znod1
     real*8 :: val1, val2, val3, val4, val5, val6, valmin
     real*8 :: coef_a, coef_b, coef_t, det_trasf

     real*8, dimension(:), allocatable :: ctgl,wwgl
     real*8, dimension(:), allocatable :: ct, ww
     real*8, dimension(nn_loc) :: xs,ys,zs
     real*8, dimension(ne_loc) :: alfa11,alfa12,alfa13
     real*8, dimension(ne_loc) :: alfa21,alfa22,alfa23
     real*8, dimension(ne_loc) :: alfa31,alfa32,alfa33
     real*8, dimension(ne_loc) :: beta11,beta12,beta13
     real*8, dimension(ne_loc) :: beta21,beta22,beta23
     real*8, dimension(ne_loc) :: beta31,beta32,beta33
     real*8, dimension(ne_loc) :: gamma1,gamma2,gamma3
     real*8, dimension(ne_loc) :: delta1,delta2,delta3

     real*8, dimension(:,:), allocatable :: dd
     real*8, dimension(:,:), allocatable :: normalxyz

      nel_dg_loc = 0     
      ned = cs_loc(0) - 1


      do im = 1,nm

         nn = sd(im) +1         

         do ie = 1,ned
            if (cs_loc(cs_loc(ie -1) + 0) .eq. tag_mat(im)) then
               
               ! face x = - 1
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(1 -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(1 -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(1 -1) +1)
                       
                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then

                  call GET_TAG_BC(con_bc, nface, &
                                   local_n_num(ne1),local_n_num(ne2), &
                                   local_n_num(ne3),local_n_num(ne4), tag_dg_el, nload_dg, tag_ind)
                  
                  if (tag_ind == 0) print*, 'ERROR IN GET_TAG_BC'
                   
                  

                  nel_dg_loc = nel_dg_loc + 1                    
                  dg_els(nel_dg_loc)%ind_el = local_el_num(ie)
                  dg_els(nel_dg_loc)%face_el = 1
                  dg_els(nel_dg_loc)%mat = tag_mat(im)
                  dg_els(nel_dg_loc)%spct_deg = nn-1
                  dg_els(nel_dg_loc)%quad_rule = nofqp
 
                  dg_els(nel_dg_loc)%proj_yn = tag_dg_yn(tag_ind)
 
                 
                  call MAKE_NORMAL(1,xs(ne1), xs(ne2), xs(ne3), xs(ne4), ys(ne1), ys(ne2), ys(ne3), ys(ne4), &
                                    zs(ne1), zs(ne2), zs(ne3), zs(ne4), normal_x, normal_y, normal_z, 0, 0)

                 dg_els(nel_dg_loc)%nx = normal_x
                 dg_els(nel_dg_loc)%ny = normal_y
                 dg_els(nel_dg_loc)%nz = normal_z    
                              
                  allocate(ct(dg_els(nel_dg_loc)%quad_rule), ww(dg_els(nel_dg_loc)%quad_rule), & 
                           dd(dg_els(nel_dg_loc)%quad_rule, dg_els(nel_dg_loc)%quad_rule))
                           
                  call MAKE_LGL_NW(dg_els(nel_dg_loc)%quad_rule, ct, ww, dd)  
                  
                  allocate(ctgl(dg_els(nel_dg_loc)%quad_rule), wwgl(dg_els(nel_dg_loc)%quad_rule))
                  
                  call MAKE_GL_NW(dg_els(nel_dg_loc)%quad_rule, ctgl, wwgl)  
                  
                  ip = 0
                  do k = 1, dg_els(nel_dg_loc)%quad_rule
                    do j = 1, dg_els(nel_dg_loc)%quad_rule
                      do i = 1,1

                        ip = ip + 1
                        scratch_dg_els(nel_dg_loc)%x_nq(ip) = alfa11(ie)*ct(i) + alfa12(ie)*ctgl(j) &
                           + alfa13(ie)*ctgl(k) + beta11(ie)*ctgl(j)*ctgl(k) &
                           + beta12(ie)*ct(i)*ctgl(k) + beta13(ie)*ct(i)*ctgl(j) &
                           + gamma1(ie)*ct(i)*ctgl(j)*ctgl(k) + delta1(ie)
                        
                        scratch_dg_els(nel_dg_loc)%y_nq(ip) = alfa21(ie)*ct(i) + alfa22(ie)*ctgl(j) &
                           + alfa23(ie)*ctgl(k) + beta21(ie)*ctgl(j)*ctgl(k) &
                           + beta22(ie)*ct(i)*ctgl(k) + beta23(ie)*ct(i)*ctgl(j) &
                           + gamma2(ie)*ct(i)*ctgl(j)*ctgl(k) + delta2(ie)
                        
                        scratch_dg_els(nel_dg_loc)%z_nq(ip) = alfa31(ie)*ct(i) + alfa32(ie)*ctgl(j) &
                           + alfa33(ie)*ctgl(k) + beta31(ie)*ctgl(j)*ctgl(k) &
                           + beta32(ie)*ct(i)*ctgl(k) + beta33(ie)*ct(i)*ctgl(j) &
                           + gamma3(ie)*ct(i)*ctgl(j)*ctgl(k) + delta3(ie)
                        
                      
                        dg_els(nel_dg_loc)%wx_pl(ip) = 1.d0
                        dg_els(nel_dg_loc)%wy_pl(ip) = wwgl(j)
                        dg_els(nel_dg_loc)%wz_pl(ip) = wwgl(k)

                        enddo
                     enddo
                   enddo
                   deallocate(ct,ww,dd)
                   deallocate(ctgl,wwgl)

               endif               
               
               ! face y = - 1
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(1 -1) +1)

                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then
 

                  call GET_TAG_BC(con_bc, nface, &
                                   local_n_num(ne1),local_n_num(ne2), &
                                   local_n_num(ne3),local_n_num(ne4), tag_dg_el, nload_dg, tag_ind)
                  
                  if (tag_ind == 0) print*, 'ERROR IN GET_TAG_BC'


                  nel_dg_loc = nel_dg_loc +1
                  dg_els(nel_dg_loc)%ind_el = local_el_num(ie)
                  dg_els(nel_dg_loc)%face_el = 2
                  dg_els(nel_dg_loc)%mat = tag_mat(im)
                  dg_els(nel_dg_loc)%spct_deg = nn-1
                  dg_els(nel_dg_loc)%quad_rule = nofqp
 
                  dg_els(nel_dg_loc)%proj_yn = tag_dg_yn(tag_ind)

                 
                  call MAKE_NORMAL(2,xs(ne1), xs(ne2), xs(ne3), xs(ne4), ys(ne1), ys(ne2), ys(ne3), ys(ne4), &
                                    zs(ne1), zs(ne2), zs(ne3), zs(ne4), normal_x, normal_y, normal_z, 0, 0)

                 dg_els(nel_dg_loc)%nx = normal_x
                 dg_els(nel_dg_loc)%ny = normal_y
                 dg_els(nel_dg_loc)%nz = normal_z                 
                 
                  allocate(ct(dg_els(nel_dg_loc)%quad_rule), ww(dg_els(nel_dg_loc)%quad_rule), &
                           dd(dg_els(nel_dg_loc)%quad_rule,dg_els(nel_dg_loc)%quad_rule))
                           
                  call MAKE_LGL_NW(dg_els(nel_dg_loc)%quad_rule,ct,ww,dd)  
                  
                  allocate(ctgl(dg_els(nel_dg_loc)%quad_rule),wwgl(dg_els(nel_dg_loc)%quad_rule))
                  
                  call MAKE_GL_NW(dg_els(nel_dg_loc)%quad_rule,ctgl,wwgl)  
                  
                  ip = 0
                  do k = 1, dg_els(nel_dg_loc)%quad_rule
                    do j = 1,1
                      do i = 1, dg_els(nel_dg_loc)%quad_rule 
                        ip = ip + 1
                        scratch_dg_els(nel_dg_loc)%x_nq(ip) = alfa11(ie)*ctgl(i) + alfa12(ie)*ct(j) &
                           + alfa13(ie)*ctgl(k) + beta11(ie)*ct(j)*ctgl(k) &
                           + beta12(ie)*ctgl(i)*ctgl(k) + beta13(ie)*ctgl(i)*ct(j) &
                           + gamma1(ie)*ctgl(i)*ct(j)*ctgl(k) + delta1(ie)
                        
                        scratch_dg_els(nel_dg_loc)%y_nq(ip) = alfa21(ie)*ctgl(i) + alfa22(ie)*ct(j) &
                           + alfa23(ie)*ctgl(k) + beta21(ie)*ct(j)*ctgl(k) &
                           + beta22(ie)*ctgl(i)*ctgl(k) + beta23(ie)*ctgl(i)*ct(j) &
                           + gamma2(ie)*ctgl(i)*ct(j)*ctgl(k) + delta2(ie)
                        
                        scratch_dg_els(nel_dg_loc)%z_nq(ip) = alfa31(ie)*ctgl(i) + alfa32(ie)*ct(j) &
                           + alfa33(ie)*ctgl(k) + beta31(ie)*ct(j)*ctgl(k) &
                           + beta32(ie)*ctgl(i)*ctgl(k) + beta33(ie)*ctgl(i)*ct(j) &
                           + gamma3(ie)*ctgl(i)*ct(j)*ctgl(k) + delta3(ie)

                               dg_els(nel_dg_loc)%wx_pl(ip) = wwgl(i)
                        dg_els(nel_dg_loc)%wy_pl(ip) = 1.d0
                        dg_els(nel_dg_loc)%wz_pl(ip) = wwgl(k)                       

                       enddo
                     enddo
                   enddo
                   deallocate(ct,ww,dd)
                   deallocate(ctgl,wwgl)
               endif               

               ! face z = - 1
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(1 -1) +1)

                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then

                  call GET_TAG_BC(con_bc, nface, &
                                   local_n_num(ne1),local_n_num(ne2), &
                                   local_n_num(ne3),local_n_num(ne4), tag_dg_el, nload_dg, tag_ind)
                  
                  if (tag_ind == 0) print*, 'ERROR IN GET_TAG_BC'


                  nel_dg_loc = nel_dg_loc +1
                  dg_els(nel_dg_loc)%ind_el = local_el_num(ie)
                  dg_els(nel_dg_loc)%face_el = 3
                  dg_els(nel_dg_loc)%mat = tag_mat(im)
                  dg_els(nel_dg_loc)%spct_deg = nn-1
                  dg_els(nel_dg_loc)%quad_rule = nofqp

                  dg_els(nel_dg_loc)%proj_yn = tag_dg_yn(tag_ind)
                  
                  call MAKE_NORMAL(3,xs(ne1), xs(ne2), xs(ne3), xs(ne4), ys(ne1), ys(ne2), ys(ne3), ys(ne4), &
                                    zs(ne1), zs(ne2), zs(ne3), zs(ne4), normal_x, normal_y, normal_z, 0, 0)

                 dg_els(nel_dg_loc)%nx = normal_x
                 dg_els(nel_dg_loc)%ny = normal_y
                 dg_els(nel_dg_loc)%nz = normal_z                 
                 
                  allocate(ct(dg_els(nel_dg_loc)%quad_rule), ww(dg_els(nel_dg_loc)%quad_rule), & 
                            dd(dg_els(nel_dg_loc)%quad_rule,dg_els(nel_dg_loc)%quad_rule))
                            
                  call MAKE_LGL_NW(dg_els(nel_dg_loc)%quad_rule, ct, ww, dd) 
                  
                  allocate(ctgl(dg_els(nel_dg_loc)%quad_rule), wwgl(dg_els(nel_dg_loc)%quad_rule))
                  
                  call MAKE_GL_NW(dg_els(nel_dg_loc)%quad_rule, ctgl, wwgl)   
                  
                  ip = 0
                  do k = 1,1
                    do j = 1, dg_els(nel_dg_loc)%quad_rule
                      do i = 1, dg_els(nel_dg_loc)%quad_rule
                        ip = ip + 1
                        scratch_dg_els(nel_dg_loc)%x_nq(ip) = alfa11(ie)*ctgl(i) + alfa12(ie)*ctgl(j) &
                           + alfa13(ie)*ct(k) + beta11(ie)*ctgl(j)*ct(k) &
                           + beta12(ie)*ctgl(i)*ct(k) + beta13(ie)*ctgl(i)*ctgl(j) &
                           + gamma1(ie)*ctgl(i)*ctgl(j)*ct(k) + delta1(ie)
                                                   
                        scratch_dg_els(nel_dg_loc)%y_nq(ip) = alfa21(ie)*ctgl(i) + alfa22(ie)*ctgl(j) &
                           + alfa23(ie)*ct(k) + beta21(ie)*ctgl(j)*ct(k) &
                           + beta22(ie)*ctgl(i)*ct(k) + beta23(ie)*ctgl(i)*ctgl(j) &
                           + gamma2(ie)*ctgl(i)*ctgl(j)*ct(k) + delta2(ie)
                        
                        scratch_dg_els(nel_dg_loc)%z_nq(ip) = alfa31(ie)*ctgl(i) + alfa32(ie)*ctgl(j) &
                           + alfa33(ie)*ct(k) + beta31(ie)*ctgl(j)*ct(k) &
                           + beta32(ie)*ctgl(i)*ct(k) + beta33(ie)*ctgl(i)*ctgl(j) &
                           + gamma3(ie)*ctgl(i)*ctgl(j)*ct(k) + delta3(ie)
                            

                        dg_els(nel_dg_loc)%wx_pl(ip) = wwgl(i)
                        dg_els(nel_dg_loc)%wy_pl(ip) = wwgl(j)
                        dg_els(nel_dg_loc)%wz_pl(ip) = 1.d0 
                                               
                       enddo
                     enddo
                   enddo
                   deallocate(ct,ww,dd)
                   deallocate(ctgl,wwgl)
             endif               


               ! face x = 1
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(nn -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(nn -1) +1)

                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then

                  call GET_TAG_BC(con_bc, nface, &
                                   local_n_num(ne1),local_n_num(ne2), &
                                   local_n_num(ne3),local_n_num(ne4), tag_dg_el, nload_dg, tag_ind)
                  
                  if (tag_ind == 0) print*, 'ERROR IN GET_TAG_BC'

                  nel_dg_loc = nel_dg_loc +1
                  dg_els(nel_dg_loc)%ind_el = local_el_num(ie)
                  dg_els(nel_dg_loc)%face_el = 4
                  dg_els(nel_dg_loc)%mat = tag_mat(im)
                  dg_els(nel_dg_loc)%spct_deg = nn-1
                  dg_els(nel_dg_loc)%quad_rule = nofqp

                  dg_els(nel_dg_loc)%proj_yn = tag_dg_yn(tag_ind)
 
                  call MAKE_NORMAL(4,xs(ne1), xs(ne2), xs(ne3), xs(ne4), ys(ne1), ys(ne2), ys(ne3), ys(ne4), &
                                    zs(ne1), zs(ne2), zs(ne3), zs(ne4), normal_x, normal_y, normal_z, 0, 0)

                 dg_els(nel_dg_loc)%nx = normal_x
                 dg_els(nel_dg_loc)%ny = normal_y
                 dg_els(nel_dg_loc)%nz = normal_z                 
                 
                  allocate(ct(dg_els(nel_dg_loc)%quad_rule),ww(dg_els(nel_dg_loc)%quad_rule), & 
                           dd(dg_els(nel_dg_loc)%quad_rule,dg_els(nel_dg_loc)%quad_rule))
                           
                  call MAKE_LGL_NW(dg_els(nel_dg_loc)%quad_rule, ct, ww, dd) 
                  
                  allocate(ctgl(dg_els(nel_dg_loc)%quad_rule),wwgl(dg_els(nel_dg_loc)%quad_rule))
                  
                  call MAKE_GL_NW(dg_els(nel_dg_loc)%quad_rule, ctgl, wwgl)   
                  
                  ip = 0
                  do k = 1, dg_els(nel_dg_loc)%quad_rule
                    do j = 1, dg_els(nel_dg_loc)%quad_rule
                      do i = dg_els(nel_dg_loc)%quad_rule, dg_els(nel_dg_loc)%quad_rule
                        ip = ip + 1
                        scratch_dg_els(nel_dg_loc)%x_nq(ip) = alfa11(ie)*ct(i) + alfa12(ie)*ctgl(j) &
                           + alfa13(ie)*ctgl(k) + beta11(ie)*ctgl(j)*ctgl(k) &
                           + beta12(ie)*ct(i)*ctgl(k) + beta13(ie)*ct(i)*ctgl(j) &
                           + gamma1(ie)*ct(i)*ctgl(j)*ctgl(k) + delta1(ie)
                        
                        scratch_dg_els(nel_dg_loc)%y_nq(ip) = alfa21(ie)*ct(i) + alfa22(ie)*ctgl(j) &
                           + alfa23(ie)*ctgl(k) + beta21(ie)*ctgl(j)*ctgl(k) &
                           + beta22(ie)*ct(i)*ctgl(k) + beta23(ie)*ct(i)*ctgl(j) &
                           + gamma2(ie)*ct(i)*ctgl(j)*ctgl(k) + delta2(ie)
                        
                        scratch_dg_els(nel_dg_loc)%z_nq(ip) = alfa31(ie)*ct(i) + alfa32(ie)*ctgl(j) &
                           + alfa33(ie)*ctgl(k) + beta31(ie)*ctgl(j)*ctgl(k) &
                           + beta32(ie)*ct(i)*ctgl(k) + beta33(ie)*ct(i)*ctgl(j) &
                           + gamma3(ie)*ct(i)*ctgl(j)*ctgl(k) + delta3(ie)
                     
                               dg_els(nel_dg_loc)%wx_pl(ip) = 1.d0
                        dg_els(nel_dg_loc)%wy_pl(ip) = wwgl(j)
                        dg_els(nel_dg_loc)%wz_pl(ip) = wwgl(k)                  

                       enddo
                     enddo
                   enddo
                   deallocate(ct,ww,dd)
                   deallocate(ctgl,wwgl)
               endif             


               ! face y = 1
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(1 -1) +1)

                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then
                 
                  call GET_TAG_BC(con_bc, nface, &
                                   local_n_num(ne1),local_n_num(ne2), &
                                   local_n_num(ne3),local_n_num(ne4), tag_dg_el, nload_dg, tag_ind)
                  
                  if (tag_ind == 0) print*, 'ERROR IN GET_TAG_BC'

                  nel_dg_loc = nel_dg_loc +1
                  dg_els(nel_dg_loc)%ind_el = local_el_num(ie)
                  dg_els(nel_dg_loc)%face_el = 5
                  dg_els(nel_dg_loc)%mat = tag_mat(im)
                  dg_els(nel_dg_loc)%spct_deg = nn-1
                  dg_els(nel_dg_loc)%quad_rule = nofqp
                  
                  dg_els(nel_dg_loc)%proj_yn = tag_dg_yn(tag_ind)

                  call MAKE_NORMAL(5,xs(ne1), xs(ne2), xs(ne3), xs(ne4), ys(ne1), ys(ne2), ys(ne3), ys(ne4), &
                                    zs(ne1), zs(ne2), zs(ne3), zs(ne4), normal_x, normal_y, normal_z, 0, 0)

                 dg_els(nel_dg_loc)%nx = normal_x
                 dg_els(nel_dg_loc)%ny = normal_y
                 dg_els(nel_dg_loc)%nz = normal_z                 
                 
               
                  allocate(ct(dg_els(nel_dg_loc)%quad_rule), ww(dg_els(nel_dg_loc)%quad_rule), & 
                           dd(dg_els(nel_dg_loc)%quad_rule,dg_els(nel_dg_loc)%quad_rule))
                           
                  call MAKE_LGL_NW(dg_els(nel_dg_loc)%quad_rule, ct, ww, dd)  
                  
                  allocate(ctgl(dg_els(nel_dg_loc)%quad_rule),wwgl(dg_els(nel_dg_loc)%quad_rule))
                  
                  call MAKE_GL_NW(dg_els(nel_dg_loc)%quad_rule, ctgl, wwgl)  
                  ip = 0
                  do k = 1, dg_els(nel_dg_loc)%quad_rule
                    do j = dg_els(nel_dg_loc)%quad_rule, dg_els(nel_dg_loc)%quad_rule
                      do i = 1, dg_els(nel_dg_loc)%quad_rule
                        ip = ip + 1
                        scratch_dg_els(nel_dg_loc)%x_nq(ip) = alfa11(ie)*ctgl(i) + alfa12(ie)*ct(j) &
                           + alfa13(ie)*ctgl(k) + beta11(ie)*ct(j)*ctgl(k) &
                           + beta12(ie)*ctgl(i)*ctgl(k) + beta13(ie)*ctgl(i)*ct(j) &
                           + gamma1(ie)*ctgl(i)*ct(j)*ctgl(k) + delta1(ie)
                        
                        scratch_dg_els(nel_dg_loc)%y_nq(ip) = alfa21(ie)*ctgl(i) + alfa22(ie)*ct(j) &
                           + alfa23(ie)*ctgl(k) + beta21(ie)*ct(j)*ctgl(k) &
                           + beta22(ie)*ctgl(i)*ctgl(k) + beta23(ie)*ctgl(i)*ct(j) &
                           + gamma2(ie)*ctgl(i)*ct(j)*ctgl(k) + delta2(ie)
                        
                        scratch_dg_els(nel_dg_loc)%z_nq(ip) = alfa31(ie)*ctgl(i) + alfa32(ie)*ct(j) &
                           + alfa33(ie)*ctgl(k) + beta31(ie)*ct(j)*ctgl(k) &
                           + beta32(ie)*ctgl(i)*ctgl(k) + beta33(ie)*ctgl(i)*ct(j) &
                           + gamma3(ie)*ctgl(i)*ct(j)*ctgl(k) + delta3(ie)

                        dg_els(nel_dg_loc)%wx_pl(ip) = wwgl(i)
                        dg_els(nel_dg_loc)%wy_pl(ip) = 1.d0
                        dg_els(nel_dg_loc)%wz_pl(ip) = wwgl(k)                                                

                       enddo
                     enddo
                   enddo
                   deallocate(ct,ww,dd)
                   deallocate(ctgl,wwgl)
               endif               


               ! face z = 1
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(1 -1) +1)

                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then
                 
                  call GET_TAG_BC(con_bc, nface, &
                                   local_n_num(ne1),local_n_num(ne2), &
                                   local_n_num(ne3),local_n_num(ne4), tag_dg_el, nload_dg, tag_ind)
                  
                  if (tag_ind == 0) print*, 'ERROR IN GET_TAG_BC'

                  nel_dg_loc = nel_dg_loc + 1                  
                  dg_els(nel_dg_loc)%ind_el = local_el_num(ie)
                  dg_els(nel_dg_loc)%face_el = 6
                  dg_els(nel_dg_loc)%mat = tag_mat(im)
                  dg_els(nel_dg_loc)%spct_deg = nn-1
                  dg_els(nel_dg_loc)%quad_rule = nofqp
                  
                  dg_els(nel_dg_loc)%proj_yn = tag_dg_yn(tag_ind)

                  call MAKE_NORMAL(6,xs(ne1), xs(ne2), xs(ne3), xs(ne4), ys(ne1), ys(ne2), ys(ne3), ys(ne4), &
                                    zs(ne1), zs(ne2), zs(ne3), zs(ne4), normal_x, normal_y, normal_z, 0, 0)

                 dg_els(nel_dg_loc)%nx = normal_x
                 dg_els(nel_dg_loc)%ny = normal_y
                 dg_els(nel_dg_loc)%nz = normal_z                 
                 
                   allocate(ct(dg_els(nel_dg_loc)%quad_rule),ww(dg_els(nel_dg_loc)%quad_rule), & 
                            dd(dg_els(nel_dg_loc)%quad_rule,dg_els(nel_dg_loc)%quad_rule))
                            
                  call MAKE_LGL_NW(dg_els(nel_dg_loc)%quad_rule, ct, ww, dd)
                  
                  allocate(ctgl(dg_els(nel_dg_loc)%quad_rule), wwgl(dg_els(nel_dg_loc)%quad_rule))
                  
                  call MAKE_GL_NW(dg_els(nel_dg_loc)%quad_rule, ctgl, wwgl)    
                  
                  ip = 0
                  do k = dg_els(nel_dg_loc)%quad_rule,dg_els(nel_dg_loc)%quad_rule
                    do j = 1,dg_els(nel_dg_loc)%quad_rule
                      do i = 1,dg_els(nel_dg_loc)%quad_rule
                        ip = ip + 1
                        scratch_dg_els(nel_dg_loc)%x_nq(ip) = alfa11(ie)*ctgl(i) + alfa12(ie)*ctgl(j) &
                           + alfa13(ie)*ct(k) + beta11(ie)*ctgl(j)*ct(k) &
                           + beta12(ie)*ctgl(i)*ct(k) + beta13(ie)*ctgl(i)*ctgl(j) &
                           + gamma1(ie)*ctgl(i)*ctgl(j)*ct(k) + delta1(ie)
                        
                        scratch_dg_els(nel_dg_loc)%y_nq(ip) = alfa21(ie)*ctgl(i) + alfa22(ie)*ctgl(j) &
                           + alfa23(ie)*ct(k) + beta21(ie)*ctgl(j)*ct(k) &
                           + beta22(ie)*ctgl(i)*ct(k) + beta23(ie)*ctgl(i)*ctgl(j) &
                           + gamma2(ie)*ctgl(i)*ctgl(j)*ct(k) + delta2(ie)
                        
                        scratch_dg_els(nel_dg_loc)%z_nq(ip) = alfa31(ie)*ctgl(i) + alfa32(ie)*ctgl(j) &
                           + alfa33(ie)*ct(k) + beta31(ie)*ctgl(j)*ct(k) &
                           + beta32(ie)*ctgl(i)*ct(k) + beta33(ie)*ctgl(i)*ctgl(j) &
                           + gamma3(ie)*ctgl(i)*ctgl(j)*ct(k) + delta3(ie)
         
                        dg_els(nel_dg_loc)%wx_pl(ip) = wwgl(i)
                        dg_els(nel_dg_loc)%wy_pl(ip) = wwgl(j)
                        dg_els(nel_dg_loc)%wz_pl(ip) = 1.d0

                       enddo
                     enddo
                   enddo
                   deallocate(ct,ww,dd)
                   deallocate(ctgl,wwgl)
               endif               
            endif
         enddo
      enddo

     
     !write in order
     !- material - element - face - normalx, normaly, normalz
     
     filempi = 'NORM000000.mpi'
     unitmpi = 40                                 

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
       write(unitmpi,*) nel_dg_loc 
       do i = 1, nel_dg_loc
          write(unitmpi,"(1I2,1X,1I12,1X,1I2,1X,3(1X,ES12.4))") &
                       dg_els(i)%mat, dg_els(i)%ind_el, dg_els(i)%face_el, &
                       dg_els(i)%nx, dg_els(i)%ny, dg_els(i)%nz 
       enddo
       close(unitmpi)     
     
       call MPI_BARRIER(mpi_comm, mpi_ierr)
     
       
       
       if(mpi_id .eq. 0) then 
         
         filename = 'NORMALL.input'
         unitname = 500
         open(unitname,file=filename)
         write(unitname,*) nel_dg_glo
         !close(unitname)
         
         allocate(mat_el_face(nel_dg_glo,3), normalxyz(nel_dg_glo,3))
         k = 1
         
         do i = 1, mpi_np
                           
              filempi = 'NORM000000.mpi'

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
             read(unitmpi,*) nel_dg_proc
                                           
             do j = 1, nel_dg_proc
                 read(unitmpi,*) mat_el_face(k,1),mat_el_face(k,2), mat_el_face(k,3), &
                                 normalxyz(k,1), normalxyz(k,2), normalxyz(k,3)
                    k=k+1
             enddo
          
             close(unitmpi)
         
          enddo
          
          do j = 1, nel_dg_glo
             write(unitname,*) mat_el_face(j,1),mat_el_face(j,2), mat_el_face(j,3), &
                            normalxyz(j,1), normalxyz(j,2), normalxyz(j,3)
          enddo
          
          deallocate(mat_el_face, normalxyz)
          close(unitname)
         
         
       endif

    
       call MPI_BARRIER(mpi_comm, mpi_ierr) 
       
       
     
      
      
      return
      
      end subroutine SETUP_DG_ELEM
