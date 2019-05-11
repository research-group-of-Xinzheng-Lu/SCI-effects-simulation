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

!> @brief Creates data structure for dg elements.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nm number of materials
!> @param[in] sd polynomial degree vector
!> @param[in] tag_mat  materials label
!> @param[in] prop_mat material properties rho, lambda, mu, gamma
!> @param[in] cs_nnz_loc length of cs_loc
!> @param[in] cs_loc local spectral connectivity vector
!> @param[in] nn_loc number of local nodes
!> @param[in] local_n_num local node numbering
!> @param[in] ne_loc number of local elements
!> @param[in] local_el_num local elment numbering
!> @param[in] xs x-coordinate of GLL nodes
!> @param[in] ys y-coordinate of GLL nodes
!> @param[in] zs z-coordinate of GLL nodes  
!> @param[in] nel_dg_loc  number of local dg element
!> @param[in] nel_dg_glo  number of global dg element
!> @param[in] i4count  i4count(i) = 1 if the i-th node belongs to a
!!                         dg surface
!> @param[in] mpi_id  mpi process id
!> @param[in] mpi_comm  mpi common world
!> @param[in] mpi_np  number of mpi process 
!> @param[in] alfa11 costant values for the bilinear map
!> @param[in] alfa12 costant values for the bilinear map
!> @param[in] alfa13 costant values for the bilinear map
!> @param[in] alfa21 costant values for the bilinear map
!> @param[in] alfa22 costant values for the bilinear map
!> @param[in] alfa23 costant values for the bilinear map
!> @param[in] alfa31 costant values for the bilinear map
!> @param[in] alfa32 costant values for the bilinear map
!> @param[in] alfa33 costant values for the bilinear map
!> @param[in] beta11 costant values for the bilinear map 
!> @param[in] beta12 costant values for the bilinear map
!> @param[in] beta13 costant values for the bilinear map
!> @param[in] beta21 costant values for the bilinear map
!> @param[in] beta22 costant values for the bilinear map
!> @param[in] beta23 costant values for the bilinear map
!> @param[in] beta31 costant values for the bilinear map
!> @param[in] beta32 costant values for the bilinear map
!> @param[in] beta33 costant values for the bilinear map
!> @param[in] gamma1 costant values for the bilinear map
!> @param[in] gamma2 costant values for the bilinear map
!> @param[in] gamma3 costant values for the bilinear map
!> @param[in] delta1 costant values for the bilinear map
!> @param[in] delta2 costant values for the bilinear map
!> @param[in] delta3 costant values for the bilinear map
!> @param[in] dg_cnst  -1 = SIPG, 0 = IIPG, 1 = NIPG
!> @param[in] penalty_c  penalty constant
!> @param[in] faces  identification of dg surface: faces(1,i) = material,
!!                                             faces(2,i) = element,
!!                                             faces(3,i) = face.
!> @param[in] area_nodes  area_nodes(1,i) area of the face i, 
!!                    area_nodes(2,i),...,area_nodes(25,i) constant for the bilinear map
!> @param[in] filename file name (DGFS.input) where DG info are written 
!> @param[in] testmode 1 if testmode is active 
!> @param[out] el_new  data structure for dg elements


   subroutine MAKE_DG_INTERFACE(nm, sd, tag_mat, prop_mat, cs_nnz_loc, cs_loc, &
                           nn_loc, local_n_num, ne_loc, local_el_num, &
                           xs, ys, zs,&
                           nel_dg_loc, nel_dg_glo, &
                           i4count, mpi_id, mpi_comm, mpi_np,&
                           alfa11, alfa12, alfa13, &
                           alfa21, alfa22, alfa23, &
                           alfa31, alfa32, alfa33, &
                           beta11, beta12, beta13, &
                           beta21, beta22, beta23, &
                           beta31, beta32, beta33, &
                           gamma1, gamma2, gamma3, &
                           delta1, delta2, delta3, dg_cnst, penalty_c, &
                           faces, area_nodes, el_new, filename, testmode)


     use max_var
     use str_mesh_after 
     use DGJUMP

     implicit none
  
     type(ELEMENT_after), dimension(:), allocatable :: dg_els
     type(el4loop), dimension(nel_dg_loc), intent(inout):: el_new                        
 
     character*70 :: filename, filename1, filename2

     integer*4 :: mpi_comm, mpi_id, mpi_ierr, mpi_np, ishift, jshift
     integer*4 :: im, nn, ie, ned, ip, mm, nnz, nnz_p, nnz_m, nnz_p_only_uv, nnz_m_only_uv
     integer*4 :: ne1, ne2, ne3, ne4, ic1, ic2, ic3, ic4
     integer*4 :: ne5, ne6, ne7, ne8, ic5, ic6, ic7, ic8
     integer*4 :: nm, cs_nnz_loc, nn_loc, ne_loc, nel_dg_loc, nel_dg_glo
     integer*4 :: n_line, ielem, iface, iene, ifacene, imne
     integer*4 :: int_trash, statuss, i, tt, ih, it, p, j,k
     integer*4 :: unit_interface, unitname, unitname1, unitname2
     integer*4 :: error
     integer*4 :: mpierror
     integer*4 :: nofne_el, ic, testmode

     integer*4, dimension(nm) :: tag_mat, sd
     integer*4, dimension(0:cs_nnz_loc) :: cs_loc
     integer*4, dimension(nn_loc) :: local_n_num, i4count
     integer*4, dimension(ne_loc) :: local_el_num
     integer*4, dimension(:), allocatable :: I4S, J4S     

     integer*4, dimension(:,:), allocatable :: con_DG
     integer*4, dimension(:,:), allocatable :: copia
     integer*4, dimension(3,nel_dg_glo) :: faces

     real*8 :: real_trash, lambda, mu, pen_p, pen_h, pen, penalty_c, dg_cnst
     real*8 :: cp_a, cp_b, cp_c, cp_d, cp_e, cp_f, cp_g, cp_n, cp_p
     real*8 :: csi, eta, zeta, normal_x, normal_y, normal_z
     real*8 :: c_alfa11, c_alfa12, c_alfa13, c_alfa21, c_alfa22, c_alfa23, c_alfa31, c_alfa32, c_alfa33
     real*8 :: c_beta11, c_beta12, c_beta13, c_beta21, c_beta22, c_beta23, c_beta31, c_beta32, c_beta33
     real*8 :: c_gamma1, c_gamma2, c_gamma3, c_delta1, c_delta2, c_delta3

     real*8, dimension(ne_loc) :: alfa11,alfa12,alfa13
     real*8, dimension(ne_loc) :: alfa21,alfa22,alfa23
     real*8, dimension(ne_loc) :: alfa31,alfa32,alfa33
     real*8, dimension(ne_loc) :: beta11,beta12,beta13
     real*8, dimension(ne_loc) :: beta21,beta22,beta23
     real*8, dimension(ne_loc) :: beta31,beta32,beta33
     real*8, dimension(ne_loc) :: gamma1,gamma2,gamma3
     real*8, dimension(ne_loc) :: delta1,delta2,delta3
     
     real*8, dimension(:), allocatable :: M4S

     real*8, dimension(:,:), allocatable :: nodes_DG
     real*8, dimension(:,:), allocatable :: JP,JM, JP_only_uv,JM_only_uv
     real*8, dimension(nn_loc) :: xs,ys,zs
     real*8, dimension(nm,4) :: prop_mat
     real*8, dimension(25,nel_dg_glo) :: area_nodes

     
     
!    Reading DGFS.input file'
     unitname = 50
     open(unitname,file=filename)
     read(unitname,*) n_line

     allocate(con_DG(6,2*n_line), nodes_DG(6,2*n_line))
     con_DG = 0
     nodes_DG = 0.d0
     
     
     do i = 1, n_line
        read(unitname,*) con_DG(1,i), con_DG(2,i), con_DG(3,i), con_DG(4,i), con_DG(5,i), con_DG(6,i), & 
                          nodes_DG(1,i), nodes_DG(2,i), nodes_DG(3,i), nodes_DG(4,i), nodes_DG(5,i), nodes_DG(6,i)
     enddo


     close(unitname)


     do i = 1, 3
        con_DG(i, n_line+1: 2*n_line) = con_DG(3+i,1:n_line)
     enddo  

     do i = 1, 3
        con_DG(3+i, n_line+1: 2*n_line) = con_DG(i,1:n_line)
     enddo  



     do i = 1, 6
        nodes_DG(i, n_line+1: 2*n_line) = nodes_DG(i,1:n_line)       
     enddo

     allocate(dg_els(nel_dg_loc))

!***************************************************************************************************************

     nel_dg_loc = 0      
     ned = cs_loc(0) - 1


      do im = 1,nm

         nn = sd(im) +1         

         do ie = 1,ned
            if (cs_loc(cs_loc(ie -1) + 0) .eq. tag_mat(im)) then
       
               ! face 1: x = - 1
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(1 -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(1 -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(1 -1) +1)
                                         
                           
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then

                  ip = 0  
                  nel_dg_loc = nel_dg_loc + 1                        
                  dg_els(nel_dg_loc)%ind_el = local_el_num(ie)
                  dg_els(nel_dg_loc)%face_el = 1
                  dg_els(nel_dg_loc)%mat = tag_mat(im)
                  dg_els(nel_dg_loc)%spct_deg = nn-1
                  dg_els(nel_dg_loc)%quad_rule = ip

                  i = 1
                  do  while (i .le. 2*n_line) 
                  
                    if(con_DG(2,i) .eq. local_el_num(ie) .and. con_DG(3,i) .eq. 1) then 
                    
                    
                      ip = ip + 1
                      dg_els(nel_dg_loc)%quad_rule = ip
                      
                      call NEWTON_RAPSON(nodes_DG(1,i), nodes_DG(2,i), nodes_DG(3,i), &
                                        alfa11(ie),alfa12(ie),alfa13(ie), &
                                        alfa21(ie),alfa22(ie),alfa23(ie), &
                                        alfa31(ie),alfa32(ie),alfa33(ie), &
                                        beta11(ie),beta12(ie),beta13(ie), &
                                        beta21(ie),beta22(ie),beta23(ie), &
                                        beta31(ie),beta32(ie),beta33(ie), & 
                                        gamma1(ie),gamma2(ie),gamma3(ie), &
                                        delta1(ie),delta2(ie),delta3(ie), & 
                                        tt, csi, eta, zeta, nofinr, mpi_id,&
                                        con_DG(2,i),con_DG(5,i), 1.d-6, 1.01d0,1)
                                        
                                        
                        dg_els(nel_dg_loc)%x_pl(ip) = -1.d0
                        dg_els(nel_dg_loc)%y_pl(ip) = eta
                        dg_els(nel_dg_loc)%z_pl(ip) = zeta
                        dg_els(nel_dg_loc)%wx_pl(ip) = nodes_DG(4,i)
                        dg_els(nel_dg_loc)%wy_pl(ip) = nodes_DG(5,i)
                        dg_els(nel_dg_loc)%wz_pl(ip) = nodes_DG(6,i) 
                        
                        dg_els(nel_dg_loc)%omega_minus(ip,0) = con_DG(4,i)
                        dg_els(nel_dg_loc)%omega_minus(ip,1) = con_DG(5,i)
                        dg_els(nel_dg_loc)%omega_minus(ip,2) = con_DG(6,i)
                        dg_els(nel_dg_loc)%omega_minus(ip,3) = 0 
                         
                         
                         ! find the corresponding neigh. elem in faces matrix
                         
                         call GET_FACE_DG(faces, nel_dg_glo, con_DG(5,i), con_DG(6,i), ih)
                         
                         call MAKE_BILINEAR_MAP(area_nodes(2:25,ih), &
                                             c_alfa11, c_alfa12, c_alfa13, &
                                             c_alfa21, c_alfa22, c_alfa23, &
                                             c_alfa31, c_alfa32, c_alfa33, &  
                                             c_beta11, c_beta12, c_beta13, & 
                                             c_beta21, c_beta22, c_beta23, & 
                                             c_beta31, c_beta32, c_beta33, &
                                             c_gamma1, c_gamma2, c_gamma3, &
                                             c_delta1, c_delta2, c_delta3)
                         
                         
                         call NEWTON_RAPSON(nodes_DG(1,i), nodes_DG(2,i), nodes_DG(3,i), &
                                        c_alfa11, c_alfa12, c_alfa13, &
                                        c_alfa21, c_alfa22, c_alfa23,  &
                                        c_alfa31, c_alfa32, c_alfa33,  &
                                        c_beta11, c_beta12, c_beta13,  &
                                        c_beta21, c_beta22, c_beta23,  &
                                        c_beta31, c_beta32, c_beta33,  & 
                                        c_gamma1, c_gamma2, c_gamma3,  &
                                        c_delta1, c_delta2, c_delta3,  & 
                                        tt, csi, eta, zeta, nofinr, mpi_id,&
                                        con_DG(2,i),con_DG(5,i), 1.d-6, 1.01d0,1)   
                                        
                         select case(con_DG(6,i))    
                            case(1)
                              csi = -1.d0
                            case(2)
                              eta = -1.d0
                            case(3)
                              zeta = -1.d0
                            case(4)
                              csi = 1.d0
                            case(5)
                              eta = 1.d0
                            case(6)
                              zeta = 1.d0
                            case default
                              write(*,*) 'error in make_interface'
                            end select             

                         dg_els(nel_dg_loc)%x_mn(ip) =  csi
                         dg_els(nel_dg_loc)%y_mn(ip) =  eta
                         dg_els(nel_dg_loc)%z_mn(ip) =  zeta    


                        endif
                        i = i + 1
                         
                  enddo
                  

                  
                 
                   call MAKE_NORMAL(1,xs(ne1), xs(ne2), xs(ne3), xs(ne4), ys(ne1), ys(ne2), ys(ne3), ys(ne4), &
                                    zs(ne1), zs(ne2), zs(ne3), zs(ne4), normal_x, normal_y, normal_z, 0, 0)

                   dg_els(nel_dg_loc)%nx = normal_x
                   dg_els(nel_dg_loc)%ny = normal_y
                   dg_els(nel_dg_loc)%nz = normal_z                 




               endif  
               
                            
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(1 -1) +1)

                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then
 
                  !face 2: y = -1
                  ip = 0
                  nel_dg_loc = nel_dg_loc +1
                  dg_els(nel_dg_loc)%ind_el = local_el_num(ie)
                  dg_els(nel_dg_loc)%face_el = 2
                  dg_els(nel_dg_loc)%mat = tag_mat(im)
                  dg_els(nel_dg_loc)%spct_deg = nn-1
                  dg_els(nel_dg_loc)%quad_rule = ip
                  i = 1


                  do  while (i .le. 2*n_line) 
                  
                    if(con_DG(2,i) .eq. local_el_num(ie) .and. con_DG(3,i) .eq. 2) then 

                      ip = ip + 1
                      dg_els(nel_dg_loc)%quad_rule = ip
                      
                      call NEWTON_RAPSON(nodes_DG(1,i), nodes_DG(2,i), nodes_DG(3,i), &
                                        alfa11(ie),alfa12(ie),alfa13(ie), &
                                        alfa21(ie),alfa22(ie),alfa23(ie), &
                                        alfa31(ie),alfa32(ie),alfa33(ie), &
                                        beta11(ie),beta12(ie),beta13(ie), &
                                        beta21(ie),beta22(ie),beta23(ie), &
                                        beta31(ie),beta32(ie),beta33(ie), & 
                                        gamma1(ie),gamma2(ie),gamma3(ie), &
                                        delta1(ie),delta2(ie),delta3(ie), & 
                                        tt, csi, eta, zeta, nofinr, mpi_id,&
                                        con_DG(2,i),con_DG(5,i), 1.d-6, 1.01d0,1)
                                        
                                        
                                        
                        dg_els(nel_dg_loc)%x_pl(ip) = csi
                        dg_els(nel_dg_loc)%y_pl(ip) = -1.d0                       
                        dg_els(nel_dg_loc)%z_pl(ip) = zeta
                        dg_els(nel_dg_loc)%wx_pl(ip) = nodes_DG(4,i)
                        dg_els(nel_dg_loc)%wy_pl(ip) = nodes_DG(5,i)
                        dg_els(nel_dg_loc)%wz_pl(ip) = nodes_DG(6,i) 
                        
                        dg_els(nel_dg_loc)%omega_minus(ip,0) = con_DG(4,i)
                        dg_els(nel_dg_loc)%omega_minus(ip,1) = con_DG(5,i)
                        dg_els(nel_dg_loc)%omega_minus(ip,2) = con_DG(6,i)
                        dg_els(nel_dg_loc)%omega_minus(ip,3) = 0 
                         
                         
                         ! find the corresponding neigh. elem in faces matrix
                         
                         call GET_FACE_DG(faces, nel_dg_glo, con_DG(5,i), con_DG(6,i), ih)
                         
                         call MAKE_BILINEAR_MAP(area_nodes(2:25,ih), &
                                             c_alfa11, c_alfa12, c_alfa13, &
                                             c_alfa21, c_alfa22, c_alfa23, &
                                             c_alfa31, c_alfa32, c_alfa33, &  
                                             c_beta11, c_beta12, c_beta13, & 
                                             c_beta21, c_beta22, c_beta23, & 
                                             c_beta31, c_beta32, c_beta33, &
                                             c_gamma1, c_gamma2, c_gamma3, &
                                             c_delta1, c_delta2, c_delta3)
                         
                         
                         call NEWTON_RAPSON(nodes_DG(1,i), nodes_DG(2,i), nodes_DG(3,i), &
                                        c_alfa11, c_alfa12, c_alfa13, &
                                        c_alfa21, c_alfa22, c_alfa23,  &
                                        c_alfa31, c_alfa32, c_alfa33,  &
                                        c_beta11, c_beta12, c_beta13,  &
                                        c_beta21, c_beta22, c_beta23,  &
                                        c_beta31, c_beta32, c_beta33,  & 
                                        c_gamma1, c_gamma2, c_gamma3,  &
                                        c_delta1, c_delta2, c_delta3,  & 
                                        tt, csi, eta, zeta, nofinr, mpi_id,&
                                        con_DG(2,i),con_DG(5,i), 1.d-6, 1.01d0,1)   
                                        
                         select case(con_DG(6,i))    
                            case(1)
                              csi = -1.d0
                            case(2)
                              eta = -1.d0
                            case(3)
                              zeta = -1.d0
                            case(4)
                              csi = 1.d0
                            case(5)
                              eta = 1.d0
                            case(6)
                              zeta = 1.d0
                            case default
                              write(*,*) 'error in make_interface'
                            end select             

                         dg_els(nel_dg_loc)%x_mn(ip) =  csi
                         dg_els(nel_dg_loc)%y_mn(ip) =  eta
                         dg_els(nel_dg_loc)%z_mn(ip) =  zeta    
                      
                     endif 
                       
                         i = i+1                    
                         
                  enddo
 
                   call MAKE_NORMAL(2,xs(ne1), xs(ne2), xs(ne3), xs(ne4), ys(ne1), ys(ne2), ys(ne3), ys(ne4), &
                                    zs(ne1), zs(ne2), zs(ne3), zs(ne4), normal_x, normal_y, normal_z, 0, 0)

                 dg_els(nel_dg_loc)%nx = normal_x
                 dg_els(nel_dg_loc)%ny = normal_y
                 dg_els(nel_dg_loc)%nz = normal_z                 

               endif  
               
                            
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(1 -1) +1)
      
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then

                  !face 3: z = -1
                  ip = 0
                  nel_dg_loc = nel_dg_loc +1
                  dg_els(nel_dg_loc)%ind_el = local_el_num(ie)
                  dg_els(nel_dg_loc)%face_el = 3
                  dg_els(nel_dg_loc)%mat = tag_mat(im)
                  dg_els(nel_dg_loc)%spct_deg = nn-1
                  dg_els(nel_dg_loc)%quad_rule = ip
                  
                  i = 1

                  do  while (i .le. 2*n_line) 
                  
                    if(con_DG(2,i) .eq. local_el_num(ie) .and. con_DG(3,i) .eq. 3) then 


                      ip = ip + 1
                      dg_els(nel_dg_loc)%quad_rule = ip
                      
                      call NEWTON_RAPSON(nodes_DG(1,i), nodes_DG(2,i), nodes_DG(3,i), &
                                        alfa11(ie),alfa12(ie),alfa13(ie), &
                                        alfa21(ie),alfa22(ie),alfa23(ie), &
                                        alfa31(ie),alfa32(ie),alfa33(ie), &
                                        beta11(ie),beta12(ie),beta13(ie), &
                                        beta21(ie),beta22(ie),beta23(ie), &
                                        beta31(ie),beta32(ie),beta33(ie), & 
                                        gamma1(ie),gamma2(ie),gamma3(ie), &
                                        delta1(ie),delta2(ie),delta3(ie), & 
                                        tt, csi, eta, zeta, nofinr, mpi_id,&
                                        con_DG(2,i),con_DG(5,i), 1.d-6, 1.01d0,1)
                                        
                                        
                                        
                        dg_els(nel_dg_loc)%x_pl(ip) = csi
                        dg_els(nel_dg_loc)%y_pl(ip) = eta
                        dg_els(nel_dg_loc)%z_pl(ip) = -1.d0                       
                        dg_els(nel_dg_loc)%wx_pl(ip) = nodes_DG(4,i)
                        dg_els(nel_dg_loc)%wy_pl(ip) = nodes_DG(5,i)
                        dg_els(nel_dg_loc)%wz_pl(ip) = nodes_DG(6,i) 
                                             
                        dg_els(nel_dg_loc)%omega_minus(ip,0) = con_DG(4,i)
                        dg_els(nel_dg_loc)%omega_minus(ip,1) = con_DG(5,i)
                        dg_els(nel_dg_loc)%omega_minus(ip,2) = con_DG(6,i)
                        dg_els(nel_dg_loc)%omega_minus(ip,3) = 0 
                         
                         
                        ! find the corresponding neigh. elem in faces matrix
                                                 
                         call GET_FACE_DG(faces, nel_dg_glo, con_DG(5,i), con_DG(6,i), ih)
                         
                         call MAKE_BILINEAR_MAP(area_nodes(2:25,ih), &
                                             c_alfa11, c_alfa12, c_alfa13, &
                                             c_alfa21, c_alfa22, c_alfa23, &
                                             c_alfa31, c_alfa32, c_alfa33, &  
                                             c_beta11, c_beta12, c_beta13, & 
                                             c_beta21, c_beta22, c_beta23, & 
                                             c_beta31, c_beta32, c_beta33, &
                                             c_gamma1, c_gamma2, c_gamma3, &
                                             c_delta1, c_delta2, c_delta3)
                         
                         
                         call NEWTON_RAPSON(nodes_DG(1,i), nodes_DG(2,i), nodes_DG(3,i), &
                                        c_alfa11, c_alfa12, c_alfa13, &
                                        c_alfa21, c_alfa22, c_alfa23,  &
                                        c_alfa31, c_alfa32, c_alfa33,  &
                                        c_beta11, c_beta12, c_beta13,  &
                                        c_beta21, c_beta22, c_beta23,  &
                                        c_beta31, c_beta32, c_beta33,  & 
                                        c_gamma1, c_gamma2, c_gamma3,  &
                                        c_delta1, c_delta2, c_delta3,  & 
                                        tt, csi, eta, zeta, nofinr, mpi_id,&
                                        con_DG(2,i),con_DG(5,i), 1.d-6, 1.01d0,1)   
                                        
                         select case(con_DG(6,i))    
                            case(1)
                              csi = -1.d0
                            case(2)
                              eta = -1.d0
                            case(3)
                              zeta = -1.d0
                            case(4)
                              csi = 1.d0
                            case(5)
                              eta = 1.d0
                            case(6)
                              zeta = 1.d0
                            case default
                              write(*,*) 'error in make_interface'
                            end select             

                         dg_els(nel_dg_loc)%x_mn(ip) =  csi
                         dg_els(nel_dg_loc)%y_mn(ip) =  eta
                         dg_els(nel_dg_loc)%z_mn(ip) =  zeta    
                     
                       endif 
                       i = i +1
                                             
                         
                  enddo

                   call MAKE_NORMAL(3,xs(ne1), xs(ne2), xs(ne3), xs(ne4), ys(ne1), ys(ne2), ys(ne3), ys(ne4), &
                                    zs(ne1), zs(ne2), zs(ne3), zs(ne4), normal_x, normal_y, normal_z, 0, 0)

                 dg_els(nel_dg_loc)%nx = normal_x
                 dg_els(nel_dg_loc)%ny = normal_y
                 dg_els(nel_dg_loc)%nz = normal_z                 

             endif  
                          
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(nn -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(nn -1) +1)

                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then
                  
                  !face 4 : x = 1
                  ip = 0
                  nel_dg_loc = nel_dg_loc +1             
                  dg_els(nel_dg_loc)%ind_el = local_el_num(ie)
                  dg_els(nel_dg_loc)%face_el = 4
                  dg_els(nel_dg_loc)%mat = tag_mat(im)
                  dg_els(nel_dg_loc)%spct_deg = nn-1
                  dg_els(nel_dg_loc)%quad_rule = ip
                  i = 1

                  do  while (i .le. 2*n_line) 
                  
                    if(con_DG(2,i) .eq. local_el_num(ie) .and. con_DG(3,i) .eq. 4) then 


                      ip = ip + 1
                      dg_els(nel_dg_loc)%quad_rule = ip
                      
                      call NEWTON_RAPSON(nodes_DG(1,i), nodes_DG(2,i), nodes_DG(3,i), &
                                        alfa11(ie),alfa12(ie),alfa13(ie), &
                                        alfa21(ie),alfa22(ie),alfa23(ie), &
                                        alfa31(ie),alfa32(ie),alfa33(ie), &
                                        beta11(ie),beta12(ie),beta13(ie), &
                                        beta21(ie),beta22(ie),beta23(ie), &
                                        beta31(ie),beta32(ie),beta33(ie), & 
                                        gamma1(ie),gamma2(ie),gamma3(ie), &
                                        delta1(ie),delta2(ie),delta3(ie), & 
                                        tt, csi, eta, zeta, nofinr, mpi_id,&
                                        con_DG(2,i),con_DG(5,i), 1.d-6, 1.01d0,1)
                                        
                        dg_els(nel_dg_loc)%x_pl(ip) = 1.d0                        
                        dg_els(nel_dg_loc)%y_pl(ip) = eta
                        dg_els(nel_dg_loc)%z_pl(ip) = zeta
                        dg_els(nel_dg_loc)%wx_pl(ip) = nodes_DG(4,i)
                        dg_els(nel_dg_loc)%wy_pl(ip) = nodes_DG(5,i)
                        dg_els(nel_dg_loc)%wz_pl(ip) = nodes_DG(6,i) 
                        
                        dg_els(nel_dg_loc)%omega_minus(ip,0) = con_DG(4,i) 
                        dg_els(nel_dg_loc)%omega_minus(ip,1) = con_DG(5,i)
                        dg_els(nel_dg_loc)%omega_minus(ip,2) = con_DG(6,i)
                        dg_els(nel_dg_loc)%omega_minus(ip,3) = 0 
                         
                         
                        ! find the corresponding neigh. elem in faces matrix
                         
                         call GET_FACE_DG(faces, nel_dg_glo, con_DG(5,i), con_DG(6,i), ih)
                         
                         call MAKE_BILINEAR_MAP(area_nodes(2:25,ih), &
                                             c_alfa11, c_alfa12, c_alfa13, &
                                             c_alfa21, c_alfa22, c_alfa23, &
                                             c_alfa31, c_alfa32, c_alfa33, &  
                                             c_beta11, c_beta12, c_beta13, & 
                                             c_beta21, c_beta22, c_beta23, & 
                                             c_beta31, c_beta32, c_beta33, &
                                             c_gamma1, c_gamma2, c_gamma3, &
                                             c_delta1, c_delta2, c_delta3)
                         
                         
                         call NEWTON_RAPSON(nodes_DG(1,i), nodes_DG(2,i), nodes_DG(3,i), &
                                        c_alfa11, c_alfa12, c_alfa13, &
                                        c_alfa21, c_alfa22, c_alfa23,  &
                                        c_alfa31, c_alfa32, c_alfa33,  &
                                        c_beta11, c_beta12, c_beta13,  &
                                        c_beta21, c_beta22, c_beta23,  &
                                        c_beta31, c_beta32, c_beta33,  & 
                                        c_gamma1, c_gamma2, c_gamma3,  &
                                        c_delta1, c_delta2, c_delta3,  & 
                                        tt, csi, eta, zeta, nofinr, mpi_id,&
                                        con_DG(2,i),con_DG(5,i), 1.d-6, 1.01d0,1)   
                                        
                         select case(con_DG(6,i))    
                            case(1)
                              csi = -1.d0
                            case(2)
                              eta = -1.d0
                            case(3)
                              zeta = -1.d0
                            case(4)
                              csi = 1.d0
                            case(5)
                              eta = 1.d0
                            case(6)
                              zeta = 1.d0
                            case default
                              write(*,*) 'error in make_interface'
                            end select             

                         dg_els(nel_dg_loc)%x_mn(ip) =  csi
                         dg_els(nel_dg_loc)%y_mn(ip) =  eta
                         dg_els(nel_dg_loc)%z_mn(ip) =  zeta    

                        endif                        
                            i = i+1                 
                       
                  enddo
 
                   call MAKE_NORMAL(4,xs(ne1), xs(ne2), xs(ne3), xs(ne4), ys(ne1), ys(ne2), ys(ne3), ys(ne4), &
                                    zs(ne1), zs(ne2), zs(ne3), zs(ne4), normal_x, normal_y, normal_z, 0, 0)

                 dg_els(nel_dg_loc)%nx = normal_x
                 dg_els(nel_dg_loc)%ny = normal_y
                 dg_els(nel_dg_loc)%nz = normal_z                 

               endif   
               
                         
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(1 -1) +1)

      
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then
                  
                  !face 5: y = 1
                  ip = 0
                  nel_dg_loc = nel_dg_loc +1
                  dg_els(nel_dg_loc)%ind_el = local_el_num(ie)
                  dg_els(nel_dg_loc)%face_el = 5
                  dg_els(nel_dg_loc)%mat = tag_mat(im)
                  dg_els(nel_dg_loc)%spct_deg = nn-1
                  dg_els(nel_dg_loc)%quad_rule = ip
                  i = 1
 
                  do  while (i .le. 2*n_line) 


                    if(con_DG(2,i) .eq. local_el_num(ie) .and. con_DG(3,i) .eq. 5) then 
                                                              
                      ip = ip + 1
                      dg_els(nel_dg_loc)%quad_rule = ip
                       
                      call NEWTON_RAPSON(nodes_DG(1,i), nodes_DG(2,i), nodes_DG(3,i), &
                                        alfa11(ie),alfa12(ie),alfa13(ie), &
                                        alfa21(ie),alfa22(ie),alfa23(ie), &
                                        alfa31(ie),alfa32(ie),alfa33(ie), &
                                        beta11(ie),beta12(ie),beta13(ie), &
                                        beta21(ie),beta22(ie),beta23(ie), &
                                        beta31(ie),beta32(ie),beta33(ie), & 
                                        gamma1(ie),gamma2(ie),gamma3(ie), &
                                        delta1(ie),delta2(ie),delta3(ie), & 
                                        tt, csi, eta, zeta, nofinr, mpi_id,&
                                        con_DG(2,i),con_DG(5,i), 1.d-6, 1.01d0,1)
                                        
                                        
                                        
                        dg_els(nel_dg_loc)%x_pl(ip) = csi
                        dg_els(nel_dg_loc)%y_pl(ip) = 1.d0
                        dg_els(nel_dg_loc)%z_pl(ip) = zeta
                        dg_els(nel_dg_loc)%wx_pl(ip) = nodes_DG(4,i)
                        dg_els(nel_dg_loc)%wy_pl(ip) = nodes_DG(5,i)
                        dg_els(nel_dg_loc)%wz_pl(ip) = nodes_DG(6,i) 
                        
                        dg_els(nel_dg_loc)%omega_minus(ip,0) = con_DG(4,i) 
                        dg_els(nel_dg_loc)%omega_minus(ip,1) = con_DG(5,i)
                        dg_els(nel_dg_loc)%omega_minus(ip,2) = con_DG(6,i)
                        dg_els(nel_dg_loc)%omega_minus(ip,3) = 0 
                         
                         
                        ! find the corresponding neigh. elem in faces matrix
                         
                         call GET_FACE_DG(faces, nel_dg_glo, con_DG(5,i), con_DG(6,i), ih)
                         
                         call MAKE_BILINEAR_MAP(area_nodes(2:25,ih), &
                                             c_alfa11, c_alfa12, c_alfa13, &
                                             c_alfa21, c_alfa22, c_alfa23, &
                                             c_alfa31, c_alfa32, c_alfa33, &  
                                             c_beta11, c_beta12, c_beta13, & 
                                             c_beta21, c_beta22, c_beta23, & 
                                             c_beta31, c_beta32, c_beta33, &
                                             c_gamma1, c_gamma2, c_gamma3, &
                                             c_delta1, c_delta2, c_delta3)
                         
                         
                         call NEWTON_RAPSON(nodes_DG(1,i), nodes_DG(2,i), nodes_DG(3,i), &
                                        c_alfa11, c_alfa12, c_alfa13, &
                                        c_alfa21, c_alfa22, c_alfa23,  &
                                        c_alfa31, c_alfa32, c_alfa33,  &
                                        c_beta11, c_beta12, c_beta13,  &
                                        c_beta21, c_beta22, c_beta23,  &
                                        c_beta31, c_beta32, c_beta33,  & 
                                        c_gamma1, c_gamma2, c_gamma3,  &
                                        c_delta1, c_delta2, c_delta3,  & 
                                        tt, csi, eta, zeta, nofinr, mpi_id,&
                                        con_DG(2,i),con_DG(5,i), 1.d-6, 1.01d0,1)   
                                        
                         select case(con_DG(6,i))    
                            case(1)
                              csi = -1.d0
                            case(2)
                              eta = -1.d0
                            case(3)
                              zeta = -1.d0
                            case(4)
                              csi = 1.d0
                            case(5)
                              eta = 1.d0
                            case(6)
                              zeta = 1.d0
                            case default
                              write(*,*) 'error in make_interface'
                            end select             

                         dg_els(nel_dg_loc)%x_mn(ip) =  csi
                         dg_els(nel_dg_loc)%y_mn(ip) =  eta
                         dg_els(nel_dg_loc)%z_mn(ip) =  zeta    

                       endif 

                      i = i+1

                         
                  enddo

                   call MAKE_NORMAL(5,xs(ne1), xs(ne2), xs(ne3), xs(ne4), ys(ne1), ys(ne2), ys(ne3), ys(ne4), &
                                    zs(ne1), zs(ne2), zs(ne3), zs(ne4), normal_x, normal_y, normal_z, 0, 0)

                 dg_els(nel_dg_loc)%nx = normal_x
                 dg_els(nel_dg_loc)%ny = normal_y
                 dg_els(nel_dg_loc)%nz = normal_z                 
                 
 
               endif       
                       
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(1 -1) +1)

                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then
  
                  !face 6: z = 1
                  ip = 0
                  nel_dg_loc = nel_dg_loc + 1
                  dg_els(nel_dg_loc)%ind_el = local_el_num(ie)
                  dg_els(nel_dg_loc)%face_el = 6
                  dg_els(nel_dg_loc)%mat = tag_mat(im)
                  dg_els(nel_dg_loc)%spct_deg = nn-1
                  dg_els(nel_dg_loc)%quad_rule = ip
                  i = 1

                  do  while (i .le. 2*n_line) 
                  
                    if(con_DG(2,i) .eq. local_el_num(ie) .and. con_DG(3,i) .eq. 6) then 

                      ip = ip + 1
                      dg_els(nel_dg_loc)%quad_rule = ip
                      
                      call NEWTON_RAPSON(nodes_DG(1,i), nodes_DG(2,i), nodes_DG(3,i), &
                                        alfa11(ie),alfa12(ie),alfa13(ie), &
                                        alfa21(ie),alfa22(ie),alfa23(ie), &
                                        alfa31(ie),alfa32(ie),alfa33(ie), &
                                        beta11(ie),beta12(ie),beta13(ie), &
                                        beta21(ie),beta22(ie),beta23(ie), &
                                        beta31(ie),beta32(ie),beta33(ie), & 
                                        gamma1(ie),gamma2(ie),gamma3(ie), &
                                        delta1(ie),delta2(ie),delta3(ie), & 
                                        tt, csi, eta, zeta, nofinr, mpi_id,&
                                        con_DG(2,i),con_DG(5,i), 1.d-6, 1.01d0,1)
                                        
                                        
                                        
                        dg_els(nel_dg_loc)%x_pl(ip) = csi
                        dg_els(nel_dg_loc)%y_pl(ip) = eta
                        dg_els(nel_dg_loc)%z_pl(ip) = 1.d0                        
                        dg_els(nel_dg_loc)%wx_pl(ip) = nodes_DG(4,i)
                        dg_els(nel_dg_loc)%wy_pl(ip) = nodes_DG(5,i)
                        dg_els(nel_dg_loc)%wz_pl(ip) = nodes_DG(6,i) 
                        


                         dg_els(nel_dg_loc)%omega_minus(ip,0) = con_DG(4,i)
                         dg_els(nel_dg_loc)%omega_minus(ip,1) = con_DG(5,i)
                         dg_els(nel_dg_loc)%omega_minus(ip,2) = con_DG(6,i)
                         dg_els(nel_dg_loc)%omega_minus(ip,3) = 0 
                         
                         
                        ! find the corresponding neigh. elem in faces matrix
                         
                         call GET_FACE_DG(faces, nel_dg_glo, con_DG(5,i), con_DG(6,i), ih)
                         
                         call MAKE_BILINEAR_MAP(area_nodes(2:25,ih), &
                                             c_alfa11, c_alfa12, c_alfa13, &
                                             c_alfa21, c_alfa22, c_alfa23, &
                                             c_alfa31, c_alfa32, c_alfa33, &  
                                             c_beta11, c_beta12, c_beta13, & 
                                             c_beta21, c_beta22, c_beta23, & 
                                             c_beta31, c_beta32, c_beta33, &
                                             c_gamma1, c_gamma2, c_gamma3, &
                                             c_delta1, c_delta2, c_delta3)
                         
                         
                         call NEWTON_RAPSON(nodes_DG(1,i), nodes_DG(2,i), nodes_DG(3,i), &
                                        c_alfa11, c_alfa12, c_alfa13, &
                                        c_alfa21, c_alfa22, c_alfa23,  &
                                        c_alfa31, c_alfa32, c_alfa33,  &
                                        c_beta11, c_beta12, c_beta13,  &
                                        c_beta21, c_beta22, c_beta23,  &
                                        c_beta31, c_beta32, c_beta33,  & 
                                        c_gamma1, c_gamma2, c_gamma3,  &
                                        c_delta1, c_delta2, c_delta3,  & 
                                        tt, csi, eta, zeta, nofinr, mpi_id,&
                                        con_DG(2,i),con_DG(5,i), 1.d-6, 1.01d0,1)   
                                        
                         select case(con_DG(6,i))    
                            case(1)
                              csi = -1.d0
                            case(2)
                              eta = -1.d0
                            case(3)
                              zeta = -1.d0
                            case(4)
                              csi = 1.d0
                            case(5)
                              eta = 1.d0
                            case(6)
                              zeta = 1.d0
                            case default
                              write(*,*) 'error in make_interface'
                            end select            

                         dg_els(nel_dg_loc)%x_mn(ip) =  csi
                         dg_els(nel_dg_loc)%y_mn(ip) =  eta
                         dg_els(nel_dg_loc)%z_mn(ip) =  zeta    

                     endif 
                         
                       i = i+1
                  enddo 
                 
                   call MAKE_NORMAL(6,xs(ne1), xs(ne2), xs(ne3), xs(ne4), ys(ne1), ys(ne2), ys(ne3), ys(ne4), &
                                    zs(ne1), zs(ne2), zs(ne3), zs(ne4), normal_x, normal_y, normal_z, 0, 0)

                 dg_els(nel_dg_loc)%nx = normal_x
                 dg_els(nel_dg_loc)%ny = normal_y
                 dg_els(nel_dg_loc)%nz = normal_z                 
                 

               endif               
            endif
         enddo
      enddo
      
      
      deallocate(con_DG, nodes_DG)


!   Finding/ordering neigh. elements

     do it = 1, nel_dg_loc
            allocate(copia(max_quad_points,0:3))    
            copia = dg_els(it)%omega_minus
            
            call GET_NEIGHBOUR_ELEM(copia, dg_els(it)%quad_rule, nofne_el, mpi_id, max_quad_points)
            
            dg_els(it)%nofne = nofne_el                
            dg_els(it)%omega_minus(:,3) = copia(:,3)    

            deallocate(copia)
     enddo
    
    
!CHECK IF ALL DG EL HAVE A NEIGHBOUR   
     do it = 1, nel_dg_loc
        if(dg_els(it)%nofne == 0) then 
          write(*,*) 'ATTENTION : element ', dg_els(it)%ind_el, ' has 0 neighbouring elements'
        endif   
     enddo
    
     
     do it = 1, nel_dg_loc

        ic = 1
        do while(ic .le. dg_els(it)%nofne)
            do ip = 1, (dg_els(it)%quad_rule)  
               if(dg_els(it)%omega_minus(ip,3) .eq. ic) then
                  dg_els(it)%conf(ic,0) = dg_els(it)%omega_minus(ip,0)
                  dg_els(it)%conf(ic,1) = dg_els(it)%omega_minus(ip,1)
                  dg_els(it)%conf(ic,2) = dg_els(it)%omega_minus(ip,2)
                                   
               endif
            enddo
            ic = ic + 1
        enddo
     enddo


!*************************************************************************************************************
!                                   JUMP FOR DG INTERFACES
!*************************************************************************************************************        
      
      !Needs OpenMP for speed-up
      !Remove OpenMP from MAKE_LOC_MATRIX_DG.f90

!!$OMP PARALLEL &
!!$OMP PRIVATE(ie, nnz_p, nnz_m, ielem, iface, nn, im, ic1, ne1, ic) &
!!$OMP PRIVATE(imne, iene, ifacene, mm, i, lambda, mu, ic2, JP, JM, pen_h, pen_p, pen) &
!!$OMP PRIVATE(cp_a, cp_b, cp_c, cp_d, cp_e, cp_f, cp_g, cp_n, cp_p, mpi_id, j) &
!!$OMP PRIVATE(c_alfa11, c_alfa12, c_alfa13, c_alfa21, c_alfa22, c_alfa23) &
!!$OMP PRIVATE(c_alfa31, c_alfa32, c_alfa33, c_beta11, c_beta12, c_beta13) &
!!$OMP PRIVATE(c_beta21, c_beta22, c_beta23, c_beta31, c_beta32, c_beta33) &
!!$OMP PRIVATE(c_gamma1, c_gamma2, c_gamma3, c_delta1, c_delta2, c_delta3) &
!!$OMP PRIVATE(k, I4S, ishift, jshift, J4S, M4S)

!!$OMP DO
      
      do ie = 1, nel_dg_loc
         nnz_p = 0
         nnz_m = 0   
         el_new(ie)%nnz_col = 0     
         if(testmode .eq. 1) then
              el_new(ie)%nnz_col_only_uv = 0; nnz_p_only_uv = 0; nnz_m_only_uv = 0   
         endif
         
         ielem = dg_els(ie)%ind_el
         iface = dg_els(ie)%face_el
         nn = dg_els(ie)%spct_deg + 1
         im = dg_els(ie)%mat
                  
         el_new(ie)%ind = ielem
         el_new(ie)%face = iface
         el_new(ie)%deg = nn   
         el_new(ie)%mate = im

         el_new(ie)%num_of_ne = dg_els(ie)%nofne

         call GET_FACE_DG(faces, nel_dg_glo, ielem, iface, ic1)
         call GET_INDLOC_FROM_INDGLO(local_el_num, ne_loc, ielem, ne1)
               
         allocate(el_new(ie)%matP(3*(nn**3),3*(nn**3)));  el_new(ie)%matP = 0.d0
         allocate(el_new(ie)%matM(dg_els(ie)%nofne))
         !--------------------------------------------------------------------------------------------
         if (testmode .eq. 1) then       
               allocate(el_new(ie)%matP_only_uv(3*(nn**3),3*(nn**3)));  el_new(ie)%matP_only_uv = 0.d0
               allocate(el_new(ie)%matM_only_uv(dg_els(ie)%nofne))
         endif      
         !--------------------------------------------------------------------------------------------
    
         do ic = 1, dg_els(ie)%nofne                   

            imne = dg_els(ie)%conf(ic,0)              
            iene = dg_els(ie)%conf(ic,1)
            ifacene = dg_els(ie)%conf(ic,2)
            
            mm = 2

            do i = 1, nm 
              if(tag_mat(i) .eq. imne ) then
                   mm = sd(i) +1
                   lambda = 2.d0*prop_mat(im,2)*prop_mat(i,2)/(prop_mat(im,2) + prop_mat(i,2))
                   mu = 2.d0*prop_mat(im,3)*prop_mat(i,3)/(prop_mat(im,3) + prop_mat(i,3))                                     
              endif
            enddo   
            
            el_new(ie)%el_conf(ic,0) = imne           
            el_new(ie)%el_conf(ic,1) = iene
            el_new(ie)%el_conf(ic,2) = ifacene
                
            call GET_FACE_DG(faces, nel_dg_glo, iene, ifacene, ic2)
            
            allocate(JP(3*(nn**3),3*(nn**3)),el_new(ie)%matM(ic)%MJUMP(3*(nn**3),3*(mm**3)),JM(3*(nn**3),3*(mm**3)))
            !---------------------------------------------------------------------------------------------------
            if(testmode .eq. 1) then 
              allocate(JP_only_uv(3*(nn**3),3*(nn**3)),&
                       el_new(ie)%matM(ic)%MJUMP_only_uv(3*(nn**3),3*(mm**3)),JM_only_uv(3*(nn**3),3*(mm**3)))           
            endif 
            !---------------------------------------------------------------------------------------------------
            
            pen_h = min(area_nodes(1,ic1), area_nodes(1,ic2))                      
            pen_p = max(nn-1,mm-1)                     

            pen = penalty_c * (lambda + 2.d0*mu) * pen_p**2.d0 / pen_h
                       
            cp_a = 0.5d0*(lambda+2.d0*mu)*dg_els(ie)%nx
            cp_b = 0.5d0*lambda*dg_els(ie)%nx
            cp_c = 0.5d0*mu*dg_els(ie)%ny
            cp_d = 0.5d0*mu*dg_els(ie)%nz
            cp_e = 0.5d0*mu*dg_els(ie)%nx
            cp_f = 0.5d0*(lambda+2.d0*mu)*dg_els(ie)%ny
            cp_g = 0.5d0*lambda*dg_els(ie)%ny
            cp_n = 0.5d0*(lambda+2.d0*mu)*dg_els(ie)%nz
            cp_p = 0.5d0*lambda*dg_els(ie)%nz

            call MAKE_BILINEAR_MAP(area_nodes(2:25,ic2), &
                                    c_alfa11, c_alfa12, c_alfa13, &
                                    c_alfa21, c_alfa22, c_alfa23, &
                                    c_alfa31, c_alfa32, c_alfa33, &  
                                    c_beta11, c_beta12, c_beta13, & 
                                    c_beta21, c_beta22, c_beta23, & 
                                    c_beta31, c_beta32, c_beta33, &
                                    c_gamma1, c_gamma2, c_gamma3, &
                                    c_delta1, c_delta2, c_delta3)



            call MAKE_LOC_MATRIX_DG(dg_els(ie)%x_pl, dg_els(ie)%y_pl, dg_els(ie)%z_pl, &
                       dg_els(ie)%wx_pl, dg_els(ie)%wy_pl, dg_els(ie)%wz_pl,&
                       dg_els(ie)%x_mn, dg_els(ie)%y_mn, dg_els(ie)%z_mn, &
                       dg_els(ie)%quad_rule, nn, mm, &
                       dg_els(ie)%omega_minus(1:dg_els(ie)%quad_rule,0:3), &
                       alfa11(ne1),alfa12(ne1),alfa13(ne1),&
                       alfa21(ne1),alfa22(ne1),alfa23(ne1),&
                       alfa31(ne1),alfa32(ne1),alfa33(ne1),&
                       beta11(ne1),beta12(ne1),beta13(ne1),&
                       beta21(ne1),beta22(ne1),beta23(ne1),&
                       beta31(ne1),beta32(ne1),beta33(ne1),&
                       gamma1(ne1),gamma2(ne1),gamma3(ne1),&
                       delta1(ne1),delta2(ne1),delta3(ne1),&
                       iene, c_alfa11, c_alfa12, c_alfa13,&
                       c_alfa21,c_alfa22,c_alfa23,&
                       c_alfa31,c_alfa32,c_alfa33,&
                       c_beta11,c_beta12,c_beta13,&
                       c_beta21,c_beta22,c_beta23,&
                       c_beta31,c_beta32,c_beta33,&
                       c_gamma1,c_gamma2,c_gamma3,&
                       c_delta1,c_delta2,c_delta3,&
                       cp_a,cp_b,cp_c,cp_d,cp_e,cp_f,cp_g,cp_n,cp_p, &
                       pen, dg_cnst, JP, JM, mpi_id, pen_h, testmode,&
                       JP_only_uv, JM_only_uv)

                       el_new(ie)%matP = el_new(ie)%matP + JP
                       el_new(ie)%matM(ic)%MJUMP = JM

                       !--------------------------------------------------------------------
                       if(testmode .eq. 1) then
                           el_new(ie)%matP_only_uv = el_new(ie)%matP_only_uv + JP_only_uv
                           el_new(ie)%matM(ic)%MJUMP_only_uv = JM_only_uv
                       endif       
                       !--------------------------------------------------------------------


                        do i = 1, 3*nn**3
                           do j = 1, 3*mm**3
                                   if(el_new(ie)%matM(ic)%MJUMP(i,j) .ne. 0.d0) nnz_m = nnz_m + 1
                                !------------------------------------------------------------------------------
                                   if(testmode .eq. 1) then
                                     if (el_new(ie)%matM(ic)%MJUMP_only_uv(i,j) .ne. 0.d0) &
                                      nnz_m_only_uv = nnz_m_only_uv + 1
                                   endif   
                                !------------------------------------------------------------------------------
                           enddo
                        enddo       
                        
                        el_new(ie)%nnz_col = el_new(ie)%nnz_col + 3*mm**3 
                        
                        !--------------------------------------------------------------------
                        if(testmode .eq. 1) then
                           el_new(ie)%nnz_col_only_uv = el_new(ie)%nnz_col_only_uv + 3*mm**3
                           deallocate(JM_only_uv,JP_only_uv)
                        endif
                        !--------------------------------------------------------------------
                           
                        deallocate(JM,JP)
        enddo


        do i = 1, 3*nn**3
           do j = 1, 3*nn**3
              if(el_new(ie)%matP(i,j) .ne. 0.d0) nnz_p = nnz_p + 1
           enddo
        enddo       
                    
        el_new(ie)%nnz_minus = nnz_m;  el_new(ie)%nnz_plus = nnz_p
        
        !--------------------------------------------------------------------------------------------            
        if(testmode .eq. 1) then
           do i = 1, 3*nn**3
              do j = 1, 3*nn**3
                 if(el_new(ie)%matP_only_uv(i,j) .ne. 0.d0) nnz_p_only_uv = nnz_p_only_uv + 1
              enddo
           enddo       
           el_new(ie)%nnz_minus_only_uv = nnz_m_only_uv;  el_new(ie)%nnz_plus_only_uv = nnz_p_only_uv
        endif                   
        !---------------------------------------------------------------------------------------------

!
! STORING MATRIX IN A SPARSE FORMAT
!

         nn = el_new(ie)%deg 
         allocate(el_new(ie)%IPlus(0:3*nn**3), el_new(ie)%IMin(0:3*nn**3))
         allocate(el_new(ie)%JPlus(el_new(ie)%nnz_plus),el_new(ie)%matPlus(el_new(ie)%nnz_plus))
         allocate(el_new(ie)%JMin(el_new(ie)%nnz_minus),el_new(ie)%matMin(el_new(ie)%nnz_minus))
         allocate(J4S(el_new(ie)%nnz_minus),M4S(el_new(ie)%nnz_minus))
         allocate(I4S(el_new(ie)%nnz_plus))

         k = 1
         do i = 1, 3*nn**3
            do j = 1, 3*nn**3
               if( el_new(ie)%matP(i,j) .ne. 0.d0 ) then
                  I4S(k) = i
                  el_new(ie)%JPlus(k) = j
                  el_new(ie)%matPlus(k) =  el_new(ie)%matP(i,j)
                  k = k + 1
                endif
             enddo
         enddo        

         deallocate(el_new(ie)%matP)

         el_new(ie)%IPlus= 0
         do i = 1, el_new(ie)%nnz_plus
            el_new(ie)%IPlus(I4S(i)) = el_new(ie)%IPlus(I4S(i)) + 1
         enddo
         do i = 1, 3*nn**3
            el_new(ie)%IPlus(i) = el_new(ie)%IPlus(i) + el_new(ie)%IPlus(i-1)
         enddo

         deallocate(I4S)
         allocate(I4S(el_new(ie)%nnz_minus))

         k = 1
         ishift = 0
         jshift = 0
         
         do ic = 1, el_new(ie)%num_of_ne

            mm = 2
            do i = 1, nm 
              if(tag_mat(i) .eq. el_new(ie)%el_conf(ic,0) ) then
                   mm = sd(i) +1
              endif
            enddo   
            
          
            do i = 1, 3*nn**3
               do j = 1, 3*mm**3
                  if( el_new(ie)%matM(ic)%MJUMP(i,j) .ne. 0.d0 ) then
                     I4S(k) = i 
                     J4S(k) = j + jshift
                     M4S(k) = el_new(ie)%matM(ic)%MJUMP(i,j)
                     
                     k = k + 1
                   endif

                enddo
            enddo        
            

            jshift = jshift + 3*mm**3

            deallocate(el_new(ie)%matM(ic)%MJUMP)
            
         
         enddo

         k = 1
         do i = 1, 3*nn**3
            do j = 1, el_new(ie)%nnz_minus
              if(I4S(j) .eq. i) then          
                 el_new(ie)%JMin(k) = J4S(j) 
                 el_new(ie)%matMin(k) = M4S(j)
                 k = k + 1
               endif
             enddo    
         enddo

           
         deallocate(J4S, M4S)  

         el_new(ie)%IMin = 0
         do i = 1, el_new(ie)%nnz_minus
            el_new(ie)%IMin(I4S(i)) = el_new(ie)%IMin(I4S(i)) + 1
         enddo
         
         do i = 1, 3*nn**3
            el_new(ie)%IMin(i) = el_new(ie)%IMin(i) + el_new(ie)%IMin(i-1)
         enddo

         deallocate(I4S)
                    
         !-----------------------------------------------------------------------------------------
         if(testmode .eq. 1) then  

           nn = el_new(ie)%deg 
           allocate(el_new(ie)%IPlus_only_uv(0:3*nn**3),el_new(ie)%IMin_only_uv(0:3*nn**3))
           allocate(el_new(ie)%JPlus_only_uv(el_new(ie)%nnz_plus_only_uv))
           allocate(el_new(ie)%matPlus_only_uv(el_new(ie)%nnz_plus_only_uv))
           allocate(el_new(ie)%JMin_only_uv(el_new(ie)%nnz_minus_only_uv))
           allocate(el_new(ie)%matMin_only_uv(el_new(ie)%nnz_minus_only_uv))
           allocate(J4S(el_new(ie)%nnz_minus_only_uv),M4S(el_new(ie)%nnz_minus_only_uv))
           allocate(I4S(el_new(ie)%nnz_plus_only_uv))

           k = 1
           do i = 1, 3*nn**3
              do j = 1, 3*nn**3
                 if( el_new(ie)%matP_only_uv(i,j) .ne. 0.d0 ) then
                    I4S(k) = i
                    el_new(ie)%JPlus_only_uv(k) = j
                    el_new(ie)%matPlus_only_uv(k) =  el_new(ie)%matP_only_uv(i,j)
                    k = k + 1
                 endif
               enddo
           enddo        

           deallocate(el_new(ie)%matP_only_uv)

           el_new(ie)%IPlus_only_uv= 0
           do i = 1, el_new(ie)%nnz_plus_only_uv
              el_new(ie)%IPlus_only_uv(I4S(i)) = el_new(ie)%IPlus_only_uv(I4S(i)) + 1
           enddo
           do i = 1, 3*nn**3
              el_new(ie)%IPlus_only_uv(i) = el_new(ie)%IPlus_only_uv(i) + el_new(ie)%IPlus_only_uv(i-1)
           enddo

           deallocate(I4S)
           allocate(I4S(el_new(ie)%nnz_minus_only_uv))

           k = 1
           ishift = 0
           jshift = 0
         
           do ic = 1, el_new(ie)%num_of_ne

              mm = 2
              do i = 1, nm 
                if(tag_mat(i) .eq. el_new(ie)%el_conf(ic,0) ) then
                     mm = sd(i) +1
                endif
              enddo   
            
          
              do i = 1, 3*nn**3
                 do j = 1, 3*mm**3
                    if( el_new(ie)%matM(ic)%MJUMP_only_uv(i,j) .ne. 0.d0 ) then
                       I4S(k) = i 
                       J4S(k) = j + jshift
                       M4S(k) = el_new(ie)%matM(ic)%MJUMP_only_uv(i,j)
                     
                       k = k + 1
                     endif
                 enddo
               enddo        
   
               jshift = jshift + 3*mm**3

               deallocate(el_new(ie)%matM(ic)%MJUMP_only_uv)
           enddo

           k = 1
           do i = 1, 3*nn**3
              do j = 1, el_new(ie)%nnz_minus_only_uv
                if(I4S(j) .eq. i) then          
                   el_new(ie)%JMin_only_uv(k) = J4S(j) 
                   el_new(ie)%matMin_only_uv(k) = M4S(j)
                   k = k + 1
                endif
              enddo    
           enddo

           deallocate(J4S, M4S)  

           el_new(ie)%IMin_only_uv = 0
           do i = 1, el_new(ie)%nnz_minus_only_uv
              el_new(ie)%IMin_only_uv(I4S(i)) = el_new(ie)%IMin_only_uv(I4S(i)) + 1
           enddo
         
           do i = 1, 3*nn**3
              el_new(ie)%IMin_only_uv(i) = el_new(ie)%IMin_only_uv(i) + el_new(ie)%IMin_only_uv(i-1)
           enddo

           deallocate(I4S)

 
          endif
          !-----------------------------------------------------------------------------------------

      enddo

!!$OMP END DO
!!$OMP END PARALLEL


     deallocate(dg_els)
         


    return

    end subroutine MAKE_DG_INTERFACE


