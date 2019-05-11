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

!> @brief Setup for DG faces. Stores data on arrays faces and area_nodes. 
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] nm  number of materials
!> @param[in] sd  polynomial degree vector
!> @param[in] tag_mat labels for materials
!> @param[in] cs_nnz_loc length of cs_loc
!> @param[in] cs_loc local spectral connectivity vector
!> @param[in] nn_loc number of local nodes
!> @param[in] local_n_num local node numeration
!> @param[in] ne_loc number of local elements
!> @param[in] local_el_num local element numeration
!> @param[in] xs x-coord of spectral nodes
!> @param[in] ys y-coord of spectral nodes
!> @param[in] zs z-coord of spectral nodes 
!> @param[in] nel_dg_loc  number of local dg elements
!> @param[in] nel_dg_glo  number of global dg elements
!> @param[in] i4count vector identifying nodes where DG conditions are applied
!> @param[in] mpi_id  mpi processor identity
!> @param[in] mpi_comm  mpi communicator
!> @param[in] alfa11 costant for bilinear mapping
!> @param[in] alfa12 costant for bilinear mapping
!> @param[in] alfa13 costant for bilinear mapping
!> @param[in] alfa21 costant for bilinear mapping
!> @param[in] alfa22 costant for bilinear mapping
!> @param[in] alfa23 costant for bilinear mapping
!> @param[in] alfa31 costant for bilinear mapping
!> @param[in] alfa32 costant for bilinear mapping
!> @param[in] alfa33 costant for bilinear mapping 
!> @param[in] beta11 costant for bilinear mapping
!> @param[in] beta12 costant for bilinear mapping
!> @param[in] beta13 costant for bilinear mapping
!> @param[in] beta21 costant for bilinear mapping
!> @param[in] beta22 costant for bilinear mapping
!> @param[in] beta23 costant for bilinear mapping
!> @param[in] beta31 costant for bilinear mapping
!> @param[in] beta32 costant for bilinear mapping
!> @param[in] beta33 costant for bilinear mapping
!> @param[in] gamma1 costant for bilinear mapping
!> @param[in] gamma2 costant for bilinear mapping
!> @param[in] gamma3 costant for bilinear mapping
!> @param[in] delta1 costant for bilinear mapping
!> @param[in] delta2 costant for bilinear mapping
!> @param[in] delta3 costant for bilinear mapping
!> @param[out] faces  material, element, face for a DG face
!> @param[out] area_nodes info about DG faces.
!!      area_nodes(1,i) = area of the face, 
!!      area_nodes(2,i),...,area_nodes(25,i) -> costants for the bilinear map 

      subroutine SETUP_DG(nm, sd, tag_mat, cs_nnz_loc, cs_loc, &
                           nn_loc, local_n_num, ne_loc, local_el_num, &
                           xs,ys,zs,&
                           nel_dg_loc, nel_dg_glo, &
                           i4count, mpi_id, mpi_comm, &
                           alfa11, alfa12, alfa13, &
                           alfa21, alfa22, alfa23, &
                           alfa31, alfa32, alfa33, &
                           beta11, beta12, beta13, &
                           beta21, beta22, beta23, &
                           beta31, beta32, beta33, &
                           gamma1, gamma2, gamma3, &
                           delta1, delta2, delta3, & 
                           faces, area_nodes)


     implicit none

     include 'SPEED.MPI'

     integer*4 :: nm, cs_nnz_loc, nn_loc, ne_loc, nel_dg_loc, nel_dg_glo
     integer*4 :: im, nn, ie, ned, err_out
     integer*4 :: ne1, ne2, ne3, ne4, ic1, ic2, ic3, ic4
     integer*4 :: ne5, ne6, ne7, ne8, ic5, ic6, ic7, ic8
     integer*4 :: mpi_comm, mpi_id, mpi_ierr

     integer*4, dimension(nm) :: tag_mat, sd
     integer*4, dimension(0:cs_nnz_loc) :: cs_loc
     integer*4, dimension(nn_loc) :: local_n_num, i4count
     integer*4, dimension(ne_loc) :: local_el_num

     integer*4, dimension(3,nel_dg_loc), intent(inout) :: faces

     real*8 :: surf

     real*8, dimension(nn_loc) :: xs,ys,zs
     real*8, dimension(ne_loc) :: alfa11,alfa12,alfa13
     real*8, dimension(ne_loc) :: alfa21,alfa22,alfa23
     real*8, dimension(ne_loc) :: alfa31,alfa32,alfa33
     real*8, dimension(ne_loc) :: beta11,beta12,beta13
     real*8, dimension(ne_loc) :: beta21,beta22,beta23
     real*8, dimension(ne_loc) :: beta31,beta32,beta33
     real*8, dimension(ne_loc) :: gamma1,gamma2,gamma3
     real*8, dimension(ne_loc) :: delta1,delta2,delta3
     
     real*8, dimension(25,nel_dg_loc), intent(inout) :: area_nodes


!*****************************************************************************************
!loading of data structure faces containing the description for DG faces
!*****************************************************************************************      

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
                  nel_dg_loc = nel_dg_loc +1
                  
                  call GET_AREA_FACE(xs(ne1), xs(ne2), xs(ne3), xs(ne4), ys(ne1), ys(ne2), ys(ne3), ys(ne4), &
                                              zs(ne1), zs(ne2), zs(ne3), zs(ne4), surf, ie)
                          
                  ! face x = 1                            
                                              
                  faces(1,nel_dg_loc) = tag_mat(im)
                  faces(2,nel_dg_loc) = local_el_num(ie)
                  faces(3,nel_dg_loc) = 1

                  area_nodes(1,nel_dg_loc) = surf
                  area_nodes(2,nel_dg_loc) = alfa11(ie)
                  area_nodes(3,nel_dg_loc) = alfa12(ie)
                  area_nodes(4,nel_dg_loc) = alfa13(ie)
                  
                  area_nodes(5,nel_dg_loc) = alfa21(ie)
                  area_nodes(6,nel_dg_loc) = alfa22(ie)
                  area_nodes(7,nel_dg_loc) = alfa23(ie)
                  
                  area_nodes(8,nel_dg_loc) = alfa31(ie)
                  area_nodes(9,nel_dg_loc) = alfa32(ie)
                  area_nodes(10,nel_dg_loc) = alfa33(ie)
                  
                  area_nodes(11,nel_dg_loc) = beta11(ie)
                  area_nodes(12,nel_dg_loc) = beta12(ie)
                  area_nodes(13,nel_dg_loc) = beta13(ie)
                  
                  area_nodes(14,nel_dg_loc) = beta21(ie)
                  area_nodes(15,nel_dg_loc) = beta22(ie)
                  area_nodes(16,nel_dg_loc) = beta23(ie)
                  
                  area_nodes(17,nel_dg_loc) = beta31(ie)
                  area_nodes(18,nel_dg_loc) = beta32(ie)
                  area_nodes(19,nel_dg_loc) = beta33(ie)
                  
                  area_nodes(20,nel_dg_loc) = gamma1(ie)
                  area_nodes(21,nel_dg_loc) = gamma2(ie)
                  area_nodes(22,nel_dg_loc) = gamma3(ie)
                  
                  area_nodes(23,nel_dg_loc) = delta1(ie)
                  area_nodes(24,nel_dg_loc) = delta2(ie)
                  area_nodes(25,nel_dg_loc) = delta3(ie)                                    
                     
               endif               

               ! face y = - 1 
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(1 -1) +1)

                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then
                  nel_dg_loc = nel_dg_loc +1
                  
                  call GET_AREA_FACE(xs(ne1), xs(ne2), xs(ne3), xs(ne4), ys(ne1), ys(ne2), ys(ne3), ys(ne4), &
                                              zs(ne1), zs(ne2), zs(ne3), zs(ne4), surf, ie)

                  ! face y = 1

                  faces(1,nel_dg_loc) = tag_mat(im)
                  faces(2,nel_dg_loc) = local_el_num(ie)
                  faces(3,nel_dg_loc) = 2
                 
                  area_nodes(1,nel_dg_loc) = surf
                  area_nodes(2,nel_dg_loc) = alfa11(ie)
                  area_nodes(3,nel_dg_loc) = alfa12(ie)
                  area_nodes(4,nel_dg_loc) = alfa13(ie)
                  
                  area_nodes(5,nel_dg_loc) = alfa21(ie)
                  area_nodes(6,nel_dg_loc) = alfa22(ie)
                  area_nodes(7,nel_dg_loc) = alfa23(ie)
                  
                  area_nodes(8,nel_dg_loc) = alfa31(ie)
                  area_nodes(9,nel_dg_loc) = alfa32(ie)
                  area_nodes(10,nel_dg_loc) = alfa33(ie)
                  
                  area_nodes(11,nel_dg_loc) = beta11(ie)
                  area_nodes(12,nel_dg_loc) = beta12(ie)
                  area_nodes(13,nel_dg_loc) = beta13(ie)
                  
                  area_nodes(14,nel_dg_loc) = beta21(ie)
                  area_nodes(15,nel_dg_loc) = beta22(ie)
                  area_nodes(16,nel_dg_loc) = beta23(ie)
                  
                  area_nodes(17,nel_dg_loc) = beta31(ie)
                  area_nodes(18,nel_dg_loc) = beta32(ie)
                  area_nodes(19,nel_dg_loc) = beta33(ie)
                  
                  area_nodes(20,nel_dg_loc) = gamma1(ie)
                  area_nodes(21,nel_dg_loc) = gamma2(ie)
                  area_nodes(22,nel_dg_loc) = gamma3(ie)
                  
                  area_nodes(23,nel_dg_loc) = delta1(ie)
                  area_nodes(24,nel_dg_loc) = delta2(ie)
                  area_nodes(25,nel_dg_loc) = delta3(ie)                                    

               endif               

               !face z = - 1
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(1 -1) +1)

               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then
                  nel_dg_loc = nel_dg_loc +1
                  
                  call GET_AREA_FACE(xs(ne1), xs(ne2), xs(ne3), xs(ne4), ys(ne1), ys(ne2), ys(ne3), ys(ne4), &
                                              zs(ne1), zs(ne2), zs(ne3), zs(ne4), surf, ie)

                  ! face z = 1

                  faces(1,nel_dg_loc) = tag_mat(im)
                  faces(2,nel_dg_loc) = local_el_num(ie)
                  faces(3,nel_dg_loc) = 3
                  
                  area_nodes(1,nel_dg_loc) = surf
                  area_nodes(2,nel_dg_loc) = alfa11(ie)
                  area_nodes(3,nel_dg_loc) = alfa12(ie)
                  area_nodes(4,nel_dg_loc) = alfa13(ie)
                  
                  area_nodes(5,nel_dg_loc) = alfa21(ie)
                  area_nodes(6,nel_dg_loc) = alfa22(ie)
                  area_nodes(7,nel_dg_loc) = alfa23(ie)
                  
                  area_nodes(8,nel_dg_loc) = alfa31(ie)
                  area_nodes(9,nel_dg_loc) = alfa32(ie)
                  area_nodes(10,nel_dg_loc) = alfa33(ie)
                  
                  area_nodes(11,nel_dg_loc) = beta11(ie)
                  area_nodes(12,nel_dg_loc) = beta12(ie)
                  area_nodes(13,nel_dg_loc) = beta13(ie)
                  
                  area_nodes(14,nel_dg_loc) = beta21(ie)
                  area_nodes(15,nel_dg_loc) = beta22(ie)
                  area_nodes(16,nel_dg_loc) = beta23(ie)
                  
                  area_nodes(17,nel_dg_loc) = beta31(ie)
                  area_nodes(18,nel_dg_loc) = beta32(ie)
                  area_nodes(19,nel_dg_loc) = beta33(ie)
                  
                  area_nodes(20,nel_dg_loc) = gamma1(ie)
                  area_nodes(21,nel_dg_loc) = gamma2(ie)
                  area_nodes(22,nel_dg_loc) = gamma3(ie)
                  
                  area_nodes(23,nel_dg_loc) = delta1(ie)
                  area_nodes(24,nel_dg_loc) = delta2(ie)
                  area_nodes(25,nel_dg_loc) = delta3(ie)                                    

 
               endif               

               ! face x = 1 
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(nn -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(nn -1) +1)

                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then
                  nel_dg_loc = nel_dg_loc +1
                  
                  call GET_AREA_FACE(xs(ne1), xs(ne2), xs(ne3), xs(ne4), ys(ne1), ys(ne2), ys(ne3), ys(ne4), &
                                              zs(ne1), zs(ne2), zs(ne3), zs(ne4), surf, ie)

                  ! face x = - 1
                              
                  faces(1,nel_dg_loc) = tag_mat(im)
                  faces(2,nel_dg_loc) = local_el_num(ie)
                  faces(3,nel_dg_loc) = 4
                  
                  area_nodes(1,nel_dg_loc) = surf
                  area_nodes(2,nel_dg_loc) = alfa11(ie)
                  area_nodes(3,nel_dg_loc) = alfa12(ie)
                  area_nodes(4,nel_dg_loc) = alfa13(ie)
                  
                  area_nodes(5,nel_dg_loc) = alfa21(ie)
                  area_nodes(6,nel_dg_loc) = alfa22(ie)
                  area_nodes(7,nel_dg_loc) = alfa23(ie)
                  
                  area_nodes(8,nel_dg_loc) = alfa31(ie)
                  area_nodes(9,nel_dg_loc) = alfa32(ie)
                  area_nodes(10,nel_dg_loc) = alfa33(ie)
                  
                  area_nodes(11,nel_dg_loc) = beta11(ie)
                  area_nodes(12,nel_dg_loc) = beta12(ie)
                  area_nodes(13,nel_dg_loc) = beta13(ie)
                  
                  area_nodes(14,nel_dg_loc) = beta21(ie)
                  area_nodes(15,nel_dg_loc) = beta22(ie)
                  area_nodes(16,nel_dg_loc) = beta23(ie)
                  
                  area_nodes(17,nel_dg_loc) = beta31(ie)
                  area_nodes(18,nel_dg_loc) = beta32(ie)
                  area_nodes(19,nel_dg_loc) = beta33(ie)
                  
                  area_nodes(20,nel_dg_loc) = gamma1(ie)
                  area_nodes(21,nel_dg_loc) = gamma2(ie)
                  area_nodes(22,nel_dg_loc) = gamma3(ie)
                  
                  area_nodes(23,nel_dg_loc) = delta1(ie)
                  area_nodes(24,nel_dg_loc) = delta2(ie)
                  area_nodes(25,nel_dg_loc) = delta3(ie)                                    
          
                  
               endif               

               ! face y = 1
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(1 -1) +1)

                           
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then
                  nel_dg_loc = nel_dg_loc +1
                  
                  call GET_AREA_FACE(xs(ne1), xs(ne2), xs(ne3), xs(ne4), ys(ne1), ys(ne2), ys(ne3), ys(ne4), &
                                              zs(ne1), zs(ne2), zs(ne3), zs(ne4), surf, ie)

                  faces(1,nel_dg_loc) = tag_mat(im)
                  faces(2,nel_dg_loc) = local_el_num(ie)
                  faces(3,nel_dg_loc) = 5

                  area_nodes(1,nel_dg_loc) = surf
                  area_nodes(2,nel_dg_loc) = alfa11(ie)
                  area_nodes(3,nel_dg_loc) = alfa12(ie)
                  area_nodes(4,nel_dg_loc) = alfa13(ie)
                  
                  area_nodes(5,nel_dg_loc) = alfa21(ie)
                  area_nodes(6,nel_dg_loc) = alfa22(ie)
                  area_nodes(7,nel_dg_loc) = alfa23(ie)
                  
                  area_nodes(8,nel_dg_loc) = alfa31(ie)
                  area_nodes(9,nel_dg_loc) = alfa32(ie)
                  area_nodes(10,nel_dg_loc) = alfa33(ie)
                  
                  area_nodes(11,nel_dg_loc) = beta11(ie)
                  area_nodes(12,nel_dg_loc) = beta12(ie)
                  area_nodes(13,nel_dg_loc) = beta13(ie)
                  
                  area_nodes(14,nel_dg_loc) = beta21(ie)
                  area_nodes(15,nel_dg_loc) = beta22(ie)
                  area_nodes(16,nel_dg_loc) = beta23(ie)
                  
                  area_nodes(17,nel_dg_loc) = beta31(ie)
                  area_nodes(18,nel_dg_loc) = beta32(ie)
                  area_nodes(19,nel_dg_loc) = beta33(ie)
                  
                  area_nodes(20,nel_dg_loc) = gamma1(ie)
                  area_nodes(21,nel_dg_loc) = gamma2(ie)
                  area_nodes(22,nel_dg_loc) = gamma3(ie)
                  
                  area_nodes(23,nel_dg_loc) = delta1(ie)
                  area_nodes(24,nel_dg_loc) = delta2(ie)
                  area_nodes(25,nel_dg_loc) = delta3(ie)                                    

               endif               

               ! face z = 1
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(1 -1) +1) 
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(1 -1) +1) 

                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then
                  nel_dg_loc = nel_dg_loc +1
                  
                  call GET_AREA_FACE(xs(ne1), xs(ne2), xs(ne3), xs(ne4), ys(ne1), ys(ne2), ys(ne3), ys(ne4), &
                                              zs(ne1), zs(ne2), zs(ne3), zs(ne4), surf, ie)


                  !face z = - 1

                  faces(1,nel_dg_loc) = tag_mat(im)
                  faces(2,nel_dg_loc) = local_el_num(ie)
                  faces(3,nel_dg_loc) = 6
                  
                  area_nodes(1,nel_dg_loc) = surf
                  area_nodes(2,nel_dg_loc) = alfa11(ie)
                  area_nodes(3,nel_dg_loc) = alfa12(ie)
                  area_nodes(4,nel_dg_loc) = alfa13(ie)
                  
                  area_nodes(5,nel_dg_loc) = alfa21(ie)
                  area_nodes(6,nel_dg_loc) = alfa22(ie)
                  area_nodes(7,nel_dg_loc) = alfa23(ie)
                  
                  area_nodes(8,nel_dg_loc) = alfa31(ie)
                  area_nodes(9,nel_dg_loc) = alfa32(ie)
                  area_nodes(10,nel_dg_loc) = alfa33(ie)
                  
                  area_nodes(11,nel_dg_loc) = beta11(ie)
                  area_nodes(12,nel_dg_loc) = beta12(ie)
                  area_nodes(13,nel_dg_loc) = beta13(ie)
                  
                  area_nodes(14,nel_dg_loc) = beta21(ie)
                  area_nodes(15,nel_dg_loc) = beta22(ie)
                  area_nodes(16,nel_dg_loc) = beta23(ie)
                  
                  area_nodes(17,nel_dg_loc) = beta31(ie)
                  area_nodes(18,nel_dg_loc) = beta32(ie)
                  area_nodes(19,nel_dg_loc) = beta33(ie)
                  
                  area_nodes(20,nel_dg_loc) = gamma1(ie)
                  area_nodes(21,nel_dg_loc) = gamma2(ie)
                  area_nodes(22,nel_dg_loc) = gamma3(ie)
                  
                  area_nodes(23,nel_dg_loc) = delta1(ie)
                  area_nodes(24,nel_dg_loc) = delta2(ie)
                  area_nodes(25,nel_dg_loc) = delta3(ie)                                    

               endif               
            endif
         enddo
      enddo
      


      return

      end subroutine SETUP_DG
