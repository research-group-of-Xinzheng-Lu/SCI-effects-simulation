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


!> @author Lorenzo Gaborini
!> @date October, 2014 - Creation
!> @version 1.0
!> @brief  SPEED exit codes
!> @warning  Fortran exit statement requires gfortran
module speed_exit_codes

implicit none

integer, parameter :: EXIT_NORMAL         = 0
integer, parameter :: EXIT_CFL            = 1
integer, parameter :: EXIT_INSTAB         = 2
integer, parameter :: EXIT_ANELASTIC      = 3
integer, parameter :: EXIT_SETUP          = 4
integer, parameter :: EXIT_SINGULARMTX    = 5
integer, parameter :: EXIT_SURF_NOTFOUND  = 6
integer, parameter :: EXIT_ENERGY_ERROR   = 7
integer, parameter :: EXIT_SYNTAX_ERROR   = 8
integer, parameter :: EXIT_MISSING_FILE   = 9
integer, parameter :: EXIT_ROOT           = 10
integer, parameter :: EXIT_ELEM_ORIENT    = 11
integer, parameter :: EXIT_NO_NODES       = 12
integer, parameter :: EXIT_NO_ELEMENTS    = 13
integer, parameter :: EXIT_DAMPING_PEAK   = 14
integer, parameter :: EXIT_NO_MATERIALS   = 15
integer, parameter :: EXIT_FUNCTION_ERROR = 16

end module speed_exit_codes


!> @brief Set maximal bounds.
!! @author Ilario Mazzieri
!> @date September, 2013 - Creation
!> @version 1.0
module max_var
        integer, parameter :: max_dim           = 3      !<dimension of the problem
        integer, parameter :: max_el_conf       = 200    !<max number of neighbouring elements 
        integer, parameter :: nofqp             = 8      !<max number of 1-D quadrature point per element 
        integer, parameter :: nofinr            = 500    !<max number of newton rapson iterations
        integer, parameter :: max_quad_points   = 8000   !<max number of quadrature nodes on a DG surface
end module max_var


!> @author Ilario Mazzieri
!> @date September, 2013 - Creation
!> @version 1.0
!> @brief Contains mesh structure (scratch)
module str_mesh_scratch

  use max_var , only: nofqp !<max number of 1-D quadrature point per element

  type scratch_ELEMENT !<contains coordinates of quadrature nodes
     real*8, dimension(nofqp**2) :: x_nq !<x coordinate
     real*8, dimension(nofqp**2) :: y_nq !<y coordinate
     real*8, dimension(nofqp**2) :: z_nq !<z coordinate
  end type scratch_ELEMENT

end module str_mesh_scratch


!> @author Ilario Mazzieri
!> @date September, 2013 - Creation
!> @version 1.0
!> @brief Contains mesh structure for DG interface elements 
module str_mesh

 use max_var

 type ELEMENT  !< Interface DG Element (quad)
   
   integer*4:: mat          !< block id (material)
   integer*4:: ind_el             !< element index (global numeration)
   integer*4:: face_el      !< face index (from 1 to 6)
   integer*4:: spct_deg     !< polynomila degree inside the element 
   integer*4:: quad_rule    !< number of quadrature point in 1D
   integer*4:: nofne        !< number of neighbouring elements
   integer*4:: proj_yn      !< 1 if the element project quad nodes 0 otherwise 
   real*8   :: nx,ny,nz     !< normal to the element
   

   real*8, dimension(nofqp**2) :: wx_pl  !< weights of the quadrature rule (x)
   real*8, dimension(nofqp**2) :: wy_pl  !< weights of the quadrature rule (y)
   real*8, dimension(nofqp**2) :: wz_pl  !< weights of the quadrature rule (z)

 end type ELEMENT

END MODULE str_mesh


!> @author Ilario Mazzieri
!> @date September, 2013 - Creation
!> @version 1.0
!> @brief Contains mesh structure for DG elements after pre-processing 
module str_mesh_after

 use max_var,  only : max_quad_points, max_el_conf

 type ELEMENT_after
   
   integer*4:: mat          !< block id (material)
   integer*4:: ind_el             !< element index (global numeration)
   integer*4:: face_el      !< face index (from 1 to 6)
   integer*4:: spct_deg     !< polynomial degree 
   integer*4:: quad_rule    !< number of quadrature point in 1D
   integer*4:: nofne        !< number of neighbouring elements
   real*8   :: nx,ny,nz     !< normal to the element

   real*8, dimension(max_quad_points) :: x_pl  !< quadrature points x (+,+)
   real*8, dimension(max_quad_points) :: y_pl  !< quadrature points y (+,+)
   real*8, dimension(max_quad_points) :: z_pl  !< quadrature points z (+,+)
   real*8, dimension(max_quad_points) :: x_mn  !< quadrature points x (+,-)
   real*8, dimension(max_quad_points) :: y_mn  !< quadrature points y (+,-)
   real*8, dimension(max_quad_points) :: z_mn  !< quadrature points z (+,-)

   real*8, dimension(max_quad_points) :: wx_pl !< quadrature weights x
   real*8, dimension(max_quad_points) :: wy_pl !< quadrature weights y
   real*8, dimension(max_quad_points) :: wz_pl !< quadrature weights z
 
   integer*4, dimension(max_quad_points,0:3) :: omega_minus !< matrix containing neigh el. info (quad points)
   integer*4, dimension(max_el_conf,0:2) :: conf            !< matrix containing neigh el. info (mat,el,ind,face)


 end type ELEMENT_after


END MODULE str_mesh_after



!> @author Ilario Mazzieri
!> @date September, 2013 - Creation
!> @version 1.0
!> @brief Contains structure for jump matrices
module DGJUMP

use max_var , only: max_el_conf

type matrix 

   real*8, dimension(:,:), pointer :: MJUMP !< jump matrix   
   real*8, dimension(:,:), pointer :: MJUMP_only_uv  !<jump matrix fot testmode (only [u][v])

end type matrix


type el4loop   !< element structure for time loop (RAM saving)

   integer*4:: ind                                   !< element index (gl. num.)
   integer*4:: face                                  !< element face  (1<-->6)
   integer*4:: deg                                   !< pol.degree
   integer*4:: mate                                  !< material
   integer*4:: num_of_ne                             !< number of neigh. elements
   integer*4:: nnz_plus                              !< nonzero els. jump matrix (+,+)
   integer*4:: nnz_minus                             !< nonzero els. jump matrix (+,-) 
   integer*4:: nnz_col                               !< length of u_m vector in timeloop
   integer*4:: nnz_plus_only_uv                      !< nonzero els. jump matrix (+,+) [u][v] - testmode
   integer*4:: nnz_minus_only_uv                     !< nonzero els. jump matrix (+,-) [u][v] - testmode
   integer*4:: nnz_col_only_uv                       !< length of u_m vector in timeloop [u][v] - testmode
   
   integer*4, dimension(:), pointer :: IPlus, JPlus  !< RCS format rows
   integer*4, dimension(:), pointer :: IMin, JMin    !< RCS format columns
   integer*4, dimension(max_el_conf,0:2) :: el_conf  !< matrix for neigh. elements
   real*8, dimension(:), pointer :: matPlus, matMin  !< RCS format matrix
   real*8, dimension(:,:), pointer :: matP           !< jump matrix (+,+)
   type(matrix), dimension(:), pointer :: matM       !< jump matrix (+,-)

   integer*4, dimension(:), pointer :: IPlus_only_uv, JPlus_only_uv  !< RCS format rows [u][v] - testmode
   integer*4, dimension(:), pointer :: IMin_only_uv, JMin_only_uv    !< RCS format columns [u][v] - testmode
   real*8, dimension(:), pointer :: matPlus_only_uv, matMin_only_uv  !< RCS format matrix [u][v] - testmode
   real*8, dimension(:,:), pointer :: matP_only_uv                   !< jump matrix (+,+) [u][v] - testmode
   type(matrix), dimension(:), pointer :: matM_only_uv               !< jump matrix (+,-) [u][v] - testmode


end type 

end module


!> @author Ilario Mazzieri
!> @date September, 2013 - Creation
!> @version 1.0
!> @brief Contains  SPEED PARAMETERS used in (SPEED, READ_INPUT_FILES,
!! MAKE_PARTION_AND_MPI_FILES, MAKE_SETUP_MPI_CONFORMING, MAKE_SPX_GRID_WITH_LOC_NUMERATION
!! MAKE_SEISMIC_MOMENT_OR_EXPLOSIVE_SOURCE, MAKE_LOAD_MATRIX, FIND_MONITOR_POSITION
!! MAKE_BOUNDARY_CONDITIONS, TIME_LOOP, DEALLOCATE_VARIABLES)
module speed_par

! EXIT CODES
      use speed_exit_codes

!************************************************************************************
!                                    LOGICAL
!************************************************************************************

!'TRUE' OR 'FALSE'
      logical :: filefound  

      ! if .TRUE., fails on negative anelastic coefficients (see: damping 2)
      logical :: b_failoncoeffs

      ! if .TRUE., do not start the TIME_LOOP
      ! only setup mesh, parameters, CFL, etc. then quit
      logical :: b_setuponly

      ! if .TRUE:, quit if CFL does not hold
      logical :: b_failCFL

      ! if .TRUE., quit if simulation becomes unstable
      logical :: b_instabilitycontrol

      ! Default values
      logical, parameter :: b_failoncoeffs_default = .FALSE.
      logical, parameter :: b_setuponly_default = .FALSE.
      logical, parameter :: b_failCFL_default = .FALSE.
      logical, parameter :: b_instabilitycontrol_default = .FALSE.
      
!************************************************************************************
!                                    CHARACTERS
!************************************************************************************

! FILENAMES
      character*70 :: file_PG, file_MPGM, file_LS, file_MLST, &
                      head_file, grid_file, mat_file, bkp_file, &
                      monitor_file, monitor_file_new, file_face, file_part, & 
                      file_mpi, mpi_file, file_mpi_new, &
                      filename

! 'YES' OR 'NOT'
      character*3  :: deltat_fixed   

! RUNGE-KUTTA
      character*10 :: rk_scheme

! SLIP DISTRIBUTION 'STD' = HERRERO, 'ARC' = ARCHULETA, 'GAL' = GALLOVIC  
      character*3  :: slip_type

!************************************************************************************
!                                    INTEGERS
!************************************************************************************

! COUNTERS
      integer*4 :: im_pgm, im_lst, &
                   im, ie, i, j, k, h, in, ic, id, &
                   ip, ielem, i4t, im_nle


! MAX RANGE NUMBER  
      integer*4 :: nnod_macro, nnod, nelem, con_nnz, con_nnz_bc, &
                   con_nnz_loc, con_nnz_bc_loc, con_nnz_dg, &
                   nnz_node_weight, nface, nmat, n_case, & 
                   nload_dirX_el, nload_neuX_el, nload_poiX_el, nload_forX_el, &
                   nload_dirY_el, nload_neuY_el, nload_poiY_el, nload_forY_el, &
                   nload_dirZ_el, nload_neuZ_el, nload_poiZ_el, nload_forZ_el, &
                   nload_neuN_el, nload_sism_el, nload_expl_el, & 
                   nload_plaX_el, nload_plaY_el, nload_plaZ_el, &
                   nload_pres_el, nload_shea_el, nload_forc_el, &
                   nload_traX_el, nload_traY_el, nload_traZ_el, &
                   nload_abc_el, nload_dg_el, &
                   nfunc, nfunc_data, &
                   nmonitors_pgm, nmonitors_lst, &
                   num_pgm, num_lst, &
                   nts, nsnaps, &
                   ns, nn, nn2, nn3, &
                   nnode_dirX,nnode_dirY,nnode_dirZ, &
                   nnode_abc, nelem_abc, nnode_dg, nelem_dg, nelem_dg_glo, &
                   nnod_loc, nelem_loc, nface_loc, &
                   dime_js, dime_jr, &
                   max_num_node_sism, length_check_node_sism, &                
                   max_num_node_expl, length_check_node_expl, &
                   nsend,nrecv, nsend_jump, nrecv_jump, &
                   nnode_dom, nelem_dom, edgecut, &
                   nmat_nle, total_els, nvec, &
                   nargs, ntime_err, nb_tra_load, n_test
                   
! 0/1 INTERGERS
      integer*4 :: file_mon_pgm, file_mon_lst, &
                   find, torf, trof, make_damping_yes_or_not, &
                   mpi_ierr
      
! DAMPING
      integer*4 :: damping_type
      integer*4, parameter :: N_SLS = 3    
            
! OTHER       
      integer*4 :: trash, &
                   opt_out_form, opt_out_data, &  
                   unit_mpi, unit_part, &
                   initial_snap, &
                   mpi_id, mpi_np, mpi_comm, &
                   ncase, ndt_mon_lst, ndt_mon_pgm, &
                   val_case, tag_case, &
                   iargc, &
                   rk_order, rk_stages, testmode, debug
                   
! FIXED DIMENSION VECTORS
      integer*4, dimension(3)  :: clock  !CLOCK (OBSOLTE)
      integer*4, dimension (6) :: opt_out_var ! SELECT OUTPUT VARIABLES   

! CONNECTIVIY VECTORS
      integer*4, dimension (:), allocatable :: con_spx, con_spx_dg, con_spx_bc, &
                                               node_weight, node_pointer, &
                                               con_spx_loc, con_spx_bc_loc, &
                                               local_el_num, local_node_num, local_node_num_dg, &
                                               n_glo, el_glo, &
                                               node_domain, elem_domain, &
                                               node_domain_loc, elem_domain_loc, &
                                               count_faces, &
                                               elem_index, node_index_seq, & !OBSOLETE
                                               inode_dirX, inode_dirY, inode_dirZ

! VARIABLES DEFINED IN FILENAME.MATE                                
      integer*4, dimension (:), allocatable :: &
                 sub_tag_all, &                    !FOR MULTI-NOTHONORING (Not yet implemented)                        
                 sdeg_mat, type_mat, tag_mat, &
                 fun_dirX_el, fun_neuX_el, fun_dirY_el, fun_neuY_el, fun_dirZ_el, fun_neuZ_el, &
                 fun_neuN_el, & 
                 fun_poiX_el, fun_forX_el, fun_poiY_el, fun_forY_el, fun_poiZ_el, fun_forZ_el, &
                 fun_plaX_el, fun_plaY_el, fun_plaZ_el, &
                 fun_pres_el, fun_shea_el, fun_forc_el, &
                 fun_sism_el, fun_expl_el, &
                 fun_traX_el, fun_traY_el, fun_traZ_el, &
                 fun_test, &
                 tag_dirX_el, tag_neuX_el, tag_dirY_el, tag_neuY_el, tag_dirZ_el, tag_neuZ_el, &
                 tag_neuN_el, &
                 tag_plaY_el, tag_plaX_el, tag_plaZ_el, &
                 tag_abc_el, tag_dg_el, tag_dg_yn, &
                 tag_sism_el, tag_expl_el, &  !tag_case         
                 tag_func, func_type, func_indx, &
                 node_tra
      
      
! OTHER  (MONITOR)
      integer*4, dimension (:), allocatable :: n_monitor_pgm, el_monitor_pgm, &
                                               n_monitor_lst, el_monitor_lst, & 
                                               monit_files

! OTHER (SEISMIC MOMENT OR EXPLOSIVE SOURCE)
      integer*4, dimension (:), allocatable :: num_node_sism, num_node_expl, &
                                               sism_el_glo, expl_el_glo
         
! OTHER (MPI SETUP)      
      integer*4, dimension (:), allocatable :: mpi_stat, &
                                               node_send, node_recv, proc_send, proc_recv, &
                                               node_send_jump, node_recv_jump, &
                                               proc_send_jump, proc_recv_jump      

! OTHER
      integer*4, dimension (:), allocatable :: itersnap, vec, i4count, &
                                               type_mat_nle, tag_mat_nle



! MATRICES
      integer*4, dimension (:,:), allocatable :: con, con_bc, & !val_case         
                                                 ielem_abc, faces, &        
                                                 sour_node_sism, check_node_sism, &                
                                                 sour_node_expl, check_node_expl, & 
                                                 prop_mat_nle 

!************************************************************************************
!                                    REALS
!************************************************************************************

     real*8 :: &  
               deltat, deltat_cfl, tstart, tstop, time_in_seconds, start1, start2, &
               start, finish, &                        ! TIME VARIABLES
               fmax, fpeak, &                          ! FREQUENCY                      
               xx_macro, yy_macro, zz_macro, &         ! NODES
               depth_search_mon_pgm, rotation_angle_mon_pgm, depth_search_mon_lst, &        ! MONITORS
               dg_c, pen_c, &                          ! DG CONSTANTS
               dist, eps, r8t, tol_case                ! OTHER 
               
               
! VECTOR WITH FIXED DIMENSION                
      real*8, dimension(5) :: prova, sum_prova    ! OBSOLTE
     
     
! SPECTRAL NODES & MAPPING TRANFORMATION      
      real*8, dimension (:), allocatable ::  & 
              zs_elev, zs_all, vs_tria, thick, &        !NOT HONORING                                          
              xx_spx, yy_spx, zz_spx, & !SPECTRAL NODES
              xx_spx_loc, yy_spx_loc, zz_spx_loc, &   !SPECTRAL NODES    
              alfa11, alfa12, alfa13, alfa21, alfa22, alfa23, alfa31, alfa32, alfa33, & !MAP
              beta11, beta12, beta13, beta21, beta22, beta23, beta31, beta32, beta33, & !MAP
              gamma1, gamma2, gamma3, delta1, delta2, delta3 !MAP

! MONITORS
      real*8, dimension (:), allocatable :: &
              x_monitor_pgm, y_monitor_pgm, z_monitor_pgm, &         
              x_monitor_lst, y_monitor_lst, z_monitor_lst, &
              x_monitor_real, y_monitor_real, z_monitor_real, &
              xr_monitor_pgm, yr_monitor_pgm, zr_monitor_pgm, &
              xr_monitor_lst, yr_monitor_lst, zr_monitor_lst, & 
              dist_monitor_lst, dist_monitor_pgm, &
              xr_glo, yr_glo, zr_glo, dist_glo, &
              x_glo_real, y_glo_real, z_glo_real, &
              highest_pgm, highest_pgm_loc, &
              highest_lst, highest_lst_loc, & 
              dist_el_glo, dist_tra

! DAMPING
      real*8, dimension(:), allocatable :: QS, QP, frequency_range
      real*8, dimension(:,:), allocatable:: Y_lambda,Y_mu
      real*8, dimension(:), allocatable :: A0_ray, A1_ray


! OTHER 
      real*8, dimension (:), allocatable :: tref_mat, & !tol_case
                                            func_data, tsnap, set_initial_snap, &      
                                            Cel,KCel, &
                                            time_error

! LOAD VECTOR
     ! real*8, dimension (:,:), allocatable :: Fel
      !!!!!!!!!!!!!!!!!!modified by ty 170410
      real*8, dimension (:,:,:), allocatable :: Fel
      !!!!!!!!!!!!!!!!!!modified by ty 170410

! (MATRICES OF) VALUES DEFINED IN FILENAME.MATE
      real*8, dimension (:,:), allocatable :: &
              val_dirX_el, val_neuX_el, val_dirY_el, val_neuY_el, &
              val_dirZ_el, val_neuZ_el, val_neuN_el, &
              val_poiX_el, val_forX_el, val_poiY_el, val_forY_el, &
              val_poiZ_el, val_forZ_el, val_plaX_el, val_plaY_el, val_plaZ_el, &
              val_pres_el, val_shea_el, val_forc_el, val_sism_el, val_expl_el, &
              val_traX_el, val_traY_el, val_traZ_el, &
              val_mat_nle, prop_mat

              
! (MATRICES FOR) SEISMIC MOMENT OR EXPLOSIVE SOURCE                       
      real*8, dimension (:,:), allocatable :: &
              factor_seismic_moment, tau_seismic_moment, dist_sour_node_sism, check_dist_node_sism, &        
              factor_explosive_source, dist_sour_node_expl, check_dist_node_expl                

! AREA DG FACES
      real*8, dimension (:,:), allocatable :: area_nodes

! OTHER
      real*8, dimension (:,:), allocatable :: max_u, max_v, max_a, max_o                                    

! INSTABILITY CONTROL
      real*8 :: instability_maxval
      ! default value
      real*8, parameter :: instability_maxval_default = 1E20

end module speed_par



!> @author Ilario Mazzieri
!> @date September, 2013 - Creation
!> @version 1.0
!> @brief Contains SPEED paramters (used in  MAKE_DG_INTERFACE_CONDITIONS) 
module speed_par_dg

      use max_var
      use str_mesh 
      use str_mesh_scratch                   
      use DGJUMP
      use speed_par, only: nelem_dg
      use speed_exit_codes

      type(el4loop), dimension(:), allocatable :: el_new                        
      type(ELEMENT), dimension(:), allocatable :: dg_els
      type(scratch_ELEMENT), dimension(:), allocatable :: scratch_dg_els

end module speed_par_dg



!> @author Ilario Mazzieri
!> @date September, 2013 - Creation
!> @version 1.0
!> @brief Contains a subset of SPEED paramters (used in TIME_LOOP) 
module speed_timeloop


  use speed_par, only:  &
                        !SLIP DISTRIBUTION
                         slip_type, &
  
                        !CONNECTIVITY 
                         nnod, nnode_dom,  nnod_loc, local_node_num, xx_spx_loc, yy_spx_loc, zz_spx_loc, &
                         nelem, nelem_loc, local_el_num, con_nnz_loc, con_spx_loc,&
                         nelem_dg, nelem_dg_glo, con_spx_dg, con_nnz_dg, local_node_num_dg, &                      
                        
                        !(NLE) MATERIALS
                         nmat, tag_mat, type_mat, sdeg_mat, tref_mat, prop_mat, &
                              nmat_nle, tag_mat_nle, type_mat_nle, prop_mat_nle, val_mat_nle, fpeak, &
                        
                        !EXTERNAL LOADS
                         nload_traX_el, nload_traY_el, nload_traZ_el, &                  
                         node_tra, dist_tra, nb_tra_load, &
                         val_traX_el, val_traY_el, val_traZ_el, &
                         fun_traX_el, fun_traY_el, fun_traZ_el, &          
                         nfunc, tag_func, func_type, func_indx, func_data, nfunc_data, Fel, &
                        
                         !COORDINATE TRANFORMATION                       
                         alfa11, alfa12, alfa13, alfa21, alfa22, alfa23, &
                         alfa31, alfa32, alfa33, beta11, beta12, beta13, &
                         beta21, beta22, beta23, beta31, beta32, beta33, &
                         gamma1, gamma2, gamma3, delta1, delta2, delta3, &
                        
                         !BOUNDARY CONDITIONS   
                         nnode_dirX, inode_dirX, nnode_dirY, inode_dirY, &
                         nnode_dirZ, inode_dirZ, nelem_abc, ielem_abc, &          
                        
                         !TIME INTEGRATION
                         nts, deltat, tstart, tstop, &
                         nsnaps, itersnap, initial_snap, &
                         rk_scheme, rk_order, rk_stages, & 
                        
                         !MPI
                         mpi_np, mpi_id, mpi_comm, mpi_stat, mpi_ierr, &
                         nsend, node_send, nrecv, node_recv, proc_send, proc_recv, &
                         nsend_jump, node_send_jump, nrecv_jump, node_recv_jump, &
                         proc_send_jump, proc_recv_jump, &
                         
                          !SEISMIC MOMENT
                         check_node_sism, check_dist_node_sism, &                                                  
                         length_check_node_sism, nload_sism_el, factor_seismic_moment, tau_seismic_moment, &
                           check_node_expl, check_dist_node_expl, &        
                         length_check_node_expl,nload_expl_el,factor_explosive_source, &
                         
                          !DAMPING
                          make_damping_yes_or_not, Y_lambda,Y_mu, N_SLS, damping_type, frequency_range, &
                          A0_ray, A1_ray,fmax, &
                        
                         !CASE
                         n_case,tag_case,val_case, & 
                         zs_elev, zs_all, vs_tria, thick, sub_tag_all, &
                        
                         !MONITOR & MAP 
                         ndt_mon_pgm, rotation_angle_mon_pgm, &                                                
                         nmonitors_pgm, n_monitor_pgm, el_monitor_pgm, &                                        
                         xr_monitor_pgm, yr_monitor_pgm, zr_monitor_pgm, &                                        
                         ndt_mon_lst, nmonitors_lst, n_monitor_lst, el_monitor_lst, &                                        
                         xr_monitor_lst, yr_monitor_lst, zr_monitor_lst, &                                        
                         opt_out_var, monitor_file, bkp_file, &
                            
                         !TESTMODE
                         testmode, ntime_err, time_error, debug, &
                         
                         !INSTABILITY CONTROL
                         b_instabilitycontrol, instability_maxval, &

                         !File
                         file_mpi, mpi_file

  use speed_exit_codes
end module speed_timeloop



!> @author Ilario Mazzieri
!> @date September, 2013 - Creation
!> @version 1.0
!> @brief  Quick-sort algorithm
module qsort_c_module    

implicit none
public :: QsortC
private :: Partition

contains

recursive subroutine QsortC(A)
  integer, intent(in out), dimension(:) :: A
  integer :: iq

  if(size(A) > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1))
     call QsortC(A(iq:))
  endif
end subroutine QsortC

subroutine Partition(A, marker)
  integer, intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j
  integer :: temp
  integer :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition

end module qsort_c_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ty!!!!!!!!!!!!!!!!!!!!!!!
Module MDOF_MyData
    implicit none
	type BlgInfo
		integer*4 :: ID	
		integer*4 :: NDOF	
		integer*4 :: StrucType	
        integer*4 :: ndt
        integer*4,allocatable:: MDOFyield(:,:),MDOFIDeath(:,:)
        integer*4,allocatable :: MDOFstate(:,:)
		real*8 :: height	
		real*8 :: Floor_h	
		real*8 :: T1	
		real*8 :: T2	
        real*8 :: TN
		real*8 :: Area	
		real*8 :: dt,dt2	
        real*8 :: Calpha,Cbeta
		real*8,allocatable :: props(:,:,:),dval(:,:) 
        real*8,allocatable :: BlgM(:,:),BlgK(:,:),BlgC(:,:)
        real*8,allocatable :: BlgM_Inv(:,:)
        real*8,allocatable :: MDOFEt(:,:), MDOFstatev(:,:,:), MDOFspd(:,:)
        real*8,allocatable :: BlgMaxIDR(:,:)
        real*8,allocatable:: MDOFIDR(:,:),dMDOFIDR(:,:),MDOFItF(:,:)
        real*8,allocatable:: tempMDOFU0(:,:),tempMDOFU1(:,:),tempMDOFU(:,:),tempMDOFA1(:,:)
	end type
    type(BlgInfo),allocatable:: Blg(:)
    integer*4 :: Maxdof,NumofBlgs,MDOFBlgi,MDOFBlgti
    integer*4 :: countyi,countyj,MDOFmonitor,MDOFbf,MDOFef
    integer*4 :: MDOFmonitor1,MDOFmonitor2,MDOFmonitor3
    integer*4 :: BRoutput(6)!Disp,Acc,Damage,Input,F,M 
    !For Disp, Acc: 2-every story,1-top
    !For Damage: 2-every story, 1-total 
    integer*4,allocatable::Blgdof(:),NewOrNot(:)
	real*8 ::  MasspArea
    real*8 :: kclose,configtmp 
    real*8, dimension (:), allocatable :: utyty1,utyty2,utyty3,fty1,fty2 
    real*8,allocatable:: MDOFacl(:,:)
    character(len=100):: MDOFftest0,MDOFftest,MDOFfile,MDOFfinfo,MDOFstr 
    character(len=100):: MDOFfdisp1,MDOFfdisp2 
    character(len=100):: MDOFfacc1,MDOFfacc2
    character(len=100):: MDOFfaccip1,MDOFfaccip2 
    character(len=100):: MDOFfreat1,MDOFfreat2
    character(len=100):: MDOFfreat4,MDOFfreat5 
    character(len=100):: MDOFfdamag1,MDOFfdamag2 
    !!!!!!!!!!!!!!!!!loop!!!!!!!!!!!!!!!!!!!! 
    integer*4 :: MDOFexi,MDOFexj,MDOFMnum     
    integer*4,dimension(:), allocatable :: MDOFftranfg   
    real*8,dimension(:), allocatable :: MDOFinput,MDOFinputfk     
    real*8,dimension(:), allocatable ::  MDOFinputbuffer,MDOFinputfkbuffer    
    character(len=100) :: mdoftest,mdoftest2
    real*8,dimension(:),allocatable:: MDOFforceinputbuffer,MDOFforceinput 
    !!!!!!!!!!!!!!!!!loop!!!!!!!!!!!!!!!!!!!! 
end module MDOF_MyData
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ty!!!!!!!!!!!!!!!!!!!!!!!