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

!> @brief Deallocates variables. 
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0

      subroutine DEALLOCATE_VARIABLES()
        
              
      use max_var
      use str_mesh 
      use str_mesh_scratch                   
      use DGJUMP

      use speed_par
      use speed_par_dg

      implicit none
 
      include 'SPEED.MPI' 
      
      deallocate(node_send, node_recv, proc_send, proc_recv)
      deallocate(tag_mat,type_mat,sdeg_mat,tref_mat,prop_mat,QS,QP)
      deallocate(xx_spx_loc,yy_spx_loc,zz_spx_loc)
      deallocate(con_spx_loc)
      deallocate(alfa11,alfa12,alfa13,alfa21,alfa22,alfa23,alfa31,alfa32,alfa33)
      deallocate(beta11,beta12,beta13,beta21,beta22,beta23,beta31,beta32,beta33)
      deallocate(gamma1,gamma2,gamma3,delta1,delta2,delta3)
      deallocate(local_node_num,local_el_num)
      deallocate(Fel)


      if (nmonitors_lst.gt.0) deallocate(x_monitor_lst,y_monitor_lst,z_monitor_lst,n_monitor_lst)
      if (nmonitors_lst.gt.0) deallocate(xr_monitor_lst,yr_monitor_lst,zr_monitor_lst,el_monitor_lst)
      if (nmonitors_pgm.gt.0) deallocate(x_monitor_pgm,y_monitor_pgm,z_monitor_pgm,n_monitor_pgm)
      if (nmonitors_pgm.gt.0) deallocate(xr_monitor_pgm,yr_monitor_pgm,zr_monitor_pgm,el_monitor_pgm)
      if (nmonitors_pgm .gt.0) deallocate(max_u, max_v, max_a, max_o)
      if (nsnaps.gt.0) deallocate(tsnap,itersnap)
      if (nload_dirX_el.gt.0) deallocate (val_dirX_el,fun_dirX_el,tag_dirX_el)
      if (nload_dirY_el.gt.0) deallocate (val_dirY_el,fun_dirY_el,tag_dirY_el)
      if (nload_dirZ_el.gt.0) deallocate (val_dirZ_el,fun_dirZ_el,tag_dirZ_el)
      if (nnode_dirX.gt.0) deallocate(inode_dirX)
      if (nnode_dirY.gt.0) deallocate(inode_dirY)
      if (nnode_dirZ.gt.0) deallocate(inode_dirZ)
      if (nload_neuX_el.gt.0) deallocate (val_neuX_el,fun_neuX_el,tag_neuX_el)
      if (nload_neuY_el.gt.0) deallocate (val_neuY_el,fun_neuY_el,tag_neuY_el)
      if (nload_neuZ_el.gt.0) deallocate (val_neuZ_el,fun_neuZ_el,tag_neuZ_el)
      if (nload_poiX_el.gt.0) deallocate (val_poiX_el,fun_poiX_el)
      if (nload_poiY_el.gt.0) deallocate (val_poiY_el,fun_poiY_el)
      if (nload_poiZ_el.gt.0) deallocate (val_poiZ_el,fun_poiZ_el)
      if (nload_forX_el.gt.0) deallocate (val_forX_el,fun_forX_el)
      if (nload_forY_el.gt.0) deallocate (val_forY_el,fun_forY_el)
      if (nload_forZ_el.gt.0) deallocate (val_forZ_el,fun_forZ_el)
      if (nload_pres_el.gt.0) deallocate (val_pres_el,fun_pres_el)
      if (nload_shea_el.gt.0) deallocate (val_shea_el,fun_shea_el)
      if (nload_abc_el.gt.0) deallocate (tag_abc_el)
      if (nload_dg_el.gt.0) deallocate (tag_dg_el, tag_dg_yn)
      if (nfunc.gt.0) deallocate (func_type,func_indx,func_data)
      if (n_case.gt.0) deallocate(zs_elev, zs_all, sub_tag_all,vs_tria, thick) !tag_case, val_case, tol_case
      if (nface.gt.0) deallocate(con_spx_bc_loc,con_bc)
      if (nload_sism_el.gt.0) deallocate(num_node_sism, factor_seismic_moment, tau_seismic_moment)
      if (nload_sism_el.gt.0) deallocate(check_node_sism, check_dist_node_sism, val_sism_el, fun_sism_el, tag_sism_el)
      if (nload_expl_el.gt.0) deallocate(num_node_expl, factor_explosive_source)
      if (nload_expl_el.gt.0) deallocate(check_node_expl, check_dist_node_expl, val_expl_el, fun_expl_el, tag_expl_el)
      if (nelem_abc.gt.0) deallocate(ielem_abc)

      if (nelem_dg .gt.0) deallocate(el_new)
      if (nelem_dg_glo .gt. 0) deallocate(node_send_jump, node_recv_jump, proc_send_jump, proc_recv_jump, con_spx_dg)
      if (trof .eq. 0 .and. nelem_dg_glo .gt. 0) deallocate(local_node_num_dg)
      if (damping_type .eq. 2) deallocate(frequency_range,Y_lambda,Y_mu)
      if (damping_type .eq. 3) deallocate(A0_ray, A1_ray)


      
      
      
      
      
      
      
      
      
      end subroutine DEALLOCATE_VARIABLES
