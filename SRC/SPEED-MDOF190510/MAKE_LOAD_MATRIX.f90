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

!> @brief Makes External Forces.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0

      subroutine MAKE_LOAD_MATRIX()
      
      use max_var
      use speed_par
      

      implicit none
 
      include 'SPEED.MPI'      


     
      if (nfunc.le.0) nfunc = 1
      !!!!!!!!!!!!!!!!!!modified by ty 170410
      allocate (Fel(nfunc,3*nnod_loc,3))
      !!!!!!!!!!!!!!!!!!modified by ty 170410
      call MAKE_EXTINT_FORCES(nnod_loc,xx_spx_loc,yy_spx_loc,zz_spx_loc,local_node_num,con_nnz_loc,con_spx_loc,&
                    nmat,tag_mat,type_mat,sdeg_mat,tref_mat,prop_mat,&
                    nelem_loc,local_el_num, &
                    alfa11,alfa12,alfa13,alfa21,alfa22,alfa23,&
                    alfa31,alfa32,alfa33,beta11,beta12,beta13,&
                    beta21,beta22,beta23,beta31,beta32,beta33,&
                    gamma1,gamma2,gamma3,delta1,delta2,delta3,&
                    con_nnz_bc_loc,con_spx_bc_loc,&
                    nload_dirX_el,val_dirX_el,fun_dirX_el,tag_dirX_el,&
                    nload_dirY_el,val_dirY_el,fun_dirY_el,tag_dirY_el,&
                    nload_dirZ_el,val_dirZ_el,fun_dirZ_el,tag_dirZ_el,&
                    nload_neuX_el,val_neuX_el,fun_neuX_el,tag_neuX_el,&
                    nload_neuY_el,val_neuY_el,fun_neuY_el,tag_neuY_el,&
                    nload_neuZ_el,val_neuZ_el,fun_neuZ_el,tag_neuZ_el,&
                    nload_neuN_el,val_neuN_el,fun_neuN_el,tag_neuN_el,&
                    nload_poiX_el,val_poiX_el,fun_poiX_el,&
                    nload_poiY_el,val_poiY_el,fun_poiY_el,&
                    nload_poiZ_el,val_poiZ_el,fun_poiZ_el,&
                    nload_traX_el,val_traX_el,fun_traX_el,&
                    nload_traY_el,val_traY_el,fun_traY_el,&
                    nload_traZ_el,val_traZ_el,fun_traZ_el,&
                    nload_plaX_el,val_plaX_el,fun_plaX_el,tag_plaX_el,&
                    nload_plaY_el,val_plaY_el,fun_plaY_el,tag_plaY_el,&
                    nload_plaZ_el,val_plaZ_el,fun_plaZ_el,tag_plaZ_el,&
                    nload_sism_el,val_sism_el,fun_sism_el,tag_sism_el,&
                    nload_expl_el,val_expl_el,fun_expl_el,tag_expl_el,&
                    nload_forX_el,val_forX_el,fun_forX_el,&
                    nload_forY_el,val_forY_el,fun_forY_el,&
                    nload_forZ_el,val_forZ_el,fun_forZ_el,&
                    nload_forc_el,val_forc_el,fun_forc_el,&
                    nload_pres_el,val_pres_el,fun_pres_el,&
                    nload_shea_el,val_shea_el,fun_shea_el,&
                    n_test,fun_test,  & !val_fun_test, &
                    nfunc,tag_func,func_type,func_indx,func_data,nfunc_data,&
                    Fel,&  
                    con,nelem,& 
                    length_check_node_sism,&                                 
                    sour_node_sism,max_num_node_sism,num_node_sism,&         
                    factor_seismic_moment,node_index_seq,&                  
                    tau_seismic_moment,&                                
                    length_check_node_expl,&                                 
                    sour_node_expl,max_num_node_expl,num_node_expl,&         
                    factor_explosive_source,&                                  
                    mpi_comm, mpi_np, mpi_id, testmode)

          
          

          
          if (mpi_id.eq.0) write(*,'(A)')'Load matrix built.'
          if (mpi_id.eq.0) write(*,'(A)')
          deallocate(con) 
          
          


!**********************************************************************************************************
!                                  SEISMIC MOMENT                           
!**********************************************************************************************************
!  check_node_sism, check_dist_node_sism are local vector
!  check_node_sism(i,1) = global node index
!  check_node_sism(i,2) = time function number
!  check_node_sism(i,3) = fault number
!  check_node_sism(i,4) = global element index
!  check_dist_node_sism(i,1) = distance from hypocenter

      if (nload_sism_el.gt.0) then        
         if (mpi_id.eq.0) write(*,'(A)')'------------------Make Seismic Moment------------------'
         allocate (check_node_sism(length_check_node_sism,4))
         allocate (check_dist_node_sism(length_check_node_sism,1))
      
         call CHECK_SISM(con_nnz_loc, con_spx_loc,&                        
                         nmat, tag_mat, sdeg_mat, &                          
                         nload_sism_el,&
                         num_node_sism, max_num_node_sism,&        
                         sour_node_sism, dist_sour_node_sism, &        
                         check_node_sism, check_dist_node_sism, &        
                         length_check_node_sism,&                
                         fun_sism_el, nfunc, tag_func, val_sism_el, &
                         nnod_loc, local_node_num)
        
         deallocate (sour_node_sism)                                
         deallocate (dist_sour_node_sism)                        

      endif                                                        

!**********************************************************************************************************
!                                  EXPLOSIVE SOURCE                           
!**********************************************************************************************************
!  check_node_expl, check_dist_node_expl are local vector
!  check_node_expl(i,1) = global node index
!  check_node_expl(i,2) = time function number
!  check_node_expl(i,3) = fault number
!  check_node_expl(i,4) = global element index
!  check_dist_node_expl(i,1) = distance from hypocenter
        
      if (nload_expl_el.gt.0) then      
         if (mpi_id.eq.0) write(*,'(A)')'-----------------Make explosive source-----------------'                  
         allocate (check_node_expl(length_check_node_expl,4))                
         allocate (check_dist_node_expl(length_check_node_expl,1))        
      
         call CHECK_EXPL(con_nnz,con_spx,&                        
                         nmat,tag_mat,sdeg_mat,&                  
                         nload_expl_el,&                        
                         num_node_expl,max_num_node_expl,&        
                         sour_node_expl,dist_sour_node_expl,&        
                         check_node_expl,check_dist_node_expl,&        
                         length_check_node_expl,&                
                         fun_expl_el,nfunc,tag_func,val_expl_el, &
                         nnod_loc, local_node_num)
        
         deallocate (sour_node_expl)                                
         deallocate (dist_sour_node_expl)                        

      endif                                                        


      
      
      end subroutine MAKE_LOAD_MATRIX
