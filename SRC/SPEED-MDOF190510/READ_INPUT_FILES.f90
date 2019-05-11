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


!> @brief Reads input files and allocates memory.
!! @author Ilario Mazzieri
!> @date November, 2014 
!> @version 1.0


     subroutine READ_INPUT_FILES()
     
      use max_var
      use speed_par
      

      implicit none
 
      include 'SPEED.MPI'

      head_file = 'SPEED.input'
      
      if (mpi_id.eq.0) then
         write(*,'(A)') 
         write(*,'(A)')'------------------Reading Header File------------------'
         write(*,'(A,A36)') 'Header File: ',head_file
      endif
      
      inquire(file=head_file,exist=filefound);
      if(filefound .eqv. .FALSE.) call EXIT(EXIT_MISSING_FILE)
       
      opt_out_var= 1;  opt_out_form = 1; opt_out_data = 1
      

      call READ_DIME_HEADER(head_file,nsnaps,ntime_err,debug) 

      if (nsnaps.gt.0) allocate(tsnap(nsnaps),itersnap(nsnaps))
      if(ntime_err .ne. 0) allocate(time_error(ntime_err))
            
      tstart = 0.0
      tstop = 0.0
      
      ! Set default variables
      b_setuponly = b_setuponly_default
      b_failoncoeffs = b_failoncoeffs_default
      b_instabilitycontrol = b_instabilitycontrol_default
      b_failCFL = b_failCFL_default
      instability_maxval = instability_maxval_default
      dg_c = -100.d0;
      pen_c = -100.d0;

      call READ_HEADER(head_file,grid_file,mat_file, mpi_file, monitor_file, bkp_file, &
                       deltat,tstart,tstop,&
                       opt_out_var,& 
                       opt_out_data,opt_out_form,&
                       nsnaps,tsnap,&
                       ndt_mon_lst,&        
                       deltat_fixed,&        
                       depth_search_mon_pgm,ndt_mon_pgm,num_pgm,&
                       rotation_angle_mon_pgm,file_mon_pgm,&        
                       depth_search_mon_lst,num_lst,&        
                       file_mon_lst,dg_c,pen_c,&
                       rk_scheme, rk_order, rk_stages, testmode, ntime_err, time_error,&
                       damping_type, &
                       b_failoncoeffs, b_setuponly, b_failCFL, b_instabilitycontrol, instability_maxval)

      if (mpi_id.eq.0) write(*,'(A)')'Read.'

      ! Instability control
      if (mpi_id .eq. 0) then
          write(*,*)
          write(*,'(A,L)') 'Instability control  :   ', b_instabilitycontrol
          if (b_instabilitycontrol) then
            write(*,'(A,E12.4,A,E12.4,A)') 'Instability threshold: ', instability_maxval,' (default: ', instability_maxval_default, ')'
          endif
          write(*,*)
      endif
      
      if (mpi_id.eq.0) then
         if(dg_c .eq. 1.d0)   write(*,'(A)') 'DG METHOD    : NIPG'
         if(dg_c .eq. 0.d0)   write(*,'(A)') 'DG METHOD    : IIPG'
         if(dg_c .eq. -1.d0)  write(*,'(A)') 'DG METHOD    : SIPG'
         if(pen_c .ge. 0.d0)  write(*,'(A,E12.4)') 'PENALIZATION : ', pen_c
         if(dg_c .eq. -100)   write(*,'(A)') 'SEM METHOD    '
      endif


      mat_file = mat_file(1:len_trim(mat_file)) // '.mate'
      
      if (mpi_id.eq.0) then
         write(*,'(A)')    
         write(*,'(A,A20)')'-----------------Reading Material File-----------------'
         write(*,'(A,A20)')'Material File    : ',mat_file
      endif
      
      inquire(file=mat_file,exist=filefound);
      if(filefound .eqv. .FALSE.) call EXIT(EXIT_MISSING_FILE)
      
      call READ_DIME_FILEMATE(mat_file,nmat, &
                            nmat_nle, &                                                                        
                            nload_dirX_el,nload_dirY_el,nload_dirZ_el, &
                            nload_neuX_el,nload_neuY_el,nload_neuZ_el, &
                            nload_neuN_el, &                                                            
                            nload_poiX_el,nload_poiY_el,nload_poiZ_el, &
                            nload_traX_el,nload_traY_el,nload_traZ_el, &
                            nload_plaX_el,nload_plaY_el,nload_plaZ_el, &                                
                            nload_forX_el,nload_forY_el,nload_forZ_el, &
                            nload_forc_el,nload_pres_el,nload_shea_el, &
                            nload_abc_el,nload_dg_el,nfunc,nfunc_data,&
                            nload_sism_el,&                                                                 
                            nload_expl_el,&                                                                
                            n_case,n_test)                   
                                                                                                                            
      if(n_test.gt.0 .and. mpi_id .eq. 0)  write(*,'(A)')'*********TEST MODE*********'

      if (mpi_id.eq.0) then
         write(*,'(A,I8)') 'Materials        : ',nmat
         if(nmat_nle.gt.0) write(*,'(A,I8)')     'Materials NL-El. : ',nmat_nle        
         if(nload_dirX_el.gt.0) write(*,'(A,I8)')'Dirichlet X B.C. : ',nload_dirX_el
         if(nload_dirY_el.gt.0) write(*,'(A,I8)')'Dirichlet Y B.C. : ',nload_dirY_el
         if(nload_dirZ_el.gt.0) write(*,'(A,I8)')'Dirichlet Z B.C. : ',nload_dirZ_el
         if(nload_neuX_el.gt.0) write(*,'(A,I8)')'Neumann X B.C.   : ',nload_neuX_el
         if(nload_neuY_el.gt.0) write(*,'(A,I8)')'Neumann Y B.C.   : ',nload_neuY_el
         if(nload_neuZ_el.gt.0) write(*,'(A,I8)')'Neumann Z B.C.   : ',nload_neuZ_el
         if(nload_neuN_el.gt.0) write(*,'(A,I8)')'Neumann N B.C.   : ',nload_neuN_el              
         if(nload_poiX_el.gt.0) write(*,'(A,I8)')'Point Loads X    : ',nload_poiX_el
         if(nload_poiY_el.gt.0) write(*,'(A,I8)')'Point Loads Y    : ',nload_poiY_el
         if(nload_poiZ_el.gt.0) write(*,'(A,I8)')'Point Loads Z    : ',nload_poiZ_el
         if(nload_traX_el.gt.0) write(*,'(A,I8)')'Travel. Load X   : ',nload_traX_el
         if(nload_traY_el.gt.0) write(*,'(A,I8)')'Travel. Load Y   : ',nload_traY_el
         if(nload_traZ_el.gt.0) write(*,'(A,I8)')'Travel. Load Z   : ',nload_traZ_el         
         if(nload_plaX_el.gt.0) write(*,'(A,I8)')'Plane Loads X    : ',nload_plaX_el              
         if(nload_plaY_el.gt.0) write(*,'(A,I8)')'Plane Loads Y    : ',nload_plaY_el              
         if(nload_plaZ_el.gt.0) write(*,'(A,I8)')'Plane Loads Z    : ',nload_plaZ_el              
         if(nload_forX_el.gt.0) write(*,'(A,I8)')'Force X          : ',nload_forX_el
         if(nload_forY_el.gt.0) write(*,'(A,I8)')'Force Y          : ',nload_forY_el
         if(nload_forZ_el.gt.0) write(*,'(A,I8)')'Force Z          : ',nload_forZ_el
         if(nload_forc_el.gt.0) write(*,'(A,I8)')'Force load       : ',nload_forc_el
         if(nload_pres_el.gt.0) write(*,'(A,I8)')'Pressure Load    : ',nload_pres_el
         if(nload_shea_el.gt.0) write(*,'(A,I8)')'Shear Load       : ',nload_shea_el
         if(nload_abc_el.gt.0)  write(*,'(A,I8)')'ABSO Boundaries  : ',nload_abc_el
         if(nload_dg_el.gt.0)   write(*,'(A,I8)')'DG Interfaces    : ',nload_dg_el
         if(nfunc.gt.0)         write(*,'(A,I8)')'Functions        : ',nfunc
         if(nload_sism_el.gt.0) write(*,'(A,I8)')'Seis. Mom. Load  : ',nload_sism_el        
         if(nload_expl_el.gt.0) write(*,'(A,I8)')'Explosive Load   : ',nload_expl_el            
         if (n_case.gt.1) then
            write(*,'(A)')'CASE WARNING: More than one case defined,'        
            write(*,'(A)')'              only the 1st case will be adopted'  
         endif
         write(*,'(A,I8)') 'CASE             : ',n_case                               
      endif
      
      if (nmat.le.0) then
         write(*,*)'Error ! nmat = 0'
         call EXIT(EXIT_NO_MATERIALS)
      endif
      
      allocate (type_mat(nmat), tref_mat(nmat), prop_mat(nmat,4), sdeg_mat(nmat), tag_mat(nmat))

      if (nmat_nle.gt.0) allocate(type_mat_nle(nmat_nle), prop_mat_nle(nmat_nle,1), val_mat_nle(nmat_nle,1), tag_mat_nle(nmat_nle))
      
      if (nload_dirX_el.gt.0) allocate (val_dirX_el(nload_dirX_el,4), fun_dirX_el(nload_dirX_el), tag_dirX_el(nload_dirX_el))      
      if (nload_dirY_el.gt.0) allocate (val_dirY_el(nload_dirY_el,4), fun_dirY_el(nload_dirY_el), tag_dirY_el(nload_dirY_el))      
      if (nload_dirZ_el.gt.0) allocate (val_dirZ_el(nload_dirZ_el,4), fun_dirZ_el(nload_dirZ_el), tag_dirZ_el(nload_dirZ_el))
      !!!!!!!!!!!!!!!modified by ty 170410
      if (nload_neuX_el.gt.0) allocate (val_neuX_el(nload_neuX_el,6), fun_neuX_el(nload_neuX_el), tag_neuX_el(nload_neuX_el))      
      if (nload_neuY_el.gt.0) allocate (val_neuY_el(nload_neuY_el,6), fun_neuY_el(nload_neuY_el), tag_neuY_el(nload_neuY_el))      
      if (nload_neuZ_el.gt.0) allocate (val_neuZ_el(nload_neuZ_el,6), fun_neuZ_el(nload_neuZ_el), tag_neuZ_el(nload_neuZ_el))
      if (nload_neuN_el.gt.0) allocate (val_neuN_el(nload_neuN_el,6), fun_neuN_el(nload_neuN_el), tag_neuN_el(nload_neuN_el))
      !!!!!!!!!!!!!!!modified by ty 170410
      if (nload_poiX_el.gt.0) allocate (val_poiX_el(nload_poiX_el,6), fun_poiX_el(nload_poiX_el))    
      if (nload_poiY_el.gt.0) allocate (val_poiY_el(nload_poiY_el,6), fun_poiY_el(nload_poiY_el))   
      if (nload_poiZ_el.gt.0) allocate (val_poiZ_el(nload_poiZ_el,6), fun_poiZ_el(nload_poiZ_el))
        !!!!!!!!!!!!!!!modified by ty 170410
      if (nload_traX_el.gt.0) allocate (val_traX_el(nload_traX_el,4), fun_traX_el(nload_traX_el))    
      if (nload_traY_el.gt.0) allocate (val_traY_el(nload_traY_el,4), fun_traY_el(nload_traY_el))   
      if (nload_traZ_el.gt.0) allocate (val_traZ_el(nload_traZ_el,4), fun_traZ_el(nload_traZ_el))


      if (nload_plaX_el.gt.0) allocate (val_plaX_el(nload_plaX_el,1), fun_plaX_el(nload_plaX_el), tag_plaX_el(nload_plaX_el))
      if (nload_plaY_el.gt.0) allocate (val_plaY_el(nload_plaY_el,1), fun_plaY_el(nload_plaY_el), tag_plaY_el(nload_plaY_el))
      if (nload_plaZ_el.gt.0) allocate (val_plaZ_el(nload_plaZ_el,1), fun_plaZ_el(nload_plaZ_el), tag_plaZ_el(nload_plaZ_el))
      
      if (nload_forX_el.gt.0) allocate (val_forX_el(nload_forX_el,4), fun_forX_el(nload_forX_el))
      if (nload_forY_el.gt.0) allocate (val_forX_el(nload_forY_el,4), fun_forX_el(nload_forY_el))
      if (nload_forZ_el.gt.0) allocate (val_forX_el(nload_forZ_el,4), fun_forX_el(nload_forZ_el))
     
     
      if (n_test.gt.0) allocate (fun_test(n_test)) !val_fun_test(n_test,10),      
                                       
      if (nload_forc_el.gt.0) allocate (val_forc_el(nload_forc_el,10), fun_forc_el(nload_forc_el))     
      if (nload_pres_el.gt.0) allocate (val_pres_el(nload_pres_el,10), fun_pres_el(nload_pres_el))      
      if (nload_shea_el.gt.0) allocate (val_shea_el(nload_shea_el,10), fun_shea_el(nload_shea_el))
      
      if (nload_abc_el.gt.0) allocate (tag_abc_el(nload_abc_el))
      if (nload_dg_el.gt.0) allocate (tag_dg_el(nload_dg_el), tag_dg_yn(nload_dg_el))


      if (nload_sism_el.gt.0) allocate (val_sism_el(nload_sism_el,21), &
                                        fun_sism_el(nload_sism_el), tag_sism_el(nload_sism_el))     
      if (nload_expl_el.gt.0) allocate (val_expl_el(nload_expl_el,20), &
                                        fun_expl_el(nload_expl_el), tag_expl_el(nload_expl_el))
      
      !if (n_case.gt.0) allocate (val_case(n_case,1), tag_case(n_case), tol_case(n_case))
      
      if (nfunc.gt.0) allocate (tag_func(nfunc), func_type(nfunc), func_indx(nfunc +1), func_data(nfunc_data))
      
      
      allocate(QS(nmat), QP(nmat)); QS = 0.d0; QP = 0.d0;
            
      call READ_FILEMATE(mat_file,nmat,prop_mat,type_mat,tref_mat,tag_mat, QS, QP,&
                nmat_nle,prop_mat_nle,val_mat_nle,type_mat_nle,tag_mat_nle, &        
                nload_dirX_el,val_dirX_el,fun_dirX_el,tag_dirX_el, &
                nload_dirY_el,val_dirY_el,fun_dirY_el,tag_dirY_el, &
                nload_dirZ_el,val_dirZ_el,fun_dirZ_el,tag_dirZ_el, &
                nload_neuX_el,val_neuX_el,fun_neuX_el,tag_neuX_el, &
                nload_neuY_el,val_neuY_el,fun_neuY_el,tag_neuY_el, &
                nload_neuZ_el,val_neuZ_el,fun_neuZ_el,tag_neuZ_el, &
                nload_neuN_el,val_neuN_el,fun_neuN_el,tag_neuN_el, &            
                nload_poiX_el,val_poiX_el,fun_poiX_el,&
                nload_poiY_el,val_poiY_el,fun_poiY_el,&
                nload_poiZ_el,val_poiZ_el,fun_poiZ_el,&
                nload_traX_el,val_traX_el,fun_traX_el,&
                nload_traY_el,val_traY_el,fun_traY_el,&
                nload_traZ_el,val_traZ_el,fun_traZ_el,&               
                nload_plaX_el,val_plaX_el,fun_plaX_el,tag_plaX_el, &            
                nload_plaY_el,val_plaY_el,fun_plaY_el,tag_plaY_el, &                 
                nload_plaZ_el,val_plaZ_el,fun_plaZ_el,tag_plaZ_el, &                 
                nload_forX_el,val_forX_el,fun_forX_el,&
                nload_forY_el,val_forY_el,fun_forY_el,&
                nload_forZ_el,val_forZ_el,fun_forZ_el,&
                nload_forc_el,val_forc_el,fun_forc_el,&
                nload_pres_el,val_pres_el,fun_pres_el,&
                nload_shea_el,val_shea_el,fun_shea_el,&
                n_test,fun_test,& !val_fun_test,&
                nload_abc_el,tag_abc_el,&
                nload_dg_el,tag_dg_el,tag_dg_yn, &
                nload_sism_el,val_sism_el,fun_sism_el,tag_sism_el, &                 
                nload_expl_el,val_expl_el,fun_expl_el,tag_expl_el, &                 
                n_case,val_case,tag_case,tol_case, &                                 
                nfunc,func_type,func_indx,func_data,tag_func,nfunc_data, &
                fmax,fpeak)
      !  write(*,*)nload_dg_el,tag_dg_el(:) ,tag_dg_yn(:)         

      if(fmax .eq. 0.d0) then
         fmax = 3.d0 
         if (mpi_id .eq. 0) write(*,*) 'ATTENTION: FMAX not defined!'
         if (mpi_id .eq. 0) write(*,*) 'FMAX ASSUMED = 3'
      endif   




      if (damping_type .eq. 1) then 

         do im = 1, nmat
            if(QS(im) .eq. 0.d0) then 
               prop_mat(im,4) = 0.d0;
            else
               prop_mat(im,4) = 4.d0*datan(1.d0)*(fmax/3.d0)/QS(im)
            endif   

         enddo
         
         
       else 
       
            prop_mat(:,4) = 0.d0;

       endif   
      
      
      
!**************************************************
! Materials can have different spectral degree

      do i = 1,nmat
         sdeg_mat(i) = type_mat(i)
      enddo

!************************************************** 
           
      if (mpi_id.eq.0) write(*,'(A)')'Read.'

      if (mpi_id.eq.0) then                        
         write(*,'(A)')                                
         if (n_case.gt.0) then
                               
              write(*,'(A)') '------------------Not-Honoring case-------------------'
              write(*,'(A,I8)') 'CASE :', tag_case    !tag_case(1)
            select case (tag_case)  !(tag_case(1))
              case(1)
                write(*,'(A)')'GRENOBLE H.'                        
              case(2)
                write(*,'(A)')'GRENOBLE'                        
              case(3)
                write(*,'(A)')'GUBBIO'                        
              case(4)
                write(*,'(A)')'SULMONA'                        
              case(5)
                write(*,'(A)')'VOLVI'                        
              case(6)
                write(*,'(A)')'FRIULI'                        
              case(7)
                write(*,'(A)')'AQUILA'   
              case(8)
                write(*,'(A)')'SANTIAGO'   
              case(11)
                write(*,'(A)')'CHRISTCHURCH NEW TOPO'        
              case(12)
                write(*,'(A)')'PO PLAIN (new model)'    
              case(13)
                write(*,'(A)')'PO PLAIN-BEDROCK (new-model)'                                 
              case(14)
                write(*,'(A)')'WELLINGTON (Benites)'                                 
              case(15)
                write(*,'(A)')'MARSICA-FUCINO'   
              case(16)
                write(*,'(A)')'ISTANBUL'                     
              case(20)
                write(*,'(A)')'ATENE'                     
              case(98,99)
                write(*,'(A)')'TEST MODE'                                
              case default
                 write(*,'(A)')'ERROR.. this case was not implemented!'
             end select 
                                   
               write(*,'(A,I8)')    'MATERIAL TYPE    : ',val_case   !(val_case(1,1))        
              write(*,'(A,E12.4)') 'MATERIAL TOL.    : ',tol_case   !(tol_case(1))
            write(*,*)
         endif
      endif
      
      if (mpi_id.eq.0) then
         write(*,'(A)')
         do im = 1,nmat    
            write(*,'(A,I8)')    'MATERIAL : ',tag_mat(im)
            !write(*,'(A,3E12.4)') 'PROPERTIES : ',(prop_mat(im,i),i=1,3)
            write(*,'(A,I8)')    'DEGREE   : ',type_mat(im)
            write(*,'(A,E12.4)') 'rho      : ',prop_mat(im,1)
            write(*,'(A,E12.4)') 'Vp       : ',((prop_mat(im,2) + 2*prop_mat(im,3))/prop_mat(im,1))**0.5
            write(*,'(A,E12.4)') 'Vs       : ',(prop_mat(im,3)/prop_mat(im,1))**0.5
            !write(*,'(A,E12.4)') 'gamma    : ',prop_mat(im,4)
            write(*,'(A,E12.4)') 'Qs       : ',QS(im)
            write(*,'(A,E12.4)') 'Qp       : ',QP(im)
            write(*,*)
         enddo
      endif

      if (mpi_id.eq.0) then                                                                                
          write(*,'(A)')                                                                                
         do im_nle = 1,nmat_nle                                                                                
            write(*,'(A,I8)') 'NON LINEAR ELASTIC - MATERIAL : ',tag_mat_nle(im_nle)                        
            !write(*,'(A,I8)') 'PROPERTIES : ',(prop_mat_nle(im_nle,1))                                        
            write(*,'(A,E12.4)') 'NLE DEPTH : ',(val_mat_nle(im_nle,1))                                        
            write(*,'(A,I8)')    'DEGREE    : ',type_mat_nle(im_nle)                                                
            write(*,*)                                                                                        
         enddo                                                                                                
         if ((fpeak.eq.0).and.(nmat_nle.gt.0)) then                                                        
             write(*,'(A)')'ERROR: Peak frequency for damping not defined!'                                
             call EXIT(EXIT_DAMPING_PEAK)
         endif                                                                                                   
       endif                                                                                                
      

      grid_file = grid_file(1:len_trim(grid_file)) // '.mesh'    
      if (mpi_id.eq.0)  write(*,'(A)') '-------------------Reading Grid File-------------------'
      if (mpi_id.eq.0)  write(*,'(A,A20)') 'Grid File : ',grid_file
      
      inquire(file=grid_file,exist=filefound);
      if(filefound .eqv. .FALSE.) call EXIT(EXIT_MISSING_FILE)
      

!counting hexahedras and squares
       call READ_DIME_FILEMESH(grid_file,nmat,tag_mat,&
                     nload_dirX_el,tag_dirX_el,nload_dirY_el,tag_dirY_el,&
                     nload_dirZ_el,tag_dirZ_el,&
                     nload_neuX_el,tag_neuX_el,nload_neuY_el,tag_neuY_el,&
                     nload_neuZ_el,tag_neuZ_el,&
                     nload_neuN_el,tag_neuN_el,& 
                     nload_abc_el,tag_abc_el,&
                     nload_dg_el,tag_dg_el,&
                     nnod_macro,nelem,nface)
      
      if (mpi_id.eq.0) then
         write(*,'(A,I8)')'Nodes     : ', nnod_macro
         write(*,'(A,I8)')'Elements  : ', nelem
         write(*,'(A,I8)')'Faces     : ', nface
      endif
      
      
      if (nelem.gt.0) then
         allocate (con(nelem,9))
      else
         write(*,*)'Error ! nelem = 0'
         call EXIT(EXIT_NO_ELEMENTS)
      endif
      
      if (nface.gt.0) allocate (con_bc(nface,5))
      
            
      call READ_FILEMESH(grid_file,nmat,tag_mat,&
                        nload_dirX_el,tag_dirX_el,nload_dirY_el,tag_dirY_el,&
                        nload_dirZ_el,tag_dirZ_el,&
                        nload_neuX_el,tag_neuX_el,nload_neuY_el,tag_neuY_el,&
                        nload_neuZ_el,tag_neuZ_el,&
                        nload_neuN_el,tag_neuN_el,& 
                        nload_abc_el,tag_abc_el,&
                        nload_dg_el,tag_dg_el,&
                        nnod_macro,xx_macro, yy_macro, zz_macro,  &
                        nelem,con,nface,con_bc)

      if(mpi_id .eq. 0) write(*,'(A)') 'Read.'
      if(mpi_id .eq. 0) write(*,'(A)')
     
     end subroutine READ_INPUT_FILES
