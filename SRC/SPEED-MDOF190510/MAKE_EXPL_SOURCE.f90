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

!> @brief Computes the explosive moment tensor.
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] nn number of 1-D GLL nodes
!> @param[in] xq GLL nodes
!> @param[in] wq GLL weights  
!> @param[in] dd spectral derivatives matrix
!> @param[in] nl_expl  number of explosive loads
!> @param[in] length_cne length of check_ns 
!> @param[in] check_ne   info about sismic load 
!> @param[in] check_dist_ne distance from epicenter
!> @param[in] ielem  element index
!> @param[in] facsexpl seismic moment factor 
!> @param[in] number of fuctions
!> @param[in] func_type function types
!> @param[in] func_indx function indixes for function data
!> @param[in] func_data  function data
!> @param[in] t_stress  time 
!> @param[in] cs_nnz length of cs
!> @param[in] cs spectral connectivity vector
!> @param[in] tag_func functin label
!> @param[in] nelem_loc number of local elements 
!> @param[in] local_el_num local element numeration
!> @param[in] nn_loc number of local nodes
!> @param[in] local_n_num local node numeration
!> @param[in,out] sxx nodal values for the stress tensor
!> @param[in,out] syy nodal values for the stress tensor
!> @param[in,out] szz nodal values for the stress tensor
!> @param[in,out] syz nodal values for the stress tensor
!> @param[in,out] szx nodal values for the stress tensor
!> @param[in,out] sxy nodal values for the stress tensor


      subroutine MAKE_EXPL_SOURCE(nn,ct,ww,dd,&                                                
                               sxx,syy,szz,syz,szx,sxy,&                                        
                               nl_expl,&                                                                        
                               check_ne,check_dist_ne,&                                                
                               length_cne,ielem,facsexpl,&                                        
                               func_type,func_indx,func_data,nfdata,nf,t_stress,&        
                               cs_nnz,cs,tag_func, &
                               nelem_loc, local_el_num,&
                               nn_loc, local_n_num)        
          
     implicit none
      
      integer*4 :: nn, nelem_loc, nn_loc
      integer*4 :: i,j,k,p,q,r
      integer*4 :: length_cne,ielem,nl_expl,nfdata                                
      integer*4 :: nf                                                                        
      integer*4 :: is,in,fn, ie                                                                
      integer*4 :: cs_nnz                                                                

      integer*4, dimension(nelem_loc) :: local_el_num
      integer*4, dimension(nn_loc) :: local_n_num
      integer*4, dimension(nf) :: tag_func                                
      integer*4, dimension(nf) :: func_type                                
      integer*4, dimension(nf+1) :: func_indx                        

      integer*4, dimension(0:cs_nnz) :: cs                                
      integer*4, dimension(length_cne,4) :: check_ne        

      real*8 :: get_func_value_sism
      real*8 :: t_stress                                                                

      real*8, dimension(nfdata) :: func_data                                        
      real*8, dimension(nn) :: ct,ww

      real*8, dimension(nn,nn) :: dd
      real*8, dimension(length_cne,1) :: check_dist_ne        
      real*8, dimension(nl_expl,6) :: facsexpl                        

      real*8, dimension(nn,nn,nn) :: sxx,syy,szz,syz,szx,sxy

!     STRESS CALCULATION
      
          
                                                                                  
         if ((ielem.ge.check_ne(1,4)).and. (ielem.le.check_ne(length_cne,4))) then                                                
            
            do fn = 1,nf                                                                        

            do i = 1,length_cne                                                                        
               if (ielem.eq.check_ne(i,4)) then                                                        
               

                  do r = 1,nn                                                                        
                    do q = 1,nn                                                                        
                      do p = 1,nn                                                                        
                                                                         

                         in = local_n_num(cs(cs(ielem -1) +is))                                                
                            
                         if ((tag_func(fn).eq.(check_ne(i,2))) .and.(check_ne(i,1).eq.in)) then                                        
                              

                         !STRESS RE-CALCULATION                                                        


                            sxx(p,q,r) = sxx(p,q,r) &
                                         - get_func_value_sism(nf,func_type,func_indx,func_data, nfdata,&        
                                         fn,t_stress,check_dist_ne(i,1),0.d0) &                        
                                         * facsexpl(check_ne(i,3),1)                                
                                         
                            syy(p,q,r) = syy(p,q,r) &
                                         - get_func_value_sism(nf,func_type,func_indx,func_data,  nfdata,&        
                                         fn,t_stress,check_dist_ne(i,1),0.d0) &                        
                                         * facsexpl(check_ne(i,3),2)                                
                                         
                            szz(p,q,r) = szz(p,q,r) &
                                         - get_func_value_sism(nf,func_type,func_indx,func_data,  nfdata,&        
                                         fn,t_stress,check_dist_ne(i,1),0.d0) &                        
                                         * facsexpl(check_ne(i,3),3)                                
                                         
                            syz(p,q,r) = syz(p,q,r) &
                                         - get_func_value_sism(nf,func_type,func_indx,func_data,  nfdata,&        
                                         fn,t_stress,check_dist_ne(i,1),0.d0) &                        
                                         * facsexpl(check_ne(i,3),4)                                
                                         
                            szx(p,q,r) = szx(p,q,r) &
                                         - get_func_value_sism(nf,func_type,func_indx,func_data,  nfdata,&        
                                         fn,t_stress,check_dist_ne(i,1),0.d0) &                        
                                         * facsexpl(check_ne(i,3),5)                                
                                         
                            sxy(p,q,r) = sxy(p,q,r) &
                                         - get_func_value_sism(nf,func_type,func_indx,func_data,  nfdata,&        
                                         fn,t_stress,check_dist_ne(i,1),0.d0) &                        
                                         * facsexpl(check_ne(i,3),6)                                                                

                          endif ! if ((tag_func(fn).eq.(check_node_expl(i,2))) &...                

                      enddo                                                        
                   enddo                                                        
                enddo                                                                

             endif !ielem                                                                
          enddo !i                                                                        
       enddo !fn                                                                        

    endif  !ielem        
                                                                                                                                                                 
      
   return
      
   end subroutine MAKE_EXPL_SOURCE

