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

!> @brief Creates the seismic moment tensor.
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] nn number of 1-D GLL nodes
!> @param[in] xq GLL nodes
!> @param[in] wq GLL weights  
!> @param[in] dd spectral derivatives matrix
!> @param[in] nl_sism  number of sismic loads
!> @param[in] length_cns length of check_ns 
!> @param[in] check_ns   info about sismic load 
!> @param[in] check_dist_ns distance from epicenter
!> @param[in] ielem  element index
!> @param[in] facsmom seismic moment factor 
!> @param[in] tausmom rise time factor
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

      subroutine MAKE_SEISMIC_MOMENT(nn,ct,ww,dd,&                                                        
                               sxx,syy,szz,syz,szx,sxy,&                                        
                               nl_sism,&                                                                        
                               check_ns,check_dist_ns,&                                                
                               length_cns,ielem,facsmom,&                                        
                               tausmom,&                
                               func_type,func_indx,func_data,nfdata,nf,t_stress,&        
                               cs_nnz,cs,tag_func, &
                               nelem_loc, local_el_num,&
                               nn_loc, local_n_num)                                                        

      implicit none
      
      integer*4 :: nn, nelem_loc, nn_loc
      integer*4 :: i,j,k,p,q,r
      integer*4 :: cs_nnz                                                                
      integer*4 :: is,in,fn,ie                                                                
      integer*4 :: length_cns,ielem,nl_sism                                
      integer*4 :: nf,nfdata                                                                        
      
      integer*4, dimension(nelem_loc) :: local_el_num
      integer*4, dimension(nn_loc) :: local_n_num
      integer*4, dimension(nf) :: func_type                                
      integer*4, dimension(nf +1) :: func_indx                        
      integer*4, dimension(nf) :: tag_func                                
      integer*4, dimension(0:cs_nnz) :: cs                                                

      integer*4, dimension(length_cns,4) :: check_ns        

      real*8 :: t_stress                                                                
      real*8 :: GET_FUNC_VALUE_SISM                                                

      real*8, dimension(nn) :: ct,ww
      real*8, dimension(nfdata) :: func_data                                        

      real*8, dimension(nn,nn) :: dd
      real*8, dimension(length_cns,1) :: check_dist_ns        
      real*8, dimension(nl_sism,6) :: facsmom                        
      real*8, dimension(nl_sism,1) :: tausmom                        

      real*8, dimension(nn,nn,nn) :: sxx,syy,szz,syz,szx,sxy
 
                                                                                        
      if ((ielem .ge. check_ns(1,4)) .and. (ielem .le. check_ns(length_cns,4))) then        
      
            do i = 1,length_cns                                                                        

               if (ielem .eq. check_ns(i,4)) then                                                        
               

                  do r = 1,nn                                                                        
                    do q = 1,nn                                                                        
                      do p = 1,nn                                                                        
                         
                         is = nn*nn*(r -1) +nn*(q -1) +p

                                                 
                          in = local_n_num(cs(cs(ielem -1) +is))
                                                 
                                 do fn = 1,nf        
                            
                                    if ((tag_func(fn) .eq. (check_ns(i,2))) .and. (check_ns(i,1) .eq. in)) then 
                              

                                         !STRESS RE-CALCULATION                                                        

                                    !t_stress = current time, 
                                    !check_dist_ns(i,1) = rupture time
                                    !tausmom = rise time
                                    !facsmom = seismic moment
                                    sxx(p,q,r) = sxx(p,q,r) - GET_FUNC_VALUE_SISM(nf,func_type,func_indx,func_data, nfdata, &
                                                                     fn,t_stress,check_dist_ns(i,1), &
                                                                     tausmom(check_ns(i,3),1)) * facsmom(check_ns(i,3),1)
                                                                               
                                    syy(p,q,r) = syy(p,q,r) - GET_FUNC_VALUE_SISM(nf,func_type,func_indx,func_data, nfdata,&
                                                                     fn,t_stress,check_dist_ns(i,1), &
                                                                     tausmom(check_ns(i,3),1)) * facsmom(check_ns(i,3),2)
                                                                               
                                    szz(p,q,r) = szz(p,q,r) - GET_FUNC_VALUE_SISM(nf,func_type,func_indx,func_data, nfdata,&
                                                                     fn,t_stress,check_dist_ns(i,1), &
                                                                     tausmom(check_ns(i,3),1)) * facsmom(check_ns(i,3),3)
                                                                                 
                                    syz(p,q,r) = syz(p,q,r) - GET_FUNC_VALUE_SISM(nf,func_type,func_indx,func_data, nfdata,&
                                                                     fn,t_stress,check_dist_ns(i,1), &
                                                                     tausmom(check_ns(i,3),1)) * facsmom(check_ns(i,3),4)
                                                                                
                                    szx(p,q,r) = szx(p,q,r) - GET_FUNC_VALUE_SISM(nf,func_type,func_indx,func_data, nfdata,&
                                                                     fn,t_stress,check_dist_ns(i,1), &
                                                                     tausmom(check_ns(i,3),1)) * facsmom(check_ns(i,3),5)
                                                                               
                                    sxy(p,q,r) = sxy(p,q,r) - GET_FUNC_VALUE_SISM(nf,func_type,func_indx,func_data, nfdata,&
                                                                     fn,t_stress,check_dist_ns(i,1), &
                                                                     tausmom(check_ns(i,3),1)) * facsmom(check_ns(i,3),6)  
                                                                                         

                                    endif ! if ((tag_func(fn).eq.(check_node_sism(i,2))) &...        
                                                  
                          enddo !fn = 1,nf        

                     enddo !p                                                        
                  enddo !q                                                        
                enddo !r

            endif !ielem
                                                                                                         
         enddo !i                                                                        
                                                                          

         endif  !ielem        
                               
      return
      
      end subroutine MAKE_SEISMIC_MOMENT

