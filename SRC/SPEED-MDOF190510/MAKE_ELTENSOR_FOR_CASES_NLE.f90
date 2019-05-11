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


!> @brief Assignes material properties node by node for non linear elasticity. 
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] vcase value case for non linear block 
!> @param[in] R_el  reduction factor 
!> @param[in] nn  polynomial degree
!> @param[in] cs_nnz length of cs
!> @param[in] cs connectivity vector
!> @param[in] ielem  element index
!> @param[in] nf number of functions
!> @param[in] func_type function type
!> @param[in] func_indx indices for the data 
!> @param[in] func_data(*)  data for the calculation (depending on type) 
!> @param[in] t_stress  time
!> @param[in] tag_func labels for functions 
!> @param[in] yon is the node nonlinear elastic yes (1) or not (0)?
!> @param[in] nnod_loc  number of local nodes
!> @param[in] vs_tria  vs values for each node
!> @param[out] rho_el material density 
!> @param[out] lambda_el Lame coefficient lambda
!> @param[out] mu_el Lame coefficient mu
!> @param[out] gamma_el damping coefficient gamma

      subroutine MAKE_ELTENSOR_FOR_CASES_NLE(vcase,R_el,&                        
                                 nn,rho_el,lambda_el,mu_el,gamma_el,&                                        
                                 cs_nnz,cs,ielem,&                                                                        
                                 func_type,func_indx,func_data,nfdata,&
                                 nf,t_stress,tag_func,yon,tcase,&
                                 nnod_loc,vs_tria)        


      
      implicit none
                                                      
      integer*4 :: vcase                        
      integer*4 :: nn
      integer*4 :: p,q,r
      integer*4 :: is,in,ielem        
      integer*4 :: cs_nnz                                
      integer*4 :: nf,fn,nfdata,tcase,nnod_loc                                                                        
      
      integer*4, dimension(nf) :: tag_func                                
      integer*4, dimension(nf) :: func_type                                
      integer*4, dimension(nf +1) :: func_indx                        
      integer*4, dimension(0:cs_nnz) :: cs                

      integer*4, dimension(nn,nn,nn) :: yon
      
      real*8 :: VS,VP
      real*8 :: Vratio
      real*8 :: t_stress                                                                
      real*8 :: get_func_value        
      real*8 :: reduction_factor                                        
      
      real*8, dimension(nfdata) :: func_data    
      real*8, dimension(nnod_loc) :: vs_tria                                                                            

      real*8, dimension(nn,nn,nn) :: rho_el,lambda_el,mu_el,gamma_el
      real*8, dimension(nn,nn,nn) :: R_el
                                  

                                     
      
      do r = 1,nn
         do q = 1,nn
            do p = 1,nn

               is = nn*nn*(r -1) +nn*(q -1) +p
               in = cs(cs(ielem -1) +is)

                if (yon(p,q,r).eq.1) then

                    VS = dsqrt( mu_el(p,q,r) / rho_el(p,q,r) )
                    VP = dsqrt( (lambda_el(p,q,r) + 2*mu_el(p,q,r) ) / rho_el(p,q,r) ) 
                    Vratio = VP/VS 

                    !-------------------------------------------------------------------
                    !
                    if(tcase .ne. 16) then
                      do fn = 1,nf                        
                         if (vcase .eq. tag_func(fn)) then                                                            
                             if (func_type(fn) .eq. 60 ) then     
                 
                                reduction_factor = get_func_value(nf,func_type,func_indx,func_data, nfdata, &  
                                                               fn,R_el(p,q,r))
                                                                                
                                mu_el(p,q,r) = mu_el(p,q,r) *  reduction_factor                      
                                VS = dsqrt(mu_el(p,q,r)/rho_el(p,q,r))                   
                                VP = VS * Vratio                     
                                lambda_el(p,q,r)= rho_el(p,q,r)*(VP**2)-2*mu_el(p,q,r) 
                                                
                             endif                                                                        
                         endif                                                                          
                       enddo
                    
                     elseif (tcase .eq. 16) then
                     
                       do fn = 1,nf                        
                         if (vcase .eq. tag_func(fn)) then                                                            
                             if (func_type(fn) .eq. 60 .and. vs_tria(in) .lt. 325.d0) then     

                               reduction_factor = get_func_value(nf,func_type,func_indx,func_data, nfdata, &  
                                                                 fn,R_el(p,q,r))
                                                                                
                               mu_el(p,q,r) = mu_el(p,q,r) *  reduction_factor                      
                               VS = dsqrt(mu_el(p,q,r)/rho_el(p,q,r))                   
                               VP = VS * Vratio                     
                               lambda_el(p,q,r)= rho_el(p,q,r)*(VP**2)-2*mu_el(p,q,r) 

                             elseif (func_type(fn) .eq. 62 .and. vs_tria(in) .lt. 450.d0) then     
                               
                               reduction_factor = get_func_value(nf,func_type,func_indx,func_data, nfdata, &  
                                                                 fn,R_el(p,q,r))
                                                                                
                               mu_el(p,q,r) = mu_el(p,q,r) *  reduction_factor                      
                               VS = dsqrt(mu_el(p,q,r)/rho_el(p,q,r))                   
                               VP = VS * Vratio                     
                               lambda_el(p,q,r)= rho_el(p,q,r)*(VP**2)-2*mu_el(p,q,r) 
                               
                              
                                                                             
                             endif
                                                                                                     
                          endif                                                                          
                        enddo
                     
                     endif!tcase = 16
                      
                      
                        !
                        !-------------------------------------------------------------------

                endif !if (yon(i,j,k).eq.1) then
               
            enddo
         enddo
      enddo
      
      return
      
      end subroutine MAKE_ELTENSOR_FOR_CASES_NLE

