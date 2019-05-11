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

!> @brief Assignes material properties node by node. 
!! @author Ilario Mazzieri
!> @date September, 2014 
!> @version 1.0
!> @param[in] tcase label for not honoring case  
!> @param[in] vcase value case for non linear block 
!> @param[in] nn  polynomial degree
!> @param[in] nn_loc number of local nodes
!> @param[in] zs_elev elevation of the nodes from topography
!> @param[in] zs_all elevation of the nodes from alluvial
!> @param[in] vs_nodes vs30 of the nodes 
!> @param[in] thick_nodes thickness of sediments
!> @param[in] cs_nnz_loc length of cs_loc
!> @param[in] cs_loc connectivity vector
!> @param[in] ielem  element index
!> @param[in] nf number of functions
!> @param[in] sub_tag_all tags for multi-not honoring
!> @param[in] xs -coordinate of the nodes 
!> @param[in] ys y-coordinate of the nodes 
!> @param[in] zs z-coordinate of the nodes 
!> @param[in] mpi_id MPI process id 
!> @param[in] local_n_num local node numeration
!> @param[in] damping_tpe 1-Kosloff&Kosloff, 2-Standard Linear Solid
!> @param[in] ierr 1-debug mode, 0-standard mode
!> @param[out] rho_el material density 
!> @param[out] lambda_el Lame coefficient lambda
!> @param[out] mu_el Lame coefficient mu
!> @param[out] gamma_el damping coefficient gamma (Kosloff&Kosloff)
!> @param[out] qs quality factor for s-waves (Standard Linear Solid)
!> @param[out] qp quality factor for p-waves (Standard Linear Solid)

      subroutine MAKE_ELTENSOR_FOR_CASES(tcase,vcase,&
                                 nn,rho_el,lambda_el,mu_el,gamma_el,&
                                 nn_loc, zs_elev, zs_all, vs_nodes, thick_nodes, &
                                 cs_nnz_loc, cs_loc, ielem, &
                                 sub_tag_all, zs, mpi_id, local_n_num, &
                                 damping_type, qs, qp, &
                                 xs, ys, ierr)
 
      
      implicit none
                                                      
      integer*4 :: tcase, ierr                
      integer*4 :: vcase, mpi_id        
      integer*4 :: nn
      integer*4 :: p, q, r, ic
      integer*4 :: nn_loc
      integer*4 :: cs_nnz_loc                                
      integer*4 :: is,in,ielem , damping_type                               

      integer*4, dimension(nn_loc) :: local_n_num
      integer*4, dimension(nn_loc) :: sub_tag_all
      integer*4, dimension(0:cs_nnz_loc) :: cs_loc                

      real*8 :: Depth, Depth_real, vs_all, vp_all, thickness, vs30
      real*8 :: VS,VP,rho,lambda,mu,gamma,ni, qs,qp
      real*8 :: x1,y1,x2,y2,coef_a, coef_b, coef_c, numer, den, distance, f_distance

      real*8, dimension(nn_loc) :: zs_elev
      real*8, dimension(nn_loc) :: zs_all
      real*8, dimension(nn_loc) :: vs_nodes, thick_nodes

      real*8, dimension(nn_loc) :: zs, xs, ys                

      real*8, dimension(nn,nn,nn) :: rho_el,lambda_el,mu_el,gamma_el
          
      
!     STRESS CALCULATION
      
      if (ierr .eq. 1) then
          open(1000,file='RS_1.out',position='APPEND')
!          open(1001,file='RS_2.out',position='APPEND')
!          open(1002,file='RS_3.out',position='APPEND')
!          open(1003,file='RS_4.out',position='APPEND')
!          open(1004,file='RS_5.out',position='APPEND')
!          open(1005,file='RS_6.out',position='APPEND')
!          open(1006,file='RS_7.out',position='APPEND')
!          open(1007,file='RS_8.out',position='APPEND')
!          open(1008,file='RS_9.out',position='APPEND')
!          open(1009,file='RS_10.out',position='APPEND')
!          open(1010,file='RS_11.out',position='APPEND')
          
      endif
      
      do r = 1,nn
         do q = 1,nn
            do p = 1,nn
                  is = nn*nn*(r -1) +nn*(q -1) +p
                  in = cs_loc(cs_loc(ielem -1) +is)
                  ic = in

                  if (ic .eq. 0 ) write(*,*) 'Error in MAKE_ELTENSOR_FOR_CASES '

                     if (tcase.eq.1) then

                        !-----------------------------------------------------------------
                           ! CASE 1: GRENOBLE honoring
                                                
                        Depth = zs_elev(ic)
                        if (Depth .lt. 0.0d0) Depth = 0.0d0

                               VS  = 300.0 + 19.0 * sqrt(Depth)        !VS: S velocity in m/s
                        VP  = 1450.0 + 1.2 * Depth                !VP: P velocity in m/s 
                        rho = 2140.0 + 0.125 * Depth                !RHO: MASS DENSITY in kg/m^3
                               lambda = rho * (VP**2 - 2*VS**2)
                               mu = rho * VS**2
                        gamma = 6.2832E-02
                        
                        !-------------------------------------------------------------------

                   elseif (tcase.eq.2) then

                        !-------------------------------------------------------------------
                               ! CASE 2: GRENOBLE NOT honoring
                        !
                        ! + MATERIAL INSIDE THE ALLUVIAL BASIN
                        
                        Depth = zs_elev(ic)        
                        if ((Depth .ge. 0.0d0) .and. (zs_all(ic) .ge. 0.0d0)) then
                                VS  = 300.0 + 19.0 * sqrt(Depth)        !VS: S velocity in m/s
                                  VP  = 1450.0 + 1.2 * Depth                !VP: P velocity in m/s 
                                rho = 2140.0 + 0.125 * Depth                !RHO: MASS DENSITY in kg/m^3
                                       lambda = rho * (VP**2 - 2*VS**2)
                                       mu = rho * VS**2
                                gamma = 6.2832E-02
                
                        ! + MATERIAL INTO THE BEDROCK (FIRST LAYER) 
                        else
                                lambda = 2.9594E+10
                                mu = 2.7853E+10
                                rho = 2720.0d0
                                gamma = 0.0d0
                        endif
                
                        !-------------------------------------------------------------------

                  elseif (tcase.eq.3) then

                        !-------------------------------------------------------------------
                               ! CASE 3: GUBBIO NOT honoring
                        !
                        ! + MATERIAL INSIDE THE ALLUVIAL BASIN
                
                        Depth = zs_elev(ic)                        !D: depth in m
                        if ((Depth .ge. 0.0d0) .and. (zs_all(ic) .ge. 0.0d0)) then
                                VS  = 250.0 + 30.0 * sqrt(Depth)        !VS: S velocity in m/s
                                VP  = 1000.0 + 30.0 * sqrt(Depth)       !VP: P velocity in m/s
                                rho = 1900.0                            !RHO: MASS DENSITY in kg/m^3
                                lambda = rho * (VP**2 - 2*VS**2)
                                mu = rho * VS**2
                                gamma = 6.2832E-02
                        else
                   
                       ! + MATERIAL INTO THE BEDROCK (FIRST LAYER) 
                                      lambda = 1.2694E+10 
                                mu = 7.1280E+09 
                                rho = 2200.0 
                                gamma = 3.9270E-02 
                        endif   
                                            
                        !
                        !-------------------------------------------------------------------
                                
                  elseif (tcase.eq.4) then


                        !-------------------------------------------------------------------
                               ! CASE 4: SULMONA NOT honoring
                        !
                        ! + MATERIAL INSIDE THE ALLUVIAL BASIN
                        
                        Depth = zs_elev(ic)                        !D: depth in m
                        if ((Depth .ge. 0.0d0) .and. (zs_all(ic) .ge. 0.0d0)) then
                                                        
                                  VS  = 500.0 + 19.0 * sqrt(Depth)        !VS: S velocity in m/s
                                  VP  = 1000.0 + 1.2 * Depth                !VP: P velocity in m/s 
                                 rho = 1900.0 + 0.125 * Depth                !RHO: MASS DENSITY in kg/m^3
                                        lambda = rho * (VP**2 - 2*VS**2)
                                        mu = rho * VS**2
                                 gamma = 6.2832E-02
                        else
                                 lambda = 5.760000E+09
                                 mu = 2.880000E+09
                                 rho = 2000.00
                                 gamma = 4.188790E-02
                        endif
                
                        !
                        !-------------------------------------------------------------------
                                
                   elseif (tcase.eq.5) then

                        !-------------------------------------------------------------------
                               ! CASE 5: VOLVI CASHIMA benchmark -  NOT honoring
                        !
                        !-------------------------------------------------------------------
                        ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 1st Layer
                
                        if (sub_tag_all(ic).eq.1) then
                                rho = 2100
                                       lambda = 4.557000E+09
                                       mu = 8.400000E+07
                                gamma = 2.094395E-01
                        
                        ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 2nd Layer
                        elseif (sub_tag_all(ic).eq.2) then
                                rho = 2100
                                       lambda = 6.289500E+09
                                       mu = 2.572500E+08
                                gamma = 1.196797E-01 
                                
                        ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 3rd Layer
                        elseif (sub_tag_all(ic).eq.3) then
                                rho = 2200
                                       lambda = 1.189100E+10
                                       mu = 9.295000E+08
                                gamma = 6.444293E-02
                        
                        ! + MATERIAL INTO THE BEDROCK 
                        elseif (sub_tag_all(ic).eq.4) then
                                if (zs(ic).ge.-3000.0) then
                                        VS  = 0.4100*(-zs(ic)) + 2190.0000        !VS: S velocity in m/s
                                        VP  = 0.8100*(-zs(ic)) + 3690.0000        !VP: P velocity in m/s 
                                        rho = 0.0680*(-zs(ic)) + 2532.0000        !RHO: MASS DENSITY in kg/m^3
                                               lambda = rho * (VP**2 - 2*VS**2)
                                               mu = rho * VS**2
                                        gamma = 1.6111E-02
                                else
                                        VS  = 0.0050*(-zs(ic)) + 3405.0000        !VS: S velocity in m/s
                                        VP  = 0.0050*(-zs(ic)) + 6105.0000        !VP: P velocity in m/s 
                                        rho = 0.0040*(-zs(ic)) + 2724.0000        !RHO: MASS DENSITY in kg/m^3
                                               lambda = rho * (VP**2 - 2*VS**2)
                                               mu = rho * VS**2
                                        gamma = 1.6111E-02
                                endif
                                
                        ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 5th Layer (nu modified)
                        elseif (sub_tag_all(ic).eq.5) then
                                rho = 2100
                                       lambda = 1.260000E+08
                                       mu = 8.400000E+07
                                gamma = 2.094395E-01
                        
                        ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 6th Layer (nu modified)
                        elseif (sub_tag_all(ic).eq.6) then
                                rho = 2100
                                       lambda = 3.858750E+08
                                       mu = 2.572500E+08
                                gamma = 1.196797E-01 
                                
                        ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 7th Layer (nu modified)
                        else 
                                rho = 2200
                                       lambda = 1.394250E+09
                                       mu = 9.295000E+08
                                gamma = 6.444293E-02
                        endif
                        
                        !
                        !-------------------------------------------------------------------
                                        
                  elseif (tcase.eq.6) then
                        !-------------------------------------------------------------------
                               ! CASE 6: FRIULI NOT honoring (Tagliamento river valley) 
                        !
                        ! + MATERIAL INSIDE THE ALLUVIAL BASIN
                        
                        Depth = zs_elev(ic)                        
                        if ((Depth .ge. 0.0d0) .and. (zs_all(ic).ge.0.0d0)) then    
                                ! Option #1
                                ! VS  = 300.0d0 + 6.50d0*Depth  !VS: S velocity in m/s     
                                ! VP  =  600.0d0 + 12.0d0*Depth         !VP: P velocity in m/s                  
                                ! rho =  1900.0d0 + 1.25*Depth        !RHO: MASS DENSITY in kg/m^3            
                                
                                ! Option #2
                                 VS  = 300.0d0 + 30.0d0*(Depth)**0.67  !VS: S velocity in m/s     
                                 VP  =  600.0d0 + 30.0d0*(Depth)**0.77         !VP: P velocity in m/s                  
                                 rho =  1900.0d0 + 1.25*Depth        !RHO: MASS DENSITY in kg/m^3
                                 lambda = rho * (VP**2 - 2*VS**2)                   
                                 mu = rho * VS**2                                   
                                 gamma = 3.1416E-02 ! (Qs=100)                       
                        else                                                    
                                ! + MATERIAL INSIDE THE BEDROCK (Vs=1500m/s)            
                                 lambda = 1.0350E+10                                
                                 mu = 5.1750E+09                                    
                                 rho = 2300.0d0                                     
                                 gamma = 2.0944E-02                                 
                        endif                                                   
                                        
                                        
                  elseif (tcase.eq.7) then
                        !-------------------------------------------------------------------
                               ! CASE 7: AQUILA NOT honoring 
                        !
                        ! + MATERIAL INSIDE THE ALLUVIAL BASIN
                        
                        Depth = zs_elev(ic)                        !D: depth in m              
                        if ((Depth .ge. 0.0d0) .and. (zs_all(ic) .ge. 0.0d0)) then    
                                 ! Option #1
                                 VS  = 500.0d0 + 10.0d0*(Depth)**0.50  !VS: S velocity in m/s     
                                 VP  =  sqrt(3.0)*VS         !VP: P velocity in m/s  (nu = 0.25)        
                                 rho =  2000.0d0               !RHO: MASS DENSITY in kg/m^3                              
                                                                
                                 lambda = rho * (VP**2 - 2*VS**2)                   
                                 mu = rho * VS**2                                   
                                 gamma = 6.2832E-02 ! (Qs=Vs/10=50)                 
                        else                                                    
                                 ! + MATERIAL INSIDE THE BEDROCK (Vs=1500m/s)            
                                 lambda = 1.2189E+10                                
                                 mu = 6.0943E+09                                    
                                 rho = 2600.0d0                                     
                                 gamma = 3.1416E-02                                                                          
                        endif                                                           
                                        
    		  elseif (tcase.eq.8) then
                        !-------------------------------------------------------------------
                        ! CASE 8: SANTIAGO NOT honoring 
			!
			! + MATERIAL INSIDE THE ALLUVIAL BASIN
			
			Depth = zs_elev(ic)			!D: depth in m              
			if ((Depth .ge. 0.0d0) .and. (zs_all(ic) .ge. 0.0d0)) then    
			         VS  = 400.0d0 + 55.0d0*(Depth)**0.50   ! VS: S velocity in m/s     
				 VP  = 1730.0d0 + 60.0d0*(Depth)**0.50  ! VP: P velocity in m/s         
				 rho =  2100.0d0 + 0.15d0*(Depth)      	! RHO: MASS DENSITY in kg/m^3
				 lambda = rho * (VP**2 - 2*VS**2) 
				 mu = rho * VS**2
				 qs = 0.1*VS;           
                                 gamma = (3.1415*(2.d0/3.d0))/qs !hy: fpeak = 2/3 Hz
				 
				 
			else                                                    
				 ! + MATERIAL AT OUCROPPING BEDROCK (Vs=2400m/s)           
				 lambda = 2.5368E+10                               
				 mu = 1.3824E+10                                    
				 rho = 2400.0d0 
				 VS = 2400.d0;  
				                                   
				 qs = 0.1*VS;           
                                 gamma = (3.1415*(2.d0/3.d0))/qs !hy: fpeak = 2/3 Hz

			endif 


                   elseif (tcase.eq.10) then
                            !-------------------------------------------------------------------
                               ! CASE 10: CHRISTCHURCH INGV - NH Staircase
                        !
                        ! + MATERIAL INSIDE THE ALLUVIAL BASIN

                        Depth = zs_elev(ic)                        !D: depth in m               
                        if ((Depth .ge. 0.0d0) .and. (zs_all(ic) .ge. 0.0d0)) then     

                              if (Depth.le.300.0d0) then       ! 0 < z < 300         
                                         rho =  1700.0d0                                                             
                                         lambda = 2.9700E+08       !Vs = 300 m/s               
                                         Mu = 1.5300E+08           !Vp = 596 m/s               
                                         gamma = 2.9920E-02        ! Qs = 70 (2 Hz)            
                                                 
                                                 
                               elseif ((Depth.gt.300.0d0).and.(Depth.le.700.0d0)) then  ! 300 < z < 700     
                                         rho =  2000.0d0                                                     
                                         lambda = 3.0000E+09       !Vs = 1000 m/s      
                                         mu = 2.0000E+09          !Vp = 1871 m/s      
                                         gamma = 2.0944E-02 ! Qs = 100 (2 Hz)          
                              
                               elseif ((Depth.gt.700.0d0)) then ! 700 < z < 1500    
                                         rho =  2300.0d0                                                                 
                                         lambda = 1.1178E+10      !Vs = 1800 m/s                
                                         mu = 7.4520E+09          !Vp = 3368 m/s               
                                         gamma = 2.0944E-02 ! Qs = 100 (2 Hz)                              
                               endif                                                                     
                        else                            

                                ! + MATERIAL INSIDE THE BEDROCK (Vs=3175 m/s)           
                                         lambda = 2.6217E+10                                
                                         mu = 2.6217E+10                                    
                                         rho = 2600.0d0                                     
                                         gamma = 1.0472E-02      
                        endif        
                                
                    elseif (tcase.eq.11) then
                   !-------------------------------------------------------------------
                          ! CASE 11: New Chch 

                         Depth = zs_elev(ic)
                         
                         if ((Depth.ge.0.0d0).and.(zs_all(ic).ge.0.0d0)) then                  
                                                 
                                if (Depth.lt.15.0d0) then               
                                         VS = 270.d0                                 
                                         ni = 0.45d0                                                           
                                         VP  =  (((2.0d0*(1.0d0 - ni))/(1.0d0-2.0d0*ni))**0.5)*VS 
                                         rho = 1700.d0      
                                         lambda = rho * (VP**2 - 2*VS**2)                                                    
                                         mu = rho * VS**2                                
                                         gamma = (3.1415*(2/3))/(70.d0)                                                      
                                 elseif (Depth.lt.50.0d0) then                                                       
                                         VS = 270.d0+11.5d0*(Depth-15.d0)                                                    
                                         ni = 0.45 - 0.0025*(Depth-15.d0)                                                    
                                         VP  =  (((2.0d0*(1.0d0 - ni))/(1.0d0-2.0d0*ni))**0.5)*VS
                                         rho = 1700.d0 + 5.d0*(Depth-15.d0)                                                  
                                         lambda = rho * (VP**2 - 2*VS**2)                                                    
                                         mu = rho * VS**2  
                                         gamma = (3.1415*(2/3))/(70.d0+0.5d0*(Depth-15.d0))                                  
                                 else                                                                                
                                         VS= 270.d0 + 11.5d0*(50.d0-15.d0) + 0.7d0*(Depth-50)                                
                                         ni= 0.45d0 - 0.0025d0*(50.d0-15.d0) - 0.000075d0*(Depth-50)                         
                                         VP  =  (((2.0d0*(1.0d0 - ni))/(1.0d0-2.0d0*ni))**0.5)*VS                            
                                         rho= 1700.d0 - 5.d0*(50.d0-15.d0) +0.5d0*(Depth-50)                                 
                                         lambda = rho * (VP**2 - 2*VS**2)                                                    
                                         mu = rho * VS**2                                                       
                                         gamma = (3.1415*(2/3))/(70.d0 + 0.5d0*(50.d0-15.d0) + 0.0775d0*(Depth-50))          
                          endif                                                                                      
                                                                                                                                                     
                          else          
                             Depth_real = abs(zs(ic))
                                                                                                                                               
                                 if (Depth_real.lt.15.0d0) then                                                                 
                                         VS = 750.d0                                                                 
                                         ni = 0.30d0                                                        
                                         VP  =  (((2.0d0*(1.0d0 - ni))/(1.0d0-2.0d0*ni))**0.5)*VS                            
                                         rho = 2000.d0                                             
                                         lambda = rho * (VP**2 - 2*VS**2)                                                    
                                         mu = rho * VS**2                                
                                         gamma = (3.1415*(2/3))/(100.d0)                                                     
                                 elseif (Depth_real.lt.50.0d0) then                                                       
                                         VS = 750.d0+14.d0*(Depth_real-15.d0)              
                                         ni = 0.30d0 - 0.0005d0*(Depth_real-15.d0)  
                                         VP  =  (((2.0d0*(1.0d0 - ni))/(1.0d0-2.0d0*ni))**0.5)*VS                            
                                         rho = 2000.d0 + 6.5d0*(Depth_real-15.d0) 
                                         lambda = rho * (VP**2 - 2*VS**2)                                                    
                                         mu = rho * VS**2       
                                         gamma = (3.1415*(2/3))/(100.d0+0.8d0*(Depth_real-15.d0))  
                                 else              
                                         VS= 750.d0 + 14.d0*(50.d0-15.d0) + 1.1d0*(Depth_real-50)         
                                         ni= 0.30d0 - 0.0005d0*(50.d0-15.d0) - 0.000022d0*(Depth_real-50) 
                                         VP  =  (((2.0d0*(1.0d0 - ni))/(1.0d0-2.0d0*ni))**0.5)*VS                            
                                         rho= 2000.d0 + 6.5d0*(50.d0-15.d0) + 0.26d0*(Depth_real-50)
                                         lambda = rho * (VP**2 - 2*VS**2)                                                    
                                         mu = rho * VS**2  
                                         gamma = (3.1415*(2/3))/(100.d0 + 0.8d0*(50.d0-15.d0) + 0.049d0*(Depth_real-50))
                                 endif                                                  
                          endif                                                                                  

                    elseif (tcase.eq.12) then
                    !-------------------------------------------------------------------
                    ! CASE 12: PO PLAIN (new model)  
                    !-------------------------------------------------------------------
                        
                         Depth = zs_elev(ic)
                         !straight line passing by (x1,y1,0), (x2,y2,0)
                         x1 =  654957.002352;  y1 = 4974060.299450;
                         x2 =  688420.525202;  y2 = 4957613.600935;
                         !coefficient of the line ax + by + c = 0
                         coef_a = 1.d0/(x2-x1);
                         coef_b = 1.d0/(y1-y2);
                         coef_c = - y1/(y1-y2) + x1/(x1-x2);
                         !distance between (x,y,0) and the line ax + by + c = 0
                         numer = coef_a*xs(ic) + coef_b*ys(ic) + coef_c
                         den = dsqrt(coef_a**2 + coef_b**2)
                         distance = dabs(numer/den)
                         
                         f_distance = 150.d0 + 1850.d0/(1.d0 + dexp(-0.0012*(distance-5000.d0)));    
                         
                         if ((Depth.ge.0.0d0).and.(zs_all(ic).ge.0.0d0)) then                                    
                                        ! + MATERIAL INSIDE THE BASIN 
                                if (Depth .le. 150.0d0) then       
                                         VS = 300.d0                   
                                         VP = 1500.d0                  
                                         rho = 1800.d0                  
                                         lambda = rho * (VP**2 - 2*VS**2)
                                         mu = rho * VS**2   
                                         qs = 0.1*VS;           
                                         gamma = (3.1415*(2.d0/3.d0))/qs !hy: fpeak = 2/3 Hz

                                elseif(Depth .gt. 150.d0 .and. Depth .le. f_distance)  then          
                                         VS  =  300.d0 + 10.d0*(Depth-150.d0)**0.5 
                                         VP =  1500.d0 + 10.d0*(Depth-150.d0)**0.5
                                         rho = 1800.d0 +  6.d0*(Depth-150.d0)**0.5
                                         lambda = rho * (VP**2 - 2*VS**2)
                                         mu = rho * VS**2              
                                         qs = 0.1*vs;
                                         gamma = (3.1415*(2.d0/3.d0))/qs                

                                else                                   
                                         VS =   800.d0 + 15.d0*(Depth-f_distance)**0.5                                 
                                         VP  = 2000.d0 + 15.d0*(Depth-f_distance)**0.5                           
                                         rho = 2100.d0 +  4.d0*(Depth-f_distance)**0.5                                 
                                         lambda = rho * (VP**2 - 2*VS**2)                                                    
                                         mu = rho * VS**2                                          
                                         qs = 0.1*vs
                                         gamma = (3.1415*(2/3))/qs          

                                endif                                                                        
                                                                                                                                                     
                          else  
                          ! + MATERIAL INSIDE THE BEDROCK         
                             Depth_real = abs(zs(ic))
                             if (Depth_real .le. 1000.0d0) then                   
                                    VS = 1200.d0                                                                 
                                    VP  = 2300.d0                          
                                    rho = 2100.d0                                         
                                    lambda = rho * (VP**2 - 2*VS**2)                                                    
                                    mu = rho * VS**2                                 
                                    gamma = (3.1415*(2/3))/(150.d0)                                                     
                                    !if(damping_type .eq. 2) then
                                    !   qs=120.d0; qp=230.d0;
                                    !endif
                                 elseif (Depth_real.le.3000.0d0) then                                                       
                                         VS = 2100.d0                                                                   
                                         VP  = 3500.d0                          
                                         rho = 2200.d0     
                                         lambda = rho * (VP**2 - 2*VS**2)                                                    
                                         mu = rho * VS**2                                 
                                         gamma = (3.1415*(2/3))/(200.d0)                                
                                    !if(damping_type .eq. 2) then
                                    !   qs=200.d0; qp=400.d0;
                                    !endif

                                 elseif (Depth_real.le.6000.0d0) then
                                         VS = 2750.d0                                                                    
                                         VP  = 4750.d0                          
                                         rho = 2400.d0                                           
                                         lambda = rho * (VP**2 - 2*VS**2)                                                    
                                         mu = rho * VS**2                                 
                                         gamma = (3.1415*(2/3))/(250.d0)                                  
                                    !if(damping_type .eq. 2) then
                                    !   qs=250.d0; qp=500.d0;
                                    !endif
                                 else
                                         VS = 3670.d0                                                                   
                                         VP  = 6340.d0                          
                                         rho = 2800.d0                   
                                         lambda = rho * (VP**2 - 2*VS**2)                                                    
                                         mu = rho * VS**2                                
                                         gamma = (3.1415*(2/3))/(350.d0)                                                
                                    !if(damping_type .eq. 2) then
                                    !   qs=350.d0; qp=700.d0;
                                    !endif
                                         
                                  endif
                          endif             

            elseif (tcase.eq.13) then
       !-------------------------------------------------------------------
       ! CASE 13: PO PLAIN-BEDROCK (new model)  
       !-------------------------------------------------------------------
                        
                         Depth = zs_elev(ic)
                         
                         if ((Depth.ge.0.0d0).and.(zs_all(ic).ge.0.0d0)) then
                                ! + MATERIAL INSIDE THE BASIN 
                               if(Depth .le. 150.0d0) then       
                                         VS = 340.d0                                 
                                         VP = 1500.d0                                           
                                         rho = 1800.d0                                          
                                         lambda = rho * (VP**2 - 2*VS**2)                                                    
                                         mu = rho * VS**2                                
                                         gamma = (3.1415*(2.d0/3.d0))/(35.d0) !hy: fpeak = 2/3 Hz

                                elseif (Depth .le. 500.0d0) then                               
                                         VS = 800.d0                                                     
                                         VP = 1800.d0                                        
                                         rho = 2100.d0                                                
                                         lambda = rho * (VP**2 - 2*VS**2)                                                    
                                         mu = rho * VS**2
                                         gamma = (3.1415*(2/3))/(80.d0)                                                       
  
                                elseif(Depth .le. 1000.0d0) then                               
                                         VS = 1200.d0                                                     
                                         VP = 2300.d0                                        
                                         rho = 2100.d0                                                
                                         lambda = rho * (VP**2 - 2*VS**2)                                                    
                                         mu = rho * VS**2
                                         gamma = (3.1415*(2/3))/(250.d0)                                                       

                                elseif(Depth .le. 3000.0d0) then                               
                                         VS = 2100.d0                                                     
                                         VP = 3500.d0                                        
                                         rho = 2200.d0                                                
                                         lambda = rho * (VP**2 - 2*VS**2)                                                    
                                         mu = rho * VS**2
                                         gamma = (3.1415*(2/3))/(200.d0)                                                       

                                elseif(Depth .le. 6000.0d0) then                               
                                         VS = 2750.d0                                                     
                                         VP = 4750.d0                                        
                                         rho = 2400.d0                                                
                                         lambda = rho * (VP**2 - 2*VS**2)                                                    
                                         mu = rho * VS**2
                                         gamma = (3.1415*(2/3))/(250.d0)                                                       
                                
                                 else 
                                         VS = 3670.d0                           
                                         VP  = 6340.d0                          
                                         rho = 2800.d0
                                         lambda = rho * (VP**2 - 2*VS**2)                                                    
                                         mu = rho * VS**2                                 
                                         gamma = (3.1415*(2/3))/(350.d0)              
                                endif 
                                                                                                                                                     
                          endif             

                                
                        elseif (tcase.eq.14) then
                        !-------------------------------------------------------------------
                        ! CASE 14: Wellington NOT honoring 
                        !
                        ! simplified model, Benites et al. 2005
                        
                        Depth = zs_elev(ic)                        !D: depth in m              
                        if ((Depth .ge. 0.0d0) .and. (zs_all(ic) .ge. 0.0d0)) then    
                           
                           ! + MATERIAL INSIDE THE ALLUVIAL BASIN 
                            if (Depth .le. 100.0d0) then 
                                VS = 300.d0                                                     
                                VP = 520.d0                                        
                                rho = 2200.d0                                                
                                lambda = rho * (VP**2 - 2.d0*VS**2)                                                    
                                mu = rho * VS**2  
                                gamma = (3.1415*(3.d0/3.d0))/(30.d0)                                                       
                            
                            elseif (Depth .le. 200.0d0) then 
                                VS = 400.d0
                                VP = 700.d0                          
                                rho = 2300.d0 
                                lambda = rho * (VP**2 - 2.d0*VS**2)                                                    
                                mu = rho * VS**2
                                gamma = (3.1415*(2.d0/3.d0))/(40.d0)  
                            
                            elseif (Depth .le. 500.0d0) then                           
                               VS = 500.d0     
                               VP = 850.d0                          
                               rho = 2400.d0                                  
                               lambda = rho * (VP**2 - 2.d0*VS**2)                                                    
                               mu = rho * VS**2                          
                               gamma = (3.1415*(2.d0/3.d0))/(50.d0)
                                                 
                            else 
                               VS = 1000.d0                        
                               VP  = 1700.d0                          
                               rho = 2400.d0       
                               lambda = rho * (VP**2 - 2.d0*VS**2)                                                    
                               mu = rho * VS**2                          
                               gamma = (3.1415*(2.d0/3.d0))/(100.d0)                                                 
                                                                     
                            endif                                               
                        else                                                    
                           ! + MATERIAL INSIDE THE BEDROCK (Vs=1500m/s)            
                           VS = 1500.d0                                                     
                           VP = 2800.d0                                        
                           rho = 2400.d0                                                
                           lambda = rho * (VP**2 - 2.d0*VS**2)                                                    
                           mu = rho * VS**2                                                                                 
                           gamma = (3.1415*(3.d0/3.d0))/(150.d0)                                                               
                        endif                                   
                        
                   elseif (tcase.eq.15) then
                    !-------------------------------------------------------------------
                    ! CASE 15: MARSICA  
                    !-------------------------------------------------------------------
                        
                         Depth = zs_elev(ic)                       
                         if ((Depth .ge. 0.0d0).and.(zs_all(ic) .ge. 0.0d0)) then                                    
                              ! + MATERIAL INSIDE THE BASIN 
                              VS = 100.d0 + 10.d0 * Depth**(0.6d0) !100.d0 + 10.d0 * Depth**(0.6d0)
                              VP  = dsqrt(10.d0)*VS
                              rho = 1530.d0 + 0.1d0*Depth**(0.54d0)                               
                              lambda = rho * (VP**2 - 2*VS**2)                                                    
                              mu = rho * VS**2                                          
                              qs = 0.1*vs
                              gamma = (3.1415*(2/3))/qs          
                                              
                              if (ierr .eq. 1)   write(1000,*) xs(ic),ys(ic),zs(ic), VS, VP                
                          else  
                             ! + MATERIAL INSIDE THE BEDROCK         
                             Depth_real = zs(ic)
                             if (Depth_real .ge. -500.0d0) then
                                    VS = 1000.d0                                                                 
                                    VP  = 1800.d0                          
                                    rho = 2300.d0                                         
                                    lambda = rho * (VP**2 - 2*VS**2)                                                    
                                    mu = rho * VS**2 
                                    qs = 0.1*vs                                
                                    gamma = (3.1415*(2.d0/3.d0))/qs       
                              if (ierr .eq. 1)   write(1001,*) xs(ic),ys(ic),zs(ic), VS, VP                
                                                                                      
                             elseif (Depth_real .le. -500.d0 .and. Depth_real .ge. -1000.0d0) then                   
                                    VS = 1700.d0                                                                 
                                    VP  = 3160.d0                          
                                    rho = 2500.d0                                         
                                    lambda = rho * (VP**2 - 2*VS**2)                                                    
                                    mu = rho * VS**2 
                                    qs = 0.1*vs                                
                                    gamma = (3.1415*(2.d0/3.d0))/qs   
                              if (ierr .eq. 1)   write(1002,*) xs(ic),ys(ic),zs(ic), VS, VP                
                                                                  
                             else
                                    VS = 2600.d0;
                                    VP = 4830.d0;
                                    rho = 2840.d0;
                                    lambda = rho * (VP**2 - 2*VS**2)                                                    
                                    mu = rho * VS**2 
                                    qs = 0.1*vs                                
                                    gamma = (3.1415*(2.d0/3.d0))/qs   
                              if (ierr .eq. 1)   write(1003,*) xs(ic),ys(ic),zs(ic), VS, VP                
                                                                  
                             endif
                         endif             
                                                                
                   elseif (tcase.eq.16) then
                    !-------------------------------------------------------------------
                    ! CASE 16: ISTANBUL 
                    !-------------------------------------------------------------------
                        
                         Depth = zs_elev(ic)
                         thickness  = thick_nodes(ic);
                         vs30 = vs_nodes(ic)

                         
                         if( vs30 .lt. 325.d0) then
                             if ( Depth .le. 160.d0) then
                                VS = 250.d0 + 43.d0*Depth**(0.5d0);
                                VP  = 700.d0 + 45.d0*Depth**(0.5d0);
                                rho = 1800.d0;
                                lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                                mu = rho * VS**2.d0;
                                qs = 0.1d0*VS;
                                gamma = 4.d0*datan(1.d0)/qs;
                             elseif(Depth .le. 2000.d0) then
                                VS = 800.d0 + 37.13d0*(Depth-160d0)**(0.5)
                                VP  = VS*1.6;
                                rho = 1800 + 12.92d0*(Depth-160d0)**(0.5);
                                lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                                mu = rho * VS**2.d0;
                                qs = 0.1d0*VS;
                                gamma = 4.d0*datan(1.d0)/qs;
                            else
                                VS = 1350.d0 + 23.33*(Depth)**(0.5);
                                VP = VS*1.6;
                                rho =  2100 + 5.69*(Depth)**(0.5);
                                lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                                mu = rho * VS**2.d0;
                                qs = 0.1*VS;
                                gamma = 4.d0*datan(1.d0)/qs;
                            endif
                         
                         elseif (vs30 .lt. 500.d0) then
                         
                            if ( Depth .le. 80.d0) then
                                VS = 325.d0 + 30.74*Depth**(0.5);
                                VP  = 800.d0 + 42 *Depth**(0.5);;
                                rho = 1850.d0;
                                lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                                mu = rho * VS**2.d0;
                                qs = 0.1*VS;
                                gamma = 4.d0*datan(1.d0)/qs;               
                            elseif(Depth .le. 120.d0) then
                                VS = 600.d0 + 31.62*(Depth-80.d0)**(0.5);
                                VP  = 1175.d0 + 26.72*(Depth-80.d0)**(0.5);
                                rho =  1850.d0
                                lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                                mu = rho * VS**2.d0;
                                qs = 0.1*VS;
                                gamma = 4.d0*datan(1.d0)/qs;
                            elseif(Depth .le. 250.d0) then
                                VS = 800.d0 + 40.75*(Depth-120.d0)**(0.5);
                                VP  = VS*1.6;
                                rho = 1850.d0 + 9.64*(Depth-120.d0)**(0.5);
                                lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                                mu = rho * VS**2.d0;
                                qs = 0.1*VS;
                                gamma = 4.d0*datan(1.d0)/qs;  
                            elseif(Depth .le. 2000.d0) then
                                VS = 700.d0 + 38.14*(Depth-30.d0)**(0.5);
                                VP  = VS*1.6;
                                rho = 1960.d0 + 9.43*(Depth-250.d0)**(0.5);
                                lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                                mu = rho * VS**2.d0;
                                qs = 0.1*VS;
                                gamma = 4.d0*datan(1.d0)/qs;  
                                
                            else
                                VS = 1350.d0 + 23.33*(Depth)**(0.5);
                                VP = VS*1.6;
                                rho =  2100 + 5.69*(Depth)**(0.5);
                                lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                                mu = rho * VS**2.d0;
                                qs = 0.1*VS;
                                gamma = 4.d0*datan(1.d0)/qs;
                            endif          
                         
                         elseif (vs30 .lt. 700.d0) then
                         
                            if ( Depth .le. 50.d0) then
                                VS = 500.d0 + 42.42*Depth**(0.5);
                                VP  = 900.d0 + 42.42*Depth**(0.5);
                                rho = 1900.d0;
                                lambda = rho * (VP**2 - 2*VS**2);
                                mu = rho * VS**2;
                                qs = 0.1*VS;
                                gamma = 4.d0*datan(1.d0)/qs;
                            elseif ( Depth .le. 250.d0) then
                                VS = 800.d0 + 33.1*(Depth-50.d0)**(0.5);
                                VP  = VS*1.6;
                                rho = 1900.d0 + 4.89*(Depth-50.d0)**(0.5);
                                lambda = rho * (VP**2 - 2*VS**2);
                                mu = rho * VS**2;
                                qs = 0.1*VS;
                                gamma = 4.d0*datan(1.d0)/qs;
                            
                           elseif(Depth .le. 2000.d0) then
                                VS = 700.d0 + 38.14*(Depth-30.d0)**(0.5);
                                VP  = VS*1.6;
                                rho = 1960.d0 + 9.43*(Depth-250.d0)**(0.5);
                                lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                                mu = rho * VS**2.d0;
                                qs = 0.1*VS;
                                gamma = 4.d0*datan(1.d0)/qs;
                            else
                                VS = 1350.d0 + 23.33*(Depth)**(0.5);
                                VP = VS*1.6;
                                rho =  2100 + 5.69*(Depth)**(0.5);
                                lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                                mu = rho * VS**2.d0;
                                qs = 0.1*VS;
                                gamma = 4.d0*datan(1.d0)/qs;
                            endif

                         elseif (vs30 .lt. 900.d0) then
                         
                            if ( Depth .le. 4000.d0) then
                                VS = 700.d0 + 37.9*(Depth)**(0.5)
                                VP  = VS*1.6;
                                rho = 1960.d0 + 8.885*(Depth)**(0.5) 
                                lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                                mu = rho * VS**2;
                                qs = 0.1*VS;
                                gamma = 4.d0*datan(1.d0)/qs;
                            else
                                VS = 1350.d0 + 23.33*(Depth)**(0.5);
                                VP = VS*1.6;
                                rho =  2100 + 5.69*(Depth)**(0.5);
                                lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                                mu = rho * VS**2.d0;
                                qs = 0.1*VS;
                                gamma = 4.d0*datan(1.d0)/qs;
                            endif
                
                         elseif (vs30 .lt. 1350.d0) then
                         
                             if ( Depth .le. 2000.d0) then
                                VS = 900.d0 + 33.38 * (Depth)**(0.5);
                                VP  = VS*1.6;
                                rho = 2050.d0 + 215.1*(Depth*0.001)**(0.5);
                                lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                                mu = rho * VS**2.d0;
                                qs = 0.1*VS;
                                gamma = 4.d0*atan(1.d0)/qs;
                             else
                                VS = 1350.d0 + 23.33*(Depth)**(0.5);
                                VP = VS*1.6;
                                rho =  2100 + 5.69*(Depth)**(0.5);
                                lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                                mu = rho * VS**2.d0;
                                qs = 0.1*VS;
                                gamma = 4.d0*datan(1.d0)/qs;
                            endif
                                                                                 
                         else
                                VS = 1350.d0 + 23.33*(Depth)**(0.5);
                                VP = VS*1.6;
                                rho =  2100 + 5.69*(Depth)**(0.5);
                                lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                                mu = rho * VS**2.d0;
                                qs = 0.1*VS;
                                gamma = 4.d0*datan(1.d0)/qs;            
                              
                         endif

                         !left
                         if(dabs(xs(ic) - 576059.d0) .le. 2000.d0) then 
                             VS = 3490.d0;
                             VP = 5770.d0;
                             rho = 2600.d0;
                             lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                             mu = rho * VS**2.d0;
                             qs = 0.1*VS;
                             gamma = 4.d0*atan(1.d0)/qs;             
                         !right
                         elseif(dabs(xs(ic) - 740948.d0) .le. 2000.d0) then 
                             VS = 3490.d0;
                             VP = 5770.d0;
                             rho = 2600.d0;
                             lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                             mu = rho * VS**2.d0;
                             qs = 0.1*VS;
                             gamma = 4.d0*atan(1.d0)/qs;             
                         !up
                         elseif(dabs(ys(ic) - 4602206.d0) .le. 2000.d0) then 
                             VS = 3490.d0;
                             VP = 5770.d0;
                             rho = 2600.d0;
                             lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                             mu = rho * VS**2.d0;
                             qs = 0.1*VS;
                             gamma = 4.d0*atan(1.d0)/qs;             
                         !down
                         elseif(dabs(ys(ic) - 4502679.d0) .le. 2000.d0) then 
                             VS = 3490.d0;
                             VP = 5770.d0;
                             rho = 2600.d0;
                             lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                             mu = rho * VS**2.d0;
                             qs = 0.1*VS;
                             gamma = 4.d0*atan(1.d0)/qs;             
                         endif                          
                                                                         


                   elseif (tcase.eq.20) then
                    !-------------------------------------------------------------------
                    ! CASE 20: ATENE  
                    !-------------------------------------------------------------------
                        
                                                
                         Depth = zs_elev(ic)
                         thickness  = thick_nodes(ic);
                         vs30 = vs_nodes(ic)
                        
                             if ( Depth .le. thickness ) then
                                VS = 350;
                                VP  = 700.d0;
                                rho = 1800.d0;
                                lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                                mu = rho * VS**2.d0;
                                qs = 0.1d0*VS;
                                gamma = 4.d0*datan(1.d0)/qs;
                                if (ierr .eq. 1)   write(1000,*) xs(ic),ys(ic),zs(ic), VS, VP       
                             else
                                VS = 1500
                                VP  = 2600;
                                rho = 2700;
                                lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                                mu = rho * VS**2.d0;
                                qs = 0.1d0*VS;
                                gamma = 4.d0*datan(1.d0)/qs;
                            endif
                                     
                                                  
                         
                   elseif (tcase.eq.21) then
                    !-------------------------------------------------------------------
                    ! CASE 20: BEIJING
                    !-------------------------------------------------------------------
                    
                         Depth = zs_elev(ic)
                         thickness  = thick_nodes(ic);
                         vs30 = vs_nodes(ic)

                         if(zs(ic) .lt. -2000) then
                               VS = 2100 
                               VP = 3500
                               rho = 2200   
                               lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);                           
                               mu = rho * VS**2.d0;
                               qs = 0.1d0*VS;
                               gamma = 4.d0*datan(1.d0)/qs;
                         else
                               if( vs30 .gt. 600.d0) then
                                   VS = vs30 + 5*dsqrt(Depth)
                                   VP = VS*1.6d0
                                   rho = 1800.d0 + 5.d0*dsqrt(Depth);  
                                   lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);                           
                                   mu = rho * VS**2.d0;
                                   qs = 0.1d0*VS;
                                   gamma = 4.d0*datan(1.d0)/qs;
                                   
                               elseif( vs30 .le. 600.d0 .and. zs_all(ic) .ge. 0.d0) then 
                                   VS = vs30 + 10*dsqrt(Depth);
                                   VP = VS*1.6d0
                                   rho = 1530.d0 + 5.d0*dsqrt(Depth);                                                           
                                   mu = rho * VS**2.d0;
                                   qs = 0.1d0*VS;
                                   gamma = 4.d0*datan(1.d0)/qs;
                                   
                               elseif( vs30 .le. 600.d0 .and. zs_all(ic) .lt. 0.d0) then 
                                   VS = 800 + 10*dsqrt(Depth);
                                   VP  = 2000.d0 + 15.d0*dsqrt(Depth);                           
                                   rho = 1800.d0 +  5.d0*dsqrt(Depth);                                 
                                   mu = rho * VS**2.d0;
                                   qs = 0.1d0*VS;
                                   gamma = 4.d0*datan(1.d0)/qs;
                               endif
                          endif     

                         !left
                         if(dabs(xs(ic) - 415552) .le. 2000.d0) then 
                               VS = 2100 
                               VP = 3500
                               rho = 2200   
                               lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);                           
                               mu = rho * VS**2.d0;
                               qs = 0.1d0*VS;
                               gamma = 4.d0*datan(1.d0)/qs;
!                         !right
                         elseif(dabs(xs(ic) - 484516) .le. 2000.d0) then 
                               VS = 2100 
                               VP = 3500
                               rho = 2200   
                               lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);                           
                               mu = rho * VS**2.d0;
                               qs = 0.1d0*VS;
                               gamma = 4.d0*datan(1.d0)/qs;
!                         !up
                         elseif(dabs(ys(ic) - 4447869) .le. 2000.d0) then 
                               VS = 2100 
                               VP = 3500
                               rho = 2200   
                               lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);                           
                               mu = rho * VS**2.d0;
                               qs = 0.1d0*VS;
                               gamma = 4.d0*datan(1.d0)/qs;
!                         !down
                         elseif(dabs(ys(ic) - 4379160) .le. 2000.d0) then 
                               VS = 2100 
                               VP = 3500
                               rho = 2200   
                               lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);                           
                               mu = rho * VS**2.d0;
                               qs = 0.1d0*VS;
                               gamma = 4.d0*datan(1.d0)/qs;
                         endif                          
                         if (ierr .eq. 1)   write(1000,*) xs(ic),ys(ic),zs(ic), VS, VP       



                                                                

                   elseif (tcase.eq.50) then

                        !-------------------------------------------------------------------
                        ! CASE 50: PLANE-WAVE benchmark -  MULTI NOT honoring
                        !
                        !-------------------------------------------------------------------
                        ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 1st Layer
                
                        if (sub_tag_all(ic).eq.1) then
                             VS = 300;
                             VP  = 600;
                             rho = 2200; 
                             lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                             mu = rho * VS**2.d0;
                             qs = 0.1*VS;
                             gamma = 4.d0*atan(1.d0)/qs;    
                        
                        ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 2nd Layer
                        elseif (sub_tag_all(ic).eq.2) then
                             VS = 2000;
                             VP  = 4000;
                             rho = 2200;
                             lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                             mu = rho * VS**2.d0;
                             qs = 0.1*VS;
                             gamma = 4.d0*atan(1.d0)/qs;    
                                
                        ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 3rd Layer
                        elseif (sub_tag_all(ic).eq.3) then
                             VS = 2000;
                             VP  = 4000;
                             rho = 2200;
                             lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                             mu = rho * VS**2.d0;
                             qs = 0.1*VS;
                             gamma = 4.d0*atan(1.d0)/qs;    
                        
                        ! + MATERIAL INTO THE BEDROCK 
                        elseif (sub_tag_all(ic).eq.4) then
                             VS = 2000;
                             VP = 4000;
                             rho = 2200;
                             lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
                             mu = rho * VS**2.d0;
                             qs = 0.1*VS;
                             gamma = 4.d0*atan(1.d0)/qs;    

                        endif
                        
                        !
                        !-------------------------------------------------------------------
                        elseif (tcase.eq.98) then
                        !-------------------------------------------------------------------
                        ! CASE 98: TEST honoring (only TOPOGRAPHY)

                                VS  = 100        !VS: S velocity in m/s
                                VP  = 200        !VS: S velocity in m/s
                                rho = 2000        !RHO: MASS DENSITY in kg/m^3
                                lambda = rho * (VP**2 - 2*VS**2)
                                mu = rho * VS**2
                                gamma = 0.0d0
                                
                                !-------------------------------------------------------------------

                        elseif (tcase.eq.99) then

                                !-------------------------------------------------------------------
                                ! CASE 99: TEST honoring (TOPOGRAPHY&ALLUVIAL)
                        
                                Depth = zs_elev(ic)    !D: depth in m

                                !if ((Depth .ge. 0.0d0) .and. (zs_all(ic) .ge. 0.0d0)) then

		                    VS  = 300.d0        !VS: S velocity in m/s
		                    VP  = 600.d0        !VP: P velocity in m/s 
		                    rho = 1800.d0        !rho_el: MASS DENSITY in kg/m^3
		                    lambda = rho * (VP**2 - 2*VS**2)
		                    mu = rho * VS**2
		                    qp = 60.d0;
		                    qs = 30.d0;
                                    gamma = 4.d0*datan(1.d0)/qs;

                                !endif
                                
                                !
                                !-------------------------------------------------------------------

                        endif

               rho_el(p,q,r) = rho
               lambda_el(p,q,r) = lambda
               mu_el(p,q,r) = mu
               gamma_el(p,q,r) = gamma




               
               
              
            enddo
         enddo
      enddo
     
      if (ierr .eq. 1)  close(1000)
!      if (ierr .eq. 1)  close(1001)
!      if (ierr .eq. 1)  close(1002)
!      if (ierr .eq. 1)  close(1003)
!      if (ierr .eq. 1)  close(1004)
!      if (ierr .eq. 1)  close(1005)
!      if (ierr .eq. 1)  close(1006)
!      if (ierr .eq. 1)  close(1007)
!      if (ierr .eq. 1)  close(1008)
!      if (ierr .eq. 1)  close(1009)
!      if (ierr .eq. 1)  close(1010)

     if (damping_type .eq. 2) then
       qs = 0; qp = 0;
       vs_all = 0.d0; vp_all=0.d0;
       do r = 1,nn
         do q = 1,nn
            do p = 1,nn
               vs_all = vs_all + dsqrt(mu_el(p,q,r)/rho_el(p,q,r))
               vp_all = vp_all + dsqrt( (lambda_el(p,q,r) + 2.d0*mu_el(p,q,r))/rho_el(p,q,r) );
            enddo
         enddo
       enddo
       
       qs = 0.1d0*vs_all/nn**3;     
       qp = 0.1d0*vp_all/nn**3;

      endif

     
     
     
     


      return
      
      end subroutine MAKE_ELTENSOR_FOR_CASES

