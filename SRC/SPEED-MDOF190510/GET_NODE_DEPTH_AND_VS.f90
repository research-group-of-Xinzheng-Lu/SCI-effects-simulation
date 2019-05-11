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

!> @brief Computes elevation from topography (XYZ.out), vs30 velocity
!! and depth of sediments according to the input file VS_RS.out 
!> @note See MAKE_NOT_HONORING.f90 
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] loc_n_num local numeration vector
!> @param[in] nn_elev number nodes in the triangular grid
!> @param[in] x_elev elevation values of local nodes
!> @param[in] y_elev elevation values of local nodes
!> @param[in] z_elev elevation values of local nodes
!> @param[in] nn_elem number of triangular elements
!> @param[in] node1_elem index triangle vertex 
!> @param[in] node2_elem index triangle vertex
!> @param[in] node3_elem index triangle vertex
!> @param[in] cs_nnx_loc  length cs_loc
!> @param[in] cs_loc local connectivity vector
!> @param[in] nm number of materials
!> @param[in] tm labels for material vector
!> @param[in] sd polynomial degree vector 
!> @param[in] nn_s number of local nodes
!> @param[in] xx_s vertex x- coordinate  of local nodes
!> @param[in] yy_s vertex y- coordinate  of local nodes
!> @param[in] zz_s vertex z- coordinate  of local nodes
!> @param[in] tagmat specific material tag given in CASE option
!> @param[in] max_es  max topography spacing
!> @param[in] tol  tolerance given in CASE option
!> @param[in] zz_alluvial elevation of the nodes from alluvial soil
!> @param[out] zz_elevation elevation of the nodes from complex topography
!> @param[out] vs vs30 velocity for the nodes
!> @param[out] thickness thickenss of sediments along z-direction


      subroutine GET_NODE_DEPTH_AND_VS(loc_n_num, nn_elev, nn_elem, &
                                            xx_elev, yy_elev, zz_elev, vs_elev, sediments, &
                                            node1_elem, node2_elem, node3_elem, &
                                            cs_nnz_loc, cs_loc, nm, tm, sd, &
                                            nn_s, xx_s, yy_s, zz_s, &
                                            zz_elevation, zz_alluvial, vs, thickness, &
                                            tagmat, max_es,tol)      

      implicit none
      
      integer*4 :: im,ie,i,j,k,nn,ip,isn,ic
      integer*4 :: nn_elev, nn_elem, cs_nnz_loc, nm, ne, nn_s
      integer*4 :: h, h_sel
      integer*4 :: tagmat

      integer*4, dimension(nn_s) :: loc_n_num                                                         
      integer*4, dimension(nn_elem) :: node1_elem,node2_elem,node3_elem
      integer*4, dimension(nm) :: tm
      integer*4, dimension(nm) :: sd
      integer*4, dimension(0:cs_nnz_loc) :: cs_loc

      real*8 :: X1,Y1,Z1                                
      real*8 :: X2,Y2,Z2                                
      real*8 :: X3,Y3,Z3                                
      real*8 :: ux,uy,uz,vx,vy,vz                
      real*8 :: a,b,c                                        
      real*8 :: max_es                                        
      real*8 :: zz_interp                                
      real*8 :: v0x,v0y,v1x,v1y,v2x,v2y                        
      real*8 :: dot00,dot01,dot02,dot11,dot12        
      real*8 :: invDenom,u,v                                        
      real*8 :: d2min 
      real*8 :: zz_elev_min 
      real*8 :: dx,dy,dz,tol, dist, dist_min

      real*8, dimension(:), allocatable :: ct,ww
      real*8, dimension(nn_elev) :: xx_elev,yy_elev,zz_elev
      real*8, dimension(nn_elem) :: vs_elev,sediments 
      real*8, dimension(nn_s) :: xx_s,yy_s,zz_s
      real*8, dimension(nn_s) :: zz_elevation
      real*8, dimension(nn_s) :: zz_alluvial
      real*8, dimension(nn_s) :: vs, thickness
      
      real*8, dimension(:,:), allocatable :: dd

      d2min = (5 * max_es)**2 
          
      zz_elev_min = zz_elev(1)
      do i = 1,nn_elev
         if (zz_elev(i).lt.zz_elev_min) then
                 zz_elev_min = zz_elev(i)
          endif
      enddo
      !write(*,*) d2min, max_es

      
      ne = cs_loc(0) - 1
      
      do ie = 1,ne
         
         im = cs_loc(cs_loc(ie -1) +0)
         
         if (im .eq. tagmat) then
         
             nn = sd(im) +1
             allocate(ct(nn),ww(nn),dd(nn,nn))
             call MAKE_LGL_NW(nn,ct,ww,dd)  

                            
               do k = 1,nn
                  do j = 1,nn
                     do i = 1,nn
                        
                        ip = nn*nn*(k -1) +nn*(j -1) +i
                        isn = cs_loc(cs_loc(ie -1) + ip)
                        ic = isn
                        
                        dist_min = 1.d+30;
                        h_sel = 0;
                        do h = 1,nn_elem
                            
                            
                            X1 = xx_elev(node1_elem(h)) 
                            Y1 = yy_elev(node1_elem(h)) 
                            Z1 = zz_elev(node1_elem(h)) 
                            
                            dist = (X1 - xx_s(ic))**2 + (Y1 - yy_s(ic))**2
                            if (dist .le. dist_min) then
                               h_sel = h;
                               dist_min = dist;
                            endif          
                        enddo
                                     
                            !write(*,*) h_sel 
                            !read(*,*)                            
                            X1 = xx_elev(node1_elem(h_sel)) 
                            Y1 = yy_elev(node1_elem(h_sel)) 
                            Z1 = zz_elev(node1_elem(h_sel)) 
 
                            X2 = xx_elev(node2_elem(h_sel)) 
                            Y2 = yy_elev(node2_elem(h_sel)) 
                            Z2 = zz_elev(node2_elem(h_sel)) 

                            X3 = xx_elev(node3_elem(h_sel)) 
                            Y3 = yy_elev(node3_elem(h_sel)) 
                            Z3 = zz_elev(node3_elem(h_sel)) 

                           
                            zz_interp = (Z1 + Z2 + Z3)/3.d0
                            zz_elevation(ic) = ( zz_interp - zz_s(ic) )
                            vs(ic) = vs_elev(h_sel)
                            thickness(ic) = sediments(h_sel)
                                                        
                            if(zz_elevation(ic) .le. tol) then
                               zz_elevation(ic) = 0.0d0
                               vs(ic) = vs_elev(h_sel)
                               thickness(ic) = sediments(h_sel)
                            endif
                                                        
                                                                                


              enddo !do k = 1,nn
            enddo !do j = 1,nn
          enddo !do i = 1,nn

          deallocate(ct,ww,dd)
                    
        endif !if im.eq.tagmat
      enddo

      
      return
      
      end subroutine GET_NODE_DEPTH_AND_VS

