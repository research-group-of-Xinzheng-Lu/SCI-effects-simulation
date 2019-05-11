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

!> @brief Computes elevation from topography (XYZ.out). 
!> @note See MAKE_NOT_HONORING.f90 
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] nn_s number of local nodes
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
!> @param[in] xx_s vertex x- coordinate  of local nodes
!> @param[in] yy_s vertex y- coordinate  of local nodes
!> @param[in] zz_s vertex z- coordinate  of local nodes
!> @param[in] nm number of materials
!> @param[in] tm labels for material vector
!> @param[in] sd polynomial degree vector 
!> @param[in] tagmat specific material tag given in CASE option
!> @param[in] max_es  max topography spacing
!> @param[in] tol  tolerance given in CASE option
!> @param[in] zz_alluvial elevation of the nodes from alluvial soil
!> @param[out] zz_elevation elevation of the nodes from complex topography

      subroutine GET_NODE_DEPTH_FROM_CMPLX(loc_n_num, nn_elev, nn_elem, &
                                            xx_elev, yy_elev, zz_elev, &
                                            node1_elem, node2_elem, node3_elem, &
                                            cs_nnz_loc, cs_loc, nm, tm, sd, &
                                            nn_s, xx_s, yy_s, zz_s, &
                                            zz_elevation, zz_alluvial, &
                                            tagmat, max_es,tol)      

      implicit none
      
      integer*4 :: im,ie,i,j,k,nn,ip,isn,ic
      integer*4 :: nn_elev, nn_elem, cs_nnz_loc, nm, ne, nn_s
      integer*4 :: h
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
      real*8 :: dx,dy,dz,tol

      real*8, dimension(:), allocatable :: ct,ww
      real*8, dimension(nn_elev) :: xx_elev,yy_elev,zz_elev
      real*8, dimension(nn_s) :: xx_s,yy_s,zz_s
      real*8, dimension(nn_s) :: zz_elevation
      real*8, dimension(nn_s) :: zz_alluvial
      
      real*8, dimension(:,:), allocatable :: dd

      d2min = (5 * max_es)**2 
          
      zz_elev_min = zz_elev(1)
      do i = 1,nn_elev
         if (zz_elev(i).lt.zz_elev_min) then
                 zz_elev_min = zz_elev(i)
          endif
      enddo
      

      nn = 2
      allocate(ct(nn),ww(nn),dd(nn,nn))
      call MAKE_LGL_NW(nn,ct,ww,dd)
      
      ne = cs_loc(0) - 1
      
      do im = 1,nm
         if ((sd(im) +1).ne.nn) then
            deallocate(ct,ww,dd)
            nn = sd(im) +1
            allocate(ct(nn),ww(nn),dd(nn,nn))
            call MAKE_LGL_NW(nn,ct,ww,dd)
         endif
         
         do ie = 1,ne
            if (cs_loc(cs_loc(ie -1) +0).eq.tagmat) then
              
               do k = 1,nn
                  do j = 1,nn
                     do i = 1,nn
                        
                        ip = nn*nn*(k -1) +nn*(j -1) +i
                        isn = cs_loc(cs_loc(ie -1) + ip)
                        ic = isn
                        


                        if ((zz_elevation(ic).eq.-1.0e+30).and.(zz_alluvial(ic).ge.0.0d0)) then

                                do h = 1,nn_elem

                                        X1 = xx_elev(node1_elem(h)) 
                                        Y1 = yy_elev(node1_elem(h)) 
                                        Z1 = zz_elev(node1_elem(h)) 
                                                                
                                        if (((X1 - xx_s(ic))**2 + (Y1 - yy_s(ic))**2).le.d2min) then 

                                                X2 = xx_elev(node2_elem(h)) 
                                                Y2 = yy_elev(node2_elem(h)) 
                                                Z2 = zz_elev(node2_elem(h)) 

                                                X3 = xx_elev(node3_elem(h)) 
                                                Y3 = yy_elev(node3_elem(h)) 
                                                Z3 = zz_elev(node3_elem(h)) 

                                                !Point in triangle test

                                                ! P = (xx_s(ic) yy_s(ic))
                                                ! A = (X1 Y1)
                                                ! B = (X2 Y2)
                                                ! C = (X3 Y3)

                                                ! Compute vectors 
                                                ! v0 = C - A
                                                v0x=(X3 - X1) 
                                                v0y=(Y3 - Y1) 

                                                !v0x=(X3 - X1)/sqrt((X3 - X1)**2+(Y3 - Y1)**2) 
                                                !v0y=(Y3 - Y1)/sqrt((X3 - X1)**2+(Y3 - Y1)**2) 
                                                                                        
                                                ! v1 = B - A
                                                v1x=(X2 - X1) 
                                                v1y=(Y2 - Y1) 

                                                !v1x=(X2 - X1)/sqrt((X2 - X1)**2+(Y2 - Y1)**2) 
                                                !v1y=(Y2 - Y1)/sqrt((X2 - X1)**2+(Y2 - Y1)**2) 

                                                ! v2 = P - A
                                                v2x=(xx_s(ic) - X1) 
                                                v2y=(yy_s(ic) - Y1) 

                                                !v2x=(xx_s(ic) - X1)/sqrt((xx_s(ic) - X1)**2+(yy_s(ic) - Y1)**2) 
                                                !v2y=(yy_s(ic) - Y1)/sqrt((xx_s(ic) - X1)**2+(yy_s(ic) - Y1)**2) 
                                                                                        
                                                ! Compute dot products
                                                ! [u].[v] = ux * vx + uy * vy
                                                ! dot([u],[v])

                                                !dot00 = dot(v0, v0)
                                                dot00 = v0x * v0x + v0y * v0y

                                                !dot01 = dot(v0, v1)
                                                dot01 = v0x * v1x + v0y * v1y

                                                !dot02 = dot(v0, v2)
                                                dot02 = v0x * v2x + v0y * v2y

                                                !dot11 = dot(v1, v1)
                                                dot11 = v1x * v1x + v1y * v1y

                                                !dot12 = dot(v1, v2)
                                                dot12 = v1x * v2x + v1y * v2y

                                                ! Compute barycentric coordinates
                                                invDenom = 1 / (dot00 * dot11 - dot01 * dot01)
                                                u = (dot11 * dot02 - dot01 * dot12) * invDenom
                                                v = (dot00 * dot12 - dot01 * dot02) * invDenom

                                                !Point in triangle test 
                                                !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                                                if ( (u.ge.0.0d0).and.(v.ge.0.0d0).and.((u + v).le.1.0d0) ) then
                                                                                        
                                                        ! Build up the plane passing through the points P1, P2 and P3

                                                        ux=(X1-X2)/sqrt((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2) 
                                                        uy=(Y1-Y2)/sqrt((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2) 
                                                        uz=(Z1-Z2)/sqrt((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2) 
                                                        vx=(X3-X2)/sqrt((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2) 
                                                        vy=(Y3-Y2)/sqrt((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2) 
                                                        vz=(Z3-Z2)/sqrt((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2) 
      
                                                        a = uy * vz - uz * vy 
                                                        b = uz * vx - ux * vz 
                                                        c = ux * vy - uy * vx 
                                                                                
                                                        zz_interp = -a/c * (xx_s(ic)-X1) -b/c * (yy_s(ic)-Y1) + Z1
                                                        zz_elevation(ic) = ( zz_interp - zz_s(ic) )
                                                        
                                                        if (abs(zz_elevation(ic)).lt.tol) then
                                                                zz_elevation(ic) = 0.0d0
                                                        endif
                                                endif !if ( (u.ge.0.0d0).and.(v.ge.0.0d0).and.((u + v).le.1.0d0) ) then
                                                                                
                                                if ( (u.ge.0.0d0).and.(v.ge.0.0d0).and.((u + v).le.1.0d0) ) exit

                                endif !if (((X1 - xx_s(ic))**2 + (Y1 - yy_s(ic))**2).le.d2min) then 

                        enddo !do h = 1,nn_elem
                endif !if ((zz_elevation(ic).eq.-1.0e+30).and.(zz_alluvial(ic).ge.0.0d0)) then
              enddo !do k = 1,nn
            enddo !do j = 1,nn
          enddo !do i = 1,nn

               
            endif
         enddo
      enddo
      
      return
      
      end subroutine GET_NODE_DEPTH_FROM_CMPLX

