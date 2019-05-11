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

!> @brief Generates explosive triangular faults.
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] Xipo Hipocenter x-coordinate
!> @param[in] Yipo Hipocenter y-coordinate
!> @param[in] Zipo Hipocenter z-coordinate 
!> @param[in] X1 x-coordinate of the triangular fault 1-node
!> @param[in] Y1 y-coordinate of the triangular fault 1-node
!> @param[in] Z1 z-coordinate of the triangular fault 1-node
!> @param[in] X2 x-coordinate of the triangular fault 2-node
!> @param[in] Y2 y-coordinate of the triangular fault 2-node
!> @param[in] Z2 z-coordinate of the triangular fault 2-node
!> @param[in] X3 x-coordinate of the triangular fault 3-node
!> @param[in] Y3 y-coordinate of the triangular fault 3-node
!> @param[in] Z3 z-coordinate of the triangular fault 3-node
!> @param[in] xs x-coord local spectral nodes
!> @param[in] ys y-coord local spectral nodes
!> @param[in] zs z-coord local spectral nodes
!> @param[in] num_ne  number of explosive nodes
!> @param[in] i  column index for filling sour_ne
!> @param[in] nl_expl  number of explosive loads
!> @param[in] nn_loc number of local nodes
!> @param[in] nnod number of local nodes
!> @param[in] loc_n_num local node numeration
!> @param[in] max_num_ne max number of nodes for the specific triangular fault 
!> @param[out] sour_ne info about explosive nodes 
!> @param[out] dist_sour_ne distance of explosive nodes in sour_ne

      subroutine READ_EXPL(Xipo,Yipo,Zipo, &                       
                                 X1,Y1,Z1, &                               
                                 X2,Y2,Z2, &                               
                                 X3,Y3,Z3, &                               
                                 nnod,xs,ys,zs, &                       
                                 num_ne,sour_ne,i,&                       
                                 dist_sour_ne,nl_expl,&               
                                 max_num_ne, loc_n_num, nn_loc)        
                                                        

      implicit none
      
      integer*4 :: nnod,node_expl,num_ne,nl_expl,max_num_ne, nn_loc
      integer*4 :: i,isn
      
      integer*4, dimension(nn_loc) :: loc_n_num(nn_loc)

      integer*4, dimension(max_num_ne,nl_expl) :: sour_ne

      real*8 :: Xipo,Yipo,Zipo,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,ux,uy,uz,vx,vy,vz,tol
      real*8 :: Xmax,Xmin,Ymax,Ymin,Zmax,Zmin,a,b,c
      real*8 :: den

      real*8, dimension(nn_loc) :: xs,ys,zs

      real*8, dimension(max_num_ne,nl_expl) :: dist_sour_ne

      node_expl = 0
      tol = 5.0d0
      
      ux=(X1-X2)/sqrt((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
      uy=(Y1-Y2)/sqrt((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
      uz=(Z1-Z2)/sqrt((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
      vx=(X3-X2)/sqrt((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2)
      vy=(Y3-Y2)/sqrt((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2)
      vz=(Z3-Z2)/sqrt((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2)
      
      a = uy * vz - uz * vy 
      b = uz * vx - ux * vz
      c = ux * vy - uy * vx
          den = sqrt(a**2 + b**2 + c**2)

         !Check on the relative position of fault end points
         !In this preliminary part we check if P1 point is the lower-left fault point or not
         !If this check would fail the two points are swapped
  
         
         Xmax=max(X1,X2,X3)+tol
         Xmin=min(X1,X2,X3)-tol
         Ymax=max(Y1,Y2,Y3)+tol
         Ymin=min(Y1,Y2,Y3)-tol
         Zmax=max(Z1,Z2,Z3)+tol
         Zmin=min(Z1,Z2,Z3)-tol
         

         do isn = 1,nnod

            if ( dabs( a*(xs(isn)-X1) + b*(ys(isn)-Y1) + c*(zs(isn)-Z1)/den ).le.tol) then
                if ((xs(isn).ge.Xmin).and.(xs(isn).le.Xmax)) then
                      if ((ys(isn).ge.Ymin).and.(ys(isn).le.Ymax)) then
                        if ((zs(isn).ge.Zmin).and.(zs(isn).le.Zmax)) then

                            node_expl = node_expl + 1   
                            sour_ne(node_expl,i) = loc_n_num(isn)
                            dist_sour_ne(node_expl,i) = sqrt((Xipo - xs(isn))**2 + (Yipo - ys(isn))**2 + (Zipo - zs(isn))**2)
                      
                        endif
                      endif
                 endif
            endif
             
         enddo

     
      return
      end subroutine READ_EXPL

