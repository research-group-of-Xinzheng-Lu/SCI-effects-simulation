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

!> @brief Generates seismic triangular faults.
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
!> @param[in] num_ns  number of seismic nodes
!> @param[in] i_sism  column index for filling sour_ns
!> @param[in] nl_sism  number of seismic loads
!> @param[in] nn_loc number of local nodes
!> @param[in] nnod number of local nodes
!> @param[in] loc_n_num local node numeration
!> @param[in] max_num max number of nodes for the specific triangular fault 
!> @param[out] sour_ns info about seismic nodes 
!> @param[out] dist_sour_ns distance of seismic nodes in sour_ns

      subroutine READ_SISM(Xipo,Yipo,Zipo, &                
                                 X1,Y1,Z1, &                        
                                 X2,Y2,Z2, &                        
                                 X3,Y3,Z3, &                        
                                 nnod,xs,ys,zs, &                
                                 num_ns,sour_ns,i_sism,&        
                                 dist_sour_ns,nl_sism,&                
                                 max_num_ns,loc_n_num, nn_loc)                


      implicit none
      
      integer*4 :: i_sism,isn,i,nn_loc
      integer*4 :: nnod,node_sism,num_ns,nl_sism,max_num_ns
      
      integer*4, dimension(nn_loc) :: loc_n_num(nn_loc)

      integer*4, dimension(max_num_ns,nl_sism) :: sour_ns

      real*8 :: Xipo,Yipo,Zipo,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,ux,uy,uz,vx,vy,vz,tol
      real*8 :: p1x,p1y,p1z,p2x,p2y,p2z 
      real*8 :: EPSILON, TWOPI, anglesum, costheta,m1,m2

      real*8, dimension(nn_loc) :: xs,ys,zs
      real*8, dimension(4) :: X,Y,Z

      real*8, dimension(max_num_ns,nl_sism) :: dist_sour_ns

      node_sism = 0

      tol = 6.28
      EPSILON = 1.d-8
      TWOPI = 6.283185307179586476925287
                    
      X(1) = X1
      Y(1) = Y1
      Z(1) = Z1

      X(2) = X2
      Y(2) = Y2
      Z(2) = Z2 
      
      X(3) = X3
      Y(3) = Y3
      Z(3) = Z3

      X(4) = X1
      Y(4) = Y1
      Z(4) = Z1


      do isn = 1,nnod

          anglesum = 0.0d0
          do i = 1,3
                   p1x = X(i) - xs(isn)
                  p1y = Y(i) - ys(isn)
                  p1z = Z(i) - zs(isn)
                  p2x = X(i+1) - xs(isn)
                  p2y = Y(i+1) - ys(isn)
                  p2z = Z(i+1) - zs(isn)

                  m1 = dsqrt(p1x*p1x + p1y*p1y + p1z*p1z)
                  m2 = dsqrt(p2x*p2x + p2y*p2y + p2z*p2z)
      
                  if ((m1*m2).le.EPSILON) then
                        anglesum = TWOPI
                  else
                        costheta = (p1x*p2x + p1y*p2y + p1z*p2z)/(m1*m2)
                        anglesum = anglesum + dacos(costheta)
                  endif
          enddo

          if (anglesum.ge.tol) then

              node_sism = node_sism + 1
                sour_ns(node_sism,i_sism) = loc_n_num(isn) 
            dist_sour_ns(node_sism,i_sism) = sqrt((Xipo - xs(isn))**2 +(Yipo - ys(isn))**2 +(Zipo - zs(isn))**2)

          endif

      enddo

     
      return
      end subroutine READ_SISM

