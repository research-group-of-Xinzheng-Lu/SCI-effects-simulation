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

!> @brief Computes local number of nodes where seismic source is imposed.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] Xipo x-coordinate of Hypocenter
!> @param[in] Yipo y-coordinate of Hypocenter
!> @param[in] Zipo z-coordinate of Hypocenter
!> @param[in] X1 vertex coordinate (x) of the triangular fault
!> @param[in] Y1 vertex coordinate (y) of the triangular fault
!> @param[in] Z1 vertex coordinate (z) of the triangular fault
!> @param[in] X2 vertex coordinate (x) of the triangular fault
!> @param[in] Y2 vertex coordinate (y) of the triangular fault
!> @param[in] Z2 vertex coordinate (z) of the triangular fault
!> @param[in] X3 vertex coordinate (x) of the triangular fault
!> @param[in] Y3 vertex coordinate (y) of the triangular fault
!> @param[in] Z3 vertex coordinate (z) of the triangular fault
!> @param[in] nnod number of local nodes
!> @param[in] xs x-coordinate of spectral nodes
!> @param[in] ys y-coordinate of spectral nodes
!> @param[in] zs z-coordinate of spectral nodes
!> @param[in] nn_loc number of local nodes
!> @param[in] loc_n_num numeration of local nodes
!> @param[in] mpi_id id of MPI process
!> @param[in] ind index dummy
!> @param[out] node_sism  local number of seismic source


      subroutine GET_DIME_SISM(Xipo,Yipo,Zipo, &
                                 X1,Y1,Z1, &       
                                 X2,Y2,Z2, &       
                                 X3,Y3,Z3, &       
                                 nnod,xs,ys,zs, &       
                                 node_sism, nn_loc, mpi_id, ind, loc_n_num)				       
         
      implicit none
      
      integer*4 :: isn,i,nn_loc,mpi_id,ind
      integer*4 :: nnod,node_sism
      
      integer*4, dimension(nn_loc) :: loc_n_num

      real*8 :: Xipo,Yipo,Zipo,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,ux,uy,uz,vx,vy,vz,tol
      real*8 :: p1x,p1y,p1z,p2x,p2y,p2z
      real*8 :: EPSILON, TWOPI, anglesum, costheta,m1,m2

      real*8, dimension(4) :: X,Y,Z
      real*8, dimension(nn_loc) :: xs,ys,zs

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
        endif

      enddo


      return
      end subroutine GET_DIME_SISM

