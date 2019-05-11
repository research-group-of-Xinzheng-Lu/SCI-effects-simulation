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

!> @brief Computes area of a face (quad).
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] x1 vertex x-coordinate defining a face of the element ielem
!> @param[in] x2 vertex x-coordinate defining a face of the element ielem
!> @param[in] x3 vertex x-coordinate defining a face of the element ielem
!> @param[in] x4 vertex x-coordinate defining a face of the element ielem
!> @param[in] y1 vertex y-coordinate defining a face of the element ielem
!> @param[in] y2 vertex y-coordinate defining a face of the element ielem
!> @param[in] y3 vertex y-coordinate defining a face of the element ielem
!> @param[in] y4 vertex y-coordinate defining a face of the element ielem
!> @param[in] z1 vertex z-coordinate defining a face of the element ielem
!> @param[in] z2 vertex z-coordinate defining a face of the element ielem
!> @param[in] z3 vertex z-coordinate defining a face of the element ielem
!> @param[in] z4 vertex z-coordinate defining a face of the element ielem
!> @param[in] index of the element
!> @param[out] surf  area of a face of the hex ielem 

     subroutine GET_AREA_FACE(x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,surf,ielem)
     
     implicit none

     integer*4, intent(in) :: ielem

     real*8 :: a,b,c,vx,vy,vz,wx,wy,wz,xbar,ybar,zbar
     real*8, intent(in) :: x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
     real*8, intent(out) :: surf
     
     surf = 0.d0
     
     xbar = (x1 + x2 + x3 + x4)/4.d0
     ybar = (y1 + y2 + y3 + y4)/4.d0
     zbar = (z1 + z2 + z3 + z4)/4.d0
     
     vx = x4 - x1
     vy = y4 - y1 
     vz = z4 - z1
     
     wx = xbar - x1
     wy = ybar - y1
     wz = zbar - z1     
     
                 
     a =  vy*wz - vz*wy
     b = -vx*wz + vz*wx
     c =  vx*wy - vy*wx
     

                        
     surf = surf + 0.5d0*dsqrt(a**2.d0 + b**2.d0 + c**2.d0)  
     
     vx = x2 - x1
     vy = y2 - y1 
     vz = z2 - z1
     
     wx = xbar - x1
     wy = ybar - y1
     wz = zbar - z1     
     
                 
     a =  vy*wz - vz*wy
     b = -vx*wz + vz*wx
     c =  vx*wy - vy*wx
     
     surf = surf + 0.5d0*dsqrt(a**2.d0 + b**2.d0 + c**2.d0)  
     
     vx = x2 - x3
     vy = y2 - y3 
     vz = z2 - z3
     
     wx = xbar - x3
     wy = ybar - y3
     wz = zbar - z3     
     
                 
     a =  vy*wz - vz*wy
     b = -vx*wz + vz*wx
     c =  vx*wy - vy*wx
     
     surf = surf + 0.5d0*dsqrt(a**2.d0 + b**2.d0 + c**2.d0)  
     
     vx = x4 - x3
     vy = y4 - y3 
     vz = z4 - z3
     
     wx = xbar - x3
     wy = ybar - y3
     wz = zbar - z3     
     
                 
     a =  vy*wz - vz*wy
     b = -vx*wz + vz*wx
     c =  vx*wy - vy*wx
     
     surf = surf + 0.5d0*dsqrt(a**2.d0 + b**2.d0 + c**2.d0)  
     
     end subroutine GET_AREA_FACE
