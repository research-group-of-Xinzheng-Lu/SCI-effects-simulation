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


!> @brief Makes normal vector of a given surface.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] ind  face of the hex
!> @param[in] xs1 x-coordinate of the 1-node defining the face 
!> @param[in] xs2 x-coordinate of the 2-node defining the face 
!> @param[in] xs3 x-coordinate of the 3-node defining the face 
!> @param[in] xs4 x-coordinate of the 4-node defining the face 
!> @param[in] ys1 y-coordinate of the 1-node defining the face 
!> @param[in] ys2 y-coordinate of the 2-node defining the face 
!> @param[in] ys3 y-coordinate of the 3-node defining the face 
!> @param[in] ys4 y-coordinate of the 4-node defining the face 
!> @param[in] zs1 z-coordinate of the 1-node defining the face 
!> @param[in] zs2 z-coordinate of the 2-node defining the face 
!> @param[in] zs3 z-coordinate of the 3-node defining the face 
!> @param[in] zs4 z-coordinate of the 4-node defining the face 
!> @param[in] par  control parameter (dummy)
!> @param[in] arr  control parameter (dummy)
!> @param[out] nx x-component of the normal vector
!> @param[out] ny y-component of the normal vector
!> @param[out] nz z-component of the normal vector


    subroutine MAKE_NORMAL(ind, xs1, xs2, xs3, xs4, &
                       ys1, ys2, ys3, ys4, &
                       zs1, zs2, zs3, zs4, &
                       nx,ny,nz, par, arr)

        implicit none
        
        integer*4, intent(in) :: ind, par, arr                                    

        real*8 :: a1,a2,a3,b1,b2,b3,aplane,bplane,cplane
        real*8 :: norm1, norm2, norm3, norm4
        real*8 :: nx1, nx2, nx3, nx4
        real*8 :: ny1, ny2, ny3, ny4
        real*8 :: nz1, nz2, nz3, nz4

        real*8, intent(in)  :: xs1, xs2, xs3, xs4, ys1, ys2, ys3, ys4, zs1, zs2, zs3, zs4                                   
        real*8, intent(out) :: nx,ny,nz

         select case (ind)

         case(1,3,5) 
            a1 = xs4-xs1  
            a2 = ys4-ys1
            a3 = zs4-zs1
            b1 = xs2-xs1  
            b2 = ys2-ys1  
            b3 = zs2-zs1                  
            aplane =  a2*b3 - a3*b2
            bplane = -a1*b3 + a3*b1
            cplane =  a1*b2 - a2*b1
                 
            norm1 =  aplane**2.d0 + bplane**2.d0 + cplane**2.d0  
                nx1 = aplane/dsqrt(norm1)
            ny1 = bplane/dsqrt(norm1)
            nz1 = cplane/dsqrt(norm1)
            
            
            a1 = xs1-xs2  
            a2 = ys1-ys2
            a3 = zs1-zs2
            b1 = xs3-xs2  
            b2 = ys3-ys2  
            b3 = zs3-zs2                  
            aplane =  a2*b3 - a3*b2
            bplane = -a1*b3 + a3*b1
            cplane =  a1*b2 - a2*b1
                 
            norm2 =  aplane**2.d0 + bplane**2.d0 + cplane**2.d0  
                nx2 = aplane/dsqrt(norm2)
            ny2 = bplane/dsqrt(norm2)
            nz2 = cplane/dsqrt(norm2)
            
            
            a1 = xs2-xs3  
            a2 = ys2-ys3
            a3 = zs2-zs3
            b1 = xs4-xs3  
            b2 = ys4-ys3  
            b3 = zs4-zs3                  
            aplane =  a2*b3 - a3*b2
            bplane = -a1*b3 + a3*b1
            cplane =  a1*b2 - a2*b1
                 
            norm3 =  aplane**2.d0 + bplane**2.d0 + cplane**2.d0  
                nx3 = aplane/dsqrt(norm3)
            ny3 = bplane/dsqrt(norm3)
            nz3 = cplane/dsqrt(norm3)
            
            
            a1 = xs3-xs4  
            a2 = ys3-ys4
            a3 = zs3-zs4
            b1 = xs1-xs4  
            b2 = ys1-ys4  
            b3 = zs1-zs4                  
            aplane =  a2*b3 - a3*b2
            bplane = -a1*b3 + a3*b1
            cplane =  a1*b2 - a2*b1
                 
            norm4 =  aplane**2.d0 + bplane**2.d0 + cplane**2.d0  
                nx4 = aplane/dsqrt(norm4)
            ny4 = bplane/dsqrt(norm4)
            nz4 = cplane/dsqrt(norm4)
            
            nx = (nx1 + nx2 + nx3 + nx4)/4.d0
            ny = (ny1 + ny2 + ny3 + ny4)/4.d0
            nz = (nz1 + nz2 + nz3 + nz4)/4.d0            
            
 
         case(2,4,6)
         
           a1 = xs2-xs1  
           a2 = ys2-ys1
           a3 = zs2-zs1
                 b1 = xs4-xs1  
           b2 = ys4-ys1  
           b3 = zs4-zs1                   
           aplane =  a2*b3 - a3*b2
           bplane = -a1*b3 + a3*b1
           cplane =  a1*b2 - a2*b1
                      
           norm1 =  aplane**2.d0 + bplane**2.d0 + cplane**2.d0  
           nx1 = aplane/dsqrt(norm1)
           ny1 = bplane/dsqrt(norm1)
           nz1 = cplane/dsqrt(norm1)
           
           
           a1 = xs3-xs2  
           a2 = ys3-ys2
           a3 = zs3-zs2
                 b1 = xs1-xs2  
           b2 = ys1-ys2  
           b3 = zs1-zs2                   
           aplane =  a2*b3 - a3*b2
           bplane = -a1*b3 + a3*b1
           cplane =  a1*b2 - a2*b1
                      
           norm2 =  aplane**2.d0 + bplane**2.d0 + cplane**2.d0  
           nx2 = aplane/dsqrt(norm2)
           ny2 = bplane/dsqrt(norm2)
           nz2 = cplane/dsqrt(norm2)
           
           
           a1 = xs4-xs3  
           a2 = ys4-ys3
           a3 = zs4-zs3
                 b1 = xs2-xs3  
           b2 = ys2-ys3  
           b3 = zs2-zs3                   
           aplane =  a2*b3 - a3*b2
           bplane = -a1*b3 + a3*b1
           cplane =  a1*b2 - a2*b1
                      
           norm3 =  aplane**2.d0 + bplane**2.d0 + cplane**2.d0  
           nx3 = aplane/dsqrt(norm3)
           ny3 = bplane/dsqrt(norm3)
           nz3 = cplane/dsqrt(norm3)
           
           
           a1 = xs1-xs4  
           a2 = ys1-ys4
           a3 = zs1-zs4
                 b1 = xs3-xs4  
           b2 = ys3-ys4  
           b3 = zs3-zs4                   
           aplane =  a2*b3 - a3*b2
           bplane = -a1*b3 + a3*b1
           cplane =  a1*b2 - a2*b1
                      
           norm4 =  aplane**2.d0 + bplane**2.d0 + cplane**2.d0  
           nx4 = aplane/dsqrt(norm4)
           ny4 = bplane/dsqrt(norm4)
           nz4 = cplane/dsqrt(norm4)


            nx = (nx1 + nx2 + nx3 + nx4)/4.d0
            ny = (ny1 + ny2 + ny3 + ny4)/4.d0
            nz = (nz1 + nz2 + nz3 + nz4)/4.d0          
           
           
         end select                                    
                                            
                        
                                    
end subroutine MAKE_NORMAL
