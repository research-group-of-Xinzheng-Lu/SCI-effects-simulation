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

!> @brief Takes values from nodes array. 
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[out] nodes costant values for the bilinear map to read
!> @param[out] c_alfa11 costant values for the bilinear map
!> @param[out] c_alfa12 costant values for the bilinear map
!> @param[out] c_alfa13 costant values for the bilinear map
!> @param[out] c_alfa21 costant values for the bilinear map
!> @param[out] c_alfa22 costant values for the bilinear map
!> @param[out] c_alfa23 costant values for the bilinear map
!> @param[out] c_alfa31 costant values for the bilinear map
!> @param[out] c_alfa32 costant values for the bilinear map
!> @param[out] c_alfa33 costant values for the bilinear map
!> @param[out] c_beta11 costant values for the bilinear map
!> @param[out] c_beta12 costant values for the bilinear map
!> @param[out] c_beta13 costant values for the bilinear map
!> @param[out] c_beta21 costant values for the bilinear map
!> @param[out] c_beta22 costant values for the bilinear map
!> @param[out] c_beta23 costant values for the bilinear map
!> @param[out] c_beta31 costant values for the bilinear map
!> @param[out] c_beta32 costant values for the bilinear map
!> @param[out] c_beta33 costant values for the bilinear map
!> @param[out] c_gamma1 costant values for the bilinear map
!> @param[out] c_gamma2 costant values for the bilinear map
!> @param[out] c_gamma3 costant values for the bilinear map
!> @param[out] c_delta1 costant values for the bilinear map
!> @param[out] c_delta2 costant values for the bilinear map
!> @param[out] c_delta3 costant values for the bilinear map

     subroutine MAKE_BILINEAR_MAP(nodes, c_alfa11, c_alfa12, c_alfa13, &
                                         c_alfa21, c_alfa22, c_alfa23, &
                                         c_alfa31, c_alfa32, c_alfa33, &
                                         c_beta11, c_beta12, c_beta13, & 
                                         c_beta21, c_beta22, c_beta23, & 
                                         c_beta31, c_beta32, c_beta33, &
                                         c_gamma1, c_gamma2, c_gamma3, &
                                         c_delta1, c_delta2, c_delta3)

     implicit none
     
     real*8 :: x1,x2,x3,x4,x5,x6,x7,x8
     real*8 :: y1,y2,y3,y4,y5,y6,y7,y8
     real*8 :: z1,z2,z3,z4,z5,z6,z7,z8
     real*8, intent(out) :: c_alfa11, c_alfa12, c_alfa13, c_alfa21, c_alfa22, c_alfa23, c_alfa31, c_alfa32, c_alfa33
     real*8, intent(out) :: c_beta11, c_beta12, c_beta13, c_beta21, c_beta22, c_beta23, c_beta31, c_beta32, c_beta33
     real*8, intent(out) :: c_gamma1, c_gamma2, c_gamma3, c_delta1, c_delta2, c_delta3

     real*8, dimension(24) :: nodes
     

     c_alfa11 = nodes(1)
     c_alfa12 = nodes(2)
     c_alfa13 = nodes(3)
     
     c_alfa21 = nodes(4)
     c_alfa22 = nodes(5)
     c_alfa23 = nodes(6)

     c_alfa31 = nodes(7)
     c_alfa32 = nodes(8)
     c_alfa33 = nodes(9)

     c_beta11 = nodes(10)
     c_beta12 = nodes(11)
     c_beta13 = nodes(12)

     c_beta21 = nodes(13)
     c_beta22 = nodes(14)
     c_beta23 = nodes(15)

     c_beta31 = nodes(16)
     c_beta32 = nodes(17)
     c_beta33 = nodes(18)

     c_gamma1 = nodes(19)
     c_gamma2 = nodes(20)
     c_gamma3 = nodes(21)

     c_delta1 = nodes(22)
     c_delta2 = nodes(23)
     c_delta3 = nodes(24)

 
    end subroutine MAKE_BILINEAR_MAP
