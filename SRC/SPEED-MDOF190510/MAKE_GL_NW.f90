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

!> @brief Makes Gauss Legendre nodes and weights 
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] ngp  number of Gauss points
!> @param[out] xabsc Gauss points
!> @param[out] weig weights for the Gauss quadrature rule

 
      subroutine  MAKE_GL_NW(ngp, xabsc, weig)

      implicit none

      INTEGER, parameter :: dbp = SELECTED_REAL_KIND (15,307)
      REAL(dbp) :: newv
      REAL(dbp) :: EPS, M_PI
      PARAMETER (EPS=3.0d-15)               !EPS is the relative precision
      PARAMETER (M_PI=3.141592654d0) 
      
      INTEGER*4  i, j, m
      REAL(dbp)  p1, p2, p3, pp, z, z1
      INTEGER*4, INTENT(IN) :: ngp            ! # of Gauss Points
      REAL(dbp), INTENT(OUT) :: xabsc(ngp), weig(ngp)


           m = (ngp + 1) / 2
!* Roots are symmetric in the interval - so only need to find half of them  */

           do i = 1, m                                ! Loop over the desired roots */

                     z = cos( M_PI * (i-0.25d0) / (ngp+0.5d0) )
!*   Starting with the above approximation to the ith root,
!*          we enter the main loop of refinement by NEWTON'S method   */

100             p1 = 1.0d0
                p2 = 0.0d0
!*  Loop up the recurrence relation to get the Legendre
!*  polynomial evaluated at z                 */

                do j = 1, ngp
                   p3 = p2
                   p2 = p1
                   p1 = ((2.0d0*j-1.0d0) * z * p2 - (j-1.0d0)*p3) / j
                enddo

!* p1 is now the desired Legendre polynomial. We next compute pp,
!* its derivative, by a standard relation involving also p2, the
!* polynomial of one lower order.      */
                pp = ngp*(z*p1-p2)/(z*z-1.0d0)
                z1 = z
                z = z1 - p1/pp             ! Newton's Method  */

                if (dabs(z-z1) .gt. EPS) GOTO  100

              xabsc(i) =  - z                            ! Roots will be bewteen -1.0 & 1.0 */
              xabsc(ngp+1-i) =  + z                        ! and symmetric about the origin   */
              weig(i) = 2.0d0/((1.0d0-z*z)*pp*pp)     ! Compute the weight and its       */
              weig(ngp+1-i) = weig(i)                 ! symmetric counterpart            */

      end do     ! i loop

   End subroutine MAKE_GL_NW

