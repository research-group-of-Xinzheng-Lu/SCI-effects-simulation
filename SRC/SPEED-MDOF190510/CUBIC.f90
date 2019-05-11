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

!> @brief Find roots of the cubic equation  x**3 + a*x**2 + b*x + c = 0
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] a coefficient
!> @param[in] b coefficient
!> @param[in] c coefficient
!> @param[out] x1 first root 
!> @param[out] x2 second root
!> @param[out] x3 third root

      subroutine CUBIC(a,b,c,x1,x2,x3)       

      
      implicit none

      real*8 :: a,b,c,D,x1,x2,x3,pi,prtr,prti
      real*8 :: oot,opf,three,srth,p,q
      real*8 :: srd,tmp,u,v, cosphi, phi, cf



     !------------------------------------------------------
     ! FDLIB
     !
     ! Copyright by C. Pozrikidis, 1999
     ! All rights reserved.
     !
     ! This program is to be used only under the
     ! stipulations of the licensing agreement.
     !-------------------------------------------------------

     !-------------------------------------------------------
     ! Roots of the cubic equation:
     !
     !  x**3 + a*x**2 + b*x + c = 0
     !
     ! where a,b,c are real coefficients.
     !
     ! The roots are computed
     ! using Cardano's analytical formulae
     !
     ! See:
     !
     ! John W Harris and Horst Stocker,
     ! ``Handbook of Mathematics and Computational Science''
     ! Springer (1998).
     !
     ! john.harris@yale.edu
     !
     !-------------------------------------------------------

     !Implicit Double Precision (a-h,o-z)


     !----------
     ! constants
     !----------


     
      
      pi = 3.14159265358d0

      oot   = 1.0d0/3.0d0
      opf   = 1.5d0
      three = 3.0d0
      srth  = dsqrt(three)

      !--------
      ! prepare
      !--------

      p = (3.0D0*b-a**2)/3.0D0
      q = c+2.0D0*a**3/27.0D0-a*b/3.0D0
      D = (p/3.0D0)**3+(q/2.0D0)**2        ! discriminant

      !--------------------------------------
      ! one real, two complex conjugate roots
      !--------------------------------------

      If(D.ge.0) then 

       srd  = Dsqrt(D)
       tmp  =  -0.5D0*q+srd
       u    =  Dabs(tmp)**oot
       If(tmp.lt.0) u = -u
       tmp  =  -0.5D0*q-srd
       v    =  Dabs(tmp)**oot
       If(tmp.lt.0) v = -v

       x1   = -a/3.0D0+u+v

       prtr = -a/3.0D0-0.5D0*(u+v)
       prti =   srth*0.5D0*(u-v)

       x2 = prtr
       x3 = prtr


      !-----------------
      ! three real roots
      !-----------------

      Else

       cosphi = -0.5D0*q/(Dabs(p)/3.0D0)**opf
       phi    = Dacos(cosphi)

       cf = 2.0D0*sqrt(Dabs(p)/3.0D0)

       x1 = -a/3.0D0 + cf*Dcos(phi/3.0D0)
       x2 = -a/3.0D0 - cf*Dcos((phi-pi)/3.0D0)
       x3 = -a/3.0D0 - cf*Dcos((phi+pi)/3.0D0)


      !-----------

      EndIf




      return

     
      end subroutine CUBIC

