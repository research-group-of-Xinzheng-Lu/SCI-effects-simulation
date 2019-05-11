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

!> @brief Performs the Newton Rapson algorithm.
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] xnod x-coord of the point
!> @param[in] ynod y-coord of the point
!> @param[in] znod z-coord of the point
!> @param[in] alfa11 costant for the bilinear map
!> @param[in] alfa12 costant for the bilinear map
!> @param[in] alfa13 costant for the bilinear map
!> @param[in] alfa21 costant for the bilinear map
!> @param[in] alfa22 costant for the bilinear map
!> @param[in] alfa23 costant for the bilinear map
!> @param[in] alfa31 costant for the bilinear map
!> @param[in] alfa32 costant for the bilinear map
!> @param[in] alfa33 costant for the bilinear map
!> @param[in] beta11 costant for the bilinear map
!> @param[in] beta12 costant for the bilinear map
!> @param[in] beta13 costant for the bilinear map
!> @param[in] beta21 costant for the bilinear map
!> @param[in] beta22 costant for the bilinear map
!> @param[in] beta23 costant for the bilinear map
!> @param[in] beta31 costant for the bilinear map
!> @param[in] beta32 costant for the bilinear map
!> @param[in] beta33 costant for the bilinear map
!> @param[in] gamma1  costant for the bilinear map
!> @param[in] gamma2  costant for the bilinear map
!> @param[in] gamma3  costant for the bilinear map
!> @param[in] delta1  costant for the bilinear map
!> @param[in] delta2  costant for the bilinear map
!> @param[in] delta3  costant for the bilinear map
!> @param[in] tt   1 if the node belongs to the element, 0 otherwise
!> @param[in] nofi  number of iterations (dummy)
!> @param[in] id_mpi  mpi processor identity
!> @param[in] ipiu control parameter (dummy)
!> @param[in] imeno control parameter (dummy)
!> @param[in] toll1 fixed tolerance
!> @param[in] toll2  fixed tolerance
!> @param[in] flag  control parameter
!> @param[out] csi_s x-coordinate in the reference element of xnod,ynod,znod
!> @param[out] eta_s y-coordinate in the reference element of xnod,ynod,znod
!> @param[out] zeta_s z-coordinate in the reference element of xnod,ynod,znod

     subroutine NEWTON_RAPSON(xnod, ynod, znod, &
                           alfa11,alfa12,alfa13, &
                           alfa21,alfa22,alfa23, &
                           alfa31,alfa32,alfa33, &
                           beta11,beta12,beta13, &
                           beta21,beta22,beta23, &
                           beta31,beta32,beta33, & 
                           gamma1,gamma2,gamma3, &
                           delta1,delta2,delta3, & 
                           tt, csi_s, eta_s, zeta_s, nofi, id_mpi,&
                           ipiu, imeno, toll1, toll2 , flag)

      implicit none
                      
      integer*4 :: nofiter
      integer*4 :: per_csi, per_eta, per_zeta
      integer*4 ::  iic, il,ih
      integer*4, intent(out) :: tt 
      integer*4, intent(in) :: nofi,id_mpi,ipiu,imeno, flag
      
      real*8 :: ll,lm ,xvt1,yvt1,xvt2,yvt2,xvc1,yvc1,xvc2,yvc2
      real*8 :: a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3
      real*8 :: e1,e2,e3,f1,f2,f3,g1,g2,g3,h1,h2,h3
      real*8 :: delta_csi, delta_eta, delta_zeta
      real*8 :: A,B,C,D,E,F,G,H,P, det
      real*8 :: a11, a12, a13, a21, a22, a23, a31, a32, a33
      real*8 :: FF1, FF2, FF3
      real*8 :: alfa11,alfa12,alfa13,alfa21,alfa22,alfa23,alfa31,alfa32,alfa33
      real*8 :: beta11,beta12,beta13,beta21,beta22,beta23,beta31,beta32,beta33
      real*8 :: gamma1,gamma2,gamma3,delta1,delta2,delta3
      real*8, intent(in) :: xnod, ynod, znod, toll1, toll2
      real*8, intent(out) ::  csi_s, eta_s, zeta_s


      a1 = 8.d0*gamma1
      a2 = 8.d0*gamma2
      a3 = 8.d0*gamma3
      b1 = 8.d0*beta13
      b2 = 8.d0*beta23
      b3 = 8.d0*beta33
      c1 = 8.d0*beta11
      c2 = 8.d0*beta21
      c3 = 8.d0*beta31
      d1 = 8.d0*beta12
      d2 = 8.d0*beta22
      d3 = 8.d0*beta32
      e1 = 8.d0*alfa11
      e2 = 8.d0*alfa21
      e3 = 8.d0*alfa31
      f1 = 8.d0*alfa12
      f2 = 8.d0*alfa22
      f3 = 8.d0*alfa32
      g1 = 8.d0*alfa13
      g2 = 8.d0*alfa23
      g3 = 8.d0*alfa33
      h1 = 8.d0*delta1
      h2 = 8.d0*delta2
      h3 = 8.d0*delta3

      csi_s = 0.d0
      eta_s = 0.d0
      zeta_s = 0.d0
      delta_csi = 0.0d0
      delta_eta = 0.0d0 
      delta_zeta = 0.0d0     
      
      tt = 0
      nofiter = 2000

      per_csi = 0
      per_eta = 0
      per_zeta = 0


      do il = 1, nofiter
                
        A = a1*eta_s*zeta_s + b1*eta_s + d1*zeta_s + e1
        B = a1*csi_s*zeta_s + b1*csi_s + c1*zeta_s + f1
        C = a1*csi_s*eta_s + c1*eta_s + d1*csi_s + g1

        D = a2*eta_s*zeta_s + b2*eta_s + d2*zeta_s + e2
        E = a2*csi_s*zeta_s + b2*csi_s + c2*zeta_s + f2
        F = a2*csi_s*eta_s + c2*eta_s + d2*csi_s + g2

        G = a3*eta_s*zeta_s + b3*eta_s + d3*zeta_s + e3
        H = a3*csi_s*zeta_s + b3*csi_s + c3*zeta_s + f3 
        P = a3*csi_s*eta_s + c3*eta_s + d3*csi_s + g3

        det = A*(E*P - F*H) - B*(D*P - F*G) + C*(D*H -E*G)
        
        if(det .eq. 0.d0) then
          return
        endif

        a11 = E*P - F*H
        a12 = -D*P + F*G
        a13 = D*H - E*G
        a21 = -B*P + C*H 
        a22 = A*P - C*G
        a23 = -A*H + B*G
        a31 = B*F - C*E
        a32 = -A*F + C*D
        a33 = A*E - B*D
        
        FF1 = a1*csi_s*eta_s*zeta_s + b1*csi_s*eta_s + c1*eta_s*zeta_s + d1*csi_s*zeta_s &
                                   + e1*csi_s + f1*eta_s + g1*zeta_s + h1 -8.d0*xnod 
        
        FF2 = a2*csi_s*eta_s*zeta_s + b2*csi_s*eta_s + c2*eta_s*zeta_s + d2*csi_s*zeta_s &
                                   + e2*csi_s + f2*eta_s + g2*zeta_s + h2 -8.d0*ynod 
        
        FF3 = a3*csi_s*eta_s*zeta_s + b3*csi_s*eta_s + c3*eta_s*zeta_s + d3*csi_s*zeta_s &
                                   + e3*csi_s + f3*eta_s + g3*zeta_s + h3 -8.d0*znod 


        delta_csi = (-1.d0/det)*(a11*FF1 + a21*FF2 + a31*FF3)

        delta_eta = (-1.d0/det)*(a12*FF1 + a22*FF2 + a32*FF3)

        delta_zeta = (-1.d0/det)*(a13*FF1 + a23*FF2 + a33*FF3)


         
        if (abs(delta_csi).le. toll1 .and. abs(delta_eta).le. toll1 .and. abs(delta_zeta).le. toll1 ) then
            tt = 1
           
           if(flag.eq.1) then 
              if (dabs(csi_s) .gt. toll2 .or. dabs(eta_s) .gt. toll2 .or. dabs(zeta_s).gt. toll2 ) then
              tt = 0
              return
            endif
          endif
         
          return
        endif
        
       
        csi_s = csi_s + delta_csi
        eta_s = eta_s + delta_eta
        zeta_s = zeta_s + delta_zeta
        
        if (csi_s .gt. 1.d0 .AND. csi_s .lt. 2.d0 .AND. per_csi .eq. 0) then
            csi_s = 1.d0 
            per_csi = 1
        endif
             
        if (csi_s .lt. -1.d0 .AND. csi_s .gt. -2.d0 .AND. per_csi .eq. 0) then
            csi_s = -1.d0
            per_csi = 1
        endif
        
        
        if (eta_s .gt. 1.d0 .AND. eta_s .lt. 2.d0 .AND. per_eta .eq. 0) then
            eta_s = 1.d0 
            per_eta = 1
        endif
            
        if (eta_s .lt. -1.d0 .AND. eta_s .gt. -2.d0 .AND. per_eta .eq. 0) then
            eta_s = -1.d0
            per_eta = 1
        endif    
            
        if (zeta_s .gt. 1.d0 .AND. zeta_s .lt. 2.d0 .AND. per_zeta .eq. 0) then
            zeta_s = 1.d0 
            per_zeta = 1
        endif
            
        if (zeta_s .lt. -1.d0 .AND. zeta_s .gt. -2.d0 .AND. per_zeta .eq. 0) then
            zeta_s = -1.d0
            per_zeta = 1
        endif    


      enddo   
      
      if(il .ge. nofiter) then      
        tt = 0
        return
      endif

     end subroutine NEWTON_RAPSON
