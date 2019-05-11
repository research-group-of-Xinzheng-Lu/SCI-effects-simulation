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

!> @brief Computes DG matrices for jumps.
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] x_pl  x-coord of Gauss Legendre node for the integration in Omega+
!> @param[in] y_pl  y-coord of Gauss Legendre node for the integration in Omega+
!> @param[in] z_pl  z-coord of Gauss Legendre node for the integration in Omega+
!> @param[in] wx_pl Gauss Legendre weight
!> @param[in] wy_pl Gauss Legendre weight
!> @param[in] wz_pl Gauss Legendre weight 
!> @param[in] x_mn  x-coord. of Gauss Legendre nodes for the integration in Omega-
!> @param[in] y_mn  y-coord. of Gauss Legendre nodes for the integration in Omega-
!> @param[in] z_mn  z-coord. of Gauss Legendre nodes for the integration in Omega-
!> @param[in] nq  number of quadrature pints
!> @param[in] nn  polynomial degree in Omega+
!> @param[in] mm  polynomial degree in Omega-
!> @param[in] o_minus  matrix containing info about Omega- o_minus(0,i) = material, 
!!            o_minus(1,i) = el index, o_minus(2,i) = face, o_minus(3,1) = neighbouring element in 
!!            a local numeration
!> @param[in] alfa11 costant values for the bilinear map of Omega+
!> @param[in] alfa12 costant values for the bilinear map of Omega+
!> @param[in] alfa13 costant values for the bilinear map of Omega+
!> @param[in] alfa21 costant values for the bilinear map of Omega+
!> @param[in] alfa22 costant values for the bilinear map of Omega+
!> @param[in] alfa23 costant values for the bilinear map of Omega+
!> @param[in] alfa31 costant values for the bilinear map of Omega+
!> @param[in] alfa32 costant values for the bilinear map of Omega+
!> @param[in] alfa33 costant values for the bilinear map of Omega+
!> @param[in] beta11 costant values for the bilinear map of Omega+
!> @param[in] beta12 costant values for the bilinear map of Omega+
!> @param[in] beta13 costant values for the bilinear map of Omega+
!> @param[in] beta21 costant values for the bilinear map of Omega+
!> @param[in] beta22 costant values for the bilinear map of Omega+
!> @param[in] beta23 costant values for the bilinear map of Omega+
!> @param[in] beta31 costant values for the bilinear map of Omega+
!> @param[in] beta32 costant values for the bilinear map of Omega+
!> @param[in] beta33 costant values for the bilinear map of Omega+
!> @param[in] gamma1 costant values for the bilinear map of Omega+
!> @param[in] gamma2 costant values for the bilinear map of Omega+
!> @param[in] gamma3 costant values for the bilinear map of Omega+
!> @param[in] delta1 costant values for the bilinear map of Omega+
!> @param[in] delta2 costant values for the bilinear map of Omega+
!> @param[in] delta3 costant values for the bilinear map of Omega+   
!> @param[in] ine index for neighbouring element
!> @param[in] alfa11_mn costant values for the bilinear map of Omega-
!> @param[in] alfa12_mn costant values for the bilinear map of Omega-
!> @param[in] alfa13_mn costant values for the bilinear map of Omega-
!> @param[in] alfa21_mn costant values for the bilinear map of Omega-
!> @param[in] alfa22_mn costant values for the bilinear map of Omega-
!> @param[in] alfa23_mn costant values for the bilinear map of Omega-
!> @param[in] alfa31_mn costant values for the bilinear map of Omega-
!> @param[in] alfa32_mn costant values for the bilinear map of Omega-
!> @param[in] alfa33_mn costant values for the bilinear map of Omega-
!> @param[in] beta11_mn costant values for the bilinear map of Omega-
!> @param[in] beta12_mn costant values for the bilinear map of Omega-
!> @param[in] beta13_mn costant values for the bilinear map of Omega-
!> @param[in] beta21_mn costant values for the bilinear map of Omega-
!> @param[in] beta22_mn costant values for the bilinear map of Omega-
!> @param[in] beta23_mn costant values for the bilinear map of Omega-
!> @param[in] beta31_mn costant values for the bilinear map of Omega-
!> @param[in] beta32_mn costant values for the bilinear map of Omega-
!> @param[in] beta33_mn costant values for the bilinear map of Omega-
!> @param[in] gamma1_mn costant values for the bilinear map of Omega-
!> @param[in] gamma2_mn costant values for the bilinear map of Omega-
!> @param[in] gamma3_mn costant values for the bilinear map of Omega-
!> @param[in] delta1_mn costant values for the bilinear map of Omega-
!> @param[in] delta2_mn costant values for the bilinear map of Omega-
!> @param[in] delta3_mn costant values for the bilinear map of Omega-
!> @param[in] cp_a constant appearing in the interface integrals
!> @param[in] cp_b constant appearing in the interface integrals
!> @param[in] cp_c constant appearing in the interface integrals
!> @param[in] cp_d constant appearing in the interface integrals
!> @param[in] cp_e constant appearing in the interface integrals
!> @param[in] cp_f constant appearing in the interface integrals
!> @param[in] cp_g constant appearing in the interface integrals
!> @param[in] cp_n constant appearing in the interface integrals
!> @param[in] cp_p constant appearing in the interface integrals
!> @param[in] pen  penalty constant
!> @param[in] dg_cnst  -1 SIPG, 0 IIPG, 1 NIPG
!> @param[in] id_mpi  mpi process identity
!> @param[in] det_scal  determinant of the transformation used for computing interface integrals
!> @param[in] testmode 1 if testmode is active
!> @param[out] JP  jump matrix for the part +,+
!> @param[out] JM  jump matrix for the part +,-
!> @param[out] JP_uv  jump matrix for the part +,+ for testmode case
!> @param[out] JM_uv  jump matrix for the part +,- for testmode case


      subroutine MAKE_LOC_MATRIX_DG(x_pl, y_pl, z_pl, &
                       wx_pl,wy_pl, wz_pl,&
                       x_mn, y_mn, z_mn, &
                       nq, nn, mm, &
                       o_minus, &
                       alfa11,alfa12,alfa13,&
                       alfa21,alfa22,alfa23,&
                       alfa31,alfa32,alfa33,&
                       beta11,beta12,beta13,&
                       beta21,beta22,beta23,&
                       beta31,beta32,beta33,&
                       gamma1,gamma2,gamma3,&
                       delta1,delta2,delta3,&
                       ine,&
                       alfa11_mn,alfa12_mn,alfa13_mn,&
                       alfa21_mn,alfa22_mn,alfa23_mn,&
                       alfa31_mn,alfa32_mn,alfa33_mn,&
                       beta11_mn,beta12_mn,beta13_mn,&
                       beta21_mn,beta22_mn,beta23_mn,&
                       beta31_mn,beta32_mn,beta33_mn,&
                       gamma1_mn,gamma2_mn,gamma3_mn,&
                       delta1_mn,delta2_mn,delta3_mn,&
                       cp_a,cp_b,cp_c,cp_d,cp_e,cp_f,cp_g,cp_n,cp_p, &
                       pen, dg_cnst,JP,JM,id_mpi,det_scal, testmode,&
                       JP_uv, JM_uv)


     use max_var
     use str_mesh
     
     implicit none
     
     
     integer*4 :: i, id_mpi, testmode
     integer*4 :: ip1 , i1, j1, h1, m1, n1, p1, L1, K1
     integer*4 :: ip2 , i2, j2, h2, m2, n2, p2, L2, K2
     integer*4 :: ip3 , i3, j3, h3, m3, n3, p3, L3, K3
     !integer*4 :: nthreads,  OMP_GET_NUM_THREADS

     integer*4, intent(in) :: nn, nq, mm, ine  

     integer*4, dimension(1:nq,0:3), intent(in) :: o_minus

     real*8 :: phi, dphi
     real*8 :: dxdx1,dxdy1,dxdz1,dydx1,dydy1,dydz1,dzdx1,dzdy1,dzdz1,det_j1, det_s1
     real*8 :: dxdx2,dxdy2,dxdz2,dydx2,dydy2,dydz2,dzdx2,dzdy2,dzdz2,det_j2, det_s2
     real*8 :: dxdx3,dxdy3,dxdz3,dydx3,dydy3,dydz3,dzdx3,dzdy3,dzdz3,det_j3, det_s3
     real*8 :: dxdu1, dydu1, dzdu1, dxdv1, dydv1, dzdv1
     real*8 :: dxdu2, dydu2, dzdu2, dxdv2, dydv2, dzdv2
     real*8 :: dxdu3, dydu3, dzdu3, dxdv3, dydv3, dzdv3
     real*8 :: dxdx_m,dxdy_m,dxdz_m,dydx_m,dydy_m,dydz_m,dzdx_m,dzdy_m,dzdz_m,det_j_m
     real*8 :: dpdc1, dpde1, dpdz1, dpdx1, dpdy1, dpdzt1, plx1, pkx1
     real*8 :: dpdc2, dpde2, dpdz2, dpdx2, dpdy2, dpdzt2, plx2, pkx2
     real*8 :: dpdc3, dpde3, dpdz3, dpdx3, dpdy3, dpdzt3, plx3, pkx3
     real*8 :: tempo1, tempo2

     real*8, intent(in) :: alfa11, alfa12, alfa13, alfa21, alfa22, alfa23, alfa31, alfa32, alfa33
     real*8, intent(in) :: beta11, beta12, beta13, beta21, beta22, beta23, beta31, beta32, beta33
     real*8, intent(in) :: gamma1, gamma2, gamma3, delta1, delta2, delta3
     real*8, intent(in) :: alfa11_mn, alfa12_mn, alfa13_mn 
     real*8, intent(in) :: alfa21_mn, alfa22_mn, alfa23_mn
     real*8, intent(in) :: alfa31_mn, alfa32_mn, alfa33_mn
     real*8, intent(in) :: beta11_mn, beta12_mn, beta13_mn
     real*8, intent(in) :: beta21_mn, beta22_mn, beta23_mn 
     real*8, intent(in) :: beta31_mn, beta32_mn, beta33_mn
     real*8, intent(in) :: gamma1_mn, gamma2_mn, gamma3_mn
     real*8, intent(in) :: delta1_mn, delta2_mn, delta3_mn
     real*8, intent(in) :: cp_a, cp_b, cp_c, cp_d, cp_e, cp_f, cp_g, cp_n, cp_p
     real*8, intent(in) :: pen, dg_cnst, det_scal
     real*8, dimension(nq), intent(in) :: x_pl, y_pl, z_pl, wx_pl, wy_pl, wz_pl, x_mn, y_mn, z_mn 
     real*8, dimension(nn) :: ctp
     real*8, dimension(mm) :: ctm

     real*8, dimension(:,:), allocatable :: DD,EE,FF,GG,HH,II,QQ
     real*8, dimension(:,:), allocatable :: AA,BB,CC,PP
     real*8, dimension(3*nn**3,3*nn**3), intent(out)  :: JP, JP_uv    
     real*8, dimension(3*nn**3,3*mm**3), intent(out)  :: JM, JM_uv
    
    call MAKE_LGL_NODES(nn,ctp)
    call MAKE_LGL_NODES(mm,ctm)

!********************************************************************************************!   
!   MATRICES (+,+)                                                                           !
!********************************************************************************************!
     
     
     
  allocate(AA(nn**3,nn**3),BB(nn**3,nn**3),CC(nn**3,nn**3),PP(nn**3,nn**3))

    AA = 0.d0
    BB = 0.d0
    CC = 0.d0
    PP = 0.d0
 

    allocate(DD(nn**3,mm**3),EE(nn**3,mm**3),FF(nn**3,mm**3),QQ(nn**3,mm**3))

    DD = 0.d0
    EE = 0.d0
    FF = 0.d0
    QQ = 0.d0



    allocate(GG(nn**3,mm**3),HH(nn**3,mm**3),II(nn**3,mm**3))
    GG = 0.d0
    HH = 0.d0
    II = 0.d0 



!!$OMP PARALLEL      
!!$OMP SECTIONS

!!$OMP SECTION


     do p1 = 1, nn
        do n1 = 1, nn
          do m1 = 1, nn
             L1 = m1 + (n1-1)*nn + (p1-1)*nn*nn
     
     
             do h1 = 1, nn
                do j1 = 1, nn
                   do i1 = 1, nn
                      K1 = i1 + (j1-1)*nn + (h1-1)*nn*nn
                        
                      do  ip1 = 1, nq   
                      
                          if( o_minus(ip1,1) .eq. ine) then
       
               dxdx1 = alfa11 + beta12*z_pl(ip1) + beta13*y_pl(ip1) + gamma1*y_pl(ip1)*z_pl(ip1)
               dydx1 = alfa21 + beta22*z_pl(ip1) + beta23*y_pl(ip1) + gamma2*y_pl(ip1)*z_pl(ip1)
               dzdx1 = alfa31 + beta32*z_pl(ip1) + beta33*y_pl(ip1) + gamma3*y_pl(ip1)*z_pl(ip1)
               
               dxdy1 = alfa12 + beta11*z_pl(ip1) + beta13*x_pl(ip1) + gamma1*z_pl(ip1)*x_pl(ip1)
               dydy1 = alfa22 + beta21*z_pl(ip1) + beta23*x_pl(ip1) + gamma2*z_pl(ip1)*x_pl(ip1)
               dzdy1 = alfa32 + beta31*z_pl(ip1) + beta33*x_pl(ip1) + gamma3*z_pl(ip1)*x_pl(ip1)
               
               dxdz1 = alfa13 + beta11*y_pl(ip1) + beta12*x_pl(ip1) + gamma1*x_pl(ip1)*y_pl(ip1)
               dydz1 = alfa23 + beta21*y_pl(ip1) + beta22*x_pl(ip1) + gamma2*x_pl(ip1)*y_pl(ip1)
               dzdz1 = alfa33 + beta31*y_pl(ip1) + beta32*x_pl(ip1) + gamma3*x_pl(ip1)*y_pl(ip1)
               
               det_j1 = dxdz1 * (dydx1*dzdy1 - dzdx1*dydy1) &
                     - dydz1 * (dxdx1*dzdy1 - dzdx1*dxdy1) &
                     + dzdz1 * (dxdx1*dydy1 - dydx1*dxdy1)


              dpdc1 = dphi(nn-1, ctp(i1), x_pl(ip1)) * phi(nn-1, ctp(j1), y_pl(ip1)) * phi(nn-1, ctp(h1), z_pl(ip1))               
              dpde1 = phi(nn-1, ctp(i1), x_pl(ip1)) * dphi(nn-1, ctp(j1), y_pl(ip1)) * phi(nn-1, ctp(h1), z_pl(ip1))
              dpdz1 = phi(nn-1, ctp(i1), x_pl(ip1)) * phi(nn-1, ctp(j1), y_pl(ip1)) * dphi(nn-1, ctp(h1), z_pl(ip1))

                plx1 = phi(nn-1, ctp(m1), x_pl(ip1)) * phi(nn-1, ctp(n1), y_pl(ip1)) * phi(nn-1, ctp(p1), z_pl(ip1))
                pkx1 = phi(nn-1, ctp(i1), x_pl(ip1)) * phi(nn-1, ctp(j1), y_pl(ip1)) * phi(nn-1, ctp(h1), z_pl(ip1))


              dpdx1 = (dydx1*(dzdy1*dpdz1 - dzdz1*dpde1) - dydy1*(dzdx1*dpdz1 - dzdz1*dpdc1) + &
                       dydz1*(dzdx1*dpde1 - dzdy1*dpdc1))/det_j1
 
              dpdy1 = (dzdx1*(dxdy1*dpdz1 - dxdz1*dpde1) - dzdy1*(dxdx1*dpdz1 - dxdz1*dpdc1) + &
                      dzdz1*(dxdx1*dpde1 - dxdy1*dpdc1))/det_j1


             dpdzt1 = (dxdx1*(dydy1*dpdz1 - dydz1*dpde1) - dxdy1*(dydx1*dpdz1 - dydz1*dpdc1) +  &
                      dxdz1*(dydx1*dpde1 - dydy1*dpdc1))/det_j1
                      

            det_s1 = det_scal/4.d0
              

              AA(L1,K1) = AA(L1,K1)  +  det_s1 * wx_pl(ip1)*wy_pl(ip1)*wz_pl(ip1) * dpdx1 * plx1  
              BB(L1,K1) = BB(L1,K1)  +  det_s1 * wx_pl(ip1)*wy_pl(ip1)*wz_pl(ip1) * dpdy1 * plx1 
              CC(L1,K1) = CC(L1,K1)  +  det_s1 * wx_pl(ip1)*wy_pl(ip1)*wz_pl(ip1) * dpdzt1 * plx1  
              PP(L1,K1) = PP(L1,K1)  +  det_s1 * wx_pl(ip1)*wy_pl(ip1)*wz_pl(ip1) * pkx1 * plx1  
              
              

                           endif
               
  
                      enddo   ! end loop on ip1

                  enddo
               enddo
            enddo  !end loop on K
            
           
        enddo
     enddo
  enddo    !end loop in L
         
!********************************************************************************************!   
!   MATRICES (-,+)                                                                           !
!********************************************************************************************!
         
!!$OMP SECTION

       do p2 = 1, nn
        do n2 = 1, nn
          do m2 = 1, nn
             L2 = m2 + (n2-1)*nn + (p2-1)*nn*nn
     
     
             do h2 = 1, mm
                do j2 = 1, mm
                   do i2 = 1, mm
                      K2 = i2 + (j2-1)*mm + (h2-1)*mm*mm
                    
      
     
     
                      do  ip2 = 1, nq  

                          if( o_minus(ip2,1) .eq. ine) then
       
               dxdx2 = alfa11 + beta12*z_pl(ip2) + beta13*y_pl(ip2) + gamma1*y_pl(ip2)*z_pl(ip2)
               dydx2 = alfa21 + beta22*z_pl(ip2) + beta23*y_pl(ip2) + gamma2*y_pl(ip2)*z_pl(ip2)
               dzdx2 = alfa31 + beta32*z_pl(ip2) + beta33*y_pl(ip2) + gamma3*y_pl(ip2)*z_pl(ip2)
              
               dxdy2 = alfa12 + beta11*z_pl(ip2) + beta13*x_pl(ip2) + gamma1*z_pl(ip2)*x_pl(ip2)
               dydy2 = alfa22 + beta21*z_pl(ip2) + beta23*x_pl(ip2) + gamma2*z_pl(ip2)*x_pl(ip2)
               dzdy2 = alfa32 + beta31*z_pl(ip2) + beta33*x_pl(ip2) + gamma3*z_pl(ip2)*x_pl(ip2)
               
               dxdz2 = alfa13 + beta11*y_pl(ip2) + beta12*x_pl(ip2) + gamma1*x_pl(ip2)*y_pl(ip2)
               dydz2 = alfa23 + beta21*y_pl(ip2) + beta22*x_pl(ip2) + gamma2*x_pl(ip2)*y_pl(ip2)
               dzdz2 = alfa33 + beta31*y_pl(ip2) + beta32*x_pl(ip2) + gamma3*x_pl(ip2)*y_pl(ip2)


               dxdx_m = alfa11_mn + beta12_mn*z_mn(ip2) + beta13_mn*y_mn(ip2) + gamma1_mn*y_mn(ip2)*z_mn(ip2)
               dydx_m = alfa21_mn + beta22_mn*z_mn(ip2) + beta23_mn*y_mn(ip2) + gamma2_mn*y_mn(ip2)*z_mn(ip2)
               dzdx_m = alfa31_mn + beta32_mn*z_mn(ip2) + beta33_mn*y_mn(ip2) + gamma3_mn*y_mn(ip2)*z_mn(ip2)
               
               dxdy_m = alfa12_mn + beta11_mn*z_mn(ip2) + beta13_mn*x_mn(ip2) + gamma1_mn*z_mn(ip2)*x_mn(ip2)
               dydy_m = alfa22_mn + beta21_mn*z_mn(ip2) + beta23_mn*x_mn(ip2) + gamma2_mn*z_mn(ip2)*x_mn(ip2)
               dzdy_m = alfa32_mn + beta31_mn*z_mn(ip2) + beta33_mn*x_mn(ip2) + gamma3_mn*z_mn(ip2)*x_mn(ip2)
               
               dxdz_m = alfa13_mn + beta11_mn*y_mn(ip2) + beta12_mn*x_mn(ip2) + gamma1_mn*x_mn(ip2)*y_mn(ip2)
               dydz_m = alfa23_mn + beta21_mn*y_mn(ip2) + beta22_mn*x_mn(ip2) + gamma2_mn*x_mn(ip2)*y_mn(ip2)
               dzdz_m = alfa33_mn + beta31_mn*y_mn(ip2) + beta32_mn*x_mn(ip2) + gamma3_mn*x_mn(ip2)*y_mn(ip2)
               
               det_j_m = dxdz_m * (dydx_m*dzdy_m - dzdx_m*dydy_m) &
                       - dydz_m * (dxdx_m*dzdy_m - dzdx_m*dxdy_m) &
                       + dzdz_m * (dxdx_m*dydy_m - dydx_m*dxdy_m)



              dpdc2 = dphi(mm-1, ctm(i2), x_mn(ip2)) * phi(mm-1, ctm(j2), y_mn(ip2)) * phi(mm-1, ctm(h2), z_mn(ip2))
              dpde2 = phi(mm-1, ctm(i2), x_mn(ip2)) * dphi(mm-1, ctm(j2), y_mn(ip2)) * phi(mm-1, ctm(h2), z_mn(ip2))
              dpdz2 = phi(mm-1, ctm(i2), x_mn(ip2)) * phi(mm-1, ctm(j2), y_mn(ip2)) * dphi(mm-1, ctm(h2), z_mn(ip2))
               

                plx2 = phi(nn-1, ctp(m2), x_pl(ip2)) * phi(nn-1, ctp(n2), y_pl(ip2)) * phi(nn-1, ctp(p2), z_pl(ip2))
                pkx2 = phi(mm-1, ctm(i2), x_mn(ip2)) * phi(mm-1, ctm(j2), y_mn(ip2)) * phi(mm-1, ctm(h2), z_mn(ip2))
                

              dpdx2 = (dydx_m*(dzdy_m*dpdz2 - dzdz_m*dpde2) - dydy_m*(dzdx_m*dpdz2 - dzdz_m*dpdc2) + &
                       dydz_m*(dzdx_m*dpde2 - dzdy_m*dpdc2))/det_j_m
 
              dpdy2 = (dzdx_m*(dxdy_m*dpdz2 - dxdz_m*dpde2) - dzdy_m*(dxdx_m*dpdz2 - dxdz_m*dpdc2) + &
                      dzdz_m*(dxdx_m*dpde2 - dxdy_m*dpdc2))/det_j_m

             dpdzt2 = (dxdx_m*(dydy_m*dpdz2 - dydz_m*dpde2) - dxdy_m*(dydx_m*dpdz2 - dydz_m*dpdc2) + &
                      dxdz_m*(dydx_m*dpde2 - dydy_m*dpdc2))/det_j_m
                      

              det_s2 = det_scal/4.d0            

              DD(L2,K2) = DD(L2,K2)  +  det_s2 * wx_pl(ip2)*wy_pl(ip2)*wz_pl(ip2) * dpdx2 * plx2  
              EE(L2,K2) = EE(L2,K2)  +  det_s2 * wx_pl(ip2)*wy_pl(ip2)*wz_pl(ip2) * dpdy2 * plx2  
              FF(L2,K2) = FF(L2,K2)  +  det_s2 * wx_pl(ip2)*wy_pl(ip2)*wz_pl(ip2) * dpdzt2 * plx2  
              QQ(L2,K2) = QQ(L2,K2)  +  det_s2 * wx_pl(ip2)*wy_pl(ip2)*wz_pl(ip2) * pkx2 * plx2  
              
              

                           endif

  
                      enddo   ! end loop on ip2

                  enddo
               enddo
            enddo  !end loop on K
            
        enddo
     enddo
  enddo    !end loop on L

              
!********************************************************************************************!   
!   MATRICES (+,-)                                                                           !
!********************************************************************************************!

!!$OMP SECTION

       do p3 = 1, nn
        do n3 = 1, nn
          do m3 = 1, nn
             L3 = m3 + (n3-1)*nn + (p3-1)*nn*nn
     
     
             do h3 = 1, mm
                do j3 = 1, mm
                   do i3 = 1, mm
                      K3 = i3 + (j3-1)*mm + (h3-1)*mm*mm
                    
                      do  ip3 = 1, nq   
                      
                          if( o_minus(ip3,1) .eq. ine) then
       
               dxdx3 = alfa11 + beta12*z_pl(ip3) + beta13*y_pl(ip3) + gamma1*y_pl(ip3)*z_pl(ip3)
               dydx3 = alfa21 + beta22*z_pl(ip3) + beta23*y_pl(ip3) + gamma2*y_pl(ip3)*z_pl(ip3)
               dzdx3 = alfa31 + beta32*z_pl(ip3) + beta33*y_pl(ip3) + gamma3*y_pl(ip3)*z_pl(ip3)
               
               dxdy3 = alfa12 + beta11*z_pl(ip3) + beta13*x_pl(ip3) + gamma1*z_pl(ip3)*x_pl(ip3)
               dydy3 = alfa22 + beta21*z_pl(ip3) + beta23*x_pl(ip3) + gamma2*z_pl(ip3)*x_pl(ip3)
               dzdy3 = alfa32 + beta31*z_pl(ip3) + beta33*x_pl(ip3) + gamma3*z_pl(ip3)*x_pl(ip3)
               
               dxdz3 = alfa13 + beta11*y_pl(ip3) + beta12*x_pl(ip3) + gamma1*x_pl(ip3)*y_pl(ip3)
               dydz3 = alfa23 + beta21*y_pl(ip3) + beta22*x_pl(ip3) + gamma2*x_pl(ip3)*y_pl(ip3)
               dzdz3 = alfa33 + beta31*y_pl(ip3) + beta32*x_pl(ip3) + gamma3*x_pl(ip3)*y_pl(ip3)
               
               det_j3 = dxdz3 * (dydx3*dzdy3 - dzdx3*dydy3) &
                     - dydz3 * (dxdx3*dzdy3 - dzdx3*dxdy3) &
                     + dzdz3 * (dxdx3*dydy3 - dydx3*dxdy3)


              dpdc3 = dphi(nn-1, ctp(m3), x_pl(ip3)) * phi(nn-1, ctp(n3), y_pl(ip3)) * phi(nn-1, ctp(p3), z_pl(ip3))
              dpde3 = phi(nn-1, ctp(m3), x_pl(ip3)) * dphi(nn-1, ctp(n3), y_pl(ip3)) * phi(nn-1, ctp(p3), z_pl(ip3))
              dpdz3 = phi(nn-1, ctp(m3), x_pl(ip3)) * phi(nn-1, ctp(n3), y_pl(ip3)) * dphi(nn-1, ctp(p3), z_pl(ip3))
              
                             

                pkx3 = phi(mm-1, ctm(i3), x_mn(ip3)) * phi(mm-1, ctm(j3), y_mn(ip3)) * phi(mm-1, ctm(h3), z_mn(ip3))


              dpdx3 = (dydx3*(dzdy3*dpdz3 - dzdz3*dpde3) - dydy3*(dzdx3*dpdz3 - dzdz3*dpdc3) + &
                       dydz3*(dzdx3*dpde3 - dzdy3*dpdc3))/det_j3
 
              dpdy3 = (dzdx3*(dxdy3*dpdz3 - dxdz3*dpde3) - dzdy3*(dxdx3*dpdz3 - dxdz3*dpdc3) + &
                      dzdz3*(dxdx3*dpde3 - dxdy3*dpdc3))/det_j3


             dpdzt3 = (dxdx3*(dydy3*dpdz3 - dydz3*dpde3) - dxdy3*(dydx3*dpdz3 - dydz3*dpdc3) + &
                      dxdz3*(dydx3*dpde3 - dydy3*dpdc3))/det_j3
                      
                      

              det_s3 = det_scal/4.d0
 
              GG(L3,K3) = GG(L3,K3)  +  det_s3 * wx_pl(ip3)*wy_pl(ip3)*wz_pl(ip3) * dpdx3 * pkx3  
              HH(L3,K3) = HH(L3,K3)  +  det_s3 * wx_pl(ip3)*wy_pl(ip3)*wz_pl(ip3) * dpdy3 * pkx3  
              II(L3,K3) = II(L3,K3)  +  det_s3 * wx_pl(ip3)*wy_pl(ip3)*wz_pl(ip3) * dpdzt3 * pkx3  
              
              
              

                           endif

  
                      enddo   ! end loop on ip3

                  enddo
               enddo
            enddo  !end loop on K
            
        enddo
     enddo
  enddo    !end loop on L
                  
                  
!!$OMP END SECTIONS
!!$OMP BARRIER
!!$OMP END PARALLEL 


  JP = 0.d0          
  JP(1 : nn**3, 1 : nn**3) =  pen*PP - cp_a*AA - cp_c*BB - cp_d*CC &
               + dg_cnst*cp_a*transpose(AA) + dg_cnst*cp_c*transpose(BB) + dg_cnst*cp_d*transpose(CC)
  
  JP(1 : nn**3, nn**3+1 : 2*nn**3) = - cp_c*AA - cp_b*BB &
               + dg_cnst*cp_g*transpose(AA) + dg_cnst*cp_e*transpose(BB) 
             
  JP(1 : nn**3, 2*nn**3+1 : 3*nn**3) = - cp_d*AA - cp_b*CC &
               + dg_cnst*cp_p*transpose(AA) + dg_cnst*cp_e*transpose(CC)
               
! write(*,*) JP(1 : nn**3, 2*nn**3+1 : 3*nn**3) 

               
  JP(nn**3+1 : 2*nn**3, 1 : nn**3) = - cp_g*AA - cp_e*BB &
               + dg_cnst*cp_c*transpose(AA) + dg_cnst*cp_b*transpose(BB) 
  
  JP(nn**3+1 : 2*nn**3, nn**3+1 : 2*nn**3) = pen*PP - cp_e*AA - cp_f*BB - cp_d*CC &
               + dg_cnst*cp_e*transpose(AA) + dg_cnst*cp_f*transpose(BB) + dg_cnst*cp_d*transpose(CC)
  
  JP(nn**3+1 : 2*nn**3, 2*nn**3+1 : 3*nn**3) = - cp_d*BB - cp_g*CC &
               + dg_cnst*cp_p*transpose(BB) + dg_cnst*cp_c*transpose(CC)
               
  JP(2*nn**3+1 : 3*nn**3, 1 : nn**3) = - cp_p*AA - cp_e*CC &
               + dg_cnst*cp_d*transpose(AA) + dg_cnst*cp_b*transpose(CC)
  
  JP(2*nn**3+1 : 3*nn**3, nn**3+1 : 2*nn**3) = - cp_p*BB - cp_c*CC &
               + dg_cnst*cp_d*transpose(BB) + dg_cnst*cp_g*transpose(CC)
               
  JP(2*nn**3+1 : 3*nn**3, 2*nn**3+1 : 3*nn**3) = pen*PP - cp_e*AA - cp_c*BB - cp_n*CC &
               + dg_cnst*cp_e*transpose(AA) + dg_cnst*cp_c*transpose(BB) + dg_cnst*cp_n*transpose(CC)
  


  
  JM = 0.d0
  JM(1 : nn**3, 1 : mm**3) = - pen*QQ - cp_a*DD - cp_c*EE - cp_d*FF &
                             - dg_cnst*cp_a*GG -dg_cnst*cp_c*HH - dg_cnst*cp_d*II
  
  JM(1 : nn**3, mm**3+1 : 2*mm**3) = - cp_c*DD - cp_b*EE - dg_cnst*cp_g*GG - dg_cnst*cp_e*HH  
             
  JM(1 : nn**3, 2*mm**3+1 : 3*mm**3) = - cp_d*DD - cp_b*FF - dg_cnst*cp_p*GG - dg_cnst*cp_e*II
  
  JM(nn**3+1 : 2*nn**3, 1 : mm**3) = - cp_g*DD - cp_e*EE - dg_cnst*cp_c*GG - dg_cnst*cp_b*HH
  
  JM(nn**3+1 : 2*nn**3, mm**3+1 : 2*mm**3) = - pen*QQ - cp_e*DD - cp_f*EE - cp_d*FF &
                                             - dg_cnst*cp_e*GG - dg_cnst*cp_f*HH - dg_cnst*cp_d*II
  
  JM(nn**3+1 : 2*nn**3, 2*mm**3+1 : 3*mm**3) = - cp_d*EE - cp_g*FF - dg_cnst*cp_p*HH - dg_cnst*cp_c*II
  
  JM(2*nn**3+1 : 3*nn**3, 1 : mm**3) = - cp_p*DD - cp_e*FF - dg_cnst*cp_d*GG - dg_cnst*cp_b*II
  
  JM(2*nn**3+1 : 3*nn**3, mm**3+1 : 2*mm**3) = - cp_p*EE - cp_c*FF - dg_cnst*cp_d*HH - dg_cnst*cp_g*II 
  
  JM(2*nn**3+1 : 3*nn**3, 2*mm**3+1 : 3*mm**3) = - pen*QQ - cp_e*DD - cp_c*EE - cp_n*FF &
                                             - dg_cnst*cp_e*GG - dg_cnst*cp_c*HH - dg_cnst*cp_n*II
                

 
   if(testmode .eq. 1) then

     JP_uv = 0.d0          
     JP_uv(1 : nn**3, 1 : nn**3) =  pen*PP 
     JP_uv(nn**3+1 : 2*nn**3, nn**3+1 : 2*nn**3) = pen*PP 
     JP_uv(2*nn**3+1 : 3*nn**3, 2*nn**3+1 : 3*nn**3) = pen*PP 
  
     JM_uv = 0.d0
     JM_uv(1 : nn**3, 1 : mm**3) = - pen*QQ 
     JM_uv(nn**3+1 : 2*nn**3, mm**3+1 : 2*mm**3) = - pen*QQ 
     JM_uv(2*nn**3+1 : 3*nn**3, 2*mm**3+1 : 3*mm**3) = - pen*QQ 
     
   endif





   deallocate(AA,BB,CC,PP)            
   deallocate(DD,EE,FF,QQ)                  
   deallocate(GG,HH,II) 


   end subroutine MAKE_LOC_MATRIX_DG
