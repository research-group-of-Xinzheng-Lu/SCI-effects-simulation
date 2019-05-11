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

!> @brief Makes local spectral grid.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nn_loc number of local node
!> @param[in,out] xx_loc x-coordinate of GLL node 
!> @param[in,out] yy_loc y-coordinate of GLL node
!> @param[in,out] zz_loc z-coordinate of GLL node
!> @param[in] cs_nnz_loc length of cs_loc
!> @param[in] cs_loc  local spectral connectivity vector
!> @param[in] nm number of materials
!> @param[in] tm material labels
!> @param[in] sd polynomial degree vector
!> @param[in] loc_n_num  local node numeration
!> @param[in] ne_loc number of local elements
!> @param[out] alfa11 costant values for the bilinear map
!> @param[out] alfa12 costant values for the bilinear map
!> @param[out] alfa13 costant values for the bilinear map
!> @param[out] alfa21 costant values for the bilinear map
!> @param[out] alfa22 costant values for the bilinear map
!> @param[out] alfa23 costant values for the bilinear map
!> @param[out] alfa31 costant values for the bilinear map
!> @param[out] alfa32 costant values for the bilinear map
!> @param[out] alfa33 costant values for the bilinear map
!> @param[out] beta11 costant values for the bilinear map 
!> @param[out] beta12 costant values for the bilinear map
!> @param[out] beta13 costant values for the bilinear map
!> @param[out] beta21 costant values for the bilinear map
!> @param[out] beta22 costant values for the bilinear map
!> @param[out] beta23 costant values for the bilinear map
!> @param[out] beta31 costant values for the bilinear map
!> @param[out] beta32 costant values for the bilinear map
!> @param[out] beta33 costant values for the bilinear map
!> @param[out] gamma1 costant values for the bilinear map
!> @param[out] gamma2 costant values for the bilinear map
!> @param[out] gamma3 costant values for the bilinear map
!> @param[out] delta1 costant values for the bilinear map
!> @param[out] delta2 costant values for the bilinear map
!> @param[out] delta3 costant values for the bilinear map


      subroutine MAKE_SPX_GRID_LOC(nn_loc, xx_loc, yy_loc, zz_loc, &
                               cs_nnz_loc, cs_loc, nm, tm, sd, ne_loc, loc_n_num, &
                               alfa11,alfa12,alfa13, &
                               alfa21,alfa22,alfa23, &
                               alfa31,alfa32,alfa33, &
                               beta11,beta12,beta13, &
                               beta21,beta22,beta23, &
                               beta31,beta32,beta33, &
                               gamma1,gamma2,gamma3, &
                               delta1,delta2,delta3 )

      
      use speed_exit_codes

      implicit none
      
      integer*4 :: nn_loc, cs_nnz_loc,nm, ne_loc
      integer*4 :: n1,n2,n3,n4,n5,n6,n7,n8
      integer*4 :: ic1,ic2,ic3,ic4,ic5,ic6,ic7,ic8
      integer*4 :: im,ie,i,j,k,nn,ip, ic
      
      integer*4, dimension(nn_loc) :: loc_n_num      
      integer*4, dimension(0:cs_nnz_loc) :: cs_loc
      integer*4, dimension(nm) :: tm
      integer*4, dimension(nm) :: sd

      real*8 :: xp,yp,zp,jac
      real*8 :: x1,x2,x3,x4,x5,x6,x7,x8
      real*8 :: y1,y2,y3,y4,y5,y6,y7,y8
      real*8 :: z1,z2,z3,z4,z5,z6,z7,z8
      
      real*8, dimension(:), allocatable :: ct,ww
      real*8, dimension(ne_loc) :: alfa11,alfa12,alfa13
      real*8, dimension(ne_loc) :: alfa21,alfa22,alfa23
      real*8, dimension(ne_loc) :: alfa31,alfa32,alfa33
      real*8, dimension(ne_loc) :: beta11,beta12,beta13
      real*8, dimension(ne_loc) :: beta21,beta22,beta23
      real*8, dimension(ne_loc) :: beta31,beta32,beta33
      real*8, dimension(ne_loc) :: gamma1,gamma2,gamma3
      real*8, dimension(ne_loc) :: delta1,delta2,delta3  
      real*8, dimension(nn_loc), intent(inout) :: xx_loc,yy_loc,zz_loc
      
      real*8, dimension(:,:), allocatable :: dd
  
      nn = 2
      allocate(ct(nn),ww(nn),dd(nn,nn))
      call MAKE_LGL_NW(nn,ct,ww,dd)
      
      do im = 1,nm
         if ((sd(im) +1).ne.nn) then
            deallocate(ct,ww,dd)
            nn = sd(im) +1
            allocate(ct(nn),ww(nn),dd(nn,nn))
            call MAKE_LGL_NW(nn,ct,ww,dd)
         endif
         
         do ie = 1,ne_loc
            if (cs_loc(cs_loc(ie -1) +0).eq.tm(im)) then
               n1 = nn*nn*(1 -1) +nn*(1 -1) +1
               n2 = nn*nn*(1 -1) +nn*(1 -1) +nn
               n3 = nn*nn*(1 -1) +nn*(nn -1) +nn
               n4 = nn*nn*(1 -1) +nn*(nn -1) +1
               n5 = nn*nn*(nn -1) +nn*(1 -1) +1
               n6 = nn*nn*(nn -1) +nn*(1 -1) +nn
               n7 = nn*nn*(nn -1) +nn*(nn -1) +nn
               n8 = nn*nn*(nn -1) +nn*(nn -1) +1
               
               ic1 = cs_loc(cs_loc(ie -1) +n1)
               ic2 = cs_loc(cs_loc(ie -1) +n2)
               ic3 = cs_loc(cs_loc(ie -1) +n3)
               ic4 = cs_loc(cs_loc(ie -1) +n4)
               ic5 = cs_loc(cs_loc(ie -1) +n5)
               ic6 = cs_loc(cs_loc(ie -1) +n6)
               ic7 = cs_loc(cs_loc(ie -1) +n7)
               ic8 = cs_loc(cs_loc(ie -1) +n8)
               
                                                                     
               x1 = xx_loc(ic1)
               y1 = yy_loc(ic1)
               z1 = zz_loc(ic1)
               
               x2 = xx_loc(ic2)
               y2 = yy_loc(ic2)
               z2 = zz_loc(ic2)
               
               x3 = xx_loc(ic3)
               y3 = yy_loc(ic3)
               z3 = zz_loc(ic3)
               
               x4 = xx_loc(ic4)
               y4 = yy_loc(ic4)
               z4 = zz_loc(ic4)
               
               x5 = xx_loc(ic5)
               y5 = yy_loc(ic5)
               z5 = zz_loc(ic5)
               
               x6 = xx_loc(ic6)
               y6 = yy_loc(ic6)
               z6 = zz_loc(ic6)
               
               x7 = xx_loc(ic7)
               y7 = yy_loc(ic7)
               z7 = zz_loc(ic7)
               
               x8 = xx_loc(ic8)
               y8 = yy_loc(ic8)
               z8 = zz_loc(ic8)
                              
                              
               alfa11(ie) = 0.125d0*(-x1 +x2 +x3 -x4 -x5 +x6 +x7 -x8)  !1/8 * e1
               alfa21(ie) = 0.125d0*(-y1 +y2 +y3 -y4 -y5 +y6 +y7 -y8)  !1/8 * e2 
               alfa31(ie) = 0.125d0*(-z1 +z2 +z3 -z4 -z5 +z6 +z7 -z8)  !1/8 * e3

               alfa12(ie) = 0.125d0*(-x1 -x2 +x3 +x4 -x5 -x6 +x7 +x8)  !1/8 * f1
               alfa22(ie) = 0.125d0*(-y1 -y2 +y3 +y4 -y5 -y6 +y7 +y8)  !1/8 * f2
               alfa32(ie) = 0.125d0*(-z1 -z2 +z3 +z4 -z5 -z6 +z7 +z8)  !1/8 * f3

               alfa13(ie) = 0.125d0*(-x1 -x2 -x3 -x4 +x5 +x6 +x7 +x8)  !1/8 * g1
               alfa23(ie) = 0.125d0*(-y1 -y2 -y3 -y4 +y5 +y6 +y7 +y8)  !1/8 * g2
               alfa33(ie) = 0.125d0*(-z1 -z2 -z3 -z4 +z5 +z6 +z7 +z8)  !1/8 * g3
               
               beta11(ie) = 0.125d0*(+x1 +x2 -x3 -x4 -x5 -x6 +x7 +x8)  !1/8 * c1
               beta21(ie) = 0.125d0*(+y1 +y2 -y3 -y4 -y5 -y6 +y7 +y8)  !1/8 * c2
               beta31(ie) = 0.125d0*(+z1 +z2 -z3 -z4 -z5 -z6 +z7 +z8)  !1/8 * c3

               beta12(ie) = 0.125d0*(+x1 -x2 -x3 +x4 -x5 +x6 +x7 -x8)  !1/8 * d1
               beta22(ie) = 0.125d0*(+y1 -y2 -y3 +y4 -y5 +y6 +y7 -y8)  !1/8 * d2
               beta32(ie) = 0.125d0*(+z1 -z2 -z3 +z4 -z5 +z6 +z7 -z8)  !1/8 * d3

               beta13(ie) = 0.125d0*(+x1 -x2 +x3 -x4 +x5 -x6 +x7 -x8)  !1/8 * b1
               beta23(ie) = 0.125d0*(+y1 -y2 +y3 -y4 +y5 -y6 +y7 -y8)  !1/8 * b2
               beta33(ie) = 0.125d0*(+z1 -z2 +z3 -z4 +z5 -z6 +z7 -z8)  !1/8 * b3
               
               gamma1(ie) = 0.125d0*(-x1 +x2 -x3 +x4 +x5 -x6 +x7 -x8)  !1/8 * a1
               gamma2(ie) = 0.125d0*(-y1 +y2 -y3 +y4 +y5 -y6 +y7 -y8)  !1/8 * a2
               gamma3(ie) = 0.125d0*(-z1 +z2 -z3 +z4 +z5 -z6 +z7 -z8)  !1/8 * a3
                
               delta1(ie) = 0.125d0*(+x1 +x2 +x3 +x4 +x5 +x6 +x7 +x8)  !1/8 * h1
               delta2(ie) = 0.125d0*(+y1 +y2 +y3 +y4 +y5 +y6 +y7 +y8)  !1/8 * h2
               delta3(ie) = 0.125d0*(+z1 +z2 +z3 +z4 +z5 +z6 +z7 +z8)  !1/8 * h3
               
               
               jac = alfa31(ie)*(alfa12(ie)*alfa23(ie) -alfa13(ie)*alfa22(ie)) &
                   - alfa32(ie)*(alfa11(ie)*alfa23(ie) -alfa13(ie)*alfa21(ie)) &
                   + alfa33(ie)*(alfa11(ie)*alfa22(ie) -alfa12(ie)*alfa21(ie))
               

               if (jac.le.0.d0) then
                  write(*,*)'Error ! Orientation non-conforming !'
                  write(*,*)'(element ',ie,' is clockwise oriented)'
                  call EXIT(EXIT_ELEM_ORIENT)
               endif
               
               do k = 1,nn
                  do j = 1,nn
                     do i = 1,nn
                        xp = alfa11(ie)*ct(i) + alfa12(ie)*ct(j) &
                           + alfa13(ie)*ct(k) + beta11(ie)*ct(j)*ct(k) &
                           + beta12(ie)*ct(i)*ct(k) + beta13(ie)*ct(i)*ct(j) &
                           + gamma1(ie)*ct(i)*ct(j)*ct(k) + delta1(ie)
                        
                        yp = alfa21(ie)*ct(i) + alfa22(ie)*ct(j) &
                           + alfa23(ie)*ct(k) + beta21(ie)*ct(j)*ct(k) &
                           + beta22(ie)*ct(i)*ct(k) + beta23(ie)*ct(i)*ct(j) &
                           + gamma2(ie)*ct(i)*ct(j)*ct(k) + delta2(ie)
                        
                        zp = alfa31(ie)*ct(i) + alfa32(ie)*ct(j) &
                           + alfa33(ie)*ct(k) + beta31(ie)*ct(j)*ct(k) &
                           + beta32(ie)*ct(i)*ct(k) + beta33(ie)*ct(i)*ct(j) &
                           + gamma3(ie)*ct(i)*ct(j)*ct(k) + delta3(ie)
                        
                        ip = nn*nn*(k -1) +nn*(j -1) +i

                        xx_loc(cs_loc(cs_loc(ie -1) + ip)) = xp
                        yy_loc(cs_loc(cs_loc(ie -1) + ip)) = yp
                        zz_loc(cs_loc(cs_loc(ie -1) + ip)) = zp

                     enddo
                  enddo
               enddo
               
            endif
         enddo
      enddo
      
      return
      
      end subroutine MAKE_SPX_GRID_LOC

