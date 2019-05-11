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


!> @brief Computes max values for Peak Ground Map 
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in]  nmonitpgn number of monitored points
!> @param[in]  elem_mpgm elem_mpgm(i) cotains the index of the element
!! to which the i-th node belongs to
!> @param[in]  ne_loc number of local element
!> @param[in]  local_el_num local element numeration
!> @param[in]  nm number of materials
!> @param[in]  sdeg_mat polynomial degree vector
!> @param[in]  cs_nnz_loc lenght of cs_loc
!> @param[in]  cs_loc spectral connectivity vector
!> @param[in]  u1 displacement vector at time t_n
!> @param[in]  u0 displacement vector at time t_{n-1}
!> @param[in]  u_1 displacement vector at time t_{n-2}
!> @param[in]  omega rotation tensor at time t_n
!> @param[in]  nnod_loc number of local nodes
!> @param[in]  xr_mpgm x-coor for monitor mpgm
!> @param[in]  yr_mpgm y-coor for monitor mpgm
!> @param[in]  zr_mpgm z-coor for monitor mpgm
!> @param[out] max_u max displacement value
!> @param[out] max_v max velocity value
!> @param[out] max_a max acceleration value
!> @param[out] max_0 max rotation values 

     subroutine GET_MAX_VALUES(nmonitpgm, elem_mpgm, local_el_num, ne_loc, &
                              sdeg_mat, nm, cs_loc, cs_nnz_loc, u1, u0, u_1,& 
                              omega, nnod_loc, &
                              xr_mpgm,yr_mpgm,zr_mpgm,&
                              max_u, max_v, max_a, max_o)


      implicit none      

      integer*4 :: nnod_loc, ne_loc, cs_nnz_loc, nm
      integer*4 :: nn,  im, imon, iaz
      integer*4 :: i, j, k, is, in, id
      integer*4 :: ie, ielem
      integer*4, dimension(0:cs_nnz_loc) :: cs_loc
      integer*4, dimension(nm) :: sdeg_mat
      integer*4, dimension(ne_loc) :: local_el_num

      real*8 :: uxm, uym, uzm
      real*8 :: vxm, vym, vzm
      real*8 :: axm, aym, azm


      real*8, dimension(:), allocatable :: ct,ww
      real*8, dimension(:,:), allocatable :: dd
      real*8, dimension(:,:,:), allocatable :: ux_el, uy_el, uz_el

      real*8, dimension(3*nnod_loc) :: u_1,u0,u1

      integer*4 :: imonpgm, nmonitpgm,ndt_monitpgm
      integer*4, dimension(nmonitpgm) :: node_mpgm
      integer*4, dimension(nmonitpgm) :: elem_mpgm
      real*8 :: tmp, dt
      real*8 :: rotang_monitpgm
      real*8 :: upm,unm,vpm,vnm,apm,anm
      real*8, dimension(nmonitpgm) :: xr_mpgm,yr_mpgm,zr_mpgm
      real*8, dimension(nmonitpgm,9) :: max_u,max_v,max_a
      real*8, dimension(nmonitpgm,3) :: max_o
      
      
      real*8 :: variable1m,variable2m,variable3m
      real*8 :: variable4m,variable5m,variable6m

      real*8, dimension(3*nnod_loc) :: omega

      real*8, dimension(:,:,:), allocatable :: variable1_el,variable2_el,variable3_el
      real*8, dimension(:,:,:), allocatable :: variable4_el,variable5_el,variable6_el

                     

     upm = 0.0d0; unm = 0.0d0; vpm = 0.0d0; vnm = 0.0d0; apm = 0.0d0; anm = 0.0d0


      do imonpgm = 1,nmonitpgm
           
            ielem = elem_mpgm(imonpgm)
            call GET_INDLOC_FROM_INDGLO(local_el_num, ne_loc, ielem, ie)

            if (ie .ne. 0) then
               
               im = cs_loc(cs_loc(ie -1) +0)
               nn = sdeg_mat(im) +1
               allocate(ct(nn),ww(nn),dd(nn,nn))
               allocate(ux_el(nn,nn,nn),uy_el(nn,nn,nn),uz_el(nn,nn,nn))
               allocate(variable1_el(nn,nn,nn),variable2_el(nn,nn,nn),variable3_el(nn,nn,nn))	
               allocate(variable4_el(nn,nn,nn),variable5_el(nn,nn,nn),variable6_el(nn,nn,nn))	
               call MAKE_LGL_NW(nn,ct,ww,dd)

               do k = 1,nn
                  do j = 1,nn
                     do i = 1,nn
                        is = nn*nn*(k -1) +nn*(j -1) +i
                        in = cs_loc(cs_loc(ie -1) + is)

                        iaz = 3*(in -1) +1
                        ux_el(i,j,k) = u1(iaz)
                        iaz = 3*(in -1) +2
                        uy_el(i,j,k) = u1(iaz)
                        iaz = 3*(in -1) +3
                        uz_el(i,j,k) = u1(iaz)
                     enddo
                  enddo
               enddo

               call GET_MONITOR_VALUE(nn,ct,ux_el,&
                            xr_mpgm(imonpgm),yr_mpgm(imonpgm),zr_mpgm(imonpgm),uxm)
               call GET_MONITOR_VALUE(nn,ct,uy_el,&
                            xr_mpgm(imonpgm),yr_mpgm(imonpgm),zr_mpgm(imonpgm),uym)
               call GET_MONITOR_VALUE(nn,ct,uz_el,&
                            xr_mpgm(imonpgm),yr_mpgm(imonpgm),zr_mpgm(imonpgm),uzm)

               do k = 1,nn
                  do j = 1,nn
                     do i = 1,nn
                        is = nn*nn*(k -1) +nn*(j -1) +i
                        in = cs_loc(cs_loc(ie -1) + is)

                        iaz = 3*(in -1) +1
                        ux_el(i,j,k) = (3.0d0*u1(iaz) - 4.0d0*u0(iaz) + 1.0d0*u_1(iaz)) / (2.0d0*dt)
                        iaz = 3*(in -1) +2
                        uy_el(i,j,k) = (3.0d0*u1(iaz) - 4.0d0*u0(iaz) + 1.0d0*u_1(iaz)) / (2.0d0*dt)
                        iaz = 3*(in -1) +3
                        uz_el(i,j,k) = (3.0d0*u1(iaz) - 4.0d0*u0(iaz) + 1.0d0*u_1(iaz)) / (2.0d0*dt)
                     enddo
                  enddo
               enddo
               
               call GET_MONITOR_VALUE(nn,ct,ux_el,&
                    xr_mpgm(imonpgm),yr_mpgm(imonpgm),zr_mpgm(imonpgm),vxm)               
               call GET_MONITOR_VALUE(nn,ct,uy_el,&
                    xr_mpgm(imonpgm),yr_mpgm(imonpgm),zr_mpgm(imonpgm),vym)               
               call GET_MONITOR_VALUE(nn,ct,uz_el,&
                    xr_mpgm(imonpgm),yr_mpgm(imonpgm),zr_mpgm(imonpgm),vzm) 

               do k = 1,nn
                  do j = 1,nn
                     do i = 1,nn
                        is = nn*nn*(k -1) +nn*(j -1) +i
                        in = cs_loc(cs_loc(ie -1) + is)

                        iaz = 3*(in -1) +1
                        ux_el(i,j,k) = (u1(iaz) -2.0*u0(iaz) +u_1(iaz)) / dt**2
                        iaz = 3*(in -1) +2
                        uy_el(i,j,k) = (u1(iaz) -2.0*u0(iaz) +u_1(iaz)) / dt**2 
                        iaz = 3*(in -1) +3
                        uz_el(i,j,k) = (u1(iaz) -2.0*u0(iaz) +u_1(iaz)) / dt**2
                     enddo
                  enddo
               enddo
               
               call GET_MONITOR_VALUE(nn,ct,ux_el,&
                            xr_mpgm(imonpgm),yr_mpgm(imonpgm),zr_mpgm(imonpgm),axm)               
               call GET_MONITOR_VALUE(nn,ct,uy_el,&
                            xr_mpgm(imonpgm),yr_mpgm(imonpgm),zr_mpgm(imonpgm),aym)               
               call GET_MONITOR_VALUE(nn,ct,uz_el,&
                           xr_mpgm(imonpgm),yr_mpgm(imonpgm),zr_mpgm(imonpgm),azm)

               if (dabs(uxm).lt.1.0e-30) uxm = 0.0e+00
               if (dabs(uym).lt.1.0e-30) uym = 0.0e+00
               if (dabs(uzm).lt.1.0e-30) uzm = 0.0e+00
               if (dabs(vxm).lt.1.0e-30) vxm = 0.0e+00
               if (dabs(vym).lt.1.0e-30) vym = 0.0e+00
               if (dabs(vzm).lt.1.0e-30) vzm = 0.0e+00
               if (dabs(axm).lt.1.0e-30) axm = 0.0e+00
               if (dabs(aym).lt.1.0e-30) aym = 0.0e+00
               if (dabs(azm).lt.1.0e-30) azm = 0.0e+00

               if (rotang_monitpgm.ne.0.0d0) then
 
                  ! upm = Displacement fault Parallel Max
                  ! unm = Displacement fault Normal Max

                  upm = uxm * cos(rotang_monitpgm) + uym * sin(rotang_monitpgm)	
                  unm = -uxm * sin(rotang_monitpgm) + uym * cos(rotang_monitpgm)
                  vpm = vxm * cos(rotang_monitpgm) + vym * sin(rotang_monitpgm)	
                  vnm = -vxm * sin(rotang_monitpgm) + vym * cos(rotang_monitpgm)
                  apm = axm * cos(rotang_monitpgm) + aym * sin(rotang_monitpgm)
                  anm = -axm * sin(rotang_monitpgm) + aym * cos(rotang_monitpgm)

                  uxm = upm
                  uym = unm
                  vxm = vpm
                  vym = vnm
                  axm = apm
                  aym = anm
                endif

              if (dabs(uxm).gt.max_u(imonpgm,1)) max_u(imonpgm,1) = dabs(uxm)
              if (dabs(uym).gt.max_u(imonpgm,2)) max_u(imonpgm,2) = dabs(uym)
              if (dabs(uzm).gt.max_u(imonpgm,3)) max_u(imonpgm,3) = dabs(uzm)
    
              tmp = dsqrt(dabs(uxm*uym))
              if (tmp.gt.max_u(imonpgm,4)) max_u(imonpgm,4) = tmp
   
              tmp = ((dabs(uxm)+dabs(uym))/2)
              if (tmp.gt.max_u(imonpgm,5)) max_u(imonpgm,5) = tmp
    
              tmp = (dabs(uxm*uym*uzm))**(0.3333333)
              if (tmp.gt.max_u(imonpgm,6)) max_u(imonpgm,6) = tmp
   
              tmp = ((dabs(uxm)+dabs(uym)+dabs(uzm))/3)
              if (tmp.gt.max_u(imonpgm,7)) max_u(imonpgm,7) = tmp

              tmp = dsqrt((uxm)**2 + (uym)**2)
              if (tmp.gt.max_u(imonpgm,8)) max_u(imonpgm,8) = tmp

              tmp = dsqrt((uxm)**2 + (uym)**2 + (uzm)**2)
              if (tmp.gt.max_u(imonpgm,9)) max_u(imonpgm,9) = tmp


              if (dabs(vxm).gt.max_v(imonpgm,1)) max_v(imonpgm,1) = dabs(vxm)
              if (dabs(vym).gt.max_v(imonpgm,2)) max_v(imonpgm,2) = dabs(vym)
              if (dabs(vzm).gt.max_v(imonpgm,3)) max_v(imonpgm,3) = dabs(vzm)
   
              tmp = dsqrt(dabs(vxm*vym))
              if (tmp.gt.max_v(imonpgm,4)) max_v(imonpgm,4) = tmp
   
              tmp = ((dabs(vxm)+dabs(vym))/2)
              if (tmp.gt.max_v(imonpgm,5)) max_v(imonpgm,5) = tmp
   
              tmp = (dabs(vxm*vym*vzm))**(0.3333333)
              if (tmp.gt.max_v(imonpgm,6)) max_v(imonpgm,6) = tmp
   
              tmp = ((dabs(vxm)+dabs(vym)+dabs(vzm))/3)
              if (tmp.gt.max_v(imonpgm,7)) max_v(imonpgm,7) = tmp

              tmp = dsqrt((vxm)**2 + (vym)**2)
              if (tmp.gt.max_v(imonpgm,8)) max_v(imonpgm,8) = tmp

              tmp = dsqrt((vxm)**2 + (vym)**2 + (vzm)**2)
              if (tmp.gt.max_v(imonpgm,9)) max_v(imonpgm,9) = tmp

              if (dabs(axm).gt.max_a(imonpgm,1)) max_a(imonpgm,1) = dabs(axm)
              if (dabs(aym).gt.max_a(imonpgm,2)) max_a(imonpgm,2) = dabs(aym)
              if (dabs(azm).gt.max_a(imonpgm,3)) max_a(imonpgm,3) = dabs(azm)
   
              tmp = dsqrt(dabs(axm*aym))
              if (tmp.gt.max_a(imonpgm,4)) max_a(imonpgm,4) = tmp
    
              tmp = ((dabs(axm)+dabs(aym))/2)
              if (tmp.gt.max_a(imonpgm,5)) max_a(imonpgm,5) = tmp
    
              tmp = (dabs(axm*aym*azm))**(0.3333333)
              if (tmp.gt.max_a(imonpgm,6)) max_a(imonpgm,6) = tmp
   
              tmp = ((dabs(axm)+dabs(aym)+dabs(azm))/3)
              if (tmp.gt.max_a(imonpgm,7)) max_a(imonpgm,7) = tmp

              tmp = dsqrt((axm)**2 + (aym)**2)
              if (tmp.gt.max_a(imonpgm,8)) max_a(imonpgm,8) = tmp

              tmp = dsqrt((axm)**2 + (aym)**2 + (azm)**2)
              if (tmp.gt.max_a(imonpgm,9)) max_a(imonpgm,9) = tmp
               
!********************************************************
! OMEGA
!********************************************************

              do k = 1,nn
                  do j = 1,nn
                     do i = 1,nn
                        is = nn*nn*(k -1) +nn*(j -1) +i
                        in = cs_loc(cs_loc(ie -1) + is)

                        iaz = 3*(in -1) +1
                        variable1_el(i,j,k) = omega(iaz)
                        iaz = 3*(in -1) +2
                        variable2_el(i,j,k) = omega(iaz)
                        iaz = 3*(in -1) +3
                        variable3_el(i,j,k) = omega(iaz)
                     enddo
                  enddo
               enddo

             
               call GET_MONITOR_VALUE(nn,ct,variable1_el,&
                    xr_mpgm(imonpgm),yr_mpgm(imonpgm),zr_mpgm(imonpgm),variable1m)
               
               call GET_MONITOR_VALUE(nn,ct,variable2_el,&
                    xr_mpgm(imonpgm),yr_mpgm(imonpgm),zr_mpgm(imonpgm),variable2m)
               
               call GET_MONITOR_VALUE(nn,ct,variable3_el,&
                    xr_mpgm(imonpgm),yr_mpgm(imonpgm),zr_mpgm(imonpgm),variable3m)

       if (dabs(variable1m).lt.1.0e-30) variable1m = 0.0e+00
       if (dabs(variable2m).lt.1.0e-30) variable2m = 0.0e+00
       if (dabs(variable3m).lt.1.0e-30) variable3m = 0.0e+00
       if (dabs(variable1m).gt.max_o(imonpgm,1)) max_o(imonpgm,1) = dabs(variable1m)
       if (dabs(variable2m).gt.max_o(imonpgm,2)) max_o(imonpgm,2) = dabs(variable2m)
       if (dabs(variable3m).gt.max_o(imonpgm,3)) max_o(imonpgm,3) = dabs(variable3m)
               !-------------------------------------------------------------
               
       deallocate(ct,ww,dd)
       deallocate(ux_el,uy_el,uz_el)
       deallocate(variable1_el,variable2_el,variable3_el)
       deallocate(variable4_el,variable5_el,variable6_el)


            endif
          enddo

       end subroutine GET_MAX_VALUES
