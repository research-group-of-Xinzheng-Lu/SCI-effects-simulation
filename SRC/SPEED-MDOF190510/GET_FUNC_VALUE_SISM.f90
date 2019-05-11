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


!> @brief Computes time evolution function for seismic source.
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] nb_fnc number of functions
!> @param[in] type_fnc function type
!> @param[in] ind_fnc indices for the data 
!> @param[in] nb_data_fnc number of data for each function
!> @param[in] data_fnc data for the calculation (depending on type_fnc) 
!> @param[in] id_fnc  number of the function
!> @param[in] time  instant time
!> @param[in] t0_delay  delay for activating seismic fault
!> @param[in] tau_new  rising time
!> @param[out] GET_FUNC_VALUE_SISM value of the time function

      real*8 function GET_FUNC_VALUE_SISM(nb_fnc, type_fnc, ind_fnc,&
                                          data_fnc,nb_data_fnc, id_fnc, time, t0_delay, tau_new) 
      
      implicit none
      
      integer*4 :: nb_fnc,id_fnc,i,nb_data_fnc
      integer*4 :: np0,np1,np2,np_current
      
      integer*4, dimension(nb_fnc) :: type_fnc
      integer*4, dimension(nb_fnc +1) :: ind_fnc

      real*8 :: val,PI,t_t0,t0,t1,v0,v1
      real*8 :: TAU, scaling, HDUR 
      real*8 :: amp 
      real*8 :: tau_new
      real*8 :: time,beta2,t0_delay 
      real*8 :: dt,Tr,Te,Tp,psv,sum1,ts,svi
      
      real*8, dimension(nb_data_fnc) :: data_fnc
      
      real*8, dimension(:), allocatable :: svf, int_svf, integ_svf
      
      GET_FUNC_VALUE_SISM = 0.0d0
      pi = 4.0d0 * datan(1.0d0)
      
      select case (type_fnc(id_fnc))
         
         case(0)
           GET_FUNC_VALUE_SISM = 1.d0

         case(1)  ! - Ricker "beta" type           
           t_t0 = time - data_fnc(ind_fnc(id_fnc) +1) - t0_delay
           GET_FUNC_VALUE_SISM = (1.0d0 - 2.0d0*data_fnc(ind_fnc(id_fnc))*t_t0*t_t0) &
                                 * dexp(-1.0d0*data_fnc(ind_fnc(id_fnc))*t_t0*t_t0)
          
         case(2)
           PI = 4.0d0 * datan(1.0d0)
           t_t0 = time - data_fnc(ind_fnc(id_fnc) +1) - t0_delay
           GET_FUNC_VALUE_SISM = dcos(PI*data_fnc(ind_fnc(id_fnc))*t_t0) &
                                 * dexp(-0.5d0*data_fnc(ind_fnc(id_fnc)) &
                                 * data_fnc(ind_fnc(id_fnc))*t_t0*t_t0)          
            
         case(3)  
            do i = ind_fnc(id_fnc),ind_fnc(id_fnc+1) -3,2
               t0 = data_fnc(i) - t0_delay
               t1 = data_fnc(i +2) - t0_delay
               v0 = data_fnc(i +1)
               v1 = data_fnc(i +3)
               if ((time.ge.t0).and.(time.le.t1)) &
                  GET_FUNC_VALUE_SISM = (v1 - v0) / (t1 - t0) * (time - t0)  + v0
            enddo
            
         case(4) ! - First derivative of the Ricker ( d(Ricker(t))/dt )
           t_t0 = time - data_fnc(ind_fnc(id_fnc) +1) - t0_delay
           beta2=data_fnc(ind_fnc(id_fnc))
           GET_FUNC_VALUE_SISM = 2.0d0*beta2*t_t0 &
                                 * (-3.0d0 + 2.0d0*beta2*t_t0*t_t0) &
                                 * dexp(-beta2*t_t0*t_t0)
         
         case(6)             
           PI = 4.0d0 * datan(1.0d0)
           t_t0 = time - data_fnc(ind_fnc(id_fnc) +1) - t0_delay
           GET_FUNC_VALUE_SISM = (2.d0*PI*data_fnc(ind_fnc(id_fnc))*t_t0) &
              * dexp(-0.5d0*4.d0*PI*PI*data_fnc(ind_fnc(id_fnc)) &
              * data_fnc(ind_fnc(id_fnc))*t_t0*t_t0)
           write(33,*)time,GET_FUNC_VALUE_SISM
         
         case(13) ! - GRENOBLE BENCHMARK
           TAU = data_fnc(ind_fnc(id_fnc)) 
           scaling = data_fnc(ind_fnc(id_fnc) +1)
           t0 = 2.0 * TAU
           HDUR = TAU/2.0

           t_t0 = time - t0 - t0_delay
           GET_FUNC_VALUE_SISM = 0.5d0 * ( 1.0d0 + erf(scaling*(t_t0)/HDUR) )
         
         case(14) ! - SCEC BENCHMARK
           TAU = data_fnc(ind_fnc(id_fnc)) 
           amp = data_fnc(ind_fnc(id_fnc) +1)
           t_t0 = time - t0_delay

           !val = amp * (1 - (1 + t_t0/TAU)*exp(-t_t0/TAU))
           if (t_t0.lt.0.0d0) then
             GET_FUNC_VALUE_SISM = 0.0d0
           elseif (t_t0.eq.0.0d0) then
             GET_FUNC_VALUE_SISM = 0.5d0
           else
             GET_FUNC_VALUE_SISM = amp * (1.0d0 - (1 + t_t0/TAU)*exp(-t_t0/TAU))
           endif         
         
         case(50) ! - GRENOBLE BENCHMARK - NEW TAU
           TAU =  tau_new
           scaling = data_fnc(ind_fnc(id_fnc) +1)
           t0 = 2.0 * TAU
           HDUR = TAU/2.0
         
           t_t0 = time - t0 - t0_delay
           GET_FUNC_VALUE_SISM = 0.5d0 * ( 1.0d0 + erf(scaling*(t_t0)/HDUR) )


         case(55) ! - ARCHULETA FUNCTION
            ! time = current time
            ! t0_delay = rupture time
            ! tau_new  = rise time
           !if ( time .lt. t0_delay + tau_new .and. time .gt. t0_delay) then 
           !  write(*,*) time, t0_delay, t0_delay + tau_new
           !endif 
            
           if (time .le. t0_delay) then
               GET_FUNC_VALUE_SISM = 0.d0;
           elseif ( time .ge. t0_delay + tau_new) then 
               GET_FUNC_VALUE_SISM = 1.d0;
           else
               dt = 0.01;
               Tr = tau_new;
               Te = 0.8*Tr;
               Tp = 0.2*Tr;
               
               np0 = int(Tp/dt+1.0)
               np1 = int(Te/dt+1.0)
               np2 = int(Tr/dt+1.0)
               np_current = int((time-t0_delay)/dt+1.0)
               
               if (np_current .gt. np2) then 
                  write(*,*) 'Error in Archuleta function'
                  write(*,*) np2, np_current
                  stop
               endif
               
               psv = sqrt(1.+100./(np0*dt))
               sum1 =0.0
               allocate(svf(1:np2)); 
              
               do i = 1, np0
                  ts = (i-1)*dt
                  svi = ts*psv/Tp*sin(0.5*pi/Tp*ts)
                  svf(i) = svi
                  sum1 = sum1+svi
               enddo
               do i = np0 + 1, np1
                  ts = (i-1)*dt
                  svi = sqrt(1.+100./ts)
                  svf(i) = svi
                  sum1 = sum1 + svi
               enddo
               do i = np1 + 1, np2
                  ts = (i-1)*dt
                  svi = sqrt(1.+100./ts)*sin((np2-i)*dt*pi*0.5/(Tr-Te))
                  svf(i) = svi
                  sum1 = sum1 + svi
              enddo
              sum1 = sum1 * dt
 
              !do i = 1, np2
              !   svf(i) = svf(i)/sum1
              !enddo
              svf = svf/sum1;
              
              allocate(int_svf(1:np2)); int_svf = 0.d0;  int_svf(1) = 0;
              allocate(integ_svf(1:np2)); integ_svf = 0.d0;  integ_svf(1) = 0;
              
              do i = 2, np2   
                  int_svf(i) =  0.5*(svf(i)+svf(i-1))*dt;
              !enddo
              !do i = 2, np2   
                integ_svf(i) =  sum(int_svf(2:i));
              enddo
              
              GET_FUNC_VALUE_SISM = integ_svf(np_current)
               if (GET_FUNC_VALUE_SISM .gt. 1.05d0) then 
                  write(*,*) 'Error in Archuleta slip function '
                  write(*,*) GET_FUNC_VALUE_SISM
                  stop
               endif
              deallocate(svf,int_svf,integ_svf)
              !write(*,*) GET_FUNC_VALUE_SISM
              !read(*,*)

           endif
           
              

         case(99) ! - CASHIMA 1 
           TAU = data_fnc(ind_fnc(id_fnc))
           amp = data_fnc(ind_fnc(id_fnc) +1) 
           t_t0 = time
           GET_FUNC_VALUE_SISM = ( exp( - (((2.0d0 * 2.0d0*dasin(1.0d0))*1.5d0) & 
                                 *(t_t0 - TAU)/amp)**2.0d0) * cos(((2.0d0 * 2.0d0*dasin(1.0d0))*1.5d0) &
                                 *(t_t0 - TAU) + 2.0d0*dasin(1.0d0)/2.0d0))         
                                 
         case default
           GET_FUNC_VALUE_SISM = 0.d0
      
      end select
      
      
      return
      
      end function GET_FUNC_VALUE_SISM

