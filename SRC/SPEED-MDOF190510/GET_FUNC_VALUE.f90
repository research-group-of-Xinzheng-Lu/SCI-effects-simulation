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

!> @brief Computes time evolution function.
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
!> @param[in] dist  distance form source point (for travelling load)
!> @param[in] vel  (constant) velocity of the travelling load
!> @param[out] GET_FUNC_VALUE value of the time function


      real*8 function GET_FUNC_VALUE(nb_fnc, type_fnc, ind_fnc, &
                                     data_fnc, nb_data_fnc, id_fnc, time, dist,vel, &
                                     Mdoftnum,MDOFforcetran,MDOFbid,MDOFbdirct,tagty)
            



      implicit none
      
      integer*4 :: nb_fnc,id_fnc,i,nb_data_fnc,Mdoftnum
      
      integer*4, dimension(nb_fnc) :: type_fnc
      integer*4, dimension(nb_fnc +1) :: ind_fnc
      
      real*8 :: PI,t_t0,t0,t1,v0,v1
      real*8 :: TAU,scaling,HDUR 
      real*8 :: amp, ps0, tplus, alpha,time,beta2,dist,vel
      
      real*8, dimension(nb_data_fnc) :: data_fnc
      real*8, dimension(6*Mdoftnum) :: MDOFforcetran  !Modified on 26th, Mar. to add moment.
      real*8 :: MDOFbid,MDOFbdirct !Modified on 10th, Apr. to reduce memory.
      integer*4 :: MDOFbidint,MDOFbdirctint !Modified on 10th, Apr. to reduce memory.
      integer*4,intent(in)::tagty
      
      GET_FUNC_VALUE = 0.0d0

      PI = 4.0d0 * datan(1.0d0)

      select case (type_fnc(id_fnc))
         
         case(0)
           GET_FUNC_VALUE = 1.0d0

         case(1)
           t_t0 = time - data_fnc(ind_fnc(id_fnc) +1) 
           GET_FUNC_VALUE = (1.0d0 - 2.0d0*data_fnc(ind_fnc(id_fnc))*t_t0*t_t0) &
                 * dexp(-1.0d0*data_fnc(ind_fnc(id_fnc))*t_t0*t_t0)
         
         case(2)
           PI = 4.0d0 * datan(1.0d0)
           t_t0 = time - data_fnc(ind_fnc(id_fnc) +1)
           GET_FUNC_VALUE = dcos(PI*data_fnc(ind_fnc(id_fnc))*t_t0) &
                 * dexp(-0.5d0*data_fnc(ind_fnc(id_fnc)) &
                 * data_fnc(ind_fnc(id_fnc))*t_t0*t_t0)
      
         case(3)
           do i = ind_fnc(id_fnc),ind_fnc(id_fnc+1) -3,2
              t0 = data_fnc(i);    t1 = data_fnc(i +2)
              v0 = data_fnc(i +1);  v1 = data_fnc(i +3)
              if ((time.ge.t0) .and. (time.le.t1))  then
                  GET_FUNC_VALUE = (v1 - v0) / (t1 - t0) * (time - t0)  + v0
                  return
              endif    
           enddo
         
         ! TIME SERIES ADDED BY TY  !!!!!!!!!!!!!ty!!!!!!!!!!!!!
         case(777)
         MDOFbidint=nint(MDOFbid)
         MDOFbdirctint=nint(MDOFbdirct)
         if(MDOFbidint.gt.0) then
         GET_FUNC_VALUE =MDOFforcetran(6*(MDOFbidint-1)+MDOFbdirctint)
         else
         GET_FUNC_VALUE=0
         end if
         ! TIME SERIES ADDED BY TY  !!!!!!!!!!!!!ty!!!!!!!!!!!!!
         case(773)
                  GET_FUNC_VALUE = data_fnc(ind_fnc(id_fnc) + tagty)

        !!!!!!!!!!!!!!!!!!!
         case(4)
           t_t0 = time - data_fnc(ind_fnc(id_fnc) +1)
           beta2 = data_fnc(ind_fnc(id_fnc))
           GET_FUNC_VALUE = 2.0d0*beta2*t_t0 &
                 * (-3.0d0 + 2.0d0*beta2*t_t0*t_t0) &
                 * dexp(-beta2*t_t0*t_t0)
        
         case(6)
           PI = 4.0d0 * datan(1.0d0)
           t_t0 = time - data_fnc(ind_fnc(id_fnc) +1)
           GET_FUNC_VALUE = (2.d0*PI*data_fnc(ind_fnc(id_fnc))*t_t0) &
                 * dexp(-0.5d0*4.d0*PI*PI*data_fnc(ind_fnc(id_fnc)) &
                 * data_fnc(ind_fnc(id_fnc))*t_t0*t_t0)

         case(12)
           !------------------------------------------------
           ! 12 - sigmf(t,[a c]) = amp*(1/(1+exp(-a*(t-c))))
           ! 
           !
           !   
           !  |                                                                                                
           !  |............*************************......amp      
           !  |          ** : 
           !  |     a   *   :
           !  |        *    :
           !  |      **     :   
           !  0******----------------------------------> Time
           !       |        |
           !       |____c___|                          |
           !       |        |                          |
           !                                       n*(1/f)+t0
           ! 
           t_t0 = time - data_fnc(ind_fnc(id_fnc) +2)
           GET_FUNC_VALUE = data_fnc(ind_fnc(id_fnc)) &
                           * (1/(1+exp(-data_fnc(ind_fnc(id_fnc) +1)*(t_t0))))

         case(13)
           ! - GRENOBLE BENCHMARK
           TAU = data_fnc(ind_fnc(id_fnc));  scaling = data_fnc(ind_fnc(id_fnc) +1)
           t0 = 2.0 * TAU;  HDUR = TAU/2.0;  t_t0 = time - t0
           GET_FUNC_VALUE = 0.5d0 * ( 1.0d0 + erf(scaling*(t_t0)/HDUR) )
 
         case(14)
           ! - SCEC BENCHMARK
           TAU = data_fnc(ind_fnc(id_fnc));  amp = data_fnc(ind_fnc(id_fnc) +1)
           t_t0 = time
           if (t_t0 .lt. 0.0d0) then
             GET_FUNC_VALUE = 0.0d0
           elseif (t_t0 .eq. 0.0d0) then
             GET_FUNC_VALUE = 0.5d0
           else
             GET_FUNC_VALUE = amp * (1.0d0 - (1 + t_t0/TAU)*exp(-t_t0/TAU))
           endif

         case(15)
           ! - EXPLOSION
           ps0 = data_fnc(ind_fnc(id_fnc)) 
           tplus = data_fnc(ind_fnc(id_fnc) +1); alpha = data_fnc(ind_fnc(id_fnc) +2)
           t_t0 = time - data_fnc(ind_fnc(id_fnc) +3)
           if (t_t0 .lt. 0.0d0) then
             GET_FUNC_VALUE = 0.0d0
           else
             GET_FUNC_VALUE = ps0 * (1 - (1 + time/tplus)*exp(-alpha*time/tplus))
           endif

         case(60,62)
           ! FUNCTION FOR G/G0
           do i = ind_fnc(id_fnc),ind_fnc(id_fnc+1) -3,2     
           t0 = data_fnc(i);     t1 = data_fnc(i +2)                    
              v0 = data_fnc(i +1);  v1 = data_fnc(i +3)                
    
              if (abs(time).le.data_fnc(ind_fnc(id_fnc))) then          
                GET_FUNC_VALUE = data_fnc(ind_fnc(id_fnc)+1)            
              elseif ((abs(time).ge.t0).and.(abs(time).le.t1)) then     
                GET_FUNC_VALUE = (v1 - v0) / (t1 - t0) * (abs(time) - t0)  + v0
             elseif (abs(time).ge.data_fnc(ind_fnc(id_fnc+1)-2)) then   
                GET_FUNC_VALUE = data_fnc(ind_fnc(id_fnc+1)-1)          
              endif                                                 
           enddo

         case(61,63)
           ! FUNCTION FOR DAMPING
           do i = ind_fnc(id_fnc),ind_fnc(id_fnc+1) -3,2      
              t0 = data_fnc(i);    t1 = data_fnc(i +2)                  
              v0 = data_fnc(i +1); v1 = data_fnc(i +3)                  

              if (abs(time).le.data_fnc(ind_fnc(id_fnc))) then          
                GET_FUNC_VALUE = data_fnc(ind_fnc(id_fnc)+1)            
              elseif ((abs(time).ge.t0).and.(abs(time).le.t1)) then    
                GET_FUNC_VALUE = (v1 - v0) / (t1 - t0) * (abs(time) - t0)  + v0 
              elseif (abs(time).ge.data_fnc(ind_fnc(id_fnc+1)-2)) then 
                GET_FUNC_VALUE = data_fnc(ind_fnc(id_fnc+1)-1)          
              endif                                                 
           enddo

         case(99)
           ! - CASHIMA1 BENCHMARK                                           
           TAU = data_fnc(ind_fnc(id_fnc));  amp = data_fnc(ind_fnc(id_fnc) +1)
           t_t0 = time                                                  
           GET_FUNC_VALUE = ( exp( - (((2.0d0 * 2.0d0*dasin(1.0d0))*1.5d0)&         
                            *(t_t0 - TAU)/amp)**2.0d0)&                
                            * cos(((2.0d0 * 2.0d0*dasin(1.0d0))*1.5d0)& 
                            *(t_t0 - TAU) + 2.0d0*dasin(1.0d0)/2.0d0))  
      
         case(100)
          ! - TESTMODE
          amp = data_fnc(ind_fnc(id_fnc))       
          GET_FUNC_VALUE = dsin(amp*PI*time)
          
          
         case(101) 
          GET_FUNC_VALUE = time
               
         
         case default
           GET_FUNC_VALUE = 0.d0
      
      end select

     
      return
      
      end function GET_FUNC_VALUE

