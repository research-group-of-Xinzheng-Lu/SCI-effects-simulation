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

!> @brief Update the damping tensor with the Heun method. 
!! @author Ilario Mazzieri
!> @date November, 2014 
!> @version 1.0
!> @param[in] nnod_loc number of local nodes
!> @param[in] N_SLS number of standard linear solids
!> @param[in] frequency_range  frequency range for SLS model
!> @param[in] strain strain tensor
!> @param[in] deltat time step
!> @param[out] strain_visc viscoelastic tensor


      subroutine UPDATE_DAMPING_TENSOR(nnod_loc, N_SLS, frequency_range, &
                                  strain_visc, strain, deltat)
 
      
      implicit none
      						
 
      integer*4 :: nnod_loc, N_SLS, i
                         
      
      real*8 :: deltat, frequency_range(N_SLS), strain_visc(6*nnod_loc,N_SLS), &
                strain(6*nnod_loc), temp(6*nnod_loc, N_SLS)
 

      
      do i = 1, N_SLS
         !EXPLICIT EULER
         !temp(:,i) = strain_visc(:,i) + deltat*frequency_range(i)*(strain - strain_visc(:,i))
         !EXPLICIT HEUN
         temp(:,i) = strain_visc(:,i) + &
                     (deltat*frequency_range(i) - 0.5*(deltat*frequency_range(i))**2)*(strain - strain_visc(:,i))
         

      enddo
      
      strain_visc = temp;
 
 
      return
      
      end subroutine UPDATE_DAMPING_TENSOR
