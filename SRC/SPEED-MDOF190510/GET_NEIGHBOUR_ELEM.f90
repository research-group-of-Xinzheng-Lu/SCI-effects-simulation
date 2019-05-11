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

!> @brief Computes neighbouring element. 
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] omega_m omega_m(:,0)  = block identification for the neighbouring element,
!!                         omega_m(:,1)  = neighbouring element number,
!!                         omega_m(:,2)  = neighbouring element face
!> @param[in] npt number of quadrature points lying on the neighbouring element
!> @param[in] id_mpi number of process mpi
!> @param[in] n_qp  number of quadrature points
!> @param[out] omega_m(:,3)  label identifying the neighbournig element for each quadrature point 
!> @param[out] n  total neighbouring elements

     subroutine GET_NEIGHBOUR_ELEM(omega_m, npt, n, id_mpi, n_qp)
          
     implicit none
     
     integer*4 :: i, lett, j     
     integer*4, intent(in) :: npt, id_mpi, n_qp
     integer*4, intent(out)  :: n

     integer*4, dimension(n_qp) :: omega_m_copy
     integer*4, dimension(n_qp,0:3), intent(inout) :: omega_m 
     
     
     n = 0
     i = 1
     omega_m_copy = omega_m(:,1)
     do while (i .le. npt)
        lett = omega_m_copy(i)
        
        if(lett .ne. 0) then
           n = n + 1
           omega_m(i,3) = n
        endif
        
        do j = i + 1, npt
           
           if (omega_m(j,1) .eq. lett) then
             omega_m(j,3) = n
             omega_m_copy(j) = 0
           endif
           
        enddo
        
        i = i + 1
        
     enddo
     
   
     
     
     
     end subroutine GET_NEIGHBOUR_ELEM
