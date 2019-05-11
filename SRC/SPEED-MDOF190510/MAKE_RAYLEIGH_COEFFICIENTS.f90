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

!> @brief Compute the A0, A1 coefficient for the Rayleigh damping matrix C = A0 M + A1 K,
!! M = mass matrix, K= stiffness matrix
!! @author Ilario Mazzieri
!> @date November, 2014 
!> @version 1.0
!> @param[in] namt number of material
!> @param[in] prop_mat material properties
!> @param[in] QS quality factors for the s-waves
!> @param[in] f_val value for the reference frequency
!> @param[in] mpi_id id of the mpi-process
!> @param[out] A0_ray Rayleigh coefficients
!> @param[out] A1_ray Rayleigh coefficients



      subroutine MAKE_RAYLEIGH_COEFFICIENTS(nmat, A0_ray, A1_ray, prop_mat, QS, f_val, mpi_id)
 
      
      implicit none
      						
 
      integer*4 :: nmat, mpi_id, im, i, j, k
                   
                   
      real*8 :: f_val, fmin, fmax, esp1, esp2, deltax, &
                rho, lambda, mu, vp2, vs2, pi
      
      real*8 :: prop_mat(nmat,4), QS(nmat), csi(nmat),&
                A0_ray(nmat), A1_ray(nmat)
                
      pi = 4.d0*datan(1.d0);
 
      csi =1./(2*QS);
      
      
      if (f_val .le. 10) then
         fmin = 0.05d0; !fmax = 5.d0;
      elseif(f_val .le. 100) then
         fmin = 5.d0; !fmax = 100.d0;   
      elseif(f_val .le. 1000) then
         fmin = 50.d0; !fmax = 1000.d0;   
      elseif(f_val .le. 10000) then
         fmin = 500.d0; !fmax = 10000.d0;   
      elseif(f_val .le. 100000) then
         fmin = 5000.d0; !fmax = 100000.d0;
      endif
      
      fmax=100*fmin;
      
      if(mpi_id .eq. 0) then
         write(*,*) 'Frequency range where Quality Factor is assumed to be constant is:'
         write(*,*) 'FMIN =' ,fmin, '  FMAX =', fmax;
      endif
             
 
 
      do im = 1, nmat
         
         
         A0_ray(im) = csi(im)* 4.*pi*fmin*fmax/(fmin+fmax);
         A1_ray(im) = csi(im)/(pi*(fmin+fmax));
         
         
         if(mpi_id .eq. 0) then
           write(*,*) 'RAYLEIGH COEFFICIENTS FOR MATERIAL ', im, '  ARE '
           write(*,*) 'A0 = ', A0_ray(im)
           write(*,*) 'A1 = ', A1_ray(im)
           write(*,*) '----------------------------------------------------'
         endif
     
       enddo 

       return
      
       end subroutine MAKE_RAYLEIGH_COEFFICIENTS
