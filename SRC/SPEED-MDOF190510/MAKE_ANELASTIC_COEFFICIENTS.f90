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

!> @brief Compute anelastic coefficients for Standard Linear Solid model
!! USE THE LIBRARY QR_SOLVE!
!! @author Ilario Mazzieri
!> @date November, 2014 
!> @version 1.0
!> @param[in] nmat number of materials 
!> @param[in] N_SLS number of subdivision in the standard linear solid (dafult N_SLS=3)
!!            possible values are N_SLS=3,4,5 (see MODULES.f90)
!> @param[in] prop_mat material properties
!> @param[in] QS quality factors for s-waves
!> @param[in] QP quality factors for p-waves
!> @param[in] f_val value for the reference frequency
!> @param[in] mpi_id id of the mpi-process 
!> @param[out] Y_lambda anelastic coefficient (Lame 1st modulus)
!> @param[out] Y_mu anelastic coefficient (Lame 2nd modulus)
!> @param[out] frequency_range frequency range where quality factors are assumed to be constant


      subroutine MAKE_ANELASTIC_COEFFICIENTS(nmat, N_SLS, prop_mat, QS, QP, f_val, &
                                         Y_lambda, Y_mu, frequency_range, mpi_id)
 
      
      use speed_exit_codes                  

      implicit none
 
      integer*4 :: nmat,mpi_id, im, N_SLS, i, j, k, &
                   PIVOT_QP(N_SLS),PIVOT_QS(N_SLS)
                   
      real*8 :: f_val, fmin, fmax, esp1, esp2, deltax, &
                rho, lambda, mu, vp2, vs2, pi
      
      real*8 :: prop_mat(nmat,4), QS(nmat), QP(nmat), &
                Y_lambda(nmat,N_SLS),&
                Y_mu(nmat,N_SLS), &
                frequency_range(N_SLS),&
                frequency_range_sampling(2*N_SLS-1),&
                YS(N_SLS),YP(N_SLS), &
                A_QP(2*N_SLS-1,N_SLS),RHS_QP(2*N_SLS-1),A_QS(2*N_SLS-1,N_SLS),RHS_QS(2*N_SLS-1)

      pi = 4.d0*datan(1.d0);
      
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
            
     esp1 = LOG10(fmin); esp2 = LOG10(fmax);
     deltax = (esp2-esp1)/(N_SLS-1);
     

     if (N_SLS .eq. 3) then

        frequency_range(1) = 2.*pi*fmin;
        frequency_range(2) = 2.*pi*10*fmin;
        frequency_range(3) = 2.*pi*fmax;
        frequency_range_sampling(1) = frequency_range(1);
        frequency_range_sampling(2) = 0.5*(frequency_range(1)+frequency_range(2))
        frequency_range_sampling(3) = frequency_range(2);
        frequency_range_sampling(4) = 0.5*(frequency_range(2)+frequency_range(3))
        frequency_range_sampling(5) = frequency_range(3);
                

     elseif (N_SLS .eq. 4) then

        frequency_range(1) = 2.*pi*fmin;
        frequency_range(2) = 2.*pi*10**(esp1+deltax)
        frequency_range(3) = 2.*pi*10**(esp1+2*deltax)
        frequency_range(4) = 2.*pi*fmax;
        frequency_range_sampling(1) = frequency_range(1);
        frequency_range_sampling(2) = 0.5*(frequency_range(1)+frequency_range(2))
        frequency_range_sampling(3) = frequency_range(2);
        frequency_range_sampling(4) = 0.5*(frequency_range(2)+frequency_range(3))
        frequency_range_sampling(5) = frequency_range(3);
        frequency_range_sampling(6) = 0.5*(frequency_range(3)+frequency_range(4))
        frequency_range_sampling(7) = frequency_range(4);


      elseif (N_SLS .eq. 5) then

        frequency_range(1) = 2.*pi*fmin;
        frequency_range(2) = 2.*pi*10**(esp1+deltax)
        frequency_range(3) = 2.*pi*10**(esp1+2*deltax)
        frequency_range(4) = 2.*pi*10**(esp1+3*deltax)
        frequency_range(5) = 2.*pi*fmax;
        frequency_range_sampling(1) = frequency_range(1);
        frequency_range_sampling(2) = 0.5*(frequency_range(1)+frequency_range(2))
        frequency_range_sampling(3) = frequency_range(2);
        frequency_range_sampling(4) = 0.5*(frequency_range(2)+frequency_range(3))
        frequency_range_sampling(5) = frequency_range(3);
        frequency_range_sampling(6) = 0.5*(frequency_range(3)+frequency_range(4))
        frequency_range_sampling(7) = frequency_range(4);
        frequency_range_sampling(6) = 0.5*(frequency_range(4)+frequency_range(5))
        frequency_range_sampling(7) = frequency_range(5);


      endif           

      if(mpi_id .eq. 0) then
         write(*,*) 'Sampled frequencies:'
         write(*,*) 'FREQ =' ,frequency_range_sampling
      endif
 
 
 
      do im = 1, nmat
         
         rho = prop_mat(im,1); lambda=prop_mat(im,2); mu=prop_mat(im,3); 
         vp2 = (lambda + 2.d0*mu)/rho;  vs2 = mu/rho;
  
         Y_mu(im,:) = 0.d0; Y_lambda(im,:) = 0.d0;
         if (QP(im) .ne. 0.d0 .and. QS(im) .ne. 0.d0) then
 
           do i = 1, 2*N_SLS-1
              do j = 1, N_SLS

                A_QP(i,j) = (frequency_range(j)*frequency_range_sampling(i) + frequency_range(j)**2*(1.d0/QP(im))) &
                       / (frequency_range(j)**2 + frequency_range_sampling(i)**2);

                A_QS(i,j) = (frequency_range(j)*frequency_range_sampling(i) + frequency_range(j)**2*(1.d0/QS(im))) &
                      / (frequency_range(j)**2 + frequency_range_sampling(i)**2);
                      
              enddo
            
              RHS_QP(i) = 1.d0/QP(im);
              RHS_QS(i) = 1.d0/QS(im);
            
           enddo
         
           !call FACTORIZE_MATRIX(A_QP,N_SLS,PIVOT_QP)
           !call FACTORIZE_MATRIX(A_QS,N_SLS,PIVOT_QS)
         
           !call DIRECT_LU_SOLVER(A_QP,RHS_QP,N_SLS,YP,PIVOT_QP)
           !call DIRECT_LU_SOLVER(A_QS,RHS_QS,N_SLS,YS,PIVOT_QS)
            
           call QR_SOLVE(2*N_SLS-1,N_SLS, A_QP, RHS_QP, YP) 
           call QR_SOLVE(2*N_SLS-1,N_SLS, A_QS, RHS_QS, YS) 
                  
           Y_mu(im,:) = YS;
           Y_lambda(im,:) = (vp2*YP - 2.d0*vs2*YS)/(vp2-2.d0*vs2);
      
         endif
     
!         do i = 1, N_SLS
!            if (mpi_id .eq. 0 .and.Y_lambda(im,i) .lt. 0.d0 .or. Y_mu(im,i) .lt. 0) then
!              write(*,*) 'ERROR! ANELASTIC COEFFICIENTS HAVE TO BE POSITIVE!'
!              write(*,*) 'ANELASTIC COEFFICIENTS FOR MATERIAL ', im, '  ARE '
!              write(*,*) 'Y_MU = ', Y_mu(im,:)
!              write(*,*) 'Y_LAMBDA = ', Y_lambda(im,:)
!              write(*,*) 'Y_P = ', YP
!              write(*,*) '----------------------------------------------------'              
!            endif
!         enddo
         
         if(mpi_id .eq. 0) then
           write(*,*) 'ANELASTIC COEFFICIENTS FOR MATERIAL ', im, '  ARE '
           write(*,*) 'Y_MU = ', Y_mu(im,:)
           write(*,*) 'Y_LAMBDA = ', Y_lambda(im,:)
           write(*,*) 'Y_P = ', YP
           write(*,*) '----------------------------------------------------'
         endif
     
       enddo 

       return
      
       end subroutine MAKE_ANELASTIC_COEFFICIENTS
 
 
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
 
      subroutine FACTORIZE_MATRIX(A,ndim,pivot)

      use speed_exit_codes

      implicit none

      real*8, dimension(ndim,ndim), intent(inout) :: A
      integer*4, intent(in) :: ndim 
      integer*4, dimension(ndim), intent(out) :: pivot

 
      real*8, dimension(ndim) :: dd
      integer*4 :: ii,ij,ik,iflag,ip,is
      real*8 :: aw,cmax,rmax,temp,sun

      pivot(:) = 0
      iflag = 1
      do ii=1,ndim
	  pivot(ii) = ii
	  rmax = 0.d0
          do ij=1,ndim
	    rmax = max(rmax, ABS(A(ii, ij)))
          enddo
      
	  if (rmax.eq.0) then
		   iflag = 0
		   rmax = 1.d0
	  endif
	  dd(ii) = rmax
      enddo

      if (ndim.le.1) return

      do ik=1,ndim-1

	  cmax = ABS(A(ik, ik)) / dd(ik)
          is=ik
          do ii=ik+1,ndim
	    aw = ABS(A(ii,ik)) / dd(ii)
	    if (aw.gt.cmax) then
		cmax = aw
		is=ii
	    endif
          enddo

	  if (cmax.eq.0) then
		iflag = 0
	  else
            if (is.gt.ik) then
		   iflag = -iflag
		   ii=pivot(is)
		   pivot(is)=pivot(ik)
		   pivot(ik) = ii
		   temp=dd(is)
		   dd(is)=dd(ik)
		   dd(ik) = temp
		    do ij=1,ndim
		       temp =A(is,ij)
		       A(is,ij)=A(ik,ij)
		       A(ik, ij) = temp
        enddo

	     endif

	       do ii=ik+1,ndim
		 A(ii, ik) = A(ii, ik) / A(ik, ik)
	         do ij=ik+1,ndim
		   A(ii, ij) = A(ii, ij) - A(ii, ik) * A(ik, ij)
                 enddo
               enddo
 
	   endif
      enddo
      if (A(ndim,ndim) .eq. 0.d0) then 
          iflag = 0
      endif

      if (iflag .eq. 0 )  then
          write(*,*) 'SINGULAR MATRIX'
          call EXIT(EXIT_SINGULARMTX)
      endif
 

      end subroutine FACTORIZE_MATRIX

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


      subroutine DIRECT_LU_SOLVER(A,tn,ndim,x_sol,pivot)


      integer*4, intent(in) :: ndim 
      real*8, dimension(ndim), intent(inout) :: x_sol
      real*8, dimension(ndim,ndim), intent(inout) :: A
      real*8, dimension(ndim), intent(inout) :: tn
      integer*4, dimension(ndim), intent(in) :: pivot

      integer*4 :: ii,ij,ik,iflag,ip
      real*8, dimension(ndim) ::y_sol
      real*8, dimension(ndim) :: dd
      real*8 :: aw,cmax,rmax,temp,sun

      if (ndim .le. 1) then
	x_sol(1) = tn(1) / A(1, 1)
        return
      endif
  
      ip = pivot(1)
      y_sol(1) = tn(ip)
      do ii=2,ndim
        sun = 0.d0
        do ij=1,ii-1
	  sun = A(ii, ij) * y_sol(ij) + sun
        enddo
	ip = pivot(ii)
	y_sol(ii) = tn(ip) - sun
      enddo
  
      x_sol(ndim) = y_sol(ndim) / A(ndim, ndim)
  
      do ii=ndim-1,1,-1
        sun = 0.d0
        do ij=ii+1,ndim
          sun = A(ii, ij) * x_sol(ij) + sun
        enddo
          x_sol(ii) = (y_sol(ii) - sun) / A(ii, ii)
      enddo

      return
      end subroutine DIRECT_LU_SOLVER
 
 
 
 
     
 
 

