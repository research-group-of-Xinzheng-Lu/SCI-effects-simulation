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

!> @brief Matrix-vector multiplication for sparse matrices (RCS format).
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] nnz number of nonzero elements
!> @param[in] As sparse matrix
!> @param[in] Jsp pointer for sparsity pattern
!> @param[in] Isp pointer for sparsity pattern
!> @param[in] n  row of the matrix
!> @param[in] m  column of the matrix
!> @param[in] vet_in input vector for the multiplication As*vet_in
!> @param[in] error control parameter (dummy)                      
!> @param[out] vet_out result of the multiplication As*vet_in

        subroutine MATMUL_SPARSE(As, nnz, Jsp, Isp, vet_out, n, vet_in, m, error)        


        implicit none
        
        integer*4 :: nnz, i, j, n, m, error   
        integer*4 :: Jsp(nnz), Isp(0:n)  
  
        real*8 :: As(nnz), vet_in(m), vet_out(n)

        vet_out = 0.d0
 
              
        do i = 1, n 
           do j = Isp(i-1) + 1, Isp(i)
              vet_out(i) = vet_out(i) + As(j)*vet_in(Jsp(j))
           enddo
        enddo
        

        return
        end subroutine MATMUL_SPARSE
