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

!> @brief Makes connectivity vector for DG elements.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nnz_dg length of cs_dg
!> @param[in] tot_el  total dg elements
!> @param[in] np number of mpi processes 
!> @param[in] mpi_file folder where mpi file are stored
!> @param[in,out] cs_dg dg connectivity vector  


     subroutine MAKE_SPX_CON_DG(cs_dg, nnz_dg, tot_el, np, mpi_file)
     
     implicit none
     
     character*70 :: file_mpi, mpi_file, file_mpi_new

     integer*4 :: nnz_dg, np
     integer*4 :: i, j, unit_mpi, nnz_mpi, tot_el, mpi_np, ishift     

     integer*4, dimension(:), allocatable :: cs_mpi
     integer*4, dimension(0:nnz_dg), intent(inout) :: cs_dg
    
     cs_dg(0) = tot_el
     ishift = 0
     
     do i = 1, np
             file_mpi = 'DGCS000000.mpi'
             unit_mpi = 40                                 

        if (i-1 .lt. 10) then                                        
            write(file_mpi(10:10),'(i1)') i-1                
        else if (i-1 .lt. 100) then                                
            write(file_mpi(9:10),'(i2)') i-1                
        else if (i-1 .lt. 1000) then                                
            write(file_mpi(8:10),'(i3)') i-1                
        else if (i-1 .lt. 10000) then                                
            write(file_mpi(7:10),'(i4)') i-1                
        else if (i-1 .lt. 100000) then        
            write(file_mpi(6:10),'(i5)') i-1                
        else if (i-1 .lt. 1000000) then                                
            write(file_mpi(5:10),'(i6)') i-1                
         endif
                                                              
        if(len_trim(mpi_file) .ne. 70) then                                                                                  
           file_mpi_new = mpi_file(1:len_trim(mpi_file)) // '/' // file_mpi
        else
           file_mpi_new = file_mpi
        endif

             open(unit_mpi,file=file_mpi_new)
                                     
             read(unit_mpi,*) nnz_mpi
             allocate(cs_mpi(0:nnz_mpi))          
        do j = 0, nnz_mpi
           read(unit_mpi,*) cs_mpi(j)                        
        enddo
        close(unit_mpi)        

        do j = 1, cs_mpi(0) - 1
        
           cs_dg(ishift + j) = cs_dg(ishift + j-1) + (cs_mpi(j) - cs_mpi(j-1))
           cs_dg(cs_dg(ishift + j-1) : cs_dg(ishift + j) - 1 ) = cs_mpi(cs_mpi(j-1) : cs_mpi(j) -1) 
           
        enddo
        ishift = ishift + cs_mpi(0)-1
        deallocate(cs_mpi)        
             

      enddo      
     

     
     end subroutine MAKE_SPX_CON_DG
