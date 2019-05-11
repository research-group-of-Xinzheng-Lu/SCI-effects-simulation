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

!> @brief Writes DGCSXXXXXX.mpi files for DG connectivity. 
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] el_new structure for DG faces
!> @param[in] nel_dg_loc  number of local dg elements
!> @param[in] nel_dg_glo  number of global dg elements
!> @param[in] mpi_id  mpi processor identity
!> @param[in] mpi_comm  mpi communicator
!> @param[in] mpi_np  number of mpi processes
!> @param[in] cs_nnz_loc length of cs_loc
!> @param[in] cs_loc local spectral connectivity vector
!> @param[in] nn_loc number of local nodes
!> @param[in] local_n_num local node numeration
!> @param[in] local_el_num local element numeration
!> @param[in] mpi_file folder where to write *.mpi file
!> @param[out] DGCSXXXXXX.mpi DG connectivity vector for MPI routines (XXXXXX -> processor number)
!> @param[out] total_cs_nnz_mpi  total number of DG nodes for mpi routines
!> @param[out] total_el  total number of DG elements for mpi routines

    subroutine SETUP_DG4MPI(el_new, nel_dg_loc, nel_dg_glo, &
                          mpi_id, mpi_comm, mpi_np, cs_nnz_loc, cs_loc, &
                          ne_loc, local_el_num, &
                          local_n_num, nn_loc, total_cs_nnz_mpi, total_el, mpi_file)
                                                    
     use max_var
     use DGJUMP

     implicit none

                    
     include 'SPEED.MPI'

     type(el4loop), dimension(nel_dg_loc), intent(in) :: el_new 
  
     character*70 :: file_mpi, mpi_file, file_mpi_new
     character*70  :: filename

     integer*4 :: nel_dg_loc, cs_nnz_loc, ne_loc, nel_dg_glo, mpi_id, mpi_comm, mpi_np, nn_loc
     integer*4 :: mpierror, cs_nnz_mpi, nofel_mpi, el_dg
     integer*4 :: i, j, k, ie, max_dime, dime_mpi, unit_mpi  
     integer*4, intent(out) :: total_cs_nnz_mpi, total_el

     integer*4, dimension(:), allocatable :: cs_mpi
     integer*4, dimension(0:cs_nnz_loc) :: cs_loc
     integer*4, dimension(ne_loc) :: local_el_num
     integer*4, dimension(nn_loc) :: local_n_num

     cs_nnz_mpi = 0 
     nofel_mpi = 1
     
     do i = 1, nel_dg_loc
         if (i .eq. 1) then
            cs_nnz_mpi = cs_nnz_mpi + el_new(i)%deg**3 + 1 + 1
            nofel_mpi = nofel_mpi + 1

         elseif(el_new(i)%ind .ne. el_new(i-1)%ind) then   
            cs_nnz_mpi = cs_nnz_mpi + el_new(i)%deg**3 + 1 + 1                
            nofel_mpi = nofel_mpi + 1            
         endif               
      enddo     
         
     allocate(cs_mpi(0:cs_nnz_mpi))
     cs_mpi =  0
     cs_mpi(0) = nofel_mpi
     k = 0
    
    
     do i = 1, nel_dg_loc

         if (i .eq. 1) then
         
            
            call GET_INDLOC_FROM_INDGLO(local_el_num, ne_loc, el_new(i)%ind, ie)                                    
            k = 1
            cs_mpi(k) = cs_mpi(k-1) + el_new(i)%deg**3 + 1
            cs_mpi(cs_mpi(k-1)) = el_new(i)%ind                 
             
            do j= 1, el_new(i)%deg**3
                cs_mpi(cs_mpi(k-1) + j) = local_n_num(cs_loc(cs_loc(ie-1) + j))
            enddo   
                 
         elseif(el_new(i)%ind .ne. el_new(i-1)%ind) then   

            call GET_INDLOC_FROM_INDGLO(local_el_num, ne_loc, el_new(i)%ind, ie)          
            k = k + 1

            cs_mpi(k) = cs_mpi(k-1) + el_new(i)%deg**3 + 1
            cs_mpi(cs_mpi(k-1)) = el_new(i)%ind                 
             
             do j= 1, el_new(i)%deg**3
                cs_mpi(cs_mpi(k-1) + j) = local_n_num(cs_loc(cs_loc(ie-1) + j))
             enddo   
                                 
         endif                       
      enddo


!      write(*,*)  'dopo', cs_mpi


     file_mpi = 'DGCS000000.mpi'
     unit_mpi = 40                 
     
     if (mpi_id .lt. 10) then                                        
         write(file_mpi(10:10),'(i1)') mpi_id                
     else if (mpi_id .lt. 100) then                                
         write(file_mpi(9:10),'(i2)') mpi_id
     else if (mpi_id .lt. 1000) then                                
         write(file_mpi(8:10),'(i3)') mpi_id
     else if (mpi_id .lt. 10000) then                                
         write(file_mpi(7:10),'(i4)') mpi_id
     else if (mpi_id .lt. 100000) then                                
         write(file_mpi(6:10),'(i5)') mpi_id
     else if (mpi_id .lt. 1000000) then                                
         write(file_mpi(5:10),'(i6)') mpi_id
     endif                                                              
     
     if(len_trim(mpi_file) .ne. 70) then                                                                                  
         file_mpi_new = mpi_file(1:len_trim(mpi_file)) // '/' // file_mpi
     else 
         file_mpi_new = file_mpi     
     endif
          
     open(unit_mpi,file=file_mpi_new)                        
     write(unit_mpi,*) cs_nnz_mpi                
     do i = 0, cs_nnz_mpi
         write(unit_mpi,*) cs_mpi(i)               
     enddo
     close(unit_mpi)        


     deallocate(cs_mpi)

     call MPI_BARRIER(mpi_comm, mpierror)    
     
     total_cs_nnz_mpi = 0
     total_el = 0
     
     do i = 1, mpi_np
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
             read(unit_mpi,*) cs_nnz_mpi
             read(unit_mpi,*) el_dg
             
             total_cs_nnz_mpi = total_cs_nnz_mpi + cs_nnz_mpi                 
             total_el = total_el + el_dg -1
             close(unit_mpi)        
     
      enddo 
      
     total_el = total_el + 1

     return                     
                          
                          
    end subroutine SETUP_DG4MPI                          
