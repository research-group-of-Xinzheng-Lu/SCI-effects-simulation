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

!> @brief Generates a set of files where output results will be written.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nmonitlist number of monitored points
!> @param[in] elem_mlst list of monitored elements
!> @param[in] local_el_num local element numbering
!> @param[in] ne_loc number of local elements
!> @param[in] count_monitor number of monitors for the mpi process 
!> @param[in] monitor_file folder where MONITOR files are stored.
!> @param[in] mpi_id id for MPI process
!> @param[in] option_out_var options for output
!> @param[in,out] filename file name 
!> @param[out] ... creation of files MONITORXXXXXX.D/.V/.A/.E/.S/.O/.INFO

      subroutine  MAKE_MONITOR_FILES(nmonitlst,elem_mlst,local_el_num, ne_loc,&
                                     count_monitor, filename, monitor_file, mpi_id, &
                                     option_out_var)

      implicit none
      
      integer*4, intent(inout) :: count_monitor
      integer*4 :: imon, ielem, ne_loc, ie, mpi_id, unit_monitor, i
      integer*4 :: nmonitlst
 
      integer*4, dimension(nmonitlst) :: elem_mlst
      integer*4, dimension(ne_loc) :: local_el_num
      integer*4, dimension (6) :: option_out_var 
      integer*4, dimension(:), allocatable :: monitor_index
      
      character*70 :: filename, monitor_file, filename_new
      character*5 :: filesuffix
     
!      write(*,*) 'Making MONITOR files'

      count_monitor = 0      
      do imon = 1, nmonitlst
         ielem = elem_mlst(imon)     !ie = index hexahedra containing monitor
         call GET_INDLOC_FROM_INDGLO(local_el_num, ne_loc, ielem, ie)                                             
         if (ie .ne. 0) count_monitor = count_monitor + 1
      enddo       
                        
      allocate(monitor_index(count_monitor));  count_monitor = 0      
            
      do imon = 1, nmonitlst
            ielem = elem_mlst(imon)     !ie = index hexahedra containing monitor
         call GET_INDLOC_FROM_INDGLO(local_el_num, ne_loc, ielem, ie)                                             
         if (ie .ne. 0) then 
           count_monitor = count_monitor + 1
           monitor_index(count_monitor) = imon !node_mlst(imon)
        endif   
      enddo       

      unit_monitor = 40 + mpi_id        
           
      if (count_monitor .ne. 0) then
                                               
      ! Loop over variables                                         
      do i = 1,6

          filesuffix = ''
          if (i .eq. 1 .and. option_out_var(i) .eq. 1 ) filesuffix = '.D'
          if (i .eq. 2 .and. option_out_var(i) .eq. 1 ) filesuffix = '.V'
          if (i .eq. 3 .and. option_out_var(i) .eq. 1 ) filesuffix = '.A'
          if (i .eq. 4 .and. option_out_var(i) .eq. 1 ) filesuffix = '.S'
          if (i .eq. 5 .and. option_out_var(i) .eq. 1 ) filesuffix = '.E'

          ! Skip writing if necessary
          if (option_out_var(i) .eq. 1) then
            ! Create output filenames for data output
           
!          	write(*,'(A,I0,A)') 'mpi_id ', mpi_id, ' suffix ' // filesuffix
            write(filename, '(A,I5.5,A5)') 'MONITOR', mpi_id, filesuffix
!            write(*,*) 'Variable = ', i
!            write(*,*) 'Output filename = ', filename
!            write(*,*) 'Monitor directory: ', monitor_file
!            write(*,*) 'length ', len_trim(monitor_file)

            ! Prepend directory name
!            if(len_trim(monitor_file) + len_trim(filename) + 1 .le. 255) then                                                                                  
      	    if(len_trim(monitor_file) .ne. 70) then  
               filename_new = monitor_file(1:len_trim(monitor_file)) // '/' // filename
            else
               filename_new = filename
            endif

            write(*,*) 'Creating file ', trim(filename_new)
            open(unit_monitor,file=trim(filename_new))        !open MONITORXXXXX.D, etc...                
            close(unit_monitor)
          endif
                        
        enddo        

        write(filename, '(A,I5.5,A5)') 'MONITOR', mpi_id, '.INFO'
!        write(*,*) 'INFO file: ', filename
                        
        if(len_trim(monitor_file) .ne. 70) then                                                                                  
!        if(len_trim(monitor_file) + len_trim(filename) + 1 .le. 255) then                                                                                  
          filename_new = monitor_file(1:len_trim(monitor_file)) // '/' // filename
        else
          filename_new = filename  
        endif

        write(*,*) 'Creating file ', trim(filename_new)
        open(unit_monitor,file=trim(filename_new))                        
        write(unit_monitor,*) count_monitor, monitor_index
        close(unit_monitor)

     endif

     write(*,*)     
     deallocate(monitor_index)
    
     end subroutine MAKE_MONITOR_FILES  
            
            
            
            
