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

!> @brief Writes file FACSXXXXX.mpi and stores info about DG
!! surfaces.
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] mpi_file directory where saving the output
!> @param[in] filename file where writing the output
!> @param[in] faces array contining info about DG surface elements
!> @param[in] nodes constant values for mapping hexes 
!> @param[in] n number of local DG faces (quads)
!> @param[in] ng  number of global DG faces (quads)
!> @param[in] mpi_id  mpi processor id
!> @param[in] mpi_comm  mpi communicator
!> @param[in] cfaces pointer for reading DG faces
!> @param[in] np number of mpi process
!> @param[out] --- FACSXXXXXX.mpi  file containig faces info for the XXXXXX-mpi process

 
     subroutine WRITE_FILE_FACES(mpi_file,filename, faces, nodes, n, ng, mpi_id, mpi_com, cfaces, np)


     implicit none

     include 'SPEED.MPI'

     CHARACTER*70 :: filename, file_mpi, mpi_file, file_mpi_new
     CHARACTER(LEN=1000) :: Format

     INTEGER*4 :: i,j
     INTEGER*4 :: mpierror, unit_mpi, dim2, jstart, unitname       

     INTEGER*4, INTENT(IN) ::  n, mpi_id, mpi_com, np, ng
     INTEGER*4, DIMENSION(np), INTENT(INOUT) :: cfaces

     INTEGER*4, DIMENSION(3,n) :: faces 
     INTEGER*4, DIMENSION(:,:), ALLOCATABLE :: faces_glo 

     REAL*8, DIMENSION(25,n) :: nodes
     REAL*8, DIMENSION(:,:), ALLOCATABLE :: nodes_glo       



     file_mpi = 'FACS000000.mpi'
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
     write(unit_mpi,*) n                
     do i = 1, n
         write(unit_mpi,"(1I2,1X,1I12,1X,1I2,25(2X,ES16.9))") &
                 faces(1,i), faces(2,i), faces(3,i), &
                 nodes(1,i), nodes(2,i), nodes(3,i), & 
                 nodes(4,i), nodes(5,i), nodes(6,i), & 
                 nodes(7,i), nodes(8,i), nodes(9,i), & 
                 nodes(10,i), nodes(11,i), nodes(12,i), & 
                 nodes(13,i), nodes(14,i), nodes(15,i), & 
                 nodes(16,i), nodes(17,i), nodes(18,i), & 
                 nodes(19,i), nodes(20,i), nodes(21,i), & 
                 nodes(22,i), nodes(23,i), nodes(24,i), & 
                 nodes(25,i) 
     enddo
     close(unit_mpi)        

!     write(*,*)  cfaces


     call MPI_BARRIER(mpi_com, mpierror)                  
     



    if(mpi_id .eq. 0) then 
    
       allocate(faces_glo(3,ng), nodes_glo(25,ng))

       do i = 1, np
       
          if (cfaces(i) .ne. 0) then
          
              file_mpi = 'FACS000000.mpi'
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
              read(unit_mpi,*) dim2
              
              if(i.eq. 1) then
                 jstart = 0
              else
                 jstart = sum(cfaces(1:i-1))
              endif
                                    
              do j = 1, dim2

                 read(unit_mpi,*) &
                         faces_glo(1,j+jstart), faces_glo(2,j+jstart), faces_glo(3,j+jstart), &
                         nodes_glo(1,j+jstart), nodes_glo(2,j+jstart), nodes_glo(3,j+jstart), & 
                         nodes_glo(4,j+jstart), nodes_glo(5,j+jstart), nodes_glo(6,j+jstart), & 
                         nodes_glo(7,j+jstart), nodes_glo(8,j+jstart), nodes_glo(9,j+jstart), & 
                         nodes_glo(10,j+jstart), nodes_glo(11,j+jstart), nodes_glo(12,j+jstart), & 
                         nodes_glo(13,j+jstart), nodes_glo(14,j+jstart), nodes_glo(15,j+jstart), & 
                         nodes_glo(16,j+jstart), nodes_glo(17,j+jstart), nodes_glo(18,j+jstart), & 
                         nodes_glo(19,j+jstart), nodes_glo(20,j+jstart), nodes_glo(21,j+jstart), & 
                         nodes_glo(22,j+jstart), nodes_glo(23,j+jstart), nodes_glo(24,j+jstart), & 
                         nodes_glo(25,j+jstart) 



              enddo
          
              close(unit_mpi)

          endif
      
       enddo


     unitname = 400
     open(unitname,file=filename)

      do j = 1, ng

        write(unitname,"(1I2,1X,1I12,1X,1I2,25(2X,ES16.9))") &
                         faces_glo(1,j), faces_glo(2,j), faces_glo(3,j), &
                         nodes_glo(1,j), nodes_glo(2,j), nodes_glo(3,j), & 
                         nodes_glo(4,j), nodes_glo(5,j), nodes_glo(6,j), & 
                         nodes_glo(7,j), nodes_glo(8,j), nodes_glo(9,j), & 
                         nodes_glo(10,j), nodes_glo(11,j), nodes_glo(12,j), & 
                         nodes_glo(13,j), nodes_glo(14,j), nodes_glo(15,j), & 
                         nodes_glo(16,j), nodes_glo(17,j), nodes_glo(18,j), & 
                         nodes_glo(19,j), nodes_glo(20,j), nodes_glo(21,j), & 
                         nodes_glo(22,j), nodes_glo(23,j), nodes_glo(24,j), & 
                         nodes_glo(25,j) 



      enddo

      close(unitname)
      deallocate(faces_glo, nodes_glo)


     endif

 
     call MPI_BARRIER(mpi_com, mpierror)   
 
        
     return

     end subroutine WRITE_FILE_FACES
