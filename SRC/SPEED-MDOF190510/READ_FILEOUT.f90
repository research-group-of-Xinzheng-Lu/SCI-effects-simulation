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

!> @brief Reads files *.out for the restart of the simulation.
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] filename  file name
!> @param[in] counter  index for snapshot
!> @param[in] procs  MPI process id 
!> @param[in] nb_vec  number of values
!> @param[out]  vec  values read

      subroutine READ_FILEOUT(filename,counter,procs,nb_vec,vec)
      
      implicit none
      
      integer*4 :: counter,procs,nb_vec
      real*8, dimension(nb_vec) :: vec
      character*70 :: filename
      
      character*70 :: out_file
      integer*4 :: i,lname
      
      lname = len_trim(filename)
      out_file = filename(1:lname) // '000000_0000.out'
      
      if (procs .lt. 10) then
         write(out_file(lname+6:lname+6),'(i1)') procs
      else if (procs .lt. 100) then
         write(out_file(lname+5:lname+6),'(i2)') procs
      else if (procs .lt. 1000) then
         write(out_file(lname+4:lname+6),'(i3)') procs     
      else if (procs .lt. 10000) then
         write(out_file(lname+3:lname+6),'(i4)') procs      
      else if (procs .lt. 100000) then
         write(out_file(lname+2:lname+6),'(i5)') procs
      else
         write(out_file(lname+1:lname+6),'(i6)') procs
      endif
      
      if (counter .lt. 10) then
         write(out_file(lname+11:lname+11),'(i1)') counter
      else if (procs .lt. 100) then
         write(out_file(lname+10:lname+11),'(i2)') counter
      else if (procs .lt. 1000) then
         write(out_file(lname+9:lname+11),'(i3)') counter    
      else
         write(out_file(lname+8:lname+11),'(i4)') counter
      endif
      
      open(20+procs, file=out_file)
      
      do i = 1,nb_vec
         read(20+procs,*) vec(i)
      enddo
      
      close(20+procs)
      
      return

      
      return
      
      end subroutine READ_FILEOUT
      

