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

!> @brief Writes output results for Restart.
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in,out] file_name directory where saving files
!> @param[in] count index for snapshot
!> @param[in] proc mpi process id
!> @param[in] nv number of values to print
!> @param[in] vec values to print

      subroutine WRITE_FILEOUT(file_name,count,proc,nv,vec)
      
      
      character*70 :: file_name

      integer*4 :: count,proc,nv
      real*8, dimension(nv) :: vec
      
      character*70 :: out_file
      integer*4 :: i,lname
      
      lname = len_trim(file_name)
      out_file = file_name(1:lname) // '000000_0000.out'
      
      if (proc .lt. 10) then
         write(out_file(lname+6:lname+6),'(i1)') proc
      else if (proc .lt. 100) then
         write(out_file(lname+5:lname+6),'(i2)') proc
      else if (proc .lt. 1000) then
         write(out_file(lname+4:lname+6),'(i3)') proc     
      else if (proc .lt. 10000) then
         write(out_file(lname+3:lname+6),'(i4)') proc      
      else if (proc .lt. 100000) then
         write(out_file(lname+2:lname+6),'(i5)') proc
      else
         write(out_file(lname+1:lname+6),'(i6)') proc
      endif
      
      if (count .lt. 10) then
         write(out_file(lname+11:lname+11),'(i1)') count
      else if (proc .lt. 100) then
         write(out_file(lname+10:lname+11),'(i2)') count
      else if (proc .lt. 1000) then
         write(out_file(lname+9:lname+11),'(i3)') count    
      else
         write(out_file(lname+8:lname+11),'(i4)') count
      endif
               
      open(20+proc, file=out_file)
      
      do i = 1,nv
         write(20+proc,*) vec(i)
      enddo
      
      close(20+proc)
      
      return
      
      end subroutine WRITE_FILEOUT

