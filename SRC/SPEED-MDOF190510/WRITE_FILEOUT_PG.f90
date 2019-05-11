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

!> @brief Writes output results for Peak Ground Map.
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] monitor_file directory where saving files
!> @param[in] proc mpi process id
!> @param[in] nv number of values to print
!> @param[in] vec values to print
!> @param[out] file_name  file containig results

      subroutine WRITE_FILEOUT_PG(monitor_file,file_name,proc,nv,vec,nvec)
      
      character*70 :: file_name,monitor_file
      character*70 :: out_file, out_file_new

      integer*4 :: proc,nv,count,nvec
      integer*4 :: i,lname

      real*8 :: val
      real*8, dimension(nv,nvec) :: vec
      
      count = 0
      lname = len_trim(file_name)
      out_file = file_name(1:lname) // '_000_00000.MAP'
      
      if (count.lt.10) then
         write(out_file(lname+4:lname+4),'(i1)')count
      else if (count.le.99) then
         write(out_file(lname+3:lname+4),'(i2)')count
      else if (count.le.999) then
         write(out_file(lname+2:lname+4),'(i3)')count
      endif
      
      if (proc.lt.10) then
         write(out_file(lname+10:lname+10),'(i1)')proc
      else if (proc.le.99) then
         write(out_file(lname+9:lname+10),'(i2)')proc
      else if (proc.le.999) then
         write(out_file(lname+8:lname+10),'(i3)')proc
      else if (proc.le.9999) then
         write(out_file(lname+7:lname+10),'(i4)')proc
      else if (proc.le.99999) then
         write(out_file(lname+6:lname+10),'(i5)')proc
      endif
      
      if(len_trim(monitor_file) .ne. 70) then                                                                                  
           out_file_new = monitor_file(1:len_trim(monitor_file)) // '/' // out_file
      else 
           out_file_new = out_file
      endif
      
      
      open(20,file=out_file_new,form='formatted')
      
      do i = 1,nv
                 if (nvec.eq.9) then 
                        write(20,'(1E14.8,1X,1E14.8,1X,1E14.8,1X,1E14.8,1X,1E14.8,1X,1E14.8,1X,1E14.8,1X,1E14.8,1X,1E14.8)') &
                                vec(i,1),vec(i,2),vec(i,3),vec(i,4),vec(i,5), &
                                vec(i,6),vec(i,7),vec(i,8),vec(i,9)
                 endif

                 if (nvec.eq.3) then 
                        write(20,'(1E14.8,1X,1E14.8,1X,1E14.8)') vec(i,1),vec(i,2),vec(i,3)
                 endif
      enddo
      
      close(20)
      
      return
      
      end subroutine WRITE_FILEOUT_PG
