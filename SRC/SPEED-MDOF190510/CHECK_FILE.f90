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

!> @brief Checks if a file is present or not
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] filename name of the file
!> @param[out] err_out  0 if filename exists, 101 if filename does not exist

     subroutine CHECK_FILE(filename,err_out)
      
      character*70 :: filename

      logical*4 :: f_ex

      integer*4 :: err_out     
      integer*4 :: i,file_a,file_b
      
      do i=1,70
        if (filename(i:i).ne.' ') exit
      enddo
      
      file_a=i
      do i=70,1,-1
        if (filename(i:i).ne.' ') exit
      enddo
      file_b=i
      
      f_ex=.false. 
      err_out=0
      inquire (file=filename,exist=f_ex)
      if (f_ex) then
        continue
      else
        err_out = 101

      endif
      
      return
      
      end subroutine CHECK_FILE

