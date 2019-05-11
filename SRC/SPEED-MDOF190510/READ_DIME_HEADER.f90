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

!> @brief Reads dimension in SPEED.input file.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] file_head name for header file (SPEED.input)
!> @param[out] nsnapshot  number of snapshots
!> @param[out] ntime_err  time instant for coputing the norm of the error, for test case

      subroutine READ_DIME_HEADER(file_head,nsnapshots,ntime_err,debug)
         
      implicit none
      
      character*70 :: file_head
      character*80 :: input_line
      character*8 :: keyword
      
      integer*4 :: status
      integer*4 :: nmonitors,nsnapshots,ntime_err,debug
      
      
      nmonitors = 0
      nsnapshots = 0
      ntime_err = 0 
      
      open(20,file=file_head)
      
      do
         read(20,'(A)',IOSTAT = status) input_line
         
         if (status.ne.0) exit
         
         keyword = input_line(1:8)

         if (keyword(1:8).eq.'SNAPSHOT') then
            nsnapshots = nsnapshots + 1
         endif
         if (keyword(1:7).eq.'TIMEERR') then
            ntime_err = ntime_err + 1
         endif
         if (keyword(1:5).eq.'DEBUG') then
            debug = 1
         endif

         
         
      enddo
      
      close(20)
      
      return
      end subroutine READ_DIME_HEADER

