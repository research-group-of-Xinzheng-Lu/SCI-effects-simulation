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

!> @brief Reads files such as MLST.input, MLST.position or MPGM.input. 
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] filec file name to read
!> @param[in] nmonitors_pgm number of monitors 
!> @param[out] n_monitor_pgm id monitor number (id spectral node)
!> @param[out] el_monitor_pgm id monitor element (id grid element)
!> @param[out] xr_monitor_pgm x-coordinate for monitor
!> @param[out] yr_monitor_pgm y-coordinate for monitor
!> @param[out] zr_monitor_pgm z-coordinate for monitor

      subroutine READ_FILE_MPGM(filec,&
                        nmonitors_pgm,&
                        n_monitor_pgm,&
                        el_monitor_pgm,&
                        xr_monitor_pgm,&
                        yr_monitor_pgm,&
                        zr_monitor_pgm)

      implicit none

      character*70 :: filec
      character*100000 :: input_line

      integer*4 :: nmonitors_pgm
      integer*4 :: i
      integer*4 :: trash
      integer*4 :: ileft,iright
      integer*4 :: status

      integer*4,dimension(nmonitors_pgm) :: n_monitor_pgm
      integer*4,dimension(nmonitors_pgm) :: el_monitor_pgm

      real*8,dimension(nmonitors_pgm) :: xr_monitor_pgm
      real*8,dimension(nmonitors_pgm) :: yr_monitor_pgm
      real*8,dimension(nmonitors_pgm) :: zr_monitor_pgm


              open(20,file=filec)
                read(20,'(A)',IOSTAT = status) input_line
                ileft = 1
                iright = len(input_line)
                read(input_line(ileft:iright),*)nmonitors_pgm

                do i = 1,nmonitors_pgm
                        read(20,'(A)',IOSTAT = status) input_line
                        if (status.ne.0) exit

                        ileft = 1
                        iright = len(input_line)
                        read(input_line(ileft:iright),*)trash,n_monitor_pgm(i),el_monitor_pgm(i),&
                                                xr_monitor_pgm(i),yr_monitor_pgm(i),zr_monitor_pgm(i)

                enddo
                    
        close(20)

      return

      end subroutine READ_FILE_MPGM

