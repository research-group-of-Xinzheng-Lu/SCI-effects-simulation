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

!> @brief Reads dimensions in gridfile (*.mesh)
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] gridfile  file name (*.mesh)
!> @param[in] nb_mate  number of blocks (materials)
!> @param[in] lab_mate  tag for materials
!> @param[in] nb_diriX  number of Dirichlet b.c. (x-dir)
!> @param[in] nb_diriY  number of Dirichlet b.c. (y-dir)
!> @param[in] nb_diriZ  number of Dirichlet b.c. (z-dir)
!> @param[in] lab_diriX  label for Dirichlet b.c. (x-dir)
!> @param[in] lab_diriY  label for Dirichlet b.c. (y-dir)
!> @param[in] lab_diriZ  label for Dirichlet b.c. (z-dir)
!> @param[in] nb_neuX  number of Neumann boundary loads (x-dir)
!> @param[in] nb_neuY  number of Neumann boundary loads (y-dir)
!> @param[in] nb_neuZ  number of Neumann boundary loads (z-dir)
!> @param[in] lab_neuX   label for Neumann boundary loads (x-dir)
!> @param[in] lab_neuY   label for Neumann boundary loads (y-dir)
!> @param[in] lab_neuZ   label for Neumann boundary loads (z-dir)
!> @param[in] nb_neuN  number of Neumann boundary loads (normal direction)
!> @param[in] lab_neuN   label for Neumann boundary loads
!> @param[in] nnb_abc  number of absorbing boundary conditions
!> @param[in] lab_abc   label for absorbing boundary condition
!> @param[in] nnb_dg  number of discontinuos interfaces (where DG scheme is applied)
!> @param[in] lab_dg   label for DG interface conditions
!> @param[out] nb_node number of mesh node
!> @param[out] nb_hexa number of hexahedral elements
!> @param[out] nb_quad number of quadrilateral elements



      subroutine READ_DIME_FILEMESH(gridfile,nb_mate,lab_mate,&
                           nb_diriX,lab_diriX,nb_diriY,lab_diriY,nb_diriZ,lab_diriZ,&
                           nb_neuX,lab_neuX,nb_neuY,lab_neuY,nb_neuZ,lab_neuZ,&
                           nb_neuN,lab_neuN, & 
                           nb_abc,lab_abc,&
                           nb_dg,lab_dg,&
                           nb_node,nb_hexa,nb_quad)
      

      implicit none
      
      character*70 :: gridfile
      character*80 :: inline
      character*20 :: el_code

      integer*4 :: nb_mate,nb_diriX,nb_diriY,nb_diriZ,nb_neuX,nb_neuY,nb_neuZ,nb_abc,nb_dg
      integer*4 :: nb_neuN 
      integer*4 :: nb_elem,ie,i,mat_code,ileft,iright,sl,trash,status,control

      integer*4, dimension(nb_mate) :: lab_mate

      integer*4, dimension(nb_diriX) :: lab_diriX
      integer*4, dimension(nb_diriY) :: lab_diriY
      integer*4, dimension(nb_diriZ) :: lab_diriZ
      integer*4, dimension(nb_neuX) :: lab_neuX
      integer*4, dimension(nb_neuY) :: lab_neuY
      integer*4, dimension(nb_neuZ) :: lab_neuZ
      integer*4, dimension(nb_neuN) :: lab_neuN
      integer*4, dimension(nb_abc) :: lab_abc
      integer*4, dimension(nb_dg) :: lab_dg
      integer*4 :: nb_node,nb_hexa,nb_quad
      
      
      nb_node = 0;  nb_hexa = 0;  nb_quad = 0
      
      open(40,file=gridfile)
      
      do 
         read(40,'(A)') inline
         if (inline(1:1) .ne. '#') exit
      enddo
      
      read(inline,*)nb_node,nb_elem
      
      do i = 1,nb_node
         read(40,'(A)')inline
      enddo
      
      do ie = 1,nb_elem
         read(40,'(A)')inline
         
         sl = len(inline)
         ileft = 0
         iright = 0 
         do i = 1,sl
            if (inline(i:i).ge.'A') exit
         enddo
         ileft = i
         do i = ileft,sl
            if (inline(i:i).lt.'A') exit
         enddo
         iright = i
         
         el_code = inline(ileft:iright)
         
         read(inline(1:ileft),*)trash,mat_code
         
         if ((el_code.eq.'hex').or.(el_code.eq.'HEX')) then
            control = 0
            do i = 1,nb_mate
               if (lab_mate(i).eq.mat_code) control = 1
            enddo
            !
            if (control.ne.0) nb_hexa = nb_hexa +1
         !
         elseif ((el_code.eq.'quad').or.(el_code.eq.'QUAD')) then
            control = 0
            do i = 1,nb_diriX
               if (lab_diriX(i).eq.mat_code) control = 1
            enddo
            do i = 1,nb_diriY
               if (lab_diriY(i).eq.mat_code) control = 1
            enddo
            do i = 1,nb_diriZ
               if (lab_diriZ(i).eq.mat_code) control = 1
            enddo
            do i = 1,nb_neuX
               if (lab_neuX(i).eq.mat_code) control = 1
            enddo
            do i = 1,nb_neuY
               if (lab_neuY(i).eq.mat_code) control = 1
            enddo
            do i = 1,nb_neuZ
               if (lab_neuZ(i).eq.mat_code) control = 1
            enddo
            do i = 1,nb_neuN                                                        
               if (lab_neuN(i).eq.mat_code) control = 1                                
            enddo                                                                
            do i = 1,nb_abc
               if (lab_abc(i).eq.mat_code) control = 1
            enddo
            do i = 1,nb_dg
               if (lab_dg(i).eq.mat_code) control = 1
            enddo
            

            if (control.ne.0) nb_quad = nb_quad +1
         endif
      enddo
      
      close(40)
      
      return
      end subroutine READ_DIME_FILEMESH

