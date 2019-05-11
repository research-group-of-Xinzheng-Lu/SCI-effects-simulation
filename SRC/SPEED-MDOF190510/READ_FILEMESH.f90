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
!> @param[in] nb_mat  number of blocks (materials)
!> @param[in] lab_mat(nb_mat)  label for the materials
!> @param[in] nb_diriX  number of Dirichlet boundary conditions (x-direction)
!> @param[in] lab_dirX label for Dirichelet boundary conditions (x-direction)
!> @param[in] nb_diriY number of Dirichlet boundary conditions (x-direction)
!> @param[in] nb_diriZ number of Dirichlet boundary conditions (x-direction)
!> @param[in]  lab_dirY label for Dirichelet boundary conditions (y-direction)
!> @param[in]  lab_dirZ label for Dirichelet boundary conditions (z-direction) 
!> @param[in]  nb_neuX  number of Neumann boundary conditions (x-direction)
!> @param[in] lab_neuX  label for Neumann boundary conditions (x-direction)
!> @param[in]  nb_neuY  number of Neumann boundary conditions (y-direction)
!> @param[in] lab_neuY  label for Neumann boundary conditions (y-direction)
!> @param[in]  nb_neuZ  number of Neumann boundary conditions (z-direction)
!> @param[in] lab_neuZ  label for Neumann boundary conditions (z-direction)
!> @param[in] nb_neuN   number of Neumann boundary conditions (N-direction)
!> @param[in] lab_neuN  label for Neumann boundary conditions (z-direction)
!> @param[in] nb_abc  number of absorbing boundary conditions
!> @param[in] lab_abc label for absorbing boundary conditions
!> @param[in] nb_dg  number of Discontinuous Galerkin interface conditions
!> @param[in] lab_dg  label for Discontinuous Galerkin interface conditions
!> @param[in] nb_node  nunmber of grid nodes
!> @param[out]  xx,yy,zz  real numbers (dummy)
!> @param[in] nb_hexa  number of hexes in the mesh
!> @param[out] con_hexa  matrix connectivity for the hexes of the mesh (1-el, 2-,...,9-nodes)
!> @param[in] nb_quad  number of quads in the mesh
!> @param[out] con_quad  matrix connectivity for the quads of the mesh (1-el, 2-,...,5-nodes)

      subroutine READ_FILEMESH(gridfile,nb_mat,lab_mat,&
                              nb_diriX,lab_dirX,nb_diriY,lab_dirY,nb_diriZ,lab_dirZ,&
                              nb_neuX,lab_neuX,nb_neuY,lab_neuY,nb_neuZ,lab_neuZ,&
                              nb_neuN,lab_neuN,&                  
                              nb_abc,lab_abc,&
                              nb_dg,lab_dg,&
                              nb_node,xx,yy,zz,nb_hexa,con_hexa,nb_quad,con_quad)
      
      implicit none
      
      character*110 :: inline
      character*20 :: el_code
      character*70 :: gridfile

      integer*4 :: nb_mat,nb_diriX,nb_diriY,nb_diriZ,nb_neuX,nb_neuY,nb_neuZ,nb_abc,nb_dg
      integer*4 :: nb_neuN
      integer*4 :: nb_node,nb_hexa,nb_quad
      integer*4 :: inode,ihexa,iquad
      integer*4 :: nb_elem,ie,i,j,mat_code,ileft,iright,sl,trash,status,control

      integer*4, dimension(nb_mat) :: lab_mat
      integer*4, dimension(nb_diriX) :: lab_dirX
      integer*4, dimension(nb_diriY) :: lab_dirY
      integer*4, dimension(nb_diriZ) :: lab_dirZ
      integer*4, dimension(nb_neuX) :: lab_neuX
      integer*4, dimension(nb_neuY) :: lab_neuY
      integer*4, dimension(nb_neuZ) :: lab_neuZ
      integer*4, dimension(nb_neuN) :: lab_neuN
      integer*4, dimension(nb_abc) :: lab_abc
      integer*4, dimension(nb_dg) :: lab_dg

      integer*4, dimension(nb_hexa,9) :: con_hexa
      integer*4, dimension(nb_quad,5) :: con_quad
      
      real*8 :: xx,yy,zz
      
      inode = 0
      ihexa = 0
      iquad = 0
      
      status = 0 
      
      open(40,file=gridfile)
      
      do 
         read(40,'(A)') inline
         if (inline(1:1) .ne. '#') exit
      enddo
      
      read(inline,*)nb_node,nb_elem
      
      do i = 1,nb_node
        read(40,*)inode,xx,yy,zz

        if (inode.ne.i) then
          status = 1
        endif
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
            do i = 1,nb_mat
               if (lab_mat(i).eq.mat_code) control = 1
            enddo
            !
            if (control.ne.0) then
               ihexa = ihexa + 1
               con_hexa(ihexa,1) = mat_code
               read(inline(iright:sl),*)(con_hexa(ihexa,j),j=2,9)
            endif
        !
         elseif ((el_code.eq.'quad').or.(el_code.eq.'QUAD')) then
            control = 0
            do i = 1,nb_diriX
               if (lab_dirX(i).eq.mat_code) control = 1
            enddo
            do i = 1,nb_diriY
               if (lab_dirY(i).eq.mat_code) control = 1
            enddo
            do i = 1,nb_diriZ
               if (lab_dirZ(i).eq.mat_code) control = 1
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

         
            if (control .ne. 0) then
               iquad = iquad +1
               con_quad(iquad,1) = mat_code
               read(inline(iright:sl),*)(con_quad(iquad,j),j=2,5)
            endif
         
         endif
      enddo
      
      close(40)
      
      return
      
      end subroutine READ_FILEMESH

