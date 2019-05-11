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


!> @brief Set boundary conditions.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0

      subroutine MAKE_BOUNDARY_CONDITIONS() 
        
              
      use max_var
      use speed_par
      

      implicit none
 
      include 'SPEED.MPI' 

!*****************************************************************************************
!                         SETUP DIRICHELET BOUNDARY CONDITIONS    
!*****************************************************************************************
     
       if (mpi_id.eq.0) write(*,'(A)')
       if (mpi_id.eq.0) write(*,'(A)')'---------------Setup boundary conditions---------------'                  

       allocate(i4count(nnod_loc))
       i4count = 0
       nnode_dirX = 0
      
       call GET_NODE_FROM_FACE(nnod_loc,con_nnz_bc_loc, con_spx_bc_loc, nload_dirX_el,tag_dirX_el,&
                          nnode_dirX, i4count, local_node_num)

       if(nnode_dirX .gt. 0) then
         allocate(inode_dirX(nnode_dirX))
         i = 1         
         do in = 1, nnod_loc
            if (i4count(in) .ne. 0) then
               inode_dirX(i) = in
               i = i + 1
            endif
         enddo
      endif            
      
       i4count = 0
       nnode_dirY = 0
      
       call GET_NODE_FROM_FACE(nnod_loc, con_nnz_bc_loc, con_spx_bc_loc, nload_dirY_el,tag_dirY_el,&
                          nnode_dirY, i4count, local_node_num)

       if(nnode_dirY .gt. 0) then
         allocate(inode_dirY(nnode_dirY))
                i = 1         
         do in = 1, nnod_loc
            if (i4count(in) .ne. 0) then
               inode_dirY(i) = in
               i = i + 1
            endif
         enddo
      endif            


       i4count = 0
       nnode_dirZ = 0
      
       call GET_NODE_FROM_FACE(nnod_loc,con_nnz_bc_loc, con_spx_bc_loc, nload_dirZ_el, tag_dirZ_el,&
                           nnode_dirZ, i4count, local_node_num)

       if(nnode_dirZ .gt. 0) then
               allocate(inode_dirZ(nnode_dirZ))
                i = 1         
         do in = 1, nnod_loc
            if (i4count(in) .ne. 0) then
               inode_dirZ(i) = in
               i = i + 1
            endif
         enddo
      endif            

      deallocate(i4count)     


!*****************************************************************************************
!                            SETUP ABC BOUNDARY CONDITIONS
!*****************************************************************************************

      nelem_abc = 0     
      nnode_abc = 0
      allocate(i4count(nnod_loc))
      i4count = 0


      call GET_NODE_FROM_FACE(nnod_loc, con_nnz_bc_loc, con_spx_bc_loc, nload_abc_el, tag_abc_el,&
                          nnode_abc, i4count, local_node_num)

              
      call GET_DIME_ABC(nmat, sdeg_mat, tag_mat, con_nnz_loc, con_spx_loc, &
                    nnod_loc, local_node_num, &
                    nelem_abc, i4count) 
 
     if(nelem_abc .gt. 0) then
        allocate(ielem_abc(nelem_abc,7)) 

        call SETUP_ABC(nmat, sdeg_mat, tag_mat, con_nnz_loc, con_spx_loc, &
                       nnod_loc, local_node_num, nelem_loc, local_el_num, &
                       nelem_abc, ielem_abc, i4count)

     endif
              

      deallocate(i4count)
      
      if (mpi_id.eq.0) write(*,'(A)')'Made.'
      if (mpi_id.eq.0) write(*,'(A)')                  




      end subroutine MAKE_BOUNDARY_CONDITIONS
