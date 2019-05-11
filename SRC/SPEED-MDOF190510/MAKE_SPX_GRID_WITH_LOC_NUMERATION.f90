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

!> @brief Makes local spectral grids with nodes numbered according
!! to the new local numeration. 
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0

     subroutine MAKE_SPX_GRID_WITH_LOC_NUMERATION()
     
      use max_var
      use speed_par

      implicit none
 
      include 'SPEED.MPI'

!***********************************************************************************
!     LOCAL NODES RENUMBERING
!***********************************************************************************
                   
     call MAKE_LOC_NODES_NUMERATION(nnod_loc, local_node_num, nrecv, node_recv)            
    

     call MAKE_SPX_LOC_NUMERATION(nnod_loc, local_node_num, &
                            con_nnz_loc, con_spx_loc, &
                            con_nnz_bc_loc, con_spx_bc_loc)

     if (nnod_loc.gt.0) then
          allocate (xx_spx_loc(nnod_loc), yy_spx_loc(nnod_loc), zz_spx_loc(nnod_loc))
     else 
         write(*,'(A)')'Error ! nnod_loc = 0'
         call EXIT(EXIT_NO_ELEMENTS)
     endif

      
     call READ_FILEMESH_LOC(grid_file, nnod_loc, xx_spx_loc, yy_spx_loc, zz_spx_loc, &
                            local_node_num)
      
     if(mpi_id .eq.0) write(*,'(A)')
     if(mpi_id.eq.0) write(*,'(A,I8)') 'Total Spectral Nodes : ',nnod

     
     allocate(alfa11(nelem_loc),alfa12(nelem_loc),alfa13(nelem_loc))
     allocate(alfa21(nelem_loc),alfa22(nelem_loc),alfa23(nelem_loc))
     allocate(alfa31(nelem_loc),alfa32(nelem_loc),alfa33(nelem_loc))
     allocate(beta11(nelem_loc),beta12(nelem_loc),beta13(nelem_loc))
     allocate(beta21(nelem_loc),beta22(nelem_loc),beta23(nelem_loc))
     allocate(beta31(nelem_loc),beta32(nelem_loc),beta33(nelem_loc))
     allocate(gamma1(nelem_loc),gamma2(nelem_loc),gamma3(nelem_loc))
     allocate(delta1(nelem_loc),delta2(nelem_loc),delta3(nelem_loc))
           
     
     call MAKE_SPX_GRID_LOC(nnod_loc,xx_spx_loc,yy_spx_loc,zz_spx_loc,&
                            con_nnz_loc, con_spx_loc, nmat, tag_mat, sdeg_mat, nelem_loc, &
                            local_node_num, &
                            alfa11,alfa12,alfa13, &
                            alfa21,alfa22,alfa23, &
                            alfa31,alfa32,alfa33, &
                            beta11,beta12,beta13, &
                            beta21,beta22,beta23, &
                            beta31,beta32,beta33, &
                            gamma1,gamma2,gamma3, &
                            delta1,delta2,delta3 )


     
     end subroutine MAKE_SPX_GRID_WITH_LOC_NUMERATION
