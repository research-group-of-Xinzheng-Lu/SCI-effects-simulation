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

!> @brief Counts number of DG elements (local and global)
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nm = number of materials
!> @param[in] sd(nm) = spectral degree vector
!> @param[in] tag_mat(nm) = tag for materials
!> @param[in] cs_nnz_loc length of cs_loc
!> @param[in] cs_loc local spectral connectivity vector
!> @param[in] nn_loc number of local nodes 
!> @param[in] local_n_num local node numeration
!> @param[in] ne_loc number of local elements
!> @param[in] local_el_num local element numeration
!> @param[in] xs x-coord local spectral nodes 
!> @param[in] ys y-coord local spectral nodes 
!> @param[in] zs z-coord local spectral nodes 
!> @param[in] i4count  vector identifying nodes where DG conditions are applied
!> @param[in] mpi_id  MPI processor identity
!> @param[in] mpi_comm  MPI common world
!> @param[out] nel_dg_loc  number of local dg elements
!> @param[out] nel_dg_glo  number of global dg elements

      subroutine GET_DIME_DG(nm, sd, tag_mat, cs_nnz_loc, cs_loc, &
                           nn_loc, local_n_num, ne_loc, local_el_num, &
                           xs,ys,zs,&
                           nel_dg_loc, nel_dg_glo, &
                           i4count, mpi_id, mpi_comm)

     implicit none

     include 'SPEED.MPI'

     integer*4 :: nm, cs_nnz_loc, nn_loc, ne_loc, nel_dg_loc, nel_dg_glo
     integer*4 :: im, nn, ie, ned, err_out
     integer*4 :: ne1, ne2, ne3, ne4, ic1, ic2, ic3, ic4
     integer*4 :: ne5, ne6, ne7, ne8, ic5, ic6, ic7, ic8
     integer*4 :: mpi_comm, mpi_id, mpi_ierr

     integer*4, dimension(nm) :: tag_mat, sd
     integer*4, dimension(0:cs_nnz_loc) :: cs_loc
     integer*4, dimension(nn_loc) :: local_n_num, i4count
     integer*4, dimension(ne_loc) :: local_el_num

     real*8, dimension(nn_loc) :: xs,ys,zs
     
!************************************************************************************************
!                        conto il numero di elemnti dg 
!************************************************************************************************   

     nel_dg_loc = 0
     ned = cs_loc(0) - 1

      do im = 1,nm
      
         nn = sd(im) +1         

         do ie = 1,ned
            if (cs_loc(cs_loc(ie -1) + 0) .eq. tag_mat(im)) then
               
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(1 -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(1 -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(1 -1) +1)
                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then
                  nel_dg_loc = nel_dg_loc +1
               endif
               
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(1 -1) +1)
                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then
                  nel_dg_loc = nel_dg_loc +1
               endif
               
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(1 -1) +1)
                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then
                  nel_dg_loc = nel_dg_loc +1
               endif
               
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(nn -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(nn -1) +1)
                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then
                  nel_dg_loc = nel_dg_loc +1
               endif
               
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(1 -1) +1)

                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then
                  nel_dg_loc = nel_dg_loc +1
               endif
               
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(1 -1) +1)

                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) .and. (i4count(ne3).ne.0) .and. (i4count(ne4).ne.0)) then
                  nel_dg_loc = nel_dg_loc +1
               endif
               
            endif
         enddo
      enddo    
      
      
      call MPI_ALLREDUCE(nel_dg_loc, nel_dg_glo, 1, SPEED_INTEGER, MPI_SUM, mpi_comm, mpi_ierr)
      
      


      end subroutine GET_DIME_DG

