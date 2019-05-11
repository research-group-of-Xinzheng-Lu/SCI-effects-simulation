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

!> @brief Makes mesh partitioning using METIS.
!> @date September, 2013 
!> @version 1.0
!> @warning Mesh partitioning is made using metis-4.0.3 library. At the moment
!! we have problems in using metis-5.0.2 library. 
!> @todo Switch to other mesh partitioner or to ParMetis library.

!> @param[in] mpi_file folder name given in SPEED.input (MPIFILE keyword) 
!> @param[in] nelem  number of elements
!> @param[in] nnode  number of mesh nodes
!> @param[in] nparts  number of processors (number of parts of the partition)
!> @param[in]  conn  connectivity matrix: conn(i,1) = number of element i
!!                                  conn(i,2),...,conn(i,9) = nodes of the hex
!> @param[in]  w  control parameter 1 - different weights for the elements 
!! 0 - same weight for all elements (for metis-5.0.2 version)
!> @param[out] elemdomain.mpi file containing elements partitioning

      subroutine MESH_PARTITIONING(mpi_file,nelem,nnode,nparts,conn,w_yn)

 
      implicit none
      

      character*70 :: mpi_file,u_name
      integer*4 :: ie,i,ic,u_mpi, w_yn, trash
      integer*4 :: ncommon, objval, edgecut

      integer*4, intent(inout) :: nelem, nnode, nparts
      integer*4, dimension(nelem,9), intent(in) :: conn

      integer*4, dimension(nnode) :: npart
      integer*4, dimension(nelem) :: epart
      integer*4, dimension(nelem) :: vwgt

      integer*4, pointer :: vsize(:) => NULL() 

      integer*4, dimension(nelem+1) :: eptr
      integer*4, dimension(nelem*8) :: eind
            
      integer*4 :: options(40) ! see METIS_NOPTIONS in metis.h
      integer*4, parameter :: METIS_OPTION_PTYPE = 1, METIS_OPTION_OBJTYPE = 2, METIS_OPTION_CTYPE = 3, &
                                METIS_OPTION_IPTYPE = 4, METIS_OPTION_RTYPE = 5, METIS_OPTION_DBGLVL = 6, &
                              METIS_OPTION_NITER = 7, METIS_OPTION_NCUTS = 8, METIS_OPTION_SEED = 9, &
                              METIS_OPTION_MINCONN = 10, METIS_OPTION_CONTIG = 11, METIS_OPTION_COMPRESS = 12, &
                              METIS_OPTION_CCORDER = 13, METIS_OPTION_PFACTOR = 14, METIS_OPTION_NSEPS = 15, &
                              METIS_OPTION_UFACTOR = 16, METIS_OPTION_NUMBERING = 17

      real*8, pointer :: tpwgts(:) => NULL() 

      eptr(1) = 1
      ic = 1
      do ie = 1, nelem
         eptr(ie+1) = eptr(ie) + 8 
         do i = 2, 9
            eind(ic) = conn(ie,i)
            ic = ic +1
         enddo
      enddo

     eptr = eptr -1
     eind = eind - 1

     ncommon = 4 

!     METIS-5.0.2
!     if(w_yn .eq. 0) then
!        vwgt = 1
!     else
!        u_name = 'element.wgt'
!        u_mpi = 40                                 
!        open(u_mpi,file=u_name)        
!        read(u_mpi,*) nelem
!        
!        do i = 1, nelem
!               read(u_mpi,*) trash, vwgt(i)
!        enddo
!        close(u_mpi)
!     endif     

!     METIS-5.0.2
!     Setting the Fortran numbering scheme
!     call METIS_SetDefaultOptions(options)  !INPUT
     
                              
!     call METIS_PartMeshDual(nelem, nnode, eptr, eind, vwgt, vsize, &
!                              ncommon, nparts, tpwgts, options, objval, &
!                              epart, npart)                        


!     METIS-4.0.3
      eind = eind +1  
      call METIS_PartMeshDual(nelem,nnode,eind,3,1,nparts,edgecut,&
                              epart,npart)
                        
      epart = epart - 1             



      if(w_yn .eq. 0) then
         if(len_trim(mpi_file) .ne. 70) then                                                                                  
            u_name = mpi_file(1:len_trim(mpi_file)) // '/elemdomain.mpi'
         else 
            u_name = 'elemdomain.mpi'     
         endif      
   
         u_mpi = 400                                 
         open(u_mpi,file=u_name)        
         write(u_mpi,*) nelem
         do i = 1, nelem
               write(u_mpi,*) i, epart(i)
         enddo
         close(u_mpi)

!       else
!
!         if(len_trim(mpi_file) .ne. 70) then                                                                                  
!            u_name = mpi_file(1:len_trim(mpi_file)) // '/final.part'
!         else 
!            u_name = 'final.part'     
!         endif
!         
!         u_mpi = 400                                 
!         open(u_mpi,file=u_name)        
!         write(u_mpi,*) nelem
!         do i = 1, nelem
!               write(u_mpi,*) i, epart(i)
!         enddo
!         close(u_mpi)

       endif       
      
      return
      
      end subroutine MESH_PARTITIONING
