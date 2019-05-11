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


!> @brief Find non-linear nodes.
!! @author Ilario Mazzieri
!> @date September, 2015
!> @version 1.0
!> @param[in] Depth depth for searching nodes
!> @param[in] nn polynomial degree
!> @param[in] nnt number of nodes
!> @param[in]  zs_elv nodes elevation from surface
!> @param[in] cs_nnz length of cs
!> @param[in] cs connectivity vector
!> @param[in] ielem  element index 
!> @param[in] ncase  id for case
!> @param[in] tcase  label for case  
!> @param[in] zs_all nodes elevation from alluvial
!> @param[in] vs_tria vs30 values
!> @param[out] yon 1 if the node is non-linear elastic, 0 otherwise
!
!************************************************************************************************** 

      subroutine MAKE_NLE(Depth,nn,&                                                                 
                          yon,&                                                                
                          nnt,zs_elev,&                                                
                          cs_nnz,cs,ielem,ncase,&                
                          tcase,zs_all,vs_tria)                          
      
   
     
      implicit none
                                                              
      integer*4 :: nn
      integer*4 :: p,q,r
      integer*4 :: nnt
      integer*4 :: cs_nnz                                
      integer*4 :: ncase
      integer*4 :: tcase                                                                        
      integer*4 :: is,in,ielem

      integer*4, dimension(0:cs_nnz) :: cs                

      integer*4, dimension(nn,nn,nn) :: yon

      real*8 :: Depth
          
      real*8, dimension(nnt) :: zs_elev
      real*8, dimension(nnt) :: zs_all
      real*8, dimension(nnt) :: vs_tria                 
      
!     STRESS CALCULATION
      
      do r = 1,nn
         do q = 1,nn
            do p = 1,nn

                yon(p,q,r) = 1

                is = nn*nn*(r -1) +nn*(q -1) +p
               in = cs(cs(ielem -1) +is)
               if (ncase.gt.0 .and. tcase .ne. 16) then
                        if (Depth.lt.zs_elev(in) .or. zs_elev(in) .lt. 0.d0 ) yon(p,q,r) = 0
                        if ((tcase .gt. 0) .and. (zs_all(in) .lt. 0.0d0)) yon(p,q,r) = 0
               elseif (ncase.gt.0 .and. tcase .eq. 16) then
                        if (Depth .lt. zs_elev(in) .or. zs_elev(in) .lt. 0.d0 ) yon(p,q,r) = 0
                        if (vs_tria(in) .ge. 450.0d0) yon(p,q,r) = 0
               endif        

            enddo
         enddo
      enddo
      
      return
      
      end subroutine MAKE_NLE

