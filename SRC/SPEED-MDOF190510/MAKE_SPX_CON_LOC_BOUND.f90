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

!> @brief Makes spectral connectivity vector for the boundary.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] cs_nnz length of cs
!> @param[in] cs spectral connectivity vector
!> @param[in] ne_bc number of boundary element
!> @param[in] cm_bc connectivity matrix for boundary conditions
!> @param[in] nm number of materials  
!> @param[in] tm material labels
!> @param[in] sd polynomial degree vector
!> @param[in] Ennz length of Ebin
!> @param[in] Ebin pointer for connectivity (see MAKE_GRID_NODES.f90)
!> @param[in] nel_loc number of local elements
!> @param[in] local_el_num local element numeration
!> @param[in] nf_loc dummy
!> @param[in] cs_nnz_bc length of cs_bc 
!> @param[out] cs_bc local spectral boundary connectivity vector 


      subroutine MAKE_SPX_CON_LOC_BOUND(cs_nnz,cs,ne_bc,cm_bc,nm,tm,sd,&
                                        Ennz,Ebin,cs_nnz_bc,cs_bc,loc_el_num, nel_loc, &
                                        nf_loc)
      
      implicit none
      
      integer*4 :: cs_nnz,ne_bc,nm,Ennz,cs_nnz_bc,nel_loc,nf_loc
      integer*4 :: im,ibc,ielem,nn, icur, iloc
      integer*4 :: i,j,k,an,bn,cn,ia,ib,ic,ja,jb,jc,ka,kb,kc
      integer*4 :: i1,i2,j1,j2,k1,k2,is,js

      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(nm) :: tm 
      integer*4, dimension(nm) :: sd 
      integer*4, dimension(0:Ennz) :: Ebin
      integer*4, dimension(0:cs_nnz_bc) :: cs_bc
      integer*4, dimension(nel_loc) :: loc_el_num
      
      integer*4, dimension(ne_bc,5) :: cm_bc
      
      
      iloc = 0
      
      do ibc = 1,ne_bc
         an = cm_bc(ibc,1 +1)
         bn = cm_bc(ibc,2 +1)
         cn = cm_bc(ibc,4 +1)
         
         call GET_ELEM_FROM_FACE(Ennz,Ebin,an,bn,cn,ielem)         
         call GET_INDLOC_FROM_INDGLO(loc_el_num, nel_loc, ielem, icur)
         
         
         if (icur .ne. 0) then
         iloc = iloc + 1
         do im = 1,nm
            if (tm(im).eq.cs(cs(icur -1))) then
                nn = sd(im) +1
                cs_bc(iloc) = cs_bc(iloc -1) +nn*nn +1
                cs_bc(cs_bc(iloc -1)) = cm_bc(ibc,1)         
            endif
         enddo
         endif
      enddo
      

      
      iloc = 0
      do ibc = 1,ne_bc
         an = cm_bc(ibc,1 +1)
         bn = cm_bc(ibc,2 +1)
         cn = cm_bc(ibc,4 +1)
         
       
         call GET_ELEM_FROM_FACE(Ennz,Ebin,an,bn,cn,ielem)
         call GET_INDLOC_FROM_INDGLO(loc_el_num, nel_loc, ielem, icur)


         if(icur.ne.0) then
           iloc = iloc + 1
           do im = 1,nm
              if (tm(im).eq.cs(cs(icur -1) +0)) nn = sd(im) +1
           enddo
         
           ia = 0
           ja = 0
           ka = 0
           ib = 0
           jb = 0
           kb = 0
           ic = 0
           jc = 0
           kc = 0
         
           do k = 1,nn
              do j = 1,nn
                 do i = 1,nn
                    js = nn*nn*(k -1) +nn*(j -1) +i
                    
                    if (cs(cs(icur -1) +js).eq.an) then
                       ia = i
                       ja = j
                       ka = k
                    endif
                  
                    if (cs(cs(icur -1) +js).eq.bn) then
                       ib = i
                       jb = j
                       kb = k
                    endif
                  
                    if (cs(cs(icur -1) +js).eq.cn) then
                       ic = i
                       jc = j
                       kc = k
                    endif
                 enddo
              enddo
           enddo
         
           i1 = min(ia,ib,ic)
           i2 = max(ia,ib,ic)
           j1 = min(ja,jb,jc)
           j2 = max(ja,jb,jc)
           k1 = min(ka,kb,kc)
           k2 = max(ka,kb,kc)
         
           is = 1
           do k = k1,k2
              do j = j1,j2
                 do i = i1,i2
                    js = nn*nn*(k -1) +nn*(j -1) +i
                     
                    cs_bc(cs_bc(iloc -1) +is) = cs(cs(icur -1) +js)
                    is = is +1
                 enddo
              enddo
            enddo
           endif        

       enddo
        
        return
       
      end subroutine MAKE_SPX_CON_LOC_BOUND

