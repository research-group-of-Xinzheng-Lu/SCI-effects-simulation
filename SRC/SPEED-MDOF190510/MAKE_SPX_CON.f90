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

!> @brief Makes spectral connectivity vector.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nelem     number of elements
!> @param[in] con_mac   connectivity matrix for spectral nodes  
!> @param[in] tag_mat   material label
!> @param[in] sdeg      polynomial degree vector 
!> @param[in] node_pntr pointer for connectivity (see MAKE_GRID_NODES.f90) 
!> @param[out] con_spx spectral connectivity vector 
!> @param[out] nnode  number of spectral nodes


      subroutine MAKE_SPX_CON(nelem,con_mac,nmat,tag_mat,sdeg,&
                                            nnz_pntr,node_pntr,con_nnz,con_spx,nnode)
 
      implicit none
      
      integer*4 :: nelem,nmat,nnz_pntr,con_nnz,nnode
      integer*4, dimension(nelem,9) :: con_mac
      integer*4, dimension(nmat) :: tag_mat
      integer*4, dimension(nmat) :: sdeg
      integer*4, dimension(0:nnz_pntr) :: node_pntr
      integer*4, dimension(0:con_nnz) :: con_spx
      
      integer*4 :: nnode_mac,nn,imat,ie,je
      integer*4 :: i,j,k,an,bn,cn,ia,ib,ic,ja,jb,jc,ka,kb,kc
      integer*4 :: jsa,dljs,dmjs,l,m,is,js
      
      
      do i = 0,con_nnz
         con_spx(i) = 0
      enddo
      
      con_spx(0) = nelem +1
      do ie = 1,nelem
         do imat = 1,nmat
            if (con_mac(ie,1).eq.tag_mat(imat)) then
               nn = sdeg(imat) +1
               con_spx(ie) = con_spx(ie -1) + nn*nn*nn +1
            endif
         enddo
      enddo
      
      nnode_mac = node_pntr(0) -1
      
      nnode = nnode_mac
      
      do ie = 1,nelem
         do imat = 1,nmat
            if (con_mac(ie,1).eq.tag_mat(imat)) then
               nn = sdeg(imat) +1
               
! First put material index
               con_spx(con_spx(ie -1) +0) = con_mac(ie,1)
               
! Then put vertices
               con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +1) = &
                    con_mac(ie,2)
               con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +nn) = &
                    con_mac(ie,3)
               con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +nn) = &
                    con_mac(ie,4)
               con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +1) = &
                    con_mac(ie,5)
               con_spx(con_spx(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +1) = &
                    con_mac(ie,6)
               con_spx(con_spx(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +nn) = &
                    con_mac(ie,7)
               con_spx(con_spx(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +nn) = &
                    con_mac(ie,8)
               con_spx(con_spx(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +1) = &
                    con_mac(ie,9)
               
                           
! Then construct edge connectivity
               
! First edge down
               
               an = con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +1)
               bn = con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +nn)
               
               je = ie
               do i = node_pntr(an -1),node_pntr(an) -1
                  do j = node_pntr(bn -1),node_pntr(bn) -1
                     if (node_pntr(i).eq.node_pntr(j)) then
                        if ((node_pntr(i).lt.ie).and.((nn*nn*nn +1)&
                        .eq.(con_spx(node_pntr(i)) - con_spx(node_pntr(i) -1)))) then
                           je = node_pntr(i)
                        endif
                     endif
                  enddo
               enddo
               
               if (je.ne.ie) then
                  ia = 0
                  ja = 0
                  ka = 0
                  ib = 0
                  jb = 0
                  kb = 0
                  
                  do k = 1,nn
                     do j = 1,nn
                        do i = 1,nn
                           js = nn*nn*(k -1) +nn*(j -1) +i
                           
                           if (con_spx(con_spx(je -1) +js).eq.an) then
                              ia = i
                              ja = j
                              ka = k
                           endif
                           
                           if (con_spx(con_spx(je -1) +js).eq.bn) then
                              ib = i
                              jb = j
                              kb = k
                           endif
                        enddo
                     enddo
                  enddo

                  jsa = nn*nn*(ka -1) +nn*(ja -1) +(ia -1) +1
                  dljs = (nn*nn*(kb -ka) +nn*(jb -ja) +(ib -ia)) /(nn -1)
                  
                  do l = 1,nn
                     is = nn*nn*(1 -1) +nn*(1 -1) +l
                     js = jsa +dljs*(l -1)
                     con_spx(con_spx(ie -1) +is) = &
                          con_spx(con_spx(je -1) +js)
                  enddo
                  
               else
                  do i = 2,(nn -1)
                     nnode = nnode +1
                     is = nn*nn*(1 -1) +nn*(1 -1) +i
                     con_spx(con_spx(ie -1) + is) = nnode
                  enddo
               endif
               
               
! Second edge down
               
               an = con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +nn)
               bn = con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +nn)
               
               je = ie
               do i = node_pntr(an -1),node_pntr(an) -1
                  do j = node_pntr(bn -1),node_pntr(bn) -1
                     if (node_pntr(i).eq.node_pntr(j)) then
                        if ((node_pntr(i).lt.ie).and.((nn*nn*nn +1)&
                        .eq.(con_spx(node_pntr(i)) - con_spx(node_pntr(i) -1)))) then
                           je = node_pntr(i)
                        endif
                     endif
                  enddo
               enddo
               
               if (je.ne.ie) then
                  ia = 0
                  ja = 0
                  ka = 0
                  ib = 0
                  jb = 0
                  kb = 0
                  
                  do k = 1,nn
                     do j = 1,nn
                        do i = 1,nn
                           js = nn*nn*(k -1) +nn*(j -1) +i
                           
                           if (con_spx(con_spx(je -1) +js).eq.an) then
                              ia = i
                              ja = j
                              ka = k
                           endif
                           
                           if (con_spx(con_spx(je -1) +js).eq.bn) then
                              ib = i
                              jb = j
                              kb = k
                           endif
                        enddo
                     enddo
                  enddo
                  
                  jsa = nn*nn*(ka -1) +nn*(ja -1) +(ia -1) +1
                  dljs = (nn*nn*(kb -ka) +nn*(jb -ja) +(ib -ia)) /(nn -1)
                  
                  do l = 1,nn
                     is = nn*nn*(1 -1) +nn*(l -1) +nn
                     js = jsa +dljs*(l -1)
                     con_spx(con_spx(ie -1) +is) = &
                          con_spx(con_spx(je -1) +js)
                  enddo
                  
               else
                  do j = 2,(nn -1)
                     nnode = nnode +1
                     is = nn*nn*(1 -1) +nn*(j -1) +nn
                     con_spx(con_spx(ie -1) + is) = nnode
                  enddo
               endif
               
               
! Third edge down
               
               an = con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +1)
               bn = con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +nn)
               
               je = ie
               do i = node_pntr(an -1),node_pntr(an) -1
                  do j = node_pntr(bn -1),node_pntr(bn) -1
                     if (node_pntr(i).eq.node_pntr(j)) then
                        if ((node_pntr(i).lt.ie).and.((nn*nn*nn +1)&
                        .eq.(con_spx(node_pntr(i)) - con_spx(node_pntr(i) -1)))) then
                           je = node_pntr(i)
                        endif
                     endif
                  enddo
               enddo
               
              
               if (je.ne.ie) then
                  ia = 0
                  ja = 0
                  ka = 0
                  ib = 0
                  jb = 0
                  kb = 0
                  
                  do k = 1,nn
                     do j = 1,nn
                        do i = 1,nn
                           js = nn*nn*(k -1) +nn*(j -1) +i
                           
                           if (con_spx(con_spx(je -1) +js).eq.an) then
                              ia = i
                              ja = j
                              ka = k
                           endif
                           
                           if (con_spx(con_spx(je -1) +js).eq.bn) then
                              ib = i
                              jb = j
                              kb = k
                           endif
                        enddo
                     enddo
                  enddo
                  
                  jsa = nn*nn*(ka -1) +nn*(ja -1) +(ia -1) +1
                  dljs = (nn*nn*(kb -ka) +nn*(jb -ja) +(ib -ia)) /(nn -1)
                  
                  do l = 1,nn
                     is = nn*nn*(1 -1) +nn*(nn -1) +l
                     js = jsa +dljs*(l -1)
                     con_spx(con_spx(ie -1) +is) = &
                          con_spx(con_spx(je -1) +js)
                  enddo
                  
               else
                  do i = 2,(nn -1)
                     nnode = nnode +1
                     is = nn*nn*(1 -1) +nn*(nn -1) +i
                     con_spx(con_spx(ie -1) + is) = nnode
                  enddo
               endif
               
               
! Fourth edge down
               
               an = con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +1)
               bn = con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +1)
               
               je = ie
               do i = node_pntr(an -1),node_pntr(an) -1
                  do j = node_pntr(bn -1),node_pntr(bn) -1
                     if (node_pntr(i).eq.node_pntr(j)) then
                        if ((node_pntr(i).lt.ie).and.((nn*nn*nn +1)&
                        .eq.(con_spx(node_pntr(i)) - con_spx(node_pntr(i) -1)))) then
                           je = node_pntr(i)
                        endif
                     endif
                  enddo
               enddo
                              
               if (je.ne.ie) then
                  ia = 0
                  ja = 0
                  ka = 0
                  ib = 0
                  jb = 0
                  kb = 0
                  
                  do k = 1,nn
                     do j = 1,nn
                        do i = 1,nn
                           js = nn*nn*(k -1) +nn*(j -1) +i
                           
                           if (con_spx(con_spx(je -1) +js).eq.an) then
                              ia = i
                              ja = j
                              ka = k
                           endif
                           
                           if (con_spx(con_spx(je -1) +js).eq.bn) then
                              ib = i
                              jb = j
                              kb = k
                           endif
                        enddo
                     enddo
                  enddo
                  
                  jsa = nn*nn*(ka -1) +nn*(ja -1) +(ia -1) +1
                  dljs = (nn*nn*(kb -ka) +nn*(jb -ja) +(ib -ia)) /(nn -1)
                  
                  do l = 1,nn
                     is = nn*nn*(1 -1) +nn*(l -1) +1
                     js = jsa +dljs*(l -1)
                     con_spx(con_spx(ie -1) +is) = &
                          con_spx(con_spx(je -1) +js)
                  enddo
                  
               else
                  do j = 2,(nn -1)
                     nnode = nnode +1
                     is = nn*nn*(1 -1) +nn*(j -1) +1
                     con_spx(con_spx(ie -1) + is) = nnode
                  enddo
               endif
               
               
! First edge side
               
               an = con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +1)
               bn = con_spx(con_spx(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +1)
               
               je = ie
               do i = node_pntr(an -1),node_pntr(an) -1
                  do j = node_pntr(bn -1),node_pntr(bn) -1
                     if (node_pntr(i).eq.node_pntr(j)) then
                        if ((node_pntr(i).lt.ie).and.((nn*nn*nn +1)&
                        .eq.(con_spx(node_pntr(i)) - con_spx(node_pntr(i) -1)))) then
                           je = node_pntr(i)
                        endif
                     endif
                  enddo
               enddo
               
               if (je.ne.ie) then
                  ia = 0
                  ja = 0
                  ka = 0
                  ib = 0
                  jb = 0
                  kb = 0
                  
                  do k = 1,nn
                     do j = 1,nn
                        do i = 1,nn
                           js = nn*nn*(k -1) +nn*(j -1) +i
                           
                           if (con_spx(con_spx(je -1) +js).eq.an) then
                              ia = i
                              ja = j
                              ka = k
                           endif
                           
                           if (con_spx(con_spx(je -1) +js).eq.bn) then
                              ib = i
                              jb = j
                              kb = k
                           endif
                        enddo
                     enddo
                  enddo
                  
                  jsa = nn*nn*(ka -1) +nn*(ja -1) +(ia -1) +1
                  dljs = (nn*nn*(kb -ka) +nn*(jb -ja) +(ib -ia)) /(nn -1)
                  
                  do l = 1,nn
                     is = nn*nn*(l -1) +nn*(1 -1) +1
                     js = jsa +dljs*(l -1)
                     con_spx(con_spx(ie -1) +is) = &
                          con_spx(con_spx(je -1) +js)
                  enddo
                  
               else
                  do k = 2,(nn -1)
                     nnode = nnode +1
                     is = nn*nn*(k -1) +nn*(1 -1) +1
                     con_spx(con_spx(ie -1) + is) = nnode                    
                  enddo
               endif
               
               
! Second edge side
               
               an = con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +nn)
               bn = con_spx(con_spx(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +nn)
               
               je = ie
               do i = node_pntr(an -1),node_pntr(an) -1
                  do j = node_pntr(bn -1),node_pntr(bn) -1
                     if (node_pntr(i).eq.node_pntr(j)) then
                        if ((node_pntr(i).lt.ie).and.((nn*nn*nn +1)&
                        .eq.(con_spx(node_pntr(i)) - con_spx(node_pntr(i) -1)))) then
                           je = node_pntr(i)
                        endif
                     endif
                  enddo
               enddo
               
               if (je.ne.ie) then
                  ia = 0
                  ja = 0
                  ka = 0
                  ib = 0
                  jb = 0
                  kb = 0
                  
                  do k = 1,nn
                     do j = 1,nn
                        do i = 1,nn
                           js = nn*nn*(k -1) +nn*(j -1) +i
                           
                           if (con_spx(con_spx(je -1) +js).eq.an) then
                              ia = i
                              ja = j
                              ka = k
                           endif
                           
                           if (con_spx(con_spx(je -1) +js).eq.bn) then
                              ib = i
                              jb = j
                              kb = k
                           endif
                        enddo
                     enddo
                  enddo
                  
                  jsa = nn*nn*(ka -1) +nn*(ja -1) +(ia -1) +1
                  dljs = (nn*nn*(kb -ka) +nn*(jb -ja) +(ib -ia)) /(nn -1)
                  
                  do l = 1,nn
                     is = nn*nn*(l -1) +nn*(1 -1) +nn
                     js = jsa +dljs*(l -1)
                     con_spx(con_spx(ie -1) +is) = &
                          con_spx(con_spx(je -1) +js)
                  enddo
                  
               else
                  do k = 2,(nn -1)
                     nnode = nnode +1
                     is = nn*nn*(k -1) +nn*(1 -1) +nn
                     con_spx(con_spx(ie -1) + is) = nnode
                  enddo
               endif
               
               
! Third edge side
               
               an = con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +nn)
               bn = con_spx(con_spx(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +nn)
               
               je = ie
               do i = node_pntr(an -1),node_pntr(an) -1
                  do j = node_pntr(bn -1),node_pntr(bn) -1
                     if (node_pntr(i).eq.node_pntr(j)) then
                        if ((node_pntr(i).lt.ie).and.((nn*nn*nn +1)&
                        .eq.(con_spx(node_pntr(i)) - con_spx(node_pntr(i) -1)))) then
                           je = node_pntr(i)
                        endif
                     endif
                  enddo
               enddo
               
               if (je.ne.ie) then
                  ia = 0
                  ja = 0
                  ka = 0
                  ib = 0
                  jb = 0
                  kb = 0
                  
                  do k = 1,nn
                     do j = 1,nn
                        do i = 1,nn
                           js = nn*nn*(k -1) +nn*(j -1) +i
                           
                           if (con_spx(con_spx(je -1) +js).eq.an) then
                              ia = i
                              ja = j
                              ka = k
                           endif
                           
                           if (con_spx(con_spx(je -1) +js).eq.bn) then
                              ib = i
                              jb = j
                              kb = k
                           endif
                        enddo
                     enddo
                  enddo
                  
                  jsa = nn*nn*(ka -1) +nn*(ja -1) +(ia -1) +1
                  dljs = (nn*nn*(kb -ka) +nn*(jb -ja) +(ib -ia)) /(nn -1)
                  
                  do l = 1,nn
                     is = nn*nn*(l -1) +nn*(nn -1) +nn
                     js = jsa +dljs*(l -1)
                     con_spx(con_spx(ie -1) +is) = &
                          con_spx(con_spx(je -1) +js)
                  enddo
                  
               else
                  do k = 2,(nn -1)
                     nnode = nnode +1
                     is = nn*nn*(k -1) +nn*(nn -1) +nn
                     con_spx(con_spx(ie -1) + is) = nnode
                  enddo
               endif
               
               
! Fourth edge side
               
               an = con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +1)
               bn = con_spx(con_spx(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +1)
               
               je = ie
               do i = node_pntr(an -1),node_pntr(an) -1
                  do j = node_pntr(bn -1),node_pntr(bn) -1
                     if (node_pntr(i).eq.node_pntr(j)) then
                        if ((node_pntr(i).lt.ie).and.((nn*nn*nn +1)&
                        .eq.(con_spx(node_pntr(i)) - con_spx(node_pntr(i) -1)))) then
                           je = node_pntr(i)
                        endif
                     endif
                  enddo
               enddo
               
               if (je.ne.ie) then
                  ia = 0
                  ja = 0
                  ka = 0
                  
                  ib = 0
                  jb = 0
                  kb = 0
                  
                  do k = 1,nn
                     do j = 1,nn
                        do i = 1,nn
                           js = nn*nn*(k -1) +nn*(j -1) +i
                           
                           if (con_spx(con_spx(je -1) +js).eq.an) then
                              ia = i
                              ja = j
                              ka = k
                           endif
                           
                           if (con_spx(con_spx(je -1) +js).eq.bn) then
                              ib = i
                              jb = j
                              kb = k
                           endif
                        enddo
                     enddo
                  enddo
                  
                  jsa = nn*nn*(ka -1) +nn*(ja -1) +(ia -1) +1
                  dljs = (nn*nn*(kb -ka) +nn*(jb -ja) +(ib -ia)) /(nn -1)
                  
                  do l = 1,nn
                     is = nn*nn*(l -1) +nn*(nn -1) +1
                     js = jsa +dljs*(l -1)
                     con_spx(con_spx(ie -1) +is) = &
                          con_spx(con_spx(je -1) +js)
                  enddo
                  
               else
                  do k = 2,(nn -1)
                     nnode = nnode +1
                     is = nn*nn*(k -1) +nn*(nn -1) +1
                     con_spx(con_spx(ie -1) + is) = nnode
                  enddo
               endif
               
               
! First edge up
               
               an = con_spx(con_spx(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +1)
               bn = con_spx(con_spx(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +nn)
               
               je = ie             
               do i = node_pntr(an -1),node_pntr(an) -1
                  do j = node_pntr(bn -1),node_pntr(bn) -1
                     if (node_pntr(i).eq.node_pntr(j)) then
                        if ((node_pntr(i).lt.ie).and.((nn*nn*nn +1)&
                        .eq.(con_spx(node_pntr(i)) - con_spx(node_pntr(i) -1)))) then
                           je = node_pntr(i)
                        endif
                     endif
                  enddo
               enddo
               
               if (je.ne.ie) then
                  ia = 0
                  ja = 0
                  ka = 0
                  
                  ib = 0
                  jb = 0
                  kb = 0
                  
                  do k = 1,nn
                     do j = 1,nn
                        do i = 1,nn
                           js = nn*nn*(k -1) +nn*(j -1) +i                         
                           if (con_spx(con_spx(je -1) +js).eq.an) then
                              ia = i
                              ja = j
                              ka = k
                           endif
                                                                    
                           if (con_spx(con_spx(je -1) +js).eq.bn) then
                              ib = i
                              jb = j
                              kb = k
                           endif
                        enddo
                     enddo
                  enddo
                  
                  jsa = nn*nn*(ka -1) +nn*(ja -1) +(ia -1) +1
                  dljs = (nn*nn*(kb -ka) +nn*(jb -ja) +(ib -ia)) /(nn -1)
                  
                  do l = 1,nn
                     is = nn*nn*(nn -1) +nn*(1 -1) +l
                     js = jsa +dljs*(l -1)
                     con_spx(con_spx(ie -1) +is) = &
                          con_spx(con_spx(je -1) +js)
                  enddo
                  
               else
                  do i = 2,(nn -1)
                     nnode = nnode +1
                     is = nn*nn*(nn -1) +nn*(1 -1) +i
                     con_spx(con_spx(ie -1) + is) = nnode
                  enddo
               endif
               
               
! Second edge up
               
               an = con_spx(con_spx(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +nn)
               bn = con_spx(con_spx(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +nn)
               
               je = ie
               do i = node_pntr(an -1),node_pntr(an) -1
                  do j = node_pntr(bn -1),node_pntr(bn) -1
                     if (node_pntr(i).eq.node_pntr(j)) then
                        if ((node_pntr(i).lt.ie).and.((nn*nn*nn +1)&
                        .eq.(con_spx(node_pntr(i)) - con_spx(node_pntr(i) -1)))) then
                           je = node_pntr(i)
                        endif
                     endif
                  enddo
               enddo

               if (je.ne.ie) then
                  ia = 0
                  ja = 0
                  ka = 0
                  ib = 0
                  jb = 0
                  kb = 0
                  
                  do k = 1,nn
                     do j = 1,nn
                        do i = 1,nn
                           js = nn*nn*(k -1) +nn*(j -1) +i
                           
                           if (con_spx(con_spx(je -1) +js).eq.an) then
                              ia = i
                              ja = j
                              ka = k
                           endif
                           
                           if (con_spx(con_spx(je -1) +js).eq.bn) then
                              ib = i
                              jb = j
                              kb = k
                           endif
                        enddo
                     enddo
                  enddo
                  
                  jsa = nn*nn*(ka -1) +nn*(ja -1) +(ia -1) +1
                  dljs = (nn*nn*(kb -ka) +nn*(jb -ja) +(ib -ia)) /(nn -1)
                  
                  do l = 1,nn
                     is = nn*nn*(nn -1) +nn*(l -1) +nn
                     js = jsa +dljs*(l -1)
                     con_spx(con_spx(ie -1) +is) = &
                          con_spx(con_spx(je -1) +js)
                  enddo
                  
               else
                  do j = 2,(nn -1)
                     nnode = nnode +1
                     is = nn*nn*(nn -1) +nn*(j -1) +nn
                     con_spx(con_spx(ie -1) + is) = nnode
                  enddo
               endif
               
               
! Third edge up
               
               an = con_spx(con_spx(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +1)
               bn = con_spx(con_spx(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +nn)
               
               je = ie
               do i = node_pntr(an -1),node_pntr(an) -1
                  do j = node_pntr(bn -1),node_pntr(bn) -1
                     if (node_pntr(i).eq.node_pntr(j)) then
                        if ((node_pntr(i).lt.ie).and.((nn*nn*nn +1)&
                        .eq.(con_spx(node_pntr(i)) - con_spx(node_pntr(i) -1)))) then
                           je = node_pntr(i)
                        endif
                     endif
                  enddo
               enddo
               

               if (je.ne.ie) then
                  ia = 0
                  ja = 0
                  ka = 0
                  ib = 0
                  jb = 0
                  kb = 0
                  
                  do k = 1,nn
                     do j = 1,nn
                        do i = 1,nn
                           js = nn*nn*(k -1) +nn*(j -1) +i
                           
                           if (con_spx(con_spx(je -1) +js).eq.an) then
                              ia = i
                              ja = j
                              ka = k
                           endif
                           
                           if (con_spx(con_spx(je -1) +js).eq.bn) then
                              ib = i
                              jb = j
                              kb = k
                           endif
                        enddo
                     enddo
                  enddo
                  
                  jsa = nn*nn*(ka -1) +nn*(ja -1) +(ia -1) +1
                  dljs = (nn*nn*(kb -ka) +nn*(jb -ja) +(ib -ia)) /(nn -1)
                  
                  do l = 1,nn
                     is = nn*nn*(nn -1) +nn*(nn -1) +l
                     js = jsa +dljs*(l -1)
                     con_spx(con_spx(ie -1) +is) = &
                          con_spx(con_spx(je -1) +js)
                  enddo
                  
               else
                  do i = 2,(nn -1)
                     nnode = nnode +1
                     is = nn*nn*(nn -1) +nn*(nn -1) +i
                     con_spx(con_spx(ie -1) + is) = nnode
                  enddo
               endif
               
               
! Fourth edge up
               
               an = con_spx(con_spx(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +1)
               bn = con_spx(con_spx(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +1)
               
               je = ie
               do i = node_pntr(an -1),node_pntr(an) -1
                  do j = node_pntr(bn -1),node_pntr(bn) -1
                     if (node_pntr(i).eq.node_pntr(j)) then
                        if ((node_pntr(i).lt.ie).and.((nn*nn*nn +1)&
                        .eq.(con_spx(node_pntr(i)) - con_spx(node_pntr(i) -1)))) then
                           je = node_pntr(i)
                        endif
                     endif
                  enddo
               enddo
               

               if (je.ne.ie) then
                  ia = 0
                  ja = 0
                  ka = 0
                  ib = 0
                  jb = 0
                  kb = 0
                  
                  do k = 1,nn
                     do j = 1,nn
                        do i = 1,nn
                           js = nn*nn*(k -1) +nn*(j -1) +i
                           
                           if (con_spx(con_spx(je -1) +js).eq.an) then
                              ia = i
                              ja = j
                              ka = k
                           endif
                           
                           if (con_spx(con_spx(je -1) +js).eq.bn) then
                              ib = i
                              jb = j
                              kb = k
                           endif
                        enddo
                     enddo
                  enddo
                  
                  jsa = nn*nn*(ka -1) +nn*(ja -1) +(ia -1) +1
                  dljs = (nn*nn*(kb -ka) +nn*(jb -ja) +(ib -ia)) /(nn -1)
                  
                  do l = 1,nn
                     is = nn*nn*(nn -1) +nn*(l -1) +1
                     js = jsa +dljs*(l -1)
                     con_spx(con_spx(ie -1) +is) = &
                          con_spx(con_spx(je -1) +js)
                  enddo
                  
               else
                  do j = 2,(nn -1)
                     nnode = nnode +1
                     is = nn*nn*(nn -1) +nn*(j -1) +1
                     con_spx(con_spx(ie -1) + is) = nnode
                  enddo
               endif
               
               
! Then construct face connectivity
               
! Face down
               
               an = con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +1)
               bn = con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +nn)
               cn = con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +1)
               
               je = ie
               do i = node_pntr(an -1),node_pntr(an) -1
                  do j = node_pntr(bn -1),node_pntr(bn) -1
                     do k = node_pntr(cn -1),node_pntr(cn) -1
                        if ((node_pntr(i).eq.node_pntr(j))&
                             .and.(node_pntr(i).eq.node_pntr(k))) then
                           if ((node_pntr(i).lt.ie).and.((nn*nn*nn +1)&
                           .eq.(con_spx(node_pntr(i)) - con_spx(node_pntr(i) -1)))) then
                              je = node_pntr(i)
                           endif
                        endif
                     enddo
                  enddo
               enddo
               
               if (je.ne.ie) then
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
                           
                           if (con_spx(con_spx(je -1) +js).eq.an) then
                              ia = i
                              ja = j
                              ka = k
                           endif
                           
                           if (con_spx(con_spx(je -1) +js).eq.bn) then
                              ib = i
                              jb = j
                              kb = k
                           endif
                           
                           if (con_spx(con_spx(je -1) +js).eq.cn) then
                              ic = i
                              jc = j
                              kc = k
                           endif
                        enddo
                     enddo
                  enddo
                  
                  jsa = nn*nn*(ka -1) +nn*(ja -1) +(ia -1) +1
                  dljs = (nn*nn*(kb -ka) +nn*(jb -ja) +(ib -ia)) /(nn -1)
                  dmjs = (nn*nn*(kc -ka) +nn*(jc -ja) +(ic -ia)) /(nn -1)
                  
                  do m = 1,nn
                     do l = 1,nn
                        is = nn*nn*(1 -1) +nn*(m -1) +l
                        js = jsa +dmjs*(m -1) +dljs*(l -1)
                        con_spx(con_spx(ie -1) +is) = &
                             con_spx(con_spx(je -1) +js)
                     enddo
                  enddo
                  
               else
                  do m = 2,(nn -1)
                     do l = 2,(nn -1)
                        nnode = nnode +1
                        is = nn*nn*(1 -1) +nn*(m -1) +l
                        con_spx(con_spx(ie -1) + is) = nnode
                     enddo
                  enddo
               endif
               
               
! First face side
               
               an = con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +1)
               bn = con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +nn)
               cn = con_spx(con_spx(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +1)
               
               je = ie
               do i = node_pntr(an -1),node_pntr(an) -1
                  do j = node_pntr(bn -1),node_pntr(bn) -1
                     do k = node_pntr(cn -1),node_pntr(cn) -1
                        if ((node_pntr(i).eq.node_pntr(j))&
                             .and.(node_pntr(i).eq.node_pntr(k))) then
                           if ((node_pntr(i).lt.ie).and.((nn*nn*nn +1)&
                           .eq.(con_spx(node_pntr(i)) - con_spx(node_pntr(i) -1)))) then
                              je = node_pntr(i)
                           endif
                        endif
                     enddo
                  enddo
               enddo
               
               if (je.ne.ie) then
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
                           
                           if (con_spx(con_spx(je -1) +js).eq.an) then
                              ia = i
                              ja = j
                              ka = k
                           endif
                           
                           if (con_spx(con_spx(je -1) +js).eq.bn) then
                              ib = i
                              jb = j
                              kb = k
                           endif
                           
                           if (con_spx(con_spx(je -1) +js).eq.cn) then
                              ic = i
                              jc = j
                              kc = k
                           endif
                        enddo
                     enddo
                  enddo
                  
                  jsa = nn*nn*(ka -1) +nn*(ja -1) +(ia -1) +1
                  dljs = (nn*nn*(kb -ka) +nn*(jb -ja) +(ib -ia)) /(nn -1)
                  dmjs = (nn*nn*(kc -ka) +nn*(jc -ja) +(ic -ia)) /(nn -1)
                  
                  do m = 1,nn
                     do l = 1,nn
                        is = nn*nn*(m -1) +nn*(1 -1) +l
                        js = jsa +dmjs*(m -1) +dljs*(l -1)
                        con_spx(con_spx(ie -1) +is) = &
                             con_spx(con_spx(je -1) +js)
                     enddo
                  enddo
                  
               else
                  do m = 2,(nn -1)
                     do l = 2,(nn -1)
                        nnode = nnode +1
                        is = nn*nn*(m -1) +nn*(1 -1) +l
                        con_spx(con_spx(ie -1) + is) = nnode
                     enddo
                  enddo
               endif
               
               
! Second face side
               
               an = con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +nn)
               bn = con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +nn)
               cn = con_spx(con_spx(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +nn)
               
               je = ie
               do i = node_pntr(an -1),node_pntr(an) -1
                  do j = node_pntr(bn -1),node_pntr(bn) -1
                     do k = node_pntr(cn -1),node_pntr(cn) -1
                        if ((node_pntr(i).eq.node_pntr(j))&
                             .and.(node_pntr(i).eq.node_pntr(k))) then
                           if ((node_pntr(i).lt.ie).and.((nn*nn*nn +1)&
                           .eq.(con_spx(node_pntr(i)) - con_spx(node_pntr(i) -1)))) then
                              je = node_pntr(i)
                           endif
                        endif
                     enddo
                  enddo
               enddo
               
               if (je.ne.ie) then
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
                           
                           if (con_spx(con_spx(je -1) +js).eq.an) then
                              ia = i
                              ja = j
                              ka = k
                           endif
                           
                           if (con_spx(con_spx(je -1) +js).eq.bn) then
                              ib = i
                              jb = j
                              kb = k
                           endif
                           
                           if (con_spx(con_spx(je -1) +js).eq.cn) then
                              ic = i
                              jc = j
                              kc = k
                           endif
                        enddo
                     enddo
                  enddo
                  
                  jsa = nn*nn*(ka -1) +nn*(ja -1) +(ia -1) +1
                  dljs = (nn*nn*(kb -ka) +nn*(jb -ja) +(ib -ia)) /(nn -1)
                  dmjs = (nn*nn*(kc -ka) +nn*(jc -ja) +(ic -ia)) /(nn -1)
                  
                  do m = 1,nn
                     do l = 1,nn
                        is = nn*nn*(m -1) +nn*(l -1) +nn
                        js = jsa +dmjs*(m -1) +dljs*(l -1)
                        con_spx(con_spx(ie -1) +is) = &
                             con_spx(con_spx(je -1) +js)
                     enddo
                  enddo
                  
               else
                  do m = 2,(nn -1)
                     do l = 2,(nn -1)
                        nnode = nnode +1
                        is = nn*nn*(m -1) +nn*(l -1) +nn
                        con_spx(con_spx(ie -1) + is) = nnode
                     enddo
                  enddo
               endif
               
               
! Third face side
               
               an = con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +1)
               bn = con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +nn)
               cn = con_spx(con_spx(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +1)
               
               je = ie
               do i = node_pntr(an -1),node_pntr(an) -1
                  do j = node_pntr(bn -1),node_pntr(bn) -1
                     do k = node_pntr(cn -1),node_pntr(cn) -1
                        if ((node_pntr(i).eq.node_pntr(j))&
                             .and.(node_pntr(i).eq.node_pntr(k))) then
                           if ((node_pntr(i).lt.ie).and.((nn*nn*nn +1)&
                           .eq.(con_spx(node_pntr(i)) - con_spx(node_pntr(i) -1)))) then
                              je = node_pntr(i)
                           endif
                        endif
                     enddo
                  enddo
               enddo
               
               if (je.ne.ie) then
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
                           
                           if (con_spx(con_spx(je -1) +js).eq.an) then
                              ia = i
                              ja = j
                              ka = k
                           endif
                           
                           if (con_spx(con_spx(je -1) +js).eq.bn) then
                              ib = i
                              jb = j
                              kb = k
                           endif
                           
                           if (con_spx(con_spx(je -1) +js).eq.cn) then
                              ic = i
                              jc = j
                              kc = k
                           endif
                        enddo
                     enddo
                  enddo
                  
                  jsa = nn*nn*(ka -1) +nn*(ja -1) +(ia -1) +1
                  dljs = (nn*nn*(kb -ka) +nn*(jb -ja) +(ib -ia)) /(nn -1)
                  dmjs = (nn*nn*(kc -ka) +nn*(jc -ja) +(ic -ia)) /(nn -1)
                  
                  do m = 1,nn
                     do l = 1,nn
                        is = nn*nn*(m -1) +nn*(nn -1) +l
                        js = jsa +dmjs*(m -1) +dljs*(l -1)
                        con_spx(con_spx(ie -1) +is) = &
                             con_spx(con_spx(je -1) +js)
                     enddo
                  enddo
                  
               else
                  do m = 2,(nn -1)
                     do l = 2,(nn -1)
                        nnode = nnode +1
                        is = nn*nn*(m -1) +nn*(nn -1) +l
                        con_spx(con_spx(ie -1) + is) = nnode
                     enddo
                  enddo
               endif
               
               
! Fourth face side
               
               an = con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +1)
               bn = con_spx(con_spx(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +1)
               cn = con_spx(con_spx(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +1)
               
               je = ie
               do i = node_pntr(an -1),node_pntr(an) -1
                  do j = node_pntr(bn -1),node_pntr(bn) -1
                     do k = node_pntr(cn -1),node_pntr(cn) -1
                        if ((node_pntr(i).eq.node_pntr(j))&
                             .and.(node_pntr(i).eq.node_pntr(k))) then
                           if ((node_pntr(i).lt.ie).and.((nn*nn*nn +1)&
                           .eq.(con_spx(node_pntr(i)) - con_spx(node_pntr(i) -1)))) then
                              je = node_pntr(i)
                           endif
                        endif
                     enddo
                  enddo
               enddo
               
               if (je.ne.ie) then
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
                           
                           if (con_spx(con_spx(je -1) +js).eq.an) then
                              ia = i
                              ja = j
                              ka = k
                           endif
                           
                           if (con_spx(con_spx(je -1) +js).eq.bn) then
                              ib = i
                              jb = j
                              kb = k
                           endif
                           
                           if (con_spx(con_spx(je -1) +js).eq.cn) then
                              ic = i
                              jc = j
                              kc = k
                           endif
                        enddo
                     enddo
                  enddo
                  
                  jsa = nn*nn*(ka -1) +nn*(ja -1) +(ia -1) +1
                  dljs = (nn*nn*(kb -ka) +nn*(jb -ja) +(ib -ia)) /(nn -1)
                  dmjs = (nn*nn*(kc -ka) +nn*(jc -ja) +(ic -ia)) /(nn -1)
                  
                  do m = 1,nn
                     do l = 1,nn
                        is = nn*nn*(m -1) +nn*(l -1) +1
                        js = jsa +dmjs*(m -1) +dljs*(l -1)
                        con_spx(con_spx(ie -1) +is) = &
                             con_spx(con_spx(je -1) +js)
                     enddo
                  enddo
                  
               else
                  do m = 2,(nn -1)
                     do l = 2,(nn -1)
                        nnode = nnode +1
                        is = nn*nn*(m -1) +nn*(l -1) +1
                        con_spx(con_spx(ie -1) + is) = nnode
                     enddo
                  enddo
               endif
               
               
! Face up
               
               an = con_spx(con_spx(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +1)
               bn = con_spx(con_spx(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +nn)
               cn = con_spx(con_spx(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +1)
               
               je = ie
               do i = node_pntr(an -1),node_pntr(an) -1
                  do j = node_pntr(bn -1),node_pntr(bn) -1
                     do k = node_pntr(cn -1),node_pntr(cn) -1
                        if ((node_pntr(i).eq.node_pntr(j))&
                             .and.(node_pntr(i).eq.node_pntr(k))) then
                           if ((node_pntr(i).lt.ie).and.((nn*nn*nn +1)&
                           .eq.(con_spx(node_pntr(i)) - con_spx(node_pntr(i) -1)))) then
                              je = node_pntr(i)
                           endif
                        endif
                     enddo
                  enddo
               enddo
               
               if (je.ne.ie) then
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
                           
                           if (con_spx(con_spx(je -1) +js).eq.an) then
                              ia = i
                              ja = j
                              ka = k
                           endif
                           
                           if (con_spx(con_spx(je -1) +js).eq.bn) then
                              ib = i
                              jb = j
                              kb = k
                           endif
                           
                           if (con_spx(con_spx(je -1) +js).eq.cn) then
                              ic = i
                              jc = j
                              kc = k
                           endif
                        enddo
                     enddo
                  enddo
                  
                  jsa = nn*nn*(ka -1) +nn*(ja -1) +(ia -1) +1
                  dljs = (nn*nn*(kb -ka) +nn*(jb -ja) +(ib -ia)) /(nn -1)
                  dmjs = (nn*nn*(kc -ka) +nn*(jc -ja) +(ic -ia)) /(nn -1)
                  
                  do m = 1,nn
                     do l = 1,nn
                        is = nn*nn*(nn -1) +nn*(m -1) +l
                        js = jsa +dmjs*(m -1) +dljs*(l -1)
                        con_spx(con_spx(ie -1) +is) = &
                             con_spx(con_spx(je -1) +js)
                     enddo
                  enddo
                  
               else
                  do m = 2,(nn -1)
                     do l = 2,(nn -1)
                        nnode = nnode +1
                        is = nn*nn*(nn -1) +nn*(m -1) +l
                        con_spx(con_spx(ie -1) + is) = nnode
                     enddo
                  enddo
               endif
               
               
! End of edge and face connectivity
               
! Finally fill inside the elements
               
               do k = 2,(nn -1)
                  do j = 2,(nn -1)
                     do i = 2,(nn -1)
                        nnode = nnode +1
                        is = nn*nn*(k -1) +nn*(j -1) +i
                        con_spx(con_spx(ie -1) + is) = nnode
                     enddo
                  enddo
               enddo
               
            endif
         enddo
      enddo
      
      return
      
      end subroutine MAKE_SPX_CON

