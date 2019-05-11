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

!> @brief Generates extend seismic or explosive fault.
!> @note Faults can degenerate to a single point for special test case. 
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0

      subroutine MAKE_SEISMIC_MOMENT_OR_EXPLOSIVE_SOURCE()

      use max_var
      use speed_par
      

      implicit none
 
      include 'SPEED.MPI'      

!*****************************************************************************************
!                                    SEISMIC MOMENT
!*****************************************************************************************
                                                                      
      if (mpi_id.eq.0 .and. nload_sism_el.gt.0) & 
           write(*,'(A)') '-----------Building the Seismic Moment vector----------'


!     Dimensioning vector 'num_node_sism'(nodes number along each faults) - begin
    
      if (nload_sism_el.gt.0) then
       
         allocate (num_node_sism(nload_sism_el))
         num_node_sism = 0
      
         do i = 1,nload_sism_el
            if ((((val_sism_el(i,1).eq.val_sism_el(i,4)).and.(val_sism_el(i,4).eq.val_sism_el(i,7))) .and. &
                (val_sism_el(i,7).eq.val_sism_el(i,10))).and.(((val_sism_el(i,2).eq.val_sism_el(i,5)).and. & 
                (val_sism_el(i,5).eq.val_sism_el(i,8))) .and. (val_sism_el(i,8).eq.val_sism_el(i,11))) &
               .and.(((val_sism_el(i,3).eq.val_sism_el(i,6)).and.(val_sism_el(i,6).eq.val_sism_el(i,9))) & 
               .and. (val_sism_el(i,9).eq.val_sism_el(i,12)))) then
               
               num_node_sism(i) = 1

            else
         
               call GET_DIME_SISM(val_sism_el(i,1),val_sism_el(i,2),val_sism_el(i,3),&
                                    val_sism_el(i,4),val_sism_el(i,5),val_sism_el(i,6),&
                                    val_sism_el(i,7),val_sism_el(i,8),val_sism_el(i,9),&
                                    val_sism_el(i,10),val_sism_el(i,11),val_sism_el(i,12),&
                                    nnode_dom,xx_spx_loc,yy_spx_loc,zz_spx_loc,&  ! nnod_loc
                                    num_node_sism(i),nnod_loc,mpi_id,i,local_node_num)


            endif
        enddo !i = 1,nload_sism_el
            
       allocate(sism_el_glo(nload_sism_el*mpi_np), vec(mpi_np))
       call MPI_BARRIER(mpi_comm, mpi_ierr)                  
       call MPI_ALLGATHER(num_node_sism, nload_sism_el, SPEED_INTEGER, sism_el_glo, &
                          nload_sism_el, SPEED_INTEGER, mpi_comm, mpi_ierr)


       if(nload_sism_el .eq. 1) then 
         num_node_sism(1) = 1
       else  
         do i = 1, nload_sism_el
            do j = 1, mpi_np
               vec(j) = sism_el_glo(i + (j-1)*nload_sism_el)          
            enddo
            num_node_sism(i) = sum(vec) 
         enddo
       endif
            
       deallocate(vec, sism_el_glo)

      !Checking the maximum number of nodes for single seismic faults

         max_num_node_sism = num_node_sism(1)

         do i = 1,nload_sism_el
            if (num_node_sism(i).gt.max_num_node_sism)  max_num_node_sism = num_node_sism(i)
         enddo
     
      
         allocate (sour_node_sism(max_num_node_sism,nload_sism_el))
         allocate (dist_sour_node_sism(max_num_node_sism,nload_sism_el))
         sour_node_sism = 0
         dist_sour_node_sism = 0.d0
      
         !Searching the node 'id' in the global numeration for each seismic faults.
         !sour_node_sism = node id (global numeration) generating the 'i'th seismic fault
         do i = 1,nload_sism_el
            if (num_node_sism(i) .eq. 1 .and. nload_sism_el .eq. 1) then
               
               
                call GET_NEAREST_NODE(nnod_loc,xx_spx_loc,yy_spx_loc,zz_spx_loc, &
                                       val_sism_el(i,1),val_sism_el(i,2),val_sism_el(i,3), &
                                       sour_node_sism(1,i),dist_sour_node_sism(1,i))
                                          
                 
                !global numeration 
                sour_node_sism(1,i) = local_node_num(sour_node_sism(1,i))

                allocate(sism_el_glo(mpi_np), dist_el_glo(mpi_np))
                call MPI_BARRIER(mpi_comm, mpi_ierr)                  
                call MPI_ALLGATHER(sour_node_sism(1,i), 1, SPEED_INTEGER, sism_el_glo, 1, SPEED_INTEGER, mpi_comm, mpi_ierr)
                call MPI_ALLGATHER(dist_sour_node_sism(1,i), 1, SPEED_DOUBLE, dist_el_glo, 1, SPEED_DOUBLE, mpi_comm, mpi_ierr)
                

                               
                call GET_MINVALUES(sism_el_glo, dist_el_glo, mpi_np, sour_node_sism(1,i), 1, mpi_np)               
                dist_sour_node_sism(1,i) = 0 
                deallocate(sism_el_glo, dist_el_glo) 
                                   

                                   

            else   !if(num_node_sism(i) .gt. 1) then

                 call READ_SISM(val_sism_el(i,1),val_sism_el(i,2),val_sism_el(i,3),&
                                      val_sism_el(i,4),val_sism_el(i,5),val_sism_el(i,6),&
                                      val_sism_el(i,7),val_sism_el(i,8),val_sism_el(i,9),&
                                      val_sism_el(i,10),val_sism_el(i,11),val_sism_el(i,12),&
                                      nnode_dom,xx_spx_loc,yy_spx_loc,zz_spx_loc,&    
                                      num_node_sism(i),sour_node_sism,i,&
                                      dist_sour_node_sism,nload_sism_el,&
                                      max_num_node_sism,local_node_num,nnod_loc)

                 allocate(sism_el_glo(max_num_node_sism*mpi_np), dist_el_glo(max_num_node_sism*mpi_np))
                 sism_el_glo = 0
                 dist_el_glo = 0.d0
                 
                 call MPI_BARRIER(mpi_comm, mpi_ierr)                  
                 call MPI_ALLGATHER(sour_node_sism(:,i), max_num_node_sism, SPEED_INTEGER, &
                                    sism_el_glo, max_num_node_sism, SPEED_INTEGER, mpi_comm, mpi_ierr)
                                    
                 call MPI_ALLGATHER(dist_sour_node_sism(:,i), max_num_node_sism, SPEED_DOUBLE, &
                                    dist_el_glo, max_num_node_sism, SPEED_DOUBLE, mpi_comm, mpi_ierr)



                 j = 0                       
                 do k = 1, mpi_np*max_num_node_sism
                       if(sism_el_glo(k).ne. 0) then 
                       
                          j = j + 1 
                          sour_node_sism(j,i) = sism_el_glo(k)
                          dist_sour_node_sism(j,i) = dist_el_glo(k)

                        endif
                     enddo   
                  deallocate(sism_el_glo, dist_el_glo)
             endif

             if (mpi_id.eq.0 .and. i .eq. 1) write(*,'(A)')
             if (mpi_id.eq.0 .and. i .eq. 1) write(*,'(A)')'Seismic faults'
             if (mpi_id.eq.0 .and. i .eq. 1) write(*,'(A)')
             if (mpi_id.eq.0) write(*,'(A,I6,A,I6,A)')'Seismic Faults - ',i, ' is generated by ',num_node_sism(i),' nodes'
             if (mpi_id.eq.0) write(*,'(I10,I12)')(j,sour_node_sism(j,i), j=1,num_node_sism(i))                               
             if (mpi_id.eq.0) write(*,'(A)')
         enddo !i = 1,nload_sism_el


      endif !if (nload_sism_el.gt.0) then
      
      ! Dimensioning vector 'num_node_sism'(nodes number along each faults) - end  
          
    
      if (nfunc.le.0) nfunc = 1
      
      if (nload_sism_el.gt.0) then                         
         allocate (factor_seismic_moment(nload_sism_el,6)) 
         allocate (tau_seismic_moment(nload_sism_el,1)) 
      endif                                                 

       if ((mpi_id.eq.0).and.(nload_sism_el.gt.0)) write(*,'(A)')'Seismic Moment vector built'


!*****************************************************************************************
!                                EXPLOSIVE SOURCE
!*****************************************************************************************

       if (mpi_id.eq.0 .and. nload_expl_el.gt.0 ) write(*,'(A)')'----------Building the Explosive Source vector---------'

      ! Dimensioning vector 'num_node_expl'(nodes number along each explosion) - begin
      !
      if (nload_expl_el.gt.0) then
       
                allocate (num_node_expl(nload_expl_el))
      
         do i = 1,nload_expl_el
            if ((((val_expl_el(i,1).eq.val_expl_el(i,4)).and.(val_expl_el(i,4).eq.val_expl_el(i,7))) &
               .and. (val_expl_el(i,7).eq.val_expl_el(i,10))).and.(((val_expl_el(i,2).eq.val_expl_el(i,5)) & 
               .and.(val_expl_el(i,5).eq.val_expl_el(i,8))) .and. (val_expl_el(i,8).eq.val_expl_el(i,11))) &
               .and.(((val_expl_el(i,3).eq.val_expl_el(i,6)).and.(val_expl_el(i,6).eq.val_expl_el(i,9))) &
               .and. (val_expl_el(i,9).eq.val_expl_el(i,12)))) then
               
               num_node_expl(i)=1
            
            else  
               call GET_DIME_EXPL(val_expl_el(i,1),val_expl_el(i,2),val_expl_el(i,3),&
                                    val_expl_el(i,4),val_expl_el(i,5),val_expl_el(i,6),&
                                    val_expl_el(i,7),val_expl_el(i,8),val_expl_el(i,9),&
                                    val_expl_el(i,10),val_expl_el(i,11),val_expl_el(i,12),&
                                    nnode_dom,xx_spx_loc,yy_spx_loc,zz_spx_loc,& 
                                    num_node_expl(i),nnod_loc)
                endif
         enddo !i = 1,nload_expl_el


         allocate(expl_el_glo(nload_expl_el*mpi_np), vec(mpi_np))
         call MPI_BARRIER(mpi_comm, mpi_ierr)                  
         call MPI_ALLGATHER(num_node_expl, nload_expl_el, SPEED_INTEGER, expl_el_glo, &
                            nload_expl_el, SPEED_INTEGER, mpi_comm, mpi_ierr)
 
       if(nload_expl_el .eq. 1) then 
          num_node_expl(1) = 1
       else   
          do i = 1, nload_expl_el
             do j = 1, mpi_np
                vec(j) = expl_el_glo(i + (j-1)*nload_expl_el)          
             enddo
           
             num_node_expl(i) = sum(vec) 
          enddo
       endif
 
       
         deallocate(vec, expl_el_glo)
       
         !Checking the maximum number of nodes for single seismic faults

         max_num_node_expl = num_node_expl(1)

         do i = 1,nload_expl_el
                if (num_node_expl(i).gt.max_num_node_expl) then
                        max_num_node_expl = num_node_expl(i)
                endif
         enddo
      
         allocate (sour_node_expl(max_num_node_expl,nload_expl_el))
         allocate (dist_sour_node_expl(max_num_node_expl,nload_expl_el))
         sour_node_expl = 0
         dist_sour_node_expl = 0.d0
      
         !Searching the node 'id' in the global numeration for each Explosive Source.
         !sour_node_expl = node id (global numeration) generating the 'i'th Explosion
         do i = 1,nload_expl_el
            if (num_node_expl(i).eq.1 .and.  nload_expl_el .eq. 1 ) then
               
                call GET_NEAREST_NODE(nnod_loc,xx_spx_loc,yy_spx_loc,zz_spx_loc, &   
                                       val_expl_el(i,1),val_expl_el(i,2),val_expl_el(i,3), &
                                       sour_node_expl(1,i),dist_sour_node_expl(1,i))
                                       
                                       
                sour_node_expl(1,i) = local_node_num(sour_node_expl(1,i))                       

                allocate(expl_el_glo(mpi_np), dist_el_glo(mpi_np))
                call MPI_BARRIER(mpi_comm, mpi_ierr)                  
                call MPI_ALLGATHER(sour_node_expl(1,i), 1, SPEED_INTEGER, expl_el_glo, 1, SPEED_INTEGER, mpi_comm, mpi_ierr)
                call MPI_ALLGATHER(dist_sour_node_expl(1,i), 1, SPEED_DOUBLE, dist_el_glo, 1, SPEED_DOUBLE, mpi_comm, mpi_ierr)
                
                call GET_MINVALUES(expl_el_glo, dist_el_glo, mpi_np, sour_node_expl(1,i), 1, mpi_np)                            
                dist_sour_node_expl(1,i) = 0
                
                deallocate(expl_el_glo, dist_el_glo) 
                                   
                
            else
                call READ_EXPL(val_expl_el(i,1),val_expl_el(i,2),val_expl_el(i,3),&
                                     val_expl_el(i,4),val_expl_el(i,5),val_expl_el(i,6),&
                                     val_expl_el(i,7),val_expl_el(i,8),val_expl_el(i,9),&
                                     val_expl_el(i,10),val_expl_el(i,11),val_expl_el(i,12),&
                                     nnode_dom,xx_spx_loc,yy_spx_loc,zz_spx_loc,&   !nnod_loc
                                     num_node_expl(i),sour_node_expl,i,&
                                     dist_sour_node_expl,nload_expl_el,&
                                     max_num_node_expl,local_node_num,nnod_loc)


                allocate(expl_el_glo(max_num_node_expl*mpi_np), dist_el_glo(max_num_node_expl*mpi_np))
                call MPI_BARRIER(mpi_comm, mpi_ierr)                  
                call MPI_ALLGATHER(sour_node_expl(:,i), max_num_node_expl, SPEED_INTEGER, &
                                   expl_el_glo, max_num_node_expl, SPEED_INTEGER, mpi_comm, mpi_ierr)
                call MPI_ALLGATHER(dist_sour_node_expl(:,i), max_num_node_expl, SPEED_DOUBLE, &
                                   dist_el_glo, max_num_node_expl, SPEED_DOUBLE, mpi_comm, mpi_ierr)
                 j = 0                       
                 do k = 1, mpi_np*max_num_node_expl
                       if(expl_el_glo(k).ne. 0) then 
                       
                          j = j + 1 
                          sour_node_expl(j,i) = expl_el_glo(k)
                          dist_sour_node_expl(j,i) = dist_el_glo(k)

                        endif
                     enddo   
                  deallocate(expl_el_glo, dist_el_glo)
                  
             endif
             if (mpi_id.eq.0 .and. i .eq. 1) write(*,'(A)')
             if (mpi_id.eq.0 .and. i .eq. 1) write(*,'(A)')'Explosive Source'
             if (mpi_id.eq.0 .and. i .eq. 1) write(*,'(A)')
             if (mpi_id.eq.0) write(*,'(A,I6,A,I6,A)')'Explosive Source - ',i, ' is generated by ',num_node_expl(i),' nodes'
             if (mpi_id.eq.0) write(*,'(I6,I6)')(j,sour_node_expl(j,i), j=1,num_node_expl(i))
             if (mpi_id.eq.0 .and. i .eq. 1) write(*,'(A)')
         enddo !i = 1,nload_expl_el


      endif !if (nload_expl_el.gt.0) then
      
      ! Dimensioning vector 'num_node_expl'(nodes number along each explosion area) - end  
       if (nfunc.le.0) nfunc = 1
     
      if (nload_expl_el.gt.0) allocate (factor_explosive_source(nload_expl_el,6))
      if ((mpi_id.eq.0).and.(nload_expl_el.gt.0)) write(*,'(A)')'Explosive Source vector built'



      end subroutine MAKE_SEISMIC_MOMENT_OR_EXPLOSIVE_SOURCE

