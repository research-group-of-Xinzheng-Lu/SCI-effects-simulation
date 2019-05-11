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


!> @brief Makes invariants for the strain tensor and computes the main strain.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nn  polynomial  degree
!> @param[in] duxdx nodal values for spatial derivatives of the displacement 
!> @param[in] duydx nodal values for spatial derivatives of the displacement
!> @param[in] duzdx nodal values for spatial derivatives of the displacement
!> @param[in] duxdy nodal values for spatial derivatives of the displacement
!> @param[in] duydy nodal values for spatial derivatives of the displacement
!> @param[in] duzdy nodal values for spatial derivatives of the displacement
!> @param[in] duxdz nodal values for spatial derivatives of the displacement
!> @param[in] duydz nodal values for spatial derivatives of the displacement      
!> @param[in] duzdz nodal values for spatial derivatives of the displacement
!> @param[in] yon 1 if the point is non-linear elastic, 0 otherwise
!> @param[out] R averaged main shear strain


      subroutine MAKE_INVARIANTS_AND_MAIN_STRAIN(nn,&                                         
                                                duxdx,duydx,duzdx,&         
                                                duxdy,duydy,duzdy,&         
                                                duxdz,duydz,duzdz,&         
                                                R,yon)
                                            
      implicit none


      integer*4 :: nn
      integer*4 :: i,j,k

      integer*4, dimension(nn,nn,nn) :: yon        

      real*8 :: exx,eyy,ezz,exy,eyz,ezx
      real*8 :: R1,R2,R3                        
      real*8 :: invar1,invar2,invar3
      real*8 :: strain_I,strain_II,strain_III
      real*8 :: R_average_elem

      real*8, dimension(nn,nn,nn) :: duxdx,duxdy,duxdz,duydx,duydy,duydz,duzdx,duzdy,duzdz
      real*8, dimension(nn,nn,nn) :: R


      R_average_elem = 0.0d0
      
      do i = 1,nn
         do j = 1,nn
            do k = 1,nn

                if (yon(i,j,k).eq.1) then
                
                        !+-----------------------------------------
                        !! Invariants of the strain tensor sigma_ij

                        exx = duxdx(i,j,k)
                        eyy = duydy(i,j,k)
                        ezz = duzdz(i,j,k)

                        exy = .5 * (duxdy(i,j,k) + duydx(i,j,k))
                        eyz = .5 * (duydz(i,j,k) + duzdy(i,j,k))
                        ezx = .5 * (duzdx(i,j,k) + duxdz(i,j,k))

                        invar1 = exx + eyy + ezz

                        invar2 = exx * eyy + &
                                 eyy * ezz + &
                                 ezz* exx - &
                                 exy**2 - &
                                 eyz**2 - &
                                 ezx**2
                                                                  
                        invar3 = exx * eyy * ezz + &
                                 2 * exy * eyz * (-ezx) - &  
                                 exx * eyz**2 - &
                                 eyy * ezx**2 - &
                                 ezz * exy**2 

                        !+-----------------------------------------
                        !! Principal strain computations

                        call CUBIC(-invar1,invar2,-invar3,&
                                  strain_I,strain_II,strain_III)

                        !call ZerosPolyO3(-invar1,invar2,-invar3,&
                        !          strain_I,strain_II,strain_III)


                        !+-----------------------------------------
                        !! Max shear strain computation

                        R(i,j,k) = (max(strain_I,strain_II,strain_III) - &
                                        min(strain_I,strain_II,strain_III) ) / 2

                        R_average_elem = R_average_elem + R(i,j,k)
                                        
                        !! 
                        !+-----------------------------------------

                endif !if (yon(i,j,k).eq.1) then
                                          
            enddo
         enddo
      enddo

      R_average_elem = R_average_elem / (nn*nn*nn)

      do i = 1,nn
        do j = 1,nn
           do k = 1,nn

                R(i,j,k) = R_average_elem

           enddo
        enddo
      enddo

      return
      
      end subroutine MAKE_INVARIANTS_AND_MAIN_STRAIN

