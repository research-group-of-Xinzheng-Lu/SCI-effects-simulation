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

!> @brief Reads and stores data from filemate (*.mate)
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] filemate  file name (*.mate)
!> @param[in] nb_mate  number of materials
!> @param[in] nb_mate_nle  number of non-linear materials
!> @param[in] nb_diriX  number of Dirichlet boundary conditions (x-dir)
!> @param[in] nb_diriY  number of Dirichlet boundary conditions (y-dir)
!> @param[in] nb_diriZ  number of Dirichlet boundary conditions (z-dir)
!> @param[in] nb_neuX  number of Neumann boundary conditions (x-dir)
!> @param[in] nb_neuY  number of Neumann boundary conditions (y-dir)
!> @param[in] nb_neuZ  number of Neumann boundary conditions (z-dir)
!> @param[in] nb_neuN  number of Neumann boundary conditions (N-dir)
!> @param[in] nb_poiX  number of point load volume force (x-dir) 
!> @param[in] nb_poiY  number of point load volume force (y-dir) 
!> @param[in] nb_poiZ  number of point load volume force (z-dir) 
!> @param[in] ntX  number of travelling point load (x-dir)
!> @param[in] ntY  number of travelling point load (y-dir)
!> @param[in] ntZ  number of travelling point load (z-dir)
!> @param[in] nb_plaX  number of plane wave loads (x-dir)
!> @param[in] nb_plaY  number of plane wave loads (y-dir)
!> @param[in] nb_plaZ  number of plane wave loads (z-dir)
!> @param[in] nb_forX  number of volume forces (x-dir)
!> @param[in] nb_forY  number of volume forces (y-dir)
!> @param[in] nb_forZ  number of volume forces (z-dir)
!> @param[in] nb_forF  number of external forces
!> @param[in] nb_pre  number of pressure loads
!> @param[in] nb_she  number of shear loads
!> @param[in] ntest  test case mode
!> @param[in] nb_abc  number of abc boundary conditions
!> @param[in] nb_dg  number od DG interface conditions
!> @param[in] nb_sism  number of sismic loads
!> @param[in] nb_expl  number of explosive loads
!> @param[in] nb_case  number of not honoring cases
!> @param[in] nb_fnc  number of functions
!> @param[out] char_mat  material properties (rho,lambda,mu,gamma)
!> @param[out] type_mat spectral degree use within the material
!> @param[out] trefm obsolete
!> @param[out] lab_mat label for the material
!> @param[out] char_mat_nle material properties (rho,lambda,mu,gamma)
!> @param[out] val_mat_nle spectral degree use within the material
!> @param[out] type_mat_nle obsolete
!> @param[out] lab_mat_nle label for the material
!> @param[out] val_diriX  val_ue of the displacement the Dirichlet boundary (x-dir)
!> @param[out] fnc_diriX  number of the time function for the Dir. b.c. (x-dir)
!> @param[out] lab_diriX  label of the Dir. b.c. (x-dir)
!> @param[out] val_diriY  val_ue of the displacement the Dirichlet boundary (y-dir)
!> @param[out] fnc_diriY  number of the time function for the Dir. b.c. (y-dir)
!> @param[out] lab_diriY  label of the Dir. b.c. (y-dir)
!> @param[out] val_diriZ  val_ue of the displacement the Dirichlet boundary (z-dir)
!> @param[out] fnc_diriZ  number of the time function for the Dir. b.c. (z-dir)
!> @param[out] lab_diriZ  label of the Dir. b.c. (z-dir)
!> @param[out] val_neuX amplitude of the Neumann load (x-dir)
!> @param[out] fnc_neuX  number of the time function for the Neu. b.c. (x-dir)
!> @param[out] lab_neuX  label of the Neu. b.c. (x-dir)
!> @param[out] val_neuY amplitude of the Neumann load (y-dir)
!> @param[out] fnc_neuY  number of the time function for the Neu. b.c. (y-dir)
!> @param[out] lab_neuY  label of the Neu. b.c. (y-dir)
!> @param[out] val_neuZ amplitude of the Neumann load (z-dir)
!> @param[out] fnc_neuZ  number of the time function for the Neu. b.c. (z-dir)
!> @param[out] lab_neuZ  label of the Neu. b.c. (z-dir)
!> @param[out] val_neuN amplitude of the Neumann load (N-dir)
!> @param[out] fnc_neuN  number of the time function for the Neu. b.c. (N-dir)
!> @param[out] lab_neuN  label of the Neu. b.c. (N-dir)
!> @param[out] val_poiX  amplitude of point load volume force (x-dir)
!> @param[out] fnc_poiX  number of the time function for point load volume force (x-dir) 
!> @param[out] val_poiY  amplitude of point load volume force (y-dir)
!> @param[out] fnc_poiY  number of the time function for point load volume force (y-dir) 
!> @param[out] val_poiZ  amplitude of point load volume force (z-dir)
!> @param[out] fnc_poiZ  number of the time function for point load volume force (z-dir) 
!> @param[out] valtX     to be completed
!> @param[out] ftX       to be completed
!> @param[out] valtY     to be completed
!> @param[out] ftY       to be completed
!> @param[out] valtZ     to be completed
!> @param[out] ftz       to be completed
!> @param[out] val_plaX  amplitude for plane wave loads (x-dir)
!> @param[out] fnc_plaX  number of the time function for plane wave loads (x-dir)
!> @param[out] lab_plaX  label of the plane wave load (x-dir)
!> @param[out] val_plaY  amplitude for plane wave loads (y-dir)
!> @param[out] fnc_plaY  number of the time function for plane wave loads (y-dir)
!> @param[out] lab_plaY  label of the plane wave load (y-dir)
!> @param[out] val_plaZ  amplitude for plane wave loads (z-dir)
!> @param[out] fnc_plaZ  number of the time function for plane wave loads (z-dir)
!> @param[out] lab_plaZ  label of the plane wave load (z-dir)
!> @param[out] val_forX  amplitude of volume forces (x-dir)
!> @param[out] fnc_forX  number of the time function for volume forces (x-dir)
!> @param[out] val_forY  amplitude of volume forces (y-dir)
!> @param[out] fnc_forY  number of the time function for volume forces (y-dir)
!> @param[out] val_forZ  amplitude of volume forces (z-dir)
!> @param[out] fnc_forZ  number of the time function for volume forces (z-dir)
!> @param[out] val_forF  val_ues for external forces
!> @param[out] fnc_forF  time function for external forces
!> @param[out] val_pre  amplitude of pressure loads
!> @param[out] fnc_pre  number of the time function for pressure loads
!> @param[out] val_she  amplitude of shear loads 
!> @param[out] fnc_she  number of the time function for shear loads
!> @param[out] ftest  function for test mode 
!> @param[out] lab_abc  label for abc boundary conditions
!> @param[out] lab_dg  label for DG interface conditions
!> @param[out] val_sism  val_ues for seismic loads
!> @param[out] fnc_sism  time function for seismic loads
!> @param[out] lab_sism  label for seismic loads
!> @param[out] val_expl  val_ues for explosive loads
!> @param[out] fnc_expl  time function for explosive loads
!> @param[out] lab_expl  label for explosive loads
!> @param[out] val_case  val_ue of not honoring case
!> @param[out] lab_case  label for not honoring cases
!> @param[out] tol_case  tolerance for not honoring routine
!> @param[out] type_fnc  function type 
!> @param[out] ind_fnc  pointers for data functions
!> @param[out] dat_fnc  data functions
!> @param[out] lab_fnc  label for time functions
!> @param[out] fmax  max frequency of the plane wave load 
!> @param[out] fpeak peak frequency for non lienar damping


      subroutine READ_FILEMATE(filemate,nb_mate,char_mat,type_mat,trefm,lab_mat,Qua_S,Qua_P, &
                                  nb_mate_nle,char_mat_nle,val_mat_nle,type_mat_nle,lab_mat_nle, & 
                                  nb_diriX,val_diriX,fnc_diriX,lab_diriX, &
                                  nb_diriY,val_diriY,fnc_diriY,lab_diriY, &
                                  nb_diriZ,val_diriZ,fnc_diriZ,lab_diriZ, &
                                  nb_neuX,val_neuX,fnc_neuX,lab_neuX, &
                                  nb_neuY,val_neuY,fnc_neuY,lab_neuY, &
                                  nb_neuZ,val_neuZ,fnc_neuZ,lab_neuZ, &
                                  nb_neuN,val_neuN,fnc_neuN,lab_neuN, & 
                                  nb_poiX,val_poiX,fnc_poiX, &
                                  nb_poiY,val_poiY,fnc_poiY, &
                                  nb_poiZ,val_poiZ,fnc_poiZ, &
                                  ntX,valtX,ftX, &
                                  ntY,valtY,ftY, &
                                  ntZ,valtZ,ftZ, &
                                  nb_plaX,val_plaX,fnc_plaX,lab_plaX, & 
                                  nb_plaY,val_plaY,fnc_plaY,lab_plaY, & 
                                  nb_plaZ,val_plaZ,fnc_plaZ,lab_plaZ, & 
                                  nb_forX,val_forX,fnc_forX, &
                                  nb_forY,val_forY,fnc_forY, &
                                  nb_forZ,val_forZ,fnc_forZ, &
                                  nb_forF,val_forF,fnc_forF, &
                                  nb_pre,val_pre,fnc_pre, &
                                  nb_she,val_she,fnc_she, &
                                  ntest,ftest, & !valftest,&
                                  nb_abc,lab_abc, &
                                  nb_dg,lab_dg,lab_dg_yn, &
                                  nb_sism,val_sism,fnc_sism,lab_sism, & 
                                  nb_expl,val_expl,fnc_expl,lab_expl, & 
                                  nb_case,val_case,lab_case,tol_case, & 
                                  nb_fnc,type_fnc,ind_fnc,dat_fnc,lab_fnc,nb_fnc_data, &
                                  fmax,fpeak)
                                  

      use speed_exit_codes
!      use speed_par, only: slip_type
      
      implicit none
      
      character*70 :: filemate,fileinput
      character*100000 :: inline
      character*4 :: keyword

      integer*4 :: nb_mate,nb_fnc, nb_fnc_data
      integer*4 :: nb_mate_nle        
      integer*4 :: nb_diriX,nb_diriY,nb_diriZ,nb_neuX,nb_neuY,nb_neuZ
      integer*4 :: nb_neuN 
      integer*4 :: nb_poiX,nb_poiY,nb_poiZ,nb_forX,nb_forY,nb_forZ,nb_forF
      integer*4 :: nb_plaX,nb_plaY,nb_plaZ,ntX,ntY,ntZ 
      integer*4 :: nb_abc, nb_dg
      integer*4 :: nb_pre,nb_she
      integer*4 :: isism,nb_sism                                        
      integer*4 :: iexpl,nb_expl                                        
      integer*4 :: icase,nb_case                                        
      
      integer*4 :: im,ifunc,idf,ndat_fnc,file_nd
      integer*4 :: im_nle        
      integer*4 :: idX,idY,idZ,inX,inY,inZ, itX, itY, itZ
      integer*4 :: inN, itest 
      integer*4 :: ipX,ipY,ipZ,ifX,ifY,ifZ,iff
      integer*4 :: iplX,iplY,iplZ 
      integer*4 :: ipr,ish,iabc,idg
      integer*4 :: ileft,iright
      integer*4 :: i,j,dummy,status
      integer*4 :: val_case, lab_case, ntest

      integer*4, dimension(nb_diriX) :: fnc_diriX,lab_diriX
      integer*4, dimension(nb_diriY) :: fnc_diriY,lab_diriY
      integer*4, dimension(nb_diriZ) :: fnc_diriZ,lab_diriZ
      integer*4, dimension(nb_neuX) :: fnc_neuX,lab_neuX
      integer*4, dimension(nb_neuY) :: fnc_neuY,lab_neuY
      integer*4, dimension(nb_neuZ) :: fnc_neuZ,lab_neuZ
      integer*4, dimension(nb_neuN) :: fnc_neuN,lab_neuN
      integer*4, dimension(ntest) :: ftest              
      integer*4, dimension(nb_forX) :: fnc_forX
      integer*4, dimension(nb_forY) :: fnc_forY
      integer*4, dimension(nb_forZ) :: fnc_forZ
      integer*4, dimension(nb_forF) :: fnc_forF
      integer*4, dimension(nb_pre) :: fnc_pre
      integer*4, dimension(nb_she) :: fnc_she     
      integer*4, dimension(nb_abc) :: lab_abc
      integer*4, dimension(nb_dg) :: lab_dg, lab_dg_yn
      integer*4, dimension(nb_sism) :: fnc_sism,lab_sism                  
      integer*4, dimension(nb_expl) :: fnc_expl,lab_expl                  
      !integer*4, dimension(*) :: lab_case                        
      integer*4, dimension(nb_poiX) :: fnc_poiX
      integer*4, dimension(nb_poiY) :: fnc_poiY
      integer*4, dimension(nb_poiZ) :: fnc_poiZ
      integer*4, dimension(ntX) :: ftX
      integer*4, dimension(ntY) :: ftY
      integer*4, dimension(ntZ) :: ftZ
      integer*4, dimension(nb_plaX) :: fnc_plaX,lab_plaX 
      integer*4, dimension(nb_plaY) :: fnc_plaY,lab_plaY 
      integer*4, dimension(nb_plaZ) :: fnc_plaZ,lab_plaZ 
      integer*4, dimension(nb_mate) :: type_mat
      integer*4, dimension(nb_mate) :: lab_mat
      integer*4, dimension(nb_fnc) :: type_fnc
      integer*4, dimension(nb_fnc) :: lab_fnc
      integer*4, dimension(nb_fnc +1) :: ind_fnc
      integer*4, dimension(nb_mate_nle) :: lab_mat_nle                
      
      integer*4, dimension(nb_mate_nle) :: type_mat_nle                
      integer*4, dimension(nb_mate_nle,1) :: char_mat_nle        

      real*8 :: fmax                                                                
      real*8 :: fpeak                                                                
      real*8 :: tol_case, rho, VS, VP
      
      real*8, dimension(nb_fnc_data) :: dat_fnc
      real*8, dimension(nb_mate) :: trefm, Qua_S, Qua_P
      
      real*8, dimension(nb_diriX,4) :: val_diriX
      real*8, dimension(nb_diriY,4) :: val_diriY
      real*8, dimension(nb_diriZ,4) :: val_diriZ
      !!!!!!!!!!!!!!!!!!modified by ty 170410
      real*8, dimension(nb_neuX,6) :: val_neuX
      real*8, dimension(nb_neuY,6) :: val_neuY
      real*8, dimension(nb_neuZ,6) :: val_neuZ
      real*8, dimension(nb_neuN,6) :: val_neuN                      
      real*8, dimension(nb_poiX,6) :: val_poiX
      real*8, dimension(nb_poiY,6) :: val_poiY
      real*8, dimension(nb_poiZ,6) :: val_poiZ
      !!!!!!!!!!!!!!!!!!modified by ty 170410
      real*8, dimension(ntX,4) :: valtX
      real*8, dimension(ntY,4) :: valtY
      real*8, dimension(ntZ,4) :: valtZ      
      real*8, dimension(nb_plaX,1) :: val_plaX 
      real*8, dimension(nb_plaY,1) :: val_plaY 
      real*8, dimension(nb_plaZ,1) :: val_plaZ 
      real*8, dimension(nb_forX,4) :: val_forX
      real*8, dimension(nb_forY,4) :: val_forY
      real*8, dimension(nb_forZ,4) :: val_forZ
      real*8, dimension(nb_forF,10) :: val_forF
      real*8, dimension(nb_pre,10) :: val_pre
      real*8, dimension(nb_she,10) :: val_she
      real*8, dimension(nb_sism,21) :: val_sism                        
      real*8, dimension(nb_expl,20) :: val_expl                        
      real*8, dimension(nb_mate_nle,1) :: val_mat_nle                
      real*8, dimension(nb_mate,4) :: char_mat


      open(40,file=filemate)
      

      
      im = 0;       im_nle = 0;          icase = 0
      idX = 0;      idY = 0;      idZ = 0
      inX = 0;      inY = 0;      inZ = 0;      inN = 0;
      ipX = 0;      ipY = 0;      ipZ = 0
      iplX = 0;     iplY = 0;     iplZ = 0 
      ifX = 0;      ifY = 0;      ifZ = 0
      iff = 0;      ipr = 0;      ish = 0
      iabc = 0;     idg = 0
      isism = 0;    iexpl = 0
      itest = 0 
      
      ifunc = 0
      fmax = 0.d0
  
      
      
      if (nb_fnc.gt.0) ind_fnc(1) = 1
      
      
      do 
         read(40,'(A)',IOSTAT = status) inline
         
         if (status.ne.0) exit
         
         keyword = inline(1:4)
         ileft = 0
         iright = len(inline)
         do i = 1,iright
            if (inline(i:i).eq.' ') exit
         enddo
         ileft = i
         
         select case (keyword)
         
           case('MATE')
            im = im + 1
            read(inline(ileft:iright),*) lab_mat(im),type_mat(im),&
                 rho, VS, VP, & !char_mat(im,4)
                 Qua_S(im), Qua_P(im)

                 
                 char_mat(im,1) = rho
                 !lambda
                 char_mat(im,2) = rho * (VP**2 - 2*VS**2) 
                 !mu 
                 char_mat(im,3) = rho * VS**2
                                  
           case('MATN')                                                                        
            im_nle = im_nle + 1                                                                                        
            read(inline(ileft:iright),*) lab_mat_nle(im_nle),type_mat_nle(im_nle),&                                
                 char_mat_nle(im_nle,1),val_mat_nle(im_nle,1)                         
         
           case('DIRX')
            idX = idX + 1
            read(inline(ileft:iright),*) lab_diriX(idX),fnc_diriX(idX),&
                 val_diriX(idX,1),val_diriX(idX,2),val_diriX(idX,3),val_diriX(idX,4)

           case('DIRY')
            idY = idY + 1
            read(inline(ileft:iright),*) lab_diriY(idY),fnc_diriY(idY),&
                 val_diriY(idY,1),val_diriY(idY,2),val_diriY(idY,3),val_diriY(idY,4)

           case('DIRZ')
            idZ = idZ + 1
            read(inline(ileft:iright),*) lab_diriZ(idZ),fnc_diriZ(idZ),&
                 val_diriZ(idZ,1),val_diriZ(idZ,2),val_diriZ(idZ,3),val_diriZ(idZ,4)
            !!!!!!!!!!!!!!!!!!modified by ty 170410
           case('NEUX')
            inX = inX + 1
            read(inline(ileft:iright),*) lab_neuX(inX),fnc_neuX(inX),&
                 val_neuX(inX,1),val_neuX(inX,2),val_neuX(inX,3),val_neuX(inX,4),&
                 val_neuX(inX,5),val_neuX(inX,6)

           case('NEUY')
            inY = inY + 1
            read(inline(ileft:iright),*) lab_neuY(inY),fnc_neuY(inY),&
                 val_neuY(inY,1),val_neuY(inY,2),val_neuY(inY,3),val_neuY(inY,4),&
                 val_neuY(inY,5),val_neuY(inY,6)

           case('NEUZ')
            inZ = inZ + 1
            read(inline(ileft:iright),*) lab_neuZ(inZ),fnc_neuZ(inZ),&
                 val_neuZ(inZ,1),val_neuZ(inZ,2),val_neuZ(inZ,3),val_neuZ(inZ,4),&
                 val_neuZ(inZ,5),val_neuZ(inZ,6)

           case('NEUN')                                         
            inN = inN + 1                                                         
            read(inline(ileft:iright),*) lab_neuN(inN),fnc_neuN(inN),&                 
                 val_neuN(inN,1),val_neuN(inN,2),val_neuN(inN,3),val_neuN(inN,4),&
                 val_neuN(inN,5),val_neuN(inN,6)             

           case('PLOX')
            ipX = ipX + 1
            read(inline(ileft:iright),*) fnc_poiX(ipX),&
                 val_poiX(ipX,1),val_poiX(ipX,2),val_poiX(ipX,3),val_poiX(ipX,4),&
                 val_poiX(ipX,5),val_poiX(ipX,6)
           case('PLOY')
            ipY = ipY + 1
            read(inline(ileft:iright),*) fnc_poiY(ipY),&
                 val_poiY(ipY,1),val_poiY(ipY,2),val_poiY(ipY,3),val_poiY(ipY,4),&
                 val_poiY(ipY,5),val_poiY(ipY,6)

           case('PLOZ')
            ipZ = ipZ + 1
            read(inline(ileft:iright),*) fnc_poiZ(ipZ),&
                 val_poiZ(ipZ,1),val_poiZ(ipZ,2),val_poiZ(ipZ,3),val_poiZ(ipZ,4),&
                 val_poiZ(ipZ,5),val_poiZ(ipZ,6)
            !!!!!!!!!!!!!!!!!!modified by ty 170410
           case('TLOX') 
            itX = itX + 1
            read(inline(ileft:iright),*) ftX(itX),valtX(itX,1),valtX(itX,4)

           case('TLOY') 
            itY = itY + 1
            read(inline(ileft:iright),*) ftY(itY),valtY(itY,1),valtY(itY,4)

           case('TLOZ') 
            itZ = itZ + 1
            read(inline(ileft:iright),*) ftZ(itZ),valtZ(itZ,1),valtZ(itZ,4)
 
           case('PLAX')                                        
            iplX = iplX + 1                                                        
            read(inline(ileft:iright),*)fnc_plaX(iplX),&                        
                 lab_plaX(iplX),val_plaX(iplX,1)                                        
 
           case('PLAY')                                        
            iplY = iplY + 1                                                        
            read(inline(ileft:iright),*)fnc_plaY(iplY),&                        
                 lab_plaY(iplY),val_plaY(iplY,1)                                        
 
           case('PLAZ')
            iplZ = iplZ + 1
            read(inline(ileft:iright),*)fnc_plaZ(iplZ),&
                 lab_plaZ(iplZ),val_plaZ(iplZ,1)

           case('FORX')
            ifX = ifX + 1
            read(inline(ileft:iright),*) fnc_forX(ifX),&
                 val_forX(ifX,1),val_forX(ifX,2),val_forX(ifX,3),val_forX(ifX,4)

           case('FORY')
            ifY = ifY + 1
            read(inline(ileft:iright),*) fnc_forY(ifY),&
                 val_forY(ifY,1),val_forY(ifY,2),val_forY(ifY,3),val_forY(ifY,4)

           case('FORZ')
            ifZ = ifZ + 1
            read(inline(ileft:iright),*) fnc_forZ(ifZ),&
                 val_forZ(ifZ,1),val_forZ(ifZ,2),val_forZ(ifZ,3),val_forZ(ifZ,4)

           case('TEST')                                        
            itest = itest + 1                                                        
            read(inline(ileft:iright),*) ftest(itest)

           case('FORC')
            iff = iff + 1
            read(inline(ileft:iright),*) fnc_forF(iff),&
                 val_forF(iff,1),val_forF(iff,2),val_forF(iff,3),val_forF(iff,4),&
                 val_forF(iff,5),val_forF(iff,6),val_forF(iff,7),&
                 val_forF(iff,8),val_forF(iff,9),val_forF(iff,10)

           case('PRES')
            ipr = ipr + 1
            read(inline(ileft:iright),*) fnc_pre(ipr),&
                 val_pre(ipr,1),val_pre(ipr,2),val_pre(ipr,3),val_pre(ipr,4),&
                 val_pre(ipr,5),val_pre(ipr,6),val_pre(ipr,7),&
                 val_pre(ipr,8),val_pre(ipr,9),val_pre(ipr,10)

           case('SHEA')
            ish = ish + 1
            read(inline(ileft:iright),*) fnc_she(ish),&
                 val_she(ish,1),val_she(ish,2),val_she(ish,3),val_she(ish,4),&
                 val_she(ish,5),val_she(ish,6),val_she(ish,7),&
                 val_she(ish,8),val_she(ish,9),val_she(ish,10)

           case('ABSO')
            iabc = iabc + 1
            read(inline(ileft:iright),*) lab_abc(iabc)

           case('DGIC')
            idg = idg + 1
            read(inline(ileft:iright),*) lab_dg(idg), lab_dg_yn(idg)
           case('SISM')        
              isism = isism + 1                                                                
              read(inline(ileft:iright),*) fnc_sism(isism),&                                
                         lab_sism(isism),val_sism(isism,1),val_sism(isism,2),&                        
                         val_sism(isism,3),val_sism(isism,4),val_sism(isism,5),&                
                         val_sism(isism,6),val_sism(isism,7),val_sism(isism,8),&                
                         val_sism(isism,9),val_sism(isism,10),val_sism(isism,11),&                
                         val_sism(isism,12),val_sism(isism,13),val_sism(isism,14),&        
                         val_sism(isism,15),val_sism(isism,16),val_sism(isism,17),&        
                         val_sism(isism,18),val_sism(isism,19),val_sism(isism,20),&        
                         val_sism(isism,21)
         

           case('EXPL')                                                        
              iexpl = iexpl + 1                                                                
              read(inline(ileft:iright),*) fnc_expl(iexpl),&                                
                         lab_expl(iexpl),val_expl(iexpl,1),val_expl(iexpl,2),&                        
                         val_expl(iexpl,3),val_expl(iexpl,4),val_expl(iexpl,5),&                
                         val_expl(iexpl,6),val_expl(iexpl,7),val_expl(iexpl,8),&                
                         val_expl(iexpl,9),val_expl(iexpl,10),val_expl(iexpl,11),&                
                         val_expl(iexpl,12),val_expl(iexpl,13),val_expl(iexpl,14),&        
                         val_expl(iexpl,15),val_expl(iexpl,16),val_expl(iexpl,17),&        
                         val_expl(iexpl,18),val_expl(iexpl,19),val_expl(iexpl,20)                

           case('CASE')                                                        
              icase = icase + 1                                                                
              read(inline(ileft:iright),*) &                                                
                 lab_case,val_case,tol_case

           case('FMAX') 
              read(inline(ileft:iright),*) fmax                                        

           case('FPEK') 
              read(inline(ileft:iright),*) fpeak        
           
           case('FUNC')
            ifunc = ifunc + 1
            read(inline(ileft:iright),*) lab_fnc(ifunc),type_fnc(ifunc)
            
            select case (type_fnc(ifunc))
            
               case(0)
                 ind_fnc(ifunc +1) = ind_fnc(ifunc) + 0 
               
               case(1)
                 ind_fnc(ifunc +1) = ind_fnc(ifunc) + 2
                 read(inline(ileft:iright),*)dummy,dummy,&
                     (dat_fnc(j), j = ind_fnc(ifunc),ind_fnc(ifunc +1) -1)
               
               case(2)
                 ind_fnc(ifunc +1) = ind_fnc(ifunc) + 2
                 read(inline(ileft:iright),*)dummy,dummy,&
                    (dat_fnc(j), j = ind_fnc(ifunc),ind_fnc(ifunc +1) -1)
                    
               case(3)
                 read(inline(ileft:iright),*)dummy,dummy,ndat_fnc,fileinput
                 ind_fnc(ifunc +1) = ind_fnc(ifunc) + 2*ndat_fnc

                 open(24,file=fileinput)
                 read(24,*) file_nd
                 if (ndat_fnc .ne. file_nd) then
                    write(*,*) 'Error reading function ! ndat_fnc .ne. file_nd !'
                    write(*,*) 'Error reading function from file ', trim(fileinput), '!'
                    write(*,*) 'Line numbers not consistent with material file.'
                    call EXIT(EXIT_FUNCTION_ERROR)
                 endif
                 do j = 1, file_nd
                    i = ind_fnc(ifunc) + 2*(j -1)
                    read(24,*) dat_fnc(i), dat_fnc(i +1)
                 enddo
                 close(24)
              
              !!!!!!!!!!!!!ty!!!!!!!!!!!!!!!!!!!
              case(777)
                  read(inline(ileft:iright),*)dummy,dummy,ndat_fnc,dat_fnc(ind_fnc(ifunc)),dat_fnc(ind_fnc(ifunc)+1) !
                  ind_fnc(ifunc +1) = ind_fnc(ifunc) + 2*ndat_fnc
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!
              case(773)
                 read(inline(ileft:iright),*)dummy,dummy,ndat_fnc,fileinput
                 ind_fnc(ifunc +1) = ind_fnc(ifunc) + ndat_fnc

                 open(24,file=fileinput)
                 read(24,*) file_nd
                 if (ndat_fnc .ne. file_nd) then
                    write(*,*) 'Error reading function ! ndat_fnc .ne. file_nd !'
                    write(*,*) 'Error reading function from file ', trim(fileinput), '!'
                    write(*,*) 'Line numbers not consistent with material file.'
                    call EXIT(EXIT_FUNCTION_ERROR)
                 endif
                 do j = 1, file_nd
                    i = ind_fnc(ifunc) + (j -1)
                    read(24,*) dat_fnc(i)
                 enddo
                 close(24)
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               case(4)
                 ind_fnc(ifunc +1) = ind_fnc(ifunc) + 2
                 read(inline(ileft:iright),*)dummy,dummy,&
                    (dat_fnc(j), j = ind_fnc(ifunc),ind_fnc(ifunc +1) -1)
            
               case(6)
                 ind_fnc(ifunc +1) = ind_fnc(ifunc) + 2
                 read(inline(ileft:iright),*)dummy,dummy,&
                    (dat_fnc(j), j = ind_fnc(ifunc),ind_fnc(ifunc +1) -1)
                        
               !SIGMOIDAL FUNCTION
               case(12)
                 ind_fnc(ifunc +1) = ind_fnc(ifunc) + 3
                 read(inline(ileft:iright),*)dummy,dummy,&
                    (dat_fnc(j), j = ind_fnc(ifunc),ind_fnc(ifunc +1) -1)

               !GRENOBLE BENCHMARK
                  case(13)
                 ind_fnc(ifunc +1) = ind_fnc(ifunc) + 2
                 read(inline(ileft:iright),*)dummy,dummy,&
                    (dat_fnc(j), j = ind_fnc(ifunc),ind_fnc(ifunc +1) -1)

               !SCEC BENCHMARK
               case(14)
                 ind_fnc(ifunc +1) = ind_fnc(ifunc) + 2
                 read(inline(ileft:iright),*)dummy,dummy,&
                    (dat_fnc(j), j = ind_fnc(ifunc),ind_fnc(ifunc +1) -1)
            
               !ERF FUNCTION
               case(15)
                 ind_fnc(ifunc +1) = ind_fnc(ifunc) + 4
                 read(inline(ileft:iright),*)dummy,dummy,&
                    (dat_fnc(j), j = ind_fnc(ifunc),ind_fnc(ifunc +1) -1)
               
               !ADD ARCHULETA FUNC     
               case(21)
                 ind_fnc(ifunc +1) = ind_fnc(ifunc) + 0
                    
                    
               case(50,55)
                 ind_fnc(ifunc +1) = ind_fnc(ifunc) + 2
                 read(inline(ileft:iright),*)dummy,dummy,&
                    (dat_fnc(j), j = ind_fnc(ifunc),ind_fnc(ifunc +1) -1)                    

               case(60,62)                                   
                 read(inline(ileft:iright),*)dummy,dummy,ndat_fnc                 
                 ind_fnc(ifunc +1) = ind_fnc(ifunc) + 2*ndat_fnc                  
                 read(inline(ileft:iright),*)dummy,dummy,dummy,&                 
                    (dat_fnc(j), j = ind_fnc(ifunc),ind_fnc(ifunc +1) -1)    
               
               case(61,63)                                   
                 read(inline(ileft:iright),*)dummy,dummy,ndat_fnc                 
                 ind_fnc(ifunc +1) = ind_fnc(ifunc) + 2*ndat_fnc                  
                 read(inline(ileft:iright),*)dummy,dummy,dummy,&                 
                    (dat_fnc(j), j = ind_fnc(ifunc),ind_fnc(ifunc +1) -1)     
               
                  case(99)
                 ind_fnc(ifunc +1) = ind_fnc(ifunc) + 2
                 read(inline(ileft:iright),*)dummy,dummy,&
                    (dat_fnc(j), j = ind_fnc(ifunc),ind_fnc(ifunc +1) -1)
                    
                  case(100)
                    ind_fnc(ifunc +1) = ind_fnc(ifunc) + 1
                 read(inline(ileft:iright),*) dummy,dummy,&
                    (dat_fnc(j), j = ind_fnc(ifunc),ind_fnc(ifunc +1) -1)                    
            
            end select
                                           
         end select

      enddo
      
      close(40)
      return
      
      end subroutine READ_FILEMATE
