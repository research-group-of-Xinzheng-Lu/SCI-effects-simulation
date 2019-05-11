!> @brief Output acc for buildings. by Y TIAN

subroutine WRITE_OUTPUT_MDOF(nmonitlst, mpi_id, elem_mlst, local_el_num, ne_loc, &
                        cs_loc, cs_nnz_loc, sdeg_mat, nmat, &
                        u2, u1, u0, v1, nnod_loc, &
                        xr_mlst, yr_mlst, zr_mlst, & 
                        tt1, dt, dt2, &
                        count_mon, &
                        MDOFinputty,mpi_np, &
                        uty1, uty2, uty3,nmonity,finternal,ftyty1,ftyty2,MDOFinputfkty)
                              
     
implicit none
     
integer*4 :: imon,ielem,ie,im,nn,k,j,i,is,in,iaz,count_mon,ishift,jshift,kshift,dbg
integer*4 :: nmonitlst, mpi_id, ne_loc, nmat, nnod_loc, cs_nnz_loc,mpi_np
integer*4, dimension(0:cs_nnz_loc) :: cs_loc
integer*4, dimension(nmonitlst) :: elem_mlst
integer*4, dimension(nmat) :: sdeg_mat
integer*4, dimension(ne_loc) :: local_el_num
integer*4 :: mpi_ierr
real*8 :: tt1, dt, dt2
real*8 :: uxm,uym,uzm
real*8 :: vxm,vym,vzm
real*8 :: axm,aym,azm
real*8 :: fxm,fym,fzm
real*8, dimension(:), allocatable :: ct, ww
real*8, dimension(:,:), allocatable :: dd
real*8, dimension(:,:,:), allocatable :: fkx_el,fky_el,fkz_el !ux_el, uy_el, uz_el, 
real*8, dimension(3*nnod_loc),intent(in) :: finternal      
real*8, dimension(3*nnod_loc) :: u1, u0, u2, v1 !, omega
real*8, dimension(nmonitlst) ::xr_mlst, yr_mlst, zr_mlst

!!!ty!!!
integer*4 :: tynpmpi,nmonity
real*8, dimension(nmonity*3) :: MDOFinputty,MDOFinputfkty
real*8, dimension(3*count_mon) :: uty1,uty2,uty3,ftyty1,ftyty2 !!!ty ftyty2:n-1,ftyty1:n
real*8, dimension(:,:,:), allocatable :: ux_el2, uy_el2, uz_el2
!!!ty!!!
uty3(1:3*count_mon)=uty2(1:3*count_mon)
uty2(1:3*count_mon)=uty1(1:3*count_mon)
ftyty2(1:3*count_mon)=ftyty1(1:3*count_mon)


ishift = 0
     
do imon = 1,nmonity ! nmonitlst
        
ielem = elem_mlst(imon)     !ie = index hexahedra containing monitor
call GET_INDLOC_FROM_INDGLO(local_el_num, ne_loc, ielem, ie)                                             
if (ie .ne. 0) then
            
    im = cs_loc(cs_loc(ie -1) +0)
    nn = sdeg_mat(im) +1

    allocate(ct(nn),ww(nn),dd(nn,nn))
    !!!ty
    allocate(fkx_el(nn,nn,nn),fky_el(nn,nn,nn),fkz_el(nn,nn,nn))
    allocate(ux_el2(nn,nn,nn),uy_el2(nn,nn,nn),uz_el2(nn,nn,nn))
    !!!ty
    call MAKE_LGL_NW(nn,ct,ww,dd)

    do k = 1,nn
        do j = 1,nn
            do i = 1,nn
                is = nn*nn*(k -1) +nn*(j -1) +i
                in = cs_loc(cs_loc(ie -1) + is)

                iaz = 3*(in -1) +1
                !!!ty
                ux_el2(i,j,k) = u2(iaz)
                fkx_el(i,j,k) = finternal(iaz)
                !!!ty
                iaz = 3*(in -1) +2
                !!!ty
                uy_el2(i,j,k) = u2(iaz)
                fky_el(i,j,k) = finternal(iaz)
                !!!ty
                iaz = 3*(in -1) +3
                !!!ty
                uz_el2(i,j,k) = u2(iaz)
                fkz_el(i,j,k) = finternal(iaz)
                !!!ty
            enddo
        enddo
    enddo
       
                        
    !-------------------------------------------------------------
    ! ACCELERATION
    !-------------------------------------------------------------
        call GET_MONITOR_VALUE(nn,ct,ux_el2,&
                        xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),uxm)
              
        call GET_MONITOR_VALUE(nn,ct,uy_el2,&
                        xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),uym)
               
        call GET_MONITOR_VALUE(nn,ct,uz_el2,&
                        xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),uzm)

        call GET_MONITOR_VALUE(nn,ct,fkx_el,&
                        xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),fxm)
              
        call GET_MONITOR_VALUE(nn,ct,fky_el,&
                        xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),fym)
               
        call GET_MONITOR_VALUE(nn,ct,fkz_el,&
                        xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),fzm)
        if (dabs(uxm).lt.1.0e-30) uxm = 0.0e+00
        if (dabs(uym).lt.1.0e-30) uym = 0.0e+00
        if (dabs(uzm).lt.1.0e-30) uzm = 0.0e+00

        uty1(ishift+1) = uxm
        uty1(ishift+2) = uym
        uty1(ishift+3) = uzm
        ftyty1(ishift+1) = fxm
        ftyty1(ishift+2) = fym
        ftyty1(ishift+3) = fzm
            axm = (uty1(ishift+1) -2.0d0*uty2(ishift+1) +uty3(ishift+1)) / dt2
            aym = (uty1(ishift+2) -2.0d0*uty2(ishift+2) +uty3(ishift+2)) / dt2
            azm = (uty1(ishift+3) -2.0d0*uty2(ishift+3) +uty3(ishift+3)) / dt2
            if (dabs(axm).lt.1.0e-30) axm = 0.0e+00
            if (dabs(aym).lt.1.0e-30) aym = 0.0e+00
            if (dabs(azm).lt.1.0e-30) azm = 0.0e+00

    !-------------------------------------------------------------
        !!!ty!!!
        if(imon.le.nmonity) then
        do tynpmpi=1,mpi_np
            if (mpi_id .eq. (tynpmpi-1)) then
            MDOFinputty(3*imon-2)=axm
            MDOFinputty(3*imon-1)=aym
            MDOFinputty(3*imon)=azm
            MDOFinputfkty(3*imon-2)=(ftyty1(ishift+1)-ftyty2(ishift+1))/dt
            MDOFinputfkty(3*imon-1)=(ftyty1(ishift+2)-ftyty2(ishift+2))/dt
            MDOFinputfkty(3*imon)=(ftyty1(ishift+3)-ftyty2(ishift+3))/dt
            end if
        end do
        end if
        !!!ty!!!
    deallocate(ct,ww,dd)
    !!!ty
    deallocate(fkx_el,fky_el,fkz_el)
    deallocate(ux_el2,uy_el2,uz_el2)
    !!!ty
    ishift = ishift + 3


endif
enddo
        
end subroutine WRITE_OUTPUT_MDOF
