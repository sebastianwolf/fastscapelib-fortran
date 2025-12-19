#include "Error.fpp"
subroutine Advect3d_lag (ierr)

  use FastScapeContext
  use omp_lib

  ! routine to advect the topography (h), basement (hb) and total erosion(etot)
  ! by a velocity field vx, vy, vz known at the same locations on a rectangular
  ! grid of nx by ny points separated by dx, dy over a time step dt

  implicit none

  integer ierr
  double precision dx,dy,advect_dt
  double precision :: dtime_in, dtime_out

  dx=xl/(nx-1)
  dy=yl/(ny-1)
  !advect_dt=advect_every_step*dt
  ! Useful if dt is variable
  advect_dt=totaltime-totaltime_before_advection

  if (step == 0) call cloud_setup(ierr)

  if ( (mod(step+1,advect_every_step)==0) ) then

   write(*,'(a,x,i6,x,f9.3,x,f13.3,x,f13.3)') 'Advection 3d step,advect_dt',step,advect_dt,totaltime,totaltime_before_advection
  totaltime_before_advection=totaltime

  write(*,*) '================= start advection ======================'
  dtime_in = omp_get_wtime()
  call advect_cloud(advect_dt,ierr)
  dtime_out = omp_get_wtime()
  write(*,*) 'advect_cloud ',dtime_out-dtime_in,' s'

  dtime_in = omp_get_wtime()
  call locate_cloud(ierr)
  dtime_out = omp_get_wtime()
  write(*,*) 'locate_cloud ',dtime_out-dtime_in,' s'

  dtime_in = omp_get_wtime()
  call update_cloud(advect_dt,ierr)
  dtime_out = omp_get_wtime()
  write(*,*) 'update_cloud ',dtime_out-dtime_in,' s'

!  dtime_in = omp_get_wtime()
!  call locate_cloud(ierr)
!  dtime_out = omp_get_wtime()
!  write(*,*) 'locate_cloud ',dtime_out-dtime_in,' s'

  dtime_in = omp_get_wtime()
  call cloud_to_eul(ierr)
  dtime_out = omp_get_wtime()
  write(*,*) 'cloud_to_eul ',dtime_out-dtime_in,' s'

  write(*,*) '================= end advection ======================'

  !call write_cell(2410)

  endif ! if (mod(step,advect_every_step)==0)

  return  

end subroutine Advect3d_lag

!--------------------------------------------------------------------------

subroutine write_cell(ipg)

  use FastScapeContext

  implicit none

  integer ipg, i, j, k, ip, ic, ncell
  integer, dimension(4) :: icell
  character(len=25) :: filename,filename2 

  j        = ceiling(real(ipg)/real(nx))
  i        = ipg - (j-1)*nx
  icell(1) = (j-2)*(nx-1) + i - 1 
  icell(2) = icell(1) + 1
  icell(3) = (j-1)*(nx-1) + i
  icell(4) = icell(3) - 1 
  ncell    = (nx-1)*(ny-1)

  write(filename,'(a,I5.5,a)') 'grid_',step,'.txt' 
  write(filename2,'(a,I5.5,a)') 'cloud_',step,'.txt'

  open (1, file = filename, status = 'new')
  open (2, file = filename2, status = 'new')
  do i=1,4
    ic = icell(i)
    if (ic>0 .and. ic <= ncell) then
      do j=1,4
        write(1,*) grid%x(grid%icon(j,ic)),grid%y(grid%icon(j,ic))
      enddo
      do k=1,grid%nn(ic)
        ip = grid%pair(k,ic)
        write(2,*) cl%x(ip), cl%y(ip)
      enddo
    endif
  enddo
  close(1) ; close(2)

end subroutine write_cell

!--------------------------------------------------------------------------

subroutine cloud_setup (ierr)

  use FastScapeContext

  implicit none
  
  integer, intent(out) :: ierr
  integer i, j, ic, ncell, np_per_cell, icy, icx, ipx, ipy, counter, np
  integer multiplicator, ndim, mpe, ncellx, ncelly
  double precision dx, dy, dxx, dyy
  double precision deltafx, deltafy, deltafxy
  double precision deltafx_b, deltafy_b, deltafxy_b
  double precision deltafx_etot, deltafy_etot, deltafxy_etot

  ierr = 0

  ncellx        = nx-1
  ncelly        = ny-1
  ncell         = ncellx*ncelly
  multiplicator = 3 
  np_per_cell   = (multiplicator+1)**2
  dx            = xl/(nx-1)
  dy            = yl/(ny-1)
  np            = ncell*np_per_cell
  ndim          = 2
  mpe           = 2**ndim
  grid%nmax     = int(38) !int(48)
  grid%nmin     = int(9)  !np_per_cell

  write(*,*) 'np_per_cell = ', np_per_cell

  allocate(grid%icon(mpe,ncell))
  allocate(grid%x(nn))
  allocate(grid%y(nn))
  allocate(grid%nn(ncell))
  allocate(grid%pair(grid%nmax,ncell))

  do counter=1,nn
    j = ceiling(real(counter)/real(nx))
    i = counter - (j-1)*nx
    grid%x(counter)=(i-1)*dx
    grid%y(counter)=(j-1)*dy
  enddo

  !Build cell to global node mapping
  counter=0
  do j=1,ncelly
    do i=1,ncellx
       counter=counter+1
       grid%icon(1,counter)=i   + (j-1)*(ncellx+1)
       grid%icon(2,counter)=i+1 + (j-1)*(ncellx+1)
       grid%icon(3,counter)=i+1 + j    *(ncellx+1)
       grid%icon(4,counter)=i   + j    *(ncellx+1)
    end do
  end do

  grid%nn = 0

  allocate (cl%x(np),cl%y(np),cl%h(np),cl%b(np),cl%etot(np),cl%erate(np))
  allocate (cl%icx(np),cl%icy(np),cl%active(np),cl%cell(np))
  allocate (cl%closest_node(np))
  cl%closest_node = -1
 
  counter = 0
  do ic=1,ncell
    icy = ceiling(real(ic)/real(nx-1))
    icx = ic - (icy-1)*(nx-1)
    deltafx  = h2(icx+1,icy) - h2(icx,icy)
    deltafy  = h2(icx,icy+1) - h2(icx,icy)
    deltafxy = h2(icx,icy) + h2(icx+1,icy+1) - h2(icx+1,icy) - h2(icx,icy+1)

    deltafx_b  = b2(icx+1,icy) - b2(icx,icy)
    deltafy_b  = b2(icx,icy+1) - b2(icx,icy)
    deltafxy_b = b2(icx,icy) + b2(icx+1,icy+1) - b2(icx+1,icy) - b2(icx,icy+1)

    deltafx_etot  = etot2(icx+1,icy) - etot2(icx,icy)
    deltafy_etot  = etot2(icx,icy+1) - etot2(icx,icy)
    deltafxy_etot = etot2(icx,icy) + etot2(icx+1,icy+1) - etot2(icx+1,icy) - etot2(icx,icy+1)

    do ipy=1,multiplicator+1
      do ipx=1,multiplicator+1
        counter     = counter + 1
        grid%nn(ic) = grid%nn(ic) + 1
        grid%pair(grid%nn(ic),ic) = counter

        if (ipx==1) then
          cl%x(counter) = (icx-1)*dx + dx/(10.d0*multiplicator)
        elseif (ipx==multiplicator+1) then
          cl%x(counter) = icx*dx - dx/(10.d0*multiplicator)
        else
          cl%x(counter) = (icx-1)*dx + (ipx-1)*(dx/multiplicator)
        endif

        if (ipy==1) then
          cl%y(counter) = (icy-1)*dy + dy/(10.d0*multiplicator)
        elseif (ipy==multiplicator+1) then
          cl%y(counter) = icy*dy - dy/(10.d0*multiplicator)
        else
          cl%y(counter) = (icy-1)*dy + (ipy-1)*(dy/multiplicator)
        endif

        ! bilinear interpolation
        dxx           = cl%x(counter) - (icx-1)*dx
        dyy           = cl%y(counter) - (icy-1)*dy 
        cl%h(counter)    = deltafx*(dxx/dx)      + deltafy*(dyy/dy)      + deltafxy*(dxx*dyy/(dx*dy)) + h2(icx,icy) 
        cl%b(counter)    = deltafx_b*(dxx/dx)    + deltafy_b*(dyy/dy)    + deltafxy_b*(dxx*dyy/(dx*dy)) + b2(icx,icy) 
        cl%etot(counter) = deltafx_etot*(dxx/dx) + deltafy_etot*(dyy/dy) + deltafxy_etot*(dxx*dyy/(dx*dy)) + etot2(icx,icy) 

        cl%icx(counter)  = icx 
        cl%icy(counter)  = icy
        cl%cell(counter) = ic
      enddo
    enddo
  enddo

  cl%npcl   = counter
  cl%active = .true.
  !cl%b      = cl%h ! not correct in case of restart
  !cl%etot   = 0.d0
  cl%erate  = 0.d0
  
  return

end subroutine cloud_setup

!--------------------------------------------------------------------------

subroutine advect_cloud (advect_dt,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out) :: ierr
  double precision, intent(in) :: advect_dt
  integer ic, i, ip, ncell, inode1, inode2, inode3, inode4
  double precision xnode1, r, s, ymin, ymax, dx, dy, xip
  double precision N1, N2, N3, N4

  ierr  = 0
  ncell = (nx-1)*(ny-1)
  dx    = xl/(nx-1)
  dy    = yl/(ny-1)

  if (allocated(cl%ip_leaving)) deallocate(cl%ip_leaving)
  allocate(cl%ip_leaving(cl%npcl))

  cl%nleaving   = 0
  cl%ip_leaving = 0

  !$omp parallel do shared(ncell,grid,cl,vx,vy,u) private(ic,inode1,inode2,inode3,inode4,xnode1,i,ip,xip,r,ymin,ymax,s,N1,N2,N3,N4)
  do ic=1,ncell
    inode1=grid%icon(1,ic)
    inode2=grid%icon(2,ic)
    inode3=grid%icon(3,ic)
    inode4=grid%icon(4,ic)
    xnode1=grid%x(inode1)
    do i=1,grid%nn(ic)
       ip=grid%pair(i,ic)
       xip=cl%x(ip)
       r=((xip-xnode1)/dx -0.5d0 ) *2.d0
       ymin=(grid%y(inode2)-grid%y(inode1))/dx*(xip-xnode1) + grid%y(inode1)
       ymax=(grid%y(inode3)-grid%y(inode4))/dx*(xip-xnode1) + grid%y(inode4)
       s=((cl%y(ip)-ymin)/(ymax-ymin) -0.5d0 ) *2.d0
       N1=0.25d0*(1.d0-r)*(1.d0-s) 
       N2=0.25d0*(1.d0+r)*(1.d0-s) 
       N3=0.25d0*(1.d0+r)*(1.d0+s) 
       N4=0.25d0*(1.d0-r)*(1.d0+s) 
       cl%x(ip)=cl%x(ip)+advect_dt*(N1 * vx(inode1) + N2 * vx(inode2) + N3 * vx(inode3) + N4 * vx(inode4))
       cl%y(ip)=cl%y(ip)+advect_dt*(N1 * vy(inode1) + N2 * vy(inode2) + N3 * vy(inode3) + N4 * vy(inode4))
       cl%h(ip)=cl%h(ip)+advect_dt*(N1 * u(inode1)  + N2 * u(inode2)  + N3 * u(inode3)  + N4 * u(inode4))

       !Advect deep basement surface with surface velocities (not optimal but identical to other advect routines)
       cl%b(ip)=cl%b(ip)+advect_dt*(N1 * u(inode1)  + N2 * u(inode2)  + N3 * u(inode3)  + N4 * u(inode4))

       if (cl%x(ip)>grid%x(grid%icon(2,ic)) .or. cl%x(ip)<grid%x(grid%icon(1,ic)) .or. &
           cl%y(ip)>grid%y(grid%icon(3,ic)) .or. cl%y(ip)<grid%y(grid%icon(1,ic)) ) then
              !$omp critical
              cl%nleaving                = cl%nleaving + 1
              cl%ip_leaving(cl%nleaving) = ip
              !$omp end critical
       endif

    end do
  end do
  !$omp end parallel do

  return

end subroutine advect_cloud

!--------------------------------------------------------------------------

subroutine locate_cloud (ierr)

  use FastScapeContext

  implicit none

  integer, intent(out) :: ierr
  integer ncell, ip, ic

  integer i, icx, icy, ncellx, ncelly
  double precision dx, dy

  ierr = 0

  if ( cl%nleaving > 0 ) then

  ncellx = nx-1
  ncelly = ny-1
  ncell  = ncellx*ncelly

!  if (allocated(cl%cell)  ) deallocate(cl%cell  )
  if (allocated(cl%active)) deallocate(cl%active)
  
!  allocate(cl%cell  (cl%npcl))
  allocate(cl%active(cl%npcl))
  cl%active=.true.
  
  !call search_cell2D (cl%npcl,cl%x,cl%y,cl%cell,cl%icy,cl%active,ierr)
  !TT update AAA 
  dx     = xl/(nx-1)
  dy     = yl/(ny-1)
  !write(*,'(a,2i15)')    'before sum(grid%nn), cl%npcl :', sum(grid%nn), cl%npcl
  !$omp parallel do shared(cl,grid,dx,dy) private(i,ip,icx,icy,ic)
  do i=1,cl%nleaving
     ip = cl%ip_leaving(i)
     icx=floor(cl%x(ip)/dx)+1
     icy=floor(cl%y(ip)/dy)+1
     if (icx>ncellx .or. icx<1 .or. icy>ncelly .or. icy<1) then
       cl%active(ip)        = .false.
       grid%nn(cl%cell(ip)) = grid%nn(cl%cell(ip)) - 1
       cycle
     endif
     ic = icx + (icy-1)*ncellx
     !$omp critical
     grid%nn(cl%cell(ip)) = grid%nn(cl%cell(ip)) - 1     
     grid%nn(ic)          = grid%nn(ic) + 1
     !$omp end critical
     cl%icy(ip)           = icy
     cl%cell(ip)          = ic
  enddo
  !$omp end parallel do

  write(*,'(a,1i15)')    'nb points changing cell :', cl%nleaving
  write(*,'(a,1i15)')    'nb points outside       :', cl%npcl-count(cl%active)
  
  if (count(cl%active)<cl%npcl) then
    call cloud_trim (ierr)
  endif
  
  !TT This loops is not needed anymore after the update AAA
  !build cell to particle mapping
  !grid%nn = 0
  !do ip=1,cl%npcl
  !  ic          = cl%cell(ip)
  !  grid%nn(ic) = grid%nn(ic) + 1
  !enddo
  !write(*,'(a,2i15)')    'after  sum(grid%nn), cl%npcl :', sum(grid%nn), cl%npcl
  if (sum(grid%nn)/=cl%npcl) stop 'pb in locate_cloud_points2D'
  
  if (allocated(grid%pair)) deallocate(grid%pair)
  !allocate(grid%pair(grid%nmax,ncell))
  allocate(grid%pair(maxval(grid%nn),ncell))

  !write(*,'(a,2i15)') 'grid%nmax',grid%nmax
  
  !if (maxval(grid%nn)>grid%nmax) then
  !  print*,'some cells have too many particles than grid%pair'
  !  stop 'stop build pair'
  !endif
  
  !now it is safe to build grid%pair
  grid%nn = 0
  do ip=1,cl%npcl
    ic                            = cl%cell(ip)
    grid%nn(ic)                   = grid%nn(ic) + 1
    grid%pair(grid%nn(ic),ic)     = ip
  enddo

  write(*,'(a,2i15)')    'located grid%nn min/max',minval(grid%nn),maxval(grid%nn)

  endif ! cl%nleaving > 0
    
  return

end subroutine locate_cloud

!--------------------------------------------------------------------------

!subroutine search_cell2D (np,partx,party,cellid,icelly,valid,ierr)
!
!  use FastScapeContext
!
!  implicit none
!
!  integer                        :: np
!  double precision,dimension(np) :: partx,party
!  integer,dimension(np)          :: cellid,icelly
!  logical,dimension(np)          :: valid
!  integer, intent(out) :: ierr
!
!  integer ip,icx,icy,ic,ncellx, ncelly
!  double precision dx, dy
!
!  ierr   = 0
!  dx     = xl/(nx-1)
!  dy     = yl/(ny-1)
!  ncellx = nx-1
!  ncelly = ny-1
!
!  !$omp parallel do shared(np,partx,party,dx,dy,ncellx,ncelly,valid,icelly,cellid) private(ip,icx,icy,ic)
!  do ip=1,np
!    icx=floor(partx(ip)/dx)+1
!    icy=floor(party(ip)/dy)+1
!
!    if (icx>ncellx .or. icx<1 .or. icy>ncelly .or. icy<1) then
!      valid(ip) = .false.
!      cycle
!    endif
!
!    ic = icx + (icy-1)*ncellx
!
!    icelly(ip) = icy
!    cellid(ip) = ic
!  end do
!  !$omp end parallel do
!
!  return
!
!end subroutine search_cell2D
!
!--------------------------------------------------------------------------

subroutine cloud_trim (ierr)

  use FastScapeContext

  implicit none

  integer, intent(out) :: ierr
  logical,dimension(:),allocatable :: mask
  integer nout,npcl_new

  ierr = 0

  nout=cl%npcl-count(cl%active)
  
  write(*,'(a,i5)') 'cloud pts advected outside:',nout 
  
  npcl_new = cl%npcl-nout
  
  !===========================================
  
  if (nout > 0) then 
  
    allocate(mask(cl%npcl))
    mask = cl%active
  
    call array_cut(cl%x        ,mask)
    call array_cut(cl%y        ,mask)
    call array_cut(cl%h        ,mask)
    call array_cut(cl%b        ,mask)
    call array_cut(cl%etot     ,mask)
    call array_cut(cl%erate    ,mask)
    call array_cut(cl%icx      ,mask)
    call array_cut(cl%icy      ,mask)
    call array_cut(cl%cell     ,mask)
    call array_cut(cl%active   ,mask)
    call array_cut(cl%closest_node ,mask)
  
  end if !end if nout > 0
  
  deallocate(mask)
  
  cl%npcl = npcl_new

  return

end subroutine cloud_trim

!--------------------------------------------------------------------------

subroutine update_cloud (advect_dt,ierr)

  use FastScapeContext

  implicit none

  double precision, intent(in) :: advect_dt
  integer, intent(out) :: ierr
  integer ic, jcell, icell, ncellx, ncelly, ncell, nnic, npcl_max, npcl_min
  integer ninject, nremove, counter_inject
  integer ip
  integer, dimension(:), allocatable :: ninject_per_cell
  integer, dimension(:), allocatable :: nremove_per_cell
  type (cloud) :: clinject
  double precision ymax, ymin
  double precision xip, yip, xmin, xmax, distmin, dist, x, y
  logical, dimension(:), allocatable :: remove
  logical, dimension(:), allocatable :: not_remove
  logical, dimension(:), allocatable :: cell_injected
  integer ichoice, k, k2, ii, k1, ndist, npcl_new, nntot

  double precision distmin_per_quadrant(4), dx, dy, refdist
  integer ipart(4)
  integer npart, jcell_k, icell_k, ic_k, jl, il, nperline
  integer counter_interpolated, counter_averaged, counter_closest
  double precision dlx, dly

  !TT velocity dependent injection
  double precision inode1, inode2, inode3, inode4
  integer ntemp, ik, jk, pair(16)
  double precision xtemp(16), ytemp(16), htemp(16), Nnneighbours(16), rcut
  double precision btemp(16), etottemp(16), eratetemp(16)
  double precision N1,N2,N3,N4,r,s
  double precision xnode1,ynode1,xnode2,ynode2,xnode3,ynode3,xnode4,ynode4
  double precision hnode1,hnode2,hnode3,hnode4,highestnode, lowestnode

  ierr    = 0

  if (any(grid%nn<grid%nmin) .or. any(grid%nn>2*grid%nmax)) then

  ncelly  = ny-1
  ncellx  = nx-1
  ncell   = ncellx*ncelly

  dx      = xl/(nx-1)
  dy      = yl/(ny-1)
  refdist = sqrt(dx**2+dy**2)
  rcut    = 2.3d0*dx

  ninject = 0
  nremove = 0

  allocate(ninject_per_cell(ncell))
  allocate(nremove_per_cell(ncell))
  ninject_per_cell = 0
  nremove_per_cell = 0
  npcl_max         = grid%nmax
  npcl_min         = grid%nmin

  ! diagnosis
  !$omp parallel do shared(ncell,grid,ninject_per_cell,nremove_per_cell,npcl_max,npcl_min) private(ic,nnic) reduction(+:ninject,nremove)
  do ic=1,ncell
     nnic     = grid%nn(ic)
  
     if (nnic < npcl_min) then
       ninject_per_cell(ic) = floor(dble(npcl_min-nnic)/9.d0) * 9 + 9
     endif
     if (nnic > npcl_max) then
       nremove_per_cell(ic) = nnic-npcl_max
     endif
!     write(*,'(5I10)') ic, ninject_per_cell(ic), nremove_per_cell(ic), npcl_min, npcl_max

     ninject = ninject + ninject_per_cell(ic)
     nremove = nremove + nremove_per_cell(ic)

  end do
  !$omp end parallel do
  !end do

  write(*,'(a,i6,a,i6)') 'injecting ',ninject,' cloud points, maxpercell',maxval(ninject_per_cell)
  write(*,'(a,i6,a,i6)') 'removing  ',nremove,' cloud points, maxpercell',maxval(nremove_per_cell)

  npcl_new = cl%npcl + ninject-nremove

  write(*,'(a,2i9)') 'npcl new/old =',npcl_new,cl%npcl


  !injection ---------------------------------------------------------------------------
  if (ninject>0) then
    allocate(cell_injected(ncell))
    cell_injected = .false.

    allocate(clinject%x(ninject))
    allocate(clinject%y(ninject))
    allocate(clinject%h(ninject))
    allocate(clinject%b(ninject))
    allocate(clinject%etot(ninject))
    allocate(clinject%erate(ninject))
    allocate(clinject%cell(ninject))
    allocate(clinject%active(ninject))
    allocate(clinject%icy(ninject))
    allocate(clinject%icx(ninject))
    allocate(clinject%closest_node(ninject))
    clinject%closest_node = -1
    counter_inject       = 0
    counter_interpolated = 0
    counter_averaged     = 0
    counter_closest      = 0
    dly                  = dy/(3+2-1)

    !!! omp parallel do shared(ncell,ncellx,ncelly,grid,cl,ninject_per_cell,nremove_per_cell,dly,dx,clinject,refdist) private(ic,nnic,ii,il,jl,dlx,nperline,distmin_per_quadrant,ipart,k,dist,distmin,npart,icell,jcell,icell_k,jcell_k,ic_k,ip,xip,yip,xmin,xmax,r,ymin,ymax,s,N1,N2,N3,N4) reduction(+:counter_inject,counter_interpolated,counter_averaged)
    do ic=1,ncell
       jcell    = ceiling(real(ic)/real(ncellx))
       icell    = ic - (jcell-1)*ncellx
       nnic     = grid%nn(ic)
       nperline = ninject_per_cell(ic)/3
       dlx      = dx/(nperline+2-1)
       if (ninject_per_cell(ic) > 0) then
          cell_injected(ic) = .true.

          if ((icell==1).or.(icell==ncellx).or.(jcell==1).or.(jcell==ncelly)) then
             inode1=grid%icon(1,ic)
             inode2=grid%icon(2,ic)
             inode3=grid%icon(3,ic)
             inode4=grid%icon(4,ic)
             xnode1=grid%x(inode1) + advect_dt*vx(inode1)
             ynode1=grid%y(inode1) + advect_dt*vy(inode1)
             hnode1=h(inode1) + advect_dt*u(inode1)
             highestnode = hnode1
             lowestnode  = hnode1
             xnode2=grid%x(inode2) + advect_dt*vx(inode2)
             ynode2=grid%y(inode2) + advect_dt*vy(inode2)
             hnode2=h(inode2) + advect_dt*u(inode2)
             if (hnode2>highestnode) highestnode=hnode2
             if (hnode2<lowestnode) lowestnode=hnode2
             xnode3=grid%x(inode3) + advect_dt*vx(inode3)
             ynode3=grid%y(inode3) + advect_dt*vy(inode3)
             hnode3=h(inode3) + advect_dt*u(inode3)
             if (hnode3>highestnode) highestnode=hnode3
             if (hnode3<lowestnode) lowestnode=hnode3
             xnode4=grid%x(inode4) + advect_dt*vx(inode4)
             ynode4=grid%y(inode4) + advect_dt*vy(inode4)
             hnode4=h(inode4) + advect_dt*u(inode4)
             if (hnode4>highestnode) highestnode=hnode4
             if (hnode4<lowestnode) lowestnode=hnode4
     
             do ii=1,ninject_per_cell(ic)
                counter_inject             = counter_inject+1
                jl                         = ceiling(real(ii)/real(nperline))
                il                         = ii - (jl-1)*nperline
                clinject%x(counter_inject) = grid%x(grid%icon(1,ic)) + il*dlx
                clinject%y(counter_inject) = grid%y(grid%icon(1,ic)) + jl*dly
   
                counter_interpolated = counter_interpolated + 1
                xip=clinject%x(counter_inject)
                yip=clinject%y(counter_inject)
   
                xmin=(xnode4-xnode1)/(ynode4-ynode1)*(yip-ynode1) + xnode1
                xmax=(xnode3-xnode2)/(ynode3-ynode2)*(yip-ynode2) + xnode2
                r=((xip-xmin)/(xmax-xmin) -0.5d0 ) *2.d0
   
                ymin=(ynode2-ynode1)/(xnode2-xnode1)*(xip-xnode1) + ynode1
                ymax=(ynode3-ynode4)/(xnode3-xnode4)*(xip-xnode4) + ynode4
                s=((yip-ymin)/(ymax-ymin) -0.5d0 ) *2.d0
                N1=0.25d0*(1.d0-r)*(1.d0-s) 
                N2=0.25d0*(1.d0+r)*(1.d0-s) 
                N3=0.25d0*(1.d0+r)*(1.d0+s) 
                N4=0.25d0*(1.d0-r)*(1.d0+s) 
                clinject%h(counter_inject)     = N1 * hnode1     + N2 * hnode2     + N3 * hnode3     + N4 * hnode4
                clinject%b(counter_inject)     = N1 * b(grid%icon(1,ic))     + N2 * b(grid%icon(2,ic))     + N3 * b(grid%icon(3,ic))     + N4 * b(grid%icon(4,ic))
                clinject%etot(counter_inject)  = N1 * etot(grid%icon(1,ic))  + N2 * etot(grid%icon(2,ic))  + N3 * etot(grid%icon(3,ic))  + N4 * etot(grid%icon(4,ic))
                clinject%erate(counter_inject) = N1 * erate(grid%icon(1,ic)) + N2 * erate(grid%icon(2,ic)) + N3 * erate(grid%icon(3,ic)) + N4 * erate(grid%icon(4,ic))
    
                clinject%b(counter_inject)      = min(clinject%b(counter_inject),clinject%h(counter_inject))
                clinject%cell(counter_inject)   = ic
                clinject%active(counter_inject) = .true.
                clinject%icy(counter_inject)    = jcell
                clinject%icx(counter_inject)    = icell

                if (clinject%h(counter_inject)>highestnode) then
                    !write(*,'(a,1F10.3)') 'WARNING: something weird with elevation of new injected particle above surrounding nodes ', clinject%h(counter_inject)
                    !write(*,*) '4 cell nodes elevation are :'
                    !write(*,'(a,1F10.3,3I10)') ' node 1 ',hnode1, ic, icell, jcell
                    !write(*,*) ' node 2 ',hnode2
                    !write(*,*) ' node 3 ',hnode3
                    !write(*,*) ' node 4 ',hnode4
                    clinject%h(counter_inject) = highestnode
                    clinject%b(counter_inject) = min(clinject%b(counter_inject),clinject%h(counter_inject))
                endif
                if (clinject%h(counter_inject)<lowestnode) then
                    !write(*,'(a,1F10.3)') 'WARNING: something weird with elevation of new injected particle below surrounding nodes ', clinject%h(counter_inject)
                    !write(*,*) '4 cell nodes elevation are :'
                    !write(*,'(a,1F10.3,3I10)') ' node 1 ',hnode1, ic, icell, jcell
                    !write(*,*) ' node 2 ',hnode2
                    !write(*,*) ' node 3 ',hnode3
                    !write(*,*) ' node 4 ',hnode4
                    clinject%h(counter_inject) = lowestnode
                    clinject%b(counter_inject) = min(clinject%b(counter_inject),clinject%h(counter_inject))
                endif


   
             enddo !ii

             !TT we smooth all particles elevation in this cell to avoid
             !   unresolved small wavelength topographic features
             do k=1,nnic
                xip=cl%x(grid%pair(k,ic))
                yip=cl%y(grid%pair(k,ic))

                xmin=(xnode4-xnode1)/(ynode4-ynode1)*(yip-ynode1) + xnode1
                xmax=(xnode3-xnode2)/(ynode3-ynode2)*(yip-ynode2) + xnode2
                r=((xip-xmin)/(xmax-xmin) -0.5d0 ) *2.d0
   
                ymin=(ynode2-ynode1)/(xnode2-xnode1)*(xip-xnode1) + ynode1
                ymax=(ynode3-ynode4)/(xnode3-xnode4)*(xip-xnode4) + ynode4
                s=((yip-ymin)/(ymax-ymin) -0.5d0 ) *2.d0
                N1=0.25d0*(1.d0-r)*(1.d0-s) 
                N2=0.25d0*(1.d0+r)*(1.d0-s) 
                N3=0.25d0*(1.d0+r)*(1.d0+s) 
                N4=0.25d0*(1.d0-r)*(1.d0+s) 
                cl%h(grid%pair(k,ic))     = N1 * hnode1     + N2 * hnode2     + N3 * hnode3     + N4 * hnode4
                cl%b(grid%pair(k,ic))     = N1 * b(grid%icon(1,ic))     + N2 * b(grid%icon(2,ic))     + N3 * b(grid%icon(3,ic))     + N4 * b(grid%icon(4,ic))
                cl%etot(grid%pair(k,ic))  = N1 * etot(grid%icon(1,ic))  + N2 * etot(grid%icon(2,ic))  + N3 * etot(grid%icon(3,ic))  + N4 * etot(grid%icon(4,ic))
                cl%erate(grid%pair(k,ic)) = N1 * erate(grid%icon(1,ic)) + N2 * erate(grid%icon(2,ic)) + N3 * erate(grid%icon(3,ic)) + N4 * erate(grid%icon(4,ic))
                cl%b(grid%pair(k,ic))     = min(cl%b(grid%pair(k,ic)),cl%h(grid%pair(k,ic)))
                if (cl%h(grid%pair(k,ic))>highestnode) then
                    !write(*,'(a,1F10.3)') 'WARNING: something weird with elevation update of particle above surrounding nodes ', cl%h(grid%pair(k,ic))
                    !write(*,*) '4 cell nodes elevation are :'
                    !write(*,'(a,1F10.3,3I10)') ' node 1 ',hnode1,ic,icell,jcell
                    !write(*,*) ' node 2 ',hnode2
                    !write(*,*) ' node 3 ',hnode3
                    !write(*,*) ' node 4 ',hnode4
                    cl%h(grid%pair(k,ic)) = highestnode
                    cl%b(grid%pair(k,ic)) = min(cl%b(grid%pair(k,ic)),cl%h(grid%pair(k,ic)))
                endif
                if (cl%h(grid%pair(k,ic))<lowestnode) then
                    !write(*,'(a,1F10.3)') 'WARNING: something weird with elevation update of particle below surrounding nodes ', cl%h(grid%pair(k,ic))
                    !write(*,*) '4 cell nodes elevation are :'
                    !write(*,'(a,1F10.3,3I10)') ' node 1 ',hnode1,ic,icell,jcell
                    !write(*,*) ' node 2 ',hnode2
                    !write(*,*) ' node 3 ',hnode3
                    !write(*,*) ' node 4 ',hnode4
                    cl%h(grid%pair(k,ic)) = lowestnode
                    cl%b(grid%pair(k,ic)) = min(cl%b(grid%pair(k,ic)),cl%h(grid%pair(k,ic)))
                endif

             enddo

          else
             ntemp = 0
             highestnode = -9.e9
             lowestnode  = 9.e9
             do ik=1,3
                do jk=1,3
                    icell_k = icell + ik - 2
                    jcell_k = jcell + jk - 2
                    ic_k    = (jcell_k-1)*ncellx + icell_k
                    if ((icell_k<1).or.(jcell_k<1).or.(icell_k>ncellx).or.(jcell_k>ncelly)) cycle
   
                    inode1       = grid%icon(1,ic_k)
                    ntemp        = ntemp + 1
                    pair(ntemp)  = ntemp
                    xtemp(ntemp) = grid%x(inode1) + advect_dt*vx(inode1)
                    ytemp(ntemp) = grid%y(inode1) + advect_dt*vy(inode1)
                    htemp(ntemp) = h(inode1) + advect_dt*u(inode1)
                    if (htemp(ntemp)>highestnode) highestnode=htemp(ntemp)
                    if (htemp(ntemp)<lowestnode) lowestnode=htemp(ntemp)
                    btemp(ntemp) = b(inode1)
                    etottemp(ntemp)  = etot(inode1)
                    eratetemp(ntemp) = erate(inode1)
   
                    if (jk==3) then
                       inode4       = grid%icon(4,ic_k)
                       ntemp        = ntemp + 1
                       pair(ntemp)  = ntemp
                       xtemp(ntemp) = grid%x(inode4) + advect_dt*vx(inode4)
                       ytemp(ntemp) = grid%y(inode4) + advect_dt*vy(inode4)
                       htemp(ntemp) = h(inode4) + advect_dt*u(inode4)
                       if (htemp(ntemp)>highestnode) highestnode=htemp(ntemp)
                       if (htemp(ntemp)<lowestnode) lowestnode=htemp(ntemp)
                       btemp(ntemp) = b(inode4)
                       etottemp(ntemp)  = etot(inode4)
                       eratetemp(ntemp) = erate(inode4)
                    endif
   
                    if (ik==3) then
                       inode2       = grid%icon(2,ic_k)
                       ntemp        = ntemp + 1
                       pair(ntemp)  = ntemp
                       xtemp(ntemp) = grid%x(inode2) + advect_dt*vx(inode2)
                       ytemp(ntemp) = grid%y(inode2) + advect_dt*vy(inode2)
                       htemp(ntemp) = h(inode2) + advect_dt*u(inode2) 
                       if (htemp(ntemp)>highestnode) highestnode=htemp(ntemp)
                       if (htemp(ntemp)<lowestnode) lowestnode=htemp(ntemp)
                       btemp(ntemp) = b(inode2)
                       etottemp(ntemp)  = etot(inode2)
                       eratetemp(ntemp) = erate(inode2)
                       if (jk==3) then
                          inode3       = grid%icon(3,ic_k)
                          ntemp        = ntemp + 1
                          pair(ntemp)  = ntemp
                          xtemp(ntemp) = grid%x(inode3) + advect_dt*vx(inode3)
                          ytemp(ntemp) = grid%y(inode3) + advect_dt*vy(inode3)
                          htemp(ntemp) = h(inode3) + advect_dt*u(inode3) 
                          if (htemp(ntemp)>highestnode) highestnode=htemp(ntemp)
                          if (htemp(ntemp)<lowestnode) lowestnode=htemp(ntemp)
                          btemp(ntemp) = b(inode3)
                          etottemp(ntemp)  = etot(inode3)
                          eratetemp(ntemp) = erate(inode3)
                       endif
                    endif
                enddo
             enddo

             ! not possible
             if (ntemp>16) write(*,*) 'ERROR: ntemp is > 16 = ',ntemp
   
             do ii=1,ninject_per_cell(ic)
                counter_inject             = counter_inject+1
                jl                         = ceiling(real(ii)/real(nperline))
                il                         = ii - (jl-1)*nperline
                clinject%x(counter_inject) = grid%x(grid%icon(1,ic)) + il*dlx
                clinject%y(counter_inject) = grid%y(grid%icon(1,ic)) + jl*dly
   
                counter_interpolated = counter_interpolated + 1
                xip=clinject%x(counter_inject)
                yip=clinject%y(counter_inject)
                call compute_SF3 (ntemp,pair(1:ntemp),ntemp,xtemp,ytemp,xip,yip,rcut,Nnneighbours(1:ntemp)) 
                clinject%h(counter_inject)     = dot_product(Nnneighbours(1:ntemp),htemp(1:ntemp))
                clinject%b(counter_inject)     = dot_product(Nnneighbours(1:ntemp),btemp(1:ntemp))
                clinject%etot(counter_inject)  = dot_product(Nnneighbours(1:ntemp),etottemp(1:ntemp))
                clinject%erate(counter_inject) = dot_product(Nnneighbours(1:ntemp),eratetemp(1:ntemp))
 
                clinject%b(counter_inject)      = min(clinject%b(counter_inject),clinject%h(counter_inject))
                clinject%cell(counter_inject)   = ic
                clinject%active(counter_inject) = .true.
                clinject%icy(counter_inject)    = jcell
                clinject%icx(counter_inject)    = icell
                if (clinject%h(counter_inject)>highestnode) then
                    !write(*,'(a,1F10.3)') 'WARNING: something weird with elevation of new injected particle above surrounding nodes (SF3)', clinject%h(counter_inject)
                    !write(*,'(a,1F10.3,3I10)') 'highest node is  :', highestnode, ic, icell, jcell
                    clinject%h(counter_inject)     = highestnode
                    clinject%b(counter_inject)     = min(clinject%b(counter_inject),clinject%h(counter_inject))
                endif
                if (clinject%h(counter_inject)<lowestnode) then
                    !write(*,'(a,1F10.3)') 'WARNING: something weird with elevation of new injected particle below surrounding nodes (SF3)', clinject%h(counter_inject)
                    !write(*,'(a,1F10.3,3I10)') 'lowest node is  :', lowestnode, ic, icell, jcell
                    clinject%h(counter_inject)     = lowestnode
                    clinject%b(counter_inject)     = min(clinject%b(counter_inject),clinject%h(counter_inject))
                endif


             enddo !ii

             !TT we smooth all particles elevation in this cell to avoid
             !   unresolved small wavelength topographic features
             do k=1,nnic
                xip=cl%x(grid%pair(k,ic))
                yip=cl%y(grid%pair(k,ic))
                if (ntemp<=16) then
                   call compute_SF3 (ntemp,pair(1:ntemp),ntemp,xtemp,ytemp,xip,yip,rcut,Nnneighbours(1:ntemp))
                   cl%h(grid%pair(k,ic))     = dot_product(Nnneighbours(1:ntemp),htemp(1:ntemp))
                   cl%b(grid%pair(k,ic))     = dot_product(Nnneighbours(1:ntemp),btemp(1:ntemp))
                   cl%etot(grid%pair(k,ic))  = dot_product(Nnneighbours(1:ntemp),etottemp(1:ntemp))
                   cl%erate(grid%pair(k,ic)) = dot_product(Nnneighbours(1:ntemp),eratetemp(1:ntemp))
                   cl%b(grid%pair(k,ic))     = min(cl%b(grid%pair(k,ic)),cl%h(grid%pair(k,ic)))  
                endif 
                if (cl%h(grid%pair(k,ic))>highestnode) then
                    !write(*,'(a,1F10.3)') 'WARNING: something weird with elevation update of particle above surrounding nodes (SF3)', cl%h(grid%pair(k,ic))
                    !write(*,'(a,1F10.3,3I10)') 'highest node is  :', highestnode, ic, icell, jcell
                    cl%h(grid%pair(k,ic))     = highestnode
                    cl%b(grid%pair(k,ic))     = min(cl%b(grid%pair(k,ic)),cl%h(grid%pair(k,ic)))
                endif
                if (cl%h(grid%pair(k,ic))<lowestnode) then
                    !write(*,'(a,1F10.3)') 'WARNING: something weird with elevation update of particle below surrounding nodes (SF3)', cl%h(grid%pair(k,ic))
                    !write(*,'(a,1F10.3,3I10)') 'lowest node is  :', lowestnode,ic,icell,jcell
                    cl%h(grid%pair(k,ic))     = lowestnode
                    cl%b(grid%pair(k,ic))     = min(cl%b(grid%pair(k,ic)),cl%h(grid%pair(k,ic)))
                endif

             enddo

          endif !boundaries, interpolation SF3 does not work there

          grid%nn(ic) = grid%nn(ic) + ninject_per_cell(ic)
       endif    !ninject_cell>0
    enddo       !ic
    !!! omp end parallel do

    if (counter_inject /= ninject) stop 'pb inject in update_cloud_structure'
    write(*,'(a,F7.3,a,F7.3,a,F7.3,a)') 'Method for interpolation during injection: shape function = ', (float(counter_interpolated)/float(ninject))*100.d0, &
                                        ' %, average = ', (float(counter_averaged)/float(ninject))*100.d0, &
                                        ' %, closest = ', (float(counter_closest)/float(ninject))*100.d0,' %'

    deallocate(cell_injected)
  endif !if (ninject>0)
    
  !removal ---------------------------------------------------------------------------
  allocate(remove(cl%npcl))
  remove=.false.

  if (nremove>0) then

    !$omp parallel do shared(ncell,grid,cl,nremove_per_cell,remove) private(ic,nnic,distmin,ii,k1,dist,ndist,x,y,k2,k,ichoice)
    do ic=1,ncell
       
       if (nremove_per_cell(ic) > 0) then
          nnic     = grid%nn(ic)
          do ii=1,nremove_per_cell(ic)
             distmin=1.d30
             !loop through points in cell
             do k1=1,nnic
                if (.not.remove(grid%pair(k1,ic))) then
                   dist=0.d0
                   ndist=0
                   x=cl%x(grid%pair(k1,ic))
                   y=cl%y(grid%pair(k1,ic))
                   !calculate summed distance to all other points in cell that
                   !haven't been chosen to be removed yet
                   do k2=1,nnic
                      if (.not.remove(grid%pair(k2,ic)) .and. k1/=k2) then
                         dist=dist+(x-cl%x(grid%pair(k2,ic)))**2 + (y-cl%y(grid%pair(k2,ic)))**2  
                         ndist=ndist+1
                      end if
                   end do 
                   !add distances to cell corner nodes
                   do k=1,4
                      dist=dist+(x-grid%x(grid%icon(k,ic)))**2+(y-grid%y(grid%icon(k,ic)))**2
                      ndist=ndist+1
                   end do
                   !find smallest average distance to other points and corner nodes
                   dist=dist/ndist
                   if (dist<distmin) then
                      distmin=dist
                      ichoice=grid%pair(k1,ic)
                   end if
                end if
             end do
             !remove point with smallest distance
             remove(ichoice)=.true.
          end do ! do ii=1,nremove_per_cell(ic)
          grid%nn(ic) = grid%nn(ic) - nremove_per_cell(ic)
       end if    !if (nremove_per_cell(ic) > 0)

    end do
    !$omp end parallel do
    !end loop elements
   
    if(count(remove)/=nremove) stop 'pb inject in update_cloud_structure'

  endif !endif nremove>0

  allocate(not_remove(cl%npcl))
  not_remove = .not. remove

  deallocate(ninject_per_cell)
  deallocate(nremove_per_cell)


  !update arrays ---------------------------------------------------------------------------
  if (nremove>0) then
    call array_cut(cl%x                     ,not_remove)
    call array_cut(cl%y                     ,not_remove)
    call array_cut(cl%h                     ,not_remove)
    call array_cut(cl%b                     ,not_remove)
    call array_cut(cl%etot                  ,not_remove)
    call array_cut(cl%erate                 ,not_remove)
    call array_cut(cl%icy                   ,not_remove)
    call array_cut(cl%cell                  ,not_remove)
    call array_cut(cl%active                ,not_remove)
    call array_cut(cl%icx                   ,not_remove)
    call array_cut(cl%closest_node          ,not_remove)
  endif !endif nremove>0
  if (ninject>0) then
    call array_affix(cl%x                   ,clinject%x)
    call array_affix(cl%y                   ,clinject%y)
    call array_affix(cl%h                   ,clinject%h)
    call array_affix(cl%b                   ,clinject%b)
    call array_affix(cl%etot                ,clinject%etot)
    call array_affix(cl%erate               ,clinject%erate)
    call array_affix(cl%icy                 ,clinject%icy)
    call array_affix(cl%cell                ,clinject%cell)
    call array_affix(cl%active              ,clinject%active)
    call array_affix(cl%icx                 ,clinject%icx)
    call array_affix(cl%closest_node        ,clinject%closest_node)
  endif !endif ninject>0
  !free memory
  if (ninject>0) then
      deallocate(clinject%x)
      deallocate(clinject%y)
      deallocate(clinject%h)
      deallocate(clinject%b)
      deallocate(clinject%etot)
      deallocate(clinject%erate)
      deallocate(clinject%cell)
      deallocate(clinject%active)
      deallocate(clinject%icy)
      deallocate(clinject%icx)
      deallocate(clinject%closest_node)
  endif !endif ninject>0
  deallocate(not_remove)
  deallocate(remove)

  if (ninject>0 .or. nremove>0) then

  cl%npcl=npcl_new

  if (sum(grid%nn)/=cl%npcl) stop 'pb grid%nn in update_cloud (1)'
 
  if (allocated(grid%pair)) deallocate(grid%pair)
  allocate(grid%pair(maxval(grid%nn),ncell))

  grid%nn = 0
  do ip=1,cl%npcl
    ic                            = cl%cell(ip)
    grid%nn(ic)                   = grid%nn(ic) + 1
    grid%pair(grid%nn(ic),ic)     = ip
  enddo

  if (sum(grid%nn)/=cl%npcl) stop 'pb grid%nn in update_cloud (2)'

  write(*,'(a,2i15)')    'located grid%nn min/max',minval(grid%nn),maxval(grid%nn)

  endif !if (ninject>0 .or. nremove>0) then

  endif !any nnic <grid%nmin or nnic>2*grid%nmax

  return

end subroutine update_cloud

!--------------------------------------------------------------------------


!--------------------------------------------------------------------------

subroutine update_cloud_v3 (advect_dt,ierr)

  use FastScapeContext

  implicit none

  double precision, intent(in) :: advect_dt
  integer, intent(out) :: ierr
  integer ic, jcell, icell, ncellx, ncelly, ncell, nnic, npcl_max, npcl_min
  integer ninject, nremove, counter_inject
  integer ip
  integer, dimension(:), allocatable :: ninject_per_cell
  integer, dimension(:), allocatable :: nremove_per_cell
  type (cloud) :: clinject
  double precision ymax, ymin
  double precision xip, yip, xmin, xmax, distmin, dist, x, y
  logical, dimension(:), allocatable :: remove
  logical, dimension(:), allocatable :: not_remove
  logical, dimension(:), allocatable :: cell_injected
  integer ichoice, k, k2, ii, k1, ndist, npcl_new, nntot

  double precision distmin_per_quadrant(4), dx, dy, refdist
  integer ipart(4)
  integer npart, jcell_k, icell_k, ic_k, jl, il, nperline
  integer counter_interpolated, counter_averaged, counter_closest
  double precision N1, N2, N3, N4, r, s, dlx, dly

  !TT velocity dependent injection
  double precision xnode1, inode1, inode2, inode3, inode4
  double precision xnode2,xnode3,xnode4,ynode1,ynode2,ynode3,ynode4
  double precision hnode1,hnode2,hnode3,hnode4
  double precision vxinject, vyinject, vzinject, vxcl, vycl, vzcl
  double precision alpha, alphacl, alphainject, PI
  integer closest

  ierr    = 0

  if (any(grid%nn<grid%nmin) .or. any(grid%nn>2*grid%nmax)) then

  PI = 4.d0*datan(1.d0)

  ncelly  = ny-1
  ncellx  = nx-1
  ncell   = ncellx*ncelly

  dx      = xl/(nx-1)
  dy      = yl/(ny-1)
  refdist = sqrt(dx**2+dy**2)

  ninject = 0
  nremove = 0

  allocate(ninject_per_cell(ncell))
  allocate(nremove_per_cell(ncell))
  ninject_per_cell = 0
  nremove_per_cell = 0
  npcl_max         = grid%nmax
  npcl_min         = grid%nmin

  ! diagnosis
  !$omp parallel do shared(ncell,grid,ninject_per_cell,nremove_per_cell,npcl_max,npcl_min) private(ic,nnic) reduction(+:ninject,nremove)
  do ic=1,ncell
     nnic     = grid%nn(ic)
  
     if (nnic < npcl_min) then
       ninject_per_cell(ic) = floor(dble(npcl_min-nnic)/9.d0) * 9 + 9
     endif
     if (nnic > npcl_max) then
       nremove_per_cell(ic) = nnic-npcl_max
     endif
!     write(*,'(5I10)') ic, ninject_per_cell(ic), nremove_per_cell(ic), npcl_min, npcl_max

     ninject = ninject + ninject_per_cell(ic)
     nremove = nremove + nremove_per_cell(ic)

  end do
  !$omp end parallel do
  !end do

  write(*,'(a,i6,a,i6)') 'injecting ',ninject,' cloud points, maxpercell',maxval(ninject_per_cell)
  write(*,'(a,i6,a,i6)') 'removing  ',nremove,' cloud points, maxpercell',maxval(nremove_per_cell)

  npcl_new = cl%npcl + ninject-nremove

  write(*,'(a,2i9)') 'npcl new/old =',npcl_new,cl%npcl


  !injection ---------------------------------------------------------------------------
  if (ninject>0) then
    allocate(cell_injected(ncell))
    cell_injected = .false.

    allocate(clinject%x(ninject))
    allocate(clinject%y(ninject))
    allocate(clinject%h(ninject))
    allocate(clinject%b(ninject))
    allocate(clinject%etot(ninject))
    allocate(clinject%erate(ninject))
    allocate(clinject%cell(ninject))
    allocate(clinject%active(ninject))
    allocate(clinject%icy(ninject))
    allocate(clinject%icx(ninject))
    allocate(clinject%closest_node(ninject))
    clinject%closest_node = -1
    counter_inject       = 0
    counter_interpolated = 0
    counter_averaged     = 0
    counter_closest      = 0
    dly                  = dy/(3+2-1)

    !!! omp parallel do shared(ncell,ncellx,ncelly,grid,cl,ninject_per_cell,nremove_per_cell,dly,dx,clinject,refdist) private(ic,nnic,ii,il,jl,dlx,nperline,distmin_per_quadrant,ipart,k,dist,distmin,npart,icell,jcell,icell_k,jcell_k,ic_k,ip,xip,yip,xmin,xmax,r,ymin,ymax,s,N1,N2,N3,N4) reduction(+:counter_inject,counter_interpolated,counter_averaged)
    do ic=1,ncell
       jcell    = ceiling(real(ic)/real(ncellx))
       icell    = ic - (jcell-1)*ncellx
       nnic     = grid%nn(ic)
       nperline = ninject_per_cell(ic)/3
       dlx      = dx/(nperline+2-1)
       if (ninject_per_cell(ic) > 0) then
          cell_injected(ic) = .true.

          !TT
          inode1=grid%icon(1,ic)
          inode2=grid%icon(2,ic)
          inode3=grid%icon(3,ic)
          inode4=grid%icon(4,ic)
          xnode1=grid%x(inode1) + advect_dt*vx(inode1)
          ynode1=grid%y(inode1) + advect_dt*vy(inode1)
          hnode1=h(inode1) + advect_dt*u(inode1)
          xnode2=grid%x(inode2) + advect_dt*vx(inode2)
          ynode2=grid%y(inode2) + advect_dt*vy(inode2)
          hnode2=h(inode2) + advect_dt*u(inode2)
          xnode3=grid%x(inode3) + advect_dt*vx(inode3)
          ynode3=grid%y(inode3) + advect_dt*vy(inode3)
          hnode3=h(inode3) + advect_dt*u(inode3)
          xnode4=grid%x(inode4) + advect_dt*vx(inode4)
          ynode4=grid%y(inode4) + advect_dt*vy(inode4)
          hnode4=h(inode4) + advect_dt*u(inode4)
 
          do ii=1,ninject_per_cell(ic)
             counter_inject             = counter_inject+1
             jl                         = ceiling(real(ii)/real(nperline))
             il                         = ii - (jl-1)*nperline
             clinject%x(counter_inject) = grid%x(grid%icon(1,ic)) + il*dlx
             clinject%y(counter_inject) = grid%y(grid%icon(1,ic)) + jl*dly
 
!             !find closest cloud point to new point
!             if (nnic>0) then
!                distmin_per_quadrant = 1.d30
!                distmin   = 1.d30
!                ipart     = -1
!                do k=1,nnic
!
!                   !TT velocity dependent injection
!                   xip=cl%x(grid%pair(k,ic))
!                   yip=cl%y(grid%pair(k,ic))
!                   r=((xip-xnode1)/dx -0.5d0 ) *2.d0
!                   ymin=(grid%y(inode2)-grid%y(inode1))/dx*(xip-xnode1) + grid%y(inode1)
!                   ymax=(grid%y(inode3)-grid%y(inode4))/dx*(xip-xnode1) + grid%y(inode4)
!                   s=((yip-ymin)/(ymax-ymin) -0.5d0 ) *2.d0
!                   N1=0.25d0*(1.d0-r)*(1.d0-s) 
!                   N2=0.25d0*(1.d0+r)*(1.d0-s) 
!                   N3=0.25d0*(1.d0+r)*(1.d0+s) 
!                   N4=0.25d0*(1.d0-r)*(1.d0+s) 
!                   vxcl = N1 * vx(inode1) + N2 * vx(inode2) + N3 * vx(inode3) + N4 * vx(inode4) 
!                   vycl = N1 * vy(inode1) + N2 * vy(inode2) + N3 * vy(inode3) + N4 * vy(inode4)
!                   !vzcl = N1 * u(inode1)  + N2 * u(inode2)  + N3 * u(inode3)  + N4 * u(inode4)
!                   alphacl = atan(dabs(vycl/vxcl))
!                   if (vycl>0 .and. vxcl<0) then
!                      alphacl = alphacl + PI/2
!                   elseif (vycl<0 .and. vxcl<0) then
!                      alphacl = alphacl + PI
!                   elseif (vycl<0 .and. vxcl>0) then
!                      alphacl = alphacl + 3*PI/2
!                   endif
!                   alpha = dabs(alphainject - alphacl)
!
!                   dist=(clinject%x(counter_inject)-cl%x(grid%pair(k,ic)))**2 + &
!                        (clinject%y(counter_inject)-cl%y(grid%pair(k,ic)))**2 
!                   if (dist < distmin) then
!                        distmin = dist
!                        closest = grid%pair(k,ic)
!                   endif
!                   if ( cl%x(grid%pair(k,ic)) < clinject%x(counter_inject) .and. &
!                        cl%y(grid%pair(k,ic)) < clinject%y(counter_inject) .and. &
!                        dist < distmin_per_quadrant(1) .and. &
!                        ( alpha < 20 .or. alpha > 340 ) ) then
!                      distmin_per_quadrant(1) = dist
!                      ipart(1)                = grid%pair(k,ic)
!                      cycle
!                   endif
!                   if ( cl%x(grid%pair(k,ic)) > clinject%x(counter_inject) .and. &
!                        cl%y(grid%pair(k,ic)) < clinject%y(counter_inject) .and. &
!                        dist < distmin_per_quadrant(2) .and. &
!                        ( alpha < 20 .or. alpha > 340 ) ) then
!                      distmin_per_quadrant(2) = dist
!                      ipart(2)                = grid%pair(k,ic)
!                      cycle
!                   endif
!                   if ( cl%x(grid%pair(k,ic)) > clinject%x(counter_inject) .and. &
!                        cl%y(grid%pair(k,ic)) > clinject%y(counter_inject) .and. &
!                        dist < distmin_per_quadrant(3) .and. &
!                        ( alpha < 20 .or. alpha > 340 ) ) then
!                      distmin_per_quadrant(3) = dist
!                      ipart(3)                = grid%pair(k,ic)
!                      cycle
!                   endif
!                   if ( cl%x(grid%pair(k,ic)) < clinject%x(counter_inject) .and. &
!                        cl%y(grid%pair(k,ic)) > clinject%y(counter_inject) .and. &
!                        dist < distmin_per_quadrant(4) .and. &
!                        ( alpha < 20 .or. alpha > 340 ) ) then
!                      distmin_per_quadrant(4) = dist
!                      ipart(4)                = grid%pair(k,ic)
!                      cycle
!                   endif
!                end do
!                npart = 0
!                do k=1,4
!                   if (ipart(k)>0) npart = npart + 1
!                enddo
!             else
!                ipart   = -1
!                npart   = 0
!                write(*,*) "WARNING: nnic = 0 (update_cloud)"
!                !stop 'nnic = 0 | not possible'
!             endif !nnic>0
! 
!             if (npart<4) then
!               ! if no particle in a quadrant we find the closest particle in the
!               ! adjacent cell
!               do k=1,4
!                 if (ipart(k)<0) then
!                   if (k==1) then
!                     icell_k = icell - 1
!                     jcell_k = jcell - 1
!                     if ((icell_k<1).or.(jcell_k<1)) cycle
!                   else if (k==2) then
!                     icell_k = icell + 1
!                     jcell_k = jcell - 1
!                     if ((icell_k>ncellx).or.(jcell_k<1)) cycle
!                   else if (k==3) then
!                     icell_k = icell + 1
!                     jcell_k = jcell + 1
!                     if ((icell_k>ncellx).or.(jcell_k>ncelly)) cycle
!                   else if (k==4) then
!                     icell_k = icell - 1
!                     jcell_k = jcell + 1
!                     if ((icell_k<1).or.(jcell_k>ncelly)) cycle
!                   endif
!                   ic_k      = (jcell_k-1)*ncellx+icell_k
!                   distmin   = 1.d30
!                   if (cell_injected(ic_k)) then
!                     nntot = grid%nn(ic_k) - ninject_per_cell(ic_k)
!                   else
!                     nntot = grid%nn(ic_k)
!                   endif
!                   if (nntot==0) then
!                     if (k==1 .or. k==4) then
!                       do while (icell_k>0 .and. nntot <= 0) 
!                         icell_k = icell_k - 1
!                         ic_k     = (jcell_k-1)*ncellx+icell_k
!                         if (cell_injected(ic_k)) then
!                           nntot = grid%nn(ic_k) - ninject_per_cell(ic_k)
!                         else
!                           nntot = grid%nn(ic_k)
!                         endif
!                       enddo
!                     endif
!                     if (k==2 .or. k==3) then
!                       do while (icell_k<=ncellx .and. nntot <= 0)
!                         icell_k = icell_k + 1
!                         ic_k     = (jcell_k-1)*ncellx+icell_k
!                         if (cell_injected(ic_k)) then
!                           nntot = grid%nn(ic_k) - ninject_per_cell(ic_k)
!                         else
!                           nntot = grid%nn(ic_k)
!                         endif
!                       enddo
!                     endif
!                   endif
!                   !write(*,'(6I10)') icell,jcell,k,ic_k, nntot,ninject_per_cell(ic_k)
!                   do ip=1,nntot
!                     !we cannot use new injected
!                     !if (grid%pair(ip,ic_k)>0.and.grid%pair(ip,ic_k)<=cl%npcl) then
!                       !write(*,*) ip,ic_k,grid%pair(ip,ic_k)
!                       dist=(clinject%x(counter_inject)-cl%x(grid%pair(ip,ic_k)))**2 + &
!                            (clinject%y(counter_inject)-cl%y(grid%pair(ip,ic_k)))**2
!                       if (dist<distmin) then
!                          distmin   = dist
!                          ipart(k)  = grid%pair(ip,ic_k)
!                       end if
!                     !endif
!                   enddo
!                   if (nntot>0) then
!                     npart = npart + 1
!                   else
!                     write(*,'(a,I10,a,I10,a,I10,a,I10)') 'WARNING: The cell ic=',ic,' icell=',icell,' and jcell=',jcell,' did not find closest particle in quadrant ', k
!                   endif
!                 endif
!               enddo
!             endif
 
              counter_interpolated = counter_interpolated + 1
              xip=clinject%x(counter_inject)
              yip=clinject%y(counter_inject)
 
              xmin=(xnode4-xnode1)/(ynode4-ynode1)*(yip-ynode1) + xnode1
              xmax=(xnode3-xnode2)/(ynode3-ynode2)*(yip-ynode2) + xnode2
              r=((xip-xmin)/(xmax-xmin) -0.5d0 ) *2.d0
 
              ymin=(ynode2-ynode1)/(xnode2-xnode1)*(xip-xnode1) + ynode1
              ymax=(ynode3-ynode4)/(xnode3-xnode4)*(xip-xnode4) + ynode4
              s=((yip-ymin)/(ymax-ymin) -0.5d0 ) *2.d0
              N1=0.25d0*(1.d0-r)*(1.d0-s) 
              N2=0.25d0*(1.d0+r)*(1.d0-s) 
              N3=0.25d0*(1.d0+r)*(1.d0+s) 
              N4=0.25d0*(1.d0-r)*(1.d0+s) 
              clinject%h(counter_inject)     = N1 * hnode1     + N2 * hnode2     + N3 * hnode3     + N4 * hnode4
              clinject%b(counter_inject)     = N1 * b(grid%icon(1,ic))     + N2 * b(grid%icon(2,ic))     + N3 * b(grid%icon(3,ic))     + N4 * b(grid%icon(4,ic))
              clinject%etot(counter_inject)  = N1 * etot(grid%icon(1,ic))  + N2 * etot(grid%icon(2,ic))  + N3 * etot(grid%icon(3,ic))  + N4 * etot(grid%icon(4,ic))
              clinject%erate(counter_inject) = N1 * erate(grid%icon(1,ic)) + N2 * erate(grid%icon(2,ic)) + N3 * erate(grid%icon(3,ic)) + N4 * erate(grid%icon(4,ic))

!             if (npart==4) then
!                counter_interpolated = counter_interpolated + 1
!                xip=clinject%x(counter_inject)
!                yip=clinject%y(counter_inject)
! 
!                xmin=(cl%x(ipart(4))-cl%x(ipart(1)))/(cl%y(ipart(4))-cl%y(ipart(1)))*(yip-cl%y(ipart(1))) + cl%x(ipart(1))
!                xmax=(cl%x(ipart(3))-cl%x(ipart(2)))/(cl%y(ipart(3))-cl%y(ipart(2)))*(yip-cl%y(ipart(2))) + cl%x(ipart(2))
!                r=((xip-xmin)/(xmax-xmin) -0.5d0 ) *2.d0
! 
!                ymin=(cl%y(ipart(2))-cl%y(ipart(1)))/(cl%x(ipart(2))-cl%x(ipart(1)))*(xip-cl%x(ipart(1))) + cl%y(ipart(1))
!                ymax=(cl%y(ipart(3))-cl%y(ipart(4)))/(cl%x(ipart(3))-cl%x(ipart(4)))*(xip-cl%x(ipart(4))) + cl%y(ipart(4))
!                s=((yip-ymin)/(ymax-ymin) -0.5d0 ) *2.d0
!                N1=0.25d0*(1.d0-r)*(1.d0-s) 
!                N2=0.25d0*(1.d0+r)*(1.d0-s) 
!                N3=0.25d0*(1.d0+r)*(1.d0+s) 
!                N4=0.25d0*(1.d0-r)*(1.d0+s) 
!                clinject%h(counter_inject)     = N1 * cl%h(ipart(1))     + N2 * cl%h(ipart(2))     + N3 * cl%h(ipart(3))     + N4 * cl%h(ipart(4))
!                clinject%b(counter_inject)     = N1 * cl%b(ipart(1))     + N2 * cl%b(ipart(2))     + N3 * cl%b(ipart(3))     + N4 * cl%b(ipart(4))
!                clinject%etot(counter_inject)  = N1 * cl%etot(ipart(1))  + N2 * cl%etot(ipart(2))  + N3 * cl%etot(ipart(3))  + N4 * cl%etot(ipart(4))
!                clinject%erate(counter_inject) = N1 * cl%erate(ipart(1)) + N2 * cl%erate(ipart(2)) + N3 * cl%erate(ipart(3)) + N4 * cl%erate(ipart(4))
!             else
!                counter_averaged = counter_averaged + 1
!                clinject%h(counter_inject)     = 0.d0
!                clinject%b(counter_inject)     = 0.d0
!                clinject%etot(counter_inject)  = 0.d0
!                clinject%erate(counter_inject) = 0.d0
!                N1 = 0.d0
!                if (npart>0) then
!                  do k=1,4
!                     if (ipart(k)>0) then
!                       clinject%h(counter_inject) = clinject%h(counter_inject) + cl%h(ipart(k))*refdist/distmin_per_quadrant(k)
!                       clinject%b(counter_inject) = clinject%b(counter_inject) + cl%b(ipart(k))*refdist/distmin_per_quadrant(k)
!                       clinject%etot(counter_inject) = clinject%etot(counter_inject) + cl%etot(ipart(k))*refdist/distmin_per_quadrant(k)
!                       clinject%erate(counter_inject) = clinject%erate(counter_inject) + cl%erate(ipart(k))*refdist/distmin_per_quadrant(k)
!                       N1 = N1 + (refdist/distmin_per_quadrant(k))
!                     endif
!                  enddo
!                  clinject%h(counter_inject)     = clinject%h(counter_inject)/N1
!                  clinject%b(counter_inject)     = clinject%b(counter_inject)/N1
!                  clinject%etot(counter_inject)  = clinject%etot(counter_inject)/N1
!                  clinject%erate(counter_inject) = clinject%erate(counter_inject)/N1
!                else
!                  counter_closest = counter_closest + 1
!                  write(*,'(a,I10,a,I10,a,I10,a,I10)') 'WARNING: The cell ic=',ic,' icell=',icell,' and jcell=',jcell,' did not find any closest particle with velocity conditions'
!                  clinject%h(counter_inject)     = cl%h(closest)
!                  clinject%b(counter_inject)     = cl%b(closest)
!                  clinject%etot(counter_inject)  = cl%etot(closest)
!                  clinject%erate(counter_inject) = cl%erate(closest)
!                endif
! 
             !endif
 
             clinject%b(counter_inject)      = min(clinject%b(counter_inject),clinject%h(counter_inject))
             clinject%cell(counter_inject)   = ic
             clinject%active(counter_inject) = .true.
             clinject%icy(counter_inject)    = jcell
             clinject%icx(counter_inject)    = icell

          enddo !ii
          grid%nn(ic) = grid%nn(ic) + ninject_per_cell(ic)
       endif    !ninject_cell>0
    enddo       !ic
    !!! omp end parallel do

    if (counter_inject /= ninject) stop 'pb inject in update_cloud_structure'
    write(*,'(a,F7.3,a,F7.3,a,F7.3,a)') 'Method for interpolation during injection: shape function = ', (float(counter_interpolated)/float(ninject))*100.d0, &
                                        ' %, average = ', (float(counter_averaged)/float(ninject))*100.d0, &
                                        ' %, closest = ', (float(counter_closest)/float(ninject))*100.d0,' %'

    deallocate(cell_injected)
  endif !if (ninject>0)
    
  !removal ---------------------------------------------------------------------------
  allocate(remove(cl%npcl))
  remove=.false.

  if (nremove>0) then

    !$omp parallel do shared(ncell,grid,cl,nremove_per_cell,remove) private(ic,nnic,distmin,ii,k1,dist,ndist,x,y,k2,k,ichoice)
    do ic=1,ncell
       
       if (nremove_per_cell(ic) > 0) then
          nnic     = grid%nn(ic)
          do ii=1,nremove_per_cell(ic)
             distmin=1.d30
             !loop through points in cell
             do k1=1,nnic
                if (.not.remove(grid%pair(k1,ic))) then
                   dist=0.d0
                   ndist=0
                   x=cl%x(grid%pair(k1,ic))
                   y=cl%y(grid%pair(k1,ic))
                   !calculate summed distance to all other points in cell that
                   !haven't been chosen to be removed yet
                   do k2=1,nnic
                      if (.not.remove(grid%pair(k2,ic)) .and. k1/=k2) then
                         dist=dist+(x-cl%x(grid%pair(k2,ic)))**2 + (y-cl%y(grid%pair(k2,ic)))**2  
                         ndist=ndist+1
                      end if
                   end do 
                   !add distances to cell corner nodes
                   do k=1,4
                      dist=dist+(x-grid%x(grid%icon(k,ic)))**2+(y-grid%y(grid%icon(k,ic)))**2
                      ndist=ndist+1
                   end do
                   !find smallest average distance to other points and corner nodes
                   dist=dist/ndist
                   if (dist<distmin) then
                      distmin=dist
                      ichoice=grid%pair(k1,ic)
                   end if
                end if
             end do
             !remove point with smallest distance
             remove(ichoice)=.true.
          end do ! do ii=1,nremove_per_cell(ic)
          grid%nn(ic) = grid%nn(ic) - nremove_per_cell(ic)
       end if    !if (nremove_per_cell(ic) > 0)

    end do
    !$omp end parallel do
    !end loop elements
   
    if(count(remove)/=nremove) stop 'pb inject in update_cloud_structure'

  endif !endif nremove>0

  allocate(not_remove(cl%npcl))
  not_remove = .not. remove

  deallocate(ninject_per_cell)
  deallocate(nremove_per_cell)


  !update arrays ---------------------------------------------------------------------------
  if (nremove>0) then
    call array_cut(cl%x                     ,not_remove)
    call array_cut(cl%y                     ,not_remove)
    call array_cut(cl%h                     ,not_remove)
    call array_cut(cl%b                     ,not_remove)
    call array_cut(cl%etot                  ,not_remove)
    call array_cut(cl%erate                 ,not_remove)
    call array_cut(cl%icy                   ,not_remove)
    call array_cut(cl%cell                  ,not_remove)
    call array_cut(cl%active                ,not_remove)
    call array_cut(cl%icx                   ,not_remove)
    call array_cut(cl%closest_node          ,not_remove)
  endif !endif nremove>0
  if (ninject>0) then
    call array_affix(cl%x                   ,clinject%x)
    call array_affix(cl%y                   ,clinject%y)
    call array_affix(cl%h                   ,clinject%h)
    call array_affix(cl%b                   ,clinject%b)
    call array_affix(cl%etot                ,clinject%etot)
    call array_affix(cl%erate               ,clinject%erate)
    call array_affix(cl%icy                 ,clinject%icy)
    call array_affix(cl%cell                ,clinject%cell)
    call array_affix(cl%active              ,clinject%active)
    call array_affix(cl%icx                 ,clinject%icx)
    call array_affix(cl%closest_node        ,clinject%closest_node)
  endif !endif ninject>0
  !free memory
  if (ninject>0) then
      deallocate(clinject%x)
      deallocate(clinject%y)
      deallocate(clinject%h)
      deallocate(clinject%b)
      deallocate(clinject%etot)
      deallocate(clinject%erate)
      deallocate(clinject%cell)
      deallocate(clinject%active)
      deallocate(clinject%icy)
      deallocate(clinject%icx)
      deallocate(clinject%closest_node)
  endif !endif ninject>0
  deallocate(not_remove)
  deallocate(remove)

  if (ninject>0 .or. nremove>0) then

  cl%npcl=npcl_new

  if (sum(grid%nn)/=cl%npcl) stop 'pb grid%nn in update_cloud (1)'
 
  if (allocated(grid%pair)) deallocate(grid%pair)
  allocate(grid%pair(maxval(grid%nn),ncell))

  grid%nn = 0
  do ip=1,cl%npcl
    ic                            = cl%cell(ip)
    grid%nn(ic)                   = grid%nn(ic) + 1
    grid%pair(grid%nn(ic),ic)     = ip
  enddo

  if (sum(grid%nn)/=cl%npcl) stop 'pb grid%nn in update_cloud (2)'

  write(*,'(a,2i15)')    'located grid%nn min/max',minval(grid%nn),maxval(grid%nn)

  endif !if (ninject>0 .or. nremove>0) then

  endif !any nnic <grid%nmin or nnic>2*grid%nmax

  return

end subroutine update_cloud_v3

!--------------------------------------------------------------------------


!--------------------------------------------------------------------------

subroutine update_cloud_v2 (ierr)

  use FastScapeContext

  implicit none

  integer, intent(out) :: ierr
  integer ic, jcell, icell, ncellx, ncelly, ncell, nnic, npcl_max, npcl_min
  integer ninject, nremove, counter_inject
  integer ip
  integer, dimension(:), allocatable :: ninject_per_cell
  integer, dimension(:), allocatable :: nremove_per_cell
  type (cloud) :: clinject
  double precision ymax, ymin
  double precision xip, yip, xmin, xmax, distmin, dist, x, y
  logical, dimension(:), allocatable :: remove
  logical, dimension(:), allocatable :: not_remove
  logical, dimension(:), allocatable :: cell_injected
  integer ichoice, k, k2, ii, k1, ndist, npcl_new, nntot

  double precision distmin_per_quadrant(4), dx, dy, refdist
  integer ipart(4)
  integer npart, jcell_k, icell_k, ic_k, jl, il, nperline
  integer counter_interpolated, counter_averaged, counter_closest
  double precision N1, N2, N3, N4, r, s, dlx, dly

  !TT velocity dependent injection
  double precision xnode1, inode1, inode2, inode3, inode4
  double precision vxinject, vyinject, vzinject, vxcl, vycl, vzcl
  double precision alpha, alphacl, alphainject, PI
  integer closest

  ierr    = 0

  if (any(grid%nn<grid%nmin) .or. any(grid%nn>2*grid%nmax)) then

  PI = 4.d0*datan(1.d0)

  ncelly  = ny-1
  ncellx  = nx-1
  ncell   = ncellx*ncelly

  dx      = xl/(nx-1)
  dy      = yl/(ny-1)
  refdist = sqrt(dx**2+dy**2)

  ninject = 0
  nremove = 0

  allocate(ninject_per_cell(ncell))
  allocate(nremove_per_cell(ncell))
  ninject_per_cell = 0
  nremove_per_cell = 0
  npcl_max         = grid%nmax
  npcl_min         = grid%nmin

  ! diagnosis
  !$omp parallel do shared(ncell,grid,ninject_per_cell,nremove_per_cell,npcl_max,npcl_min) private(ic,nnic) reduction(+:ninject,nremove)
  do ic=1,ncell
     nnic     = grid%nn(ic)
  
     if (nnic < npcl_min) then
       ninject_per_cell(ic) = floor(dble(npcl_min-nnic)/9.d0) * 9 + 9
     endif
     if (nnic > npcl_max) then
       nremove_per_cell(ic) = nnic-npcl_max
     endif
!     write(*,'(5I10)') ic, ninject_per_cell(ic), nremove_per_cell(ic), npcl_min, npcl_max

     ninject = ninject + ninject_per_cell(ic)
     nremove = nremove + nremove_per_cell(ic)

  end do
  !$omp end parallel do
  !end do

  write(*,'(a,i6,a,i6)') 'injecting ',ninject,' cloud points, maxpercell',maxval(ninject_per_cell)
  write(*,'(a,i6,a,i6)') 'removing  ',nremove,' cloud points, maxpercell',maxval(nremove_per_cell)

  npcl_new = cl%npcl + ninject-nremove

  write(*,'(a,2i9)') 'npcl new/old =',npcl_new,cl%npcl


  !injection ---------------------------------------------------------------------------
  if (ninject>0) then
    allocate(cell_injected(ncell))
    cell_injected = .false.

    allocate(clinject%x(ninject))
    allocate(clinject%y(ninject))
    allocate(clinject%h(ninject))
    allocate(clinject%b(ninject))
    allocate(clinject%etot(ninject))
    allocate(clinject%erate(ninject))
    allocate(clinject%cell(ninject))
    allocate(clinject%active(ninject))
    allocate(clinject%icy(ninject))
    allocate(clinject%icx(ninject))
    allocate(clinject%closest_node(ninject))
    clinject%closest_node = -1
    counter_inject       = 0
    counter_interpolated = 0
    counter_averaged     = 0
    counter_closest      = 0
    dly                  = dy/(3+2-1)

    !!! omp parallel do shared(ncell,ncellx,ncelly,grid,cl,ninject_per_cell,nremove_per_cell,dly,dx,clinject,refdist) private(ic,nnic,ii,il,jl,dlx,nperline,distmin_per_quadrant,ipart,k,dist,distmin,npart,icell,jcell,icell_k,jcell_k,ic_k,ip,xip,yip,xmin,xmax,r,ymin,ymax,s,N1,N2,N3,N4) reduction(+:counter_inject,counter_interpolated,counter_averaged)
    do ic=1,ncell
       jcell    = ceiling(real(ic)/real(ncellx))
       icell    = ic - (jcell-1)*ncellx
       nnic     = grid%nn(ic)
       nperline = ninject_per_cell(ic)/3
       dlx      = dx/(nperline+2-1)
       if (ninject_per_cell(ic) > 0) then
          cell_injected(ic) = .true.

          !TT velocity dependent injection
          inode1=grid%icon(1,ic)
          inode2=grid%icon(2,ic)
          inode3=grid%icon(3,ic)
          inode4=grid%icon(4,ic)
          xnode1=grid%x(inode1)

          do ii=1,ninject_per_cell(ic)
             counter_inject             = counter_inject+1
             jl                         = ceiling(real(ii)/real(nperline))
             il                         = ii - (jl-1)*nperline
             clinject%x(counter_inject) = grid%x(grid%icon(1,ic)) + il*dlx
             clinject%y(counter_inject) = grid%y(grid%icon(1,ic)) + jl*dly
 
             !TT velocity dependent injection
             xip=clinject%x(counter_inject)
             yip=clinject%y(counter_inject)
             r=((xip-xnode1)/dx -0.5d0 ) *2.d0
             ymin=(grid%y(inode2)-grid%y(inode1))/dx*(xip-xnode1) + grid%y(inode1)
             ymax=(grid%y(inode3)-grid%y(inode4))/dx*(xip-xnode1) + grid%y(inode4)
             s=((yip-ymin)/(ymax-ymin) -0.5d0 ) *2.d0
             N1=0.25d0*(1.d0-r)*(1.d0-s) 
             N2=0.25d0*(1.d0+r)*(1.d0-s) 
             N3=0.25d0*(1.d0+r)*(1.d0+s) 
             N4=0.25d0*(1.d0-r)*(1.d0+s) 
             vxinject = N1 * vx(inode1) + N2 * vx(inode2) + N3 * vx(inode3) + N4 * vx(inode4) 
             vyinject = N1 * vy(inode1) + N2 * vy(inode2) + N3 * vy(inode3) + N4 * vy(inode4)
             !vzinject = N1 * u(inode1)  + N2 * u(inode2)  + N3 * u(inode3)  + N4 * u(inode4)
             alphainject = atan(dabs(vyinject/vxinject))
             if (vyinject>0 .and. vxinject<0) then
                alphainject = alphainject + PI/2
             elseif (vyinject<0 .and. vxinject<0) then
                alphainject = alphainject + PI
             elseif (vyinject<0 .and. vxinject>0) then
                alphainject = alphainject + 3*PI/2
             endif

             !find closest cloud point to new point
             if (nnic>0) then
                distmin_per_quadrant = 1.d30
                distmin   = 1.d30
                ipart     = -1
                do k=1,nnic

                   !TT velocity dependent injection
                   xip=cl%x(grid%pair(k,ic))
                   yip=cl%y(grid%pair(k,ic))
                   r=((xip-xnode1)/dx -0.5d0 ) *2.d0
                   ymin=(grid%y(inode2)-grid%y(inode1))/dx*(xip-xnode1) + grid%y(inode1)
                   ymax=(grid%y(inode3)-grid%y(inode4))/dx*(xip-xnode1) + grid%y(inode4)
                   s=((yip-ymin)/(ymax-ymin) -0.5d0 ) *2.d0
                   N1=0.25d0*(1.d0-r)*(1.d0-s) 
                   N2=0.25d0*(1.d0+r)*(1.d0-s) 
                   N3=0.25d0*(1.d0+r)*(1.d0+s) 
                   N4=0.25d0*(1.d0-r)*(1.d0+s) 
                   vxcl = N1 * vx(inode1) + N2 * vx(inode2) + N3 * vx(inode3) + N4 * vx(inode4) 
                   vycl = N1 * vy(inode1) + N2 * vy(inode2) + N3 * vy(inode3) + N4 * vy(inode4)
                   !vzcl = N1 * u(inode1)  + N2 * u(inode2)  + N3 * u(inode3)  + N4 * u(inode4)
                   alphacl = atan(dabs(vycl/vxcl))
                   if (vycl>0 .and. vxcl<0) then
                      alphacl = alphacl + PI/2
                   elseif (vycl<0 .and. vxcl<0) then
                      alphacl = alphacl + PI
                   elseif (vycl<0 .and. vxcl>0) then
                      alphacl = alphacl + 3*PI/2
                   endif
                   alpha = dabs(alphainject - alphacl)

                   dist=(clinject%x(counter_inject)-cl%x(grid%pair(k,ic)))**2 + &
                        (clinject%y(counter_inject)-cl%y(grid%pair(k,ic)))**2 
                   if (dist < distmin) then
                        distmin = dist
                        closest = grid%pair(k,ic)
                   endif
                   if ( cl%x(grid%pair(k,ic)) < clinject%x(counter_inject) .and. &
                        cl%y(grid%pair(k,ic)) < clinject%y(counter_inject) .and. &
                        dist < distmin_per_quadrant(1) .and. &
                        ( alpha < 20 .or. alpha > 340 ) ) then
                      distmin_per_quadrant(1) = dist
                      ipart(1)                = grid%pair(k,ic)
                      cycle
                   endif
                   if ( cl%x(grid%pair(k,ic)) > clinject%x(counter_inject) .and. &
                        cl%y(grid%pair(k,ic)) < clinject%y(counter_inject) .and. &
                        dist < distmin_per_quadrant(2) .and. &
                        ( alpha < 20 .or. alpha > 340 ) ) then
                      distmin_per_quadrant(2) = dist
                      ipart(2)                = grid%pair(k,ic)
                      cycle
                   endif
                   if ( cl%x(grid%pair(k,ic)) > clinject%x(counter_inject) .and. &
                        cl%y(grid%pair(k,ic)) > clinject%y(counter_inject) .and. &
                        dist < distmin_per_quadrant(3) .and. &
                        ( alpha < 20 .or. alpha > 340 ) ) then
                      distmin_per_quadrant(3) = dist
                      ipart(3)                = grid%pair(k,ic)
                      cycle
                   endif
                   if ( cl%x(grid%pair(k,ic)) < clinject%x(counter_inject) .and. &
                        cl%y(grid%pair(k,ic)) > clinject%y(counter_inject) .and. &
                        dist < distmin_per_quadrant(4) .and. &
                        ( alpha < 20 .or. alpha > 340 ) ) then
                      distmin_per_quadrant(4) = dist
                      ipart(4)                = grid%pair(k,ic)
                      cycle
                   endif
                end do
                npart = 0
                do k=1,4
                   if (ipart(k)>0) npart = npart + 1
                enddo
             else
                ipart   = -1
                npart   = 0
                write(*,*) "WARNING: nnic = 0 (update_cloud)"
                !stop 'nnic = 0 | not possible'
             endif !nnic>0
 
!             if (npart<4) then
!               ! if no particle in a quadrant we find the closest particle in the
!               ! adjacent cell
!               do k=1,4
!                 if (ipart(k)<0) then
!                   if (k==1) then
!                     icell_k = icell - 1
!                     jcell_k = jcell - 1
!                     if ((icell_k<1).or.(jcell_k<1)) cycle
!                   else if (k==2) then
!                     icell_k = icell + 1
!                     jcell_k = jcell - 1
!                     if ((icell_k>ncellx).or.(jcell_k<1)) cycle
!                   else if (k==3) then
!                     icell_k = icell + 1
!                     jcell_k = jcell + 1
!                     if ((icell_k>ncellx).or.(jcell_k>ncelly)) cycle
!                   else if (k==4) then
!                     icell_k = icell - 1
!                     jcell_k = jcell + 1
!                     if ((icell_k<1).or.(jcell_k>ncelly)) cycle
!                   endif
!                   ic_k      = (jcell_k-1)*ncellx+icell_k
!                   distmin   = 1.d30
!                   if (cell_injected(ic_k)) then
!                     nntot = grid%nn(ic_k) - ninject_per_cell(ic_k)
!                   else
!                     nntot = grid%nn(ic_k)
!                   endif
!                   if (nntot==0) then
!                     if (k==1 .or. k==4) then
!                       do while (icell_k>0 .and. nntot <= 0) 
!                         icell_k = icell_k - 1
!                         ic_k     = (jcell_k-1)*ncellx+icell_k
!                         if (cell_injected(ic_k)) then
!                           nntot = grid%nn(ic_k) - ninject_per_cell(ic_k)
!                         else
!                           nntot = grid%nn(ic_k)
!                         endif
!                       enddo
!                     endif
!                     if (k==2 .or. k==3) then
!                       do while (icell_k<=ncellx .and. nntot <= 0)
!                         icell_k = icell_k + 1
!                         ic_k     = (jcell_k-1)*ncellx+icell_k
!                         if (cell_injected(ic_k)) then
!                           nntot = grid%nn(ic_k) - ninject_per_cell(ic_k)
!                         else
!                           nntot = grid%nn(ic_k)
!                         endif
!                       enddo
!                     endif
!                   endif
!                   !write(*,'(6I10)') icell,jcell,k,ic_k, nntot,ninject_per_cell(ic_k)
!                   do ip=1,nntot
!                     !we cannot use new injected
!                     !if (grid%pair(ip,ic_k)>0.and.grid%pair(ip,ic_k)<=cl%npcl) then
!                       !write(*,*) ip,ic_k,grid%pair(ip,ic_k)
!                       dist=(clinject%x(counter_inject)-cl%x(grid%pair(ip,ic_k)))**2 + &
!                            (clinject%y(counter_inject)-cl%y(grid%pair(ip,ic_k)))**2
!                       if (dist<distmin) then
!                          distmin   = dist
!                          ipart(k)  = grid%pair(ip,ic_k)
!                       end if
!                     !endif
!                   enddo
!                   if (nntot>0) then
!                     npart = npart + 1
!                   else
!                     write(*,'(a,I10,a,I10,a,I10,a,I10)') 'WARNING: The cell ic=',ic,' icell=',icell,' and jcell=',jcell,' did not find closest particle in quadrant ', k
!                   endif
!                 endif
!               enddo
!             endif
 
!             if (npart==4) then
!                counter_interpolated = counter_interpolated + 1
!                xip=clinject%x(counter_inject)
!                yip=clinject%y(counter_inject)
! 
!                xmin=(cl%x(ipart(4))-cl%x(ipart(1)))/(cl%y(ipart(4))-cl%y(ipart(1)))*(yip-cl%y(ipart(1))) + cl%x(ipart(1))
!                xmax=(cl%x(ipart(3))-cl%x(ipart(2)))/(cl%y(ipart(3))-cl%y(ipart(2)))*(yip-cl%y(ipart(2))) + cl%x(ipart(2))
!                r=((xip-xmin)/(xmax-xmin) -0.5d0 ) *2.d0
! 
!                ymin=(cl%y(ipart(2))-cl%y(ipart(1)))/(cl%x(ipart(2))-cl%x(ipart(1)))*(xip-cl%x(ipart(1))) + cl%y(ipart(1))
!                ymax=(cl%y(ipart(3))-cl%y(ipart(4)))/(cl%x(ipart(3))-cl%x(ipart(4)))*(xip-cl%x(ipart(4))) + cl%y(ipart(4))
!                s=((yip-ymin)/(ymax-ymin) -0.5d0 ) *2.d0
!                N1=0.25d0*(1.d0-r)*(1.d0-s) 
!                N2=0.25d0*(1.d0+r)*(1.d0-s) 
!                N3=0.25d0*(1.d0+r)*(1.d0+s) 
!                N4=0.25d0*(1.d0-r)*(1.d0+s) 
!                clinject%h(counter_inject)     = N1 * cl%h(ipart(1))     + N2 * cl%h(ipart(2))     + N3 * cl%h(ipart(3))     + N4 * cl%h(ipart(4))
!                clinject%b(counter_inject)     = N1 * cl%b(ipart(1))     + N2 * cl%b(ipart(2))     + N3 * cl%b(ipart(3))     + N4 * cl%b(ipart(4))
!                clinject%etot(counter_inject)  = N1 * cl%etot(ipart(1))  + N2 * cl%etot(ipart(2))  + N3 * cl%etot(ipart(3))  + N4 * cl%etot(ipart(4))
!                clinject%erate(counter_inject) = N1 * cl%erate(ipart(1)) + N2 * cl%erate(ipart(2)) + N3 * cl%erate(ipart(3)) + N4 * cl%erate(ipart(4))
!             else
                counter_averaged = counter_averaged + 1
                clinject%h(counter_inject)     = 0.d0
                clinject%b(counter_inject)     = 0.d0
                clinject%etot(counter_inject)  = 0.d0
                clinject%erate(counter_inject) = 0.d0
                N1 = 0.d0
                if (npart>0) then
                  do k=1,4
                     if (ipart(k)>0) then
                       clinject%h(counter_inject) = clinject%h(counter_inject) + cl%h(ipart(k))*refdist/distmin_per_quadrant(k)
                       clinject%b(counter_inject) = clinject%b(counter_inject) + cl%b(ipart(k))*refdist/distmin_per_quadrant(k)
                       clinject%etot(counter_inject) = clinject%etot(counter_inject) + cl%etot(ipart(k))*refdist/distmin_per_quadrant(k)
                       clinject%erate(counter_inject) = clinject%erate(counter_inject) + cl%erate(ipart(k))*refdist/distmin_per_quadrant(k)
                       N1 = N1 + (refdist/distmin_per_quadrant(k))
                     endif
                  enddo
                  clinject%h(counter_inject)     = clinject%h(counter_inject)/N1
                  clinject%b(counter_inject)     = clinject%b(counter_inject)/N1
                  clinject%etot(counter_inject)  = clinject%etot(counter_inject)/N1
                  clinject%erate(counter_inject) = clinject%erate(counter_inject)/N1
                else
                  counter_closest = counter_closest + 1
                  write(*,'(a,I10,a,I10,a,I10,a,I10)') 'WARNING: The cell ic=',ic,' icell=',icell,' and jcell=',jcell,' did not find any closest particle with velocity conditions'
                  clinject%h(counter_inject)     = cl%h(closest)
                  clinject%b(counter_inject)     = cl%b(closest)
                  clinject%etot(counter_inject)  = cl%etot(closest)
                  clinject%erate(counter_inject) = cl%erate(closest)
                endif
 
             !endif
 
             clinject%b(counter_inject)      = min(clinject%b(counter_inject),clinject%h(counter_inject))
             clinject%cell(counter_inject)   = ic
             clinject%active(counter_inject) = .true.
             clinject%icy(counter_inject)    = jcell
             clinject%icx(counter_inject)    = icell

          enddo !ii
          grid%nn(ic) = grid%nn(ic) + ninject_per_cell(ic)
       endif    !ninject_cell>0
    enddo       !ic
    !!! omp end parallel do

    if (counter_inject /= ninject) stop 'pb inject in update_cloud_structure'
    write(*,'(a,F7.3,a,F7.3,a,F7.3,a)') 'Method for interpolation during injection: shape function = ', (float(counter_interpolated)/float(ninject))*100.d0, &
                                        ' %, average = ', (float(counter_averaged)/float(ninject))*100.d0, &
                                        ' %, closest = ', (float(counter_closest)/float(ninject))*100.d0,' %'

    deallocate(cell_injected)
  endif !if (ninject>0)
    
  !removal ---------------------------------------------------------------------------
  allocate(remove(cl%npcl))
  remove=.false.

  if (nremove>0) then

    !$omp parallel do shared(ncell,grid,cl,nremove_per_cell,remove) private(ic,nnic,distmin,ii,k1,dist,ndist,x,y,k2,k,ichoice)
    do ic=1,ncell
       
       if (nremove_per_cell(ic) > 0) then
          nnic     = grid%nn(ic)
          do ii=1,nremove_per_cell(ic)
             distmin=1.d30
             !loop through points in cell
             do k1=1,nnic
                if (.not.remove(grid%pair(k1,ic))) then
                   dist=0.d0
                   ndist=0
                   x=cl%x(grid%pair(k1,ic))
                   y=cl%y(grid%pair(k1,ic))
                   !calculate summed distance to all other points in cell that
                   !haven't been chosen to be removed yet
                   do k2=1,nnic
                      if (.not.remove(grid%pair(k2,ic)) .and. k1/=k2) then
                         dist=dist+(x-cl%x(grid%pair(k2,ic)))**2 + (y-cl%y(grid%pair(k2,ic)))**2  
                         ndist=ndist+1
                      end if
                   end do 
                   !add distances to cell corner nodes
                   do k=1,4
                      dist=dist+(x-grid%x(grid%icon(k,ic)))**2+(y-grid%y(grid%icon(k,ic)))**2
                      ndist=ndist+1
                   end do
                   !find smallest average distance to other points and corner nodes
                   dist=dist/ndist
                   if (dist<distmin) then
                      distmin=dist
                      ichoice=grid%pair(k1,ic)
                   end if
                end if
             end do
             !remove point with smallest distance
             remove(ichoice)=.true.
          end do ! do ii=1,nremove_per_cell(ic)
          grid%nn(ic) = grid%nn(ic) - nremove_per_cell(ic)
       end if    !if (nremove_per_cell(ic) > 0)

    end do
    !$omp end parallel do
    !end loop elements
   
    if(count(remove)/=nremove) stop 'pb inject in update_cloud_structure'

  endif !endif nremove>0

  allocate(not_remove(cl%npcl))
  not_remove = .not. remove

  deallocate(ninject_per_cell)
  deallocate(nremove_per_cell)


  !update arrays ---------------------------------------------------------------------------
  if (nremove>0) then
    call array_cut(cl%x                     ,not_remove)
    call array_cut(cl%y                     ,not_remove)
    call array_cut(cl%h                     ,not_remove)
    call array_cut(cl%b                     ,not_remove)
    call array_cut(cl%etot                  ,not_remove)
    call array_cut(cl%erate                 ,not_remove)
    call array_cut(cl%icy                   ,not_remove)
    call array_cut(cl%cell                  ,not_remove)
    call array_cut(cl%active                ,not_remove)
    call array_cut(cl%icx                   ,not_remove)
    call array_cut(cl%closest_node          ,not_remove)
  endif !endif nremove>0
  if (ninject>0) then
    call array_affix(cl%x                   ,clinject%x)
    call array_affix(cl%y                   ,clinject%y)
    call array_affix(cl%h                   ,clinject%h)
    call array_affix(cl%b                   ,clinject%b)
    call array_affix(cl%etot                ,clinject%etot)
    call array_affix(cl%erate               ,clinject%erate)
    call array_affix(cl%icy                 ,clinject%icy)
    call array_affix(cl%cell                ,clinject%cell)
    call array_affix(cl%active              ,clinject%active)
    call array_affix(cl%icx                 ,clinject%icx)
    call array_affix(cl%closest_node        ,clinject%closest_node)
  endif !endif ninject>0
  !free memory
  if (ninject>0) then
      deallocate(clinject%x)
      deallocate(clinject%y)
      deallocate(clinject%h)
      deallocate(clinject%b)
      deallocate(clinject%etot)
      deallocate(clinject%erate)
      deallocate(clinject%cell)
      deallocate(clinject%active)
      deallocate(clinject%icy)
      deallocate(clinject%icx)
      deallocate(clinject%closest_node)
  endif !endif ninject>0
  deallocate(not_remove)
  deallocate(remove)

  if (ninject>0 .or. nremove>0) then

  cl%npcl=npcl_new

  if (sum(grid%nn)/=cl%npcl) stop 'pb grid%nn in update_cloud (1)'
 
  if (allocated(grid%pair)) deallocate(grid%pair)
  allocate(grid%pair(maxval(grid%nn),ncell))

  grid%nn = 0
  do ip=1,cl%npcl
    ic                            = cl%cell(ip)
    grid%nn(ic)                   = grid%nn(ic) + 1
    grid%pair(grid%nn(ic),ic)     = ip
  enddo

  if (sum(grid%nn)/=cl%npcl) stop 'pb grid%nn in update_cloud (2)'

  write(*,'(a,2i15)')    'located grid%nn min/max',minval(grid%nn),maxval(grid%nn)

  endif !if (ninject>0 .or. nremove>0) then

  endif !any nnic <grid%nmin or nnic>2*grid%nmax

  return

end subroutine update_cloud_v2

!--------------------------------------------------------------------------

subroutine update_cloud_v1 (ierr)

  use FastScapeContext

  implicit none

  integer, intent(out) :: ierr
  integer ic, jcell, icell, ncellx, ncelly, ncell, nnic, npcl_max, npcl_min
  integer ninject, nremove, counter_inject
  integer ip
  integer, dimension(:), allocatable :: ninject_per_cell
  integer, dimension(:), allocatable :: nremove_per_cell
  type (cloud) :: clinject
  double precision x1,x2,x3,x4,x5,x6,x7
  double precision y1,y2,y3,y4,y5,y7, ymax, ymin
  double precision xip, yip, xmin, xmax, distmin, dist, x, y
  logical, dimension(:), allocatable :: remove
  logical, dimension(:), allocatable :: not_remove
  logical left, top
  integer npcl_in_quadrant(4),chosen_quadrant(1), ichoice, k, k2, ii, k1, ndist, npcl_new

  double precision distmin_per_quadrant(4), dx, dy, refdist
  integer ipart(4)
  integer npart, jcell_k, icell_k, ic_k
  integer counter_interpolated, counter_averaged, counter_closest
  double precision N1, N2, N3, N4, r, s

  ierr    = 0
  ncelly  = ny-1
  ncellx  = nx-1
  ncell   = ncellx*ncelly

  dx      = xl/(nx-1)
  dy      = yl/(ny-1)
  refdist = sqrt(dx**2+dy**2)

  ninject = 0
  nremove = 0

  allocate(ninject_per_cell(ncell))
  allocate(nremove_per_cell(ncell))
  ninject_per_cell = 0
  nremove_per_cell = 0

  ! diagnosis
  ic=0
  do jcell=1,ncelly
  do icell=1,ncellx
     ic       = ic+1
     nnic     = grid%nn(ic)
     npcl_max = grid%nmax
     npcl_min = grid%nmin
  
     if (nnic < npcl_min) then
       ninject_per_cell(ic) = npcl_min-nnic
     endif
     if (nnic > npcl_max) then
       nremove_per_cell(ic) = nnic-npcl_max
     endif

     ninject = ninject + ninject_per_cell(ic)
     nremove = nremove + nremove_per_cell(ic)

  end do
  end do

  write(*,'(a,i6,a)') 'injecting ',ninject,' cloud points'
  write(*,'(a,i6,a)') 'removing  ',nremove,' cloud points'

  npcl_new = cl%npcl + ninject-nremove

  write(*,'(a,2i9)') 'npcl new/old =',npcl_new,cl%npcl


  !injection ---------------------------------------------------------------------------
  if (ninject>0) then
    allocate(clinject%x(ninject))
    allocate(clinject%y(ninject))
    allocate(clinject%h(ninject))
    allocate(clinject%b(ninject))
    allocate(clinject%etot(ninject))
    allocate(clinject%erate(ninject))
    allocate(clinject%cell(ninject))
    allocate(clinject%active(ninject))
    allocate(clinject%icy(ninject))
    allocate(clinject%icx(ninject))
    counter_inject=0
    counter_interpolated=0
    counter_averaged=0
    counter_closest=0
    !loop elements
    ic=0
    do jcell=1,ncelly
    do icell=1,ncellx
       ic=ic+1
       nnic=grid%nn(ic)
       if (ninject_per_cell(ic) > 0) then
         x1=grid%x(grid%icon(1,ic)) ; y1=grid%y(grid%icon(1,ic))
         x2=grid%x(grid%icon(2,ic)) ; y2=grid%y(grid%icon(2,ic))
         x3=grid%x(grid%icon(3,ic)) ; y3=grid%y(grid%icon(3,ic))
         x4=grid%x(grid%icon(4,ic)) ; y4=grid%y(grid%icon(4,ic))
         x5 = x1                    ; y5 = 0.5d0*(y1+y4) 
         x6 = 0.5d0*(x1+x3)
         x7 = x2                    ; y7 = 0.5d0*(y2+y3)

         npcl_in_quadrant=0
         do ip=1,nnic
            xip=cl%x(grid%pair(ip,ic))
            yip=cl%y(grid%pair(ip,ic))

            left = (xip < x6 ) 
            
            top = (yip > (y7-y5)/(x7-x5)*(xip-x5) + y5)

            if (       left  .and. (.not.top)) npcl_in_quadrant(1)=npcl_in_quadrant(1)+1
            if ( (.not.left) .and. (.not.top)) npcl_in_quadrant(2)=npcl_in_quadrant(2)+1
            if ( (.not.left) .and.       top ) npcl_in_quadrant(3)=npcl_in_quadrant(3)+1
            if (       left  .and.       top ) npcl_in_quadrant(4)=npcl_in_quadrant(4)+1
         end do

         if (sum(npcl_in_quadrant)/=nnic) stop 'youpla pb'

         do ii=1,ninject_per_cell(ic)

            counter_inject=counter_inject+1

            chosen_quadrant = minloc (npcl_in_quadrant) 

            call random_number(xip)
            call random_number(yip)

            !get coord of injected point
            if (chosen_quadrant(1)==1) then

               xmin=x1
               xmax=x6

               clinject%x(counter_inject)= xmin + xip*(xmax-xmin)

               ymin=(y2-y1)/(x2-x1)*(clinject%x(counter_inject)-x1) + y1
               ymax=(y7-y5)/(x7-x5)*(clinject%x(counter_inject)-x5) + y5

               clinject%y(counter_inject)= ymin + yip*(ymax-ymin)

               npcl_in_quadrant(1)=npcl_in_quadrant(1)+1

            elseif (chosen_quadrant(1)==2) then

               xmin=x6
               xmax=x2

               clinject%x(counter_inject)= xmin + xip*(xmax-xmin)

               ymin=(y2-y1)/(x2-x1)*(clinject%x(counter_inject)-x1) + y1
               ymax=(y7-y5)/(x7-x5)*(clinject%x(counter_inject)-x5) + y5

               clinject%y(counter_inject)= ymin + yip*(ymax-ymin)

               npcl_in_quadrant(2)=npcl_in_quadrant(2)+1

            elseif (chosen_quadrant(1)==3) then

               xmin=x6
               xmax=x2

               clinject%x(counter_inject)= xmin + xip*(xmax-xmin)

               ymin=(y7-y5)/(x7-x5)*(clinject%x(counter_inject)-x5) + y5
               ymax=(y3-y4)/(x3-x4)*(clinject%x(counter_inject)-x4) + y4

               clinject%y(counter_inject)= ymin + yip*(ymax-ymin)

               npcl_in_quadrant(3)=npcl_in_quadrant(3)+1

            else

               xmin=x1
               xmax=x6

               clinject%x(counter_inject)= xmin + xip*(xmax-xmin)

               ymin=(y7-y5)/(x7-x5)*(clinject%x(counter_inject)-x5) + y5
               ymax=(y3-y4)/(x3-x4)*(clinject%x(counter_inject)-x4) + y4

               clinject%y(counter_inject)= ymin + yip*(ymax-ymin)

               npcl_in_quadrant(4)=npcl_in_quadrant(4)+1

            end if

            if ( (clinject%y(counter_inject) > y3) .or. (clinject%y(counter_inject) > y4) .or. (clinject%y(counter_inject) < y1) .or. (clinject%y(counter_inject) < y2) .or. &
                 (clinject%x(counter_inject) > x3) .or. (clinject%x(counter_inject) > x2) .or. (clinject%x(counter_inject) < x4) .or. (clinject%x(counter_inject) < x1) ) then
              write(*,*) 'WARNING: something weird with injection new particle outside the cell ', clinject%x(counter_inject),clinject%y(counter_inject)
              write(*,*) '4 cell nodes coordinates are :'
              write(*,*) ' node 1 ',x1,y1
              write(*,*) ' node 2 ',x2,y2
              write(*,*) ' node 3 ',x3,y3
              write(*,*) ' node 4 ',x4,y4
            endif

            !find closest cloud point to new point
            if (nnic>0) then
               ichoice = grid%pair(1,ic)
               distmin = 1.d30
               distmin_per_quadrant = 1.d30
               ipart = -1
               do k=1,nnic
                  dist=(clinject%x(counter_inject)-cl%x(grid%pair(k,ic)))**2 + &
                       (clinject%y(counter_inject)-cl%y(grid%pair(k,ic)))**2 
                  if (dist<distmin) then
                     distmin=dist
                     ichoice=grid%pair(k,ic)
                  end if
                  if ( cl%x(grid%pair(k,ic)) < clinject%x(counter_inject) .and. &
                       cl%y(grid%pair(k,ic)) < clinject%y(counter_inject) .and. &
                       dist < distmin_per_quadrant(1) ) then
                     distmin_per_quadrant(1) = dist
                     ipart(1)                = grid%pair(k,ic)
                     cycle
                  endif
                  if ( cl%x(grid%pair(k,ic)) > clinject%x(counter_inject) .and. &
                       cl%y(grid%pair(k,ic)) < clinject%y(counter_inject) .and. &
                       dist < distmin_per_quadrant(2) ) then
                     distmin_per_quadrant(2) = dist
                     ipart(2)                = grid%pair(k,ic)
                     cycle
                  endif
                  if ( cl%x(grid%pair(k,ic)) > clinject%x(counter_inject) .and. &
                       cl%y(grid%pair(k,ic)) > clinject%y(counter_inject) .and. &
                       dist < distmin_per_quadrant(3) ) then
                     distmin_per_quadrant(3) = dist
                     ipart(3)                = grid%pair(k,ic)
                     cycle
                  endif
                  if ( cl%x(grid%pair(k,ic)) < clinject%x(counter_inject) .and. &
                       cl%y(grid%pair(k,ic)) > clinject%y(counter_inject) .and. &
                       dist < distmin_per_quadrant(4) ) then
                     distmin_per_quadrant(4) = dist
                     ipart(4)                = grid%pair(k,ic)
                     cycle
                  endif
               end do
               npart = 0
               do k=1,4
                  if (ipart(k)>0) npart = npart + 1
               enddo
            else
               ipart   = -1
               npart   = 0
               ichoice = -1
               write(*,*) "WARNING: nnic = 0 (update_cloud)"
               !stop 'nnic = 0 | not possible'
            endif !nnic>0

            if (npart<4) then
              ! if no particle in a quadrant we find the closest particle in the
              ! adjacent cell
              do k=1,4
                if (ipart(k)<0) then
                  if (k==1) then
                    icell_k = icell - 1
                    jcell_k = jcell - 1
                    if ((icell_k<1).or.(jcell_k<1)) cycle
                  else if (k==2) then
                    icell_k = icell + 1
                    jcell_k = jcell - 1
                    if ((icell_k>ncellx).or.(jcell_k<1)) cycle
                  else if (k==3) then
                    icell_k = icell + 1
                    jcell_k = jcell + 1
                    if ((icell_k>ncellx).or.(jcell_k>ncelly)) cycle
                  else if (k==4) then
                    icell_k = icell - 1
                    jcell_k = jcell + 1
                    if ((icell_k<1).or.(jcell_k>ncelly)) cycle
                  endif
                  ic_k      = (jcell_k-1)*ncellx+icell_k
                  distmin   = 1.d30
                  do ip=1,grid%nn(ic_k)
                    dist=(clinject%x(counter_inject)-cl%x(grid%pair(ip,ic_k)))**2 + &
                         (clinject%y(counter_inject)-cl%y(grid%pair(ip,ic_k)))**2
                    if (dist<distmin) then
                       distmin   = dist
                       ipart(k)  = grid%pair(ip,ic_k)
                    end if
                  enddo
                  npart = npart + 1
                endif
              enddo
            endif

            if (npart==4) then
               counter_interpolated = counter_interpolated + 1
               xip=clinject%x(counter_inject)
               yip=clinject%y(counter_inject)

               xmin=(cl%x(ipart(4))-cl%x(ipart(1)))/(cl%y(ipart(4))-cl%y(ipart(1)))*(yip-cl%y(ipart(1))) + cl%x(ipart(1))
               xmax=(cl%x(ipart(3))-cl%x(ipart(2)))/(cl%y(ipart(3))-cl%y(ipart(2)))*(yip-cl%y(ipart(2))) + cl%x(ipart(2))
               r=((xip-xmin)/(xmax-xmin) -0.5d0 ) *2.d0

               ymin=(cl%y(ipart(2))-cl%y(ipart(1)))/(cl%x(ipart(2))-cl%x(ipart(1)))*(xip-cl%x(ipart(1))) + cl%y(ipart(1))
               ymax=(cl%y(ipart(3))-cl%y(ipart(4)))/(cl%x(ipart(3))-cl%x(ipart(4)))*(xip-cl%x(ipart(4))) + cl%y(ipart(4))
               s=((yip-ymin)/(ymax-ymin) -0.5d0 ) *2.d0
               N1=0.25d0*(1.d0-r)*(1.d0-s) 
               N2=0.25d0*(1.d0+r)*(1.d0-s) 
               N3=0.25d0*(1.d0+r)*(1.d0+s) 
               N4=0.25d0*(1.d0-r)*(1.d0+s) 
               clinject%h(counter_inject)     = N1 * cl%h(ipart(1))     + N2 * cl%h(ipart(2))     + N3 * cl%h(ipart(3))     + N4 * cl%h(ipart(4))
               clinject%b(counter_inject)     = N1 * cl%b(ipart(1))     + N2 * cl%b(ipart(2))     + N3 * cl%b(ipart(3))     + N4 * cl%b(ipart(4))
               clinject%etot(counter_inject)  = N1 * cl%etot(ipart(1))  + N2 * cl%etot(ipart(2))  + N3 * cl%etot(ipart(3))  + N4 * cl%etot(ipart(4))
               clinject%erate(counter_inject) = N1 * cl%erate(ipart(1)) + N2 * cl%erate(ipart(2)) + N3 * cl%erate(ipart(3)) + N4 * cl%erate(ipart(4))
            else
               counter_averaged = counter_averaged + 1
               clinject%h(counter_inject)     = 0.d0
               clinject%b(counter_inject)     = 0.d0
               clinject%etot(counter_inject)  = 0.d0
               clinject%erate(counter_inject) = 0.d0
               N1 = 0.d0
               do k=1,4
                  if (ipart(k)>0) then
                    clinject%h(counter_inject) = clinject%h(counter_inject) + cl%h(ipart(k))*refdist/distmin_per_quadrant(k)
                    clinject%b(counter_inject) = clinject%b(counter_inject) + cl%b(ipart(k))*refdist/distmin_per_quadrant(k)
                    clinject%etot(counter_inject) = clinject%etot(counter_inject) + cl%etot(ipart(k))*refdist/distmin_per_quadrant(k)
                    clinject%erate(counter_inject) = clinject%erate(counter_inject) + cl%erate(ipart(k))*refdist/distmin_per_quadrant(k)
                    N1 = N1 + (refdist/distmin_per_quadrant(k))
                  endif
               enddo
               clinject%h(counter_inject)     = clinject%h(counter_inject)/N1
               clinject%b(counter_inject)     = clinject%b(counter_inject)/N1
               clinject%etot(counter_inject)  = clinject%etot(counter_inject)/N1
               clinject%erate(counter_inject) = clinject%erate(counter_inject)/N1               

!               counter_closest                = counter_closest + 1
!               clinject%h(counter_inject)     = cl%h(ichoice)
!               clinject%b(counter_inject)     = cl%b(ichoice)
!               clinject%etot(counter_inject)  = cl%etot(ichoice)
!               clinject%erate(counter_inject) = cl%erate(ichoice)
            endif

            clinject%b(counter_inject)      = min(clinject%b(counter_inject),clinject%h(counter_inject))
            !clinject%cell(counter_inject)   = cl%cell(ichoice)
            !clinject%active(counter_inject) = cl%active(ichoice)
            !clinject%icy(counter_inject)    = cl%icy(ichoice)
            !clinject%icx(counter_inject)    = cl%icx(ichoice)
            clinject%cell(counter_inject)   = ic
            clinject%active(counter_inject) = .true.
            clinject%icy(counter_inject)    = jcell
            clinject%icx(counter_inject)    = icell

         end do !end ii to ninject
       endif ! if (ninject_per_cell(ic) > 0)
    enddo
    enddo
    if (counter_inject /= ninject) stop 'pb inject in update_cloud_structure'
    write(*,'(a,F7.3,a,F7.3,a,F7.3,a)') 'Method for interpolation during injection: shape function = ', (float(counter_interpolated)/float(ninject))*100.d0, &
                                        ' %, average = ', (float(counter_averaged)/float(ninject))*100.d0, &
                                        ' %, closest = ', (float(counter_closest)/float(ninject))*100.d0,' %'


  endif !if (ninject>0)

    
  !removal ---------------------------------------------------------------------------
  allocate(remove(cl%npcl))
  remove=.false.

  if (nremove>0) then
    !loop elements
    ic=0
    do jcell=1,ncelly
    do icell=1,ncellx
       ic=ic+1
       
       npcl_max=grid%nmax
   
       nnic=grid%nn(ic)
       if (nremove_per_cell(ic) > 0) then
          do ii=1,nremove_per_cell(ic)
             distmin=1.d30
             !loop through points in cell
             do k1=1,nnic
                if (.not.remove(grid%pair(k1,ic))) then
                   dist=0.d0
                   ndist=0
                   x=cl%x(grid%pair(k1,ic))
                   y=cl%y(grid%pair(k1,ic))
                   !calculate summed distance to all other points in cell that
                   !haven't been chosen to be removed yet
                   do k2=1,nnic
                      if (.not.remove(grid%pair(k2,ic)) .and. k1/=k2) then
                         dist=dist+(x-cl%x(grid%pair(k2,ic)))**2 + (y-cl%y(grid%pair(k2,ic)))**2  
                         ndist=ndist+1
                      end if
                   end do 
                   !add distances to cell corner nodes
                   do k=1,4
                      dist=dist+(x-grid%x(grid%icon(k,ic)))**2+(y-grid%y(grid%icon(k,ic)))**2
                      ndist=ndist+1
                   end do
                   !find smallest average distance to other points and corner nodes
                   dist=dist/ndist
                   if (dist<distmin) then
                      distmin=dist
                      ichoice=grid%pair(k1,ic)
                   end if
                end if
             end do
             !remove point with smallest distance
             remove(ichoice)=.true.
          end do ! do ii=1,nremove_per_cell(ic)
       end if !if (nremove_per_cell(ic) > 0)
    end do
    end do
    !end loop elements
   
    if(count(remove)/=nremove) stop 'pb inject in update_cloud_structure'

  endif !endif nremove>0

  allocate(not_remove(cl%npcl))
  not_remove = .not. remove

  deallocate(ninject_per_cell)
  deallocate(nremove_per_cell)


  !update arrays ---------------------------------------------------------------------------
  if (nremove>0) then
    call array_cut(cl%x                     ,not_remove)
    call array_cut(cl%y                     ,not_remove)
    call array_cut(cl%h                     ,not_remove)
    call array_cut(cl%b                     ,not_remove)
    call array_cut(cl%etot                  ,not_remove)
    call array_cut(cl%erate                 ,not_remove)
    call array_cut(cl%icy                   ,not_remove)
    call array_cut(cl%cell                  ,not_remove)
    call array_cut(cl%active                ,not_remove)
    call array_cut(cl%icx                   ,not_remove)
  endif !endif nremove>0
  if (ninject>0) then
    call array_affix(cl%x                   ,clinject%x)
    call array_affix(cl%y                   ,clinject%y)
    call array_affix(cl%h                   ,clinject%h)
    call array_affix(cl%b                   ,clinject%b)
    call array_affix(cl%etot                ,clinject%etot)
    call array_affix(cl%erate               ,clinject%erate)
    call array_affix(cl%icy                 ,clinject%icy)
    call array_affix(cl%cell                ,clinject%cell)
    call array_affix(cl%active              ,clinject%active)
    call array_affix(cl%icx                 ,clinject%icx)
  endif !endif ninject>0
  !free memory
  if (ninject>0) then
      deallocate(clinject%x)
      deallocate(clinject%y)
      deallocate(clinject%h)
      deallocate(clinject%b)
      deallocate(clinject%etot)
      deallocate(clinject%erate)
      deallocate(clinject%cell)
      deallocate(clinject%active)
      deallocate(clinject%icy)
      deallocate(clinject%icx)
  endif !endif ninject>0
  deallocate(not_remove)
  deallocate(remove)

  cl%npcl=npcl_new

  return

end subroutine update_cloud_v1

!--------------------------------------------------------------------------

subroutine cloud_to_eul (ierr)

  use FastScapeContext

  implicit none

  integer, intent(out) :: ierr
  integer counter,i,j,k, ic, ncell, ncellx, ncelly,ip,ichoice,nneighbours,ichoice2
  integer, dimension(4) :: icell
  double precision, dimension(4) :: xcell, ycell
  double precision distmin, dist, distmincell, rcutlim, xip, yip
  integer, dimension(4) :: pair
  double precision ymin,xmin,ymax,xmax,r,s,N1,N2,N3,N4, dx, dy, maxelev 
  integer counter_interpolated, counter_averaged, counter_closest

  ierr   = 0
  ncellx = nx-1
  ncelly = ny-1
  ncell  = ncellx*ncelly
  dx     = xl/(nx-1)
  dy     = yl/(ny-1)
  rcutlim = sqrt((dx/3)**2+(dy/3)**2)/2.d0   ! 3 = multiplicator
  pair   = -1
  counter_interpolated=0
  counter_averaged=0
  counter_closest=0
  cl%closest_node = -1

  !$omp parallel do shared(nn,nx,ncell,grid,cl,rcutlim,h,b,etot,erate) private(counter,i,j,icell,distmin,nneighbours,k,ic,distmincell,ip,dist,pair,ichoice2,xcell,ycell,ichoice,xip,yip,xmin,xmax,r,ymin,ymax,s,N1,N2,N3,N4) reduction(+:counter_averaged,counter_closest,counter_interpolated)
  do counter=1,nn
    j        = ceiling(real(counter)/real(nx))
    i        = counter - (j-1)*nx
    icell(1) = (j-2)*(nx-1) + i - 1 
    icell(2) = icell(1) + 1
    icell(3) = (j-1)*(nx-1) + i
    icell(4) = icell(3) - 1 
    if (i==nx) then
      icell(3) = -1 ; icell(2) = -1
    endif
    if (i==1) then
      icell(4) = -1 ; icell(1) = -1
    endif
    if (j==1) then
      icell(1) = -1 ; icell(2) = -1
    endif
    if (j==ny) then
      icell(3) = -1 ; icell(4) = -1
    endif
    ! select closest particules in the surrounding cells
!    distmin     = 1.d30
    maxelev     = -1.d30
    nneighbours = 0
    do k=1,4
      ic = icell(k)
      if (ic>0 .and. ic<=ncell) then
        distmincell = 1.d30
        do ip=1,grid%nn(ic)
          dist=(grid%x(counter)-cl%x(grid%pair(ip,ic)))**2 + &
               (grid%y(counter)-cl%y(grid%pair(ip,ic)))**2
          if (dist<distmincell) then
            distmincell = dist
            ichoice2    = grid%pair(ip,ic)
          endif
          if (cl%h(grid%pair(ip,ic))>maxelev) then
            maxelev = cl%h(grid%pair(ip,ic))
          endif
        enddo
        nneighbours       = nneighbours+1
        pair(nneighbours) = ichoice2
        xcell(k)          = cl%x(ichoice2)
        ycell(k)          = cl%y(ichoice2)
        cl%closest_node(ichoice2) = counter 
!        if (distmincell<distmin) then
!          distmin = distmincell
!          ichoice = ichoice2
!        endif
      endif
    enddo

!    if (sqrt(distmin)>rcutlim) then
!      !interpolation
!      !write(*,'(a,3I10)') 'WARNING: no sufficiently close particle for this eul pt ',counter,i,j
!      !write(*,*) 'using nneighbours particles for interpolation ', nneighbours
!      !write(*,'(2F15.6)') sqrt(distmin), rcutlim
!      !stop 'distmin>rcutlim'
      if (nneighbours==4) then
        counter_interpolated = counter_interpolated + 1
        ! shape function
        xip=grid%x(counter)
        yip=grid%y(counter)
        xmin=(xcell(4)-xcell(1))/(ycell(4)-ycell(1))*(yip-ycell(1)) + xcell(1)
        xmax=(xcell(3)-xcell(2))/(ycell(3)-ycell(2))*(yip-ycell(2)) + xcell(2)
        r=((xip-xmin)/(xmax-xmin) -0.5d0 ) *2.d0
        ymin=(ycell(2)-ycell(1))/(xcell(2)-xcell(1))*(xip-xcell(1)) + ycell(1)
        ymax=(ycell(3)-ycell(4))/(xcell(3)-xcell(4))*(xip-xcell(4)) + ycell(4)
        s=((yip-ymin)/(ymax-ymin) -0.5d0 ) *2.d0
        N1=0.25d0*(1.d0-r)*(1.d0-s) 
        N2=0.25d0*(1.d0+r)*(1.d0-s) 
        N3=0.25d0*(1.d0+r)*(1.d0+s) 
        N4=0.25d0*(1.d0-r)*(1.d0+s) 
        h(counter)     = N1 * cl%h(pair(1))     + N2 * cl%h(pair(2))     + N3 * cl%h(pair(3))     + N4 * cl%h(pair(4))
        if (h(counter)>maxelev+epsilon(maxelev)) then
          write(*,*) 'WARNING: something weird with shape function h(counter) > maxelev ', h(counter), maxelev
          write(*,*) '4 closest particles coordinates are :'
          do k=1,4
            write(*,*) 'xcell ',k,xcell(k)
            write(*,*) 'ycell ',k,ycell(k)
            write(*,*) 'h     ',k,cl%h(pair(k))
          enddo
        endif
        b(counter)     = N1 * cl%b(pair(1))     + N2 * cl%b(pair(2))     + N3 * cl%b(pair(3))     + N4 * cl%b(pair(4))
        etot(counter)  = N1 * cl%etot(pair(1))  + N2 * cl%etot(pair(2))  + N3 * cl%etot(pair(3))  + N4 * cl%etot(pair(4))
        erate(counter) = N1 * cl%erate(pair(1)) + N2 * cl%erate(pair(2)) + N3 * cl%erate(pair(3)) + N4 * cl%erate(pair(4))
      else
        counter_averaged = counter_averaged + 1
        !model boundaries
        h(counter)     = 0.d0
        b(counter)     = 0.d0
        etot(counter)  = 0.d0
        erate(counter) = 0.d0
        do k=1,nneighbours
          h(counter)     = h(counter)     + cl%h(pair(k))
          b(counter)     = b(counter)     + cl%b(pair(k))
          etot(counter)  = etot(counter)  + cl%etot(pair(k))
          erate(counter) = erate(counter) + cl%erate(pair(k))
        enddo
        h(counter)     = h(counter)/nneighbours
        b(counter)     = b(counter)/nneighbours
        etot(counter)  = etot(counter)/nneighbours
        erate(counter) = erate(counter)/nneighbours
      endif
!    else
!      counter_closest = counter_closest + 1
!      h(counter)     = cl%h(ichoice)
!      b(counter)     = cl%b(ichoice)
!      erate(counter) = cl%erate(ichoice)
!      etot(counter)  = cl%etot(ichoice)
!    endif

  enddo
  !$omp end parallel do

  write(*,'(a,F7.3,a,F7.3,a,F7.3,a)') 'Method for interpolation: shape function = ', (float(counter_interpolated)/float(nn))*100.d0, &
                                      ' %, average = ', (float(counter_averaged)/float(nn))*100.d0, &
                                      ' %, closest = ', (float(counter_closest)/float(nn))*100.d0,' %'

  b=min(b,h)

  return

end subroutine cloud_to_eul

!--------------------------------------------------------------------------

subroutine EulToLag (h_before_sp,b_before_sp,etot_before_sp,erate_before_sp,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out) :: ierr
  double precision, intent(in), dimension(nn) :: h_before_sp,b_before_sp,etot_before_sp,erate_before_sp
  double precision, dimension(:), allocatable :: dh, db, detot, derate
  integer ic, i, ip, ncell, inode1, inode2, inode3, inode4
  double precision xnode1, ynode1, r, s, yip, dx, dy, xip
  !double precision ymin, ymax
  double precision N1, N2, N3, N4, VERYSMALL,dhp
  logical DeltaOrSea, useBilinear, useClosestNodeOption

  !Bilinear interpolation
  double precision deltafx,deltafy,deltafxy,deltafx_b,deltafy_b,deltafxy_b,deltafx_etot,deltafy_etot,deltafxy_etot,deltafx_erate,deltafy_erate,deltafxy_erate 
  double precision dxx,dyy

  VERYSMALL            = 1.d-3
  useBilinear          = .true.
  useClosestNodeOption = .false.

  ierr  = 0
  ncell = (nx-1)*(ny-1)
  dx    = xl/(nx-1)
  dy    = yl/(ny-1)
  allocate(dh(nn), db(nn), detot(nn), derate(nn))
  dh = h-h_before_sp
  db = b-b_before_sp
  detot = etot-etot_before_sp
  derate = erate-erate_before_sp

  if (useBilinear) then

  !$omp parallel do shared(ncell,grid,cl,dx,dy,dh,db,detot,derate,h,sealevel,VERYSMALL) private(ic,ip,inode1,inode2,inode3,inode4,xnode1,ynode1,deltafx,deltafy,deltafxy,deltafx_b,deltafy_b,deltafxy_b,deltafx_etot,deltafy_etot,deltafxy_etot,deltafx_erate,deltafy_erate,deltafxy_erate,i,dxx,dyy,DeltaOrSea,dhp)
  do ic=1,ncell
    inode1=grid%icon(1,ic)
    inode2=grid%icon(2,ic)
    inode3=grid%icon(3,ic)
    inode4=grid%icon(4,ic)
    deltafx  = dh(inode2) - dh(inode1)
    deltafy  = dh(inode4) - dh(inode1)
    deltafxy = dh(inode1) + dh(inode3) - dh(inode2) - dh(inode4)

    deltafx_b  = db(inode2) - db(inode1)
    deltafy_b  = db(inode4) - db(inode1)
    deltafxy_b = db(inode1) + db(inode3) - db(inode2) - db(inode4)

    deltafx_etot  = detot(inode2) - detot(inode1)
    deltafy_etot  = detot(inode4) - detot(inode1)
    deltafxy_etot = detot(inode1) + detot(inode3) - detot(inode2) - detot(inode4)
 
    deltafx_erate  = derate(inode2) - derate(inode1)
    deltafy_erate  = derate(inode4) - derate(inode1)
    deltafxy_erate = derate(inode1) + derate(inode3) - derate(inode2) - derate(inode4)
  
    xnode1=grid%x(inode1)
    ynode1=grid%y(inode1)

    DeltaOrSea = .false.
    if ( h(inode1) <= sealevel .or. h(inode2) <= sealevel .or. h(inode3) <= sealevel .or. h(inode4) <= sealevel ) then
       DeltaOrSea = .true.
    endif

    do i=1,grid%nn(ic)
       ip=grid%pair(i,ic)
       dxx = cl%x(ip) - xnode1
       dyy = cl%y(ip) - ynode1
       dhp = deltafx*(dxx/dx) + deltafy*(dyy/dy) + deltafxy*(dxx*dyy/(dx*dy)) + dh(inode1)
       if (DeltaOrSea) then
         if (cl%h(ip)< sealevel .and. cl%h(ip)+dhp>sealevel) dhp = sealevel - cl%h(ip)
         if (cl%h(ip)>=sealevel .and. dhp>0)                 dhp = 0.d0                            ! this particle close to the shoreline could be eroded or have deposition but we cannot know
                                                                                                   ! but for sure we want to avoid to transfer sediments from offshore to onshore because of interpolation
         if (cl%h(ip)>=sealevel .and. cl%h(ip)+dhp<sealevel) dhp = sealevel + VERYSMALL - cl%h(ip)
       else
         if (cl%h(ip)+dhp< sealevel) dhp = sealevel + VERYSMALL - cl%h(ip)
       endif
       
       if (cl%closest_node(ip)>0 .and. useClosestNodeOption) then
         cl%h(ip)     = cl%h(ip) + dh(cl%closest_node(ip))
         cl%b(ip)     = cl%b(ip) + db(cl%closest_node(ip))
         cl%b(ip)     = min(cl%b(ip),cl%h(ip))
         !if (cl%h(ip)< sealevel) then
         !  !cl%etot(ip)  = 0.d0  we keep tracking total continental erosion below sea level
         !  cl%erate(ip) = 0.d0
         !else
         cl%etot(ip)     = cl%etot(ip)  + detot(cl%closest_node(ip))
         cl%erate(ip)    = cl%erate(ip) + derate(cl%closest_node(ip))
         if (cl%h(ip)< sealevel) then
           cl%erate(ip) = 0.d0
         endif
       else
         cl%h(ip)     = cl%h(ip) + dhp
         cl%b(ip)     = cl%b(ip) + deltafx_b*(dxx/dx) + deltafy_b*(dyy/dy) + deltafxy_b*(dxx*dyy/(dx*dy)) + db(inode1)
         cl%b(ip)     = min(cl%b(ip),cl%h(ip))
         cl%etot(ip)  = cl%etot(ip)  + deltafx_etot*(dxx/dx) + deltafy_etot*(dyy/dy) + deltafxy_etot*(dxx*dyy/(dx*dy)) + detot(inode1)
         cl%erate(ip) = cl%erate(ip) + deltafx_erate*(dxx/dx) + deltafy_erate*(dyy/dy) + deltafxy_erate*(dxx*dyy/(dx*dy)) + derate(inode1)
         if (cl%h(ip)< sealevel) then
           !cl%etot(ip)  = 0.d0 we keep tracking total continental erosion below sea level
           cl%erate(ip) = 0.d0
         endif
         !else
         !  cl%etot(ip)  = cl%etot(ip)  + deltafx_etot*(dxx/dx) + deltafy_etot*(dyy/dy) + deltafxy_etot*(dxx*dyy/(dx*dy)) + detot(inode1)
         !  cl%erate(ip) = cl%erate(ip) + deltafx_erate*(dxx/dx) + deltafy_erate*(dyy/dy) + deltafxy_erate*(dxx*dyy/(dx*dy)) + derate(inode1)
         !endif
       endif !cl%closest_node(ip)>0
    enddo

  end do
  !$omp end parallel do

  else

  !$omp parallel do shared(ncell,grid,cl,dx,dy,dh,db,detot,derate,h,sealevel,VERYSMALL) private(ic,ip,inode1,inode2,inode3,inode4,xnode1,ynode1,i,xip,yip,r,s,N1,N2,N3,N4,DeltaOrSea,dhp)
  do ic=1,ncell
    inode1=grid%icon(1,ic)
    inode2=grid%icon(2,ic)
    inode3=grid%icon(3,ic)
    inode4=grid%icon(4,ic)
    xnode1=grid%x(inode1)
    ynode1=grid%y(inode1)

    DeltaOrSea = .false.
    if ( h(inode1) <= sealevel .or. h(inode2) <= sealevel .or. h(inode3) <= sealevel .or. h(inode4) <= sealevel ) then
       DeltaOrSea = .true.
    endif

    do i=1,grid%nn(ic)
       ip=grid%pair(i,ic)
       xip=cl%x(ip)
       yip=cl%y(ip)
       r=((xip-xnode1)/dx -0.5d0 ) *2.d0
       s=((yip-ynode1)/dy -0.5d0 ) *2.d0
       !ymin=(grid%y(inode2)-grid%y(inode1))/dx*(xip-xnode1) + grid%y(inode1)
       !ymax=(grid%y(inode3)-grid%y(inode4))/dx*(xip-xnode1) + grid%y(inode4)
       !s=((cl%y(ip)-ymin)/(ymax-ymin) -0.5d0 ) *2.d0
       N1=0.25d0*(1.d0-r)*(1.d0-s) 
       N2=0.25d0*(1.d0+r)*(1.d0-s) 
       N3=0.25d0*(1.d0+r)*(1.d0+s) 
       N4=0.25d0*(1.d0-r)*(1.d0+s)
       dhp = N1 * dh(inode1)     + N2 * dh(inode2)     + N3 * dh(inode3)     + N4 * dh(inode4)
       if (DeltaOrSea) then
         if (cl%h(ip)< sealevel .and. cl%h(ip)+dhp>sealevel) dhp = sealevel - cl%h(ip)
         if (cl%h(ip)>=sealevel .and. dhp>0)                 dhp = 0.d0                            ! this particle close to the shoreline could be eroded or have deposition but we cannot know
                                                                                                   ! but for sure we want to avoid to transfer sediments from offshore to onshore because of interpolation
         if (cl%h(ip)>=sealevel .and. cl%h(ip)+dhp<sealevel) dhp = sealevel + VERYSMALL - cl%h(ip)
       else
         if (cl%h(ip)+dhp< sealevel) dhp = sealevel + VERYSMALL - cl%h(ip)
       endif
       if (cl%closest_node(ip)>0 .and. useClosestNodeOption) then
         cl%h(ip)     = cl%h(ip) + dh(cl%closest_node(ip))
         cl%b(ip)     = cl%b(ip) + db(cl%closest_node(ip))
         cl%b(ip)     = min(cl%b(ip),cl%h(ip))
         cl%etot(ip)     = cl%etot(ip)  + detot(cl%closest_node(ip))
         cl%erate(ip)    = cl%erate(ip) + derate(cl%closest_node(ip))
         if (cl%h(ip)< sealevel) then
         !  cl%etot(ip)  = 0.d0
           cl%erate(ip) = 0.d0
         endif
         !else
         !  cl%etot(ip)     = cl%etot(ip)  + detot(cl%closest_node(ip))
         !  cl%erate(ip)    = cl%erate(ip) + derate(cl%closest_node(ip))
         !endif
       else
         cl%h(ip)     = cl%h(ip) + dhp
         cl%b(ip)     = cl%b(ip)     + N1 * db(inode1)     + N2 * db(inode2)     + N3 * db(inode3)     + N4 * db(inode4)
         cl%b(ip)     = min(cl%b(ip),cl%h(ip))
         cl%etot(ip)  = cl%etot(ip)  + N1 * detot(inode1)  + N2 * detot(inode2)  + N3 * detot(inode3)  + N4 * detot(inode4)
         cl%erate(ip) = cl%erate(ip) + N1 * derate(inode1) + N2 * derate(inode2) + N3 * derate(inode3) + N4 * derate(inode4)
         if (cl%h(ip)< sealevel) then
           !cl%etot(ip)  = 0.d0
           cl%erate(ip) = 0.d0
         endif
         !else
         !  cl%etot(ip)  = cl%etot(ip)  + N1 * detot(inode1)  + N2 * detot(inode2)  + N3 * detot(inode3)  + N4 * detot(inode4)
         !  cl%erate(ip) = cl%erate(ip) + N1 * derate(inode1) + N2 * derate(inode2) + N3 * derate(inode3) + N4 * derate(inode4)
         !endif
       endif
    end do
  end do
  !$omp end parallel do

  endif

  deallocate(dh,db,detot,derate)
  return

end subroutine EulToLag


