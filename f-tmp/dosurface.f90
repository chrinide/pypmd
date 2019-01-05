! Determine the surface for all atoms or specified atoms
subroutine dosurface (filename)

  !$ use omp_lib, only: omp_get_wtime
  use mod_prec, only: rp, i4=>ip
  use mod_memory, only: alloc, free
  use mod_param, only: pi
  use mod_io, only: fourchar, uout, mline, string, flush_unit
  use mod_surf, only: neqsurf, insurf, rlimsurf, nlimsurf, minter, epsilon, &
                      xnuc, inuc, xyzrho, nangleb, rotgrid, sphere, rmaxsurf, &
                      allocate_space_for_surface, deallocate_space_for_surface
  implicit none
 
  character(len=*), intent(in) :: filename
 
  integer(kind=i4) :: ntheta, nphi, npang, iqudt
  real(kind=rp), allocatable, dimension(:) :: ct, st, cp, sp, angw
  real(kind=rp), allocatable, dimension(:) :: tp, tw
  integer(kind=i4) :: i, ii, ip, isave, it, j, k, nsurf
  integer(kind=i4) :: lsu, ltxt
  real(kind=rp) :: thang, time1, time2, phi
  real(kind=rp) :: rsurf(minter,2), delphi, rmaxs, rmins
  character(len=mline) :: files
  character(len=4) :: d4

  interface
    subroutine surf (cot,sit,cop,sip,rsurf,nsurf)
      import rp, i4, minter
      integer(kind=i4), intent(out) :: nsurf
      real(kind=rp), intent(in) :: cot, sit, cop, sip
      real(kind=rp), intent(out) :: rsurf(minter,2) 
    end subroutine
    subroutine weightheta (iq,x,w,n)
      import rp, i4
      integer(kind=i4), intent(in) :: iq, n
      real(kind=rp), intent(out) :: x(n), w(n)
    end subroutine
    subroutine lebgrid (ct,st,cp,sp,w,npoints)
      import rp, i4
      integer(kind=i4), intent(in) :: npoints
      real(kind=rp), intent(out) :: ct(npoints)
      real(kind=rp), intent(out) :: st(npoints)
      real(kind=rp), intent(out) :: cp(npoints)
      real(kind=rp), intent(out) :: sp(npoints)
      real(kind=rp), intent(out) :: w(npoints)
    end subroutine
    subroutine rotagrid (ct,st,cp,sp,npang)
      import rp, i4
      integer(kind=i4), intent(in) :: npang
      real(kind=rp), intent(inout) :: ct(npang)
      real(kind=rp), intent(inout) :: st(npang)
      real(kind=rp), intent(inout) :: cp(npang)
      real(kind=rp), intent(inout) :: sp(npang)
    end subroutine
  end interface
      
  ! Init
  ltxt = 999
  lsu = 998
  files = trim(filename)//".surf"
  
  ! Find nucleus
  call findnuc ()

  ! Order the indices
  do i = 1,neqsurf
    do j = i+1,neqsurf
      if (insurf(j).lt.insurf(i)) then
        isave = insurf(i)
        insurf(i) = insurf(j)
        insurf(j) = isave
      end if
    end do
  end do

  ! Begin
  call flush_unit (uout)
  do ii = 1,neqsurf
    inuc = insurf(ii)
    d4 = fourchar(inuc)
    xnuc(:) = xyzrho(inuc,:)
    call flush_unit (uout)
    write (uout,'(1x,a,1x,i0)') string('# Computing SURFACE for atom'), inuc
    if (nangleb(inuc,1).eq.1_i4) then
      npang = nangleb(inuc,2)
      call good_lebedev (npang)
      call alloc ('dosurface', 'ct', ct, npang)
      call alloc ('dosurface', 'st', st, npang)
      call alloc ('dosurface', 'cp', cp, npang)
      call alloc ('dosurface', 'sp', sp, npang)
      call alloc ('dosurface', 'angw', angw, npang)
      call lebgrid (ct,st,cp,sp,angw,npang)
    else
      ntheta = nangleb(inuc,2)
      nphi = nangleb(inuc,3)
      iqudt = nangleb(inuc,4)
      npang = ntheta*nphi  
      call alloc ('dosurface', 'tp', tp, ntheta)
      call alloc ('dosurface', 'tw', tw, ntheta)
      call alloc ('dosurface', 'ct', ct, npang)
      call alloc ('dosurface', 'st', st, npang)
      call alloc ('dosurface', 'cp', cp, npang)
      call alloc ('dosurface', 'sp', sp, npang)
      call alloc ('dosurface', 'angw', angw, npang)
      call weightheta (iqudt,tp,tw,ntheta)
      delphi = 2.0_rp*pi/nphi
      i = 0_i4
      do ip = 0,nphi-1
        phi = ip*delphi
        do it = 1,ntheta
          i = i + 1_i4
          thang = tp(it)
          ct(i) = thang
          st(i) = sqrt(1.0_rp-thang*thang)
          cp(i) = cos(phi)
          sp(i) = sin(phi)
          angw(i) = tw(it)*delphi
        end do
      end do
      call free ('dosurface', 'tp', tp)
      call free ('dosurface', 'tw', tw)
    end if
    !if (rotgrid) then
    !  call rotagrid (ct,st,cp,sp,npang)
    !nd if
    call allocate_space_for_surface (npang,minter)
    call cpu_time (time1)
    !$ time1 = omp_get_wtime()
    if (.not.sphere) then
      !$omp parallel default(none) &
      !$omp private(j,nsurf,rsurf) &
      !$omp shared(npang,ct,st,cp,sp,epsilon,rlimsurf,nlimsurf)
      !$omp do schedule(dynamic)
      do j = 1,npang
        call surf (ct(j),st(j),cp(j),sp(j),rsurf,nsurf)
        do k = 1,nsurf
          rlimsurf(j,k) = rsurf(k,2)
        end do
        nlimsurf(j) = nsurf
      end do
      !$omp end do nowait
      !$omp end parallel
    end if
    call cpu_time (time2)
    !$ time2 = omp_get_wtime()
    write (uout,'(1x,a,1x,f12.5)') string('# Elapsed seconds :'), time2-time1
    open (ltxt,file=trim(files)//"-txt"//d4)
    open (lsu,file=trim(files)//d4,form='unformatted')
    write (ltxt,1111) npang, inuc
    write (lsu) npang, inuc
    if (.not.sphere) then
      write (ltxt,3333) (nlimsurf(j),j=1,npang)
      write (ltxt,1090)
      write (lsu) (nlimsurf(j),j=1,npang)
      rmins = 1000_rp
      rmaxs = 0.0_rp
      do j = 1,npang
        nsurf = nlimsurf(j)
        rmins = min(rmins,rlimsurf(j,1))
        rmaxs = max(rmaxs,rlimsurf(j,nsurf))
        write (ltxt,2222) ct(j),st(j),cp(j),sp(j),angw(j),(rlimsurf(j,k),k=1,nsurf)
        write (lsu) ct(j),st(j),cp(j),sp(j),angw(j),(rlimsurf(j,k),k=1,nsurf)
      end do
    else
      write (ltxt,3333) (1,j=1,npang)
      write (ltxt,1090)
      write (lsu)(1,j=1,npang)
      rmins = 0.0_rp
      rmaxs = rmaxsurf(inuc)
      do j = 1,npang
        write (ltxt,2222) ct(j),st(j),cp(j),sp(j),angw(j),rmaxs
        write (lsu) ct(j),st(j),cp(j),sp(j),angw(j),rmaxs
      end do
    end if
    write (ltxt,2222) rmins,rmaxs
    write (lsu) rmins,rmaxs
    close (ltxt)
    close (lsu)
    call deallocate_space_for_surface ()
    call free ('dosurface', 'ct', ct)
    call free ('dosurface', 'ct', st)
    call free ('dosurface', 'ct', cp)
    call free ('dosurface', 'ct', sp)
    call free ('dosurface', 'angw', angw)
  end do

1090 format (9x,'cos(theta)',13x,'sin(theta)',13x,'cos(phi)',15x,'sin(phi)',15x,'weight')
1111 format (2(1x,i5),' <--- (Angular points & Atom)')
3333 format (20(1x,i2),4x,'(Surface intersections)')
2222 format (15(1x,F22.15))
 
end subroutine
