module mod_grids

  use mod_prec, only: rp, ip
  implicit none
  private

  integer(kind=ip), parameter :: ibecke = 3

  public :: gen_becke_mesh

contains

  !> The q-mesh is uniform on the interval (-1,+1). 
  !> Transformation is
  !>            ( 1 + q ) 
  !> r  =  rmid ---------
  !>            ( 1 - q )
  subroutine rmeshb(n,rmid,r,wintr)
    
    use mod_param, only: pi  
    implicit none
    integer(kind=ip), intent(in) :: n
    real(kind=rp), intent(in) :: rmid
    real(kind=rp), intent(out), dimension(n) :: r, wintr
    
    real(kind=rp), parameter :: fourpi = 4.0_rp*pi

    real(kind=rp) :: h, q
    integer(kind=ip) :: i

    h = 1.0_rp/real(n+1,rp)
    do i = 1,n
      q = cos(h*i*pi)
      r(i) = ((1.0_rp+q)/(1.0_rp-q))*rmid
      wintr(i) = 2.0_rp*pi*h*rmid**3*(1.0_rp+q)**2.5_rp/(1.0_rp-q)**3.5_rp*fourpi
    end do

  end subroutine

  !> Original scheme
  subroutine gen_becke_mesh(sizes, radii, npoints, weights, points)

    use mod_memory, only: free, alloc
    use mod_mole, only: ncent_, coords_
    implicit none

    integer(kind=ip), intent(in) :: npoints
    integer(kind=ip), intent(in), dimension(2,ncent_) :: sizes
    real(kind=rp), intent(in), dimension(ncent_) :: radii
    real(kind=rp), intent(out), dimension(npoints) :: weights
    real(kind=rp), intent(out), dimension(3,npoints) :: points
    
    real(kind=rp) :: rr(ncent_,ncent_), smat(ncent_,ncent_), pvec(ncent_) 
    real(kind=rp) :: rmid, r, r1, r2, rmiu, rmiu2, tmps
    real(kind=rp) :: beckew, chi, uij, aij, x(3)
    integer(kind=ip) :: kk, ii, i, j, k, nrad, nang, ir, il, ib 
    real(kind=rp), allocatable :: rads(:), wrads(:)
    real(kind=rp), allocatable :: xang(:), yang(:), zang(:), wang(:)

    !> Compute interatomic distances
    rr = 0.0_rp
    do i = 1,ncent_
      do j = i+1,ncent_
        rr(i,j) = sqrt((coords_(1,i) - coords_(1,j))**2 + &
                       (coords_(2,i) - coords_(2,j))**2 + &
                       (coords_(3,i) - coords_(3,j))**2)
        rr(j,i) = rr(i,j)
      end do
    end do

    !> Gen the mesh
    kk = 0_ip
    do i = 1,ncent_
      nrad = sizes(1,i)
      nang = sizes(2,i)
      call alloc ('mod_grid', 'rads', rads, nrad)
      call alloc ('mod_grid', 'wrads', wrads, nrad)
      rmid = radii(i)
      call rmeshb (nrad,rmid,rads,wrads)
      call alloc ('mod_grid', 'xang', xang, nang)
      call alloc ('mod_grid', 'yang', yang, nang)
      call alloc ('mod_grid', 'zang', zang, nang)
      call alloc ('mod_grid', 'wang', wang, nang)
      if (nang == 6) then
        call ld0006(xang,yang,zang,wang,nang)
      else if (nang == 14) then
        call ld0014(xang,yang,zang,wang,nang)
      else if (nang == 26) then
        call ld0026(xang,yang,zang,wang,nang)
      else if (nang == 38) then
        call ld0038(xang,yang,zang,wang,nang)
      else if (nang == 50) then
        call ld0050(xang,yang,zang,wang,nang)
      else if (nang == 74) then
        call ld0074(xang,yang,zang,wang,nang)
      else if (nang == 86) then
        call ld0086(xang,yang,zang,wang,nang)
      else if (nang == 110) then
        call ld0110(xang,yang,zang,wang,nang)
      else if (nang == 146) then
        call ld0146(xang,yang,zang,wang,nang)
      else if (nang == 170) then
        call ld0170(xang,yang,zang,wang,nang)
      else if (nang == 194) then
        call ld0194(xang,yang,zang,wang,nang)
      else if (nang == 230) then
        call ld0230(xang,yang,zang,wang,nang)
      else if (nang == 266) then
        call ld0266(xang,yang,zang,wang,nang)
      else if (nang == 302) then
        call ld0302(xang,yang,zang,wang,nang)
      else if (nang == 350) then
        call ld0350(xang,yang,zang,wang,nang)
      else if (nang == 434) then
        call ld0434(xang,yang,zang,wang,nang)
      else if (nang == 590) then
        call ld0590(xang,yang,zang,wang,nang)
      else if (nang == 770) then
        call ld0770(xang,yang,zang,wang,nang)
      else if (nang == 974) then
        call ld0974(xang,yang,zang,wang,nang)
      else if (nang == 1202) then
        call ld1202(xang,yang,zang,wang,nang)
      else if (nang == 1454) then
        call ld1454(xang,yang,zang,wang,nang)
      else if (nang == 1730) then
        call ld1730(xang,yang,zang,wang,nang)
      else if (nang == 2030) then
        call ld2030(xang,yang,zang,wang,nang)
      else if (nang == 2354) then
        call ld2354(xang,yang,zang,wang,nang)
      else if (nang == 2702) then
        call ld2702(xang,yang,zang,wang,nang)
      else if (nang == 3074) then
        call ld3074(xang,yang,zang,wang,nang)
      else if (nang == 3470) then
        call ld3470(xang,yang,zang,wang,nang)
      else if (nang == 3890) then
        call ld3890(xang,yang,zang,wang,nang)
      else if (nang == 4334) then
        call ld4334(xang,yang,zang,wang,nang)
      else if (nang == 4802) then
        call ld4802(xang,yang,zang,wang,nang)
      else if (nang == 5294) then
        call ld5294(xang,yang,zang,wang,nang)
      else if (nang == 5810) then
        call ld5810(xang,yang,zang,wang,nang)
      end if

      !> Set origin and mix angular and radial part
      do ir = 1,nrad
        r = rads(ir)
        do il = 1,nang
          x = coords_(:,i) + r*(/xang(il),yang(il),zang(il)/)
          smat = 1.0_rp
          do j = 1,ncent_
            r1 = sqrt((x(1)-coords_(1,j))**2 + &
                      (x(2)-coords_(2,j))**2 + &
                      (x(3)-coords_(3,j))**2)
            do k = 1,ncent_
              if (j==k) cycle
              r2 = sqrt((x(1)-coords_(1,k))**2 + &
                        (x(2)-coords_(2,k))**2 + &
                        (x(3)-coords_(3,k))**2)
              rmiu = (r1-r2)/rr(j,k)
              chi = radii(j)/radii(k)
              uij = (chi-1.0_rp)/(chi+1.0_rp)
              aij = uij/(uij**2-1.0_rp)
              if (aij >  0.5_rp) aij =  0.5_rp
              if (aij < -0.5_rp) aij = -0.5_rp
              rmiu2 = rmiu+aij*(1.0_rp-rmiu**2)
              tmps = 1.5_rp*rmiu2-0.5_rp*rmiu2**3
              do ib = 2,ibecke
                tmps = 1.5_rp*(tmps)-0.5_rp*(tmps)**3
              end do
              smat(j,k) = 0.5_rp*(1.0_rp-tmps)
            end do
          end do
          pvec = 1.0_rp
          do ii = 1,ncent_
            pvec = pvec*smat(:,ii)
          end do
          beckew = pvec(i)/sum(pvec)
          kk = kk + 1
          weights(kk) = beckew * wrads(ir) * wang(il)
          points(:,kk) = x
        end do
      end do
      call free ('mod_grid', 'rads', rads)
      call free ('mod_grid', 'wrads', wrads)
      call free ('mod_grid', 'xang', xang)
      call free ('mod_grid', 'yang', yang)
      call free ('mod_grid', 'zang', zang)
      call free ('mod_grid', 'wang', wang)
    end do

  end subroutine 

end module mod_grids
