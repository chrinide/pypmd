! Given the cartesian coordinates of the n atoms of a molecule, this
! routine determines the connectivity matrix, diagonalizes it, 
! determines the characteristic polynomiun, the first n powers of the 
! connectivity matrix, and the so-called distance matrix.
! TODO: put again forall
subroutine connect (covx)

  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: rp, ip
  use mod_io, only: string, faterr, ferror
  use mod_wfn, only: xyz, ncent, charge, atnam
  use mod_datatm, only: covr
  use mod_linalg, only: jacobi
  use mod_futils, only: iqcksort
  implicit none
      
  real(kind=rp), intent(in) :: covx

  integer(kind=ip), parameter :: maxcoord = 20
  integer(kind=ip), parameter :: infinity = -1

  logical :: connected
  real(kind=rp) :: c(0:ncent), ankl
  character(len=2) :: this, oth(maxcoord)
  integer(kind=ip) :: wh(ncent,maxcoord), nu, iclus(ncent), iclaux(ncent)
  real(kind=rp) :: cnx(ncent,ncent), catom(ncent), v(ncent,ncent)  
  integer(kind=ip) :: coord(ncent), bonded(ncent), iord(ncent)
  real(kind=rp) :: xdis(3), rbond, covk, covm, dis2, d(ncent)
  integer(kind=ip) :: k, m, ichm, ichk, i, j, nwh, nbonds, nclaux, cdis(ncent)
  integer(kind=ip) :: p(0:ncent), nrot, madis(ncent,ncent), nclus
  character(len=1) :: labdis(ncent,ncent), digs(-1:18) 
  data digs / '-','-','1','2','3','4','5','6','7','8','9', &
              'a','b','c','d','e','f','g','h','i'/

  ! Compute connectivity matrix. Two atoms are bonded if the dis-
  ! tance between them is smaller the sum of their covalent radius 
  ! multiplied by a factor covx
  cnx = 0.0_rp
  do k = 1,ncent-1
    do m = k+1,ncent
      xdis(:) = xyz(k,:)-xyz(m,:)
      dis2 = xdis(1)**2+xdis(2)**2+xdis(3)**2
      ichk = int(charge(k))
      ichm = int(charge(m))
      covk = covr(ichk)
      covm = covr(ichm)
      rbond = (covk+covm)*covx
      if (dis2 .lt. rbond*rbond) then 
        cnx(k,m) = 1.0_rp
        cnx(m,k) = 1.0_rp
      end if
    end do
  end do

  ! Compute the coordinations of all the atoms.
  nbonds = 0 
  do i = 1,ncent
    coord(i) = 0
    do j = 1,ncent
      nbonds = nbonds + nint(cnx(i,j))
      coord(i) = coord(i) + nint(cnx(i,j))
    end do
  end do
  nbonds = nbonds/2

  ! Write coordination indices and connectivity matrix.
  write (uout,132) covx !TODO option
  write (uout,'(1x,a)') string('# Coordination indices and Connectivity Matrix')
  write (uout,'(1x,a)') string('# --------------------------------------------')
  do k = 1,ncent
    nwh = 0
    do m = 1,ncent
      if (nint(cnx(k,m)).eq.1) then
        nwh = nwh + 1
        if (nwh.gt.MAXCOORD) then
          call ferror ('connect', 'increase the value of maxcoord', faterr)
        end if
        wh(k,nwh) = m
      end if
    end do
    this(1:2) = atnam(k)(1:2)
    do m = 1,nwh
      oth(m)(1:2) = atnam(wh(k,m))(1:2)
    end do
    if (ncent.lt.10) then
      write (uout,2201) this,k,coord(k),(oth(m)(1:2),wh(k,m),m=1,nwh)
    else if (ncent.lt.100) then
      write (uout,2202) this,k,coord(k),(oth(m)(1:2),wh(k,m),m=1,nwh)
    else if (ncent.lt.1000) then
      write (uout,2202) this,k,coord(k),(oth(m)(1:2),wh(k,m),m=1,nwh)
    else
      write (uout,2204) this,k,coord(k),(oth(m)(1:2),wh(k,m),m=1,nwh)
    end if
  end do
  bonded = coord
  write (uout,10) nbonds

  ! Order the coordinations of all the atoms.
  forall (i=1:ncent) iord(i) = i
  call iqcksort (coord,iord,1,ncent)
  forall (i=1:ncent) catom(i) = coord(iord(i))

  ! Diagonalize the connectivity matrix.
  call jacobi (cnx,d,v,nrot)
  if (nrot .lt. 0) call ferror ('connect','fail diag connectivity', faterr)
  ! Detemine the characteristic polynomium.
  call polich (d,ncent,ncent,c)
  forall (i=0:ncent) p(i) = nint(c(i))
 
  ! Compute the distance matrix. Algorithm: See Chemical Graph Theory.
  ! Chapter 2 by O. E. Polansky, Section 2.9, Eq. (51).
  connected = .true.
  do i = 1,ncent-1
    madis(i,i) = 0
    do j = i+1,ncent
      do nu = 1,ncent
        ankl = 0.0_rp
        do k = 1,ncent
          ankl = ankl + v(i,k)*v(j,k)*d(k)**nu
        end do
        if (nint(ankl).gt.0) then
          madis(i,j) = nu
          madis(j,i) = nu
          go to 4
        end if
      end do
      ! This cluster is a disconnected graph.
      connected = .false.
      madis(i,j) = infinity
      madis(j,i) = infinity
 4    continue
    end do
  end do
  madis(ncent,ncent) = 0

  ! Compute the sum of distance indices for all the atoms and order them.
  forall (i=1:ncent) coord(i) = sum(madis(i,:))
 
  ! Write the sum of indices of distances and the distance matrix.
  if (connected) then
    write (uout,'(1x,a)') '# Connected graph'
  else
    write (uout,'(1x,a)') '# Non-Connected graph'
  end if
  write (uout,'(1x,a)') 'Distance Matrix: "-" means non connected atoms'
  if (ncent.le.100) then
    labdis(1:ncent,1:ncent) = digs(0)
    do i = 1,ncent
      do j = 1,ncent
        labdis(i,j) = digs(madis(i,j))
      end do
      write (uout,'(100a4)') (labdis(i,j),j=1,ncent)
    end do
  end if
 
  ! Order the sums of indices of distances.
  forall (i=1:ncent) iord(i) = i
  call iqcksort (coord,iord,1,ncent)
  forall (i=1:ncent) cdis(i) = coord(iord(i))
 
  ! Determine and identify the different non-connected clusters.
  forall (i=1:ncent) iclus(i) = 0
  nclus = 0
  do i = 1,ncent
    if (iclus(i).eq.0) then
      nclus = nclus + 1
      iclus(i) = nclus
      do j = i+1,ncent
        if (madis(i,j).ne.infinity) iclus(j) = nclus
      end do
    end if
  end do
  if (nclus.gt.1) then
    write (uout,7) nclus
    do i = 1,nclus
      nclaux = 0
      do j = 1,ncent
        if (iclus(j).eq.i) then 
          nclaux = nclaux + 1
          iclaux(nclaux) = j
        end if
      end do
      write (uout,8) i, nclaux
      write (uout,9) (iclaux(j),j=1,nclaux)
    end do
  end if

  ! Formats
 10   format (i4,' bonds')
 132  format (1x,'# Bonding Criterion = ',f12.5,' x Sum of covalent radii')
 2201 format (1x,'# (',a2,i1,') Coor = ',i1,' -->',30(1x,'(',a2,i1,')'))
 2202 format (1x,'# (',a2,i2,') Coor = ',i2,' -->',30(1x,'(',a2,i2,')'))
 2204 format (1x,'# (',a2,i4,') Coor = ',i4,' -->',30(1x,'(',a2,i4,')'))
 7    format (1x,'Molecule made of ',I2,' non-connected fragmets')
 8    format (1x,'Fragment ',I4,' contains ', I4,' atoms :')
 9    format (20(1x,I4))

end subroutine connect
