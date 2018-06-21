! multiply all symmetry operators by the xleft matrix
! from the left and by the xright matrix from the right. This is
! intended to produce a similarity transformation.
subroutine transformsym (xleft, xright)

  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: rp, ip
  use mod_param, only: pi, debug
  use mod_io, only: faterr, ferror, warning
  use mod_math, only: inv, euler
  use mod_sym, only: opsym, nopsym, tolisint, opeuler, opproper, &
                     opt_identity, opangle, opaxis, optype, oporder, &
                     opt_inversion, opm, opsymbol, tolnull, opt_rotation, &
                     opt_sigma, opt_imp_rotation
  implicit none

  real(kind=rp), dimension(3,3), intent(in) :: xleft, xright

  integer(kind=ip), parameter :: MOP = 20

  integer(kind=ip) :: i, j, k, l, iop, ierr, ii, ibest
  real(kind=rp), dimension(3,3) :: xmat, xinv, xvect
  real(kind=rp), dimension(3) :: xeul, xval
  real(kind=rp) :: xtrace, xdet, xnm, ang1, angle, err, besterr
  character(len=40) :: fmt
  logical :: found

  ! Similarity transformation
  do iop = 1, nopsym
    do i = 1,3
      do j = 1,3
        xmat(i,j) = 0.0_rp
        do k = 1,3
          do l = 1,3
            xmat(i,j) = xmat(i,j) + xleft(i,k)*opsym(iop,k,l)*xright(l,j)
          end do
        end do
      end do
    end do
    opsym(iop,:,:) = xmat
  end do

  ! Recalculate the properties of the transformed matrices:
  do iop = 1,nopsym
    xmat = opsym(iop,:,:)
    xtrace = 0.0_rp
    do i = 1,3
      xtrace = xtrace + xmat(i,i)
    end do
    call inv (xmat, xinv, xdet, ierr)
    if (ierr .lt. 0) then
      call ferror ('transformsym','operator has a singular matrix', faterr)
    end if
    ! Get the Euler angles for the rotation 
    call euler (xmat, xeul)
    opeuler(iop,:) = xeul(:)
    ! Is this a proper or improper rotation?
    if (abs(abs(xdet)-1.0_rp) .gt. TOLisint) then
      call ferror ('transformsym', 'determinant is not +-1', faterr)
    else if (xdet .gt. 0.0_rp) then
      opproper(iop) = .true.
    else
      opproper(iop) = .false.
    end if
    ! Get the rotation axis except for E and i, for which it is
    ! fully degenerated else Get the rotation angle and the n 
    ! and m of the C_n^m or S_n^m symmetry operator 
    if (abs(abs(xtrace)-3.0_rp) .le. TOLisint) then
      if (xtrace .gt. 0.0_rp) then
        opsymbol(iop) = 'E'
        opaxis(iop,:) = 0.0_rp
        oporder(iop) = 1_ip
        opm(iop) = 1_ip
        opangle(iop) = pi + pi
        optype(iop) = opt_identity
      else
        opsymbol(iop) = 'i'
        opaxis(iop,1) = 0.0_rp
        opaxis(iop,2) = 0.0_rp
        opaxis(iop,3) = 1.0_rp
        oporder(iop) = 2
        opm(iop) = 1
        opangle(iop) = pi
        optype(iop) = opt_inversion
      end if
    else
      ang1 = 0.5_rp*(xtrace-xdet)
      if (abs(ang1) .gt. 1.0_rp) ang1 = sign(1.0_rp, ang1)
      angle = acos(ang1)
      if (abs(angle) .le. TOLnull) then
        xnm = 1.0_rp
        angle = pi + pi
      else
        xnm = 2.0_rp*pi/angle
      end if
      opangle(iop) = angle
      ii = 0_ip
      found = .false.
      do while(.not. found .and. ii.le.MOP)
        ii = ii + 1_ip
        found = abs(xnm*ii - nint(xnm*ii)) .le. TOLisint
      end do
      if (found) then
        oporder(iop) = nint(xnm*ii)
        opm(iop) = ii
      else
        call ferror ('transformsym', 'using best approx order', warning)
        ibest = 1_ip
        besterr = abs(xnm - nint(xnm))
        do ii = 2,MOP
          err = abs(xnm*ii - nint(xnm*ii))
          if (err .le. besterr) then
            besterr = err
            ibest = ii
          end if
        end do
        oporder(iop) = nint(xnm*ibest)
        opm(iop) = ibest
      end if
      write (fmt,200) 1+int(log10(1.0*oporder(iop))),1+int(log10(1.0*opm(iop)))
      call eigensym (xmat, xval, xvect, xdet, ierr)
      if (ierr .ne. 0) then
        call ferror ('transformsym', 'error finding the rotation axis', faterr)
      end if
      opaxis(iop,1) = xvect(1,3)
      opaxis(iop,2) = xvect(2,3)
      opaxis(iop,3) = xvect(3,3)
      if(xdet .gt. 0d0) then
        write (opsymbol(iop),fmt) 'C', oporder(iop), opm(iop)
        optype(iop) = opt_rotation
      else if(abs(xtrace - 1d0) .le. TOLisint) then
        write (opsymbol(iop),'(a)') 'sigma'
        optype(iop) = opt_sigma
      else if (oporder(iop).eq.1 .and. opm(iop).eq.1) then
        write (opsymbol(iop),'(a)') 'sigma'
        optype(iop) = opt_sigma
      else
        write (opsymbol(iop),fmt) 'S', oporder(iop), opm(iop)
        optype(iop) = opt_imp_rotation
      end if
    end if
    if (debug) write (uout,700) iop, opsymbol(iop), fmt, optype(iop), xdet, &
                      xtrace, oporder(iop), opm(iop), ierr, ((xmat(i,j),j=1,3),i=1,3)             
  end do

200 format ('(a, "_", i', i2, ', "^", i', i2, ')')
700 format (/                                                       &
      1x, '++SYMTRANSFORM:'/                                        &
      1x, 'DBG(symtransform): Trial symmetry matrix number ', i6/   &
      1x, '   Symbol  ', a/                                         &
      1x, '   Format  ', a/                                         &
      1x, '   Op type ', i6/                                        &
      1x, '   Det and trace ', 2f15.9/                              &
      1x, '   Order   ', 2i6/                                       &
      1x, '   Ierror  ', i6/                                        &
      1x, '   Matrix: '/                                            &
      (1x, 3f15.9))

end subroutine transformsym
