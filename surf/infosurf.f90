subroutine info ()

  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: ip
  use mod_io, only: string, mline
  use mod_surf, only: rmaxsurf, epsiscp, epsilon, inuc, &
                      ntrial, rprimer, nangleb, steeper
  implicit none

  character(len=mline), dimension(5) :: rqudstr
  character(len=mline), dimension(5) :: ssteeper

  rqudstr(1) = 'Gauss-Legendre'
  rqudstr(2) = 'Clenshaw-Curtis'
  rqudstr(3) = 'Gauss-Chebychev 1st kind'
  rqudstr(4) = 'Gauss-Chebychev 2nd kind'
  rqudstr(5) = 'Perez-Jorda (Gauss-Chebychev) 2nd kind'

  ssteeper(1) = 'Runge-Kutta-Cash-Karp'
  ssteeper(2) = 'Calvo-Montijano-Randez'
  ssteeper(3) = 'Dormand-Prince method'

  write (uout,'(1x,a,1x,e13.6)') string('# Rmaxsur ='), rmaxsurf 
  write (uout,'(1x,a,1x,a)') string('# Steeper ='), string(ssteeper(steeper))
  write (uout,'(1x,a,1x,e13.6)') string('# Surface precision ='), epsilon
  write (uout,'(1x,a,1x,e13.6)') string('# EPSISCP parameter ='), epsiscp
  write (uout,'(1x,a,1x,i0)') string('# Ntrial ='), ntrial
  write (uout,'(1x,a,1x,e13.6)') string('# Rprimer ='), rprimer
! logical, allocatable, dimension(:), public :: lstart
! integer(kind=ip), allocatable, dimension(:), public :: nrsearch
! real(kind=rp), allocatable, dimension(:,:), public :: rstart
  write (uout,'(1x,a)') string('# Angular Quadratures')
  if (nangleb(1).eq.1) then
    write (uout,'(1x,a,1x,i0,1x,a,1x,i0)') &
    string('# Atom'), inuc, string('lebedev points'), nangleb(2)
  else if (nangleb(1).eq.0) then
    write (uout,'(1x,a)') string('# Phi quadrature is always trapezoidal')
    write (uout,'(1x,a,1x,i0,1x,a,1x,i0,1x,i0,1x,a)') &
    string('# Atom'), inuc, string('(ntheta,nphi,iqudt'), nangleb(2), nangleb(3), &
    string(rqudstr(nangleb(4)))
  end if

end subroutine info
