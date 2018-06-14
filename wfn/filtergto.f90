! Determine the maximum distance at which it is necessary to compute
! a shell from the value of the 'cuttz' variable
subroutine filtergto ()

  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: rp, ip
  use mod_param, only: verbose
  use mod_wfn, only: ncent, numshells, cuttz, nlm, oexp, ityp, &
                     nden, rcutte, iden, rint, nshell, ishell, &
                     atcenter, nuexp, ngroup, nzexp, rmaxatom
  implicit none

  real(kind=rp) :: dis, x1, xmax, zz
  integer(kind=ip) :: ic, i, j, jc, k, lsum, m
  character(len=1) :: lb(0:5), lbl
  data (lb(i),i=0,5) /'S','P','D','F','G','H'/
  logical :: wrout, okcen
  integer(kind=ip), parameter :: ncentw = 100

  if (ncent.gt.ncentw) then
    wrout = .false.
  else
    wrout = .true.
  end if

  ! Maximum distance at which it is necessary to compute a shell.
  if (verbose) write (uout,223) rmaxatom
  if (verbose) write (uout,222) cuttz
  do ic = 1,ncent
    if (wrout.and.verbose) write (uout,210) ic
    do m = 1,ngroup(ic)
      i = nuexp(ic,m,1)
      lsum = nlm(ityp(i),1)+nlm(ityp(i),2)+nlm(ityp(i),3)
      zz = oexp(i)
      x1 = 0.1_rp
      do 
        if (x1**lsum*exp(-zz*x1*x1).le.abs(cuttz)) exit
        x1 = x1 + 0.1_rp
      end do
      rcutte(ic,m) = x1
      if (wrout.and.verbose) then
        lbl = lb(lsum)
        write (uout,613) lbl,zz,x1,(nuexp(ic,m,k),k=1,nzexp(ic,m))
      end if
    end do
  end do

  nden(1:ncent) = 0_ip
  nshell(1:ncent) = 0_ip
  if (wrout.and.verbose) then
    write (uout,'(1x,a,i0)') '# Total number of shells = ',numshells
  end if
  do ic = 1,ncent
    xmax = rmaxatom
    do jc = 1,ncent
      dis = rint(ic,jc)
      okcen = .false.
      ! The shell 'm' of center 'jc' contributes to the density, orbitals, 
      ! orbital products, etc, at any point inside the center 'ic' if the 
      ! following condition holds.
      do m = 1,ngroup(jc)
        if (dis.lt.xmax+rcutte(jc,m)) then
          nshell(ic) = nshell(ic) + 1_ip
          atcenter(ic,nshell(ic)) = jc
          ishell(ic,nshell(ic)) = m
          okcen = .true.
        end if
      end do
      if (okcen) then
        nden(ic) = nden(ic) + 1_ip
        iden(ic,nden(ic)) = jc
      end if
    end do
  end do
  if (wrout.and.verbose) then
    do ic = 1,ncent
      write (uout,300) nshell(ic),ic
      write (uout,301) (ishell(ic,j),atcenter(ic,j),j=1,nshell(ic))
    end do
  end if

  do ic = 1,ncent
    do m = 1,ngroup(ic)
      rcutte(ic,m) = rcutte(ic,m)*rcutte(ic,m)
    end do
  end do

 222  format (1x,'# CUTOFF for GTOs, eps = ',e17.10)
 223  format (1x,'# RMAXATOM for Shells = ',e17.10)
 300  format (' # ',i0,' shells contribute to the basin of center ',i0, &
       /,' # [ shell(atom) means shell number "shell" of atom "atom" ]')
 301  format (1000(1x,8(I6,'(',I4,')'),/))
 210  format (1x,'# CENTER ',i0)
 613  format (1x,'# ',a,' Shell   Exp = ',e16.8,4x,'Cutoff = ',f13.6, &
              4x,'Primitives : ',21(1x,i0))
 
end subroutine
