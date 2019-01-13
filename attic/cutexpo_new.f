      subroutine cutexpo (lw,largwr)
c
c.....Determine the maximum distance at which it is necessary to compute
c     a shell from the value of the 'cuttz' variable defined in 'param.inc'
c
      use mod_prec, only: dp
      use mod_wfn, only: ncent, numshells
      use mod_param, only: cuttz, nlm
      use space_for_wfnbasis, only: oexp, ityp, nden, rcutte, iden, rint
      use space_for_primgto, only: nuexp, ngroup, nzexp
      use space_for_rmax, only: rmaxatom
      use space_for_shells, only: nshell, ishell, atcenter
      implicit none
c
      integer, intent(in) :: lw
      logical, intent(in) :: largwr
c
      real(kind=dp) :: expon, xeps
      real(kind=dp) :: dis, x1, xmax, zz, der0, der1, der2, x0
      integer :: ic, i, j, jc, k, lsum, m, iter, ii, maxcycles
      character(len=1) :: lb(0:5), lbl
      data (lb(i),i=0,5) /'S','P','D','F','G','H'/
      logical :: wrout, okcen, converge, greater
      integer, parameter :: ncentw = 100
c
c.....Maximum distance at which it is necessary to compute a shell.
c
      if (ncent.gt.ncentw) then
        wrout = .false.
      else
        wrout = .true.
      endif
c
      write (lw,222) cuttz
 222  format (1x,'# CUTOFF for GTOs, eps = ',E17.10)
      do ic=1,ncent
        if (wrout.and.largwr) write (lw,210) ic
        do m       = 1,ngroup(ic)
          i        = nuexp(ic,m,1)
          zz       = oexp(i)
          ii       = ityp(i)
          lsum     = nlm(ii,1) + nlm(ii,2) + nlm(ii,3)
          converge = .false.
          iter     = 0
          x0       = 1d0
          do while ((.not.converge).and.iter.le.maxcycles)
            iter   = iter + 1
            x1     = sqrt((lsum*log(x0)-log(cuttz))/zz)
            if (abs(x1-x0).lt.xeps) then
               converge=.true.
            else
              x0=x1
            endif
          enddo
          if (.not.converge) write (lw,614) zz
c
          greater=.true.
          iter=0
          do while (greater)
            iter=iter+1
            expon=exp(-zz*x1*x1)
            der0=x1**lsum*expon
            if (lsum.eq.0) then
              der1=abs(2d0*zz*x1*der0)
              der2=abs((4*zz**2*x1**2-2*zz)*der0)
            elseif (lsum.eq.1) then
              der1=abs(expon*(1d0-2d0*zz*x1*x1))
              der2=abs(expon*(4*zz**2*x1**3-6d0*zz*x1))
            else
              der1=x1**(lsum-1)*(lsum-2*zz*x1*x1)*expon
              der1=abs(der1)
              der2=lsum*(lsum-1)
              der2=der2-(2*zz*(lsum-1)+4*zz+2*zz*lsum)*x1*x1
              der2=der2+4*zz*zz*x1**4
              der2=abs(der2*expon*x1**(lsum-2))
            endif
            if (der0.le.cuttz.and.der1.le.cuttz.and.der2.le.cuttz) then
              greater=.false.
            else
              x1=x1+0.01d0
            endif
          enddo
          rcutte(ic,m)=x1
          if (wrout.and.largwr) then
            lbl=lb(lsum)
            write (lw,613) lbl,zz,x1,(nuexp(ic,m,k),k=1,nzexp(ic,m))
          endif
        enddo
      enddo
c
      nden(1:ncent)   = 0
      nshell(1:ncent) = 0
      if (wrout) write (lw,*) '# Total number of shells = ',numshells
      do ic=1,ncent
        xmax=rmaxatom(ic)
        do jc=1,ncent
          dis   = rint(ic,jc)
          okcen = .false.
          do m=1,ngroup(jc)
c
c-------------The shell 'm' of center 'jc' contributes to the density, 
c             orbitals, orbital products, etc, at any point inside the 
c             center 'ic' if the following condition holds.
c
            if (dis.lt.xmax+rcutte(jc,m)) then
              nshell(ic)=nshell(ic)+1
              atcenter(ic,nshell(ic))=jc
              ishell(ic,nshell(ic))=m
              okcen=.true.
            endif
          enddo
          if (okcen) then
            nden(ic)=nden(ic)+1
            iden(ic,nden(ic))=jc
          endif
        enddo
      enddo
      if (wrout.and.largwr) then
        do ic=1,ncent
          write (lw,300) nshell(ic),ic
          write (lw,301) (ishell(ic,j),atcenter(ic,j),j=1,nshell(ic))
        enddo
      endif
c
c.....Formats
c
 300  format (' # ',I3,' shells contribute to the basin of center ',I4,
     & /,' # [ shell(atom) means shell number "shell" of atom "atom" ]')
 301  format (1000(1x,8(I6,'(',I4,')'),/))
 210  format (1x,'# CENTER ',I3)
 211  format (1x,'# ',a,' Shell (Z=',e17.10,') : ',30I4)
 212  format (1x,'# This seems to be a [ ',
     &   I2,'s | ',I2,'p |',I2,'d | ',I2,'f | ',I2,'g ] basis')
 613  format (1x,'# ',a,' Shell   Exp = ',e16.8,4x,'Cutoff = ',F13.6,
     &   4x,'Primitives : ',15I5)
 614  format (1x,'# !!!! WARNING : Cutoff for exponent ',e16.8,
     &  ' does not converge !!!!')
c
      end subroutine
