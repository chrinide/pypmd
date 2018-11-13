      subroutine info(filedat)
c
c.....initialize some variables
c
      use mod_prec, only: dp
      use mod_io, only: uout, faterr, ferror
      use mod_datatm, only: namatm, wgatm
      use mod_param, only: corr, rhf, rohf, uhf, nlm
      use mod_wfn, only: ncent, nmo, nprims, npor
      use mod_corr, only: nact, nalpha, nbeta, ncore, ndets, ndou,
     &                    nelact, nsin, totel 
      use space_for_wfnbasis, only: xyz, atnam, charge, ityp
      use space_for_corr, only: ialpha, ibeta
      implicit none
c
      character(len=*), intent(in) :: filedat
c
      real(kind=dp) :: wmol
      integer :: nform(0:106), nn, i, ia, ib, ier, ii, j
      integer :: maxlprim, ni, nl, nsal, nsbe, nz
      integer, allocatable, dimension(:) :: ialps, ibets
c
c.....write some info about molecule
c
      nl=131
      do while (filedat(nl:nl).eq.' ')
        nl=nl-1
      enddo
      write (uout,10) filedat(1:nl), ncent
      call cmcq ()
      do i=0,106
        nform(i)=0
      enddo
      wmol=0d0
c
      write (uout,11, advance='no')
      do i=1,ncent
        nz=int(charge(i))
        if (nz.le.103) wmol=wmol+wgatm(nz)
        nform(nz)=nform(nz)+1
      enddo
      nn=0
      do i=0,106
        nn=nform(i)
        if (nn.ne.0) then 
          write (uout,12,advance='no') namatm(i)(1:2),'(',nn,')'     
        end if
        nn=0
      enddo
      write (uout,13) wmol
c
      write (uout,20)
      do i=1,ncent
        write (uout,21) i, atnam(i)(1:4),(xyz(i,j),j=1,3)
      enddo
c
      maxlprim=0
      do i=1,nprims
        ii=ityp(i)
        ni=nlm(ii,1)+nlm(ii,2)+nlm(ii,3)
        maxlprim=max(ni,maxlprim)
      enddo
      write (uout,30) npor, nprims, nmo, maxlprim
c
      allocate (ialps(nmo),stat=ier)
      if (ier.ne.0) then  
        call ferror('info', 'cannot allocate ialps', faterr)
      endif
      allocate (ibets(nmo),stat=ier)
      if (ier.ne.0) then  
        call ferror('info', 'cannot allocate ibets', faterr)
      endif
c
      if (corr) then
        write (uout,41) 
        write (uout,31) totel, ndets, nelact, nact, ncore 
      else
        if (rhf) then
          write (uout,42) 
        elseif (rohf) then
          write (uout,43) 
        elseif (uhf) then
          write (uout,44) 
        else
          write (uout,45) 
        endif
        write (uout,46) totel,nsin,ndou,nalpha-ndou,nbeta-ndou
        if (uhf.or.rohf) then
          nsal=0
          do i=1,nalpha
            ia=ialpha(i)
            do j=1,nbeta
              if (ia.eq.ibeta(j)) goto 49
            enddo
            nsal=nsal+1
            ialps(nsal)=ia
 49         continue
          enddo
          write (uout,47) 'ALPHA',(ialps(i),i=1,nsal)
          nsbe=0
          do i=1,nbeta
            ib=ibeta(i)
            do j=1,nalpha
              if (ib.eq.ialpha(j)) goto 59
            enddo
            nsbe=nsbe+1
            ibets(nsbe)=ib
 59         continue
          enddo
          write (uout,47) 'BETA ',(ibets(i),i=1,nsbe)
        endif
      endif
c
      deallocate (ialps,stat=ier)
      if (ier.ne.0) then  
        call ferror('info', 'cannot deallocate ialps', faterr)
      endif
      deallocate (ibets,stat=ier)
      if (ier.ne.0) then  
        call ferror('info', 'cannot deallocate ibets', faterr)
      endif
c
c.....Formats
c
 10   format (//
     &        1x, '# WAVEFUNCTION ACCOUNT:',/,
     &        1x, '# File: ', a,/
     &        1x, '# Number of Centers:' ,i4)
 11   format (1x, '# Molecule: ')
 13   format (/1x,'# Molecular Weight:', f14.6)
 12   format (a2,a1,i4,a1,1x)
 20   format (1x, '# Cartesian coordinates of centers:'/)
 21   format (6x,i3,1x,a4,1x,3f14.8)
 30   format (/,1x,'# Number of Primitives         : ',i4,4x,
     &             'reduced to ',I4,/
     &        1x, '# Number of Molecular Orbitals : ', i4,/
     &        1x, '# Maximum l in Basis Set       : ', i4,/
     &        1x, '#',72('-'))
 41   format (1x, '# MULTIDETERMINANTAL WAVEFUNCTION')
 42   format (1x, '# RESTRICTED CLOSED SHELL HF WAVEFUNCTION (RHF)')
 43   format (1x, '# GENERAL MONODETERMINANTAL WAVEFUNCTION')
 44   format (1x, '# UNRESTRICTED HF WAVEFUNCTION (UHF)')
 45   format (1x, '# WAVEFUNCTION IS NOT RHF, ROHF, NOR UHF')
 31   format (1x, '# Number of Electrons          :', F16.10,/
     &        1x, '# Number of Determinants       :', i9,/
     &        1x, '# Number of Active Electrons   :', i6,/
     &        1x, '# Number of Active Orbitals    :', i6,/
     &        1x, '# Number of Core Orbitals      :', i6,/
     &        1x, '#',72('-'))
 46   format (1x, '# Number of Electrons           :',  F16.10,/
     &        1x, '# Number of Singly occupied MOs :',  i6,/
     &        1x, '# Number of Double occupied MOs :',  i6,/
     &        1x, '# Open-Shell (ALPHA/BETA) MOs   :', 2i6,/
     &        1x, '#',72('-'))
 47   format (1x, '# Open Shell ',a,' MOs : ',50I3)
c
      end subroutine
      DO 130 I = 1,NPRIMS
       if (ityp(I).GT.maxtype) GOTO 999
       if (ityp(i).gt.10) morethand=.true.
 130  CONTINUE
C
C-----Compute here the total number of electrons
C     (Evelio Francisco, Nov-2003)
C     ... and ALPHA and BETA electrons
C     ((Evelio Francisco, April-2007)
C
      do 120 i = 1,nmo
        read (iwfn,106) occ(i),eorb(i) 
        totel = totel + occ(i)
c
        if (occ(i).gt.occmax) occmax=occ(i)
        if (occ(i).lt.occmin) occmin=occ(i)
        if (occ(i).lt.0d0) then
          if (abs(abs(occ(i))-1d0).lt.1d-6) then
            nbeta=nbeta+1
            ibeta(nbeta)=i
            nsin=nsin+1
          else
            wfnat = .true.
          endif
        else
          if (abs(occ(i)-1d0).lt.1d-6) then
            !nalpha=nalpha+1
            !ialpha(nalpha)=i
            nbeta=nbeta+1
            ibeta(nbeta)=i
            nsin=nsin+1
          else if (abs(occ(i)-2d0).lt.1d-6) then
            nalpha=nalpha+1
            nbeta=nbeta+1
            ialpha(nalpha)=i
            ibeta(nbeta)=i
            ndou=ndou+1
          else
            wfnat = .true.
          endif
        endif
        read (iwfn,107) (coef(i,j),j=1,nprims)
        coef(i+nmo,1:nprims)=coef(i,1:nprims)
 120  continue
c
c-----Read the coefficients (virtual MOs also) from UDATT
c
      if (wfntfile) then
        do i=1,nmofull
          read (udatt,106) dummy1,dummy2
          read (udatt,107) (coeft(i,j),j=1,nprims)
        enddo
        close (udatt)
      endif
c
      if (abs(occmin-2d0).lt.1d-6) then
        rhf = .true.
      elseif (abs(abs(occmin)-1d0).lt.1d-6  .and.
     &        abs(occmax-1d0).lt.1d-6) then
        uhf = .true.
        spinpol = .true.
      elseif (abs(occmin-1d0).lt.1d-6 .and.
     &        abs(occmax-2d0).lt.1d-6) then
        rohf = .true.
        spinpol = .true.
      else
*       write (uout,312) occmin,occmax
*       write (uout,313) 
      endif
c
c.....If UHF in case map orbitals here
c
      if (uhf) then
        nalpha = iualpha
        nbeta = iubeta      
        do i = 1,nalpha
          ialpha(i)=i
*         write(*,*) "alpha", ialpha(i)
        end do
        do i = 1,nbeta
          ibeta(i)=i+nalpha
*         write(*,*) "beta", ibeta(i)
        end do
      end if
c
 309  format (' # ',a,' WAVEFUNCTION')
 312  format (1x,'#',/,
     & 1x,'# rdwdn: !!! This WFN is not RHF, nor UHF, nor ROHF !!!',/,
     & 1x,'# rdwdn: (OccMin,OccMax) = ',2(2x,F12.6))
 313  format (1x,'#',/,
     & 1x,'# rdwdn: !!! Only computations involving exclusively the',/,
     & 1x,'#        !!! natural molecular orbitals have any sense')
C
        multiplicity=1
c
c.....Determine the maximum distance at which is necessary
c     to compute a primitive.
c
      write (uout,*)'# '
      write (uout,222) cuttz
 222  format (1x,'# Cutoff for primitives, eps = ',E16.10)
      write (uout,*)'# '
      do i=1,nprims
        zz=oexp(i)
        isu=nlm(ityp(i),1)+nlm(ityp(i),2)+nlm(ityp(i),3)
        converge=.false.
        iter=0
        x0=1d0
        do while ((.not.converge).and.iter.le.50)
          iter=iter+1
          x1=sqrt((isu*log(x0)-log(cuttz))/zz)
          if (abs(x1-x0).lt.1d-5) then
             converge=.true.
          else
            x0=x1
          endif
        enddo
        if (.not.converge) write (uout,614) zz
c
        greater=.true.
        iter=0
        do while (greater)
          iter=iter+1
          expon=exp(-zz*x1*x1)
          der0=x1**isu*expon
          if (isu.eq.0) then
            der1=abs(2d0*zz*x1*der0)
            der2=abs((4*zz**2*x1**2-2*zz)*der0)
          elseif (isu.eq.1) then
            der1=abs(expon*(1d0-2d0*zz*x1*x1))
            der2=abs(expon*(4*zz**2*x1**3-6d0*zz*x1))
          else
            der1=x1**(isu-1)*(isu-2*zz*x1*x1)*expon
            der1=abs(der1)
            der2=isu*(isu-1)
            der2=der2-(2*zz*(isu-1)+4*zz+2*zz*isu)*x1*x1
            der2=der2+4*zz*zz*x1**4
            der2=abs(der2*expon*x1**(isu-2))
          endif
c       
          if (der0.le.cuttz.and.der1.le.cuttz.and.der2.le.cuttz) then
            greater=.false.
          else
            x1=x1+0.01d0
          endif
        enddo
        rcutte(i)=x1*x1 
      enddo
