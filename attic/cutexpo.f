      subroutine cutexpo ()
c
c.....Determine the maximum distance at which is necessary
c.....to compute a primitive from the value of the 'cuttz '
c
      use mod_prec, only: dp
      use mod_io, only: uout, ferror, warning, string
      use mod_param, only: cuttz, nlm, largwr
      use mod_wfn, only: ncent
      use space_for_wfnbasis, only: oexp, rcutte, ityp
      use space_for_primgto, only: nzexp, nuexp, ngroup
      implicit none

      ! Local vars
      real(kind=dp) :: der0, der1, der2, expon, x0, x1, zz
      integer :: ic, ii, isu, iter, k, m
      logical :: converge, greater
      character(len=1) :: iotipe(0:5)
      data iotipe(0) /'S'/, iotipe(1) /'P'/, iotipe(2) /'D'/,
     &     iotipe(3) /'F'/, iotipe(4) /'G'/, iotipe(5) /'H'/

      ! Determine the maximum distance at which is necessary
      ! to compute a primitive.
      if (largwr) then
        write (uout,*)
        write (uout,222) cuttz
      end if
      do ic = 1,ncent
        if (largwr) write (uout,210) ic
        do m = 1,ngroup(ic)
          zz = oexp(nuexp(ic,m,1))
          ii = nuexp(ic,m,1)
          isu = nlm(ityp(ii),1) + nlm(ityp(ii),2) + nlm(ityp(ii),3)
          converge = .false.
          iter = 0
          x0 = 1d0
          do while ((.not.converge).and.iter.le.50)
            iter = iter + 1
            x1 = sqrt((isu*log(x0)-log(cuttz))/zz)
            if (abs(x1-x0).lt.1d-5) then
              converge = .true.
            else
              x0 = x1
            end if
          end do
          if (.not.converge) then
            call ferror('cutexpo', 'cuttof for exponenet '  
     &        //string(zz,'e')//' does not converge', warning)
          end if
          greater = .true.
          iter = 0
          do while (greater)
            iter = iter + 1
            expon = exp(-zz*x1*x1)
            der0 = x1**isu*expon
            if (isu.eq.0) then
              der1 = abs(2d0*zz*x1*der0)
              der2 = abs((4*zz**2*x1**2-2*zz)*der0)
            else if (isu.eq.1) then
              der1 = abs(expon*(1d0-2d0*zz*x1*x1))
              der2 = abs(expon*(4*zz**2*x1**3-6d0*zz*x1))
            else
              der1 = x1**(isu-1)*(isu-2*zz*x1*x1)*expon
              der1 = abs(der1)
              der2 = isu*(isu-1)
              der2 = der2-(2*zz*(isu-1)+4*zz+2*zz*isu)*x1*x1
              der2 = der2+4*zz*zz*x1**4
              der2 = abs(der2*expon*x1**(isu-2))
            endif
            if (der0.le.cuttz.and.der1.le.cuttz.and.der2.le.cuttz) then
              greater = .false.
            else
              x1 = x1 + 0.01d0
            end if
          end do
          ! Store the square of the cutoff instead of the cutoff itself.
          do k = 1,nzexp(ic,m)
            rcutte(nuexp(ic,m,k)) = x1*x1 
          end do
          if (largwr) write (uout,613) iotipe(isu),zz,x1, 
     &                    (nuexp(ic,m,k),k=1,nzexp(ic,m))
        end do
      end do
   
      ! Formats
  210 format (1x,'# CENTER ',I3)
  222 format (1x,'# Re-defining CUTOFF for GTOs, eps = ',E16.10)
  613 format (1x,'# ',a,' Shell   Exp = ',e16.8,4x,'Cutoff = ',F12.6,  
     &        4x,'Primitives : ',15I5)
   
      end subroutine
