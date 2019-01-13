      subroutine computerdm (det,ndet,coef,epsdet,nmo,nint,sa,sb)
!
      use mod_prec, only: ip, size_t, rp
      use mod_io, only: faterr, ferror
      use mod_parallel, only: ichunk
      use mod_wfn, only: c1ea, c1eb, c2e
      implicit none
!
      integer, intent(in) :: ndet, nint, nmo
      integer(kind=size_t), intent(in) :: det(nint,2,ndet)
      real(kind=rp), intent(in) :: coef(ndet)
      real(kind=rp), intent(in) :: epsdet
      real(kind=rp), intent(out) :: sa(nmo), sb(nmo)
!
!.....Declare auxiliary arrays as openmp doesnt feel confortable with
!.....allocatable (from calling routine) arrays in the reduction clause, 
!.....at final stage a transfer is made.
!
      integer(kind=ip) :: nmopair
      real(kind=rp) :: tmpc1ea(nmo,nmo)
      real(kind=rp) :: tmpc1eb(nmo,nmo)
      real(kind=rp) :: tmpsa(nmo), tmpsb(nmo)
      real(kind=rp), dimension(:,:), allocatable :: tmpc2e
!
      integer(kind=size_t) :: det1(Nint,2),det2(Nint,2)
      integer(kind=size_t) :: buffer1, buffer2
      real(kind=rp) :: cd1, cd2, cd, cdt, phase,c,ct
      integer :: i1, i2, k, l, il, ll, lm, in1, in2, im, jm, jn
      integer :: inm, in, ji, jl, km, kn, mn, ishift1, ishift2
      integer :: deg, ispin, ispin1, ispin2, exc(0:2,2,2), ier
!
      interface
        integer function nexcitations (det1,det2,nint)
          import :: size_t
          implicit none
          integer(kind=size_t), intent(in) :: det1(nint,2)
          integer(kind=size_t), intent(in) :: det2(nint,2)
          integer, intent(in) :: nint
        end function nexcitations
      end interface
!
      nmopair = (nmo*(nmo+1))/2
!
!     exc(i,j,k)
!     k = 1,2   ==> alpha,beta
!     i = 1,2   ==> 1st excitation, 2nd excitation
!     exc(1,2,k) = 1st created particle in det2 with spin=k
!     exc(1,1,k) = 1st created hole     in det1 with spin=k
!     exc(2,2,k) = 2nd created particle in det2 with spin=k
!     exc(2,1,k) = 2nd created hole     in det1 with spin=k
!     exc(0,j,k) = Number of created particles in det2 (j=2) or
!                  number of created holes     in det1 (j=1) with spin=k
!
      c1ea = 0.0_rp
      c1eb = 0.0_rp
      c2e  = 0.0_rp
      sa   = 0.0_rp
      sb   = 0.0_rp

!$omp parallel &
!$omp default(none) &
!$omp shared(ndet,coef,epsdet,det,nint,c1ea,sa,c1eb,sb,c2e,nmo,ichunk,nmopair) &  
!$omp private(cd1,cd,cdt,det1,ishift1,buffer1,jn,in,ishift2,buffer2, &
!$omp jm,im,i1,i2,inm,cd2,deg,det2,exc,phase,c,ct,ispin1,ispin2,ispin,ji,km, &
!$omp jl,ll,lm,il,kn,mn,tmpc1ea,tmpsa,tmpc1eb,tmpsb,tmpc2e,ier)  
      allocate(tmpc2e(nmopair,nmopair),stat=ier)
      if (ier.ne.0) then
        call ferror('computerdm', 'cannot allocate tmpc2e', faterr)
      end if
      tmpc1ea(:,:) = 0d0
      tmpc1eb(:,:) = 0d0
      tmpc2e(:,:)  = 0d0
      tmpsa(:)     = 0d0
      tmpsb(:)     = 0d0
!$omp do &
!$omp schedule(static,ichunk)  
      do k=1,Ndet
        cd1 = coef(k)
        cd  = cd1 * cd1
        cdt = cd + cd
!
!       if (abs(cd) <= abs(epsdet)) goto 1000 ! goto not allowed in parallel
!
!       Diagonal terms
!
        det1(:,:) = det(:,:,k)
!
!       Alpha RDMs
!
        ishift1 = 0
        do in1=1,Nint
          buffer1 = det1(in1,1)
          do while (buffer1 /= 0_8)
            jn = trailz(buffer1) + ishift1 + 1
            in = jn * (jn + 1) / 2
            tmpc1ea(jn,jn) = tmpc1ea(jn,jn) + cd
            tmpsa(jn)      = tmpsa(jn) + cd
            ishift2 = 0
            do in2=1,Nint
              buffer2 = det1(in2,1)
              do while (buffer2 /= 0_8)
                jm = trailz(buffer2) + ishift2 + 1
                if (jm < jn) then
                  im = jm * (jm + 1) / 2
                  tmpc2e(in,im) = tmpc2e(in,im) + cd
                  tmpc2e(im,in) = tmpc2e(im,in) + cd
                  i1 = max(jn,jm) 
                  i2 = min(jn,jm)
                  inm = i1 * ( i1 - 1 ) / 2 + i2
                  tmpc2e(inm,inm) = tmpc2e(inm,inm) - cdt
                endif
                buffer2 = iand(buffer2,buffer2-1_8)
              enddo
              ishift2 = ishift2 + 63
            enddo
!
!           Alpha-Beta RDM coupling terms
!
            ishift2 = 0
            do in2=1,Nint
              buffer2 = det1(in2,2)
              do while (buffer2 /= 0_8)
                jm = trailz(buffer2) + ishift2 + 1
                im = jm * (jm + 1) / 2
                tmpc2e(in,im) = tmpc2e(in,im) + cd
                tmpc2e(im,in) = tmpc2e(im,in) + cd
                buffer2 = iand(buffer2,buffer2-1_8)
              enddo
              ishift2 = ishift2 + 63
            enddo
            buffer1 = iand(buffer1,buffer1-1_8)
          end do
          ishift1 = ishift1 + 63
        end do
!
!       Beta RDMs
!
        ishift1 = 0
        do in1=1,Nint
          buffer1 = det1(in1,2)
          do while (buffer1 /= 0_8)
            jn = trailz(buffer1) + ishift1 + 1
            in = jn * (jn + 1) / 2
            tmpc1eb(jn,jn) = tmpc1eb(jn,jn) + cd
            tmpsb(jn)      = tmpsb(jn) + cd
            ishift2 = 0
            do in2=1,Nint
              buffer2 = det1(in2,2)
              do while (buffer2 /= 0_8)
                jm = trailz(buffer2) + ishift2 + 1
                if (jm < jn) then
                  im = jm * (jm + 1) / 2
                  tmpc2e(in,im) = tmpc2e(in,im) + cd
                  tmpc2e(im,in) = tmpc2e(im,in) + cd
                  i1 = max(jn,jm) 
                  i2 = min(jn,jm)
                  inm = i1 * ( i1 - 1 ) / 2 + i2
                  tmpc2e(inm,inm) = tmpc2e(inm,inm) - cdt
                endif
                buffer2 = iand(buffer2,buffer2-1_8)
              enddo
              ishift2 = ishift2 + 63
            enddo
            buffer1 = iand(buffer1,buffer1-1_8)
          end do
          ishift1 = ishift1 + 63
        end do
!1000   continue
!
!       Non-Diagonal terms
!
        do l=1,k-1
!
          cd2 = coef(l)
          cd  = cd1 * cd2
          cdt = cd + cd
          if (abs(cd) <= abs(epsdet)) cycle
          det2(:,:) = det(:,:,l)
          deg = nexcitations (det1,det2,Nint)
          select case (deg)
          case (0)
!
!           This cannot happen since we are in a Non-Diagonal term
!
            cycle
          case (1)
!
!           Single excitation
!
            call getexcitation (det1,det2,exc,deg,phase,Nint)
            c  = phase * cd
            ct = phase * cdt
!
!           The excitation corresponds to a ALPHA elecron ...
!
            if (exc(0,1,1) == 1) then
              ispin1 = 1
              ispin2 = 2
              ji = exc(1,1,ispin1)
              km = exc(1,2,ispin1)
              tmpc1ea(km,ji) = tmpc1ea(km,ji) + c
              tmpc1ea(ji,km) = tmpc1ea(ji,km) + c
            else
              ispin1 = 2
              ispin2 = 1
              ji = exc(1,1,ispin1)
              km = exc(1,2,ispin1)
              tmpc1eb(km,ji) = tmpc1eb(km,ji) + c
              tmpc1eb(ji,km) = tmpc1eb(ji,km) + c
            endif
            i1 = max(ji,km)
            i2 = min(ji,km)
            im = i1 * ( i1 - 1 ) / 2 + i2
            ishift1 = 0
            do in1=1,Nint
              buffer1 = det1(in1,ispin1)
              do while (buffer1 /= 0_8)
                jl = trailz(buffer1) + ishift1 + 1
                if (jl /= ji ) then
                  ll = jl * ( jl + 1 ) / 2
                  tmpc2e(im,ll) = tmpc2e(im,ll) + ct 
                  tmpc2e(ll,im) = tmpc2e(ll,im) + ct 
                  i1 = max(km,jl)
                  i2 = min(km,jl)
                  lm  = i1 * ( i1 - 1 ) / 2 + i2
                  i1 = max(ji,jl)
                  i2 = min(ji,jl)
                  il  = i1 * ( i1 - 1 ) / 2 + i2
                  tmpc2e(lm,il) = tmpc2e(lm,il) - ct 
                  tmpc2e(il,lm) = tmpc2e(il,lm) - ct 
                endif
                buffer1 = iand(buffer1,buffer1-1_8)
              enddo
              ishift1 = ishift1 + 63
            enddo
            ishift1 = 0
            do in1=1,Nint
              buffer2 = det1(in1,ispin2)
              do while (buffer2 /= 0_8)
                jl = trailz(buffer2) + ishift1 + 1
                ll = jl * ( jl + 1 ) / 2
                tmpc2e(im,ll) = tmpc2e(im,ll) + ct
                tmpc2e(ll,im) = tmpc2e(ll,im) + ct
                buffer2 = iand(buffer2,buffer2-1_8)
              enddo
              ishift1 = ishift1 + 63
            enddo
          case (2)
!
!           Double excitation
!
            call getexcitation (det1,det2,exc,deg,phase,Nint)
            c  = phase * cd
            ct = phase * cdt

            if (exc(0,1,1) == 2 .or. exc(0,1,2) == 2) then
!
!             Both excitations correspond to ALPHA electrons (ispin=1) or
!             both excitations correspond to BETA  electrons (ispin=2).
!
              if (exc(0,1,1) == 2) ispin=1
              if (exc(0,1,2) == 2) ispin=2
              ji = exc(1,1,ispin)
              jl = exc(2,1,ispin)
              km = exc(1,2,ispin)
              kn = exc(2,2,ispin)
              i1 = max(ji,km)
              i2 = min(ji,km)
              il = i1 * ( i1 - 1 ) / 2 + i2
              i1 = max(jl,kn)
              i2 = min(jl,kn)
              mn = i1 * ( i1 - 1 ) / 2 + i2
              tmpc2e(il,mn) = tmpc2e(il,mn) + ct
              tmpc2e(mn,il) = tmpc2e(mn,il) + ct
              i1 = max(ji,kn)
              i2 = min(ji,kn)
              in = i1 * ( i1 - 1 ) / 2 + i2
              i1 = max(jl,km)
              i2 = min(jl,km)
              lm = i1 * ( i1 - 1 ) / 2 + i2
              tmpc2e(in,lm) = tmpc2e(in,lm) - ct
              tmpc2e(lm,in) = tmpc2e(lm,in) - ct
            elseif (exc(0,1,1) == 1) then
!
!             One excitation is ALPHA and the other is BETA
!
              ji = exc(1,1,1)
              km = exc(1,2,1)
              jl = exc(1,1,2)
              kn = exc(1,2,2)      
              i1 = max(ji,km)
              i2 = min(ji,km)
              il = i1 * ( i1 - 1 ) / 2 + i2
              i1 = max(jl,kn)
              i2 = min(jl,kn)
              mn = i1 * ( i1 - 1 ) / 2 + i2
              tmpc2e(il,mn) = tmpc2e(il,mn) + ct
              tmpc2e(mn,il) = tmpc2e(mn,il) + ct
            endif
          case(3:)
!
!           Triple excitations do not contribute to the RDMs
!
            cycle
          end select
        end do
      end do
!$omp end do nowait
!$omp critical (reduction)
      c1ea = c1ea + tmpc1ea
      c1eb = c1eb + tmpc1eb
      c2e  = c2e  + tmpc2e 
      sa   = sa   + tmpsa  
      sb   = sb   + tmpsb  
!$omp end critical (reduction)
      deallocate(tmpc2e,stat=ier)
!$omp end parallel
!
      end subroutine computerdm
