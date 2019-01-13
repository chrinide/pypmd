      subroutine getexcitation (det1,det2,exc,degree,phase,nint)
      use mod_prec, only: size_t, rp
      implicit none
      integer, intent(in) :: nint
      integer(kind=size_t), intent(in) :: det1(nint,2), det2(nint,2)
      integer, intent(out) :: exc(0:2,2,2)
      integer, intent(out) :: degree
      real(kind=rp), intent(out) :: phase

      interface
        integer function nexcitations (det1,det2,nint)
          import :: size_t
          implicit none
          integer(kind=size_t), intent(in) :: det1(nint,2)
          integer(kind=size_t), intent(in) :: det2(nint,2)
          integer, intent(in) :: nint
        end function nexcitations
      end interface

      degree = nexcitations (det1,det2,Nint)
      select case (degree)
      case (3:)
        degree = -1
        return
      case (2)
        call doublexcitation (det1,det2,exc,phase,Nint)
        return
      case (1)
        call singlexcitation (det1,det2,exc,phase,Nint)
        return
      case(0)
        return
      end select
      end subroutine
!
      integer function nexcitations (det1,det2,nint)
      use mod_prec, only: size_t
      implicit none
      integer, intent(in) :: nint
      integer(kind=size_t), intent(in) :: det1(nint,2), det2(nint,2)

      integer(kind=size_t) :: d1a,d2a,d1b,d2b
      integer :: l

      d1a    = det1(1,1)
      d2a    = det2(1,1)
      d1b    = det1(1,2)
      d2b    = det2(1,2)
      nexcitations = popcnt(xor(d1a,d2a)) + popcnt(xor(d1b,d2b))
      do l=2,Nint
        d1a  = det1(l,1)
        d2a  = det2(l,1)
        d1b  = det1(l,2)
        d2b  = det2(l,2)
        nexcitations = nexcitations + popcnt(xor(d1a,d2a)) &
                                    + popcnt(xor(d1b,d2b))
      end do
      nexcitations = ishft(nexcitations,-1)
      end function
!
      subroutine singlexcitation (det1,det2,exc,phase,nint)
      use mod_prec, only: size_t, rp
      implicit none
      integer, intent(in) :: nint
      integer(kind=size_t), intent(in) :: det1(nint,2)
      integer(kind=size_t), intent(in) :: det2(nint,2)
      integer, intent(out) :: exc(0:2,2,2)
      real(kind=rp), intent(out) :: phase

      integer :: tz,l,ispin,ishift,nperm
      integer :: i,j,k,m,n,high,low
      integer(kind=size_t) :: hole, particle, tmp
      real(kind=rp), parameter :: phase_dble(0:1) = (/1.d0,-1.d0/)

      exc(0,1,1) = 0
      exc(0,2,1) = 0
      exc(0,1,2) = 0
      exc(0,2,2) = 0
      do ispin = 1,2
        ishift = 0
        do l=1,Nint
          if (det1(l,ispin) == det2(l,ispin)) cycle
          tmp = xor( det1(l,ispin), det2(l,ispin) )
          particle = iand(tmp, det2(l,ispin))
          hole = iand(tmp, det1(l,ispin))
          if (particle /= 0_8) then
            tz = trailz(particle)
            exc(0,2,ispin) = 1
            exc(1,2,ispin) = tz+ishift+1
          end if
          if (hole /= 0_8) then
            tz = trailz(hole)
            exc(0,1,ispin) = 1
            exc(1,1,ispin) = tz+ishift+1
          end if
          if ( iand(exc(0,1,ispin),exc(0,2,ispin)) == 1 ) then
            low = min(exc(1,1,ispin),exc(1,2,ispin))
            high = max(exc(1,1,ispin),exc(1,2,ispin))
            j = ishft(low-1,-6)+1
            n = iand(low,63)
            k = ishft(high-1,-6)+1
            m = iand(high,63)
            if (j==k) then
              nperm = popcnt(iand(det1(j,ispin), &
              iand( ibset(0_8,m-1)-1_8, ibclr(-1_8,n)+1_8 ) ))
            else
              nperm = popcnt(iand(det1(k,ispin), ibset(0_8,m-1)-1_8)) &
                    + popcnt(iand(det1(j,ispin), ibclr(-1_8,n) +1_8))
              do i=j+1,k-1
                nperm = nperm + popcnt(det1(i,ispin))
              end do
            end if
            phase = phase_dble(iand(nperm,1))
            return
          end if
          ishift = ishift + 63
        end do
      end do
      end subroutine
!
      subroutine doublexcitation (det1,det2,exc,phase,nint)
      use mod_prec, only: size_t, rp
      implicit none
      integer, intent(in) :: nint
      integer(kind=size_t), intent(in) :: det1(nint,2), det2(nint,2)
      integer, intent(out) :: exc(0:2,2,2)
      real(kind=rp), intent(out) :: phase

      integer :: l,ispin,idx_hole
      integer :: idx_particle,ishift
      integer :: i,j,k,m,n,high,low,a,b,c,d
      integer :: nperm,tz,nexc
      integer(kind=size_t) :: hole,particle,tmp
      real(kind=rp), parameter :: phase_dble(0:1) = (/1.d0,-1.d0/)

      exc(0,1,1) = 0
      exc(0,2,1) = 0
      exc(0,1,2) = 0
      exc(0,2,2) = 0
      nexc       = 0
      nperm      = 0
      do ispin = 1,2
        idx_particle = 0
        idx_hole     = 0
        ishift       = 0
        do l=1,Nint
          if (det1(l,ispin) == det2(l,ispin)) then
            cycle
          end if
          tmp      = xor( det1(l,ispin), det2(l,ispin) )
          particle = iand(tmp, det2(l,ispin))
          hole     = iand(tmp, det1(l,ispin))
          do while (particle /= 0_8)
            tz             = trailz(particle)
            nexc           = nexc+1
            idx_particle   = idx_particle + 1
            exc(0,2,ispin) = exc(0,2,ispin) + 1
            exc(idx_particle,2,ispin) = tz+ishift+1
            particle       = iand(particle,particle-1_8)
          end do
          do while (hole /= 0_8)
            tz                    = trailz(hole)
            nexc                  = nexc+1
            idx_hole              = idx_hole + 1
            exc(0,1,ispin)        = exc(0,1,ispin) + 1
            exc(idx_hole,1,ispin) = tz+ishift+1
            hole                  = iand(hole,hole-1_8)
          end do
          if (nexc == 4) exit
          ishift = ishift + 63
        end do
        do i=1,exc(0,1,ispin)
          low  = min(exc(i,1,ispin),exc(i,2,ispin))
          high = max(exc(i,1,ispin),exc(i,2,ispin))
          j    = ishft(low-1,-6)+1
          n    = iand(low,63)
          k    = ishft(high-1,-6)+1
          m    = iand(high,63)
          if (j==k) then
            nperm = nperm + popcnt(iand(det1(j,ispin), &
                    iand( ibset(0_8,m-1)-1_8, ibclr(-1_8,n)+1_8 ) ))
          else
            nperm = nperm + popcnt(iand(det1(k,ispin), &
                    ibset(0_8,m-1)-1_8))               &
            + popcnt(iand(det1(j,ispin),               &
            ibclr(-1_8,n) +1_8))
            do l=j+1,k-1
              nperm = nperm + popcnt(det1(l,ispin))
            end do
          end if
        end do
        if (exc(0,1,ispin) == 2) then
          a = min(exc(1,1,ispin), exc(1,2,ispin))
          b = max(exc(1,1,ispin), exc(1,2,ispin))
          c = min(exc(2,1,ispin), exc(2,2,ispin))
          d = max(exc(2,1,ispin), exc(2,2,ispin))
          if (c>a .and. c<b .and. d>b) nperm = nperm + 1
          exit
        end if
      end do
      phase = phase_dble(iand(nperm,1))
      end subroutine
!
      subroutine get_particles (det1,det2,particles,nint,nmo)
      use mod_prec, only: size_t
      implicit none
      integer, intent(in) :: nint,nmo
      integer(kind=size_t), intent(in) :: det1(nint,2), det2(nint,2)
      integer(kind=size_t), intent(out) :: particles(nmo,2)

      integer(kind=size_t) :: p,position
      integer :: l,k,ispin

      do ispin = 1,2
        k = 1
        do l = 1, nint
          p = and ( xor(det1(l,1),det2(l,1)), det2(l,1))
          do while ( p /= 0)
            position = trailz (p)
            particles (k,ispin) = 1 + 64 * (l-1) + position
            p = ibclr(p,position)
            k = k + 1
          enddo
        enddo
      enddo
      end subroutine
!
      subroutine get_holes (det1,det2,holes,nint,nmo)
      use mod_prec, only: size_t
      implicit none
      integer, intent(in) :: nint,nmo
      integer(kind=size_t), intent(in) :: det1(nint,2), det2(nint,2)
      integer(kind=size_t), intent(out) :: holes(nmo,2)

      integer(kind=size_t) :: h,position
      integer :: l,k,ispin

      do ispin = 1,2
        k = 1
        do l = 1, nint
          h = and ( xor(det1(l,1),det2(l,1)), det1(l,1))
          do while ( h /= 0)
            position = trailz (h)
            holes (k,ispin) = 1 + 64 * (l-1) + position
            h = ibclr(h,position)
            k = k + 1
          enddo
        enddo
      enddo
      end subroutine
