      subroutine findnuc ()

      use mod_prec, only: rp, ip
      use mod_io, only: string, faterr, uout, ferror
      use mod_surf, only: inuc, xyzrho
      use mod_wfn, only: xyz, charge, ncent
      implicit none

      integer(kind=ip) :: i
      real(kind=rp) :: rho, grad(3), gradmod, p(3)
      logical :: inf

      interface
        subroutine gradrho (xpoint,hi,iup,inf)
          import rp, ip
          logical, intent(inout) :: inf
          integer(kind=ip), intent(in) :: iup
          real(kind=rp), intent(in) :: hi
          real(kind=rp), intent(inout) :: xpoint(3)
        end subroutine
        subroutine pointr1 (p,rho,grad,gradmod)
          import rp
          real(kind=rp), intent(in) :: p(3)
          real(kind=rp), intent(out) :: grad(3)
          real(kind=rp), intent(out) :: rho
          real(kind=rp), intent(out) :: gradmod
        end subroutine
        subroutine pointshells (p,rho,grad,gradmod)
          import rp
          real(kind=rp), intent(in) :: p(3)
          real(kind=rp), intent(out) :: grad(3)
          real(kind=rp), intent(out) :: rho
          real(kind=rp), intent(out) :: gradmod
        end subroutine
      end interface
  
      do i = 1,ncent
        inuc = i
        p(:) = xyz(i,:)
        call gradrho (p,0.05d0,1,inf)
*       call pointr1 (p,rho,grad,gradmod)
        call pointshells (p,rho,grad,gradmod)
        if (gradmod.gt.1d-4) then
          if (charge(i).gt.2.0d0) then
            write (uout,321) i
            xyzrho(i,:) = xyz(i,:)
          else
            call ferror('findnuc', 
     &                  'failed finding nucleus '//string(i), faterr)
          end if
        else
          xyzrho(i,:) = xyz(i,:)
        end if
      end do
 
 321  format (1x,'# Assuming nuclei ',i0,' position: Check!')

      end subroutine
