module mod_futils

  implicit none
  private

  public :: timer

contains

  !> Timer interface
  !> accumulates and prints out the elapsed times of a series of up to
  !> MPROCESS processes.
  !> Input parameters are:
  !> key ........... 0 = reset all time tables.
  !>                 1 = reset pid entry and begin counting for it.
  !>                 2 = continue the cont for pid process. Do not
  !>                     reset previous entry times.
  !>                 3 = end of pid process. Free pid entry.
  !>                 4 = end of pid process. Do not free entry.
  !>                 5 = end of run. Print out all time tables.
  !>                 6 = close all processes and print out tables.
  !> pid ........... process identification number (1..MPROCESS).
  !> name .......... process name (used only in the print out).
  !> lw ............ printer logical unit. Negative if print out is not
  !>                 desired.
  !>
  !>.key controls what data are used in the run, as the following table
  !> resumes:
  !>
  !> key value    pid      name       lw
  !> ---------  -------   -------   -------
  !>    0       ignored   ignored   ignored
  !>   1,2      output    input     ignored
  !>   3,4      input     ignored   input
  !>   5,6      ignored   ignored   input
  subroutine timer (key,pid,name,lw)

!$  use omp_lib
    use mod_prec, only: rp, ip
    use mod_io, only: ferror, warning
    implicit none

    integer(kind=ip), intent(in) :: key, lw
    integer(kind=ip), intent(inout) :: pid
    character(len=10), intent(in) :: name

    integer(kind=ip), parameter :: mprocess = 500
    character(len=10) ::  pname(mprocess)
    real(kind=rp) :: time(mprocess), cumtime(mprocess), timedum
    logical :: popen(mprocess), pused(mprocess), firsttime
    integer(kind=ip) :: pcalls(mprocess), i

    save time, cumtime, pname, popen, pused, pcalls, firsttime
    data firsttime /.true./
 
    if (key.eq.0 .or. firsttime) then ! initiallize all entries: 
      firsttime = .false.
      time(:) = 0.0_rp
      cumtime(:) = 0.0_rp
      popen(:) = .false.
      pused(:) = .false.
      pname(:) = '          '
      pcalls(:) = 0_ip
    else if (key.eq.1 .or. key.eq.2) then ! begin pid count:
      call cpu_time(timedum)
!$    timedum = omp_get_wtime()
      i = 1
 15   if (pused(i)) then
        if (pname(i)(1:8).ne.name(1:8)) then
          i = i + 1
          if (i.gt.mprocess) then
            call ferror('mod_futils/timer', 'pid out of bonds', warning)
            return
          end if
          goto 15
        end if
      end if
      pid = i
      if (key.eq.1) then
        cumtime(pid) = 0.0_rp
        pcalls(pid) = 1_ip
      else
        pcalls(pid) = pcalls(pid) + 1_ip
      end if
        time(pid) = timedum
        popen(pid) = .true.
        pused(pid) = .true.
        pname(pid) = name
    else if (key.eq.3 .or. key.eq.4) then  ! end pid accounting: 
      if (pid.le.0 .or. pid.gt.mprocess) then
         call ferror('mod_futils/timer', 'pid out of bonds', warning)
      else if (.not.popen(pid)) then
         call ferror('mod_futils/timer', 'pid unused or closed', warning)
      else
        call cpu_time(timedum)
!$      timedum = omp_get_wtime()
        time(pid) = timedum - time(pid)
        cumtime(pid) = cumtime(pid) + time(pid)
        popen(pid) = .false.
        if (lw.gt.0) write (lw,100) pname(pid),time(pid)
        if (key.eq.3) then
          pused(pid) = .false.
          pcalls(pid) = 0_ip
          cumtime(pid) = 0.0_rp
          time(pid) = 0.0_rp
        end if
      end if
    else if (key.eq.5 .or. key.eq.6) then ! print out the time tables: 
      write (lw,105)
      call cpu_time(timedum)
!$    timedum = omp_get_wtime()
      do i = 1,mprocess
        if (pused(i)) then
          if (popen(i)) then
            time(i) = timedum - time(i)
            cumtime(i) = cumtime(i) + time(i)
            if (key.eq.6) popen(i)=.false.
          end if
          write (lw,110) i,pname(i),cumtime(i),pcalls(i),popen(i)
        end if
      end do
      write (lw,115)
    else
      call ferror('mod_futils/timer', 'key value not recognized', warning)
    end if
 
100 format (/' #------------------------------------',/'#*** timer:'/ &
      ' #    process name:',a10,5x,'elapsed time (sec):',f10.3/)
105 format (/' #',53('-'),/,' #    timer:'/,' #   '/                  &
      ' # -pid--------name----------cumtime--------pcalls--popen-')
110 format (' # ',i3,6x,a10,1x,f14.6,1x,i10,6x,l2)
115 format (' #'/)

  end subroutine timer

end module mod_futils
