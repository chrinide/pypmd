! Copyright (c) 2018 
! Jose Luis Casals Sainz <jluiscasalssainz@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and 
! Evelio Francisco Miguelez <evelio@uniovi.es>. 
!
! promolden is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
! 
! promolden is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
module mod_futils

  implicit none
  private

  public :: timer, iqcksort

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

  !> Same as qcksort, for integer arrays.
  !> WARNING! It is not stable!?!?
  subroutine iqcksort(iarr, iord, first, last)

    use mod_prec, only: rp, ip
    use mod_io, only: ferror, faterr
    integer(kind=ip), dimension(:), intent(in) :: iarr !< Array to be sorted
    integer(kind=ip), dimension(:), intent(inout) :: iord !< Permutation array
    integer(kind=ip), intent(in) :: first !< First index
    integer(kind=ip), intent(in) :: last  !< Last index

    integer(kind=ip), parameter :: m = 7, nstack = 100
    real(kind=rp), parameter :: fm = 7875d0, fa = 211d0
    real(kind=rp), parameter :: fc = 1663d0, fmi = 1.2698413d-4
    real(kind=rp) :: fx
    integer(kind=ip) :: istack(nstack), jstack, l, ir, j, na, a, i, iq 

    jstack=0
    l=first
    ir=last
    fx=0
10  if(ir-l.lt.m)then
       do j=l+1,ir
          na=iord(j)
          a=iarr(na)
          do i=j-1,first,-1
             if(iarr(iord(i)).le.a) go to 12
             iord(i+1)=iord(i)
          enddo
          i=first-1
12        iord(i+1)=na
       enddo
       if(jstack.eq.0)return
       ir=istack(jstack)
       l=istack(jstack-1)
       jstack=jstack-2
    else
       i=l
       j=ir
       fx=mod(fx*fa+fc,fm)
       iq=l+int((ir-l+1)*(fx*fmi))
       na=iord(iq)
       a=iarr(na)
       iord(iq)=iord(l)
20     continue
21     if (j.ge.first) then 
          if (a.lt.iarr(iord(j))) then
             j=j-1
             goto 21
          endif
       endif
       if(j.le.i)then
          iord(i)=na
          go to 30
       endif
       iord(i)=iord(j)
       i=i+1
22     if (i.le.last) then 
          if (a.gt.iarr(iord(i))) then
             i=i+1
             goto 22
          endif
       endif
       if(j.le.i)then
          iord(j)=na
          i=j
          go to 30
       endif
       iord(j)=iord(i)
       j=j-1
       go to 20
30     jstack=jstack+2
       if(jstack.gt.nstack) call ferror ('qcksort','Increase nstack',faterr) 
       if(ir-i.ge.i-l)then
          istack(jstack)=ir
          istack(jstack-1)=i+1
          ir=i-1
       else
          istack(jstack)=i-1
          istack(jstack-1)=l
          l=i+1
       endif
    endif
    go to 10

  end subroutine iqcksort

end module mod_futils
