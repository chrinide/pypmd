      subroutine colord (o1,o2,n,maxn,nper,nsig,ndif,idif)
c
c.....Given two integer arrays o1(1:n) and o2(1:n) with their 
c     elements being also integer numbers (all of them different within
c     each of these two arrays) and provided that < o1(i) | o2(j) > = 1 
c     if o1(i).eq.o2(j) and < o1(i) | o2(j) > = 0 if o1(i).ne.o2(j), 
c     this routine re-orders the elements of the second array (o2) in 
c     such a way that the overlap matrix is diagonal. 
c
c.....There can not be repeated values in o1 and o2 arrays. However,
c     it is a responsability of the program that calls this routine
c     to test this point. 
c
c.....The input  of the routine is o1(i) and o2(i) arrays (i=1,n)
c     The output of the routine is the new   o2(i) array  (i=1,n),
c     nper = Number of permutations required to put back the final 
c            o2() values to the orignal orderging, and
c     nsig = +1 if nper is even
c            -1 if nper is odd.
c     ndif = Number of zeroes in the diagonal after re-ordering o2().
c     idif = indices of the NDIF elements.
c
c.....Example
c
c     Let o1() = ( 1  2  4  6  8  9 10 12 13 15 ) and 
c         o2() = ( 2  3  4  6  8 10 11 12 14 18 ) 
c
c     The original overlap matrix is
c
c     o1\o2 -->    2  3  4  6  8 10 11 12 14 18 
c              --------------------------------
c              1   0  0  0  0  0  0  0  0  0  0
c              2   1  0  0  0  0  0  0  0  0  0
c              4   0  0  1  0  0  0  0  0  0  0
c              6   0  0  0  1  0  0  0  0  0  0
c              8   0  0  0  0  1  0  0  0  0  0
c              9   0  0  0  0  0  0  0  0  0  0
c              10  0  0  0  0  0  1  0  0  0  0
c              12  0  0  0  0  0  0  0  1  0  0
c              13  0  0  0  0  0  0  0  0  0  0
c              15  0  0  0  0  0  0  0  0  0  0
c
c      Output o2 array is o2() = ( 3  2  4  6  8 11 10 12 14 18 )
c
c      that gives the following output overlap matrix.
c
c     o1\o2 -->    3  2  4  6  8 11 10 12 14 18 
c              --------------------------------
c              1   0  0  0  0  0  0  0  0  0  0
c              2   0  1  0  0  0  0  0  0  0  0
c              4   0  0  1  0  0  0  0  0  0  0
c              6   0  0  0  1  0  0  0  0  0  0
c              8   0  0  0  0  1  0  0  0  0  0
c              9   0  0  0  0  0  0  0  0  0  0
c              10  0  0  0  0  0  0  1  0  0  0
c              12  0  0  0  0  0  0  0  1  0  0
c              13  0  0  0  0  0  0  0  0  0  0
c              15  0  0  0  0  0  0  0  0  0  0
c-----------------------------------------------------------------------
c
      use mod_prec, only: ip
      implicit none
      integer(kind=ip) :: n,maxn,nper,nsig,k,j,i,ndif
      integer(kind=ip) :: o1(maxn),o2(maxn),idif(maxn)
c
      nper = 0
      nsig = 1
      do i = 1,n
        do j = 1,n
          if (i.ne.j.and.o2(j).eq.o1(i)) then
            nper = nper + 1
            nsig = -nsig
            k = o2(i)
            o2(i) = o2(j)
            o2(j) = k
            goto 1
          end if
        end do
 1    end do
      ndif = 0
      do i = 1,n
        if (o1(i).ne.o2(i)) then
          ndif = ndif + 1
          idif(ndif) = i
        end if
      end do
c
      end subroutine
