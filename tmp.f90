
      allocate (f1(0:ndata),a1(0:ndata),b1(0:ndata),c1(0:ndata), stat=istat)
      if (istat /= 0) call error('atomin','could not allocate memory for f1, a1, b1, c1',faterr)
      h = 1d0 / (ndata+1)
      f1(0) = 0d0
      do j = 1,ndata
        q = h*j
        rdata = rmid*q/(1.d0-q)
        f1(j) = ftot(j,z(i))
      end do
      call spline(h,f1,a1,b1,c1,ndata,0.d0)
      atrho = 0d0
      do kk = 1, tpoints
        x = meshx(:,kk) - xyz(i,:)
        r = sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
        q = r/(r + rmid)
        intq = int((ndata+1)*q)
        dq = q - intq*h
        arho = abs((f1(intq)+dq*(a1(intq)+dq*(b1(intq)+dq*c1(intq)))))/(r*r)
        promol(kk) = promol(kk) + arho
        atrho = atrho + meshw(kk)*arho
      end do
      write (uout,'(1x,a,i3,1x,f12.8)') "*** Promolecular den for atom ", i, atrho
      deallocate (f1,a1,b1,c1)
    end do

    ! final stuff
    qpro = sum(meshw*promol)
