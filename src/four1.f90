subroutine four1(cdata,nn,isign)
  INTEGER isign,nn
  COMPLEX cdata(nn)
  REAL rdata(2*nn)
  INTEGER i,istep,j,m,mmax,n
  REAL tempi,tempr
  DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
  n=2*nn
!
! Pack complex data into real array
!      
  do i=1,n,2
    rdata(i)=real(cdata((i+1)/2))
    rdata(i+1)=imag(cdata((i+1)/2))
  end do
!
! Fourier transform
!
 
  j=1
  do i=1,n,2
    if (j.gt.i) then
      tempr=rdata(j)
      tempi=rdata(j+1)
      rdata(j)=rdata(i)
      rdata(j+1)=rdata(i+1)
      rdata(i)=tempr
      rdata(i+1)=tempi
    end if
    m=n/2
1   if ((m.ge.2).and.(j.gt.m)) then
      j=j-m
      m=m/2
      goto 1
    end if
    j=j+m
  end do

  mmax=2
2 if (n.gt.mmax) then
    istep=2*mmax
    theta=6.28318530717959d0/(isign*mmax)
    wpr=-2.d0*sin(0.5d0*theta)**2
    wpi=sin(theta)
    wr=1.d0
    wi=0.d0
    do m=1,mmax,2
      do i=m,n,istep
        j=i+mmax
        tempr=sngl(wr)*rdata(j)-sngl(wi)*rdata(j+1)
        tempi=sngl(wr)*rdata(j+1)+sngl(wi)*rdata(j)
        rdata(j)=rdata(i)-tempr
        rdata(j+1)=rdata(i+1)-tempi
        rdata(i)=rdata(i)+tempr
        rdata(i+1)=rdata(i+1)+tempi
      end do
      wtemp=wr
      wr=wr*wpr-wi*wpi+wr
      wi=wi*wpr+wtemp*wpi+wi
    end do
    mmax=istep
   goto 2
  end if
!
! Pack real data back into complex array
!      
  do i=1,nn
    cdata(i)=rdata(2*i-1) + (0.,1.)*rdata(2*i)
  end do
  return
end subroutine four1
