subroutine four1 (cdata,nn,isign)
!
!-----------------------------------------------------------------------
!
! ****** FOUR1.  Fourier transform.
!
!-----------------------------------------------------------------------
!
  use iso_fortran_env
!
!-----------------------------------------------------------------------
!
  implicit none
!
!-----------------------------------------------------------------------
!
  integer, parameter :: r_typ = REAL64
!
!-----------------------------------------------------------------------
!
  integer :: isign,nn
  complex(r_typ) :: cdata(nn)
  real(r_typ) :: rdata(2*nn)
  integer i,istep,j,m,mmax,n
  real(r_typ) :: tempi,tempr
  real(r_typ) :: theta,wi,wpi,wpr,wr,wtemp
!
!-----------------------------------------------------------------------
!  
  n = 2*nn
!
! ****** Pack complex data into real array
!
  do i=1,n,2
    rdata(i)   = real(cdata((i+1)/2),r_typ)
    rdata(i+1) = aimag(cdata((i+1)/2))
  enddo
!
! ****** Fourier transform
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
  enddo

  mmax = 2
!  
2 if (n.gt.mmax) then
    istep=2*mmax
    theta=6.28318530717959_r_typ/(isign*mmax)
    wpr=-2.0_r_typ*sin(0.5_r_typ*theta)**2
    wpi=sin(theta)
    wr=1.0_r_typ
    wi=0.
    do m=1,mmax,2
      do i=m,n,istep
        j=i+mmax
        tempr=real(wr,r_typ)*rdata(j)-real(wi,r_typ)*rdata(j+1)
        tempi=real(wr,r_typ)*rdata(j+1)+real(wi,r_typ)*rdata(j)
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
! ****** Pack real data back into complex array
!
  do i=1,nn
    cdata(i) = rdata(2*i-1) + (0.,1.0_r_typ)*rdata(2*i)
  enddo
!  
end subroutine four1
