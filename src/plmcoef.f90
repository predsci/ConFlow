subroutine plmcoef(lmaxp,coef)
!***********************************************************************
!
!  calculate coefficient matrix for generation of legendre polynomials
!  of degree l and order m in subroutine plm.
!
!***********************************************************************
  real*4 coef(lmaxp,lmaxp)

  lmax=lmaxp-1
  m=0
  x1=0.5
  coef(m+1,m+1)=sqrt(x1)
  do l=m+1,lmax
    coef(l+1,m+1)=sqrt(((2.*l+1.)/(l+m))*((2.*l-1.)/(l-m)))
    coef(m+1,l+1)=sqrt(((2.*l+1.)/(l+m))*((l+m-1.)/(l-m))*((l-m-1.)/(2.*l-3)))
  end do
  do m=1,lmax-1
    x1=x1*(2.*m+1.)/(2.*m)
    coef(m+1,m+1)=sqrt(x1)
    do l=m+1,lmax
      coef(l+1,m+1)=sqrt(((2.*l+1.)/(l+m))*((2.*l-1.)/(l-m)))
      coef(m+1,l+1)=sqrt(((2.*l+1.)/(l+m))*((l+m-1.)/(l-m))*((l-m-1.)/(2.*l-3)))
    enddo
  enddo
  m=lmax
  x1=x1*(2.*m+1.)/(2.*m)
  coef(m+1,m+1)=sqrt(x1)
  return
end subroutine plmcoef
