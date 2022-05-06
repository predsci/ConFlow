subroutine plm(m,x,lmaxp,coef,p)
!***********************************************************************
!
!  This subroutine calculates the associated legendre polynomial
!  using the recurrance relation for increasing degree l for given m
!  The coefficients for the relation are calculated once in
!  S/R PLMCOEF and passed as arguement coef.
!  The results, p(l,m,x), are returned in array p(l+1).
!
!***********************************************************************
  real*4 p(lmaxp),coef(lmaxp,lmaxp)
  integer lls
!***********************************************************************
!
!  lls is the index to the first non-zero p.                           *
!   mm .le. lls .le. lmaxp+1                                           *
!   If all p's are zero; i.e., .lt. 10**minlp, then lls=lmaxp+1        *
!
!***********************************************************************
  integer l,ll,mm,lp,lp1
  integer :: minlp = -20
  real    x2,alp,pa,pm1,pm2
  real, dimension(20) :: rlpa
  real :: r10p10 = 1.0e+10
  real :: r10m10 = 1.0e-10

  rlpa = (/1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, &
           1.0e-6, 1.0e-7, 1.0e-8, 1.0e-9, 1.0e-10,&
           1.0e-11, 1.0e-12, 1.0e-13, 1.0e-14, 1.0e-15,&
           1.0e-16, 1.0e-17, 1.0e-18, 1.0e-19, 1.0e-20/)

  mm=m+1
  lls=mm
  lmax=lmaxp-1

  do l=0,lmax
    p(l+1)=0.
  end do

  if (abs(x) .gt. 1.) then
    write(*,*) 'Legendre polynomial argument out of range.'
    return
  end if

  if (mm .gt. lmaxp) then
    write(*,*) 'Order of Legendre polynomial out of range.'
    return
  end if

  if (x .eq. 1. .and. m .eq. 0) then
    do l=0,lmax
      p(l+1)=sqrt(l+0.5)
    end do
    return
  end if
  if (x .eq. -1. .and. m .eq. 0) then
    do l=0,lmax
      isgn=1-2*mod(l,2)
      p(l+1)=isgn*sqrt(l+0.5)
    end do
    return
  end if

  if (abs(x) .eq. 1.)return

  x2=sqrt(1.-x*x)
  alp=alog10(coef(mm,mm))+float(m)*alog10(x2)
  lp=alp
  if (lp .lt. minlp) then
    alp=alp-float(lp)
    pa=10.0**alp
    lls=mm+1
  else
    pa=coef(mm,mm)*(x2**m)
    p(mm)=pa
  end if

  pm2=0.0
  pm1=pa
  if (lp .lt. minlp) then
    do ll=mm+1,lmaxp
      pa=coef(ll,mm)*x*pm1-coef(mm,ll)*pm2
      pm2=pm1
      pm1=pa
      if (pa .gt. r10p10) then
        pm2=pm2*r10m10
        pm1=pm1*r10m10
        lp=lp+10
        lp1=lp
        if (lp1 .ge. minlp) then
          pm2=pm2*rlpa(-lp1)
          pm1=pm1*rlpa(-lp1)
          lp=0
          go to 100
        end if
      end if
      lls=ll+1
    end do
  end if

100 continue

  if (lls .ge. lmaxp) return

  do ll=lls+1,lmaxp
    pa=coef(ll,mm)*x*pm1-coef(mm,ll)*pm2
    pm2=pm1
    pm1=pa
    p(ll)=pm1
  end do

  return
end subroutine plm
