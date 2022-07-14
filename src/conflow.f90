!#######################################################################
!
!       ______            _______ _
!      / _____)          (_______) |
!     | /      ___  ____  _____  | | ___  _ _ _
!     | |     / _ \|  _ \|  ___) | |/ _ \| | | |
!     | \____| |_| | | | | |     | | |_| | | | |
!      \______)___/|_| |_|_|     |_|\___/ \____|
!
!
! ****** ConFlow: Super Granular Convective Flow Generator
!
!     This program creates 1024 by 512 velocity maps of supergranules
!     It constructs a spectrum of spherical harmonic amplitudes
!     s(l,m) and calculates the components of the velocity field at
!     an array of points in theta=colatitude and phi=longitude.
!
!     The velocity vector (u,v) is in the (phi,theta) direction.
!
!     The data arrays are written to direct access disk files.
!
!     Differential rotation is simulated by evolving
!     the spectral coefficients using 3-coefficient fits
!     (North-South symmetric)
!
!     Meridional flow is simulated by evolving
!     the spectral coefficients using 6-coefficient fits
!
! ****** Authors:  Raphael Attie
!                  Ronald M. Caplan
!                  David H. Hathaway
!
!     Code originally derived from SynchronicMapsSupergranulesV6.f90
!     by David H. Hathaway
!
!#######################################################################
! Copyright 2022.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
! implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!#######################################################################
!
!#######################################################################
module ident
!
!-----------------------------------------------------------------------
! ****** Set the name, version, and date of code.
!-----------------------------------------------------------------------
!
      character(*), parameter :: cname='Conflow'
      character(*), parameter :: cvers='0.4.0'
      character(*), parameter :: cdate='07/13/2022'
!
end module
!#######################################################################
module timing
!
      use number_types
!
      implicit none
!
      real(r_typ) :: wtime_tmp = 0.
!
      real(r_typ) :: wtime_io = 0.
!
      real(r_typ) :: wtime_total = 0.
!
end module
!#######################################################################
program conflow
!
!-----------------------------------------------------------------------
!
  use number_types
  use timing
  use ident
!
!-----------------------------------------------------------------------
!
  implicit none
!
!-----------------------------------------------------------------------
!
  integer,parameter :: nl=512
  integer,parameter :: nphi=1024
!
  character(8) :: fname
  character(5) :: ext
  character(100) :: path
  character(6) :: velFileNum
!
  integer :: i,j,l,m,l1,m1,l2,m2,jn,js,ierr
  integer :: lmax,lenPath,nlhalf,ifile,nfiles,itime
!
  real(r_typ) :: pi,root3,root5,root7,root9
  real(r_typ) :: rSun,daySec,deltaT,dphi,dtheta,theta,sintheta,x,rst,eo,v1,v2
  real(r_typ) :: tmax,curr_time
  real(r_typ) :: s0,s1,s2,s3,s4,s5,t0,t1,t2,t3,t4,el,em,amp2,amp3,randamp,phase,taper,xlifetime
  real(r_typ) :: elfunc,elfuncMF,DR0,DR1,DR2,DR3,DR4,MF0,MF1,MF2,MF3,MF4,MF5
  real(r_typ) :: BAALM4,BAALM3,BAALM2,BAALM1,BAALP0,BAALP1,BAALP2,BAALP3,BAALP4,BAALP5
  real(r_typ) :: BAA0,BAA1,BAA2,BAA3,BAA4,BAA5,BAA6
  real(r_typ) :: u(nphi,nl),v(nphi,nl)
  real(r_typ) :: coef(nl,nl),p(nl)
  real(r_typ) :: Omega(nl,nl)
  real(r_typ) :: OmegaM1(nl,nl),OmegaM2(nl,nl),OmegaM3(nl,nl),OmegaM4(nl,nl)
  real(r_typ) :: OmegaP1(nl,nl),OmegaP2(nl,nl),OmegaP3(nl,nl),OmegaP4(nl,nl)
  real(r_typ) :: MFlow0(nl,nl)
  real(r_typ) :: MFlowM1(nl,nl),MFlowM2(nl,nl),MFlowM3(nl,nl),MFlowM4(nl,nl)
  real(r_typ) :: MFlowP1(nl,nl),MFlowP2(nl,nl),MFlowP3(nl,nl),MFlowP4(nl,nl)
  real(r_typ) :: MFlowP5(nl,nl),MFlowM5(nl,nl),MFlowP6(nl,nl),MFlowM6(nl,nl)
  real(r_typ) :: coefDR(5),coefMF(6)
!
  real(r_typ), dimension(:), allocatable :: tvec,pvec
  real(r_typ) :: wtime,wt1
  real(r_typ) :: rnd
!
  complex(r_typ) :: s(nl,nl),ds(nl,nl,4)
  complex(r_typ) :: t(nl,nl),dt(nl,nl,4)
  complex(r_typ) :: xi,arg,sum1,sum2,sum3,sum4
  complex(r_typ) :: unorth(nphi),usouth(nphi)
  complex(r_typ) :: vnorth(nphi),vsouth(nphi)
!
  integer :: seed_n
  integer, allocatable, dimension(:) :: seed_old,seed_new
!
!-----------------------------------------------------------------------
!
! ****** Start wall clock timer.
!
  wtime_tmp = wtime()
!
!-----------------------------------------------------------------------
!
  write(*,*)
  write(*,*) cname,' v',cvers,' ',cdate
  write(*,*)
!
  path= './'
  fname=' '
  ext='.data'
  lenPath=len(trim(path))
!
  lmax=nl-1
  nlhalf=nl/2
!
! ****** Initialize random number generator.
! ****** Set a specific seed for reproducibility.
! ****** Later this should be made into an option!!!! [RMC]
!
  call RANDOM_INIT (.true., .true.)
  call RANDOM_SEED (size=seed_n)
  allocate (seed_old(seed_n))
  allocate (seed_new(seed_n))
  call RANDOM_SEED(get=seed_old)
  seed_new(:)=12345
  call RANDOM_SEED(put=seed_new)
!
!c***********************************************************************
!c                                                                      *
!c  Numeric constants
!c                                                                      *
!c***********************************************************************
  pi=4.0_r_typ*atan(1.0_r_typ)
  root3=sqrt(3.0_r_typ)
  root5=sqrt(5.0_r_typ)
  root7=sqrt(7.0_r_typ)
  root9=sqrt(9.0_r_typ)
  xi=(0.,1.0_r_typ)
!c***********************************************************************
!c                                                                      *
!c  Physical constants
!c                                                                      *
!c***********************************************************************
  rSun = 6.96e+08
  daySec = 86400.0_r_typ
!c***********************************************************************
!c                                                                      *
!c  Time (seconds) and space (radians) steps
!c                                                                      *
!c***********************************************************************
  deltaT = 15.0_r_typ*60.0_r_typ
  tmax = 28*24*3600
  curr_time = 0.
  nfiles = int(tmax/deltaT)
  dphi = 2.0_r_typ*pi/nphi
  dtheta = pi/nl
!
! ****** Set up theta (co-lat) and phi 1D scale arrays.
!
  allocate(tvec(nl))
  allocate(pvec(nphi))

  do i=1,nl
    tvec(i)=dtheta*0.5_r_typ+(i-1)*dtheta
  enddo
!
  do i=1,nphi
    pvec(i)=dphi*0.5_r_typ+(i-1)*dphi
  enddo
!
!***********************************************************************
!c                                                                      *
!c  Differential Rotation coefficients m/s relative to Carrington
!c                                                                      *
!c***********************************************************************
2 format(5(4x,f8.3))
3 format(6(4x,f8.3))
  open(unit=2,file=path(1:lenPath) // 'conflow_flow_parameters_v6.txt',status='old')
    read(2,2) t0,t1,t2,t3,t4
    write(*,*)
    write(*,*) 'Differential rotation coeffs (m/s) t0,t1,t2,t3,t4:'
    write(*,*) t0,t1,t2,t3,t4
    coefDR(1)=t0*deltaT/rSun
    coefDR(2)=t1*deltaT/rSun
    coefDR(3)=t2*deltaT/rSun
    coefDR(4)=t3*deltaT/rSun
    coefDR(5)=t4*deltaT/rSun
!c***********************************************************************
!c
!c Meridional flow speed of supergranules in m/s poleward
!c
!c***********************************************************************
    read(2,3) s0,s1,s2,s3,s4,s5
    write(*,*)
    write(*,*) 'Differential rotation coeffs (m/s) s0,s1,s2,s3,s4,s5:'
    write(*,*) s0,s1,s2,s3,s4,s5
    coefMF(1)=s0*deltaT/rSun
    coefMF(2)=s1*deltaT/rSun
    coefMF(3)=s2*deltaT/rSun
    coefMF(4)=s3*deltaT/rSun
    coefMF(5)=s4*deltaT/rSun
    coefMF(6)=s5*deltaT/rSun
  close(2)
!c***********************************************************************
!c
!c  Create coefficients for spherical harmonics.
!c
!c***********************************************************************
  call plmcoef(nl,coef)
!
  write(*,*)
  write(*,*) 'Legendre recurrance coefficients calculated.'
!c***********************************************************************
!c                                                                      *
!c  Construct the convection spectrum.                                  *
!c                                                                      *
!c***********************************************************************
  do l=1,lmax
    l1=l+1
    el=real(l,r_typ)
    amp2 = 0.08*(1. - tanh(el/165.)) + 0.0024*(1. - tanh(el/2000.))
    amp3 = 1.5*(1. - 0.5*sqrt(el/1000.))/el
    taper=1.0
    if (l .gt. 384) taper=0.5_r_typ*(1.0_r_typ + cos(pi*(l-384.0_r_typ)/(512.0_r_typ-384.0_r_typ)))
    amp2 = taper*amp2
    amp3 = taper*amp3
    do m=1,l
      m1=m+1
      call RANDOM_NUMBER (rnd)
      phase=2.0_r_typ*pi*rnd
      arg=cos(phase)+xi*sin(phase)
      call RANDOM_NUMBER (rnd)
      randamp=1.8*rnd
      s(l1,m1)=randamp*amp2*arg
      call RANDOM_NUMBER (rnd)
      phase=2.*pi*rnd
      arg=cos(phase)+xi*sin(phase)
      call RANDOM_NUMBER (rnd)
      randamp=1.8*rnd
      t(l1,m1)=randamp*amp3*arg
    end do
  end do
  write(*,*)
  write(*,*) 'Velocity spectrum calculated'
!c***********************************************************************
!c
!c Coupling Coefficients
!c
!c***********************************************************************
  do m=1,lmax
    m1=m+1
    em=real(m,r_typ)
    do l=m,lmax
      l1=l+1
      el=real(l,r_typ)
!c
!c Modify differential rotation and meridional flow with l (depth)
!c
      elfunc=1.0
      elfuncMF=1.0
      DR0=elfunc*coefDR(1)
      DR1=elfunc*coefDR(2)
      DR2=elfunc*coefDR(3)
      DR3=elfunc*coefDR(4)
      DR4=elfunc*coefDR(5)
      MF0=elfuncMF*coefMF(1)
      MF1=elfuncMF*coefMF(2)
      MF2=elfuncMF*coefMF(3)
      MF3=elfuncMF*coefMF(4)
      MF4=elfuncMF*coefMF(5)
      MF5=elfuncMF*coefMF(6)
!
! Common factors in calculations - A(l,m)^-2
!
      BAALM4=(el+em-4.)*(el-em-4.)/((2.*el-7.)*(2.*el-9.))
      BAALM3=(el+em-3.)*(el-em-3.)/((2.*el-5.)*(2.*el-7.))
      BAALM2=(el+em-2.)*(el-em-2.)/((2.*el-3.)*(2.*el-5.))
      BAALM1=(el+em-1.)*(el-em-1.)/((2.*el-1.)*(2.*el-3.))
      BAALP0=(el+em+0.)*(el-em+0.)/((2.*el+1.)*(2.*el-1.))
      BAALP1=(el+em+1.)*(el-em+1.)/((2.*el+3.)*(2.*el+1.))
      BAALP2=(el+em+2.)*(el-em+2.)/((2.*el+5.)*(2.*el+3.))
      BAALP3=(el+em+3.)*(el-em+3.)/((2.*el+7.)*(2.*el+5.))
      BAALP4=(el+em+4.)*(el-em+4.)/((2.*el+9.)*(2.*el+7.))
      BAALP5=(el+em+5.)*(el-em+5.)/((2.*el+11.)*(2.*el+9.))
!
!  Coupling with l+6 component
!
      if ((l+6) .lt. lmax) then
        BAA0=1./(coef(l1+6,m1)*coef(l1+5,m1)*coef(l1+4,m1)*coef(l1+3,m1)*coef(l1+2,m1)*coef(l1+1,m1))
        MFlowP6(l1,m1)=-(el+7.)*BAA0*MF5
      end if
!
!  Coupling with l+5 component
!
      if ((l+5) .lt. lmax) then
        BAA0=1./(coef(l1+5,m1)*coef(l1+4,m1)*coef(l1+3,m1)*coef(l1+2,m1)*coef(l1+1,m1))
        MFlowP5(l1,m1)=-(el+6.)*BAA0*MF4
      end if
!
!  Coupling with l+4 component
!
      if ((l+4) .lt. lmax) then
        BAA0=1./(coef(l1+4,m1)*coef(l1+3,m1)*coef(l1+2,m1)*coef(l1+1,m1))
        OmegaP4(l1,m1)=em*BAA0*DR4

        BAA1=BAA0*BAALP5
        BAA2=BAA0*(BAALP4 + BAALP3 + BAALP2 + BAALP1 + BAALP0)
        MFlowP4(l1,m1)=-(el+5.)*BAA0*MF3 + ((el+4.)*BAA1-(el+5.)*BAA2)*MF5
      end if
!
!  Coupling with l+3 component
!
      if ((l+3) .lt. lmax) then
        BAA0=1./(coef(l1+3,m1)*coef(l1+2,m1)*coef(l1+1,m1))
        OmegaP3(l1,m1)=em*BAA0*DR3

        BAA1=BAA0*BAALP4
        BAA2=BAA0*(BAALP3 + BAALP2 + BAALP1 + BAALP0)
        MFlowP3(l1,m1)=-(el+4.)*BAA0*MF2 + ((el+3.)*BAA1-(el+4.)*BAA2)*MF4
      end if
!
!  Coupling with l+2 component
!
      if ((l+2) .lt. lmax) then
        BAA0=1./(coef(l1+2,m1)*coef(l1+1,m1))
        BAA1=BAA0*(BAALP3 + BAALP2 + BAALP1 + BAALP0)
        OmegaP2(l1,m1)=em*(BAA0*DR2+BAA1*DR4)

        BAA1=BAA0*BAALP3
        BAA2=BAA0*(BAALP2 + BAALP1 + BAALP0)
        BAA3=BAA0*BAALP3*(BAALP4+BAALP3+BAALP2+BAALP1+BAALP0)
        BAA4=BAA0*(BAALP2*(BAALP3+BAALP2+BAALP1) + BAALP1*(BAALP2+BAALP1+BAALP0) + BAALP0*(BAALP2+BAALP1+BAALP0+BAALM1))
        MFlowP2(l1,m1)=-(el+3.)*BAA0*MF1 + ((el+2.)*BAA1-(el+3.)*BAA2)*MF3 + ((el+2.)*BAA3-(el+3.)*BAA4)*MF5
      end if
!
!  Coupling with l+1 component
!
      if ((l+1) .lt. lmax) then
        BAA0=1./coef(l1+1,m1)
        BAA1=BAA0*(BAALP2 + BAALP1 + BAALP0)
        OmegaP1(l1,m1)=em*(BAA0*DR1+BAA1*DR3)

!        BAA1=BAA0*BAALP0  ! This seem to be in error
        BAA1=BAA0*BAALP2
        BAA2=BAA0*(BAALP1+BAALP0)
        BAA3=BAA0*BAALP2*(BAALP3+BAALP2+BAALP1+BAALP0)
        BAA4=BAA0*(BAALP1*(BAALP2+BAALP1+BAALP0) + BAALP0*(BAALP1+BAALP0+BAALM1))
        MFlowP1(l1,m1)=-(el+2.)*BAA0*MF0 + ((el+1.)*BAA1-(el+2.)*BAA2)*MF2 + ((el+1.)*BAA3-(el+2.)*BAA4)*MF4
      end if
!c
!c Coupling with l component
!c
      BAA2=BAALP1+BAALP0
      BAA4=BAALP1*(BAALP2+BAALP1+BAALP0) + BAALP0*(BAALP1+BAALP0+BAALM1)
      Omega(l1,m1)=em*(DR0+BAA2*DR2+BAA4*DR4)

      BAA1=BAALP1
      BAA2=BAALP0
      BAA3=BAALP1*(BAALP2+BAALP1+BAALP0)
      BAA4=BAALP0*(BAALP1+BAALP0+BAALM1)
      BAA5=BAALP1*(BAALP2*(BAALP3+BAALP2+BAALP1) + BAALP1*(BAALP2+BAALP1+BAALP0) + BAALP0*(BAALP2+BAALP1+BAALP0+BAALM1))
      BAA6=BAALP0*(BAALP1*(BAALP2+BAALP1+BAALP0+BAALM1) + BAALP0*(BAALP1+BAALP0+BAALM1) + BAALM1*(BAALP0+BAALM1+BAALM2))
      MFlow0(l1,m1)=(el*BAA1-(el+1.)*BAA2)*MF1 + (el*BAA3-(el+1.)*BAA4)*MF3 + (el*BAA5-(el+1.)*BAA6)*MF5
!c
!c  Coupling with l-1 component
!c
      if ((l-1) .ge. m) then
        BAA0=1./coef(l1,m1)
        BAA1=BAA0*(BAALP1 + BAALP0 + BAALM1)
        OmegaM1(l1,m1)=em*(BAA0*DR1+BAA1*DR3)

        BAA1=BAA0*(BAALP1 + BAALP0)
        BAA2=BAA0*BAALM1
        BAA3=BAA0*(BAALP1*(BAALP2+BAALP1+BAALP0) + BAALP0*(BAALP1+BAALP0+BAALM1))
        BAA4=BAA0*BAALM1*(BAALP1+BAALP0+BAALM1+BAALM2)
        MFlowM1(l1,m1)=(el-1.)*BAA0*MF0 + ((el-1.)*BAA1-el*BAA2)*MF2 + ((el-1.)*BAA3-el*BAA4)*MF4
      end if
!c
!c  Coupling with l-2 component
!c
      if ((l-2) .ge. m) then
        BAA0=1./(coef(l1,m1)*coef(l1-1,m1))
        BAA1=BAA0*(BAALP1+BAALP0+BAALM1+BAALM2)
        OmegaM2(l1,m1)=em*(BAA0*DR2+BAA1*DR4)

        BAA2=BAA0*(BAALP1+BAALP0+BAALM1)
        BAA3=BAA0*BAALM2
        BAA4=BAA0*(BAALP1*(BAALP2+BAALP1+BAALP0+BAALM1) + BAALP0*(BAALP1+BAALP0+BAALM1) + BAALM1*(BAALP0+BAALM1+BAALM2))
        BAA5=BAA0*BAALM2*(BAALP1+BAALP0+BAALM1+BAALM2+BAALM3)
        MFlowM2(l1,m1)=(el-2.)*BAA0*MF1 + ((el-2.)*BAA2-(el-1.)*BAA3)*MF3 + ((el-2.)*BAA4-(el-1.)*BAA5)*MF5
      end if
!c
!c  Coupling with l-3 component
!c
      if ((l-3) .ge. m) then
        BAA0=1./(coef(l1,m1)*coef(l1-1,m1)*coef(l1-2,m1))
        OmegaM3(l1,m1)=em*BAA0*DR3

        BAA2=BAA0*(BAALP1 + BAALP0 + BAALM1 + BAALM2)
        BAA3=BAA0*BAALM3
        MFlowM3(l1,m1)=(el-3.)*BAA0*MF2 + ((el-3.)*BAA2-(el-2.)*BAA3)*MF4
      end if
!c
!c  Coupling with l-4 component
!c
      if ((l-4) .ge. m) then
        BAA0=1./(coef(l1,m1)*coef(l1-1,m1)*coef(l1-2,m1)*coef(l1-3,m1))
        OmegaM4(l1,m1)=em*BAA0*DR4

        BAA1=BAA0*(BAALP1 + BAALP0 + BAALM1 + BAALM2)
        BAA2=BAA0*BAALM4
        MFlowM4(l1,m1)=(el-4.)*BAA0*MF3 + ((el-4.)*BAA1-(el-3.)*BAA2)*MF5
      end if
!c
!c  Coupling with l-5 component
!c
      if ((l-5) .ge. m) then
        BAA0=1./(coef(l1,m1)*coef(l1-1,m1)*coef(l1-2,m1)*coef(l1-3,m1)*coef(l1-4,m1))
        MFlowM5(l1,m1)=(el-5.)*BAA0*MF4
      end if
!c
!c  Coupling with l-6 component
!c
      if ((l-6) .ge. m) then
        BAA0=1./(coef(l1,m1)*coef(l1-1,m1)*coef(l1-2,m1)*coef(l1-3,m1)*coef(l1-4,m1)*coef(l1-5,m1))
        MFlowM6(l1,m1)=(el-6.)*BAA0*MF5
      end if
    end do ! End of m-loop
  end do ! End of l=loop
!c***********************************************************************
!c
!c  Construct series of velocity maps at deltaT (in seconds) intervals
!c
!c***********************************************************************
!
  write(*,*)
  write(*,*) 'A time step of ', deltaT, 'seconds will produce ', &
             nfiles,' vt and vp flow files.'
  write(*,*)
  flush(OUTPUT_UNIT)
!
! ****** Open text file to list output maps.
!
  call ffopen (12,'flow_output_list.csv','rw',ierr)
  write (12,'(A11,A1,A11,A1,A11)') &
        'TIME(JD)',',','VTFILENAME',',','VPFILENAME'
  close(12)
!
!****** START MAIN LOOP. ******
!
  do itime=1,nfiles
    ifile=itime
!
!c***********************************************************************
!c
!c  Evolve spectral coefficients using 4th order Runga-Kutta
!c
!c Step 1
!c
!c***********************************************************************
    do m=1,lmax-1
      m1=m+1
      do l=m,lmax-1
        l1=l+1
        ds(l1,m1,1)=0.
        dt(l1,m1,1)=0.
!c
!c  contribution from l+6 component
!c
        if ((l+6) .lt. lmax-1) then
          ds(l1,m1,1)=ds(l1,m1,1)+MFlowP6(l1,m1)*s(l1+6,m1)
          dt(l1,m1,1)=dt(l1,m1,1)+MFlowP6(l1,m1)*t(l1+6,m1)
        end if
!c
!c  contribution from l+5 component
!c
        if ((l+5) .lt. lmax-1) then
          ds(l1,m1,1)=ds(l1,m1,1)+MFlowP5(l1,m1)*s(l1+5,m1)
          dt(l1,m1,1)=dt(l1,m1,1)+MFlowP5(l1,m1)*t(l1+5,m1)
        end if
!c
!c  contribution from l+4 component
!c
        if ((l+4) .lt. lmax-1) then
          ds(l1,m1,1)=ds(l1,m1,1)+(xi*OmegaP4(l1,m1)+MFlowP4(l1,m1))*s(l1+4,m1)
          dt(l1,m1,1)=dt(l1,m1,1)+(xi*OmegaP4(l1,m1)+MFlowP4(l1,m1))*t(l1+4,m1)
         end if
!c
!c  contribution from l+3 component
!c
        if ((l+3) .lt. lmax-1) then
          ds(l1,m1,1)=ds(l1,m1,1)+(xi*OmegaP3(l1,m1)+MFlowP3(l1,m1))*s(l1+3,m1)
          dt(l1,m1,1)=dt(l1,m1,1)+(xi*OmegaP3(l1,m1)+MFlowP3(l1,m1))*t(l1+3,m1)
        end if
!c
!c  contribution from l+2 component
!c
        if ((l+2) .lt. lmax-1) then
         ds(l1,m1,1)=ds(l1,m1,1)+(xi*OmegaP2(l1,m1)+MFlowP2(l1,m1))*s(l1+2,m1)
         dt(l1,m1,1)=dt(l1,m1,1)+(xi*OmegaP2(l1,m1)+MFlowP2(l1,m1))*t(l1+2,m1)
        end if
!c
!c  contribution from l+1 component
!c
        if ((l+1) .lt. lmax-1) then
          ds(l1,m1,1)=ds(l1,m1,1)+(xi*OmegaP1(l1,m1)+MFlowP1(l1,m1))*s(l1+1,m1)
          dt(l1,m1,1)=dt(l1,m1,1)+(xi*OmegaP1(l1,m1)+MFlowP1(l1,m1))*t(l1+1,m1)
        end if
!c
!c  contribution from l component
!c
        ds(l1,m1,1)=ds(l1,m1,1)+(0.0*xi*Omega(l1,m1)+MFlow0(l1,m1))*s(l1,m1) ! No Omega(l1,m1) term?
        dt(l1,m1,1)=dt(l1,m1,1)+(0.0*xi*Omega(l1,m1)+MFlow0(l1,m1))*t(l1,m1) ! No Omega(l1,m1) term?
!c
!c  contribution from l-1 component
!c
        if ((l-1) .ge. m) then
          ds(l1,m1,1)=ds(l1,m1,1)+(xi*OmegaM1(l1,m1)+MFlowM1(l1,m1))*s(l1-1,m1)
          dt(l1,m1,1)=dt(l1,m1,1)+(xi*OmegaM1(l1,m1)+MFlowM1(l1,m1))*t(l1-1,m1)
        end if
!c
!c  contribution from l-2 component
!c
        if ((l-2) .ge. m) then
          ds(l1,m1,1)=ds(l1,m1,1)+(xi*OmegaM2(l1,m1)+MFlowM2(l1,m1))*s(l1-2,m1)
          dt(l1,m1,1)=dt(l1,m1,1)+(xi*OmegaM2(l1,m1)+MFlowM2(l1,m1))*t(l1-2,m1)
        end if
!c
!c  contribution from l-3 component
!c
        if ((l-3) .ge. m) then
          ds(l1,m1,1)=ds(l1,m1,1)+(xi*OmegaM3(l1,m1)+MFlowM3(l1,m1))*s(l1-3,m1)
          dt(l1,m1,1)=dt(l1,m1,1)+(xi*OmegaM3(l1,m1)+MFlowM3(l1,m1))*t(l1-3,m1)
        end if
!c
!c  contribution from l-4 component
!c
        if ((l-4) .ge. m) then
          ds(l1,m1,1)=ds(l1,m1,1)+(xi*OmegaM4(l1,m1)+MFlowM4(l1,m1))*s(l1-4,m1)
          dt(l1,m1,1)=dt(l1,m1,1)+(xi*OmegaM4(l1,m1)+MFlowM4(l1,m1))*t(l1-4,m1)
        end if
!c
!c  contribution from l-5 component
!c
        if ((l-5) .ge. m) then
          ds(l1,m1,1)=ds(l1,m1,1)+MFlowM5(l1,m1)*s(l1-5,m1)
          dt(l1,m1,1)=dt(l1,m1,1)+MFlowM5(l1,m1)*t(l1-5,m1)
        end if
!c
!c  contribution from l-6 component
!c
        if ((l-6) .ge. m) then
          ds(l1,m1,1)=ds(l1,m1,1)+MFlowM6(l1,m1)*s(l1-6,m1)
          dt(l1,m1,1)=dt(l1,m1,1)+MFlowM6(l1,m1)*t(l1-6,m1)
        end if
      end do ! End l-loop step1
    end do ! End m-loop step1

!c***********************************************************************
!c
!c Step 2
!c
!c***********************************************************************
    do m=1,lmax-1
      m1=m+1
      do l=m,lmax-1
        l1=l+1
        ds(l1,m1,2)=0.
        dt(l1,m1,2)=0.
!c
!c  contribution from l+6 component
!c
        if ((l+6) .lt. lmax-1) then
          ds(l1,m1,2)=ds(l1,m1,2)+MFlowP6(l1,m1)*(s(l1+6,m1)+ds(l1+6,m1,1)/2.)*exp(xi*Omega(l1+6,m1)/2.)
          dt(l1,m1,2)=dt(l1,m1,2)+MFlowP6(l1,m1)*(t(l1+6,m1)+dt(l1+6,m1,1)/2.)*exp(xi*Omega(l1+6,m1)/2.)
        end if
!c
!c  contribution from l+5 component
!c
        if ((l+5) .lt. lmax-1) then
          ds(l1,m1,2)=ds(l1,m1,2)+MFlowP5(l1,m1)*(s(l1+5,m1)+ds(l1+5,m1,1)/2.)*exp(xi*Omega(l1+5,m1)/2.)
          dt(l1,m1,2)=dt(l1,m1,2)+MFlowP5(l1,m1)*(t(l1+5,m1)+dt(l1+5,m1,1)/2.)*exp(xi*Omega(l1+5,m1)/2.)
        end if
!c
!c  contribution from l+4 component
!c
        if ((l+4) .lt. lmax-1) then
          ds(l1,m1,2)=ds(l1,m1,2)+(xi*OmegaP4(l1,m1)+MFlowP4(l1,m1))*(s(l1+4,m1)+ds(l1+4,m1,1)/2.)*exp(xi*Omega(l1+4,m1)/2.)
          dt(l1,m1,2)=dt(l1,m1,2)+(xi*OmegaP4(l1,m1)+MFlowP4(l1,m1))*(t(l1+4,m1)+dt(l1+4,m1,1)/2.)*exp(xi*Omega(l1+4,m1)/2.)
        end if
!c
!c  contribution from l+3 component
!c
        if ((l+3) .lt. lmax-1) then
          ds(l1,m1,2)=ds(l1,m1,2)+(xi*OmegaP3(l1,m1)+MFlowP3(l1,m1))*(s(l1+3,m1)+ds(l1+3,m1,1)/2.)*exp(xi*Omega(l1+3,m1)/2.)
          dt(l1,m1,2)=dt(l1,m1,2)+(xi*OmegaP3(l1,m1)+MFlowP3(l1,m1))*(t(l1+3,m1)+dt(l1+3,m1,1)/2.)*exp(xi*Omega(l1+3,m1)/2.)
        end if
!c
!c  contribution from l+2 component
!c
        if ((l+2) .lt. lmax-1) then
          ds(l1,m1,2)=ds(l1,m1,2)+(xi*OmegaP2(l1,m1)+MFlowP2(l1,m1))*(s(l1+2,m1)+ds(l1+2,m1,1)/2.)*exp(xi*Omega(l1+2,m1)/2.)
          dt(l1,m1,2)=dt(l1,m1,2)+(xi*OmegaP2(l1,m1)+MFlowP2(l1,m1))*(t(l1+2,m1)+dt(l1+2,m1,1)/2.)*exp(xi*Omega(l1+2,m1)/2.)
        end if
!c
!c  contribution from l+1 component
!c
        if ((l+1) .lt. lmax-1) then
          ds(l1,m1,2)=ds(l1,m1,2)+(xi*OmegaP1(l1,m1)+MFlowP1(l1,m1))*(s(l1+1,m1)+ds(l1+1,m1,1)/2.)*exp(xi*Omega(l1+1,m1)/2.)
          dt(l1,m1,2)=dt(l1,m1,2)+(xi*OmegaP1(l1,m1)+MFlowP1(l1,m1))*(t(l1+1,m1)+dt(l1+1,m1,1)/2.)*exp(xi*Omega(l1+1,m1)/2.)
        end if
!c
!c  contribution from l component
!c
        ds(l1,m1,2)=ds(l1,m1,2)+MFlow0(l1,m1)*(s(l1,m1)+ds(l1,m1,1)/2.)*exp(xi*Omega(l1,m1)/2.)
        dt(l1,m1,2)=dt(l1,m1,2)+MFlow0(l1,m1)*(t(l1,m1)+dt(l1,m1,1)/2.)*exp(xi*Omega(l1,m1)/2.)
!c
!c  contribution from l-1 component
!c
        if ((l-1) .ge. m) then
          ds(l1,m1,2)=ds(l1,m1,2)+(xi*OmegaM1(l1,m1)+MFlowM1(l1,m1))*(s(l1-1,m1)+ds(l1-1,m1,1)/2.)*exp(xi*Omega(l1-1,m1)/2.)
          dt(l1,m1,2)=dt(l1,m1,2)+(xi*OmegaM1(l1,m1)+MFlowM1(l1,m1))*(t(l1-1,m1)+dt(l1-1,m1,1)/2.)*exp(xi*Omega(l1-1,m1)/2.)
        end if
!c
!c  contribution from l-2 component
!c
        if ((l-2) .ge. m) then
          ds(l1,m1,2)=ds(l1,m1,2)+(xi*OmegaM2(l1,m1)+MFlowM2(l1,m1))*(s(l1-2,m1)+ds(l1-2,m1,1)/2.)*exp(xi*Omega(l1-2,m1)/2.)
          dt(l1,m1,2)=dt(l1,m1,2)+(xi*OmegaM2(l1,m1)+MFlowM2(l1,m1))*(t(l1-2,m1)+dt(l1-2,m1,1)/2.)*exp(xi*Omega(l1-2,m1)/2.)
        end if
!c
!c  contribution from l-3 component
!c
        if ((l-3) .ge. m) then
          ds(l1,m1,2)=ds(l1,m1,2)+(xi*OmegaM3(l1,m1)+MFlowM3(l1,m1))*(s(l1-3,m1)+ds(l1-3,m1,1)/2.)*exp(xi*Omega(l1-3,m1)/2.)
          dt(l1,m1,2)=dt(l1,m1,2)+(xi*OmegaM3(l1,m1)+MFlowM3(l1,m1))*(t(l1-3,m1)+dt(l1-3,m1,1)/2.)*exp(xi*Omega(l1-3,m1)/2.)
        end if
!c
!c  contribution from l-4 component
!c
        if ((l-4) .ge. m) then
          ds(l1,m1,2)=ds(l1,m1,2)+(xi*OmegaM4(l1,m1)+MFlowM4(l1,m1))*(s(l1-4,m1)+ds(l1-4,m1,1)/2.)*exp(xi*Omega(l1-4,m1)/2.)
          dt(l1,m1,2)=dt(l1,m1,2)+(xi*OmegaM4(l1,m1)+MFlowM4(l1,m1))*(t(l1-4,m1)+dt(l1-4,m1,1)/2.)*exp(xi*Omega(l1-4,m1)/2.)
        end if
!c
!c  contribution from l-5 component
!c
        if ((l-5) .ge. m) then
          ds(l1,m1,2)=ds(l1,m1,2)+MFlowM5(l1,m1)*(s(l1-5,m1)+ds(l1-5,m1,1)/2.)*exp(xi*Omega(l1-5,m1)/2.)
          dt(l1,m1,2)=dt(l1,m1,2)+MFlowM5(l1,m1)*(t(l1-5,m1)+dt(l1-5,m1,1)/2.)*exp(xi*Omega(l1-5,m1)/2.)
        end if
!c
!c  contribution from l-6 component
!c
        if ((l-6) .ge. m) then
          ds(l1,m1,2)=ds(l1,m1,2)+MFlowM6(l1,m1)*(s(l1-6,m1)+ds(l1-6,m1,1)/2.)*exp(xi*Omega(l1-6,m1)/2.)
          dt(l1,m1,2)=dt(l1,m1,2)+MFlowM6(l1,m1)*(t(l1-6,m1)+dt(l1-6,m1,1)/2.)*exp(xi*Omega(l1-6,m1)/2.)
        end if
      end do ! End l-loop Step2
    end do ! End m-loop Step2

!c***********************************************************************
!c
!c Step 3
!c
!c***********************************************************************
    do m=1,lmax-1
      m1=m+1
      do l=m,lmax-1
        l1=l+1
        ds(l1,m1,3)=0.
        dt(l1,m1,3)=0.
!c
!c  contribution from l+6 component
!c
        if ((l+6) .lt. lmax-1) then
          ds(l1,m1,3)=ds(l1,m1,3)+MFlowP6(l1,m1)*(s(l1+6,m1)+ds(l1+6,m1,2)/2.)*exp(xi*Omega(l1+6,m1)/2.)
          dt(l1,m1,3)=dt(l1,m1,3)+MFlowP6(l1,m1)*(t(l1+6,m1)+dt(l1+6,m1,2)/2.)*exp(xi*Omega(l1+6,m1)/2.)
        end if
!c
!c  contribution from l+5 component
!c
        if ((l+5) .lt. lmax-1) then
         ds(l1,m1,3)=ds(l1,m1,3)+MFlowP5(l1,m1)*(s(l1+5,m1)+ds(l1+5,m1,2)/2.)*exp(xi*Omega(l1+5,m1)/2.)
         dt(l1,m1,3)=dt(l1,m1,3)+MFlowP5(l1,m1)*(t(l1+5,m1)+dt(l1+5,m1,2)/2.)*exp(xi*Omega(l1+5,m1)/2.)
        end if
!c
!c  contribution from l+4 component
!c
        if ((l+4) .lt. lmax-1) then
          ds(l1,m1,3)=ds(l1,m1,3)+(xi*OmegaP4(l1,m1)+MFlowP4(l1,m1))*(s(l1+4,m1)+ds(l1+4,m1,2)/2.)*exp(xi*Omega(l1+4,m1)/2.)
          dt(l1,m1,3)=dt(l1,m1,3)+(xi*OmegaP4(l1,m1)+MFlowP4(l1,m1))*(t(l1+4,m1)+dt(l1+4,m1,2)/2.)*exp(xi*Omega(l1+4,m1)/2.)
        end if
!c
!c  contribution from l+3 component
!c
        if ((l+3) .lt. lmax-1) then
          ds(l1,m1,3)=ds(l1,m1,3)+(xi*OmegaP3(l1,m1)+MFlowP3(l1,m1))*(s(l1+3,m1)+ds(l1+3,m1,2)/2.)*exp(xi*Omega(l1+3,m1)/2.)
          dt(l1,m1,3)=dt(l1,m1,3)+(xi*OmegaP3(l1,m1)+MFlowP3(l1,m1))*(t(l1+3,m1)+dt(l1+3,m1,2)/2.)*exp(xi*Omega(l1+3,m1)/2.)
        end if
!c
!c  contribution from l+2 component
!c
        if ((l+2) .lt. lmax-1) then
          ds(l1,m1,3)=ds(l1,m1,3)+(xi*OmegaP2(l1,m1)+MFlowP2(l1,m1))*(s(l1+2,m1)+ds(l1+2,m1,2)/2.)*exp(xi*Omega(l1+2,m1)/2.)
          dt(l1,m1,3)=dt(l1,m1,3)+(xi*OmegaP2(l1,m1)+MFlowP2(l1,m1))*(t(l1+2,m1)+dt(l1+2,m1,2)/2.)*exp(xi*Omega(l1+2,m1)/2.)
        end if
!c
!c  contribution from l+1 component
!c
        if ((l+1) .lt. lmax-1) then
          ds(l1,m1,3)=ds(l1,m1,3)+(xi*OmegaP1(l1,m1)+MFlowP1(l1,m1))*(s(l1+1,m1)+ds(l1+1,m1,2)/2.)*exp(xi*Omega(l1+1,m1)/2.)
          dt(l1,m1,3)=dt(l1,m1,3)+(xi*OmegaP1(l1,m1)+MFlowP1(l1,m1))*(t(l1+1,m1)+dt(l1+1,m1,2)/2.)*exp(xi*Omega(l1+1,m1)/2.)
        end if
!c
!c  contribution from l component
!c
        ds(l1,m1,3)=ds(l1,m1,3)+MFlow0(l1,m1)*(s(l1,m1)+ds(l1,m1,2)/2.)*exp(xi*Omega(l1,m1)/2.)
        dt(l1,m1,3)=dt(l1,m1,3)+MFlow0(l1,m1)*(t(l1,m1)+dt(l1,m1,2)/2.)*exp(xi*Omega(l1,m1)/2.)
!c
!c  contribution from l-1 component
!c
        if ((l-1) .ge. m) then
          ds(l1,m1,3)=ds(l1,m1,3)+(xi*OmegaM1(l1,m1)+MFlowM1(l1,m1))*(s(l1-1,m1)+ds(l1-1,m1,2)/2.)*exp(xi*Omega(l1-1,m1)/2.)
          dt(l1,m1,3)=dt(l1,m1,3)+(xi*OmegaM1(l1,m1)+MFlowM1(l1,m1))*(t(l1-1,m1)+dt(l1-1,m1,2)/2.)*exp(xi*Omega(l1-1,m1)/2.)
        end if
!c
!c  contribution from l-2 component
!c
        if ((l-2) .ge. m) then
          ds(l1,m1,3)=ds(l1,m1,3)+(xi*OmegaM2(l1,m1)+MFlowM2(l1,m1))*(s(l1-2,m1)+ds(l1-2,m1,2)/2.)*exp(xi*Omega(l1-2,m1)/2.)
          dt(l1,m1,3)=dt(l1,m1,3)+(xi*OmegaM2(l1,m1)+MFlowM2(l1,m1))*(t(l1-2,m1)+dt(l1-2,m1,2)/2.)*exp(xi*Omega(l1-2,m1)/2.)
        end if
!c
!c  contribution from l-3 component
!c
        if ((l-3) .ge. m) then
          ds(l1,m1,3)=ds(l1,m1,3)+(xi*OmegaM3(l1,m1)+MFlowM3(l1,m1))*(s(l1-3,m1)+ds(l1-3,m1,2)/2.)*exp(xi*Omega(l1-3,m1)/2.)
          dt(l1,m1,3)=dt(l1,m1,3)+(xi*OmegaM3(l1,m1)+MFlowM3(l1,m1))*(t(l1-3,m1)+dt(l1-3,m1,2)/2.)*exp(xi*Omega(l1-3,m1)/2.)
        end if
!c
!c  contribution from l-4 component
!c
        if ((l-4) .ge. m) then
          ds(l1,m1,3)=ds(l1,m1,3)+(xi*OmegaM4(l1,m1)+MFlowM4(l1,m1))*(s(l1-4,m1)+ds(l1-4,m1,2)/2.)*exp(xi*Omega(l1-4,m1)/2.)
          dt(l1,m1,3)=dt(l1,m1,3)+(xi*OmegaM4(l1,m1)+MFlowM4(l1,m1))*(t(l1-4,m1)+dt(l1-4,m1,2)/2.)*exp(xi*Omega(l1-4,m1)/2.)
        end if
!c
!c  contribution from l-5 component
!c
        if ((l-5) .ge. m) then
          ds(l1,m1,3)=ds(l1,m1,3)+MFlowM5(l1,m1)*(s(l1-5,m1)+ds(l1-5,m1,2)/2.)*exp(xi*Omega(l1-5,m1)/2.)
          dt(l1,m1,3)=dt(l1,m1,3)+MFlowM5(l1,m1)*(t(l1-5,m1)+dt(l1-5,m1,2)/2.)*exp(xi*Omega(l1-5,m1)/2.)
        end if
!c
!c  contribution from l-6 component
!c
        if ((l-6) .ge. m) then
          ds(l1,m1,3)=ds(l1,m1,3)+MFlowM6(l1,m1)*(s(l1-6,m1)+ds(l1-6,m1,2)/2.)*exp(xi*Omega(l1-6,m1)/2.)
          dt(l1,m1,3)=dt(l1,m1,3)+MFlowM6(l1,m1)*(t(l1-6,m1)+dt(l1-6,m1,2)/2.)*exp(xi*Omega(l1-6,m1)/2.)
        end if
      end do ! End l-loop Step3
    end do ! End m-loop Step3

!c***********************************************************************
!c
!c Step 4
!c
!c***********************************************************************
    do m=1,lmax-1
      m1=m+1
      do l=m,lmax-1
        l1=l+1
        ds(l1,m1,4)=0.
        dt(l1,m1,4)=0.
!c
!c  contribution from l+6 component
!c
        if ((l+6) .lt. lmax-1) then
          ds(l1,m1,4)=ds(l1,m1,4)+MFlowP6(l1,m1)*(s(l1+6,m1)+ds(l1+6,m1,3))*exp(xi*Omega(l1+6,m1))
          dt(l1,m1,4)=dt(l1,m1,4)+MFlowP6(l1,m1)*(t(l1+6,m1)+dt(l1+6,m1,3))*exp(xi*Omega(l1+6,m1))
        end if
!c
!c  contribution from l+5 component
!c
        if ((l+5) .lt. lmax-1) then
          ds(l1,m1,4)=ds(l1,m1,4)+MFlowP5(l1,m1)*(s(l1+5,m1)+ds(l1+5,m1,3))*exp(xi*Omega(l1+5,m1))
          dt(l1,m1,4)=dt(l1,m1,4)+MFlowP5(l1,m1)*(t(l1+5,m1)+dt(l1+5,m1,3))*exp(xi*Omega(l1+5,m1))
        end if
!c
!c  contribution from l+4 component
!c
        if ((l+4) .lt. lmax-1) then
          ds(l1,m1,4)=ds(l1,m1,4)+(xi*OmegaP4(l1,m1)+MFlowP4(l1,m1))*(s(l1+4,m1)+ds(l1+4,m1,3))*exp(xi*Omega(l1+4,m1))
          dt(l1,m1,4)=dt(l1,m1,4)+(xi*OmegaP4(l1,m1)+MFlowP4(l1,m1))*(t(l1+4,m1)+dt(l1+4,m1,3))*exp(xi*Omega(l1+4,m1))
        end if
!c
!c  contribution from l+3 component
!c
        if ((l+3) .lt. lmax-1) then
          ds(l1,m1,4)=ds(l1,m1,4)+(xi*OmegaP3(l1,m1)+MFlowP3(l1,m1))*(s(l1+3,m1)+ds(l1+3,m1,3))*exp(xi*Omega(l1+3,m1))
          dt(l1,m1,4)=dt(l1,m1,4)+(xi*OmegaP3(l1,m1)+MFlowP3(l1,m1))*(t(l1+3,m1)+dt(l1+3,m1,3))*exp(xi*Omega(l1+3,m1))
        end if
!c
!c  contribution from l+2 component
!c
        if ((l+2) .lt. lmax-1) then
          ds(l1,m1,4)=ds(l1,m1,4)+(xi*OmegaP2(l1,m1)+MFlowP2(l1,m1))*(s(l1+2,m1)+ds(l1+2,m1,3))*exp(xi*Omega(l1+2,m1))
          dt(l1,m1,4)=dt(l1,m1,4)+(xi*OmegaP2(l1,m1)+MFlowP2(l1,m1))*(t(l1+2,m1)+dt(l1+2,m1,3))*exp(xi*Omega(l1+2,m1))
        end if
!c
!c  contribution from l+1 component
!c
        if ((l+1) .lt. lmax-1) then
          ds(l1,m1,4)=ds(l1,m1,4)+(xi*OmegaP1(l1,m1)+MFlowP1(l1,m1))*(s(l1+1,m1)+ds(l1+1,m1,3))*exp(xi*Omega(l1+1,m1))
          dt(l1,m1,4)=dt(l1,m1,4)+(xi*OmegaP1(l1,m1)+MFlowP1(l1,m1))*(t(l1+1,m1)+dt(l1+1,m1,3))*exp(xi*Omega(l1+1,m1))
        end if
!c
!c  contribution from l component
!c
        ds(l1,m1,4)=ds(l1,m1,4)+MFlow0(l1,m1)*(s(l1,m1)+ds(l1,m1,3))*exp(xi*Omega(l1,m1))
        dt(l1,m1,4)=dt(l1,m1,4)+MFlow0(l1,m1)*(t(l1,m1)+dt(l1,m1,3))*exp(xi*Omega(l1,m1))
!c
!c  contribution from l-1 component
!c
        if ((l-1) .ge. m) then
          ds(l1,m1,4)=ds(l1,m1,4)+(xi*OmegaM1(l1,m1)+MFlowM1(l1,m1))*(s(l1-1,m1)+ds(l1-1,m1,3))*exp(xi*Omega(l1-1,m1))
          dt(l1,m1,4)=dt(l1,m1,4)+(xi*OmegaM1(l1,m1)+MFlowM1(l1,m1))*(t(l1-1,m1)+dt(l1-1,m1,3))*exp(xi*Omega(l1-1,m1))
        end if
!c
!c  contribution from l-2 component
!c
        if ((l-2) .ge. m) then
          ds(l1,m1,4)=ds(l1,m1,4)+(xi*OmegaM2(l1,m1)+MFlowM2(l1,m1))*(s(l1-2,m1)+ds(l1-2,m1,3))*exp(xi*Omega(l1-2,m1))
          dt(l1,m1,4)=dt(l1,m1,4)+(xi*OmegaM2(l1,m1)+MFlowM2(l1,m1))*(t(l1-2,m1)+dt(l1-2,m1,3))*exp(xi*Omega(l1-2,m1))
        end if
!c
!c  contribution from l-3 component
!c
        if ((l-3) .ge. m) then
          ds(l1,m1,4)=ds(l1,m1,4)+(xi*OmegaM3(l1,m1)+MFlowM3(l1,m1))*(s(l1-3,m1)+ds(l1-3,m1,3))*exp(xi*Omega(l1-3,m1))
          dt(l1,m1,4)=dt(l1,m1,4)+(xi*OmegaM3(l1,m1)+MFlowM3(l1,m1))*(t(l1-3,m1)+dt(l1-3,m1,3))*exp(xi*Omega(l1-3,m1))
        end if
!c
!c  contribution from l-4 component
!c
        if ((l-4) .ge. m) then
          ds(l1,m1,4)=ds(l1,m1,4)+(xi*OmegaM4(l1,m1)+MFlowM4(l1,m1))*(s(l1-4,m1)+ds(l1-4,m1,3))*exp(xi*Omega(l1-4,m1))
          dt(l1,m1,4)=dt(l1,m1,4)+(xi*OmegaM4(l1,m1)+MFlowM4(l1,m1))*(t(l1-4,m1)+dt(l1-4,m1,3))*exp(xi*Omega(l1-4,m1))
        end if
!c
!c  contribution from l-5 component
!c
        if ((l-5) .ge. m) then
          ds(l1,m1,4)=ds(l1,m1,4)+MFlowM5(l1,m1)*(s(l1-5,m1)+ds(l1-5,m1,3))*exp(xi*Omega(l1-5,m1))
          dt(l1,m1,4)=dt(l1,m1,4)+MFlowM5(l1,m1)*(t(l1-5,m1)+dt(l1-5,m1,3))*exp(xi*Omega(l1-5,m1))
        end if
!c
!c  contribution from l-6 component
!c
        if ((l-6) .ge. m) then
          ds(l1,m1,4)=ds(l1,m1,4)+MFlowM6(l1,m1)*(s(l1-6,m1)+ds(l1-6,m1,3))*exp(xi*Omega(l1-6,m1))
          dt(l1,m1,4)=dt(l1,m1,4)+MFlowM6(l1,m1)*(t(l1-6,m1)+dt(l1-6,m1,3))*exp(xi*Omega(l1-6,m1))
        end if
      end do ! End l-loop Step4
    end do ! End m-loop Step4
!
! ****** Add changes and evolve
!
    do m=1,lmax-1
      m1=m+1
      do l=m,lmax-1
        l1=l+1
        el=real(l,r_typ)
        xlifetime=10.*3600.*(100./(el+0.1))**2.0
        !        xlifetime=xlifetime/2.
        call RANDOM_NUMBER (rnd)
        phase=Omega(l1,m1)+2.*(rnd-0.5)*sqrt(deltaT/xlifetime)
        arg=cos(phase)+xi*sin(phase)
        s(l1,m1)=(s(l1,m1) + ds(l1,m1,1)/6. + ds(l1,m1,2)/3. + ds(l1,m1,3)/3. + ds(l1,m1,4)/6.)*arg

        call RANDOM_NUMBER (rnd)
        phase=Omega(l1,m1)+2.*(rnd-0.5)*sqrt(deltaT/xlifetime)
        arg=cos(phase)+xi*sin(phase)
        t(l1,m1)=(t(l1,m1) + dt(l1,m1,1)/6. + dt(l1,m1,2)/3. + dt(l1,m1,3)/3. + dt(l1,m1,4)/6.)*arg
      end do
    end do
    !write(*,*) 'Spectral Coefficients updated'
!c***********************************************************************
!c                                                                      *
!c  Calculate the vector velocity components at each latitude.
!c  Equator is at jj=nx/2 + 2.5 due to wrap-around border
!c  j=1 is centered 0.5*dtheta above and below the equator
!c  j=nxhalf is centered 0.5*dtheta inside each pole
!c                                                                      *
!c***********************************************************************
    do j=1,nlhalf
      jn=nlhalf+j
      js=nlhalf+1-j
      theta=0.5*pi-(j-0.5)*dtheta
      x=cos(theta)
      sintheta=sin(theta)
      rst=1.0/sintheta
!c***********************************************************************
!c                                                                      *
!c  calculate the spectral coefficients for wavenumber m at all x.      *
!c                                                                      *
!c***********************************************************************
      m=0
      m1=m+1
      unorth(m1)=0.
      usouth(m1)=0.
      vnorth(m1)=0.
      vsouth(m1)=0.
!c***********************************************************************
!c                                                                      *
!c  Split non-axisymmetric signal into equal positive and               *
!c  negative freqencies.                                                *
!c                                                                      *
!c***********************************************************************
      do m=1,lmax-1
        m1=m+1
        m2=nphi+1-m
        call plm(m,x,nl,coef,p)
        sum1=(0.,0.)
        sum2=(0.,0.)
        sum3=(0.,0.)
        sum4=(0.,0.)
        do l=m,lmax-1
          l1=l+1
          l2=l+2
!c***********************************************************************
!c
!c  Construct velocity functions
!c
!c***********************************************************************
          eo=1.-2.*mod(l-m,2)
          v1=l*p(l2)/coef(l2,m1)-(l+1.)*p(l)/coef(l1,m1)
          v2=-m*p(l1)
          sum1 = sum1 + xi*s(l1,m1)*v2 - t(l1,m1)*v1
          sum2 = sum2 + eo*(xi*s(l1,m1)*v2 + t(l1,m1)*v1)
          sum3 = sum3 + s(l1,m1)*v1 + xi*t(l1,m1)*v2
          sum4 = sum4 - eo*(s(l1,m1)*v1 - xi*t(l1,m1)*v2)
        end do ! End l-loop
        unorth(m1)=0.5*sum1*rst
        usouth(m1)=0.5*sum2*rst
        unorth(m2)=0.5*conjg(sum1)*rst
        usouth(m2)=0.5*conjg(sum2)*rst
        vnorth(m1)=0.5*sum3*rst
        vsouth(m1)=0.5*sum4*rst
        vnorth(m2)=0.5*conjg(sum3)*rst
        vsouth(m2)=0.5*conjg(sum4)*rst
      end do ! End m-loop
      do m=lmax,nphi-lmax
        m1=m+1
        unorth(m1)=0.
        usouth(m1)=0.
        vnorth(m1)=0.
        vsouth(m1)=0.
      end do
!***********************************************************************
!                                                                      *
!  Calculate the vector velocity components at all phi positions.      *
!                                                                      *
!***********************************************************************
      call four1(unorth,nphi,-1)
      call four1(usouth,nphi,-1)
      call four1(vnorth,nphi,-1)
      call four1(vsouth,nphi,-1)
!
      do i=1,nphi
!
! ****** Get real part of complex values from FFT.
!
        u(i,jn)=real(unorth(i),r_typ)
        u(i,js)=real(usouth(i),r_typ)
        v(i,jn)=real(vnorth(i),r_typ)
        v(i,js)=real(vsouth(i),r_typ)
      end do ! End phi-loop
    end do ! End latitude-loop
!***********************************************************************
!                                                                      *
!  Write longitudinal and co-latitudinal velocity to disk file.        *
!                                                                      *
!***********************************************************************

!  [RMC] INSERTED PSI HDF5 OUTPUT HERE (SET 1D SCALES ABOVE TIME LOOP)
!   NOTE!  For now, I have disabled old binary output.
!          It should be put back in as an option?

    wt1 = wtime()
    write(velFileNum, '(i0.6)') ifile
!    open(unit=1,file=path(1:lenPath) // fname // ext,access='direct', &
!         status='unknown',recl=8*nphi)
!      do j=1,nl
!        write(1,rec=j) (u(i,j),i=1,nphi)
!      end do
!    close(1)
    wtime_io = wtime_io + (wtime() - wt1)
!
    call write_2d_file ('vp'//velFileNum//'.h5',nphi,nl,u,pvec,tvec,ierr)
!
    wt1 = wtime()
!    open(unit=1,file=path(1:lenPath) // fname // ext,access='direct', &
!         status='unknown',recl=8*nphi)
!      do j=1,nl
!        write(1,rec=j) (v(i,j),i=1,nphi)
!      enddo
!    close(1)
    wtime_io = wtime_io + (wtime() - wt1)
!
    call write_2d_file ('vt'//velFileNum//'.h5',nphi,nl,v,pvec,tvec,ierr)
!
    call ffopen (12,'flow_output_list.csv','a',ierr)
    write (12,'(F11.5,A1,A11,A1,A11)') curr_time/(3600*24),',', &
                              trim('vt'//velFileNum//'.h5'),',', &
                              trim('vp'//velFileNum//'.h5')
    close(12)
!
    write(*,*) 'Completed step ',ifile,' of ', nfiles
    write(*,*) '    max(vt)=',maxval(v),' min(vt)=',minval(v)
    write(*,*) '    max(vp)=',maxval(u),' min(vp)=',minval(u)
    flush(OUTPUT_UNIT)
!
! ***** Update time.
!
    curr_time = curr_time + deltaT
!
  enddo ! End time-loop
!
! ****** Get wall-clock time.
!
  wtime_total = wtime() - wtime_tmp

  call write_timing
!
end program conflow
!#######################################################################
subroutine write_2d_file (fname,ln1,ln2,f,s1,s2,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Write 2D data to H5 file FNAME.
!
!-----------------------------------------------------------------------
!
      use number_types
      use ds_def
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      real(r_typ), dimension(ln1,ln2) :: f
      real(r_typ), dimension(ln1) :: s1
      real(r_typ), dimension(ln2) :: s2
      integer :: ln1,ln2
      real(r_typ) :: t1,wtime
!
!-----------------------------------------------------------------------
!
      type(ds) :: s
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      t1 = wtime()
!
! ****** Set the structure components.
!
      s%ndim = 2
      s%dims(1) = ln1
      s%dims(2) = ln2
      s%dims(3) = 1
      s%scale = .true.
      s%hdf32 = .true.
!
      allocate (s%scales(1)%f(ln1))
      allocate (s%scales(2)%f(ln2))
      allocate (s%f(ln1,ln2,1))
!
      s%scales(1)%f(:) = s1(:)
      s%scales(2)%f(:) = s2(:)
      s%f(:,:,1) = f(:,:)
!
! ****** Write the data set.
!
      call wrh5 (fname,s,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRITE_2D_FILE:'
        write (*,*) '### Could not write the 2D data set.'
        write (*,*) 'File name: ',trim(fname)
        return
      end if
!
! ****** Free up memory.
!
      deallocate (s%scales(1)%f)
      deallocate (s%scales(2)%f)
      deallocate (s%f)
!
      wtime_io = wtime_io + (wtime() - t1)
!
end subroutine
!#######################################################################
subroutine write_timing
!
!-----------------------------------------------------------------------
!
! ****** Write out timing.
!
!-----------------------------------------------------------------------
!
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*), parameter :: FMT = '(a20,f20.2)'
!
!-----------------------------------------------------------------------
!
      write(*,*)
      write(*,"(a40)") repeat("-", 40)
      write(*,FMT) "Wall clock time:   ",wtime_total
      write(*,"(a40)") repeat("-", 40)
      write(*,FMT) "--> I/O:           ",wtime_io
      write(*,"(a40)") repeat("-", 40)
!
      write(*,*)
      flush(OUTPUT_UNIT)
!
end subroutine
!#######################################################################
function wtime ()
!
!*********************************************************************72
!
!  WTIME returns a reading of the wall clock time.
!
!  Discussion:
!    To get the elapsed wall clock time, call WTIME before and after
!    a given operation, and subtract the first reading from the second.
!  Licensing:
!    This code is distributed under the GNU LGPL license.
!  Author:
!    John Burkardt
!  Parameters:
!    Output, r_typ WTIME, the wall clock reading, in seconds.
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer*8 :: clock_max
      integer*8 :: clock_rate
      integer*8 :: clock_reading
!
      real(r_typ) :: wtime
!
!-----------------------------------------------------------------------
!
      call SYSTEM_CLOCK (clock_reading,clock_rate,clock_max)
!
      wtime = real(clock_reading,r_typ)/real(clock_rate,r_typ)
!
      return
end function
!#######################################################################
!
!-----------------------------------------------------------------------
!
! ****** Update log:
!
! 05/20/2022, RC, Version 0.1.2:
!   - Merged updates from Raphael's version.
!
! 06/07/2022, RC, Version 0.2.0:
!   - Added initial hdf5 output.
!
! 06/20/2022, RC, Version 0.3.0:
!   - Added timers, fixed implicit variables.
!   - Converted code to double precision.
!   - Temporarily disabled old binary output until it is made
!     as an option.
!
! 07/14/2022, RC, Version 0.4.0:
!   - Cleaned up some code
!   - Changes output file names to v[tp]######.h5
!   - Created output csv file listing output flow files and times.
!
!-----------------------------------------------------------------------
!
!#######################################################################
