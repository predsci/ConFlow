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
      character(*), parameter :: cvers='0.5.1'
      character(*), parameter :: cdate='07/19/2022'
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
module constants
!
!-----------------------------------------------------------------------
! ****** Constants in r_typ precision for use throughout the code.
! ****** Used for simplicity and readability.
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: zero = 0.0_r_typ
      real(r_typ), parameter :: one = 1.0_r_typ
      integer*8,   parameter :: one_int = 1
      real(r_typ), parameter :: two = 2.0_r_typ
      integer*8,   parameter :: two_int = 2
      real(r_typ), parameter :: three = 3._r_typ
      integer*8,   parameter :: three_int = 3
      real(r_typ), parameter :: four = 4._r_typ
      integer*8,   parameter :: four_int = 4
      real(r_typ), parameter :: five = 5._r_typ
      real(r_typ), parameter :: six = 6._r_typ
      real(r_typ), parameter :: seven = 7._r_typ
      real(r_typ), parameter :: nine = 9._r_typ
      real(r_typ), parameter :: ten = 10._r_typ
      real(r_typ), parameter :: fifteen = 15._r_typ
      real(r_typ), parameter :: sixteen = 16._r_typ
      real(r_typ), parameter :: half = 0.5_r_typ
      real(r_typ), parameter :: quarter = 0.25_r_typ
      real(r_typ), parameter :: twentyfour = 24.0_r_typ
      real(r_typ), parameter :: twentyfive = 25.0_r_typ
      real(r_typ), parameter :: fivehundred_i = 0.002_r_typ
      real(r_typ), parameter :: three_quarter = 0.75_r_typ
      real(r_typ), parameter :: two_third = 0.66666666666666666_r_typ
      real(r_typ), parameter :: third = 0.33333333333333333_r_typ
!
      real(r_typ), parameter :: pi = 3.1415926535897932_r_typ
      real(r_typ), parameter :: pi_two = 1.5707963267948966_r_typ
      real(r_typ), parameter :: pi_i = 0.3183098861837907_r_typ
      real(r_typ), parameter :: twopi = 6.2831853071795864_r_typ
      real(r_typ), parameter :: twopi_i = 0.15915494309189535_r_typ
      real(r_typ), parameter :: threepi_two = 4.71238898038469_r_typ
      real(r_typ), parameter :: threepi_four = 2.356194490192345_r_typ
!
      real(r_typ), parameter :: d2r = 0.017453292519943295_r_typ
      real(r_typ), parameter :: r2d = 57.29577951308232_r_typ
!
      real(r_typ), parameter :: rsun_cm = 6.96e10_r_typ
      real(r_typ), parameter :: rsun_cm2 = 4.84416e21_r_typ
!
      real(r_typ), parameter :: &
                         diff_km2_s_to_rs2_s = 2.0643413925221e-12_r_typ
      real(r_typ), parameter :: &
                         diff_km2_s_to_rs2_hr = 7.43162901307967e-09_r_typ
      real(r_typ), parameter :: &
                         km_s_to_rs_hr = 0.005172413793103448_r_typ
      real(r_typ), parameter :: &
                         m_s_to_rs_hr = 5.172413793103448e-06_r_typ
      real(r_typ), parameter :: &
                         km_s_to_rs_s = 1.4367816091954023e-06_r_typ
      real(r_typ), parameter :: output_flux_fac = 1.0e-21_r_typ
!
      real(r_typ), parameter :: small_value = tiny(one)
      real(r_typ), parameter :: large_value = huge(one)
      real(r_typ), parameter :: safety = 0.95_r_typ
!
end module
!#######################################################################
program conflow
!
!-----------------------------------------------------------------------
!
  use number_types
  use timing
  use constants
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
  integer :: i,j,l,m,l1,m1,l2,m2,jn,js,ierr,ii,jj
  integer :: lmax,lenPath,nlhalf,ifile,nfiles,itime
!
  real(r_typ) :: root3,root5,root7,root9
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
  real(r_typ), dimension(:), allocatable :: p_main,p_half,t_main,t_half
  real(r_typ), dimension(:), allocatable :: p_ext,t_ext
  real(r_typ), dimension(:,:), allocatable :: vt_psi,vp_psi,v_ext
  real(r_typ) ::  dphi_psi,dtheta_psi,pole

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
!-----------------------------------------------------------------------
!                                                                     
!  Numeric constants
!                                                                      
!-----------------------------------------------------------------------
  root3=sqrt(three)
  root5=sqrt(five)
  root7=sqrt(seven)
  root9=sqrt(nine)
  xi=(0.,one)
!-----------------------------------------------------------------------
!                                                                      
!  Physical constants
!                                                                      
!-----------------------------------------------------------------------
  rSun = 6.96e+08_r_typ   ! (meters)
  daySec = 86400.0_r_typ
!-----------------------------------------------------------------------
!                                                                      
!  Time (seconds) and space (radians) steps
!                                                                      
!-----------------------------------------------------------------------
  deltaT = fifteen*60.0_r_typ
  tmax = 28*24*3600
  curr_time = 0.
  nfiles = int(tmax/deltaT)
  dphi = twopi/nphi
  dtheta = pi/nl
!
!-----------------------------------------------------------------------
!
!****** Set up output scales and arrays.
!
! Currently, the VT and VP are on the same inner-half-half mesh
! with 1024x512 points.
!
! For HipFT input, we need VT to be 1024x513 on a PSI main-half mesh
! and                      VP to be 1025x512 on a PSI half-main mesh
!
! Until this can be done directly/correctly in this code,
! we add here interpolation to those grids before outputting the flows.
! This definetly should be replaced with modifying the above code to
! set these meshes directly to avoid interpolation.
!
! ****** Set up extended grid for current flows so that half-mesh
!        dimensions can be interpolated properly.
  allocate (v_ext(nphi+2,nl+2))
  allocate (p_ext(nphi+2))
  allocate (t_ext(nl+2))
!
! ****** Set extended ConFlow meshes.
!
  do i=1,nphi+2
    p_ext(i) = -half*dphi + (i-1)*dphi
  enddo
  do i=1,nl+2
    t_ext(i) = -half*dtheta + (i-1)*dtheta
  enddo
  print*,p_ext(1),t_ext(1),p_ext(nphi+2),t_ext(nl+2)
!
! ***** Allocate psi output arrays.
!
  allocate (vt_psi(nphi,nl+1))
  allocate (vp_psi(nphi+1,nl))
!
! ****** Set PSI meshes.
!
  allocate (p_main(nphi))
  allocate (p_half(nphi+1))
  allocate (t_main(nl))
  allocate (t_half(nl+1))
!
  dphi_psi = twopi/(nphi-1)
  dtheta_psi = pi/(nl-1)
!
  do i=1,nphi
    p_main(i) = (i-1)*dphi_psi
  enddo
  do i=1,nphi+1
    p_half(i) = -half*dphi_psi+(i-1)*dphi_psi
  enddo
  do i=1,nl
    t_main(i) = (i-1)*dtheta_psi
  enddo
  do i=1,nl+1
    t_half(i) = -half*dtheta_psi + (i-1)*dtheta_psi
  enddo
  print*,p_main(1),t_main(1),p_main(nphi),t_main(nl)
  print*,p_half(1),t_half(1),p_half(nphi+1),t_half(nl+1)
!
!-----------------------------------------------------------------------
!                                                                      
!  Differential Rotation coefficients m/s relative to Carrington
!                                                                      
!-----------------------------------------------------------------------
2 format(5(4x,f8.3))
3 format(6(4x,f8.3))
  open(unit=2,file=path(1:lenPath) // 'conflow_flow_parameters_aft_v1.txt', &
       status='old')
    read(2,2) t0,t1,t2,t3,t4
    write(*,*)
    write(*,*) 'Differential rotation coeffs (m/s) t0,t1,t2,t3,t4:'
    write(*,*) t0,t1,t2,t3,t4
    coefDR(1)=t0*deltaT/rSun
    coefDR(2)=t1*deltaT/rSun
    coefDR(3)=t2*deltaT/rSun
    coefDR(4)=t3*deltaT/rSun
    coefDR(5)=t4*deltaT/rSun
!-----------------------------------------------------------------------
!
! Meridional flow speed of supergranules in m/s poleward
!
!-----------------------------------------------------------------------
    read(2,3) s0,s1,s2,s3,s4,s5
    write(*,*)
    write(*,*) 'Meridional flow coeffs (m/s) s0,s1,s2,s3,s4,s5:'
    write(*,*) s0,s1,s2,s3,s4,s5
    coefMF(1)=s0*deltaT/rSun
    coefMF(2)=s1*deltaT/rSun
    coefMF(3)=s2*deltaT/rSun
    coefMF(4)=s3*deltaT/rSun
    coefMF(5)=s4*deltaT/rSun
    coefMF(6)=s5*deltaT/rSun
  close(2)
!-----------------------------------------------------------------------
!
!  Create coefficients for spherical harmonics.
!
!-----------------------------------------------------------------------
  call plmcoef(nl,coef)
!
  write(*,*)
  write(*,*) 'Legendre recurrance coefficients calculated.'
!-----------------------------------------------------------------------
!                                                                     
!  Construct the convection spectrum.                                 
!                                                                     
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
!
! Coupling Coefficients
!
!-----------------------------------------------------------------------
  do m=1,lmax
    m1=m+1
    em=real(m,r_typ)
    do l=m,lmax
      l1=l+1
      el=real(l,r_typ)
!
! Modify differential rotation and meridional flow with l (depth)
!
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
! -----------------------------------------------------------------------
! Raphael: It is the couplings below that take into account
! the advection equation satisfied by the velocity in presence of 
! the Meridional Flow (MF) and Differential Rotation (DR). 
!
! Here, "Coupling" means that the spectral coefficient for the velocity  
! components at a given l and m value (S_lm and T_lm in Hathaway et al.) 
! depend on spectral coefficient at different l values, up to l+/-6. 
!
! These couplings come from the projection and orthoganlity of the 
! associated Legendre Polynomials "P_lm". E.g. their recursion relation 
! eliminate the powers of cos(theta) of the formula of the DR, which is 
! why we do not see them anywhere in the algorithm.  
! Without advection by MF and DR, these coupling terms would not exist. 
!
! See Appendix of Hathaway et al. (2010). However, only the details for 
! the radial component (R_lm) of the velocity are explicited. 
! So far, there is no reference showing even a partial derivation of the 
! couplings for the toroidal and poloidal components (S_lm and T_lm). 
!  
!
! TODO: partial or full derivation of the coupling terms for 
! the toroidal and poloidal components (S_lm and T_lm)
! that is not documented in Dave's papers. 
!
! Outcome: Per Dave's explanations, S_lm and T_lm are not following the 
! the analytical development set by equation (2) and (3) for v_theta
! and v_phi. Instead, they are assumed to follow the same development
! as equation (1) for R_lm to keep things simple. 
! This choice was made as (2) and (3) may not be analytically developed.
! This sets a possible discrepancy in both
! the math and the physics of the evolution of horizontal flow field
! with respect to what equation (2) and (3) otherwise imposes. 
! The unrealistic network may be a consequence of this simplification. 
! 
! -----------------------------------------------------------------------
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
!
! Coupling with l component
!
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
!
!  Coupling with l-1 component
!
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
!
!  Coupling with l-2 component
!
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
!
!  Coupling with l-3 component
!
      if ((l-3) .ge. m) then
        BAA0=1./(coef(l1,m1)*coef(l1-1,m1)*coef(l1-2,m1))
        OmegaM3(l1,m1)=em*BAA0*DR3

        BAA2=BAA0*(BAALP1 + BAALP0 + BAALM1 + BAALM2)
        BAA3=BAA0*BAALM3
        MFlowM3(l1,m1)=(el-3.)*BAA0*MF2 + ((el-3.)*BAA2-(el-2.)*BAA3)*MF4
      end if
!
!  Coupling with l-4 component
!
      if ((l-4) .ge. m) then
        BAA0=1./(coef(l1,m1)*coef(l1-1,m1)*coef(l1-2,m1)*coef(l1-3,m1))
        OmegaM4(l1,m1)=em*BAA0*DR4

        BAA1=BAA0*(BAALP1 + BAALP0 + BAALM1 + BAALM2)
        BAA2=BAA0*BAALM4
        MFlowM4(l1,m1)=(el-4.)*BAA0*MF3 + ((el-4.)*BAA1-(el-3.)*BAA2)*MF5
      end if
!
!  Coupling with l-5 component
!
      if ((l-5) .ge. m) then
        BAA0=1./(coef(l1,m1)*coef(l1-1,m1)*coef(l1-2,m1)*coef(l1-3,m1)*coef(l1-4,m1))
        MFlowM5(l1,m1)=(el-5.)*BAA0*MF4
      end if
!
!  Coupling with l-6 component
!
      if ((l-6) .ge. m) then
        BAA0=1./(coef(l1,m1)*coef(l1-1,m1)*coef(l1-2,m1)*coef(l1-3,m1)*coef(l1-4,m1)*coef(l1-5,m1))
        MFlowM6(l1,m1)=(el-6.)*BAA0*MF5
      end if
    end do ! End of m-loop
  end do ! End of l=loop
!-----------------------------------------------------------------------
!
!  Construct series of velocity maps at deltaT (in seconds) intervals
!
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
!
!  Evolve spectral coefficients using 4th order Runga-Kutta
!
!  Raphael: In what follows, s and t are the same as S_lm and T_lm 
!  in Hathaway 2010. 
!
! Step 1
!
!-----------------------------------------------------------------------
    do m=1,lmax-1
      m1=m+1
      do l=m,lmax-1
        l1=l+1
        ds(l1,m1,1)=0.
        dt(l1,m1,1)=0.
!
!  contribution from l+6 component
!
        if ((l+6) .lt. lmax-1) then
          ds(l1,m1,1)=ds(l1,m1,1)+MFlowP6(l1,m1)*s(l1+6,m1)
          dt(l1,m1,1)=dt(l1,m1,1)+MFlowP6(l1,m1)*t(l1+6,m1)
        end if
!
!  contribution from l+5 component
!
        if ((l+5) .lt. lmax-1) then
          ds(l1,m1,1)=ds(l1,m1,1)+MFlowP5(l1,m1)*s(l1+5,m1)
          dt(l1,m1,1)=dt(l1,m1,1)+MFlowP5(l1,m1)*t(l1+5,m1)
        end if
!
!  contribution from l+4 component
!
        if ((l+4) .lt. lmax-1) then
          ds(l1,m1,1)=ds(l1,m1,1)+(xi*OmegaP4(l1,m1)+MFlowP4(l1,m1))*s(l1+4,m1)
          dt(l1,m1,1)=dt(l1,m1,1)+(xi*OmegaP4(l1,m1)+MFlowP4(l1,m1))*t(l1+4,m1)
         end if
!
!  contribution from l+3 component
!
        if ((l+3) .lt. lmax-1) then
          ds(l1,m1,1)=ds(l1,m1,1)+(xi*OmegaP3(l1,m1)+MFlowP3(l1,m1))*s(l1+3,m1)
          dt(l1,m1,1)=dt(l1,m1,1)+(xi*OmegaP3(l1,m1)+MFlowP3(l1,m1))*t(l1+3,m1)
        end if
!
!  contribution from l+2 component
!
        if ((l+2) .lt. lmax-1) then
         ds(l1,m1,1)=ds(l1,m1,1)+(xi*OmegaP2(l1,m1)+MFlowP2(l1,m1))*s(l1+2,m1)
         dt(l1,m1,1)=dt(l1,m1,1)+(xi*OmegaP2(l1,m1)+MFlowP2(l1,m1))*t(l1+2,m1)
        end if
!
!  contribution from l+1 component
!
        if ((l+1) .lt. lmax-1) then
          ds(l1,m1,1)=ds(l1,m1,1)+(xi*OmegaP1(l1,m1)+MFlowP1(l1,m1))*s(l1+1,m1)
          dt(l1,m1,1)=dt(l1,m1,1)+(xi*OmegaP1(l1,m1)+MFlowP1(l1,m1))*t(l1+1,m1)
        end if
!
!  contribution from l component
!
        ds(l1,m1,1)=ds(l1,m1,1)+(0.0*xi*Omega(l1,m1)+MFlow0(l1,m1))*s(l1,m1) ! No Omega(l1,m1) term?
        dt(l1,m1,1)=dt(l1,m1,1)+(0.0*xi*Omega(l1,m1)+MFlow0(l1,m1))*t(l1,m1) ! No Omega(l1,m1) term?
!
!  contribution from l-1 component
!
        if ((l-1) .ge. m) then
          ds(l1,m1,1)=ds(l1,m1,1)+(xi*OmegaM1(l1,m1)+MFlowM1(l1,m1))*s(l1-1,m1)
          dt(l1,m1,1)=dt(l1,m1,1)+(xi*OmegaM1(l1,m1)+MFlowM1(l1,m1))*t(l1-1,m1)
        end if
!
!  contribution from l-2 component
!
        if ((l-2) .ge. m) then
          ds(l1,m1,1)=ds(l1,m1,1)+(xi*OmegaM2(l1,m1)+MFlowM2(l1,m1))*s(l1-2,m1)
          dt(l1,m1,1)=dt(l1,m1,1)+(xi*OmegaM2(l1,m1)+MFlowM2(l1,m1))*t(l1-2,m1)
        end if
!
!  contribution from l-3 component
!
        if ((l-3) .ge. m) then
          ds(l1,m1,1)=ds(l1,m1,1)+(xi*OmegaM3(l1,m1)+MFlowM3(l1,m1))*s(l1-3,m1)
          dt(l1,m1,1)=dt(l1,m1,1)+(xi*OmegaM3(l1,m1)+MFlowM3(l1,m1))*t(l1-3,m1)
        end if
!
!  contribution from l-4 component
!
        if ((l-4) .ge. m) then
          ds(l1,m1,1)=ds(l1,m1,1)+(xi*OmegaM4(l1,m1)+MFlowM4(l1,m1))*s(l1-4,m1)
          dt(l1,m1,1)=dt(l1,m1,1)+(xi*OmegaM4(l1,m1)+MFlowM4(l1,m1))*t(l1-4,m1)
        end if
!
!  contribution from l-5 component
!
        if ((l-5) .ge. m) then
          ds(l1,m1,1)=ds(l1,m1,1)+MFlowM5(l1,m1)*s(l1-5,m1)
          dt(l1,m1,1)=dt(l1,m1,1)+MFlowM5(l1,m1)*t(l1-5,m1)
        end if
!
!  contribution from l-6 component
!
        if ((l-6) .ge. m) then
          ds(l1,m1,1)=ds(l1,m1,1)+MFlowM6(l1,m1)*s(l1-6,m1)
          dt(l1,m1,1)=dt(l1,m1,1)+MFlowM6(l1,m1)*t(l1-6,m1)
        end if
      end do ! End l-loop step1
    end do ! End m-loop step1
!-----------------------------------------------------------------------
!
! Step 2
!
!-----------------------------------------------------------------------
    do m=1,lmax-1
      m1=m+1
      do l=m,lmax-1
        l1=l+1
        ds(l1,m1,2)=0.
        dt(l1,m1,2)=0.
!
!  contribution from l+6 component
!
        if ((l+6) .lt. lmax-1) then
          ds(l1,m1,2)=ds(l1,m1,2)+MFlowP6(l1,m1)*(s(l1+6,m1)+ds(l1+6,m1,1)/2.)*exp(xi*Omega(l1+6,m1)/2.)
          dt(l1,m1,2)=dt(l1,m1,2)+MFlowP6(l1,m1)*(t(l1+6,m1)+dt(l1+6,m1,1)/2.)*exp(xi*Omega(l1+6,m1)/2.)
        end if
!
!  contribution from l+5 component
!
        if ((l+5) .lt. lmax-1) then
          ds(l1,m1,2)=ds(l1,m1,2)+MFlowP5(l1,m1)*(s(l1+5,m1)+ds(l1+5,m1,1)/2.)*exp(xi*Omega(l1+5,m1)/2.)
          dt(l1,m1,2)=dt(l1,m1,2)+MFlowP5(l1,m1)*(t(l1+5,m1)+dt(l1+5,m1,1)/2.)*exp(xi*Omega(l1+5,m1)/2.)
        end if
!
!  contribution from l+4 component
!
        if ((l+4) .lt. lmax-1) then
          ds(l1,m1,2)=ds(l1,m1,2)+(xi*OmegaP4(l1,m1)+MFlowP4(l1,m1))*(s(l1+4,m1)+ds(l1+4,m1,1)/2.)*exp(xi*Omega(l1+4,m1)/2.)
          dt(l1,m1,2)=dt(l1,m1,2)+(xi*OmegaP4(l1,m1)+MFlowP4(l1,m1))*(t(l1+4,m1)+dt(l1+4,m1,1)/2.)*exp(xi*Omega(l1+4,m1)/2.)
        end if
!
!  contribution from l+3 component
!
        if ((l+3) .lt. lmax-1) then
          ds(l1,m1,2)=ds(l1,m1,2)+(xi*OmegaP3(l1,m1)+MFlowP3(l1,m1))*(s(l1+3,m1)+ds(l1+3,m1,1)/2.)*exp(xi*Omega(l1+3,m1)/2.)
          dt(l1,m1,2)=dt(l1,m1,2)+(xi*OmegaP3(l1,m1)+MFlowP3(l1,m1))*(t(l1+3,m1)+dt(l1+3,m1,1)/2.)*exp(xi*Omega(l1+3,m1)/2.)
        end if
!
!  contribution from l+2 component
!
        if ((l+2) .lt. lmax-1) then
          ds(l1,m1,2)=ds(l1,m1,2)+(xi*OmegaP2(l1,m1)+MFlowP2(l1,m1))*(s(l1+2,m1)+ds(l1+2,m1,1)/2.)*exp(xi*Omega(l1+2,m1)/2.)
          dt(l1,m1,2)=dt(l1,m1,2)+(xi*OmegaP2(l1,m1)+MFlowP2(l1,m1))*(t(l1+2,m1)+dt(l1+2,m1,1)/2.)*exp(xi*Omega(l1+2,m1)/2.)
        end if
!
!  contribution from l+1 component
!
        if ((l+1) .lt. lmax-1) then
          ds(l1,m1,2)=ds(l1,m1,2)+(xi*OmegaP1(l1,m1)+MFlowP1(l1,m1))*(s(l1+1,m1)+ds(l1+1,m1,1)/2.)*exp(xi*Omega(l1+1,m1)/2.)
          dt(l1,m1,2)=dt(l1,m1,2)+(xi*OmegaP1(l1,m1)+MFlowP1(l1,m1))*(t(l1+1,m1)+dt(l1+1,m1,1)/2.)*exp(xi*Omega(l1+1,m1)/2.)
        end if
!
!  contribution from l component
!
        ds(l1,m1,2)=ds(l1,m1,2)+MFlow0(l1,m1)*(s(l1,m1)+ds(l1,m1,1)/2.)*exp(xi*Omega(l1,m1)/2.)
        dt(l1,m1,2)=dt(l1,m1,2)+MFlow0(l1,m1)*(t(l1,m1)+dt(l1,m1,1)/2.)*exp(xi*Omega(l1,m1)/2.)
!
!  contribution from l-1 component
!
        if ((l-1) .ge. m) then
          ds(l1,m1,2)=ds(l1,m1,2)+(xi*OmegaM1(l1,m1)+MFlowM1(l1,m1))*(s(l1-1,m1)+ds(l1-1,m1,1)/2.)*exp(xi*Omega(l1-1,m1)/2.)
          dt(l1,m1,2)=dt(l1,m1,2)+(xi*OmegaM1(l1,m1)+MFlowM1(l1,m1))*(t(l1-1,m1)+dt(l1-1,m1,1)/2.)*exp(xi*Omega(l1-1,m1)/2.)
        end if
!
!  contribution from l-2 component
!
        if ((l-2) .ge. m) then
          ds(l1,m1,2)=ds(l1,m1,2)+(xi*OmegaM2(l1,m1)+MFlowM2(l1,m1))*(s(l1-2,m1)+ds(l1-2,m1,1)/2.)*exp(xi*Omega(l1-2,m1)/2.)
          dt(l1,m1,2)=dt(l1,m1,2)+(xi*OmegaM2(l1,m1)+MFlowM2(l1,m1))*(t(l1-2,m1)+dt(l1-2,m1,1)/2.)*exp(xi*Omega(l1-2,m1)/2.)
        end if
!
!  contribution from l-3 component
!
        if ((l-3) .ge. m) then
          ds(l1,m1,2)=ds(l1,m1,2)+(xi*OmegaM3(l1,m1)+MFlowM3(l1,m1))*(s(l1-3,m1)+ds(l1-3,m1,1)/2.)*exp(xi*Omega(l1-3,m1)/2.)
          dt(l1,m1,2)=dt(l1,m1,2)+(xi*OmegaM3(l1,m1)+MFlowM3(l1,m1))*(t(l1-3,m1)+dt(l1-3,m1,1)/2.)*exp(xi*Omega(l1-3,m1)/2.)
        end if
!
!  contribution from l-4 component
!
        if ((l-4) .ge. m) then
          ds(l1,m1,2)=ds(l1,m1,2)+(xi*OmegaM4(l1,m1)+MFlowM4(l1,m1))*(s(l1-4,m1)+ds(l1-4,m1,1)/2.)*exp(xi*Omega(l1-4,m1)/2.)
          dt(l1,m1,2)=dt(l1,m1,2)+(xi*OmegaM4(l1,m1)+MFlowM4(l1,m1))*(t(l1-4,m1)+dt(l1-4,m1,1)/2.)*exp(xi*Omega(l1-4,m1)/2.)
        end if
!
!  contribution from l-5 component
!
        if ((l-5) .ge. m) then
          ds(l1,m1,2)=ds(l1,m1,2)+MFlowM5(l1,m1)*(s(l1-5,m1)+ds(l1-5,m1,1)/2.)*exp(xi*Omega(l1-5,m1)/2.)
          dt(l1,m1,2)=dt(l1,m1,2)+MFlowM5(l1,m1)*(t(l1-5,m1)+dt(l1-5,m1,1)/2.)*exp(xi*Omega(l1-5,m1)/2.)
        end if
!
!  contribution from l-6 component
!
        if ((l-6) .ge. m) then
          ds(l1,m1,2)=ds(l1,m1,2)+MFlowM6(l1,m1)*(s(l1-6,m1)+ds(l1-6,m1,1)/2.)*exp(xi*Omega(l1-6,m1)/2.)
          dt(l1,m1,2)=dt(l1,m1,2)+MFlowM6(l1,m1)*(t(l1-6,m1)+dt(l1-6,m1,1)/2.)*exp(xi*Omega(l1-6,m1)/2.)
        end if
      end do ! End l-loop Step2
    end do ! End m-loop Step2
!-----------------------------------------------------------------------
!
! Step 3
!
!-----------------------------------------------------------------------
    do m=1,lmax-1
      m1=m+1
      do l=m,lmax-1
        l1=l+1
        ds(l1,m1,3)=0.
        dt(l1,m1,3)=0.
!
!  contribution from l+6 component
!
        if ((l+6) .lt. lmax-1) then
          ds(l1,m1,3)=ds(l1,m1,3)+MFlowP6(l1,m1)*(s(l1+6,m1)+ds(l1+6,m1,2)/2.)*exp(xi*Omega(l1+6,m1)/2.)
          dt(l1,m1,3)=dt(l1,m1,3)+MFlowP6(l1,m1)*(t(l1+6,m1)+dt(l1+6,m1,2)/2.)*exp(xi*Omega(l1+6,m1)/2.)
        end if
!
!  contribution from l+5 component
!
        if ((l+5) .lt. lmax-1) then
         ds(l1,m1,3)=ds(l1,m1,3)+MFlowP5(l1,m1)*(s(l1+5,m1)+ds(l1+5,m1,2)/2.)*exp(xi*Omega(l1+5,m1)/2.)
         dt(l1,m1,3)=dt(l1,m1,3)+MFlowP5(l1,m1)*(t(l1+5,m1)+dt(l1+5,m1,2)/2.)*exp(xi*Omega(l1+5,m1)/2.)
        end if
!
!  contribution from l+4 component
!
        if ((l+4) .lt. lmax-1) then
          ds(l1,m1,3)=ds(l1,m1,3)+(xi*OmegaP4(l1,m1)+MFlowP4(l1,m1))*(s(l1+4,m1)+ds(l1+4,m1,2)/2.)*exp(xi*Omega(l1+4,m1)/2.)
          dt(l1,m1,3)=dt(l1,m1,3)+(xi*OmegaP4(l1,m1)+MFlowP4(l1,m1))*(t(l1+4,m1)+dt(l1+4,m1,2)/2.)*exp(xi*Omega(l1+4,m1)/2.)
        end if
!
!  contribution from l+3 component
!
        if ((l+3) .lt. lmax-1) then
          ds(l1,m1,3)=ds(l1,m1,3)+(xi*OmegaP3(l1,m1)+MFlowP3(l1,m1))*(s(l1+3,m1)+ds(l1+3,m1,2)/2.)*exp(xi*Omega(l1+3,m1)/2.)
          dt(l1,m1,3)=dt(l1,m1,3)+(xi*OmegaP3(l1,m1)+MFlowP3(l1,m1))*(t(l1+3,m1)+dt(l1+3,m1,2)/2.)*exp(xi*Omega(l1+3,m1)/2.)
        end if
!
!  contribution from l+2 component
!
        if ((l+2) .lt. lmax-1) then
          ds(l1,m1,3)=ds(l1,m1,3)+(xi*OmegaP2(l1,m1)+MFlowP2(l1,m1))*(s(l1+2,m1)+ds(l1+2,m1,2)/2.)*exp(xi*Omega(l1+2,m1)/2.)
          dt(l1,m1,3)=dt(l1,m1,3)+(xi*OmegaP2(l1,m1)+MFlowP2(l1,m1))*(t(l1+2,m1)+dt(l1+2,m1,2)/2.)*exp(xi*Omega(l1+2,m1)/2.)
        end if
!
!  contribution from l+1 component
!
        if ((l+1) .lt. lmax-1) then
          ds(l1,m1,3)=ds(l1,m1,3)+(xi*OmegaP1(l1,m1)+MFlowP1(l1,m1))*(s(l1+1,m1)+ds(l1+1,m1,2)/2.)*exp(xi*Omega(l1+1,m1)/2.)
          dt(l1,m1,3)=dt(l1,m1,3)+(xi*OmegaP1(l1,m1)+MFlowP1(l1,m1))*(t(l1+1,m1)+dt(l1+1,m1,2)/2.)*exp(xi*Omega(l1+1,m1)/2.)
        end if
!
!  contribution from l component
!
        ds(l1,m1,3)=ds(l1,m1,3)+MFlow0(l1,m1)*(s(l1,m1)+ds(l1,m1,2)/2.)*exp(xi*Omega(l1,m1)/2.)
        dt(l1,m1,3)=dt(l1,m1,3)+MFlow0(l1,m1)*(t(l1,m1)+dt(l1,m1,2)/2.)*exp(xi*Omega(l1,m1)/2.)
!
!  contribution from l-1 component
!
        if ((l-1) .ge. m) then
          ds(l1,m1,3)=ds(l1,m1,3)+(xi*OmegaM1(l1,m1)+MFlowM1(l1,m1))*(s(l1-1,m1)+ds(l1-1,m1,2)/2.)*exp(xi*Omega(l1-1,m1)/2.)
          dt(l1,m1,3)=dt(l1,m1,3)+(xi*OmegaM1(l1,m1)+MFlowM1(l1,m1))*(t(l1-1,m1)+dt(l1-1,m1,2)/2.)*exp(xi*Omega(l1-1,m1)/2.)
        end if
!
!  contribution from l-2 component
!
        if ((l-2) .ge. m) then
          ds(l1,m1,3)=ds(l1,m1,3)+(xi*OmegaM2(l1,m1)+MFlowM2(l1,m1))*(s(l1-2,m1)+ds(l1-2,m1,2)/2.)*exp(xi*Omega(l1-2,m1)/2.)
          dt(l1,m1,3)=dt(l1,m1,3)+(xi*OmegaM2(l1,m1)+MFlowM2(l1,m1))*(t(l1-2,m1)+dt(l1-2,m1,2)/2.)*exp(xi*Omega(l1-2,m1)/2.)
        end if
!
!  contribution from l-3 component
!
        if ((l-3) .ge. m) then
          ds(l1,m1,3)=ds(l1,m1,3)+(xi*OmegaM3(l1,m1)+MFlowM3(l1,m1))*(s(l1-3,m1)+ds(l1-3,m1,2)/2.)*exp(xi*Omega(l1-3,m1)/2.)
          dt(l1,m1,3)=dt(l1,m1,3)+(xi*OmegaM3(l1,m1)+MFlowM3(l1,m1))*(t(l1-3,m1)+dt(l1-3,m1,2)/2.)*exp(xi*Omega(l1-3,m1)/2.)
        end if
!
!  contribution from l-4 component
!
        if ((l-4) .ge. m) then
          ds(l1,m1,3)=ds(l1,m1,3)+(xi*OmegaM4(l1,m1)+MFlowM4(l1,m1))*(s(l1-4,m1)+ds(l1-4,m1,2)/2.)*exp(xi*Omega(l1-4,m1)/2.)
          dt(l1,m1,3)=dt(l1,m1,3)+(xi*OmegaM4(l1,m1)+MFlowM4(l1,m1))*(t(l1-4,m1)+dt(l1-4,m1,2)/2.)*exp(xi*Omega(l1-4,m1)/2.)
        end if
!
!  contribution from l-5 component
!
        if ((l-5) .ge. m) then
          ds(l1,m1,3)=ds(l1,m1,3)+MFlowM5(l1,m1)*(s(l1-5,m1)+ds(l1-5,m1,2)/2.)*exp(xi*Omega(l1-5,m1)/2.)
          dt(l1,m1,3)=dt(l1,m1,3)+MFlowM5(l1,m1)*(t(l1-5,m1)+dt(l1-5,m1,2)/2.)*exp(xi*Omega(l1-5,m1)/2.)
        end if
!
!  contribution from l-6 component
!
        if ((l-6) .ge. m) then
          ds(l1,m1,3)=ds(l1,m1,3)+MFlowM6(l1,m1)*(s(l1-6,m1)+ds(l1-6,m1,2)/2.)*exp(xi*Omega(l1-6,m1)/2.)
          dt(l1,m1,3)=dt(l1,m1,3)+MFlowM6(l1,m1)*(t(l1-6,m1)+dt(l1-6,m1,2)/2.)*exp(xi*Omega(l1-6,m1)/2.)
        end if
      end do ! End l-loop Step3
    end do ! End m-loop Step3
!-----------------------------------------------------------------------
!
! Step 4
!
!-----------------------------------------------------------------------
    do m=1,lmax-1
      m1=m+1
      do l=m,lmax-1
        l1=l+1
        ds(l1,m1,4)=0.
        dt(l1,m1,4)=0.
!
!  contribution from l+6 component
!
        if ((l+6) .lt. lmax-1) then
          ds(l1,m1,4)=ds(l1,m1,4)+MFlowP6(l1,m1)*(s(l1+6,m1)+ds(l1+6,m1,3))*exp(xi*Omega(l1+6,m1))
          dt(l1,m1,4)=dt(l1,m1,4)+MFlowP6(l1,m1)*(t(l1+6,m1)+dt(l1+6,m1,3))*exp(xi*Omega(l1+6,m1))
        end if
!
!  contribution from l+5 component
!
        if ((l+5) .lt. lmax-1) then
          ds(l1,m1,4)=ds(l1,m1,4)+MFlowP5(l1,m1)*(s(l1+5,m1)+ds(l1+5,m1,3))*exp(xi*Omega(l1+5,m1))
          dt(l1,m1,4)=dt(l1,m1,4)+MFlowP5(l1,m1)*(t(l1+5,m1)+dt(l1+5,m1,3))*exp(xi*Omega(l1+5,m1))
        end if
!
!  contribution from l+4 component
!
        if ((l+4) .lt. lmax-1) then
          ds(l1,m1,4)=ds(l1,m1,4)+(xi*OmegaP4(l1,m1)+MFlowP4(l1,m1))*(s(l1+4,m1)+ds(l1+4,m1,3))*exp(xi*Omega(l1+4,m1))
          dt(l1,m1,4)=dt(l1,m1,4)+(xi*OmegaP4(l1,m1)+MFlowP4(l1,m1))*(t(l1+4,m1)+dt(l1+4,m1,3))*exp(xi*Omega(l1+4,m1))
        end if
!
!  contribution from l+3 component
!
        if ((l+3) .lt. lmax-1) then
          ds(l1,m1,4)=ds(l1,m1,4)+(xi*OmegaP3(l1,m1)+MFlowP3(l1,m1))*(s(l1+3,m1)+ds(l1+3,m1,3))*exp(xi*Omega(l1+3,m1))
          dt(l1,m1,4)=dt(l1,m1,4)+(xi*OmegaP3(l1,m1)+MFlowP3(l1,m1))*(t(l1+3,m1)+dt(l1+3,m1,3))*exp(xi*Omega(l1+3,m1))
        end if
!
!  contribution from l+2 component
!
        if ((l+2) .lt. lmax-1) then
          ds(l1,m1,4)=ds(l1,m1,4)+(xi*OmegaP2(l1,m1)+MFlowP2(l1,m1))*(s(l1+2,m1)+ds(l1+2,m1,3))*exp(xi*Omega(l1+2,m1))
          dt(l1,m1,4)=dt(l1,m1,4)+(xi*OmegaP2(l1,m1)+MFlowP2(l1,m1))*(t(l1+2,m1)+dt(l1+2,m1,3))*exp(xi*Omega(l1+2,m1))
        end if
!
!  contribution from l+1 component
!
        if ((l+1) .lt. lmax-1) then
          ds(l1,m1,4)=ds(l1,m1,4)+(xi*OmegaP1(l1,m1)+MFlowP1(l1,m1))*(s(l1+1,m1)+ds(l1+1,m1,3))*exp(xi*Omega(l1+1,m1))
          dt(l1,m1,4)=dt(l1,m1,4)+(xi*OmegaP1(l1,m1)+MFlowP1(l1,m1))*(t(l1+1,m1)+dt(l1+1,m1,3))*exp(xi*Omega(l1+1,m1))
        end if
!
!  contribution from l component
!
        ds(l1,m1,4)=ds(l1,m1,4)+MFlow0(l1,m1)*(s(l1,m1)+ds(l1,m1,3))*exp(xi*Omega(l1,m1))
        dt(l1,m1,4)=dt(l1,m1,4)+MFlow0(l1,m1)*(t(l1,m1)+dt(l1,m1,3))*exp(xi*Omega(l1,m1))
!
!  contribution from l-1 component
!
        if ((l-1) .ge. m) then
          ds(l1,m1,4)=ds(l1,m1,4)+(xi*OmegaM1(l1,m1)+MFlowM1(l1,m1))*(s(l1-1,m1)+ds(l1-1,m1,3))*exp(xi*Omega(l1-1,m1))
          dt(l1,m1,4)=dt(l1,m1,4)+(xi*OmegaM1(l1,m1)+MFlowM1(l1,m1))*(t(l1-1,m1)+dt(l1-1,m1,3))*exp(xi*Omega(l1-1,m1))
        end if
!
!  contribution from l-2 component
!
        if ((l-2) .ge. m) then
          ds(l1,m1,4)=ds(l1,m1,4)+(xi*OmegaM2(l1,m1)+MFlowM2(l1,m1))*(s(l1-2,m1)+ds(l1-2,m1,3))*exp(xi*Omega(l1-2,m1))
          dt(l1,m1,4)=dt(l1,m1,4)+(xi*OmegaM2(l1,m1)+MFlowM2(l1,m1))*(t(l1-2,m1)+dt(l1-2,m1,3))*exp(xi*Omega(l1-2,m1))
        end if
!
!  contribution from l-3 component
!
        if ((l-3) .ge. m) then
          ds(l1,m1,4)=ds(l1,m1,4)+(xi*OmegaM3(l1,m1)+MFlowM3(l1,m1))*(s(l1-3,m1)+ds(l1-3,m1,3))*exp(xi*Omega(l1-3,m1))
          dt(l1,m1,4)=dt(l1,m1,4)+(xi*OmegaM3(l1,m1)+MFlowM3(l1,m1))*(t(l1-3,m1)+dt(l1-3,m1,3))*exp(xi*Omega(l1-3,m1))
        end if
!
!  contribution from l-4 component
!
        if ((l-4) .ge. m) then
          ds(l1,m1,4)=ds(l1,m1,4)+(xi*OmegaM4(l1,m1)+MFlowM4(l1,m1))*(s(l1-4,m1)+ds(l1-4,m1,3))*exp(xi*Omega(l1-4,m1))
          dt(l1,m1,4)=dt(l1,m1,4)+(xi*OmegaM4(l1,m1)+MFlowM4(l1,m1))*(t(l1-4,m1)+dt(l1-4,m1,3))*exp(xi*Omega(l1-4,m1))
        end if
!
!  contribution from l-5 component
!
        if ((l-5) .ge. m) then
          ds(l1,m1,4)=ds(l1,m1,4)+MFlowM5(l1,m1)*(s(l1-5,m1)+ds(l1-5,m1,3))*exp(xi*Omega(l1-5,m1))
          dt(l1,m1,4)=dt(l1,m1,4)+MFlowM5(l1,m1)*(t(l1-5,m1)+dt(l1-5,m1,3))*exp(xi*Omega(l1-5,m1))
        end if
!
!  contribution from l-6 component
!
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
!-----------------------------------------------------------------------
!                                                                     
!  Calculate the vector velocity components at each latitude.
!  Equator is at jj=nx/2 + 2.5 due to wrap-around border
!  j=1 is centered 0.5*dtheta above and below the equator
!  j=nxhalf is centered 0.5*dtheta inside each pole
!                                                                     
!-----------------------------------------------------------------------
    do j=1,nlhalf
      jn=nlhalf+j
      js=nlhalf+1-j
      theta=0.5*pi-(j-0.5)*dtheta
      x=cos(theta)
      sintheta=sin(theta)
      rst=1.0/sintheta
!-----------------------------------------------------------------------
!                                                                     
!  calculate the spectral coefficients for wavenumber m at all x.     
!                                                                     
!-----------------------------------------------------------------------
      m=0
      m1=m+1
      unorth(m1)=0.
      usouth(m1)=0.
      vnorth(m1)=0.
      vsouth(m1)=0.
!-----------------------------------------------------------------------
!                                                                     
!  Split non-axisymmetric signal into equal positive and              
!  negative freqencies.                                               
!                                                                     
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
!
!  Construct velocity functions
!
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
!                                                                     
!  Calculate the vector velocity components at all phi positions.     
!                                                                     
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
!                                                                     
!  Write longitudinal and co-latitudinal velocity to disk file.       
!                                                                     
!-----------------------------------------------------------------------
!
! ****** Interpolate ConFlow flows into PSI grid flow arrays.
!
    wt1 = wtime()
!
! ****** VT (NT,NPM) ******
!
! ****** Set up v_ext for vt.
! ****** Since "v" is in latitude, need to flip theta axis in the process.
!
    do ii=2,nphi+1
      do jj=2,nl+1
        v_ext(ii,jj) = v(ii-1,nl+1-(jj-1))
      enddo
    enddo
!
! ****** Set "past the pole" values:
!
    pole = SUM(v_ext(2:nphi+1,2))/nphi
    v_ext(:,1) = two*pole - v_ext(:,2)
    pole = SUM(v_ext(2:nphi+1,nl+1))/nphi
    v_ext(:,nl+2) = two*pole - v_ext(:,nl+1)
!
! ****** Set periodicity:
!
    v_ext(1,:) = v_ext(nphi+1,:)
    v_ext(nphi+2,:) = v_ext(2,:)
!
! ****** Now interpolate into PSI grid (only inner grid for half-mesh in theta):
!
    call interp2d (nphi+2,nl+2,p_ext,t_ext,v_ext,                      &
                   nphi,nl-1,p_main,t_half(2:nl),vt_psi(:,2:nl),ierr)
!
! ****** Now set poles for PSI vt:
!
    pole = SUM(vt_psi(1:nphi-1,2))/(nphi-1)
    vt_psi(:,1) = two*pole - vt_psi(:,2)
    pole = SUM(vt_psi(1:nphi-1,nl))/(nphi-1)
    vt_psi(:,nl+1) = two*pole - vt_psi(:,nl)
!
! ****** VP (NTM,NP) ******
!
! ****** Set up v_ext for vp:
! ****** Since "u" is in latitude, need to flip theta axis in the process.
!
    do ii=2,nphi+1
      do jj=2,nl+1
        v_ext(ii,jj) = u(ii-1,nl+1-(jj-1))
      enddo
    enddo
!
! ****** Set "past the pole" values:
!
    pole = SUM(v_ext(2:nphi+1,2))/nphi
    v_ext(:,1) = two*pole - v_ext(:,2)
    pole = SUM(v_ext(2:nphi+1,nl+1))/nphi
    v_ext(:,nl+2) = two*pole - v_ext(:,nl+1)
!
! ****** Set periodicity:
!
    v_ext(1,:) = v_ext(nphi+1,:)
    v_ext(nphi+2,:) = v_ext(2,:)
!
! ****** Now interpolate into PSI grid (only inner grid for half-mesh in phi):
!
    call interp2d (nphi+2,nl+2,p_ext,t_ext,v_ext,                      &
                   nphi-1,nl,p_half(2:nphi),t_main,vp_psi(2:nphi,:),ierr)
!
! ****** Now set periodicity for vp:
!
    vp_psi(1,:) = vp_psi(nphi,:)
    vp_psi(nphi+1,:) = vp_psi(2,:)
!
    write(velFileNum, '(i0.6)') ifile
!
    call write_2d_file ('vt'//velFileNum//'.h5',nphi,nl+1,vt_psi,p_main,t_half,ierr)
    call write_2d_file ('vp'//velFileNum//'.h5',nphi+1,nl,vp_psi,p_half,t_main,ierr)
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
    wtime_io = wtime_io + (wtime() - wt1)
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
subroutine interp2d (nxi,nyi,xi,yi,fi,nx,ny,x,y,f,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Interpolate a 2D field from array FI(NXI,NYI), defined
! ****** on the mesh XI(NXI) x YI(NYI), into the array F(NX,NY),
! ****** defined on the mesh X(NX) x Y(NY).
!
! ****** Zero values are returned at data points outside the
! ****** bounds of the XI x YI mesh.
!
!-----------------------------------------------------------------------
!
      use number_types
      use constants
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: nxi,nyi
      real(r_typ), dimension(nxi) :: xi
      real(r_typ), dimension(nyi) :: yi
      real(r_typ), dimension(nxi,nyi) :: fi
      integer :: nx,ny
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      real(r_typ), dimension(nx,ny) :: f
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      integer :: i,j,iip1,jjp1
      integer :: ii=0,jj=0
      real(r_typ) :: dummy,ax,ay,xv,yv
!
!-----------------------------------------------------------------------
!
      real(r_typ), external :: flint
!
!-----------------------------------------------------------------------
!
      ierr=0
!
! ****** Check that the scales XI and YI are monotonic.
!
      dummy=flint(.true.,zero,nxi,xi,xi,ierr)
      dummy=flint(.true.,zero,nyi,yi,yi,ierr)
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in INTRP2D:'
        write (*,*) '### Scales are not monotonically increasing.'
        ierr=1
        return
      end if
!
! ****** Interpolate the data.
!
      do j=1,ny
        yv=y(j)
        if (yv.lt.yi(1).or.yv.gt.yi(nyi)) then
          f(:,j)=0.
          cycle
        else
          call interp (nyi,yi,yv,jj,jjp1,ay,ierr)
          if (ierr.ne.0) then
            f(:,j)=0.
            cycle
          end if
        end if
        do i=1,nx
          xv=x(i)
          if (xv.lt.xi(1).or.xv.gt.xi(nxi)) then
            f(i,j)=0.
            cycle
          else
            call interp (nxi,xi,xv,ii,iip1,ax,ierr)
            if (ierr.ne.0) then
              f(i,j)=0.
              cycle
            end if
          end if
          f(i,j)=(one-ax)*((one-ay)*fi(ii  ,jj  )+ay*fi(ii  ,jjp1)) &
                 +ax *((one-ay)*fi(iip1,jj  )+ay*fi(iip1,jjp1))
        enddo
      enddo
!
end subroutine
!#######################################################################
function flint (check,x,n,xn,fn,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Interpolate a function linearly.
!
!-----------------------------------------------------------------------
!
! ****** The funcion is defined at N nodes, with values given by
! ****** FN(N) at positions XN(N).  The function value returned is
! ****** the linear interpolant at X.
!
! ****** Note that if X.lt.XN(1), the function value returned
! ****** is FN(1), and if X.gt.XN(N), the function value returned
! ****** is FN(N).
!
! ****** Call once with CHECK=.true. to check that the values
! ****** in XN(N) are monotonically increasing.  In this mode
! ****** the array XN(N) is checked, and X and FN(N) are not
! ****** accessed.  If the check is passed, IERR=0 is returned.
! ****** Otherwise, IERR=1 is returned.
!
!-----------------------------------------------------------------------
!
      use number_types
      use constants
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      logical :: check
      real(r_typ) :: x
      integer :: n
      real(r_typ), dimension(n) :: xn,fn
      integer :: ierr
      real(r_typ) :: flint
!
!-----------------------------------------------------------------------
!
      integer :: i
      real(r_typ) :: x1,x2,alpha
!
!-----------------------------------------------------------------------
!
      ierr=0
      flint=0.
!
! ****** If CHECK=.true., check the abscissa table.
!
      if (check) then
        if (n.le.0) then
          write (*,*)
          write (*,*) '### ERROR in FLINT:'
          write (*,*) '### Invalid abscissa table dimension.'
          write (*,*) 'N = ',n
          ierr=1
          return
        end if
        do i=1,n-1
          if (xn(i+1).le.xn(i)) then
            write (*,*)
            write (*,*) '### ERROR in FLINT:'
            write (*,*) '### Abscissa table values are not'// &
                       ' monotonically increasing.'
            write (*,*) 'N = ',n
            write (*,*) 'XN = ',xn
            ierr=1
            return
          end if
        enddo
        return
      end if
!
! ****** Get the interpolated value.
!
      if (x.le.xn(1)) then
        flint=fn(1)
      else if (x.gt.xn(n)) then
        flint=fn(n)
      else
        do i=1,n-1
          if (x.ge.xn(i).and.x.lt.xn(i+1)) exit
        enddo
        x1=xn(i)
        x2=xn(i+1)
        alpha=(x-x1)/(x2-x1)
        flint=fn(i)*(one-alpha)+fn(i+1)*alpha
      end if
!
      return
end function
!#######################################################################
subroutine interp (n,x,xv,i,ip1,a,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Get the interpolation factor at XV from the table X(N).
!
!-----------------------------------------------------------------------
!
! ****** This routine does not do the actual interpolation.  Use the
! ****** returned values of I, IP1, and A to get the interpolated
! ****** value.
!
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: n
      real(r_typ), dimension(n) :: x
      real(r_typ) :: xv
      integer :: i,ip1
      real(r_typ) :: a
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      ierr=0
!
! ****** Check if the x-scale has only one point.
!
      if (n.eq.1.and.xv.eq.x(1)) then
        ip1=1
        a=0.
        return
      end if
!
! ****** Find the interval and compute the interpolation factor.
!
      do i=1,n-1
        if (xv.ge.x(i).and.xv.le.x(i+1)) then
          ip1=i+1
          if (x(i).eq.x(i+1)) then
            a=0.
          else
            a=(xv-x(i))/(x(i+1)-x(i))
          end if
          return
        end if
      enddo
!
! ****** ERROR: the value was not found.
!
      ierr=1
!
end subroutine
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
! 07/15/2022, RC, Version 0.5.0:
!   - Changed output to have vt and vp on HipFT staggered meshes.
!     This should be done directly, but for now, the conflow
!     flows are interpolated to the HipFT grids before being written out.
!   - Added constants module.
!
! 07/19/2022, RC, Version 0.5.1:
!   - BUG FIX: Flipped the theta axis in the output arrays since u and v
!              are in latitude (south to north), and vt/vp is output in
!              colatitude (north to south).
!   - BUG FIX: u and v were being written to vt and vp instead of vp and vt.
!
!-----------------------------------------------------------------------
!
!#######################################################################
