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
!     This program creates maps of supergranule flows.
!     It constructs a spectrum of spherical harmonic amplitudes
!     and calculates the components of the velocity field at
!     an array of points in theta (colatitude) and phi (longitude).
!
!     The velocity vector (u,v) is in the (phi,theta) direction.
!
!     The data arrays are written to files.
!
!     Differential rotation is simulated by evolving
!     the spectral coefficients using 5-coefficient fits
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
      character(*), parameter :: cvers='0.8.2'
      character(*), parameter :: cdate='07/10/2023'
!
end module
!#######################################################################
module timing
!
      use number_types
!
      implicit none
!
      real(r_typ) :: wtime_tmp   = 0.
      real(r_typ) :: wtime_tmp2  = 0.
      real(r_typ) :: wtime_io    = 0.
      real(r_typ) :: wtime_vcalc = 0.
      real(r_typ) :: wtime_scalc = 0.
      real(r_typ) :: wtime_setup = 0.
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
      real(r_typ), parameter :: rSun = 6.96e+08_r_typ
      real(r_typ), parameter :: day_to_sec = 86400.0_r_typ
!
      real(r_typ), parameter :: root3 = sqrt(three)
      real(r_typ), parameter :: root5 = sqrt(five)
      real(r_typ), parameter :: root7 = sqrt(seven)
      real(r_typ), parameter :: root9 = sqrt(nine)
!
      real(r_typ), parameter :: small_value = TINY(one)
      real(r_typ), parameter :: large_value = HUGE(one)
!
end module
!#######################################################################
module input_parameters
!
!-----------------------------------------------------------------------
! ****** Input parameters.
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
      logical :: verbose = .false.
!
! ****** Output options.
!
      logical :: output_flows_on_uniform_grid = .true.
      logical :: output_flows_on_staggered_grid = .true.
      logical :: output_spectrum = .true.
      character(512) :: output_directory = '.'
!
! ****** Resolution.
!
      integer :: n_lat = 512
      integer :: n_long = 1024
!
! ****** Time paramters (seconds).
!
      real(r_typ) :: tmax = 24*28*3600
      real(r_typ) :: dtime = fifteen*60.0_r_typ
!
! ****** Differential rotation coefficients (m/s).
!
      real(r_typ) :: flow_dr_t0 = 46.0_r_typ
      real(r_typ) :: flow_dr_t1 = zero
      real(r_typ) :: flow_dr_t2 = -262.0_r_typ
      real(r_typ) :: flow_dr_t3 = zero
      real(r_typ) :: flow_dr_t4 = -379.0_r_typ
!
! ****** Meridional flow coefficients (m/s).
!
      real(r_typ) :: flow_mf_s0 = zero
      real(r_typ) :: flow_mf_s1 = 22.0_r_typ
      real(r_typ) :: flow_mf_s2 = zero
      real(r_typ) :: flow_mf_s3 = 11.0_r_typ
      real(r_typ) :: flow_mf_s4 = zero
      real(r_typ) :: flow_mf_s5 = -28.0_r_typ
!
! ****** Random seed options.
!
      logical :: set_random_seed = .true.
      integer :: random_seed_value = 12345
!
! ****** Spectrum taper options.
!
      integer :: spectrum_taper_model = 1
      real(r_typ) :: spectrum_taper_val1 = 180
      real(r_typ) :: spectrum_taper_val2 = 204
!           (1) Original taper
!           (2) Raphael tamper
!           (3) Ron 1 
!           (4) Ron 2
!           (5) Original with cut-off
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
  use input_parameters
!
!-----------------------------------------------------------------------
!
  implicit none
!
!-----------------------------------------------------------------------
!
! ****** Main arrays.
!
  real(r_typ) :: coefDR(5),coefMF(6)
!
  real(r_typ),    dimension(:,:),   allocatable :: vt_psi,vp_psi,v_ext
  real(r_typ),    dimension(:,:),   allocatable :: u,v,coef,Omega,MFlow0
  real(r_typ),    dimension(:,:),   allocatable :: OmegaM1,OmegaM2
  real(r_typ),    dimension(:,:),   allocatable :: OmegaM4,OmegaM3
  real(r_typ),    dimension(:,:),   allocatable :: OmegaP1,OmegaP2
  real(r_typ),    dimension(:,:),   allocatable :: OmegaP4,OmegaP3
  real(r_typ),    dimension(:,:),   allocatable :: MFlowM1,MFlowM2,MFlowM3
  real(r_typ),    dimension(:,:),   allocatable :: MFlowM4,MFlowM5,MFlowM6
  real(r_typ),    dimension(:,:),   allocatable :: MFlowP1,MFlowP2,MFlowP3
  real(r_typ),    dimension(:,:),   allocatable :: MFlowP4,MFlowP5,MFlowP6
  real(r_typ),    dimension(:,:),   allocatable :: s_abs
  real(r_typ),    dimension(:),     allocatable :: latscale
  complex(r_typ), dimension(:,:),   allocatable :: s,t
  complex(r_typ), dimension(:,:,:), allocatable :: ds,dt
  complex(r_typ), dimension(:),     allocatable :: unorth,usouth
  complex(r_typ), dimension(:),     allocatable :: vnorth,vsouth
  integer,        dimension(:),     allocatable :: seed_old,seed_new
  real(r_typ),    dimension(:),     allocatable :: p_main,p_half,p
  real(r_typ),    dimension(:),     allocatable :: t_main,t_half
  real(r_typ),    dimension(:),     allocatable :: p_ext,t_ext
!
  character(512) :: fname = ' '
  character(  6) :: velFileNum = '000000'
!
  integer :: seed_n
  integer :: i,j,l,m,l1,m1,l2,m2,jn,js,ierr,ii,jj
  integer :: n_lat_2x, n_long_2x
  integer :: lmax,nlhalf,ifile,nfiles,itime
  complex(r_typ) :: xi,arg,sum1,sum2,sum3,sum4
  real(r_typ) :: n_lat_2x_real
  real(r_typ) :: dphi,dtheta,theta,sintheta,x,rst,eo,v1,v2
  real(r_typ) :: el,em,randamp,phase,xlifetime
  real(r_typ) :: elfunc,elfuncMF,DR0,DR1,DR2,DR3,DR4,MF0,MF1,MF2,MF3,MF4
  real(r_typ) :: BAALM4,BAALM3,BAALM2,BAALM1,BAALP0,BAALP1,BAALP2,BAALP3
  real(r_typ) :: BAA0,BAA1,BAA2,BAA3,BAA4,BAA5,BAA6,BAALP4,BAALP5,MF5
  real(r_typ) :: dphi_psi,dtheta_psi,pole
  real(r_typ) :: wtime,wt1,rnd,curr_time
  real(r_typ) :: ampS,ampT,taper_l0,taper_l1,taper
!
!-----------------------------------------------------------------------
!
! ****** Start wall clock timer.
!
  wtime_tmp = wtime()
!
!-----------------------------------------------------------------------
!
! ****** READ INPUT FILE.
!
  call read_input_file
!
  call write_welcome_message
!
! ****** Set grid resolutions.
!
  n_long_2x     = 2*n_long
  n_lat_2x_real = real(n_lat_2x,r_typ)
  n_lat_2x      = 2*n_lat
  nlhalf        = n_lat_2x/2
  lmax          = n_lat_2x-1
  dphi          = twopi/n_long_2x
  dtheta        = pi/n_lat_2x
! ****** Allocate arrays.
!
  allocate (               p(n_lat_2x))
  allocate (     u(n_long_2x,n_lat_2x))
  allocate (     v(n_long_2x,n_lat_2x))
  allocate (   coef(n_lat_2x,n_lat_2x))
!
  allocate (  Omega(n_lat_2x,n_lat_2x))
  allocate (OmegaM1(n_lat_2x,n_lat_2x))
  allocate (OmegaM2(n_lat_2x,n_lat_2x))
  allocate (OmegaM3(n_lat_2x,n_lat_2x))
  allocate (OmegaM4(n_lat_2x,n_lat_2x))
!
  allocate (OmegaP1(n_lat_2x,n_lat_2x))
  allocate (OmegaP2(n_lat_2x,n_lat_2x))
  allocate (OmegaP3(n_lat_2x,n_lat_2x))
  allocate (OmegaP4(n_lat_2x,n_lat_2x))
!
  allocate ( MFlow0(n_lat_2x,n_lat_2x))
  allocate (MFlowM1(n_lat_2x,n_lat_2x))
  allocate (MFlowM2(n_lat_2x,n_lat_2x))
  allocate (MFlowM3(n_lat_2x,n_lat_2x))
  allocate (MFlowM4(n_lat_2x,n_lat_2x))
  allocate (MFlowM5(n_lat_2x,n_lat_2x))
  allocate (MFlowM6(n_lat_2x,n_lat_2x))
!
  allocate (MFlowP1(n_lat_2x,n_lat_2x))
  allocate (MFlowP2(n_lat_2x,n_lat_2x))
  allocate (MFlowP3(n_lat_2x,n_lat_2x))
  allocate (MFlowP4(n_lat_2x,n_lat_2x))
  allocate (MFlowP5(n_lat_2x,n_lat_2x))
  allocate (MFlowP6(n_lat_2x,n_lat_2x))
!
  allocate (s(n_lat_2x,n_lat_2x))
  allocate (ds(n_lat_2x,n_lat_2x,4))
  allocate (t(n_lat_2x,n_lat_2x))
  allocate (dt(n_lat_2x,n_lat_2x,4))
!
  allocate (unorth(n_long_2x))
  allocate (usouth(n_long_2x))
  allocate (vnorth(n_long_2x))
  allocate (vsouth(n_long_2x))
!
! ****** Initialize random number generator.
!
  call RANDOM_INIT (.true., .true.)

  if (set_random_seed) then
    call RANDOM_SEED (size=seed_n)
    allocate (seed_old(seed_n))
    allocate (seed_new(seed_n))
    call RANDOM_SEED(get=seed_old)
    seed_new(:)=random_seed_value
    call RANDOM_SEED(put=seed_new)
  end if
!
!-----------------------------------------------------------------------
!
!  Numeric constants
!
!-----------------------------------------------------------------------
  xi=(0.,one)
!-----------------------------------------------------------------------
!
!  Time (seconds) and space (radians) steps
!
!-----------------------------------------------------------------------
  curr_time = 0.
  nfiles = int(tmax/dtime)
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
!
  allocate (v_ext(n_long_2x+2,n_lat_2x+2))
  allocate (p_ext(n_long_2x+2))
  allocate (t_ext(n_lat_2x+2))
!
! ****** Set extended ConFlow meshes.
!
  do i=1,n_long_2x+2
    p_ext(i) = -half*dphi + (i-1)*dphi
  enddo
  do i=1,n_lat_2x+2
    t_ext(i) = -half*dtheta + (i-1)*dtheta
  enddo
!  print*,p_ext(1),t_ext(1),p_ext(n_long_2x+2),t_ext(n_lat_2x+2)
!
! ***** Allocate psi output arrays.
!
  allocate (vt_psi(n_long,n_lat+1))
  allocate (vp_psi(n_long+1,n_lat))
!
! ****** Set PSI meshes.
!
  allocate (p_main(n_long))
  allocate (p_half(n_long+1))
  allocate (t_main(n_lat))
  allocate (t_half(n_lat+1))
!
  dphi_psi = twopi/(n_long-1)
  dtheta_psi = pi/(n_lat-1)
!
  do i=1,n_long
    p_main(i) = (i-1)*dphi_psi
  enddo
  do i=1,n_long+1
    p_half(i) = -half*dphi_psi+(i-1)*dphi_psi
  enddo
  do i=1,n_lat
    t_main(i) = (i-1)*dtheta_psi
  enddo
  do i=1,n_lat+1
    t_half(i) = -half*dtheta_psi + (i-1)*dtheta_psi
  enddo
  print*,p_main(1),t_main(1),p_main(n_long),t_main(n_lat)
  print*,p_half(1),t_half(1),p_half(n_long+1),t_half(n_lat+1)
!
!-----------------------------------------------------------------------
!
!  Differential Rotation coefficients m/s relative to Carrington
!
!-----------------------------------------------------------------------
  write(*,*)
  write(*,*) 'Differential rotation coeffs (m/s) flow_dr_t0,flow_dr_t1,flow_dr_t2,flow_dr_t3,flow_dr_t4:'
  write(*,*) flow_dr_t0,flow_dr_t1,flow_dr_t2,flow_dr_t3,flow_dr_t4
  coefDR(1)=flow_dr_t0*dtime/rSun
  coefDR(2)=flow_dr_t1*dtime/rSun
  coefDR(3)=flow_dr_t2*dtime/rSun
  coefDR(4)=flow_dr_t3*dtime/rSun
  coefDR(5)=flow_dr_t4*dtime/rSun
!-----------------------------------------------------------------------
!
! Meridional flow speed of supergranules in m/s poleward
!
!-----------------------------------------------------------------------
  write(*,*)
  write(*,*) 'Meridional flow coeffs (m/s) flow_mf_s0,flow_mf_s1,flow_mf_s2,flow_mf_s3,flow_mf_s4,flow_mf_s5:'
  write(*,*) flow_mf_s0,flow_mf_s1,flow_mf_s2,flow_mf_s3,flow_mf_s4,flow_mf_s5
  coefMF(1)=flow_mf_s0*dtime/rSun
  coefMF(2)=flow_mf_s1*dtime/rSun
  coefMF(3)=flow_mf_s2*dtime/rSun
  coefMF(4)=flow_mf_s3*dtime/rSun
  coefMF(5)=flow_mf_s4*dtime/rSun
  coefMF(6)=flow_mf_s5*dtime/rSun
!-----------------------------------------------------------------------
!
!  Create coefficients for spherical harmonics.
!
!-----------------------------------------------------------------------
  call plmcoef(n_lat_2x,coef)
!
  write(*,*)
  write(*,*) 'Legendre recurrence coefficients calculated.'
!-----------------------------------------------------------------------
!
!  Construct the convection spectrum.
!
!-----------------------------------------------------------------------
!
  do l=1,lmax
    l1=l+1
    el=real(l,r_typ)
!
    ! From Raphael: Hathaway's Spectrum inconsistent
    !               with Hathaway et al. 2010 and earlier.
!
    ampS = 0.08_r_typ*(one - tanh(el/165.0_r_typ)) + 0.0024_r_typ*(one - tanh(el/2000._r_typ))
    ampT = 1.5_r_typ*(one - half*sqrt(el/1000.0_r_typ))/el
!
    select case (spectrum_taper_model) 
      case (1)   ! Original ConFlow
        taper_l0 = 384.0_r_typ
        taper_l1 = 512.0_r_typ
        taper = one
        if (el .gt. taper_l0) taper = half*(one + cos(pi*(el - taper_l0)/(taper_l1 - taper_l0)))  
      case (2)   ! Raphel v1
        taper_l0 = 200.0_r_typ
        taper_l1 = 200.0_r_typ
        ampS = 0.08_r_typ*(one - tanh(el/300.0_r_typ))
        taper = half*(one + cos(pi*el/taper_l1))
      case (3)  ! Ron v1
        taper_l0 = spectrum_taper_val1
        taper_l1 = spectrum_taper_val2
        taper = one
        if (el .gt. taper_l0) taper = (half*(one + cos(pi*(el - taper_l0)/(taper_l1 - taper_l0))))**(0.2_r_typ)
      case (4)  ! Ron v2
        taper_l0 = spectrum_taper_val1
        taper_l1 = spectrum_taper_val2
        taper = (half*(one + cos(pi*el/taper_l1)))**(0.05_r_typ)
      case (5)  ! Original ConFlow with hard cut-off
        taper_l0 = zero
        taper_l1 = spectrum_taper_val2
        taper = one
      case default  ! Original ConFlow 
        taper_l0 = 384.0_r_typ
        taper_l1 = 512.0_r_typ
        taper = one
        if (el .gt. taper_l0) taper = half*(one + cos(pi*(el - taper_l0)/(taper_l1 - taper_l0)))  
    end select
!
    if (el .gt. taper_l1) taper = zero  
!
    ampS = taper*ampS
    ampT = taper*ampT
!
    do m=1,l
      m1=m+1
      call RANDOM_NUMBER (rnd)
      phase = two*pi*rnd
      arg = cos(phase) + xi*sin(phase)
      call RANDOM_NUMBER (rnd)
      randamp = 1.8*rnd
      s(l1,m1) = randamp*ampS*arg
      call RANDOM_NUMBER (rnd)
      phase = two*pi*rnd
      arg = cos(phase) + xi*sin(phase)
      call RANDOM_NUMBER (rnd)
      randamp = 1.8*rnd
      t(l1,m1) = randamp*ampT*arg
    end do
  end do
!
! ****** Write spectrum to file if desired.
!
  if (output_spectrum) then
!
    wt1 = wtime()
!
    allocate (s_abs(n_lat_2x,n_lat_2x))
    allocate (latscale(n_lat_2x))
!
    do concurrent (l=1:n_lat_2x)
      latscale(l) = real(l,r_typ)
    enddo
!
    do concurrent (l=1:n_lat_2x,m=1:n_lat_2x)
      s_abs(l,m) = ABS(s(l,m))
    enddo
!
    fname = TRIM(output_directory)//'/'//'spectrum.h5'
    call write_2d_file (TRIM(fname),n_lat_2x,n_lat_2x,   &
                        s_abs,latscale,latscale,ierr)
    print*, "Wrote spectrum to: "//fname
!
    deallocate (s_abs)
    deallocate (latscale)
!
    wtime_io = wtime_io + (wtime() - wt1)
!
  end if
!
  write(*,*) 'Velocity spectrum calculated.'
!
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
!
  wtime_setup = wtime() - wtime_tmp
!
!-----------------------------------------------------------------------
!
!  Construct series of velocity maps at dtime (in seconds) intervals
!
!-----------------------------------------------------------------------
!
  write(*,*)
  write(*,*) 'A time step of ', dtime, 'seconds will produce ', &
             nfiles,' vt and vp flow files.'
  write(*,*)
  flush(OUTPUT_UNIT)
!
! ****** Open text file to list output maps.
!
  call ffopen (12,TRIM(output_directory)//'/flow_output_list.csv','rw',ierr)
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
!           in Hathaway 2010.
!
! Step 1
!
!-----------------------------------------------------------------------
    wt1 = wtime()
!
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
    wtime_scalc = wtime_scalc + (wtime() - wt1)
!
! ****** Add changes and evolve
!
    wt1 = wtime()
!
    do m=1,lmax-1
      m1=m+1
      do l=m,lmax-1
        l1=l+1
        el=real(l,r_typ)
        xlifetime=10.*3600.*(100./(el+0.1))**2.0
        !        xlifetime=xlifetime/2.
        call RANDOM_NUMBER (rnd)
        phase=Omega(l1,m1)+2.*(rnd-0.5)*sqrt(dtime/xlifetime)
        arg=cos(phase)+xi*sin(phase)
        s(l1,m1)=(s(l1,m1) + ds(l1,m1,1)/6. + ds(l1,m1,2)/3. + ds(l1,m1,3)/3. + ds(l1,m1,4)/6.)*arg

        call RANDOM_NUMBER (rnd)
        phase=Omega(l1,m1)+2.*(rnd-0.5)*sqrt(dtime/xlifetime)
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
!  negative frequencies.
!
!-----------------------------------------------------------------------
      do m=1,lmax-1
        m1=m+1
        m2=n_long_2x+1-m
        call plm(m,x,n_lat_2x,coef,p)
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
!
      do concurrent (m=lmax:n_long_2x-lmax)
        unorth(m+1)=0.
        usouth(m+1)=0.
        vnorth(m+1)=0.
        vsouth(m+1)=0.
      enddo
!-----------------------------------------------------------------------
!
!  Calculate the vector velocity components at all phi positions.
!
!-----------------------------------------------------------------------
      call four1(unorth,n_long_2x,-1)
      call four1(usouth,n_long_2x,-1)
      call four1(vnorth,n_long_2x,-1)
      call four1(vsouth,n_long_2x,-1)
!
! ****** Get real part of complex values from FFT.
!
      do concurrent (i=1:n_long_2x)
        u(i,jn)=real(unorth(i),r_typ)
        u(i,js)=real(usouth(i),r_typ)
        v(i,jn)=real(vnorth(i),r_typ)
        v(i,js)=real(vsouth(i),r_typ)
      end do
!
    end do ! End latitude-loop
!
    wtime_vcalc = wtime_vcalc + (wtime() - wt1)
!-----------------------------------------------------------------------
!
!  Write longitudinal and co-latitudinal velocity to disk file.
!
!-----------------------------------------------------------------------
!
    wt1 = wtime()
!
    if (output_flows_on_uniform_grid) then
      write(velFileNum, '(i0.6)') ifile
      fname='vp' // velFileNum
      open(unit=1,file=TRIM(output_directory)//'/'//TRIM(fname)//'.data',access='direct', &
           status='unknown',recl=8*n_long_2x)
        do j=1,n_lat_2x
          write(1,rec=j) (u(i,j),i=1,n_long_2x)
        end do
      close(1)
!
      fname='vt' // velFileNum
      open(unit=1,file=TRIM(output_directory)//'/'//TRIM(fname)//'.data',access='direct', &
           status='unknown',recl=8*n_long_2x)
        do j=1,n_lat_2x
          write(1,rec=j) (v(i,j),i=1,n_long_2x)
        enddo
      close(1)
    end if ! end of Dave's i/o block

    if (output_flows_on_staggered_grid) then
!
! ****** Interpolate ConFlow flows into PSI grid flow arrays.
!
! ****** VT (NT,NPM) ******
!
! ****** Set up v_ext for vt.
! ****** Since "v" is in latitude, need to flip theta axis in the process.
!
      do ii=2,n_long_2x+1
        do jj=2,n_lat_2x+1
          v_ext(ii,jj) = v(ii-1,n_lat_2x+1-(jj-1))
        enddo
      enddo
!
! ****** Set "past the pole" values:
!
      pole = SUM(v_ext(2:n_long_2x+1,2))/n_long_2x
      v_ext(:,1) = two*pole - v_ext(:,2)
!
      pole = SUM(v_ext(2:n_long_2x+1,n_lat_2x+1))/n_long_2x
      v_ext(:,n_lat_2x+2) = two*pole - v_ext(:,n_lat_2x+1)
!
! ****** Set periodicity:
!
      v_ext(1,:) = v_ext(n_long_2x+1,:)
      v_ext(n_long_2x+2,:) = v_ext(2,:)
!
! ****** Now select points on PSI grid (only inner grid for half-mesh in theta):
!
      vt_psi(1:n_long, 2:n_lat) = v_ext(1:n_long_2x:2,2:n_lat_2x-2:2)
!
! ****** Now set poles for PSI vt:
!
      pole = SUM(vt_psi(1:n_long-1,2))/(n_long-1)
      vt_psi(:,1) = two*pole - vt_psi(:,2)
!
      pole = SUM(vt_psi(1:n_long-1,n_lat))/(n_long-1)
      vt_psi(:,n_lat+1) = two*pole - vt_psi(:,n_lat)
!
! ****** VP (NTM,NP) ******
!
! ****** Set up v_ext for vp:
! ****** Since "u" is in latitude, need to flip theta axis in the process.
!
      do ii=2,n_long_2x+1
        do jj=2,n_lat_2x+1
          v_ext(ii,jj) = u(ii-1,n_lat_2x+1-(jj-1))
        enddo
      enddo
!
! ****** Set "past the pole" values:
!
      pole = SUM(v_ext(2:n_long_2x+1,2))/n_long_2x
      v_ext(:,1) = two*pole - v_ext(:,2)
      pole = SUM(v_ext(2:n_long_2x+1,n_lat_2x+1))/n_long_2x
      v_ext(:,n_lat_2x+2) = two*pole - v_ext(:,n_lat_2x+1)
!
! ****** Set periodicity:
!
      v_ext(1,:) = v_ext(n_long_2x+1,:)
      v_ext(n_long_2x+2,:) = v_ext(2,:)
!
! ****** Now select points on PSI grid (only inner grid for half-mesh in phi):
!
      vp_psi(2:n_long,1:n_lat) = v_ext(2:n_long_2x-2:2,1:n_lat_2x:2)
!
! ****** Now set periodicity for vp:
!
      vp_psi(1,:) = vp_psi(n_long,:)
      vp_psi(n_long+1,:) = vp_psi(2,:)
!
      write(velFileNum, '(i0.6)') ifile
!
      call write_2d_file (TRIM(output_directory)//'/'//'vt'//velFileNum//'.h5',n_long,n_lat+1,vt_psi,p_main,t_half,ierr)
      call write_2d_file (TRIM(output_directory)//'/'//'vp'//velFileNum//'.h5',n_long+1,n_lat,vp_psi,p_half,t_main,ierr)
!
      call ffopen (12,TRIM(output_directory)//'/'//'flow_output_list.csv','a',ierr)
      write (12,'(F11.5,A1,A11,A1,A11)') curr_time/(3600*24),',',     &
                                trim('vt'//velFileNum//'.h5'),',',    &
                                trim('vp'//velFileNum//'.h5')
      close(12)
!
      write(*,*) 'Completed step ',ifile,' of ', nfiles,              &
                 ' max(vt)=',maxval(v),' min(vt)=',minval(v),         &
                 ' max(vp)=',maxval(u),' min(vp)=',minval(u)
      flush(OUTPUT_UNIT)
    end if ! end of PSI i/o block
!
    wtime_io = wtime_io + (wtime() - wt1)
!
! ***** Update time.
!
    curr_time = curr_time + dtime
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
subroutine write_welcome_message
!
!-----------------------------------------------------------------------
!
! ****** Write welcome message
!
!-----------------------------------------------------------------------
!
      use ident
      use input_parameters
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      write (*,*) ''
      write (*,*) '       ______            _______ _'
      write (*,*) '      / _____)          (_______) |'
      write (*,*) '     | /      ___  ____  _____  | | ___  _ _ _'
      write (*,*) '     | |     / _ \|  _ \|  ___) | |/ _ \| | | |'
      write (*,*) '     | \____| |_| | | | | |     | | |_| | | | |'
      write (*,*) '      \______)___/|_| |_|_|     |_|\___/ \____|'
      write (*,*) ''
      write (*,*) ''
      write (*,*) ' ****** ConFlow: Super Granular Convective Flow Generator'
      write (*,*) ''
      write (*,*) '     Version: ',cvers,' of ',cdate
      write (*,*) ''
      write (*,*) '     Authors:  Raphael Attie'
      write (*,*) '               Ronald M. Caplan'
      write (*,*) '               David H. Hathaway'
      write (*,*) ''
!
! ****** [RMC] Update this with relevent info in a nice way.
!
      write (*,*) ''
!
      FLUSH(OUTPUT_UNIT)
!
end subroutine
!#######################################################################
subroutine write_2d_file (fname,ln1,ln2,f,flow_mf_s1,flow_mf_s2,ierr)
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
      real(r_typ), dimension(ln1) :: flow_mf_s1
      real(r_typ), dimension(ln2) :: flow_mf_s2
      integer :: ln1,ln2
!
!-----------------------------------------------------------------------
!
      type(ds) :: s
      integer :: ierr
!
!-----------------------------------------------------------------------
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
      s%scales(1)%f(:) = flow_mf_s1(:)
      s%scales(2)%f(:) = flow_mf_s2(:)
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
      write(*,FMT) "Wall clock time:      ",wtime_total
      write(*,FMT) "--> SETUP:         ",wtime_setup
      write(*,FMT) "--> COMP Spectrum: ",wtime_scalc
      write(*,FMT) "--> COMP Velocity: ",wtime_vcalc
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
subroutine read_input_file
!
!-----------------------------------------------------------------------
!
! ****** Read the input file.
!
!-----------------------------------------------------------------------
!
      use input_parameters
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      namelist /conflow_input_parameters/                             &
               verbose, output_flows_on_uniform_grid,                 &
               output_flows_on_staggered_grid, output_spectrum,       &
               flow_dr_t0, flow_dr_t1, flow_dr_t2, flow_dr_t3,        &
               flow_dr_t4, flow_mf_s0, flow_mf_s1, flow_mf_s2,        &
               flow_mf_s3, flow_mf_s4, flow_mf_s5, set_random_seed,   &
               random_seed_value, output_directory, n_lat, n_long,    &
               tmax, dtime, spectrum_taper_model, spectrum_taper_val1,&
               spectrum_taper_val2
!
!-----------------------------------------------------------------------
!
      integer :: ierr
      character(80) :: infile='conflow.dat'
!
!-----------------------------------------------------------------------
!
! ****** Read the input file.
!
      call ffopen (8,TRIM(infile),'r',ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READ_INPUT_FILE:'
        write (*,*) '### Could not open the input file.'
        write (*,*) 'File name: ',trim(infile)
        STOP
      end if
!
      read (8,conflow_input_parameters)
      close (8)
!
! ****** Add input parameter checks here.
!
end subroutine
!#######################################################################
subroutine ffopen (iun,fname,mode,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Open file FNAME and link it to unit IUN.
!
! ****** If there is an error, this routine returns IERR.ne.0.
!
!-----------------------------------------------------------------------
!
! ****** When MODE='r', the file must exist.
! ****** When MODE='w', the file is created.
! ****** When MODE='rw', the file must exist, but can be overwritten.
! ****** When MODE='a', the file is created if it does not exist,
! ******                otherwise, it is appended.
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: iun
      character(*) :: fname
      character(*) :: mode
      integer :: ierr
      logical :: ex
!
!-----------------------------------------------------------------------
!
      ierr=0
!
      if (mode.eq.'r') then
        open (iun,file=fname,form="FORMATTED",status='old',err=900)
      else if (mode.eq.'rw') then
        open (iun,file=fname,form="FORMATTED",status='replace',err=900)
      else if (mode.eq.'w') then
        open (iun,file=fname,form="FORMATTED",status='new',err=900)
      elseif (mode.eq.'a') then
        inquire(file=fname, exist=ex)
        if (ex) then
          open (iun,file=fname,form="FORMATTED",position='append',err=900)
        else
          open (iun,file=fname,form="FORMATTED",status='new',err=900)
        end if
      else
        write (*,*)
        write (*,*) '### ERROR in FFOPEN:'
        write (*,*) '### Invalid MODE requested.'
        write (*,*) 'MODE = ',mode
        write (*,*) 'File name: ',trim(fname)
        ierr=2
        return
      end if
!
      return
!
  900 continue
!
      write (*,*)
      write (*,*) '### ERROR in FFOPEN:'
      write (*,*) '### Error while opening the requested file.'
      write (*,*) 'File name: ',trim(fname)
      write (*,*) 'MODE = ',mode
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
! 03/17/2023, RA/RC, Version 0.6.0:
!   - Added command line argument to indicate whether to
!     output flows to the original AFT-style grid (nonstag),
!     the PSI HipFT grid (stag), or both.
!     Set the input argument to: stag|nonstag|both
!   - Changed spectrum in order to have flows
!     resolved at the 512x1024 resoltion with at least 5
!     grid cells for the flow blobs.
!
! 03/20/2023, RC, Version 0.7.0:
!   - Lots of refactoring.
!   - Added namelist to specify all input parameters.
!   - Changed output of spectra to hdf5.
!
! 03/24/2023, RC, Version 0.7.1:
!   - Added new timers
!   - Some cosmetic changes.
!
! 04/26/2023, RC, Version 0.8.0:
!   - Added spectrum_taper_model input parameter to select type of 
!     spectrum model to use (5 options currently).
!
! 06/02/2023, RC, Version 0.8.1:
!   - Changed spectra cutoff in option 5 to 170 instead of 200.
!
! 07/10/2023, RC, Version 0.8.2:
!   - Added some taper parameters to input file.
!     Now onw can set the cutoff for method 5.
!
!-----------------------------------------------------------------------
!
!#######################################################################
