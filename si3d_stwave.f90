!************************************************************************
                          MODULE si3d_stwave
!************************************************************************
!
!  Purpose: Procedures that implement routines related to the modelling
!           of ecological processees
!
!-------------------------------------------------------------------------

  USE si3d_types
  USE omp_lib

  IMPLICIT NONE
  SAVE

CONTAINS

!        d1_stw   = constant depth along one side of the grid when BC types are 0,1,3.
!  df_const_stw   = delta_stw frequency constant.
!        dx   = grid spacing (rectangular) in meters in x direction.
!        dy   = grid spacing (rectangular) in meters in y direction.
!      dadd_stw   = water level adjustment
! etot_stw(:,:)   = hm0
!        f0_stw   = Spectral frequencies begin at f0_stw (minimum value).
!        fm_stw   = peak frequency
! fma_stw(:, :)   = Spectral Peak frequency redefined for local wind growth cases.
!  fp(:, :)   = Spectral Peak frequency for each grid location (i_stw, j_stw).
! Tmm1_stw(:,:)   = Tm-01, mean period.
!         g   = gravitational acceleration
! i_bc_stw(1-4)   = Boundary Condition Type options for grid sides 1 through 4:
!             = 0 - constant spectrum set equal to zero
!             = 1 - constant TMA spectrum with H_stw, T, and dir specified
!             = 2 - constant spectrum read from a .eng file (file code 15)
!             = 3 - 1D transformed spectrum (adjacent boundary must be 0, 1, or 2)
!     ibnd_stw    = 0  - read single spectrum for each time_stw step                                
!             = 1  - read spectra from offshore output & linearly interpolate to boundary   
!             = 2  - read spectra from offshore output & morphicly interpolate to boundary  
!    ibreak_stw   = 0  - don't write breaker indices 
!             = 1  - write breaker indices to file break
!      icur_stw   = 0  - no current
!             = 1  - currents
!             = 2  - current field constant for all input spectra
!     ifric_stw   = 1 (JONSWAP constant) bottom friction
!             = 2 (JONSWAP variable) bottom friction
!             = 3 (Manning constant) bottom friction
!             = 4 (Manning variable) bottom friction
!    isurge_stw   = 0  - read constant depth adjustment from spectral input file
!             = 1  - read spatially variable depth adjustment from unit 14 (fort.14)
!     iwind_stw   = 0  - read constant wind speed and direction from spectral input file
!             = 1  - read spatially variable depth adjustment from unit 20 (fort.20)
!       idd_stw   = integer identification label (such as a date-time_stw code)
!  idep_opt_stw   = Read options for bathymetry:
!             = 0  - read bathymetry from DepFile (default), 
!             = 1  - plane slope. Read grid parameters from OptsFile i_side_stw =1,4 (deep end of grid), 
!   iout_stw(:)   = i_stw-location of nnth special output point.
!   jout_stw(:)   = j_stw-location of nnth special output point.
!    iplane_stw   = 0  - half-plane version
!             = 1  - full-plane version
!    i_side_stw   = Grid sides are number 1-4 to assist in grid sweep direction.
!             = 1 = West side of grid.       3 = East side of grid.
!             = 2 = South side of grid.      4 = North side of grid.          
!i_windside_stw   = which side of grid the wind input is on.
!      iprp_stw   = 0  - source terms
!             = 1  - no source terms
!       irs_stw   = 0  - don't calculate radiation stresses
!             = 1  - calculate radiation stresses and write to radstress
!  iter_sect_stw   = number of times to process wave propagation sweeps. Until it converges, 1-10 iter_stw.
!   jout_stw(:)   = j_stw-location of nnth special output point.
!     kkmax   = for each grid location (i_stw, j_stw), kkmax identifies the spectral
!             peak frequency (freq_stw(kkmax) (integrated over all angles).
!   kmax_stw(:)   = kmax_stw(1, nj_stw) stores kkmax for each element in the grid
!             column of the current grid row(i_stw). kmax_stw(:) is redefined
!             for subsequent grid rows. wkpeak_stw(:, :, :) requires kmax_stw(:).
!   m1_stw - m4_stw   = Identifies quadrant of grid being swept (m1_stw=first quadrant, etc).
!    m_wind_stw   = sweep order of grid sides including wind direction
!    n_wind_stw   = the angle band of wind direction.
!        na_stw   = number of angles in spectrum (=35 by default)
!      nest_stw   = number of output nesting points
!   nest_in   = number of nesting input points
!     nfreq_stw   = number of frequencies in spectrum
!        ni_stw   = number of columns in grid
!        nj_stw   = number of rows in grid
!    nselct_stw   = number of special output points 
!        pi   = 0.5e0*twopi_stw or (approx 3.141592654)
!    radfac_stw   = degees to radian conversion factor
!     slope   = slope input when BC types are 0,1,3.  
!     twopi_stw   = 8.0e0*atan(1.0e0) or (approx 6.283185307)       
!         uairStwave     = wind speed in m_stw/sec
! uairdirStwave = wind direction (vector heading) in degrees (stwave coord sys)
! wangle_stw(:,:)   = average wave angle
!___________________________________________________________________________________________________
!ADDED BY PATRICIO MORENO (04/2010)
!___________________________________________________________________________________________________
!  tau_stwave1(:,:)  = Shield stress magnitude at the bottom bed, value calculated per gridcell
!  taux_stwave(:,:)  = Shield stress component in the x-direction, at the bottom bed, calculated per gridcell
!  tauy_stwave(:,:)  = Shield stress component in the y-direction, at the bottom bed, calculated per gridcell

! ********************************************************************
SUBROUTINE stwave_input(n)
! ********************************************************************
!
! Purpose: Estimate the mean for wind speed and direction. This is 
!          then used in stwave for bottom shear stress and wave 
!          height calculations.
!
! --------------------------------------------------------------------
  ! Subroutine arguments
  integer, intent(in) :: n
  integer :: i,j,l,iteration,liter

  real, dimension(jm1,im1) :: u_stwave
  real, dimension(jm1,im1) :: udir_stwave
  real, dimension(jm1,im1) :: tau_stwave1
  real, dimension(jm1,im1) :: uair_stwave
  real, dimension(jm1,im1) :: vair_stwave

  ! Initialize the wind arrays
  uair_stwave(:, :) = 0.0
  vair_stwave(:, :) = 0.0

  select case (ifsurfbc)
    case (0)
      uair_stwave = -wa * sin(pi*phi/180)
      vair_stwave = -wa * cos(pi*phi/180)
    case (1:3,20)
      uair_stwave(:, :) = uair(1)
      vair_stwave(:, :) = vair(1)
    case (10:11)
      do liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
        l = id_column(liter)
        i = l2i(l)
        j = l2j(l)
        uair_stwave(j, i) = uair(l)
        vair_stwave(j, i) = vair(l)
      end do
  end select

  !$omp barrier

  do liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
    l = id_column(liter)

    if (l == 1) then
      if (n == 1) then
        uair_tmp = 0.0
        udir_tmp = 0.0
        iteration = 1
      elseif (mod(n,int(Ti_4_stwave*3600/dt)+1) == 0) then 
        iteration = 1
      elseif (mod(n,int(Ti_4_stwave*3600/dt)) == 0) then
        iteration = mod(n,int(Ti_4_stwave*3600/dt)+1)
      else
        iteration = mod(n,int(Ti_4_stwave*3600/dt))
      end if

      do i = 1, im1
        do j = 1, jm1
          uair_tmp(j,i,iteration) = sqrt(uair_stwave(j,i)**2 + vair_stwave(j,i)**2)
          udir_tmp(j,i,iteration) = mod(270. - atan2d(vair_stwave(j,i),uair_stwave(j,i)), 360.)
          if (udir_tmp(j,i,iteration) .le. 270) then
            udir_tmp(j,i,iteration) = -(90 + udir_tmp(j,i,iteration))
          else
            udir_tmp(j,i,iteration) = 270 - udir_tmp(j,i,iteration)
          end if
          if (mod(n,int(Ti_4_stwave*3600/dt)) == 0) then
            u_stwave(j,i) = sum(uair_tmp(j,i,:)) / size(uair_tmp(j,i,:))
            udir_stwave(j,i) = sum(udir_tmp(j,i,:)) / size(udir_tmp(j,i,:))
          end if 
        end do
      end do

      if (mod(n,int(Ti_4_stwave*3600/dt)) == 0) then
        tau_stwave1 = 0.0
        call stwave(im1,jm1,u_stwave,udir_stwave,tau_stwave1)
        tau_stwave = transpose(tau_stwave1)
      end if
    end if
  end do

END SUBROUTINE stwave_input

! ********************************************************************
SUBROUTINE stwave(im_stw,jm_stw,u_stwave,udir_stwave,tau_stwave1)
! ********************************************************************
!
! Purpose: If suspended sediments is modeled, this subroutine
!          calculates the source and sink terms that depend on 
!          the suspended sediments. 
!
! --------------------------------------------------------------------

  integer, intent(in) :: im_stw
  integer, intent(in) :: jm_stw
  real, intent(in), dimension(jm_stw,im_stw):: u_stwave
  real, intent(in), dimension(jm_stw,im_stw):: udir_stwave
  real, intent(inout), dimension(jm_stw,im_stw):: tau_stwave1
  real, allocatable :: angav_stw(:,:)
  real, allocatable:: angfac_stw(:)
  real, allocatable:: anglz_stw(:)
  real, allocatable::  c_stwave(:,:,:)
  real, allocatable:: cg_stwave(:,:,:)
  real, allocatable:: delf_stw(:)
  real, allocatable:: dep_side_stw(:)
  real, allocatable:: ddepdx_stw(:,:,:)
  real, allocatable:: ddepdy_stw(:,:,:)
  real, allocatable:: delr_stw(:,:)
  real, allocatable:: delr_m_stw(:,:)
  real, allocatable:: delr_p_stw(:,:)
  real, allocatable:: e_stw(:,:,:,:)
  real, allocatable:: e_bc_stw(:,:,:)
  real, allocatable:: e_bc_sea_stw(:,:,:)
  real, allocatable:: etot_stw(:,:)
  real, allocatable:: fm_sea_stw(:,:)
  real, allocatable:: fma_stw(:,:)
  real, allocatable:: freq_stw(:)
  real, allocatable:: H_stw(:,:)
  real, allocatable:: sea_stw(:,:,:,:)
  real, allocatable:: wangle_stw(:,:)
  real, allocatable:: wk_stw(:,:,:)
  real, allocatable:: wt1_stw(:,:)
  real, allocatable:: wt2_stw(:,:)
  real, allocatable:: sxx_stw(:,:)
  real, allocatable:: syy_stw(:,:)
  real, allocatable:: sxy_stw(:,:)
  real, allocatable:: sxxx_stw(:,:)
  real, allocatable:: sxyx_stw(:,:)
  real, allocatable:: sxyy_stw(:,:)
  real, allocatable:: syyy_stw(:,:)
  real, allocatable:: wxrs_stw(:,:)
  real, allocatable:: wyrs_stw(:,:)
  real, allocatable:: dadd_f_stw(:,:)
  real, allocatable:: cf_stw(:,:)
  real, allocatable:: Tmm1_stw(:,:)
  real, allocatable:: taux_stwave(:,:)
  real, allocatable:: tauy_stwave(:,:)
  real, allocatable:: ubottom_stw(:, :)
  real, allocatable:: uairStwave(:,:)
  real, allocatable:: uairdirStwave(:,:)
  integer, allocatable:: i1_stw(:,:)
  integer, allocatable:: i2_stw(:,:)
  integer, allocatable:: iout_stw(:)
  integer, allocatable:: isweep_order_stw(:)
  integer, allocatable:: j1_stw(:,:)
  integer, allocatable:: j2_stw(:,:)
  integer, allocatable:: jout_stw(:)
  integer, allocatable:: l_sweep_stw(:,:)
  integer, allocatable:: inest_stw(:)
  integer, allocatable:: jnest_stw(:)
  integer, allocatable:: ibr_stw(:,:)


  real, dimension(4) :: dep_bc_stw
  real, dimension(4) :: fm_sea_bc_stw
  integer, dimension(4) :: i_st_stw
  integer, dimension(4) :: i_en_stw
  integer, dimension(4) :: i_inc_stw
  integer,  dimension(4) :: j_st_stw
  integer,  dimension(4) :: j_en_stw
  integer,  dimension(4) :: j_inc_stw
  integer,  dimension(4) :: m_sweep_stw
  integer,  dimension(4) :: i_ws_stw
  integer,  dimension(4) :: j_ws_stw

  integer, dimension(4):: i_bc_stw
  integer :: surge_cnt_stw
  integer :: i_side_stw

  real :: epsd_stw
  real :: dth_stw
  real :: epse_stw
  integer :: nfreq_stw
  real :: twopi_stw
  real :: pi2_stw
  real :: radfac_stw
  real :: x0_stw
  real :: y0_stw
  real :: azimuth_stw
  integer :: nest_stw
  real :: dadd_stw
  integer :: m_wind_stw
  integer :: n_wind_stw
  integer :: i_windside_stw


  ! PARAMETERS OF INPUT FILE 
  integer :: iplane_stw
  integer :: iprp_stw
  integer :: icur_stw
  integer :: ibreak_stw
  integer :: irs_stw
  integer :: nselct_stw
  integer :: ibnd_stw
  integer :: iter_sect_stw
  integer :: ifric_stw
  integer :: isurge_stw
  integer :: iwind_stw
  integer :: idep_opt_stw

  integer :: na_4_stw
  integer :: na_stw
  real :: f0_stw
  real :: df_const_stw
  integer :: k_stw

  integer :: ni_stw
  integer :: nj_stw
  integer :: i_stw
  integer :: j_stw
  integer :: l_stw
  integer :: idd_stw
  integer :: i_mid_stw
  integer :: j_mid_stw

  integer :: i_case_stw
  real :: dadd_old_stw

  real :: ang_calc_stw
  real :: angdif_stw

  integer :: m4_stw
  integer :: m3_stw
  integer :: m2_stw
  integer :: m1_stw

  integer :: kmax_stw
  integer :: j_bc3_stw

  real :: sum_stw 
  integer :: i_bc3_stw
  integer :: iter_stw
  real :: xcmp_stw
  real :: ycmp_stw
  real :: emax_stw
  real :: sumfm1_stw
  real :: e1sum_stw
  real :: e_all_stw
  real :: local_wind_stw

! -----------------------------------------------------------------------
  print*,'  Running Stwave'

  ! Initialize constant parameters.
  epsd_stw = 1.0e-2
  epse_stw = 1.0e-8
  twopi_stw = 2.0e0*pi
  pi2_stw = pi/2.0e0
  radfac_stw = pi/180.0e0
  nest_stw = 0
  x0_stw   = 0.0e0      
  y0_stw   = 0.0e0      
  azimuth_stw = 0.0e0  
  dadd_stw = 0.0e0

  ! Set sweep direction defaults for propagation only case.
  m_wind_stw = 1
  n_wind_stw = 1
  i_windside_stw = 1

  ! -----------------------------------------------------------------------
  ! DEFINING THE INPUT PARAMETERS 
  ! FIRST LINE OF _JMS.STD
  iplane_stw = 1
  iprp_stw = 0
  icur_stw = 0
  ibreak_stw = 0
  irs_stw = 0
  nselct_stw = 0
  ibnd_stw = 0
  iter_sect_stw = 3
  ifric_stw = 3
  isurge_stw = 0
  if ((ifsurfbc .le. 3) .or. (ifsurfbc == 20)) then
    iwind_stw = 0
  elseif ((ifsurfbc .ge. 10) .and. (ifsurfbc .lt. 20)) then
    iwind_stw = 1
  end if

  !SECOND LINE _JMS.STD
  idep_opt_stw = 0
  !THIRD LINE _JMS.STD
  i_bc_stw(:) = 0
  !---Set up grid sweep direction information
  !
  ! Boundary Condition types (e_bc_stw(freq_stw,dir,side)):
  !   0 = constant spectrum set equal to zero
  !   1 = constant TMA spectrum with H_stw, T, and dir specified
  !   2 = constant spectrum read from spectral input file (WIS format?)
  !   3 = 1D transformed spectrum (adjacent boundary must be 0, 1, or 2)
  !
  ! Boundary definition:
  !                       ^ y
  !                       |
  !     |  side 4 (n_stw)
  !      -------------
  !     |             |
  !     |             |
  ! side 1 (W)|             |  side 3 (e_stw)
  !     |             |
  !     |             |
  !   0,0  -------------  --> x
  !       side 2 (S)
  !
  ! i_bc_stw  = 0 - constant spectrum set equal to zero
  !     = 1 - constant TMA spectrum with H_stw, T, and dir specified
  !     = 2 - constant spectrum read from a .eng file (file code 8)
  !     = 3 - 1D transformed spectrum (adjacent boundary must be 0, 1, or 2)
  !FOURTH LINE _JMS.STD
  nfreq_stw = 30
  na_stw = 36
  f0_stw = 0.2
  df_const_stw = 0.06
  ! FIFTH LINE OF _JMS.STD
  ! This is the wind idd_stw, the wind speed, wind direction and a dadd_stw
  ! The wind speed and wind direction are assigned below as u_stwave
  ! and udir_stwave. These values come from the surface boundary condition from Si3D
  idd_stw = 1
  dadd_stw = 0

  if (idep_opt_stw .eq. 1) then
    print*, 'STOPPING PROGRAM, NEEDS TO BE IMPLEMENTED'
    stop
  else
    ni_stw = im_stw
    nj_stw = jm_stw
  end if 

  ! set up locations for checking input windsea
  i_mid_stw = ni_stw/2
  j_mid_stw = nj_stw/2
  i_ws_stw(1) = 1
  i_ws_stw(2) = i_mid_stw
  i_ws_stw(3) = ni_stw
  i_ws_stw(4) = i_mid_stw
  j_ws_stw(1) = j_mid_stw
  j_ws_stw(2) = 1
  j_ws_stw(3) = j_mid_stw
  j_ws_stw(4) = nj_stw

  if (i_bc_stw(1) .eq. 2  .or.  i_bc_stw(2) .eq. 2 .or. i_bc_stw(3) .eq. 2 .or. i_bc_stw(4) .eq. 2) then
    ! Read constant spectrum from spectral input file when any 
    ! boundary option is 2.
    print*,'STOPPING PROGRAM, SECTION NEEDS TO BE IMPLEMENTED'
    stop
    if (ibnd_stw.eq.0) then
      ! Coarse offshore grid input is from a Single point on the boundary.
      ! read (8, *)  nfreq_stw, na_stw
    else
      ! Nested nearshore grid input is from:
      ! ibnd_stw=1  a single point on the boundary. No interpolation required.
      ! ibnd_stw=2  user-defined multiple points coinciding with offshore boundary of nearshore grid.
      ! (additional parameters: number of nesting input points and coarse grid cell
      ! location and azimuth_stw in degrees counter clockwise from East)                            
      ! read (8, *)  nfreq_stw, na_stw, nest_in, azimuth_nest                    
    endif
    !allocate (freq_stw(nfreq_stw))
    !read (8, *) (freq_stw(k_stw),k_stw=1,nfreq_stw)
  else
    !FOURTH LINE _JMS.STD
    nfreq_stw = nfreq_stw
    na_stw = na_stw
    f0_stw = f0_stw
    df_const_stw = df_const_stw
    allocate (freq_stw(nfreq_stw))
    do k_stw = 1,nfreq_stw
      freq_stw(k_stw) = f0_stw + (k_stw-1) * df_const_stw
    end do
  endif

  na_4_stw = int(0.25 * na_stw)

  !   Allocate arrays.
  allocate ( angav_stw(nj_stw,ni_stw), angfac_stw(na_stw), anglz_stw(na_stw), c_stwave(nfreq_stw,nj_stw,ni_stw), &
  cg_stwave(nfreq_stw,nj_stw,ni_stw), delf_stw(nfreq_stw),      &
  dep_side_stw(4), ddepdx_stw(na_stw,nj_stw,ni_stw), ddepdy_stw(na_stw,nj_stw,ni_stw),            &
  delr_stw(na_4_stw,4), delr_m_stw(na_4_stw,4), delr_p_stw(na_4_stw,4),               &
  e_stw(nfreq_stw,na_stw,nj_stw,ni_stw), e_bc_stw(nfreq_stw,na_stw,4), e_bc_sea_stw(nfreq_stw,na_stw,4),  &
  etot_stw(nj_stw,ni_stw), fm_sea_stw(nj_stw,ni_stw), fma_stw(nj_stw,ni_stw), H_stw(nj_stw,ni_stw), uairStwave(nj_stw,ni_stw), &
  uairdirStwave(nj_stw,ni_stw), sea_stw(nfreq_stw,na_stw,nj_stw,ni_stw), wangle_stw(nj_stw,ni_stw), wk_stw(nfreq_stw,nj_stw,ni_stw),        &
  wt1_stw(na_4_stw,4), wt2_stw(na_4_stw,4), dadd_f_stw(nj_stw, ni_stw),Tmm1_stw(nj_stw,ni_stw), &
  taux_stwave(nj_stw,ni_stw), tauy_stwave(nj_stw,ni_stw), ubottom_stw(nj_stw,ni_stw))
  allocate (i1_stw(na_4_stw,4), i2_stw(na_4_stw,4), iout_stw(nselct_stw), isweep_order_stw(5), &
   j1_stw(na_4_stw,4), j2_stw(na_4_stw,4), jout_stw(nselct_stw), l_sweep_stw(na_4_stw,4), &
   ibr_stw(nj_stw,ni_stw) )

  if (ibreak_stw .eq. 1) then
    PRINT*, 'STOPPING PROGRAM SECTION NEES TO BE DEVELOPED'
    stop
    ! open (22, file = BreakFile, status = 'unknown')
    ! write (22, *)ni_stw, nj_stw, dx, dy
  endif
  if (irs_stw .eq. 1) then
    PRINT*, 'STOPPING PROGRAM SECTION NEES TO BE DEVELOPED'
    stop
    ! allocate (sxx_stw(nj_stw, ni_stw), syy_stw(nj_stw, ni_stw), sxy_stw(nj_stw, ni_stw), &
    ! sxxx_stw(nj_stw, ni_stw), sxyx_stw(nj_stw, ni_stw), sxyy_stw(nj_stw, ni_stw), &
    ! syyy_stw(nj_stw, ni_stw), wxrs_stw(nj_stw, ni_stw), wyrs_stw(nj_stw, ni_stw))
    ! sxx_stw(1:nj_stw, 1:ni_stw) = 0.0e0
    ! syy_stw(1:nj_stw, 1:ni_stw) = 0.0e0
    ! sxy_stw(1:nj_stw, 1:ni_stw) = 0.0e0
    ! sxxx_stw(1:nj_stw, 1:ni_stw) = 0.0e0
    ! sxyx_stw(1:nj_stw, 1:ni_stw) = 0.0e0
    ! sxyy_stw(1:nj_stw, 1:ni_stw) = 0.0e0
    ! syyy_stw(1:nj_stw, 1:ni_stw) = 0.0e0
    ! wxrs_stw(1:nj_stw, 1:ni_stw) = 0.0e0
    ! wyrs_stw(1:nj_stw, 1:ni_stw) = 0.0e0
    ! open (18, file = RadsFile, status = 'unknown')
    ! write (18, *) ni_stw, nj_stw, dx, dy
    !        open (21, file = 'radstress.21',status = 'unknown')
    !        write (21, *) ni_stw, nj_stw, dx, dy
  endif

  if (ifric_stw .ge. 1) then
    allocate (cf_stw(nj_stw, ni_stw))
    ! open (25, file = FricFile, status = 'old')
    if (ifric_stw .eq. 1 .or. ifric_stw .eq. 3) then
      ! read (25, *)cf_const_stw
      cf_stw(1:nj_stw, 1:ni_stw) = cd
    else
      PRINT*, 'STOPPING PROGRAM SECTION NEES TO BE DEVELOPED'
      stop
      ! read(25,*) nif, njf
      ! if ((nif .ne. ni_stw) .or. (njf .ne. nj_stw)) then
      !   print *,'Friction coefficent field size      ',nif, njf
      !   print *,'does not equal Depth grid dimensions',ni_stw, nj_stw
      !   stop
      ! endif
      ! do j_stw = nj_stw, 1, -1
      !   read (25, *) (cf_stw(j_stw, i_stw), i_stw = 1, ni_stw)
      ! end do
    end if
  endif

  if(isurge_stw .eq. 1) then
    PRINT*, 'STOPPING PROGRAM SECTION NEES TO BE DEVELOPED'
    stop
    ! open(17,file = SurgeFile,  status = 'old')
    ! read (17, *) ni_s, nj_s, dx_s, dy_s  
  endif

  if(iwind_stw .eq. 1)then
    ! PRINT*, 'STOPPING PROGRAM SECTION NEES TO BE DEVELOPED'
    ! stop
    ! open(20,file = WindFile,  status = 'old')
    ! read(20,*) ni_w, nj_w, dx_w, dy_w
    ! print *,ni_w,nj_w,dx_w,dy_w,ni_stw,nj_stw,dx,dy
    ! if(ni_w.ne.ni_stw.or.nj_w.ne.nj_stw)then
    ! stop
    ! endif
  endif

  delf_stw(1) = freq_stw(2) - freq_stw(1)
  delf_stw(nfreq_stw) = freq_stw(nfreq_stw) - freq_stw(nfreq_stw-1)
  do k_stw = 2, nfreq_stw-1
    delf_stw(k_stw) = 0.5 * (freq_stw(k_stw+1) - freq_stw(k_stw-1))
  enddo

  !  Read in special output points.
  if (nselct_stw .ge. 1) then
    PRINT*, 'STOPPING PROGRAM SECTION NEES TO BE DEVELOPED'
    stop
    ! do nn_stw = 1, nselct_stw
    ! read (11, *) iout_stw(nn_stw), jout_stw(nn_stw)
    ! enddo
  endif

  !  read in nesting output points                       
  ! read (11,*,end=10) nest_stw    
  if (nest_stw .gt. 0) then   
    PRINT*, 'STOPPING PROGRAM SECTION NEES TO BE DEVELOPED'
    stop                    
    ! allocate (inest_stw(nest_stw), jnest_stw(nest_stw))               
    ! do nn_stw = 1, nest_stw                                   
    ! read (11, *) inest_stw(nn_stw), jnest_stw(nn_stw)              
    ! enddo    
    ! open (16, file = NestFile, status = 'unknown')    
    ! !       azimuth_stw read from .sim file in Subroutine STWfiles
    ! write (16,*) nfreq_stw, na_stw, nest_stw, azimuth_stw             
    ! write (16,*) (freq_stw(k_stw), k_stw = 1, nfreq_stw)                                            
  endif

  ! Read bathymetry or calculate plane beach

  ! determine average depth on each side
  dep_bc_stw(:) = 0.0e0
  do j_stw=1,nj_stw
    dep_bc_stw(1) = dep_bc_stw(1) + dep_stwave(j_stw,1)
    dep_bc_stw(3) = dep_bc_stw(3) + dep_stwave(j_stw,ni_stw)
  enddo
  dep_bc_stw(1) = dep_bc_stw(1)/nj_stw
  dep_bc_stw(3) = dep_bc_stw(3)/nj_stw
  do i_stw=1,ni_stw
    dep_bc_stw(2) = dep_bc_stw(2) + dep_stwave(1,i_stw)
    dep_bc_stw(4) = dep_bc_stw(4) + dep_stwave(nj_stw,i_stw)
  enddo
  dep_bc_stw(2) = dep_bc_stw(2)/ni_stw
  dep_bc_stw(4) = dep_bc_stw(4)/ni_stw

  dth_stw = 360. / na_stw
  dth_stw = dth_stw * radfac_stw
  do l_stw = 1, na_stw
    anglz_stw(l_stw) = (l_stw-1) * dth_stw
  end do

  ! Set up sweep method for no current and constant depth.
  call sweeps (na_4_stw, na_stw, dth_stw, nj_stw, ni_stw, anglz_stw, l_sweep_stw, i1_stw, i2_stw, j1_stw, j2_stw, wt1_stw, wt2_stw, delr_stw, delr_m_stw, delr_p_stw)
  call setsweeps (i_st_stw, i_en_stw, i_inc_stw, j_st_stw, j_en_stw, j_inc_stw, ni_stw, nj_stw)

  ! Generate or read boundary spectrum for each side.
  ! Do type 1 and 2 boundaries first.
  i_case_stw = 0
  dadd_old_stw = 0.0e0

400 continue

  if(isurge_stw.eq.1)then
    PRINT*, 'STOPPING PROGRAM SECTION NEES TO BE DEVELOPED'
    stop 
    ! read (17,*,end=421)idd_adcirc
    ! do j_stw = nj_stw, 1, -1
    !   read (17,*) (dadd_f_stw(j_stw, i_stw), i_stw = 1, ni_stw)
    ! enddo
  else
    dadd_f_stw(:,:) = dadd_stw
  endif

  ! if(iwind_stw.eq.1)then
  !   ! read (20,*,end=422)idd_wind
  !   ! do j_stw = nj_stw, 1, -1
  !   !   read (20,*) (uairStwave(j_stw,i_stw), i_stw=1,ni_stw)
  !   ! enddo
  !   ! do j_stw = nj_stw, 1, -1
  !   !   read (20,*) (uairdirStwave(j_stw,i_stw), i_stw=1,ni_stw)
  !   ! enddo
  ! else
  uairStwave = u_stwave
  uairdirStwave = udir_stwave
  ! endif

  if((isurge_stw .eq. 1) .and. (iwind_stw .eq. 1)) then 
    ! print *,idd_adcirc,idd_wind
  end if

  ! Initalize arrays (BC = 3 requires no additional action).  
  fma_stw(:,:)       = freq_stw(nfreq_stw)
  fm_sea_stw(:,:)    = freq_stw(nfreq_stw)
  fm_sea_bc_stw(:)   = freq_stw(nfreq_stw)
  e_stw(:,:,:,:)     = 0.0e0
  sea_stw(:,:,:,:)   = 0.0e0
  e_bc_stw(:,:,:)    = 0.0e0
  e_bc_sea_stw(:,:,:)= 0.0e0
  Tmm1_stw(:,:)      = 0.0e0
  
  ! Deterine which side of the grid the wind input is on.
  uairdirStwave(:,:) = uairdirStwave(:,:)*radfac_stw

  do j_stw = 1, nj_stw
    do i_stw = 1, ni_stw
      if (uairdirStwave(j_stw,i_stw) .lt. 0.0) then 
        uairdirStwave(j_stw,i_stw) = uairdirStwave(j_stw,i_stw) + twopi_stw ! resolve negative direction 
      elseif (uairdirStwave(j_stw,i_stw) .gt. twopi_stw) then
        uairdirStwave(j_stw,i_stw) = uairdirStwave(j_stw,i_stw) - twopi_stw
      end if
    enddo
  enddo
  i_windside_stw = (uairdirStwave(j_mid_stw,i_mid_stw) + 0.5*pi2_stw)/pi2_stw + 1.0e0
  if (i_windside_stw .gt. 4) then
    i_windside_stw = 1
  elseif (iwind_stw .eq. 1) then
    i_windside_stw=1
  end if

  ! if there is no wind, set the last sweep to input boundary
  if ( uairStwave(j_mid_stw,i_mid_stw) .eq. 0.0) then
    if ( (i_bc_stw(1) .gt. 0) .and. (i_bc_stw(1) .lt. 3)) then
      i_windside_stw = 1
      m_wind_stw = 1
    else if ((i_bc_stw(2) .gt. 0) .and. (i_bc_stw(2) .lt. 3)) then
      i_windside_stw = 2
      m_wind_stw = 2
    else if ((i_bc_stw(3) .gt. 0) .and. (i_bc_stw(3) .lt. 3)) then
      i_windside_stw = 3
      m_wind_stw = 3
    else
      i_windside_stw = 4
      m_wind_stw = 4
    endif
  endif

  ! Add depth adjustment for this case for all points on the grid.
  do i_stw = 1, ni_stw
    do j_stw = 1, nj_stw
      dep_stwave(j_stw,i_stw) = dep_stwave(j_stw,i_stw) + dadd_f_stw(j_stw,i_stw)
    enddo
  enddo

  call depgrad (nj_stw, ni_stw, na_stw, na_4_stw, l_sweep_stw, i1_stw, i2_stw, j1_stw, j2_stw, anglz_stw,   &
  delr_stw, ddepdx_stw, ddepdy_stw)

  i_case_stw = i_case_stw + 1

  !if (i_case_stw .eq. 1  .or.  dadd_stw .ne. dadd_old_stw) then
  call celerity (nfreq_stw, ni_stw, nj_stw, wk_stw, c_stwave, cg_stwave, freq_stw)
  ! endif
  dadd_old_stw = dadd_stw

  ! Set type 1 and type 2 boundaries (single value on each side, for now).
  ! 1 - constant TMA spectrum with H_stw, T, and dir specified in OptsFile
  ! 2 - constant spectrum read from EngFile (WIS format?)
  do i_side_stw = 1, 4
    if(i_bc_stw(i_side_stw) .gt. 0 .and. i_bc_stw(i_side_stw) .lt. 3)then
      if (i_bc_stw(i_side_stw) .eq. 1) then
        ! read (11, *) H_in, Tp_in, wvang_in
        ! call specgen (dep_bc_stw(i_side_stw) + dadd_stw, e_bc_stw, delf_stw, freq_stw, anglz_stw, &
        ! H_in, Tp_in, wvang_in, i_side_stw)
      else if (i_bc_stw(i_side_stw) .eq. 2) then
        ! if (ibnd_stw.eq.0) then                                  
        !   !       Single point boundary input for a coarse offshore grid.
        !   read (8, *, end = 420) idd_stw
        !   !       read boundary spectra
        !   read (8, *) ( (e_bc_stw(k_stw,l_stw,i_side_stw), l_stw=1,na_stw), k_stw=1,nfreq_stw)
        !   ! If depth location has water, assign boundary array values.
        !   select case (i_side_stw)
        !     case(1)
        !       do j_stw = 1, nj_stw
        !         if (dep_stwave(j_stw,1) .gt. 0.0e0) then
        !           e_stw(:,:,j_stw,1)   = e_bc_stw(:,:,1)
        !         endif
        !       enddo
        !     case(2)  
        !       do i_stw = 1, ni_stw
        !         if (dep_stwave(1,i_stw) .gt. 0.0e0) then
        !           e_stw(:,:,1,i_stw)   = e_bc_stw(:,:,2)
        !         endif
        !       enddo
        !     case(3)
        !       do j_stw = 1, nj_stw
        !         if( dep_stwave(j_stw,ni_stw) .gt. 0.0e0) then
        !           e_stw(:,:,j_stw,ni_stw)   = e_bc_stw(:,:,3)
        !         endif
        !       enddo
        !     case(4)
        !       do i_stw = 1, ni_stw
        !         if (dep_stwave(nj_stw,i_stw) .gt. 0.0e0) then
        !           e_stw(:,:,nj_stw,i_stw)   = e_bc_stw(:,:,4)
        !         endif      
        !       enddo
        !     end select
        ! else
        !   !       ibnd_stw >= 1 when specifying spectra for each offshore cell (j_stw = 1, nj_stw)        
        !   call nest_interp(ni_stw,nj_stw,nest_in,ibnd_stw,azimuth_nest,idd_stw,uairStwave,uairdirStwave,fm_stw,dadd_stw,e_stw)  
        ! endif 
      endif
    endif
  enddo

  ! Determine wind sea_stw component for the wind input side of the grid (i_windside_stw).
  if (iprp_stw .eq. 0) then
    m_wind_stw = int(uairdirStwave(j_mid_stw,i_mid_stw)/pi2_stw + 1.00001)
    !   angle band of wind direction
    n_wind_stw = (uairdirStwave(j_mid_stw,i_mid_stw)-(m_wind_stw-1)*pi2_stw)/dth_stw + 1.00001         
  end if

  do l_stw = 1, na_stw
    ang_calc_stw = cos(uairdirStwave(j_mid_stw,i_mid_stw)) * cos(anglz_stw(l_stw)) + sin(uairdirStwave(j_mid_stw,i_mid_stw)) * sin(anglz_stw(l_stw))
    if (ang_calc_stw .gt. 1.0) ang_calc_stw = 1.0
    if (ang_calc_stw.lt.-1.0) ang_calc_stw = -1.0
    angdif_stw = acos(ang_calc_stw)
    if (abs(angdif_stw) .gt. pi2_stw) angdif_stw = pi2_stw
    angfac_stw(l_stw) = cos(angdif_stw)
    if (angfac_stw(l_stw) .lt. 1.0e-2) angfac_stw(l_stw) = 0.0e0
  end do

  ! for variable winds, propogate all swell first
  if(iwind_stw.eq.1)then
    i_windside_stw = 1
    m_wind_stw = 1
    n_wind_stw = 1
    m4_stw = 1
    m3_stw = 4
    m2_stw = 3
    m1_stw = 2
    ibr_stw(:,:) = 0
    ! Sweeps 1 and 2
    ! 1D BC on upwind end
    if(i_bc_stw(m2_stw) .eq. 3) then
      if(i_bc_stw(m1_stw) .eq. 1 .or. i_bc_stw(m1_stw) .eq. 2) then
        i_stw = i_en_stw(m_wind_stw)+i_inc_stw(m_wind_stw)
        do j_stw = j_st_stw(m_wind_stw), j_en_stw(m_wind_stw), j_inc_stw(m_wind_stw)
          call one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw,delr_stw,delr_m_stw,delr_p_stw,wk_stw, &
              c_stwave,cg_stwave,freq_stw,i_stw,j_stw,i_stw,j_stw-j_inc_stw(m_wind_stw),m1_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
          ! Check for breaking.
          call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
          fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
        end do
      end if
      if(i_bc_stw(m3_stw) .eq. 1 .or. i_bc_stw(m3_stw) .eq. 2) then
        i_stw = i_en_stw(m_wind_stw)+i_inc_stw(m_wind_stw)
        do j_stw = j_en_stw(m_wind_stw), j_st_stw(m_wind_stw), -j_inc_stw(m_wind_stw)
          call one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw,delr_stw,delr_m_stw,delr_p_stw,wk_stw, &
            c_stwave,cg_stwave,freq_stw,i_stw,j_stw,i_stw,j_stw+j_inc_stw(m_wind_stw),m2_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
          ! Check for breaking.
          call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
          fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
        end do
      end if
    end if
    !   Sweep j_stw as inner loop (i_windside_stw = 1 or 3).
    do i_stw = i_en_stw(m_wind_stw), i_st_stw(m_wind_stw), -i_inc_stw(m_wind_stw)
      ! 1D BC on lateral side
      if(i_bc_stw(m1_stw) .eq. 3) then
        j_bc3_stw = j_inc_stw(m_wind_stw)
        j_stw = j_st_stw(m_wind_stw)-j_inc_stw(m_wind_stw)
        call one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw,delr_stw,delr_m_stw,delr_p_stw,wk_stw, &
        c_stwave,cg_stwave,freq_stw,i_stw,j_stw,i_stw+i_inc_stw(m_wind_stw),j_stw,m1_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
        ! Check for breaking.
        call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
        fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
      end if
      ! First sweep
      do j_stw = j_st_stw(m_wind_stw), j_en_stw(m_wind_stw), j_inc_stw(m_wind_stw)
        call prop (nj_stw, ni_stw, na_stw, na_4_stw, anglz_stw, l_sweep_stw, i1_stw, i2_stw, j1_stw, j2_stw, wt1_stw, wt2_stw, e_stw, sea_stw, dth_stw,& 
          delr_stw, delr_m_stw,delr_p_stw, ddepdx_stw, ddepdy_stw, wk_stw, c_stwave, cg_stwave, freq_stw, i_stw, j_stw, m1_stw, epsd_stw, pi2_stw, nfreq_stw, twopi_stw, epse_stw)
        !________________________________________________________________________________
        !   Modified by Patricio Moreno 04/02/2010 (added tau_stwave1 to call friction)
        !________________________________________________________________________________
        if(ifric_stw .ge. 1) then 
          call friction (i_stw, j_stw, k_stw, ni_stw, nj_stw, na_4_stw, e_stw, sea_stw,freq_stw, delf_stw, wk_stw, cg_stwave, cf_stw, ifric_stw, &
            delr_stw, l_sweep_stw, m1_stw, tau_stwave1, ubottom_stw, dth_stw,epsd_stw,nfreq_stw,twopi_stw,pi2_stw,epse_stw,na_stw)
        end if 
        ! Check for breaking.
        call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
        fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
      end do
      j_bc3_stw = 0
      ! 1D BC on lateral side
      if(i_bc_stw(m3_stw) .eq. 3) then
        j_bc3_stw = -j_inc_stw(m_wind_stw)
        j_stw = j_en_stw(m_wind_stw)+j_inc_stw(m_wind_stw)
        call one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw, delr_stw,delr_m_stw,delr_p_stw,wk_stw, & 
          c_stwave,cg_stwave,freq_stw,i_stw,j_stw,i_stw+i_inc_stw(m_wind_stw),j_stw,m2_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
      end if
      ! Sweep 2 (added by Patricio Moreno 07/2010 assuming this is sweep 2)
      do j_stw = j_en_stw(m_wind_stw), j_st_stw(m_wind_stw), -j_inc_stw(m_wind_stw)
        local_wind_stw = int(uairdirStwave(j_stw,i_stw)/pi2_stw + 1.00001)
        call prop (nj_stw, ni_stw, na_stw, na_4_stw, anglz_stw, l_sweep_stw, i1_stw, i2_stw, j1_stw, j2_stw,wt1_stw, wt2_stw, e_stw, sea_stw, dth_stw,&
            delr_stw, delr_m_stw, delr_p_stw, ddepdx_stw, ddepdy_stw, wk_stw, c_stwave, cg_stwave, freq_stw, i_stw, j_stw, m2_stw, epsd_stw, pi2_stw, nfreq_stw, twopi_stw, epse_stw)
        !________________________________________________________________________________
        !   Modified by Patricio Moreno 04/02/2010 (added tau_stwave1 to call friction)
        !________________________________________________________________________________
        if(ifric_stw .ge. 1) then 
          call friction (i_stw, j_stw, k_stw, ni_stw, nj_stw, na_4_stw, e_stw, sea_stw, freq_stw, delf_stw, wk_stw, cg_stwave, cf_stw, ifric_stw,  &
            delr_stw, l_sweep_stw, m2_stw, tau_stwave1, ubottom_stw, dth_stw,epsd_stw,nfreq_stw,twopi_stw,pi2_stw,epse_stw,na_stw)
        end if 
        ! Check for breaking.
        call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
        fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
      end do
    end do
    print *,' Wave ', idd_stw, ';  50% complete iteration', iter_stw
    ! Sweep 3 and 4.
    ! 1D BC on wind side
    if(i_bc_stw(m4_stw) .eq. 3) then
      if(i_bc_stw(m3_stw) .eq. 1 .or. i_bc_stw(m3_stw) .eq. 2) then
        i_stw = i_st_stw(m_wind_stw)-i_inc_stw(m_wind_stw)
        do j_stw = j_en_stw(m_wind_stw), j_st_stw(m_wind_stw), -j_inc_stw(m_wind_stw)
          call one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw,delr_stw,delr_m_stw,delr_p_stw,wk_stw, &
              c_stwave,cg_stwave,freq_stw,i_stw,j_stw,i_stw,j_stw+j_inc_stw(m_wind_stw),m3_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
          ! Check for breaking.
          call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
          fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
        end do
      end if
      if(i_bc_stw(m1_stw) .eq. 1 .or. i_bc_stw(m1_stw) .eq. 3) then
        i_stw = i_st_stw(m_wind_stw)-i_inc_stw(m_wind_stw)
        do j_stw = j_st_stw(m_wind_stw), j_en_stw(m_wind_stw), j_inc_stw(m_wind_stw)
          call one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw,delr_stw,delr_m_stw,delr_p_stw,wk_stw, &
              c_stwave,cg_stwave,freq_stw,i_stw,j_stw,i_stw,j_stw-j_inc_stw(m_wind_stw),m4_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
          ! Check for breaking.
          call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
          fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
        end do
      end if
    end if
    do i_stw = i_st_stw(m_wind_stw), i_en_stw(m_wind_stw), i_inc_stw(m_wind_stw)
      j_bc3_stw = 0
      ! 1D BC on lateral side
      if(i_bc_stw(m3_stw) .eq. 3) then
        j_bc3_stw = -j_inc_stw(m_wind_stw)
        j_stw = j_en_stw(m_wind_stw)+j_inc_stw(m_wind_stw)
        call one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw,delr_stw,delr_m_stw,delr_p_stw,wk_stw, & 
            c_stwave,cg_stwave,freq_stw,i_stw,j_stw,i_stw-i_inc_stw(m_wind_stw),j_stw,m3_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
        ! Check for breaking.
        call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
        fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
      end if
      ! Sweep 3
      do j_stw = j_en_stw(m_wind_stw) ,j_st_stw(m_wind_stw), -j_inc_stw(m_wind_stw)
        call prop (nj_stw, ni_stw, na_stw, na_4_stw, anglz_stw, l_sweep_stw, i1_stw, i2_stw, j1_stw, j2_stw,wt1_stw, wt2_stw, e_stw, sea_stw, dth_stw,& 
            delr_stw, delr_m_stw, delr_p_stw, ddepdx_stw, ddepdy_stw, wk_stw, c_stwave, cg_stwave, freq_stw, i_stw, j_stw, m3_stw, epsd_stw, pi2_stw, nfreq_stw, twopi_stw, epse_stw)
        !_______________________________________________________________________________
        !   Modified by Patricio Moreno 04/02/2010 (added tau_stwave1 to call friction)
        !________________________________________________________________________________
        if(ifric_stw.ge.1)call friction (i_stw, j_stw, k_stw, ni_stw, nj_stw, na_4_stw, e_stw, sea_stw,freq_stw, delf_stw, wk_stw, &
          cg_stwave, cf_stw, ifric_stw, delr_stw, l_sweep_stw, m3_stw, tau_stwave1, ubottom_stw, dth_stw,epsd_stw,nfreq_stw,twopi_stw,pi2_stw,epse_stw,na_stw)
        ! Check for breaking.
        call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
        fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
      end do
      j_bc3_stw = 0
      if(i_bc_stw(m1_stw) .eq. 3) then
        ! 1D BC on lateral side
        j_bc3_stw = j_inc_stw(m_wind_stw)
        j_stw  = j_st_stw(m_wind_stw)-j_inc_stw(m_wind_stw)
        call one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw,delr_stw,delr_m_stw,delr_p_stw,wk_stw, &
            c_stwave,cg_stwave,freq_stw,i_stw,j_stw,i_stw-i_inc_stw(m_wind_stw),j_stw,m4_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
        ! Check for breaking.
        call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
        fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
      end if
      ! Sweep 4
      do j_stw = j_st_stw(m_wind_stw), j_en_stw(m_wind_stw), j_inc_stw(m_wind_stw)
        call prop (nj_stw, ni_stw, na_stw, na_4_stw, anglz_stw, l_sweep_stw, i1_stw, i2_stw, j1_stw, j2_stw,wt1_stw, wt2_stw, e_stw, sea_stw, dth_stw,&
          delr_stw, delr_m_stw, delr_p_stw, ddepdx_stw, ddepdy_stw, wk_stw, c_stwave, cg_stwave, freq_stw, i_stw, j_stw, m4_stw, epsd_stw, pi2_stw, nfreq_stw, twopi_stw, epse_stw)
        !________________________________________________________________________________
        !   Modified by Patricio Moreno 04/02/2010 (added tau_stwave1 to call friction)
        !________________________________________________________________________________
        if(ifric_stw .ge. 1) then 
          call friction (i_stw, j_stw, k_stw, ni_stw, nj_stw, na_4_stw, e_stw, sea_stw,freq_stw, delf_stw, wk_stw, cg_stwave, cf_stw, &
            ifric_stw, delr_stw, l_sweep_stw, m4_stw, tau_stwave1, ubottom_stw, dth_stw,epsd_stw,nfreq_stw,twopi_stw,pi2_stw,epse_stw,na_stw)
        end if 
        ! Check for breaking.
        call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
        fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
      end do
    end do 
  print *,' Wave ', idd_stw, ';  Complete prop only step for spatially variable winds'
  endif

  !-----------------------------------------------------------------------------------
  ! Proprogate waves in 4 sweeps.  If there is no wind input, sweeps sides 1 - 4.  
  ! If there IS wind input (iprp_stw=0), then sweep the wind side last and calculate
  ! generation during the last sweep.  
  ! m_wind_stw is sweep order of grid sides including wind direction
  ! n_wind_stw is the angle band of wind direction
  ! i_windside_stw is side with wind input
  !
  ! Determine sweep order, where m_stw is m1_stw is first sweep and m4_stw is last sweep, 
  ! value indicate quadrant being swept

  m4_stw = m_wind_stw
  if (i_windside_stw .eq. m_wind_stw) then
    m3_stw = m_wind_stw - 1
    if (m3_stw .lt. 1) then
      m3_stw = m3_stw + 4
    end if
    m2_stw = m3_stw - 1
    if (m2_stw .lt. 1) then 
      m2_stw = m2_stw + 4
    end if 
    m1_stw = m2_stw - 1
    if(m1_stw .lt. 1) then 
      m1_stw = m1_stw + 4
    end if
  else
    m3_stw = m_wind_stw + 1
    if (m3_stw .gt. 4) then
      m3_stw = m3_stw - 4
    end if 
    m2_stw = m3_stw + 1
    if (m2_stw .gt. 4) then 
      m2_stw = m2_stw - 4
    end if 
    m1_stw = m2_stw + 1
    if (m1_stw .gt. 4) then 
      m1_stw = m1_stw - 4
    end if 
  endif

  !---Iterations.
  do iter_stw = 1 ,iter_sect_stw
    ibr_stw(:,:) = 0
    ! Sweeps 1 and 2  
    if(i_windside_stw .eq. 1  .or.  i_windside_stw.eq.3) then
      ! 1D BC on upwind end
      if(i_bc_stw(m2_stw) .eq. 3 .and. iter_stw .eq. 1) then
        if(i_bc_stw(m1_stw) .eq. 1 .or. i_bc_stw(m1_stw) .eq. 2) then
          i_stw = i_en_stw(m_wind_stw)+i_inc_stw(m_wind_stw)
          do j_stw = j_st_stw(m_wind_stw), j_en_stw(m_wind_stw), j_inc_stw(m_wind_stw)
            call one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw,delr_stw,delr_m_stw,delr_p_stw,wk_stw, &
                c_stwave,cg_stwave,freq_stw,i_stw,j_stw,i_stw,j_stw-j_inc_stw(m_wind_stw),m1_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
            ! Check for breaking.
            call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
            fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
          end do
        end if
        if(i_bc_stw(m3_stw) .eq. 1 .or. i_bc_stw(m3_stw) .eq. 2) then
          i_stw = i_en_stw(m_wind_stw)+i_inc_stw(m_wind_stw)
          do j_stw = j_en_stw(m_wind_stw), j_st_stw(m_wind_stw), -j_inc_stw(m_wind_stw)
            call one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw,delr_stw,delr_m_stw,delr_p_stw,wk_stw, &
                c_stwave,cg_stwave,freq_stw,i_stw,j_stw,i_stw,j_stw+j_inc_stw(m_wind_stw),m2_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
            ! Check for breaking.
            call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
            fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
          end do
        end if
      end if
      !   Sweep j_stw as inner loop (i_windside_stw = 1 or 3).
      do i_stw = i_en_stw(m_wind_stw), i_st_stw(m_wind_stw), -i_inc_stw(m_wind_stw)
        ! 1D BC on lateral side
        if(i_bc_stw(m1_stw) .eq. 3 .and. iter_stw .eq. 1) then
          j_bc3_stw = j_inc_stw(m_wind_stw)
          j_stw = j_st_stw(m_wind_stw)-j_inc_stw(m_wind_stw)
          call one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw,delr_stw,delr_m_stw,delr_p_stw,wk_stw, &
            c_stwave,cg_stwave,freq_stw,i_stw,j_stw,i_stw+i_inc_stw(m_wind_stw),j_stw,m1_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
          ! Check for breaking.
          call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
          fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
        end if
        ! First sweep
        do j_stw = j_st_stw(m_wind_stw), j_en_stw(m_wind_stw), j_inc_stw(m_wind_stw)
          call prop (nj_stw, ni_stw, na_stw, na_4_stw, anglz_stw, l_sweep_stw, i1_stw, i2_stw, j1_stw, j2_stw, wt1_stw, wt2_stw, e_stw, sea_stw, dth_stw,& 
            delr_stw, delr_m_stw,delr_p_stw,ddepdx_stw, ddepdy_stw, wk_stw, c_stwave, cg_stwave, freq_stw, i_stw, j_stw, m1_stw, epsd_stw, pi2_stw, nfreq_stw, twopi_stw, epse_stw)
          if (iprp_stw .eq. 0) then
            local_wind_stw = int(uairdirStwave(j_stw,i_stw)/pi2_stw + 1.00001)
            ! Solve energy source terms in arbitrary depth water.
            if(local_wind_stw .eq. m1_stw) then
              call gen (nj_stw, ni_stw, na_4_stw, anglz_stw, e_stw, sea_stw, delf_stw, wk_stw, c_stwave, cg_stwave, &
                freq_stw, etot_stw, angav_stw, fma_stw, fm_sea_stw,i_stw, j_stw, m_wind_stw, n_wind_stw, wt1_stw, wt2_stw, i1_stw, &
                i2_stw, j1_stw, j2_stw, idd_stw, uairStwave, uairdirStwave ,na_stw, epsd_stw, nfreq_stw, dth_stw,pi2_stw, twopi_stw, epse_stw,angfac_stw)  
            end if
          endif
          !________________________________________________________________________________
          !   Modified by Patricio Moreno 04/02/2010 (added tau_stwave1 to call friction)
          !________________________________________________________________________________
          if(ifric_stw .ge. 1) then
            call friction (i_stw, j_stw, k_stw, ni_stw, nj_stw, na_4_stw, e_stw, sea_stw,freq_stw, delf_stw, wk_stw, cg_stwave, &
              cf_stw, ifric_stw, delr_stw, l_sweep_stw, m1_stw, tau_stwave1, ubottom_stw, dth_stw,epsd_stw,nfreq_stw,twopi_stw,pi2_stw,epse_stw,na_stw)
          end if

          ! Check for breaking.
          call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
          fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
        end do
         
        j_bc3_stw = 0
        ! 1D BC on lateral side
        if(i_bc_stw(m3_stw) .eq. 3) then
          j_bc3_stw = -j_inc_stw(m_wind_stw)
          j_stw = j_en_stw(m_wind_stw)+j_inc_stw(m_wind_stw)
          call one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw,delr_stw,delr_m_stw,delr_p_stw,wk_stw, &
            c_stwave,cg_stwave,freq_stw,i_stw,j_stw,i_stw+i_inc_stw(m_wind_stw),j_stw,m2_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
        endif
        do j_stw = j_en_stw(m_wind_stw), j_st_stw(m_wind_stw), -j_inc_stw(m_wind_stw)
          call prop (nj_stw, ni_stw, na_stw, na_4_stw, anglz_stw, l_sweep_stw, i1_stw, i2_stw, j1_stw, j2_stw, wt1_stw, wt2_stw, e_stw, sea_stw, dth_stw,&
            delr_stw, delr_m_stw, delr_p_stw, ddepdx_stw, ddepdy_stw, wk_stw, c_stwave, cg_stwave, freq_stw, i_stw, j_stw, m2_stw, epsd_stw, pi2_stw, nfreq_stw, twopi_stw, epse_stw)
          if (iprp_stw .eq. 0) then
            local_wind_stw = int(uairdirStwave(j_stw,i_stw)/pi2_stw + 1.00001)
            !Solve energy source terms in arbitrary depth water.
            if(local_wind_stw .eq. m2_stw) then
              call gen (nj_stw, ni_stw, na_4_stw, anglz_stw, e_stw, sea_stw, delf_stw, wk_stw, c_stwave, cg_stwave, & 
                freq_stw, etot_stw, angav_stw, fma_stw, fm_sea_stw, i_stw, j_stw, m_wind_stw, n_wind_stw, wt1_stw, wt2_stw, i1_stw, &
                i2_stw, j1_stw, j2_stw, idd_stw, uairStwave, uairdirStwave,na_stw, epsd_stw, nfreq_stw, dth_stw,pi2_stw, twopi_stw, epse_stw,angfac_stw)
            end if
          end if
          !________________________________________________________________________________
          !   Modified by Patricio Moreno 04/02/2010 (added tau_stwave1 to call friction)
          !________________________________________________________________________________
          if(ifric_stw .ge. 1) then
            call friction (i_stw, j_stw, k_stw, ni_stw, nj_stw, na_4_stw, e_stw, sea_stw, freq_stw, delf_stw, wk_stw, &
              cg_stwave, cf_stw, ifric_stw, delr_stw, l_sweep_stw, m2_stw, tau_stwave1, ubottom_stw, dth_stw,epsd_stw,nfreq_stw,twopi_stw,pi2_stw,epse_stw,na_stw)
          end if 
          !     Check for breaking.
          call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
          fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
        end do
      end do

      print *,' Wave ', idd_stw, ';  50% complete iteration', iter_stw
      ! Sweep 3 and 4.
      ! 1D BC on wind side
      if(i_bc_stw(m4_stw) .eq. 3) then
        if(i_bc_stw(m3_stw) .eq. 1 .or. i_bc_stw(m3_stw) .eq. 2) then
          i_stw = i_st_stw(m_wind_stw)-i_inc_stw(m_wind_stw)
          do j_stw = j_en_stw(m_wind_stw), j_st_stw(m_wind_stw), -j_inc_stw(m_wind_stw)
            call one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw,delr_stw,delr_m_stw,delr_p_stw,wk_stw, &
              c_stwave,cg_stwave,freq_stw,i_stw,j_stw, i_stw,j_stw+j_inc_stw(m_wind_stw),m3_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
            ! Check for breaking.
            call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
            fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
          end do
        end if
        if(i_bc_stw(m1_stw) .eq. 1 .or. i_bc_stw(m1_stw) .eq. 3) then
          i_stw = i_st_stw(m_wind_stw)-i_inc_stw(m_wind_stw)
          do j_stw = j_st_stw(m_wind_stw), j_en_stw(m_wind_stw), j_inc_stw(m_wind_stw)
            call one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw,delr_stw,delr_m_stw,delr_p_stw,wk_stw, &
              c_stwave,cg_stwave,freq_stw,i_stw,j_stw,i_stw,j_stw-j_inc_stw(m_wind_stw),m4_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
            ! Check for breaking.
            call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
            fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
          end do
        end if
      end if

      do i_stw = i_st_stw(m_wind_stw), i_en_stw(m_wind_stw), i_inc_stw(m_wind_stw)
        j_bc3_stw = 0
        ! 1D BC on lateral side
        if(i_bc_stw(m3_stw).eq.3)then
          j_bc3_stw = -j_inc_stw(m_wind_stw)
          j_stw = j_en_stw(m_wind_stw)+j_inc_stw(m_wind_stw)
          call one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw,delr_stw,delr_m_stw,delr_p_stw,wk_stw, &
            c_stwave,cg_stwave,freq_stw,i_stw,j_stw,i_stw-i_inc_stw(m_wind_stw),j_stw,m3_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
            !     Check for breaking.
            call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
            fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
        end if
        ! Sweep 3
        do j_stw = j_en_stw(m_wind_stw) ,j_st_stw(m_wind_stw), -j_inc_stw(m_wind_stw)
          call prop (nj_stw, ni_stw, na_stw, na_4_stw, anglz_stw, l_sweep_stw, i1_stw, i2_stw, j1_stw, j2_stw, wt1_stw, wt2_stw, e_stw, sea_stw, dth_stw,&
            delr_stw, delr_m_stw, delr_p_stw, ddepdx_stw, ddepdy_stw, wk_stw, c_stwave, cg_stwave, freq_stw, i_stw, j_stw, m3_stw, epsd_stw, pi2_stw, nfreq_stw, twopi_stw, epse_stw)
          if (iprp_stw .eq. 0) then
            local_wind_stw = int(uairdirStwave(j_stw,i_stw)/pi2_stw + 1.00001)
            ! Solve energy source terms in arbitrary depth water.
            if(local_wind_stw .eq. m3_stw) then
              call gen (nj_stw, ni_stw, na_4_stw, anglz_stw, e_stw, sea_stw, delf_stw, wk_stw, c_stwave, cg_stwave, freq_stw, & 
                etot_stw, angav_stw, fma_stw, fm_sea_stw, i_stw, j_stw, m_wind_stw, n_wind_stw, wt1_stw, wt2_stw, i1_stw,i2_stw, j1_stw, j2_stw, idd_stw, uairStwave, uairdirStwave,na_stw, epsd_stw, nfreq_stw, dth_stw,pi2_stw, twopi_stw, epse_stw,angfac_stw)
            end if
          end if
          !________________________________________________________________________________
          !   Modified by Patricio Moreno 04/02/2010 (added tau_stwave1 to call friction)
          !________________________________________________________________________________
          if(ifric_stw.ge.1) then
            call friction (i_stw, j_stw, k_stw, ni_stw, nj_stw, na_4_stw, e_stw, sea_stw, freq_stw, delf_stw, wk_stw, cg_stwave, cf_stw, ifric_stw, &
              delr_stw, l_sweep_stw, m3_stw, tau_stwave1, ubottom_stw, dth_stw,epsd_stw,nfreq_stw,twopi_stw,pi2_stw,epse_stw,na_stw)
          end if 
          ! Check for breaking.
          call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
          fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
        end do
        j_bc3_stw = 0 
        if(i_bc_stw(m1_stw) .eq. 3) then
          ! 1D BC on lateral side
          j_bc3_stw = j_inc_stw(m_wind_stw)
          j_stw  = j_st_stw(m_wind_stw)-j_inc_stw(m_wind_stw)
          call one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw,delr_stw,delr_m_stw,delr_p_stw,wk_stw, & 
            c_stwave,cg_stwave,freq_stw,i_stw,j_stw,i_stw-i_inc_stw(m_wind_stw),j_stw,m4_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
          ! Check for breaking.
          call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
          fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
        end if
        ! Sweep 4      
        do j_stw = j_st_stw(m_wind_stw), j_en_stw(m_wind_stw), j_inc_stw(m_wind_stw)
          call prop (nj_stw, ni_stw, na_stw, na_4_stw, anglz_stw, l_sweep_stw, i1_stw, i2_stw, j1_stw, j2_stw, wt1_stw, wt2_stw, e_stw, sea_stw, dth_stw,& 
            delr_stw, delr_m_stw, delr_p_stw,ddepdx_stw, ddepdy_stw, wk_stw, c_stwave, cg_stwave, freq_stw, i_stw, j_stw, m4_stw, epsd_stw, pi2_stw, nfreq_stw, twopi_stw, epse_stw)
          if (iprp_stw .eq. 0) then
            local_wind_stw = int(uairdirStwave(j_stw,i_stw)/pi2_stw + 1.00001)
            ! Solve energy source terms in arbitrary depth water.
            if(local_wind_stw .eq. m4_stw) then 
              call gen (nj_stw, ni_stw, na_4_stw, anglz_stw, e_stw, sea_stw, delf_stw, wk_stw, c_stwave, cg_stwave, freq_stw, &
                etot_stw, angav_stw, fma_stw, fm_sea_stw,i_stw, j_stw, m_wind_stw, n_wind_stw, wt1_stw, wt2_stw, i1_stw,i2_stw, j1_stw, j2_stw, idd_stw, uairStwave, uairdirStwave,na_stw, epsd_stw, nfreq_stw, dth_stw,pi2_stw, twopi_stw, epse_stw,angfac_stw)  
            end if 
          end if
          !________________________________________________________________________________
          !   Modified by Patricio Moreno 04/02/2010 (added tau_stwave1 to call friction)
          !________________________________________________________________________________
          if(ifric_stw .ge. 1) then
            call friction (i_stw, j_stw, k_stw, ni_stw, nj_stw, na_4_stw, e_stw, sea_stw,freq_stw, delf_stw, wk_stw, cg_stwave, cf_stw, ifric_stw, &
            delr_stw, l_sweep_stw, m4_stw, tau_stwave1, ubottom_stw, dth_stw,epsd_stw,nfreq_stw,twopi_stw,pi2_stw,epse_stw,na_stw)
          end if
          !     Check for breaking.
          call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
          fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
        end do
      end do
      print *,' Wave ', idd_stw, ';  Complete iteration', iter_stw

    else
      ! Sweep i_stw as inner loops (i_windside_stw = 2 or 4).
      ! 1D BC on upwind end
      if(i_bc_stw(m2_stw) .eq. 3) then
        if(i_bc_stw(m1_stw) .eq. 1 .or. i_bc_stw(m1_stw) .eq. 2) then
          j_stw = j_en_stw(m_wind_stw) + j_inc_stw(m_wind_stw)
          do i_stw = i_st_stw(m_wind_stw), i_en_stw(m_wind_stw), i_inc_stw(m_wind_stw)
            call one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw,delr_stw,delr_m_stw,delr_p_stw,wk_stw, &
              c_stwave,cg_stwave,freq_stw,i_stw,j_stw,i_stw-i_inc_stw(m_wind_stw),j_stw,m1_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
            ! Check for breaking.
            call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
            fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
          end do
        end if
        if(i_bc_stw(m3_stw) .eq. 1 .or. i_bc_stw(m3_stw) .eq. 3) then
          j_stw = j_en_stw(m_wind_stw) + j_inc_stw(m_wind_stw)
          do i_stw = i_en_stw(m_wind_stw), i_st_stw(m_wind_stw), -i_inc_stw(m_wind_stw)
            call one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw,delr_stw,delr_m_stw,delr_p_stw,wk_stw, &
              c_stwave,cg_stwave,freq_stw,i_stw,j_stw,i_stw+i_inc_stw(m_wind_stw),j_stw,m2_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
            ! Check for breaking.
            call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
            fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
          end do
        end if
      end if
      ! Sweep j_stw as outloop and i_stw as inner loop

      do j_stw = j_en_stw(m_wind_stw), j_st_stw(m_wind_stw), - j_inc_stw(m_wind_stw)
        i_bc3_stw = 0
        ! 1D BC on lateral side
        if(i_bc_stw(m1_stw).eq.3)then 
          ! i_bc3_stw = i_inc_stw(m_wind_stw)
          i_stw = i_st_stw(m_wind_stw)-i_inc_stw(m_wind_stw)
          call one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw,delr_stw,delr_m_stw,delr_p_stw,wk_stw, &
            c_stwave,cg_stwave,freq_stw,i_stw,j_stw,i_stw,j_stw+j_inc_stw(m_wind_stw),m1_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
          ! Check for breaking.
          call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
          fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
        end if
        ! First sweep
        do i_stw = i_st_stw(m_wind_stw), i_en_stw(m_wind_stw)+i_bc3_stw, i_inc_stw(m_wind_stw)
          call prop(nj_stw, ni_stw, na_stw, na_4_stw, anglz_stw, l_sweep_stw, i1_stw, i2_stw, j1_stw, j2_stw,wt1_stw, wt2_stw, e_stw, sea_stw, dth_stw,&
            delr_stw, delr_m_stw, delr_p_stw,ddepdx_stw, ddepdy_stw, wk_stw, c_stwave, cg_stwave, freq_stw, i_stw, j_stw, m1_stw, epsd_stw, pi2_stw, nfreq_stw, twopi_stw, epse_stw)
          if (iprp_stw .eq. 0) then
            local_wind_stw = int(uairdirStwave(j_stw,i_stw)/pi2_stw + 1.00001)
            ! Solve energy source terms in arbitrary depth water.
            if(local_wind_stw .eq. m1_stw) then
              call gen (nj_stw, ni_stw, na_4_stw, anglz_stw, e_stw, sea_stw, delf_stw, wk_stw, c_stwave, cg_stwave, &
                freq_stw, etot_stw, angav_stw, fma_stw, fm_sea_stw, i_stw, j_stw, m_wind_stw, n_wind_stw, wt1_stw, wt2_stw, i1_stw,   &
                i2_stw, j1_stw, j2_stw, idd_stw, uairStwave, uairdirStwave,na_stw, epsd_stw, nfreq_stw, dth_stw,pi2_stw, twopi_stw, epse_stw,angfac_stw)  
            end if
          end if
          !________________________________________________________________________________
          !   Modified by Patricio Moreno 04/02/2010 (added tau_stwave1 to call friction)
          !________________________________________________________________________________
          if(ifric_stw .ge. 1) then 
            call friction (i_stw, j_stw, k_stw, ni_stw, nj_stw, na_4_stw, e_stw, sea_stw, freq_stw, delf_stw, wk_stw, &
              cg_stwave, cf_stw, ifric_stw, delr_stw, l_sweep_stw, m1_stw, tau_stwave1, ubottom_stw, dth_stw,epsd_stw,nfreq_stw,twopi_stw,pi2_stw,epse_stw,na_stw)
          end if 
          ! Check for breaking.
          call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
          fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
        end do
        j_bc3_stw = 0
        ! 1D BC on lateral side
        if(i_bc_stw(m3_stw) .eq. 3) then
          ! i_bc3_stw = -i_inc_stw(m_wind_stw)
          i_stw = i_en_stw(m_wind_stw)+i_inc_stw(m_wind_stw)
          call one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw,delr_stw,delr_m_stw,delr_p_stw,wk_stw, &
            c_stwave,cg_stwave,freq_stw,i_stw,j_stw,i_stw,j_stw+j_inc_stw(m_wind_stw),m2_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
          ! Check for breaking.
          call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
          fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
        endif
        ! Second Sweep
        do i_stw = i_en_stw(m_wind_stw), i_st_stw(m_wind_stw)+i_bc3_stw, - i_inc_stw(m_wind_stw)
          call prop (nj_stw, ni_stw, na_stw, na_4_stw, anglz_stw, l_sweep_stw, i1_stw, i2_stw, j1_stw, j2_stw, wt1_stw, wt2_stw, e_stw, sea_stw, dth_stw,&
            delr_stw, delr_m_stw, delr_p_stw, ddepdx_stw, ddepdy_stw, wk_stw, c_stwave, cg_stwave, freq_stw, i_stw, j_stw, m2_stw, epsd_stw, pi2_stw, nfreq_stw, twopi_stw, epse_stw)
          if (iprp_stw .eq. 0) then
            local_wind_stw = int(uairdirStwave(j_stw,i_stw)/pi2_stw + 1.00001)
            ! Solve energy source terms in arbitrary depth water.
            if(local_wind_stw .eq. m2_stw) then 
              call gen (nj_stw, ni_stw, na_4_stw, anglz_stw, e_stw, sea_stw, delf_stw, wk_stw, c_stwave, cg_stwave, & 
                freq_stw, etot_stw, angav_stw, fma_stw, fm_sea_stw,i_stw, j_stw, m_wind_stw, n_wind_stw, wt1_stw, wt2_stw, i1_stw,   &
                i2_stw, j1_stw, j2_stw, idd_stw, uairStwave, uairdirStwave,na_stw, epsd_stw, nfreq_stw, dth_stw,pi2_stw, twopi_stw, epse_stw,angfac_stw)  
            end if
          end if
          !________________________________________________________________________________
          !   Modified by Patricio Moreno 04/02/2010 (added tau_stwave1 to call friction)
          !________________________________________________________________________________
          if(ifric_stw.ge.1) then 
            call friction (i_stw, j_stw, k_stw, ni_stw, nj_stw, na_4_stw, e_stw, sea_stw, freq_stw, delf_stw, wk_stw, &
              cg_stwave, cf_stw, ifric_stw, delr_stw, l_sweep_stw, m2_stw, tau_stwave1, ubottom_stw, dth_stw,epsd_stw,nfreq_stw,twopi_stw,pi2_stw,epse_stw,na_stw)
          end if 
          !     Check for breaking.
          call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
          fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
        end do
      end do
      print *,' Wave ', idd_stw, ';  50% complete iteration', iter_stw
      ! Sweep 3 and 4
      ! 1D BC on wind side
      if(i_bc_stw(m4_stw) .eq. 3) then
        if(i_bc_stw(m3_stw) .eq. 1 .or. i_bc_stw(m3_stw) .eq. 2) then
          j_stw = j_st_stw(m_wind_stw)-j_inc_stw(m_wind_stw)
          do i_stw = i_en_stw(m_wind_stw), i_st_stw(m_wind_stw), -i_inc_stw(m_wind_stw)
            call one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw,delr_stw,delr_m_stw,delr_p_stw,wk_stw, &
              c_stwave,cg_stwave,freq_stw,i_stw,j_stw, i_stw+i_inc_stw(m_wind_stw),j_stw,m3_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
            ! Check for breaking.
            call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
            fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
          end do
        end if 
        if(i_bc_stw(m1_stw).eq.1.or.i_bc_stw(m1_stw).eq.3)then
          j_stw = j_st_stw(m_wind_stw)-j_inc_stw(m_wind_stw)
          do i_stw = i_st_stw(m_wind_stw), i_en_stw(m_wind_stw), i_inc_stw(m_wind_stw)
            call one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw,delr_stw,delr_m_stw,delr_p_stw,wk_stw, &
              c_stwave,cg_stwave,freq_stw,i_stw,j_stw,i_stw-i_inc_stw(m_wind_stw),j_stw,m4_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
            !     Check for breaking.
            call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
            fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
          end do
        end if
      end if
      ! Sweep 3
      do j_stw = j_st_stw(m_wind_stw), j_en_stw(m_wind_stw), j_inc_stw(m_wind_stw)
        i_bc3_stw = 0
        ! 1D BC on lateral side
        if(i_bc_stw(m3_stw) .eq. 3) then
          ! i_bc3_stw = -i_inc_stw(m_wind_stw)
          i_stw = i_en_stw(m_wind_stw)+i_inc_stw(m_wind_stw)
          call one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw,delr_stw,delr_m_stw,delr_p_stw,wk_stw, &
            c_stwave,cg_stwave,freq_stw,i_stw,j_stw,i_stw,j_stw-j_inc_stw(m_wind_stw),m3_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
          ! Check for breaking.
          call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
          fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
        end if
        do i_stw = i_en_stw(m_wind_stw), i_st_stw(m_wind_stw)+i_bc3_stw, - i_inc_stw(m_wind_stw)
          call prop (nj_stw, ni_stw, na_stw, na_4_stw, anglz_stw, l_sweep_stw, i1_stw, i2_stw, j1_stw, j2_stw,wt1_stw, wt2_stw, e_stw, sea_stw, dth_stw,&
            delr_stw, delr_m_stw, delr_p_stw,ddepdx_stw, ddepdy_stw, wk_stw, c_stwave, cg_stwave, freq_stw, i_stw, j_stw, m3_stw, epsd_stw, pi2_stw, nfreq_stw, twopi_stw, epse_stw)
          if (iprp_stw .eq. 0) then
            local_wind_stw = int(uairdirStwave(j_stw,i_stw)/pi2_stw + 1.00001)
            ! Solve energy source terms in arbitrary depth water.
            if(local_wind_stw .eq. m3_stw) then 
              call gen (nj_stw, ni_stw, na_4_stw, anglz_stw, e_stw, sea_stw, delf_stw, wk_stw, c_stwave, cg_stwave, &
                freq_stw, etot_stw, angav_stw, fma_stw, fm_sea_stw,i_stw, j_stw, m_wind_stw, n_wind_stw, wt1_stw, wt2_stw, i1_stw,   &
                i2_stw, j1_stw, j2_stw, idd_stw, uairStwave, uairdirStwave,na_stw, epsd_stw, nfreq_stw, dth_stw,pi2_stw, twopi_stw, epse_stw,angfac_stw)
            end if 
          end if
          !________________________________________________________________________________
          !   Modified by Patricio Moreno 04/02/2010 (added tau_stwave1 to call friction)
          !________________________________________________________________________________
          if(ifric_stw.ge.1) then
            call friction (i_stw, j_stw, k_stw, ni_stw, nj_stw, na_4_stw, e_stw, sea_stw,freq_stw, delf_stw, wk_stw, cg_stwave, cf_stw, ifric_stw, &
              delr_stw, l_sweep_stw, m3_stw, tau_stwave1, ubottom_stw, dth_stw,epsd_stw,nfreq_stw,twopi_stw,pi2_stw,epse_stw,na_stw)
          end if 
          ! Check for breaking.
          call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
          fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
        end do

        j_bc3_stw = 0 
        if(i_bc_stw(m1_stw) .eq. 3) then
          ! 1D BC on lateral side
          !     i_bc3_stw = i_inc_stw(m_wind_stw)
          i_stw = i_st_stw(m_wind_stw)-i_inc_stw(m_wind_stw)
          call one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw,delr_stw,delr_m_stw,delr_p_stw,wk_stw, &
            c_stwave,cg_stwave,freq_stw,i_stw,j_stw, i_stw,j_stw-j_inc_stw(m_wind_stw),m4_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
          ! Check for breaking.
          call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
          fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
        end if
        ! Sweep 4
        do i_stw = i_st_stw(m_wind_stw), i_en_stw(m_wind_stw)+i_bc3_stw, i_inc_stw(m_wind_stw)
          call prop(nj_stw, ni_stw, na_stw, na_4_stw, anglz_stw, l_sweep_stw, i1_stw, i2_stw, j1_stw, j2_stw,wt1_stw, wt2_stw, e_stw, sea_stw, dth_stw,&
            delr_stw, delr_m_stw, delr_p_stw, ddepdx_stw, ddepdy_stw, wk_stw, c_stwave, cg_stwave, freq_stw, i_stw, j_stw, m4_stw, epsd_stw, pi2_stw, nfreq_stw, twopi_stw, epse_stw)
          if (iprp_stw .eq. 0) then
            local_wind_stw = int(uairdirStwave(j_stw,i_stw)/pi2_stw + 1.00001)
            ! Solve energy source terms in arbitrary depth water.
            if(local_wind_stw .eq. m4_stw) then 
              call gen (nj_stw, ni_stw, na_4_stw, anglz_stw, e_stw, sea_stw, delf_stw, wk_stw,c_stwave, cg_stwave, freq_stw, &
                etot_stw, angav_stw, fma_stw, fm_sea_stw,i_stw, j_stw, m_wind_stw, n_wind_stw, wt1_stw, wt2_stw, i1_stw,i2_stw, j1_stw, j2_stw, idd_stw, uairStwave, uairdirStwave,na_stw, epsd_stw, nfreq_stw, dth_stw,pi2_stw, twopi_stw, epse_stw,angfac_stw)
            end if   
          end if
          !________________________________________________________________________________
          !   Modified by Patricio Moreno 04/02/2010 (added tau_stwave1 to call friction)
          !________________________________________________________________________________
          if(ifric_stw .ge. 1) then 
            call friction (i_stw, j_stw, k_stw, ni_stw, nj_stw, na_4_stw, e_stw, sea_stw,freq_stw, delf_stw, wk_stw, cg_stwave, &
              cf_stw, ifric_stw, delr_stw, l_sweep_stw, m4_stw, tau_stwave1, ubottom_stw, dth_stw,epsd_stw,nfreq_stw,twopi_stw,pi2_stw,epse_stw,na_stw)
          end if 
          ! Check for breaking.
          call break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
          fma_stw(j_stw,i_stw)=freq_stw(kmax_stw)
        end do
      end do
      print *,' Wave ', idd_stw, ';  Complete iteration', iter_stw
    end if
    continue 
  end do

  do j_stw = 1, nj_stw
    do i_stw = 1,ni_stw
      H_stw(j_stw,i_stw) = 0.0e0
      wangle_stw(j_stw,i_stw) = 0.0e0
      sum_stw = 0.0e0
      xcmp_stw = 0.0e0
      ycmp_stw = 0.0e0
      kmax_stw = nfreq_stw
      emax_stw = 0.0e0
      sumfm1_stw = 0.0e0
      do k_stw = 1, nfreq_stw
        e1sum_stw = 0.0e0
        do l_stw = 1, na_stw
          e_all_stw = e_stw(k_stw,l_stw,j_stw,i_stw) + sea_stw(k_stw,l_stw,j_stw,i_stw)
          if (e_all_stw .lt. epse_stw) then
            e_stw(k_stw, l_stw, j_stw, i_stw) = 0.0e0
            sea_stw(k_stw, l_stw, j_stw, i_stw) = 0.0e0
            e_all_stw = 0.0e0
          end if
          e1sum_stw = e1sum_stw + e_all_stw
          sumfm1_stw = sumfm1_stw + e_all_stw/freq_stw(k_stw)*delf_stw(k_stw)*dth_stw
          if (e1sum_stw .gt. emax_stw) then
            emax_stw = e1sum_stw
            kmax_stw = k_stw
          end if
          xcmp_stw = xcmp_stw + (e_all_stw)* cos(anglz_stw(l_stw)) * delf_stw(k_stw)
          ycmp_stw = ycmp_stw + (e_all_stw)* sin(anglz_stw(l_stw)) * delf_stw(k_stw)
        end do
        e1sum_stw = e1sum_stw * dth_stw
        sum_stw = sum_stw + e1sum_stw * delf_stw(k_stw)
        fma_stw(j_stw,i_stw) = freq_stw(kmax_stw)
      end do
      !  H_stw contains H_m0.
      !  wangle_stw contains average wave angle
      if(sum_stw .le. 0.0e0) then 
        sum_stw=0.0e0
      end if 
      H_stw(j_stw,i_stw) = 4.0e0 * sqrt(sum_stw)
      if(sum_stw .gt. 0.0e0) then 
        Tmm1_stw(j_stw,i_stw) = sumfm1_stw/sum_stw
      end if  
      if ((abs(xcmp_stw) + abs(ycmp_stw)) .lt. 1.0e-10) then
        wangle_stw(j_stw, i_stw) = 0.0e0
        !________________________________________________________________________________
        !   Modified by Patricio Moreno 09/2010 (calculating taux_stwave and tauy_stwave per cell)
        !________________________________________________________________________________           
        taux_stwave(j_stw,i_stw)=tau_stwave1(j_stw,i_stw)*cos(wangle_stw(j_stw,i_stw))
        tauy_stwave(j_stw,i_stw)=tau_stwave1(j_stw,i_stw)*sin(wangle_stw(j_stw,i_stw))
      else
        wangle_stw(j_stw,i_stw) = atan2 (ycmp_stw, xcmp_stw) * 180. / pi
        if (wangle_stw(j_stw,i_stw) .lt. -0.05) then 
          wangle_stw(j_stw,i_stw) = wangle_stw(j_stw,i_stw) + 360.
        end if 
        !________________________________________________________________________________
        !   Modified by Patricio Moreno 09/2010 (calculating taux_stwave and tauy_stwave per cell)
        !________________________________________________________________________________             
        taux_stwave(j_stw,i_stw)=tau_stwave1(j_stw,i_stw)*cos(wangle_stw(j_stw,i_stw))
        tauy_stwave(j_stw,i_stw)=tau_stwave1(j_stw,i_stw)*sin(wangle_stw(j_stw,i_stw))
      end if
    end do
  end do

END SUBROUTINE stwave

! ********************************************************************
SUBROUTINE gen(nj_stw,ni_stw,na_4_stw,anglz_stw,e_stw,sea_stw,delf_stw,wk_stw,c_stwave,cg_stwave,freq_stw,etot_stw,angav_stw,  &
                 fma_stw,fm_sea_stw,i_stw,j_stw,m_stw,n_wind_stw,wt1_stw,wt2_stw,i1_stw,i2_stw,j1_stw,j2_stw,idd_stw,uairStwave,uairdirStwave,na_stw, epsd_stw, nfreq_stw, dth_stw,pi2_stw, twopi_stw, epse_stw,angfac_stw)
! ********************************************************************
!
! Purpose: subroutine to solve energy source terms in arbitrary depth
!          water (steady state version) flxcon_stw is the constant
!          epsilon in  eq. 36 -- resio (1987) ucon_stw is the 
!          "universal" equilibrium range constant alpha with a
!          subscript of uairStwave -- see figure 2 -- resio 1988
!          sigb_stw is the jonswap spectral width pararmeter (rear face)
!          rofac_stw is the ration_stw of the density of air to water
!          pwaves_stw is the partitioning constant which defines the
!          percentage of total momemtum entering directly into the wave
!          field to the total momentum transferred from the air to the
!          water. ftchcon_stw is the jonswap coefficient in the
!          nodimensional energy vs fetch relationship (using a friction
!          velocity scaling). zconst_stw is an empirical constant used
!          in limiting how high a peak frequency can be as a function of
!          the total energy. it is totally unimportant except in
!          situations where the peak frequency is such that it no longer
!          represents the energy-containing region of the spectrum. In
!          this case the peak frequency is lowered until the steepness
!          value of hmo divided by the deep-water wavelength is less than
!          2pi/zconst_stw. 
!
! --------------------------------------------------------------------
  ! Subroutine arguments
  integer, intent(in)                                             :: nj_stw
  integer, intent(in)                                             :: ni_stw
  integer, intent(in)                                             :: na_4_stw
  integer, intent(in)                                             :: i_stw
  integer, intent(in)                                             :: j_stw
  integer, intent(in)                                             :: m_stw
  integer, intent(in)                                             :: n_wind_stw
  integer, intent(in), dimension(na_4_stw,4)                      :: i1_stw
  integer, intent(in), dimension(na_4_stw,4)                      :: i2_stw
  integer, intent(in), dimension(na_4_stw,4)                      :: j1_stw
  integer, intent(in), dimension(na_4_stw,4)                      :: j2_stw
  integer, intent(in)                                             :: idd_stw
  integer, intent(in)                                             :: na_stw
  integer, intent(in)                                             :: nfreq_stw
  real, intent(in), dimension(na_stw)                             :: anglz_stw
  real, intent(in), dimension(nfreq_stw)                          :: delf_stw
  real, intent(in), dimension(nfreq_stw,nj_stw,ni_stw)            :: wk_stw
  real, intent(in), dimension(nfreq_stw,nj_stw,ni_stw)            :: c_stwave
  real, intent(in), dimension(nfreq_stw,nj_stw,ni_stw)            :: cg_stwave
  real, intent(in), dimension(nfreq_stw)                          :: freq_stw
  real, intent(in)                                                :: epsd_stw
  real, intent(in), dimension(nj_stw,ni_stw)                      :: uairStwave
  real, intent(in), dimension(nj_stw,ni_stw)                      :: uairdirStwave
  real, intent(in), dimension(na_4_stw,4)                         :: wt1_stw
  real, intent(in), dimension(na_4_stw,4)                         :: wt2_stw  
  real, intent(in)                                                :: dth_stw
  real, intent(in)                                                :: pi2_stw
  real, intent(in)                                                :: twopi_stw
  real, intent(in)                                                :: epse_stw
  real, intent(inout), dimension(nj_stw,ni_stw)                   :: etot_stw
  real, intent(inout), dimension(nj_stw,ni_stw)                   :: angav_stw
  real, intent(inout), dimension(nj_stw,ni_stw)                   :: fma_stw
  real, intent(inout), dimension(nj_stw,ni_stw)                   :: fm_sea_stw
  real, intent(inout), dimension(nfreq_stw,na_stw,nj_stw,ni_stw)  :: e_stw
  real, intent(inout), dimension(nfreq_stw,na_stw,nj_stw,ni_stw)  :: sea_stw

  integer :: k_stw, l_stw, i_angdif_stw, iwind_band_stw, iband1_stw
  integer :: iband2_stw, l1_stw, l2_stw, l3_stw, l4_stw, mm_stw
  integer :: nn_stw, ii1_stw, ii2_stw, jj1_stw, jj2_stw
  real :: fmn_stw, flxcon_stw, ucon_stw, sigb_stw, rofac_stw, pwaves_stw
  real :: ftchcon_stw, zconst_stw, sqrtg_stw, cdrag_stw, gamma_stw
  real :: ust_stw, xlam_stw, alpstr_stw, zxx_stw, etotzz_stw, xsum_stw 
  real :: ysum_stw, emax_stw, kmax_stw
  real :: e1sum_stw, angav_sea_stw, ang_calc_stw, angdif_stw, wkdp_stw
  real :: bconst_stw, udtmp_stw, depij_stw, fm_stw, zfm_stw, tmn_stw
  real :: wkm_stw, cm_stw, tkd_stw, cgm_stw, delta_stw, time_stw, wkh_stw
  real :: flux_stw, eout_stw, uocm_stw, fmbar_stw, afact_stw, efcon_stw
  real :: btermz_stw, weight1_stw, ee2o_stw, dep_up_stw, ee2_stw, sum_stw
  real :: einpt_stw, e_upwind_stw, efin_stw, xxx_stw, equil_stw, zzz_stw
  real :: zzz2_stw, yyy_stw, equilz_stw, xxo_stw, yyo_stw, equilo_stw
  real :: dele1_stw, dele_stw, sum_test_stw, sum_test2_stw, ration_stw
  real :: ein_stw, ediff_stw, eee_stw, exx_stw, depfac_stw, ratt_stw, sume1_stw
  real, dimension(nfreq_stw) :: e1_stw
  real, dimension(nfreq_stw,na_stw) :: edens_stw
  real, dimension(nfreq_stw) :: zjac_stw
  real, dimension(nfreq_stw) :: fka_stw
  real, intent(inout), dimension(na_stw) :: angfac_stw

  ! check on zero depth
  fmn_stw = freq_stw(nfreq_stw)

  if (dep_stwave(j_stw,i_stw) .le. epsd_stw) then
    ! fmn_stw = freq_stw(nfreq_stw)
    do k_stw = 1, nfreq_stw
      e1_stw(k_stw) = 0.0e0
      do l_stw = 1, na_stw
        e_stw(k_stw, l_stw, j_stw, i_stw) = 0.0e0
        sea_stw(k_stw, l_stw, j_stw, i_stw) = 0.0e0
      end do
    end do
    etot_stw(j_stw, i_stw) = 0.0e0
    fm_sea_stw(j_stw,i_stw) = fmn_stw
    angav_stw(j_stw, i_stw) = uairdirStwave(j_stw, i_stw)
  else
    flxcon_stw = 30.0e0
    ucon_stw = 0.004e0
    sigb_stw = 0.09e0
    rofac_stw = 0.00123e0
    pwaves_stw = 0.75e0
    ftchcon_stw = 1.4e-4
    !         zconst_stw = 16.5e0
    zconst_stw = 1.0e0
    sqrtg_stw = sqrt(g)
    cdrag_stw = (1.2e0 + 0.025e0*uairStwave(j_stw, i_stw))*1.0e-3
    gamma_stw = 1.5e0
    ust_stw = sqrt(cdrag_stw)*uairStwave(j_stw, i_stw)
    xlam_stw = 1.75e0
    alpstr_stw = 0.0390e0
    zxx_stw = 0.1715e0/(alpstr_stw*xlam_stw)

    etotzz_stw = 0.0e0
    xsum_stw = 0.0e0
    ysum_stw = 0.0e0
    emax_stw = 0.0e0
    kmax_stw = nfreq_stw
    edens_stw(:,:) =  0.0e0
    ! estimate wind sea_stw freq_stw and direction
    do k_stw = 1, nfreq_stw
      e1sum_stw = 0.0e0
      do l_stw = 1, na_stw
        if (e_stw(k_stw, l_stw, j_stw, i_stw) .lt. epse_stw) then
          e_stw(k_stw, l_stw, j_stw, i_stw) = 0.0e0
        end if
        e1sum_stw = e1sum_stw + e_stw(k_stw, l_stw, j_stw, i_stw)
        xsum_stw = xsum_stw + e_stw(k_stw, l_stw, j_stw, i_stw)*delf_stw(k_stw)*cos(anglz_stw(l_stw))
        ysum_stw = ysum_stw + e_stw(k_stw, l_stw, j_stw, i_stw)*delf_stw(k_stw)*sin(anglz_stw(l_stw))
      end do
      e1sum_stw = e1sum_stw*dth_stw
      if (e1sum_stw .gt. emax_stw) then
        emax_stw = e1sum_stw
        kmax_stw = k_stw
      end if
        etotzz_stw = etotzz_stw + e1sum_stw*delf_stw(k_stw)
        e1_stw(k_stw) = e1sum_stw
    end do
    ! fm_sea_stw(j_stw, i_stw) = min(fm_sea_stw(j_stw,i_stw),freq_stw(kmax_stw))
    if (etotzz_stw .lt. 1.0e-3) then
      angav_stw(j_stw, i_stw) = uairdirStwave(j_stw, i_stw)
    else
      angav_stw(j_stw, i_stw) = atan2 (ysum_stw, xsum_stw)
    end if
    angav_sea_stw = angav_stw(j_stw, i_stw)

    ! determine if existing wave field is within +/-90 of winds
    ! if not, start growth of new wave train
    ang_calc_stw = cos(uairdirStwave(j_stw, i_stw))*cos(angav_stw(j_stw,i_stw))+sin(uairdirStwave(j_stw, i_stw))*sin(angav_stw(j_stw,i_stw))
    if(ang_calc_stw .gt. 1.0) then
      ang_calc_stw = 1.0
    end if
    if(ang_calc_stw .lt. -1.0) then 
      ang_calc_stw = -1.0
    end if 
    angdif_stw = acos(ang_calc_stw)
    i_angdif_stw = 0
    if(abs(angdif_stw) .gt. pi2_stw) then
      i_angdif_stw = 1
      ! determine the angle band with the wind and sum_stw from +/- 60 deg (pi/3) to define sea_stw
      iwind_band_stw = int(uairdirStwave(j_stw,i_stw)/dth_stw + 0.5) + 1
      iband1_stw = iwind_band_stw - int(pi/3./dth_stw)
      iband2_stw = iwind_band_stw + int(pi/3./dth_stw)
      l1_stw = iband1_stw
      l2_stw = iband2_stw
      l3_stw = 0
      l4_stw = 0
      if(iband1_stw .lt. 1) then
        l1_stw = 1
        l3_stw = iband1_stw + na_stw
        l4_stw = na_stw
      end if
      if(iband2_stw .gt. na_stw) then
        l2_stw = na_stw
        l3_stw = 1
        l4_stw = iband2_stw - na_stw
      end if
      etotzz_stw = 0.0e0
      xsum_stw = 0.0e0
      ysum_stw = 0.0e0
      emax_stw = 0.0e0
      kmax_stw = nfreq_stw
      ! estimate wind sea_stw freq_stw and direction
      do k_stw = 1, nfreq_stw
        e1sum_stw = 0.0e0
        do l_stw = l1_stw,l2_stw
          e1sum_stw = e1sum_stw + e_stw(k_stw, l_stw, j_stw, i_stw)
          xsum_stw = xsum_stw + e_stw(k_stw, l_stw, j_stw, i_stw)*delf_stw(k_stw)*cos(anglz_stw(l_stw))
          ysum_stw = ysum_stw + e_stw(k_stw, l_stw, j_stw, i_stw)*delf_stw(k_stw)*sin(anglz_stw(l_stw))
        end do
        if(l3_stw .gt. 0) then
          do l_stw = l3_stw,l4_stw
            e1sum_stw = e1sum_stw + e_stw(k_stw, l_stw, j_stw, i_stw)
            xsum_stw = xsum_stw + e_stw(k_stw, l_stw, j_stw, i_stw)*delf_stw(k_stw)*cos(anglz_stw(l_stw))
            ysum_stw = ysum_stw + e_stw(k_stw, l_stw, j_stw, i_stw)*delf_stw(k_stw)*sin(anglz_stw(l_stw))
          end do
        end if
        e1sum_stw = e1sum_stw*dth_stw
        if (e1sum_stw .gt. emax_stw) then
          emax_stw = e1sum_stw
          kmax_stw = k_stw
        end if
        etotzz_stw = etotzz_stw + e1sum_stw*delf_stw(k_stw)
        e1_stw(k_stw) = e1sum_stw
      end do
      ! fm_sea_stw(j_stw, i_stw) = freq_stw(kmax_stw)
      if (etotzz_stw .lt. 1.0e-3) then
        angav_sea_stw = uairdirStwave(j_stw, i_stw)
        fm_sea_stw(j_stw, i_stw) = freq_stw(nfreq_stw)
      else
        angav_sea_stw = atan2 (ysum_stw, xsum_stw)
      end if
      ang_calc_stw = cos(uairdirStwave(j_stw, i_stw))*cos(angav_sea_stw)+sin(uairdirStwave(j_stw, i_stw))*sin(angav_sea_stw)
      if(ang_calc_stw.gt.1.0) then 
        ang_calc_stw = 1.0
      end if 
      if(ang_calc_stw.lt.-1.0) then 
        ang_calc_stw = -1.0
      end if 
      angdif_stw = acos(ang_calc_stw)
    end if

    ! wk_stw is the wavenumber matrix
    ! zjac_stw is the finite-depth jacobian matrix
    do k_stw = 1, nfreq_stw
      wkdp_stw = (twopi_stw*freq_stw(k_stw))**2/g
      zjac_stw(k_stw) = (wkdp_stw/wk_stw(k_stw,j_stw,i_stw))**3
      if (zjac_stw(k_stw) .lt. 1.0e-6) then
        zjac_stw(k_stw) = 0.0e0
      end if 
    end do

    ! fka_stw is an array containing the deep-water equilibrium range
    ! values for each frequency
    do k_stw = 1, nfreq_stw
      fka_stw(k_stw) = g/(twopi_stw**3*freq_stw(k_stw)**4)
    end do

    ! bconst_stw is a constant which, when multiplied times nondimensional peak
    ! frequency (which is proportional to ustar/(peak phase velocity)
    ! in deep water) and ustar**2 produces a wind input term consistent
    ! with equation 2.29 in resio and perrie (1989).
    ! it is derived to produce a constant percentage of total momentum
    ! flux_stw into the wave field.
    bconst_stw = 2.0e0*rofac_stw*pwaves_stw*twopi_stw*twopi_stw/(1.75e0*ucon_stw*g)
    ! note: the sqrt(0.0013) term is added here since the fmbar_stw term is now
    ! defined in terms of friction velocity. the 0.0013 factor is an
    ! estimate of the average coefficient of drag in past wind input
    ! experiments (typically performed at about 8-12 m_stw/sec)
    bconst_stw = bconst_stw/sqrt(0.0013e0)
    udtmp_stw = uairdirStwave(j_stw, i_stw)
    depij_stw = dep_stwave(j_stw, i_stw)
    ! estimate "upstream fm_stw" for wave growth
    mm_stw = int(uairdirStwave(j_stw,i_stw)/pi2_stw + 1.00001)
    nn_stw =  (uairdirStwave(j_stw,i_stw)-(mm_stw-1)*pi2_stw)/dth_stw + 1.00001
    ! note tracing back using same interpolation as wave propagation
    ! print *,i_stw,j_stw,uairdirStwave(j_stw,i_stw),mm_stw,nn_stw
    ii1_stw = i_stw + i1_stw(nn_stw,mm_stw)
    ii2_stw = i_stw + i2_stw(nn_stw,mm_stw)
    jj1_stw = j_stw + j1_stw(nn_stw,mm_stw)
    jj2_stw = j_stw + j2_stw(nn_stw,mm_stw)
    if(wt1_stw(nn_stw,mm_stw) .lt. 0.0001) then
      ii1_stw = ii2_stw
      jj1_stw = jj2_stw
    elseif(wt2_stw(nn_stw,mm_stw) .lt. 0.0001) then
      ii2_stw = ii1_stw
      jj2_stw = jj1_stw
    end if
    ! mod here for testing
    ! fm_stw = fm_sea_stw(j_stw,i_stw)
    fm_stw = wt1_stw(nn_stw,mm_stw)*fm_sea_stw(jj1_stw,ii1_stw)+wt2_stw(nn_stw,mm_stw)*fm_sea_stw(jj2_stw,ii2_stw)
    zfm_stw = sqrt(g/(twopi_stw*zconst_stw*4.01e0*sqrt(etotzz_stw) + 1.0e-10))
    ! mod 7/13/2005 Prevent fm_stw=0.0
    if (fm_stw .eq. 0.0) then
      fm_stw = zfm_stw
    end if 
    fm_stw = min (fm_stw, zfm_stw)
    fmn_stw = fm_stw
    tmn_stw = 1.0e0/fmn_stw
    ! test on zconst_stw sensitivity
    ! see note under definition of zconst_stw to interpret zfm_stw.!
    ! solve for group velocity of sea_stw component
    wkm_stw = wkfnc_stw(fmn_stw,depij_stw)
    cm_stw = fmn_stw*twopi_stw/wkm_stw
    tkd_stw = min(2.0e0*wkm_stw*depij_stw, 20.0e0)
    cgm_stw = 0.5e0*cm_stw*(1.0e0 + (2.0e0*wkm_stw*depij_stw/sinh(tkd_stw)))
    ! note: time_stw is estimated here based on a wave travel speed of 90% of the
    ! group velocity of the spectral peak.  as discussed on page 59 of
    ! resio (1988), the "beta" factor for converting to a <mean> speed
    ! of propagation (see eqs 11, 12, and 13) typically falls between
    ! 0.8 and 1.0; hence the motivation for the choice of 0.9 here.
    delta_stw = min(abs(dx/cos(angav_sea_stw)),abs(dy/sin(angav_sea_stw)))
    time_stw = delta_stw/(0.9e0*cgm_stw)
    wkh_stw = min (wkm_stw*depij_stw, 20.0e0)
    ! note: time_stw is estimated here based on a wave travel speed of 100% of the
    ! group velocity of the spectral peak.  as discussed on page 59 of
    ! resio (1988), the "beta" factor for converting to a <mean> speed
    ! of propagation (see eqs 11, 12, and 13) typically falls between
    ! 0.8 and 1.0.  the appropriate estimate for swell is 1.

    ! flux_stw is the energy flux_stw to high frequencies (where it is presumed to
    ! be lost to wave breaking).  see eq. 36 -- resio (1987)
    ! this flux_stw is due to the higher frequency wavetrain.
    flux_stw = flxcon_stw*sqrtg_stw*etotzz_stw**3*wkm_stw**4.5e0/tanh(wkh_stw)**0.75e0

    ! flux2 is the energy flux_stw to high frequencies (analogous to flux_stw).
    ! this flux_stw is due to the lower frequency wavetrain (if one exists).
    eout_stw = flux_stw*time_stw
     
    ! section of code for evolution of fm_stw

    ! ust_stw is the friction velocity
    ! xlam_stw is an estimated "average" value for lamda in eq 2.5 -- resio and
    ! perrie (1989)
    ! alpstr_stw is the alpha parameter consistent with wave growth in dwave.
    ! zxx_stw is the fetch-growth constant consistent with wave growth in dwave.

    ! check on km*H_stw limit for fm_stw evolution
    if (wkh_stw .ge. 0.6e0) then
      ! check on p-m_stw limit for fm_stw evolution
      uocm_stw = uairStwave(j_stw, i_stw)/cm_stw
      if (uocm_stw .ge. 0.98e0) then
        tmn_stw = (1.0e0/fm_stw**3.333e0 + zxx_stw*delta_stw*    &
        ust_stw**1.333e0/g**2.333e0)**0.3e0
        fmn_stw = 1.0e0/(tmn_stw + 1.0e-10)
      end if
    end if
    fmbar_stw = fmn_stw*ust_stw/g
    if(abs(angdif_stw).gt.pi2_stw)then
      afact_stw = 0.
    else
      afact_stw = cos((angdif_stw)*0.5e0)**2
    end if
    efcon_stw = afact_stw*ucon_stw*uairStwave(j_stw,i_stw)
    btermz_stw = fmbar_stw*bconst_stw*cdrag_stw*afact_stw

    ! weight1_stw is a relaxation parameter    weight1_stw=1 for tma type adjustment
    !                                  weight1_stw=0 for no tma adjustment
    !                                  0<weight1_stw<1 for in between adjustment
    weight1_stw = min ((etotzz_stw - eout_stw)/etotzz_stw, 1.0e0)
    if (etotzz_stw .lt. 1.0e-4) then
      weight1_stw = 0.0e0
    end if 
    weight1_stw = max (weight1_stw, 0.0e0)

    ! ee2o_stw is the estimated energy in spectrum with peak frequency at fm_stw
    ee2o_stw = gamma_stw*afact_stw*ucon_stw*uairStwave(j_stw,i_stw)*g/(twopi_stw**3*fm_stw**4)
    dep_up_stw = wt1_stw(nn_stw,mm_stw)*dep_stwave(jj1_stw,ii1_stw)+ wt2_stw(nn_stw,mm_stw)*dep_stwave(jj2_stw,ii2_stw)
    if (dep_up_stw .lt. epsd_stw) then
      ee2o_stw = 0.0e0
    end if 
    ! ee2_stw is the estimated energy in spectrum with peak frequency at fmn_stw
    ee2_stw = gamma_stw*afact_stw*ucon_stw*uairStwave(j_stw,i_stw)*g/(twopi_stw**3*fmn_stw**4)
    do l_stw = 1, na_stw
      ang_calc_stw = cos(uairdirStwave(j_stw,i_stw)) * cos(anglz_stw(l_stw)) + sin(uairdirStwave(j_stw,i_stw)) * sin(anglz_stw(l_stw))
      if (abs(ang_calc_stw) .gt. 1.0) then
        angfac_stw(l_stw) = 0.0e0
      else
        angdif_stw = acos(ang_calc_stw)
        if (abs(angdif_stw) .gt. pi2_stw) then
          angdif_stw = pi2_stw
        end if 
        angfac_stw(l_stw) = cos(angdif_stw)
      end if
      if (angfac_stw(l_stw) .lt. 1.0e-2) then
        angfac_stw(l_stw) = 0.0e0
      end if
    end do

    ! start frequency loop
    do k_stw = 1, nfreq_stw
      if (freq_stw(k_stw) .ge. fmn_stw) then
        sum_stw = 0.0e0
        ! entering rear face of spectrum source (f > fm_stw)
        einpt_stw = (ee2_stw)*16./(5.*pi)
        if(i_angdif_stw.eq.0)then
          do l_stw = 1, na_stw
            ! btermz_stw is the estimated wind input source function
            e_upwind_stw = wt1_stw(nn_stw,mm_stw)*(e_stw(k_stw,l_stw,jj1_stw,ii1_stw)+sea_stw(k_stw,l_stw,jj1_stw,ii1_stw)) + wt2_stw(nn_stw,mm_stw)*(e_stw(k_stw,l_stw,jj2_stw,ii2_stw)+sea_stw(k_stw,l_stw,jj2_stw,ii2_stw))
            efin_stw = btermz_stw*0.5e0*(e_stw(k_stw, l_stw, j_stw, i_stw) + e_upwind_stw)*time_stw*angfac_stw(l_stw)
            if (freq_stw(k_stw) .lt. fm_stw) then
              efin_stw = efin_stw + einpt_stw*angfac_stw(l_stw)**6
            end if 
            edens_stw(k_stw, l_stw) = e_stw(k_stw, l_stw, j_stw, i_stw) + efin_stw
            sum_stw = sum_stw + edens_stw(k_stw, l_stw)*dth_stw
          end do
        else
          edens_stw(k_stw,:) = e_stw(k_stw, :, j_stw, i_stw)
          do l_stw = l1_stw, l2_stw
            ! btermz_stw is the estimated wind input source function
            e_upwind_stw = wt1_stw(nn_stw,mm_stw)*(e_stw(k_stw,l_stw,jj1_stw,ii1_stw)+sea_stw(k_stw,l_stw,jj1_stw,ii1_stw)) + wt2_stw(nn_stw,mm_stw)*(e_stw(k_stw,l_stw,jj2_stw,ii2_stw)+sea_stw(k_stw,l_stw,jj2_stw,ii2_stw))
            efin_stw = btermz_stw*0.5e0*(e_stw(k_stw, l_stw, j_stw, i_stw) + e_upwind_stw)*time_stw*angfac_stw(l_stw)
            if (freq_stw(k_stw) .lt. fm_stw) then
              efin_stw = efin_stw + einpt_stw*angfac_stw(l_stw)**6
            end if
            edens_stw(k_stw, l_stw) = e_stw(k_stw, l_stw, j_stw, i_stw) + efin_stw
            sum_stw = sum_stw + edens_stw(k_stw, l_stw)*dth_stw
          end do
          if(l3_stw .gt. 0) then
            do l_stw = l3_stw,l4_stw
              ! btermz_stw is the estimated wind input source function
              e_upwind_stw = wt1_stw(nn_stw,mm_stw)*(e_stw(k_stw,l_stw,jj1_stw,ii1_stw)+sea_stw(k_stw,l_stw,jj1_stw,ii1_stw)) + wt2_stw(nn_stw,mm_stw)*(e_stw(k_stw,l_stw,jj2_stw,ii2_stw)+sea_stw(k_stw,l_stw,jj2_stw,ii2_stw))
              efin_stw = btermz_stw*0.5e0*(e_stw(k_stw, l_stw, j_stw, i_stw) + e_upwind_stw)*time_stw*angfac_stw(l_stw)
              if (freq_stw(k_stw) .lt. fm_stw) then
                efin_stw = efin_stw + einpt_stw*angfac_stw(l_stw)**6
              end if 
              edens_stw(k_stw, l_stw) = e_stw(k_stw, l_stw, j_stw, i_stw) + efin_stw
              sum_stw = sum_stw + edens_stw(k_stw, l_stw)*dth_stw
            end do
          end if
        end if
        xxx_stw =  -0.5e0*((freq_stw(k_stw) - fmn_stw)/(sigb_stw*fmn_stw))**2
        xxx_stw = gamma_stw**exp(max (xxx_stw, -15.0e0))
        ! equil_stw is the equilibrium range value for this point in the spectrum,
        ! including the jacobian term "zjak"
        equil_stw = efcon_stw*fka_stw(k_stw)*xxx_stw*zjac_stw(k_stw)
        zzz_stw = weight1_stw + (1.0e0 - weight1_stw)*equil_stw/(e1_stw(k_stw) + 1.0e-6)
        e1_stw(k_stw) = zzz_stw*(e1_stw(k_stw) + 1.0e-6)
        ! zzz2_stw is normalization coeffiecient to force 2-d_stw spectral densities to
        ! integrate to 1-d_stw value
        zzz2_stw = e1_stw(k_stw)/(sum_stw+1.0e-6)
        if(i_angdif_stw .eq. 0) then
          do l_stw = 1, na_stw
            edens_stw(k_stw,l_stw) = edens_stw(k_stw,l_stw)*zzz2_stw
          end do
        else
          do l_stw = l1_stw,l2_stw
            edens_stw(k_stw,l_stw) = edens_stw(k_stw,l_stw)*zzz2_stw
          end do
          if(l3_stw .gt. 0) then
            do l_stw = l3_stw,l4_stw
              edens_stw(k_stw,l_stw) = edens_stw(k_stw,l_stw)*zzz2_stw
            end do
          end if
        end if
      else
        ! entering front face of spectrum source (f .le. fm_stw)
        xxx_stw = fmn_stw/freq_stw(k_stw)
        yyy_stw = max (1.0e0 - xxx_stw**4, -15.0e0)
        ! equilz_stw is the estimated spectral density at the new spectral peak
        ! frequency
        equilz_stw = ee2_stw*exp(yyy_stw)*zjac_stw(k_stw)
        xxo_stw = fm_stw/freq_stw(k_stw)
        yyo_stw = max (1.0e0 - xxo_stw**4, -15.0e0)
        ! equilo_stw is the estimated spectral density at the old spectral peak
        ! frequency
        equilo_stw = ee2o_stw*exp(yyo_stw)*zjac_stw(k_stw)
        ! dele1_stw is the "gain" in energy at the kth frequency (on the front
        ! portion of the spectrum)
        dele1_stw = max ((equilz_stw - equilo_stw), 0.0e0)
        ! change units for directional spectrum energy densities
        dele_stw = dele1_stw/(5./16.*pi)
        sum_stw = 0.0e0
        sum_test_stw = 0.0e0
        sum_test2_stw = 0.0e0
        if (equilz_stw .lt. 1.0e-6) then
          equilz_stw = 0.0e0
        end if 
        do l_stw = 1, na_stw
          if (e_stw(k_stw, l_stw, j_stw, i_stw) .lt. epse_stw) then
            e_stw(k_stw, l_stw, j_stw, i_stw) = 0.0e0
          end if
          ! edens_stw is the new spectral energy for the kth  frequency and lth angle
          ! component
          edens_stw(k_stw, l_stw) = e_stw(k_stw, l_stw, j_stw, i_stw) + dele_stw*angfac_stw(l_stw)**6
          sum_stw = sum_stw + edens_stw(k_stw, l_stw)
        end do
        e1_stw(k_stw) = e1_stw(k_stw) + dele1_stw
        ration_stw = 1.0e0
        do l_stw = 1, na_stw
          edens_stw(k_stw, l_stw) = edens_stw(k_stw, l_stw)*ration_stw
        end do
          ! completed front face of spectrum source
      end if
      ! end of rear portion of spectrum source
    enddo
    ! ein_stw is the total integrated energy gain expected to be retained by
    ! the waves
    ein_stw = pwaves_stw*rofac_stw*0.85e0*cm_stw*ust_stw*ust_stw/g*time_stw
    ! normalize to integral energy changes
    ! check on relative depth (depth/(deep-water wavelength))
    ! if this is greater than 1.5 we can ignore energy losses.
    ! this keeps a small amount of energy from "leaking" out of the
    ! system during propagation in deep water.
    if (depij_stw*twopi_stw*fmn_stw*fmn_stw/g .gt. 0.5e0) then
      eout_stw = 0.0e0
    end if
    ediff_stw = ein_stw - eout_stw
    ! rescale to fetch-growth case for ein_stw >eout_stw
    if (ediff_stw .gt. 0.0e0) then
      ediff_stw = delta_stw*ust_stw*ust_stw/g*ftchcon_stw
    end if
    eee_stw = max (etotzz_stw + ediff_stw, 0.0e0)
    if (eee_stw .lt. 0.0e0) then
      eee_stw = 0.0e0
    end if
    ! if eee_stw is not very different from etot_stw then there is no need for more
    ! involved analysis as performed between branch point and 2095 label
    if (eee_stw .lt. 0.9e0*etotzz_stw) then
      exx_stw = max (etotzz_stw, 1.0e-12)
      depfac_stw = eout_stw/(exx_stw**3*time_stw)
      eee_stw = etotzz_stw-(eout_stw-ein_stw)*exp(-depfac_stw*time_stw)
    end if
    ! integrate spectrum to get total energy and direction components
    sum_stw = 0.0e0
    xsum_stw = 0.0e0
    ysum_stw = 0.0e0
    kmax_stw = nfreq_stw
    do k_stw = 1, nfreq_stw
      do l_stw = 1, na_stw
        if (edens_stw(k_stw, l_stw) .lt. 1.0e-10) then
          edens_stw(k_stw, l_stw) = 0.0e0
        end if
        sum_stw = sum_stw + edens_stw(k_stw, l_stw)*delf_stw(k_stw)
        xsum_stw = xsum_stw + edens_stw(k_stw, l_stw)*cos(anglz_stw(l_stw))
        ysum_stw = ysum_stw + edens_stw(k_stw, l_stw)*sin(anglz_stw(l_stw))
      end do
    end do

    sum_stw = max (sum_stw, 1.0e-10)
    ! ratt_stw ensures that 2-d_stw energies sum_stw to integral balance
    !            ratt_stw = min (eee_stw/(sum_stw*dth_stw), 1.0e0)
    ratt_stw = 1.0e0
    sum_stw = sum_stw*ratt_stw
    sume1_stw=0.
    do k_stw = 1, nfreq_stw
      e1_stw(k_stw) = e1_stw(k_stw)*ratt_stw
      sume1_stw=sume1_stw+e1_stw(k_stw)*delf_stw(k_stw)
      if (e1_stw(k_stw) .gt. e1_stw(kmax_stw)) then
        kmax_stw = k_stw
      end if
      do l_stw = 1, na_stw
        sea_stw(k_stw, l_stw, j_stw, i_stw) = edens_stw(k_stw, l_stw)*ratt_stw - e_stw(k_stw, l_stw, j_stw, i_stw)
      end do
    end do
    sum_stw = sum_stw*dth_stw
    etot_stw(j_stw, i_stw) = sum_stw
    xsum_stw = xsum_stw + 1.0e-10
    angav_stw(j_stw, i_stw) = atan2 (ysum_stw, xsum_stw)
    if (etot_stw(j_stw, i_stw) .lt. 1.0e-3) then
      angav_stw(j_stw, i_stw) = uairdirStwave(j_stw, i_stw)
    end if
    fm_sea_stw(j_stw, i_stw) = fmn_stw
  end if

END SUBROUTINE gen

! ********************************************************************
SUBROUTINE friction (i_stw, j_stw, k_stw, ni_stw, nj_stw, na_4_stw, e_stw, sea_stw, freq_stw, delf_stw, &
                           wk_stw, cg_stwave, cf_stw, ifric_stw, delr_stw, l_sweep_stw, m_stw, tau_stwave1, ubottom_stw, dth_stw,epsd_stw,nfreq_stw,twopi_stw,pi2_stw,epse_stw,na_stw)
! ********************************************************************
!
! Purpose: To estimate the bottom friction induced by waves
!          Modified by Patricio Moreno 04/05/2010 (added tau_stwave1
!          to subroutine friction)
!
! --------------------------------------------------------------------
  ! Subroutine arguments
  integer, intent(in)                                             :: i_stw
  integer, intent(in)                                             :: j_stw
  integer, intent(in)                                             :: k_stw
  integer, intent(in)                                             :: ni_stw
  integer, intent(in)                                             :: nj_stw
  integer, intent(in)                                             :: na_4_stw
  integer, intent(in)                                             :: nfreq_stw
  integer, intent(in)                                             :: na_stw
  integer, intent(in)                                             :: ifric_stw
  integer, intent(in), dimension(na_4_stw,4)                      :: l_sweep_stw
  integer, intent(in)                                             :: m_stw
  real, intent(in)                                                :: dth_stw
  real, intent(in)                                                :: epsd_stw
  real, intent(in)                                                :: twopi_stw
  real, intent(in)                                                :: pi2_stw
  real , intent(in)                                               :: epse_stw
  real, intent(in), dimension(nfreq_stw)                          :: freq_stw
  real, intent(in), dimension(nfreq_stw)                          :: delf_stw
  real, intent(in), dimension(nfreq_stw,nj_stw,ni_stw)            :: wk_stw
  real, intent(in), dimension(nfreq_stw,nj_stw,ni_stw)            :: cg_stwave
  real, intent(in), dimension(nj_stw,ni_stw)                      :: cf_stw
  real, intent(in), dimension(na_4_stw,4)                         :: delr_stw
  real, intent(inout), dimension(nj_stw,ni_stw)                   :: tau_stwave1
  real, intent(inout), dimension(nj_stw,ni_stw)                   :: ubottom_stw
  real, intent(inout), dimension(nfreq_stw,na_stw,nj_stw,ni_stw)  :: e_stw
  real, intent(inout), dimension(nfreq_stw,na_stw,nj_stw,ni_stw)  :: sea_stw

  integer :: l_stw, n_stw, kk_stw
  real :: urms_stw, wkd_stw, sigma_stw, b_stw

  ! Calculate bottom friction.
  ! epse_stw (a small value of energy) = 1.0e-8
  if(dep_stwave(j_stw,i_stw) .gt. epsd_stw) then
    urms_stw=0.0
    do kk_stw = 1, nfreq_stw
      do l_stw = 1, na_stw
        if (e_stw(kk_stw,l_stw,j_stw,i_stw) + sea_stw(kk_stw,l_stw,j_stw,i_stw) .gt. epse_stw .and. wk_stw(kk_stw,j_stw,i_stw) .gt. 0.0) then
          wkd_stw = min(wk_stw(kk_stw,j_stw,i_stw) * dep_stwave(j_stw,i_stw),20.)
          urms_stw = urms_stw + (twopi_stw*freq_stw(kk_stw)/sinh(wkd_stw))**2 * (e_stw(kk_stw,l_stw,j_stw,i_stw)+sea_stw(kk_stw,l_stw,j_stw,i_stw))*delf_stw(kk_stw)*dth_stw
        end if
      end do
    end do
    urms_stw = sqrt(urms_stw)
    ubottom_stw(j_stw,i_stw) = urms_stw !added ubottom_stw(j_stw,i_stw). MODIFIED BY PM 07,2010

    !added by PM 07/2011: calculate tau_stwave1 for water depths greater than zero.
    if (dep_stwave(j_stw,i_stw) .le. 0) then
      tau_stwave1(j_stw,i_stw) = 0.0
    else
      tau_stwave1(j_stw,i_stw) = 1000*9.81*cf_stw(j_stw,i_stw)**2*dep_stwave(j_stw,i_stw)**(-0.3333)*ubottom_stw(j_stw,i_stw)*abs(ubottom_stw(j_stw,i_stw))
    endif
    do kk_stw = 1, nfreq_stw
      do n_stw = 1, na_4_stw
        l_stw = l_sweep_stw(n_stw,m_stw)
        if(e_stw(kk_stw,l_stw,j_stw,i_stw)+sea_stw(kk_stw,l_stw,j_stw,i_stw) .gt. epse_stw) then
          wkd_stw = min(wk_stw(kk_stw,j_stw,i_stw)*dep_stwave(j_stw,i_stw),20.0e0)
          sigma_stw = twopi_stw*freq_stw(kk_stw)
          if(ifric_stw .le. 2) then
            b_stw = -cf_stw(j_stw,i_stw)/g*(sigma_stw/sinh(wkd_stw))**2
            e_stw(kk_stw, l_stw, j_stw, i_stw) = e_stw(kk_stw, l_stw, j_stw, i_stw)*exp(b_stw*delr_stw(n_stw,m_stw)/cg_stwave(kk_stw,j_stw,i_stw))
            sea_stw(kk_stw, l_stw, j_stw, i_stw) = sea_stw(kk_stw, l_stw, j_stw, i_stw)*exp(b_stw*delr_stw(n_stw,m_stw)/cg_stwave(kk_stw,j_stw,i_stw))
          else
            b_stw = -cf_stw(j_stw,i_stw)*cf_stw(j_stw,i_stw)*dep_stwave(j_stw,i_stw)**(-0.3333)*(sigma_stw/sinh(wkd_stw))**2
            e_stw(kk_stw, l_stw, j_stw, i_stw) = e_stw(kk_stw, l_stw, j_stw, i_stw)*exp(b_stw*urms_stw*delr_stw(n_stw,m_stw)/cg_stwave(kk_stw,j_stw,i_stw))
            sea_stw(kk_stw, l_stw, j_stw, i_stw) = sea_stw(kk_stw, l_stw, j_stw, i_stw)*exp(b_stw*urms_stw*delr_stw(n_stw,m_stw)/cg_stwave(kk_stw,j_stw,i_stw))
          end if
        end if
      end do
    end do
  end if

END SUBROUTINE friction

! ********************************************************************
subroutine prop(nj_stw,ni_stw,na_stw,na_4_stw,anglz_stw,l_sweep_stw,i1_stw,i2_stw,j1_stw,j2_stw,wt1_stw,wt2_stw,e_stw,sea_stw,dth_stw,  &
          delr_stw,delr_m_stw,delr_p_stw,ddepdx_stw,ddepdy_stw,wk_stw,c_stwave,cg_stwave,freq_stw,i_stw,j_stw,m_stw, epsd_stw, pi2_stw, nfreq_stw, twopi_stw, epse_stw)
! ********************************************************************
!
! Purpose: 
!
! --------------------------------------------------------------------
  ! Subroutine arguments
  integer, intent(in)                                            :: nj_stw
  integer, intent(in)                                            :: ni_stw
  integer, intent(in)                                            :: na_stw
  integer, intent(in)                                            :: na_4_stw
  integer, intent(in), dimension(na_4_stw,4)                     :: l_sweep_stw
  integer, intent(in), dimension(na_4_stw,4)                     :: i1_stw
  integer, intent(in), dimension(na_4_stw,4)                     :: i2_stw
  integer, intent(in), dimension(na_4_stw,4)                     :: j1_stw
  integer, intent(in), dimension(na_4_stw,4)                     :: j2_stw
  integer, intent(in)                                            :: i_stw
  integer, intent(in)                                            :: j_stw
  integer, intent(in)                                            :: m_stw
  integer, intent(in)                                            :: nfreq_stw
  real, intent(in), dimension(na_stw)                            :: anglz_stw
  real, intent(in), dimension(na_4_stw,4)                        :: wt1_stw
  real, intent(in), dimension(na_4_stw,4)                        :: wt2_stw
  real, intent(in)                                               :: dth_stw
  real, intent(in), dimension(na_4_stw,4)                        :: delr_stw
  real, intent(in), dimension(na_4_stw,4)                        :: delr_m_stw
  real, intent(in), dimension(na_4_stw,4)                        :: delr_p_stw
  real, intent(in), dimension(na_stw,nj_stw,ni_stw)              :: ddepdx_stw
  real, intent(in), dimension(na_stw,nj_stw,ni_stw)              :: ddepdy_stw
  real, intent(in), dimension(nfreq_stw,nj_stw,ni_stw)           :: wk_stw
  real, intent(in), dimension(nfreq_stw,nj_stw,ni_stw)           :: c_stwave
  real, intent(in), dimension(nfreq_stw,nj_stw,ni_stw)           :: cg_stwave
  real, intent(in), dimension(nfreq_stw)                         :: freq_stw
  real, intent(in)                                               :: epsd_stw
  real, intent(in)                                               :: pi2_stw
  real, intent(in)                                               :: twopi_stw
  real, intent(in)                                               :: epse_stw
  real, intent(inout), dimension(nfreq_stw,na_stw,nj_stw,ni_stw) :: e_stw
  real, intent(inout), dimension(nfreq_stw,na_stw,nj_stw,ni_stw) :: sea_stw

  integer :: n_stw, l_stw, jj1_stw, ii1_stw, jj2_stw, ii2_stw, k_stw
  integer :: lm_stw, lp_stw, ll_stw
  real :: dth2_stw, d1_stw, dddn_stw, dddn_m_stw, dddn_p_stw, tkd_stw
  real :: ang_turn, ang1_stw, ang1_m_stw, ang1_p_stw, del_ang_stw
  real :: e1_stw, wk1_stw, cg1_stw, c1_stw, tkd1_stw, dang_m_stw
  real :: dang_p_stw, fac_stw

  ! dth2_stw is half the direction bin (represents edge of the bin)
  dth2_stw = 0.5 * dth_stw

  do n_stw = 1, na_stw/4
    l_stw = l_sweep_stw(n_stw,m_stw)
    jj1_stw = j_stw+j1_stw(n_stw,m_stw)
    ii1_stw = i_stw+i1_stw(n_stw,m_stw)
    jj2_stw = j_stw+j2_stw(n_stw,m_stw)
    ii2_stw = i_stw+i2_stw(n_stw,m_stw)
    if(wt2_stw(n_stw,m_stw) .le. 0.00001) then
      jj2_stw = jj1_stw
      ii2_stw = ii1_stw
    elseif(wt1_stw(n_stw,m_stw) .le. 0.00001) then
      jj1_stw = jj2_stw
      ii1_stw = ii2_stw
    end if
    d1_stw = wt1_stw(n_stw,m_stw) * dep_stwave(jj1_stw,ii1_stw) + wt2_stw(n_stw,m_stw) * dep_stwave(jj2_stw,ii2_stw)
    if(d1_stw .le. epsd_stw .or. dep_stwave(j_stw,i_stw) .le. epsd_stw) then
      do k_stw = 1, nfreq_stw
        e_stw(k_stw,l_stw,j_stw,i_stw) = 0.0e0
      end do
    else  
      dddn_stw = -ddepdx_stw(l_stw,j_stw,i_stw)*sin(anglz_stw(l_stw)) +ddepdy_stw(l_stw,j_stw,i_stw)*cos(anglz_stw(l_stw))
      dddn_m_stw = -ddepdx_stw(l_stw,j_stw,i_stw)*sin(anglz_stw(l_stw)-dth2_stw) + ddepdy_stw(l_stw,j_stw,i_stw)*cos(anglz_stw(l_stw)-dth2_stw)
      dddn_p_stw = -ddepdx_stw(l_stw,j_stw,i_stw)*sin(anglz_stw(l_stw)+dth2_stw) + ddepdy_stw(l_stw,j_stw,i_stw)*cos(anglz_stw(l_stw)+dth2_stw)
      do k_stw = 1, nfreq_stw
        tkd_stw = min(2.0e0*wk_stw(k_stw,j_stw,i_stw)*dep_stwave(j_stw,i_stw), 20.0e0)
        ang_turn = delr_stw(n_stw,m_stw)/cg_stwave(k_stw,j_stw,i_stw)*(-c_stwave(k_stw,j_stw,i_stw)*wk_stw(k_stw,j_stw,i_stw)/sinh(tkd_stw))*dddn_stw
        ! for cases where the bathymetry is sufficiently resolved (angle turning greater than 90 deg), zero energy
        if(abs(ang_turn) .gt. pi2_stw) then
          e_stw(k_stw,l_stw,j_stw,i_stw) = 0.0e0
          sea_stw(k_stw,l_stw,j_stw,i_stw) = 0.0e0
        else
          ang1_stw = anglz_stw(l_stw) - ang_turn
          ang1_p_stw = anglz_stw(l_stw) + dth2_stw - (delr_p_stw(n_stw,m_stw)/cg_stwave(k_stw,j_stw,i_stw))*(-c_stwave(k_stw,j_stw,i_stw)*wk_stw(k_stw,j_stw,i_stw)/sinh(tkd_stw)*dddn_p_stw)
          ang1_m_stw = anglz_stw(l_stw) - dth2_stw - (delr_m_stw(n_stw,m_stw)/cg_stwave(k_stw,j_stw,i_stw))*(-c_stwave(k_stw,j_stw,i_stw)*wk_stw(k_stw,j_stw,i_stw)/sinh(tkd_stw)*dddn_m_stw)
          del_ang_stw = ang1_p_stw - ang1_m_stw
          if(del_ang_stw .le. 0.0 .or. del_ang_stw .gt. pi) then
            e_stw(k_stw,l_stw,j_stw,i_stw)=0.0
            sea_stw(k_stw,l_stw,j_stw,i_stw)=0.0
          else
            lm_stw = int((ang1_m_stw - dth2_stw)/dth_stw + 2)
            lp_stw = int((ang1_p_stw - dth2_stw)/dth_stw + 2)
            if(lm_stw .lt. 1) then
              lm_stw = lm_stw + na_stw
            end if
            if(lm_stw .gt. na_stw) then
             lm_stw = lm_stw - na_stw
            end if
            if(lp_stw .lt. 1) then
              lp_stw = lp_stw + na_stw
            end if
            if(lp_stw .gt. na_stw) then
              lp_stw = lp_stw - na_stw
            end if
            dang_m_stw=(asin(sin(anglz_stw(lm_stw))*cos(ang1_m_stw) - cos(anglz_stw(lm_stw))*sin(ang1_m_stw)) -dth2_stw)/dth_stw
            dang_p_stw=(asin(sin(anglz_stw(lp_stw))*cos(ang1_p_stw) - cos(anglz_stw(lp_stw))*sin(ang1_p_stw)) +dth2_stw)/dth_stw
            e1_stw = wt1_stw(n_stw,m_stw)*((e_stw(k_stw,lm_stw,jj1_stw,ii1_stw))*(dang_m_stw) - (e_stw(k_stw,lp_stw,jj1_stw,ii1_stw))*(dang_p_stw))+ wt2_stw(n_stw,m_stw)*((e_stw(k_stw,lm_stw,jj2_stw,ii2_stw))*(dang_m_stw)-(e_stw(k_stw,lp_stw,jj2_stw,ii2_stw))*(dang_p_stw))+ &
                 wt1_stw(n_stw,m_stw)*((sea_stw(k_stw,lm_stw,jj1_stw,ii1_stw))*(dang_m_stw)-(sea_stw(k_stw,lp_stw,jj1_stw,ii1_stw))*(dang_p_stw))+ wt2_stw(n_stw,m_stw)*((sea_stw(k_stw,lm_stw,jj2_stw,ii2_stw))*(dang_m_stw)-(sea_stw(k_stw,lp_stw,jj2_stw,ii2_stw))*(dang_p_stw))
            if(lm_stw .le. lp_stw) then
              do ll_stw = lm_stw,lp_stw
                e1_stw = e1_stw + wt1_stw(n_stw,m_stw)*(e_stw(k_stw,ll_stw,jj1_stw,ii1_stw)) + wt2_stw(n_stw,m_stw)*(e_stw(k_stw,ll_stw,jj2_stw,ii2_stw)) + wt1_stw(n_stw,m_stw)*(sea_stw(k_stw,ll_stw,jj1_stw,ii1_stw)) + wt2_stw(n_stw,m_stw)*(sea_stw(k_stw,ll_stw,jj2_stw,ii2_stw))
              end do
            else
              do ll_stw = lm_stw,na_stw
                e1_stw = e1_stw + wt1_stw(n_stw,m_stw)*(e_stw(k_stw,ll_stw,jj1_stw,ii1_stw)) + wt2_stw(n_stw,m_stw)*(e_stw(k_stw,ll_stw,jj2_stw,ii2_stw)) + wt1_stw(n_stw,m_stw)*(sea_stw(k_stw,ll_stw,jj1_stw,ii1_stw)) + wt2_stw(n_stw,m_stw)*(sea_stw(k_stw,ll_stw,jj2_stw,ii2_stw))
              end do
              do ll_stw = 1,lp_stw
                e1_stw = e1_stw + wt1_stw(n_stw,m_stw)*(e_stw(k_stw,ll_stw,jj1_stw,ii1_stw)) + wt2_stw(n_stw,m_stw)*(e_stw(k_stw,ll_stw,jj2_stw,ii2_stw)) + wt1_stw(n_stw,m_stw)*(sea_stw(k_stw,ll_stw,jj1_stw,ii1_stw)) + wt2_stw(n_stw,m_stw)*(sea_stw(k_stw,ll_stw,jj2_stw,ii2_stw))
              end do
            end if
            wk1_stw = wkfnc_stw(freq_stw(k_stw),d1_stw)
            c1_stw = freq_stw(k_stw) * twopi_stw / wk1_stw
            tkd1_stw = min (2.0e0 * wk1_stw * d1_stw, 20.0e0)
            cg1_stw = 0.5e0 * c1_stw * (1.0e0 + (2.0e0 * wk1_stw * d1_stw / sinh(tkd1_stw)))
            e_stw(k_stw,l_stw,j_stw,i_stw) = e1_stw * c1_stw * cg1_stw / (c_stwave(k_stw,j_stw,i_stw) *cg_stwave(k_stw,j_stw,i_stw)) * (dth_stw / del_ang_stw)
            fac_stw = c1_stw * cg1_stw / (c_stwave(k_stw,j_stw,i_stw) *cg_stwave(k_stw,j_stw,i_stw)) * (dth_stw / del_ang_stw)
            if (e_stw(k_stw,l_stw,j_stw,i_stw) .le. epse_stw) then 
              e_stw(k_stw,l_stw,j_stw,i_stw) = 0.0e0
            end if 
          end if
        end if
      end do
    end if
  end do

END SUBROUTINE prop

! ********************************************************************
SUBROUTINE one_d_prop (nj_stw,ni_stw,na_4_stw,anglz_stw,l_sweep_stw,e_stw,sea_stw, &
             delr_stw,delr_m_stw,delr_p_stw,wk_stw,c_stwave,cg_stwave,freq_stw,i_stw,j_stw,i_in,j_in,m_stw,ifric_stw,cf_stw,dth_stw, epsd_stw, nfreq_stw,na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw)
! ********************************************************************
!
! Purpose: 
!
! --------------------------------------------------------------------
  ! Subroutine arguments
  integer, intent(in)                                            :: nj_stw
  integer, intent(in)                                            :: ni_stw
  integer, intent(in)                                            :: na_4_stw
  integer, intent(in)                                            :: i_stw
  integer, intent(in)                                            :: j_stw
  integer, intent(in)                                            :: i_in
  integer, intent(in)                                            :: j_in
  integer, intent(in)                                            :: m_stw
  integer, intent(in)                                            :: ifric_stw
  integer, intent(in)                                            :: nfreq_stw
  integer, intent(in)                                            :: na_stw
  integer, intent(in), dimension(na_4_stw,4)                     :: l_sweep_stw
  real, intent(in), dimension(na_stw)                            :: anglz_stw
  real, intent(in), dimension(na_4_stw,4)                        :: delr_stw
  real, intent(in), dimension(na_4_stw,4)                        :: delr_m_stw
  real, intent(in), dimension(na_4_stw,4)                        :: delr_p_stw
  real, intent(in), dimension(nfreq_stw,nj_stw,ni_stw)           :: wk_stw
  real, intent(in), dimension(nfreq_stw,nj_stw,ni_stw)           :: c_stwave
  real, intent(in), dimension(nfreq_stw,nj_stw,ni_stw)           :: cg_stwave
  real, intent(in), dimension(nfreq_stw)                         :: freq_stw
  real, intent(in), dimension(nj_stw,ni_stw)                     :: cf_stw
  real, intent(in)                                               :: dth_stw
  real, intent(in)                                               :: epsd_stw
  real, intent(in)                                               :: twopi_stw
  real, intent(in)                                               :: pi2_stw
  real, intent(in)                                               :: radfac_stw  
  real, intent(in)                                               :: epse_stw
  real, intent(inout), dimension(nfreq_stw,na_stw,nj_stw,ni_stw) :: e_stw
  real, intent(inout), dimension(nfreq_stw,na_stw,nj_stw,ni_stw) :: sea_stw

  integer :: l_stw, n_stw, k_stw, lm_stw, lp_stw, dang_m_stw, dang_p_stw, ll_stw
  real :: dth2_stw, ddepdx_stw, ddepdy_stw, dddn_stw, dddn_m_stw, dddn_p_stw
  real :: tkd_stw, ang_turn, ang1_stw, ang1_m_stw, ang1_p_stw, del_ang_stw
  real :: e1_stw, b_stw, sigma_stw

! dth2_stw is half the direction bin (represents edge of the bin)
  dth2_stw = 0.5 * dth_stw

  if(j_stw.eq.j_in)then
    ddepdx_stw = (dep_stwave(j_stw,i_stw) - dep_stwave(j_in,i_in))/((i_stw-i_in)*dx)
    ddepdy_stw = 0.0
  else
    ddepdx_stw = 0.0
    ddepdy_stw = (dep_stwave(j_stw,i_stw) - dep_stwave(j_in,i_in))/((j_stw-j_in)*dy)
  endif

  do n_stw = 1, na_stw/4
    l_stw = l_sweep_stw(n_stw,m_stw)
    if(dep_stwave(j_in,i_in) .le. epsd_stw .or. dep_stwave(j_stw,i_stw) .le. epsd_stw) then
      do k_stw = 1, nfreq_stw
        e_stw(k_stw,l_stw,j_stw,i_stw) = 0.0e0
        sea_stw(k_stw,l_stw,j_stw,i_stw) = 0.0e0
      end do
    else
      dddn_stw = -ddepdx_stw*sin(anglz_stw(l_stw)) +ddepdy_stw*cos(anglz_stw(l_stw))
      dddn_m_stw = -ddepdx_stw*sin(anglz_stw(l_stw)-dth2_stw) + ddepdy_stw*cos(anglz_stw(l_stw)-dth2_stw)
      dddn_p_stw = -ddepdx_stw*sin(anglz_stw(l_stw)+dth2_stw) + ddepdy_stw*cos(anglz_stw(l_stw)+dth2_stw)
      do k_stw = 1, nfreq_stw
        tkd_stw = min(2.0e0*wk_stw(k_stw,j_stw,i_stw)*dep_stwave(j_stw,i_stw), 20.0e0)
        ang_turn = delr_stw(n_stw,m_stw)/cg_stwave(k_stw,j_stw,i_stw)*(-c_stwave(k_stw,j_stw,i_stw)*wk_stw(k_stw,j_stw,i_stw)/sinh(tkd_stw))*dddn_stw
        ! for cases where the bathymetry is sufficiently resolved (angle turning greater than 90 deg), zero energy
        if(abs(ang_turn) .gt. pi2_stw) then
          e_stw(k_stw,l_stw,j_stw,i_stw) = 0.0e0
        else
          ang1_stw = anglz_stw(l_stw) - ang_turn
          ang1_p_stw = anglz_stw(l_stw) + dth2_stw - (delr_p_stw(n_stw,m_stw)/cg_stwave(k_stw,j_stw,i_stw))*(-c_stwave(k_stw,j_stw,i_stw)*wk_stw(k_stw,j_stw,i_stw)/sinh(tkd_stw)*dddn_p_stw)
          ang1_m_stw = anglz_stw(l_stw) - dth2_stw - (delr_m_stw(n_stw,m_stw)/cg_stwave(k_stw,j_stw,i_stw))*(-c_stwave(k_stw,j_stw,i_stw)*wk_stw(k_stw,j_stw,i_stw)/sinh(tkd_stw)*dddn_m_stw)
          del_ang_stw = ang1_p_stw - ang1_m_stw
          if(del_ang_stw .le. 0.0) then
            e_stw(k_stw,l_stw,j_stw,i_stw)=0.0
          else
            lm_stw = int((ang1_m_stw - dth2_stw)/dth_stw + 2)
            lp_stw = int((ang1_p_stw - dth2_stw)/dth_stw + 2)
            if(lm_stw .lt. 1) then
              lm_stw = lm_stw + na_stw
            end if
            if(lm_stw .gt. na_stw) then
              lm_stw = lm_stw - na_stw
            end if
            if(lp_stw.lt.1) then
              lp_stw = lp_stw + na_stw
            end if
            if(lp_stw.gt.na_stw) then
              lp_stw = lp_stw - na_stw
            end if
            dang_m_stw=(asin(sin(anglz_stw(lm_stw))*cos(ang1_m_stw)-cos(anglz_stw(lm_stw))*sin(ang1_m_stw)) - dth2_stw)/dth_stw
            dang_p_stw=(asin(sin(anglz_stw(lp_stw))*cos(ang1_p_stw)-cos(anglz_stw(lp_stw))*sin(ang1_p_stw)) + dth2_stw)/dth_stw
            e1_stw = ((e_stw(k_stw,lm_stw,j_in,i_in))*(dang_m_stw) - (e_stw(k_stw,lp_stw,j_in,i_in))*(dang_p_stw))
            if(lm_stw .le. lp_stw) then
              do ll_stw = lm_stw,lp_stw
                e1_stw = e1_stw + (e_stw(k_stw,ll_stw,j_in,i_in))
              end do
            else
              do ll_stw = lm_stw,na_stw
                e1_stw = e1_stw + (e_stw(k_stw,ll_stw,j_in,i_in))                
              end do
              do ll_stw = 1,lp_stw
                e1_stw = e1_stw + (e_stw(k_stw,ll_stw,j_in,i_in))
              end do
            end if
            e_stw(k_stw,l_stw,j_stw,i_stw) = e1_stw * c_stwave(k_stw,j_in,i_in) * cg_stwave(k_stw,j_in,i_in) / (c_stwave(k_stw,j_stw,i_stw) *cg_stwave(k_stw,j_stw,i_stw)) * (dth_stw / del_ang_stw)
            if (e_stw(k_stw,l_stw,j_stw,i_stw) .le. epse_stw) then
              e_stw(k_stw,l_stw,j_stw,i_stw) = 0.0e0
            end if 
            if(ifric_stw .ge. 1) then
              sigma_stw = twopi_stw*freq_stw(k_stw)
              if(ifric_stw .le. 2) then
                b_stw = -cf_stw(j_stw,i_stw)/g*(sigma_stw/sinh(tkd_stw))**2
              else
                b_stw = -cf_stw(j_stw,i_stw)*cf_stw(j_stw,i_stw)*dep_stwave(j_stw,i_stw)**(-0.3333)*(sigma_stw/sinh(tkd_stw))**2
              end if
              e_stw(k_stw, l_stw, j_stw, i_stw) = e_stw(k_stw, l_stw, j_stw, i_stw)*exp(b_stw*delr_stw(n_stw,m_stw)/cg_stwave(k_stw,j_stw,i_stw))
              ! sea_stw(k_stw, l_stw, j_stw, i_stw) = sea_stw(k_stw, l_stw, j_stw, i_stw)*exp(b_stw*delr_stw(n_stw,m_stw)/Cg(k_stw,j_stw,i_stw))
            end if
          end if
        end if
      end do
    end if
  end do

END SUBROUTINE one_d_prop

! ********************************************************************
SUBROUTINE break(ni_stw,nj_stw,i_stw,j_stw,e_stw,sea_stw,delf_stw,wk_stw,ibr_stw,kmax_stw, na_stw,twopi_stw,pi2_stw,radfac_stw, epse_stw, nfreq_stw, dth_stw)
! ********************************************************************
!
! Purpose: 
!
! --------------------------------------------------------------------
  ! Subroutine arguments
  integer, intent(in)                                            :: nj_stw
  integer, intent(in)                                            :: ni_stw
  integer, intent(in)                                            :: i_stw
  integer, intent(in)                                            :: j_stw
  integer, intent(in)                                            :: na_stw
  integer, intent(in)                                            :: nfreq_stw
  integer, intent(inout), dimension(nj_stw,ni_stw)               :: ibr_stw
  integer, intent(inout)                                         :: kmax_stw
  real, intent(in), dimension(nfreq_stw)                         :: delf_stw
  real, intent(in), dimension(nfreq_stw,nj_stw,ni_stw)           :: wk_stw
  real, intent(in)                                               :: twopi_stw
  real, intent(in)                                               :: pi2_stw
  real, intent(in)                                               :: radfac_stw  
  real, intent(in)                                               :: epse_stw
  real, intent(in)                                               :: dth_stw
  real, intent(inout), dimension(nfreq_stw,na_stw,nj_stw,ni_stw) :: e_stw
  real, intent(inout), dimension(nfreq_stw,na_stw,nj_stw,ni_stw) :: sea_stw
  
  integer :: k_stw, l_stw
  real :: sum_e_stw, emax_stw, sum_ef_stw, edplim_stw, wkpeak_stw, estlim_stw, ration_stw

! Check for breaking.

  sum_e_stw = 0.0e0
  kmax_stw = nfreq_stw
  emax_stw = 0.0e0
  do k_stw = 1, nfreq_stw
    sum_ef_stw = 0.0e0
    do l_stw = 1, na_stw
      sum_e_stw = sum_e_stw + (e_stw(k_stw,l_stw,j_stw,i_stw) + sea_stw(k_stw,l_stw,j_stw,i_stw)) * delf_stw(k_stw) * dth_stw
      sum_ef_stw = sum_ef_stw + (e_stw(k_stw,l_stw,j_stw,i_stw) + sea_stw(k_stw,l_stw,j_stw,i_stw)) *  delf_stw(k_stw) * dth_stw
    end do
    if (sum_ef_stw .ge. emax_stw) then
      kmax_stw = k_stw
      emax_stw = sum_ef_stw
    end if
  end do
  if (dep_stwave(j_stw,i_stw) .gt. 0.0e0 .and. sum_e_stw .gt. 0.0e0) then
    edplim_stw = 0.64e0 * dep_stwave(j_stw,i_stw)
    wkpeak_stw = wk_stw(kmax_stw,j_stw,i_stw)
    estlim_stw = 0.2e0 * pi / wkpeak_stw * tanh(wkpeak_stw * dep_stwave(j_stw,i_stw))
    edplim_stw = min (edplim_stw, estlim_stw)
    if(sum_e_stw.le.0.0e0) then
      sum_e_stw=0.0e0
    end if
    sum_e_stw = 4.0e0 * sqrt(sum_e_stw)
    if (sum_e_stw .gt. edplim_stw) then
      ibr_stw(j_stw,i_stw) = 1
      ration_stw = (edplim_stw / (sum_e_stw + 1.0e-6)) **2
      do k_stw = 1, nfreq_stw
        do l_stw = 1, na_stw
          e_stw(k_stw,l_stw,j_stw,i_stw) = ration_stw * e_stw(k_stw,l_stw,j_stw,i_stw)
          sea_stw(k_stw,l_stw,j_stw,i_stw) = ration_stw * sea_stw(k_stw,l_stw,j_stw,i_stw)
        end do
      end do
    end if
  end if

END SUBROUTINE break

! ********************************************************************
SUBROUTINE celerity(nfreq_stw,ni_stw,nj_stw,wk_stw,c_stwave,cg_stwave,freq_stw)
! ********************************************************************
!
! Purpose: 
!
! --------------------------------------------------------------------
  ! Subroutine arguments
  real, intent(inout), dimension(nfreq_stw,nj_stw,ni_stw) :: c_stwave
  real, intent(inout), dimension(nfreq_stw,nj_stw,ni_stw) :: wk_stw
  real, intent(inout), dimension(nfreq_stw,nj_stw,ni_stw) :: cg_stwave
  real, intent(in), dimension(nfreq_stw)                  :: freq_stw
  integer, intent(in)                                     :: ni_stw
  integer, intent(in)                                     :: nj_stw
  integer, intent(in)                                     :: nfreq_stw
  integer :: i_stw, j_stw, k_stw
  real :: tkd_stw
  real :: twopi_stw = 2.0e0*pi

  ! calculate wave c_stwave for entire grid for all frequencies

  do k_stw = 1, nfreq_stw
    do j_stw = 1,nj_stw
      do i_stw = 1,ni_stw
        if (dep_stwave(j_stw,i_stw) .gt. 0.0e0) then
          wk_stw(k_stw,j_stw,i_stw) = wkfnc_stw(freq_stw(k_stw),dep_stwave(j_stw,i_stw))
          c_stwave(k_stw,j_stw,i_stw) = freq_stw(k_stw)*twopi_stw/wk_stw(k_stw,j_stw,i_stw)
          tkd_stw = min(2.0e0*wk_stw(k_stw,j_stw,i_stw)*(dep_stwave(j_stw,i_stw)), 20.0e0)
          cg_stwave(k_stw,j_stw,i_stw) = 0.5e0 * c_stwave(k_stw,j_stw,i_stw) * (1.0 + (2.0 * wk_stw(k_stw,j_stw,i_stw) * dep_stwave(j_stw,i_stw) / sinh(tkd_stw)))
        else
          wk_stw(k_stw,j_stw,i_stw) = 0.0e0
          c_stwave(k_stw,j_stw,i_stw) = 0.0e0
          cg_stwave(k_stw,j_stw,i_stw) = 0.0e0
        end if
      end do
    end do
  end do

END SUBROUTINE celerity

! ********************************************************************
REAL FUNCTION wkfnc_stw(fm_stw, d_stw)
! ********************************************************************
!
! Purpose: 
!
! --------------------------------------------------------------------
  ! Subroutine arguments
  real             :: y
  real, intent(in) :: fm_stw
  real, intent(in) :: d_stw
  real             :: sum_stw
  real             :: coe1_stw = 0.66667e0
  real             :: coe2_stw = 0.35550e0
  real             :: coe3_stw = 0.16084e0
  real             :: coe4_stw = 0.06320e0
  real             :: coe5_stw = 0.02174e0
  real             :: coe6_stw = 0.00654e0
  real             :: coe7_stw = 0.00171e0
  real             :: coe8_stw = 0.00039e0
  real             :: coe9_stw = 0.00011e0
  real             :: twopi_stw = 2.0e0*pi

  y = (twopi_stw*twopi_stw*fm_stw*fm_stw)*d_stw/g

  !  SERIES APPROXIMATION:
  sum_stw = 1.0e0 + y*(coe1_stw + y*(coe2_stw + y*(coe3_stw + y*(coe4_stw + y*(coe5_stw &
  + y*(coe6_stw + y*(coe7_stw + y*(coe8_stw + y*coe9_stw))))))))
  sum_stw = 1.0e0/(y + 1.0e0/sum_stw)
  sum_stw = sqrt(g*d_stw*sum_stw)
  wkfnc_stw = twopi_stw*fm_stw/sum_stw

  return

END FUNCTION wkfnc_stw
! ********************************************************************
SUBROUTINE depgrad(nj_stw,ni_stw,na_stw,na_4_stw,l_sweep_stw,i1_stw,i2_stw,j1_stw,j2_stw,anglz_stw,delr_stw,ddepdx_stw,ddepdy_stw)
! ********************************************************************
!
! Purpose: 
!
! --------------------------------------------------------------------
  ! Subroutine arguments
  integer, intent(in)                                  :: nj_stw
  integer, intent(in)                                  :: ni_stw
  integer, intent(in)                                  :: na_4_stw
  integer, intent(in)                                  :: na_stw
  integer, intent(in), dimension(na_4_stw,4)           :: i1_stw
  integer, intent(in), dimension(na_4_stw,4)           :: i2_stw
  integer, intent(in), dimension(na_4_stw,4)           :: j1_stw
  integer, intent(in), dimension(na_4_stw,4)           :: j2_stw
  integer, intent(in), dimension(na_4_stw,4)           :: l_sweep_stw
  real, intent(inout), dimension(na_stw,nj_stw,ni_stw) :: ddepdx_stw
  real, intent(inout), dimension(na_stw,nj_stw,ni_stw) :: ddepdy_stw
  real, intent(in), dimension(na_stw)                  :: anglz_stw
  real, intent(in), dimension(na_4_stw,4)              :: delr_stw

  integer :: m_stw, n_stw, l_stw, i_stw, j_stw, ix1_stw, ix2_stw
  integer :: jx1_stw, jx2_stw, jy1_stw, jy2_stw
  real :: atan_ang_stw, delx_stw, dely_stw

  do m_stw = 1, 4
    do n_stw = 1, na_4_stw
      l_stw = l_sweep_stw(n_stw,m_stw)
      atan_ang_stw = abs(tan(anglz_stw(l_stw)))
      if (atan_ang_stw .gt. 4.0e0) then
        ix1_stw = -1
        ix2_stw = 1
        delx_stw = 2.0e0*dx
      elseif (i1_stw(n_stw,m_stw) .lt. 0) then
        ix1_stw = -1
        ix2_stw = 0
        delx_stw = dx
      else
        ix1_stw = 0
        ix2_stw = 1
        delx_stw = dx
      endif
      if (atan_ang_stw .lt. 0.25e0) then
        jy1_stw = -1
        jy2_stw = 1
        dely_stw = 2.0e0*dy
      elseif (j1_stw(n_stw,m_stw) .lt. 0) then
        jy1_stw = -1
        jy2_stw = 0
        dely_stw = dy
      else
        jy1_stw = 0
        jy2_stw = 1
        dely_stw = dy
      endif
      do i_stw = 2, ni_stw-1
        do j_stw = 2, nj_stw-1
          ddepdx_stw(l_stw,j_stw,i_stw) = (dep_stwave(j_stw,i_stw+ix2_stw)-dep_stwave(j_stw,i_stw+ix1_stw))/delx_stw
          ddepdy_stw(l_stw,j_stw,i_stw) = (dep_stwave(j_stw+jy2_stw,i_stw)-dep_stwave(j_stw+jy1_stw,i_stw))/dely_stw
        enddo
      enddo
      do j_stw = 2, nj_stw-1
        ddepdx_stw(l_stw,j_stw,1) = ddepdx_stw(l_stw,j_stw,2)
        ddepdy_stw(l_stw,j_stw,1) = ddepdy_stw(l_stw,j_stw,2)
      enddo
      do j_stw = 2, nj_stw-1
        ddepdx_stw(l_stw,j_stw,ni_stw) = ddepdx_stw(l_stw,j_stw,ni_stw-1)
        ddepdy_stw(l_stw,j_stw,ni_stw) = ddepdy_stw(l_stw,j_stw,ni_stw-1)
      enddo
      do i_stw = 2, ni_stw-1
        ddepdx_stw(l_stw,1,i_stw) = ddepdx_stw(l_stw,2,i_stw)
        ddepdy_stw(l_stw,1,i_stw) = ddepdy_stw(l_stw,2,i_stw)
      enddo
      do i_stw = 2, ni_stw-1
        ddepdx_stw(l_stw,nj_stw,i_stw) = ddepdx_stw(l_stw,nj_stw-1,i_stw)
        ddepdy_stw(l_stw,nj_stw,i_stw) = ddepdy_stw(l_stw,nj_stw-1,i_stw)
      enddo 
    enddo
  enddo

END SUBROUTINE depgrad

! ********************************************************************
SUBROUTINE sweeps(na_4_stw,na_stw,dth_stw,nj_stw,ni_stw,anglz_stw,l_sweep_stw,i1_stw,i2_stw,j1_stw,j2_stw,wt1_stw,wt2_stw,delr_stw,delr_m_stw,delr_p_stw)
! ********************************************************************
!
! Purpose: to determine order to sweep grid so the secto with the
!          wind input is swept last
!
! --------------------------------------------------------------------
  ! Subroutine arguments
  integer, intent(in)                           :: na_4_stw
  integer, intent(in)                           :: na_stw
  integer, intent(in)                           :: ni_stw
  integer, intent(in)                           :: nj_stw
  integer, intent(inout), dimension(na_4_stw,4) :: i1_stw
  integer, intent(inout), dimension(na_4_stw,4) :: i2_stw
  integer, intent(inout), dimension(na_4_stw,4) :: j1_stw
  integer, intent(inout), dimension(na_4_stw,4) :: j2_stw
  integer, intent(inout), dimension(na_4_stw,4) :: l_sweep_stw
  real, intent(inout), dimension(na_4_stw,4)    :: wt1_stw
  real, intent(inout), dimension(na_4_stw,4)    :: wt2_stw
  real, intent(inout), dimension(na_4_stw,4)    :: delr_stw
  real, intent(inout), dimension(na_4_stw,4)    :: delr_m_stw
  real, intent(inout), dimension(na_4_stw,4)    :: delr_p_stw
  real, intent(inout), dimension(na_stw)        :: anglz_stw
  real, intent(in)                              :: dth_stw

  integer :: m_stw, n_stw, l_stw
  real :: atan_ang_stw

  do m_stw = 1, 4
    l_sweep_stw(1,m_stw) = 1 + (m_stw-1) * na_4_stw
    do n_stw = 2, na_4_stw
      l_sweep_stw(n_stw,m_stw) = l_sweep_stw(n_stw-1,m_stw) + 1
    enddo
  enddo

  ! determine points to interpolate for each angle band in forward sweep

  do m_stw = 1,4
    do n_stw = 1,na_4_stw
      l_stw = l_sweep_stw(n_stw,m_stw)
      atan_ang_stw = abs (tan(anglz_stw(l_stw)))
      if (atan_ang_stw / (dy / dx) .le. 1.0e0) then

        ! interpolate in y
        if (cos(anglz_stw(l_stw)) .ge. 0.0e0) then
          i1_stw(n_stw,m_stw) = -1
        else
          i1_stw(n_stw,m_stw) = 1
        endif
        i2_stw(n_stw,m_stw) = i1_stw(n_stw,m_stw)
        if (sin(anglz_stw(l_stw)) .ge. 0.0e0) then
          j1_stw(n_stw,m_stw) = -1
        else
          j1_stw(n_stw,m_stw) = 1
        endif
        j2_stw(n_stw,m_stw) = 0

        ! determine interpolation weights

        wt1_stw(n_stw,m_stw) = dx / dy * atan_ang_stw
        wt2_stw(n_stw,m_stw) = 1.0e0 - wt1_stw(n_stw,m_stw)
        delr_stw(n_stw,m_stw) = dx / abs(cos(anglz_stw(l_stw)))
        delr_p_stw(n_stw,m_stw) = dx / abs(cos(anglz_stw(l_stw)+0.5*dth_stw))
        delr_m_stw(n_stw,m_stw) = dx / abs(cos(anglz_stw(l_stw)-0.5*dth_stw))

      else

        ! interpolate in x
        if (sin(anglz_stw(l_stw)) .gt. 0.0e0)then
          j1_stw(n_stw,m_stw) = -1
        else
          j1_stw(n_stw,m_stw) = 1
        endif
        j2_stw(n_stw,m_stw) = j1_stw(n_stw,m_stw)
        if (cos(anglz_stw(l_stw)) .gt. 0.0e0) then
          i1_stw(n_stw,m_stw) = -1
        else
          i1_stw(n_stw,m_stw) = 1
        endif
          i2_stw(n_stw,m_stw) = 0

        ! determine interpolation weights
        wt1_stw(n_stw,m_stw) = dy / atan_ang_stw / dx
        wt2_stw(n_stw,m_stw) = 1.0e0 - wt1_stw(n_stw,m_stw)
        delr_stw(n_stw,m_stw) = dy / abs(sin(anglz_stw(l_stw)))
        delr_p_stw(n_stw,m_stw) = dy / abs(sin(anglz_stw(l_stw)+0.5*dth_stw))
        delr_m_stw(n_stw,m_stw) = dy / abs(sin(anglz_stw(l_stw)-0.5*dth_stw))
      endif
    enddo
  enddo

END SUBROUTINE sweeps

! ********************************************************************
SUBROUTINE setsweeps(i_st_stw,i_en_stw,i_inc_stw,j_st_stw,j_en_stw,j_inc_stw,ni_stw,nj_stw)
! ********************************************************************
!
! Purpose: to determine order to sweep grid so the secto with the
!          wind input is swept last
!
! --------------------------------------------------------------------
  ! Subroutine arguments
  integer, intent(inout), dimension(4)  :: i_st_stw
  integer, intent(inout), dimension(4)  :: i_en_stw
  integer, intent(inout), dimension(4)  :: i_inc_stw
  integer, intent(inout), dimension(4)  :: j_st_stw
  integer, intent(inout), dimension(4)  :: j_en_stw
  integer, intent(inout), dimension(4)  :: j_inc_stw
  integer, intent(in)                   :: ni_stw
  integer, intent(in)                   :: nj_stw

  i_st_stw(1) = 2
  i_st_stw(2) = ni_stw-1
  i_st_stw(3) = ni_stw-1
  i_st_stw(4) = 2
  i_en_stw(1) = ni_stw-1
  i_en_stw(2) = 2
  i_en_stw(3) = 2
  i_en_stw(4) = ni_stw-1
  i_inc_stw(1) = 1
  i_inc_stw(2) = -1
  i_inc_stw(3) = -1
  i_inc_stw(4) = 1
  j_st_stw(1) = 2
  j_st_stw(2) = 2
  j_st_stw(3) = nj_stw-1
  j_st_stw(4) = nj_stw-1
  j_en_stw(1) = nj_stw-1
  j_en_stw(2) = nj_stw-1
  j_en_stw(3) = 2
  j_en_stw(4) = 2
  j_inc_stw(1) = 1
  j_inc_stw(2) = 1
  j_inc_stw(3) = -1
  j_inc_stw(4) = -1

END SUBROUTINE setsweeps
!************************************************************************
                        END MODULE si3d_stwave
!************************************************************************