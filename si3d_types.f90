!************************************************************************
                          MODULE si3d_types
!************************************************************************
!
!  Purpose: Global type declarations for the semi-implicit 3-D (si3d)
!           hydrodynamic model.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-------------------------------------------------------------------------

   IMPLICIT NONE
   SAVE

   !                   -----Named Constants-----

   !.....I/O file numbers.....
   INTEGER, PARAMETER :: i4=4, i5=50, i6=60, i82=82, i83=83, i93=93,    &
                       & i94=94, i96=96, i97=97, i98=98, i99 = 99
   INTEGER, PARAMETER :: i54=i6-6, i55=i6-5, i56=i6-4, i57=i6-3,        &
                       & i58=i6-2, i59=i6-1, i53=i6-7, i52=i6-8

   !.....Define kind numbers for single- and double-precision reals.....
   INTEGER, PARAMETER :: SPV = KIND(1.0)
   INTEGER, PARAMETER :: DPV = KIND(1.0D0)

   !.....Define kind number for half-precision (2-byte) integers.....
   INTEGER, PARAMETER :: INT2 = SELECTED_INT_KIND(3)

   !.....Geophysical constants.....
   REAL(SPV), PARAMETER :: g = 9.806, rhoair = 1.3, kappaS = .41
   REAL(SPV), PARAMETER :: pi = 3.1415926535897932_spv
   REAL     , PARAMETER :: HMIN = 0.0E0, ZERO = 0.0E0, SMALL = 1E-1

   !.....Initial condition constants.....
   REAL, PARAMETER :: u0 = 0.0, u0h0 = 0.0, v0 = 0.0, v0h0 = 0.0, w0 = 0.0

   !.....Maximum number of nodes where time series output is desired.....
   INTEGER, PARAMETER :: maxnodes = 20

   !.....Maximum number of open boundaries.....
   INTEGER, PARAMETER :: maxnopen = 10

   !.....Harmonic constants for tidal test problem.....
   REAL, PARAMETER :: h0 = 0.0, amp1= .9144, period1 = 12.4
   REAL, PARAMETER :: coef1 = 2.*pi/period1

   !.....Harmonic constants for salinity open boundaries.....
   REAL, PARAMETER :: sal00 = 0.0, amp2 = 5.0, period2 = 12.4
   REAL, PARAMETER :: coef2 = 2.*pi/period2

   !                      -----Variables-----

   !.....Program timing variables.....
   REAL    :: t_exmom=0, t_matmom=0, t_matcon=0, t_solver=0, t_vel=0,    &
            & t_exsal=0, t_salin=0, t_outt=0, t_trid=0, t_turb=0,        &
            & t_settrap=0, t_save=0, t_setmask=0, tot_subs=0
   INTEGER :: n_trid=0, n_turb=0, n_exmom=0

   !.....Input variables.....
   REAL    :: xl, yl, zl, tl, cd, cw, wa, phi, f, amp, cstar, ax0, ay0,  &
            & Av0, Dv0, Av00, a1, b1, a2, b2, alp, beta, theta, dzmin,   &
            & datadj, t0, sal0, tramp, begind, chi, zetainit
   REAL    :: idz, idt, idx, idy
   INTEGER :: iexplt, ivde, itrap, niter, ismooth,        &
            & iturb, ilin, ihomo, ipt, ipx, ipc, iadv, ihd, ibc, isal,   &
            & nnodes, nopen, iseich, idbg, iextrp, iyr0, imon0, iday0,   &
            & ihr0, istd, igs, ivpv, iupwind, ioutg, ipv, ipsal, ipxml,  &
            & itspf, itrsca, itrmom, ibathyf, apxml
   CHARACTER :: title*80, sal_ic_file*50, wse_file*7, flw_file*7,        &
              & hcn_file*7, sal_file*7, barrier_file*50, commentline*50

   !.....Frequently used numerical constants and coefficients.....
   REAL    :: dt, dx, dy, ddz, twodt, dtdx, dtdy, gdtdx, gdtdy, gdt2dx2, &
            & gdt2dy2, gthx, gthy, gth1x, gth1y, cwind, alp4, drho,      &
            & dsal, rho0, twodx, twody, twodxdx, twodydy, dt_min,        &
            & dxdx, dydy, fourdx, fourdy, beta2, chi1, twochi1, dxdy
   INTEGER :: im, im1, i1, jm, jm1, j1, km, km1, k1, ndx, ndy, ndz, nts, &
            & ndim, maxnz, nw, inw, ibdwd, ifirst, ilast, jfirst, jlast, &
            & lfirst,llast

   !.....Frequently used scalar variables.....
   REAL    :: thrs, tz, its
   INTEGER :: n, istep, lastiter

   !.....Open boundary variables.....
   INTEGER :: nwse_t, nwse_hc, nflw, nsal

   !.....Date variables.....
   INTEGER :: iyr, imon, iday, ihr
   REAL :: isec0,isec
   REAL    :: doy, doyp

   !           -----Explicit-Shape Array Declarations-----

   !.....Node numbers arrays where time series output is desired.....
   INTEGER, DIMENSION(maxnodes) :: inode, jnode

   !         -----Permanent Allocatable Array Declarations-----

   !.....Arrays thread.....
   INTEGER :: num_threads
   INTEGER, ALLOCATABLE, DIMENSION(:) :: lh, lh_aux, lhi, lhf, ph,iauxs,iauxe
   INTEGER, ALLOCATABLE, DIMENSION(:) :: id_column, lhiE,lhfE,lhiW,lhfW
   INTEGER, ALLOCATABLE, DIMENSION(:) :: id_columnCE, lhiECE,lhfECE,lhiWCE,lhfWCE
   INTEGER, ALLOCATABLE, DIMENSION(:) :: lhiCE,lhfCE,lhiCN,lhfCN,flag,nopth
   INTEGER, ALLOCATABLE, DIMENSION(:) :: id_columnCN, lhiECN,lhfECN,lhiWCN,lhfWCN,nopenH,nopenHH
   INTEGER, ALLOCATABLE, DIMENSION (:,:) :: isbcH, jsbcH, iebcH, jebcH,noh2no,noh2noH
   INTEGER, ALLOCATABLE, DIMENSION (:,:) :: siptNBI,eiptNBI,siptNBIH,eiptNBIH,isbcHH, iebcHH
   REAL, ALLOCATABLE, DIMENSION (:) :: areatot,cont,contNG


   !.....Arrays at u-pts.....
   REAL, ALLOCATABLE, DIMENSION(:,:) :: uh, uhp, uhpp, u, up, upp
   REAL, ALLOCATABLE, DIMENSION(:,:) :: hu, hup, hupp
   REAL, ALLOCATABLE, DIMENSION(:,:) :: ex, agx, arx, ex2, extr
   REAL, ALLOCATABLE, DIMENSION(:) :: eagx, earx

   !.....Arrays at v-pts.....
   REAL, ALLOCATABLE, DIMENSION(:,:) :: vh, vhp, vhpp, v, vp, vpp
   REAL, ALLOCATABLE, DIMENSION(:,:) :: hv, hvp, hvpp
   REAL, ALLOCATABLE, DIMENSION(:,:) :: agy, ary
   REAL, ALLOCATABLE, DIMENSION(:) :: eagy, eary

   !.....Arrays at pressure points.....
   REAL, ALLOCATABLE, DIMENSION(:,:) :: h, hp, hpp
   REAL, ALLOCATABLE, DIMENSION(:,:) :: sal, salp, salpp, rhop
   REAL, ALLOCATABLE, DIMENSION(:  ) :: uout, vout, wout, scout, &
                                        Avout, Dvout, sal1, uhout

   !.....Arrays at vertical interfaces between pressure-pts.....
   REAL, ALLOCATABLE, DIMENSION(:,:) :: Av, Dv, Dvm  !Andrea
   REAL, ALLOCATABLE, DIMENSION(:,:) :: wp

   !.....Arrays at zeta-pts in 3D space (i.e. in the x-y plane) .....
   REAL   , ALLOCATABLE, DIMENSION(:) :: s, sp, spp, sx, sy
   REAL   , ALLOCATABLE, DIMENSION(:) :: hhs, hhu, hhv
   REAL   , ALLOCATABLE, DIMENSION(:) :: dd, qq, rr
   INTEGER, ALLOCATABLE, DIMENSION(:) :: kmz, k1z,k1u,k1v

   ! ... Logical arrays for 3D space
   LOGICAL, ALLOCATABLE, DIMENSION(:,:  ) :: mask2d
   LOGICAL, ALLOCATABLE, DIMENSION(:  ) :: mask

   !.....Matrix solution arrays & constants for block solver
   INTEGER, PARAMETER :: mdim  = 3
   INTEGER, DIMENSION(mdim ) :: jcoef
   REAL   , ALLOCATABLE, DIMENSION(:,:) :: coef
   REAL   , ALLOCATABLE, DIMENSION(:  ) :: rhs, zeta, wksp, ubar, rparm
   INTEGER, ALLOCATABLE, DIMENSION(:  ) :: iwksp, iparm, ip, jp

   ! ... Matrix solution arrays & constants for sparse solver
   INTEGER, PARAMETER :: mdimA = 5
   INTEGER :: ndimA, maxnzA, nwA, inwA, ibdwdA
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jcoefA ! Storage format 1 - NorthDelta
   REAL   , ALLOCATABLE, DIMENSION(:,:) :: coeffA ! Storage format 1 - NorthDelta

   !.....Other arrays.....
   REAL   , ALLOCATABLE, DIMENSION(:,:) :: sbc1, qbc1, salbc1
   REAL   , ALLOCATABLE, DIMENSION(:  ) :: ds
   REAL   , ALLOCATABLE, DIMENSION(:  ) :: zlevel

   !         -----Transient Allocatable Array Declarations-----
   !.....Arrays for interpolations at vertical interfaces.....
   REAL, ALLOCATABLE, DIMENSION(:,:) :: haxpp, haypp, hdxpp, hdypp,haxpp2,kh
   REAL, ALLOCATABLE, DIMENSION(:,:) :: haxpptr, hayppsal, haxppsal
   REAL, ALLOCATABLE, DIMENSION(:,:) :: th, th1,haypp2,th2,th12

   !         -----Surface Temperature Boundary Conditions Vars. & Arrays ---------

   INTEGER, PARAMETER :: nvSurfbcP = 6; ! No. of variables read on pre-process mode
   INTEGER, PARAMETER :: nvSurfbcR = 9; ! No. of variables read on run-time mode
   INTEGER :: ifSurfbc          ! Flag that specifies whether surf.BC. are used
   INTEGER :: nmetstat          ! No. of met stations to be used in interpolating met variables
   REAL    :: Qsw, Qn, Qlw, eta ! Shortwave, Net, longwave fluxes and attenuation coeff.
   REAL    :: Ta , Pa, Rh , Cc  ! Air temp., pressure, relative humidity & cloud cover
   REAL    :: dtSurfbc          ! Time (secs.) between records in surfbc.txt file
   REAL, ALLOCATABLE, DIMENSION(:) :: cdw, uair, vair   ! Wind drag coeff., EW and NS wind speeds
   REAL, ALLOCATABLE, DIMENSION(:,:) :: surfbc1
   REAL, ALLOCATABLE, DIMENSION(:,:) :: HeatSource, QswFr
   REAL, ALLOCATABLE, DIMENSION(:  ) :: Qsw2d , Qlw2d    ! Shortwave flux, Longwave flux = constant
   REAL, ALLOCATABLE, DIMENSION(:  ) :: Ta2d  , RH2d     ! Air Temperature and Relative Humidity
   REAL, ALLOCATABLE, DIMENSION(:  ) :: Cc2d             ! Cloud cover
   REAL, ALLOCATABLE, DIMENSION(:  ) :: uair2d, vair2d   ! EW and NS component of air vel.
   REAL, ALLOCATABLE, DIMENSION(:  ) :: metxy
   REAL, ALLOCATABLE, DIMENSION(:,:,:) :: weightst
   REAL, ALLOCATABLE, DIMENSION(:,:  ) :: weightn  !new weightst MAC
   REAL, ALLOCATABLE, DIMENSION(:,:) :: cw_sin, cw_cos
   ! --- Multiquadric interpolation
   REAL, ALLOCATABLE, DIMENSION(:,:  ) :: mqQij  ! arrays used in the MQ interpolation method
   INTEGER, ALLOCATABLE, DIMENSION (:) :: indice ! arrays used in the MQ interpolation method
   REAL, ALLOCATABLE, DIMENSION(:    ) :: QswMQ            ! Shortwave flux, Net flux attn.coeff = constant
   REAL, ALLOCATABLE, DIMENSION(:    ) :: TaMQ             ! Air Temperature
   REAL, ALLOCATABLE, DIMENSION(:    ) :: RhMQ             ! Relative Humidity
   REAL, ALLOCATABLE, DIMENSION(:    ) :: QlwMQ            ! Longwave - Atmospheric Pressure considered as constant
   REAL, ALLOCATABLE, DIMENSION(:    ) :: uairMQ, vairMQ   ! Drag coeff. for wind, EW and NS component of air vel.

   !         ----- Mapping arrays from 2D-lk to 3D-ijk -----------------------------
   INTEGER :: lm, lm1,cm1
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ij2l
   INTEGER, ALLOCATABLE, DIMENSION(:  ) :: l2i, l2j, &
                                           lEC, lWC, &
                                           lNC, lSC,c2l,l2c

   !         ----- Arrays & vars. for 2Eq. Turb. models ----------------------------
   INTEGER, PARAMETER :: real_G1 = SELECTED_REAL_KIND(P=15,R=80)
   REAL         , PARAMETER :: AvMolecular = 1.3e-6 ! Molecular viscosity at 10oC
   REAL         , PARAMETER :: DvMolecular = 1.4e-7 ! Molecular diffusivity at 10oC
   REAL(real_G1), PARAMETER :: q2_min = 2.e-6              !(Gaspar et al , 1990)
   REAL(real_G1), PARAMETER :: q2l_min= 6.80236723501E-009 !(Ivey&Imberger, 1991)
   REAL(real_G1), PARAMETER :: Gh_max_Kc = 2.9E-2 ! Max.Gh (Kantha&Clayson, 1994)
   REAL(real_G1), PARAMETER :: Gh_min_Kc =-2.8E-1 ! Min.Gn (Kantha&Clayson, 1994)
   REAL(real_G1), PARAMETER :: A_1 = 0.92, A_2 = 0.74, B_1 = 16.6, &
                             & B_2 = 10.1, E_1 = 1.80, E_2 = 1.33, &
                             & E_3 = 0.25, C_1 = 0.08, C_2 = 0.70, &
                             & C_3 = 0.20, Sq  = 0.20
   REAL(real_G1)            :: t_1, t_2, t_3, t_4, t_5, &
                             & var_1, var_2, var_3,  &
                             & var_4, var_5, var_6
   REAL(real_G1), ALLOCATABLE, DIMENSION(:,:) :: q2 , q2p , q2pp
   REAL(real_G1), ALLOCATABLE, DIMENSION(:,:) :: q2l, q2lp, q2lpp
   REAL(real_G1), ALLOCATABLE, DIMENSION(:  ) :: dsT, sal1T
   REAL(real_G1), ALLOCATABLE, DIMENSION(:,:) :: aaT

   !         ----- Arrays & vars. for GOTM interface --------------------------------
   REAL, ALLOCATABLE, DIMENSION(:,:) :: si3dtke
   REAL, ALLOCATABLE, DIMENSION(:,:) :: si3deps
   REAL, ALLOCATABLE, DIMENSION(:,:) :: si3dlen

   !         ----- Arrays & vars. for open lateral boundary conditions ---------------
   REAL                                :: dtsecOpenBC
   INTEGER, DIMENSION(maxnopen)        :: iside, itype
   INTEGER, DIMENSION(maxnopen)        :: isbc, jsbc, iebc, jebc
   REAL, ALLOCATABLE, DIMENSION(:,:,:) :: varsOpenBC
   REAL, ALLOCATABLE, DIMENSION(:,:  ) :: uhEB  , uhWB  , huEB  , huWB  , &
                                          uhEBp , uhWBp , huEBp , huWBp , &
                                          uhEBpp, uhWBpp, huEBpp, huWBpp, nnH,nnHH
   REAL, ALLOCATABLE, DIMENSION(:,:  ) :: vhNB  , vhSB  , hvNB  , hvSB  , &
                                          vhNBp , vhSBp , hvNBp , hvSBp , &
                                          vhNBpp, vhSBpp, hvNBpp, hvSBpp

   !         ----- Arrays & vars. for nesting boundary outputs  -----------------------
   INTEGER, PARAMETER                  :: ioNBTOGGLE = 1  ! = 0: No nesting; >= 1: Nesting routines on
   ! ... Variables used in outputting Nested Grid Boundaries
   INTEGER                             :: nxNBO           ! No. of interior boundaries IB to output
   INTEGER                             :: ioNBO           ! No. of steps between consecutive outputs
   INTEGER                             :: xxNBO           ! Horizontal scaling ratio in nesting procedure
   INTEGER, DIMENSION (maxnopen)       :: isbcNBO,       &
                                       &  jsbcNBO,       &
                                       &  iebcNBO,       &
                                       &  jebcNBO         ! Cells where IB being output start and finishes
   INTEGER, DIMENSION (maxnopen)       :: nfrNBO,        &! No. of frames being output for each IB
                                       &  iptNBO,        &! No. of cells being output for each IB
                                       &  ntrNBO,        &! No. of tracers being output for each IB
                                       &  isdNBO          ! Side for each IB file being output
   !  ... Variables used to input Nested Grid Boundaries
   INTEGER, DIMENSION (maxnopen)       :: nfrNBI,        &! No. of frames being read for each IB
                                       &  iptNBI,        &! No. of cells being read for each IB
                                       &  ntrNBI,        &! No. of tracers being output for each IB
                                       &  isdNBI          ! Side for each IB file being input
   ! .... Variables used in to force the model through Nested Grid Boundaries
   REAL   ,ALLOCATABLE, DIMENSION(:,:)::  uhNGB, uhNGBp, &
                                       &  vhNGB, vhNGBp, &
                                       &  scNGB, scNGBp
   REAL ,ALLOCATABLE, DIMENSION(:,:,:)::  trNGB, trNGBp
   INTEGER,ALLOCATABLE, DIMENSION(:,:)::  kNGB, iNGB, jNGB
   REAL ,ALLOCATABLE, DIMENSION(:) ::  thrsNGB, thrsNGBp   !MAC
   INTEGER, PARAMETER :: nbiid0 = 980
   INTEGER, PARAMETER :: nboid0 = 960

   !         ----- Vars. & Arrays used in new OUTPUT routines (sections, planes
   !               isotherms and complete outputs ---------
   INTEGER :: iop, iox, ioi, ioc, n_planes, n_sections, n_isot, isot_points
   INTEGER, PARAMETER :: max_planes = 10, max_sections = 30, max_isot = 20
   INTEGER, PARAMETER :: max_section_cells = 100
   INTEGER, ALLOCATABLE, DIMENSION (:) :: interior_plane_points
   INTEGER, ALLOCATABLE, DIMENSION (:) :: interior_section_points
   INTEGER, DIMENSION (max_planes  ) :: p_out
   INTEGER, DIMENSION (max_sections) :: n_section_cells
   REAL   , DIMENSION (max_isot    ) :: isot_value
   INTEGER, DIMENSION (max_sections, max_section_cells) :: xinode, xjnode
   REAL, ALLOCATABLE, DIMENSION(:,:) :: out_array


   !          ----- Active & Non-active Scalar transport models  **********************
   INTEGER :: iotr       ! No. steps between output to file
   INTEGER :: ntr        ! No. of tracers requested for simulation
   INTEGER, PARAMETER  :: ntrmax = 15 ! Max. No. of tracers that can be simulated
   REAL   , ALLOCATABLE, DIMENSION(:,:)  :: fluxX, fluxY, fluxZ
   REAL   , ALLOCATABLE, DIMENSION(:,:)  :: fluxXtr, fluxY2, fluxZ2, fluxXsal
   REAL   , ALLOCATABLE, DIMENSION(:,:,:):: tracer, tracerpp, sourcesink
   REAL   , ALLOCATABLE, DIMENSION(:,:)  :: trout
   INTEGER, ALLOCATABLE, DIMENSION(:  )  :: trct0, trctn
   REAL   , ALLOCATABLE, DIMENSION(:  )  :: trcx0, trcy0, trcz0, &
                                            trcsx, trcsy, trcsz, trcpk

   !          ----- Point Sources & Sinks Eqs. Vars. & Arrays *********************
   ! ... Variables used specificallly to model plumes - 
   REAL    :: k4sod      ! 5.7870E-6 ! g/m2/s
   REAL(8) :: salamb     ! Vars. for generalized version of plume model
   REAL(8) :: patm       ! Vars. for generalized version of plume model
   REAL   , ALLOCATABLE, DIMENSION(:    )  :: lambda ! Plume width
   REAL   , ALLOCATABLE, DIMENSION(:    )  :: diammb ! Initial diammeter of bubbles
   INTEGER, ALLOCATABLE, DIMENSION(:    )  :: kdetr          ! k for detrainment cell
   INTEGER, ALLOCATABLE, DIMENSION(:    )  :: idetr          ! How detrainment is modelled in diffuser devices 
   REAL   , ALLOCATABLE, DIMENSION(:    )  :: dfL            ! Diffuser length (device)

   ! ... Variables for boundary conditions as point sources and sinks BCasPSS
   INTEGER :: npssdev    ! No. of devices producing point sinks and sources (i.e. diffusers)
   INTEGER :: iopss      ! No. of colums with sources/sinks
   REAL    :: dtsecpss   ! Time in seconds between consecutive records in io files  
   INTEGER, ALLOCATABLE, DIMENSION(:    )  :: iopssH
   INTEGER, ALLOCATABLE, DIMENSION(:,:  )  :: ioph2iop
   INTEGER, ALLOCATABLE, DIMENSION(:    )  :: ptype          ! Type of PSS simulated (device)
   REAL   , ALLOCATABLE, DIMENSION(:,:,:)  :: varspss        ! flows, temps. & tracers = f(time,dev)
   INTEGER, ALLOCATABLE, DIMENSION(:    )  :: pdt            ! Frequency of update (device) 
   INTEGER, ALLOCATABLE, DIMENSION(:    )  :: iodev          ! Device No. for each pss (column) 
   INTEGER, ALLOCATABLE, DIMENSION(:    )  :: ipss,jpss      ! Grid location of pss (column)
   REAL   , ALLOCATABLE, DIMENSION(:,:  )  :: Qpss	         ! Inflow/outflow rate at n+1 (column)
   REAL   , ALLOCATABLE, DIMENSION(:,:  )  :: Tpss           ! Temp. for pss at n+1 - (column)
   REAL   , ALLOCATABLE, DIMENSION(:,:,:)  :: Rpss           ! Tracers for pss at n + 1 - (column)
   REAL   , ALLOCATABLE, DIMENSION(:    )  :: uWpss          ! Enforce flow rates on EW direction (device)
   REAL   , ALLOCATABLE, DIMENSION(:    )  :: uEpss          ! Enforce flow rates on EW direction (device)
   REAL   , ALLOCATABLE, DIMENSION(:    )  :: vNpss          ! Enforce flow rates on NS direction (device)
   REAL   , ALLOCATABLE, DIMENSION(:    )  :: vSpss          ! Enforce flow rates on NS direction (device)
   REAL   , ALLOCATABLE, DIMENSION(:    )  :: flpss          ! Water inflow rate = f(type of pss) (device)
   REAL   , ALLOCATABLE, DIMENSION(:    )  :: scpss          ! Active scalar in point source (device)
   REAL   , ALLOCATABLE, DIMENSION(:,:  )  :: trpss          ! Tracer load (grams/second) (device)
   REAL   , ALLOCATABLE, DIMENSION(:    )  :: qthrs          ! Threshold to determine whether diffuser is on/off


   !        -------- ***************************************
   INTEGER, PARAMETER :: testcase = 3 ! FJR

   !        -------- Energy & Scalar balances *****************************************
   INTEGER, PARAMETER :: iobal = 1 ! IF scalar/energy balances are to be conducted
   REAL(real_G1)::  ShearProduction, BuoyancyProduction, Dissipation, TKinE, &
   &                ShearCum, BuocyCum, DisspCum

   !        -------- Interpolation procedures *****************************************
   INTEGER :: iinterp  ! Changed 12/2010 SWA
   REAL :: gammaB, delNfactor  ! Changed 12/2010 SWA

   !        -------- ECOMOD - Size Structure ******************************************
   INTEGER :: ecomod ! Type of behaviour assigned to tracers 
   REAL, ALLOCATABLE, DIMENSION (:) :: diamSZ  ! Diameter (m)
   REAL, ALLOCATABLE, DIMENSION (:) :: rmaxSZ  ! Maximum growth rate (s-1)
   REAL, ALLOCATABLE, DIMENSION (:) :: vdep    ! Settling velocity (m/s)

   !        -------- ECOMOD - Sediment Transport **************************************
   REAL, ALLOCATABLE, DIMENSION (:) :: sdenSED  ! Submerged specific density (adimensional) 
   REAL, ALLOCATABLE, DIMENSION (:) :: diamSED  ! Particle diameter (m) 
   REAL, PARAMETER :: Ased = 1.3e7; 
   !                                   vdep (uses the same array as ECOMOD = 1)
   REAL :: dthrsWIBSS, &                      ! Hours between WIBSS fields (input)
           thrsWIBSS , &                      ! Hours from start to current WIBSS field
           thrsWIBSSp, &                      ! Hours from start to past    WIBSS field
           thrsWIBSS0                         ! Hours from start to first WIBSS field
   REAL, ALLOCATABLE, DIMENSION (:,:) :: wibssIOx, wibssIOxp
   REAL, ALLOCATABLE, DIMENSION (:,:) :: wibssIOy, wibssIOyp
   REAL, ALLOCATABLE, DIMENSION (:,:) :: wibssx  , wibssy

   !         ------- ECOMOD - Water Quality   ******************************************

   ! ... Integer switches to deterimine if constituent is modeled
   INTEGER:: iARB, iDO , iPON, iDON, iNH4, iNO3, iPOP, iDOP, iPO4
   INTEGER:: iALG, iDOM, iPOM, iSOD

   ! ... Integer to determine constituent location
   INTEGER:: LARB, LDO , LPON, LDON, LNH4, LNO3, LPOP, LDOP, LPO4
   INTEGER:: LALG, LDOM, LPOM, LSOD

   ! ... Model Constants - read from input file and many used for calibration
   ! - Stochiometeric constants
   REAL   :: acc, anc, apc, roc, ron
 	
   ! - Half-saturation values and algal prefernce for NH4	
   REAL   :: KDOM, KNIT, KSN, KSP, FNH4	

   ! - Model rates
   REAL   :: k_a, k_arb, k_dn, k_DOM, k_ex, k_gr, k_hc
   REAL	  :: k_hn, k_hp, k_mn, k_mor, k_mp, k_n, k_ra
   REAL	  :: k_rs, k_set, k_setarb, k_vn, mu_max		

   ! - Temperature depependence factors
   REAL	  :: Theta_a  , Theta_arb, Theta_dn, Theta_DOM, Theta_ex
   REAL	  :: Theta_gr , Theta_hc , Theta_hn, Theta_hp , Theta_mn
   REAL	  :: Theta_mor, Theta_mp , Theta_mu, Theta_n  , Theta_PON
   REAL	  :: Theta_ra , Theta_SOD, Theta_vn						

   ! - Atmoshperic deposition rates
   REAL	  :: ATM_DON, ATM_NH4, ATM_NO3, ATM_PO4, ATM_POP

   ! Groundwater flux rates
   REAL	  :: GW_NH4, GW_NO3, GW_PO4

   ! Sediment release rates
   REAL	  :: J_NH4, J_NO3, J_PO4



!                        -----Data Dictionary-----

END MODULE si3d_types

