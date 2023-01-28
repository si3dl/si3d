!************************************************************************
                          MODULE si3d_sed
!************************************************************************
!
!  Purpose: Procedures that implement routines related to the modelling
!           of ecological processees
!
!-------------------------------------------------------------------------

  USE si3d_types

  IMPLICIT NONE
  SAVE

CONTAINS


! ********************************************************************
SUBROUTINE sourceSS(kwq,lwq)
! ********************************************************************
!
! Purpose: If suspended sediments is modeled, this subroutine
!          calculates the source and sink terms that depend on 
!          the suspended sediments. 
!
! --------------------------------------------------------------------
  ! Arguments of subroutine
  integer, intent(in)   :: kwq               !< Layer index
  integer, intent(in)   :: lwq               !< Column index
  real                  :: taub              !< (Pa) Bottom shear stress
  real                  :: ustarb            !< (m/s) Shear velocity
  real                  :: w_dens            !< (kg/m3) Water density
  real, dimension(sedNumber) :: settling_vel !< (m/s) Settling velocity of sediment
  real, dimension(sedNumber) :: Rep          !< Explicity Particle Reynolds Number
  real, dimension(sedNumber) :: tauCrt       !< Critical Shear Stress
  real, dimension(sedNumber) :: erosionFlux  !<
  real, dimension(sedNumber) :: depositionFlux  !<
  real :: conc
  integer :: pn
  integer :: i
  integer :: kms

  kms = kmz(lwq)
  
  w_dens = rhop(kwq, lwq) + 1000

  ! print*, '--------------------------'
  ! print*, 'k =',kwq,'l =',lwq

  ! Estimate properties of sediment for a given water density at bottom cell
  call get_sed_prop(settling_vel,Rep,tauCrt,sed_diameter,sed_dens,w_dens)
  ! Estimate bottom shear stress
  call tauBottom(kwq, lwq, taub, ustarb)
  ! print*, 'taub =',taub
  ! Estimate erosion flux
  pn = 0

  do i = 1, sedNumber
    ! print*, 'pn = ', pn + 1
    if (kwq .eq. kms) then 
      call erosion(erosionFlux(i), taub, ustarb, Rep(i), settling_vel(i), tauCrt(i))
    else
      erosionFlux(i) = 0.0
    end if 
    conc = tracerpp(kwq,lwq,LSS1 + pn)
    ! print*, 'conc =',conc

    call deposition(depositionFlux(i), settling_vel(i), tauCrt(i), taub, conc, kwq, kms)
    ! print*, 'erosionFlux =',erosionFlux(i)
    ! print*, 'depositionFlux =',depositionFlux(i)

    sourcesink(kwq, lwq, LSS1 + pn) = (erosionFlux(i) - depositionFlux(i)) / hp(kwq, lwq)

    ! print*,'sourcesink LSS1 = ',sourcesink(kwq,lwq,LSS1)
    ! print*,'sourcesink LSS2 = ',sourcesink(kwq,lwq,LSS2)
    pn = pn + 1
    
  end do

END SUBROUTINE sourceSS

! ********************************************************************
SUBROUTINE erosion(erosionFlux, taub, ustarb, Rep, settling_vel, tauCrt)
! ********************************************************************
!
! Purpose: To estimate the erosion caused by the flow at the bottom
!         cell. It follows the flux method by Garcia and Parker 1991,
!         1993, Reardon et al., 2014, others. 
!
! --------------------------------------------------------------------
  ! Arguments
  real, intent(in) :: taub
  real, intent(in) :: ustarb
  real, intent(in) :: tauCrt
  real, intent(in) :: Rep !< Explicity Particle Reynolds Number
  real, intent(in) :: settling_vel !< (m/s) Settling velocity of sediment
  real   :: z_u 
  real  :: E_s
  real, intent(out) :: erosionFlux

  if (Rep .lt. 3.5) then
    z_u = 0.708 * ustarb * (Rep ** 0.6) / settling_vel
  else
    z_u = ustarb * (Rep ** 0.6) / settling_vel
  end if 

  ! Sediment entrainment coefficient
  E_s = Ased * (z_u ** 5) / (1 + (z_u ** 5) * Ased/0.3)
  if (taub .gt. tauCrt) then 
    erosionFlux = E_s * settling_vel
  else
    erosionFlux = 0
  end if   

  ! print*, 'E_s',E_s

  return
END SUBROUTINE erosion


! ********************************************************************
SUBROUTINE deposition(depositionFlux, settling_vel, tauCrt, taub, conc, kwq, kms)
! ********************************************************************
!
! Purpose: To estimate the suspended sediment deposition at the
!         bottom layer for each wet column
!
! --------------------------------------------------------------------
  ! Arguments
  real, intent(in) :: settling_vel
  real, intent(in) :: taub
  real, intent(in) :: tauCrt
  real, intent(in) :: conc
  integer, intent(in) :: kwq
  integer, intent(in) :: kms
  real, intent(out) :: depositionFlux

    if ((taub .lt. tauCrt) .and. (kwq .eq. kms)) then
      depositionFlux = settling_vel * conc * (1 - taub/tauCrt)
    else
      depositionFlux = settling_vel * conc
    end if

  return
END SUBROUTINE deposition

! ********************************************************************
SUBROUTINE get_sed_prop(settling_vel,Rep,tauCrt,sed_diameter,sed_dens,w_dens)
! SUBROUTINE get_sed_prop(settling_vel,Rep,tauCrt,sed_diameter,sed_dens,w_dens)
! ********************************************************************
!
! Purpose: Estimate particle dependent parameters / properties
!
! --------------------------------------------------------------------
  implicit none
  ! Arguments
  real, dimension(sedNumber), intent(in)  :: sed_diameter     !< (um) Sediment diameter
  real, intent(in)                        :: w_dens           !< (kg/m3) water density
  real, dimension(sedNumber), intent(in)  :: sed_dens         !< (kg/m3) sediment density
  real, dimension(sedNumber), intent(out) :: settling_vel     !< (m/s) settling velocity
  real, dimension(sedNumber), intent(out) :: Rep              !< Explicit Particle Reynolds Number
  real, dimension(sedNumber), intent(out) :: tauCrt           !< (Pa) Critical shear stress 
  real, dimension(sedNumber)              :: submerged_spec_g !< Sediment submerged specific gravity
  real, dimension(sedNumber)              :: sed_diamm        !< (m) Sediment diameter
  real, dimension(sedNumber)              :: sed_spec_g       !< specific sediment gravity
  logical                                 :: ivanRijn         !< Flag for using van Rijn (1984) formula or Dietrich (1982). The default is van Rijn

  ivanRijn = .true.

  ! print*, 'w_dens',w_dens,'sed_diameter',sed_diameter
  ! To estimate sediment diamenter in m
  sed_diamm = sed_diameter * 0.000001

  sed_spec_g = sed_dens/w_dens - 1.

  call submergedSpecificGravity(submerged_spec_g, sed_dens, w_dens)

  ! print*, 'R',submerged_spec_g

  call partReynolds_Number(Rep, sed_diamm, kinematic_viscosity, submerged_spec_g)

  ! print*, 'Rep',Rep

  call settling_velocity(settling_vel, g, submerged_spec_g, sed_spec_g, &
                         Rep, sed_diamm, kinematic_viscosity, ivanRijn)

  ! print*, 'settling_vel',settling_vel

  call tauCritical(tauCrt,g,sed_diamm,submerged_spec_g, sed_dens,kinematic_viscosity, Rep)

  ! print*, 'tauCrt',tauCrt

  return
END SUBROUTINE get_sed_prop

! ********************************************************************
SUBROUTINE submergedSpecificGravity(submerged_spec_g, sed_dens, w_dens)
! ********************************************************************
!
! Purpose: To estimate the submerged specific gravity for a given
!          Particle type
!
! --------------------------------------------------------------------
  ! Arguments of subroutine
  real, dimension(sedNumber), intent(in)  :: sed_dens         !< (kg/m3) Sediment density
  real, intent(in)                        :: w_dens           !< (kg/m3) Water density
  real, dimension(sedNumber), intent(out) :: submerged_spec_g !< Sediment submerged specific gravity

  ! Estimate submerged specific gravity
  submerged_spec_g = sed_dens / w_dens - 1

  return
END SUBROUTINE submergedSpecificGravity

! ********************************************************************
SUBROUTINE partReynolds_Number(Rep, sed_diamm, kinematic_viscosity, submerged_spec_g)
! ********************************************************************
!
! Purpose: To estimate the explicit Particle Reynolds Number.
!        The equation used is from Garcia and Parker 1993.
!        Experiments on the entrainment of sediment into suspension
!        by a dense bottom current. DOI: 10.1029/92JC02404
!
! --------------------------------------------------------------------

  ! Arguments of subroutine
  real, dimension(sedNumber), intent(in)  :: sed_diamm            !< (m) Sediment diameter D50
  real, intent(in)                        :: kinematic_viscosity  !< (m2/sec) kinematic viscosity of water
  real, dimension(sedNumber), intent(in)  :: submerged_spec_g     !< Sediment submerged specific gravity
  real, dimension(sedNumber), intent(out) :: Rep                  !< Explicit Particle Reynolds Number

  Rep = sqrt(g * submerged_spec_g * sed_diamm ** 3) / kinematic_viscosity

  return
END SUBROUTINE partReynolds_Number

! ********************************************************************
SUBROUTINE settling_velocity(settling_vel, g, submerged_spec_g, sed_spec_g, &
                             Rep, sed_diamm, kinematic_viscosity, ivanRijn)
! ********************************************************************
!
! Purpose: To estimate the settling velocity for a given particle
!          size. The estimate uses 
!
! --------------------------------------------------------------------
  ! Arguments of subroutine
  real, intent(in)                        :: g                    !< (m/s2)Gravitational acceleration (m/s**2)
  real, dimension(sedNumber), intent(in)  :: submerged_spec_g     !< Sediment submerged specific gravity
  real, intent(in)                        :: kinematic_viscosity  !< (m2/s) Kinematic viscosity of water
  real, dimension(sedNumber), intent(in)  :: sed_diamm            !< (m) Sediment Diameter
  real, dimension(sedNumber), intent(in)  :: sed_spec_g           !< Sediment specific gravity
  real, dimension(sedNumber), intent(in)  :: Rep                  !< Explicit Particle Reynolds Number
  logical, optional                       :: ivanRijn             !< Flag for using van Rijn (1984) formula or Dietrich (1982). The default is van Rijn      
  real                                    :: dimless_fall_vel     !< dimensionaless fall velocity
  logical                                 :: vanRijnFlag
  integer                                 :: i
  ! Parameter for Dietrich (1982) equation
  ! Commented values are from Bombardelli and Moreno 2012 found in Reardon et al., 2014
  real                                    :: b_1 = 3.76715 ! 2.891394
  real                                    :: b_2 = 1.92944 ! 0.95296
  real                                    :: b_3 = 0.09815 ! 0.056835
  real                                    :: b_4 = 0.00575 ! 0.002892
  real                                    :: b_5 = 0.00056 ! 0.000245
  real, dimension(sedNumber), intent(out) :: settling_vel         !< (m/s) Settling

  if ( present(ivanRijn) ) then
    vanRijnFlag = ivanRijn
  end if

  do i = 1,sedNumber

    SELECT CASE (vanRijnFlag)
      CASE (.true.)
        ! Van Rijn Formula
        if (sed_diamm(i) .ge. 1.0d-3) then
          settling_vel(i) = 1.1 * sqrt(submerged_spec_g(i) * g * sed_diamm(i))
        elseif (sed_diamm(i) .gt. 1.0d-4 .and. sed_diamm(i) .le. 1.0d-3) then
          settling_vel(i) = (10 * kinematic_viscosity / sed_diamm(i)) *        & 
                         (sqrt(1 + 0.01 * (submerged_spec_g(i) * g       &
                          * sed_diamm(i) **3) / kinematic_viscosity ** 2.) - 1)
        elseif (sed_diamm(i) .gt. 6.5d-5 .and. sed_diamm(i) .le. 1.0d-4) then
          ! Stokes Law
          settling_vel(i) = (submerged_spec_g(i) * g * sed_diamm(i) ** 2.) / (18.0 * kinematic_viscosity)
        end if

      CASE (.false.)
        dimless_fall_vel = exp(-1.*b_1 + b_2 * log(Rep(i)) - b_3 * (log(Rep(i))) ** 2.0 - b_4 * (log(Rep(i))) ** 3. + b_5 * (log(Rep(i))) ** 4.)
        if ( sed_diamm(i) .lt. 1.0d-5) then
          settling_vel(i) = (submerged_spec_g(i) * g * sed_diamm(i)**2)/(18.*kinematic_viscosity)
        else
          settling_vel(i) = dimless_fall_vel * sqrt(submerged_spec_g(i) * g * sed_diamm(i))
        end if
    END SELECT
  end do
  return
END SUBROUTINE settling_velocity

! ********************************************************************
SUBROUTINE tauCritical(tauCrt,g,sed_diamm,submerged_spec_g, sed_dens,kinematic_viscosity, Rep)
! ********************************************************************
!
! Purpose: To estimate the critical shear stress for a given particle
!          size. The estimate uses the nondimensional critical shields
!          parameter from Parker et al., 2003.
!
! --------------------------------------------------------------------
  ! Arguments of subroutine
  real, intent(in)                        :: g                   !< (m/s2) Gravitational acceleration (m/s**2)
  real, dimension(sedNumber), intent(in)  :: sed_diamm           !< (m) Sediment Diameter
  real, dimension(sedNumber), intent(in)  :: submerged_spec_g    !< Sediment submerged specific gravity
  real, dimension(sedNumber), intent(in)  :: sed_dens            !< (kg/m3) Sediment density
  real, intent(in)                        :: kinematic_viscosity !< (m2/s) Kinematic viscosity of water
  real, dimension(sedNumber), intent(in)  :: Rep                 !< Explicit Particle Reynolds Number
  real, dimension(sedNumber)              :: shields_param       !< Nondimensional Critical Shields Parameter
  real, dimension(sedNumber), intent(out) :: tauCrt              !< (Pa) Critical shear stress

  ! Estimate of the nondimensional critical Shields parameter
  shields_param = 0.5 * (0.22 * Rep ** (-0.6) &
                  + 0.06 * 10 ** (-7.7 * Rep ** (-0.6)))

  ! Estimate of critical shear stress for given water and sediment properties
  tauCrt = shields_param * g * submerged_spec_g * sed_diamm * sed_dens

  return 
END SUBROUTINE tauCritical

! ********************************************************************
SUBROUTINE tauBottom(kwq,lwq,taub,ustarb)
! ********************************************************************
!
! Purpose: If suspended sediments is modeled, this subroutine
!          calculates the source and sink terms that depend on 
!          the suspended sediments. 
!
! --------------------------------------------------------------------
  ! Arguments of subroutine
  integer, intent(in)   :: kwq    !< Layer index
  integer, intent(in)   :: lwq    !< Column index
  integer               :: kmxp   !< min bottom layer btwn cell on east and point
  integer               :: kmyp   !< min bottom layer btwn cell on north and point
  integer               :: kmxm   !< min bottom layer btwn cell on west and point
  integer               :: kmym   !< min bottom layer btwn cell on south and point
  integer               :: kms    ! bottom layer
  real                  :: ubott  ! (m/s) vel u at bottom cell
  real                  :: vbott  ! (m/s) vel v at bottom cell
  real                  :: taubx  ! (Pa) bottom shear stress in x
  real                  :: tauby  ! (Pa) bottom shear stress in y
  real, intent(out)     :: taub   !< (Pa) Bottom shear stress
  real, intent(out)     :: ustarb !< (m/s) Shear velocity      

  !.....Compute bottom layer numbers  ....
  ! kwq must be the bottom
  kms = kmz(lwq)
  if (kwq .eq. kms) then
    kmyp= MIN( kmz( lNC(lwq)), kms)
    kmym= MIN( kmz( lSC(lwq)), kms)
    kmxp= MIN( kmz( lEC(lwq)), kms)
    kmxm= MIN( kmz( lWC(lwq)), kms)

    ! ....Compute Currend-Induced Bottom Shear Stress CIBSS (cd*rho*U^2)
    ubott = ( uhp( kmxp, lwq ) + uhp( kmxm, lWC(lwq)) ) / 2. / hp(kms,lwq)
    vbott = ( vhp( kmyp, lwq ) + vhp( kmym, lSC(lwq)) )/2. / hp(kms,lwq)
    taubx = cd * (rhop(kms,lwq)+1000.) * ubott**2.
    tauby = cd * (rhop(kms,lwq)+1000.) * vbott**2.

    ! ... Add Wave-Induced Bottom Shear Stress & calculate friction velocity
    ! taubx = taubx + wibssx(i,j)
    ! tauby = tauby + wibssy(i,j)
    taub  = sqrt(taubx**2.+tauby**2.)
    ustarb = sqrt(taub/1000.)

  endif

END SUBROUTINE tauBottom



















































































END MODULE si3d_sed