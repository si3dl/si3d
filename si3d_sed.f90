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
  integer, intent(in)        :: kwq               !< Layer index
  integer, intent(in)        :: lwq               !< Column index
  integer                    :: kms
  integer                    :: i
  real                       :: taub              !< (Pa) Bottom shear stress
  real                       :: ustarb            !< (m/s) Shear velocity
  real                       :: w_dens            !< (kg/m3) Water density
  real, dimension(sedNumber) :: Rep               !< Explicity Particle Reynolds Number
  real, dimension(sedNumber) :: tauCrt            !< Critical Shear Stress
  real(kind=8), dimension(sedNumber) :: resus_flux      !<
  real(kind=8), dimension(sedNumber) :: resus_sed  !<
  real(kind=8), dimension(sedNumber) :: flux_fluff
  real, dimension(sedNumber) :: deposition_flux   !<
  real, dimension(sedNumber) :: burial_flux
  real                       :: cb
  real :: c_sed

  resus_flux(:) = 0.0
  resus_sed(:) = 0.0
  deposition_flux(:) = 0.0
  burial_flux(:) = 0.0
  flux_fluff(:) = 0.0

  kms = kmz(lwq)
  if (kwq .eq. kms) then
    w_dens = (rhop(kwq, lwq) + 1000)
    ! Estimate bottom shear stress
    call tauBottom(taub, ustarb, kwq, lwq)
    call totalMass(sed_frac, kwq, lwq)

    do i = 1, sedNumber
      ! Estimate properties of sediment for a given water density at bottom cell
      cb = tracerpp(kwq,lwq,LSS1 + i - 1) ! kg/m3
      c_sed = tracerpp(kwq + 1, lwq, LSS1 + i - 1) ! kg/m3

      call get_sed_prop(settling_vel(i), Rep(i), tauCrt(i), sed_diameter(i), sed_dens(i), w_dens)
      ! Estimate erosion flux
      ! Estimate erosion flux fluff layer
      if ((sed_frac(i) .le. 0.0) .or. (c_sed .le. 0.0)) then
        flux_fluff(i) = 0.0
      else
        if ((l2i(lwq) .ge. 185) .and. (l2j(lwq) .ge. 65)) then
        ! if (i .eq. 2) then
          ! flux_fluff(i) = min(1e-3, 1.5e-6 * (sed_h * c_sed) * max(0.0, (taub / 0.02 - 1.0)))
          flux_fluff(i) = 0.0
        else
          ! flux_fluff(i) = min(1e-3, 1.5e-6 * (sed_h * c_sed) * max(0.0, (taub / 0.03 - 1.0)))
          flux_fluff(i) = 0.0
        end if
      end if
      
      if (taub .gt. tauCrt(i)) then
        if (sed_type(i) == 0) then
          call resuspension_noncohesive(resus_sed(i), ustarb, Rep(i), settling_vel(i), lwq)
          resus_flux(i) = resus_sed(i) * sed_frac(i) * c_sed
        else if (sed_type(i) == 1) then
          call resuspension_cohesive(resus_sed(i), taub, tauCrt(i), M_cohesive(i), lwq)
          resus_flux(i) = resus_sed(i) * sed_frac(i)
        end if
      else
        resus_flux(i) = 0.0
      end if

      resus_flux(i) = resus_flux(i) + flux_fluff(i)
      ! Estimate deposition flux
      if (cb .gt. 0.0) then
        if (sed_type(i) == 0) then
          call deposition_noncohesive(deposition_flux(i), settling_vel(i), tauCrt(i), taub, cb)
        else if (sed_type(i) == 1) then
          call deposition_cohesive(deposition_flux(i), settling_vel(i), tauCrt(i), taub, cb, ustarb, sed_dens(i), hp(kwq, lwq), lwq)
        end if
        ! To correct deposition flux as it can not remove more sediment than what is in the water layer on top of sediment
        if (deposition_flux(i) .gt. (cb * hp(kwq, lwq) / dt)) then
          deposition_flux(i) = cb * hp(kwq, lwq) / dt
        end if
      else
        deposition_flux(i) = 0.0
      end if

      sourcesink(kwq, lwq, LSS1 + i - 1) = resus_flux(i) - deposition_flux(i)

      call burial(burial_flux(i), resus_flux(i), deposition_flux(i))

      ! Estimate source and sink for the sediment cell.
      sourcesink(kwq + 1,lwq, LSS1 + i - 1) = deposition_flux(i) - resus_flux(i) - burial_flux(i)
    end do

    ! if ((l2i(lwq) .eq. 185) .and. (l2j(lwq) .eq. 80)) then
    !   print*, '----------------- OA04 SS Model ------------------'
    !   print*, 'vs =', settling_vel * 86400, 'm/d'
    !   print*, 'taub =',taub
    !   print*, 'tauCr =', tauCrt
    !   print*, 'sed_frac = ', sed_frac
    !   print*, 'resuspension_flux = ', resus_flux
    !   print*, 'depo_flux', deposition_flux
    !   print*, 'sed_conc = ', tracerpp(kwq + 1,lwq,LSS1:LSS3)
    ! end if
    ! if ((l2i(lwq) .eq. 83) .and. (l2j(lwq) .eq. 134)) then
    !   print*, '----------------- UA06 SS Model ------------------'
    !   print*, 'vs =', settling_vel * 86400, 'm/d'
    !   print*, 'taub =',taub
    !   print*, 'tauCr =', tauCrt
    !   print*, 'sed_frac = ', sed_frac
    !   print*, 'resuspension_flux = ', resus_flux
    !   print*, 'depo_flux', deposition_flux
    !   print*, 'sed_conc = ', tracerpp(kwq + 1,lwq,LSS1:LSS3)
    ! end if
    ! if ((l2i(lwq) .eq. 170) .and. (l2j(lwq) .eq. 47)) then
    !   print*, '----------------- LA03 SS Model ------------------'
    !   print*, 'vs =', settling_vel * 86400, 'm/d'
    !   print*, 'taub =',taub
    !   print*, 'tauCr =', tauCrt
    !   print*, 'sed_frac = ', sed_frac
    !   print*, 'resuspension_flux = ', resus_flux
    !   print*, 'depo_flux', deposition_flux
    !   print*, 'sed_conc = ', tracerpp(kwq + 1,lwq,LSS1:LSS3)
    ! end if

  else
    do i = 1, sedNumber
      sourcesink(kwq,lwq,LSS1 + i - 1) = 0.0
    end do
  end if

  fluxes_out(kwq, lwq, 33) = sum(resus_flux)
  fluxes_out(kwq, lwq, 34) = sum(deposition_flux)

  fluxes_out(kwq + 1, lwq, 33) = sum(resus_flux)
  fluxes_out(kwq + 1, lwq, 34) = sum(deposition_flux)

END SUBROUTINE sourceSS

! ********************************************************************
SUBROUTINE totalMass(sed_frc, kwq, lwq)
! ********************************************************************
!
! Purpose: To estimate the total mass in the sediment layer and
!         estimate the sediment fraction in the bed for each type
!         of sediment.
!
! --------------------------------------------------------------------
  ! Arguments
  real, intent(inout), dimension(sedNumber) :: sed_frc
  integer, intent(in)                       :: kwq
  integer, intent(in)                       :: lwq
  real , dimension(sedNumber)               :: massSed
  real                                      :: totalMass_bed
  integer                                   :: i

  totalMass_bed = 0.0
  do i = 1, sedNumber
    massSed(i) = tracerpp(kwq+1,lwq,LSS1 + i - 1) * dx * dy * h(kwq+1,lwq)
  end do

  totalMass_bed = sum(massSed)
  do i = 1, sedNumber
    if (totalMass_bed .gt. 0.0) then
      sed_frc(i) = massSed(i) / totalMass_bed
    else
      sed_frc(i) = 1 / sedNumber
    end if
  end do

END SUBROUTINE totalMass

! ********************************************************************
SUBROUTINE resuspension_noncohesive(resus_flux, ustarb, Rep, sett_vel, lwq)
! ********************************************************************
!
! Purpose: To estimate the erosion caused by the flow at the bottom
!         cell. It follows the flux method by Garcia and Parker 1991,
!         1993, Reardon et al., 2014, others. 
!
! --------------------------------------------------------------------
  ! Arguments
  real, intent(in)  :: ustarb       !< (m/s) Shear velocity at bottom
  real, intent(in)  :: Rep          !< Explicit Particle Reynolds Number
  real(kind=8), intent(in) :: sett_vel !< (m/s) Settling velocity of sediment
  real              :: z_u          !< Similarity variable for uniform sediment
  real              :: E_s          !< Dimensionless coefficient for sediment entrainment. Under quasi-equilibrium conditions (Garcia & Parker 1991)
  real(kind=8), intent(out) :: resus_flux  !< Vertical erosion flux
  integer, intent(in) :: lwq

  if ((Rep .gt. 0.1) .and. (Rep .le. 1)) then
    z_u = 1 * ustarb * (Rep ** 3.75) / sett_vel
  else if ((Rep .gt. 1.0) .and. (Rep .le. 3.5)) then
    z_u = 0.586 * ustarb * (Rep ** 1.23) / sett_vel
  else if (Rep .gt. 3.5) then
    z_u = 1 * ustarb * (Rep ** 0.6) / sett_vel
  end if

  ! Sediment entrainment coefficient
  E_s = Ased * (z_u ** 5) / (1 + (z_u ** 5) * Ased/0.3)
  resus_flux = E_s * sett_vel

  return
END SUBROUTINE resuspension_noncohesive

! ********************************************************************
SUBROUTINE resuspension_cohesive(resus_flux, taub, tauCrt, M_param, lwq)
! ********************************************************************
!
! Purpose: To estimate the erosion caused by the flow at the bottom
!         cell. It follows the erosion method for cohesive particles
!         by Raudkivi 2020.
! --------------------------------------------------------------------
  ! Arguments
  real, intent(in)  :: taub         !< (Pa) Shear stress at bottom
  real, intent(in)  :: tauCrt       !< Critical shear stress for sediment type
  real, intent(in)  :: M_param      !< Surface erosion rate 
  real              :: Beta_ss         !< Dimensionless coefficient for method
  real(kind=8), intent(out) :: resus_flux  !< Vertical erosion flux
  integer, intent(in) :: lwq

  !Beta CAN BE BETWEEN 1 AND 3.6 
  Beta_ss = 1.0
  ! Sediment entrainment flux
  ! Erosion flux for cohesive sediment
  resus_flux = M_param * max(0.0, (taub - tauCrt) / tauCrt) ** Beta_ss
  ! resus_flux = 0.0

  ! if ((l2i(lwq) .ge. 185) .and. (l2j(lwq) .ge. 65)) then
  ! ! if (i .eq. 2) then
  !   resus_flux = resus_flux + 0.3 * M_param * max(0.0, ((taub / (0.5 * tauCrt)) - 1.0))
  !   ! flux_fluff(i) = 0.0
  ! else
  !   resus_flux = resus_flux + 0.3 * M_param * max(0.0, ((taub / (0.5 * tauCrt)) - 1.0))
  !   ! flux_fluff(i) = min(1e-3, 1.5e-6 * (sed_h * c_sed) * max(0.0, (taub / 0.03 - 1.0)))
  !   ! resus_flux = resus_flux + 0.0
  ! end if

  return
END SUBROUTINE resuspension_cohesive

! ********************************************************************
SUBROUTINE deposition_noncohesive(deposition_flux, sett_vel, tauCrt, taub, cb)
! ********************************************************************
!
! Purpose: To estimate the suspended sediment deposition at the
!         bottom layer for each wet column
!
! --------------------------------------------------------------------
  ! Arguments
  real(kind=8), intent(in) :: sett_vel
  real, intent(in) :: taub
  real, intent(in) :: tauCrt
  real, intent(in) :: cb
  real, intent(out) :: deposition_flux

    if (taub .le. tauCrt) then
      deposition_flux = sett_vel * cb
    else
      deposition_flux = 0.0
    end if

  return
END SUBROUTINE deposition_noncohesive

! ********************************************************************
SUBROUTINE deposition_cohesive(deposition_flux, sett_vel, tauCrt, taub, c_ref, ustar, sed_den, z_ref, lwq)
! ********************************************************************
!
! Purpose: To estimate the suspended sediment deposition at the
!         bottom layer for each wet column
!
! --------------------------------------------------------------------
  ! Arguments
  real(kind=8), intent(in) :: sett_vel
  real, intent(in) :: taub
  real, intent(in) :: tauCrt
  real, intent(in) :: c_ref
  real, intent(in) :: ustar
  real, intent(in) :: sed_den
  real, intent(in) :: z_ref
  real, intent(out) :: deposition_flux
  integer, intent(in) :: lwq

  real :: b, zb, cb, c_gel, P

  P = min(3.5, max(2.5, sett_vel / (kappaS * ustar)))
  zb = 0.003
  c_gel = 0.04 * sed_den ! [kg/m3] Concentration at which sediment deposition is zero
  cb = min(c_gel, c_ref * (zb / (z_ref)) ** (-P))

  deposition_flux = (sett_vel * (max(0.0, (1.0 - cb / c_gel)) ** 2)) * (c_ref) * max(0.0, (1.0 - taub / (0.7 * tauCrt)))
  ! deposition_flux = (sett_vel * c_ref) * max(0.0, (1.0 - taub / (0.7 * tauCrt)))
  ! deposition_flux = 0.0

    ! if (taub .le. tauCrt) then
    !   deposition_flux = sett_vel * c_ref * (1.0 - taub/tauCrt)
    ! else
    !   deposition_flux = 0.0
    ! end if

  ! if ((l2i(lwq) .eq. 185) .and. (l2j(lwq) .eq. 80)) then
  !   print*, '-------------- SS MODEL OA04 ----------------'
  !   print*, 'Cb', cb, 'cref', c_ref
  !   print*, 'c_gel', c_gel, 'P', P
  !   print*, 'deposition_flux', deposition_flux
  ! end if
  ! if ((l2i(lwq) .eq. 83) .and. (l2j(lwq) .eq. 134)) then
  !   print*, '-------------- SS MODEL UA06 ----------------'
  !   print*, 'Cb', cb, 'cref', c_ref
  !   print*, 'c_gel', c_gel, 'P', P
  !   print*, 'deposition_flux', deposition_flux
  ! end if
  ! if ((l2i(lwq) .eq. 170) .and. (l2j(lwq) .eq. 47)) then
  !   print*, '-------------- SS MODEL LA03 ----------------'
  !   print*, 'Cb', cb, 'cref', c_ref
  !   print*, 'c_gel', c_gel, 'P', P
  !   print*, 'deposition_flux', deposition_flux
  ! end if

  return
END SUBROUTINE deposition_cohesive

! ********************************************************************
SUBROUTINE get_sed_prop(sett_vel,Rep,tauCrt,sed_d,rho_sed,w_dens)
! ********************************************************************
!
! Purpose: Estimate particle dependent parameters / properties
!
! --------------------------------------------------------------------
  implicit none
  ! Arguments
  real, intent(in)  :: sed_d            !< (m) Sediment diameter
  real, intent(in)  :: w_dens           !< (kg/m3) water density
  real, intent(in)  :: rho_sed         !< (kg/m3) sediment density
  real(kind=8), intent(out) :: sett_vel     !< (m/s) settling velocity
  real, intent(out) :: Rep              !< Explicit Particle Reynolds Number
  real, intent(out) :: tauCrt           !< (Pa) Critical shear stress 
  real              :: submerged_spec_g !< Sediment submerged specific gravity
  logical           :: ivanRijn         !< Flag for using van Rijn (1984) formula or Dietrich (1982). The default is van Rijn

  ivanRijn = .false.

  call submergedSpecificGravity(submerged_spec_g, rho_sed, w_dens)

  call partReynolds_Number(Rep, sed_d, kinematic_viscosity, submerged_spec_g)

  call settling_velocity(sett_vel, g, submerged_spec_g, Rep, sed_d, kinematic_viscosity, ivanRijn)

  call tauCritical(tauCrt, g, sed_d, submerged_spec_g, w_dens, kinematic_viscosity, Rep)

  return
END SUBROUTINE get_sed_prop

! ********************************************************************
SUBROUTINE fvs_ss(vs_ss, sed_d, rho_sed, w_dens)
! ********************************************************************
!
! Purpose: Estimate particle dependent parameters / properties
!
! --------------------------------------------------------------------
  implicit none
  ! Arguments
  real, intent(in)  :: sed_d            !< (m) Sediment diameter
  real, intent(in)  :: w_dens           !< (kg/m3) water density
  real, intent(in)  :: rho_sed          !< (kg/m3) sediment density
  real(kind=8), intent(out) :: vs_ss            !< (m/s) settling velocity
  real              :: Rep              !< Explicit Particle Reynolds Number
  real              :: submerged_spec_g !< Sediment submerged specific gravity
  logical           :: ivanRijn         !< Flag for using van Rijn (1984) formula or Dietrich (1982). The default is van Rijn

  ivanRijn = .false.

  call submergedSpecificGravity(submerged_spec_g, rho_sed, w_dens)

  call partReynolds_Number(Rep, sed_d, kinematic_viscosity, submerged_spec_g)

  call settling_velocity(vs_ss, g, submerged_spec_g, Rep, sed_d, kinematic_viscosity, ivanRijn)

  return
END SUBROUTINE fvs_ss

! ********************************************************************
SUBROUTINE submergedSpecificGravity(submerged_spec_g, rho_sed, w_dens)
! ********************************************************************
!
! Purpose: To estimate the submerged specific gravity for a given
!          Particle type
!
! --------------------------------------------------------------------
  ! Arguments of subroutine
  real, intent(in)  :: rho_sed         !< (mg/m3) Sediment density
  real, intent(in)  :: w_dens           !< (mg/m3) Water density
  real, intent(out) :: submerged_spec_g !< Sediment submerged specific gravity

  ! Estimate submerged specific gravity
  submerged_spec_g = (rho_sed / w_dens) - 1

  return
END SUBROUTINE submergedSpecificGravity

! ********************************************************************
SUBROUTINE partReynolds_Number(Rep, sed_d, ki_visc, submerged_spec_g)
! ********************************************************************
!
! Purpose: To estimate the explicit Particle Reynolds Number.
!        The equation used is from Garcia and Parker 1993.
!        Experiments on the entrainment of sediment into suspension
!        by a dense bottom current. DOI: 10.1029/92JC02404
!
! --------------------------------------------------------------------

  ! Arguments of subroutine
  real, intent(in)  :: sed_d                !< (m) Sediment diameter D50
  real, intent(in)  :: ki_visc              !< (m2/sec) kinematic viscosity of water
  real, intent(in)  :: submerged_spec_g     !< Sediment submerged specific gravity
  real, intent(out) :: Rep                  !< Explicit Particle Reynolds Number

  Rep = sqrt(g * submerged_spec_g * sed_d ** 3) / ki_visc

  return
END SUBROUTINE partReynolds_Number

! ********************************************************************
SUBROUTINE settling_velocity(sett_vel, g_ss, submerged_spec_g, Rep, sed_d, ki_visc, ivanRijn)
! ********************************************************************
!
! Purpose: To estimate the settling velocity for a given particle
!          size. The estimate uses 
!
! --------------------------------------------------------------------
  ! Arguments of subroutine
  real, intent(in)  :: g_ss                    !< (m/s2)Gravitational acceleration (m/s**2)
  real, intent(in)  :: submerged_spec_g     !< Sediment submerged specific gravity
  real, intent(in)  :: ki_visc              !< (m2/s) Kinematic viscosity of water
  real, intent(in)  :: sed_d                !< (m) Sediment Diameter
  ! real, intent(in)  :: sed_spec_g           !< Sediment specific gravity
  real, intent(in)  :: Rep                  !< Explicit Particle Reynolds Number
  logical, optional :: ivanRijn             !< Flag for using van Rijn (1984) formula or Dietrich (1982). The default is van Rijn      
  real              :: dimless_fall_vel     !< dimensionaless fall velocity
  logical           :: vanRijnFlag
  integer           :: i
  ! Parameter for Dietrich (1982) equation
  ! Values are from dsm2
  real              :: b_1ss = 3.76715
  real              :: b_2ss = 1.92944 
  real              :: b_3ss = 0.09815 
  real              :: b_4ss = 0.00575
  real              :: b_5ss = 0.00056
  ! values are from Bombardelli and Moreno 2012 found in Reardon et al., 2014
  ! real              :: b_1 = 2.891394
  ! real              :: b_2 = 0.95296 
  ! real              :: b_3 = 0.056835 
  ! real              :: b_4 = 0.002892
  ! real              :: b_5 = 0.000245 
  real(kind=8), intent(out) :: sett_vel         !< (m/s) Settling

  if ( present(ivanRijn) ) then
    vanRijnFlag = ivanRijn
  end if

  SELECT CASE (vanRijnFlag)
    CASE (.true.)
      ! Van Rijn Formula
      if (sed_d .gt. 1.0d-3) then
        sett_vel = 1.1 * sqrt(submerged_spec_g * g_ss * sed_d)
      elseif (sed_d .gt. 1.0d-4 .and. sed_d .le. 1.0d-3) then
        sett_vel = (10 * ki_visc / sed_d) *        & 
                       (sqrt(1 + 0.01 * (submerged_spec_g * g_ss       &
                        * sed_d **3) / ki_visc ** 2.) - 1)
      else
      ! Stokes Law
        sett_vel = (submerged_spec_g * g_ss * sed_d ** 2.) / (18.0 * ki_visc)
      end if

    CASE (.false.)
      dimless_fall_vel = exp(-1.*b_1ss + b_2ss * log(Rep) - b_3ss * (log(Rep)) ** 2.0 - b_4ss * (log(Rep)) ** 3. + b_5ss * (log(Rep)) ** 4.)
      ! if ( sed_diamm .lt. 1.0d-5) then
      !   sett_vel = (submerged_spec_g * g * sed_diamm**2)/(18.*ki_visc)
      ! else
      sett_vel = dimless_fall_vel * sqrt(submerged_spec_g * g_ss * sed_d)
      ! end if
  END SELECT

  return
END SUBROUTINE settling_velocity

! ********************************************************************
SUBROUTINE tauCritical(tauCrt, g_ss, sed_d, submerged_spec_g, w_dens, ki_visc, Rep)
! ********************************************************************
!
! Purpose: To estimate the critical shear stress for a given particle
!          size. The estimate uses the nondimensional critical shields
!          parameter from Parker et al., 2003.
!
! --------------------------------------------------------------------
  ! Arguments of subroutine
  real, intent(in)  :: g_ss                   !< (m/s2) Gravitational acceleration (m/s**2)
  real, intent(in)  :: sed_d               !< (m) Sediment Diameter
  real, intent(in)  :: submerged_spec_g    !< Sediment submerged specific gravity
  real, intent(in)  :: w_dens              !< (kg/m3) Water density
  real, intent(in)  :: ki_visc             !< (m2/s) Kinematic viscosity of water
  real, intent(in)  :: Rep                 !< Explicit Particle Reynolds Number
  real              :: shields_param       !< Nondimensional Critical Shields Parameter
  real, intent(out) :: tauCrt              !< (Pa) Critical shear stress

  ! Estimate of the nondimensional critical Shields parameter
  shields_param = 0.3 * (0.22 * Rep ** (-0.6) + 0.06 * 10 ** (-7.7 * Rep ** (-0.6)))

  ! Estimate of critical shear stress for given water and sediment properties
  tauCrt = shields_param * g_ss * submerged_spec_g * sed_d * (w_dens)
  ![kgm/s2/m2] = [-]     * [m/s2]  *      [-]      *   [m] * [kg/m3]
  return 
END SUBROUTINE tauCritical

! ********************************************************************
SUBROUTINE tauBottom(taub, ustarb,kwq,lwq)
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
  integer               :: kms    !< bottom layer
  real                  :: ubott  !< (m/s) vel u at bottom cell
  real                  :: vbott  !< (m/s) vel v at bottom cell
  real                  :: taubx  !< (Pa) bottom shear stress in x
  real                  :: tauby  !< (Pa) bottom shear stress in y
  ! real, intent(in)      :: tauwbx !< (Pa) bottom shear stress induced by waves
  ! real, intent(in)      :: tauwby !< (Pa) bottom shear stress induced by waves  
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

    ! ....Compute Currend-Induced Bottom Shear Stress
    ubott = ( uhp( kmxp, lwq ) + uhp( kmxm, lWC(lwq)) ) / 2. / hp(kms,lwq)
    vbott = ( vhp( kmyp, lwq ) + vhp( kmym, lSC(lwq)) )/2. / hp(kms,lwq)
    taubx = cd * (rhop(kms,lwq)+1000.) * ubott**2.
    tauby = cd * (rhop(kms,lwq)+1000.) * vbott**2.

    ! ... Add Wave-Induced Bottom Shear Stress & calculate friction velocity
    if (iSTWAVE == 1) then
      taub  = sqrt( (sqrt((taubx)**2. + (tauby)**2.)) ** 2 + tau_stwave(l2i(lwq), l2j(lwq)) ** 2)
    else
      taub  = sqrt((taubx)**2. + (tauby)**2.)
    end if
    ustarb = sqrt(taub/(rhop(kms, lwq) + 1000.))
  endif

END SUBROUTINE tauBottom

!************************************************************************
SUBROUTINE burial(burial_flux, resus_flux, deposition_flux)
!************************************************************************
!
!   Purpose: To estimate erosion of MeHg adsorbed to sediments
!
!
!------------------------------------------------------------------------

  ! Arguments
  real, intent(in)  :: deposition_flux
  real(kind=8), intent(in)  :: resus_flux
  real, intent(out) :: burial_flux

  burial_flux = deposition_flux - resus_flux

END SUBROUTINE burial

!************************************************************************
                        END MODULE si3d_sed
!************************************************************************