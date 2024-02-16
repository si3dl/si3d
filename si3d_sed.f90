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
  ! real, dimension(sedNumber) :: settling_vel      !< (m/s) Settling velocity of sediment
  real, dimension(sedNumber) :: Rep               !< Explicity Particle Reynolds Number
  real, dimension(sedNumber) :: tauCrt            !< Critical Shear Stress
  real, dimension(sedNumber) :: erosion_flux      !<
  real, dimension(sedNumber) :: deposition_flux    !<
  real                       :: cb

  kms = kmz(lwq)
  w_dens = rhop(kwq, lwq) + 1000

  if (kwq .eq. kms) then
    ! Estimate bottom shear stress
    call tauBottom(taub, ustarb, kwq, lwq)
    call totalMass(sed_frac, kwq, lwq)

    do i = 1, sedNumber
      ! Estimate properties of sediment for a given water density at bottom cell
      cb = tracerpp(kwq,lwq,LSS1 + i - 1)
      
      call get_sed_prop(settling_vel(i),Rep(i),tauCrt(i),sed_diameter(i),sed_dens(i),w_dens)
      ! Estimate erosion flux
      if (taub .ge. tauCrt(i)) then
        if (sed_type(i) == 0) then
          call erosion_noncohesive(erosion_flux(i), ustarb, Rep(i), settling_vel(i), sed_frac(i),sed_dens(i))
          if (lwq == 50) then
            print*,'Erosion Flux for noncohesive'
            print*,'diameter = ',sed_diameter(i),'fraction = ', sed_frac(i)
          end if
        else if (sed_type(i) == 1) then
          call erosion_cohesive(erosion_flux(i),taub, tauCrt(i), sed_frac(i), sed_dens(i))
          if (lwq == 50) then
            print*,'Erosion Flux for cohesive'
            print*,'diameter = ',sed_diameter(i),'fraction = ', sed_frac(i)
          end if
        end if 
      else
        erosion_flux(i) = 0.0
      end if

      ! Estimate deposition flux
      if (cb .gt. 0.0) then
        if (sed_type(i) == 0) then
          call deposition_noncohesive(deposition_flux(i), settling_vel(i), tauCrt(i), taub, cb)
          if (lwq == 50) then
            print*,'Deposition Flux for noncohesive'
            print*,'diameter = ',sed_diameter(i),'fraction = ', sed_frac(i)
          end if
        else if (sed_type(i) == 1) then
          call deposition_cohesive(deposition_flux(i), settling_vel(i), tauCrt(i), taub, cb)
          if (lwq == 50) then
            print*,'Deposition Flux for cohesive'
            print*,'diameter = ',sed_diameter(i),'fraction = ', sed_frac(i)
          end if
        end if
        ! To correct deposition flux as it can not remove more sediment than what is in the water layer on top of sediment
        if (deposition_flux(i) .gt. (cb * hp(kwq, lwq) / dt)) then
          deposition_flux(i) = cb * hp(kwq, lwq) / dt
        end if
      else
        deposition_flux(i) = 0.0
      end if

      sourcesink(kwq, lwq, LSS1 + i - 1) = erosion_flux(i) - deposition_flux(i)

      ! if ((tauCrt(i) .lt. taub)) then
      if (lwq == 50) then
        print*, '--------------------------'
        print*, 'k =',kwq,'l =',lwq,'i = ',l2i(lwq),'j = ',l2j(lwq)
        print*, 'dt = ',dt,'h = ',hp(kwq,lwq)
        print*, 'cb = ',cb
        print*, 'taub =',taub
        print*, 'tauCr =', tauCrt(i)
        print*, 'erosion_flux = ', erosion_flux(i)
        print*, 'depositionFlux = ', deposition_flux(i)
        print*, 'sourcesink = ', sourcesink(kwq, lwq, LSS1+i-1)
        print*, 'sourcesinklim = ', -cb * hp(kwq,lwq) / dt
        print*, 'Rep = ', Rep(i)
      end if

      if (sourcesink(kwq,lwq,LSS1 + i - 1) .lt. (-1*cb * hp(kwq,lwq) / dt)) then
        sourcesink(kwq, lwq, LSS1 + i - 1) = -1 * cb * hp(kwq,lwq) / dt
      end if

      ! if (tauCrt(i) .lt. taub) then
      if (lwq == 50) then
        print*, 'sourcesinkNEW = ', sourcesink(kwq, lwq, LSS1+i-1)        
        print*, 'sourcesinklim = ', -cb * hp(kwq,lwq) / dt 
        print*, 'cb in/out = ', sourcesink(kwq, lwq, LSS1+i-1) * dt / hp(kwq,lwq)
      end if

      ! Estimate source and sink for the sediment cell.
      sourcesink(kwq+1,lwq,LSS1 + i - 1) = deposition_flux(i) - erosion_flux(i)
      ! Estimate erosion for mercury
      erosion_Hgpn(i) = erosion_flux(i) / (sed_dens(i) * sed_frac(i))
    end do

  elseif (kwq .ne. kms) then 
    do i = 1, sedNumber
      sourcesink(kwq,lwq,LSS1 + i - 1) = 0.0
      erosion_Hgpn(i) = 0.0
    end do
  end if

END SUBROUTINE sourceSS

! ********************************************************************
SUBROUTINE totalMass(sed_frac, kwq, lwq)
! ********************************************************************
!
! Purpose: To estimate the total mass in the sediment layer and
!         estimate the sediment fraction in the bed for each type
!         of sediment.
!
! --------------------------------------------------------------------
  ! Arguments
  real, intent(inout), dimension(sedNumber) :: sed_frac
  integer, intent(in)                       :: kwq
  integer, intent(in)                       :: lwq
  real , dimension(sedNumber)               :: massSed
  real                                      :: totalMass_bed
  integer                                   :: i

  totalMass_bed = 0.0
  do i = 1, sedNumber
    massSed(i) = tracerpp(kwq+1,lwq,LSS1 + i - 1) * dx * dy * h(kwq,lwq)
    if (lwq == 50) then
    print*,'i = ',i,'h = ',h(kwq,lwq),'dx = ',dx,'dy = ',dy
    print*,'LSS1 + i -1 = ',LSS1 + i - 1
    print*,'tracer LSS1 + i -1 = ', tracerpp(kwq+1,lwq,LSS1 + i - 1)
    print*,'massSed = ',massSed(i)
    end if 
  end do
  totalMass_bed = sum(massSed)
  do i = 1, sedNumber
    sed_frac(i) = massSed(i) / totalMass_bed
    if (lwq == 50) then
    print*,'total_mass = ',totalMass_bed
    print*,'i = ',i,' sed_frac(i) = ',sed_frac(i)
    end if 
  end do

END SUBROUTINE totalMass

! ********************************************************************
SUBROUTINE erosion_noncohesive(erosion_flux, ustarb, Rep, settling_vel, sed_frac, sed_dens)
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
  real, intent(in)  :: settling_vel !< (m/s) Settling velocity of sediment
  real, intent(in)  :: sed_frac     !< fraction of bed of a given sediment type
  real, intent(in)  :: sed_dens     !< Sediment density to keep units consistent
  real              :: z_u          !< Similarity variable for uniform sediment
  real              :: E_s          !< Dimensionless coefficient for sediment entrainment. Under quasi-equilibrium conditions (Garcia & Parker 1991)
  real, intent(out) :: erosion_flux  !< Vertical erosion flux

  if ((Rep .gt. 0.4) .and. (Rep .le. 1)) then
    z_u = 1 * ustarb * (Rep ** 3.75) / settling_vel
  else if ((Rep .gt. 1.0) .and. (Rep .le. 3.5)) then
    z_u = 0.586 * ustarb * (Rep ** 1.23) / settling_vel
  else if (Rep .gt. 3.5) then
    z_u = 1 * ustarb * (Rep ** 0.6) / settling_vel
  end if

  ! Sediment entrainment coefficient
  E_s = Ased * (z_u ** 5) / (1 + (z_u ** 5) * Ased/0.3)

  erosion_flux = E_s * settling_vel * sed_frac * sed_dens

  return
END SUBROUTINE erosion_noncohesive

! ********************************************************************
SUBROUTINE erosion_cohesive(erosion_flux, taub, tauCrt, sed_frac, sed_dens)
! ********************************************************************
!
! Purpose: To estimate the erosion caused by the flow at the bottom
!         cell. It follows the erosion method for cohesive particles
!         by Raudkivi 2020.
! --------------------------------------------------------------------
  ! Arguments
  real, intent(in)  :: taub         !< (Pa) Shear stress at bottom
  real, intent(in)  :: tauCrt       !< Critical shear stress for sediment type
  real, intent(in)  :: sed_frac     !< fraction of bed of a given sediment type
  real, intent(in)  :: sed_dens     !< Sediment density to keep units consistent
  real              :: M_param      !< Surface erosion rate 
  real              :: Beta         !< Dimensionless coefficient for method
  real, intent(out) :: erosion_flux  !< Vertical erosion flux

  !Beta CAN BE BETWEEN 1 AND 3.6 
  M_param = 3.6e-4

  print*, 'M =', M_param
  Beta = 3
  ! Sediment entrainment flux
  erosion_flux = M_param * (taub / tauCrt - 1) ** Beta * sed_frac * sed_dens

  return
END SUBROUTINE erosion_cohesive

! ********************************************************************
SUBROUTINE deposition_noncohesive(deposition_flux, settling_vel, tauCrt, taub, cb)
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
  real, intent(in) :: cb
  real, intent(out) :: deposition_flux

    if (taub .le. tauCrt) then
      deposition_flux = settling_vel * cb * (1 - taub/tauCrt)
    else
      deposition_flux = settling_vel * cb
    end if

  return
END SUBROUTINE deposition_noncohesive

! ********************************************************************
SUBROUTINE deposition_cohesive(deposition_flux, settling_vel, tauCrt, taub, cb)
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
  real, intent(in) :: cb
  real, intent(out) :: deposition_flux

    if (taub .le. tauCrt) then
      deposition_flux = settling_vel * cb * (1 - taub/tauCrt)
    else
      deposition_flux = 0
    end if

  return
END SUBROUTINE deposition_cohesive

! ********************************************************************
SUBROUTINE get_sed_prop(settling_vel,Rep,tauCrt,sed_diameter,sed_dens,w_dens)
! ********************************************************************
!
! Purpose: Estimate particle dependent parameters / properties
!
! --------------------------------------------------------------------
  implicit none
  ! Arguments
  real, intent(in)  :: sed_diameter     !< (um) Sediment diameter
  real, intent(in)  :: w_dens           !< (kg/m3) water density
  real, intent(in)  :: sed_dens         !< (kg/m3) sediment density
  real, intent(out) :: settling_vel     !< (m/s) settling velocity
  real, intent(out) :: Rep              !< Explicit Particle Reynolds Number
  real, intent(out) :: tauCrt           !< (Pa) Critical shear stress 
  real              :: submerged_spec_g !< Sediment submerged specific gravity
  real              :: sed_diamm        !< (m) Sediment diameter
  logical           :: ivanRijn         !< Flag for using van Rijn (1984) formula or Dietrich (1982). The default is van Rijn

  ivanRijn = .false.

  ! To estimate sediment diamenter in m
  sed_diamm = sed_diameter * 0.000001

  ! sed_spec_g = sed_dens/w_dens

  call submergedSpecificGravity(submerged_spec_g, sed_dens, w_dens)

  call partReynolds_Number(Rep, sed_diamm, kinematic_viscosity, submerged_spec_g)

  call settling_velocity(settling_vel, g, submerged_spec_g, Rep, sed_diamm, kinematic_viscosity, ivanRijn)

  call tauCritical(tauCrt,g,sed_diamm,submerged_spec_g, sed_dens,kinematic_viscosity, Rep)

  return
END SUBROUTINE get_sed_prop

! ********************************************************************
SUBROUTINE fvs_ss(vs_ss,sed_diameter,sed_dens,w_dens)
! ********************************************************************
!
! Purpose: Estimate particle dependent parameters / properties
!
! --------------------------------------------------------------------
  implicit none
  ! Arguments
  real, intent(in)  :: sed_diameter     !< (um) Sediment diameter
  real, intent(in)  :: w_dens           !< (kg/m3) water density
  real, intent(in)  :: sed_dens         !< (kg/m3) sediment density
  real, intent(out) :: vs_ss            !< (m/s) settling velocity
  real              :: Rep              !< Explicit Particle Reynolds Number
  real              :: submerged_spec_g !< Sediment submerged specific gravity
  real              :: sed_diamm        !< (m) Sediment diameter
  logical           :: ivanRijn         !< Flag for using van Rijn (1984) formula or Dietrich (1982). The default is van Rijn

  ivanRijn = .true.

  ! To estimate sediment diamenter in m
  sed_diamm = sed_diameter * 0.000001

  ! sed_spec_g = sed_dens/w_dens

  call submergedSpecificGravity(submerged_spec_g, sed_dens, w_dens)

  call partReynolds_Number(Rep, sed_diamm, kinematic_viscosity, submerged_spec_g)

  call settling_velocity(vs_ss, g, submerged_spec_g, Rep, sed_diamm, kinematic_viscosity, ivanRijn)

  return
END SUBROUTINE fvs_ss

! ********************************************************************
SUBROUTINE submergedSpecificGravity(submerged_spec_g, sed_dens, w_dens)
! ********************************************************************
!
! Purpose: To estimate the submerged specific gravity for a given
!          Particle type
!
! --------------------------------------------------------------------
  ! Arguments of subroutine
  real, intent(in)  :: sed_dens         !< (kg/m3) Sediment density
  real, intent(in)  :: w_dens           !< (kg/m3) Water density
  real, intent(out) :: submerged_spec_g !< Sediment submerged specific gravity

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
  real, intent(in)  :: sed_diamm            !< (m) Sediment diameter D50
  real, intent(in)  :: kinematic_viscosity  !< (m2/sec) kinematic viscosity of water
  real, intent(in)  :: submerged_spec_g     !< Sediment submerged specific gravity
  real, intent(out) :: Rep                  !< Explicit Particle Reynolds Number

  Rep = sqrt(g * submerged_spec_g * sed_diamm ** 3) / kinematic_viscosity

  return
END SUBROUTINE partReynolds_Number

! ********************************************************************
SUBROUTINE settling_velocity(settling_vel, g, submerged_spec_g, Rep, sed_diamm, kinematic_viscosity, ivanRijn)
! ********************************************************************
!
! Purpose: To estimate the settling velocity for a given particle
!          size. The estimate uses 
!
! --------------------------------------------------------------------
  ! Arguments of subroutine
  real, intent(in)  :: g                    !< (m/s2)Gravitational acceleration (m/s**2)
  real, intent(in)  :: submerged_spec_g     !< Sediment submerged specific gravity
  real, intent(in)  :: kinematic_viscosity  !< (m2/s) Kinematic viscosity of water
  real, intent(in)  :: sed_diamm            !< (m) Sediment Diameter
  ! real, intent(in)  :: sed_spec_g           !< Sediment specific gravity
  real, intent(in)  :: Rep                  !< Explicit Particle Reynolds Number
  logical, optional :: ivanRijn             !< Flag for using van Rijn (1984) formula or Dietrich (1982). The default is van Rijn      
  real              :: dimless_fall_vel     !< dimensionaless fall velocity
  logical           :: vanRijnFlag
  integer           :: i
  ! Parameter for Dietrich (1982) equation
  ! Values are from dsm2
  real              :: b_1 = 3.76715
  real              :: b_2 = 1.92944 
  real              :: b_3 = 0.09815 
  real              :: b_4 = 0.00575
  real              :: b_5 = 0.00056
  ! values are from Bombardelli and Moreno 2012 found in Reardon et al., 2014
  ! real              :: b_1 = 2.891394
  ! real              :: b_2 = 0.95296 
  ! real              :: b_3 = 0.056835 
  ! real              :: b_4 = 0.002892
  ! real              :: b_5 = 0.000245 
  real, intent(out) :: settling_vel         !< (m/s) Settling

  if ( present(ivanRijn) ) then
    vanRijnFlag = ivanRijn
  end if

  do i = 1,sedNumber

    SELECT CASE (vanRijnFlag)
      CASE (.true.)
        ! Van Rijn Formula
        if (sed_diamm .gt. 1.0d-3) then
          settling_vel = 1.1 * sqrt(submerged_spec_g * g * sed_diamm)
        elseif (sed_diamm .gt. 1.0d-4 .and. sed_diamm .le. 1.0d-3) then
          settling_vel = (10 * kinematic_viscosity / sed_diamm) *        & 
                         (sqrt(1 + 0.01 * (submerged_spec_g * g       &
                          * sed_diamm **3) / kinematic_viscosity ** 2.) - 1)
        ! elseif (sed_diamm .gt. 6.5d-5 .and. sed_diamm .le. 1.0d-4) then
        !   ! Stokes Law
        !   settling_vel = (submerged_spec_g * g * sed_diamm ** 2.) / (18.0 * kinematic_viscosity)
        else
        ! Stokes Law
          settling_vel = (submerged_spec_g * g * sed_diamm ** 2.) / (18.0 * kinematic_viscosity)
        end if

      CASE (.false.)
        dimless_fall_vel = exp(-1.*b_1 + b_2 * log(Rep) - b_3 * (log(Rep)) ** 2.0 - b_4 * (log(Rep)) ** 3. + b_5 * (log(Rep)) ** 4.)
        ! if ( sed_diamm .lt. 1.0d-5) then
        !   settling_vel = (submerged_spec_g * g * sed_diamm**2)/(18.*kinematic_viscosity)
        ! else
        settling_vel = dimless_fall_vel * sqrt(submerged_spec_g * g * sed_diamm)
        ! end if
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
  real, intent(in)  :: g                   !< (m/s2) Gravitational acceleration (m/s**2)
  real, intent(in)  :: sed_diamm           !< (m) Sediment Diameter
  real, intent(in)  :: submerged_spec_g    !< Sediment submerged specific gravity
  real, intent(in)  :: sed_dens            !< (kg/m3) Sediment density
  real, intent(in)  :: kinematic_viscosity !< (m2/s) Kinematic viscosity of water
  real, intent(in)  :: Rep                 !< Explicit Particle Reynolds Number
  real              :: shields_param       !< Nondimensional Critical Shields Parameter
  real, intent(out) :: tauCrt              !< (Pa) Critical shear stress

  ! Estimate of the nondimensional critical Shields parameter
  shields_param = 0.5 * (0.22 * Rep ** (-0.6) + 0.06 * 10 ** (-7.7 * Rep ** (-0.6)))

  ! Estimate of critical shear stress for given water and sediment properties
  tauCrt = shields_param * g * submerged_spec_g * sed_diamm * sed_dens

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
      taub  = sqrt((taubx)**2. + (tauby)**2.) + tau_stwave(l2i(lwq),l2j(lwq))
    else
      taub  = sqrt((taubx)**2. + (tauby)**2.)
    end if
    ustarb = sqrt(taub/1000.)
  endif

END SUBROUTINE tauBottom

!************************************************************************
                        END MODULE si3d_sed
!************************************************************************