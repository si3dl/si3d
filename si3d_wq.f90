!************************************************************************
                          MODULE si3d_wq
!************************************************************************
!
!  Purpose: Procedures that implement routines related to the modelling
!           of water quality in lakes
!
!-------------------------------------------------------------------------

 USE si3d_types
 USE si3d_sed
 USE si3d_Hg
 IMPLICIT NONE
 SAVE

CONTAINS

!*********************************************************************
SUBROUTINE sourceDO(kwq,lwq) 
!********************************************************************
!
!  Purpose: if dissolved oxygen is modeled, this subroutine
!  calculates source and sink terms that depend on dissolved
!  oxygen concentrations
!  Important: Sinksource term must have flux units [M/L^2/T], or
!             [M/L^3]*L*(1/T)
!  Units of DO in the input file are mg/m3
!
!-----------------------------------------------------------------------

  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq, lwq

  !. . . Local Variables
  REAL    ::   Tk, lnOS, OS, Patm_do, ln_Pwv, Pwv, theta2, f_SOD 
  REAL    :: reaeration, sedoxydemand
  real    :: ws, R_SOD_ij
  integer :: i, j

  reaeration = 0.0
  sedoxydemand = 0.0

  ! ...Calculate DO saturation
  Tk = salp(kwq, lwq) + 273
  lnOS = -139.34410 + 1.575701*1E5 /(Tk    ) &
  &                 - 6.642308*1E7 /(Tk**2.) &
  &                 + 1.243800*1E10/(Tk**3.) &
  &                 - 8.621949*1E11/(Tk**4.)
  OS = EXP(lnos)

  ! Correct for Patmospheric (Pa - declared in si3d_types and defined in surfbc0)
  Patm_do   = Pa * 0.00000986923; ! Transform atmospheric pressure from Pa to atm
  ln_Pwv = 11.8751 - (3840.70/Tk) - (216961/(Tk**2.))
  Pwv    = EXP(ln_Pwv)
  theta2 = 0.000975 - 1.426*1E-5 * salp(kwq,lwq) + &
  &                    6.436*1E-8 * salp(kwq,lwq)**2.
  OS = OS*Patm_do*((1-Pwv/Patm_do) *(1-theta2*Patm_do))&
  &           /((1-Pwv)*(1-theta2) )
  ! Estimate of OS is in mg/L. Si3D uses mg/m3 then:
  OS = OS * 1000

  ! ...Calculate reaeration only at the lake surface
  ! for now using constant reaeration defined in wq_inp, but in future, can have
  ! alternatives for reaeration rates.
  IF (kwq .le. k1z(lwq) + 1) THEN
    ws = SQRT(uair(lwq)**2. + vair(lwq)**2.)
    ! if (ws .le. 0.6) then
    !   ws = 0.6
    ! end if
    ! if (kwq .eq. k1z(lwq)) then
      reaeration = R_reaer * (ws ** 1.64)* (OS - tracerpp(kwq,lwq,LDO))
    ! else
    !   reaeration = 1/2 * (R_reaer * (ws ** 1.64))* (OS - tracerpp(kwq,lwq,LDO))
    ! end if

    ! Units: [mg/m^2/s] = [m/s] * [mg/m^3]
  ELSE
     reaeration  = 0.0
  END IF

  ! ...Calculate the sediment oxygen demand from the sediments (only bottom cell)
  IF (kwq .eq. kmz(lwq)) THEN
    f_SOD = tracerpp(kwq,lwq,LDO) /(KSOD + tracerpp(kwq,lwq,LDO) ) ! DO inhibition of sediment oxygen demand
    i = l2i(lwq)
    j = l2j(lwq)
    if (((i >= 1) .and. (i <= 139)) .and. ((j >=1) .and. (j <= 195))) then
      R_SOD_ij = R_SOD * 0.45
    elseif ((i > 170) .and. ((j >= 1) .and. (j <= 63))) then
      R_SOD_ij = R_SOD * 0.05
    elseif (((i > 139) .and. (i <= 180)) .and. ((j > 63) .and. (j <= 70))) then
      R_SOD_ij = R_SOD * 0.05
    else
      R_SOD_ij = R_SOD
    end if
    ! Units of KSDO need to be mg/m3
    sedoxydemand = R_SOD_ij * f_SOD * (Theta_SOD ** (salp(kwq, lwq) - 20)) 
    ! Units: [mg/m^2/s] = [mg/m^2/s] * [-] * [-] 
  ELSE
     sedoxydemand = 0.0

  END IF

  ! ... Incorporate all terms to the source-sink DO term
   sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO)    &
              & +  reaeration                           &
              & -  sedoxydemand

  fluxes_out(kwq, lwq, 2) = reaeration
  fluxes_out(kwq, lwq, 3) = sedoxydemand
  fluxes_out(kwq + 1, lwq, 3) = sedoxydemand

  ! If ALG are modeled: Add photosynthetic oxygen production and respiration consumption  by phytoplankton
  ! If NH4 is modeled: Add consumption  by nitrification
  ! If DOC is modeled: Add consumption  by mineratlization of DOC (water column biological oxygen demand)

END SUBROUTINE sourceDO

!************************************************************************
SUBROUTINE sourcePON(kwq,lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on PON
!  Important: Sinksource term must have flux units [M/L^2/T], or
!             [M/L^3]*L*(1/T)
!--------------------------------------------------------------------------

  ! ... Arguments
    INTEGER, INTENT (IN) :: kwq,lwq

  !... Local variables
  REAL:: decompositionPON, f_decom, depositionPON, resuspensionPON

  depositionPON = 0.0
  resuspensionPON = 0.0
  !... Calculate decompositionN
    ! Calculate DO inhibition of decomposition
    IF (IDO == 1) THEN
      f_decom = tracerpp(kwq,lwq,LDO) /(KDECMIN + tracerpp(kwq,lwq,LDO) )  ! We assume that decomposition and mineralzation 
                                                                           ! half-saturation values are very similar
    ELSE
      f_decom = 1.0
    END IF
  decompositionPON = R_decom_pon * f_decom * (Theta_decom**(salp(kwq,lwq) - 20.0)) *tracerpp(kwq,lwq,LPON) * hpp(kwq,lwq)
  ! Units: [mg/m^2/s] = [1/s] * [-] * [-] * [mg/m^3] * [m]

  ! ... Calculate deposition of PON only in the bottom layer
  IF (kwq .eq. kmz(lwq)) THEN
    depositionPON = R_settl * tracerpp(kwq,lwq,LPON) 
    ! Units: [mg/m^2/s] =  [m/s] * [mg/m^3] 
    ! ... Calculate resusupension of PON only in the bottom layer
    resuspensionPON = R_resusp * tracerpp(kwq,lwq,LPON)
    ! Units: [mg/m^2/s] =  [m/s] * [mg/m^3]

    ! resuspensionPON = R_resusp * tracerpp(kwq,lwq,LPON) 
    ! Units: [mg/m^2/s] =  [m/s] * [mg/m^3]
    if (depositionPON .gt. (tracerpp(kwq, lwq, LPON) * hp(kwq, lwq) / dt)) then
      depositionPON = tracerpp(kwq, lwq, LPON) * hp(kwq, lwq) / dt
    end if
  END IF

  ! ... Incorporate all terms to the source-sink PON term
  sourcesink(kwq,lwq,LPON) = sourcesink(kwq,lwq,LPON)       &
                        &    - decompositionPON             &   
                        &    - depositionPON                  &                   
                        &    + resuspensionPON
                        !    + mortality  - only if IALG = 1; calcualted in sourceALG                         

  ! Add contribution of decomposition to DON concentration
  IF (iDON == 1) THEN
    sourcesink(kwq,lwq,LDON) = sourcesink(kwq,lwq,LDON) +  decompositionPON
  END IF

END SUBROUTINE sourcePON

!************************************************************************
SUBROUTINE sourceDON(kwq,lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on DON
!  Important: Sinksource term must have flux units [M/L^2/T], or
!             [M/L^3]*L*(1/T)
!--------------------------------------------------------------------------

  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq

  !. . . Local variables
  REAL:: mineralizationDON, f_miner, atmosdepositionDON, f_sedflux, sedfluxDON

  !. . Mineralization by bacteria
    ! Calculate DO inhibition of mineralization
    IF (IDO == 1) THEN
      f_miner = tracerpp(kwq,lwq,LDO) /(KDECMIN + tracerpp(kwq,lwq,LDO) )
    ELSE
      f_miner = 1.0
    END IF
  mineralizationDON = R_miner_don * f_miner * (Theta_miner**(salp(kwq,lwq) - 20.0)) *tracerpp(kwq,lwq,LDON) * hpp(kwq,lwq)
  ! Units: [mg/m^2/s] = [1/s] * [-] * [-] * [mg/m^3] * [m]

  !. . .Add contribution from atmospheric deposition to top layer
  IF (kwq .eq. k1z(lwq)) THEN
    atmosdepositionDON = ATM_DON
    ! Units: [mg/m^2/s] =  [mg/m^2/s] 
  ELSE
    atmosdepositionDON = 0.0
  END IF

  !. . .Add contribution from sediment flux to the bottom layer
  IF (kwq .eq. kmz(lwq)) THEN
      IF (IDO == 1) THEN
         ! Calculate DO inhibition of sediment flux
         f_sedflux = tracerpp(kwq,lwq,LDO) /(KSED + tracerpp(kwq,lwq,LDO) )
         ! f_sedflux = KSED / (KSED + tracerpp(kwq,lwq,LDO))
      ELSE
         f_sedflux = 1.0
      END IF
      sedfluxDON = SED_DON * f_sedflux* (Theta_sedflux**(salp(kwq,lwq) - 20))  
      ! Units: [mg/m^2/s] = [mg/m^2/s] * [-] * [-] 
  ELSE
      sedfluxDON = 0.0
  END IF

  ! ... Incorporate all terms to the source-sink DON term
  sourcesink(kwq,lwq,LDON) = sourcesink(kwq,lwq,LDON)         &
                         &        -  mineralizationDON        &
                         &        + atmosdepositionDON        &
                         &        + sedfluxDON
                        !         + decompositionPON - only if IPON = 1; caluclated in sourcePON

  ! Add mineralization of DON to NH4
  IF (INH4 == 1) THEN   
    sourcesink(kwq,lwq,LNH4) = sourcesink(kwq,lwq,LNH4) + mineralizationDON
  END IF

END SUBROUTINE sourceDON

!************************************************************************
SUBROUTINE sourceNH4(kwq,lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on NH4
!  Important: Sinksource term must have flux units [M/L^2/T], or
!             [M/L^3]*L*(1/T)
!--------------------------------------------------------------------------

  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq

  !. . . Local Variables
  REAL:: nitrification, f_nitrif, atmosdepositionNH4, sedfluxNH4, f_sedflux

  ! ... Eliminate negative values of NH4 (min detection)
  IF (tracerpp(kwq,lwq,LNH4) .lt. 0) THEN
     tracerpp(kwq,lwq,LNH4) = 1.000
  END IF

  !. . . Calculate nitrification
    ! Calculate DO inhibition of nitrification
    IF (IDO == 1) THEN
      f_nitrif = tracerpp(kwq,lwq,LDO) /(KNIT + tracerpp(kwq,lwq,LDO) )
    ELSE
      f_nitrif = 1.0
    END IF
  nitrification = R_nitrif * f_nitrif * (Theta_nitrif**(salp(kwq,lwq) - 20.0)) * tracerpp(kwq,lwq,LNH4) * hpp(kwq,lwq)
  ! Units: [mg/m^2/s] = [1/s] * [-] * [-] * [mg/m^3] * [m]

  !. . .Add contribution from atmospheric deposition to top layer
  IF (kwq .eq. k1z(lwq)) THEN
    atmosdepositionNH4 = ATM_NH4
    ! Units: [mg/m^2/s] =  [mg/m^2/s] 
  ELSE
    atmosdepositionNH4 = 0.0
  END IF

  !. . . Add contribution from sediment flux to the bottom layer
  IF (kwq .eq. kmz(lwq)) THEN
      IF (IDO == 1) THEN
         ! Calculate DO inhibition of sediment flux
         ! f_sedflux = tracerpp(kwq,lwq,LDO) /(KSED + tracerpp(kwq,lwq,LDO) )
         f_sedflux = KSED / (KSED + tracerpp(kwq,lwq,LDO))
      ELSE
         f_sedflux = 1.0
      END IF
      sedfluxNH4 = SED_NH4 * f_sedflux* (Theta_sedflux**(salp(kwq,lwq) - 20))  
      ! Units: [mg/m^2/s] = [mg/m^2/s] * [-] * [-]
  ELSE
      sedfluxNH4 = 0.0

  END IF

  ! ... Incorporate all terms to the source-sink NH4 term
  sourcesink(kwq,lwq,LNH4) = sourcesink(kwq,lwq,LNH4)         &
                         &        -  nitrification            &
                         &        + atmosdepositionNH4        &
                         &        + sedfluxNH4
                        !         + mineralizationDON - only if IDON = 1; caluclated in sourcePON
                        !         - algal uptake    - if IALG = 1; calculated in sourceALG


  !. . Add contribution from nitrification to NO3 and DO sourcesink
  IF (INO3 == 1) THEN
    sourcesink(kwq,lwq,LNO3) = sourcesink(kwq,lwq,LNO3) + nitrification
  END IF

  IF (IDO == 1) THEN
    sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) - ron*nitrification
  END IF

END SUBROUTINE sourceNH4

!************************************************************************
SUBROUTINE sourceNO3(kwq,lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on NO3
!  Important: Sinksource term must have flux units [M/L^2/T], or
!             [M/L^3]*L*(1/T)
!--------------------------------------------------------------------------

  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq

  !. . . Local Variables
  REAL:: denitrification, atmosdepositionNO3, sedfluxNO3, f_sedflux

  ! ... Eliminate negative values of NO3 (min detection)
  IF (tracerpp(kwq,lwq,LNO3) .lt. 0) THEN
     tracerpp(kwq,lwq,LNO3) = 1.000
  END IF

  ! ... Denitrification (release of N2 gas)
  denitrification = R_denit* (Theta_denit**(salp(kwq,lwq) - 20.0)) * tracerpp(kwq,lwq,LNO3) * hpp(kwq,lwq)
  ! Units: [mg/m^2/s] = [1/s] * [-] * [-] * [mg/m^3] * [m]

  !. . .Add contribution from atmospheric deposition to top layer
  IF (kwq .eq. k1z(lwq)) THEN
    atmosdepositionNO3 = ATM_NO3
    ! Units: [mg/m^2/s] =  [mg/m^2/s] 
  ELSE
    atmosdepositionNO3 = 0.0
  END IF

  !. . . Add contribution from sediment flux to the bottom layer
  IF (kwq .eq. kmz(lwq)) THEN
      IF (IDO == 1) THEN
         ! Calculate DO inhibition of sediment flux
         f_sedflux = tracerpp(kwq,lwq,LDO) /(KSED + tracerpp(kwq,lwq,LDO) )
      ELSE
         f_sedflux = 1.0
      END IF
      sedfluxNO3 = SED_NO3 * f_sedflux* (Theta_sedflux**(salp(kwq,lwq) - 20)) 
      ! Units: [mg/m^2/s] = [mg/m^2/s] * [-] * [-]
  ELSE
      sedfluxNO3 = 0.0

  END IF

  ! ... Incorporate all terms to the source-sink NO3 term
  sourcesink(kwq,lwq,LNO3) = sourcesink(kwq,lwq,LNO3)         &
                         &        - denitrification           &
                         &        + atmosdepositionNO3        &
                         &        + sedfluxNO3
                        !         + nitrification   - if INH4 = 1; calculated in sourceNH4
                        !         - algal uptake    - if IALG = 1; calculated in sourceALG

END SUBROUTINE sourceNO3

!************************************************************************
SUBROUTINE sourcePOP(kwq,lwq)
!*********************************************************************

! Purpose: To calculate sourcesink terms that depend on POP
!  Important: Sinksource term must have flux units [M/L^2/T], or
!             [M/L^3]*L*(1/T)
!--------------------------------------------------------------------------

  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq

  !... Local variables
  REAL:: decompositionPOP, f_decom, depositionPOP, resuspensionPOP

  depositionPOP = 0.0
  resuspensionPOP = 0.0
  decompositionPOP = 0.0

  !... Calculate decompositionPOP
    ! Calculate DO inhibition of decomposition
    IF (IDO == 1) THEN
      f_decom = tracerpp(kwq,lwq,LDO) /(KDECMIN + tracerpp(kwq,lwq,LDO) )  ! We assume that decomposition and mineralzation 
                                                                            ! half-saturation values are very similar
    ELSE
      f_decom = 1.0
    END IF
  decompositionPOP = R_decom_pop * f_decom * (Theta_decom**(salp(kwq,lwq) - 20.0)) *tracerpp(kwq,lwq,LPOP) * hpp(kwq,lwq)
  ! Units: [mg/m^2/s] = [1/s] * [-] * [-] * [mg/m^3] * [m]

  ! ... Calculate deposition of POP only in the bottom layer
  IF (kwq .eq. kmz(lwq)) THEN
    depositionPOP = R_settl * tracerpp(kwq, lwq, LPOP)
    ! Units: [mg/m^2/s] =  [m/s] * [mg/m^3] 
    ! ... Calculate resusupension of POP only in the bottom layer
    resuspensionPOP = R_resusp * tracerpp(kwq,lwq,LPOP)
    ! Units: [mg/m^2/s] =  [m/s] * [mg/m^3]

    ! resuspensionPOP = R_resusp * tracerpp(kwq, lwq, LPOP)
    ! Units: [mg/m^2/s] =  [m/s] * [mg/m^3]
    if (depositionPOP .gt. (tracerpp(kwq, lwq, LPOP) * hp(kwq, lwq) / dt)) then
      depositionPOP = tracerpp(kwq, lwq, LPOP) * hp(kwq, lwq) / dt
    end if
  END IF

  ! ... Incorporate all terms to the source-sink POP term
  sourcesink(kwq,lwq,LPOP) = sourcesink(kwq,lwq,LPOP)       &
                        &    - decompositionPOP             &   
                        &    - depositionPOP                  &                   
                        &    + resuspensionPOP
                        !    + mortality  - only if IALG = 1; calcualted in sourceALG
                              
  ! Add contribution of decomposition to DOP concentration
  IF (iDOP == 1) THEN
    sourcesink(kwq, lwq, LDOP) = sourcesink(kwq, lwq, LDOP) +  decompositionPOP
  END IF

END SUBROUTINE sourcePOP

!************************************************************************
SUBROUTINE sourceDOP(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on DOP
!  Important: Sinksource term must have flux units [M/L^2/T], or
!             [M/L^3]*L*(1/T)
!--------------------------------------------------------------------------

  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq

  !. . . Local variables
  REAL:: mineralizationDOP, f_miner, atmosdepositionDOP, f_sedflux, sedfluxDOP

  atmosdepositionDOP = 0.0
  sedfluxDOP = 0.0

  !. . Mineralization by bacteria
    ! Calculate DO inhibition of mineralization
    IF (IDO == 1) THEN
      f_miner = tracerpp(kwq,lwq,LDO) /(KDECMIN + tracerpp(kwq,lwq,LDO) )
    ELSE
      f_miner = 1.0
    END IF
  mineralizationDOP = R_miner_dop * f_miner * (Theta_miner**(salp(kwq,lwq) - 20.0)) *tracerpp(kwq,lwq,LDOP) * hpp(kwq,lwq)
  ! Units: [mg/m^2/s] = [1/s] * [-] * [-] * [mg/m^3] * [m]

  !. . .Add contribution from atmospheric deposition to top layer
  IF (kwq .eq. k1z(lwq)) THEN
    atmosdepositionDOP = ATM_DOP
    ! Units: [mg/m^2/s] =  [mg/m^2/s] 
  END IF

  !. . .Add contribution from sediment flux to the bottom layer
  IF (kwq .eq. kmz(lwq)) THEN
    IF (IDO == 1) THEN
      ! Calculate DO inhibition of sediment flux
      ! f_sedflux = tracerpp(kwq,lwq,LDO) /(KSED + tracerpp(kwq,lwq,LDO) )
      f_sedflux = KSED /(KSED + tracerpp(kwq,lwq,LDO) )
    ELSE
      f_sedflux = 1.0
    END IF
    sedfluxDOP = SED_DOP * f_sedflux* (Theta_sedflux**(salp(kwq,lwq) - 20))  
    ! Units: [mg/m^2/s] =  [mg/m^2/s] * [-] * [-]
  END IF

  ! ... Incorporate all terms to the source-sink DON term
  sourcesink(kwq,lwq,LDOP) = sourcesink(kwq,lwq,LDOP)         &
                         &        -  mineralizationDOP        &
                         &        + atmosdepositionDOP        &
                         &        + sedfluxDOP
                        !         + decompositionPOP - only if IPOP = 1; caluclated in sourcePON

  ! Add mineralization of DOP to PO4
  IF (IPO4 == 1) THEN   
    sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) + mineralizationDOP
  END IF

END SUBROUTINE sourceDOP

!************************************************************************
SUBROUTINE sourcePO4(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on PO4
!  Important: Sinksource term must have flux units [M/L^2/T], or
!             [M/L^3]*L*(1/T)
!--------------------------------------------------------------------------

  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq

  !. . . Local Variables
  REAL:: atmosdepositionPO4, sedfluxPO4, f_sedflux

  atmosdepositionPO4 = 0.0

  !. . .Add contribution from atmospheric deposition to top layer
  IF (kwq .eq. k1z(lwq)) THEN
    atmosdepositionPO4 = ATM_PO4
    ! Units: [mg/m^2/s] =  [mg/m^2/s] 
  END IF

  !. . . Add contribution from sediment flux to the bottom layer
  IF (kwq .eq. kmz(lwq)) THEN
    IF (IDO == 1) THEN
      ! Calculate DO inhibition of sediment flux
      ! f_sedflux = tracerpp(kwq,lwq,LDO) /(KSED + tracerpp(kwq,lwq,LDO) )
      f_sedflux = KSED / (KSED + tracerpp(kwq, lwq, LDO))
    ELSE
      f_sedflux = 1.0
    END IF
    sedfluxPO4 = SED_PO4 * f_sedflux * (Theta_sedflux ** (salp(kwq, lwq) - 20))  
    ! Units: [mg/m^2/s] =  [mg/m^2/s] * [-] * [-]
  ELSE
    sedfluxPO4 = 0.0
  END IF

  ! ... Incorporate all terms to the source-sink PO4 term
  sourcesink(kwq,lwq,LPO4) = sourcesink(kwq, lwq, LPO4)         &
                          & + atmosdepositionPO4        &
                          & + sedfluxPO4
                          ! + mineralizationDOP - if IDOP = 1; calculated in sourceDOP
                          ! - algal uptake      - if IALG = 1; calculated in sourceALG

END SUBROUTINE sourcePO4

!************************************************************************
SUBROUTINE sourcePOC (kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on POC
!  Important: Sinksource term must have flux units [M/L^2/T], or
!             [M/L^3]*L*(1/T)
!--------------------------------------------------------------------------

  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq

  !... Local variables
  REAL:: decompositionPOC, f_decom, depositionPOC, resuspensionPOC
  real :: R_resusp_ij
  real :: decomp_poc_sed
  integer :: i, j

  depositionPOC = 0.0
  resuspensionPOC = 0.0
  decompositionPOC = 0.0
  decomp_poc_sed = 0.0

  !... Calculate decompositionPOC
  ! Calculate DO inhibition of decomposition
  IF (IDO == 1) THEN
    f_decom = tracerpp(kwq, lwq, LDO) / (KDECMIN + tracerpp(kwq, lwq, LDO) )   ! We assume that decomposition and mineralzation 
                                                                          ! half-saturation values are very similar
  ELSE
    f_decom = 1.0
  END IF
  
  decompositionPOC = R_decom_poc * f_decom * (Theta_decom ** (salp(kwq, lwq) - 20.0)) * tracerpp(kwq, lwq, LPOC) * hpp(kwq, lwq)
  ! Units: [mg/m^2/s] =  [1/s] * [-] * [-] * [mg/m^3] * [m]

  ! ... Calculate deposition of POC only in the bottom layer
  IF (kwq .eq. kmz(lwq)) THEN
    depositionPOC = vspoc * tracerpp(kwq,lwq,LPOC)
    ! Units: [mg/m^2/s] =  [m/s] * [mg/m^3] 
    ! ... Calculate resusupension of POC only in the bottom layer
    resuspensionPOC = R_resusp * tracerpp(kwq,lwq,LPOC)
    ! Units: [mg/m^2/s] =  [m/s] * [mg/m^3]
    if (depositionPOC .gt. (tracerpp(kwq, lwq, LPOC) * hp(kwq, lwq) / dt)) then
      depositionPOC = tracerpp(kwq, lwq, LPOC) * hp(kwq, lwq) / dt
    end if

    decomp_poc_sed = R_decom_poc *  (Theta_decom ** (salp(kwq, lwq) - 20.0)) * tracerpp(kwq + 1, lwq, LPOC) * hpp(kwq + 1, lwq)
    sourcesink(kwq + 1, lwq, LPOC) = sourcesink(kwq + 1, lwq, LPOC) - decomp_poc_sed + depositionPOC - resuspensionPOC
  END IF

  ! ... Incorporate all terms to the source-sink POC term
  sourcesink(kwq,lwq,LPOC) = sourcesink(kwq,lwq,LPOC) - decompositionPOC  - depositionPOC + resuspensionPOC
  !    + mortality  - only if IALG = 1; calcualted in sourceALG
      
  ! Add contribution of decomposition to DOC concentration
  IF (iDOC == 1) THEN
    sourcesink(kwq,lwq,LDOC) = sourcesink(kwq,lwq,LDOC) + decompositionPOC
    if (kwq .eq. kmz(lwq)) then
      sourcesink(kwq + 1, lwq, LDOC) = sourcesink(kwq + 1, lwq, LDOC) + decomp_poc_sed
    end if
  END IF

  fluxes_out(kwq, lwq, 4) = decompositionPOC
  fluxes_out(kwq, lwq, 5) = depositionPOC
  fluxes_out(kwq, lwq, 6) = resuspensionPOC
  fluxes_out(kwq + 1, lwq, 4) = decomp_poc_sed
  fluxes_out(kwq + 1, lwq, 5) = depositionPOC
  fluxes_out(kwq + 1, lwq, 6) = resuspensionPOC

END SUBROUTINE sourcePOC

!************************************************************************
SUBROUTINE sourceDOC(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on DOC
!  Important: Sinksource term must have flux units [M/L^2/T], or
!             [M/L^3]*L*(1/T)
!--------------------------------------------------------------------------

  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq

  !. . . Local variables
  REAL:: mineralizationDOC, f_miner, atmosdepositionDOC, f_sedflux, sedfluxDOC
  integer :: i, j
  real :: doc_miner_sed

  mineralizationDOC = 0.0
  f_miner = 0.0
  atmosdepositionDOC = 0.0
  f_sedflux = 0.0
  sedfluxDOC = 0.0
  doc_miner_sed = 0.0

  !. . Mineralization of DO by bacteria into inorganic nutrients (there is biological oxygen demand)
    ! Calculate DO inhibition of mineralization
    IF (IDO == 1) THEN
      f_miner = tracerpp(kwq,lwq,LDO) /(KDECMIN + tracerpp(kwq,lwq,LDO) )
    ELSE
      f_miner = 1.0
    END IF
  mineralizationDOC = R_miner_doc * f_miner * (Theta_miner**(salp(kwq,lwq) - 20.0)) * tracerpp(kwq,lwq,LDOC) * hpp(kwq,lwq)
  ! Units: [mg/m^2/s] =  [1/s] * [-] * [-] * [mg/m^3] * [m]

  !. . .Add contribution from atmospheric deposition to top layer
  IF (kwq .eq. k1z(lwq)) THEN
    atmosdepositionDOC = ATM_DOC
    ! Units: [mg/m^2/s] =  [mg/m^2/s]
  ELSE
    atmosdepositionDOC = 0.0
  END IF

  !. . .Add contribution from sediment flux to the bottom layer
  IF (kwq .eq. kmz(lwq)) THEN
    IF (IDO == 1) THEN
       ! Calculate DO inhibition of sediment flux
       f_sedflux = tracerpp(kwq,lwq,LDO) /(KSED + tracerpp(kwq,lwq,LDO) )
       ! f_sedflux = KSED /(KSED + tracerpp(kwq,lwq,LDO) )
    ELSE
       f_sedflux = 1.0
    END IF
    sedfluxDOC = SED_DOC * f_sedflux * (Theta_sedflux**(salp(kwq,lwq) - 20))
    ! Units: [mg/m^2/s] =  [mg/m^2/s] * [-] * [-]
    doc_miner_sed = R_miner_doc * (Theta_miner ** (salp(kwq, lwq) - 20.0)) * tracerpp(kwq + 1, lwq, LDOC) * hpp(kwq + 1, lwq)
    sourcesink(kwq + 1, lwq, LDOC) = sourcesink(kwq + 1, lwq, LDOC) - doc_miner_sed
  ELSE
      sedfluxDOC = 0.0
  END IF

  ! ... Incorporate all terms to the source-sink DON term
  sourcesink(kwq,lwq,LDOC) = sourcesink(kwq,lwq,LDOC) - mineralizationDOC + atmosdepositionDOC + sedfluxDOC
                        !         + decompositionPOC - only if IPOC = 1; caluclated in sourcePOC

  ! Remove oxygen due to microbial uptake (DOC mineralization)
  IF (IDO == 1) THEN   
      sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) - roc*mineralizationDOC
  END IF

  fluxes_out(kwq, lwq, 7) = atmosdepositionDOC
  fluxes_out(kwq, lwq, 8) = mineralizationDOC
  fluxes_out(kwq, lwq, 9) = sedfluxDOC
  fluxes_out(kwq + 1, lwq, 9) = sedfluxDOC

END SUBROUTINE sourceDOC

!************************************************************************
SUBROUTINE sourceALG1(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on ALG1
!  Important: Sinksource term must have flux units [M/L^2/T], or
!             [M/L^3]*L*(1/T)
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq

  integer :: i,j

  !. . .Local Variables
  REAL::  mu1, f_L1, f_T, f_N, f_P, N_conc
  REAL::  growth1, mort1, graz1, deposi1, resus1
  REAL::  Tmax1, Tmin1

  deposi1 = 0.0
  resus1 = 0.0

  !. .  Calculate growth limiting factors
  ! Light Limitation
  f_L1 = ((Qsw*QswFr(kwq,lwq))/light_sat1) *EXP(1 -((Qsw*QswFr(kwq,lwq))/light_sat1)) ! Steele (1962) = Photoinhibited
  IF (f_L1 .lt. 0.01) THEN
    f_L1 = 0.01
  END IF

  ! temperature limitaton
  Tmax1 = 35
  Tmin1 = 5
  IF (salp(kwq,lwq) .lt. Tmin1) THEN
    f_T = 0
  ELSE IF ((salp(kwq,lwq) .gt. Tmin1) .AND. (salp(kwq,lwq) .lt. Topt1)) THEN
    f_T = (salp(kwq,lwq) - Tmin1) / (Topt1 - Tmin1)
  ELSE IF ((salp(kwq,lwq) .gt. Topt1) .AND. (salp(kwq, lwq) .lt. Tmax1)) THEN 
    f_T = (Tmax1 - salp(kwq,lwq)) / (Tmax1 - Topt1)
  ELSE IF (salp(kwq, lwq) .gt. Tmax1) THEN
    f_T = 0
  END IF

  ! nutrient limitation - but only if the nutrients are modeled
  IF ((INH4 ==1) .AND. (INO3 ==1)) THEN
    N_conc = tracerpp(kwq,lwq,LNH4) + tracerpp(kwq,lwq,LNO3)
    f_N = N_conc/(KSN + N_conc)
  ELSE IF (INH4 == 1 .AND. INO3 ==0) THEN
    N_conc = tracerpp(kwq,lwq,LNH4)
    f_N = N_conc/(KSN + N_conc)
  ELSE IF (INH4 == 0 .AND. INO3 ==1) THEN
    N_conc = tracerpp(kwq,lwq,LNO3)
    f_N = N_conc/(KSN + N_conc)
  ELSE
    f_N = 1.0
  END IF

  IF (f_N .lt. 0) THEN ! Assume low nutrient concentration when f_N<0
    f_N =5.000/(KSN + 5.000)
  END IF

  IF (IPO4 ==1) THEN
    f_P = tracerpp(kwq,lwq,LPO4) /(KSP + tracerpp(kwq,lwq,LPO4) )
  ELSE
    f_P = 1.0
  END IF

  IF (f_P .lt. 0) THEN ! Assume low nutrient concentration when f_P<0
    f_P =5.000/(KSP + 5.000)
  END IF

  !. . Calculate growth
  mu1 = mu_max1 * MIN(f_L1,f_N,f_P)

  i = l2i(lwq)
  j = l2j(lwq)
  if (((i >= 1) .and. (i <= 139)) .and. ((j >=1) .and. (j <= 195))) then
    mu1 = mu1 * 1.0
  elseif ((i > 170) .and. ((j >= 1) .and. (j <= 65))) then
    mu1 = mu1 * 1.3
  elseif (((i > 139) .and. (i <= 170)) .and. ((j > 65) .and. (j <= 70))) then
    mu1 = mu1 * 1.0
  else
    mu1 = mu1
  end if

  growth1 = mu1 * f_T * tracerpp(kwq,lwq,LALG1) * hpp(kwq,lwq)
  ! Units: [mg/m^2/s] =  [1/s] * [-] * [mg/m^3] * [m]
  IF (tracerpp(kwq,lwq,LALG1) + (growth1 * dt / hpp(kwq, lwq)) .le. 10.00) THEN
    growth1 = 0.0 ! This limits growth to a minium phytoplankton concentration. If phyto < 0.01 ug/L, then growth will be zero
  END IF

  !. . Calculate mortality, respiration & excretion
  mort1   = R_mor1 * Theta_mor**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LALG1) * hpp(kwq,lwq)
  ! Units: [mg/m^2/s] =  [1/s] * [-] * [mg/m^3] * [m]
  IF ((tracerpp(kwq,lwq,LALG1) - mort1 ) .le. 10.00) THEN
    mort1 = 0.0 ! This limits grazing to a minium phytoplankton concentration. If phyto < 0.01 ug/L, then grazing will be zero
  END IF

  !. . Calculate grazing
  graz1   = R_gr1 * Theta_gr**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LALG1) * hpp(kwq,lwq)  
  ! Units: [mg/m^2/s] =  [1/s] * [-] * [mg/m^3] * [m]
  IF (tracerpp(kwq,lwq,LALG1) - (graz1 * dt / hpp(kwq, lwq)) .le. 10.00) THEN
    graz1 = 0.0 ! This limits grazing to a minium phytoplankton concentration. If phyto < 0.01 ug/L, then grazing will be zero
  END IF

  if (kwq .eq. kmz(lwq)) then
    !. . Calculate deposition
    deposi1 = vspa * tracerpp(kwq,lwq,LALG1)
    ! Units: [mg/m^2/s] =  [m/s] * [mg/m^3] 
    if (deposi1 .gt. (tracerpp(kwq, lwq, LALG1) * hp(kwq, lwq) / dt)) then
      deposi1 = tracerpp(kwq, lwq, LALG1) * hp(kwq, lwq) / dt
    end if 
    !. . Calculate resuspension
    resus1 = R_settl * tracerpp(kwq, lwq, LALG1)
    ! Units: [mg/m^2/s] =  [m/s] * [mg/m^3]
  end if

  ! Source-Sink equation
  sourcesink(kwq,lwq,LALG1) = sourcesink(kwq,lwq,LALG1) + growth1 - mort1 - graz1 - deposi1 + resus1

  ! .... Contributions from algae to other sub-routines (nutrients)
  ! If dissolved oxygen is modeled, alter sourcesink(kwq,lwq,LDO) to reflect photosynthetic oxygen production 
  ! and respiration consumption by phytoplankton
  IF (IDO ==1) THEN
    sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) + (roc*growth1 - roc*mort1)
  END IF

  ! IF PON is modeled, alter sourcesink(kwq,lwq,LPON) to reflect mortality of algae
  IF (IPON == 1) THEN
    sourcesink(kwq,lwq,LPON) = sourcesink(kwq,lwq,LPON) + rnc*mort1
  END IF

  ! IF NH4 is modeled, alter sourcesink(kwq,lwq,LNH4) to include uptake of NH4 by algae
  IF (INH4 == 1) THEN
    sourcesink(kwq,lwq,LNH4) = sourcesink(kwq,lwq,LNH4) - (rnc*FNH4*growth1)
  END IF

  ! IF NO3 is modeled, alter sourcesink(kwq,lwq,NO3) to include uptake of NO3 by algae
  IF (INO3 == 1) THEN
    sourcesink(kwq,lwq,LNO3) = sourcesink(kwq,lwq,LNO3) - (rnc*(1-FNH4) * growth1)
  END IF

  ! If POP is modeled, alter sourcesink(kwq,lwq,POP) to include mortality of algae
  IF (IPOP == 1) THEN
    sourcesink(kwq,lwq,LPOP) = sourcesink(kwq,lwq,LPOP) + (rpc*mort1)
  END IF

  ! If PO4 is modeled, alter sourcesink(kwq,lwq,LPO4) to include uptake
  IF (IPO4 == 1) THEN
    sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) - (rpc* growth1)
  END IF

  ! If POC is modeled, alter sourcesink(kwq,lwq,LPOC) to include mortality
  IF (IPOC == 1) THEN
    sourcesink(kwq,lwq,LPOC) = sourcesink(kwq,lwq,LPOC) + mort1
  END IF

  fluxes_out(kwq, lwq, 10) = growth1
  fluxes_out(kwq, lwq, 11) = mort1
  fluxes_out(kwq, lwq, 12) = graz1
  fluxes_out(kwq, lwq, 13) = deposi1
  fluxes_out(kwq + 1, lwq, 13) = deposi1
  fluxes_out(kwq, lwq, 14) = resus1
  fluxes_out(kwq + 1, lwq, 14) = resus1

END SUBROUTINE sourceALG1

!************************************************************************
SUBROUTINE sourceALG2(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on ALG1
!  Important: Sinksource term must have flux units [M/L^2/T], or
!             [M/L^3]*L*(1/T)
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq

  !. . .Local Variables
  REAL::  mu2, f_L2, f_T, f_N, f_P, N_conc
  REAL::  growth2, mort2, graz2, sett2, resus2
  REAL::  Tmax2, Tmin2

  !. .  Calculate growth limiting factors
    ! Light Limitation
    f_L2 = ((Qsw*QswFr(kwq,lwq))/light_sat2) *EXP(1 -((Qsw*QswFr(kwq,lwq))/light_sat2)) ! Steele (1962) = Photoinhibited
    IF (f_L2 .lt. 0.01) THEN
        f_L2 = 0.01
    END IF

    ! temperature limitaton
    Tmax2 = 30
    Tmin2 = 5
    IF (salp(kwq,lwq) .lt. Tmin2) THEN
      f_T = 0
    ELSE IF ((salp(kwq,lwq) .gt. Tmin2) .AND. (salp(kwq,lwq) .lt. Topt2)) THEN
      f_T = (salp(kwq,lwq) - Tmin2) / (Topt2 - Tmin2)
    ELSE IF (salp(kwq,lwq) .gt. Topt2) THEN 
      f_T = (Tmax2 - salp(kwq,lwq)) / (Tmax2 - Topt2)
    END IF

    ! nutrient limitation - but only if the nutrients are modeled
    IF ((INH4 ==1) .AND. (INO3 ==1)) THEN
      N_conc = tracerpp(kwq,lwq,LNH4) + tracerpp(kwq,lwq,LNO3)
      f_N = N_conc/(KSN + N_conc)
    ELSE IF (INH4 == 1 .AND. INO3 ==0) THEN
      N_conc = tracerpp(kwq,lwq,LNH4)
      f_N = N_conc/(KSN + N_conc)
    ELSE IF (INH4 == 0 .AND. INO3 ==1) THEN
      N_conc = tracerpp(kwq,lwq,LNO3)
      f_N = N_conc/(KSN + N_conc)
    ELSE
      f_N = 1.0
    END IF

    IF (f_N .lt. 0) THEN ! Assume low nutrient concentration when f_N<0
      f_N =5.0/(KSN + 5.0)
    END IF

    IF (IPO4 ==1) THEN
      f_P = tracerpp(kwq,lwq,LPO4) /(KSP + tracerpp(kwq,lwq,LPO4) )
    ELSE
      f_P = 1.0
    END IF

    IF (f_P .lt. 0) THEN ! Assume low nutrient concentration when f_P<0
       f_P =5.0/(KSP + 5.0)
    END IF

  !. . Calculate growth
  mu2 = mu_max2 * MIN(f_L2,f_N,f_P)
  growth2 = mu2 * f_T * tracerpp(kwq,lwq,LALG2) * hpp(kwq,lwq)
  ! Units: [mg/m^2/s] =  [1/s] * [-] * [mg/m^3] * [m]
  IF ((tracerpp(kwq,lwq,LALG2) + growth2) .le. 10.0) THEN
    growth2 = 0.0 ! This limits growth to a minium phytoplankton concentration. If phyto < 0.01 ug/L, then growth will be zero
  END IF

  !. . Calculate mortality, respiration & excretion
  mort2   = R_mor2 * Theta_mor**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LALG2) * hpp(kwq,lwq)
  ! Units: [mg/m^2/s] =  [1/s] * [-] * [mg/m^3] * [m]
  IF ((tracerpp(kwq,lwq,LALG2)- mort2 ) .le. 10.0) THEN
    mort2 = 0.0 ! This limits grazing to a minium phytoplankton concentration. If phyto < 0.01 ug/L, then grazing will be zero
  END IF

  !. . Calculate grazing
  graz2   = R_gr2 * Theta_gr**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LALG2) * hpp(kwq,lwq)  
  ! Units: [mg/m^2/s] =  [1/s] * [-] * [mg/m^3] * [m]
  IF ((tracerpp(kwq,lwq,LALG2)- graz2) .le. 10.0) THEN
    graz2 = 0.0 ! This limits grazing to a minium phytoplankton concentration. If phyto < 0.01 ug/L, then grazing will be zero
  END IF

  !. . Calculate deposition
  sett2   = R_settl * tracerpp(kwq,lwq,LALG2)
  ! Units: [mg/m^2/s] =  [m/s] * [mg/m^3]  
  !. . Calculate resuspension
  resus2   = R_resusp * tracerpp(kwq,lwq,LALG2)
  ! Units: [mg/m^2/s] =  [m/s] * [mg/m^3] 

  ! Source-Sink equation
  sourcesink(kwq,lwq,LALG2) = sourcesink(kwq,lwq,LALG2) + growth2 - mort2 - graz2 - sett2 + resus2

  ! .... Contributions from algae to other sub-routines (nutrients)
  ! If dissolved oxygen is modeled, alter sourcesink(kwq,lwq,LDO) to reflect photosynthetic oxygen production 
  ! and respiration consumption by phytoplankton
  IF (IDO ==1) THEN
    sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) + (roc*growth2 - roc*mort2)
  END IF

  ! IF PON is modeled, alter sourcesink(kwq,lwq,LPON) to reflect mortality of algae
  IF (IPON == 1) THEN
    sourcesink(kwq,lwq,LPON) = sourcesink(kwq,lwq,LPON) + rnc*mort2
  END IF

  ! IF NH4 is modeled, alter sourcesink(kwq,lwq,LNH4) to include uptake of NH4 by algae
  IF (INH4 == 1) THEN
    sourcesink(kwq,lwq,LNH4) = sourcesink(kwq,lwq,LNH4) - (rnc*FNH4*growth2)
  END IF

  ! IF NO3 is modeled, alter sourcesink(kwq,lwq,NO3) to include uptake of NO3 by algae
  IF (INO3 == 1) THEN
    sourcesink(kwq,lwq,LNO3) = sourcesink(kwq,lwq,LNO3) - (rnc*(1-FNH4) * growth2)
  END IF

  ! If POP is modeled, alter sourcesink(kwq,lwq,POP) to include mortality of algae
  IF (IPOP == 1) THEN
    sourcesink(kwq,lwq,LPOP) = sourcesink(kwq,lwq,LPOP) + (rpc*mort2)
  END IF

  ! If PO4 is modeled, alter sourcesink(kwq,lwq,LPO4) to include uptake
  IF (IPO4 == 1) THEN
    sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) - (rpc* growth2)
  END IF

  ! If POC is modeled, alter sourcesink(kwq,lwq,LPOC) to include mortality
  IF (IPOC == 1) THEN
    sourcesink(kwq,lwq,LPOC) = sourcesink(kwq,lwq,LPOC) + mort2
  END IF

END SUBROUTINE sourceALG2

!************************************************************************
SUBROUTINE sourceALG3(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on ALG3
!  Important: Sinksource term must have flux units [M/L^2/T], or
!             [M/L^3]*L*(1/T)
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq

  !. . .Local Variables
  REAL::  mu3, f_L3, f_T, f_N, f_P, N_conc
  REAL::  growth3, mort3, graz3, sett3, resus3
  REAL::  Tmax3, Tmin3

  !. .  Calculate growth limiting factors
    ! Light Limitation
    f_L3 = ((Qsw*QswFr(kwq,lwq))/light_sat3) *EXP(1 -((Qsw*QswFr(kwq,lwq))/light_sat3)) ! Steele (1962) = Photoinhibited
    IF (f_L3 .lt. 0.01) THEN
        f_L3 = 0.01
    END IF

    ! temperature limitaton
    Tmax3 = 30
    Tmin3 = 5
    IF (salp(kwq,lwq) .lt. Tmin3) THEN
      f_T = 0
    ELSE IF ((salp(kwq,lwq) .gt. Tmin3) .AND. (salp(kwq,lwq) .lt. Topt3)) THEN
      f_T = (salp(kwq,lwq) - Tmin3) / (Topt3 - Tmin3)
    ELSE IF (salp(kwq,lwq) .gt. Topt3) THEN 
      f_T = (Tmax3 - salp(kwq,lwq)) / (Tmax3 - Topt3)
    END IF

    ! nutrient limitation - but only if the nutrients are modeled
    IF ((INH4 ==1) .AND. (INO3 ==1)) THEN
      N_conc = tracerpp(kwq,lwq,LNH4) + tracerpp(kwq,lwq,LNO3)
      f_N = N_conc/(KSN + N_conc)
    ELSE IF (INH4 == 1 .AND. INO3 ==0) THEN
      N_conc = tracerpp(kwq,lwq,LNH4)
      f_N = N_conc/(KSN + N_conc)
    ELSE IF (INH4 == 0 .AND. INO3 ==1) THEN
      N_conc = tracerpp(kwq,lwq,LNO3)
      f_N = N_conc/(KSN + N_conc)
    ELSE
      f_N = 1.0
    END IF

    IF (f_N .lt. 0) THEN ! Assume low nutrient concentration when f_N<0
      f_N =5.000/(KSN + 5.000)
    END IF

    IF (IPO4 ==1) THEN
      f_P = tracerpp(kwq,lwq,LPO4) /(KSP + tracerpp(kwq,lwq,LPO4) )
    ELSE
      f_P = 1.0
    END IF

    IF (f_P .lt. 0) THEN ! Assume low nutrient concentration when f_P<0
       f_P =5.000/(KSP + 5.000)
    END IF

  !. . Calculate growth
  mu3 = mu_max3 * MIN(f_L3,f_N,f_P)
  growth3 = mu3 * f_T * tracerpp(kwq,lwq,LALG3) * hpp(kwq,lwq)
  ! Units: [mg/m^2/s] =  [1/s] * [-] * [mg/m^3] * [m]
  IF ((tracerpp(kwq,lwq,LALG3) + growth3) .le. 10.0) THEN
    growth3 = 0.0 ! This limits growth to a minium phytoplankton concentration. If phyto < 0.01 ug/L, then growth will be zero
  END IF

  !. . Calculate mortality, respiration & excretion
  mort3   = R_mor3 * Theta_mor**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LALG3) * hpp(kwq,lwq)
  ! Units: [mg/m^2/s] =  [1/s] * [-] * [mg/m^3] * [m]
  IF ((tracerpp(kwq,lwq,LALG3)- mort3 ) .le. 10.0) THEN
    mort3 = 0.0 ! This limits grazing to a minium phytoplankton concentration. If phyto < 0.01 ug/L, then grazing will be zero
  END IF

  !. . Calculate grazing
  graz3   = R_gr3 * Theta_gr**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LALG3) * hpp(kwq,lwq)  
  ! Units: [mg/m^2/s] =  [1/s] * [-] * [mg/m^3] * [m]
  IF ((tracerpp(kwq,lwq,LALG3)- graz3) .le. 10.0) THEN
    graz3 = 0.0 ! This limits grazing to a minium phytoplankton concentration. If phyto < 0.01 ug/L, then grazing will be zero
  END IF

  !. . Calculate deposition
  sett3   = R_settl * tracerpp(kwq,lwq,LALG3)
  ! Units: [mg/m^2/s] =  [m/s] * [mg/m^3]  
  !. . Calculate resuspension
  resus3   = R_resusp * tracerpp(kwq,lwq,LALG3)
  ! Units: [mg/m^2/s] =  [m/s] * [mg/m^3] 

  ! Source-Sink equation
  sourcesink(kwq,lwq,LALG3) = sourcesink(kwq,lwq,LALG3) + growth3 - mort3 - graz3 - sett3 + resus3

  ! .... Contributions from algae to other sub-routines (nutrients)
  ! If dissolved oxygen is modeled, alter sourcesink(kwq,lwq,LDO) to reflect photosynthetic oxygen production 
  ! and respiration consumption by phytoplankton
  IF (IDO ==1) THEN
    sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) + (roc*growth3 - roc*mort3)
  END IF

  ! IF PON is modeled, alter sourcesink(kwq,lwq,LPON) to reflect mortality of algae
  IF (IPON == 1) THEN
    sourcesink(kwq,lwq,LPON) = sourcesink(kwq,lwq,LPON) + rnc*mort3
  END IF

  ! IF NH4 is modeled, alter sourcesink(kwq,lwq,LNH4) to include uptake of NH4 by algae
  IF (INH4 == 1) THEN
    sourcesink(kwq,lwq,LNH4) = sourcesink(kwq,lwq,LNH4) - (rnc*FNH4*growth3)
  END IF

  ! IF NO3 is modeled, alter sourcesink(kwq,lwq,NO3) to include uptake of NO3 by algae
  IF (INO3 == 1) THEN
    sourcesink(kwq,lwq,LNO3) = sourcesink(kwq,lwq,LNO3) - (rnc*(1-FNH4) * growth3)
  END IF

  ! If POP is modeled, alter sourcesink(kwq,lwq,POP) to include mortality of algae
  IF (IPOP == 1) THEN
    sourcesink(kwq,lwq,LPOP) = sourcesink(kwq,lwq,LPOP) + (rpc*mort3)
  END IF

  ! If PO4 is modeled, alter sourcesink(kwq,lwq,LPO4) to include uptake
  IF (IPO4 == 1) THEN
    sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) - (rpc* growth3)
  END IF

  ! If POC is modeled, alter sourcesink(kwq,lwq,LPOC) to include mortality
  IF (IPOC == 1) THEN
    sourcesink(kwq,lwq,LPOC) = sourcesink(kwq,lwq,LPOC) + mort3
  END IF

END SUBROUTINE sourceALG3

!************************************************************************
SUBROUTINE sourceALG4(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on ALG4
!  Important: Sinksource term must have flux units [M/L^2/T], or
!             [M/L^3]*L*(1/T)
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq

  !. . .Local Variables
  REAL::  mu4, f_L4, f_T, f_N, f_P, N_conc
  REAL::  growth4, mort4, graz4, sett4, resus4
  REAL::  Tmax4, Tmin4

  !. .  Calculate growth limiting factors
    ! Light Limitation
    f_L4 = ((Qsw*QswFr(kwq,lwq))/light_sat4) *EXP(1 -((Qsw*QswFr(kwq,lwq))/light_sat4)) ! Steele (1962) = Photoinhibited
    IF (f_L4 .lt. 0.01) THEN
        f_L4 = 0.01
    END IF

    ! temperature limitaton
    Tmax4 = 30
    Tmin4 = 5
    IF (salp(kwq,lwq) .lt. Tmin4) THEN
      f_T = 0
    ELSE IF ((salp(kwq,lwq) .gt. Tmin4) .AND. (salp(kwq,lwq) .lt. Topt4)) THEN
      f_T = (salp(kwq,lwq) - Tmin4) / (Topt4 - Tmin4)
    ELSE IF (salp(kwq,lwq) .gt. Topt4) THEN 
      f_T = (Tmax4 - salp(kwq,lwq)) / (Tmax4 - Topt4)
    END IF

    ! nutrient limitation - but only if the nutrients are modeled
    IF ((INH4 ==1) .AND. (INO3 ==1)) THEN
      N_conc = tracerpp(kwq,lwq,LNH4) + tracerpp(kwq,lwq,LNO3)
      f_N = N_conc/(KSN + N_conc)
    ELSE IF (INH4 == 1 .AND. INO3 ==0) THEN
      N_conc = tracerpp(kwq,lwq,LNH4)
      f_N = N_conc/(KSN + N_conc)
    ELSE IF (INH4 == 0 .AND. INO3 ==1) THEN
      N_conc = tracerpp(kwq,lwq,LNO3)
      f_N = N_conc/(KSN + N_conc)
    ELSE
      f_N = 1.0
    END IF

    IF (f_N .lt. 0) THEN ! Assume low nutrient concentration when f_N<0
      f_N =5.000/(KSN + 5.000)
    END IF

    IF (IPO4 ==1) THEN
      f_P = tracerpp(kwq,lwq,LPO4) /(KSP + tracerpp(kwq,lwq,LPO4) )
    ELSE
      f_P = 1.0
    END IF

    IF (f_P .lt. 0) THEN ! Assume low nutrient concentration when f_P<0
       f_P =5.000/(KSP + 5.000)
    END IF

  !. . Calculate growth
  mu4 = mu_max4 * MIN(f_L4,f_N,f_P)
  growth4 = mu4 * f_T * tracerpp(kwq,lwq,LALG4) * hpp(kwq,lwq)
  ! Units: [mg/m^2/s] =  [1/s] * [-] * [mg/m^3] * [m]
  IF ((tracerpp(kwq,lwq,LALG4) + growth4) .le. 10.0) THEN
    growth4 = 0.0 ! This limits growth to a minium phytoplankton concentration. If phyto < 0.01 ug/L, then growth will be zero
  END IF

  !. . Calculate mortality, respiration & excretion
  mort4   = R_mor4 * Theta_mor**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LALG4) * hpp(kwq,lwq)
  ! Units: [mg/m^2/s] =  [1/s] * [-] * [mg/m^3] * [m]
  IF ((tracerpp(kwq,lwq,LALG4)- mort4 ) .le. 10.0) THEN
    mort4 = 0.0 ! This limits grazing to a minium phytoplankton concentration. If phyto < 0.01 ug/L, then grazing will be zero
  END IF

  !. . Calculate grazing
  graz4   = R_gr4 * Theta_gr**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LALG4) * hpp(kwq,lwq)  
  ! Units: [mg/m^2/s] =  [1/s] * [-] * [mg/m^3] * [m]
  IF ((tracerpp(kwq,lwq,LALG4)- graz4) .le. 10.0) THEN
    graz4 = 0.0 ! This limits grazing to a minium phytoplankton concentration. If phyto < 0.01 ug/L, then grazing will be zero
  END IF

  !. . Calculate deposition
  sett4   = R_settl * tracerpp(kwq,lwq,LALG4)
  ! Units: [mg/m^2/s] =  [m/s] * [mg/m^3]  
  !. . Calculate resuspension
  resus4   = R_resusp * tracerpp(kwq,lwq,LALG4)
  ! Units: [mg/m^2/s] =  [m/s] * [mg/m^3] 

  ! Source-Sink equation
  sourcesink(kwq,lwq,LALG4) = sourcesink(kwq,lwq,LALG4) + growth4 - mort4 - graz4 - sett4 + resus4

  ! .... Contributions from algae to other sub-routines (nutrients)
  ! If dissolved oxygen is modeled, alter sourcesink(kwq,lwq,LDO) to reflect photosynthetic oxygen production 
  ! and respiration consumption by phytoplankton
  IF (IDO ==1) THEN
    sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) + (roc*growth4 - roc*mort4)
  END IF

  ! IF PON is modeled, alter sourcesink(kwq,lwq,LPON) to reflect mortality of algae
  IF (IPON == 1) THEN
    sourcesink(kwq,lwq,LPON) = sourcesink(kwq,lwq,LPON) + rnc*mort4
  END IF

  ! IF NH4 is modeled, alter sourcesink(kwq,lwq,LNH4) to include uptake of NH4 by algae
  IF (INH4 == 1) THEN
    sourcesink(kwq,lwq,LNH4) = sourcesink(kwq,lwq,LNH4) - (rnc*FNH4*growth4)
  END IF

  ! IF NO3 is modeled, alter sourcesink(kwq,lwq,NO3) to include uptake of NO3 by algae
  IF (INO3 == 1) THEN
    sourcesink(kwq,lwq,LNO3) = sourcesink(kwq,lwq,LNO3) - (rnc*(1-FNH4) * growth4)
  END IF

  ! If POP is modeled, alter sourcesink(kwq,lwq,POP) to include mortality of algae
  IF (IPOP == 1) THEN
    sourcesink(kwq,lwq,LPOP) = sourcesink(kwq,lwq,LPOP) + (rpc*mort4)
  END IF

  ! If PO4 is modeled, alter sourcesink(kwq,lwq,LPO4) to include uptake
  IF (IPO4 == 1) THEN
    sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) - (rpc* growth4)
  END IF

  ! If POC is modeled, alter sourcesink(kwq,lwq,LPOC) to include mortality
  IF (IPOC == 1) THEN
    sourcesink(kwq,lwq,LPOC) = sourcesink(kwq,lwq,LPOC) + mort4
  END IF

END SUBROUTINE sourceALG4

!************************************************************************
SUBROUTINE sourceALG5(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on ALG5
!  Important: Sinksource term must have flux units [M/L^2/T], or
!             [M/L^3]*L*(1/T)
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq

  !. . .Local Variables
  REAL::  mu5, f_L5, f_T, f_N, f_P, N_conc
  REAL::  growth5, mort5, graz5, sett5, resus5
  REAL::  Tmax5, Tmin5

  !. .  Calculate growth limiting factors
    ! Light Limitation
    f_L5 = ((Qsw*QswFr(kwq,lwq))/light_sat5) *EXP(1 -((Qsw*QswFr(kwq,lwq))/light_sat5)) ! Steele (1962) = Photoinhibited
    IF (f_L5 .lt. 0.01) THEN
        f_L5 = 0.01
    END IF

    ! temperature limitaton
    Tmax5 = 30
    Tmin5 = 5
    IF (salp(kwq,lwq) .lt. Tmin5) THEN
      f_T = 0
    ELSE IF ((salp(kwq,lwq) .gt. Tmin5) .AND. (salp(kwq,lwq) .lt. Topt5)) THEN
      f_T = (salp(kwq,lwq) - Tmin5) / (Topt5 - Tmin5)
    ELSE IF (salp(kwq,lwq) .gt. Topt5) THEN 
      f_T = (Tmax5 - salp(kwq,lwq)) / (Tmax5 - Topt5)
    END IF

    ! nutrient limitation - but only if the nutrients are modeled
    IF ((INH4 ==1) .AND. (INO3 ==1)) THEN
      N_conc = tracerpp(kwq,lwq,LNH4) + tracerpp(kwq,lwq,LNO3)
      f_N = N_conc/(KSN + N_conc)
    ELSE IF (INH4 == 1 .AND. INO3 ==0) THEN
      N_conc = tracerpp(kwq,lwq,LNH4)
      f_N = N_conc/(KSN + N_conc)
    ELSE IF (INH4 == 0 .AND. INO3 ==1) THEN
      N_conc = tracerpp(kwq,lwq,LNO3)
      f_N = N_conc/(KSN + N_conc)
    ELSE
      f_N = 1.0
    END IF

    IF (f_N .lt. 0) THEN ! Assume low nutrient concentration when f_N<0
      f_N =5.000/(KSN + 5.000)
    END IF

    IF (IPO4 ==1) THEN
      f_P = tracerpp(kwq,lwq,LPO4) /(KSP + tracerpp(kwq,lwq,LPO4) )
    ELSE
      f_P = 1.0
    END IF

    IF (f_P .lt. 0) THEN ! Assume low nutrient concentration when f_P<0
       f_P =5.000/(KSP + 5.000)
    END IF

  !. . Calculate growth
  mu5 = mu_max5 * MIN(f_L5,f_N,f_P)
  growth5 = mu5 * f_T * tracerpp(kwq,lwq,LALG5) * hpp(kwq,lwq)
  ! Units: [mg/m^2/s] =  [1/s] * [-] * [mg/m^3] * [m]
  IF ((tracerpp(kwq,lwq,LALG5) + growth5) .le. 10.0) THEN
    growth5 = 0.0 ! This limits growth to a minium phytoplankton concentration. If phyto < 0.01 ug/L, then growth will be zero
  END IF

  !. . Calculate mortality, respiration & excretion
  mort5   = R_mor5 * Theta_mor**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LALG5) * hpp(kwq,lwq)
  ! Units: [mg/m^2/s] =  [1/s] * [-] * [mg/m^3] * [m]
  IF ((tracerpp(kwq,lwq,LALG5)- mort5 ) .le. 10.0) THEN
    mort5 = 0.0 ! This limits grazing to a minium phytoplankton concentration. If phyto < 0.01 ug/L, then grazing will be zero
  END IF

  !. . Calculate grazing
  graz5   = R_gr5 * Theta_gr**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LALG5) * hpp(kwq,lwq)  
  ! Units: [mg/m^2/s] =  [1/s] * [-] * [mg/m^3] * [m]
  IF ((tracerpp(kwq,lwq,LALG5)- graz5) .le. 10.0) THEN
    graz5 = 0.0 ! This limits grazing to a minium phytoplankton concentration. If phyto < 0.01 ug/L, then grazing will be zero
  END IF

  !. . Calculate deposition
  sett5   = R_settl * tracerpp(kwq,lwq,LALG5)
  ! Units: [mg/m^2/s] =  [m/s] * [mg/m^3]  
  !. . Calculate resuspension
  resus5   = R_resusp * tracerpp(kwq,lwq,LALG5)
  ! Units: [mg/m^2/s] =  [m/s] * [mg/m^3] 

  ! Source-Sink equation
  sourcesink(kwq,lwq,LALG5) = sourcesink(kwq,lwq,LALG5) + growth5 - mort5 - graz5 - sett5 + resus5

  ! .... Contributions from algae to other sub-routines (nutrients)
  ! If dissolved oxygen is modeled, alter sourcesink(kwq,lwq,LDO) to reflect photosynthetic oxygen production 
  ! and respiration consumption by phytoplankton
  IF (IDO ==1) THEN
    sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) + (roc*growth5 - roc*mort5)
  END IF

  ! IF PON is modeled, alter sourcesink(kwq,lwq,LPON) to reflect mortality of algae
  IF (IPON == 1) THEN
    sourcesink(kwq,lwq,LPON) = sourcesink(kwq,lwq,LPON) + rnc*mort5
  END IF

  ! IF NH4 is modeled, alter sourcesink(kwq,lwq,LNH4) to include uptake of NH4 by algae
  IF (INH4 == 1) THEN
    sourcesink(kwq,lwq,LNH4) = sourcesink(kwq,lwq,LNH4) - (rnc*FNH4*growth5)
  END IF

  ! IF NO3 is modeled, alter sourcesink(kwq,lwq,NO3) to include uptake of NO3 by algae
  IF (INO3 == 1) THEN
    sourcesink(kwq,lwq,LNO3) = sourcesink(kwq,lwq,LNO3) - (rnc*(1-FNH4) * growth5)
  END IF

  ! If POP is modeled, alter sourcesink(kwq,lwq,POP) to include mortality of algae
  IF (IPOP == 1) THEN
    sourcesink(kwq,lwq,LPOP) = sourcesink(kwq,lwq,LPOP) + (rpc*mort5)
  END IF

  ! If PO4 is modeled, alter sourcesink(kwq,lwq,LPO4) to include uptake
  IF (IPO4 == 1) THEN
    sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) - (rpc* growth5)
  END IF

  ! If POC is modeled, alter sourcesink(kwq,lwq,LPOC) to include mortality
  IF (IPOC == 1) THEN
    sourcesink(kwq,lwq,LPOC) = sourcesink(kwq,lwq,LPOC) + mort5
  END IF

END SUBROUTINE sourceALG5

!************************************************************************
                        END MODULE si3d_wq
!************************************************************************