!************************************************************************
                          MODULE si3d_wq
!************************************************************************
!
!  Purpose: Procedures that implement routines related to the modelling
!           of water quality in lakes
!
!-------------------------------------------------------------------------

 USE si3d_types
 IMPLICIT NONE
 SAVE

CONTAINS

!*********************************************************************
SUBROUTINE sourceDO(kwq,lwq, num_threads) 
!********************************************************************
!
!  Purpose: if dissolved oxygen is modeled, this subroutine
!  calculates source and sink terms that depend on dissolved
!  oxygen concentrations
!
!-----------------------------------------------------------------------

  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq
  REAL, INTENT(IN) :: num_threads

  !. . . Local Variables
  REAL    ::   Tk, lnOS, OS, Patm, ln_Pwv, Pwv, theta2, f_SOD 
  REAL    :: reaeration, sedoxydemand

  ! ...Calculate DO saturation
  Tk = salp(kwq,lwq) + 273
  lnOS = -139.34410 + 1.575701*1E5 /(Tk    ) &
  &                 - 6.642308*1E7 /(Tk**2.) &
  &                 + 1.243800*1E10/(Tk**3.) &
  &                 - 8.621949*1E11/(Tk**4.)
  OS = EXP(lnos)

  ! Correct for Patmospheric (Pa - declared in si3d_types and defined in surfbc0)

  Patm   = Pa * 0.00000986923; ! Transform atmospheric pressure from Pa to atm
  ln_Pwv = 11.8751 - (3840.70/Tk) - (216961/(Tk**2.))
  Pwv    = EXP(ln_Pwv)
  theta2 = 0.000975 - 1.426*1E-5 * salp(kwq,lwq) + &
                      6.436*1E-8 * salp(kwq,lwq)**2.
  OS = OS*Patm*((1-Pwv/Patm) *(1-theta2*Patm))&
  &           /((1-Pwv)*(1-theta2) )



 ! ...Calculate reaeration only at the lake surface
 ! for now using constant reaeration defined in wq_inp, but in future, can have
 ! alternatives for reaeration rates. 
 IF (kwq .eq. k1z(lwq)) THEN
    reaeration  = R_reaer*(OS - tracerpp(kwq,lwq,LDO))
    reaeration  = reaeration 
 ELSE
    reaeration  = 0.0
 END IF

 ! ...Calculate the sediment oxygen demand from the sediments (only bottom cell)
 IF (kwq .eq. kmz(lwq)) THEN
    f_SOD = tracerpp(kwq,lwq,LDO) /(KSOD + tracerpp(kwq,lwq,LDO) ) ! DO inhibition of sediment oxygen demand
    sedoxydemand = R_SOD * f_SOD* (Theta_SOD**(salp(kwq,lwq) - 20))  * (1/(hpp(kwq,lwq))) 
    sedoxydemand = sedoxydemand
ELSE
    sedoxydemand = 0.0

END IF

!IF (lwq .eq. 11) THEN
!      IF (kwq .eq. k1z(lwq)) THEN
!        PRINT *, 'DOs = ', OS,', DO = ',tracerpp(kwq,lwq,LDO),', reaeration = ',reaeration
!      END IF
!END IF

! ... Incorporate all terms to the source-sink DO term
 sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO)    &
            & +  reaeration                            &
            & -  sedoxydemand

! If ALG are modeled: Add photosynthetic oxygen production and respiration consumption  by phytoplankton
! If NH4 is modeled: Add consumption  by nitrification
! If DOC is modeled: Add consumption  by mineratlization of DOC

END SUBROUTINE sourceDO

!************************************************************************
SUBROUTINE sourcePON(kwq,lwq, num_threads)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on PON
!
!--------------------------------------------------------------------------

! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq
  REAL, INTENT(IN) :: num_threads

!... Local variables
REAL:: decompositionPON, f_decom, settlingPON, resuspensionPON

!... Calculate decompositionN
  ! Calculate DO inhibition of decomposition
  IF (IDO == 1) THEN
    f_decom = tracerpp(kwq,lwq,LDO) /(KDECMIN + tracerpp(kwq,lwq,LDO) )  ! We assume that decomposition and mineralzation 
                                                                          ! half-saturation values are very similar
  ELSE
    f_decom = 1.0
  END IF
decompositionPON = R_decom_pon * f_decom * (Theta_decom**(salp(kwq,lwq) - 20.0)) *tracerpp(kwq,lwq,LPON)
decompositionPON = decompositionPON

! ... Calculate settling of PON
settlingPON = R_settl * tracerpp(kwq,lwq,LPON) / hpp(kwq,lwq)
settlingPON = settlingPON

! ... Calculate resusupension of PON
resuspensionPON = R_resusp * tracerpp(kwq,lwq,LPON) / hpp(kwq,lwq)
resuspensionPON = resuspensionPON

! ... Incorporate all terms to the source-sink PON term
sourcesink(kwq,lwq,LPON) = sourcesink(kwq,lwq,LPON)       &
                      &    - decompositionPON             &   
                      &    - settlingPON                  &                   
                      &    + resuspensionPON
                      !    + mortality  - only if IALG = 1; calcualted in sourceALG
                            

  ! Add contribution of decomposition to DON concentration
  IF (iDON == 1) THEN
    sourcesink(kwq,lwq,LDON) = sourcesink(kwq,lwq,LDON) +  decompositionPON
  END IF

END SUBROUTINE sourcePON


!************************************************************************
SUBROUTINE sourceDON(kwq,lwq, num_threads)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on DON
!
!--------------------------------------------------------------------------

! ... Arguments
INTEGER, INTENT (IN) :: kwq,lwq
REAL, INTENT(IN) :: num_threads 

!. . . Local variables
REAL:: mineralizationDON, f_miner, atmosdepositionDON, f_sedflux, sedfluxDON

!. . Mineralization by bacteria
  ! Calculate DO inhibition of mineralization
  IF (IDO == 1) THEN
    f_miner = tracerpp(kwq,lwq,LDO) /(KDECMIN + tracerpp(kwq,lwq,LDO) )
  ELSE
    f_miner = 1.0
  END IF
mineralizationDON = R_miner_don * f_miner * (Theta_miner**(salp(kwq,lwq) - 20.0)) *tracerpp(kwq,lwq,LDON)
mineralizationDON = mineralizationDON

!. . .Add contribution from atmospheric deposition to top layer
IF (kwq .eq. k1z(lwq)) THEN
  atmosdepositionDON = ATM_DON/hpp(kwq,lwq)
  atmosdepositionDON = atmosdepositionDON
ELSE
  atmosdepositionDON = 0.0
END IF

!. . .Add contribution from sediment flux to the bottom layer
 IF (kwq .eq. kmz(lwq)) THEN
    IF (IDO == 1) THEN
       ! Calculate DO inhibition of sediment flux
       f_sedflux = tracerpp(kwq,lwq,LDO) /(KSED + tracerpp(kwq,lwq,LDO) )
         ELSE
       f_sedflux = 1.0
    END IF
    sedfluxDON = SED_DON * f_sedflux* (Theta_sedflux**(salp(kwq,lwq) - 20))  * (1/(hpp(kwq,lwq))) 
    sedfluxDON = sedfluxDON
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
SUBROUTINE sourceNH4(kwq,lwq, num_threads)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on NH4
!
!--------------------------------------------------------------------------

! ... Arguments
   INTEGER, INTENT (IN) :: kwq,lwq
   REAL, INTENT(IN) :: num_threads

!. . . Local Variables
REAL:: nitrification, f_nitrif, atmosdepositionNH4, sedfluxNH4, f_sedflux

!. . . Calculate nitrification
  ! Calculate DO inhibition of nitrification
  IF (IDO == 1) THEN
    f_nitrif = tracerpp(kwq,lwq,LDO) /(KNIT + tracerpp(kwq,lwq,LDO) )
  ELSE
    f_nitrif = 1.0
  END IF

nitrification = R_nitrif * f_nitrif * (Theta_nitrif**(salp(kwq,lwq) - 20.0)) *tracerpp(kwq,lwq,LNH4)
nitrification = nitrification


!. . .Add contribution from atmospheric deposition to top layer
IF (kwq .eq. k1z(lwq)) THEN
  atmosdepositionNH4 = ATM_NH4/hpp(kwq,lwq)
  atmosdepositionNH4 = atmosdepositionNH4
ELSE
  atmosdepositionNH4 = 0.0
END IF

!. . . Add contribution from sediment flux to the bottom layer
 IF (kwq .eq. kmz(lwq)) THEN
    IF (IDO == 1) THEN
       ! Calculate DO inhibition of sediment flux
       f_sedflux = tracerpp(kwq,lwq,LDO) /(KSED + tracerpp(kwq,lwq,LDO) )
         ELSE
       f_sedflux = 1.0
    END IF
    sedfluxNH4 = SED_NH4 * f_sedflux* (Theta_sedflux**(salp(kwq,lwq) - 20))  * (1/(hpp(kwq,lwq))) 
    sedfluxNH4 = sedfluxNH4
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
SUBROUTINE sourceNO3(kwq,lwq, num_threads)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on NO3
!
!--------------------------------------------------------------------------

! ... Arguments
INTEGER, INTENT (IN) :: kwq,lwq
REAL, INTENT(IN) :: num_threads

!. . . Local Variables
REAL:: denitrification, atmosdepositionNO3, sedfluxNO3, f_sedflux

! ... Denitrification (release of N2 gas)
denitrification = R_denit* (Theta_denit**(salp(kwq,lwq) - 20.0)) *tracerpp(kwq,lwq,LNO3)
denitrification = denitrification


!. . .Add contribution from atmospheric deposition to top layer
IF (kwq .eq. k1z(lwq)) THEN
  atmosdepositionNO3 = ATM_NO3/hpp(kwq,lwq)
  atmosdepositionNO3 = atmosdepositionNO3
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
    sedfluxNO3 = SED_NO3 * f_sedflux* (Theta_sedflux**(salp(kwq,lwq) - 20))  * (1/(hpp(kwq,lwq))) 
    sedfluxNO3 = sedfluxNO3
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
SUBROUTINE sourcePOP(kwq,lwq, num_threads)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on POP
!
!--------------------------------------------------------------------------

! ... Arguments
   INTEGER, INTENT (IN) :: kwq,lwq
   REAL, INTENT(IN) :: num_threads

!... Local variables
REAL:: decompositionPOP, f_decom, settlingPOP, resuspensionPOP

!... Calculate decompositionPOP
  ! Calculate DO inhibition of decomposition
  IF (IDO == 1) THEN
    f_decom = tracerpp(kwq,lwq,LDO) /(KDECMIN + tracerpp(kwq,lwq,LDO) )  ! We assume that decomposition and mineralzation 
                                                                          ! half-saturation values are very similar
  ELSE
    f_decom = 1.0
  END IF
decompositionPOP = R_decom_pop * f_decom * (Theta_decom**(salp(kwq,lwq) - 20.0)) *tracerpp(kwq,lwq,LPOP)
decompositionPOP = decompositionPOP

! ... Calculate settling of POP
settlingPOP = R_settl * tracerpp(kwq,lwq,LPOP) / hpp(kwq,lwq)
settlingPOP = settlingPOP

! ... Calculate resusupension of POP
resuspensionPOP = R_resusp * tracerpp(kwq,lwq,LPOP) / hpp(kwq,lwq)
resuspensionPOP = resuspensionPOP

! ... Incorporate all terms to the source-sink POP term
sourcesink(kwq,lwq,LPOP) = sourcesink(kwq,lwq,LPOP)       &
                      &    - decompositionPOP             &   
                      &    - settlingPOP                  &                   
                      &    + resuspensionPOP
                      !    + mortality  - only if IALG = 1; calcualted in sourceALG
                            

  ! Add contribution of decomposition to DOP concentration
  IF (iDOP == 1) THEN
    sourcesink(kwq,lwq,LDOP) = sourcesink(kwq,lwq,LDOP) +  decompositionPOP
  END IF

END SUBROUTINE sourcePOP

!************************************************************************
SUBROUTINE sourceDOP(kwq, lwq, num_threads)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on DOP
!
!--------------------------------------------------------------------------

! ... Arguments
   INTEGER, INTENT (IN) :: kwq,lwq
   REAL, INTENT(IN) :: num_threads

!. . . Local variables
REAL:: mineralizationDOP, f_miner, atmosdepositionDOP, f_sedflux, sedfluxDOP

!. . Mineralization by bacteria
  ! Calculate DO inhibition of mineralization
  IF (IDO == 1) THEN
    f_miner = tracerpp(kwq,lwq,LDO) /(KDECMIN + tracerpp(kwq,lwq,LDO) )
  ELSE
    f_miner = 1.0
  END IF
mineralizationDOP = R_miner_dop * f_miner * (Theta_miner**(salp(kwq,lwq) - 20.0)) *tracerpp(kwq,lwq,LDOP)
mineralizationDOP = mineralizationDOP

!. . .Add contribution from atmospheric deposition to top layer
IF (kwq .eq. k1z(lwq)) THEN
  atmosdepositionDOP = ATM_DOP/hpp(kwq,lwq)
  atmosdepositionDOP = atmosdepositionDOP
ELSE
  atmosdepositionDOP = 0.0
END IF

!. . .Add contribution from sediment flux to the bottom layer
 IF (kwq .eq. kmz(lwq)) THEN
    IF (IDO == 1) THEN
       ! Calculate DO inhibition of sediment flux
       f_sedflux = tracerpp(kwq,lwq,LDO) /(KSED + tracerpp(kwq,lwq,LDO) )
         ELSE
       f_sedflux = 1.0
    END IF
    sedfluxDOP = SED_DOP * f_sedflux* (Theta_sedflux**(salp(kwq,lwq) - 20))  * (1/(hpp(kwq,lwq))) 
    sedfluxDOP = sedfluxDOP
ELSE
    sedfluxDOP = 0.0

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
SUBROUTINE sourcePO4(kwq, lwq, num_threads)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on PO4
!
!--------------------------------------------------------------------------

! ... Arguments
INTEGER, INTENT (IN) :: kwq,lwq
REAL, INTENT(IN) :: num_threads

!. . . Local Variables
REAL:: atmosdepositionPO4, sedfluxPO4, f_sedflux

!. . .Add contribution from atmospheric deposition to top layer
IF (kwq .eq. k1z(lwq)) THEN
  atmosdepositionPO4 = ATM_PO4/hpp(kwq,lwq)
  atmosdepositionPO4 = atmosdepositionPO4
ELSE
  atmosdepositionPO4 = 0.0
END IF

!. . . Add contribution from sediment flux to the bottom layer
 IF (kwq .eq. kmz(lwq)) THEN
    IF (IDO == 1) THEN
       ! Calculate DO inhibition of sediment flux
       f_sedflux = tracerpp(kwq,lwq,LDO) /(KSED + tracerpp(kwq,lwq,LDO) )
         ELSE
       f_sedflux = 1.0
    END IF
    sedfluxPO4 = SED_PO4 * f_sedflux* (Theta_sedflux**(salp(kwq,lwq) - 20))  * (1/(hpp(kwq,lwq))) 
    sedfluxPO4 = sedfluxPO4
ELSE
    sedfluxPO4 = 0.0

END IF

! ... Incorporate all terms to the source-sink PO4 term
sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4)         &
                       &        + atmosdepositionPO4        &
                       &        + sedfluxPO4
                      !         + mineralizationDOP - if IDOP = 1; calculated in sourceDOP
                      !         - algal uptake      - if IALG = 1; calculated in sourceALG

END SUBROUTINE sourcePO4

!************************************************************************
SUBROUTINE sourcePOC (kwq, lwq, num_threads)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on POC
!
!--------------------------------------------------------------------------

  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq
  REAL, INTENT(IN) :: num_threads

!... Local variables
REAL:: decompositionPOC, f_decom, settlingPOC, resuspensionPOC

!... Calculate decompositionPOC
  ! Calculate DO inhibition of decomposition
  IF (IDO == 1) THEN
    f_decom = tracerpp(kwq,lwq,LDO) /(KDECMIN + tracerpp(kwq,lwq,LDO) )  ! We assume that decomposition and mineralzation 
                                                                          ! half-saturation values are very similar
  ELSE
    f_decom = 1.0
  END IF
decompositionPOC = R_decom_poc * f_decom * (Theta_decom**(salp(kwq,lwq) - 20.0)) *tracerpp(kwq,lwq,LPOC)
decompositionPOC = decompositionPOC

! ... Calculate settling of POC
settlingPOC = R_settl * tracerpp(kwq,lwq,LPOC) / hpp(kwq,lwq)
settlingPOC = settlingPOC

! ... Calculate resusupension of POC
resuspensionPOC = R_resusp * tracerpp(kwq,lwq,LPOC) / hpp(kwq,lwq)
resuspensionPOC = resuspensionPOC

! ... Incorporate all terms to the source-sink POC term
sourcesink(kwq,lwq,LPOC) = sourcesink(kwq,lwq,LPOC)       &
                      &    - decompositionPOC             &   
                      &    - settlingPOC                  &                   
                      &    + resuspensionPOC
                      !    + mortality  - only if IALG = 1; calcualted in sourceALG
                            

  ! Add contribution of decomposition to DOC concentration
  IF (iDOC == 1) THEN
    sourcesink(kwq,lwq,LDOC) = sourcesink(kwq,lwq,LDOC) +  decompositionPOC
  END IF

END SUBROUTINE sourcePOC


!************************************************************************
SUBROUTINE sourceDOC(kwq, lwq, num_threads)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on DOC
!
!--------------------------------------------------------------------------

  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq
  REAL, INTENT(IN) :: num_threads

!. . . Local variables
REAL:: mineralizationDOC, f_miner, atmosdepositionDOC, f_sedflux, sedfluxDOC

!. . Mineralization by bacteria
  ! Calculate DO inhibition of mineralization
  IF (IDO == 1) THEN
    f_miner = tracerpp(kwq,lwq,LDO) /(KDECMIN + tracerpp(kwq,lwq,LDO) )
  ELSE
    f_miner = 1.0
  END IF
mineralizationDOC = R_miner_doc * f_miner * (Theta_miner**(salp(kwq,lwq) - 20.0)) *tracerpp(kwq,lwq,LDOC)
mineralizationDOC = mineralizationDOC

!. . .Add contribution from atmospheric deposition to top layer
IF (kwq .eq. k1z(lwq)) THEN
  atmosdepositionDOC = ATM_DOC/hpp(kwq,lwq)
  atmosdepositionDOC = atmosdepositionDOC
ELSE
  atmosdepositionDOC = 0.0
END IF

!. . .Add contribution from sediment flux to the bottom layer
 IF (kwq .eq. kmz(lwq)) THEN
    IF (IDO == 1) THEN
       ! Calculate DO inhibition of sediment flux
       f_sedflux = tracerpp(kwq,lwq,LDO) /(KSED + tracerpp(kwq,lwq,LDO) )
         ELSE
       f_sedflux = 1.0
    END IF
    sedfluxDOC = SED_DOC * f_sedflux* (Theta_sedflux**(salp(kwq,lwq) - 20))  * (1/(hpp(kwq,lwq))) 
    sedfluxDOC = sedfluxDOC
ELSE
    sedfluxDOC = 0.0

END IF

! ... Incorporate all terms to the source-sink DON term
sourcesink(kwq,lwq,LDOC) = sourcesink(kwq,lwq,LDOC)         &
                       &        -  mineralizationDOC        &
                       &        + atmosdepositionDOC        &
                       &        + sedfluxDOC
                      !         + decompositionPOC - only if IPOC = 1; caluclated in sourcePOC


  ! Remove oxygen due to microbial uptake (DOC mineralization)
  IF (IDO == 1) THEN   
    sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) - roc*mineralizationDOC
  END IF

END SUBROUTINE sourceDOC
!************************************************************************
SUBROUTINE sourceALG1(kwq, lwq, num_threads)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on ALG1
!
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq
  REAL, INTENT(IN) :: num_threads

  !. . .Local Variables
  REAL::  mu1, f_L1, f_T, f_N, f_P, N_conc
  REAL::  growth1, mort1, graz1, sett1, resus1

  !. .  Calculate growth limiting factors
    ! Light Limitation
      f_L1 = ((Qsw*QswFr(kwq,lwq)*0.47)/light_sat1) *EXP(1 -((Qsw*QswFr(kwq,lwq)*0.47)/light_sat1)) ! Steele (1962) = Photoinhibited
      IF (f_L1 .lt. 0) THEN
         f_L1 = 0
      END IF

   ! temperature limitaton
    f_T = Theta_mu**(salp(kwq,lwq) - 20.0)

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

      IF (f_N .lt. 0) THEN
         f_N =0.005/(KSN + 0.005)
      END IF

  IF (IPO4 ==1) THEN
    f_P = tracerpp(kwq,lwq,LPO4) /(KSP + tracerpp(kwq,lwq,LPO4) )
  ELSE
    f_P = 1.0
  END IF

! ... Due to the paralelization, each process gets repeated as many times as the number of threads. The following lines get rid of that redundancy 
!. . Calculate growth
    mu1 = mu_max1 * MIN(f_L1,f_N,f_P)
    growth1 = mu1 * f_T * tracerpp(kwq,lwq,LALG1)
    !growth1 = growth1/REAL(num_threads)
!. . Calculate mortality, respiration & excretion
    mort1   = R_mor1 * Theta_mor**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LALG1)
    !mort1 = mort1 /REAL(num_threads)
      IF (mort1 .lt. 0.0) THEN
         mort1 = 0.0
      END IF

!. . Calculate grazing
    graz1   = R_gr1 * Theta_gr**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LALG1) 
    !graz1 = graz1/REAL(num_threads)
    IF ((tracerpp(kwq,lwq,LALG1)- graz1) .le. 0.01) THEN
        graz1 = 0.0 ! This limits grazing to a minium phytoplankton concentration. If phyto < 0.01 ug/L, then grazing will be zero
    END IF
!. . Calculate settling
    sett1   = R_settl * tracerpp(kwq,lwq,LALG1)/ hpp(kwq,lwq)
    !sett1 = sett1 /REAL(num_threads)
!. . Calculate resuspension
    resus1   = R_resusp * tracerpp(kwq,lwq,LALG1)/ hpp(kwq,lwq)
    !resus1 = resus1 /REAL(num_threads)

! Source-Sink equation
sourcesink(kwq,lwq,LALG1) = sourcesink(kwq,lwq,LALG1) + growth1 - mort1 - graz1 - sett1 + resus1

 !  IF (lwq .eq. 300) THEN
 !     IF (kwq .eq. 5) THEN
 !      PRINT *, 'lwq = ', lwq, ', kwq = ', kwq, ', f_L1 = ', f_L1, ', f_N = ', f_N, 'f_P = ', f_P,' , sourcesink(ALG1)  = ',sourcesink(kwq,lwq,LALG1) 
 !      PRINT *, 'growth1 = ', growth1, ', mort1 = ', mort1, ', graz1 = ', graz1,' , PhytoC = ',tracerpp(kwq,lwq,LALG1)
 !     END IF
 !  END IF

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

END SUBROUTINE sourceALG1


!************************************************************************
SUBROUTINE sourceALG2(kwq, lwq, num_threads)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on ALG2
!
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq
  REAL, INTENT(IN) :: num_threads

  !. . .Local Variables
  REAL::  mu2, f_L2, f_T, f_N, f_P, N_conc
  REAL::  growth2, mort2, graz2, sett2, resus2

  !. .  Calculate growth limiting factors
    ! Light Limitation
      f_L2 = ((Qsw*QswFr(kwq,lwq)*0.47)/light_sat2) *EXP(1 -((Qsw*QswFr(kwq,lwq)*0.47)/light_sat2)) ! Steele (1962) = Photoinhibited

   ! temperature limitaton
    f_T = Theta_mu**(salp(kwq,lwq) - 20.0)

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

  IF (IPO4 ==1) THEN
    f_P = tracerpp(kwq,lwq,LPO4) /(KSP + tracerpp(kwq,lwq,LPO4) )
  ELSE
    f_P = 1.0
  END IF

! ... Due to the paralelization, each process gets repeated as many times as the number of threads. The following lines get rid of that redundancy 
!. . Calculate growth
    mu2 = mu_max2 * MIN(f_L2,f_N,f_P)
    growth2 = mu2 * f_T * tracerpp(kwq,lwq,LALG2)
    growth2 = growth2
!. . Calculate mortality & excretion
    mort2   = R_mor2 * Theta_mor**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LALG2)
    mort2 = mort2 
!. . Calculate grazing
    graz2   = R_gr2 * Theta_gr**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LALG2) 
    graz2 = graz2
    IF ((tracerpp(kwq,lwq,LALG2)- graz2) .le. 0.01) THEN
        graz2 = 0.0 ! This limits grazing to a minium phytoplankton concentration. If phyto < 0.01 ug/L, then grazing will be zero
    END IF
!. . Calculate settling
    sett2   = R_settl * tracerpp(kwq,lwq,LALG2)
    sett2 = sett2
!. . Calculate resuspension
    resus2   = R_resusp * tracerpp(kwq,lwq,LALG2)
    resus2 = resus2

! Source-Sink equation
sourcesink(kwq,lwq,LALG2) = sourcesink(kwq,lwq,LALG2) + growth2 - mort2 - graz2 - sett2 + resus2


! If dissolved oxygen is modeled, alter sourcesink(kwq,lwq,LDO) to reflect growth and resp
! of algae population
IF (IDO ==1) THEN
  sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) + roc*growth2 - roc*mort2
END IF

! IF PON is modeled, alter sourcesink(kwq,lwq,LPON) to reflect mortality of algae
IF (IPON == 1) THEN
  sourcesink(kwq,lwq,LPON) = sourcesink(kwq,lwq,LPON) + rnc*mort2
END IF

! IF NH4 is modeled, alter sourcesink(kwq,lwq,LNH4) to include uptake of NH4 by algae
IF (INH4 == 1) THEN
  sourcesink(kwq,lwq,LNH4) = sourcesink(kwq,lwq,LNH4) - rnc*FNH4*growth2
END IF

! IF NO3 is modeled, alter sourcesink(kwq,lwq,NO3) to include uptake of NO3 by algae
IF (INO3 == 1) THEN
  sourcesink(kwq,lwq,LNO3) = sourcesink(kwq,lwq,LNO3) - rnc*(1-FNH4) * growth2
END IF

! If POP is modeled, alter sourcesink(kwq,lwq,POP) to include mortality of algae
IF (IPOP == 1) THEN
  sourcesink(kwq,lwq,LPOP) = sourcesink(kwq,lwq,LPOP) + rpc*mort2
END IF

! If PO4 is modeled, alter sourcesink(kwq,lwq,LPO4) to include uptake
IF (IPO4 == 1) THEN
  sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) - rpc* growth2
END IF

! If POC is modeled, alter sourcesink(kwq,lwq,LPOC) to include mortality
IF (IPOC == 1) THEN
  sourcesink(kwq,lwq,LPOC) = sourcesink(kwq,lwq,LPOC) + mort2
END IF

END SUBROUTINE sourceALG2


!************************************************************************
SUBROUTINE sourceALG3(kwq, lwq, num_threads)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on ALG3
!
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq
  REAL, INTENT(IN) :: num_threads

  !. . .Local Variables
  REAL::  mu3, f_L3, f_T, f_N, f_P, N_conc
  REAL::  growth3, mort3, graz3, sett3, resus3

  !. .  Calculate growth limiting factors
    ! Light Limitation
      f_L3 = ((Qsw*QswFr(kwq,lwq)*0.47)/light_sat3) *EXP(1 -((Qsw*QswFr(kwq,lwq)*0.47)/light_sat3)) ! Steele (1962) = Photoinhibited

   ! temperature limitaton
    f_T = Theta_mu**(salp(kwq,lwq) - 20.0)

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

  IF (IPO4 ==1) THEN
    f_P = tracerpp(kwq,lwq,LPO4) /(KSP + tracerpp(kwq,lwq,LPO4) )
  ELSE
    f_P = 1.0
  END IF

! ... Due to the paralelization, each process gets repeated as many times as the number of threads. The following lines get rid of that redundancy 
!. . Calculate growth
    mu3 = mu_max3 * MIN(f_L3,f_N,f_P)
    growth3 = mu3 * f_T * tracerpp(kwq,lwq,LALG3)
    growth3 = growth3
!. . Calculate mortality & excretion
    mort3   = R_mor3 * Theta_mor**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LALG3)
    mort3 = mort3 
!. . Calculate grazing
    graz3   = R_gr3 * Theta_gr**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LALG3) 
    graz3 = graz3
    IF ((tracerpp(kwq,lwq,LALG3)- graz3) .le. 0.01) THEN
        graz3 = 0.0 ! This limits grazing to a minium phytoplankton concentration. If phyto < 0.01 ug/L, then grazing will be zero
    END IF
!. . Calculate settling
    sett3   = R_settl * tracerpp(kwq,lwq,LALG3)
    sett3 = sett3 
!. . Calculate resuspension
    resus3   = R_resusp * tracerpp(kwq,lwq,LALG3)
    resus3 = resus3 

! Source-Sink equation
sourcesink(kwq,lwq,LALG3) = sourcesink(kwq,lwq,LALG3) + growth3 - mort3 - graz3 - sett3 + resus3


! If dissolved oxygen is modeled, alter sourcesink(kwq,lwq,LDO) to reflect growth and resp
! of algae population
IF (IDO ==1) THEN
  sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) + roc*growth3 - roc*mort3
END IF

! IF PON is modeled, alter sourcesink(kwq,lwq,LPON) to reflect mortality of algae
IF (IPON == 1) THEN
  sourcesink(kwq,lwq,LPON) = sourcesink(kwq,lwq,LPON) + rnc*mort3
END IF

! IF NH4 is modeled, alter sourcesink(kwq,lwq,LNH4) to include uptake of NH4 by algae
IF (INH4 == 1) THEN
  sourcesink(kwq,lwq,LNH4) = sourcesink(kwq,lwq,LNH4) - rnc*FNH4*growth3
END IF

! IF NO3 is modeled, alter sourcesink(kwq,lwq,NO3) to include uptake of NO3 by algae
IF (INO3 == 1) THEN
  sourcesink(kwq,lwq,LNO3) = sourcesink(kwq,lwq,LNO3) - rnc*(1-FNH4) * growth3
END IF

! If POP is modeled, alter sourcesink(kwq,lwq,POP) to include mortality of algae
IF (IPOP == 1) THEN
  sourcesink(kwq,lwq,LPOP) = sourcesink(kwq,lwq,LPOP) + rpc*mort3
END IF

! If PO4 is modeled, alter sourcesink(kwq,lwq,LPO4) to include uptake
IF (IPO4 == 1) THEN
  sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) - rpc* growth3
END IF

! If POC is modeled, alter sourcesink(kwq,lwq,LPOC) to include mortality
IF (IPOC == 1) THEN
  sourcesink(kwq,lwq,LPOC) = sourcesink(kwq,lwq,LPOC) + mort3
END IF

END SUBROUTINE sourceALG3

!************************************************************************
SUBROUTINE sourceALG4(kwq, lwq, num_threads)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on ALG4
!
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq
  REAL, INTENT(IN) :: num_threads

  !. . .Local Variables
  REAL::  mu4, f_L4, f_T, f_N, f_P, N_conc
  REAL::  growth4, mort4, graz4, sett4, resus4

  !. .  Calculate growth limiting factors
    ! Light Limitation
      f_L4 = ((Qsw*QswFr(kwq,lwq)*0.47)/light_sat4) *EXP(1 -((Qsw*QswFr(kwq,lwq)*0.47)/light_sat4)) ! Steele (1962) = Photoinhibited

   ! temperature limitaton
    f_T = Theta_mu**(salp(kwq,lwq) - 20.0)

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

  IF (IPO4 ==1) THEN
    f_P = tracerpp(kwq,lwq,LPO4) /(KSP + tracerpp(kwq,lwq,LPO4) )
  ELSE
    f_P = 1.0
  END IF

! ... Due to the paralelization, each process gets repeated as many times as the number of threads. The following lines get rid of that redundancy 
!. . Calculate growth
    mu4 = mu_max4 * MIN(f_L4,f_N,f_P)
    growth4 = mu4 * f_T * tracerpp(kwq,lwq,LALG4)
    growth4 = growth4
!. . Calculate mortality & excretion
    mort4   = R_mor4 * Theta_mor**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LALG4)
    mort4 = mort4 
!. . Calculate grazing
    graz4   = R_gr4 * Theta_gr**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LALG4) 
    graz4 = graz4
    IF ((tracerpp(kwq,lwq,LALG4)- graz4) .le. 0.01) THEN
        graz4 = 0.0 ! This limits grazing to a minium phytoplankton concentration. If phyto < 0.01 ug/L, then grazing will be zero
    END IF
!. . Calculate settling
    sett4   = R_settl * tracerpp(kwq,lwq,LALG4)
    sett4 = sett4 
!. . Calculate resuspension
    resus4   = R_resusp * tracerpp(kwq,lwq,LALG4)
    resus4 = resus4 

! Source-Sink equation
sourcesink(kwq,lwq,LALG4) = sourcesink(kwq,lwq,LALG4) + growth4 - mort4 - graz4 - sett4 + resus4


! If dissolved oxygen is modeled, alter sourcesink(kwq,lwq,LDO) to reflect growth and resp
! of algae population
IF (IDO ==1) THEN
  sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) + roc*growth4 - roc*mort4
END IF

! IF PON is modeled, alter sourcesink(kwq,lwq,LPON) to reflect mortality of algae
IF (IPON == 1) THEN
  sourcesink(kwq,lwq,LPON) = sourcesink(kwq,lwq,LPON) + rnc*mort4
END IF

! IF NH4 is modeled, alter sourcesink(kwq,lwq,LNH4) to include uptake of NH4 by algae
IF (INH4 == 1) THEN
  sourcesink(kwq,lwq,LNH4) = sourcesink(kwq,lwq,LNH4) - rnc*FNH4*growth4
END IF

! IF NO3 is modeled, alter sourcesink(kwq,lwq,NO3) to include uptake of NO3 by algae
IF (INO3 == 1) THEN
  sourcesink(kwq,lwq,LNO3) = sourcesink(kwq,lwq,LNO3) - rnc*(1-FNH4) * growth4
END IF

! If POP is modeled, alter sourcesink(kwq,lwq,POP) to include mortality of algae
IF (IPOP == 1) THEN
  sourcesink(kwq,lwq,LPOP) = sourcesink(kwq,lwq,LPOP) + rpc*mort4
END IF

! If PO4 is modeled, alter sourcesink(kwq,lwq,LPO4) to include uptake
IF (IPO4 == 1) THEN
  sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) - rpc* growth4
END IF

! If POC is modeled, alter sourcesink(kwq,lwq,LPOC) to include mortality
IF (IPOC == 1) THEN
  sourcesink(kwq,lwq,LPOC) = sourcesink(kwq,lwq,LPOC) + mort4
END IF

END SUBROUTINE sourceALG4


!************************************************************************
!SUBROUTINE sourceFT(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate FT
!
!--------------------------------------------------------------------------
!  ! ... Arguments
!  INTEGER, INTENT (IN) :: kwq,lwq  
! ! temperature limitaton
!   
!   tracerpp(kwq,lwq,LFT) = Theta_mu**(salp(kwq,lwq) - 20.0)
!
!END SUBROUTINE sourceFT


!************************************************************************
!SUBROUTINE sourceFN(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate FN
!
!--------------------------------------------------------------------------
!  ! ... Arguments
!  INTEGER, INTENT (IN) :: kwq,lwq  
!  REAL:: N_conc
!     
! ! nutrient limitation - but only if the nutrients are modeled
! IF ((INH4 ==1) .AND. (INO3 ==1)) THEN
!   N_conc = tracerpp(kwq,lwq,LNH4) + tracerpp(kwq,lwq,LNO3)
!   tracerpp(kwq,lwq,LFN)  = N_conc/(KSN + N_conc)
!
! ELSE IF (INH4 == 1 .AND. INO3 ==0) THEN
!   N_conc = tracerpp(kwq,lwq,LNH4)
!   tracerpp(kwq,lwq,LFN) = N_conc/(KSN + N_conc)
!
! ELSE IF (INH4 == 0 .AND. INO3 ==1) THEN
!   N_conc = tracerpp(kwq,lwq,LNO3)
!   tracerpp(kwq,lwq,LFN) = N_conc/(KSN + N_conc)
!
! ELSE
!   tracerpp(kwq,lwq,LFN) = 1.0
!   
! END IF
!
!END SUBROUTINE sourceFN

!************************************************************************
!SUBROUTINE sourceFP(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate FP
!
!--------------------------------------------------------------------------
!  ! ... Arguments
!  INTEGER, INTENT (IN) :: kwq,lwq    
!   
! ! nutrient limitation - but only if the nutrients are modeled
! IF (IPO4 ==1) THEN
!   tracerpp(kwq,lwq,LFP) = tracerpp(kwq,lwq,LPO4) /(KSP + tracerpp(kwq,lwq,LPO4) )
!   
! ELSE
!   tracerpp(kwq,lwq,LFP) = 1.0
!   
! END IF
!
!END SUBROUTINE sourceFP

!************************************************************************
!SUBROUTINE sourceFL(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate FL
!
!--------------------------------------------------------------------------
!  ! ... Arguments
!  INTEGER, INTENT (IN) :: kwq,lwq  
!
! tracerpp(kwq,lwq,LFL) = ((Qsw*QswFr(kwq,lwq)*0.47)/light_sat) *EXP(1 -((Qsw*QswFr(kwq,lwq)*0.47)/light_sat)) ! Steele (1962) = Photoinhibited 
! 
!END SUBROUTINE sourceFL

!************************************************************************
                        END MODULE si3d_wq
!************************************************************************