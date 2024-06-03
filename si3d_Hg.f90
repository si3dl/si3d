!*************************************************************************
                          MODULE si3d_Hg
!*************************************************************************
!
!  Purpose: Procedures that implement routines related to the modelling
!           of ecological processees for mercury
!
!-------------------------------------------------------------------------

  USE si3d_types
  USE si3d_sed

  IMPLICIT NONE
  SAVE

CONTAINS

!*************************************************************************
SUBROUTINE sourceHg(kwq, lwq)
!*************************************************************************
!
!   Purpose:
!
!-------------------------------------------------------------------------
  ! ... Arguments
  integer, intent(in) :: kwq
  integer, intent(in) :: lwq
  integer             :: kms
  integer             :: k1s

  real                       :: fwd2
  real                       :: fwdoc2
  real                       :: fwpa2
  real                       :: fwpom2
  real                       :: fsd2
  real                       :: fsdoc2
  real                       :: fspom2
  real                       :: HgII_wd
  real                       :: HgII_wdoc
  real                       :: HgII_wpa
  real                       :: HgII_wpom
  real                       :: HgII_sd
  real                       :: HgII_sdoc
  real                       :: HgII_spom
  real, dimension(sedNumber) :: fwpn2
  real, dimension(sedNumber) :: fspn2
  real, dimension(sedNumber) :: HgII_wpn
  real, dimension(sedNumber) :: HgII_spn

  real                       :: fwd3
  real                       :: fwdoc3
  real                       :: fwpa3
  real                       :: fwpom3 
  real                       :: fsd3
  real                       :: fsdoc3
  real                       :: fspom3
  real                       :: MeHg_wd
  real                       :: MeHg_wdoc
  real                       :: MeHg_wpa
  real                       :: MeHg_wpom
  real                       :: MeHg_sd
  real                       :: MeHg_sdoc
  real                       :: MeHg_spom
  real, dimension(sedNumber) :: fwpn3
  real, dimension(sedNumber) :: fspn3
  real, dimension(sedNumber) :: MeHg_wpn
  real, dimension(sedNumber) :: MeHg_spn

  real :: HgII_wddoc
  real :: MeHg_wddoc
  real :: HgII_sddoc
  real :: MeHg_sddoc
  real :: Hg0w
  real :: HgIIw
  real :: MeHgw
  real :: Hg0s
  real :: HgIIs
  real :: MeHgs

  real :: MeHgw_diffusion
  real :: HgIIw_diffusion
  real :: Hg0w_diffusion
  real :: HgIIs_methy
  real :: MeHgs_demethy
  real :: HgIIw_atmdep
  real :: MeHgw_atmdep
  real :: Hg0w_vol
  real :: MeHgw_vol
  real :: HgIIw_reduction
  real :: Hg0w_oxidation
  real :: MeHgw_photodeg
  real :: HgIIw_methy
  real :: MeHgw_demethy
  real :: HgIIw_deposition
  real :: MeHgw_deposition
  real :: HgIIs_burial
  real :: MeHgs_burial
  real :: HgIIs_erosion
  real :: MeHgs_erosion

  MeHgw_diffusion  = 0.0
  HgIIw_diffusion  = 0.0
  Hg0w_diffusion   = 0.0
  HgIIs_methy      = 0.0
  MeHgs_demethy    = 0.0
  HgIIw_atmdep     = 0.0
  MeHgw_atmdep     = 0.0
  Hg0w_vol         = 0.0
  MeHgw_vol        = 0.0
  HgIIw_reduction  = 0.0
  Hg0w_oxidation   = 0.0
  MeHgw_photodeg   = 0.0
  HgIIw_methy      = 0.0
  MeHgw_demethy    = 0.0
  HgIIw_deposition = 0.0
  HgIIs_erosion    = 0.0
  MeHgw_deposition = 0.0
  MeHgs_erosion    = 0.0
  HgIIs_burial = 0.0
  MeHgs_burial = 0.0

  HgII_wddoc = 0.0
  MeHg_wddoc = 0.0
  HgII_sddoc = 0.0
  MeHg_sddoc = 0.0
  Hg0w = 0.0
  HgIIw = 0.0
  MeHgw = 0.0
  Hg0s = 0.0
  HgIIs = 0.0
  MeHgs = 0.0

  fwd2 = 0.0
  fwdoc2 = 0.0
  fwpa2 = 0.0
  fwpom2 = 0.0
  fsd2 = 0.0
  fsdoc2 = 0.0
  fspom2 = 0.0
  HgII_wd = 0.0
  HgII_wdoc = 0.0
  HgII_wpa = 0.0
  HgII_wpom = 0.0
  HgII_sd = 0.0
  HgII_sdoc = 0.0
  HgII_spom = 0.0
  fwpn2(:) = 0.0
  fspn2(:) = 0.0
  HgII_wpn(:) = 0.0
  HgII_spn(:) = 0.0

  fwd3 = 0.0
  fwdoc3 = 0.0
  fwpa3 = 0.0
  fwpom3 = 0.0
  fsd3 = 0.0
  fsdoc3 = 0.0
  fspom3 = 0.0
  MeHg_wd = 0.0
  MeHg_wdoc = 0.0
  MeHg_wpa = 0.0
  MeHg_wpom = 0.0
  MeHg_sd = 0.0
  MeHg_sdoc = 0.0
  MeHg_spom = 0.0
  fwpn3(:) = 0.0
  fspn3(:) = 0.0
  MeHg_wpn(:) = 0.0
  MeHg_spn(:) = 0.0

  kms = kmz(lwq)
  k1s = k1z(lwq)

  ! Assign total mercury concentrations
  if (iHg0 .eq. 1) then
    Hg0w = tracerpp(kwq, lwq, LHg0)
    if (kwq .eq. kms) then
      Hg0s = tracerpp(kwq + 1, lwq, LHg0)
    end if
  end if

  if (iHgII .eq. 1) then
    HgIIw = tracerpp(kwq, lwq, LHgII)
    if (kwq .eq. kms) then
      HgIIs = tracerpp(kwq + 1, lwq, LHgII)
    end if
    ! Define partitioning (Instantaneous) 
    call HgII_partitioning(kms, kwq, lwq, HgIIw, HgIIs, fwd2, fwdoc2, fwpa2, fwpom2, fwpn2, &
       & fsd2, fsdoc2, fspom2, fspn2, HgII_wpn, HgII_wd, HgII_wdoc, HgII_wpa, &
       & HgII_wpom, HgII_spn, HgII_sd, HgII_sdoc, HgII_spom)
  end if

  if (iMeHg .eq. 1) then
    MeHgw = tracerpp(kwq, lwq, LMeHg)
    if (kwq .eq. kms) then
      MeHgs = tracerpp(kwq + 1, lwq, LMeHg)
    end if
    ! Define partitioning (Instantaneous)
    call MeHg_partitioning(kms, kwq, lwq, MeHgw, MeHgs, fwd3, fwdoc3, fwpa3, fwpom3, fwpn3, &
          & fsd3, fsdoc3, fspom3, fspn3, MeHg_wpn, MeHg_wd, MeHg_wdoc, &
          & MeHg_wpa, MeHg_wpom, MeHg_spn, MeHg_sd, MeHg_sdoc, MeHg_spom)
  end if

  !--------------------------------------------------------------------------------
  ! Source / Sink for the Mercury processes in the water column 
  
  ! Total dissolved HgII in water. Sum of dissolved and dissolved bounded to DOC
  HgII_wddoc = HgII_wd + HgII_wdoc
  MeHg_wddoc = MeHg_wd + MeHg_wdoc

  MeHg_sddoc = MeHg_sd + MeHg_sdoc
  HgII_sddoc = HgII_sd + HgII_sdoc

  ! HgII_wddoc = HgIIw
  ! MeHg_wddoc = MeHgw
  ! MeHg_sddoc = MeHgs
  ! HgII_sddoc = HgIIs

  call HgII_reduction(HgIIw_reduction, kwq, lwq, HgII_wddoc)
  call HgIIw_methylation(HgIIw_methy, kwq, lwq, HgII_wddoc)
  call MeHgw_demethylation(MeHgw_demethy, kwq, lwq, MeHg_wddoc)
  call MeHg_photodegradation(MeHgw_photodeg, kwq, lwq, MeHg_wddoc)
  call Hg0_oxidation(Hg0w_oxidation, Hg0w, HgII_wddoc, HgIIw_reduction, lwq)

  if (kwq .eq. k1s) then
    call HgII_atm_deposition(HgIIw_atmdep, kwq, lwq)
    call MeHg_atm_deposition(MeHgw_atmdep, kwq, lwq)
    call MeHg_volatilization(MeHgw_vol, MeHg_wddoc)
    call Hg0_volatilization(Hg0w_vol, Hg0w)
  end if

  if (kwq .eq. kms) then
    call Hg0_diffusion(Hg0w_diffusion, Hg0w, Hg0s)
    call MeHg_diffusion(MeHgw_diffusion, MeHg_wddoc, MeHg_sddoc)
    call HgII_diffusion(HgIIw_diffusion, HgII_wddoc, HgII_sddoc)
    call HgII_deposition(HgIIw_deposition, HgII_wpa, HgII_wpom, HgII_wpn)
    call MeHg_deposition(MeHgw_deposition, MeHg_wpa, MeHg_wpom, MeHg_wpn)
    call HgII_erosion(HgIIs_erosion, HgII_spn)
    call MeHg_erosion(MeHgs_erosion, MeHg_spn)
    call HgIIs_methylation(HgIIs_methy, kwq + 1, lwq, HgII_sddoc)
    call MeHgs_demethylation(MeHgs_demethy, kwq + 1, lwq, MeHg_sddoc)
    call HgII_burial(HgIIs_burial, HgIIw_deposition, HgIIs_erosion)
    call MeHg_burial(MeHgs_burial, MeHgw_deposition, MeHgs_erosion)
  end if

  sourcesink(kwq, lwq, LHg0)  = HgIIw_reduction + MeHgw_photodeg + Hg0w_diffusion - Hg0w_oxidation - Hg0w_vol
  sourcesink(kwq, lwq, LHgII) = HgIIw_atmdep + HgIIw_diffusion + Hg0w_oxidation + MeHgw_demethy - HgIIw_reduction - HgIIw_methy - HgIIw_deposition + HgIIs_erosion
  sourcesink(kwq, lwq, LMeHg) = MeHgw_atmdep + MeHgw_diffusion + HgIIw_methy - MeHgw_vol - MeHgw_photodeg - MeHgw_demethy - MeHgw_deposition + MeHgs_erosion


  ! Source / Sink for Mercury processes in the sediment layer
  if (kwq .eq. kms) then
    sourcesink(kwq + 1, lwq, LHg0)  = - Hg0w_diffusion
    sourcesink(kwq + 1, lwq, LHgII) = - HgIIw_diffusion + MeHgs_demethy - HgIIs_methy + HgIIw_deposition - HgIIs_erosion  - HgIIs_burial
    sourcesink(kwq + 1, lwq, LMeHg) = - MeHgw_diffusion - MeHgs_demethy + HgIIs_methy + MeHgw_deposition - MeHgs_erosion  - MeHgs_burial
  end if

  if (lwq == 106) then
    print*,'k = ',kwq
    print*,'h = ',h(kwq + 1,lwq)
    print*,'MeHgw = ', MeHgw
    print*,'MeHg_wd = ', MeHg_wd, 'MeHg_wdoc = ', MeHg_wdoc, 'MeHg_wpa = ', MeHg_wpa, 'MeHg_wpom = ', MeHg_wpom
    print*,'MeHg_wpn = ', MeHg_wpn
    print*,'MeHgw_veri = ', MeHg_wd + MeHg_wdoc + MeHg_wpa + MeHg_wpom + sum(MeHg_wpn)
    print*,'MeHgs = ', MeHgs
    print*,'MeHg_sd = ', MeHg_sd, 'MeHg_sdoc = ', MeHg_sdoc, 'MeHg_spom = ', MeHg_spom
    print*,'MeHg_spn = ', MeHg_spn
    print*,'MeHgs_veri = ', MeHg_sd + MeHg_sdoc + MeHg_spom + sum(MeHg_spn)
    print*,'MeHg_wddoc = ',MeHg_wddoc

    print*,'HgIIw = ', HgIIw
    print*,'HgII_wd = ', HgII_wd, 'HgII_wdoc = ', HgII_wdoc, 'HgII_wpa = ', HgII_wpa, 'HgII_wpom = ', HgII_wpom
    print*,'HgII_wpn = ', HgII_wpn
    print*,'HgIIw_veri = ', HgII_wd + HgII_wdoc + HgII_wpa + HgII_wpom + sum(HgII_wpn)
    print*,'HgIIs = ', HgIIs
    print*,'HgII_sd = ', HgII_sd, 'HgII_sdoc = ', HgII_sdoc, 'HgII_spom = ', HgII_spom
    print*,'HgII_spn = ', HgII_spn
    print*,'HgIIs_veri = ', HgII_sd + HgII_sdoc + HgII_spom + sum(HgII_spn)
    print*,'HgII_wddoc = ',HgII_wddoc

    print*,'Hg0w = ', Hg0w
    print*,'Hg0s = ', Hg0s
    print*,'kws = ',kws
    print*,'diff = ', kws * (MeHgs - MeHgw)
    print*,'MeHg_atmdep = ', MeHgw_atmdep
    print*,'HgII_atmdep = ', HgIIw_atmdep
    print*,'HgII_reduction =', HgIIw_reduction
    print*,'Hg0_oxidation =', Hg0w_oxidation
    print*,'MeHg_volatilization =', MeHgw_vol
    print*,'Hg0_volatilization =', Hg0w_vol
    print*,'MeHgw_diffusion = ', MeHgw_diffusion
    print*,'HgIIw_diffusion = ', HgIIw_diffusion
    print*,'Hg0w_diffusion = ', Hg0w_diffusion
    print*,'MeHgs_demethy = ', MeHgs_demethy
    print*,'HgIIs_methy = ', HgIIs_methy
    print*,'MeHgw_demethy = ', MeHgw_demethy
    print*,'HgIIw_methy = ', HgIIw_methy
    print*, 'HgIIs_erosion = ', HgIIs_erosion
    print*, 'HgIIw_deposition = ', HgIIw_deposition
    print*, 'HgIIs_burial = ', HgIIs_burial
    print*,'MeHgs_erosion = ', MeHgs_erosion
    print*,'MeHgw_deposition = ', MeHgw_deposition
    print*,'MeHgs_burial = ', MeHgs_burial
    print*,'MeHgw_photodeg = ', MeHgw_photodeg
    print*,'sourcesink MeHg= ', sourcesink(kwq, lwq, LMeHg)
    print*,'sourcesink HgII= ', sourcesink(kwq, lwq, LHgII)
    print*,'sourcesink Hg0= ', sourcesink(kwq, lwq, LHg0)
    if (kwq == kms) then
      print*,'sourcesink_sed MeHg= ',sourcesink(kwq + 1, lwq, LMeHg)
      print*,'sourcesink_sed HgII= ',sourcesink(kwq + 1, lwq, LHgII)
      print*,'sourcesink_sed Hg0= ',sourcesink(kwq + 1, lwq, LHg0)
    end if 
  end if

END SUBROUTINE sourceHg

!*************************************************************************
SUBROUTINE HgII_partitioning(kms, kwq, lwq, HgIIw, HgIIs, fwd2, fwdoc2, fwpa2, &
          & fwpom2, fwpn2, fsd2, fsdoc2, fspom2, fspn2, HgII_wpn, HgII_wd, &
          & HgII_wdoc, HgII_wpa, HgII_wpom, HgII_spn, HgII_sd, HgII_sdoc, HgII_spom)
!*************************************************************************
!  Purpose:
!
!*************************************************************************
  ! Arguments
  integer, intent(in)                     :: kms
  integer, intent(in)                     :: kwq
  integer, intent(in)                     :: lwq
  real,    intent(in)                     :: HgIIw
  real,    intent(in)                     :: HgIIs
  integer                                 :: i
  real                                    :: DOC
  real                                    :: ALG
  real                                    :: POM
  real,              dimension(sedNumber) :: HgII_SS
  real                                    :: R_HgIIw
  real                                    :: R_HgIIs
  real, intent(out)                       :: fwd2
  real, intent(out)                       :: fwdoc2
  real, intent(out)                       :: fwpa2
  real, intent(out)                       :: fwpom2
  real, intent(out), dimension(sedNumber) :: fwpn2
  real, intent(out)                       :: fsd2
  real, intent(out)                       :: fsdoc2
  real, intent(out)                       :: fspom2
  real, intent(out), dimension(sedNumber) :: fspn2
  real, intent(out), dimension(sedNumber) :: HgII_wpn 
  real, intent(out)                       :: HgII_wd
  real, intent(out)                       :: HgII_wdoc
  real, intent(out)                       :: HgII_wpa
  real, intent(out)                       :: HgII_wpom
  real, intent(out), dimension(sedNumber) :: HgII_spn
  real, intent(out)                       :: HgII_sd
  real, intent(out)                       :: HgII_sdoc
  real, intent(out)                       :: HgII_spom

  ! WATER COLUMN
  ! DOC = tracerpp(kwq, lwq, LDOC)
  ! ALG = tracerpp(kwq, lwq, LALG1)
  ! POM = tracerpp(kwq, lwq, LPOC)

  ! Values assumed for testing of partitioning
  DOC = 1.0 * 1000 * 1000
  ALG = 0.0
  POM = 0.0


  do i = 1, sedNumber
    HgII_SS(i) = kd_wpn2(i) * tracerpp(kwq, lwq, LSS1 + i - 1)
  end do

  R_HgIIw =  1 + (kd_wdoc2 * DOC) + (kd_wpa2 * ALG) + (kd_wpom2 * POM) + sum(HgII_SS)

  fwd2 = 1 / R_HgIIw
  fwdoc2 = kd_wdoc2 * DOC / R_HgIIw
  fwpa2 = kd_wpa2 * ALG / R_HgIIw
  fwpom2 = kd_wpom2 * POM / R_HgIIw

  do i = 1,sedNumber
    fwpn2(i) = kd_wpn2(i) * tracerpp(kwq, lwq, LSS1 + i - 1) / R_HgIIw
    HgII_wpn(i) = fwpn2(i) * HgIIw
  end do

  HgII_wd = fwd2 * HgIIw
  HgII_wdoc = fwdoc2 * HgIIw
  HgII_wpa = fwpa2 * HgIIw
  HgII_wpom = fwpom2 * HgIIw

  ! SEDIMENT LAYER
  if (kwq .eq. kms) then
    ! DOC = tracerpp(kwq + 1, lwq, LDOC)
    ! POM = tracerpp(kwq + 1, lwq, LPOC)

    DOC = 1.0 * 1000 * 1000
    POM = 0.0

    do i = 1, sedNumber
      HgII_SS(i) = kd_spn2(i) * tracerpp(kwq + 1, lwq, LSS1 + i - 1)
    end do

    R_HgIIs =  1 + (kd_sdoc2 * DOC) + (kd_spom2 * POM) + sum(HgII_SS)

    fsd2 = 1 / R_HgIIs
    fsdoc2 = kd_sdoc2 * DOC / R_HgIIs
    fspom2 = kd_spom2 * POM / R_HgIIs
    do i = 1,sedNumber
      fspn2(i) = kd_spn2(i) * tracerpp(kwq + 1, lwq, LSS1 + i - 1) / R_HgIIs
      HgII_spn(i) = fspn2(i) * HgIIs
    end do

    HgII_sd = fwd2 * HgIIs
    HgII_sdoc = fwdoc2 * HgIIs
    HgII_spom = fwpom2 * HgIIs
  end if 

END SUBROUTINE HgII_partitioning

!*************************************************************************
SUBROUTINE MeHg_partitioning(kms, kwq, lwq, MeHgw, MeHgs, fwd3, fwdoc3, fwpa3, fwpom3, fwpn3, &
          & fsd3, fsdoc3, fspom3, fspn3, MeHg_wpn, MeHg_wd, MeHg_wdoc, &
          & MeHg_wpa, MeHg_wpom, MeHg_spn, MeHg_sd, MeHg_sdoc, MeHg_spom)
!*************************************************************************
! Purposes: 
!
!
!*************************************************************************
  ! Arguments
  integer, intent(in)                     :: kms
  integer, intent(in)                     :: kwq
  integer, intent(in)                     :: lwq
  real,    intent(in)                     :: MeHgw
  real,    intent(in)                     :: MeHgs
  integer                                 :: i
  real                                    :: DOC
  real                                    :: ALG
  real                                    :: POM
  real,              dimension(sedNumber) :: MeHg_SS
  real                                    :: R_MeHgw
  real                                    :: R_MeHgs
  real, intent(out)                       :: fwd3
  real, intent(out)                       :: fwdoc3
  real, intent(out)                       :: fwpa3
  real, intent(out)                       :: fwpom3
  real, intent(out)                       :: fsd3
  real, intent(out)                       :: fsdoc3
  real, intent(out)                       :: fspom3
  real, intent(out)                       :: MeHg_wd
  real, intent(out)                       :: MeHg_wdoc
  real, intent(out)                       :: MeHg_wpa
  real, intent(out)                       :: MeHg_wpom
  real, intent(out)                       :: MeHg_sd
  real, intent(out)                       :: MeHg_sdoc
  real, intent(out)                       :: MeHg_spom
  real, intent(out), dimension(sedNumber) :: fwpn3
  real, intent(out), dimension(sedNumber) :: fspn3
  real, intent(out), dimension(sedNumber) :: MeHg_wpn
  real, intent(out), dimension(sedNumber) :: MeHg_spn

  ! WATER COLUMN
  ! DOC = tracerpp(kwq, lwq, LDOC)
  ! ALG = tracerpp(kwq, lwq, LALG1)
  ! POM = tracerpp(kwq, lwq, LPOC)

  ! Values assumed for testing of partitioning
  DOC = 1.0 * 1000 * 1000 ! mg/L in ng/L
  ALG = 0.0
  POM = 0.0

  do i = 1, sedNumber
    MeHg_SS(i) = kd_wpn3(i) * tracerpp(kwq, lwq, LSS1 + i - 1)
  end do

  R_MeHgw =  1 + kd_wdoc3 * DOC + kd_wpa3 * ALG + kd_wpom3 * POM + sum(MeHg_SS)
  
  fwd3 = 1 / R_MeHgw
  fwdoc3 = kd_wdoc3 * DOC / R_MeHgw
  fwpa3 = kd_wpa3 * ALG / R_MeHgw
  fwpom3 = kd_wpom3 * POM / R_MeHgw
  do i = 1,sedNumber
    fwpn3(i) = kd_wpn3(i) * tracerpp(kwq, lwq, LSS1 + i - 1) / R_MeHgw
    MeHg_wpn(i) = fwpn3(i) * MeHgw
  end do

  MeHg_wd = fwd3 * MeHgw
  MeHg_wdoc = fwdoc3 * MeHgw
  MeHg_wpa = fwpa3 * MeHgw
  MeHg_wpom = fwpom3 * MeHgw

  if (kwq .eq. kms) then
    ! DOC = tracerpp(kwq + 1, lwq, LDOC)
    ! POM = tracerpp(kwq + 1, lwq, LPOC)

    DOC = 1.0 * 1000 * 1000
    POM = 0.0

    do i = 1, sedNumber
      MeHg_SS(i) = kd_spn3(i) * tracerpp(kwq + 1, lwq, LSS1 + i - 1)
    end do

    R_MeHgs =  1 + kd_sdoc3 * DOC + kd_spom3 * POM + sum(MeHg_SS)

    fsd3 = 1 / R_MeHgs
    fsdoc3 = kd_sdoc3 * DOC / R_MeHgs
    fspom3 = kd_spom3 * POM / R_MeHgs
    do i = 1,sedNumber
      fspn3(i) = kd_spn3(i) * tracerpp(kwq + 1, lwq, LSS1 + i - 1) / R_MeHgs
      MeHg_spn(i) = fspn3(i) * MeHgs
    end do

    MeHg_sd = fsd3 * MeHgs
    MeHg_sdoc = fsdoc3 * MeHgs
    MeHg_spom = fspom3 * MeHgs
  end if 

END SUBROUTINE MeHg_partitioning

!*************************************************************************
SUBROUTINE HgII_reduction(HgIIw_reduction, kwq, lwq, HgII)
!*************************************************************************
!
!  Purpose: To estimate the reduction of HgII_w. This plays as a source
!           for Hg0_w and a sink for HgII_w. The process applies to both
!           the dissolved phase (_wd) and the adsorbed to DOC (_wdoc)
!
!-------------------------------------------------------------------------

  !Arguments
  integer, intent(in)  :: kwq
  integer, intent(in)  :: lwq
  real,    intent(in)  :: HgII          !< Concentration HgII in water (sum of dissolved and DOC adsorbed)
  real,    intent(out) :: HgIIw_reduction !< Source/Sink term for the balance of Hg0 and HgII
  real                 :: Io              !< Solar radiation at the surface

  ! Case select for different type of surface boundary conditions
  ! Current state of the model considers constant surface shortwave radiation over the whole surface regardless of the case selected
  
  select case (ifSurfbc)
  case (1, 2, 3)
    Io = Qsw
  case (10, 11)
    Io = Qsw
  end select

  HgIIw_reduction = h(kwq, lwq) * (Io * QswFr(kwq, lwq)) * kw21 * HgII

END SUBROUTINE HgII_reduction

!*************************************************************************
SUBROUTINE Hg0_oxidation(Hg0w_oxidation, Hg0w, HgIIwd, HgIIw_reduction, lwq)
!*************************************************************************
!
!  Purpose: To estimate the oxidation of Hg0_w. This plays as a source
!           for HgII_w and a sink for Hg0_w. The process is estimated as
!           instantaneous oxidation and it is estimated using a desired 
!           ratio following data for DOC and PH.
!
!-------------------------------------------------------------------------

  ! Arguments
  real, intent(in)    :: Hg0w
  real, intent(in)    :: HgIIwd
  real, intent(in)    :: HgIIw_reduction
  real, intent(out)   :: Hg0w_oxidation
  real                :: DGMr
  real                :: DGMra
  real                :: DOC
  real                :: PH
  integer, intent(in) :: lwq

  DGMr = Hg0w / HgIIwd
  DOC = 10.0 
  PH = 7.0
  DGMra = 0.00188 * DOC ** (-0.73385) * PH ** (2.0409)

  if (DGMr > DGMra) then
    Hg0w_oxidation = HgIIw_reduction * (DGMr / DGMra)
  else if (DGMr .le. DGMra) then
    Hg0w_oxidation = 0
  end if

  END SUBROUTINE Hg0_oxidation

!*************************************************************************
SUBROUTINE MeHg_photodegradation(MeHgw_photodeg, kwq, lwq, MeHg)
!*************************************************************************
!
!  Purpose: This estimates the photodegradation of MeHg_w. This plays as
!           a sink for MeHg and as a source for Hg0_w. The process applies
!           to both the dissolved (_wd) and DOC adsorbed (_wdoc) phases.
!
!-------------------------------------------------------------------------
  !Arguments
  integer, intent(in)  :: kwq
  integer, intent(in)  :: lwq
  real,    intent(in)  :: MeHg            !< Concentration MeHg in water (sum of dissolved and DOC adsorbed)
  real,    intent(out) :: MeHgw_photodeg  !< Source/Sink term for the balance of Hg0 and HgII
  real                 :: Io              !< Solar radiation at the surface
  
  ! Case select for different type of surface boundary conditions
  ! Current state of the model considers constant surface shortwave radiation over the whole surface regardless of the case selected
  select case (ifSurfbc)
  case (1, 2, 3)
    Io = Qsw
  case (10, 11)
    Io = Qsw
  end select

  MeHgw_photodeg = h(kwq, lwq) * (Io * QswFr(kwq, lwq)) * kw31 * MeHg

  ! if (lwq == 50) then
  !   print*,'kw31 = ',kw31,'MeHg_reduction = ',MeHgw_photodeg,'h = ',h(kwq, lwq)
  ! end if

END SUBROUTINE MeHg_photodegradation

!***********************************************************************
SUBROUTINE MeHg_atm_deposition(MeHgw_atmdep, kwq, lwq)
!***********************************************************************
!
!  Purpose: This estimates the atmospheric deposition of MeHg. This 
!           plays as a sink for MeHg_w and as a source for Hg0_w. The 
!           mercury enters as dissolved phase
!-----------------------------------------------------------------------
  ! Arguments
  integer, intent(in)  :: kwq
  integer, intent(in)  :: lwq
  real,    intent(out) :: MeHgw_atmdep !< Source/Sink term for the balance of MeHg dissolved in water

  MeHgw_atmdep = (dx * dy) * atm_MeHg * h(kwq, lwq)

END SUBROUTINE MeHg_atm_deposition

!***********************************************************************
SUBROUTINE HgII_atm_deposition(HgIIw_atmdep, kwq, lwq)
!***********************************************************************
!
!  Purpose: This estimates the atmospheric deposition of MeHg. This 
!           plays as a sink for MeHg_w and as a source for Hg0_w. The 
!           mercury enters as dissolved phase
!-----------------------------------------------------------------------
  ! Arguments
  integer, intent(in)  :: kwq
  integer, intent(in)  :: lwq
  real,    intent(out) :: HgIIw_atmdep !< Source/Sink term for the balance of MeHg dissolved in water

  HgIIw_atmdep = (dx * dy) * atm_HgII * h(kwq, lwq)

END SUBROUTINE HgII_atm_deposition

!***********************************************************************
SUBROUTINE Hg0_volatilization(Hg0w_vol, Hg0)
!***********************************************************************
!
!  Purpose: This estimates the atmospheric deposition of MeHg. This 
!           plays as a sink for MeHg_w and as a source for Hg0_w. The 
!           mercury enters as dissolved phase
!-----------------------------------------------------------------------
  ! Arguments
  real, intent(in)  :: Hg0
  real, intent(out) :: Hg0w_vol

  ! if (Hg0 - (Hg0atm / K_H_Hg0w) .lt. 0.0) then
  !   Hg0w_vol = 0.0
  ! else
    Hg0w_vol = (1 / ((1 / k_Hg0w) + (1 / (k_Hg0atm * K_H_Hg0w)))) * (Hg0 - (Hg0atm / K_H_Hg0w))
  ! end if

END SUBROUTINE Hg0_volatilization

!***********************************************************************
SUBROUTINE MeHg_volatilization(MeHgw_vol, MeHg)
!***********************************************************************
!
!  Purpose: This estimates the atmospheric deposition of MeHg. This 
!           plays as a sink for MeHg_w and as a source for Hg0_w. The 
!           mercury enters as dissolved phase
!-----------------------------------------------------------------------
  ! Arguments
  real, intent(in)  :: MeHg
  real, intent(out) :: MeHgw_vol

  ! if (MeHg - (MeHgatm / K_H_MeHgw) .lt. 0.0) then
  !   MeHgw_vol = 0.0
  ! else
    MeHgw_vol = (1 / ((1 / k_MeHgw) + (1 / (k_MeHgatm * K_H_MeHgw)))) * (MeHg - (MeHgatm / K_H_MeHgw))
  ! end if

END SUBROUTINE MeHg_volatilization

!***********************************************************************
SUBROUTINE HgIIw_methylation(HgIIw_methy, kwq, lwq, HgII)
!***********************************************************************
!
!  Purpose: This estimates the methylation of HgII in the water column 
!           and plays as a sink for HgII_w and as a source for MeHg_w.
!           It only applies to the dissolved phase and for both
!           dissolved and absorbed to DOC = HgII_w 
!-----------------------------------------------------------------------
  ! Arguments
  integer, intent(in)  :: kwq
  integer, intent(in)  :: lwq
  real,    intent(in)  :: HgII
  real,    intent(out) :: HgIIw_methy
  real                 :: DO_w
  real                 :: Canox
  real                 :: Q10
  real                 :: Q10_methyl
  real                 :: T_methyl

  
  Q10_methyl = 2.0 ! Source Reed
  T_methyl = 15.0
  Q10 = Q10_methyl ** ((salp(kwq, lwq) - T_methyl) / 10)

  if (kwq .lt. 9) then
    DO_w = 8.0 ! mg/L because KDO is input as mg/L it needs consistency for now. Pay attention how the variable is modeled. 
  else
    DO_w = 0.0
  end if
    
  ! DO_w = tracerpp(kwq, lwq, LDO)

  if (DO_w .lt. 4) then
    Canox = KDO / (DO_w + KDO)
  else
    Canox = 0.0
  end if

  HgIIw_methy = h(kwq, lwq) * kw23 * Q10 * Canox * HgII

END SUBROUTINE HgIIw_methylation

!***********************************************************************
SUBROUTINE MeHgw_demethylation(MeHgw_demethy, kwq, lwq, MeHg)
!***********************************************************************
!
!  Purpose: This estimates the demethylation of MeHg in the water column 
!           and plays as a sink for MeHg_w and as a source for HgII_w.
!           It only applies to the dissolved phase and for both
!           dissolved and absorbed to DOC = MeHg_w 
!-----------------------------------------------------------------------
  ! Arguments
  integer, intent(in)  :: kwq
  integer, intent(in)  :: lwq
  real,    intent(in)  :: MeHg
  real,    intent(out) :: MeHgw_demethy
  real                 :: DO_w
  real                 :: Canox
  
  ! DO_w = tracerpp(kwq, lwq, LDO)
  Canox = 1.0

  MeHgw_demethy = h(kwq, lwq) * kw32 * Canox * MeHg

END SUBROUTINE MeHgw_demethylation

!***********************************************************************
SUBROUTINE MeHg_diffusion(MeHgw_diffusion, MeHgw, MeHgs)
!***********************************************************************
!
!  Purpose: This estimates the diffusion of MeHg in the water column 
!           and plays as a sink for MeHg_s and as a source for MeHg_w.
!           It means that for a possitive flux MeHg is moved from 
!           sediment layer into the water column 
!-----------------------------------------------------------------------
  ! Arguments
  real, intent(in)  :: MeHgw
  real, intent(in)  :: MeHgs
  real, intent(out) :: MeHgw_diffusion

  MeHgw_diffusion = kws * (MeHgs - MeHgw)

END SUBROUTINE MeHg_diffusion

!***********************************************************************
SUBROUTINE HgII_diffusion(HgIIw_diffusion, HgIIw, HgIIs)
!***********************************************************************
!
!  Purpose: This estimates the diffusion of HgII in the water column 
!           and plays as a sink for HgII_s and as a source for HgII_w.
!           It means that for a possitive flux HgII is moved from 
!           sediment layer into the water column 
!-----------------------------------------------------------------------
  ! Arguments
  real, intent(in)  :: HgIIw
  real, intent(in)  :: HgIIs
  real, intent(out) :: HgIIw_diffusion

  HgIIw_diffusion = kws * (HgIIs - HgIIw)

END SUBROUTINE HgII_diffusion

!***********************************************************************
SUBROUTINE Hg0_diffusion(Hg0w_diffusion, Hg0w, Hg0s)
!***********************************************************************
!
!  Purpose: This estimates the diffusion of Hg0 in the water column 
!           and plays as a sink for Hg0_s and as a source for Hg0_w.
!           It means that for a possitive flux Hg0 is moved from 
!           sediment layer into the water column 
!-----------------------------------------------------------------------
  ! Arguments
  real, intent(in)  :: Hg0w
  real, intent(in)  :: Hg0s
  real, intent(out) :: Hg0w_diffusion

  Hg0w_diffusion = kws * (Hg0s - Hg0w)

END SUBROUTINE Hg0_diffusion

! ********** FROM NOW ON ARE THE PROCESSES IN SEDIMENT LAYER *************

!***********************************************************************
SUBROUTINE HgIIs_methylation(HgIIs_methy, kwq, lwq, HgII)
!***********************************************************************
!
!  Purpose: This estimates the methylation of HgII in the sediment layer 
!           and plays as a sink for HgII_s and as a source for MeHg_s.
!           It only applies to the dissolved phase and for both
!           dissolved and absorbed to DOC = HgII_s 
!-----------------------------------------------------------------------
  ! Arguments
  integer, intent(in)  :: kwq
  integer, intent(in)  :: lwq
  real,    intent(in)  :: HgII
  real,    intent(out) :: HgIIs_methy
  real                 :: Q10
  real                 :: Q10_methyl
  real                 :: T_methyl

  Q10_methyl = 2.6 ! From Urban, N. R., Brezonik, P. L., Baker, L. A., & Sherman, L. A. (1994). Sulfate reduction and diffusion in sediments of Little Rock Lake, Wisconsin. Limnology and Oceanography, 39(4), 797–815. https://doi.org/10.4319/lo.1994.39.4.0797

  T_methyl = 15

  Q10 = Q10_methyl ** ((salp(kwq,lwq) - T_methyl) / 10)
  HgIIs_methy = h(kwq, lwq) * ks23 * Q10 * (1 + miu_so4 * (SO4 / (KSO4 + SO4)) ) * HgII

END SUBROUTINE HgIIs_methylation

!***********************************************************************
SUBROUTINE MeHgs_demethylation(MeHgs_demethy, kwq, lwq, MeHg)
!***********************************************************************
!
!  Purpose: This estimates the demethylation of MeHg in the sediment layer 
!           and plays as a sink for MeHg_s and as a source for HgII_s.
!           It only applies to the dissolved phase and for both
!           dissolved and absorbed to DOC = MeHg_s 
!-----------------------------------------------------------------------
  ! Arguments
  integer, intent(in)  :: kwq
  integer, intent(in)  :: lwq
  real,    intent(in)  :: MeHg
  real,    intent(out) :: MeHgs_demethy
  real                 :: Q10
  real                 :: Q10_demethyl
  real                 :: T_demethyl

  Q10_demethyl = 2.6 ! From Urban, N. R., Brezonik, P. L., Baker, L. A., & Sherman, L. A. (1994). Sulfate reduction and diffusion in sediments of Little Rock Lake, Wisconsin. Limnology and Oceanography, 39(4), 797–815. https://doi.org/10.4319/lo.1994.39.4.0797
  T_demethyl = 15

  Q10 = Q10_demethyl ** ((salp(kwq,lwq) - T_demethyl) / 10)

  MeHgs_demethy = h(kwq, lwq) * ks32 * Q10 * MeHg

END SUBROUTINE MeHgs_demethylation

!***********************************************************************
SUBROUTINE HgII_deposition(HgIIw_deposition, HgII_wpa, HgII_wpom, HgII_wpn)
!***********************************************************************
!
!  Purpose: This estimates the deposition of HgII from the water column 
!           and plays as a sink for HgII_wp and as a source for MeHg_sp.
!           It only applies to the particulate phase
!-----------------------------------------------------------------------
  ! Arguments
  real, intent(in)                       :: HgII_wpa
  real, intent(in)                       :: HgII_wpom
  real, intent(in), dimension(sedNumber) :: HgII_wpn
  real, intent(out)                      :: HgIIw_deposition
  real,             dimension(sedNumber) :: settling_pn
  integer                                :: i

  if (inst_eq .eq. 1) then
    do i = 1, sedNumber
      settling_pn(i) = settling_vel(i) * HgII_wpn(i)
    end do
  end if

  HgIIw_deposition = vspa * HgII_wpa + vspoc * HgII_wpom + sum(settling_pn)

END SUBROUTINE HgII_deposition

!***********************************************************************
SUBROUTINE MeHg_deposition(MeHgw_deposition, MeHg_wpa, MeHg_wpom, MeHg_wpn)
!***********************************************************************
!
!  Purpose: This estimates the deposition of MeHg from the water column 
!           and plays as a sink for MeHg_wp and as a source for MeHg_sp.
!           It only applies to the particulate phase
!-----------------------------------------------------------------------
  ! Arguments
  real, intent(in)                       :: MeHg_wpa
  real, intent(in)                       :: MeHg_wpom
  real, intent(in), dimension(sedNumber) :: MeHg_wpn
  real, intent(out)                      :: MeHgw_deposition
  real,             dimension(sedNumber) :: settling_pn
  integer                                :: i

  if (inst_eq .eq. 1) then
    do i = 1, sedNumber
      settling_pn(i) = settling_vel(i) * MeHg_wpn(i)
    end do
  end if

  MeHgw_deposition = vspa * MeHg_wpa + vspoc * MeHg_wpom + sum(settling_pn)

END SUBROUTINE MeHg_deposition

!************************************************************************
SUBROUTINE HgII_erosion(HgIIs_erosion, HgII_spn)
!************************************************************************
!
!   Purpose: To estimate the erosion of HgII adsorbed to sediments.
!
!
!------------------------------------------------------------------------

  ! Arguments
  real, intent(in), dimension(sedNumber) :: HgII_spn
  real, intent(out)                      :: HgIIs_erosion
  real, dimension(sedNumber)             :: flux
  integer                                :: i

  do i = 1, sedNumber
    flux(i) = erosion_Hgpn(i) * HgII_spn(i)
  end do

  HgIIs_erosion = sum(flux)

END SUBROUTINE HgII_erosion

!************************************************************************
SUBROUTINE MeHg_erosion(MeHgs_erosion, MeHg_spn)
!************************************************************************
!
!   Purpose: To estimate erosion of MeHg adsorbed to sediments
!
!
!------------------------------------------------------------------------

  ! Arguments
  real, intent(in), dimension(sedNumber) :: MeHg_spn
  real, intent(out)                      :: MeHgs_erosion
  real, dimension(sedNumber)             :: flux
  integer                                :: i

  do i = 1, sedNumber
    flux(i) = erosion_Hgpn(i) * MeHg_spn(i)
  end do

  MeHgs_erosion = sum(flux)

END SUBROUTINE MeHg_erosion

!************************************************************************
SUBROUTINE MeHg_burial(MeHgs_burial, MeHgs_erosion, MeHgw_deposition)
!************************************************************************
!
!   Purpose: To estimate erosion of MeHg adsorbed to sediments
!
!
!------------------------------------------------------------------------

  ! Arguments
  real, intent(in)  :: MeHgw_deposition
  real, intent(in)  :: MeHgs_erosion
  real, intent(out) :: MeHgs_burial

  MeHgs_burial = MeHgw_deposition - MeHgs_erosion

END SUBROUTINE MeHg_burial

!************************************************************************
SUBROUTINE HgII_burial(HgIIs_burial, HgIIs_erosion, HgIIw_deposition)
!************************************************************************
!
!   Purpose: To estimate erosion of MeHg adsorbed to sediments
!
!
!------------------------------------------------------------------------

  ! Arguments
  real, intent(in)  :: HgIIw_deposition
  real, intent(in)  :: HgIIs_erosion
  real, intent(out) :: HgIIs_burial

  HgIIs_burial = HgIIw_deposition - HgIIs_erosion

END SUBROUTINE HgII_burial

!************************************************************************
                        END MODULE si3d_Hg
!************************************************************************





































! !***********************************************************************
! SUBROUTINE HgIIw_ads_des(kwq, lwq, HgIIw_ad_de)
! !***********************************************************************
! !
! !  Purpose: This estimates the adsorption/desorption of HgII in the 
! !           water column. It includes adsorption to algae, pom, and SS 
! !           Plays as a sink for HgII_w and a source for HgII_wp.
! !-----------------------------------------------------------------------  

!   ! Arguments
!   integer, intent(in) :: kwq
!   integer, intent(in) :: lwq
!   real, intent(out)   :: HgIIw_ad_de
!   real                :: HgII_w
!   real                :: HgII_wd
!   real                :: HgII_wdoc
!   real                :: HgII_wpa
!   real                :: HgII_wpom
!   real                :: HgII_wpn
!   real                :: HgII_wp


!   if (inst_eq .eq. 1) then
!     HgII_w = tracerpp(kwq, lwq, LHgII)
!     HgII_wd = HgII_w / (1 + kd_wdoc2)
!     HgII_wdoc = kd_wdoc2 * HgII_wd
!     HgII_wpa = kd_wpa2 * HgII_wd
!     HgII_wpom = kd_wpom2 * HgII_wd
!     do i = 1, sedNumber
!       HgII_wpn(i) = kd_wpn2(i) * HgII_wd
!     end do
!     HgII_wp = sum(HgII_wpn)
!   end if
!   HgIIw_ad_de = HgII_wpa + HgII_wpom + HgII_wp

! END SUBROUTINE HgIIw_ads_des

! !***********************************************************************
! SUBROUTINE MeHgw_ads_des(kwq, lwq, MeHgw_ad_de)
! !***********************************************************************
! !
! !  Purpose: This estimates the adsorption/desorption of HgII in the 
! !           water column. It includes adsorption to algae, pom, and SS 
! !           Plays as a sink for HgII_w and a source for HgII_wp.
! !-----------------------------------------------------------------------  

!   ! Arguments
!   integer, intent(in) :: kwq
!   integer, intent(in) :: lwq
!   real, intent(out)   :: MeHgw_ad_de
!   real                :: MeHg_w
!   real                :: MeHg_wd
!   real                :: MeHg_wdoc
!   real                :: MeHg_wpa
!   real                :: MeHg_wpom
!   real                :: MeHg_wpn
!   real                :: MeHg_wp


!   if (inst_eq .eq. 1) then
!     MeHg_w = tracerpp(kwq, lwq, LMeHg)
!     MeHg_wd = MeHg_w / (1 + kd_wdoc3)
!     MeHg_wdoc = kd_wdoc3 * MeHg_wd
!     MeHg_wpa = kd_wpa3 * MeHg_wd
!     MeHg_wpom = kd_wpom3 * MeHg_wd
!     do i = 1, sedNumber
!       MeHg_wpn(i) = kd_wpn3(i) * MeHg_wd
!     end do
!     MeHg_wp = sum(MeHg_wpn)
!   end if
!   MeHgw_ad_de = MeHg_wpa + MeHg_wpom + MeHg_wp 

! END SUBROUTINE MeHgw_ads_des
























! !***********************************************************************
! SUBROUTINE HgIIs_ads_des(kwq, lwq, HgIIs_ad_de)
! !***********************************************************************
! !
! !  Purpose: This estimates the adsorption/desorption of HgII in the 
! !           sediment layer. It includes adsorption to algae, pom, and SS 
! !           Plays as a sink for HgII_w and a source for HgII_wp.
! !-----------------------------------------------------------------------  

!   ! Arguments
!   integer, intent(in) :: kwq
!   integer, intent(in) :: lwq
!   real, intent(out)   :: HgIIs_ad_de
!   real                :: HgII_s
!   real                :: HgII_sd
!   real                :: HgII_sdoc
!   real                :: HgII_spa
!   real                :: HgII_spom
!   real                :: HgII_spn
!   real                :: HgII_sp


!   if (inst_eq .eq. 1) then
!     HgII_s = tracerpp(kwq+1, lwq, LHgII)
!     HgII_sd = HgII_s / (1 + kd_sdoc2)
!     HgII_sdoc = kd_sdoc2 * HgII_sd
!     HgII_spa = kd_spa2 * HgII_sd
!     HgII_spom = kd_spom2 * HgII_sd
!     do i = 1, sedNumber
!       HgII_spn(i) = kd_spn2(i) * HgII_sd
!     end do
!     HgII_sp = sum(HgII_spn)
!   end if
!   HgIIs_ad_de = HgII_spa + HgII_spom + HgII_sp

! END SUBROUTINE HgIIs_ads_des

! !***********************************************************************
! SUBROUTINE MeHgs_ads_des(kwq, lwq, MeHgs_ad_de)
! !***********************************************************************
! !
! !  Purpose: This estimates the adsorption/desorption of HgII in the 
! !           sediment layer. It includes adsorption to algae, pom, and SS 
! !           Plays as a sink for HgII_w and a source for HgII_wp.
! !-----------------------------------------------------------------------  

!   ! Arguments
!   integer, intent(in) :: kwq
!   integer, intent(in) :: lwq
!   real, intent(out)   :: MeHgs_ad_de
!   real                :: MeHg_s
!   real                :: MeHg_sd
!   real                :: MeHg_sdoc
!   real                :: MeHg_spa
!   real                :: MeHg_spom
!   real                :: MeHg_spn
!   real                :: MeHg_sp


!   if (inst_eq .eq. 1) then
!     MeHg_s = tracerpp(kwq, lwq, LMeHg)
!     MeHg_sd = MeHg_s / (1 + kd_sdoc3)
!     MeHg_sdoc = kd_sdoc3 * MeHg_sd
!     MeHg_spa = kd_spa3 * MeHg_sd
!     MeHg_spom = kd_spom3 * MeHg_sd
!     do i = 1, sedNumber
!       MeHg_spn(i) = kd_spn3(i) * MeHg_sd
!     end do
!     MeHg_sp = sum(MeHg_spn)
!   end if
!   MeHgs_ad_de = MeHg_spa + MeHg_spom + MeHg_sp 

! END SUBROUTINE MeHgs_ads_des






















































! !*************************************************************************
! SUBROUTINE sourceHg0(kwq, lwq)
! !*************************************************************************
! !
! !   Purpose: This subroutine calls the different processes implemented
! !           for the modeling of elemental mercury Hg0, and calculates 
! !           the sources and sinks for Hg0
! !-------------------------------------------------------------------------

!   ! ... Arguments
!   integer, intent(in) :: kwq
!   integer, intent(in) :: lwq
!   integer             :: kms
!   integer             :: k1s
!   real                :: HgIIw_reduction = 0.0
!   real                :: MeHgw_photodeg = 0.0
!   real                :: Hg0w_diffusion = 0.0

!   kms = kmz(lwq)
!   k1s = k1z(lwq)

!   HgIIw_reduction = 0.0
!   MeHgw_photodeg = 0.0
!   Hg0w_diffusion = 0.0
 
!   ! call HgII_reduction(kwq, lwq, HgIIw_reduction)
!   ! call Hg0_oxidation 
!   ! call MeHg_photodegradation(kwq, lwq, MeHgw_photodeg)

!   if (kwq .eq. k1s) then
!     ! call Hg0_volatilization
!   end if 
!   ! call HgII_atm_deposition
!   ! call MeHg_atm_deposition
!   ! call MeHg_volatilization
!   ! call Hg0_volatilization
!   ! call HgII_methylation
!   ! call MeHg_demethylation
!   ! call HgII_adsorption_desorption
!   ! call MeHg_adsorption_desorption
!   ! call MeHg_erosion
!   ! call MeHg_deposition
!   ! call HgII_erosion
!   ! call HgII_deposition
  
!   ! HgIIw_reduction + MeHgw_photodeg

!   if (kwq .eq. kms) then
!     call Hg0_diffusion(kwq, lwq, Hg0w_diffusion)
!     ! Section to call on the subroutines that calculate sources and sink of bottom sediment layer
!     sourcesink(kwq + 1, lwq, LHg0) = -Hg0w_diffusion
!   end if 

!   if (lwq == 50) then
!     print*,'Hg0'
!     print*,'kwq = ',kwq,'kms = ',kms
!     print*,'erosion_flux',erosion_flux(:)
!     print*,'settling_vel',settling_vel(:)
!     print*,'diffusion_flux',Hg0w_diffusion
!   end if
!   sourcesink(kwq, lwq, LHg0) = HgIIw_reduction + MeHgw_photodeg + Hg0w_diffusion

! END SUBROUTINE sourceHg0

! !*************************************************************************
! SUBROUTINE sourceHgII(kwq,lwq)
! !*************************************************************************
! !
! !  Purpose: This subroutine calls the different processes implemented
! !           for the modeling of divalent mercury HgII, and calculates 
! !           the sources and sinks for HgII
! !-------------------------------------------------------------------------

!   ! ... Arguments
!   integer, intent(in)        :: kwq
!   integer, intent(in)        :: lwq
!   integer                    :: kms
!   integer                    :: k1s
!   real                       :: HgIIw_reduction
!   real                       :: HgIIw_atmdep
!   real                       :: HgIIw_diffusion

!   real                       :: fwd2
!   real                       :: fwdoc2
!   real                       :: fwpa2
!   real                       :: fwpom2
  
!   real                       :: fsd2
!   real                       :: fsdoc2
!   real                       :: fspom2
!   real                       :: HgII_wd
!   real                       :: HgII_wdoc
!   real                       :: HgII_wpa
!   real                       :: HgII_wpom
!   real                       :: HgII_sd
!   real                       :: HgII_sdoc
!   real                       :: HgII_spom
!   real, dimension(sedNumber) :: fwpn2
!   real, dimension(sedNumber) :: fspn2
!   real, dimension(sedNumber) :: HgII_wpn
!   real, dimension(sedNumber) :: HgII_spn

!   real                       :: HgII !< Estimated variable used in subroutines. Changes based on process to represent.

!   kms = kmz(lwq)
!   k1s = k1z(lwq)

!   HgIIw_reduction = 0.0
!   HgIIw_atmdep = 0.0
!   HgIIw_diffusion = 0.0


!   call HgII_partitioning(kms, kwq, lwq, fwd2, fwdoc2, fwpa2, fwpom2, fwpn2, &
!        & fsd2, fsdoc2, fspom2, fspn2, HgII_wpn, HgII_wd, HgII_wdoc, HgII_wpa, &
!        & HgII_wpom, HgII_spn, HgII_sd, HgII_sdoc, HgII_spom)

!   ! Estimate Dissolved mercury in water
!   HgII = HgII_wd + HgII_wdoc

!   ! call HgII_reduction(kwq, lwq, HgII, HgIIw_reduction)
!   ! call Hg0_oxidation(kwq, lwq, Hg0w_to_HgIIw) 
!   ! call MeHg_photodegradation

!   if (kwq .eq. k1s) then
!     call HgII_atm_deposition(kwq, lwq, HgIIw_atmdep)
!   end if 

!   ! call MeHg_volatilization
!   ! call Hg0_volatilization
!   ! call HgII_methylation
!   ! call MeHg_demethylation
!   ! call HgII_adsorption_desorption
!   ! call MeHg_adsorption_desorption
!   ! call MeHg_erosion
!   ! call MeHg_deposition
!   ! call HgII_erosion
!   ! call HgII_deposition

!   if (kwq .eq. kms) then
!     call HgII_diffusion(kwq, lwq, HgIIw_diffusion)
!     ! Section to call on the subroutines that calculate sources and sink of bottom sediment layer
!     sourcesink(kwq + 1, lwq, LHgII) = -HgIIw_diffusion
!   end if

!   if (lwq == 50) then
!     print*,'HgII'
!     print*,'kwq = ',kwq,'kms = ',kms
!     print*,'erosion_flux',erosion_flux(:)
!     print*,'settling_vel',settling_vel(:)
!     print*,'diffusion_flux',HgIIw_diffusion
!     print*,'atm_dep',HgIIw_atmdep
!   end if
!   sourcesink(kwq, lwq, LHgII) = -HgIIw_reduction + HgIIw_atmdep + HgIIw_diffusion
  

! END SUBROUTINE sourceHgII

! !*************************************************************************
! SUBROUTINE sourceMeHg(kwq,lwq)
! !*************************************************************************
! !
! !  Purpose: This subroutine calls the different processes implemented
! !           for the modeling of divalent mercury HgII, and calculates 
! !           the sources and sinks for HgII
! !-------------------------------------------------------------------------

!   ! ... Arguments
!   integer, intent(in)        :: kwq
!   integer, intent(in)        :: lwq
!   integer                    :: kms
!   integer                    :: k1s
!   real                       :: MeHgw_photodeg
!   real                       :: MeHgw_atmdep
!   real                       :: MeHgw_diffusion
!   real                       :: fwd3
!   real                       :: fwdoc3
!   real                       :: fwpa3
!   real                       :: fwpom3 
!   real                       :: fsd3
!   real                       :: fsdoc3
!   real                       :: fspom3
!   real                       :: MeHg_wd
!   real                       :: MeHg_wdoc
!   real                       :: MeHg_wpa
!   real                       :: MeHg_wpom
!   real                       :: MeHg_sd
!   real                       :: MeHg_sdoc
!   real                       :: MeHg_spom
!   real, dimension(sedNumber) :: fwpn3
!   real, dimension(sedNumber) :: fspn3
!   real, dimension(sedNumber) :: MeHg_wpn
!   real, dimension(sedNumber) :: MeHg_spn

!   kms = kmz(lwq)
!   k1s = k1z(lwq)

!   MeHgw_photodeg = 0.0
!   MeHgw_diffusion = 0.0
!   MeHgw_atmdep = 0.0

!   call MeHg_partitioning(kms, kwq, lwq, fwd3, fwdoc3, fwpa3, fwpom3, fwpn3, &
!           & fsd3, fsdoc3, fspom3, fspn3, MeHg_wpn, MeHg_wd, MeHg_wdoc, &
!           & MeHg_wpa, MeHg_wpom, MeHg_spn, MeHg_sd, MeHg_sdoc, MeHg_spom)

!   ! call HgII_reduction(kwq, lwq, HgIIw_to_HgOw)
!   ! call Hg0_oxidation 
!   ! call MeHg_photodegradation(kwq, lwq, MeHgw_photodeg)

!   if (kwq .eq. k1s) then
!     call MeHg_atm_deposition(kwq, lwq, MeHgw_atmdep)
!   end if 
!   ! call MeHg_volatilization
!   ! call Hg0_volatilization
!   ! call HgII_methylation
!   ! call MeHg_demethylation
!   ! call HgII_adsorption_desorption
!   ! call MeHg_adsorption_desorption
!   ! call MeHg_erosion
!   ! call MeHg_deposition
!   ! call HgII_erosion
!   ! call HgII_deposition

!   if (kwq .eq. kms) then
!     ! Section to call on the subroutines that calculate sources and sink of bottom sediment layer
!     ! call MeHgp_settling(kwq, lwq, MeHgs_settling)
!     call MeHg_diffusion(kwq, lwq, MeHgw_diffusion)
!     sourcesink(kwq + 1, lwq, LMeHg) = -MeHgw_diffusion
!   end if

!   if (lwq == 50) then
!       print*,'MeHg'
!       print*,'kwq = ',kwq,'kms = ',kms,'k1s =',k1s
!       print*,'erosion_flux',erosion_flux(:)
!       print*,'settling_vel',settling_vel(:)
!       print*,'diffusion_flux',MeHgw_diffusion
!       print*,'MeHgw_atmdep',MeHgw_atmdep
!   end if

!   sourcesink(kwq, lwq, LMeHg) = -MeHgw_photodeg + MeHgw_atmdep + MeHgw_diffusion

! END SUBROUTINE sourceMeHg