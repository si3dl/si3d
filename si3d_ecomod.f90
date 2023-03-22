!************************************************************************
                          MODULE si3d_ecomod
!************************************************************************
!
!  Purpose: Procedures that implement routines related to the modelling
!           of ecological processees
!
!-------------------------------------------------------------------------

   USE omp_lib
   USE si3d_types
   USE si3d_sed
   USE si3d_stwave

   IMPLICIT NONE
   SAVE

CONTAINS


!***********************************************************************
SUBROUTINE srcsnk00
!***********************************************************************
!
!  Purpose: Set sources and sinks to zero - used when conservative
!           tracers (ecomod == 0) are modelled
!
!-----------------------------------------------------------------------

  sourcesink = 0.0E0

END SUBROUTINE srcsnk00

!***********************************************************************
SUBROUTINE TRCinput
!***********************************************************************
!
!  Purpose: To read in values of peak concentrations, location of peaks
!           and dispersion of clouds -
!
!-----------------------------------------------------------------------

  !. . . .Local Variables
  CHARACTER(LEN=12) :: trciofile="si3d_tr.txt"
  CHARACTER(LEN=18) :: iotrcfmt
  INTEGER:: ios, itr ,idep, istat
  INTEGER, PARAMETER :: niotrc = 9

  !.....Open input parameter file.....
  OPEN (UNIT=i99, FILE=trciofile, STATUS="OLD",  FORM="FORMATTED", IOSTAT=ios)
  IF (ios /= 0) CALL open_error ( "Error opening "//trciofile, ios )

  ! ... Allocate space for arrays holding information on size groups

  ALLOCATE ( trct0(ntr), trcpk(ntr), trctn(ntr), &
             trcx0(ntr), trcy0(ntr), trcz0(ntr), &
             trcsx(ntr), trcsy(ntr), trcsz(ntr), STAT=istat)

  IF (istat /= 0) CALL allocate_error ( istat, 121 )

  ! Skip over first five header records in SIZE file
  READ (UNIT=i99, FMT='(/////)', IOSTAT=ios)
  IF (ios /= 0) CALL input_error ( ios, 90 )

  ! Write the format of the data records into an internal file
  WRITE (UNIT=iotrcfmt, FMT='("(10X,",I3,"G11.2)")') niotrc

  ! Read information for each size group (diameter, growth rate & deposition)
  DO itr = 1, ntr
    READ (UNIT=i99, FMT=iotrcfmt, IOSTAT=ios) &
         trct0(itr), trcpk(itr), trctn(itr),  &
         trcx0(itr), trcy0(itr), trcz0(itr),  &
         trcsx(itr), trcsy(itr), trcsz(itr)
    PRINT *, itr, trct0(itr), trcpk(itr), trctn(itr),trcx0(itr), trcy0(itr), trcz0(itr),trcsx(itr), trcsy(itr), trcsz(itr)
    IF (ios /= 0) CALL input_error ( ios, 91 )
  END DO

  CLOSE (i99)

END SUBROUTINE TRCinput

!***********************************************************************
SUBROUTINE InitTracerCloud
!***********************************************************************

  INTEGER:: itr, i, j, k, l
  REAL   :: x, y, z

  DO itr = 1, ntr

    DO l = 1, lm
      i = l2i(l);
      j = l2j(l);
      DO k = k1,kmz(l)
        x = ( REAL ( i - i1 ) + 0.5 ) * dx - trcx0 (itr)
        y = ( REAL ( j - j1 ) + 0.5 ) * dy - trcy0 (itr)
        z = ( REAL ( k - k1 ) + 0.5 ) * ddz- trcz0 (itr)
        tracer  (k,l,itr) = trcpk (itr) *                        &
          & EXP( - x ** 2. / ( 2 * trcsx(itr) ** 2. ) ) *        &
          & EXP( - y ** 2. / ( 2 * trcsy(itr) ** 2. ) ) *        &
          & EXP( - z ** 2. / ( 2 * trcsz(itr) ** 2. ) )
      ENDDO
    ENDDO

  ENDDO

END SUBROUTINE InitTracerCloud

!***********************************************************************
SUBROUTINE SDinput
!***********************************************************************
!
!  Purpose: To read in values of parameters for sediment transport
!           model E = alpha*(tau-taucr/taucr)^m & vdep
!
!-----------------------------------------------------------------------

  !. . . .Local Variables
  CHARACTER(LEN=12) :: sdiofile="si3d_sed.txt"
  CHARACTER(LEN=18) :: iosdfmt
  INTEGER:: ios, itr , idep, istat
  INTEGER, PARAMETER :: niosd = 3

  !.....Open input parameter file.....
  OPEN (UNIT=i99, FILE=sdiofile, STATUS="OLD",  FORM="FORMATTED", IOSTAT=ios)
  IF (ios /= 0) CALL open_error ( "Error opening "//sdiofile, ios )

  ! ... Allocate space for arrays holding information for sediments
  ALLOCATE ( sdenSED(ntr), diamSED(ntr), vdep(ntr), STAT=istat)
  IF (istat /= 0) CALL allocate_error ( istat, 99 )

  ! ....Skip over first five header records in SEDIMENT file
  READ (UNIT=i99, FMT='(/////)', IOSTAT=ios)
  IF (ios /= 0) CALL input_error ( ios, 90 )

  ! ....Read hours between consecutive WIBSS fields
  READ (UNIT=i99, FMT='(10X,G11.2)', IOSTAT=ios) dthrsWIBSS
  IF (ios /= 0) CALL input_error ( ios, 90 )

  ! ....Read hours from start to first WIBSS field
  READ (UNIT=i99, FMT='(10X,G11.2)', IOSTAT=ios) thrsWIBSS0
  IF (ios /= 0) CALL input_error ( ios, 90 )

  ! ....Skip over one header records in SEDIMENT file
  READ (UNIT=i99, FMT='(//)', IOSTAT=ios)
  IF (ios /= 0) CALL input_error ( ios, 90 )

  ! ....Write the format of the data records into an internal file
  WRITE (UNIT=iosdfmt, FMT='("(10X,",I3,"G11.2)")') niosd

  ! ....Read information for each sediment category
  DO itr = 1, ntr
    READ (UNIT=i99, FMT=iosdfmt, IOSTAT=ios) &
         diamSED(itr), sdenSED(itr), vdep(itr)
    PRINT *, itr, sdenSED(itr), diamSED(itr), vdep(itr)
    IF (ios /= 0) CALL input_error ( ios, 91 )
  END DO

  CLOSE (i99)

  ! ... Set thrsWIBSS to zero
  thrsWIBSSp = 0.0
  thrsWIBSS  = 0.0


END SUBROUTINE SDinput

!***********************************************************************
SUBROUTINE srcsnkSD
!***********************************************************************
!
!  Purpose: To evaluate the source-sink terms in the reactive-transport
!           equations for sediments - the only sources considered
!           occur at the bottom cell due to resuspension E
!
!-----------------------------------------------------------------------

   ! ... Local variables
   INTEGER :: i, j, k, l, k1s, kms, itr, kmyp, kmym, kmxp, kmxm
   REAL    :: Rep, fRep, Zu5, ubott, vbott, ustar, taubx, tauby, taub

   ! ... Initialize to zero on first call
   IF (n .LE. 1) sourcesink = 0.0E0

   ! ... Calculate Wave-Induced Bottom Shear Stress WIBSS fields
   !     (computed with stwave in pre-process mode)
   CALL CalculateWIBSS

   ! ... Loop over tracers
   DO itr = 1, ntr

     ! ... Calculate particle Reynolds number
     Rep  = sqrt(g*sdenSED(itr)*(diamSED(itr)**3.))/1.E-6;
     IF      ( Rep .GT. 1. ) THEN
       fRep = 0.586*(Rep**1.23);
     ELSE IF ( Rep .LE. 1. ) THEN
       fRep = Rep**3.75;
     ENDIF
     IF      ( Rep .GT. 3. ) THEN
       PRINT *, '!!!!!!!!!!!! WARNING !!!!!!!!!!!'
       PRINT *, 'Rep is not within range of      '
       PRINT *, 'validity of sediment model      '
       PRINT *, 'Size Class = ', itr
       PRINT *, 'Rep = ', Rep
     ENDIF

     ! ... Loop over wet columns
     DO l = 1, lm

       ! ... 3D-(i,j) indexes for l ...........
       i = l2i(l); j = l2j(l);

       !.....Compute bottom layer numbers  ....
       kms = kmz(l);
       kmyp= MIN(kmz(lNC(l)), kmz(l))
       kmym= MIN(kmz(lSC(l)), kmz(l))
       kmxp= MIN(kmz(lEC(l)), kmz(l))
       kmxm= MIN(kmz(lWC(l)), kmz(l))

       ! ....Compute Currend-Induced Bottom Shear Stress CIBSS (cd*rho*U^2)
       ubott = (uhp(kmxp,l)+uhp(kmxm,lWC(l)))/2./hp(kms,l)
       vbott = (vhp(kmyp,l)+vhp(kmym,lSC(l)))/2./hp(kms,l)
       taubx = cd * (rhop(kms,l)+1000.) * ubott**2.
       tauby = cd * (rhop(kms,l)+1000.) * vbott**2.

       ! ... Add Wave-Induced Bottom Shear Stress & calculate friction velocity
       taubx = taubx + wibssx(i,j)
       tauby = tauby + wibssy(i,j)
       taub  = sqrt(taubx**2.+tauby**2.)
       ustar = sqrt(taub/1000.)

       ! ... Compute Erosion flux E (kg/m2/s) at bottom cell
       Zu5  = ( (ustar/vdep(itr)) * fRep )**5.;
       sourcesink(kms,l,itr)= Ased*(Zu5) / ( 1. + (Ased/0.3)*(Zu5))

     ENDDO ! End loop over wet columns

   ENDDO ! End loop over tracers

END SUBROUTINE srcsnkSD

!***********************************************************************
SUBROUTINE CalculateWIBSS
!***********************************************************************

   !.....Local variables.....
   INTEGER           :: i, j, k, l, ios, istat
   INTEGER, SAVE     :: nn = 0
   REAL              :: weight
   CHARACTER(LEN=12) :: wibssfmt, wibssfile

   ! ... Allocate space & initialize on first call
   IF ( nn .EQ. 0) THEN
     ALLOCATE ( wibssIOx(im1,jm1), wibssIOxp(im1,jm1),         &
              & wibssIOy(im1,jm1), wibssIOyp(im1,jm1),         &
              & wibssx  (im1,jm1), wibssy   (im1,jm1), STAT=istat)
     IF (istat /= 0) CALL allocate_error ( istat, 90 )
     wibssIOx = 0.0E0; wibssIOxp = 0.0E0
     wibssIOy = 0.0E0; wibssIOyp = 0.0E0
     wibssx   = 0.0E0; wibssx    = 0.0E0
   ENDIF

   ! Save and read new frame if thrs > thrsWIBSS.......
   IF (thrs > thrsWIBSS) THEN

     ! ... Update times
     nn = nn + 1
     thrsWIBSSp = thrsWIBSS
     thrsWIBSS  = thrsWIBSS + dthrsWIBSS

     ! ... Set initial time to thrsWIBSS0 for nn = 1
     IF (nn .EQ. 1) thrsWIBSS = thrsWIBSS0

     ! ... Update wibssIOxp & wibssIOyp
     wibssIOxp = wibssIOx
     wibssIOyp = wibssIOy

     ! ... Input file name
     wibssfile = "wibss00 .txt"
     IF (                nn <  10) WRITE ( wibssfile(8:8), FMT='(I1)' ) nn
     IF ( nn >= 10 .AND. nn < 100) WRITE ( wibssfile(7:8), FMT='(I2)' ) nn
     IF ( nn >= 100              ) WRITE ( wibssfile(6:8), FMT='(I3)' ) nn

     !.....Open wibss file.....
     OPEN (UNIT=i99, FILE=wibssfile, STATUS="OLD",  FORM="FORMATTED", IOSTAT=ios)
     IF (ios /= 0) CALL open_error ( "Error opening "//wibssfile, ios )

     ! ... Read wibssIOx & wibssIOy
     WRITE (UNIT=wibssfmt, FMT='("(", I4, "G9.0)")') im-i1+1
     DO j = jm, j1, -1
       READ (UNIT=i99, FMT=wibssfmt, IOSTAT=ios) (wibssIOx(i,j), i = i1, im)
       IF (ios /= 0) CALL input_error ( ios, 92 )
     END DO
     DO j = jm, j1, -1
       READ (UNIT=i99, FMT=wibssfmt, IOSTAT=ios) (wibssIOy(i,j), i = i1, im)
       IF (ios /= 0) CALL input_error ( ios, 92 )
     END DO
     CLOSE(i99)

   ENDIF

   ! ... Estimate wibss at present time step by interpolation
   weight  = (thrs - thrsWIBSSp)/(thrsWIBSS-thrsWIBSSp)
   DO l = 1, lm
     i = l2i(l);
	 j = l2j(l);
     wibssx (i,j) = wibssIOx (i,j)*(   weight) + &
                    wibssIOxp(i,j)*(1.-weight)
     wibssy (i,j) = wibssIOy (i,j)*(   weight) + &
                    wibssIOyp(i,j)*(1.-weight)
   ENDDO

END SUBROUTINE CalculateWIBSS

!***********************************************************************
SUBROUTINE SZinput
!***********************************************************************
!
!  Purpose: To read in values for diameter, max growth rate & settling
!           velocity for each size group (up to ntr)
!
!-----------------------------------------------------------------------

  !. . . .Local Variables
  CHARACTER(LEN=12) :: sziofile="si3d_sz.txt"
  CHARACTER(LEN=18) :: ioszfmt
  INTEGER:: ios, itr ,idep, istat
  INTEGER, PARAMETER :: niosize = 4

  !.....Open input parameter file.....
  OPEN (UNIT=i99, FILE=sziofile, STATUS="OLD",  FORM="FORMATTED", IOSTAT=ios)
  IF (ios /= 0) CALL open_error ( "Error opening "//sziofile, ios )

  ! ... Allocate space for arrays holding information on size groups
  ALLOCATE ( diamsz(ntr), rmaxsz(ntr), vdep(ntr), STAT=istat)
  IF (istat /= 0) CALL allocate_error ( istat, 99 )

  ! Skip over first five header records in SIZE file
  READ (UNIT=i99, FMT='(/////)', IOSTAT=ios)
  IF (ios /= 0) CALL input_error ( ios, 90 )

  ! Write the format of the data records into an internal file
  WRITE (UNIT=ioszfmt, FMT='("(10X,",I3,"G11.2)")') niosize

  ! Read information for each size group (diameter, growth rate & deposition)
  DO itr = 1, ntr
    READ (UNIT=i99, FMT=ioszfmt, IOSTAT=ios) &
         diamsz(itr), rmaxsz(itr), idep, vdep(itr)
    PRINT *, itr, diamsz(itr), rmaxsz(itr), idep, vdep(itr)
    IF (ios /= 0) CALL input_error ( ios, 91 )
    IF (idep > 0) THEN ! Stokes Law
      vdep(itr) = 10. ** (1.112 * LOG(diamsz(itr))/LOG(10.) - 1.647); ! Speed in m/day
      vdep(itr) = vdep(itr) / 86400. ! Speed in m/s
    ENDIF
  END DO

  CLOSE (i99)

END SUBROUTINE SZinput

!***********************************************************************
SUBROUTINE srcsnkSZ
!***********************************************************************
!
!  Purpose: To evaluate the source-sink terms in the reactive-transport
!           equations, used in the analysis of size structures
!
!-----------------------------------------------------------------------

   ! ... Local variables
   INTEGER :: i, j, k, l, k1s, kms, itr
   REAL    :: r, zt

   ! ... Loop over tracers
   DO itr = 1, ntr

     ! ... Loop over wet columns
     DO l = 1, lm

       ! ... 3D-(i,j) indexes for l
       i = l2i(l); j = l2j(l);

       !.....Compute top & bottom layer numbers & No. of layers ....
       kms = kmz(l);
       k1s = k1z(l);

       ! ... Growth rate (express it in s^-1)
       k = k1s;
       zt = hpp(k,l)/2.
       sourcesink(k,l,itr)= rmaxsz(itr)*tracerpp(k,l,itr)* &
                            EXP(-0.3*zt)*hpp(k,l) / 86400.
       DO k = k1s+1, kms
         zt = zt + (hpp(k-1,l)+hpp(k,l))/2.
         sourcesink(k,l,itr)= rmaxsz(itr)*tracerpp(k,l,itr)* &
                              EXP(-0.3*zt)*hpp(k,l) / 86400.
       ENDDO

     ENDDO ! End loop over wet columns

   ENDDO ! End loop over tracers

END SUBROUTINE srcsnkSZ

!************************************************************
SUBROUTINE WQinput
!***********************************************************
!
!            Purpose: To read wq_inp.txt
!
!------------------------------------------------------------


  !. . . .Local Variables
  CHARACTER(LEN=15):: wq_input_file="si3d_wq_inp.txt"
  INTEGER::   ios, nn

  !.....Open input parameter file.....
  OPEN (UNIT=i99, FILE=wq_input_file, STATUS="OLD", IOSTAT=ios)
  IF (ios /= 0) CALL open_error ( "Error opening "//wq_input_file, ios )

  !.....Read header record containing comments about run................
  READ (UNIT=i99, FMT='(/(A))', IOSTAT=ios) title
  IF (ios /= 0) CALL input_error ( ios, 91)

  !. . Read list of tracerpps modeled: Dissolved oxygen, N forms, P forms, C forms 

  READ (UNIT=i99,FMT='(///(18X,I20))',IOSTAT=ios) iDO,  &
      iPON, iDON, iNH4, iNO3, iPOP, iDOP, iPO4, iPOC,   &
      iDOC, iALG1, iALG2, iALG3, iALG4, iMeHg, iHg2,    &
      iHg0, iSS 
  IF (ios /= 0) CALL input_error ( ios, 92)

  !. . . Read model stochiometeric constants and other constants
  READ (UNIT=i99,FMT='(///(18X,G20.3))',IOSTAT=ios) rnc, rpc, roc, ron, &
  &     KSOD, KDECMIN, KSED, KNIT, KSN, KSP, FNH4,  &
  &     light_sat1, light_sat2, light_sat3, light_sat4
  IF (ios /= 0) CALL input_error ( ios, 93)

  !. . . Read model rates
  READ (UNIT=i99,FMT='(///(18X,G20.2))',IOSTAT=ios) R_reaer, R_SOD, &
  &    mu_max1, R_mor1, R_gr1, mu_max2, R_mor2, R_gr2, mu_max3, R_mor3, R_gr3, mu_max4, R_mor4, R_gr4, &
  &    R_decom_pon, R_miner_don, R_nitrif, R_denit, &
  &    R_decom_pop, R_miner_dop, R_decom_poc, R_miner_doc, &
  &    R_settl, R_resusp
  IF (ios /= 0) CALL input_error ( ios, 94)

  !. . . Read model temperature rates
  READ (UNIT=i99,FMT='(///(18X,G20.2))',IOSTAT=ios) Theta_SOD, Theta_mu, Theta_mor, Theta_gr, &
  &     Theta_decom, Theta_miner, Theta_sedflux, Theta_nitrif , Theta_denit 

  IF (ios /= 0) CALL input_error ( ios, 95)

  !. . . Read miscillaneous rates
  READ (UNIT=i99,FMT='(///(18X,G20.2))',IOSTAT=ios) ATM_DON, ATM_NH4,    &
  &    ATM_NO3, ATM_DOP, ATM_PO4,  ATM_DOC, SED_DON, SED_NH4, SED_NO3, &
  &    SED_DOP, SED_PO4, SED_DOC
  IF (ios /= 0) CALL input_error ( ios, 96)

  !. . . Read Sediment Parameters
  READ (UNIT=i99,FMT='(///(18X,I20))',IOSTAT=ios) sedNumber
  IF (ios /= 0) CALL input_error ( ios, 97 )

  IF (sedNumber > 0) THEN
    if (sedNumber .gt. sedMax) then
      print*,('************************************************')
      print*,('                  WARNING                       ')
      print*,('       Number of Sediment Particles greater     ')
      print*,('  Than the possible number (i.e., sedNumber> 3) ')
      print*,('                EXITING MODEL                   ')
      print*,('****************************************')
      STOP 
    end if
   ALLOCATE(sed_diameter(sedNumber),sed_dens(sedNumber),sed_frac(sedNumber))
   READ (UNIT=i99, FMT='(18X,5F)', IOSTAT=ios) (sed_diameter(nn), nn = 1, sedNumber)
   IF (ios /= 0) CALL input_error ( ios, 98 )
   READ (UNIT=i99, FMT='(18X,5F)', IOSTAT=ios) (sed_dens(nn), nn = 1, sedNumber)
   IF (ios /= 0) CALL input_error ( ios, 99 )
   READ (UNIT=i99, FMT='(18X,5F)', IOSTAT=ios) (sed_frac(nn), nn = 1, sedNumber)
   IF (ios /= 0) CALL input_error ( ios, 100 )
  ELSE IF (sedNumber == 0) THEN
   READ (UNIT=i5, FMT='(18X,5F)', IOSTAT=ios)
   READ (UNIT=i5, FMT='(18X,5F)', IOSTAT=ios)
   READ (UNIT=i5, FMT='(18X,5F)', IOSTAT=ios)
   IF (ios /= 0) CALL input_error ( ios, 101 )
  END IF

  !... Convert model rates [1/s]. Input file has 1/day values. WQ module is run every hour
  ! DO
  R_reaer   =  R_reaer/86400.0
  R_SOD     =  R_SOD/86400.0
  ! ALG 
  mu_max1 =  mu_max1/86400.0
  R_mor1 =  R_mor1/86400.0
  R_gr1 =  R_gr1/86400.0
  mu_max2 =  mu_max2/86400.0
  R_mor2 =  R_mor2/86400.0
  R_gr2 =  R_gr2/86400.0
  mu_max3 =  mu_max3/86400.0
  R_mor3 =  R_mor3/86400.0
  R_gr3 =  R_gr3/86400.0
  mu_max4 =  mu_max4/86400.0
  R_mor4 =  R_mor4/86400.0
  R_gr4 =  R_gr4/86400.0
  ! Nutrients
  R_decom_pon =  R_decom_pon/86400.0
  R_miner_don =  R_miner_don/86400.0
  R_nitrif =  R_nitrif/86400.0
  R_denit =  R_denit/86400.0
  R_decom_pop =  R_decom_pop/86400.0
  R_miner_dop =  R_miner_dop/86400.0
  R_decom_poc =  R_decom_poc/86400.0
  R_miner_doc =  R_miner_doc/86400.0
  R_settl =  R_settl/86400.0
  R_resusp =  R_resusp/86400.0

  ATM_DON = ATM_DON/86400.0
  ATM_NH4 = ATM_NH4/86400.0
  ATM_NO3 = ATM_NO3/86400.0
  ATM_DOP = ATM_DOP/86400.0
  ATM_PO4 = ATM_PO4/86400.0
  ATM_DOC = ATM_DOC/86400.0
  SED_DON = SED_DON/86400.0
  SED_NH4 = SED_NH4/86400.0
  SED_NO3 = SED_NO3/86400.0
  SED_DOP = SED_DOP/86400.0
  SED_PO4 = SED_PO4/86400.0
  SED_DOC = SED_DOC/86400.0

  IF (idbg == 1) THEN
    PRINT*, "iDO  = ", iDO , "iPOC = ", iPOC, "iDOC = ", iDOC
    PRINT*, "iPON = ", iPON, "iDON = ", iDON, "iNH4 = ", iNH4, "iNO3 = ", iNO3
    PRINT*, "iPOP = ", iPOP, "iDOP = ", iDOP, "iPO4 = ", iPO4
    PRINT*, "iALG1 = ", iALG1, "iALG2 = ", iALG2, "iALG3 = ", iALG3, "iALG4 = ", iALG4
    PRINT*, 'iMeHg = ', iMeHg, 'iHg(II) = ',iHg2, 'iHg(0) = ', iHg0, 'iSS = ', iSS
    PRINT*, 'sed_diam',sed_diameter,'sed_dens',sed_dens 
  END IF

  CALL WQinit !ACortes 09/24/2021

END SUBROUTINE WQinput

!************************************************************************
SUBROUTINE WQinit
!************************************************************************
!
! Purpose: initializes water quality variables
!
!-----------------------------------------------------------------------

  !. . . Local Variables
  INTEGER, DIMENSION(ntrmax):: tracerpplocal
  INTEGER:: i,j, sumtr

  ! ... Initialize tracerpp local
  tracerpplocal = 0

  !. . Initialize constituent locations
  LDO =0;
  LPON=0; LDON=0; LNH4=0; LNO3=0;
  LPOP=0; LDOP=0; LPO4=0
  LALG1=0; LALG2=0; LALG3=0; LALG4=0;
  LDOC=0; LPOC=0;
  LSS1 = 0 
  LSS2 = 0
  LSS3 = 0
  
  !. . Assign Lxx to each constituent modeled
  !. . .. first need to define intermediate array tracerpplocal
  i=1

  IF (iDO == 1) THEN
    tracerpplocal(i) = 1
    i = i+1
  END IF

  IF (iPON == 1) THEN
    tracerpplocal(i) = 2
    i = i+1
  END IF

  IF (iDON == 1) THEN
    tracerpplocal(i) = 3
    i = i+1
  END IF

  IF (iNH4 == 1) THEN
    tracerpplocal(i) = 4
    i = i+1
  END IF

  IF (iNO3 == 1) THEN
    tracerpplocal(i) = 5
    i = i+1
  END IF

  IF (iPOP == 1) THEN
    tracerpplocal(i) = 6
    i = i+1
  END IF

  IF (iDOP == 1) THEN
    tracerpplocal(i) = 7
    i = i+1
  END IF

  IF (iPO4 == 1) THEN
    tracerpplocal(i) = 8
    i = i+1
  END IF

  IF (iPOC == 1) THEN
    tracerpplocal(i) = 9
    i = i+1
  END IF

  IF (iDOC == 1) THEN
    tracerpplocal(i) = 10
    i = i+1
  END IF

  IF (iALG1 == 1) THEN
    tracerpplocal(i) = 11
    i = i+1
  END IF

    IF (iALG2 == 1) THEN
    tracerpplocal(i) = 12
    i = i+1
  END IF
  
  IF (iALG3 == 1) THEN
    tracerpplocal(i) = 13
    i = i+1
  END IF

  IF (iALG4 == 1) THEN
    tracerpplocal(i) = 14
    i = i+1
  END IF

  IF (iMeHg == 1) THEN
    tracerpplocal(i) = 15
    i = i+1
  END IF

  IF (iHg2 == 1) THEN
    tracerpplocal(i) = 16
    i = i+1
  END IF

  IF (iHg0 == 1) THEN
    tracerpplocal(i) = 17
    i = i+1
  END IF

  IF (iSS == 1) THEN
    do j = 1, sedNumber
      tracerpplocal(i) = 18 + (j-1)
      i = i+1
    end do 
  END IF

  PRINT*, 'i',i
  PRINT*, 'tracerpplocal', tracerpplocal

  sumtr = 0
  DO j = 1, ntrmax
    print*,'j',j
    IF (tracerpplocal(j)>0) sumtr = sumtr + 1
  END DO

  IF (idbg == 1) THEN
    PRINT*, "ntr = ", ntr
    PRINT*, "sumtr = ", sumtr
  END IF

  IF (sumtr .ne. ntr) THEN
    print*,('****************************************')
    print*,('               WARNING                  ')
    print*,('    NTR AND NO. tracerpp DO NOT MATCH     ')
    print*,('           EXITING MODEL                ')
    print*,('****************************************')
    STOP
  END IF

  !. . Next define Lxx
  DO i = 1, ntr
  IF (tracerpplocal(i) == 1) THEN
    LDO  = i
  ELSEIF (tracerpplocal(i) == 2) THEN
    LPON = i
  ELSEIF (tracerpplocal(i) == 3) THEN
    LDON = i
  ELSEIF (tracerpplocal(i) == 4) THEN
    LNH4 = i
  ELSEIF (tracerpplocal(i) == 5) THEN
    LNO3 = i
  ELSEIF (tracerpplocal(i) == 6) THEN
    LPOP = i
  ELSEIF (tracerpplocal(i) == 7) THEN
    LDOP = i
  ELSEIF (tracerpplocal(i) == 8) THEN
    LPO4 = i
  ELSEIF (tracerpplocal(i) == 9) THEN
    LPOC = i
  ELSEIF (tracerpplocal(i) == 10) THEN
    LDOC = i
  ELSEIF (tracerpplocal(i) == 11) THEN
    LALG1 = i
  ELSEIF (tracerpplocal(i) == 12) THEN
    LALG2 = i   
  ELSEIF (tracerpplocal(i) == 13) THEN
    LALG3 = i 
  ELSEIF (tracerpplocal(i) == 14) THEN
    LALG4 = i 
  ELSEIF (tracerpplocal(i) == 15) THEN
    LMeHg = i
  ELSEIF (tracerpplocal(i) == 16) THEN
    LHg2 = i
  ELSEIF (tracerpplocal(i) == 17) THEN
    LHg0 = i
  ELSEIF (tracerpplocal(i) == 18) THEN
    LSS1 = i
  ELSEIF (tracerpplocal(i) == 19) THEN
    LSS2 = i
  ELSEIF (tracerpplocal(i) == 20) THEN
    LSS3 = i
  END IF
  END DO

  IF (idbg == 1) THEN
    PRINT*, "LDO  = ", LDO, "LPOC = ", LPOC, "LDOC = ", LDOC
    PRINT*, "LPON = ", LPON, "LDON = ", LDON, "LNH4 = ", LNH4, "LNO3 = ", LNO3
    PRINT*, "LPOP = ", LPOP, "LDOP = ", LDOP, "LPO4 = ", LPO4
    PRINT*, "LALG1 = ", LALG1, "LALG2 = ", LALG2, "LALG3 = ", LALG3, "LALG4 = ", LALG4
    PRINT*, "LMeHg = ", LMeHg, "LHg2 = ", LHg2, "LHg0 = ", LHg0
    PRINT*, "LSS1 = ", LSS1, 'LSS2',LSS2,'LSS3 = ', LSS3
  END IF

END SUBROUTINE WQinit

!************************************************************
SUBROUTINE srcsnkWQ(n)
!***********************************************************
!
!   Purpose: to call all source subroutines for all cells
!
!------------------------------------------------------------

  !... Local variables
  INTEGER:: i, j, k, l, liter, k1s, kms, iteration
  integer, intent(in) :: n 
  REAL :: thrs1

  ! reset soursesink = 0
  sourcesink = 0.0
  thrs1 = n * dt / 3600

  ! STWAVE controlling section. (SergioValbuena 03-11-2023) 
  if ((iSS == 1) .and. (iSTWAVE == 1)) then
    call stwave_input(n)
    print*,'tau'
    print*,transpose(tau_stwave(2:im1-1,2:jm1-1))
  end if 

  !$omp barrier  

  DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
    l = id_column(liter)
    ! ... Retrieve top & bottom wet sal-pts .................
    kms = kmz(l)
    k1s = k1z(l)

    DO k = k1s, kms;

      IF (iDO == 1) THEN
        CALL sourceDO(k,l,thrs1)
      END IF
      IF (iPON == 1) THEN
        CALL sourcePON(k,l,thrs1)
      END IF
      IF (iDON == 1) THEN
        CALL sourceDON(k,l,thrs1)
      END IF
      IF (iNH4 == 1) THEN
        CALL sourceNH4(k,l,thrs1)
      END IF
      IF (iNO3 == 1) THEN
        CALL sourceNO3(k,l,thrs1)
      END IF
      IF (iPOP == 1) THEN
        CALL sourcePOP(k,l,thrs1)
      END IF
      IF (iDOP == 1) THEN
        CALL sourceDOP(k,l,thrs1)
      END IF
      IF (iPO4 == 1) THEN
        CALL sourcePO4(k,l,thrs1)
      END IF
      IF (iDOC == 1) THEN
        CALL sourceDOC(k,l,thrs1)
      END IF
      IF (iPOC == 1) THEN
        CALL sourcePOC(k,l,thrs1)
      END IF
      IF (iALG1 == 1) THEN
        CALL sourceALG1(k,l,thrs1)
      END IF
      IF (iALG2 == 1) THEN
        CALL sourceALG2(k,l,thrs1)
      END IF
      IF (iALG3 == 1) THEN
        CALL sourceALG3(k,l,thrs1)
      END IF
      IF (iALG4 == 1) THEN
        CALL sourceALG4(k,l,thrs1)
      END IF
      IF (iMeHg == 1) THEN
        ! CALL sourceALG4(k,l,thrs1)
      END IF
      IF (iHg2 == 1) THEN
        ! CALL sourceALG4(k,l,thrs1)
      END IF
      IF (iHg0 == 1) THEN
        ! CALL sourceALG4(k,l,thrs1)
      END IF
      IF (iSS == 1) THEN
        call sourceSS(k, l)
      END IF
    END DO
  END DO

END SUBROUTINE srcsnkWQ

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

!***********************************************************************
FUNCTION parabn ( frstpt, x, fx, dx )
!***********************************************************************
!
!  Purpose: To interpolate parabolically between the functional values
!           within the array fx. May not be used
!
!-----------------------------------------------------------------------

   REAL, DIMENSION(:), INTENT(IN) :: fx      ! Assumed-shape array
   REAL, INTENT(IN) :: frstpt, x, dx
   REAL :: parabn
   REAL :: om, theta
   INTEGER :: m

   m = (x - frstpt)/dx
   om = m
   theta = (x - frstpt - om*dx) / dx
   IF (m == 0) THEN
      m = 2
      theta = theta - 1.0
   ELSE
      m = m + 1
   END IF
   parabn=fx(m)+0.5*theta*(fx(m+1)-fx(m-1)+theta*(fx(m+1)+fx(m-1)-2.0*fx(m)))

END FUNCTION parabn

!************************************************************************
SUBROUTINE allocate_error ( istat, ierror_code )
!************************************************************************
!
!  Purpose: Prints messages regarding any errors encountered during
!           allocation of space for model arrays. Stops program after
!           messages.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!------------------------------------------------------------------------

   !.....Arguments.....
   INTEGER, INTENT(IN) :: istat, ierror_code

   !.....Print error messages....
   PRINT *, " Program could not allocate space for arrays"
   PRINT '(" istat= ", I5, "  error code=", I3)', istat, ierror_code
   PRINT *, "  "
   PRINT *, "  "
   PRINT *, " ****STOPPING si3d due to allocate error"
   STOP

END SUBROUTINE allocate_error


!************************************************************************
SUBROUTINE open_error ( string, iostat )
!************************************************************************
!
!  Purpose: Prints 'string' regarding any error encountered during the
!           opening of a file. Also prints the iostat error code. The
!           program is stopped after printing messages.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!------------------------------------------------------------------------

   !.....Arguments.....
   CHARACTER(LEN=*), INTENT(IN) :: string
   INTEGER, INTENT(IN) :: iostat

   !.....Print error messages....
   PRINT '(A)', string
   PRINT '(" The iostat error number is ", I5)', iostat
   PRINT *, "  "
   PRINT *, "  "
   PRINT *, " ****STOPPING si3d due to OPEN error"
   STOP

END SUBROUTINE open_error


!************************************************************************
SUBROUTINE input_error ( ios, ierror_code )
!************************************************************************
!
!  Purpose: Prints messages regarding any errors encountered during
!           reading of the input file. Stops program after messages.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!------------------------------------------------------------------------

   !.....Arguments.....
   INTEGER, INTENT(IN) :: ios, ierror_code

   !.....Determine type of read error.....
   SELECT CASE (ios)

   !.....End of file (ios=-1).....
   CASE (:-1)
      PRINT *, " Unexpected end-of-file encountered while reading file"

   !.....Error during reading (ios>0).....
   CASE (1:)
      PRINT *, " Error during reading data"
   END SELECT

   !.....Print ios and error statement number.....
   PRINT '(" The iostat error number is ", I5)',     ios
   PRINT '(" Error on read statement number ", I3)', ierror_code
   PRINT *, "  "
   PRINT *, "  "
   PRINT *, " ****STOPPING si3d due to read error"
   STOP

END SUBROUTINE input_error

END MODULE si3d_ecomod