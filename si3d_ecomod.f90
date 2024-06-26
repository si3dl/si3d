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
   USE si3d_wq
   USE si3d_sed
   USE si3d_stwave
   USE si3d_Hg

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
  READ (UNIT=i99, FMT='(/(A))', IOSTAT=ios)
  IF (ios /= 0) CALL input_error ( ios, 91)

  !. . Read list of tracerpps modeled: Dissolved oxygen, N forms, P forms, C forms 

  READ (UNIT=i99,FMT='(///(18X,I20))',IOSTAT=ios) iDO,  &
      iPON, iDON, iNH4, iNO3, iPOP, iDOP, iPO4, iPOC,   &
      iDOC, iALG1, iALG2, iALG3, iALG4, iALG5, iMeHg, iHgII,    &
      iHg0, iSS 
  IF (ios /= 0) CALL input_error ( ios, 92)

  !. . . Read model stochiometeric constants and other constants
  READ (UNIT=i99,FMT='(///(18X,G20.3))',IOSTAT=ios) rnc, rpc, roc, ron, &
  &     KSOD, KDECMIN, KSED, KNIT, KSN, KSP, FNH4,  &
  &     light_sat1, light_sat2, light_sat3, light_sat4, light_sat5
  IF (ios /= 0) CALL input_error ( ios, 93)

  !. . . Read model rates
  READ (UNIT=i99,FMT='(///(18X,G20.2))',IOSTAT=ios)  mu_max1, R_mor1, R_gr1,  &
  &    mu_max2, R_mor2, R_gr2, mu_max3, R_mor3, R_gr3, mu_max4, R_mor4, R_gr4, mu_max5, R_mor5, R_gr5, &
  &    R_decom_pon, R_miner_don, R_nitrif, R_denit, &
  &    R_decom_pop, R_miner_dop, R_decom_poc, R_miner_doc, &
  &    R_reaer, R_settl, R_resusp, vspa, vspoc
  IF (ios /= 0) CALL input_error ( ios, 94)

  !. . . Read model temperature correction factors [-]
  READ (UNIT=i99,FMT='(///(18X,G20.2))',IOSTAT=ios) Topt1, Topt2, Topt3, Topt4, Topt5, &
  &     Theta_SOD, Theta_mor, Theta_gr, &
  &     Theta_decom, Theta_miner, Theta_sedflux, Theta_nitrif , Theta_denit

  IF (ios /= 0) CALL input_error ( ios, 95)

  !. . . Read miscillaneous fluxes in mg/m2/d
  READ (UNIT=i99,FMT='(///(18X,G20.2))',IOSTAT=ios) R_SOD, ATM_DON, ATM_NH4,    &
  &    ATM_NO3, ATM_DOP, ATM_PO4,  ATM_DOC, SED_DON, SED_NH4, SED_NO3, &
  &    SED_DOP, SED_PO4, SED_DOC
  IF (ios /= 0) CALL input_error ( ios, 96)

  !. . . Read Sediment Parameters
  READ (UNIT=i99,FMT='(///(18X,I20))',IOSTAT=ios) sedNumber
  IF (ios /= 0) CALL input_error ( ios, 97 )

  IF ((iSS == 1) .AND. (sedNumber > 0)) THEN
    if (sedNumber .gt. sedMax) then
      print*,('************************************************')
      print*,('                  WARNING                       ')
      print*,('       Number of Sediment Particles greater     ')
      print*,('  Than the possible number (i.e., sedNumber> 3) ')
      print*,('                EXITING MODEL                   ')
      print*,('************************************************')
      STOP
    end if
    ALLOCATE(sed_diameter(sedNumber), sed_dens(sedNumber), sed_frac(sedNumber), sed_type(sedNumber))
    READ (UNIT=i99,FMT='(18X,I20)',IOSTAT=ios) iSTWAVE
    IF (ios /= 0) CALL input_error ( ios, 98)
    READ (UNIT=i99,FMT='(18X,G20.2)',IOSTAT=ios) sed_h
    IF (ios /= 0) CALL input_error ( ios, 99)
    READ (UNIT=i99, FMT='(18X,5F)', IOSTAT=ios) (sed_diameter(nn), nn = 1, sedNumber)
    IF (ios /= 0) CALL input_error ( ios, 100)
    READ (UNIT=i99, FMT='(18X,5F)', IOSTAT=ios) (sed_dens(nn), nn = 1, sedNumber)
    IF (ios /= 0) CALL input_error ( ios, 101)
    READ (UNIT=i99, FMT='(18X,5F)', IOSTAT=ios) (sed_frac(nn), nn = 1, sedNumber)
    IF (ios /= 0) CALL input_error ( ios, 102)
    READ (UNIT=i99,FMT='(18X,5I)', IOSTAT=ios) (sed_type(nn), nn = 1, sedNumber)
    IF (ios /= 0) CALL input_error ( ios, 103)
  ELSE IF (sedNumber == 0) THEN
    READ (UNIT=i99, FMT='(18X,I20)', IOSTAT=ios)
    READ (UNIT=i99, FMT='(18X,G20.2)', IOSTAT=ios) sed_h
    READ (UNIT=i99, FMT='(18X,5F)', IOSTAT=ios)
    READ (UNIT=i99, FMT='(18X,5F)', IOSTAT=ios)
    READ (UNIT=i99, FMT='(18X,5F)', IOSTAT=ios)
    READ (UNIT=i99, FMT='(18X,5I)', IOSTAT=ios)
    IF (ios /= 0) CALL input_error ( ios, 104)
  END IF

  if (iMeHg == 1) then
    allocate(kd_wpn3(sedNumber), kd_spn3(sedNumber))
    READ (UNIT=i99, FMT='(///(18X,G20.5))', IOSTAT=ios) inst_eq
    IF (ios /= 0) CALL input_error ( ios, 105)
    READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) kw31
    IF (ios /= 0) CALL input_error ( ios, 106)
    READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) atm_MeHg
    IF (ios /= 0) CALL input_error ( ios, 107)
    READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) k_MeHgw
    IF (ios /= 0) CALL input_error ( ios, 108)
    READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) k_MeHgatm
    IF (ios /= 0) CALL input_error ( ios, 109)
    READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) MeHgatm
    IF (ios /= 0) CALL input_error ( ios, 110)
    READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) K_H_MeHgw
    IF (ios /= 0) CALL input_error ( ios, 111)
    READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) kw32
    IF (ios /= 0) CALL input_error ( ios, 112)
    READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) ks32
    IF (ios /= 0) CALL input_error ( ios, 113)
    READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) kws
    IF (ios /= 0) CALL input_error ( ios, 114)

    ! if (iDOC == 1) then
      READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) kd_wdoc3
      IF (ios /= 0) CALL input_error ( ios, 115)
      READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) kd_sdoc3
      IF (ios /= 0) CALL input_error ( ios, 115)
    ! else
      ! READ (UNIT=i99, FMT='(18X,G20.5)', IOSTAT=ios)
      ! IF (ios /= 0) CALL input_error ( ios, 115)
      ! READ (UNIT=i99, FMT='(18X,G20.5)', IOSTAT=ios)
      ! IF (ios /= 0) CALL input_error ( ios, 115)
    ! end if

    ! if (iALG1 == 1) then
      READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) kd_wpa3
      IF (ios /= 0) CALL input_error ( ios, 116)
    ! else
      ! READ (UNIT=i99, FMT='(18X,G20.5)', IOSTAT=ios)
      ! IF (ios /= 0) CALL input_error ( ios, 116)
    ! end if

    ! if (iPOC == 1) then
      READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) kd_wpom3
      IF (ios /= 0) CALL input_error ( ios, 117)
      READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) kd_spom3
      IF (ios /= 0) CALL input_error ( ios, 117)
    ! else
      ! READ (UNIT=i99, FMT='(18X,G20.5)', IOSTAT=ios)
      ! IF (ios /= 0) CALL input_error ( ios, 117)
      ! READ (UNIT=i99, FMT='(18X,G20.5)', IOSTAT=ios)
      ! IF (ios /= 0) CALL input_error ( ios, 117)
    ! end if

    ! if (iSS == 1) then
      READ (UNIT=i99, FMT='(18X,5F)', IOSTAT=ios) (kd_wpn3(nn), nn = 1, sedNumber)
      IF (ios /= 0) CALL input_error ( ios, 119)
      READ (UNIT=i99, FMT='(18X,5F)', IOSTAT=ios) (kd_spn3(nn), nn = 1, sedNumber)
      IF (ios /= 0) CALL input_error ( ios, 120)
    ! else
    !   READ (UNIT=i99, FMT='(18X,5F)', IOSTAT=ios)
    !   IF (ios /= 0) CALL input_error ( ios, 119)
    !   READ (UNIT=i99, FMT='(18X,5F)', IOSTAT=ios)
    !   IF (ios /= 0) CALL input_error ( ios, 120)
    ! end if

    READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) MeHg_sed
    IF (ios /= 0) CALL input_error ( ios, 121)

  else if (iMeHg == 0) then
    READ (UNIT=i99, FMT='(////////////////////(18X,G20.5))', IOSTAT=ios)
    IF (ios /= 0) CALL input_error ( ios, 122)
  end if

  if (iHgII == 1) then
    allocate(kd_wpn2(sedNumber), kd_spn2(sedNumber))
    READ (UNIT=i99, FMT='(///(18X,G20.5))', IOSTAT=ios) kw21
    IF (ios /= 0) CALL input_error ( ios, 123)
    READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) atm_HgII
    IF (ios /= 0) CALL input_error ( ios, 124)
    READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) kw23
    IF (ios /= 0) CALL input_error ( ios, 125)
    READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) ks23
    IF (ios /= 0) CALL input_error ( ios, 126)
    READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) KDO
    IF (ios /= 0) CALL input_error ( ios, 127)
    READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) KSO4
    IF (ios /= 0) CALL input_error ( ios, 128)
    READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) SO4
    IF (ios /= 0) CALL input_error ( ios, 129)
    READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) miu_so4
    IF (ios /= 0) CALL input_error ( ios, 130)

    ! if (iDOC == 1) then
      READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) kd_wdoc2
      IF (ios /= 0) CALL input_error ( ios, 131)
      READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) kd_sdoc2
      IF (ios /= 0) CALL input_error ( ios, 131)
    ! else
      ! READ (UNIT=i99, FMT='(18X,G20.5)', IOSTAT=ios)
      ! IF (ios /= 0) CALL input_error ( ios, 131)
      ! READ (UNIT=i99, FMT='(18X,G20.5)', IOSTAT=ios)
      ! IF (ios /= 0) CALL input_error ( ios, 131)
    ! end if

    ! if (iALG1 == 1) then
      READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) kd_wpa2
      IF (ios /= 0) CALL input_error ( ios, 132)
    ! else
      ! READ (UNIT=i99, FMT='(18X,G20.5)', IOSTAT=ios)
      ! IF (ios /= 0) CALL input_error ( ios, 132)
    ! end if

    ! if (iPOC == 1) then
      READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) kd_wpom2
      IF (ios /= 0) CALL input_error ( ios, 133)
      READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) kd_spom2
      IF (ios /= 0) CALL input_error ( ios, 133)
    ! else
      ! READ (UNIT=i99, FMT='(18X,G20.5)', IOSTAT=ios)
      ! IF (ios /= 0) CALL input_error ( ios, 133)
      ! READ (UNIT=i99, FMT='(18X,G20.5)', IOSTAT=ios)
      ! IF (ios /= 0) CALL input_error ( ios, 133)
    ! end if

    ! if (iSS == 1) then
      READ (UNIT=i99, FMT='(18X,5F)', IOSTAT=ios) (kd_wpn2(nn), nn = 1, sedNumber)
      IF (ios /= 0) CALL input_error ( ios, 134)
      READ (UNIT=i99, FMT='(18X,5F)', IOSTAT=ios) (kd_spn2(nn), nn = 1, sedNumber)
      IF (ios /= 0) CALL input_error ( ios, 134)
    ! else
    !   READ (UNIT=i99, FMT='(18X,5F)', IOSTAT=ios)
    !   IF (ios /= 0) CALL input_error ( ios, 134)
    !   READ (UNIT=i99, FMT='(18X,5F)', IOSTAT=ios)
    !   IF (ios /= 0) CALL input_error ( ios, 134)
    ! end if
    READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) HgII_sed
    IF (ios /= 0) CALL input_error ( ios, 135)
  else if (iHgII == 0) then
    READ (UNIT=i99, FMT='(//////////////////(18X,G20.5))', IOSTAT=ios)
  end if

  if (iHg0 == 1) then
    READ (UNIT=i99, FMT='(///(18X,G20.5))', IOSTAT=ios) DGMra
    READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) k_Hg0w
    READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) k_Hg0atm
    READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) Hg0atm
    READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) K_H_Hg0w
    READ (UNIT=i99, FMT='((18X,G20.5))', IOSTAT=ios) Hg0_sed
  else if (iHg0 == 0) then
    READ (UNIT=i99, FMT='(////////(18X,G20.5))', IOSTAT=ios)
  end if

  !.....Close wq input file.....
   CLOSE (UNIT=i99)

  !... Convert model rates [1/s] and [m/s]. Input file has 1/day values. WQ module is run every hour
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
  mu_max5 =  mu_max5/86400.0
  R_mor5 =  R_mor5/86400.0
  R_gr5 =  R_gr5/86400.0

  ! Nutrients
  R_decom_pon =  R_decom_pon/86400.0
  R_miner_don =  R_miner_don/86400.0
  R_nitrif =  R_nitrif/86400.0
  R_denit =  R_denit/86400.0
  R_decom_pop =  R_decom_pop/86400.0
  R_miner_dop =  R_miner_dop/86400.0
  R_decom_poc =  R_decom_poc/86400.0
  R_miner_doc =  R_miner_doc/86400.0
  R_settl =  R_settl/86400.0          ! [m/s]
  R_resusp =  R_resusp/86400.0        ! [m/s]

  ! DO
  R_reaer   =  R_reaer/86400.0        ! [m/s]

  R_SOD   = R_SOD/86400.0
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


  ! Tranformation to seconds (/86400) and m3 (*1000) of input constants for Hg Model
   
  kw31 = kw31 / 86400.0
  atm_MeHg = atm_MeHg / 86400.0
  k_MeHgw = k_MeHgw / 86400.0
  k_MeHgatm = k_MeHgatm / 86400.0
  ! print*, 'K_H_MeHgw =', K_H_MeHgw
  ! print*, 'kw32 =', kw32
  ks32 = ks32 / 86400.0
  kws = kws / 86400.0
  kd_wdoc3 = kd_wdoc3 / 1000
  kd_sdoc3 = kd_sdoc3 / 1000
  kd_wpa3 =  kd_wpa3 / 1000
  kd_wpom3 = kd_wpom3 / 1000
  kd_spom3 = kd_spom3 / 1000
  kd_wpn3 = kd_wpn3 / 1000
  kd_spn3 = kd_spn3 / 1000

  kw21 = kw21 / 86400.0
  atm_HgII = atm_HgII / 86400.0
  ! print*, 'kw23 =', kw23
  ks23 = ks23 / 86400.0
  ! print*, 'KDO =', KDO
  ! print*, 'KSO4 =', KSO4
  ! print*, 'SO4 =', SO4
  ! print*, 'miu_so4 =', miu_so4
  kd_wdoc2 = kd_wdoc2 / 1000 !to go from L/ug to L/ng
  kd_sdoc2 = kd_sdoc2 / 1000
  kd_wpa2 = kd_wpa2 / 1000
  kd_wpom2 = kd_wpom2 / 1000
  kd_spom2 = kd_spom2 / 1000
  kd_wpn2 = kd_wpn2 / 1000
  kd_spn2 = kd_spn2 / 1000

  ! print*, 'DGMra =', DGMra
  k_Hg0w = k_Hg0w / 86400.0
  k_Hg0atm = k_Hg0atm / 86400.0
  ! print*, 'K_H_Hg0w =', K_H_Hg0w

  IF (idbg == 1) THEN
    PRINT*, 'iDO  = ', iDO , 'iPOC = ', iPOC, 'iDOC = ', iDOC
    PRINT*, 'iPON = ', iPON, 'iDON = ', iDON, 'iNH4 = ', iNH4
    PRINT*, 'iPOP = ', iPOP, 'iDOP = ', iDOP, 'iPO4 = ', iPO4
    PRINT*, 'iNO3 = ', iNO3, 'iALG1 = ', iALG1, 'iALG2 = ', iALG2
    PRINT*, 'iALG3 = ', iALG3, 'iALG4 = ', iALG4, 'iALG5 = ', iALG5
    PRINT*, 'iMeHg = ', iMeHg, 'iHgII = ',iHgII, 'iHg0 = ', iHg0, 'iSS = ', iSS
    print*, 'mu_max1 =', mu_max1, 'R_mor1 =', R_mor1, 'R_gr1 =', R_gr1
    print*, 'mu_max2 =', mu_max2, 'R_mor2 = ', R_mor2, 'R_gr2 =', R_gr2
    print*, 'mu_max3 =', mu_max3, 'R_mor3 = ', R_mor3, 'R_gr3 =', R_gr3
    print*, 'mu_max4 =', mu_max4, 'R_mor4 = ', R_mor4, 'R_gr4 =', R_gr4
    print*, 'mu_max5 =', mu_max5, 'R_mor5 = ', R_mor5, 'R_gr5 =', R_gr5
    print*, 'R_decom_pon =' , R_decom_pon, 'R_miner_don = ', R_miner_don 
    print*, 'R_nitrif = ', R_nitrif, 'R_denit = ', R_denit
    print*, 'R_decom_pop =' , R_decom_pop, 'R_miner_dop = ', R_miner_dop
    print*, 'R_decom_poc =' , R_decom_poc, 'R_miner_doc = ', R_miner_doc
    print*, 'R_settl = ', R_settl, 'R_resusp = ', R_resusp
    print*, 'R_reaer = ', R_reaer, 'R_SOD = ', R_SOD
    print*, 'ATM_DOC = ', ATM_DOC, 'ATM_DON = ', ATM_DON, 'ATM_DOP = ', ATM_DOP
    print*, 'ATM_NH4 = ', ATM_NH4, 'ATM_NO3 = ', ATM_NO3, 'ATM_PO4 = ', ATM_NO3
    print*, 'SED_DOC = ', SED_DOC, 'SED_DON = ', SED_DON, 'SED_DOP = ', SED_DOP
    print*, 'SED_NH4 = ', SED_NH4, 'SED_NO3 = ', SED_NO3, 'SED_PO4 = ', SED_PO4
    PRINT*, 'sed_h = ', sed_h
    PRINT*, 'sed_diam = ',sed_diameter
    PRINT*,'sed_dens = ',sed_dens
    PRINT*,'sed_frac = ',sed_frac
    PRINT*, 'sed_type = ',sed_type
    print*, 'inst_eq =',inst_eq
    
    print*, '----- Settings MeHg -----'
    print*, 'kw31 =',kw31
    print*, 'atm_MeHg =', atm_MeHg
    print*, 'k_MeHgw =', k_MeHgw
    print*, 'k_MeHgatm =', k_MeHgatm
    print*, 'MeHgatm =', MeHgatm
    print*, 'K_H_MeHgw =', K_H_MeHgw
    print*, 'kw32 =', kw32
    print*, 'ks32 =', ks32
    print*, 'kws =', kws
    print*, 'kd_wdoc3 =', kd_wdoc3
    print*, 'kd_sdoc3 =', kd_sdoc3
    print*, 'kd_wpa3 =', kd_wpa3
    print*, 'kd_wpom3 =', kd_wpom3
    print*, 'kd_spom3 =', kd_spom3
    print*, 'kd_wpn3 =', kd_wpn3
    print*, 'kd_spn3 =', kd_spn3

    print*, '----- Settings HgII -----'
    print*, 'kw21 =',kw21
    print*, 'atm_HgII =', atm_HgII
    print*, 'kw23 =', kw23
    print*, 'ks23 =', ks23
    print*, 'KDO =', KDO
    print*, 'KSO4 =', KSO4
    print*, 'SO4 =', SO4
    print*, 'miu_so4 =', miu_so4
    print*, 'kd_wdoc2 =', kd_wdoc2
    print*, 'kd_sdoc2 =', kd_sdoc2
    print*, 'kd_wpa2 =', kd_wpa2
    print*, 'kd_wpom2 =', kd_wpom2
    print*, 'kd_spom2 =', kd_spom2
    print*, 'kd_wpn2 =', kd_wpn2
    print*, 'kd_spn2 =', kd_spn2

    print*, '----- Settings Hg0 -----'
    print*, 'DGMra =', DGMra
    print*, 'k_Hg0w =', k_Hg0w
    print*, 'k_Hg0atm =', k_Hg0atm
    print*, 'Hg0atm =', Hg0atm
    print*, 'K_H_Hg0w =', K_H_Hg0w

  END IF

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
  INTEGER:: i,j, ios, ver
  character(11), allocatable, dimension(:) :: tracer_list
  

  allocate ( tracer_list (ntr), stat= ios)
  if (ios /= 0) then
    print*, 'Error alloc. init. arrays'
    stop
  end if

  ! ... Initialize tracerpp local
  tracerpplocal = 0

  !. . Initialize constituent locations
  LDO =0; LPON=0; LDON=0; LNH4=0; LNO3=0; LPOP=0; LDOP=0; LPO4=0
  LALG1=0; LALG2=0; LALG3=0; LALG4=0; LALG5=0; LDOC=0; LPOC=0;
  LSS1 = 0; LSS2 = 0; LSS3 = 0; LHg0 = 0; LMeHg = 0; LHgII = 0

  ! Read tracer line from si3d_init file 
  !.....Open initial condition file.....
  sal_ic_file = 'si3d_init.txt'
  OPEN (UNIT=i4, FILE='si3d_init.txt', STATUS="OLD", FORM="FORMATTED", IOSTAT=ios)
  IF(ios /= 0) CALL open_error ( "Error opening "//sal_ic_file, ios )

  ! Skip over first five header records in open boundary condition file
  READ (UNIT=i4, FMT='(////)', IOSTAT=ios)
  IF (ios /= 0) CALL input_error ( ios, 13 )

  READ (UNIT=i4, FMT='(A)', IOSTAT=ios) tracer_line
  CLOSE(i4)
  
  !. . Assign Lxx to each constituent modeled
  !. . .. first need to define intermediate array tracerpplocal

  tracer_line = tracer_line(22:)
  do i = 1, size(tracer_list)
      read(tracer_line, '(A11, A)', iostat=ios) tracer_list(i), tracer_line
      tracer_list(i) = adjustl(tracer_list(i))
  end do

  if (iSS == 1) then
    allocate(settling_vel(sedNumber), erosion_Hgpn(sedNumber))
    settling_vel(:) = 0.0
    erosion_Hgpn(:) = 0.0
  end if

  print*, 'Constituents to model:'
  print*, tracer_list(:)

  ver = iDO + iPON + iDON + iNH4 + iNO3 + iPOP + iDOP + iPO4 + iPOC + iDOC + iALG1 + iALG2 + iALG3 + iALG4 + iALG5 + (sedNumber * iSS) + iMeHg + iHg0 + iHgII
  if (ver .ne. ntr) then
    print*, 'ERROR - Number of tracers to model does '
    print*, ' not match list of tracer simulated in  '
    print*, ' si3d_wq_inp.txt                        '
    stop
  end if

  do i = 1, ntr
    if ((tracer_list(i) == 'DO') .and. (iDO == 1)) then
      LDO = i
    else if ((tracer_list(i) == 'PON') .and.(iPON == 1)) then
      LPON = i
    else if ((tracer_list(i) == 'DON') .and. (iDON == 1)) then
      LDON = i
    else if ((tracer_list(i) == 'NH4') .and. (iNH4 == 1)) then
      LNH4 = i
    else if ((tracer_list(i) == 'NO3') .and. (iNO3 == 1)) then
      LNO3 = i
    else if ((tracer_list(i) == 'POP') .and. (iPOP == 1)) then
      LPOP = i
    else if ((tracer_list(i) == 'DOP') .and. (iDOP == 1)) then
      LDOP = i
    else if ((tracer_list(i) == 'PO4') .and. (iPO4 == 1)) then
      LPO4 = i
    else if ((tracer_list(i) == 'POC') .and. (iPOC == 1)) then
      LPOC = i
    else if ((tracer_list(i) == 'DOC') .and. (iDOC == 1)) then
      LDOC = i
    else if ((tracer_list(i) == 'ALG1') .and. (iALG1 == 1)) then
      LALG1 = i
    else if ((tracer_list(i) == 'ALG2') .and. (iALG2 == 1)) then
      LALG2 = i
    else if ((tracer_list(i) == 'ALG3') .and. (iALG3 == 1)) then
      LALG3 = i
    else if ((tracer_list(i) == 'ALG4') .and. (iALG4 == 1)) then
      LALG4 = i
    else if ((tracer_list(i) == 'ALG5') .and. (iALG5 == 1)) then
      LALG5 = i
    else if ((tracer_list(i) == 'SS1') .and. (iSS == 1)) then
      LSS1 = i
    else if ((tracer_list(i) == 'SS2') .and. (iSS == 1)) then
      LSS2 = i
    else if ((tracer_list(i) == 'SS3') .and. (iSS == 1)) then
      LSS3 = i
    else if ((tracer_list(i) == 'Hg0') .and. (iHg0 == 1)) then
      LHg0 = i
    else if ((tracer_list(i) == 'HgII') .and. (iHgII == 1)) then
      LHgII = i
    else if ((tracer_list(i) == 'MeHg') .and. (iMeHg == 1)) then
      LMeHg = i
    else
      print*, 'ERROR - Name of water quality parameter'
      print*, 'is not set to be modeled - STOPPING CODE'
      print*, 'Revise si3d_init and si3d_wq_inp'
      stop
    end if
  end do

  IF (idbg == 1) THEN
    PRINT*,'LDO  = ', LDO
    PRINT*,'LPON = ', LPON
    PRINT*,'LDON = ', LDON
    PRINT*,'LNH4 = ', LNH4
    PRINT*,'LNO3 = ', LNO3
    PRINT*,'LPOP = ', LPOP
    PRINT*,'LDOP = ', LDOP
    PRINT*,'LPO4 = ', LPO4
    PRINT*,'LPOC = ', LPOC
    PRINT*,'LDOC = ', LDOC    
    PRINT*,'LALG1 = ', LALG1
    PRINT*,'LALG2 = ', LALG2
    PRINT*,'LALG3 = ', LALG3
    PRINT*,'LALG4 = ', LALG4
    PRINT*,'LALG5 = ', LALG5
    PRINT*,'LHg0 = ', LHg0
    PRINT*,'LHgII = ', LHgII
    PRINT*,'LMeHg = ', LMeHg
    PRINT*,'LSS1 = ', LSS1
    PRINT*,'LSS2 = ', LSS2
    PRINT*,'LSS3 = ', LSS3
  END IF

  if (iSS .eq. 1) then
    if (sedNumber .gt. 1) then
      if ((LSS2 - LSS1 .ne.  1) .or. (LSS3 - LSS2 .ne.  1 )) then
      print*,'ERROR - SS1, SS2, and SS3 must be consecutive within si3d_init'
      stop
      end if
    else
      if ((LSS1 .gt. LHg0) .and. (LHg0 .eq. 1))then
        print*,'ERROR - SS1 must be before Hg variables within si3d_init'
        stop
      end if
    end if
  end if 

END SUBROUTINE WQinit

!************************************************************
SUBROUTINE srcsnkWQ(n)
!***********************************************************
!
!   Purpose: to call all source subroutines for all cells
!
!------------------------------------------------------------

  !... Local variables
  INTEGER:: i, k, l, liter, k1s, kms, iteration
  integer, intent(in) :: n 

  ! reset soursesink = 0
  sourcesink = 0.0

  ! STWAVE controlling section. (SergioValbuena 03-11-2023) 
  if ((iSS == 1) .and. (iSTWAVE == 1)) then
    call stwave_input(n)
  end if 

  !$omp barrier  

  DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
    l = id_column(liter)
    ! ... Retrieve top & bottom wet sal-pts .................
    kms = kmz(l)
    k1s = k1z(l)

    DO k = k1s, kms;

      IF (iDO == 1) THEN
        CALL sourceDO(k,l)
      END IF
      IF (iPON == 1) THEN
        CALL sourcePON(k,l)
      END IF
      IF (iDON == 1) THEN
        CALL sourceDON(k,l)
      END IF
      IF (iNH4 == 1) THEN
        CALL sourceNH4(k,l)
      END IF
      IF (iNO3 == 1) THEN
        CALL sourceNO3(k,l)
      END IF
      IF (iPOP == 1) THEN
        CALL sourcePOP(k,l)
      END IF
      IF (iDOP == 1) THEN
        CALL sourceDOP(k,l)
      END IF
      IF (iPO4 == 1) THEN
        CALL sourcePO4(k,l)
      END IF
      IF (iDOC == 1) THEN
        CALL sourceDOC(k,l)
      END IF
      IF (iPOC == 1) THEN
        CALL sourcePOC(k,l)
      END IF
      IF (iALG1 == 1) THEN
        CALL sourceALG1(k,l)
      END IF
      IF (iALG2 == 1) THEN
        CALL sourceALG2(k,l)
      END IF
      IF (iALG3 == 1) THEN
        CALL sourceALG3(k,l)
      END IF
      IF (iALG4 == 1) THEN
        CALL sourceALG4(k,l)
      END IF
      IF (iALG5 == 1) THEN
        CALL sourceALG5(k,l)
      END IF
      IF (iSS == 1) THEN
        call sourceSS(k, l)
      END IF
      IF ((iHg0 == 1) .OR. (iHgII == 1) .OR. (iMeHg == 1)) THEN
        CALL sourceHg(k, l)
      END IF
    END DO
  END DO

END SUBROUTINE srcsnkWQ

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