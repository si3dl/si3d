!************************************************************************
                          MODULE si3d_ecomod
!************************************************************************
!
!  Purpose: Procedures that implement routines related to the modelling
!           of ecological processees
!
!-------------------------------------------------------------------------

   !USE si3d_BoundaryConditions
   USE si3d_types

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
  CHARACTER(LEN=12):: wq_input_file="si3d_wq.txt"
  INTEGER::		ios

  !.....Open input parameter file.....
  OPEN (UNIT=i99, FILE=wq_input_file, STATUS="OLD", IOSTAT=ios)
  IF (ios /= 0) CALL open_error ( "Error opening "//wq_input_file, ios )

  !.....Read header record containing comments about run................
  READ (UNIT=i99, FMT='(/(A))', IOSTAT=ios) title
  IF (ios /= 0) CALL input_error ( ios, 91)

  !. . Read list of tracerpps modeled: Dissolved oxygen, N forms, P forms, C forms and 5 algae groups
  ! ...... Representative Taxon and size range (or group)
  ! ...... 1) Picoplankton: < 2 um
  ! ...... 2) Cyclotella sp: 2-6 um
  ! ...... 3) Cryptomonas: 6-30 um
  ! ...... 4) Synedra: >30 um
  ! ...... 5) Mycrocystis (cyanobacteria)

  READ (UNIT=i99,FMT='(///(14X,I20))',IOSTAT=ios) iDO,  &
  &    iPON, iDON, iNH4, iNO3, &
  &    iPOP, iDOP, iPO4, &
  &    iPOC, iDOC, &
  &    iALG1, iALG2, iALG3, iALG4, iALG5, &
  &    iFT, iFN, iFP, iFL1, iFL2, iFL3, iFL4, iFL5
  IF (ios /= 0) CALL input_error ( ios, 92)

  !. . . Read model stochiometeric constants and other constants
  READ (UNIT=i99,FMT='(///(14X,G20.3))',IOSTAT=ios) rnc, rpc, roc, ron, &
  &     KNIT, KSN, KSP, FNH4, KDOC, SOD, KSOD, &
  &     light_sat1, light_sat2, light_sat3, light_sat4, light_sat5, &
  &     Topt1, Topt2, Topt3, Topt4, Topt5
  IF (ios /= 0) CALL input_error ( ios, 93)

  !. . . Read model rates
  READ (UNIT=i99,FMT='(///(14X,G20.2))',IOSTAT=ios) k_a,  &
  &    mu_max1, k_mor1, k_ex1, k_gr1, &
  &    mu_max2, k_mor2, k_ex2, k_gr2, &
  &    mu_max3, k_mor3, k_ex3, k_gr3, &
  &    mu_max4, k_mor4, k_ex4, k_gr4, &
  &    mu_max5, k_mor5, k_ex5, k_gr5, &
  &    k_dcn, k_mn, k_n, k_dn, &
  &    k_dcp, k_mp, k_dcc, &
  &    k_set, k_rs
  IF (ios /= 0) CALL input_error ( ios, 94)


  !... Convert model rates [1/s]. Input file has 1/day values
  ! DO
  k_a   =  k_a*idt/86400.0
  SOD   =  SOD*idt/86400.0
  ! ALG1 : Pycoplankton
  mu_max1 =  mu_max1*idt/86400.0
  k_mor1 =  k_mor1*idt/86400.0
  k_ex1 =  k_ex1*idt/86400.0
  k_gr1 =  k_gr1*idt/86400.0
  ! ALG2 : Cyclotella
  mu_max2 =  mu_max2*idt/86400.0
  k_mor2 =  k_mor2*idt/86400.0
  k_ex2 =  k_ex2*idt/86400.0
  k_gr2 =  k_gr2*idt/86400.0
  ! ALG3: Cryptomonas
  mu_max3 =  mu_max3*idt/86400.0
  k_mor3 =  k_mor3*idt/86400.0
  k_ex3 =  k_ex3*idt/86400.0
  k_gr3 =  k_gr3*idt/86400.0
  ! ALG4: Synedra
  mu_max4 =  mu_max4*idt/86400.0
  k_mor4 =  k_mor4*idt/86400.0
  k_ex4 =  k_ex4*idt/86400.0
  k_gr4 =  k_gr4*idt/86400.0
  ! ALG5: Microcystis
  mu_max5 =  mu_max5/86400.0
  k_mor5 =  k_mor5*idt/86400.0
  k_ex5 =  k_ex5*idt/86400.0
  k_gr5 =  k_gr5*idt/86400.0
  ! Nutrients
  k_dcn =  k_dcn*idt/86400.0
  k_mn =  k_mn*idt/86400.0
  k_n =  k_n*idt/86400.0
  k_dn =  k_dn*idt/86400.0
  k_dcp =  k_dcp*idt/86400.0
  k_mp =  k_mp*idt/86400.0
  k_set =  k_set*idt/86400.0
  k_rs =  k_rs*idt/86400.0
  k_dcc =  k_dcc*idt/86400.0

  !. . . Read model temperature rates
  READ (UNIT=i99,FMT='(///(14X,G20.2))',IOSTAT=ios) Theta_a, Theta_sod, Theta_mu, Theta_mor, Theta_exc, Theta_gr, &
  &     Theta_dcn, Theta_mn, Theta_n , Theta_dn , &
  &     Theta_dcp , Theta_mp , Theta_dcc , Theta_DOC

  IF (ios /= 0) CALL input_error ( ios, 95)

  !. . . Read miscillaneous rates
  READ (UNIT=i99,FMT='(///(14X,G20.2))',IOSTAT=ios) ATM_DON, ATM_NH4,    &
  &    ATM_NO3, ATM_PO4, ATM_DOP, ATM_DOC, GW_NH4, GW_NO3,  GW_PO4, GW_DOC, J_NH4, J_NO3, &
  &    J_PO4, J_DOC
  IF (ios /= 0) CALL input_error ( ios, 96)

  IF (idbg == 1) THEN
    PRINT*, "iDO  = ", iDO , "iPOC = ", iPOC, "iDOC = ", iDOC
    PRINT*, "iPON = ", iPON, "iDON = ", iDON, "iNH4 = ", iNH4, "iNO3 = ", iNO3
    PRINT*, "iPOP = ", iPOP, "iDOP = ", iDOP, "iPO4 = ", iPO4
	PRINT*, "iALG1 = ", iALG1, "iALG2 = ", iALG2, "iALG3 = ", iALG3, "iALG4 = ", iALG4, "iALG5 = ", iALG5
    PRINT*, "iFT = ", iFT, "iFN = ", iFN, "iFP = ", iFP, "iFL1 = ", iFL1, "iFL2 = ", iFL2, "iFL3 = ", iFL3, "iFL4 = ", iFL4, "iFL5 = ", iFL5
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
  LALG1=0; LALG2=0; LALG3=0; LALG4=0; LALG5=0;
  LDOC=0; LPOC=0;
  LFT=0; LFN=0; LFP=0; LFL1=0; LFL2=0; LFL3=0; LFL4=0; LFL5=0;

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

  IF (iALG5 == 1) THEN
    tracerpplocal(i) = 15
    i = i+1
  END IF

  IF (iFT == 1) THEN
    tracerpplocal(i) = 16
    i = i+1
  END IF
  
  IF (iFN == 1) THEN
    tracerpplocal(i) = 17
    i = i+1
  END IF

  IF (iFP == 1) THEN
    tracerpplocal(i) = 18
    i = i+1
  END IF

  IF (iFL1 == 1) THEN
    tracerpplocal(i) = 19
    i = i+1
  END IF

  IF (iFL2 == 1) THEN
    tracerpplocal(i) = 20
    i = i+1
  END IF

  IF (iFL3 == 1) THEN
    tracerpplocal(i) = 21
    i = i+1
  END IF

  IF (iFL4 == 1) THEN
    tracerpplocal(i) = 22
    i = i+1
  END IF

  IF (iFL5 == 1) THEN
    tracerpplocal(i) = 23
    i = i+1
  END IF
  
  sumtr = 0
  DO j = 1, ntrmax
    IF (tracerpplocal(j)>0) sumtr = sumtr + 1
  END DO

  IF (sumtr .ne. ntr) THEN
    print*,('****************************************')
    print*,('               WARNING                  ')
    print*,('    NTR AND NO. tracerpp DO NOT MATCH     ')
    print*,('           EXITING MODEL                ')
    print*,('****************************************')
    STOP
  END IF

IF (idbg == 1) THEN
  PRINT*, "ntr = ", ntr
  PRINT*, "sumtr = ", sumtr
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
    LALG5 = i
	ELSEIF (tracerpplocal(i) == 16) THEN
		LFT = i		
	ELSEIF (tracerpplocal(i) == 17) THEN
		LFN = i	
	ELSEIF (tracerpplocal(i) == 18) THEN
		LFP = i	
	ELSEIF (tracerpplocal(i) == 19) THEN
		LFL1 = i	
	ELSEIF (tracerpplocal(i) == 20) THEN
		LFL2 = i	
	ELSEIF (tracerpplocal(i) == 21) THEN
		LFL3 = i	
	ELSEIF (tracerpplocal(i) == 22) THEN
		LFL4 = i	
	ELSEIF (tracerpplocal(i) == 23) THEN
		LFL5 = i	
	END IF
  END DO

  IF (idbg == 1) THEN
    PRINT*, "LDO  = ", LDO, "LPOC = ", LPOC, "LDOC = ", LDOC
    PRINT*, "LPON = ", LPON, "LDON = ", LDON, "LNH4 = ", LNH4, "LNO3 = ", LNO3
    PRINT*, "LPOP = ", LPOP, "LDOP = ", LDOP
    PRINT*, "LALG1 = ", LALG1, "LALG2 = ", LALG2, "LALG3 = ", LALG3, "LALG4 = ", LALG4, "LALG5 = ", LALG5
    PRINT*, "LFT = ", LFT, "LFN = ", LFN, "LFP = ", LFP, "LFL1 = ", LFL1, "LFL2 = ", LFL2, "LFL3 = ", LFL3, "LFL4 = ", LFL4, "LFL5 = ", LFL5
  END IF

END SUBROUTINE WQinit

!************************************************************
SUBROUTINE srcsnkWQ(thrs)
!***********************************************************
!
!   Purpose: to call all source subroutines for all cells
!
!------------------------------------------------------------


  !... Local variables
  INTEGER:: i,j,k,l, k1s, kms,itr
  REAL, INTENT(IN) :: thrs

  ! reset soursesink = 0
  sourcesink = 0;

  DO l = 1, lm;

    ! ... Map l- into (i,j)-indexes .........................
    !i = l2i(l); j = l2j(l);

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
      IF (iALG1 == 1) THEN
        CALL sourceALG1(k,l,thrs)
      END IF
      IF (iALG2 == 1) THEN
        CALL sourceALG2(k,l,thrs)
      END IF
      IF (iALG3 == 1) THEN
        CALL sourceALG3(k,l,thrs)
      END IF
      IF (iALG4 == 1) THEN
        CALL sourceALG4(k,l,thrs)
      END IF
      IF (iALG5 == 1) THEN
        CALL sourceALG5(k,l)
      END IF
      IF (iDOC == 1) THEN
        CALL sourceDOC(k,l)
      END IF
      IF (iPOC == 1) THEN
        CALL sourcePOC(k,l)
      END IF
      IF (iFT == 1) THEN
        CALL sourceFT(k,l)
  	  END IF
  	  IF (iFN == 1) THEN
          CALL sourceFN(k,l)
  	  END IF
  	  IF (iFP == 1) THEN
          CALL sourceFP(k,l)
  	  END IF
  	  IF (iFL1 == 1) THEN
          CALL sourceFL1(k,l)
  	  END IF
  	  IF (iFL2 == 1) THEN
          CALL sourceFL2(k,l)
  	  END IF
  	  IF (iFL3 == 1) THEN
          CALL sourceFL3(k,l)
  	  END IF
  	  IF (iFL4 == 1) THEN
          CALL sourceFL4(k,l)
  	  END IF
  	  IF (iFL5 == 1) THEN
          CALL sourceFL5(k,l)
      END IF

    END DO
  END DO

END SUBROUTINE srcsnkWQ

!*********************************************************************
SUBROUTINE sourceDO(kwq,lwq)
!********************************************************************
!
!  Purpose: if dissolved oxygen is modeled, this subroutine
!  calculates source and sink terms that depend on dissolved
!  oxygen concentrations
!
!-----------------------------------------------------------------------

! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq

  !. . . Local Variables
!  REAL		::	 Tk, lnOS, OS, ln_Pwv, Pwv, theta2, Patm
   REAL     ::   OS, f_SOD 

!  ! Calculate DO saturation
!  Tk = salp(kwq,lwq) + 273
!  lnos = -139.34410 + 1.575701*1E5 /(Tk    ) &
!  &                 - 6.642308*1E7 /(Tk**2.) &
!  &	                + 1.243800*1E10/(Tk**3.) &
!  &                 - 8.621949*1E11/(Tk**4.)
!  os = EXP(lnos)

!  ! Correct for Patmospheric (Pa - declared in si3d_types and defined in surfbc0)
!  Patm   = Pa * 0.00000986923; ! Transform atmospheric pressure from Pa to atm
!  ln_Pwv = 11.8751 - 3840.70/Tk - 216961/Tk
!  Pwv    = EXP(ln_Pwv)
!  theta2 = 0.000975 - 1.426*1E-5 * salp(kwq,lwq) + &
!                      6.436*1E-8 * salp(kwq,lwq)**2
!  os = os*Pa*((1-Pwv/Pa) *(1-theta2*Pa))&
!  &           /((1-Pwv)*(1-theta2) )

   OS = 10

  ! Calculate reaertaion
  ! for now using constant rearation defined in wq_inp, but in future, can have
  ! alternatives for calcualting reaeration.
sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO)

 IF (kwq .eq. k1z(lwq)) THEN
 sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO)       &
         &       +    k_a*(OS - tracerpp(kwq,lwq,LDO))   !reaeration

 END IF


!. . Add contribution from sediment release into bottom cells

IF (kwq .eq. kmz(lwq)) THEN

 f_SOD = tracerpp(kwq,lwq,LDO) /(KSOD + tracerpp(kwq,lwq,LDO) )

 sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO)    &
						& +	 SOD * Theta_sod**(salp(kwq,lwq) - 20) * f_SOD * (1/(hpp(kwq,lwq)))    ! sediment oxygen demand
END IF

!PRINT*, "k = ", kwq,", l = ", lwq,", DO = ", sourcesink(kwq,lwq,LDO) 

END SUBROUTINE sourceDO
!************************************************************************
SUBROUTINE sourcePON(kwq,lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on PON
!
!--------------------------------------------------------------------------

! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq

  !... Local variables
  REAL:: decompositionn

  !. Calculate hydrolysis
  decompositionn = k_dcn * Theta_dcn**(salp(kwq,lwq) - 20.0) *tracerpp(kwq,lwq,LPON)

  sourcesink(kwq,lwq,LPON) = sourcesink(kwq,lwq,LPON)               &
  &                          - decompositionn			                  &   ! decomposition
  &                          - k_set* tracerpp(kwq,lwq,LPON)	      &	! settling
  &                          + k_rs * tracerpp(kwq,lwq,LPON)			! resusupension
						                ! + mortality	- only if IALG = 1; calcualted in sourceALG
                            

  ! Add contribution of mineralization to DON concentration
  IF (iDON == 1) THEN
    sourcesink(kwq,lwq,LDON) = sourcesink(kwq,lwq,LDON)	+  decompositionn
  END IF

END SUBROUTINE sourcePON


!************************************************************************
SUBROUTINE sourceDON(kwq,lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on DON
!
!--------------------------------------------------------------------------

  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq

  !. . . Local variables
  REAL:: minrln

  !. . mineralization
  minrln = k_mn * Theta_mn**(salp(kwq,lwq) - 20.0) *tracerpp(kwq,lwq,LDON)

  sourcesink(kwq,lwq,LDON) = sourcesink(kwq,lwq,LDON)					&
  &                         -  minrln
							! + decomposition	- only if IPON = 1; caluclated in sourcePON
							! + Excr	    - only if IAlG = 1; calculated in sourceALG

!. . .Add contribution from atmospheric deposition to top layer
IF (kwq .eq. k1z(lwq)) THEN

	sourcesink(kwq,lwq, LDON) = sourcesink(kwq, lwq, LDON) + &
								& (hpp(kwq,lwq))*ATM_DON		! Atmoshperic deposition
END IF

  IF (INH4 == 1) THEN		! add mineralization of DON to NH4
    sourcesink(kwq,lwq,LNH4) = sourcesink(kwq,lwq,LNH4) + minrln
  END IF

END SUBROUTINE sourceDON

!************************************************************************
SUBROUTINE sourceNH4(kwq,lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on NH4
!
!--------------------------------------------------------------------------

! ... Arguments
   INTEGER, INTENT (IN) :: kwq,lwq

!. . . Local Variables
REAL:: nitrif, f_DO

!. . . Calculate nitrification
  ! Calculate DO inhibition of nitrification
	IF (IDO == 1) THEN
		f_DO = tracerpp(kwq,lwq,LDO) /(KNIT + tracerpp(kwq,lwq,LDO) )
	ELSE
		f_DO = 1.0
	END IF

	nitrif = k_n* Theta_n**(salp(kwq,lwq) - 20.0) *tracerpp(kwq,lwq,LNH4)

sourcesink(kwq,lwq,LNH4) = sourcesink(kwq,lwq,LNH4)			&
						&  -  nitrif			          ! nitrificatin

						! - Algal Uptake		- if IALG = 1; calculated in sourceALG
						! + mineralization		- if IDON = 1; calculated in sourceDON

!. . Add contribution from sediment release and GW flux into bottom cells
IF (kwq .eq. kmz(lwq)) THEN
 sourcesink(kwq,lwq,LNH4) = sourcesink(kwq,lwq,LNH4) +   &
							& (hpp(kwq,lwq))*J_NH4     +	 &	! sediment release
							& (hpp(kwq,lwq))*GW_NH4	     	! GW flux
END IF

!. . Add contribution from atmoshperic deposition into top cells
IF (kwq .eq. k1z(lwq)) THEN
	sourcesink(kwq, lwq, LNH4) = sourcesink(kwq, lwq, LNH4) +(hpp(kwq, lwq)) * ATM_NH4	 
END IF


!. . Add contribution from nitrification to nitrate sourcesink

IF (INO3 == 1) THEN
	sourcesink(kwq,lwq,LNO3) = sourcesink(kwq,lwq,LNO3) + nitrif
END IF

IF (IDO == 1) THEN
	sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) - ron*nitrif
END IF

END SUBROUTINE sourceNH4

!************************************************************************
SUBROUTINE sourceNO3(kwq,lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on NO3
!
!--------------------------------------------------------------------------

! ... Arguments
   INTEGER, INTENT (IN) :: kwq,lwq


sourcesink(kwq,lwq,LNO3) = sourcesink(kwq,lwq,LNO3)		-	&
						&  k_dn* Theta_dn**(salp(kwq,lwq) - 20.0) *tracerpp(kwq,lwq,LNO3)	! denitrification
						! + nitrification		- if INH4 = 1; calculated in sourceNH4
						! - algal uptake		- if IALG = 1; calculated in sourceALG

!. . Add contribution from sediment release and GW flux into bottom cells
IF (kwq .eq. kmz(lwq)) THEN
 sourcesink(kwq,lwq,LNO3) = sourcesink(kwq,lwq,LNO3) +   &
							& (hpp(kwq,lwq))*J_NO3     +	 &	! sediment release
							& (hpp(kwq,lwq))*GW_NO3	     	! GW flux
END IF

!. . Add contribution from atmoshperic deposition into top cells
IF (kwq .eq. k1z(lwq)) THEN
	sourcesink(kwq, lwq, LNO3) = sourcesink(kwq, lwq, LNO3) + &
								& (hpp(kwq, lwq)) * ATM_NO3	 ! ATM deposition
END IF

END SUBROUTINE sourceNO3

!************************************************************************
SUBROUTINE sourcePOP(kwq,lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on POP
!
!--------------------------------------------------------------------------

! ... Arguments
   INTEGER, INTENT (IN) :: kwq,lwq

!. . . Local Variables
REAL:: decomposition

! Calculate decomposition

decomposition = k_dcp* Theta_dcp**(salp(kwq,lwq) - 20.0) *tracerpp(kwq,lwq,LPOP)

sourcesink(kwq,lwq,LPOP) = sourcesink(kwq,lwq,LPOP)		    &
						&   - decomposition                           & ! decomposition
						&   - k_set * tracerpp(kwq, lwq, LPOP)        &  ! settling
						&   + k_rs*tracerpp(kwq,lwq,LPOP)			          ! Resuspension
						! + mortality	- only if IALG = 1; calcualted in sourceALG

! Caclulate decomposition contribution to DOP
IF (IDOP == 1) THEN
	sourcesink(kwq,lwq,LDOP) = sourcesink(kwq,lwq,LDOP) + decomposition
END IF

END SUBROUTINE sourcePOP

!************************************************************************
SUBROUTINE sourceDOP(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on DOP
!
!--------------------------------------------------------------------------

! ... Arguments
   INTEGER, INTENT (IN) :: kwq,lwq

!. . .Local Variables
REAL:: minrl

! Calculate mineralization

minrl = k_mp* Theta_mp**(salp(kwq,lwq) - 20.0) *tracerpp(kwq,lwq,LDOP)

sourcesink(kwq,lwq,LDOP) = sourcesink(kwq,lwq,LDOP)		&
						&   -  minrl				                  ! decomposition
						! + excretion	      - if IAlG = 1; calculated in sourceALG
						! + decomposition		- if iPOP = 1; calculated in sourcePOP

!... Add atmoshperic deposition to top layer
 IF (kwq .eq. k1z(lwq)) THEN
	sourcesink(kwq,lwq,LDOP) = sourcesink(kwq,lwq,LDOP)				+	&
						&  (h(kwq,lwq))*ATM_DOP					! Atmospheric deposition
END IF

! Caclulate mineralization contribution to PO4
IF (IPO4 == 1) THEN
	sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) + minrl
END IF

END SUBROUTINE sourceDOP

!************************************************************************
SUBROUTINE sourcePO4(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on PO4
!
!--------------------------------------------------------------------------

   ! ... Arguments
   INTEGER, INTENT (IN) :: kwq,lwq

   sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4)
 ! + mineralization		- if IDOP = 1; calculated in sourceDOP
   ! + algal exc			- if IALG = 1; calculated in sourceALG
   ! - algal uptake		- if IALG = 1; calculated in sourceALG

!. . Add contribution from sediment release and GW flux into bottom cells
IF (kwq .eq. kmz(lwq)) THEN
 sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) +   &
							& (hpp(kwq,lwq))*J_PO4     +	 &	! sediment release
							& (hpp(kwq,lwq))*GW_PO4	     	! GW flux
END IF

!. . Add contribution from atmoshperic deposition into top cells
IF (kwq .eq. k1z(lwq)) THEN
	sourcesink(kwq, lwq, LPO4) = sourcesink(kwq, lwq, LPO4) + &
								& (hpp(kwq, lwq)) * ATM_PO4	 ! ATM deposition
END IF

END SUBROUTINE sourcePO4

!************************************************************************
SUBROUTINE sourcePOC (kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on POC
!
!--------------------------------------------------------------------------

  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq

  !. . Local Variables
  REAL:: decomposition

  ! Calculate decomposition
  decomposition = k_dcc*Theta_dcc**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LPOC)

  sourcesink(kwq,lwq,LPOC) = sourcesink(kwq,lwq,LPOC)	&
						&  -   decomposition						          &              ! decomposition
						&  -   k_set * tracerpp(kwq,lwq,LPOC) 	  &	             ! settling
						&  +   k_rs * tracerpp(kwq,lwq,LPOC)			                 ! Resuspension
						! + mortality		                        - if IALG = 1; calculated in sourceALG

  ! Caclulate decomposition contribution to DOC
  IF (IDOC == 1) THEN
    sourcesink(kwq,lwq,LDOC) = sourcesink(kwq,lwq,LDOC) + decomposition
   END IF

END SUBROUTINE sourcePOC

!************************************************************************
SUBROUTINE sourceDOC(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on DOC
!
!--------------------------------------------------------------------------

  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq

  !. . .Local Variables
  REAL:: F_DO, oxid

  !. . Calculate how DO concn impedes oxidation of DOM
  IF (IDO ==1) THEN
    F_DO = (tracerpp(kwq,lwq,LDO))/(KDOC + tracerpp(kwq,lwq,LDO) )
  END IF

  ! Calculate Microbial Uptake (oxidation)
  oxid = KDOC * Theta_DOC**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LDOC) * F_DO

  sourcesink(kwq,lwq,LDOC) = sourcesink(kwq,lwq,LDOC)  &
                     &    - oxid                         ! microbial uptake (oxidation)
						! + decomposition	                -if IPOC = 1; calculated in sourcePOC
						! + algal excretion             - if IALG = 1; calculated in sourceALG


!. . Add contribution from sediment release and GW flux into bottom cells
 IF (kwq .eq. kmz(lwq)) THEN
 sourcesink(kwq,lwq,LDOC) = sourcesink(kwq,lwq,LDOC) +   &
							& (hpp(kwq,lwq))*J_DOC     +	 &	! sediment release
							& (hpp(kwq,lwq))*GW_DOC	     	! GW flux
END IF

!. . Add contribution from atmoshperic deposition into top cells
 IF (kwq .eq. k1z(lwq)) THEN
	sourcesink(kwq, lwq, LDOC) = sourcesink(kwq, lwq, LDOC) + &
								& (hpp(kwq, lwq)) * ATM_DOC	 ! ATM deposition
END IF

  IF (IDO == 1) THEN
    sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) - roc*oxid
  END IF

END SUBROUTINE sourceDOC

!************************************************************************
SUBROUTINE sourceALG1(kwq, lwq, thrs)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on ALG1
!
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq
  REAL, INTENT(IN) :: thrs

  !. . .Local Variables
  REAL::	mu1, f_L1, f_T, f_N, f_P, N_conc
  REAL::	growth1, mort1, excr1,  graz1, sett1, resus1
  REAL::  dthrs_zoop

  ! Calculate mu, growth rate

  !. .  Calculate growth limiting factors
    ! Light Limitation - by Steele equation (Jassby and Platt, 1976)
   	  f_L1 = ((Qsw*QswFr(kwq,lwq)*0.47)/light_sat1) *EXP(1 -((Qsw*QswFr(kwq,lwq)*0.47)/light_sat1)) ! Steele (1962) = Photoinhibited
      !f_L1 = 1 - EXP(-(Qsw*QswFr(kwq,lwq)*0.47)/light_k1) ! Webb et al (1974) in absence of photoinhibition 

	! temperature limitaton
		f_T = Theta_mu**(-0.01*(salp(kwq,lwq) - Topt1)**2)

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
	  
	  
! Growth rate considering limiting factors (light and nutrients)
mu1 = mu_max1 * MIN(f_L1,f_N,f_P)

! Read Zooplankton. Daily values between April 13 and June 15, 2018
    rotifers =  (/192.54,  196.37,  200.26,  204.19,  208.16,  212.17,  216.21,  220.29,  224.39,  228.51,  232.65,  236.80,  240.96,  245.12,  249.28,  253.44,  257.58,  &
  &    261.72,  265.83,  269.93,  273.99,  278.03,  282.03,  285.99,  289.91,  293.78,  297.60,  301.36,  305.06,  308.69,  312.26,  315.75,  319.16,  322.49,  325.73,  328.88,  &
  &    331.94,  334.89,  337.74,  340.49,  343.12,  345.63,  348.12,  350.68,  353.30,  355.97,  358.69,  361.45,  364.24,  367.05,  369.87,  372.69,  375.52,  378.33,  381.12,  383.89,  &
  &    386.61,  389.30,  391.93,  394.51,  397.01,  399.44,  401.79,  404.04/)
    rotifers = rotifers/1000

    bosmina = (/0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
    & 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
    & 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
    & 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
    & 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00/)
    bosmina = bosmina/1000

! Interpolate Zooplankton to the simulation time interval (thrs)
dthrs_zoop = 24.
rotifersnew = parabn(0.,thrs,rotifers,dthrs_zoop)
bosminanew = parabn(0.,thrs,bosmina,dthrs_zoop)

!. . Calculate growth
		growth1 = mu1 * f_T * tracerpp(kwq,lwq,LALG1)
!. . Calculate mortality
    mort1   = k_mor1 * Theta_mor**(-0.01*(salp(kwq,lwq) - Topt1)**2) * tracerpp(kwq,lwq,LALG1)
!. . Calculate excretion
		excr1   = k_ex1 * Theta_exc**(-0.01*(salp(kwq,lwq) - Topt1)**2) * tracerpp(kwq,lwq,LALG1)
!. . Calculate grazing
		!graz1   = k_gr1 * Theta_gr**(-0.01*(salp(kwq,lwq) - Topt1)**2) * tracerpp(kwq,lwq,LALG1)
    graz1 = (rotifersnew*0.0003936/4.0 + bosminanew*0.00495/3.0)*(86400.0)! LT2
!. . Calculate settling
		sett1   = k_set * tracerpp(kwq,lwq,LALG1)
!. . Calculate resuspension
    resus1   = k_rs * tracerpp(kwq,lwq,LALG1)

sourcesink(kwq,lwq,LALG1) = sourcesink(kwq,lwq,LALG1)	+ growth1 -  mort1  - graz1 - sett1 + resus1


! If dissolved oxygen is modeled, alter sourcesink(kwq,lwq,LDO) to reflect growth and resp
! of algae population
IF (IDO ==1) THEN
	sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) + roc*(growth1)
END IF

! IF PON is modeled, alter sourcesink(kwq,lwq,LPON) to reflect mortality of algae
IF (IPON == 1) THEN
	sourcesink(kwq,lwq,LPON) = sourcesink(kwq,lwq,LPON) + rnc*mort1
END IF

! If DON is modeled, alter soucesink(kwq,lwq,LDON) to reflect excretion of algae (used mortality)
IF (IDON ==1) THEN
	sourcesink(kwq,lwq,LDON) = sourcesink(kwq,lwq,LDON)  + rnc*mort1
END IF

! IF NH4 is modeled, alter sourcesink(kwq,lwq,LNH4) to include uptake of NH4 by algae
IF (INH4 == 1) THEN
	sourcesink(kwq,lwq,LNH4) = sourcesink(kwq,lwq,LNH4) - rnc*FNH4*growth1
END IF

! IF NO3 is modeled, alter sourcesink(kwq,lwq,NO3) to include uptake of NO3 by algae
IF (INO3 == 1) THEN
	sourcesink(kwq,lwq,LNO3) = sourcesink(kwq,lwq,LNO3) - rnc*(1-FNH4) * growth1
END IF

! If POP is modeled, alter sourcesink(kwq,lwq,POP) to include mortality of algae
IF (IPOP == 1) THEN
	sourcesink(kwq,lwq,LPOP) = sourcesink(kwq,lwq,LPOP) + rpc*mort1
END IF

! If DOP is modeled, alter sourcesink(kwq,lwq,LDOP) to include excretion (used mortality)
IF (IDOP == 1) THEN
	sourcesink(kwq,lwq,LDOP) = sourcesink(kwq,lwq,LDOP) + rpc*mort1
END IF

! If PO4 is modeled, alter sourcesink(kwq,lwq,LPO4) to include excretion and uptake
IF (IPO4 == 1) THEN
	sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) - rpc* growth1 ! excr from phyto removed
END IF

! If POC is modeled, alter sourcesink(kwq,lwq,LPOC) to include mortality
IF (IPOC == 1) THEN
	sourcesink(kwq,lwq,LPOC) = sourcesink(kwq,lwq,LPOC) + mort1
END IF

! If DOC is modeled, alter sourcesink(kwq,lwq,LDOC) to include excretion (used mortality)
IF (IDOC == 1) THEN
	sourcesink(kwq,lwq,LDOC) = sourcesink(kwq,lwq,LDOC) + mort1
END IF


END SUBROUTINE sourceALG1

!************************************************************************
SUBROUTINE sourceALG2(kwq, lwq, thrs)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on ALG2
!
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq
  REAL, INTENT(IN) :: thrs

  !. . .Local Variables
  REAL::  mu2, f_L2, f_T, f_N, f_P, N_conc
  REAL::  growth2, mort2, excr2,  graz2, sett2, resus2
  REAL::  dthrs_zoop, sources2, sinks2, var2
  REAL::  minimum_nut2, alg_min

  ! Calculate mu, growth rate

  !. .  Calculate growth limiting factors
    ! Light Limitation - by Steele equation (Jassby and Platt, 1976)
      f_L2 = ((Qsw*QswFr(kwq,lwq)*0.47)/light_sat2) *EXP(1 -((Qsw*QswFr(kwq,lwq)*0.47)/light_sat2)) ! Steele (1962) = Photoinhibited
      !f_L2 = 1 - EXP(-(Qsw*QswFr(kwq,lwq)*0.47)/light_k2) ! Webb et al (1974) in absence of photoinhibition 

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
    
    
! Growth rate considering limiting factors (light and nutrients)
mu2 = mu_max2 * MIN(f_L2,f_N,f_P)

! Read Zooplankton. Daily values between April 13 and June 15, 2018
    rotifers =  (/192.54,  196.37,  200.26,  204.19,  208.16,  212.17,  216.21,  220.29,  224.39,  228.51,  232.65,  236.80,  240.96,  245.12,  249.28,  253.44,  257.58,  &
  &    261.72,  265.83,  269.93,  273.99,  278.03,  282.03,  285.99,  289.91,  293.78,  297.60,  301.36,  305.06,  308.69,  312.26,  315.75,  319.16,  322.49,  325.73,  328.88,  &
  &    331.94,  334.89,  337.74,  340.49,  343.12,  345.63,  348.12,  350.68,  353.30,  355.97,  358.69,  361.45,  364.24,  367.05,  369.87,  372.69,  375.52,  378.33,  381.12,  383.89,  &
  &    386.61,  389.30,  391.93,  394.51,  397.01,  399.44,  401.79,  404.04/)
    rotifers = rotifers/1000.0

    bosmina = (/0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
    & 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
    & 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
    & 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
    & 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00/)
    bosmina = bosmina/1000.0

    codotisnapuli = (/147.54,  150.47,  153.49,  156.60,  159.80,  163.10,  166.48,  169.96,  173.52,  177.18,  180.92,  184.76,  188.69,  192.70,  196.81,  201.00, 205.29,  209.67,  214.13, 218.68, &  
    &  223.33,  228.06,  232.88,  237.79,  242.79,  247.88,  253.06,  258.33,  263.68,  269.12,  274.65,  280.27,  285.98,  291.77,  297.66,  303.63,  309.68,  315.83,  322.06, &
    &  328.38,  334.79,  341.28,  349.34,  360.36,  374.17,  390.63,  409.58,  430.86,  454.34,  479.85,  507.23,  536.35,  567.04,  599.16,  632.54,  667.04,  702.51,  738.78,  775.71, &
    &  813.15,  850.94,  888.94,  926.97,  964.90/)
    codotisnapuli = codotisnapuli/1000.0

! Interpolate Zooplankton to the simulation time interval (thrs)
dthrs_zoop = 24.
rotifersnew = parabn(0.,thrs,rotifers,dthrs_zoop)
bosminanew = parabn(0.,thrs,bosmina,dthrs_zoop)
codotisnapulinew = parabn(0.,thrs,codotisnapuli,dthrs_zoop)

!. . Calculate growth
    growth2 = mu2 * f_T * tracerpp(kwq,lwq,LALG2)
!. . Calculate mortality
    mort2   = k_mor2 * Theta_mor**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LALG2)
!. . Calculate excretion
    excr2   = k_ex2 * Theta_exc**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LALG2)
!. . Calculate grazing
    !!graz2   = k_gr2 * Theta_gr**(-0.01*(salp(kwq,lwq) - Topt2)) * tracerpp(kwq,lwq,LALG2) ! NOT USED
    IF (kwq .lt. 37) THEN
        graz2 = k_gr2 * f_T * (((rotifersnew*0.0003936/4.0) + (bosminanew*0.00495/3.0) + (codotisnapulinew*0.185/2.0))*(idt/86400.0))! 2T6
    ELSE IF (kwq .gt. 36) THEN
        graz2 = 0
    END IF

!. . Calculate settling
    sett2   = k_set * tracerpp(kwq,lwq,LALG2)
!. . Calculate resuspension
    resus2   = k_rs * tracerpp(kwq,lwq,LALG2)


 ! Estimate the minimum alloable nutrient value and recalculate growth if nutrient < 0
 
 !minimum_nut2 = MIN(sourcesink(kwq,lwq,LPO4)-(rpc*growth2),sourcesink(kwq,lwq,LNO3)-(rnc*growth2*(1-FNH4)),sourcesink(kwq,lwq,LNH4)-(rnc*growth2*FNH4))
 !IF (minimum_nut2 .lt. 0.0) THEN
 !   growth2 = 0.9*(MIN(sourcesink(kwq,lwq,LPO4)/rpc, sourcesink(kwq,lwq,LNO3)/(rnc*(1-FNH4)), sourcesink(kwq,lwq,LNH4)/(rnc*FNH4)))
 !      IF (growth2 .lt. 0.0) THEN
 !            growth2 = 0
 !        END IF
 !END IF

! To account for paralelization
growth2 = growth2/8.0
sourcesink(kwq,lwq,LALG2) = sourcesink(kwq,lwq,LALG2) + growth2

 ! Limit the minimum algae value (0.01) for mortality
! alg_min = 0.01/86400
! IF ((sourcesink(kwq,lwq,LALG2) - mort2) .le. alg_min) THEN
!      mort2 = 0.9 * (sourcesink(kwq,lwq,LALG2) - alg_min) !  In run 12C the following IF loop is commented and 0.9 factor, Uncommented in 12D. In 12B states term =0 if this calculation is negative (org)
!        IF (mort2 .lt. 0.0) THEN
!              mort2 = 0
!         END IF
! END IF

!! Estimate the minimum alloable nutrient value and recalculate mortality if nutrient < 0
!! minimum_nut2 = 0.0
!! minimum_nut2 = MIN(sourcesink(kwq,lwq,LALG2)-mort2,sourcesink(kwq,lwq,LALG2)*rnc + growth2*rnc - mort2*rnc, sourcesink(kwq,lwq,LALG2)*rpc 
!+ growth2*rpc - mort2*rpc)
!! IF (minimum_nut2 .lt. 0.0) THEN
!!     mort2 = 0.9*(MIN(sourcesink(kwq,lwq,LALG2), (sourcesink(kwq,lwq,LALG2)*rnc + growth2*rnc)/rnc,(sourcesink(kwq,lwq,LALG2)*rpc + growth2*rpc)/rpc))
!!         IF (mort2 .lt. 0.0) THEN
!!             mort2 = 0.0
!!         END IF
!!END IF

! To account for paralelization
mort2 = mort2 /8.0
sourcesink(kwq,lwq,LALG2) = sourcesink(kwq,lwq,LALG2) - mort2

 ! Limit the minimum algae value (0.01) for grazing
! alg_min = 0.01/86400
! IF ((sourcesink(kwq,lwq,LALG2) - graz2) .le. alg_min) THEN
!      graz2 = 0.1 * (sourcesink(kwq,lwq,LALG2) - alg_min) !  In run 12C the following IF loop is commented and 0.9 factor,
! Uncommented in 12D. In 12B states term =0 if this calculation is negative (org)
 !        IF (graz2 .lt. 0.0) THEN
 !              graz2 = 0
 !        END IF
 !END IF

 ! To account for paralelization
graz2 = graz2/8.0
IF ((tracerpp(kwq,lwq,LALG2)- graz2) .le. 0.01) THEN
    graz2 = 0.0
END IF
sourcesink(kwq,lwq,LALG2) = sourcesink(kwq,lwq,LALG2) - graz2

!!! Estimate the minimum alloable nutrient value and recalculate grazing if nutrient < 0
!! minimum_nut2 = 0.0
!! minimum_nut2 = MIN(sourcesink(kwq,lwq,LALG2)-graz2,sourcesink(kwq,lwq,LALG2)*rnc + growth2*rnc - mort2*rnc - graz2*rnc, sourcesink(kwq,lwq,LALG2)*rpc 
!+ growth2*rpc - mort2*rpc - graz2*rpc)
!! IF (minimum_nut2 .lt. 0.0) THEN
!!    graz2 = 0.9*(MIN(sourcesink(kwq,lwq,LALG2), (sourcesink(kwq,lwq,LALG2)*rnc + growth2*rnc - mort2*rnc)/rnc,(sourcesink(kwq,lwq,LALG2)*rpc 
!+ growth2*rpc - mort2*rpc)/rpc)
!!        IF (graz2 .lt. 0.0) THEN
!!             graz2 = 0.0
!!         END IF
!!END IF


!sourcesink(kwq,lwq,LALG2) = sourcesink(kwq,lwq,LALG2) - sett2 + resus2
!!sourcesink(kwq,lwq,LALG2) = sourcesink(kwq,lwq,LALG2) + growth2 - mort2  graz2 - sett2 + resus2
var2 = sourcesink(kwq,lwq,LALG2)
   IF (lwq .eq. 300) THEN
      IF (kwq .eq. 5) THEN
       !PRINT *, 'lwq2 = ', lwq, ', kwq2 = ', kwq, ', growth2 = ', growth2, ', mort2 = ', mort2, 'grazing2 = ', graz2, ', SS2 = ', var2, ', tracer2 = ',tracerpp(kwq,lwq,LALG2)
       !PRINT *, 'lwq2 = ', lwq, ', kwq2 = ', kwq, ', growth2 = ', growth2, ', SS2 = ', var2, ', tracer2 = ',tracerpp(kwq,lwq,LALG2), ', f_L2 = ', f_L2, ', f_N2 = ', f_N, ', f_P2 = ', f_P, ', !PARpenetrates2 = ',Qsw*QswFr(kwq,lwq)*0.47
      END IF
   END IF


! If dissolved oxygen is modeled, alter sourcesink(kwq,lwq,LDO) to reflect growth and resp
! of algae population
IF (IDO ==1) THEN
  sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) + roc*(growth2)
END IF

! IF PON is modeled, alter sourcesink(kwq,lwq,LPON) to reflect mortality of algae
IF (IPON == 1) THEN
  sourcesink(kwq,lwq,LPON) = sourcesink(kwq,lwq,LPON) + rnc*mort2
END IF

! If DON is modeled, alter soucesink(kwq,lwq,LDON) to reflect excretion of algae (used motality as a sink)
IF (IDON ==1) THEN
  sourcesink(kwq,lwq,LDON) = sourcesink(kwq,lwq,LDON)  + rnc*mort2
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

! If DOP is modeled, alter sourcesink(kwq,lwq,LDOP) to include excretion (used motality as a sink)
IF (IDOP == 1) THEN
  sourcesink(kwq,lwq,LDOP) = sourcesink(kwq,lwq,LDOP) + rpc*mort2
END IF

! If PO4 is modeled, alter sourcesink(kwq,lwq,LPO4) to include excretion and uptake
IF (IPO4 == 1) THEN
  sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) - rpc* growth2 ! excr from phyto removed
END IF

! If POC is modeled, alter sourcesink(kwq,lwq,LPOC) to include mortality
IF (IPOC == 1) THEN
  sourcesink(kwq,lwq,LPOC) = sourcesink(kwq,lwq,LPOC) + mort2
END IF

! If DOC is modeled, alter sourcesink(kwq,lwq,LDOC) to include excretion (used motality as a sink)
IF (IDOC == 1) THEN
  sourcesink(kwq,lwq,LDOC) = sourcesink(kwq,lwq,LDOC) + mort2
END IF


END SUBROUTINE sourceALG2


!************************************************************************
SUBROUTINE sourceALG3(kwq, lwq, thrs)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on ALG3
!
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq
  REAL, INTENT(IN) :: thrs

  !. . .Local Variables
  REAL::  mu3, f_L3, f_T, f_N, f_P, N_conc
  REAL::  growth3, mort3, excr3,  graz3, sett3, resus3
  REAL::  dthrs_zoop, sources3, sinks3, var3
  REAL::  minimum_nut3, alg_min


  ! Calculate mu, growth rate

  !. .  Calculate growth limiting factors
    ! Light Limitation - by Steele equation (Jassby and Platt, 1976)
      f_L3 = ((Qsw*QswFr(kwq,lwq)*0.47)/light_sat3) *EXP(1 -((Qsw*QswFr(kwq,lwq)*0.47)/light_sat3)) ! Steele (1962) = Photoinhibited
      !f_L3 = 1 - EXP(-(Qsw*QswFr(kwq,lwq)*0.47)/light_k3) ! Webb et al (1974) in absence of photoinhibition 

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
    
    
! Growth rate considering limiting factors (light and nutrients)
mu3 = mu_max3 * MIN(f_L3,f_N,f_P)

! Read Zooplankton. Daily values between April 13 and June 15, 2018
    rotifers =  (/192.54,  196.37,  200.26,  204.19,  208.16,  212.17,  216.21,  220.29,  224.39,  228.51,  232.65,  236.80,  240.96,  245.12,  249.28,  253.44,  257.58,  &
  &    261.72,  265.83,  269.93,  273.99,  278.03,  282.03,  285.99,  289.91,  293.78,  297.60,  301.36,  305.06,  308.69,  312.26,  315.75,  319.16,  322.49,  325.73,  328.88,  &
  &    331.94,  334.89,  337.74,  340.49,  343.12,  345.63,  348.12,  350.68,  353.30,  355.97,  358.69,  361.45,  364.24,  367.05,  369.87,  372.69,  375.52,  378.33,  381.12,  383.89,  &
  &    386.61,  389.30,  391.93,  394.51,  397.01,  399.44,  401.79,  404.04/)
    !rotifers = rotifers/1000.0

    diatomus =(/91.52,  96.11,  100.74,  105.42,  110.13,  114.87,  119.63,  124.41,  129.21,  134.00, 138.80,  143.59,  148.37,  153.13,  157.86,  162.57,  167.24,  171.86,  176.44,  180.96, &
    &  185.42, 189.81,  194.13,  198.37,  202.53,  206.59,  210.56,  214.42,  218.18,  221.82,  225.34,  228.73,  231.99,  235.11,  238.08,  240.91,  243.57,  246.07,  248.41,  250.57,  252.54, &
    &  254.33,  255.99,  257.60,  259.15,  260.65,  262.11,  263.51,  264.87,  266.19,  267.48,  268.73,  269.94,  271.13,  272.29,  273.42,  274.53,  275.62,  276.70,  277.76,  278.82,  279.86, &       
    &  280.90,281.93/)
    diatomus = diatomus/1000.0

    bosmina = (/0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
    & 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
    & 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
    & 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
    & 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00/)
    bosmina = bosmina/1000.0

    daphnia =(/0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
    & 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
    & 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
    & 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
    & 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00/)
    daphnia = daphnia/1000.0

    epischura = (/317.73,  333.08,  348.53,  364.05,  379.62,  395.21,  410.80,  426.37,  441.89,  457.34,  472.69,  487.92,  503.01,  517.92,  532.65,  547.16,  561.43,  575.43,  589.14,  &
    &  602.55, 615.61,  628.31,  640.63,  652.54,  664.02,  675.04,  685.58,  695.61,  705.11,  714.06,  722.43,  730.20,  737.35,  743.84,  749.66,  754.78,  759.18,  762.84,  765.72,  767.81, &
    &  769.08,  769.51, 768.93,  767.23,  764.49,  760.75,  756.10,  750.59,  744.30,  737.28,  729.60,  721.34,  712.54,  703.29,  693.64,  683.67,  673.43,  663.00, 652.43,  641.80,  631.17, &
    &  620.60,  610.17,  599.93/)
    epischura = epischura/1000.0

    mysis = (/0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
    & 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
    & 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
    & 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
    & 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00/)
    mysis = mysis/1000.0

! Interpolate Zooplankton to the simulation time interval (thrs)
dthrs_zoop = 24.0
rotifersnew = parabn(0.,thrs,rotifers,dthrs_zoop)
diatomusnew = parabn(0.,thrs,diatomus,dthrs_zoop)
bosminanew = parabn(0.,thrs,bosmina,dthrs_zoop)
daphnianew = parabn(0.,thrs,daphnia,dthrs_zoop)
epischuranew = parabn(0.,thrs,epischura,dthrs_zoop)
mysisnew = parabn(0.,thrs,mysis,dthrs_zoop)

!. . Calculate growth
    growth3 = mu3 * f_T * tracerpp(kwq,lwq,LALG3)
!. . Calculate mortality
    mort3   = k_mor3 * Theta_mor**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LALG3)
!. . Calculate excretion
    excr3   = k_ex3 * Theta_exc**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LALG3)
!. . Calculate grazing
    !!graz3   = k_gr3 * Theta_gr**(salp(kwq,lwq) - 20.0) * tracerpp(kwq,lwq,LALG3)
    IF (kwq .lt. 37) THEN
    graz3 = k_gr3 * f_T *(((rotifersnew*0.002035/4.0)+(diatomusnew*0.2267)+(bosminanew*0.00495/3.0)+(daphnianew*0.04752)+(epischuranew*0.078)+(mysisnew*0.0))*(idt/86400.0))!
    ELSE IF (kwq .gt. 36) THEN
        graz3 = 0
    END IF
!. . Calculate settling
    sett3   = k_set * tracerpp(kwq,lwq,LALG3)
!. . Calculate resuspension
    resus3   = k_rs * tracerpp(kwq,lwq,LALG3)

 ! Estimate the minimum alloable nutrient value and recalculate growth if nutrient < 0 . Uncomment in v12C
 !minimum_nut3 = MIN(sourcesink(kwq,lwq,LPO4)-(rpc*growth3),sourcesink(kwq,lwq,LNO3)-(rnc*growth3*(1-FNH4)),sourcesink(kwq,lwq,LNH4)-(rnc*growth3*FNH4))
 !IF (minimum_nut3 .lt. 0.0) THEN
 !   growth3 = 0.9*(MIN(sourcesink(kwq,lwq,LPO4)/rpc, sourcesink(kwq,lwq,LNO3)/(rnc*(1-FNH4)), sourcesink(kwq,lwq,LNH4)/(rnc*FNH4)))
 !       IF (growth3 .lt. 0.0) THEN ! In v12c, next line equal to zero
 !            growth3 = 0
 !      END IF
 !END IF

! To account for paralelization
growth3 = growth3 /8.0 ! comment in V12C
sourcesink(kwq,lwq,LALG3) = sourcesink(kwq,lwq,LALG3) + growth3

 ! Limit the minimum algae value (0.01) for mortality
 !alg_min = 0.01/86400
 !IF ((sourcesink(kwq,lwq,LALG3) - mort3) .le. alg_min) THEN ! Uncomment next 1 lines in v12C
 !     mort3 = 0.9 * (sourcesink(kwq,lwq,LALG3) - alg_min) ! 
 !        IF (mort3 .lt. 0.0) THEN ! Loop commented in v12C
 !             mort3 = 0
 !        END IF
 !END IF

!! Estimate the minimum alloable nutrient value and recalculate mortality if nutrient < 0. Uncomment in v12C
!! minimum_nut3 = 0.0
!! minimum_nut3 = MIN(sourcesink(kwq,lwq,LALG3)-mort3,sourcesink(kwq,lwq,LALG3)*rnc + growth3*rnc - mort3*rnc, sourcesink(kwq,lwq,LALG3)*rpc + growth3*rpc - mort3*rpc)
!! IF (minimum_nut3 .lt. 0.0) THEN
!!    mort3 = 0.9*(MIN(sourcesink(kwq,lwq,LALG3), (sourcesink(kwq,lwq,LALG3)*rnc + growth3*rnc)/rnc,(sourcesink(kwq,lwq,LALG3)*rpc + growth3*rpc)/rpc))
!!         IF (mort3 .lt. 0.0) THEN ! Keep commented in v12c
!!             mort3 = 0.0
!!         END IF
!!END IF

! To account for paralelization
mort3 = mort3 /8.0 ! Comment in v12C
sourcesink(kwq,lwq,LALG3) = sourcesink(kwq,lwq,LALG3) - mort3

 ! Limit the minimum algae value (0.01) for grazing
! alg_min = 0.01/86400 
! IF ((sourcesink(kwq,lwq,LALG3) - graz3) .le. alg_min) THEN   ! Uncomment next 1 lines in v12C
!      graz3 = 0.9 * (sourcesink(kwq,lwq,LALG3) - alg_min) !  In run 12C the following IF loop is commented and 0.9 factor, Uncommented in 12D. In 12B states term =0 if this calculation is negative (org)
!         IF (graz3 .lt. 0.0) THEN  ! keep commented in v12C
!               graz3 = 0
!         END IF
! END IF

 ! To account for paralelization
graz3 = graz3/8.0
IF ((tracerpp(kwq,lwq,LALG3)- graz3) .le. 0.01) THEN
    graz3 = 0.0
END IF
sourcesink(kwq,lwq,LALG3) = sourcesink(kwq,lwq,LALG3) - graz3

!! Estimate the minimum alloable nutrient value and recalculate grazing if nutrient < 0
!! minimum_nut3 = 0.0
!! minimum_nut3 = MIN(sourcesink(kwq,lwq,LALG3)-graz3,sourcesink(kwq,lwq,LALG3)*rnc + growth3*rnc - mort3*rnc - graz3*rnc, sourcesink(kwq,lwq,LALG3)*rpc + growth3*rpc - mort3*rpc - graz3*rpc)
!! IF (minimum_nut3 .lt. 0.0) THEN
!!    graz3 = 0.9*(MIN(sourcesink(kwq,lwq,LALG3), (sourcesink(kwq,lwq,LALG3)*rnc + growth3*rnc - mort3*rnc)/rnc,(sourcesink(kwq,lwq,LALG3)*rpc + growth3*rpc - mort3*rpc)/rpc)
!!        IF (graz3 .lt. 0.0) THEN
!!             graz3 = 0.0
!!         END IF
!!END IF

!sourcesink(kwq,lwq,LALG3) = sourcesink(kwq,lwq,LALG3) - sett3 + resus3
!!sourcesink(kwq,lwq,LALG3) = sourcesink(kwq,lwq,LALG3) + growth3 - mort3  graz3 - sett3 + resus3

var3 = sourcesink(kwq,lwq,LALG3)
   IF (lwq .eq. 300) THEN
      IF (kwq .eq. 5) THEN
       !PRINT *, 'growth3 = ', growth3, ', mort3 = ', mort3, 'grazing3 = ', graz3, ', SS3 = ', var3, ', tracer3 = ',tracerpp(kwq,lwq,LALG3)
       !PRINT *, 'epischuranew = ', epischuranew, ', codotisnapulinew = ', codotisnapulinew, 'diatomusnew = ', diatomusnew
       !PRINT *, 'lwq3 = ', lwq, ', kwq3 = ', kwq, ', growth3 = ', growth3, ', SS3 = ', var3, ', tracer3 = ',tracerpp(kwq,lwq,LALG3), ', f_L3 = ', f_L3, ', f_N3 = ', f_N, ', f_P3 = ', f_P, ', !PARpenetrates3 = ',Qsw*QswFr(kwq,lwq)*0.47
      END IF
   END IF

! If dissolved oxygen is modeled, alter sourcesink(kwq,lwq,LDO) to reflect growth and resp
! of algae population
IF (IDO ==1) THEN
  sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) + roc*(growth3)
END IF

! IF PON is modeled, alter sourcesink(kwq,lwq,LPON) to reflect mortality of algae
IF (IPON == 1) THEN
  sourcesink(kwq,lwq,LPON) = sourcesink(kwq,lwq,LPON) + rnc*mort3
END IF

! If DON is modeled, alter soucesink(kwq,lwq,LDON) to reflect excretion of algae (used motality as a sink)
IF (IDON ==1) THEN
  sourcesink(kwq,lwq,LDON) = sourcesink(kwq,lwq,LDON)  + rnc*mort3
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

! If DOP is modeled, alter sourcesink(kwq,lwq,LDOP) to include excretion (used motality as a sink)
IF (IDOP == 1) THEN
  sourcesink(kwq,lwq,LDOP) = sourcesink(kwq,lwq,LDOP) + rpc*mort3
END IF

! If PO4 is modeled, alter sourcesink(kwq,lwq,LPO4) to include excretion and uptake
IF (IPO4 == 1) THEN
  sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) - rpc* growth3 ! excr from phyto removed
END IF

! If POC is modeled, alter sourcesink(kwq,lwq,LPOC) to include mortality
IF (IPOC == 1) THEN
  sourcesink(kwq,lwq,LPOC) = sourcesink(kwq,lwq,LPOC) + mort3
END IF

! If DOC is modeled, alter sourcesink(kwq,lwq,LDOC) to include excretion (used motality as a sink)
IF (IDOC == 1) THEN
  sourcesink(kwq,lwq,LDOC) = sourcesink(kwq,lwq,LDOC) + mort3
END IF


END SUBROUTINE sourceALG3

!************************************************************************
SUBROUTINE sourceALG4(kwq, lwq, thrs)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on ALG4
!
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq
  REAL, INTENT(IN) :: thrs

  !. . .Local Variables
  REAL::  mu4, f_L4, f_T, f_N, f_P, N_conc
  REAL::  growth4, mort4, excr4,  graz4, sett4, resus4, varprint
  REAL::  dthrs_zoop

  ! Calculate mu, growth rate

  !. .  Calculate growth limiting factors
    ! Light Limitation - by Steele equation (Jassby and Platt, 1976)
      f_L4 = ((Qsw*QswFr(kwq,lwq)*0.47)/light_sat4) *EXP(1 -((Qsw*QswFr(kwq,lwq)*0.47)/light_sat4)) ! Steele (1962) = Photoinhibited
      !f_L4 = 1 - EXP(-(Qsw*QswFr(kwq,lwq)*0.47)/light_k4) ! Webb et al (1974) in absence of photoinhibition 

  ! temperature limitaton
    f_T = Theta_mu**(-0.01*(salp(kwq,lwq) - Topt1)**2)

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
    
    
! Growth rate considering limiting factors (light and nutrients)
mu4 = mu_max4 * MIN(f_L4,f_N,f_P)

! Read Zooplankton. Daily values between April 13 and June 15, 2018
    rotifers =  (/192.54,  196.37,  200.26,  204.19,  208.16,  212.17,  216.21,  220.29,  224.39,  228.51,  232.65,  236.80,  240.96,  245.12,  249.28,  253.44,  257.58,  &
  &    261.72,  265.83,  269.93,  273.99,  278.03,  282.03,  285.99,  289.91,  293.78,  297.60,  301.36,  305.06,  308.69,  312.26,  315.75,  319.16,  322.49,  325.73,  328.88,  &
  &    331.94,  334.89,  337.74,  340.49,  343.12,  345.63,  348.12,  350.68,  353.30,  355.97,  358.69,  361.45,  364.24,  367.05,  369.87,  372.69,  375.52,  378.33,  381.12,  383.89,  &
  &    386.61,  389.30,  391.93,  394.51,  397.01,  399.44,  401.79,  404.04/)
    rotifers = rotifers/1000

! Interpolate Zooplankton to the simulation time interval (thrs)
dthrs_zoop = 24.
rotifersnew = parabn(0.,thrs,rotifers,dthrs_zoop)

!. . Calculate growth
    growth4 = mu4 * f_T * tracerpp(kwq,lwq,LALG4)
!. . Calculate mortality
    mort4   = k_mor4 * Theta_mor**(-0.01*(salp(kwq,lwq) - Topt4)**2) * tracerpp(kwq,lwq,LALG4)
!. . Calculate excretion
    excr4   = k_ex4 * Theta_exc**(-0.01*(salp(kwq,lwq) - Topt4)**2) * tracerpp(kwq,lwq,LALG4)
!. . Calculate grazing
    !graz4   = k_gr4 * Theta_gr**(-0.01*(salp(kwq,lwq) - Topt4)**2) * tracerpp(kwq,lwq,LALG4)
     graz4 = (rotifersnew*0.0003936/4.0)*(86400.0)! G30

!. . Calculate settling
    sett4   = k_set * tracerpp(kwq,lwq,LALG4)
!. . Calculate resuspension
    resus4   = k_rs * tracerpp(kwq,lwq,LALG4)

sourcesink(kwq,lwq,LALG4) = sourcesink(kwq,lwq,LALG4) + growth4 - excr4 - mort4  - graz4 - sett4 + resus4



! If dissolved oxygen is modeled, alter sourcesink(kwq,lwq,LDO) to reflect growth and resp
! of algae population
IF (IDO ==1) THEN
  sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) + roc*(growth4)
END IF

! IF PON is modeled, alter sourcesink(kwq,lwq,LPON) to reflect mortality of algae
IF (IPON == 1) THEN
  sourcesink(kwq,lwq,LPON) = sourcesink(kwq,lwq,LPON) + rnc*mort4
END IF

! If DON is modeled, alter soucesink(kwq,lwq,LDON) to reflect excretion of algae
IF (IDON ==1) THEN
  sourcesink(kwq,lwq,LDON) = sourcesink(kwq,lwq,LDON)  + rnc*excr4
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

! If DOP is modeled, alter sourcesink(kwq,lwq,LDOP) to include excretion
IF (IDOP == 1) THEN
  sourcesink(kwq,lwq,LDOP) = sourcesink(kwq,lwq,LDOP) + rpc*excr4
END IF

! If PO4 is modeled, alter sourcesink(kwq,lwq,LPO4) to include excretion and uptake
IF (IPO4 == 1) THEN
  sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) - rpc* growth4 ! excr from phyto removed
END IF

! If POC is modeled, alter sourcesink(kwq,lwq,LPOC) to include mortality
IF (IPOC == 1) THEN
  sourcesink(kwq,lwq,LPOC) = sourcesink(kwq,lwq,LPOC) + mort4
END IF

! If DOC is modeled, alter sourcesink(kwq,lwq,LDOC) to include excretion
IF (IDOC == 1) THEN
  sourcesink(kwq,lwq,LDOC) = sourcesink(kwq,lwq,LDOC) + excr4
END IF


END SUBROUTINE sourceALG4

!************************************************************************
SUBROUTINE sourceALG5(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on ALG5
!
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq

  !. . .Local Variables
  REAL::  mu5, f_L5, f_T, f_N, f_P, N_conc
  REAL::  growth5, mort5, excr5,  graz5, sett5, resus5, varprint

  ! Calculate mu, growth rate

  !. .  Calculate growth limiting factors
    ! Light Limitation - by Steele equation (Jassby and Platt, 1976)
      f_L5 = ((Qsw*QswFr(kwq,lwq)*0.47)/light_sat5) *EXP(1 -((Qsw*QswFr(kwq,lwq)*0.47)/light_sat5)) ! Steele (1962) = Photoinhibited
      !f_L5 = 1 - EXP(-(Qsw*QswFr(kwq,lwq)*0.47)/light_k5) ! Webb et al (1974) in absence of photoinhibition 

  ! temperature limitaton
    f_T = Theta_mu**(-0.01*(salp(kwq,lwq) - Topt1)**2)

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
    
    
! Growth rate considering limiting factors (light and nutrients)
mu5 = mu_max5 * MIN(f_L5,f_N,f_P)

!. . Calculate growth
    growth5 = mu5 * f_T * tracerpp(kwq,lwq,LALG5)
!. . Calculate mortality
    mort5   = k_mor5 * Theta_mor**(-0.01*(salp(kwq,lwq) - Topt5)**2) * tracerpp(kwq,lwq,LALG5)
!. . Calculate excretion
    excr5   = k_ex5 * Theta_exc**(-0.01*(salp(kwq,lwq) - Topt5)**2) * tracerpp(kwq,lwq,LALG5)
!. . Calculate grazing
    graz5   = k_gr5 * Theta_gr**(-0.01*(salp(kwq,lwq) - Topt5)**2) * tracerpp(kwq,lwq,LALG5)
!. . Calculate settling
    sett5   = k_set * tracerpp(kwq,lwq,LALG5)
!. . Calculate resuspension
    resus5   = k_rs * tracerpp(kwq,lwq,LALG5)

sourcesink(kwq,lwq,LALG5) = sourcesink(kwq,lwq,LALG5) + growth5 - excr5 - mort5  - graz5 - sett5 + resus5


! If dissolved oxygen is modeled, alter sourcesink(kwq,lwq,LDO) to reflect growth and resp
! of algae population
IF (IDO ==1) THEN
  sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) + roc*(growth5)
END IF

! IF PON is modeled, alter sourcesink(kwq,lwq,LPON) to reflect mortality of algae
IF (IPON == 1) THEN
  sourcesink(kwq,lwq,LPON) = sourcesink(kwq,lwq,LPON) + rnc*mort5
END IF

! If DON is modeled, alter soucesink(kwq,lwq,LDON) to reflect excretion of algae
IF (IDON ==1) THEN
  sourcesink(kwq,lwq,LDON) = sourcesink(kwq,lwq,LDON)  + rnc*excr5
END IF

! IF NH4 is modeled, alter sourcesink(kwq,lwq,LNH4) to include uptake of NH4 by algae
IF (INH4 == 1) THEN
  sourcesink(kwq,lwq,LNH4) = sourcesink(kwq,lwq,LNH4) - rnc*FNH4*growth5
END IF

! IF NO3 is modeled, alter sourcesink(kwq,lwq,NO3) to include uptake of NO3 by algae
IF (INO3 == 1) THEN
  sourcesink(kwq,lwq,LNO3) = sourcesink(kwq,lwq,LNO3) - rnc*(1-FNH4) * growth5
END IF

! If POP is modeled, alter sourcesink(kwq,lwq,POP) to include mortality of algae
IF (IPOP == 1) THEN
  sourcesink(kwq,lwq,LPOP) = sourcesink(kwq,lwq,LPOP) + rpc*mort5
END IF

! If DOP is modeled, alter sourcesink(kwq,lwq,LDOP) to include excretion
IF (IDOP == 1) THEN
  sourcesink(kwq,lwq,LDOP) = sourcesink(kwq,lwq,LDOP) + rpc*excr5
END IF

! If PO4 is modeled, alter sourcesink(kwq,lwq,LPO4) to include excretion and uptake
IF (IPO4 == 1) THEN
  sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) - rpc* growth5 ! excr from phyto removed
END IF

! If POC is modeled, alter sourcesink(kwq,lwq,LPOC) to include mortality
IF (IPOC == 1) THEN
  sourcesink(kwq,lwq,LPOC) = sourcesink(kwq,lwq,LPOC) + mort5
END IF

! If DOC is modeled, alter sourcesink(kwq,lwq,LDOC) to include excretion
IF (IDOC == 1) THEN
  sourcesink(kwq,lwq,LDOC) = sourcesink(kwq,lwq,LDOC) + excr5
END IF


END SUBROUTINE sourceALG5

!************************************************************************
SUBROUTINE sourceFT(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate FT
!
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq	
	! temperature limitaton
		
		tracerpp(kwq,lwq,LFT) = Theta_mu**(-0.01*(salp(kwq,lwq) - Topt3)**2)

END SUBROUTINE sourceFT


!************************************************************************
SUBROUTINE sourceFN(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate FN
!
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq	
  REAL::	N_conc
		
		
	! nutrient limitation - but only if the nutrients are modeled
	IF ((INH4 ==1) .AND. (INO3 ==1)) THEN
		N_conc = tracerpp(kwq,lwq,LNH4) + tracerpp(kwq,lwq,LNO3)
		tracerpp(kwq,lwq,LFN)  = N_conc/(KSN + N_conc)

	ELSE IF (INH4 == 1 .AND. INO3 ==0) THEN
		N_conc = tracerpp(kwq,lwq,LNH4)
		tracerpp(kwq,lwq,LFN) = N_conc/(KSN + N_conc)

	ELSE IF (INH4 == 0 .AND. INO3 ==1) THEN
		N_conc = tracerpp(kwq,lwq,LNO3)
		tracerpp(kwq,lwq,LFN) = N_conc/(KSN + N_conc)

	ELSE
		tracerpp(kwq,lwq,LFN) = 1.0
		
	END IF

END SUBROUTINE sourceFN

!************************************************************************
SUBROUTINE sourceFP(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate FP
!
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq	

		
		
	! nutrient limitation - but only if the nutrients are modeled
	IF (IPO4 ==1) THEN
		tracerpp(kwq,lwq,LFP) = tracerpp(kwq,lwq,LPO4) /(KSP + tracerpp(kwq,lwq,LPO4) )
		
	ELSE
		tracerpp(kwq,lwq,LFP) = 1.0
		
	END IF

END SUBROUTINE sourceFP

!************************************************************************
SUBROUTINE sourceFL1(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate FL1
!
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq	

	tracerpp(kwq,lwq,LFL1) = ((Qsw*QswFr(kwq,lwq)*0.47)/light_sat1) *EXP(1 -((Qsw*QswFr(kwq,lwq)*0.47)/light_sat1)) ! Steele (1962) = Photoinhibited	
	
END SUBROUTINE sourceFL1

!************************************************************************
SUBROUTINE sourceFL2(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate FL2
!
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq	

	tracerpp(kwq,lwq,LFL2) = ((Qsw*QswFr(kwq,lwq)*0.47)/light_sat2) *EXP(1 -((Qsw*QswFr(kwq,lwq)*0.47)/light_sat2)) ! Steele (1962) = Photoinhibited	
		
END SUBROUTINE sourceFL2

!************************************************************************
SUBROUTINE sourceFL3(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate FL3
!
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq	

	tracerpp(kwq,lwq,LFL3) = ((Qsw*QswFr(kwq,lwq)*0.47)/light_sat3) *EXP(1 -((Qsw*QswFr(kwq,lwq)*0.47)/light_sat3)) ! Steele (1962) = Photoinhibited	
		
END SUBROUTINE sourceFL3

!************************************************************************
SUBROUTINE sourceFL4(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate FL4
!
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq	

	tracerpp(kwq,lwq,LFL4) = ((Qsw*QswFr(kwq,lwq)*0.47)/light_sat4) *EXP(1 -((Qsw*QswFr(kwq,lwq)*0.47)/light_sat4)) ! Steele (1962) = Photoinhibited	
		
END SUBROUTINE sourceFL4

!************************************************************************
SUBROUTINE sourceFL5(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate FL5
!
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq	

	tracerpp(kwq,lwq,LFL5) = ((Qsw*QswFr(kwq,lwq)*0.47)/light_sat5) *EXP(1 -((Qsw*QswFr(kwq,lwq)*0.47)/light_sat5)) ! Steele (1962) = Photoinhibited	
		
END SUBROUTINE sourceFL5

!***********************************************************************
FUNCTION parabn ( frstpt, x, fx, dx )
!***********************************************************************
!
!  Purpose: To interpolate parabolically between the functional values
!           within the array fx.
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
