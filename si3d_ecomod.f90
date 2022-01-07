!************************************************************************
                          MODULE si3d_ecomod
!************************************************************************
!
!  Purpose: Procedures that implement routines related to the modelling
!           of ecological processees
!
!-------------------------------------------------------------------------

   !USE si3d_ecomod
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
  &    iDOC, iPOC, &
  &    iALG1, iALG2, iALG3, iALG4, iALG5, &
  &    iZOO
  IF (ios /= 0) CALL input_error ( ios, 92)

  !. . . Read model stochiometeric constants and other constants
  READ (UNIT=i99,FMT='(///(14X,G20.3))',IOSTAT=ios) rnc, rpc, roc, ron, &
  &     KNIT, KSN, KSP, FNH4, KDOC, &
  &     light_sat1, light_sat2, light_sat3, light_sat4, light_sat5, BacteriaC
  IF (ios /= 0) CALL input_error ( ios, 93)

  !. . . Read model rates
  READ (UNIT=i99,FMT='(///(14X,G20.2))',IOSTAT=ios) k_a,  &
  &    mu_max1, k_mor1, k_ex1, k_res1, k_gr1, &
  &    mu_max2, k_mor2, k_ex2, k_res2, k_gr2, &
  &    mu_max3, k_mor3, k_ex3, k_res3, k_gr3, &
  &    mu_max4, k_mor4, k_ex4, k_res4, k_gr4, &
  &    mu_max5, k_mor5, k_ex5, k_res5, k_gr5, &
  &    k_dcn, k_mn, k_n, k_dn, &
  &    k_dcp, k_mp, k_dcc, &
  &    k_set, k_rs, &
  &    k_morz, k_exz, k_resz, k_grdet, k_grbac
  IF (ios /= 0) CALL input_error ( ios, 94)


  !... Convert model rates to time of [1/timestep]. Input file has 1/day values
  ! DO
  k_a   = idt * k_a/86400.0
  ! ALG1 : Pycoplankton
  mu_max1 = idt * mu_max1/86400.0
  k_mor1 = idt * k_mor1/86400.0
  k_ex1 = idt * k_ex1/86400.0
  k_res1 = idt * k_res1/86400.0
  k_gr1 = idt * k_gr1/86400.0
  ! ALG2 : Cyclotella
  mu_max2 = idt * mu_max2/86400.0
  k_mor2 = idt * k_mor2/86400.0
  k_ex2 = idt * k_ex2/86400.0
  k_res2 = idt * k_res2/86400.0
  k_gr2 = idt * k_gr2/86400.0
  ! ALG3: Cryptomonas
  mu_max3 = idt * mu_max3/86400.0
  k_mor3 = idt * k_mor3/86400.0
  k_ex3 = idt * k_ex3/86400.0
  k_res3 = idt * k_res3/86400.0
  k_gr3 = idt * k_gr3/86400.0
  ! ALG4: Synedra
  mu_max4 = idt * mu_max4/86400.0
  k_mor4 = idt * k_mor4/86400.0
  k_ex4 = idt * k_ex4/86400.0
  k_res4 = idt * k_res4/86400.0
  k_gr4 = idt * k_gr4/86400.0
  ! ALG5: Microcystis
  mu_max5 = idt * mu_max5/86400.0
  k_mor5 = idt * k_mor5/86400.0
  k_ex5 = idt * k_ex5/86400.0
  k_res5 = idt * k_res5/86400.0
  k_gr5 = idt * k_gr5/86400.0
  ! Nutrients
  k_dcn = idt * k_dcn/86400.0
  k_mn = idt * k_mn/86400.0
  k_n = idt * k_n/86400.0
  k_dn = idt * k_dn/86400.0
  k_dcp = idt * k_dcp/86400.0
  k_mp = idt * k_mp/86400.0
  k_set = idt * k_set/86400.0
  k_rs = idt * k_rs/86400.0
  k_dcc = idt * k_dcc/86400.0
  ! Zooplankton
  k_morz = idt * k_morz/86400.0
  k_exz = idt * k_exz/86400.0
  k_resz = idt * k_resz/86400.0
  k_grdet = idt * k_grdet/86400.0
  k_grbac = idt * k_grbac/86400.0

  !. . . Read model temperature rates
  READ (UNIT=i99,FMT='(///(14X,G20.2))',IOSTAT=ios) Theta_a, Theta_mu, Theta_mor, Theta_res, Theta_gr, &
  &     Theta_dcn, Theta_mn, Theta_n , Theta_dn , &
  &     Theta_dcp , Theta_mp , Theta_dcc , Theta_DOC, Theta_morz, Theta_resz

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
    PRINT*, "iZOO = ", iZOO
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
  LZOO=0;

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

  IF (iALG1 == 1) THEN
    tracerpplocal(i) = 9
    i = i+1
  END IF

  IF (iALG2 == 1) THEN
    tracerpplocal(i) = 10
    i = i+1
  END IF

  IF (iALG3 == 1) THEN
    tracerpplocal(i) = 11
    i = i+1
  END IF

  IF (iALG4 == 1) THEN
    tracerpplocal(i) = 12
    i = i+1
  END IF

  IF (iALG5 == 1) THEN
    tracerpplocal(i) = 13
    i = i+1
  END IF

  IF (iDOC == 1) THEN
    tracerpplocal(i) = 14
    i = i+1
  END IF

  IF (iPOC == 1) THEN
    tracerpplocal(i) = 15
    i = i+1
  END IF

  IF (iZOO == 1) THEN
    tracerpplocal(i) = 16
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
		LALG1 = i
  ELSEIF (tracerpplocal(i) == 10) THEN
    LALG2 = i
  ELSEIF (tracerpplocal(i) == 11) THEN
    LALG3 = i
  ELSEIF (tracerpplocal(i) == 12) THEN
 		LALG4 = i
  ELSEIF (tracerpplocal(i) == 13) THEN
    LALG5 = i
	ELSEIF (tracerpplocal(i) == 14) THEN
		LDOC = i
	ELSEIF (tracerpplocal(i) == 15) THEN
		LPOC = i
  ELSEIF (tracerpplocal(i) == 16) THEN
		LZOO = i
	END IF
  END DO

  IF (idbg == 1) THEN
    PRINT*, "LDO  = ", LDO, "LPOC = ", LPOC, "LDOC = ", LDOC
    PRINT*, "LPON = ", LPON, "LDON = ", LDON, "LNH4 = ", LNH4, "LNO3 = ", LNO3
    PRINT*, "LPOP = ", LPOP, "LDOP = ", LDOP
    PRINT*, "LALG1 = ", LALG1, "LALG2 = ", LALG2, "LALG3 = ", LALG3, "LALG4 = ", LALG4, "LALG5 = ", LALG5
    PRINT*, "LZOO = ", LZOO
  END IF

END SUBROUTINE WQinit

!************************************************************
SUBROUTINE srcsnkWQ
!***********************************************************
!
!   Purpose: to call all source subroutines for all cells
!
!------------------------------------------------------------


  !... Local variables
  INTEGER:: i,j,k,l, k1s, kms,itr

  ! reset soursesink = 0
  !sourcesink = 0;

  DO l = 1, lm;

    ! ... Map l- into (i,j)-indexes .........................
    !i = l2i(l); j = l2j(l);

    ! ... Retrieve top & bottom wet sal-pts .................
    kms = kmz(j)
    k1s = k1z(j)

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
      IF (iDOC == 1) THEN
        CALL sourceDOC(k,l)
      END IF
      IF (iPOC == 1) THEN
        CALL sourcePOC(k,l)
      END IF
      IF (iZOO == 1) THEN
        CALL sourceZOO(k,l)
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
  REAL		::	 Tk, lnOS, OS, ln_Pwv, Pwv, theta2, Patm

  ! Calculate DO saturation
  Tk = salp(kwq,lwq) + 273
  lnos = -139.34410 + 1.575701*1E5 /(Tk    ) &
  &                 - 6.642308*1E7 /(Tk**2.) &
  &	                + 1.243800*1E10/(Tk**3.) &
  &                 - 8.621949*1E11/(Tk**4.)
  os = EXP(lnos)

  ! Correct for Patmospheric (Pa - declared in si3d_types and defined in surfbc0)
  Patm   = Pa * 0.00000986923; ! Transform atmospheric pressure from Pa to atm
  ln_Pwv = 11.8751 - 3840.70/Tk - 216961/Tk
  Pwv    = EXP(ln_Pwv)
  theta2 = 0.000975 - 1.426*1E-5 * salp(kwq,lwq) + &
                      6.436*1E-8 * salp(kwq,lwq)**2
  os = os*Pa*((1-Pwv/Pa) *(1-theta2*Pa))&
  &           /((1-Pwv)*(1-theta2) )

   os = 10

  ! Calculate reaertaion
  ! for now using constant rearation defined in wq_inp, but in future, can have
  ! alternatives for calcualting reaeration.

  IF (kwq .eq. k1z(ij2l(l2i(lwq),l2j(lwq)))) THEN
  sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO)              &
          &       +    k_a*(OS - tracerpp(kwq,lwq,LDO)) 		!reaeration
					     ! + photosynthesis	- only if IALG = 1; calculated in sourceALG
					     ! - Respiration		- only if IALG = 1; calculated in sourceALG
               ! - RespirationZ		- only if IZOO = 1; calculated in sourceZOO
					     ! - Nitrification	- only if INH4 = 1; calculated in sourceNH4

 END IF

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
  REAL:: decomposition

  !. Calculate hydrolysis
  decomposition = k_dcn *Theta_dcn**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LPON)

  sourcesink(kwq,lwq,LPON) = sourcesink(kwq,lwq,LPON)               &
  &                          - decomposition			                  &   ! decomposition
  &                          - k_set* tracerpp(kwq,lwq,LPON)	      &	! settling
  &                          + k_rs * tracerpp(kwq,lwq,LPON)			! resusupension
						                ! + mortality	- only if IALG = 1; calcualted in sourceALG
                            ! + mortalityZ	- only if IZOO = 1; calcualted in sourceZOO

  ! Add contribution of mineralization to DON concentration
  IF (iDON == 1) THEN
    sourcesink(kwq,lwq,LDON) = sourcesink(kwq,lwq,LDON)	+  decomposition
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
  REAL:: minrl

  !. . mineralization
  minrl = k_mn * Theta_mn**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LDON)

  sourcesink(kwq,lwq,LDON) = sourcesink(kwq,lwq,LDON)					&
  &                         -  minrl
							! + decomposition	- only if IPON = 1; caluclated in sourcePON
							! + Excr	    - only if IAlG = 1; calculated in sourceALG

!. . .Add contribution from atmospheric deposition to top layer
IF (kwq .eq. k1z(ij2l(l2i(lwq),l2j(lwq)))) THEN
	sourcesink(kwq,lwq, LDON) = sourcesink(kwq, lwq, LDON) + &
								& (hpp(kwq,lwq)+2)*ATM_DON		! Atmoshperic deposition
END IF

  IF (INH4 == 1) THEN		! add mineralization of DON to NH4
    sourcesink(kwq,lwq,LNH4) = sourcesink(kwq,lwq,LNH4) + minrl
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

	nitrif = k_n*Theta_n**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LNH4) * f_DO

sourcesink(kwq,lwq,LNH4) = sourcesink(kwq,lwq,LNH4)			&
						&  -  nitrif			          ! nitrificatin

						! - Algal Uptake		- if IALG = 1; calculated in sourceALG
						! + mineralization		- if IDON = 1; calculated in sourceDON
            ! + excretion         - if IZOO = 1; calculated in source ZOO

!. . Add contribution from sediment release and GW flux into bottom cells
IF (kwq .eq. kmz(ij2l(l2i(lwq),l2j(lwq)))) THEN
 sourcesink(kwq,lwq,LNH4) = sourcesink(kwq,lwq,LNH4) +   &
							& (hpp(kwq,lwq)-1)*J_NH4     +	 &	! sediment release
							& (hpp(kwq,lwq)-1)*GW_NH4	     	! GW flux
END IF

!. . Add contribution from atmoshperic deposition into top cells
IF (kwq .eq. k1z(ij2l(l2i(lwq),l2j(lwq)))) THEN
	sourcesink(kwq, lwq, LNH4) = sourcesink(kwq, lwq, LNH4) + &
								& (hpp(kwq, lwq)+2) * ATM_NH4	 ! ATM deposition
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
						&  k_dn*Theta_dn**(salp(kwq,lwq) - 20)* tracerpp(kwq,lwq,LNO3)	! denitrification
						! + nitrification		- if INH4 = 1; calculated in sourceNH4
						! - algal uptake		- if IALG = 1; calculated in sourceALG

!. . Add contribution from sediment release and GW flux into bottom cells
IF (kwq .eq. kmz(ij2l(l2i(lwq),l2j(lwq)))) THEN
 sourcesink(kwq,lwq,LNO3) = sourcesink(kwq,lwq,LNO3) +   &
							& (hpp(kwq,lwq)-1)*J_NO3     +	 &	! sediment release
							& (hpp(kwq,lwq)-1)*GW_NO3	     	! GW flux
END IF

!. . Add contribution from atmoshperic deposition into top cells
IF (kwq .eq. k1z(ij2l(l2i(lwq),l2j(lwq)))) THEN
	sourcesink(kwq, lwq, LNO3) = sourcesink(kwq, lwq, LNO3) + &
								& (hpp(kwq, lwq)+2) * ATM_NO3	 ! ATM deposition
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

decomposition = k_dcp*Theta_dcp**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LPOP)

sourcesink(kwq,lwq,LPOP) = sourcesink(kwq,lwq,LPOP)		    &
						&   - decomposition                           & ! decomposition
						&   - k_set * tracerpp(kwq, lwq, LPOP)        &  ! settling
						&   + k_rs*tracerpp(kwq,lwq,LPOP)			          ! Resuspension
						! + mortality	- only if IALG = 1; calcualted in sourceALG
            ! + mortalityZ	- only if IZOO = 1; calcualted in sourceZOO

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

minrl = k_mp*Theta_mp**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LDOP)

sourcesink(kwq,lwq,LDOP) = sourcesink(kwq,lwq,LDOP)		&
						&   -  minrl				                  ! decomposition
						! + excretion	      - if IAlG = 1; calculated in sourceALG
						! + decomposition		- if iPOP = 1; calculated in sourcePOP

!... Add atmoshperic deposition to top layer
IF (kwq == k1z(ij2l(l2i(lwq),l2j(lwq)))) THEN
	sourcesink(kwq,lwq,LDOP) = sourcesink(kwq,lwq,LDOP)				+	&
						&  (h(kwq,lwq)+2)*ATM_DOP					! Atmospheric deposition
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
   ! + zoo exc			- if IZOO = 1; calculated in sourceZOO

!. . Add contribution from sediment release and GW flux into bottom cells
IF (kwq .eq. kmz(ij2l(l2i(lwq),l2j(lwq)))) THEN
 sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) +   &
							& (hpp(kwq,lwq)-1)*J_PO4     +	 &	! sediment release
							& (hpp(kwq,lwq)-1)*GW_PO4	     	! GW flux
END IF

!. . Add contribution from atmoshperic deposition into top cells
IF (kwq .eq. k1z(ij2l(l2i(lwq),l2j(lwq)))) THEN
	sourcesink(kwq, lwq, LPO4) = sourcesink(kwq, lwq, LPO4) + &
								& (hpp(kwq, lwq)+2) * ATM_PO4	 ! ATM deposition
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
  REAL:: decomposition, grazingdet

  ! Calculate decomposition
  decomposition = k_dcc*Theta_dcc**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LPOC)
  ! Calculate grazing detritus
  grazingdet = k_grdet * Theta_gr**(salp(kwq,lwq) - 20) *tracerpp(kwq,lwq,LPOC)

  sourcesink(kwq,lwq,LPOC) = sourcesink(kwq,lwq,LPOC)	&
						&  -   decomposition						          &              ! decomposition
            &  -   grazingdet                         &              ! grazing of detritus
						&  -   k_set * tracerpp(kwq,lwq,LPOC) 	  &	             ! settling
						&  +   k_rs * tracerpp(kwq,lwq,LPOC)			                 ! Resuspension
						! + mortality		                        - if IALG = 1; calculated in sourceALG

  ! Caclulate decomposition contribution to DOC
  IF (IDOC == 1) THEN
    sourcesink(kwq,lwq,LDOC) = sourcesink(kwq,lwq,LDOC) + decomposition
   END IF

   ! Calculate ZOO change due to detritus grazing
   IF (IZOO == 1) THEN
     sourcesink(kwq,lwq,LZOO) = sourcesink(kwq,lwq,LZOO) + grazingdet
   END IF

   ! Calculate PON change due to detritus grazing
   IF (IPON == 1) THEN
     sourcesink(kwq,lwq,LPON) = sourcesink(kwq,lwq,LPON) - rnc * grazingdet
   END IF

   ! Calculate POP change due to detritus grazing
   IF (IPOP == 1) THEN
     sourcesink(kwq,lwq,LPOP) = sourcesink(kwq,lwq,LPOP)  - rpc * grazingdet
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
  oxid = KDOC * Theta_DOC**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LDOC) * F_DO

  sourcesink(kwq,lwq,LDOC) = sourcesink(kwq,lwq,LDOC)  &
                     &    - oxid                         ! microbial uptake (oxidation)
						! + decomposition	                -if IPOC = 1; calculated in sourcePOC
						! + algal excretion             - if IALG = 1; calculated in sourceALG


!. . Add contribution from sediment release and GW flux into bottom cells
IF (kwq .eq. kmz(ij2l(l2i(lwq),l2j(lwq)))) THEN
 sourcesink(kwq,lwq,LDOC) = sourcesink(kwq,lwq,LDOC) +   &
							& (hpp(kwq,lwq)-1)*J_DOC     +	 &	! sediment release
							& (hpp(kwq,lwq)-1)*GW_DOC	     	! GW flux
END IF

!. . Add contribution from atmoshperic deposition into top cells
IF (kwq .eq. k1z(ij2l(l2i(lwq),l2j(lwq)))) THEN
	sourcesink(kwq, lwq, LDOC) = sourcesink(kwq, lwq, LDOC) + &
								& (hpp(kwq, lwq)+2) * ATM_DOC	 ! ATM deposition
END IF

  IF (IDO == 1) THEN
    sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) - roc*oxid
  END IF

END SUBROUTINE sourceDOC

!************************************************************************
SUBROUTINE sourceALG1(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on ALG1
!
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq

  !. . .Local Variables
  REAL::	mu1, f_L1, f_T, f_N, f_P, N_conc
  REAL::	growth1, resp1, excr1, mort1, graz1, sett1, resus1

  ! Calculate mu, growth rate

  !. .  Calculate growth limiting factors
    ! Light Limitation - by Steele equation (Jassby and Platt, 1976)
   		f_L1 = ((Qsw*QswFr(kwq,lwq)*0.47)/light_sat1) *EXP(1 -((Qsw*QswFr(kwq,lwq)*0.47)/light_sat1))
		IF (f_L1 == 0) THEN
		   f_L1 = 1
		END IF

	! temperature limitaton
		f_T = Theta_mu**(salp(kwq,lwq) -20)

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

!. . Calculate growth
		growth1 = mu1 * f_T * tracerpp(kwq,lwq,LALG1)
!. . Calculate respiration
		resp1   = k_res1 * Theta_res**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LALG1)
!. . Calculate excretion
		excr1   = k_ex1 * Theta_res**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LALG1)
!. . Calculate mortality
		mort1   = k_mor1 * Theta_mor**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LALG1)
!. . Calculate grazing
		graz1   = k_gr1 * Theta_gr**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LALG1)
!. . Calculate settling
		sett1   = k_set * tracerpp(kwq,lwq,LALG1)
!. . Calculate resuspension
    resus1   = k_rs * tracerpp(kwq,lwq,LALG1)

sourcesink(kwq,lwq,LALG1) = sourcesink(kwq,lwq,LALG1)	+ growth1 - resp1	- excr1	- mort1  - graz1 - sett1 + resus1

! If dissolved oxygen is modeled, alter sourcesink(kwq,lwq,LDO) to reflect growth and resp
! of algae population
IF (IDO ==1) THEN
	sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) + roc*(growth1-resp1)
END IF

! IF PON is modeled, alter sourcesink(kwq,lwq,LPON) to reflect mortality of algae
IF (IPON == 1) THEN
	sourcesink(kwq,lwq,LPON) = sourcesink(kwq,lwq,LPON) + rnc*mort1
END IF

! If DON is modeled, alter soucesink(kwq,lwq,LDON) to reflect excretion of algae
IF (IDON ==1) THEN
	sourcesink(kwq,lwq,LDON) = sourcesink(kwq,lwq,LDON)  + rnc*excr1
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

! If DOP is modeled, alter sourcesink(kwq,lwq,LDOP) to include excretion
IF (IDOP == 1) THEN
	sourcesink(kwq,lwq,LDOP) = sourcesink(kwq,lwq,LDOP) + rpc*excr1
END IF

! If PO4 is modeled, alter sourcesink(kwq,lwq,LPO4) to include excretion and uptake
IF (IPO4 == 1) THEN
	sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) + rpc*(excr1 - growth1)
END IF

! If POC is modeled, alter sourcesink(kwq,lwq,LPOC) to include mortality
IF (IPOC == 1) THEN
	sourcesink(kwq,lwq,LPOC) = sourcesink(kwq,lwq,LPOC) + mort1
END IF

! If DOC is modeled, alter sourcesink(kwq,lwq,LDOC) to include excretion
IF (IDOC == 1) THEN
	sourcesink(kwq,lwq,LDOC) = sourcesink(kwq,lwq,LDOC) + excr1
END IF

! If ZOO is modeled, alter sourcesink(kwq,lwq,LZOO) to include grazing
IF (IZOO == 1) THEN
	sourcesink(kwq,lwq,LZOO) = sourcesink(kwq,lwq,LZOO) + graz1
END IF

END SUBROUTINE sourceALG1

!************************************************************************
SUBROUTINE sourceALG2(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on ALG2
!
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq

  !. . .Local Variables
  REAL::	mu2, f_L2, f_T, f_N, f_P, N_conc
  REAL::	growth2, resp2, excr2, mort2, graz2, sett2, resus2

  ! Calculate mu, growth rate

  !. .  Calculate growth limiting factors
    ! Light Limitation - by Steele equation (Jassby and Platt, 1976)
   		f_L2 = ((Qsw*QswFr(kwq,lwq)*0.47)/light_sat2) *EXP(1 -((Qsw*QswFr(kwq,lwq)*0.47)/light_sat2))
		IF (f_L2 == 0) THEN
		   f_L2 = 1
		END IF

	! temperature limitaton
		f_T = Theta_mu**(salp(kwq,lwq) -20)

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

!. . Calculate growth
		growth2 = mu2 * f_T * tracerpp(kwq,lwq,LALG2)
!. . Calculate respiration
		resp2   = k_res2 * Theta_res**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LALG2)
!. . Calculate excretion
		excr2   = k_ex2 * Theta_res**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LALG2)
!. . Calculate mortality
		mort2   = k_mor2 * Theta_mor**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LALG2)
!. . Calculate grazing
		graz2   = k_gr2 * Theta_gr**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LALG2)
!. . Calculate settling
		sett2   = k_set * tracerpp(kwq,lwq,LALG2)
!. . Calculate resuspension
    resus2   = k_rs * tracerpp(kwq,lwq,LALG2)

sourcesink(kwq,lwq,LALG2) = sourcesink(kwq,lwq,LALG2)	+ growth2 - resp2	- excr2	- mort2  - graz2 - sett2 + resus2

! If dissolved oxygen is modeled, alter sourcesink(kwq,lwq,LDO) to reflect growth and resp
! of algae population
IF (IDO ==1) THEN
	sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) + roc*(growth2-resp2)
END IF

! IF PON is modeled, alter sourcesink(kwq,lwq,LPON) to reflect mortality of algae
IF (IPON == 1) THEN
	sourcesink(kwq,lwq,LPON) = sourcesink(kwq,lwq,LPON) + rnc*mort2
END IF

! If DON is modeled, alter soucesink(kwq,lwq,LDON) to reflect excretion of algae
IF (IDON ==1) THEN
	sourcesink(kwq,lwq,LDON) = sourcesink(kwq,lwq,LDON)  + rnc*excr2
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

! If DOP is modeled, alter sourcesink(kwq,lwq,LDOP) to include excretion
IF (IDOP == 1) THEN
	sourcesink(kwq,lwq,LDOP) = sourcesink(kwq,lwq,LDOP) + rpc*excr2
END IF

! If PO4 is modeled, alter sourcesink(kwq,lwq,LPO4) to include excretion and uptake
IF (IPO4 == 1) THEN
	sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) + rpc*(excr2 - growth2)
END IF

! If POC is modeled, alter sourcesink(kwq,lwq,LPOC) to include mortality
IF (IPOC == 1) THEN
	sourcesink(kwq,lwq,LPOC) = sourcesink(kwq,lwq,LPOC) + mort2
END IF

! If DOC is modeled, alter sourcesink(kwq,lwq,LDOC) to include excretion
IF (IDOC == 1) THEN
	sourcesink(kwq,lwq,LDOC) = sourcesink(kwq,lwq,LDOC) + excr2
END IF

! If ZOO is modeled, alter sourcesink(kwq,lwq,LZOO) to include grazing
IF (IZOO == 1) THEN
	sourcesink(kwq,lwq,LZOO) = sourcesink(kwq,lwq,LZOO) + graz2
END IF

END SUBROUTINE sourceALG2

!************************************************************************
SUBROUTINE sourceALG3(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on ALG3
!
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq

  !. . .Local Variables
  REAL::	mu3, f_L3, f_T, f_N, f_P, N_conc
  REAL::	growth3, resp3, excr3, mort3, graz3, sett3, resus3

  ! Calculate mu, growth rate

  !. .  Calculate growth limiting factors
    ! Light Limitation - by Steele equation (Jassby and Platt, 1976)
   		f_L3 = ((Qsw*QswFr(kwq,lwq)*0.47)/light_sat3) *EXP(1 -((Qsw*QswFr(kwq,lwq)*0.47)/light_sat3))
		IF (f_L3 == 0) THEN
		   f_L3 = 1
		END IF

	! temperature limitaton
		f_T = Theta_mu**(salp(kwq,lwq) -20)

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

!. . Calculate growth
		growth3 = mu3 * f_T * tracerpp(kwq,lwq,LALG3)
!. . Calculate respiration
		resp3   = k_res3 * Theta_res**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LALG3)
!. . Calculate excretion
		excr3   = k_ex3 * Theta_res**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LALG3)
!. . Calculate mortality
		mort3   = k_mor3 * Theta_mor**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LALG3)
!. . Calculate grazing
		graz3   = k_gr3 * Theta_gr**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LALG3)
!. . Calculate settling
		sett3   = k_set * tracerpp(kwq,lwq,LALG3)
!. . Calculate resuspension
    resus3   = k_rs * tracerpp(kwq,lwq,LALG3)

sourcesink(kwq,lwq,LALG3) = sourcesink(kwq,lwq,LALG3)	+ growth3 - resp3	- excr3	- mort3  - graz3 - sett3 + resus3

! If dissolved oxygen is modeled, alter sourcesink(kwq,lwq,LDO) to reflect growth and resp
! of algae population
IF (IDO ==1) THEN
	sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) + roc*(growth3-resp3)
END IF

! IF PON is modeled, alter sourcesink(kwq,lwq,LPON) to reflect mortality of algae
IF (IPON == 1) THEN
	sourcesink(kwq,lwq,LPON) = sourcesink(kwq,lwq,LPON) + rnc*mort3
END IF

! If DON is modeled, alter soucesink(kwq,lwq,LDON) to reflect excretion of algae
IF (IDON ==1) THEN
	sourcesink(kwq,lwq,LDON) = sourcesink(kwq,lwq,LDON)  + rnc*excr3
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

! If DOP is modeled, alter sourcesink(kwq,lwq,LDOP) to include excretion
IF (IDOP == 1) THEN
	sourcesink(kwq,lwq,LDOP) = sourcesink(kwq,lwq,LDOP) + rpc*excr3
END IF

! If PO4 is modeled, alter sourcesink(kwq,lwq,LPO4) to include excretion and uptake
IF (IPO4 == 1) THEN
	sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) + rpc*(excr3 - growth3)
END IF

! If POC is modeled, alter sourcesink(kwq,lwq,LPOC) to include mortality
IF (IPOC == 1) THEN
	sourcesink(kwq,lwq,LPOC) = sourcesink(kwq,lwq,LPOC) + mort3
END IF

! If DOC is modeled, alter sourcesink(kwq,lwq,LDOC) to include excretion
IF (IDOC == 1) THEN
	sourcesink(kwq,lwq,LDOC) = sourcesink(kwq,lwq,LDOC) + excr3
END IF

! If ZOO is modeled, alter sourcesink(kwq,lwq,LZOO) to include grazing
IF (IZOO == 1) THEN
	sourcesink(kwq,lwq,LZOO) = sourcesink(kwq,lwq,LZOO) + graz3
END IF

END SUBROUTINE sourceALG3

!************************************************************************
SUBROUTINE sourceALG4(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on ALG4
!
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq

  !. . .Local Variables
  REAL::	mu4, f_L4, f_T, f_N, f_P, N_conc
  REAL::	growth4, resp4, excr4, mort4, graz4, sett4, resus4

  ! Calculate mu, growth rate

  !. .  Calculate growth limiting factors
    ! Light Limitation - by Steele equation (Jassby and Platt, 1976)
   		f_L4 = ((Qsw*QswFr(kwq,lwq)*0.47)/light_sat4) *EXP(1 -((Qsw*QswFr(kwq,lwq)*0.47)/light_sat4))
		IF (f_L4 == 0) THEN
		   f_L4 = 1
		END IF

	! temperature limitaton
		f_T = Theta_mu**(salp(kwq,lwq) -20)

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

!. . Calculate growth
		growth4 = mu4 * f_T * tracerpp(kwq,lwq,LALG4)
!. . Calculate respiration
		resp4   = k_res4 * Theta_res**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LALG4)
!. . Calculate excretion
		excr4   = k_ex4 * Theta_res**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LALG4)
!. . Calculate mortality
		mort4   = k_mor4 * Theta_mor**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LALG4)
!. . Calculate grazing
		graz4   = k_gr4 * Theta_gr**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LALG4)
!. . Calculate settling
		sett4   = k_set * tracerpp(kwq,lwq,LALG4)
!. . Calculate resuspension
    resus4   = k_rs * tracerpp(kwq,lwq,LALG4)

sourcesink(kwq,lwq,LALG4) = sourcesink(kwq,lwq,LALG4)	+ growth4 - resp4	- excr4	- mort4  - graz4 - sett4 + resus4

! If dissolved oxygen is modeled, alter sourcesink(kwq,lwq,LDO) to reflect growth and resp
! of algae population
IF (IDO ==1) THEN
	sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) + roc*(growth4-resp4)
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
	sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) + rpc*(excr4 - growth4)
END IF

! If POC is modeled, alter sourcesink(kwq,lwq,LPOC) to include mortality
IF (IPOC == 1) THEN
	sourcesink(kwq,lwq,LPOC) = sourcesink(kwq,lwq,LPOC) + mort4
END IF

! If DOC is modeled, alter sourcesink(kwq,lwq,LDOC) to include excretion
IF (IDOC == 1) THEN
	sourcesink(kwq,lwq,LDOC) = sourcesink(kwq,lwq,LDOC) + excr4
END IF

! If ZOO is modeled, alter sourcesink(kwq,lwq,LZOO) to include grazing
IF (IZOO == 1) THEN
	sourcesink(kwq,lwq,LZOO) = sourcesink(kwq,lwq,LZOO) + graz4
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
  REAL::	mu5, f_L5, f_T, f_N, f_P, N_conc
  REAL::	growth5, resp5, excr5, mort5, graz5, sett5, resus5

  ! Calculate mu, growth rate

  !. .  Calculate growth limiting factors
    ! Light Limitation - by Steele equation (Jassby and Platt, 1976)
   		f_L5 = ((Qsw*QswFr(kwq,lwq)*0.47)/light_sat5) *EXP(1 -((Qsw*QswFr(kwq,lwq)*0.47)/light_sat5))
		IF (f_L5 == 0) THEN
		   f_L5 = 1
		END IF

	! temperature limitaton
		f_T = Theta_mu**(salp(kwq,lwq) -20)

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
!. . Calculate respiration
		resp5   = k_res5 * Theta_res**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LALG5)
!. . Calculate excretion
		excr5   = k_ex5 * Theta_res**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LALG5)
!. . Calculate mortality
		mort5   = k_mor5 * Theta_mor**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LALG5)
!. . Calculate grazing
		graz5   = k_gr5 * Theta_gr**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LALG5)
!. . Calculate settling
		sett5   = k_set * tracerpp(kwq,lwq,LALG5)
!. . Calculate resuspension
    resus5   = k_rs * tracerpp(kwq,lwq,LALG5)

sourcesink(kwq,lwq,LALG5) = sourcesink(kwq,lwq,LALG5)	+ growth5 - resp5	- excr5	- mort5  - graz5 - sett5 + resus5

! If dissolved oxygen is modeled, alter sourcesink(kwq,lwq,LDO) to reflect growth and resp
! of algae population
IF (IDO ==1) THEN
	sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) + roc*(growth5-resp5)
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
	sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) + rpc*(excr5 - growth5)
END IF

! If POC is modeled, alter sourcesink(kwq,lwq,LPOC) to include mortality
IF (IPOC == 1) THEN
	sourcesink(kwq,lwq,LPOC) = sourcesink(kwq,lwq,LPOC) + mort5
END IF

! If DOC is modeled, alter sourcesink(kwq,lwq,LDOC) to include excretion
IF (IDOC == 1) THEN
	sourcesink(kwq,lwq,LDOC) = sourcesink(kwq,lwq,LDOC) + excr5
END IF

! If ZOO is modeled, alter sourcesink(kwq,lwq,LZOO) to include grazing
IF (IZOO == 1) THEN
	sourcesink(kwq,lwq,LZOO) = sourcesink(kwq,lwq,LZOO) + graz5
END IF

END SUBROUTINE sourceALG5

!************************************************************************
SUBROUTINE sourceZOO(kwq, lwq)
!*********************************************************************
!
! Purpose: To calculate sourcesink terms that depend on ZOO
!
!--------------------------------------------------------------------------
  ! ... Arguments
  INTEGER, INTENT (IN) :: kwq,lwq

  !. . .Local Variables
  REAL :: grazingBacteria, morz, resz, exz

  !... Calculate Grazing by Bacteria (nominal concentration)
  grazingBacteria = k_grbac * Theta_gr**(salp(kwq,lwq) - 20) * BacteriaC
  ! ... Calculate Mortality of Zooplankton
  morz = k_morz * Theta_morz **(salp(kwq,lwq) - 20) * sourcesink(kwq,lwq,LZOO)
  ! ... Calculate Respiration of Zooplankton
  resz = k_resz * Theta_resz **(salp(kwq,lwq) - 20) * sourcesink(kwq,lwq,LZOO)
  ! ... Calculate Excretion of Zooplankton
  exz = k_exz * Theta_resz **(salp(kwq,lwq) - 20) * sourcesink(kwq,lwq,LZOO)

    sourcesink(kwq,lwq,LZOO) = sourcesink(kwq,lwq,LZOO)  &
                &   + grazingBacteria                    &  ! grazing of bacteria
                &   - morz                               &  ! mortality of Zooplankton
                &   - resz                               &  ! respiration of zooplankton
                &   - exz                                   ! excretion of zooplankton
                  ! + grazing detritus	       -if IPOC = 1; calculated in sourcePOC
                  ! + grazing algae            - if IALG = 1; calculated in sourceALG


! If dissolved oxygen is modeled, alter sourcesink(kwq,lwq,LDO) to reflect resp
IF (IDO ==1) THEN
   sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) - roc*resz
END IF

! If PON is modeled, alter sourcesink(kwq,lwq,LPON) to include mortality
IF (IPON == 1) THEN
    sourcesink(kwq,lwq,LPON) = sourcesink(kwq,lwq,LPON) + rnc * morz
END IF

! If POP is modeled, alter sourcesink(kwq,lwq,LPOP) to include mortality
IF (IPOP == 1) THEN
    sourcesink(kwq,lwq,LPOP) = sourcesink(kwq,lwq,LPOP) + rpc * morz
END IF

! If POC is modeled, alter sourcesink(kwq,lwq,LPOC) to include mortality
IF (IPOC == 1) THEN
    sourcesink(kwq,lwq,LPOC) = sourcesink(kwq,lwq,LPOC) + morz
END IF

! If NH4 is modeled, alter sourcesink(kwq,lwq,LNH4) to include excretion
IF (INH4 == 1) THEN
    sourcesink(kwq,lwq,LNH4) = sourcesink(kwq,lwq,LNH4) + rnc * exz
END IF

! If PO4 is modeled, alter sourcesink(kwq,lwq,LPO4) to include excretion
IF (INH4 == 1) THEN
    sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) + rpc * exz
END IF

END SUBROUTINE sourceZOO

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
