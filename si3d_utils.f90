!************************************************************************
  MODULE si3d_utils
!************************************************************************
!
!  Purpose: Include routines which are useful for other modules (such
!           as computing dates, or input-output routines, allocating
!           space for variables,
!
!-------------------------------------------------------------------------

  USE omp_lib
  USE si3d_Types
  USE si3d_ecomod

  IMPLICIT NONE
  SAVE

  CONTAINS




!************************************************************************
SUBROUTINE input
!************************************************************************
!
!  Purpose: To read si3d model input parameters.
!
!------------------------------------------------------------------------

   !.....Local variables.................................................
   CHARACTER(LEN=12) :: input_file = "si3d_inp.txt"
   INTEGER :: ios, nn, i, j, istat

   !.....Open input parameter file.....
   OPEN (UNIT=i5, FILE=input_file, STATUS="OLD", IOSTAT=ios)
   IF (ios /= 0) CALL open_error ( "Error opening "//input_file, ios )


   !.....Read header record containing comments about run................
   READ (UNIT=i5, FMT='(/(A))', IOSTAT=ios) title
   IF (ios /= 0) CALL input_error ( ios, 1)

   !.....Read & define start date of the run.............................
   READ (UNIT=i5,FMT='(///(14X,G20.2))',IOSTAT=ios) iyr0,imon0,iday0,ihr0
   IF (ios /= 0) CALL input_error ( ios, 2 )
   CALL compute_date (0.0)

   !.....Read space-time domains, cell size & time step .................
   READ (UNIT=i5,FMT='(///(14X,G20.2))',IOSTAT=ios) xl,yl,zl,tl,idx,idy, &
   & idz,dzmin, datadj, zetainit,idt, ibathyf
   IF (ios /= 0) CALL input_error ( ios, 3 )

   ! ... Read parameters controlling solution algorithm .................
   READ (UNIT=i5,FMT='(///(14X,G20.2))',IOSTAT=ios) itrap,niter,ismooth, &
       & beta, iturb, Av0, Dv0, iadv, itrmom, ihd, Ax0, Ay0, f, theta,   &
       & ibc,isal, itrsca, cd, ifsurfbc, dtsurfbc, cw, wa, phi, idbg, num_threads
   IF (ios /= 0) CALL input_error ( ios, 4 )

   !.....Read node numbers for time series output .......................
   READ (UNIT=i5, FMT='(///(14X,I20))', IOSTAT=ios) ipt
   IF (ios /= 0) CALL input_error ( ios, 5 )
   IF (ipt > 0) THEN
     READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) nnodes
     IF (ios /= 0) CALL input_error ( ios, 5 )
     READ (UNIT=i5, FMT='(14X,10I5)', IOSTAT=ios) (inode(nn), nn = 1, nnodes) ! Changed to I5 12/2010 SWA
     IF (ios /= 0) CALL input_error ( ios, 5 )
     READ (UNIT=i5, FMT='(14X,10I5)', IOSTAT=ios) (jnode(nn), nn = 1, nnodes) ! Changed to I5 12/2010 SWA
     IF (ios /= 0) CALL input_error ( ios, 5 )
   ENDIF

   ! ... Read nodes for horizontal plane output ...........................
   READ (UNIT=i5, FMT='(///(14X,I20))', IOSTAT=ios) iop
   IF (ios /= 0) CALL input_error ( ios, 6 )
   IF (iop /= 0) THEN
     READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) n_planes
     IF (ios /= 0) CALL input_error ( ios, 6 )
     IF (n_planes > max_planes) THEN
       PRINT *,'ERROR: # of planes requested > maximum allowed'
       STOP
     END IF
     DO j = 1, n_planes
       ! ... Read number of cells in X-section j
       READ (UNIT=i5, FMT='(14X,I20)' , IOSTAT=ios) p_out(j)
     ENDDO
   ENDIF

   ! ... Read nodes for vertical plane output ...........................
   READ (UNIT=i5, FMT='(///(14X,I20))', IOSTAT=ios) iox
   IF (ios /= 0) CALL input_error ( ios, 7 )
   IF (iox /= 0) THEN
     READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) n_sections
     IF (ios /= 0) CALL input_error ( ios, 7 )
     IF (n_sections > max_sections) THEN
       PRINT *,'ERROR: # of sections requested > maximum allowed'
       STOP
     END IF
     DO j = 1, n_sections
       ! ... Read number of cells in X-section j
       READ (UNIT=i5, FMT='(/14X,I20)' , IOSTAT=ios) n_section_cells(j)
       IF (ios /= 0) CALL input_error ( ios, 7 )
       ! ... Read i coordinates for cells in X-section j
       READ (UNIT=i5, FMT='(14X,10I5)', IOSTAT=ios)                       &  ! Changed to I5 12/2010 SWA
       &    (xinode(j, nn), nn = 1, n_section_cells(j) )
       IF (ios /= 0) CALL input_error ( ios, 7 )
       ! ... Read j coordinates for cells in X-section j
       READ (UNIT=i5, FMT='(14X,10I5)', IOSTAT=ios)                       &  ! Changed to I5 12/2010 SWA
       &    (xjnode(j, nn),nn=1,n_section_cells(j))
       IF (ios /= 0) CALL input_error ( ios, 7 )
     END DO
   ENDIF

   !! ... Read toggles for 3D output  .....................................
   READ (UNIT=i5, FMT='(///(14X,I20))', IOSTAT=ios) ipxml
   IF (ios /= 0) CALL input_error ( ios, 7 )
   READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) itspf
   IF (ios /= 0) CALL input_error ( ios, 7 )

   !.....Read info on open boundaries.....................................
   READ (UNIT=i5, FMT='(///14X,I20)', IOSTAT=ios) nopen
   IF (ios /= 0) CALL input_error ( ios, 8 )
   IF (nopen > 0) THEN
      READ (UNIT=i5, FMT='(14X,G20.2)', IOSTAT=ios) dtsecopenbc
      IF (ios /= 0) CALL input_error ( ios, 8 )
      DO nn = 1, nopen
         READ (UNIT=i5, FMT='(/(14X,I20))', IOSTAT=ios) iside(nn)
         IF (ios /= 0) CALL input_error ( ios, 8 )
         READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) itype(nn)
         IF (ios /= 0) CALL input_error ( ios, 8 )
         READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) isbc(nn)
         IF (ios /= 0) CALL input_error ( ios, 8 )
         READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) jsbc(nn)
         IF (ios /= 0) CALL input_error ( ios, 8 )
         READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) iebc(nn)
         IF (ios /= 0) CALL input_error ( ios, 8 )
         READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) jebc(nn)
         IF (ios /= 0) CALL input_error ( ios, 8 )
      END DO
   END IF

   !.....Read info to output nested grid boundaries ......................
   IF (ioNBTOGGLE > 0 ) THEN
     READ (UNIT=i5,  FMT='(///14X,I20)', IOSTAT=ios) nxNBO
     IF (ios /= 0) CALL input_error ( ios, 8 )
     IF (nxNBO > 0) THEN
       READ (UNIT=i5, FMT='(14X,G20.2)', IOSTAT=ios) ioNBO
       IF (ios /= 0) CALL input_error ( ios, 8 )
       READ (UNIT=i5, FMT='(14X,G20.2)', IOSTAT=ios) xxNBO
       IF (ios /= 0) CALL input_error ( ios, 8 )
       DO nn = 1, nxNBO
         READ (UNIT=i5, FMT='(/(14X,I20))', IOSTAT=ios) isdNBO(nn)
         IF (ios /= 0) CALL input_error ( ios, 8 )
         READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) isbcNBO(nn)
         IF (ios /= 0) CALL input_error ( ios, 8 )
         READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) jsbcNBO(nn)
         IF (ios /= 0) CALL input_error ( ios, 8 )
         READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) iebcNBO(nn)
         IF (ios /= 0) CALL input_error ( ios, 8 )
         READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) jebcNBO(nn)
         IF (ios /= 0) CALL input_error ( ios, 8 )


       END DO
     END IF
   END IF

   !.... Read info for tracers ...........................................
   READ (UNIT=i5, FMT='(///14X,I20)', IOSTAT=ios) ntr
   IF (ios  /= 0) CALL input_error ( ios, 9 )
   iotr = 0; ! Default value
   IF (ntr > 0) THEN
     READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) ecomod
     READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) iotr
     IF (ios  /= 0) CALL input_error ( ios, 9 )
   ENDIF
    print *,"aqui interp"
   !.... Read info for plume models & oxygenation systems ................
   READ (UNIT=i5, FMT='(///14X,I20)', IOSTAT=ios) iopss
   IF (ios  /= 0) CALL input_error ( ios, 10 )
   IF (iopss > 0) THEN
     ! ... Read no. of devices, each with its own inflow/outflow rate
     READ (UNIT=i5, FMT='( 14X,I20)', IOSTAT=ios) npssdev
     IF (ios  /= 0) CALL input_error ( ios, 10 )
     IF ( npssdev > iopss .OR. npssdev < 1) THEN
       PRINT *, '*** ERROR ***'
       PRINT *, 'No. of devices creating point sources & sinks cannot be > '
       PRINT *, 'No. of water columns with point sources & sinks'
       STOP
     ENDIF
     ! Read time in seconds between consecutive records from time files
     READ (UNIT=i5, FMT='(14X,G20.2)', IOSTAT=ios) dtsecpss

     ! ... Allocate space for arrays holding location
     !     of point sources and sinks and devices
     ALLOCATE ( ipss  (iopss), jpss   (iopss), &
                iodev (iopss), STAT=istat)
     IF (istat /= 0) CALL allocate_error ( istat, 100 )

     ! ... Read in locations & characteristics of diffusers
     !     At this point, they are pressumed constants in time
     READ (UNIT=i5, FMT='(A)', IOSTAT=ios) commentline
     i = 0;
     DO j = 1, iopss
       IF (ios/= 0) CALL input_error ( ios, 10)
       READ (UNIT=i5, FMT='(9X,3I4)', IOSTAT=ios)  &
            ipss(j), jpss(j), iodev(j)
       IF (ios/= 0) CALL input_error ( ios, 10 )
       IF (iodev(j).NE.i) i = iodev(j)
     END DO
     IF (i .NE. npssdev) THEN
       PRINT *, '*** ERROR ***'
       PRINT *, 'No. of devices causing point sources'
       STOP
     ENDIF

     ! ... Read sources & sinks specifications -
     CALL PointSourceSinkInput

   ENDIF
   print *,"antes interp"

   ! ... Input instructions & parameters controlling the solution of
   !     the tracer transport equations.
   IF (ntr > 0) THEN

     ! ... Allocate space for some arrays - they are initialized
	 !     only if ecomod < 0, but used allways in determining
	 !     if the tracer transport equation is used or not in subroutine fd.
     !ALLOCATE ( trct0(ntr), trcpk(ntr), trctn(ntr), &
     !          trcx0(ntr), trcy0(ntr), trcz0(ntr), &
     !           trcsx(ntr), trcsy(ntr), trcsz(ntr), STAT=istat)
     !IF (istat /= 0) CALL allocate_error ( istat, 121 )

     ! .... Initialize trct0 and trctn to default values
	 trct0 = 1E7;
	 trctn =   0;

     ! .... Define other input variables
     SELECT CASE (ecomod)
     CASE (-1) ! Tracer Cloud Releases
        PRINT *, 'Tracer Cloud Modelling activated'
        CALL trcinput
     CASE (1) ! Water quality routines
        PRINT *, 'Water Quality Model activated'
        CALL wqinput
     CASE (2) ! Size Structure distribution
        PRINT *, 'Size Structure Model activated'
        CALL szinput
     CASE (3) ! Sediment transport routines
        PRINT *, 'Sediment transport activated'
        CALL sdinput
     END SELECT

   ENDIF

   ! ... Read info for interpolation method (Added 12/2010 by SWA)
   IF (ifsurfbc >=10) THEN
     READ (UNIT=i5, FMT='(///14X,I20)', IOSTAT=ios) iinterp
     IF (ios  /= 0) CALL input_error ( ios, 10 )
     IF (iinterp==2) THEN
       READ (UNIT=i5, FMT='(14X,G20.2)', IOSTAT=ios) gammaB
       IF (ios  /= 0) CALL input_error ( ios, 10 )
       READ (UNIT=i5, FMT='(14X,G20.2)', IOSTAT=ios) delNfactor
       IF (ios  /= 0) CALL input_error ( ios, 10 )
     ENDIF
   ENDIF

   !.....Close input file.....
   CLOSE (UNIT=i5)

   !.....Define frequently used numerical constants and coefficients.....
   dt=idt; dx=idx; dy=idy; ddz=idz; twodt=2.*dt; dtdx=dt/dx; dtdy=dt/dy
   gdtdx=g*dtdx; gdtdy=g*dtdy; gdt2dx2=gdtdx*dtdx; gdt2dy2=gdtdy*dtdy
   !gthx=gdtdx*theta; gthy=gdtdy*theta; gth1x=gdtdx*2.*(1.-theta)
   !gth1y=gdtdy*2.*(1.-theta);
   cwind=2.*dt*cw*rhoair*wa*wa; alp4=(1.-alp)/4.
   twodx=2.*dx; twody=2.*dy; fourdx=2.*twodx; fourdy=2.*twody
   dxdx=dx*dx; dydy=dy*dy; twodxdx=2.*dxdx; twodydy=2.*dydy; dxdy=dx*dy ! Changed 12/2010 SWA
   beta2=beta/2.; chi1=1.-chi; twochi1=2.*chi1
   im=nint(xl/dx)+1; im1=im+1; i1=2; ndx=im1-i1
   jm=nint(yl/dy)+1; jm1=jm+1; j1=2; ndy=jm1-j1
   nts=tl/dt+.5; apxml = ABS(ipxml)
   isec0=REAL(ihr0)/100.*3600.; !dt_min=idt/60.

   ! ... Generate grid dimensions in Z-direction
   CALL ZGridDimensions

END SUBROUTINE input

!************************************************************************
SUBROUTINE ZGridDimensions
!************************************************************************
!
!  Purpose: To read in depths of levels between consecutive layers
!
!------------------------------------------------------------------------

   !.....Local variables.................................................
   CHARACTER(LEN=14) :: input_file = "si3d_layer.txt"
   INTEGER :: ios, k, istat

   ! ... Variable layer thickness
   IF (ibathyf < 0 ) THEN

     !.....Open input file.....
     OPEN (UNIT=i5, FILE=input_file, STATUS="OLD", IOSTAT=ios)
     IF (ios /= 0) CALL open_error ( "Error opening "//input_file, ios )

     !.....Read information ...
     ! Skip over first line (header)
     READ (UNIT=i5, FMT='(//)', IOSTAT=ios)
     IF (ios /= 0) CALL input_error ( ios, 40 )

     ! Read number of layers from second header record
     READ (UNIT=i5, FMT='(10X,I11)', IOSTAT=ios) km1
     IF (ios /= 0) CALL input_error ( ios, 41 )

     !..... Allocate space for zlevel array
     ALLOCATE (zlevel(1:km1), STAT=istat )
     IF (istat /= 0) CALL allocate_error ( istat, 0 )

     ! .... Read array with levels to layer interfaces
     DO k = 1, km1
       READ (UNIT=i5, FMT='(10X,G11.2)', IOSTAT=ios) zlevel(k)
       IF (ios /= 0) CALL input_error ( ios, 42 )
     END DO

     ! ... Generate grid dimensions in z-direction
     km  = km1 - 1  ;
     k1  = 2        ;
     ndz = km1 - k1 ;

   ! ... Constant layer thickness (original si3d)
   ELSE

     ! ... Generate grid dimensions in z-direction
     !km=CEILING((zl-dzmin)/ddz)+1
     km=CEILING((zl)/ddz)+1
     km1=km+1; k1=2; ndz=km1-k1

     !..... Allocate space for zlevel array
     ALLOCATE (zlevel(1:km1), STAT=istat )
     IF (istat /= 0) CALL allocate_error ( istat, 0 )

     !.....Initialize arrays with levels to layer interfaces
     zlevel(k1) = 0.0;
     DO k = k1+1, km1
       zlevel(k)=zlevel(k-1)+ddz
     END DO
     zlevel(k1) = -100.
     zlevel(1 ) = -100.

   ENDIF


END SUBROUTINE ZGridDimensions

!************************************************************************
SUBROUTINE AllocateSpace
!************************************************************************
!
!  Purpose: To allocate space for model arrays at a size determined during
!           execution. These arrays are 'permanently'
!           allocated for the entire duration of a model run.
!
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!------------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: istat, one=1



   !..... Allocate geometry arrays
   ALLOCATE (kmz   (lm1), k1z(lm1), &
             k1u   (lm1), k1v(lm1), STAT=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 2 )

   !.....Allocate matrix solution arrays.....
   ALLOCATE ( ubar(1), rparm(16), iparm(25), ip(1), jp(1), STAT=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 3 )

   !.....Allocate other miscellaneous arrays.....Map 2D-l into 3D-(i,j) indexes
   ALLOCATE ( ds(km), STAT=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 4 )

   !.....Allocate arrays with variables at zeta-pts in 3D space &
   !     arrays used in solution procedures....
   ALLOCATE ( s   (lm1), sp  (lm1), &
            & spp (lm1),                &
            & sx  (lm1), sy  (lm1), &
            & dd  (lm1), qq  (lm1), &
            & eagx(lm1), earx(lm1), &
            & eagy(lm1), eary(lm1), &
            & rr  (lm1), hhs (lm1), &
            & hhu (lm1), hhv (lm1), &
            & uair(lm1), vair(lm1), &
            & cdw (lm1), STAT=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 5 )

   ! .... Allocate arrays used in model output
   ALLOCATE(  uout (km1) , vout (km1) , wout(km1),  &
            & Avout(km1) , Dvout(km1) , sal1(ndz),  &
            & uhout(km1) , scout(km1),  trout(km1,ntrmax), STAT=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 7 )

   ! ...  Allocate space for output routines
   ALLOCATE ( interior_plane_points   ( n_planes  ),  &
            & interior_section_points ( n_sections), STAT = istat  )
   IF (istat /= 0) CALL allocate_error ( istat, 8 )

   ! .... Allocate arrays used to store velocity boundary conditions
   IF (nopen > 0) THEN
     ALLOCATE(  uhWB  (km1,jm1), uhEB  (km1,jm1),   &
                huWB  (km1,jm1), huEB  (km1,jm1),   &
                vhSB  (km1,im1), vhNB  (km1,im1),   &
                hvSB  (km1,im1), hvNB  (km1,im1),   &
                uhWBpp(km1,jm1), uhEBpp(km1,jm1),   &
                huWBpp(km1,jm1), huEBpp(km1,jm1),   &
                vhSBpp(km1,im1), vhNBpp(km1,im1),   &
                hvSBpp(km1,im1), hvNBpp(km1,im1),   &
                uhWBp (km1,jm1), uhEBp (km1,jm1),   &
                huWBp (km1,jm1), huEBp (km1,jm1),   &
                vhSBp (km1,im1), vhNBp (km1,im1),   &
                hvSBp (km1,im1), hvNBp (km1,im1), STAT=istat )
     IF (istat /= 0) CALL allocate_error ( istat, 9 )
   ENDIF

END SUBROUTINE AllocateSpace

!************************************************************************
SUBROUTINE AllocateSpace2
!************************************************************************
!
!  Purpose: To allocate space for model arrays at a size determined during
!           execution. These arrays are 'permanently'
!           allocated for the entire duration of a model run.
!
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!------------------------------------------------------------------------

   ! ... Local variables
   INTEGER:: istat

   ALLOCATE ( ij2l(im1,jm1), STAT=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 10 )

   ALLOCATE ( l2i(lm1), l2j(lm1), &
   &          lEC(lm1), lWC(lm1), &
   &          lNC(lm1), lSC(lm1), STAT=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 11 )

   ALLOCATE ( lh(num_threads), lh_aux(num_threads), &
   &          lhi(num_threads), lhf(num_threads), &
   &          id_column(lm1+((jm*num_threads*2)-(jm*2)) ),   &
   &		  lhiE(num_threads), lhfE(num_threads), &
   &          lhiW(num_threads), lhfW(num_threads), &
   &          iauxs(num_threads), iauxe(num_threads), &
   &          ph(num_threads))

   ALLOCATE ( lhiCE(num_threads), lhfCE(num_threads), &
   &          lhiCN(num_threads), lhfCN(num_threads), &
   &          id_columnCE(lm1+((jm*num_threads*2)-(jm*2)) ),   &
   &          id_columnCN(lm1+((jm*num_threads*2)-(jm*2)) ),   &
   &		  lhiECE(num_threads), lhfECE(num_threads), &
   &          lhiWCE(num_threads), lhfWCE(num_threads), &
   &          lhiWCN(num_threads), lhfWCN(num_threads), &
   &          lhiECN(num_threads), lhfECN(num_threads))

   ALLOCATE ( coeffA(lm,5), jcoefA(lm,5), &
   &          rhs(lm), zeta(lm))

   ALLOCATE ( nnH(num_threads,maxnopen), nnHH(num_threads,maxnopen) )

   ALLOCATE (isbcH(maxnopen,num_threads),iebcH(maxnopen,num_threads),jsbcH(maxnopen,num_threads), nopenH(num_threads), &
    & jebcH(maxnopen,num_threads),noh2no(maxnopen,num_threads),eiptNBI(maxnopen,num_threads),siptNBI(maxnopen,num_threads), thrsNGBp(maxnopen), &  !MAC
    & eiptNBIH(maxnopen,num_threads),siptNBIH(maxnopen,num_threads),isbcHH(maxnopen,num_threads),iebcHH(maxnopen,num_threads), thrsNGB(maxnopen), & !MAC
    & areatot(maxnopen),flag(maxnopen),nopth(maxnopen),cont(maxnopen),noh2noH(maxnopen,num_threads),nopenHH(num_threads),contNG(maxnopen))

   !.....Allocate arrays at u-pts.....
   ALLOCATE ( uh(km1,lm1), uhp(km1,lm1), uhpp(km1,lm1), kh(km1,lm1),      &
            & u (km1,lm1), up (km1,lm1), upp (km1,lm1),      &
            & ex(km1,lm1), agx(km1,lm1), arx (km1,lm1),ex2(km1,lm1),      &
            & hu(km1,lm1), hup(km1,lm1), hupp(km1,lm1),extr(km1,lm1),STAT=istat)
   IF (istat /= 0) CALL allocate_error ( istat, 12 )

   !.....Allocate arrays at v-pts.....
   ALLOCATE ( vh(km1,lm1), vhp(km1,lm1), vhpp(km1,lm1),      &
            & v (km1,lm1), vp (km1,lm1), vpp (km1,lm1),      &
            &              agy(km1,lm1), ary (km1,lm1),      &
            & hv(km1,lm1), hvp(km1,lm1), hvpp(km1,lm1),STAT=istat)
   IF (istat /= 0) CALL allocate_error ( istat, 13 )

   !.....Allocate arrays at pressure-pts.....
   ALLOCATE ( h   (km1,lm1), hp  (km1,lm1), hpp   (km1,lm1), &
            & sal (km1,lm1), salp(km1,lm1), salpp (km1,lm1), &
            & rhop(km1,lm1), STAT=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 14 )

   !.....Allocate arrays at vertical interfaces between pressure-pts.....
   ALLOCATE ( Av  (km1,lm1), Dv  (km1,lm1), Dvm(km1,lm1),             &!Andrea PT
              wp  (km1,lm1), STAT=istat    )
   IF (istat /= 0) CALL allocate_error ( istat, 15 )

   ! ....Allocate arrays at pressure points
   ALLOCATE ( QswFr     (km1,lm1), &
           &  HeatSource(km1,lm1), STAT=istat)
   IF (istat /= 0) CALL allocate_error ( istat, 16 )

   !.....Allocate horizontal diffusion, th, and th1 arrays.....
   ALLOCATE (haypp(km1,lm1),haxpp(km1,lm1), &
             hdypp(km1,lm1),hdxpp(km1,lm1), &
             th   (km1,lm1),th1  (km1,lm1),haxpp2(km1,lm1) , &
             haypp2(km1,lm1),th2(km1,lm1),th12(km1,lm1),  &
             hayppsal(km1,lm1),haxpptr(km1,lm1),haxppsal(km1,lm1),STAT=istat)
   IF (istat /=0) CALL allocate_error ( istat, 17 )

   !.... Allocate space for Higher-Order turbulence models......
   IF (iturb > 0) THEN
     ALLOCATE ( q2 (km1,lm1), q2p (km1,lm1 ), q2pp (km1,lm1),   &
                q2l(km1,lm1), q2lp(km1,lm1 ), q2lpp(km1,lm1),   &
                dsT(km1)    , aaT (3, km1+1), sal1T(ndz+1), STAT = istat )
     IF (istat /= 0) CALL allocate_error ( istat, 18 )
   ELSE IF (iturb < 0) THEN
     ALLOCATE ( si3dtke (km1,lm1), &
                si3deps (km1,lm1), &
                si3dlen (km1,lm1), STAT = istat )
     IF (istat /= 0) CALL allocate_error ( istat, 18 )
   ENDIF

   ! .... Allocate arrays used in advection for tracers & plumes &
   !      when scalar transport is done with flux-limiters ..............
   ALLOCATE(  fluxX (km1,lm1),fluxY(km1,lm1),fluxZ(km1,lm1),  &
              fluxXtr (km1,lm1),fluxY2(km1,lm1),fluxZ2(km1,lm1),  &
              fluxXsal(km1,lm1), STAT = istat )
   IF (istat /= 0) CALL allocate_error ( istat, 19 )

   ! .... Allocate arrays for tracers ...................................
   IF (ntr > 0 ) THEN
     ALLOCATE(  tracer     (km1,lm1,ntr),  &
                tracerpp   (km1,lm1,ntr),  &
                STAT = istat )
     IF (istat /= 0) CALL allocate_error ( istat, 20 )
   ENDIF

   ! .... Allocate arrays used in oxygenation simulations ...............

   ALLOCATE(Qpss(km1,iopss), Tpss(km1,iopss), iopssH(num_threads), &
            Rpss(km1,iopss,ntr), ioph2iop(iopss,num_threads),STAT=istat)
   IF (istat /= 0) CALL allocate_error ( istat, 21 )


END SUBROUTINE AllocateSpace2

!************************************************************************
SUBROUTINE bathy
!************************************************************************
!
!  Purpose: To read the file with bathymetry for the basin. The depths
!           are assumed to be defined at the corners of each computational
!           cell and stored in the array 'h4'. The average depths along
!           each cell face are computed and used to define the 3-D layer
!           thickness array hhs of bottom depths from datum at zeta points.
!           2-D logical mask arrays are defined with .TRUE. values for all
!           wet cells and .FALSE. values for all dry cells. The 2-D arrays,
!           hhu & hhv, of depths at u- & v-points are defined from hhs.
!           The depths are in meters below a datum.
!  13/11/08 Input depths at zeta point
!
!-------------------------------------------------------------------------

   !.....Local variables.....
   CHARACTER(LEN=50) :: bathymetry_file = "h"
   INTEGER :: i, j, k, l, c, ios, istat, kb, is, ie, js, je
   INTEGER :: imm, jmm, ncols, ncols1, nc, nn, ia, ib
   REAL :: hs1, udepth, vdepth, ztop
   CHARACTER(LEN=14) :: fmt
   REAL, DIMENSION(:,:), POINTER :: h4

   !.....Open bathymetry file.....
   OPEN (UNIT=i5, FILE=bathymetry_file, STATUS="OLD", IOSTAT=ios)
   IF(ios /= 0) CALL open_error ( "Error opening "//bathymetry_file, ios )

   !.....Read header information.....
   READ (UNIT=i5, FMT='(37X,I5,6X,I5,8X,I5//)', IOSTAT=ios) imm, jmm, ncols  ! Changed to I5 12/2010 SWA
   IF (ios /= 0) CALL input_error ( ios, 11 )

   !.....Check grid dimensions against input parameters.....
   IF ((im /= imm+1) .OR. (jm /= jmm+1)) THEN
      PRINT *, " ****ERROR -- Grid size computed from input file does not"
      PRINT *, "              agree with the header in the bathymetry file"
      PRINT '(4(A,I5))', " im=", im, " imm+1=", imm+1, " jm=", jm, " jmm+1=", jmm+1
      PRINT '(A/)', " "
      PRINT *, " ****STOPPING si3d in SUBROUTINE bathy"
      STOP
   END IF

   !.....Write data format to an internal file.....
   IF (ibathyf == 1) THEN ! Stockton
     WRITE (UNIT=fmt, FMT='("(5X,", I4, "G5.0)")') ncols
   ELSE                   !General case (Deeper lakes > 100 m)
     WRITE (UNIT=fmt, FMT='("(5X,", I4, "G5.0)")') ncols
   ENDIF

   !.....Allocate space for bathymetry array.....
   ALLOCATE ( h4(1:im1, 1:jm1), STAT=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 15 )

   !.....Allocate logical mask arrays.....
   ALLOCATE ( mask2d(im1,jm1    ), STAT=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 1 )

   !.....Allocate logical mask arrays 1D.....
   ALLOCATE ( mask(im1*jm1    ), STAT=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 1 )


   !.....Read bathymetry.....
   ncols1 = ncols     - 1;
   nc     = imm/ncols + 1;
   ia     = -ncols1   + 1;
   DO nn = 1, nc
      ia = ia + ncols
      ib = ia + ncols1
      IF ( ia > im ) EXIT
      IF ( ib > im ) ib = im
      DO j = jm, j1, -1 ! change 1 to j1
         !print *,"jjjj:",j
         READ (UNIT=i5, FMT=fmt, IOSTAT=ios) (h4(i,j), i = ia, ib)
         IF (ios /= 0) CALL input_error ( ios, 12 )
      END DO
   END DO

   !.....Close bathymetry file.....
   CLOSE (UNIT=i5)

   !.....Change units (from dm to m) and adjust bathymetry datum.....
   h4 = h4*0.1 + datadj

   !.....Be sure mask is false in the fictitious
   !     row/column around the outer edge of the grid.....
   h4(  1,1:jm1) = 0.0; h4(1:im1,jm1) = 0.0
   h4(im1,1:jm1) = 0.0; h4(1:im1,  1) = 0.0

   !.....Define mask2d array.....
   DO j = 1, jm1; DO i = 1, im1
     IF ( h4(i,j) > 0.0 ) THEN
        mask2d(i,j) = .TRUE.
     ELSE
        mask2d(i,j) = .FALSE.
     END IF
   END DO; END DO

   !.....Define mask array.....
   cm1=0.0
   DO i = 1, im1; DO j = 1, jm1
   cm1=cm1+1
     IF ( h4(i,j) > 0.0 ) THEN
        mask(cm1) = .TRUE.
     ELSE
        mask(cm1) = .FALSE.
     END IF
   END DO; END DO


   !.....Add fictitious row/column of depths around grid.....
   h4(1,j1:jm  ) = h4(2,j1:jm )         ! west side
   h4(i1:im,jm1) = h4(i1:im,jm)         ! north side
   h4(im1,j1:jm) = h4(im,j1:jm)         ! east side
   h4(i1:im,1  ) = h4(i1:im,2 )         ! south side
   ! Take care of corners
   h4(  1,  1) = h4( 2, 2); h4(  1,jm1) = h4( 2,jm)
   h4(im1,jm1) = h4(im,jm); h4(im1,  1) = h4(im, 2)

   !.....Compute the first and last column and row of grid
   !     with wet points. Variables used in subr. solver
   ifirst = im; ilast = i1; jfirst = jm; jlast = j1
   DO j = j1, jm; DO i = i1, im
     IF ( .NOT. mask2d(i,j) ) CYCLE   ! Ignore dry points
     IF ( i < ifirst ) ifirst = i
     IF ( i > ilast  ) ilast  = i
     IF ( j < jfirst ) jfirst = j
     IF ( j > jlast  ) jlast  = j
   END DO; END DO

   !.....Compute the first and last column and row of grid
   !     with wet points. Variables used in subr. solver
   lfirst = cm1; llast = 1;
   DO c = 1, cm1
     IF ( .NOT. mask(c) ) CYCLE   ! Ignore dry points
     IF ( c < lfirst ) lfirst = c
     IF ( c > llast  ) llast  = c
   END DO

   ! ... Find dimensions for 2D-lk arrays
   l = 0
   DO c = 1, cm1
     IF ( .NOT. mask(c) ) CYCLE
     l = l + 1
   END DO
   lm = l ; lm1 = lm + 1;

   !.....Allocate
   ALLOCATE ( c2l(lm1    ), STAT=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 1 )

   !.....Allocate
   ALLOCATE ( l2c(im1*jm1    ), STAT=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 1 )

!   l = 0
!   c2l = cm1
!   l2c = lm1
!   DO c = 1, cm1
!     IF ( .NOT. mask(c) ) CYCLE
!     l = l + 1
!     c2l(l)=c
!     l2c(c)=l
!   END DO


   ! ... Allocate space for arrays in 2D-lk coordinates
   CALL AllocateSpace2

   ! ... Mapping functions from/to 3D-(i,j)- to/from 2D-l-indexes
   l = 0
   c2l = cm1
   l2c = lm1
   c = 0
   ij2l = lm1; ! For dry(i,j) the map function ij2l will yield lm1
   DO i = 1, im1; DO j = 1, jm1;
   c = c + 1
     IF(mask2d(i,j)) THEN
     l = l + 1
     l2i (l  ) = i ; ! Goes from ipl to index i in the i,j plane
     l2j (l  ) = j ; ! Goes from ipl to index j in the i,j plane
     ij2l(i,j) = l ; ! Goes from i,j to the ipl index in the ipl line
!     print *,"i",i,"j",j,"l",l
     c2l(l)=c
     l2c(c)=l
     END IF
   END DO; END DO;



   ! ... Assign E,W,N,S colums for each l-column in the 2D-l space
   l = 0
   DO i = i1, im; DO j = j1, jm;
     IF(.NOT. mask2d(i,j)) CYCLE
     l = l + 1
     lEC(l) = ij2l(i+1,j); ! Defines water column East  of l
     lWC(l) = ij2l(i-1,j); ! Defines water column West  of l
     lNC(l) = ij2l(i,j+1); ! Defines water column North of l
     lSC(l) = ij2l(i,j-1); ! Defines water column South of l
   END DO; END DO;
!.....Allocate space for arrays.....
   CALL AllocateSpace

   !.....Define layer No. for bottom cell (kmz) &
   !     bottom depth from datum at zeta-points (hhs).....
   hhs = ZERO;
   DO j = 1,jm1; DO i = 1,im1

     SELECT CASE (mask2d(i,j))

     CASE (.FALSE.)      ! Cells with all dry layers

!       kmz(i,j) = km1

     CASE (.TRUE.)      ! Cells with wet layers

!            i = l2i(l);
!           j = l2j(l);
            l=ij2l(i,j)
       ! ... Take depth at zeta-point as given in h file
       hs1 = h4(i,j)

       ! ... Compute No. of wet layers & depths at zeta-points
       ! --- A. constant layer thickness (ibathyf >= 0)
       IF (ibathyf >= 0) THEN
         kmz(l) = FLOOR(hs1/ddz)
         IF( (hs1-kmz(l)*ddz) > dzmin ) THEN
           kmz(l) = kmz(l) + 1;
           hhs(l) = hs1
         ELSE
           hhs(l) = kmz(l)*ddz
         ENDIF
         ! Add the fictitious first layer above the water surface
         kmz(l) = kmz(l) + 1

       ! --- B. Variable layer thickness (ibathyf < 0)
       ELSE
         DO k = k1, km
           IF (zlevel(k+1)>=hs1) THEN
             ! ... Option 1 - it works
             hhs(l) = zlevel(k)+MAX(dzmin,(hs1-zlevel(k)));
             kmz(l) = k
             EXIT
             !! ... Option 2 - preferable from a theoretical stand point
             !IF ( hs1-zlevel(k) > dzmin) THEN
             !  hhs(i,j) = hs1
             !  kmz(i,j) = k
             !ELSE
             !  hhs(i,j) = zlevel(k)
             !  kmz(i,j) = k - 1
             !  IF ( kmz(i,j) == 1 ) THEN
             !     hhs(i,j) = dzmin
             !     kmz(i,j) = k1
             !  ENDIF
             !ENDIF
             !EXIT
           ENDIF
         ENDDO
       ENDIF

     END SELECT

   END DO;END DO
   kmz(lm1) = km1
   !.....Define bottom depths from datum at u- and v- points
   hhu = ZERO;
   hhv = ZERO;

   DO j = 1,jm; DO i = 1,im;
      IF(mask2d(i+1,j) .AND. mask2d(i,j)) THEN
        l=ij2l(i,j)
        hhu(l) = MIN(hhs(lEC(l)), hhs(l))
      ENDIF
      IF(mask2d(i,j+1) .AND. mask2d(i,j)) THEN
        l=ij2l(i,j)
        hhv(l) = MIN(hhs(lNC(l)), hhs(l))
      ENDIF

   END DO;END DO

   !.....Process and output the bathymetry needed for
   !     graphics and particle tracking if ioutg=1.....
   IF ( ipxml > 0 ) CALL outg ( h4 )

   !.....Deallocate h4 pointer array.....
   DEALLOCATE ( h4 )

END SUBROUTINE bathy

!***********************************************************************
FUNCTION parab ( frstpt, x, fx, dx )
!***********************************************************************
!
!  Purpose: To interpolate parabolically between the functional values
!           within the array fx.
!
!-----------------------------------------------------------------------

   REAL, DIMENSION(:), INTENT(IN) :: fx      ! Assumed-shape array
   REAL, INTENT(IN) :: frstpt, x, dx
   REAL :: parab
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
   parab=fx(m)+0.5*theta*(fx(m+1)-fx(m-1)+theta*(fx(m+1)+fx(m-1)-2.0*fx(m)))

END FUNCTION parab

!***********************************************************************
FUNCTION linear ( frstpt, x, fx, dx )
!***********************************************************************
!
!  Purpose: To interpolate linearly between the functional values
!           within the array fx.
!


!-----------------------------------------------------------------------

   REAL, DIMENSION(:), INTENT(IN) :: fx      ! Assumed-shape array
   REAL, INTENT(IN) :: frstpt, x, dx
   REAL :: linear
   REAL :: om, theta
   INTEGER :: m

   m = FLOOR((x - frstpt)/dx)
   theta = (x - frstpt - m*dx) / dx
   linear=(1-theta)*fx(m+1)+theta*fx(m+2)


END FUNCTION linear

!***********************************************************************
SUBROUTINE outr
!***********************************************************************
!
!  Purpose: To output model run parameters & performance measures to file.
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   CHARACTER(LEN=12) :: output_file="si3d_out.txt"
   CHARACTER :: date*8, time*10, zone*5
   INTEGER, DIMENSION(8) :: values
   INTEGER :: ios

   !.....Open output file.....
   OPEN (UNIT=i6, FILE=output_file, IOSTAT=ios)
   IF(ios /= 0) CALL open_error ( "Error opening "//output_file, ios )

   !.....Get date and time of run.....
   CALL date_and_time ( date, time, zone, values )

   !.....Output run parameters.....
   WRITE (UNIT=i6, FMT='(A)') title
   WRITE (UNIT=i6, FMT='("Run number = ",A8,A4,",  Start date of run:  ",I2, &
                       & "/",I2,"/",I4," at ",I4.4," hours")')date,time(1:4),&
                       & imon,iday,iyr,ihr

   !.....Output space-time domains, cell size & time step .................
   WRITE (UNIT=i6,FMT=1) xl,yl,zl,idx,idy,idz,tl,idt, dzmin,zetainit

   ! ... Read parameters controlling solution algorithm .................
   WRITE (UNIT=i6,FMT=2) itrap,niter,ismooth,                            &
       & beta, iturb, Av0, Dv0, iadv, itrmom, ihd, Ax0, Ay0, f, theta,   &
       & ibc,isal, itrsca, nopen, cd, ifsurfbc, dtsurfbc, cw, wa, phi, idbg

 1 FORMAT (/              &
     &' xl=   '  , F9.1,/ &
     &' yl=   '  , F9.1,/ &
     &' zl=   '  , F9.3,/ & ! idt real
     &' idx=  '  , F9.3,/ & ! idt real
     &' idy=  '  , F9.3,/ & ! idt real
     &' idz=  '  , F9.3,/ & ! idt real
     &' tl=   '  , G9.2,/ &
     &' idt=  '  , G9.2,/ & ! idt real
     &' dzmin='  , F9.2,/ &
     &' zeta0='  , F9.2 / )

 2 FORMAT(/              &
     &' itrap= ' ,  I9,/ &
     &' niter= ' ,  I9,/ &
     &' smooth=' ,  I9,/ &
     &' beta=  ' ,F9.4,/ &
     &' iturb= ' ,  I9,/ &
     &' Av0=   ' ,F9.4,/ &
     &' Dv0=   ' ,F9.4,/ &
     &' iadv=  ' ,  I9,/ &
     &' trmom= ' ,  I9,/ &
     &' ihd=   ' ,  I9,/ &
     &' Ax0=   ' ,F9.4,/ &
     &' Ay0=   ' ,F9.4,/ &
     &' f=     ' ,F9.4,/ &
     &' theta= ' ,F9.4,/ &
     &' ibc=   ' ,  I9,/ &
     &' isal=  ' ,  I9,/ &
     &' trsal= ' ,  I9,/ &
     &' nopen= ' ,  I9,/ &
     &' cd=    ' ,F9.4,/ &
     &' isbc=  ' ,  I9,/ &
     &' tsbc=  ' ,F9.2,/ &
     &' cw=    ' ,F9.2,/ &
     &' wa=    ' ,F9.2,/ &
     &' phi=   ' ,F9.2,/ &
     &' idbg=  ' ,  I9 / )

END SUBROUTINE outr

!***********************************************************************
SUBROUTINE outt(n,thrs)
!***********************************************************************
!
!  Purpose: To write output to timefile(s). A separate timefile is
!           opened for each node where output is requested.
!
!-----------------------------------------------------------------------

   INTEGER,INTENT(IN) :: n
   REAL, INTENT(IN) :: thrs

   !.....Local variables.....
   CHARACTER :: date*8, time*10, zone*5
   CHARACTER(LEN=9)  :: nodeno    ="         "
   CHARACTER(LEN=15) :: filenm    ="               "
   REAL :: qu, stidal, tdays
   INTEGER, DIMENSION(8)     :: values
   INTEGER :: nn, i, j, k, l, kkk, itdays, ios, nchar, it, laux
   INTEGER, SAVE :: i10, i30, i60
   LOGICAL, SAVE :: first_entry = .TRUE.

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)

   !.....Open timefiles on first entry into the subroutine.....
   IF( first_entry ) THEN
      first_entry = .FALSE.
      DO nn = 1, nnodes

         i = inode(nn)
         j = jnode(nn)

         ! Convert node numbers to a character variable
         CALL nodech ( i, j, nodeno, nchar )

         ! Name the standard si3d timefile
         filenm = 'tf'//nodeno(1:nchar)//'.txt'

         ! Open the timefile
         i60 = i6 + nn    ! Use file numbers 61-80
         OPEN ( UNIT=i60, FILE=filenm, IOSTAT=ios )
         IF(ios /= 0) CALL open_error ( "Error opening "//filenm, ios )

         !.....Get date and time of run.....
         CALL date_and_time ( date, time, zone, values )

         !.....Output run title and column headings for standard si3d format.....
         WRITE (UNIT=i60, FMT='(A)') title
         WRITE (UNIT=i60, FMT='("Run number = ", A8, A4,                      &
                 & ",  Start date of run:  ",I2,                              &
                 & "/",I2,"/",I4," at ",I4.4," hours")') date,time(1:4),      &
                 & imon,iday,iyr,ihr
         IF (idt .GE. 0.01 .AND. ddz .GE. 0.01) THEN ! idt real
           WRITE (UNIT=i60, FMT=1) i,j,(kmz(ij2l(i,j))-k1+1),idt,hhs(ij2l(i,j)),ddz,         &
                                    & iexplt,itrap,cd,ismooth,beta,niter,iextrp, &
                                    & f,tramp,iupwind
         ELSE ! idt real
           WRITE (UNIT=i60, FMT=8) i,j,(kmz(ij2l(i,j))-k1+1),idt,hhs(ij2l(i,j)),ddz,         &
                                    & iexplt,itrap,cd,ismooth,beta,niter,iextrp, &
                                    & f,tramp,iupwind
         ENDIF ! idt real
       1 FORMAT( "i =",I4,"  j =",I4, "  km =", I4,"   dt =",F5.2," sec",     & ! idt real
                  & "  hhs =", F6.3," m","   dz =", F5.2," m"/                & ! idt real
                  & "iexplt =",I2, "  itrap =", I2, "  cd = ", F7.4, 2X,      &
                  & "ismooth =", I2, "  beta =", F6.3," niter =", I2/         &
                  & "iextrp =",  I2, "  f =", F7.4, "  tramp=", F9.1, 2X,     &
                  & "iupwind =", I2 )
       8 FORMAT( "i =",I4,"  j =",I4, "  km =", I4,"   dt =",F5.4," sec",     & ! idt real
                  & "  hhs =", F6.3," m","   dz =", F5.4," m"/                & ! idt real
                  & "iexplt =",I2, "  itrap =", I2, "  cd = ", F7.4, 2X,      &
                  & "ismooth =", I2, "  beta =", F6.3," niter =", I2/         &
                  & "iextrp =",  I2, "  f =", F7.4, "  tramp=", F9.1, 2X,     &
                  & "iupwind =", I2 )
         WRITE (UNIT=i60, FMT=2)
       2 FORMAT( 1X,"   time     ","  step     ","  zeta ","   depth   " &
                    "    u       ","  v      "," w       ",              &
                    "   Av       ","        Dv      ","  scalar    "," Tracers-> " )
         WRITE (UNIT=i60, FMT=3)
       3 FORMAT( 1X,"    hrs     ","   no      ","   cm       "," m    " &
                    "   cm/s   "  ,"   cm/s   " ,"  cm/s      " ,      &
                    " cm2/s      ","     cm2/s     ","  oC        ","  g/l   -> " )

      END DO
   END IF

   !.....Output values at time step  n  to timefile(s).....
   DO nn = 1, nnodes

      i = inode(nn);
      j = jnode(nn);
      i60 = i6 + nn;

      ! ... Map (i,j)- into l-index
      l = ij2l(i,j)

      ! ... Initialize output variables to -99.
      uout  = -99.0E-2 ! 10-2 since the output is cm /s
      vout  = -99.0E-2
      wout  = -99.0E-2
      Avout = -99.0E-4 ! 10-4 since the output is cm2/s
      Dvout = -99.0E-4
      uhout = -99.0
      scout = -99.0
      trout = -99.0

      DO k  = k1, kmz(l)
        IF (h(k,l)<=ZERO) CYCLE
        !uout(k)  = 0.5 * (u  (k,l) + u  (k,lWC(l)))
        !vout(k)  = 0.5 * (v  (k,l) + v  (k,lSC(l)))
        !wout(k)  = 0.5 * (wp (k,l) + wp (k+1,l   ))
        uout(k)  = u (k, l)
        vout(k)  = v (k, l)
        wout(k)  = wp(k, l)
        Avout(k) = 0.5 * (Av (k,l) + Av (k+1,l   ))
        Dvout(k) = 0.5 * (Dv (k,l) + Dv (k+1,l   ))
        uhout(k) = 0.5 * (uh (k,l) + uh (k,lWC(l)))
        scout(k) = sal(k,l)
        IF (ntr>0) THEN
          DO it = 1, ntr
            trout(k,it) = tracer(k,l,it)
          ENDDO
        ENDIF
     END DO

     ! ... Write variables to output file
     IF (ntr <= 0) THEN
       WRITE (UNIT=i60, FMT=4) thrs, n, s(l), (zlevel(k+1),     &
           & uout (k), vout(k), wout(k), Avout(k), Dvout(k),      &
           & scout(k), k =k1,kmz(l))
     ELSE
       WRITE (UNIT=i60, FMT=5) thrs, n, s(l), (zlevel(k+1) ,     &
            & uout (k), vout(k) , wout(k), Avout(k), Dvout(k),     &
            & scout(k  ),                                          &
            & trout(k,1),                                          &
            & trout(k,2),                                          &
            & trout(k,3),                                          &
            & trout(k,4),                                          &
            & trout(k,5),                                          &
            & trout(k,6),                                          &
            & trout(k,7),                                          &
            & trout(k,8),                                          &
            & trout(k,9),                                          &
            & trout(k,10),                                         &
            & trout(k,11),                                         &
            & trout(k,12),                                         &
            & trout(k,13),                                         &
            & trout(k,14),                                         &
            & trout(k,15),                                         &
            & k = k1,kmz(l))
     ENDIF

   4 FORMAT(1X,F10.4,I10,2PF9.2,0PF9.2,2(2PF10.2),2PF9.4,2(4PF15.7),   0PE15.7 / &
	                & ( 30X,0PF9.2,2(2PF10.2),2PF9.4,2(4PF15.7),   0PE15.7 ))
   5 FORMAT(1X,F10.4,I10,2PF9.2,0PF9.2,2(2PF10.2),2PF9.4,2(4PF15.7),16(0PF15.7)/ &
	                & ( 30X,0PF9.2,2(2PF10.2),2PF9.4,2(4PF15.7),16(0PF15.7)))
   END DO

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_outt = t_outt + (etime - btime)

END SUBROUTINE outt



!***********************************************************************
SUBROUTINE outw(n)
!***********************************************************************
!
!  Purpose: To write the wind field provided to the model as boundary
!           condition
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: n

  !.....Local variables.....
  CHARACTER (LEN = 10) :: wind_file
  INTEGER, PARAMETER   :: wind_id0 = 420
  INTEGER :: wind_id
  INTEGER :: i, j, k, ios, istat, n_frames, k_out, m1,m2, kp, kb, c
  INTEGER, SAVE :: ipoints
!  REAL, ALLOCATABLE, DIMENSION(:,:) :: out_array
  INTEGER :: year_out, day_out, mon_out
  REAL    :: hour_out

  ! ... Return if meteorological field is uniform -
  IF ( ifsurfbc < 10 ) RETURN

  !..... Do only on first entry....
  IF( n == 0 ) THEN

    ! Open spacefile on first entry and print specifications for reading the file
    n_frames = nts/MAX(iop,1)
    wind_id = wind_id0
    wind_file = 'wfield.bnr'
    OPEN(unit=wind_id,file=wind_file,FORM='UNFORMATTED',IOSTAT=ios)
    IF(ios /= 0) THEN; PRINT *, "Error opening wind file ", ios; STOP; ENDIF
    !... Write number of time slices to output file
    WRITE(wind_id) n_frames
    ! ... Find out the number of surface points & store it
    ipoints = 0
    DO c=1,cm1
      ! Ignore dry cells
      IF (.NOT. mask(c)  ) CYCLE
      ipoints = ipoints + 1
    END DO
    WRITE(wind_id) ipoints

    ! ... Time stamp
    year_out = iyr
    mon_out  = imon
    day_out  = iday
    hour_out = ihr

    ! ... Output sheets
!    ALLOCATE( out_array ( ipoints, 4 ), STAT=istat )
!    IF (istat /= 0) THEN; PRINT *, 'ERROR allocating space for out_array'; STOP; ENDIF
    k_out = 0;
    DO c=1,cm1
      ! Ignore dry cells
      IF (.NOT. mask(c)) CYCLE
      ! Update counter
      k_out = k_out + 1;
      ! Save wind stress components into output variable
      out_array(k_out,1) = FLOAT(i)
      out_array(k_out,2) = FLOAT(j)
      out_array(k_out,3) = 0.0 ! Wind stress component in the EW direction
      out_array(k_out,4) = 0.0 ! Wind stress component in the NS direction
    END DO

    ! ... Id # for plane file
    wind_id = wind_id0
    ! ... Print time stamp followed by the records
    WRITE(wind_id) n,year_out,mon_out,day_out,hour_out,          &
    &            ((out_array(m1,m2),m2=1,4),m1=1,k_out)
!    DEALLOCATE (out_array)

  ! ... On successive time steps
  ELSE

    ! ... Time stamp
    year_out = iyr
    mon_out  = imon
    day_out  = iday
    hour_out = ihr

    ! ... Output sheets
!    ALLOCATE( out_array ( ipoints, 2 ), STAT=istat )
!    IF (istat /= 0) THEN; PRINT *, 'ERROR allocating space for out_array'; STOP; ENDIF

    k_out = 0;
    DO c=1,cm1
      ! Ignore dry cells
      IF (.NOT. mask(c)) CYCLE
      ! Update counter
      k_out = k_out + 1;
      ! Save wind stress components into output variable
      out_array(k_out,1) = uair(l2c(c)) ! wind component in the EW direction
      out_array(k_out,2) = vair(l2c(c)) ! wind component in the NS direction
    END DO

    ! ... Id # for plane file
    wind_id = wind_id0
    ! ... Print time stamp followed by the records
    WRITE(wind_id) n,year_out,mon_out,day_out,hour_out,          &
    &            ((out_array(m1,m2),m2=1,2),m1=1,k_out)
!    DEALLOCATE (out_array)

END IF

END SUBROUTINE outw



!***********************************************************************
SUBROUTINE outv(n)
!***********************************************************************
!
!  Purpose: To write output @ a cross section to a binary file
!
!-----------------------------------------------------------------------

   INTEGER, INTENT(IN) :: n

   !.....Local variables.....
   CHARACTER (LEN = 12) :: section_file
   INTEGER, PARAMETER   :: section_id0 = 900
   INTEGER :: section_id, year_out, mon_out, day_out, i, j, k, l
   INTEGER :: m1, m2, ios, istat, n_frames, ipoints, k_out
   REAL    :: hour_out
!   REAL, ALLOCATABLE, DIMENSION(:,:) :: out_array

   !.....Open on first entry and print initial conditions ....
   IF( n == 0 ) THEN
   n_frames = nts/MAX(iox,1)
   DO j = 1, n_sections
      section_id = section_id0 + j
      section_file = "section_    "
      IF ( j < 10 ) WRITE ( section_file(9:12), FMT='(I1,"   ")' ) j
      IF(( j >= 10 ) .AND. ( j < 100)) &
      &             WRITE ( section_file(9:12), FMT='(I2,"  " )' ) j
      IF ( j > 100) WRITE ( section_file(9:12), FMT='(I3," "  )' ) j
      OPEN(unit=section_id,file=section_file,FORM='UNFORMATTED',IOSTAT=ios)
      IF(ios /= 0) THEN
        PRINT *, "Error opening section file = ", j, ios; STOP
      ENDIF
      !... Write number of time slices to output file
      WRITE(section_id) n_frames
      ipoints = 0
      CALL CountSectionCells ( ipoints, j )
      interior_section_points (j) = ipoints
      WRITE(section_id) ipoints, km1
    END DO

   ! ... Time stamp
   year_out = iyr
   mon_out  = imon
   day_out  = iday
   hour_out = ihr


   ! ... Output sections
   DO k = 1, n_sections
     ipoints = interior_section_points (k)

!     ALLOCATE( out_array ( ipoints, 10 ), STAT=istat )

!     IF (istat /= 0) THEN
!       PRINT *, 'ERROR in out_V_plane allocating space for out_array'
!     STOP
!     ENDIF
     k_out = 0
     DO m1 = 1, n_section_cells (k)
       i = xinode (k,m1)
       j = xjnode (k,m1)
       l = ij2l(i,j)
       IF ( .NOT. mask(c2l(l))) CYCLE
       DO m2 = k1, kmz(l)
         k_out = k_out + 1
         out_array(k_out,1) = FLOAT(i )
         out_array(k_out,2) = FLOAT(j )
         out_array(k_out,3) = FLOAT(m1)
         out_array(k_out,4) = zlevel(m2+1)-hp(m2,l) ! FLOAT(m1)
         out_array(k_out,5) = uhp (m2,l) ! Changed 12/2010 SWA
         out_array(k_out,6) = vhp (m2,l) ! Changed 12/2010 SWA
         out_array(k_out,7) = wp  (m2,l)
         out_array(k_out,8) = salp(m2,l);
           IF (hp(m2,l)<= 0.) &
           out_array(k_out,8) = -99.;
         IF (ntr > 0) THEN
           out_array(k_out,9) = tracerpp(m2,l,ntr);
             IF (hp(m2,l)<= 0.) &
             out_array(k_out,9) = -99.;
         ELSE
           out_array(k_out,9) = 0.5 * (Av (m2,l) + Av (m2+1,l))
         ENDIF
         out_array(k_out,10) = 0.5 * (Dv (m2,l) + Dv (m2+1,l))
       END DO

     END DO

     ! ... Id # for plane file
     section_id = section_id0 + k

     ! ... Print time stamp followed by the records
     WRITE(section_id) n,year_out,mon_out,day_out,hour_out,        &
    &            ((out_array(m1,m2),m2=1,10),m1=1,ipoints)

!     DEALLOCATE (out_array)

   END DO

!   DEALLOCATE (out_array)
   ELSE

   ! ... Time stamp
   year_out = iyr
   mon_out  = imon
   day_out  = iday
   hour_out = ihr

   ! ... Output sections
   DO k = 1, n_sections
     ipoints = interior_section_points (k)
!     ALLOCATE( out_array ( ipoints, 6 ), STAT=istat )
!     IF (istat /= 0) THEN
!       PRINT *, 'ERROR in out_V_plane allocating space for out_array'
!       STOP
!     ENDIF
     k_out = 0
     DO m1 = 1, n_section_cells(k)
       i = xinode (k,m1)
       j = xjnode (k,m1)
       l = ij2l(i,j)
       IF ( .NOT. mask(c2l(l))) CYCLE
       DO m2 = k1, kmz(l)
         k_out = k_out + 1
         out_array(k_out,1) = uhp (m2,l) ! Changed 12/2010 SWA
         out_array(k_out,2) = vhp (m2,l) ! Changed 12/2010 SWA
         out_array(k_out,3) = wp  (m2,l)
         out_array(k_out,4) = salp(m2,l);
         IF (hp(m2,l)<= 0.) &
           out_array(k_out,4) = -99.;
         IF (ntr > 0) THEN
           out_array(k_out,5) = tracer(m2,l,ntr);
           IF (hp(m2,l)<= 0.) &
             out_array(k_out,5) = -99.;
         ELSE
           out_array(k_out,5) = 0.5 * (Av (m2,l) + Av (m2+1,l))
         ENDIF
         out_array(k_out,6) = 0.5 * (Dv (m2,l) + Dv (m2+1,l))
       END DO
     END DO
     ! ... Id # for plane file
     section_id = section_id0 + k
     ! ... Print time stamp followed by the records
     WRITE(section_id) n,year_out,mon_out, day_out,hour_out,          &
     &            ((out_array(m1,m2),m2=1,6),m1=1,ipoints)
!     DEALLOCATE (out_array)
   END DO

 END IF

END SUBROUTINE outv



!***********************************************************************
SUBROUTINE CountSectionCells ( i1, ksection )
!***********************************************************************
!
!  Purpose: To count wet points in a given X-section (ksection)
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!  02/22/00          F.J. Rueda        Original f90 code
!  03/27/00          F.J. Rueda        Check for dry cells (not counted)
!
!-----------------------------------------------------------------------
!

 INTEGER, INTENT (inout) :: i1
 INTEGER, INTENT (in)    :: ksection

 ! ... Local variables
 INTEGER :: count, i, j, k, m, l

 count = 0
 DO  m = 1, n_section_cells (ksection)
    i = xinode(ksection,m)
    j = xjnode(ksection,m)
    IF( .NOT. mask2D(i,j) ) THEN
      PRINT *, '****************** WARNING **********************'
      PRINT *, 'Section # ', ksection, ': DRY CELL (',i,j,')'
      PRINT *, '*************************************************'
      CYCLE
    END IF
    count = count + kmz(ij2l(i,j)) - 1
 END DO
 i1 = count

END SUBROUTINE CountSectionCells

!***********************************************************************
SUBROUTINE outNB(n,thrs)
!***********************************************************************
!
!  Purpose: To write transport and scalar values at a cross section to a
!           binary file to be used as boundary conditions in nesting
!           procedure.
!
!-----------------------------------------------------------------------

   INTEGER, INTENT(IN) :: n
   REAL, INTENT(IN) :: thrs

   !.....Local variables.....
   CHARACTER (LEN = 12) :: xfile, Ifile
   INTEGER :: nboid, m1, m2, ios, istat, ipts, ko, iboid
   INTEGER :: i , j , k , l, kmx, kmy, js, je, is, ie, icl, lv
   INTEGER :: ii, jj, kk, iis, iie, jjs, jje, kks, kke, nn
   REAL    :: dzi, uflow, uflowi, vflow, vflowi, dflow,dzi2
   REAL, ALLOCATABLE, DIMENSION(:,:) :: outvar
   INTEGER, PARAMETER:: iboid0 = 20

   ! ... Return if interior boundaries are not requested
   IF (nxNBO <= 0 .OR. ioNBO<= 0 .OR. ioNBTOGGLE <= 0) RETURN

   !.....Open files on first entry & write dims. of the problem
   IF( n == 0 ) THEN
     nfrNBO = nts/MAX(ioNBO,1)
     ntrNBO = ntr
     DO nn = 1, nxNBO
       ! ... Find No. cells in nested boundary and set value of iptNBO
       CALL FindCellsNBO(ipts, nn)
       iptNBO(nn) = ipts
       ! ... Open files
       nboid = nboid0 + nn
       iboid = iboid0 + nn
       xfile = "nbofilex0    "
       Ifile = "nbofilei0    "
       IF ( nn < 10 ) WRITE ( xfile(10:12), FMT='(I1,"  ")' ) nn
       IF ( nn>= 10 ) WRITE ( xfile( 9:12), FMT='(I2,"  ")' ) nn
       IF ( nn < 10 ) WRITE ( Ifile(10:12), FMT='(I1,"  ")' ) nn
       IF ( nn>= 10 ) WRITE ( Ifile( 9:12), FMT='(I2,"  ")' ) nn
       OPEN(unit=nboid,file=xfile,FORM='UNFORMATTED',IOSTAT=ios)
       OPEN(unit=iboid,file=Ifile,FORM='FORMATTED'  ,IOSTAT=ios)
       IF(ios /= 0) THEN
         PRINT *, "Error opening xfile = ", nn, ios; STOP
       ENDIF
       !... Write side of nested boundary (i.e. N, S, E or W)
       WRITE(nboid) isdNBO(nn)
       !... Write time information (no. of frames and time between frames)
       WRITE(nboid) nfrNBO(nn)
       !... Write spatial information (no. of cells) & no. of tracers
       WRITE(nboid) iptNBO(nn), ntrNBO(nn)
     END DO
   ENDIF

   ! ... Assing output variables and write to output files
   DO nn = 1, nxNBO

     ! ... Allocate space for output variables
     ipts = iptNBO(nn)
     ALLOCATE( outvar ( ipts, 5+ntr ), STAT=istat )
     IF (istat /= 0) CALL allocate_error (istat,30)

     SELECT CASE (isdNBO(nn))

       ! East or West boundary
       CASE (1,3)

         ! Get i-, j- indexes for bdry. point
         i  = isbcNBO(nn);
         js = jsbcNBO(nn);
         je = jebcNBO(nn)

         ! ... Assign variable values to cells within nested grid
         uflow = 0.0E0
         icl = 0
         DO j = js,je
           ! ... Get l- from (i,j)
           l   = ij2l(i,j)
           lv  = l; IF (isdNBO(nn)== 1) THEN; lv = ij2l(i-1,j); ENDIF
           ! ... Indexes for fine-grid cells in coarse-grid grid cell
           ii  = ((i-i1)+1)*xxNBO+i1-1
           jjs = ((j-j1)  )*xxNBO+j1
           jje = ((j-j1)+1)*xxNBO+j1-1
           kmx = kmz(l)
           DO jj = jjs, jje
             DO k = k1, kmx
               icl = icl + 1
               IF (icl > ipts) THEN
                 PRINT *, ' STOP - Counters in SUBROUTINE outNB do not agree'
                 STOP
               ENDIF
               outvar(icl,1) = FLOAT(ii)
               outvar(icl,2) = FLOAT(jj)
               outvar(icl,3) = FLOAT(k )
               outvar(icl,4) = uh (k,lv)
               outvar(icl,5) = sal(k,l )
               IF (ntr > 0) THEN
                 outvar(icl,6:5+ntr) = tracer(k,l,1:ntr)
               ENDIF
               uflow = uflow + uh(k,lv)
             ENDDO
           ENDDO
         ENDDO
         dzi = 0.0E0
         DO i = i1,isbcNBO(nn);
           DO j = j1, jm
             IF (.NOT. mask2d(i,j)) CYCLE
             ! ... Get l- from (i,j)
             l   = ij2l(i,j)
             ! ... Add water surface displacements
             dzi = dzi + s(l)
           ENDDO
         ENDDO

       ! North or South boundaries
       CASE (2,4)

         ! Get i-, j- indexes for bdry. point
         j  = jsbcNBO(nn);
         is = isbcNBO(nn);
         ie = iebcNBO(nn)

         ! ... Assign variable values to cells within nested grid
         icl = 0; uflow = 0.0E0
         DO i = is,ie
           ! ... Get l- from (i,j)
           l   = ij2l(i,j)
           lv  = l; IF (isdNBO(nn)== 4) THEN; lv = ij2l(i,j-1); ENDIF
           ! ... Indexes for fine-grid cells within coarse-grid grid cell
           jj  = ((j-j1)+1)*xxNBO+j1-1
           iis = ((i-i1)  )*xxNBO+i1
           iie = ((i-i1)+1)*xxNBO+i1-1
           kmy = kmz(l)
           DO ii = iis, iie
             DO k = k1, kmy
               icl = icl + 1
               IF (icl > ipts) THEN
                 PRINT *, ' STOP - Counters in SUBROUTINE outNB do not agree'
                 STOP
               ENDIF
               outvar(icl,1) = FLOAT(ii)
               outvar(icl,2) = FLOAT(jj)
               outvar(icl,3) = FLOAT(k )
               outvar(icl,4) = vh (k,lv)
               outvar(icl,5) = sal(k,l )
               IF (ntr > 0) THEN
                 outvar(icl,6:5+ntr) = tracer(k,l,1:ntr)
               ENDIF
               uflow = uflow + vh(k,lv)
             ENDDO
           ENDDO
         ENDDO

         ! ... Water Surface Elevation --> mass             MAC BEZNAR
         dzi = 0.0E0
         DO j = j1,jsbcNBO(nn);
           DO i = i1, im
             IF (.NOT. mask2d(i,j)) CYCLE
             ! ... Get l- from (i,j)
             l   = ij2l(i,j)
             ! ... Add water surface displacements
             dzi = dzi + s(l)
           ENDDO
         ENDDO

     END SELECT


     ! ... Double check whether icl = ipts
     IF ( icl .NE. ipts) THEN
       PRINT *, ' STOP - Counters in SUBROUTINE outNB do not agree'
       STOP
     ENDIF

     ! ... Id # for nesting boundary file
     nboid = nboid0 + nn
     iboid = iboid0 + nn

     ! ... Write variables to nesting boundary files
     IF((MOD(n,MAX(ioNBO,1)) == 0)) THEN
       IF (n == 0 ) THEN
         WRITE(nboid) thrs,((outvar(m1,m2),m2=1,5+ntr),m1=1,ipts)
       ELSE
         WRITE(nboid) thrs,((outvar(m1,m2),m2=4,5+ntr),m1=1,ipts)
       END IF
     ENDIF

     DEALLOCATE (outvar)


     ! ... Mass Balance Check
     WRITE (UNIT=iboid, FMT='(4E20.11)') thrs, uflow, dzi*dx*dy


   ENDDO


END SUBROUTINE outNB

!***********************************************************************
SUBROUTINE outNBold(n,thrs)
!***********************************************************************
!
!  Purpose: To write transport and scalar values at a cross section to a
!           binary file to be used as boundary conditions in nesting
!           procedure.
!
!-----------------------------------------------------------------------

   INTEGER, INTENT(IN) :: n
   REAL, INTENT(IN) :: thrs

   !.....Local variables.....
   CHARACTER (LEN = 12) :: xfile, Ifile
   INTEGER :: nboid, m1, m2, ios, istat, ipts, ko, iboid
   INTEGER :: i , j , k , l, kmx, kmy, js, je, is, ie, icl
   INTEGER :: ii, jj, kk, iis, iie, jjs, jje, kks, kke, nn
   REAL    :: dzi, uflow, uflowi, vflow, vflowi, dflow,dzi2
   REAL, ALLOCATABLE, DIMENSION(:,:) :: outvar
   INTEGER, PARAMETER:: iboid0 = 20

   ! ... Return if interior boundaries are not requested
   IF (nxNBO <= 0 .OR. ioNBO<= 0 .OR. ioNBTOGGLE <= 0) RETURN

   !.....Open files on first entry & write dims. of the problem
   IF( n == 0 ) THEN
     nfrNBO = nts/MAX(ioNBO,1)
     ntrNBO = ntr
     DO nn = 1, nxNBO
       ! ... Find No. cells in nested boundary and set value of iptNBO
       CALL FindCellsNBO(ipts, nn)
       iptNBO(nn) = ipts
       ! ... Open files
       nboid = nboid0 + nn
       iboid = iboid0 + nn
       xfile = "nbofilex0    "
       Ifile = "nbofilei0    "
       IF ( nn < 10 ) WRITE ( xfile(10:12), FMT='(I1,"  ")' ) nn
       IF ( nn>= 10 ) WRITE ( xfile( 9:12), FMT='(I2,"  ")' ) nn
       IF ( nn < 10 ) WRITE ( Ifile(10:12), FMT='(I1,"  ")' ) nn
       IF ( nn>= 10 ) WRITE ( Ifile( 9:12), FMT='(I2,"  ")' ) nn
       OPEN(unit=nboid,file=xfile,FORM='UNFORMATTED',IOSTAT=ios)
       OPEN(unit=iboid,file=Ifile,FORM='FORMATTED',IOSTAT=ios)

       IF(ios /= 0) THEN
         PRINT *, "Error opening xfile = ", nn, ios; STOP
       ENDIF
       !... Write side of nested boundary (i.e. N, S, E or W)
       WRITE(nboid) isdNBO(nn)
       !... Write time information (no. of frames and time between frames)
       WRITE(nboid) nfrNBO(nn)
       !... Write spatial information (no. of cells) & no. of tracers
       WRITE(nboid) iptNBO(nn), ntrNBO(nn)
     END DO
   ENDIF

   ! ... Assing output variables and write to output files
   DO nn = 1, nxNBO

     ! ... Allocate space for output variables
     ipts = iptNBO(nn)
     ALLOCATE( outvar ( ipts, 5+ntr ), STAT=istat )
     IF (istat /= 0) CALL allocate_error (istat,30)

     SELECT CASE (isdNBO(nn))

       ! East or West boundary
       CASE (1,3)

         ! Get i-, j- indexes for bdry. point
         i  = isbcNBO(nn);
         js = jsbcNBO(nn);
         je = jebcNBO(nn)

         ! ... Assign variable values to cells within nested grid
         uflow = 0.0E0
         icl = 0
         DO j = js,je
           ! ... Get l- from (i,j)
           l   = ij2l(i,j)
           ! ... Indexes for fine-grid cells in coarse-grid grid cell
           ii  = ((i-i1)+1)*xxNBO+i1-1
           jjs = ((j-j1)  )*xxNBO+j1
           jje = ((j-j1)+1)*xxNBO+j1-1
           kmx = kmz(l)
           DO jj = jjs, jje
             DO k = k1, kmx
               icl = icl + 1
               IF (icl > ipts) THEN
                 PRINT *, ' STOP - Counters in SUBROUTINE outNB do not agree'
                 STOP
               ENDIF
               outvar(icl,1) = FLOAT(ii)
               outvar(icl,2) = FLOAT(jj)
               outvar(icl,3) = FLOAT(k )
               outvar(icl,4) = uh (k,l )
               outvar(icl,5) = sal(k,l )
               IF (ntr > 0) THEN
                 outvar(icl,6:5+ntr) = tracer(k,l,1:ntr)
               ENDIF
               uflow = uflow + uh(k,l)
             ENDDO
           ENDDO
         ENDDO
         dzi = 0.0E0
         DO i = i1,isbcNBO(nn);
           DO j = j1, jm
             IF (.NOT. mask2d(i,j)) CYCLE
             ! ... Get l- from (i,j)
             l   = ij2l(i,j)
             ! ... Add water surface displacements
             dzi = dzi + s(l)
           ENDDO
         ENDDO


       ! North or South boundaries
       CASE (2,4)

         ! Get i-, j- indexes for bdry. point
         j  = jsbcNBO(nn);
         is = isbcNBO(nn);
         ie = iebcNBO(nn)

         ! ... Assign variable values to cells within nested grid
         icl = 0; uflow = 0.0E0
         DO i = is,ie
           ! ... Get l- from (i,j)
           l   = ij2l(i,j)
           ! ... Indexes for fine-grid cells within coarse-grid grid cell
           jj  = ((j-j1)+1)*xxNBO+j1-1
           iis = ((i-i1)  )*xxNBO+i1
           iie = ((i-i1)+1)*xxNBO+i1-1
           kmy = kmz(l)
           DO ii = iis, iie
             DO k = k1, kmy
               icl = icl + 1
               IF (icl > ipts) THEN
                 PRINT *, ' STOP - Counters in SUBROUTINE outNB do not agree'
                 STOP
               ENDIF
               outvar(icl,1) = FLOAT(ii)
               outvar(icl,2) = FLOAT(jj)
               outvar(icl,3) = FLOAT(k )
               outvar(icl,4) = vh (k,l )
               outvar(icl,5) = sal(k,l )
               IF (ntr > 0) THEN
                 outvar(icl,6:5+ntr) = tracer(k,l,1:ntr)
               ENDIF
               uflow = uflow + vh(k,l)
             ENDDO
           ENDDO
         ENDDO

         ! ... Water Surface Elevation --> mass
         dzi = 0.0E0
         DO j = j1,jsbcNBO(nn);
           DO i = i1, im
             IF (.NOT. mask2d(i,j)) CYCLE
             ! ... Get l- from (i,j)
             l   = ij2l(i,j)
             ! ... Add water surface displacements
             dzi = dzi + s(l)
           ENDDO
         ENDDO

     END SELECT

         dzi2 = 0.0E0
         DO i = 6,27;
           DO j = 30, 96
             IF (.NOT. mask2d(i,j)) CYCLE
             ! ... Get l- from (i,j)
             l   = ij2l(i,j)
             ! ... Add water surface displacements
             dzi2 = dzi2 + s(l)
           ENDDO
         ENDDO

     ! ... Double check whether icl = ipts
     IF ( icl .NE. ipts) THEN
       PRINT *, ' STOP - Counters in SUBROUTINE outNB do not agree'
       STOP
     ENDIF

     ! ... Id # for nesting boundary file
     nboid = nboid0 + nn
     iboid = iboid0 + nn

     ! ... Write variables to nesting boundary files
     IF((MOD(n,MAX(ioNBO,1)) == 0)) THEN
       IF (n == 0 ) THEN
         WRITE(nboid) thrs,((outvar(m1,m2),m2=1,5+ntr),m1=1,ipts)
       ELSE
         WRITE(nboid) thrs,((outvar(m1,m2),m2=4,5+ntr),m1=1,ipts)
       END IF
     ENDIF

     DEALLOCATE (outvar)

     ! ... Mass Balance Check
     WRITE (UNIT=iboid, FMT='(5E20.11)') thrs, uflow, dzi*dx*dy, dzi2*dx*dy


   ENDDO


END SUBROUTINE outNBold

!***********************************************************************
SUBROUTINE outcheckMass(n,thrs)
!***********************************************************************
!
!  Purpose: To write transport and scalar values at a cross section to a
!           binary file to be used as boundary conditions in nesting
!           procedure.
!
!-----------------------------------------------------------------------

   INTEGER, INTENT(IN) :: n
   REAL, INTENT(IN) :: thrs

   !.....Local variables.....
   CHARACTER (LEN = 12) :: xfile, Ifile
   INTEGER :: nboid, m1, m2, ios, istat, ipts, ko, iboid
   INTEGER :: i , j , k , l, kmx, kmy, js, je, is, ie, icl
   INTEGER :: ii, jj, kk, iis, iie, jjs, jje, kks, kke, nn
   REAL    :: dzi, uflow, uflowi, vflow, vflowi, dflow
   REAL, ALLOCATABLE, DIMENSION(:,:) :: outvar
   INTEGER, PARAMETER:: iboid0 = 20

   ! ... Return if interior boundaries are not requested
   IF (nopen <= 0) RETURN

   !.....Open files on first entry & write dims. of the problem
   IF( n == 0 ) THEN

     DO nn = 1, nopen



       ! ... Open files

       iboid = iboid0 + nn

       Ifile = "nbofilei0    "

       IF ( nn < 10 ) WRITE ( Ifile(10:12), FMT='(I1,"  ")' ) nn
       IF ( nn>= 10 ) WRITE ( Ifile( 9:12), FMT='(I2,"  ")' ) nn
       !print *,"antes de formatear"
       OPEN(unit=iboid,file=Ifile,FORM='FORMATTED',IOSTAT=ios)
       !print *,"despues formt"
       IF(ios /= 0) THEN
         PRINT *, "Error opening xfile = ", nn, ios; STOP
       ENDIF
       !print *,"despues formt2"
       !... Write side of nested boundary (i.e. N, S, E or W)
       !WRITE(iboid) isdNBI(nn)
       !print *,"despues formt3"
       !... Write time information (no. of frames and time between frames)
       !WRITE(iboid) nfrNBI(nn)
       !... Write spatial information (no. of cells) & no. of tracers
       !WRITE(iboid) iptNBI(nn), ntrNBI(nn)
     END DO
   ENDIF
   !print *,"despues de formatear iboid"
   ! ... Assing output variables and write to output files
   DO nn = 1, nopen



     ! ... Allocate space for output variables
     ipts = iptNBI(nn)
     ALLOCATE( outvar ( ipts, 5+ntr ), STAT=istat )
     IF (istat /= 0) CALL allocate_error (istat,30)

     SELECT CASE (isdNBI(nn))

       ! East or West boundary
       CASE (1,3)

         ! Get i-, j- indexes for bdry. point
         i  = isbc(nn);
         js = jsbc(nn);
         je = jebc(nn)

         ! ... Assign variable values to cells within nested grid
         uflow = 0.0E0
         icl = 0
         DO j = js,je
           ! ... Get l- from (i,j)
           l   = ij2l(i,j)
           ! ... Indexes for fine-grid cells in coarse-grid grid cell

           kmx = kmz(l)

             DO k = k1, kmx

               uflow = uflow + uh(k,l)
             ENDDO

         ENDDO
         dzi = 0.0E0
         DO i = i1,im;
           DO j = j1, jm
             IF (.NOT. mask2d(i,j)) CYCLE
             ! ... Get l- from (i,j)
             l   = ij2l(i,j)
             ! ... Add water surface displacements
             dzi = dzi + s(l)
           ENDDO
         ENDDO

       ! North or South boundaries
       CASE (2,4)

         ! Get i-, j- indexes for bdry. point
         j  = jsbc(nn);
         is = isbc(nn);
         ie = iebc(nn)

         ! ... Assign variable values to cells within nested grid
         icl = 0; uflow = 0.0E0
         DO i = is,ie

           ! ... Get l- from (i,j)
           l   = ij2l(i,j)
           ! ... Indexes for fine-grid cells within coarse-grid grid cell

           kmy = kmz(l)

             DO k = k1, kmy

               uflow = uflow + vh(k,l)
             ENDDO

         ENDDO

         ! ... Water Surface Elevation --> mass
         dzi = 0.0E0
         DO j = j1,jm;
           DO i = i1, im
             IF (.NOT. mask2d(i,j)) CYCLE
             ! ... Get l- from (i,j)
             l   = ij2l(i,j)
             ! ... Add water surface displacements
             dzi = dzi + s(l)

           ENDDO
         ENDDO

     END SELECT



     ! ... Id # for nesting boundary file

     iboid = iboid0 + nn



     DEALLOCATE (outvar)

     ! ... Mass Balance Check
     WRITE (UNIT=iboid, FMT='(3E20.11)') thrs, uflow, dzi*dx*dy


   ENDDO


END SUBROUTINE outcheckMass



!***********************************************************************
SUBROUTINE FindCellsNBO ( ix, nn )
!***********************************************************************
!
!  Purpose: To count wett cells ix in a nested grid within a given
!           X-section nn
!
!-----------------------------------------------------------------------

 ! ... Arguments
 INTEGER, INTENT (inout) :: ix
 INTEGER, INTENT (in)    :: nn

 ! ... Local variables
 INTEGER :: count, i, j, k, l, kmx, kmy, js, je, is, ie

 SELECT CASE (isdNBO(nn))

   ! East or West boundary
   CASE (1,3)

     ! Get i-, j- indexes for bdry. point
     i  = isbcNBO(nn);
     js = jsbcNBO(nn);
     je = jebcNBO(nn);
     print *,"i:",i,"js:",js,"je",je
     ! Make sure i,j locations are wett cells
     DO j = js, je
       IF ( .NOT. mask2d(i,j) ) THEN
         PRINT *, "  "
         PRINT *, " ****STOP - EW nesting bdry. DRY"
         STOP
       END IF
     ENDDO

     ! ... Count of No.cells at boundary of nested grid
     count = 0;
     DO j = js,je
       kmx = kmz(ij2l(i,j))
       print *,"i:",i,"j:",j,"l:",ij2l(i,j),"kmx:",kmx
       count = count + (kmx-k1+1)
       print *,"count:",count
     ENDDO
     count = count * xxNBO
     print *,"xxNBO:",xxNBO
   ! North or South boundary
   CASE (2,4)

     ! ... Get i-, j- indexes for bdry. point
     j  = jsbcNBO(nn);
     is = isbcNBO(nn);
     ie = iebcNBO(nn);
     print *,"lalalalala"
     ! ... Make sure i,j locations are wett cells
     DO i = is, ie
       IF ( .NOT. mask2d(i,j) ) THEN
         PRINT *, "  "
         PRINT *, " ****STOP - NS nesting bdry. DRY "
         STOP
       END IF
     ENDDO

     ! ... Count of No.cells at boundary of nested grid
     count = 0;
     DO i = is, ie
       kmy = kmz(ij2l(i,j))
       count = count + (kmy-k1+1)
     ENDDO
     count = count * xxNBO

 END SELECT
 ix = count


END SUBROUTINE FindCellsNBO

!***********************************************************************
SUBROUTINE outh(n)
!***********************************************************************
!
!  Purpose: To write output at a specific layer in binary format
!
!-----------------------------------------------------------------------

   INTEGER, INTENT(IN) :: n

   !.....Local variables.....
   CHARACTER (LEN = 10) :: plane_file
   INTEGER, PARAMETER   :: plane_id0 = 800
   INTEGER :: plane_id
   INTEGER :: i, j, k, l, ios, istat, n_frames, ipoints , k_out, m1,m2, kp, kb, c
!   REAL, ALLOCATABLE, DIMENSION(:,:) :: out_array
   INTEGER :: year_out, day_out, mon_out
   REAL    :: hour_out
   print *,"a,"
   !.....Open spacefile on first entry and print initial conditions ....
   IF( n == 0 ) THEN

     ! ... Determine No. of time slices to output
     n_frames = nts/MAX(iop,1)
     print *,"aa,"
     DO j = 1, n_planes
       ! ... Check that plane no. is below 2
       IF ( p_out(j) < k1 ) THEN
          PRINT *, 'ERROR: Output plane is not in DOMAIN'
          STOP
       END IF
       ! ... Open file to output solution
       plane_id = plane_id0 + j
       plane_file = "plane_    "
       IF ( p_out(j) < 10 ) WRITE ( plane_file(7:10), FMT='(I1,"   ")' ) p_out(j)
       IF(( p_out(j) >= 10 ) .AND. ( p_out(j) < 100)) &
       &                    WRITE ( plane_file(7:10), FMT='(I2,"  " )' ) p_out(j)
       IF ( p_out(j) > 100) WRITE ( plane_file(7:10), FMT='(I3," "  )' ) p_out(j)
       OPEN(unit=plane_id,file=plane_file,FORM='UNFORMATTED',IOSTAT=ios)
       IF(ios /= 0) PRINT *, "Error opening plane file = ", p_out(j), ios
       !... Write number of time slices to output file
       WRITE(plane_id) n_frames
       !... Find & write number of wett cells in the plane
       ipoints = 0;
       CALL CountPlaneCells ( ipoints, p_out(j) )
       interior_plane_points (j) = ipoints
       WRITE(plane_id) ipoints
     END DO
     print *,"bb,"
     ! ... Time stamp
     year_out = iyr
     mon_out  = imon
     day_out  = iday
     hour_out = ihr

     ! ... Output planes
     DO k = 1, n_planes
       print *,"ff"
       ipoints = interior_plane_points (k)
!       ALLOCATE( out_array ( ipoints, 8 ), STAT=istat )
!       IF (istat /= 0) THEN
!         PRINT *, 'ERROR allocating space for out_array'
!         STOP
!       ENDIF
       kp = p_out(k)
       k_out = 0
       print *,"fff"
       ! ... Bottom cell plane
       IF ( kp > km ) THEN
         DO c=1,cm1
           ! ... Ignore dry cells
           IF (.NOT. mask(c)) CYCLE
           ! ... Determine k- & l- indexes for output cell

           l  = l2c(c)
           kb = kmz (l)
           ! ... Define output variables
           k_out = k_out + 1
           out_array(k_out,1) = FLOAT(l2i(l))
           out_array(k_out,2) = FLOAT(l2j(l))
           out_array(k_out,3) = up (kb,l)
           out_array(k_out,4) = vp (kb,l)
           out_array(k_out,5) = wp (kb,l)
           out_array(k_out,6) = salp(kb,l)
           IF(hp(kb,l)<=ZERO) out_array(k_out,6) = -99.
           out_array(k_out,7) = 0.5*(Av(kb,l)+Av(kb+1,l))
           !out_array(k_out,8) = 0.5*(Dv(kb,l)+Dv(kb+1,l))
           out_array(k_out,8) = s(l) ! Changed 12/2010 SWA
           !out_array(k_out,8) = hhs(i,j)
         END DO
       print *,"bc,"
       ! ... Interior plane
       ELSE
         DO c=1,cm1
           ! Ignore dry cells
           IF (.NOT. mask(c)) CYCLE
           IF ( kmz(l2c(c)) < kp   ) CYCLE
           ! ... Determine k- & l- indexes for output cell
           kb = kp
           l  = l2c(c)
           ! ... Define output variables
           k_out = k_out + 1
           out_array(k_out,1) = FLOAT(l2i(l))
           out_array(k_out,2) = FLOAT(l2j(l))
           out_array(k_out,3) = up (kb,l)
           out_array(k_out,4) = vp (kb,l)
           out_array(k_out,5) = wp (kb,l)
           out_array(k_out,6) = salp(kb,l)
           IF(hp(kb,l)<=ZERO) out_array(k_out,6) = -99.
           out_array(k_out,7) = 0.5*(Av(kb,l)+Av(kb+1,l))
           !out_array(k_out,8) = 0.5*(Dv(kb,l)+Dv(kb+1,l))
           out_array(k_out,8) = s(l) ! Changed 12/2010 SWA
           !out_array(k_out,8) = hhs(i,j)
         END DO
       END IF

       ! ... Id # for plane file
       plane_id = plane_id0 + k
        print *,"bh,"
       ! ... Print time stamp followed by the records
       WRITE(plane_id) n,year_out,mon_out,day_out,hour_out,          &
       &            ((out_array(m1,m2),m2=1,8),m1=1,ipoints)
!       DEALLOCATE (out_array)
        print *,"j"
     END DO
     print *,"c,"
   ELSE

     ! ... Time stamp
     year_out = iyr
     mon_out  = imon
     day_out  = iday
     hour_out = ihr

     ! ... Output planes
     DO k = 1, n_planes

       ipoints = interior_plane_points (k)
!       ALLOCATE( out_array ( ipoints, 8 ), STAT=istat )
!       IF (istat /= 0) THEN
!         PRINT *, 'ERROR allocating space for out_array'
!         STOP
!       ENDIF
       kp = p_out(k)
       k_out = 0

       ! ... Bottom cell plane
       IF ( kp > km ) THEN
         DO c=1,cm1
           ! ... Ignore dry cells
           IF (.NOT. mask(c)) CYCLE
           ! ... Determine k- & l- indexes for output cell

           l  = l2c(c)
           kb = kmz (l)
           ! ... Define output variables
           k_out = k_out + 1
           out_array(k_out,3) = (up (kb,l) + up(kb,lWC(l)))/2.
           out_array(k_out,4) = (vp (kb,l) + vp(kb,lSC(l)))/2.
           out_array(k_out,5) = (wp (kb,l) + wp(kb+1,  l ))/2.
           out_array(k_out,6) = salp(kb,l)
           IF(hp(kb,l)<=ZERO) out_array(k_out,6) = -99.
           out_array(k_out,7) = 0.5*(Av(kb,l)+Av(kb+1,l))
           !out_array(k_out,8) = 0.5*(Dv(kb,l)+Dv(kb+1,l))
           out_array(k_out,8) = s(l) ! Changed 12/2010 SWA
           !out_array(k_out,8) = hhs(i,j)
         END DO
       ! ... Interior plane
       ELSE
         DO c=1,cm1
           ! Ignore dry cells
           IF (.NOT. mask(c)) CYCLE
           IF ( kmz(l2c(c)) < kp   ) CYCLE
           ! ... Determine k- & l- indexes for output cell
           kb = kp
           l  = l2c(c)
           ! ... Define output variables
           k_out = k_out + 1
           out_array(k_out,3) = (up (kb,l) + up(kb,lWC(l)))/2.
           out_array(k_out,4) = (vp (kb,l) + vp(kb,lSC(l)))/2.
           out_array(k_out,5) = (wp (kb,l) + wp(kb+1,  l ))/2.
           out_array(k_out,6) = salp(kb,l);
           IF(hp(kb,l)<=ZERO) out_array(k_out,6) = -99.
           out_array(k_out,7) = 0.5*(Av(kb,l)+Av(kb+1,l))
           !out_array(k_out,8) = 0.5*(Dv(kb,l)+Dv(kb+1,l))
           out_array(k_out,8) = s(l) ! Changed 12/2010 SWA
           !out_array(k_out,8) = hhs(i,j)
         END DO
       END IF

       ! ... Id # for plane file
       plane_id = plane_id0 + k
       ! ... Print time stamp followed by the records
       WRITE(plane_id) n,year_out,mon_out,day_out, hour_out,         &
       &            ((out_array(m1,m2),m2=3,8),m1=1,ipoints)
!       DEALLOCATE (out_array)
     END DO
   END IF

END SUBROUTINE outh



!***********************************************************************
SUBROUTINE CountPlaneCells ( count, kplane )
!***********************************************************************
!
!  Purpose: To count wet points in a given kplane layer
!
!-----------------------------------------------------------------------

 ! ... Arguments
 INTEGER, INTENT (inout) :: count
 INTEGER, INTENT (in)    :: kplane

 ! ... Local variables
 INTEGER :: i, j, c

 ! ... Code
 count = 0
 IF ( km < kplane ) THEN
   DO c=1,cm1
     ! Ignore dry cells
     IF (.NOT. mask(c)  ) CYCLE
       count = count + 1
   END DO
 ELSE
   DO c=1,cm1
     ! Ignore dry cells
     IF (.NOT. mask(c)  ) CYCLE
     IF ( kmz(l2c(c)) < kplane ) CYCLE;
       count = count + 1
   END DO
 END IF

END SUBROUTINE CountPlaneCells



!***********************************************************************
SUBROUTINE outz(n)
!***********************************************************************
!
!  Purpose: To write tracer concentrations in computational domain
!
!-----------------------------------------------------------------------

   INTEGER, INTENT(IN) :: n

   !.....Local variables.....
   CHARACTER (LEN = 11) :: tracer_file
   CHARACTER (LEN = 13) :: tracerbc_file
   INTEGER, PARAMETER   :: tracer_id0 = 1200
   INTEGER, PARAMETER   :: tracerbc_id0 = 1400
   INTEGER :: tracer_id, tracerbc_id
   INTEGER :: i, j, k, l, k_t,ios, istat, n_frames,  k_out, m1,m2, c
   INTEGER, SAVE :: ipoints
   !REAL, ALLOCATABLE, DIMENSION(:,:) :: out_arrayZ
   INTEGER :: year_out, day_out, mon_out
   REAL    :: hour_out
   INTEGER:: jj
   CHARACTER(LEN=24) :: fmt

   IF ( ntr <= 0 ) RETURN

   !.....Open spacefile on first entry and print initial conditions ....
   IF( n == 0 ) THEN

     ipoints = 0
     DO l = 1, lm
        DO k = k1, kmz(l)
           ipoints = ipoints + 1;
        ENDDO
     ENDDO

     DO j = 1, ntr
       n_frames = nts/MAX(iotr,1)
       tracer_id = tracer_id0 + j
       tracer_file = "tracer_    "
       IF ( j < 10 ) WRITE ( tracer_file(8:11), FMT='(I1,"   ")' ) j
       IF ((j >=10 ) .AND. (j < 100)) WRITE ( tracer_file(8:11), FMT='(I2,"  ")' ) j
       IF ( j >100 ) WRITE ( tracer_file(8:11), FMT='(I3," ")' ) j
       OPEN(unit=tracer_id,file=tracer_file,FORM='UNFORMATTED',IOSTAT=ios)
       IF(ios /= 0) THEN; PRINT *, "Error opening tracer file = ", j, ios;STOP;ENDIF
       !... Write number of time slices to output file
       WRITE(tracer_id) n_frames
       WRITE(tracer_id) ipoints
     END DO

     ! ... Time stamp
     year_out = iyr
     mon_out  = imon
     day_out  = iday
     hour_out = ihr

     !ALLOCATE( out_arrayZ ( ipoints, 4 ), STAT=istat )
     !IF (istat /= 0) THEN;
     !  PRINT *, 'ERROR allocating space in output_tracer'
     !  STOP
     !ENDIF

     ! ... Output tracer concentrations
     DO k_t = 1, ntr
        k_out = 0
        ! ... Assign values to the output array
        DO c=1,cm1
             IF (.NOT. mask(c)) CYCLE
             l = l2c(c)
             DO k = k1, kmz(l)
             k_out = k_out + 1
             out_array(k_out,1) = FLOAT(l2i(l))
             out_array(k_out,2) = FLOAT(l2j(l))
             out_array(k_out,3) = FLOAT(k)
             out_array(k_out,4) = tracer(k,l,k_t) * h(k,l)
             out_array(k_out,5) = tracer(k,l,k_t)! cintia_trazador

             END DO;
        END DO
        ! ... Id # for plane file
        tracer_id = tracer_id0 + k_t
        ! ... Print time stamp followed by the records
        WRITE(tracer_id) n,year_out,mon_out, day_out,hour_out,  &
       &            ((out_array(m1,m2),m2=1,5),m1=1,ipoints)
     END DO
     !DEALLOCATE (out_arrayZ)

   ELSE

     ! ... Time stamp
     year_out = iyr
     mon_out  = imon
     day_out  = iday
     hour_out = ihr

     ! ... Allocate space
     !ALLOCATE( out_arrayZ ( ipoints, 1 ), STAT=istat )
     !IF (istat /= 0) THEN;
     !  PRINT *, 'ERROR allocating space in output_tracer'
     !  STOP
     !ENDIF

     ! ... Output tracer concentrations
     DO k_t = 1, ntr
        k_out = 0
        ! ... Assign values to the output array
        DO c=1,cm1
             IF (.NOT. mask(c)) CYCLE
             l = l2c(c)
             DO k = k1, kmz(l)
             k_out = k_out + 1
             out_array(k_out,1) = tracer(k,l,k_t) * h(k,l)
             out_array(k_out,2) = tracer(k,l,k_t)! cintia_trazador
             END DO;
        END DO
        ! ... Id # for plane file
        tracer_id = tracer_id0 + k_t
        ! ... Print time stamp followed by the records
        WRITE(tracer_id) n,year_out,mon_out, day_out,hour_out,   &
        &            ((out_array(m1,m2),m2=1,2),m1=1,ipoints)
     END DO
     !DEALLOCATE (out_arrayZ)

   END IF

   10 FORMAT( 1X,2I7,F10.3,F40.10)

END SUBROUTINE outz



!***********************************************************************
SUBROUTINE outp(n)
!***********************************************************************
!
!  Purpose: To write the complete solution in the computational domain
!           The resulting file is used to drive particle tracking
!           simulations with PTRACK-TOOL.
!-----------------------------------------------------------------------

   INTEGER, INTENT(IN) :: n

   !.....Local variables.....
   CHARACTER (LEN = 16) :: ptrack_file
   INTEGER, PARAMETER   :: ptrack_id = 1002
   INTEGER              :: i, j, k, l, ios, istat, Noframes, m1, m2,  &
                           kout, year_out, day_out, mon_out, c
   REAL                 :: hour_out
   INTEGER, SAVE        :: ipoints
!   REAL, ALLOCATABLE, DIMENSION(:,:) :: out_array

   IF( n == 0 ) THEN

     ! ... Determine No. of frames to output & No. of interior points
     Noframes = nts/MAX(apxml,1)
     ipoints = 0
     DO l = 1, lm
        i = l2i(l); j = l2j(l)
        DO k = k1, kmz(l)
           ipoints = ipoints + 1;
        ENDDO
     ENDDO

     ! ... Open output file & print data & initial conditions
     ptrack_file = "ptrack_hydro.bnr"
     OPEN(unit=ptrack_id,file=ptrack_file,FORM='UNFORMATTED',IOSTAT=ios)
     IF(ios /= 0) THEN
       PRINT *, "Error opening ptrack hydro file = ", ios
       STOP
     ENDIF
     !... Write number of time slices to output file
     WRITE(ptrack_id) Noframes
     WRITE(ptrack_id) ipoints

     ! ... Time stamp
     year_out = iyr
     mon_out  = imon
     day_out  = iday
     hour_out = ihr

!     ALLOCATE( out_array ( ipoints, 9 ), STAT=istat )
!     IF (istat /= 0) THEN;
!       PRINT *, 'ERROR allocating space in output routine for PTRACK'
!       STOP
!     ENDIF

     ! ... Output tracer concentrations
     kout = 0
     ! ... Assign values to the output array
     DO c=1,cm1
       IF (.NOT. mask(c)) CYCLE
       l = l2c(c)
       i = l2i(l); j = l2j(l); !Andrea PT
       DO k = k1, kmz(l)
         kout = kout + 1

         out_array(kout,1) = FLOAT(i)
         out_array(kout,2) = FLOAT(j)
         out_array(kout,3) = FLOAT(k)
         out_array(kout,4) = hp (k,l)
         out_array(kout,5) = up (k,l)
         out_array(kout,6) = vp (k,l)
         out_array(kout,7) = wp (k,l)
         out_array(kout,8) = Dv (k,l)
         !out_array(kout,8) = Dvm (k,l) !Andrea Ptrack.Cambio Dv por Dvm
         out_array(kout,9) = sal (k,l)
         out_array(kout,10) = q2p (k,l)
         out_array(kout,11) = q2lp (k,l)
         out_array(kout,12) = kh (k,l)
         out_array(kout,13) = Av (k,l)
       END DO;
     END DO
     ! ... Print time stamp followed by the records
     WRITE(ptrack_id) n,year_out,mon_out, day_out,hour_out,  &
     &              ((out_array(m1,m2),m2=1,13),m1=1,ipoints)
!     DEALLOCATE (out_array)

 ELSE

     ! ... Time stamp
     year_out = iyr
     mon_out  = imon
     day_out  = iday
     hour_out = ihr

!     ALLOCATE( out_array ( ipoints, 10 ), STAT=istat )
!     IF (istat /= 0) THEN;
!       PRINT *, 'ERROR allocating space in output routine for PTRACK'
!       STOP
!     ENDIF

     ! ... Output tracer concentrations
     kout = 0
     ! ... Assign values to the output array
     DO c=1,cm1
       IF (.NOT. mask(c)) CYCLE
       l = l2c(c)
       DO k = k1, kmz(l)
         kout = kout + 1
         out_array(kout,1) = hp(k,l)
         out_array(kout,2) = up(k,l)
         out_array(kout,3) = vp(k,l)
         out_array(kout,4) = wp(k,l)
         out_array(kout,5) = Dv(k,l)
         !out_array(kout,5) = Dvm(k,l) !Andrea Ptrack
         out_array(kout,6) = sal (k,l)
         out_array(kout,7) = q2p (k,l)
         out_array(kout,8) = q2lp (k,l)
         out_array(kout,9) = kh (k,l)
         out_array(kout,10) = Av (k,l)
       END DO;
     END DO
     ! ... Print time stamp followed by the records
     WRITE(ptrack_id) n,year_out,mon_out, day_out,hour_out,  &
     &              ((out_array(m1,m2),m2=1,10),m1=1,ipoints)
     !DEALLOCATE (out_array)
     !Dvm = 0.0   !Andrea Ptrack
 END IF

END SUBROUTINE outp

!***********************************************************************
SUBROUTINE AverageOutp(n)   ! subroutine  added by ABH   marked 4 CINTIA
!***********************************************************************
!
!  Purpose: To write ascii output in xml format to a file used for
!           velocity and particle-tracking animations with the Gr
!           application. The file, called 'spacefile.xml', is
!           essentially a header file for the sequential binary files
!           (spacefile3d.bin  and  spacefile2d.bin) written out in
!           SUB outs_bin. The binary files contain the 2d and 3d data
!           from all the wet spatial nodes in the solution at
!           snapshots in time.
!
!-----------------------------------------------------------------------

   INTEGER, INTENT(IN) :: n

   Dvm = Dvm + Dv/100

  ! PRINT*,'Dvm',Dvm(20,10000)



END SUBROUTINE AverageOutp



!***********************************************************************
SUBROUTINE outs(n,its)
!***********************************************************************
!
!  Purpose: To write ascii output in xml format to a file used for
!           velocity and particle-tracking animations with the Gr
!           application. The file, called 'spacefile.xml', is
!           essentially a header file for the sequential binary files
!           (spacefile3d.bin  and  spacefile2d.bin) written out in
!           SUB outs_bin. The binary files contain the 2d and 3d data
!           from all the wet spatial nodes in the solution at
!           snapshots in time.
!
!-----------------------------------------------------------------------

   INTEGER, INTENT(IN) :: n
   REAL, INTENT(IN) :: its


   !.....Local variables.....
   INTEGER :: i, j, k, ios, numdim_2d, numdim_3d, idt2, ihr2, imin2, isec2
   INTEGER :: int_zeta, int_u, int_v, int_uh, int_vh, int_w, int_dz
   INTEGER, SAVE :: itspf1
   REAL(KIND=DPV) :: thrs_dpv, tdays_dpv
   CHARACTER, SAVE :: date*8, time*10, zone*5
   CHARACTER(LEN=10) :: ch_idt, ch_n, ch_int_zeta, ch_int_u, ch_int_v
   CHARACTER(LEN=10) :: ch_int_w, ch_int_dz, ch_numdim_2d, ch_numdim_3d
   CHARACTER(LEN=20) :: ch_thrs_dpv, ch_tdays_dpv
   CHARACTER(LEN=13) :: output_file_xml  = "spacefile.xml"
   CHARACTER(LEN=15) :: output_file2d_bin= "spacefile2d.bin"
   CHARACTER(LEN=15) :: output_file3d_bin= "spacefile3d.bin"
   INTEGER, DIMENSION(8), SAVE :: values
   LOGICAL, SAVE :: first_entry = .TRUE.

   !               -----First entry into subroutine-----

    IF ( first_entry ) THEN

      ! ... Only allowed if layers are of uniform thickness
      IF ( ibathyf < 0 ) THEN
        PRINT *, 'STOP - Output to xml files only allowed for ibathyf > 0'
        PRINT *, '*******  Please revise your input files   *************'
        STOP
      ENDIF

      !.....Open the xml spacefile.....
      first_entry = .FALSE.
      OPEN ( UNIT=i96, FILE=output_file_xml, IOSTAT=ios )
      IF(ios /= 0) CALL open_error ( "Error opening "//output_file_xml, ios )

      !.....Open the binary spacefiles.....
      OPEN ( UNIT=i97, FILE=output_file2d_bin, FORM='UNFORMATTED',        &
	      &   ACCESS='SEQUENTIAL', IOSTAT=ios )
      IF(ios /= 0) CALL open_error ( "Error opening "//output_file2d_bin, &
         & ios )
      OPEN ( UNIT=i98, FILE=output_file3d_bin, FORM='UNFORMATTED',        &
	      &   ACCESS='SEQUENTIAL', IOSTAT=ios )
      IF(ios /= 0) CALL open_error ( "Error opening "//output_file3d_bin, &
         & ios )

      !.....Get date and time of run for insertion into spacefile.xml.....
      CALL date_and_time ( date, time, zone, values )

      !.....Output titles and run parameters into header records of the xml spacefile.....
      WRITE (UNIT=i96, FMT='("<?xml version=""1.0""?>")' )
      WRITE (UNIT=i96, FMT='(1X)' )
      WRITE (UNIT=i96, FMT='("<!-- Si3D Output for Gr -->")')
      WRITE (UNIT=i96, FMT='("<!-- xml header file for binary spacefiles -->")')
      WRITE (UNIT=i96, FMT=1) idt,ddz,im,jm,itrsca,itrap,cd,ismooth,beta,    &
         &                    niter,itrmom,f,wa
      ! idt,ddz,im,jm,iexplt,itrap,cd,ismooth,beta,    &
      !   &                    niter,iextrp,f,tramp  --->
      ! idt,ddz,im,jm,itrsca,itrap,cd,ismooth,beta,    &
      !   &                    niter,itrmom,f,wa
    !1 FORMAT ( "<!-- ", "Run Parameters:  dt =",I5," s", "  dz =",F5.2, " m",& ! idt real
    1 FORMAT("<!-- ", "Run Parameters:  dt =",F5.2," s", "  dz =",F5.2, " m",& ! idt real
         &     "  im =", I5, "  jm =", I5,/"     itrsca =", I2, "  itrap =", &
         &     I2, "  cd = ", F7.4, "  ismooth =", I2, "  beta =", F6.3,     &
         &     "  niter =", I2/"     itrmom =", I2, "  f =", F7.4,           &
         &     "  wa =", F8.1, " -->" )
      WRITE (UNIT=i96, FMT='(1X)' )
      WRITE (UNIT=i96, FMT='("<gov.usgs.gr>")' )
      WRITE (UNIT=i96, FMT='("   <obj class=""gov.usgs.sfhydro.si3d.Si3dOutput""")')
      WRITE (UNIT=i96, FMT='("        runTitle=""",      A, """" )' ) TRIM(title)
      WRITE (UNIT=i96, FMT='("        runNumber=""", A8,A4, """" )' ) date,time(1:4)
      ! Compute the date and time at which data is first written to the spacefiles
      idt2 = ipxml*idt
      IF (itspf == 0 ) THEN
         itspf1 = itspf
         !CALL compute_date ( 0 ) ! idt real
         CALL compute_date (0.0) ! idt real
      ELSE
         itspf1 = INT(itspf/idt2) * idt2
         IF (MOD(itspf,idt2) /= 0) itspf1 = itspf1 + idt2
         !CALL compute_date ( itspf1 ) ! idt real
         CALL compute_date (FLOAT(itspf1)) ! idt real
      END IF
      ihr2  = INT(isec/3600); isec2 = isec - (ihr2*3600)
      imin2 = INT(isec2/60)
      isec2 = isec2 - (imin2*60)
      ch_idt = int_to_char ((idt2),0)
      WRITE (UNIT=i96, FMT='("        startTime=""",I4,"/",I2.2,"/",I2.2," ",I2.2,    &
         &  ":",I2.2,":",I2.2,""""," timeStepInSeconds=",A,">")' ) iyr, imon, iday,   &
         &  ihr2, imin2, isec2, TRIM(ch_idt)
      CALL compute_date ( 0.0 )     ! reset date to start time of run - idt real
      !CALL compute_date ( 0 )     ! reset date to start time of run ! idt real

      !.....Output header records for 2-D variables into the xml spacefile.....
      WRITE (UNIT=i96, FMT='(1X)' )
      WRITE (UNIT=i96, FMT='("      <obj class=""gov.usgs.gr.GrObjectContainer""", &
         &                          " title=""2-D x-y Variables"">")' )
      WRITE (UNIT=i96, FMT='("         <dataseries class=""gov.usgs.gr.data.RafData&
         &Series""")' )
      WRITE (UNIT=i96, FMT='("            file=""", A, """")' ) output_file2d_bin
      numdim_2d = 1   ! For the time being, hardwire for one 2D variable only (Zeta)
      ch_numdim_2d = int_to_char (numdim_2d,0)
      WRITE (UNIT=i96, FMT='("            title=""2-D Variables""",                &
         &                                " numDimensions=", A, ">")' )            &
         &                                TRIM(ch_numdim_2d)
      WRITE (UNIT=i96, FMT='("            <dim num=""0"" title=""Zeta"" units=""met&
         &ers"" unitScale=""0.0001""/>")' )
      WRITE (UNIT=i96, FMT='(1X)' )
      WRITE (UNIT=i96, FMT='("         </dataseries>")' )
      WRITE (UNIT=i96, FMT='("      </obj>")' )

      !.....Output header records for 3-D variables into the xml spacefile.....
      WRITE (UNIT=i96, FMT='(1X)' )
      WRITE (UNIT=i96, FMT='("      <obj class=""gov.usgs.gr.GrObjectContainer""", &
         &                          " title=""3-D Variables"">")' )
      WRITE (UNIT=i96, FMT='("         <dataseries class=""gov.usgs.gr.data.RafData&
         &Series""")' )
      WRITE (UNIT=i96, FMT='("            file=""", A, """")' ) output_file3d_bin
      numdim_3d = 4   ! For the time being, hardwire for four 3D variables only
      ch_numdim_3d = int_to_char (numdim_3d,0)
      WRITE (UNIT=i96, FMT='("            title=""3-D Variables""",                &
         &                                " numDimensions=", A, ">")' )            &
         &                                TRIM(ch_numdim_3d)
      WRITE(UNIT=i96, FMT='("            <dim num=""0"" title=""u"" ",  &
         & "units=""m/s"" unitScale=""0.0001""/>"    )' )
      WRITE(UNIT=i96, FMT='("            <dim num=""1"" title=""v"" ",  &
         & "units=""m/s"" unitScale=""0.0001""/>"    )' )
      WRITE (UNIT=i96, FMT='("            <dim num=""2"" title=""w"" ",  &
         & "units=""m/s"" unitScale=""0.00001""/>"  )' )
      WRITE (UNIT=i96, FMT='("            <dim num=""3"" title=""dz"" ", &
         & "units=""m**2/s"" unitScale=""0.0001""/>" )' )
      WRITE (UNIT=i96, FMT='(1X)' )

      !.....Output closing xml element end tags to spacefile.xml.....
      WRITE (UNIT=i96, FMT='("         </dataseries>")' )
      WRITE (UNIT=i96, FMT='("      </obj>")' )
      WRITE (UNIT=i96, FMT='("   </obj>")' )
      WRITE (UNIT=i96, FMT='("</gov.usgs.gr>")' )
      WRITE (UNIT=i96, FMT='(1X)' )

   END IF

END SUBROUTINE outs




!***********************************************************************
SUBROUTINE outs_bin
!***********************************************************************
!
!  Purpose: To write output in binary format to two spacefiles for velocity
!           and particle-tracking animations with the Gr application. The
!           two files are: spacefile3d.bin and spacefile2d.bin.  3-D
!           variables are stored in the '3d' file and 2-d variables in
!           the '2d' file.  The ordering of nodes written to the file
!           must agree with the  si3d_bathy.xml  file written out in
!           SUB outg.  The velocities written to the file are taken from
!           the faces of the grid cell control volumes.  The two unformatted
!           files written to in this subroutine are sequential. They are
!           written as 'streaming' data (no record breaks in the file) using
!           the 'little-endian' bit order for the storage layout. The data
!           values themselves are stored as 2-byte integers. Each
!           call to the subroutine outputs one time step of data to
!           the spacefiles.  The CALL and file OPEN statements for this
!           subroutine are within SUB outs_xml.
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: i, j, k, k1s, kms, l, c
   INTEGER(KIND=INT2) :: int2_zeta, int2_u, int2_v, int2_w, int2_dz

   ! ... Only allowed if layers are of uniform thickness
   IF ( ibathyf < 0 ) THEN
     PRINT *, 'STOP - Output to xml files only allowed for ibathyf > 0'
     PRINT *, '*******  Please revise your input files   *************'
     STOP
   ENDIF

   !.....Output zeta values at time step  n  to spacefile2d.bin.....
   ! Loop over all nodes
   DO c=1,cm1

      ! Ignore dry cells
      IF (.NOT. mask(c)) CYCLE
      ! Convert zeta to 2-byte integer in meters*10000.
      int2_zeta = NINT(s(l2c(c))*1.E4, int2)
      WRITE (UNIT=i97) int2_zeta

   END DO

   !.....Output velocity values at time step  n  to spacefile3d.bin.....
   ! Loop over all nodes
   DO c=1,cm1

         ! Ignore dry cells
         IF (.NOT. mask(c)) CYCLE
         l   = l2c(c)
         k1s = k1z (l)
         kms = kmz (l)
         DO k = k1s, kms
            ! Convert u,v to 2-byte integers in cm/s*100
            int2_u  = NINT(up (k,l)*1.E4, int2)
            int2_v  = NINT(vp (k,l)*1.E4, int2)
            ! Convert w to a 2-byte integer in cm/s*1000
            int2_w  = NINT(wp (k,l)*1.E5, int2)
            ! Convert dz to a 2-byte integer in cm**2/s
            int2_dz = NINT(Dv(k,l)*1.E4, int2)
            WRITE (UNIT=i98) int2_u, int2_v, int2_w, int2_dz

      END DO
   END DO

END SUBROUTINE outs_bin

!***********************************************************************
SUBROUTINE outg ( h4 )
!***********************************************************************
!
!  Purpose: To process and output bathymetry to a file for use in graphics
!           post-processing and particle tracking. The bathymetry output by
!           this subroutine is processed from the original bathymetry
!           read into the program in SUB bathy.  The processed bathymetry
!           is similar to that processed in SUB bathy for the actual
!           model calculations except that depths are defined at the
!           corners of each grid cell rather than the mid-sides. (For
!           graphics and particle tracking it is necessary to have
!           the bathymetry defined at corners.)  The processed bathymetry
!           is stored in the pointer array h44(1:im1,1:jm1,4).  The corner
!           points of each (i,j) cell are stored in h44 with the following
!           numbering system:
!                               4      3
!                                *----*
!                                |    |
!                                *----*
!                               1      2
!           This subroutine is called from within SUB bathy.
!
!  Dummy argument:
!  h4 = Bathymetry pointer array to be processed. A single depth is defined
!       at each cell corner. Prior to being passed into this subroutine,
!       the array h4 has already been modified in SUB bathy (from the
!       original model input bathymetry array) by correcting the datum
!       and adding a fictitious row/column of depths around the grid. The
!       dimensions of h4 are (0:im1,0:jm1).
!
!-----------------------------------------------------------------------

   !.....Argument.....
   REAL, DIMENSION(:,:), POINTER :: h4

   !.....Local variables.....
   REAL, DIMENSION(:,:,:), POINTER :: h44
   CHARACTER(LEN=14) :: output_file="              "
   CHARACTER :: date*8, time*10, zone*5, ch_itype*5, ch_ikind*3
   CHARACTER(LEN=10) :: ch_idx, ch_idy, ch_idz, ch_nopen
   !CHARACTER(LEN=10):: ch_nbarr
   CHARACTER(LEN=10) :: ch_isbc, ch_jsbc, ch_iebc, ch_jebc
   !CHARACTER(LEN=10):: ch_isbarr, ch_jsbarr, ch_iebarr, ch_jebarr
   CHARACTER(LEN=14) :: ch_xglobal, ch_yglobal, ch_zglobal, ch_rotation
   INTEGER, DIMENSION(8) :: values
   INTEGER :: ios, istat, i, j, k, kb, niters, niterations, nn, ixml
   INTEGER, SAVE :: i90
   REAL :: ztop, hmin, hmax, hmax_e, hmax_n, rotation, xglobal, yglobal, &
         & zglobal

   !.....Allocate space for the processed bathymetry array h44.....
   ALLOCATE ( h44(1:im1, 1:jm1, 4), stat=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 27 )

   ! ... Only allowed if layers are of uniform thickness
   IF ( ibathyf < 0 ) THEN
     PRINT *, 'STOP - Output to xml files only allowed for ibathyf > 0'
     PRINT *, '*******  Please revise your input files   *************'
     STOP
   ENDIF

   !.....Choose output file as either .xml or .txt.....
   output_file = "si3d_bathy.xml"

   !.....Define global (UTM) coordinates of the corner
   !     grid cell and the rotation angle of grid.....
   ! Hardwire these variables (for the time being)
   xglobal  = 0.0  ! Corner grid cell location is taken as the center of
   yglobal  = 0.0  !    fictitious cell (1,1)
   zglobal  = 0.0  ! zglobal is just a vertical datum adjustment
   rotation = 0.0  ! measured in degrees counterclockwise from North

   !.....Convert variables for .xml header records to characters.....
   ch_xglobal  = real_to_char(xglobal,1,0)
   ch_yglobal  = real_to_char(yglobal,1,0)
   ch_zglobal  = real_to_char(zglobal,1,0)
   ch_rotation = real_to_char(rotation,1,0)
!   ch_idx   = int_to_char(idx,0)
!  ch_idy   = int_to_char(idy,0)
   ch_idx   = real_to_char(idx,1,0) ! idt real
   ch_idy   = real_to_char(idy,1,0) ! idt real
   ch_idz   = real_to_char(idz,1,0)
   ch_nopen = int_to_char(nopen,0)
   !ch_nbarr= int_to_char(nbarr,0)

   !.....Adjust all four corners of each cell into the bottom layer.....
   DO j = 1, jm1
      DO i = 1, im1
         ! Ignore dry cells
         IF ( .NOT. mask2d(i,j) ) CYCLE
         kb = kmz(ij2l(i,j))           ! Store the bottom layer number as kb
         ztop = (kb-k1)*ddz      ! depth to top of bottom layer
         hmin = ztop + dzmin     ! minimum depth allowable
         hmax = ztop + ddz       ! maximum depth allowable
         ! Initialize four corner depths of cell in h44 array
         h44(i,j,1) = h4(i-1,j-1); h44(i,j,2) = h4(i,j-1)
         h44(i,j,3) = h4(i  ,j  ); h44(i,j,4) = h4(i-1,j)
         ! Make sure the corner depths are not
         ! less than hmin or greater then hmax
         DO k = 1, 4
            IF (h44(i,j,k) < hmin) h44(i,j,k) = hmin
            IF (h44(i,j,k) > hmax) h44(i,j,k) = hmax
         END DO
      END DO
   END DO

   !.....Adjust corner depths so that vertical faces in
   !     the bottom profile only begin at layer boundaries.....
   niterations = 2
   DO niters = 1, niterations
      ! Loop over all cells in the grid
      DO j = 1, jm1
         DO i = 1, im1
            ! Ignore dry cells
            IF ( .NOT. mask2d(i,j) ) CYCLE

            kb = kmz(ij2l(i,j))
            ztop = (kb-k1)*ddz
            hmin = ztop + dzmin
            hmax = ztop + ddz
            ! Work on the east side of cell
            IF ( i /= im1 ) THEN
               ! Check if the adjacent cell to the east has
               ! more layers (is deeper) than the present cell
               IF ( kmz(ij2l(i+1,j)) > kb ) THEN
                  ! Change depths in present cell to hmax
                  h44(i,j,2:3) = hmax
               END IF
               ! Check if the adjacent cell to the east has fewer
               ! layers (is shallower) than the present cell
               IF ( kmz(ij2l(i+1,j)) < kb ) THEN
                  hmax_e = (kmz(ij2l(i+1,j))-k1+1)*ddz
                  ! Change depths in eastern cell to hmax_e
                  h44(i+1,j,1) = hmax_e;  h44(i+1,j,4) = hmax_e
               END IF
               ! Check if the adjacent cell to the east has the
               ! same number of layers as the present cell
               IF ( kmz(ij2l(i+1,j)) == kb ) THEN
                  ! Equate depths along common cell boundary
                  h44(i+1,j,1) = h44(i,j,2);  h44(i+1,j,4) = h44(i,j,3)
               END IF
            END IF
            ! Work on the north side of cell
            IF ( j /= jm1 ) THEN
               ! Check if the adjacent cell to the north has
               ! more layers (is deeper) than the present cell
               IF ( kmz(ij2l(i,j+1)) > kb ) THEN
                  ! Change depths in present cell to hmax
                  h44(i,j,3:4) = hmax
               END IF
               ! Check if the adjacent cell to the north has fewer
               ! layers (is shallower) than the present cell
               IF ( kmz(ij2l(i,j+1)) < kb ) THEN
                  hmax_n = (kmz(ij2l(i,j+1))-k1+1)*ddz
                  ! Change depths in northern cell to hmax_n
                  h44(i,j+1,1:2) = hmax_n
               END IF
               ! Check if the adjacent cell to the north has the
               ! same number of layers as the present cell
               IF ( kmz(ij2l(i,j+1)) == kb ) THEN
                  ! Equate depths along common cell boundary
                  h44(i,j+1,1) = h44(i,j,4);  h44(i,j+1,2) = h44(i,j,3)
               END IF
            END IF
         END DO
      END DO
   END DO

   !.....Open output bathymetry file.....
   i90 = i6 + 30
   OPEN (UNIT=i90, FILE=output_file, IOSTAT=ios)
   IF(ios /= 0) CALL open_error ( "Error opening "//output_file, ios )

   !.....Get date and time of run for insertion into header.....
   CALL date_and_time ( date, time, zone, values )

   !.....Output header records to .xml bathymetry file.....
   WRITE (UNIT=i90, FMT='("<?xml version=""1.0""?>")' )
   WRITE (UNIT=i90, FMT='(1X)' )
   WRITE (UNIT=i90, FMT='("<!-- Si3D Processed Bathymetry for Gr -->")')
   WRITE (UNIT=i90, FMT='(1X)' )
   WRITE (UNIT=i90, FMT='("<gov.usgs.gr>")' )
   WRITE (UNIT=i90, FMT='("  <obj class=""gov.usgs.sfhydro.si3d.Bathymetry""")')
   WRITE (UNIT=i90, FMT='("       runTitle=""",      A, """" )' ) TRIM(title)
   WRITE (UNIT=i90, FMT='("       runNumber=""", A8,A4, """" )' ) date,time(1:4)
   WRITE (UNIT=i90, FMT='("       xSpacing=",A," ySpacing=",A," zSpacing=",A, &
      &  " rotation=", A, ">")' ) TRIM(ch_idx),                               &
      &                           TRIM(ch_idy),                               &
      &                           TRIM(ch_idz),                               &
      &                           TRIM(ch_rotation)
   WRITE (UNIT=i90, FMT='(1X)' )
   WRITE (UNIT=i90, FMT='("<!-- Define global coordinates of grid corner -->")')
   WRITE (UNIT=i90, FMT='("  <obj class=""gov.usgs.gr.data.Location"""  &
      &                           " x=",A," y=",A," z=",A,"/>" )' )     &
      &                           TRIM(ch_xglobal),                     &
      &                           TRIM(ch_yglobal),                     &
      &                           TRIM(ch_zglobal)
   WRITE (UNIT=i90, FMT='(1X)' )

   ! Output open boundary information into the file if nopen>0
   IF (nopen > 0) THEN
      WRITE (UNIT=i90, FMT='("<!-- Define node numbers (i,j) for start and end&
         & points of open boundaries -->")' )
      WRITE (UNIT=i90, FMT='("  <obj class=""gov.usgs.gr.GrObjectContainer""")')
      WRITE (UNIT=i90, FMT='("       title=""Open Boundaries"" nopen=",A,">")')&
         &                           TRIM(ch_nopen)
      DO nn = 1, nopen
         IF ( itype(nn) == 1) ch_itype = '"wse"'
         IF ( itype(nn) == 2) ch_itype = '"flw"'
         WRITE (UNIT=i90, FMT='("     <dataseries title=""Boundary", I2, """", &
            & " numDimensions=""2"" itype=", A, ">")' ) nn, ch_itype
         ! Convert node numbers from integers to character variables
         ch_isbc = int_to_char(isbc(nn),1); ch_iebc = int_to_char(iebc(nn),1)
         ch_jsbc = int_to_char(jsbc(nn),1); ch_jebc = int_to_char(jebc(nn),1)
         ! Output start and end node numbers of each open boundary
         WRITE (UNIT=i90, FMT='(A,T6,A)' ) TRIM(ch_isbc), TRIM(ch_jsbc)
         WRITE (UNIT=i90, FMT='(A,T6,A)' ) TRIM(ch_iebc), TRIM(ch_jebc)
         WRITE (UNIT=i90, FMT='("     </dataseries>")' )
      END DO
      WRITE (UNIT=i90, FMT='("  </obj>")' )   ! end open boundary object
      WRITE (UNIT=i90, FMT='(1X)' )
   END IF

   ! Output header information for bathymetry dataseries
   WRITE (UNIT=i90, FMT='("<!-- Bathymetry data series -->")' )
   WRITE (UNIT=i90, FMT='("  <dataseries title=""Bathymetry""",  &
      &                     " numDimensions=""7"">")' )
   WRITE (UNIT=i90, FMT='("     <dim num=""0"" title=""I""/>")' )
   WRITE (UNIT=i90, FMT='("     <dim num=""1"" title=""J""/>")' )
   WRITE (UNIT=i90, FMT='("     <dim num=""2"" title=""Number of Wet Layers""/>")' )
   WRITE (UNIT=i90, FMT='("     <dim num=""3"" title=""SW depth""/>")' )
   WRITE (UNIT=i90, FMT='("     <dim num=""4"" title=""SE depth""/>")' )
   WRITE (UNIT=i90, FMT='("     <dim num=""5"" title=""NE depth""/>")' )
   WRITE (UNIT=i90, FMT='("     <dim num=""6"" title=""NW depth""/>")' )
   WRITE (UNIT=i90, FMT='(1X)' )
   WRITE (UNIT=i90, FMT='("<!--i    j nwlayers hsw     hse     hne     hnw -->")')

   !.....Output depths at the corners of each grid cell in meters.....
   DO j = 1, jm1
      DO i = 1, im1
         ! Ignore dry cells
         IF ( .NOT. mask2d(i,j) ) CYCLE
         WRITE (UNIT=i90, FMT='(3I5,4F8.1)' ) i,j,(kmz(ij2l(i,j))-1),(h44(i,j,k), k=1,4)
      END DO
   END DO

   !.....Output closing xml element end tags.....
   WRITE (UNIT=i90, FMT='("   </dataseries>")' )
   WRITE (UNIT=i90, FMT='(1X)' )
   WRITE (UNIT=i90, FMT='("  </obj>")' )   ! End bathymetry object
   WRITE (UNIT=i90, FMT='("</gov.usgs.gr>")' )

   !.....Deallocate h44 pointer array.....
   DEALLOCATE ( h44 )

END SUBROUTINE outg

!************************************************************************
FUNCTION int_to_char ( ivalue, noquote )
!************************************************************************
!
!  Purpose: To convert an integer value into a 10-character string. Quotes
!           around the number are added if noquote=0.  The returned string
!           has no leading blanks. Typically the string will have trailing
!           blanks.
!
!  More details:
!     The quotemarks are placed around the value for use in xml formatted
!     output files. The character string that is returned from this routine
!     can easily be processed with the TRIM function to remove any
!     trailing blanks.  The length of the returned sting (10 characters)
!     was chosen somewhat arbitrarily.
!
!  Dummy argument:
!  ivalue = An integer number less than 10**8.
!  noquote = Parameter indicating whether quotes are placed around the
!            number within the character string or not. ( 1--> no quotes,
!            0--> quotes )
!
!  Result variable:
!  int_to_char = A 10-character string variable with the character
!                representation of the integer number. The non-blank
!                characters in the string are left-justified.
!
!------------------------------------------------------------------------

   ! .....Argument.....
   INTEGER,INTENT(IN) :: ivalue, noquote
   CHARACTER(LEN=10) :: int_to_char

   !.....Local variables.....
   CHARACTER :: q = '"'
   CHARACTER :: blank = " "
   CHARACTER(LEN=9) :: string = "         "

   !.....Write the integer number to an internal file
   !     and add a trailing quotemark if requested.....
   IF (noquote == 1) THEN
      WRITE ( string, FMT='(I8,A)' ) ivalue, blank
   ELSE
      WRITE ( string, FMT='(I8,A)' ) ivalue, q
   END IF

   !.....Add a leading quotemark if requested
   !     and left-justify the result variable.....
   IF (noquote == 1) THEN
      int_to_char = ADJUSTL(string)//blank
   ELSE
      int_to_char = q//ADJUSTL(string)
   END IF

END FUNCTION int_to_char


!************************************************************************
FUNCTION real_to_char ( value, nd, noquote )
!************************************************************************
!
!  Purpose: To convert a real value into a 14-character string. Quotes
!           around the number are added if noquote=0. The value is placed
!           in as readable a format as possible considering its range.  An
!           exponential format is used for very large or very small numbers.
!           Otherwise, the routine will use "nd" decimal places of accuracy.
!           The returned string has no leading blanks. Typically the string
!           will have trailing blanks.  The routine prints out the number
!           according to the following rules:
!              1. value > 9999999.                 ES12.5
!              2. value < -999999.                 ES12.5
!              3. 0.    < ABS(value) < 0.01        ES12.5
!              4. value = 0.0                      F12.nd
!              5. Otherwise                        F12.nd

!
!  More details:
!     The quotemarks are placed around the value for use in xml formatted
!     output files. The character string that is returned from this routine
!     can easily be processed with the TRIM function to remove any
!     trailing blanks.  The length of the returned sting (14 characters)
!     was chosen somewhat arbitrarily.
!
!  Dummy argument:
!  value = A real number.
!  nd = An integer number specifying the number of decimal places to use in
!       the character representation of the real number. (Must be in the
!       range: 0<nd<10)
!  noquote = Parameter indicating whether quotes are placed around the
!            number within the character string or not. ( 1--> no quotes,
!            0--> quotes )
!
!  Result variable:
!  real_to_char = A 14-character string variable with the character
!                 representation of the real number enclosed in
!                 quotemarks. The non-blank characters in the string
!                 are left-justified.
!
!------------------------------------------------------------------------

   ! .....Argument.....
   INTEGER,INTENT(IN) :: nd, noquote
   REAL,INTENT(IN) :: value
   CHARACTER(LEN=14) :: real_to_char

   !.....Local variables.....
   CHARACTER :: q = '"'
   CHARACTER :: blank = " "
   CHARACTER(LEN=13) :: string = "             "
   CHARACTER(LEN=13) :: fmt    = "             "

   !.....Select proper format and include a trailing quotemark.....
   IF ( noquote == 0 ) THEN
      IF ( value > 9999999. ) THEN
         fmt = '(ES12.5,"""")'
      ELSE IF ( value < -999999. ) THEN
         fmt = '(ES12.5,"""")'
      ELSE IF ( value == 0.) THEN
         WRITE ( fmt, FMT='( "(" , "F12." , I1 , "," , """""""""", ") " )' ) nd
      ELSE IF ( ABS(value) < 0.01 ) THEN
         fmt = '(ES12.5,"""")'
      ELSE
         WRITE ( fmt, FMT='( "(" , "F12." , I1 , "," , """""""""", ") " )' ) nd
      END IF
   END IF

   !.....Select proper format and do not include a trailing quotemark.....
   IF ( noquote /= 0 ) THEN
      IF ( value > 9999999. ) THEN
         fmt = '(ES12.5," ")'
      ELSE IF ( value < -999999. ) THEN
         fmt = '(ES12.5," ")'
      ELSE IF ( value == 0.) THEN
         WRITE ( fmt, FMT='( "(" , "F12." , I1 , "," , """ """, ") " )' ) nd
      ELSE IF ( ABS(value) < 0.01 ) THEN
         fmt = '(ES12.5," ")'
      ELSE
         WRITE ( fmt, FMT='( "(" , "F12." , I1 , "," , """ """, ") " )' ) nd
      END IF
   END IF

   !.....Write the real number to an internal file.....
   WRITE ( string, fmt ) value

   !.....Add a leading quotemark if requested
   !     and left-justify the result variable.....
   IF (noquote == 1) THEN
      real_to_char = ADJUSTL(string)//blank
   ELSE
      real_to_char = q//ADJUSTL(string)
   END IF

END FUNCTION real_to_char


!************************************************************************
FUNCTION double_to_char ( value_dpv, nd, noquote )
!************************************************************************
!
!  Purpose: To convert a double precision value into a 20-character string.
!           Quotes around the number are added if noquote=0. The value is
!           placed in as readable a format as possible considering its range.
!           An exponential format is used for very large or very small numbers.
!           Otherwise, the routine will use "nd" decimal places of accuracy.

!           The returned string has no leading blanks. Typically the string
!           will have trailing blanks.  The routine prints out the number
!           according to the following rules:
!              1. value_dpv > 9999999.                 ES18.11
!              2. value_dpv < -999999.                 ES18.11
!              3. 0.    < ABS(value_dpv) < 0.01        ES18.11
!              4. value_dpv = 0.0                      F18.nd
!              5. Otherwise                            F18.nd

!
!  More details:
!     The quotemarks are placed around the value for use in xml formatted
!     output files. The character string that is returned from this routine
!     can easily be processed with the TRIM function to remove any
!     trailing blanks.  The length of the returned sting (20 characters)
!     was chosen somewhat arbitrarily.
!
!  Dummy argument:
!  value_dpv = A real number of KIND(1.0D0).
!  nd = An integer number specifying the number of decimal places to use in
!       the character representation of the real number. (Must be in the
!       range: 0<nd<10)
!  noquote = Parameter indicating whether quotes are placed around the
!            number within the character string or not. ( 1--> no quotes,
!            0--> quotes )
!
!  Result variable:
!  double_to_char = A 20-character string variable with the character
!                   representation of the double precision number enclosed
!                   in quotemarks. The non-blank characters in the string
!                   are left-justified.
!
!
!------------------------------------------------------------------------

   ! .....Argument.....
   INTEGER,INTENT(IN) :: nd, noquote
   REAL(DPV),INTENT(IN) :: value_dpv
   CHARACTER(LEN=20) :: double_to_char

   !.....Local variables.....
   CHARACTER :: q = '"'
   CHARACTER :: blank = " "
   CHARACTER(LEN=19) :: string = "                   "
   CHARACTER(LEN=14) :: fmt    = "              "

   !.....Select proper format and include a trailing quotemark.....
   IF ( noquote == 0 ) THEN
      IF ( value_dpv > 9999999. ) THEN
         fmt = '(ES18.11,"""")'
      ELSE IF ( value_dpv < -999999. ) THEN
         fmt = '(ES18.11,"""")'
      ELSE IF ( value_dpv == 0._dpv) THEN
         WRITE ( fmt, FMT='( "(" , "F18." , I1 , "," , """""""""", ")  " )' ) nd
      ELSE IF ( ABS(value_dpv) < 0.01 ) THEN
         fmt = '(ES18.11,"""")'
      ELSE
         WRITE ( fmt, FMT='( "(" , "F18." , I1 , "," , """""""""", ")  " )' ) nd
      END IF
   END IF

   !.....Select proper format and do not include a trailing quotemark.....
   IF ( noquote /= 0 ) THEN
      IF ( value_dpv > 9999999. ) THEN
         fmt = '(ES18.11," ")'
      ELSE IF ( value_dpv < -999999. ) THEN
         fmt = '(ES18.11," ")'
      ELSE IF ( value_dpv == 0.) THEN
         WRITE ( fmt, FMT='( "(" , "F18." , I1 , "," , """ """, ") "  )' ) nd
      ELSE IF ( ABS(value_dpv) < 0.01 ) THEN
         fmt = '(ES18.11," ")'
      ELSE
         WRITE ( fmt, FMT='( "(" , "F18." , I1 , "," , """ """, ") "  )' ) nd
      END IF
   END IF

   !.....Write the double precision number to an internal file.....
   WRITE ( string, fmt ) value_dpv

   !.....Add a leading quotemark if requested
   !     and left-justify the result variable.....
   IF (noquote == 1) THEN
      double_to_char = ADJUSTL(string)//blank
   ELSE
      double_to_char = q//ADJUSTL(string)
   END IF

END FUNCTION double_to_char

!***********************************************************************
SUBROUTINE nodech ( i, j, nodeno, nchar )
!***********************************************************************
!
!  Purpose: To convert the node number where model time series output is
!           desired to a character variable.
!
!  Dummy arguments:
!  i,j    = x and y node mumbers expressed as integers (maximum allowable
!           value for i or j is 9999)
!  nodeno = node number expressed as a character variable
!  nchar  = number of characters in nodeno (excluding trailing blanks)
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!  6/18/01           P.E. Smith        Original f90 code
!
!-----------------------------------------------------------------------

   !.....Arguments.....
   INTEGER, INTENT(IN)  :: i, j
   INTEGER, INTENT(OUT) :: nchar
   CHARACTER(LEN=9), INTENT(OUT) :: nodeno

   !.....Convert node numbers to character variables.....
   IF ( i < 10) THEN
      WRITE ( nodeno(1:2), FMT='(I1,"_")' ) i
      IF ( j < 10 ) THEN
         WRITE ( nodeno(3:9), FMT='(I1,"      ")' ) j
         nchar = 3
      END IF
      IF(( j >= 10 ) .AND. ( j < 100 )) THEN
         WRITE ( nodeno(3:9), FMT='(I2,"     " )' ) j
         nchar = 4
      END IF
      IF ((j >= 100) .AND. ( j < 1000)) THEN
         WRITE ( nodeno(3:9), FMT='(I3,"    "  )' ) j
         nchar = 5
      END IF
      IF ( j >= 1000) THEN
         WRITE ( nodeno(3:9), FMT='(I4,"   "   )' ) j
         nchar = 6
      END IF
   END IF
   IF (( i >= 10 ) .AND. ( i < 100 )) THEN
      WRITE ( nodeno(1:3), FMT='(I2,"_")' ) i
      IF ( j < 10 ) THEN
         WRITE ( nodeno(4:9), FMT='(I1,"     ")' ) j
         nchar = 4
      END IF
      IF(( j >= 10 ) .AND. ( j < 100 )) THEN
         WRITE ( nodeno(4:9), FMT='(I2,"    " )' ) j
         nchar = 5
      END IF
      IF ((j >= 100) .AND. ( j < 1000)) THEN
         WRITE ( nodeno(4:9), FMT='(I3,"   "  )' ) j
         nchar = 6
      END IF
      IF ( j >= 1000) THEN
         WRITE ( nodeno(4:9), FMT='(I4,"  "   )' ) j
         nchar = 7
      END IF
   END IF
   IF (( i >= 100 ) .AND. ( i < 1000)) THEN
      WRITE ( nodeno(1:4), FMT='(I3,"_")' ) i
      IF ( j < 10 ) THEN
         WRITE ( nodeno(5:9), FMT='(I1,"    ")' ) j
         nchar = 5
      END IF
      IF(( j >= 10 ) .AND. ( j < 100 )) THEN
         WRITE ( nodeno(5:9), FMT='(I2,"   " )' ) j
         nchar = 6
      END IF
      IF (( j >= 100) .AND. (j < 1000)) THEN
         WRITE ( nodeno(5:9), FMT='(I3,"  "  )' ) j
         nchar = 7
      END IF
      IF ( j >= 1000) THEN
         WRITE ( nodeno(5:9), FMT='(I4," "   )' ) j
         nchar = 8
      END IF
   END IF
   IF ( i >= 1000 ) THEN
      WRITE ( nodeno(1:5), FMT='(I4,"_")' ) i
      IF ( j < 10 ) THEN
         WRITE ( nodeno(6:9), FMT='(I1,"   ")' ) j
         nchar = 6
      END IF
      IF(( j >= 10 ) .AND. ( j < 100 )) THEN
         WRITE ( nodeno(6:9), FMT='(I2,"  " )' ) j
         nchar = 7
      END IF
      IF (( j >= 100) .AND. (j < 1000)) THEN
         WRITE ( nodeno(6:9), FMT='(I3," "  )' ) j
         nchar = 8
      END IF
      IF ( j >= 1000) THEN
         WRITE ( nodeno(6:9), FMT='(I4      )' ) j
         nchar = 9
      END IF
   END IF

END SUBROUTINE nodech


!***********************************************************************
SUBROUTINE outlog(thrs)
!***********************************************************************
!
!  Purpose: To write the current simulation time and time step number
!           to a log file that can be used to externally monitor the
!           progress of the run.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-----------------------------------------------------------------------

   REAL, INTENT(IN) :: thrs

   !.....Local variables.....
   INTEGER :: ios

   !.....Open the log file.....
   OPEN ( UNIT=i82, FILE="si3d_log.txt", IOSTAT = ios )
   IF ( ios /= 0 ) CALL open_error ( "Error opening si3d_log.txt", ios )

   !.....Write time in hours and time step number to file.....
   WRITE (UNIT=i82, FMT='(" time=", F10.4, " hours", "    n=", I6)' ) thrs, n

   !.....Close the log file during each subroutine call to cause
   !     any output held in a buffer to be written to the file.....
   IF ( n < nts ) THEN
      CLOSE (UNIT=i82, STATUS="KEEP")
   ELSE
      CLOSE (UNIT=i82, STATUS="DELETE")
   END IF

END SUBROUTINE outlog

!***********************************************************************
SUBROUTINE check_stopfile(n)
!***********************************************************************
!
!  Purpose: To check whether the execution of the program is to be
!           stopped, and then, if requested, to stop the program.  The
!           routine reads a variable 'istop' that can be set manually
!           during the execution of the program.  The value of 'istop'
!           is initialized to zero.  If during the execution of the
!           program the value of 'istop' is manually reset to 1, the
!           program will close all files and immediately terminate
!           execution of the simulation.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!-----------------------------------------------------------------------

   INTEGER,INTENT(IN) :: n

   !.....Local variables.....
   INTEGER :: ios, istop
   INTEGER, SAVE :: nfirst = 0

   !.....Open stop file.....
   OPEN (UNIT=i83,FILE="si3d_stop.txt",ACTION="READWRITE",IOSTAT=ios)
   IF ( ios /= 0 ) CALL open_error ( "Error opening si3d_stop.txt", ios )

   !.....Write istop=0 on the first entry into the subroutine.....
   IF ( nfirst == 0 ) THEN
      nfirst = 1
      WRITE (UNIT=i83, FMT='(" istop = 0")')
      ENDFILE i83
      CLOSE (UNIT=i83)
      RETURN
   END IF

   !.....Read the value of istop.....
   READ (UNIT=i83, FMT='(8X, I2)' ) istop

   !.....Stop execution of program if istop=1.....
   IF ( istop == 1 ) THEN
      CLOSE (UNIT=i83, STATUS="DELETE")  ! Close and delete stop file
      PRINT *, " ****STOPPING si3d after istop was manually set to 1"
      WRITE (UNIT=i6, FMT='(" ")')
      WRITE (UNIT=i6, FMT='(" ****STOPPING si3d after istop was manually &
         &set to 1")')
      WRITE (UNIT=i6, FMT='(" n = ", I7)') n
      STOP
   END IF

   !.....Close the stop file.....
   IF ( n <  nts ) THEN
      CLOSE (UNIT=i83)   ! Close stop file every time step for PC
   ELSE
      CLOSE (UNIT=i83, STATUS="DELETE")   ! Delete file when n=nts
   END IF

END SUBROUTINE check_stopfile

!***********************************************************************
SUBROUTINE compute_date ( itime_sec )
!***********************************************************************
!
!  Purpose: To compute the date and time given the time in seconds from
!           a previous date and time.  The previous date and time are
!           saved from the last call to the subroutine.
!
!  Variables:
!     iyr0,imon0,iday0,ihr0,isec0 -- Start date and time of the run defined
!                                    in the input data file. (Used only
!                                    in the first call to the subroutine.)
!     iyrp,imonp,idayp,ihrp,isecp -- Previous date and time.
!                                    (Saved during the previous call to
!                                    the subroutine.)
!     iyr ,imon ,iday ,ihr ,isec  -- New date and time.
!                                    (Computed in the subroutine and
!                                    returned to the calling program.)
!     doy, doyp                   -- New & previous times (fractional julian day)
!
!  More details:
!     iyr  = year  (4 digit integer)
!     imon = month (2 digit integer)
!     iday = day   (2 digit integer)
!     hrs  = time in decimal hours from beginning of iday (real)
!     ihr  = time in decimal hours multiplied by 100 and rounded to the
!            nearest integer (4 digit integer) (ihr=NINT(hrs*100.))
!     isec = time in seconds from beginning of iday (5 digit integer)
!
!  Dummy argument:
!  itime_sec = time in seconds from the previous date and time.
!
!-----------------------------------------------------------------------

   !.....Argument.....
   REAL, INTENT(IN) :: itime_sec

   !.....Local variables.....
   INTEGER, DIMENSION(12) :: month=(/31,28,31,30,31,30,31,31,30,31,30,31/)
   INTEGER, SAVE :: iyrp, imonp, idayp, ihrp  ! idt real
   REAL   , SAVE :: isecp                     ! idt real
   INTEGER :: max_possible_months, jmonths
   REAL :: hrs

   !.....Initialize date on first entry into subroutine
   !     with the starting date of the simulation.....
   IF (itime_sec == 0) THEN

      ! . Year, month, day, hour, and seconds
      iyrp=iyr0; imonp=imon0; idayp=iday0; ihrp=ihr0; isecp=isec0
      iyr =iyr0; imon =imon0; iday =iday0; ihr =ihr0; isec =isec0

      ! . Julian day
      IF (leap_year(iyrp)) month(2) = 29
      IF (imon0 == 1) THEN
        doyp = iday0+FLOAT(ihr)/2400.
        doy  = iday0+FLOAT(ihr)/2400.
      ELSE
        doyp = SUM(month(1:imon0-1))+iday0+FLOAT(ihr)/2400.
        doy  = SUM(month(1:imon0-1))+iday0+FLOAT(ihr)/2400.
      ENDIF
      RETURN

   END IF

   ! ... Store previous values for day of year doy
   doyp = doy;

   !.....Compute new date.....
   isecp = isecp + itime_sec
   IF (isecp < 86400) THEN
      isec = isecp;
      hrs  = isec/3600.;
      ihrp = NINT(hrs*100.);
      ihr  = ihrp
      !doy  = doyp + FLOAT(itime_sec)/86400. ! idt real
      doy  = doyp +       (itime_sec)/86400. ! idt real
      RETURN
   END IF
   idayp = idayp + isecp/86400
   !isecp = MOD(isecp,86400); isec = isecp  ! idt real
   isecp = MOD(isecp,86400.); isec = isecp  ! idt real
   hrs   = isec/3600.
   ihrp  = NINT(hrs*100.);
   ihr = ihrp
   ! Take care of leap year
   IF (imonp == 2) THEN
      month(2) = 28
      IF (leap_year(iyrp)) month(2) = 29
   END IF
   IF (idayp <= month(imonp)) THEN
      iday = idayp
      IF (imon0 == 1) THEN
        !doy = iday+FLOAT(isec)/86400. ! idt real
        doy = iday+      (isec)/86400. ! idt real
      ELSE
        !doy = SUM(month(1:imon0-1))+iday0+FLOAT(isec)/86400. ! idt real
        doy = SUM(month(1:imon0-1))+iday0+      (isec)/86400. ! idt real
      ENDIF
      RETURN
   END IF
   ! Take care of case where 'itime_sec' is very large (> 1 month)
   max_possible_months = (itime_sec/2500000) + 1
   DO jmonths = 1, max_possible_months
      idayp = idayp - month(imonp)
      imonp = imonp + 1
      IF (imonp == 13) THEN
         imonp = 1; iyrp = iyrp + 1
         IF (leap_year(iyrp)) month(2) = 29
         IF (.NOT. leap_year(iyrp)) month(2) = 28
      END IF
      IF (idayp <= month(imonp)) THEN
         iday = idayp
         EXIT
      END IF
   END DO
   IF (imon0 == 1) THEN
     !doy  = iday+FLOAT(isec)/86400. ! idt real
     doy  = iday+      (isec)/86400. ! idt real
   ELSE
     !doy  = SUM(month(1:imon0-1))+iday0+FLOAT(isec)/86400.  ! idt real
     doy  = SUM(month(1:imon0-1))+iday0+      (isec)/86400.  ! idt real
   ENDIF
   imon = imonp
   iyr  = iyrp

END SUBROUTINE compute_date

!***********************************************************************
LOGICAL FUNCTION leap_year ( iyear )
!***********************************************************************
!
!  Purpose: To identify whether 'iyear' is a leap year or not. The
!           function returns .TRUE. if 'iyear' is a leap year.
!
!  (Note: Formally, a leap year must be divisible by 4 but not by 100,
!         or it must be divisible by 400.  For this subroutine, I only
!         check if a year is divisible by 4.)
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-----------------------------------------------------------------------

   !.....Argument.....
   INTEGER, INTENT(IN) :: iyear

   leap_year  = MOD(iyear,4) == 0

END FUNCTION leap_year

!************************************************************************
  SUBROUTINE OutScalarEnergyBalance(n,thrs)
!************************************************************************
!
!  Purpose: To output mass & energy existing at the basin
!           scale.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!---------------------------------------------------------------

   INTEGER, INTENT(IN) :: n
   REAL, INTENT(IN) :: thrs

   ! ... Local Variables
   CHARACTER(LEN=25):: flux_file="ScalarBalance.txt"
   INTEGER, SAVE	:: i899 = 899
   INTEGER              :: ios, i,j,k,l,k1s, kms
   REAL                 :: elev, uijk, vijk, wijk, rijk
   REAL(real_G1)	:: HeatB, OxygB, PotE, KinE
   REAL, PARAMETER	:: SpecificHeat = 4181.6

   ! ... Open output file on first call
   IF ( n == 0 ) THEN
      OPEN ( UNIT=i899, FILE=flux_file, IOSTAT=ios)
      IF(ios /= 0) PRINT *, "Error opening "//flux_file
      ShearCum = 0.0
      BuocyCum = 0.0
      DisspCum = 0.0
   END IF

   ! ... Heat & Oxygen in the domain
   HeatB = 0.0E0
   OxygB = 0.0E0
   DO l = 1, lm;
     i = l2i(l); j = l2j(l);
     kms = kmz(l)
     k1s = k1z(l)
     DO k = k1s, kms
       HeatB = HeatB+sal   (k,l    )*h(k,l)
       IF (ntr > 0) &
       OxygB = OxygB+tracer(k,l,ntr)*h(k,l)
      END DO
   END DO;

   ! ... Basin scale potential energy
   PotE = 0.0
   DO l = 1, lm;
     i = l2i(l); j = l2j(l);
     kms = kmz(l)
     k1s = k1z(l)
     elev = hhs(l)-h(kms,l)/2.
     rijk = densty_s(salp(kms,l),0.0)+1000.
     PotE = PotE + rijk*g*elev*h(kms,l)
     IF (k1s == kms) CYCLE
     DO k = kms-1, k1s, -1
       elev = elev + (hp(k+1,l)+hp(k,l))/2.
       rijk = densty_s(salp(k,l),0.0)+1000.
       PotE = PotE + rijk*g*elev*hp(k,l)
     END DO
   END DO;
   PotE = PotE*dx*dy

   ! ... Basin scale kinetic energy
   KinE = 0.0
   DO l = 1, lm;
     i = l2i(l); j = l2j(l);
     kms = kmz(l)
     k1s = k1z(l)
     ! ... Basin scale potential energy
     uijk = (up(kms,l) + up(kms  ,lWC(l)))/2.
     vijk = (vp(kms,l) + vp(kms  ,lSC(l)))/2.
     wijk = (wp(kms,l) + wp(kms+1,    l ))/2.
     rijk = densty_s(salp(kms,l),0.0)+1000.
     KinE = KinE + 0.5*rijk*(uijk**2.+vijk**2.+wijk**2.)*h(kms,l)
     IF (k1s == kms) CYCLE
     DO k = kms-1, k1s, -1
       uijk = (up(k,l) + up(k  ,lWC(l)))/2.
       vijk = (vp(k,l) + vp(k  ,lSC(l)))/2.
       wijk = (wp(k,l) + wp(k+1,    l ))/2.
       rijk = densty_s(salp(k,l),0.0)+1000.
       KinE = KinE + 0.5*rijk*(uijk**2.+vijk**2.+wijk**2.)*h(k,l)
     END DO
   END DO;
   KinE = KinE*dx*dy;

   WRITE (UNIT=i899, FMT='(9E20.11)') thrs, HeatB, OxygB, &
                                      PotE, KinE, TKinE,  &
                                      ShearCum, BuocyCum, DisspCum

   ! ... Set to zero integral variables
   ShearCum = 0.0
   BuocyCum = 0.0
   DisspCum = 0.0

  END SUBROUTINE OutScalarEnergyBalance




!***********************************************************************
SUBROUTINE cputimes(t_exmom2,t_matmom2,t_matcon2)
!***********************************************************************
!
!  Purpose: To output CPU times for various subroutines in the model.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-----------------------------------------------------------------------

   REAL, INTENT(IN) :: t_exmom2,t_matmom2,t_matcon2

   !.....Calculate total time in subroutines.....
   tot_subs=t_exmom2+t_matmom2+t_matcon2+t_solver+t_vel+t_exsal+t_salin &
      &    +t_outt+t_turb+t_settrap+t_save

   !.....Print CPU times.....
   PRINT '(A, I10)', " Number of calls to turb  =", n_turb
   PRINT '(A, I10)', " Number of calls to exmom =", n_exmom
   PRINT '(A)' , " "
   PRINT '(A)', "        --CPU Times--"
   PRINT '(A, F10.3, A)', " t_turb   =", t_turb ,  " seconds"
   print *,"hilo:",omp_get_thread_num(),"t_exmom =",t_exmom2, " seconds"
   print *,"hilo:",omp_get_thread_num(),"t_matmom =",t_matmom2, " seconds"
   print *,"hilo:",omp_get_thread_num(),"t_matcon =",t_matcon2, " seconds"
 !  PRINT '(A, F10.3, A)', " t_exmom  =", t_exmom,  " seconds"
 !  PRINT '(A, F10.3, A)', " t_matmom =", t_matmom, " seconds"
  !PRINT '(A, F10.3, A)', "   t_trid =", t_trid,   " seconds"
 !  PRINT '(A, F10.3, A)', " t_matcon =", t_matcon, " seconds"
   PRINT '(A, F10.3, A)', " t_solver =", t_solver, " seconds"
   PRINT '(A, F10.3, A)', " t_vel    =", t_vel,    " seconds"
   PRINT '(A, F10.3, A)', " t_exsal  =", t_exsal,  " seconds"
   PRINT '(A, F10.3, A)', " t_salin  =", t_salin,  " seconds"
   PRINT '(A, F10.3, A)', " t_settrap=", t_settrap," seconds"
   PRINT '(A, F10.3, A)', " t_outt   =", t_outt ,  " seconds"
   PRINT '(A, F10.3, A)', " t_save   =", t_save ,  " seconds"
   PRINT '(A, F10.3, A)', " t_setmask=", t_setmask," seconds"
   PRINT '(A, F10.3, A)', " tot_subs =", tot_subs, " seconds"
   PRINT '(A)' , " "

END SUBROUTINE cputimes


!***********************************************************************
PURE FUNCTION densty_s ( temperature, salinity )
!***********************************************************************
!
!  Purpose: To compute density (in kg/m**3) from active scalars
!           It uses UNESCO Eq.of state for density of freshwater
!           taken from Gill(1982) - Atmosphere-Ocean Dynamics, Appendix 3
!           However, at this point pressure (depth) effects are not
!           included in the calculation of water density.
!           This function is based on the original function written by
!           P.E. Smith in which the first arg. was salinity and the 2nd temp.
!           Here, we use temp. as first argument, as it is the first arg. whose
!           transport equation is solved in the code. These changes
!           were made as temp. is the main active scalar in Stockton Channel.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-----------------------------------------------------------------------

    ! ... Io variables
	REAL, INTENT(IN) :: temperature, salinity
	REAL             :: densty_s

    densty_s =999.842594                &
      +6.793952e-2*temperature          &
      -9.095290e-3*temperature**2.      &
      +1.001685e-4*temperature**3.      &
      -1.120083e-6*temperature**4.      &
      +6.536332e-9*temperature**5.

END FUNCTION densty_s


!***********************************************************************
SUBROUTINE PointSourceSinkSolve(n,istep,thrs)
!***********************************************************************
!
!  Purpose: Interface between si3d and other models requiring
!           sources and sinks (including VT-plume model).
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!   29-may-09        FJRUeda           Modifies call to rwps to be consisten
!                                      with latest version of VTPlume
!   21-jul-10        FJRueda           Version that includes boundary conditions
!                                      & sources sinks of WATER
!
!-----------------------------------------------------------------------

   INTEGER,INTENT(IN) :: n,istep
   REAL,INTENT(IN) :: thrs

   INTEGER :: nn, inn, i, j, k, kk, l, k1s, kms, nwl, ktop, ksrc, plmdim, itr,innH
   REAL    :: areatot, Tsource, Rsource
   REAL(8) :: wselev, dfelev, dflgth, hcell, rjulday, lambnot, diamm
   REAL(8) :: elevt, qwt, tplumt, comgpt, tplum0, comgp0, qscfm, frconot
   REAL(8), DIMENSION (1:km1) :: zamb, Tamb, DOamb ! B.C. for plume model
   REAL(8), DIMENSION (1:km1) :: qwd               ! Outflow rate for plume
   LOGICAL, SAVE :: DiffON

   ! ... Return if no points sources/sinks are specified
   IF (iopss <= 0) RETURN

   ! ... Define Qpss, Tpss, & Rpss for columns in each device
   DO nn = 1, npssdev

     SELECT CASE (ptype(nn))

     ! ************ Section Boundary Conditions ***********************
     CASE (-2)

       ! ... Only do computations when there is flow
       IF (ABS(flpss(nn)) < qthrs(nn)) CYCLE

       ! ... Only determine boundary conditions on leapfrog iterations
       IF( (n > 1) .AND. (istep > 1)) CYCLE

       ! ... Interpolate forcing variables (flpss, scpss, trpss)
       !     to present time using time series input -
       CALL PointSourceSinkForcing (nn,thrs)

       ! ... Find total volume of cells holding boundary conditions
       areatot = 0.0
       DO innH = 1, iopssH(omp_get_thread_num ( )+1)
            inn = ioph2iop(innH,omp_get_thread_num ( )+1)
         IF (iodev(inn) .NE. nn) CYCLE
         ! ... Define i,j,l indexes
         i = ipss(inn);
         j = jpss(inn);
         l = ij2l(i,j);
         ! ... Define top and bottom cells & no. of layers
         !     The diffuser is set to be one cell above the bottom
         k1s = k1z(l) ;
         kms = kmz(l) ;
         DO k = k1s, kms
           areatot = areatot + hp(k,l)
         ENDDO
       ENDDO

	   ! ... Determine flow rate, temp. and tracer for each cell
	   !     Inflow rate is pressumed uniform in space
       DO innH = 1, iopssH(omp_get_thread_num ( )+1)
            inn = ioph2iop(innH,omp_get_thread_num ( )+1)

         IF (iodev(inn) .NE. nn) CYCLE

         ! ... Define i,j,l indexes
         i = ipss(inn);
         j = jpss(inn);
         l = ij2l(i,j);

         ! ... Define k- indexes
         k1s = k1z(l) ;
         kms = kmz(l) ;

         ! ... Loop over cells in the water column
         DO k = k1s, kms
           Qpss(k,inn) = flpss(nn) * hp(k,l) / areatot
           IF (scpss(nn)<0.0 .OR. flpss(nn)<=0.0) THEN
             Tpss(k,inn) = salp (k,l)
           ELSE
             Tpss(k,inn) = scpss(nn )
           ENDIF
           IF (ntr > 0) THEN
             DO itr = 1, ntr
               IF (trpss(nn,itr)<0.0 .OR. flpss(nn)<=0.0) THEN
                 Rpss(k,inn,itr) = tracerpp(k,l,itr)
               ELSE
                 Rpss(k,inn,itr) = trpss(nn,itr)
               ENDIF
             ENDDO
           ENDIF
         ENDDO
       ENDDO


     ! ************ Bottom cell Boundary Conditions  ********************
     CASE (-1)

       ! ... Only do computations when there is flow
       IF (ABS(flpss(nn)) < qthrs(nn)) CYCLE

       ! ... Only determine boundary conditions on leapfrog iterations
       IF( (n > 1) .AND. (istep > 1)) CYCLE

       ! ... Interpolate forcing variables (flpss, scpss, trpss) for each device
       !     to present time,  using time series input -
       CALL PointSourceSinkForcing (nn,thrs)

       DO innH = 1, iopssH(omp_get_thread_num ( )+1)
            inn = ioph2iop(innH,omp_get_thread_num ( )+1)

         IF (iodev(inn) .NE. nn) CYCLE

         ! ... Define i,j,l indexes
         i = ipss(inn);
         j = jpss(inn);
         l = ij2l(i,j);

         ! ... Define k- indexes
         kms = kmz(l) ;

         ! ... Loop over cells in the water column
		 Qpss(:,inn)   = 0.0;
		 Tpss(:,inn)   = 0.0;
		 Rpss(:,inn,:) = 0.0;
	     Qpss(kms,inn  ) = flpss(nn)
		 !PRINT *, Qpss(kms,inn), flpss(nn)
	     IF (scpss(nn)<0.0 .OR. flpss(nn)<=0.0) THEN
  	       Tpss(kms,inn) = salp (kms,l)
	     ELSE
  	       Tpss(kms,inn) = scpss(nn   )
		 !PRINT *, scpss(nn)
         ENDIF
	     IF (ntr > 0) THEN
	       DO itr = 1, ntr
             IF (trpss(nn,itr)<0.0 .OR. flpss(nn) <= 0.0) THEN
  	           Rpss(kms,inn,itr) = tracerpp(kms,l,itr)
             ELSE
               Rpss(kms,inn,itr) = trpss(nn,itr)
	         ENDIF
           ENDDO
         ENDIF
       ENDDO ! Loop over columns in device


     ! ************ Water pumped inflow ****************************
     CASE (0)

        PRINT *, '***************** ERROR *****************'
        PRINT *, 'Water pumped inflow still NOT incorporated'
        PRINT *, '***************** ERROR *****************'
	    STOP


     ! ************ Oxygen-gas diffuser ****************************
     CASE (1:)


       ! ... Calculate en-detrainment flows induced by diffuser
       IF (ABS(flpss(nn)) < qthrs(nn)) THEN  ! Diffuser OFF
         DiffON = .FALSE.
         DO inn = 1, iopss
           IF (iodev(inn) .NE. nn) CYCLE
           Qpss (inn,:) = 0.0E0
		   kdetr(inn  ) = km1
         ENDDO
       ELSE                                  ! Update diffuser FLOWS
         IF &
         ( (DiffON == .FALSE.)       .OR.  & ! Diffuser is TURNED ON
         (  istep  ==  1            .AND.  & ! Update on first iterations
         (MOD(n,MAX(pdt(nn),1))==0))) THEN   ! Update every pdt time steps

           DiffON = .TRUE.                     ! Diffuser remains ON
           DO innH = 1, iopssH(omp_get_thread_num ( )+1)
            inn = ioph2iop(innH,omp_get_thread_num ( )+1)
             IF (iodev(inn) .NE. nn) CYCLE

             ! ... Define i,j,l indexes
             i = ipss(inn);
             j = jpss(inn);
             l = ij2l(i,j);

             ! ... Define k- indexes
             k1s = k1z(l) ;
             kms = kmz(l) ;
             nwl = kms-k1s+1;

             ! ... Interpolate forcing variables (flpss, scpss, trpss) for each device
             !     to present time,  using time series input -
             CALL PointSourceSinkForcing (nn,thrs)

             ! ... Define ambient temperatures
             Tamb(k1s:kms) = salpp(k1s:kms,l)
             Tamb(kms+1  ) = Tamb(kms);
             Tamb(1      ) = Tamb(k1s);

             ! ... Define DO concentrations
             DOamb(k1s:kms) = tracerpp(k1s:kms,l,ntr)
             DOamb(kms+1  ) = DOamb(kms);
             DOamb(1      ) = DOamb(k1s);

             ! ... Depths for cells in plume column from datum
             zamb(k1s  ) = hp(k1s,l)/2.
             DO k = k1s+1, kms
               zamb(k) = zamb(k-1) + (hp(k-1,l)+hp(k,l))/2.
             END DO
             zamb(kms+1) =  zamb(kms)+hp(kms,l)
             zamb(1    ) = -zamb(k1s)

             ! ... Inputs for plume model
             dfLgth  = dfL(nn)         ;       ! Length of diffuser
             rjulday = doy             ;       ! Julian day (arbitrary)
             wselev  = 0.0000          ;       ! Elevation of free surface
             ksrc    = kms-1           ;       ! Layer No. where diffuser is located
             dfelev  = -zamb(ksrc)     ;       ! Elevation of diffuser
             hcell   = ddz             ;       ! Pressumed constant - thickess of cells
             qwd     = 0.0E0           ;       ! Initialize qwd
             qscfm   = flpss(nn)       ;       ! Air flow rate
             frconot = 0.21            ;       ! Fraction of O2 in air (not used?)
			 lambnot = lambda(nn)      ;       ! Half-width
			 diamm   = diammb(nn)      ;       ! Initial bubble diameter
             IF (ptype(nn) <= 2) THEN
               plmdim = 1 ! Linear Plume
             ELSE
               plmdim = 2 ! Circular Plume
             ENDIF

             ! ... Run plume model
             CALL lineplu_v1(iyr,rjulday,wselev,dfelev,kms,dfLgth, &
                               lambnot,salamb,patm, diamm, plmdim, &
                               qscfm, frconot,                     & ! boundary conditions
                               ksrc, hcell,                        &
                               zamb (1:kms),                       &
                               Tamb (1:kms),                       &
                               DOamb(1:kms),                       &
                               elevt, qwt , tplumt, comgpt,        &
                               qwd(1:ksrc), ktop)

             ! ....Detrainment cell - save it for later
             kdetr(inn) = ktop

             ! ... Define flow at cells
             Qpss(:,inn)   = 0.0;

             ! ... Define flow at entrainment cells
             DO k = ktop+1,ksrc
                Qpss(k,inn) = -qwd(k)*dy/dfLgth
             ENDDO

             ! ... Define flow at detrainment cell to force volume conservation
             Qpss(ktop,inn) = -SUM(Qpss(ktop+1:ksrc,inn));

             PRINT *, '****************************************************'
             PRINT *, '----OUTPUT FROM PLUME ROUTINES----------------------'
             PRINT *, '****************************************************'
             PRINT *, 'Results of Plume Model elev', elevt, ktop
             PRINT *, 'Resutls of Plume Model qwt ', qwt
             PRINT *, 'Results of Plume Model T&O ', tplumt, comgpt
             PRINT *, '****************************************************'
           ENDDO
         ENDIF
       ENDIF

       ! ... Define temperature of entrained and detrained water in plume
       !     based on existing values of Qpss & temperatures
       DO innH = 1, iopssH(omp_get_thread_num ( )+1)
            inn = ioph2iop(innH,omp_get_thread_num ( )+1)
         IF (iodev(inn) .NE. nn) CYCLE

         ! ... Define i,j,l indexes
         i = ipss(inn);
         j = jpss(inn);
         l = ij2l(i,j);

         ! ... Define k- indexes
         k1s = k1z(l) ;
         kms = kmz(l) ;
         nwl = kms-k1s+1;
         Tpss(:,inn) = salpp(:,l)
         k = kdetr(inn);
         IF (k < kms) THEN
           DO kk = k+1,kms
             Tsource  = Tsource + salpp(kk,l)*Qpss(kk,inn)
           ENDDO
           Tsource = Tsource / SUM(Qpss(k+1:kms,inn))
         ELSE
           Tsource = salpp(k,l)
         ENDIF
         Tpss(k,inn) = Tsource

         ! ... Define tracer concentrations in entrained & detrained water
         !     in the plume (assumes dx = dy) based on existing values of
         !     Qpss, tracer concs. and location of detrainment cell
         IF (ntr > 0) THEN
           DO itr = 1, ntr
             Rpss(:,inn,itr) = tracerpp(:,l,itr)
             k = kdetr(inn);
             IF (k < kms) THEN
               Rsource = 0.0
               DO kk = k+1,kms
                 Rsource  = Rsource + tracerpp(kk,l,itr) * Qpss(kk,inn)
               ENDDO
               Rsource = Rsource/SUM(Qpss(k+1:kms,inn))+  &
                         trpss(nn,itr)*dy/dfL(nn)/Qpss(k,inn)
             ELSE
               Rsource = trpss(nn,itr)*dy/dfL(nn)/Qpss(k,inn)
             ENDIF
             Rpss(k,inn,itr) = Rsource
           ENDDO
         ENDIF
       ENDDO

     END SELECT

   ENDDO

END SUBROUTINE PointSourceSinkSolve

!************************************************************************
SUBROUTINE PointSourceSinkInput
!************************************************************************
!
!  Purpose: This routine is called at the beginning of the program
!           to open any files with data for point sources/sinks - it will
!           read the time series data and assign the initial values at time t=0.
!
!------------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: i, j, k, l, nn, itr, ios, istat, nptspss, ncdev
   CHARACTER(LEN=14) :: pssfmt, pssfile

   ! ... Allocate space for device characteristics
   ALLOCATE (  ptype (npssdev), &               ! Type of device simulated
               dfL   (npssdev), &               ! Length of diffuser (not allways used)
			   pdt   (npssdev), &               ! Update frequency of forcing variables
               diammb(npssdev), &               ! Initial diameter of bubbles
			   lambda(npssdev), &			    ! Half-width of the plume
               idetr (npssdev), &
			   STAT = istat)
   IF (istat /= 0) CALL allocate_error ( istat, 22 )

   ! ... Allocate space for input information on forcing variables pss -
   ALLOCATE (  flpss (npssdev), scpss (npssdev), &
               uEpss (npssdev), uWpss (npssdev), &
               vNpss (npssdev), vSpss (npssdev), &
               qthrs (npssdev), STAT = istat)
   IF (istat /= 0) CALL allocate_error ( istat, 22 )
   IF (ntr > 0) THEN
      ALLOCATE (trpss (npssdev,ntr), STAT=istat)
      IF (istat /= 0) CALL allocate_error ( istat, 23 )
   ENDIF

   ! ... Allocate space for variables holding column information
   ALLOCATE (  kdetr (iopss), STAT = istat)
   IF (istat /= 0) CALL allocate_error ( istat, 24 )

   !               -----Read files with pss data-----

   ! ----Loop over npssdev to Open & read files with pss data-----
   DO nn = 1, npssdev

      ! ... Construct file name (beware that nopen cannot be > 99)
      pssfile = "pss0 .txt"
      IF ( nn < 10  ) WRITE ( pssfile(5:5), FMT='(I1)' ) nn
      IF ( nn >= 10 ) WRITE ( pssfile(4:5), FMT='(I2)' ) nn

      ! ... Open IO unit
      OPEN (UNIT=i52, FILE=pssfile, STATUS="old", IOSTAT=ios)
      IF (ios /= 0) CALL open_error ( "Error opening "//pssfile, ios )

      ! Skip over first five header records in open boundary condition file
      READ (UNIT=i52, FMT='(//////)', IOSTAT=ios)
      IF (ios /= 0) CALL input_error ( ios, 48 )

      ! Read type of source-sink simulated (plumes, boundary conditions, pumped inflows)
      READ (UNIT=i52, FMT='(10X,G7.2)', IOSTAT=ios) ptype(nn)
      IF (ios /= 0) CALL input_error ( ios, 48 )

      ! Read number of points in file (it has to be equal in all arrays)
      READ (UNIT=i52, FMT='(10X,I7)', IOSTAT=ios) nptspss
      IF (ios /= 0) CALL input_error ( ios, 48 )

      ! Allocate space for the array of data first time in loop
      IF ( nn == 1) THEN
        ALLOCATE ( varspss(npssdev, ntr+2, nptspss), STAT=istat )
        IF (istat /= 0) CALL allocate_error ( istat, 24 )
      ENDIF

      ! Write the format of the data records into an internal file
      WRITE (UNIT=pssfmt, FMT='("(10X,",I3,"G11.2)")') ntr+2

      SELECT CASE (ptype(nn))

      ! ************** Section boundary inflows **************
      CASE (-2)

        ! Read how water is detrained at east face
        READ (UNIT=i52, FMT='(A)', IOSTAT=ios) commentline
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) uEpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 49 )

        ! Read how water is detrained at west face
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) uWpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 49 )

        ! Read how water is detrained at north face
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) vNpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 49 )

        ! Read how water is detrained at south face
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) vSpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 49 )

        ! Read flow threshold - to determine when it works and when not
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) qthrs(nn)
        IF (ios /= 0) CALL input_error ( ios, 49 )

        ! Read data array
        READ (UNIT=i52, FMT='(A)', IOSTAT=ios) commentline
        DO j = 1, nptspss
         READ (UNIT=i52, FMT=pssfmt, IOSTAT=ios) &
              (varspss(nn,i,j), i=1,ntr+2)
         IF (ios /= 0) CALL input_error ( ios, 49 )
        END DO

        ! ... Assign initial values to forcing variables
        flpss (nn) = varspss(nn,1,1) ! Flow rate
        scpss (nn) = varspss(nn,2,1) ! Active scalar concentration (temp.)
        IF (ntr > 0) THEN
          DO itr = 1, ntr
            trpss(nn,itr) = varspss(nn,2+itr,1) ! Tracer load
          ENDDO
        ENDIF

        ! ... Set idetr to 1 - (default value)
		idetr(nn) = 1;

        ! *************** Cell inflows ***********************
        CASE (-1)

        ! Read how water is detrained at east face
        READ (UNIT=i52, FMT='(A)', IOSTAT=ios) commentline
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) uEpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 50 )

        ! Read how water is detrained at west face
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) uWpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 50 )

        ! Read how water is detrained at north face
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) vNpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 50 )

        ! Read how water is detrained at south face
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) vSpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 50 )

        ! Read flow threshold - to determine when it works and when not
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) qthrs(nn)
        IF (ios /= 0) CALL input_error ( ios, 50 )

        ! Read data array
        READ (UNIT=i52, FMT='(A)', IOSTAT=ios) commentline
        DO j = 1, nptspss
         READ (UNIT=i52, FMT=pssfmt, IOSTAT=ios) &
              (varspss(nn,i,j), i=1,ntr+2)
         IF (ios /= 0) CALL input_error ( ios, 50 )
        END DO

        ! ... Assign initial values to forcing variables
        flpss (nn) = varspss(nn,1,1) ! Flow rate
        scpss (nn) = varspss(nn,2,1) ! Active scalar concentration (temp.)
        IF (ntr > 0) THEN
          DO itr = 1, ntr
            trpss(nn,itr) = varspss(nn,2+itr,1) ! Tracer load
          ENDDO
        ENDIF

        ! ... Set idetr to 1 - (default value)
        idetr(nn) = 1;

      ! ************* Water Pumped inflows *****************
      CASE (0)

        PRINT *, '***************** ERROR *****************'
        PRINT *, 'Water pumped inflow still NOT incorporated'
        PRINT *, '***************** ERROR *****************'
        STOP

      ! ******** Plume oxygenation models *******************
      CASE (1:)

        ! Make sure that, at least, one tracer is simulated
        ! which corresponds to oxygen
        IF (ntr < 1) THEN
          PRINT *, '***************** ERROR *****************'
          PRINT *, 'No. of tracers should be > 1 to simulate '
          PRINT *, 'the addition of oxygen through the plumes'
          PRINT *, '***************** ERROR *****************'
          STOP
        ENDIF

        ! Read how water is detrained at west face
        READ (UNIT=i52, FMT='(A)', IOSTAT=ios) commentline
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) uWpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read how water is detrained at north face
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) vNpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read how water is detrained at east face
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) uEpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read how water is detrained at south face
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) vSpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read flow threshold - to determine when it works and when not
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) qthrs(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read how water velocity at detrainment cell is calculated
        READ (UNIT=i52, FMT='(A)', IOSTAT=ios) commentline
        READ (UNIT=i52, FMT='(10X,I11)', IOSTAT=ios) idetr(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read half diffuser length (m)
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) lambda(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read initial bubble diameter (mm)
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) diammb(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read ambient salinity (constant, uS/cm)
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) salamb
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read atmospheric pressure (Pascals)
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) patm
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read constant sediment oxygen demand
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) k4sod
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read frequency of update
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) pdt(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read data array
        READ (UNIT=i52, FMT='(A)', IOSTAT=ios) commentline
        DO j = 1, nptspss
           READ (UNIT=i52, FMT=pssfmt, IOSTAT=ios) &
                (varspss(nn,i,j), i=1,ntr+2)
           IF (ios /= 0) CALL input_error ( ios, 51 )
        END DO

        ! ... Assign values to forcing variables
        !     (ntr should be >= 1 otherwise the model issues an error message)
        flpss (nn) = varspss(nn,1,1) ! Flow rate
        scpss (nn) = varspss(nn,2,1) ! Active scalar concentration (temp.) - Not used here
        DO itr = 1, ntr              ! Tracer loads -
          trpss(nn,itr) = varspss(nn,2+itr,1) ! Tracer load
        ENDDO

        ! ... Calculate diffuser length & set kdetr to default values
        ncdev = 0
        DO j = 1, iopss
          IF (iodev(j) == nn) THEN
		    kdetr(j) = km1
            ncdev    = ncdev + 1
          ENDIF
        ENDDO
        dfL(nn) = ncdev * idx

     END SELECT

     ! ... Close IO unit
     CLOSE (i52)

   ENDDO

END SUBROUTINE PointSourceSinkInput

!************************************************************************
SUBROUTINE PointSourceSinkForcing (nn,thrs)
!************************************************************************
!
!  Purpose: To assign values of flow, temperature and tracer loads
!           in sources and sinks devices
!
!------------------------------------------------------------------------

   !.....Local variables.....
   INTEGER, INTENT(IN) :: nn
   REAL, INTENT(IN) :: thrs
   REAL    :: dthrs_pss
   INTEGER :: i, j, k, l, itr

   ! ... Time (hrs) between consecutive records in files
   dthrs_pss = dtsecpss/3600.

   ! ... Get new values (NewOpenBC)
   flpss(nn)  = linear(0.,thrs,varspss(nn,1,:),dthrs_pss)
   scpss(nn)  = linear(0.,thrs,varspss(nn,2,:),dthrs_pss)
   IF (ntr > 0) THEN
     DO itr = 1, ntr
       trpss(nn,itr) = linear(0.,thrs,varspss(nn,2+itr,:),dthrs_pss)
      ENDDO
   ENDIF

END SUBROUTINE PointSourceSinkForcing

END MODULE si3d_utils
