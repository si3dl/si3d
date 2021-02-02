!************************************************************************
  MODULE si3d_BoundaryConditions
!************************************************************************
!
!  Purpose: Define Boundary Conditions to use in the solution of NS eqns.
!           Includes openbc routines, originally included in si3d_proc and
!           surfbc_heat routines created to represent heat fluxes through
!           the free surface.
!
!-------------------------------------------------------------------------

  USE omp_lib
  USE si3d_Types
  USE si3d_Utils
  USE si3d_ecomod

  IMPLICIT NONE
  SAVE

  CONTAINS

!************************************************************************
SUBROUTINE openbc0
!************************************************************************
!
!  Purpose: This routine is called at the beginning of the program
!           to open any files with boundary condition data, to
!           read the boundary condition data, and to 
!           assign the initial boundary values at time t=0.  
!
!------------------------------------------------------------------------

   !.....Local variables.....
   REAL :: areatot, uavg, vavg, sbc, qbc, hbc
   INTEGER :: i, j, k, l, nn, is, ie, js, je, ks, ke, & 
              kmx, kmy, nwlayers, ios, istat,laux
   INTEGER :: nbiid, noNGB, maxiptsNB, niNGB, icl, m1, m2
   INTEGER :: npts_wse, npts_flw, npts_sal, nptsOpenBC
   CHARACTER(LEN=14) :: openbcfmt, openbcfile
   REAL, ALLOCATABLE, DIMENSION(:,:) :: inputvar 

   ! ... Return if no open boundaries exist
   IF (nopen <= 0) RETURN

   !               -----Read files with bc data-----

   ! ---- Initialize No. of nested grid boundaries to zero
   noNGB = 0; maxiptsNB = -1;

   ! ----Loop over nopen to Open & read files with bc data-----  
   DO nn = 1, nopen

      SELECT CASE (itype(nn))

      CASE(1:3) ! Observed values at boundaries provided at txt files .....

      ! ... Construct file name (beware that nopen cannot be > 99)
      openbcfile = "openbc0 .txt"
      IF ( nn <  10 ) WRITE ( openbcfile(8:8), FMT='(I1)' ) nn
      IF ( nn >= 10 ) WRITE ( openbcfile(7:8), FMT='(I2)' ) nn

      ! ... Open IO unit
      OPEN (UNIT=i52, FILE=openbcfile, STATUS="old", IOSTAT=ios)
      IF (ios /= 0) CALL open_error ( "Error opening "//openbcfile, ios )

      ! Skip over first five header records in open boundary condition file
      READ (UNIT=i52, FMT='(////)', IOSTAT=ios)
      IF (ios /= 0) CALL input_error ( ios, 20 )

      ! Read number of points in file from seventh header record
      READ (UNIT=i52, FMT='(10X,I7)', IOSTAT=ios) nptsOpenBC
      IF (ios /= 0) CALL input_error ( ios, 21 )

      ! Allocate space for the array of data first time in loop
      IF ( nn == 1) THEN
        ALLOCATE ( varsOpenBC(nopen, ntr+2, nptsOpenBC), STAT=istat )
        IF (istat /= 0) CALL allocate_error ( istat, 16 )
      ENDIF

      ! Write the format of the data records into an internal file
      WRITE (UNIT=openbcfmt, FMT='("(10X,",I3,"G11.2)")') ntr+2
 
      ! Read data array
      DO j = 1, nptsOpenBC
         READ (UNIT=i52, FMT=openbcfmt, IOSTAT=ios) &
              (varsOpenBC(nn,i,j), i=1,ntr+2)
         IF (ios /= 0) CALL input_error ( ios, 22 )
      END DO

      ! ... Close IO unit
      CLOSE (i52)

      CASE(4:) ! Nested boundary conditions constructed in previous runs on coarser grid

        ! ... Open file
        openbcfile = "nbifilex0    "
        IF ( nn <  10) WRITE ( openbcfile(10:12), FMT='(I1,"  ")' ) nn
        IF ( nn >= 10) WRITE ( openbcfile( 9:12), FMT='(I2,"  ")' ) nn
        nbiid = nbiid0 + nn;
        OPEN(UNIT=nbiid,file=openbcfile,FORM='UNFORMATTED',IOSTAT=ios)
        IF (ios /= 0) CALL open_error ( "Error opening "//openbcfile, ios )

        !... Read type of boundary (needs to agree with input file)
        READ (nbiid) isdNBI(nn)
        !... Read time information (no. of frames) 
        READ (nbiid) nfrNBI(nn)
        !... Read spatial information (no. of cells)  
        READ (nbiid) iptNBI(nn), ntrNBI(nn)

        ! ... Check whether the length of simulations in the fine & coarse
        !     grids
        IF ( nfrNBI(nn) * dtsecOpenBC < tl) THEN
          PRINT * , 'Length of Coarse & Fine grid runs DISAGREE'
          PRINT * , 'Nested Boundary File no. = ', nn
          PRINT * , 'nframes in NB file = ', nfrNBI(nn)
          PRINT * , 'Time (s) between frames  = ', dtsecOpenBC
          PRINT * , 'Length of time simulated = ', tl 
          STOP
        ENDIF

        ! ... Increase the no. of nested grid boundaries by 1, update 
        !     max. no. of points within nested grid boundaries, and
        !     check for consistency between information from the coarse 
        !     and information required by the fine grids
        noNGB = noNGB + 1
        maxiptsNB = MAX(maxiptsNB, iptNBI(nn))
        IF ( ntrNBI(nn) .NE. ntr) THEN
          PRINT *, '******************************************************'
          PRINT *, 'STOP - No. of tracers different in coarse & fine grids'
          PRINT *, 'Please CHECK Coarse No. Tracers(',nn,')=', ntrNBI(nn)
          PRINT *, '             Fine   No. Tracers(',nn,')=', ntr
          PRINT *, '******************************************************'
          STOP
        ENDIF
        IF ( isdNBI(nn) .NE. iside(nn)) THEN
          PRINT *, '******************************************************'
          PRINT *, 'STOP - Boundary sides in coarse & fine grids dis-agree'
          PRINT *, 'Please CHECK Coarse Grid SIDE(',nn,')=', isdNBI(nn)
          PRINT *, '             Fine   Grid SIDE(',nn,')=', iside (nn)
          PRINT *, '******************************************************'
          STOP
        ENDIF

      END SELECT

   ENDDO

   ! -------------- Nested Grid Boundaries READ & CHECK -------------------
   IF ( noNGB > 0) THEN

     ! ... Allocate space for arrays holding nested grid boundaries 
     ALLOCATE ( uhNGB(maxiptsNB,noNGB), uhNGBp(maxiptsNB,noNGB), &
                vhNGB(maxiptsNB,noNGB), vhNGBp(maxiptsNB,noNGB), &
                scNGB(maxiptsNB,noNGB), scNGBp(maxiptsNB,noNGB), &
                STAT=istat)
     IF (istat /= 0) CALL allocate_error ( istat, 31 )
     ALLOCATE ( kNGB (maxiptsNB,noNGB), & 
                iNGB (maxiptsNB,noNGB), &
                jNGB (maxiptsNB,noNGB), STAT=istat)
     IF (istat /= 0) CALL allocate_error ( istat, 32 )

     ! ... Initialize arrays to zero
     uhNGB  = 0.0E0
     vhNGB  = 0.0E0
     scNGB  = 0.0E0
     uhNGBp = 0.0E0
     vhNGBp = 0.0E0
     scNGBp = 0.0E0
     kNGB   = 0

     ! ... Allocate space & initialize arrays if tracers are modelled
     IF (ntr > 0) THEN
       ALLOCATE ( trNGB (maxiptsNB,noNGB,ntr), &
                  trNGBp(maxiptsNB,noNGB,ntr), STAT=istat)
       IF (istat /= 0) CALL allocate_error ( istat, 33 )
       trNGB  = 0.0E0
       trNGBp = 0.0E0
     ENDIF

     ! Read frames 1 & 2 for nested grid boundaries
     niNGB = 0
     DO nn = 1, nopen

       IF ( itype(nn) < 4) CYCLE
       nbiid = nbiid0 + nn;
       niNGB = niNGB  + 1 ;

       ! ... Allocate space for temporary input variable array 
       ALLOCATE( inputvar ( iptNBI(nn), 5+ntr ), STAT=istat )
       IF (istat /= 0) CALL allocate_error (istat,34)

       ! ... Read FIRST FRAME & store variables ...............................
       READ(nbiid) thrsNGBp(nn), &
           ((inputvar(m1,m2),m2=1,5+ntr),m1=1,iptNBI(nn))

       SELECT CASE (iside(nn))
       CASE(1,3)
         DO icl = 1, iptNBI(nn)
           iNGB  (icl,niNGB) = inputvar(icl,1)
           jNGB  (icl,niNGB) = inputvar(icl,2)
           kNGB  (icl,niNGB) = inputvar(icl,3)
           uhNGBp(icl,niNGB) = inputvar(icl,4)
           scNGBp(icl,niNGB) = inputvar(icl,5) 
           IF (ntr > 0) THEN
             trNGBp(icl,1:ntr,niNGB) = inputvar(icl,6:5+ntr) 
           ENDIF
         ENDDO
       CASE(2,4)
         DO icl = 1, iptNBI(nn)
           iNGB  (icl,niNGB) = inputvar(icl,1)
           jNGB  (icl,niNGB) = inputvar(icl,2)
           kNGB  (icl,niNGB) = inputvar(icl,3)
           vhNGBp(icl,niNGB) = inputvar(icl,4)
           scNGBp(icl,niNGB) = inputvar(icl,5) 
           IF (ntr > 0) THEN
             trNGBp(icl,1:ntr,niNGB) = inputvar(icl,6:5+ntr) 
           ENDIF
         ENDDO
       END SELECT
!       print *,"inGB:",sum(inGB(:,:)),"jnGB:",sum(jnGB(:,:)),"KnGB:",sum(knGB(:,:))
!       print *,"uhNGBp:",sum(uhNGBp(:,:)),"vhNGBp:",sum(vhNGBp(:,:)),"scNGBp:",sum(scNGBp(:,:))
       ! ... Check grid size for consistency between coarse & fine grid .........
       SELECT CASE (iside(nn))
       CASE(1,3)
         i   = isbc(nn); 
         js  = jsbc(nn); 
         je  = jebc(nn);
         icl = 0

         
 
         DO j = js, je
           DO k = k1, kmz(ij2l(i,j))
             icl = icl + 1      
             
             IF(kNGB(icl,niNGB) .NE. k) THEN
                PRINT *, '******************************************************'
                PRINT *, 'STOOOP - Coarse & fine grids numbering NOT CONSISTENT'
                PRINT *, 'Boundary No. ', nn
                PRINT *, 'Please CHECK Coarse Grid (    k) =', iNGB(icl,niNGB), &
                &                                              jNGB(icl,niNGB), &
                &                                              kNGB(icl,niNGB)
                PRINT *, '             Fine   Grid (i,j,k) =', i,j,k
                PRINT *, '******************************************************'
                STOP
             ENDIF
           ENDDO
         ENDDO
         IF (icl .NE. iptNBI(nn) ) THEN
           PRINT *, '******************************************************'
           PRINT *, 'STOP - Coarse & fine grids numbering NOT CONSISTENT'
           PRINT *, 'Boundary No. ', nn
           PRINT *, 'Cells in output nested boundary = ', iptNBI(nn)
           PRINT *, 'Cells in fine grid at  boundary = ', icl
           PRINT *, '******************************************************'
           STOP
         ENDIF

       CASE(2,4)
         j   = jsbc(nn); 
         is  = isbc(nn); 
         ie  = iebc(nn);
         icl = 0
         DO i = is, ie
           DO k = k1, kmz(ij2l(i,j))
             icl = icl + 1      
             IF(kNGB(icl,niNGB) .NE. k) THEN
                PRINT *, '******************************************************'
                PRINT *, 'STOP - Coarse & fine grids numbering NOT CONSISTENT'
                PRINT *, 'Boundary No. ', nn
                PRINT *, 'Please CHECK Coarse Grid (i,j,k) =', iNGB(icl,niNGB), &
                &                                              jNGB(icl,niNGB), &
                &                                              kNGB(icl,niNGB)
                PRINT *, '             Fine   Grid (i,j,k) =', i,j,k
                PRINT *, '******************************************************'
                STOP
             ENDIF
           ENDDO
         ENDDO
         IF (icl .NE. iptNBI(nn) ) THEN
           PRINT *, '******************************************************'
           PRINT *, 'STOP - Coarse & fine grids numbering NOT CONSISTENT'
           PRINT *, 'Boundary No. ', nn
           PRINT *, 'Cells in output nested boundary = ', iptNBI(nn)
           PRINT *, 'Cells in fine grid at  boundary = ', icl
           PRINT *, '******************************************************'
           STOP
         ENDIF

       END SELECT      
       print *,"tttantes:",thrsNGB(nn)
       ! ... Read SECOND FRAME & store variables ...............................
       READ(nbiid) thrsNGB(nn), & 
                 ((inputvar(m1,m2),m2=4,5+ntr),m1=1,iptNBI(nn))
       print *,"ttt:",thrsNGB(nn)
       SELECT CASE (iside(nn))
       CASE(1,3)
         DO icl = 1, iptNBI(nn)
           uhNGB(icl,niNGB) = inputvar(icl,4)
           scNGB(icl,niNGB) = inputvar(icl,5) 
           IF (ntr > 0) THEN
             trNGB(icl,1:ntr,niNGB) = inputvar(icl,6:5+ntr) 
           ENDIF
         ENDDO
       CASE(2,4)
         DO icl = 1, iptNBI(nn)
           vhNGB(icl,niNGB) = inputvar(icl,4)
           scNGB(icl,niNGB) = inputvar(icl,5) 
           IF (ntr > 0) THEN
             trNGB(icl,1:ntr,niNGB) = inputvar(icl,6:5+ntr) 
           ENDIF
         ENDDO
       END SELECT

       DEALLOCATE (inputvar)
!       print *,"uhNGB:",sum(uhNGB(:,:)),"vhNGB:",sum(vhNGB(:,:)),"scNGB:",sum(scNGB(:,:))
     ENDDO

   ENDIF         

   !           -----Assign boundary values at time t=0.0-----

   !.....Initialize to zero boundary flows & thickness at time n+1
   uhEB   = 0.0  ; huEB   = 0.0  ;
   uhWB   = 0.0  ; huWB   = 0.0  ;
   vhNB   = 0.0  ; hvNB   = 0.0  ;
   vhSB   = 0.0  ; hvSB   = 0.0  ;
   niNGB  = 0    ;

   !.....Loop over open boundaries..........................................
   DO nn = 1, nopen

      SELECT CASE ( itype(nn) )
 
      ! ..... CASE 1 -- wse specified ......................................
      CASE (1)

         ! Get first boundary value of zeta & make sure the i,j locations
         ! are wett boundary cells (Assume first value applies for t=0.0. 
         ! It should be consistent with the initial condition for zeta 
         ! defined in SUBROUTINE init)
         sbc = varsOpenBC(nn,1,1);
        
         ! Identify wse boundary as on the west, north, east, or south
         SELECT CASE ( iside(nn) )

         ! West boundary
         CASE (1)

            ! ... Retrieve (i,j) index where BC is specified
            i = isbc(nn); js = jsbc(nn); je = jebc(nn)
            
            
            ! ... Make sure i,j locations are wett boundary cells
            DO j = js, je   
              laux=ij2l(i,j)
              ! ... Assign boundary condition
              sp(laux) = sbc
              IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i-1,j) ) THEN
                PRINT *, "  "
                PRINT *, " ****STOPPING - West bdry. not a bdry. point"
                STOP 
              END IF

              nwlayers = (kmz(laux) - k1z(laux)) + 1
              IF(nwlayers <= 0) THEN  ! dry point on boundary is not allowed
                PRINT *, " ERROR--dry point on west boundary" 
                PRINT *, "  "
                PRINT *, "  "
                PRINT *, " ****STOPPING si3d for boundary condition error"
                STOP 
              END IF
            END DO

         ! North boundary
         CASE (2)

            ! ... Retrieve (i,j) index where BC is specified
            j = jsbc(nn); is = isbc(nn); ie = iebc(nn)
            
            ! ... Make sure i,j locations are wett boundary cells
            DO i = is, ie   
              laux=ij2l(i,j)
              ! ... Assign boundary condition
              sp(laux) = sbc
              IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i,j+1) ) THEN
                PRINT *, "  "
                PRINT *, " ****STOPPING - North bdry. not a bdry. point"
                STOP 
              ENDIF
              nwlayers = (kmz(laux) - k1z(laux)) + 1
              IF (nwlayers < 1) THEN  ! dry point on boundary is not allowed
                PRINT *, " ERROR--dry point on north boundary" 
                PRINT *, "  "
                PRINT *, "  "
                PRINT *, " ****STOPPING si3d for boundary condition error"
                STOP 
              ENDIF
            END DO

         ! East boundary
         CASE (3)

            ! ... Retrieve (i,j) index where BC is specified
            i = isbc(nn); js = jsbc(nn); je = jebc(nn)
            
            ! ... Make sure i,j locations are wett boundary cells
            DO j = js, je   
              laux=ij2l(i,j)
              ! ... Assign boundary condition
              sp(laux) = sbc
              IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i+1,j) ) THEN
                PRINT *, "  "
                PRINT *, " ****STOPPING - East bdry. not a bdry. point"
                STOP 
              ENDIF
              nwlayers = (kmz(laux) - k1z(laux)) + 1
              IF ( nwlayers < 1) THEN ! dry point on boundary is not allowed

                PRINT *, " ERROR--dry point on east boundary" 




                PRINT *, "  "
                PRINT *, "  "
                PRINT *, " ****STOPPING si3d for boundary condition error"
                STOP 
              ENDIF
            ENDDO

         ! South boundary
         CASE (4)

            ! ... Retrieve (i,j) index where BC is specified
            j = jsbc(nn); is = isbc(nn); ie = iebc(nn)
            
            ! ... Make sure i,j locations are wett boundary cells
            DO i = is, ie 
              laux=ij2l(i,j)
              ! ... Assign boundary condition
              sp(laux) = sbc
              IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i,j-1) ) THEN
                PRINT *, "  "
                PRINT *, " ****STOPPING - South bdry. not a bdry. point"
                STOP 
              ENDIF
              nwlayers = (kmz(laux) - k1z(laux)) + 1
              IF (nwlayers < 1) THEN ! dry point on boundary is not allowed
                PRINT *, " ERROR--dry point on south boundary" 
                PRINT *, "  "
                PRINT *, "  "
                PRINT *, " ****STOPPING si3d for boundary condition error"
                STOP 
              ENDIF
            END DO 

         END SELECT

      !..... CASE 2 -- Free surface flow specified .........................
      CASE (2)
      
         ! Get first boundary value of flow (in units of m**3/sec)
         ! (Assume first value applies for t=0.0. It should be consistent
         !  with the initial condition for uh, vh, u,and v)
         qbc = varsOpenBC(nn,1,1);

         ! Identify boundary as on the west, north, east, or south
         SELECT CASE ( iside(nn) )

         ! West boundary
         CASE (1)
 
            ! Get i-, j- indexes for bdry. point
            i = isbc(nn); js = jsbc(nn); je = jebc(nn)

            ! Make sure i,j locations are wett boundary cells
            DO j = js, je   
              IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i-1,j) ) THEN
                 PRINT *, "  "
                 PRINT *, " ****STOPPING - West bdry. not a bdry. point"
                 STOP
              END IF
            ENDDO

            ! Get max. depth in section with bdry. values specified
            hbc = -1.e6; 
            DO j = js,je
              laux=ij2l(i,j)
              IF(hhs(laux)>hbc) hbc = hhs(laux)
            ENDDO

            ! ... Define free surface location at bdry.
            laux=ij2l(i,j)
            sbc = MAX(s(laux), -hbc+dzmin)

            ! ... Define thickness of bdry. wet cells & total area
            areatot = 0.0; huWB(:,js:je) = ZERO
            DO j = js,je
              laux=ij2l(i,j)
              kmx = kmz(laux)
              DO k = k1, kmx
                huWB (k,j)=AMIN1(zlevel(k+1),hhs(laux)) -        &
                &          AMAX1(zlevel(  k),-sbc)
                IF(huWB(k,j) <= HMIN) huWB(k,j) = ZERO;
                areatot = areatot + huWB(k,j) * dy
              ENDDO
            ENDDO

            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            uavg = qbc/areatot; uhWB(:,js:je) = 0.0E0
            DO j = js, je
              laux=ij2l(i,j)
              kmx = kmz(laux)
              DO k = k1, kmx
                uhWB(k,j) = uavg * huWB(k,j)
              END DO
            ENDDO

         ! East boundary
         CASE (3)
 
            ! ... Get i-, j- indexes for bdry. point
            i = isbc(nn); js = jsbc(nn); je = jebc(nn)
            ! ... Make sure i,j locations are wett boundary cells
            DO j = js, je   
              IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i+1,j) ) THEN
                PRINT *, "  "
                PRINT *, " ****STOPPING - East bdry. not a bdry. point"
                STOP
              END IF
            ENDDO
            ! ... Get max. depth in section with bdry. values specified
            hbc = -1.e6; 
            DO j = js,je
              laux=ij2l(i,j)
              IF(hhs(laux)>hbc) hbc = hhs(laux)
            ENDDO
            ! ... Define free surface location at bdry.
            sbc = MAX(s(ij2l(i,j)), -hbc+dzmin)
            ! ... Define thickness of bdry. wet cells & total area
            areatot = 0.0; huEB(:,js:je) = ZERO
            DO j = js,je
              kmx = kmz(ij2l(i,j))
              DO k = k1, kmx
                huEB (k,j)=AMIN1(zlevel(k+1),hhs(ij2l(i,j))) -        &
                &          AMAX1(zlevel(  k),-sbc)
                IF(huEB(k,j) <= HMIN) huEB(k,j) = ZERO;
                areatot = areatot + huEB(k,j) * dy
              ENDDO
            ENDDO
            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            uavg = qbc/areatot; uhEB(:,js:je) = 0.0E0
            DO j = js, je
              kmx = kmz(ij2l(i,j))
              DO k = k1, kmx
                 uhEB(k,j) = uavg * huEB(k,j)
              END DO
            ENDDO
            

         ! North boundary
         CASE (2)
 
            ! ... Get i-, j- indexes for bdry. point
            j = jsbc(nn); is = isbc(nn); ie = iebc(nn)
            ! ... Make sure i,j locations are wett boundary cells
            DO i = is, ie   
              IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i,j+1) ) THEN
                PRINT *, "  "
                PRINT *, " ****STOPPING - North bdry. not a bdry. point"
                STOP 
              END IF
            ENDDO
            ! ... Get max. depth in section with bdry. values specified
            hbc = -1.e6; 
            DO i = is,ie
              laux=ij2l(i,j)
              IF(hhs(laux)>hbc) hbc = hhs(laux)
            ENDDO
            ! ... Define free surface location at bdry.
            sbc = MAX(s(ij2l(i,j)), -hbc+dzmin)
            ! ... Define thickness of bdry. wet cells & total area
            areatot = 0.0; hvNB(:,is:ie) = ZERO
            DO i = is,ie
              kmy = kmz(ij2l(i,j))
              DO k = k1, kmy
                hvNB (k,i)=AMIN1(zlevel(k+1),hhs(ij2l(i,j))) -        &
                &          AMAX1(zlevel(  k),-sbc)
                IF(hvNB(k,i) <= HMIN) hvNB(k,i) = ZERO;
                areatot = areatot + hvNB(k,i) * dy
              ENDDO
            ENDDO
            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            vavg = qbc/areatot; vhNB(:,is:ie) = 0.0E0
            DO i = is, ie
              laux=ij2l(i,j)
              kmy = kmz(laux)
              DO k = k1z(laux), kmy
                 vhNB(k,i) = vavg * hvNB(k,i)
              END DO
            ENDDO

         ! South boundary
         CASE (4)
 
            ! ... Get i-, j- indexes for bdry. point
            j = jsbc(nn); is = isbc(nn); ie = iebc(nn)
            ! ... Make sure i,j locations are wett boundary cells
            DO i = is, ie   
              IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i,j-1) ) THEN
                PRINT *, "  "
                PRINT *, " ****STOPPING - South bdry. not a bdry. point"
                STOP 
              END IF
            ENDDO
            ! ... Get max. depth in section with bdry. values specified
            hbc = -1.e6; 
            DO i = is,ie
              laux=ij2l(i,j)
              IF(hhs(laux)>hbc) hbc = hhs(laux)
            ENDDO
            ! ... Define free surface location at bdry.
            sbc = MAX(s(ij2l(i,j)), -hbc+dzmin)
            ! ... Define thickness of bdry. wet cells & total area
            areatot = 0.0; hvSB(:,is:ie) = ZERO
            DO i = is,ie
              kmy = kmz(ij2l(i,j))
              DO k = k1, kmy
                hvSB (k,i)=AMIN1(zlevel(k+1),hhs(ij2l(i,j))) -        &
                &          AMAX1(zlevel(  k),-sbc)
                IF(hvSB(k,i) <= HMIN) hvSB(k,i) = ZERO;
                areatot = areatot + hvSB(k,i) * dy
              ENDDO
            ENDDO
            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            vavg = qbc/areatot; vhSB(:,is:ie) = 0.0E0
            DO i = is, ie
              kmy = kmz(ij2l(i,j))
              DO k = k1, kmy
                 vhSB(k,i) = vavg * hvSB(k,i)
              END DO
            ENDDO

         END SELECT         

      !..... CASE 3 -- Submerged flow specified.............................
      CASE (3)
      
         ! Get first boundary value of flow (in units of m**3/sec)
         ! (Assume first value applies for t=0.0. It should be consistent
         !  with the initial condition for uh, vh, u,and v)
         qbc = varsOpenBC(nn,1,1);

         ! Identify boundary as on the west, north, east, or south
         SELECT CASE ( iside(nn) )

         ! West boundary
         CASE (1)
 
            ! ... Get i-, j- indexes for bdry. point
            i  = isbc(nn); 
            j  = jsbc(nn); 
            l  = ij2l(i,j)
            ks = iebc(nn); 
            ke = jebc(nn); 
            ! ... Make sure i,j locations are wett boundary cells
            IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i-1,j) ) THEN
              PRINT *, "  "
              PRINT *, " ****STOPPING - West bdry. not a bdry. point"
              STOP
            END IF
            ! ... Make sure k- location is wett, define thickness of bdry. cell
            !     & total area for outflow-inflow section
            areatot = 0.0
            DO k = ks, ke
              huWB (k,j)= hp(k,l) 
              IF ( huWB(k,j) <= ZERO ) THEN 
                PRINT *, "  "
                PRINT *, " ****STOPPING - Spec. West bdry. below free surface"
                STOP
              ENDIF
              areatot = areatot + huWB(k,j) * dy
            ENDDO
            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            uavg = qbc/areatot; 
            DO k = ks, ke
              uhWB(k,j) = uavg * huWB(k,j)
            ENDDO

         ! East boundary
         CASE (3)
 
            ! ... Get i-, j- indexes for bdry. point
            i  = isbc(nn); 
            j  = jsbc(nn); 
            l  = ij2l(i,j)
            ks = iebc(nn); 
            ke = jebc(nn); 

            ! ... Make sure i,j locations are wett boundary cells
            IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i+1,j) ) THEN
              PRINT *, "  "
              PRINT *, " ****STOPPING - East bdry. not a bdry. point"
              STOP
            END IF
            ! ... Make sure k- location is wett, define thickness of bdry. cell
            !     & total area for outflow-inflow section
            areatot = 0.0
            DO k = ks, ke
              huEB (k,j)= hp(k,l) 
              IF ( huEB(k,j) <= ZERO ) THEN 
                PRINT *, "  "
                PRINT *, " ****STOPPING - East bdry. below free surface"
                STOP
              ENDIF
              areatot = areatot + huEB(k,j) * dy
            ENDDO
            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            uavg = qbc/areatot; 
            DO k = ks, ke
              uhEB(k,j) = uavg * huEB(k,j)
            ENDDO

         ! North boundary
         CASE (2)
 
            ! ... Get i-, j- indexes for bdry. point
            j  = jsbc(nn); 
            i  = isbc(nn);
            l  = ij2l(i,j) 
            ks = iebc(nn);
            ke = jebc(nn);
            ! ... Make sure i,j locations are wett boundary cells
            IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i,j+1) ) THEN
              PRINT *, "  "
              PRINT *, " ****STOPPING - North bdry. not a bdry. point"
              STOP 
            END IF
            ! ... Make sure k- location is wett, define thickness of bdry. cell
            !     & total area for outflow-inflow section
            areatot = 0.0
            DO k = ks, ke
              hvNB (k,i)= hp(k,l) 
              IF ( hvNB(k,i) <= ZERO ) THEN 
                PRINT *, "  "
                PRINT *, " ****STOPPING - Submerged North bdry. on a DRY CELL"
                STOP
              ENDIF
              areatot = areatot + hvNB(k,i) * dx
            ENDDO
            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            vavg = qbc/areatot; 
            DO k = ks, ke
              vhNB(k,i) = vavg * hvNB(k,i)
            ENDDO

         ! South boundary
         CASE (4)
 
            ! ... Get i-, j- indexes for bdry. point
            j  = jsbc(nn); 
            i  = isbc(nn);
            l  = ij2l(i,j) 
            ks = iebc(nn);
            ke = jebc(nn);
            ! ... Make sure i,j locations are wett boundary cells
            IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i,j-1) ) THEN
              PRINT *, "  "
              PRINT *, " ****STOPPING - South bdry. not a bdry. point"
              STOP 
            END IF
            ! ... Make sure k- location is wett, define thickness of bdry. cell
            !     & total area for outflow-inflow section
            areatot = 0.0
            DO k = ks, ke
              hvSB (k,i)= hp(k,l) 
              IF ( hvSB(k,i) <= ZERO ) THEN 
                PRINT *, "  "
                PRINT *, " ****STOPPING - Submerged South bdry. on a DRY CELL"
                STOP
              ENDIF
              areatot = areatot + hvSB(k,i) * dx
            ENDDO

            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            vavg = qbc/areatot; 
            DO k = ks, ke
              vhSB(k,i) = vavg * hvSB(k,i)
            ENDDO

         END SELECT

      !..... CASE 4 -- Nested grid boundaries ..............................
      CASE (4)

        niNGB = niNGB + 1;
        SELECT CASE (iside(nn))
        CASE(1)
          i   = isbc(nn); 
          js  = jsbc(nn); 
          je  = jebc(nn);
          icl = 0
          DO j = js, je; 
            l   = ij2l(i,j)
            IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i-1,j) ) THEN
              PRINT *, "  "
              PRINT *, " ****STOPPING - West bdry. not a bdry. point"
              STOP
            END IF
            DO k = k1, kmz(l)
              icl = icl + 1      
              uhWB(k,j) = uhNGBp(icl,niNGB)
              huWB(k,j) = hp(k,l)
            ENDDO 
          ENDDO 
        CASE(3)
          i   = isbc(nn); 
          js  = jsbc(nn); 
          je  = jebc(nn);
          icl = 0
          DO j = js, je; 
            l  = ij2l(i,j)
            IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i+1,j) ) THEN
              PRINT *, "  "
              PRINT *, " ****STOPPING - East bdry. not a bdry. point"
              STOP
            END IF
            DO k = k1, kmz(l)
              icl = icl + 1      
              uhEB(k,j) = uhNGBp(icl,niNGB)
              huEB(k,j) = hp(k,l)
            ENDDO
          ENDDO
        CASE(2)
          j   = jsbc(nn); 
          is  = isbc(nn); 
          ie  = iebc(nn);
          icl = 0
          DO i = is, ie; 
            l  = ij2l(i,j)
            IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i,j+1) ) THEN
              PRINT *, "  "
              PRINT *, " ****STOPPING - North bdry. not a bdry. point"
              STOP
            END IF
            DO k = k1, kmz(l)
              icl = icl + 1      
              vhNB(k,i) = vhNGBp(icl,niNGB)
              hvNB(k,i) = hp(k,l)
            ENDDO
          ENDDO
        CASE(4)
          j   = jsbc(nn); 
          is  = isbc(nn); 
          ie  = iebc(nn);
          icl = 0
          DO i = is, ie; 
            l  = ij2l(i,j)
            IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i,j-1) ) THEN
              PRINT *, "  "
              PRINT *, " ****STOPPING - South bdry. not a bdry. point"
              STOP
            END IF


            DO k = k1, kmz(l)
              icl = icl + 1       
              vhSB(k,i) = vhNGBp(icl,niNGB)
              hvSB(k,i) = hp(k,l)
            ENDDO
          ENDDO
        END SELECT   
          
      END SELECT
!      print *,"uhWB:",sum(uhWB(:,:)),"huWB:",sum(huWB(:,:))
   END DO

   ! ... Initialize flow variables at time n-1 & n
   uhEBp  = uhEB ; huEBp  = huEB ;
   uhWBp  = uhWB ; huWBp  = huWB ;
   vhNBp  = vhNB ; hvNBp  = hvNB ;
   vhSBp  = vhSB ; hvSBp  = hvSB ;
   uhEBpp = uhEBp; huEBpp = huEBp;
   uhWBpp = uhWBp; huWBpp = huWBp;
   vhNBpp = vhNBp; hvNBpp = hvNBp;
   vhSBpp = vhSBp; hvSBpp = hvSBp;

END SUBROUTINE openbc0

!************************************************************************
SUBROUTINE openbcUVH(thrs)
!************************************************************************
!
!  Purpose: To assign values of water surface elevation or velocity
!           along open boundaries at the new (n+1) time level. 
!
!------------------------------------------------------------------------

   REAL, INTENT(IN) :: thrs

   !.....Local variables.....
   REAL    :: areatotal, uavg, vavg, sbc, qbc, hbc, dthrs_wse, dthrs_flw
   INTEGER :: i, j, k, l, ios, nn, is, ie, js, je, kmx, kmy, ks, ke
   INTEGER :: icl, m1, m2, niNGB, nbid, istat,iaux,jaux,laux,no,ide_t
   REAL    :: weight


   ide_t = omp_get_thread_num()+1
   ! ... Return if no open boundaries exist
   IF (nopenH(ide_t) <= 0) RETURN
   
   ! ... Initialize to zero counter for nested grid boundaries
   niNGB = 0;

   !.....Loop over open boundaries.....
   DO nn = 1, nopenHH(ide_t)
      no = noh2noH(nn,ide_t)
      SELECT CASE ( itype(no) )
 
      !.....Case 1 -- wse specified.......................................
      CASE (1)

         ! Get new boundary value of zeta (NewOpenBC)
         dthrs_wse = dtsecOpenBC/3600.
         sbc = parab(0.,thrs,varsOpenBC(no,1,:),dthrs_wse)

         ! Identify wse boundary as on the west, north, east, or south
         SELECT CASE ( iside(no) )
 
         ! West boundary
         CASE (1)
            i  = isbcHH(no,ide_t); 
            js = jsbcH(no,ide_t); 
            je = jebcH(no,ide_t)
            DO jaux=js,je
            	s(ij2l(i,jaux)) = sbc
            END DO
            
         ! North boundary
         CASE (2)
            j  = jsbcH(no,ide_t); 
            is = isbcHH(no,ide_t); 
            ie = iebcHH(no,ide_t)
            DO iaux=is,ie
            	s(ij2l(iaux,j)) = sbc
            END DO
            
         ! East boundary
         CASE (3)
            i  = isbcHH(no,ide_t); 
            js = jsbcH(no,ide_t); 
            je = jebcH(no,ide_t)
            DO jaux=js,je
            	s(ij2l(i,jaux)) = sbc
            END DO
            
         ! South boundary              
         CASE (4)
            j  = jsbcH(no,ide_t); 
            is = isbcHH(no,ide_t);  
            ie = iebcHH(no,ide_t)          
            DO iaux=is,ie
            	s(ij2l(iaux,j)) = sbc
            END DO
            
         END SELECT

      !.....Case 2 -- Free surface flow specified........................
      CASE (2)
      
         ! Get new boundary value of flow (in units of m**3/sec) (NewOpenBC)
         dthrs_flw = dtsecOpenBC/3600.
         qbc = parab(0.,thrs,varsOpenBC(no,1,:),dthrs_flw)

         ! Identify boundary as on the west, north, east, or south
         SELECT CASE ( iside(no) )

         ! West boundary
         CASE (1)
 
            ! Get i-, j- indexes for bdry. point
            i  = isbcH(no,ide_t); 
            js = jsbcH(no,ide_t); 
            je = jebcH(no,ide_t)

            ! Get max. depth in section with bdry. values specified
            hbc = -1.e6; 
            DO j = js,je
              laux = ij2l(i,j)
              IF(hhs(laux)>hbc) hbc = hhs(laux)
            ENDDO

            ! ... Define free surface location at bdry.
            sbc = MAX(s(ij2l(i,j)), -hbc+dzmin)

            ! ... Define thickness of bdry. wet cells & total area
            areatot(no) = 0.0

            huWB(:,js:je) = ZERO
            DO j = js,je
              laux = ij2l(i,j)
              kmx = kmz(laux)
              DO k = k1, kmx
                huWB (k,j)=AMIN1(zlevel(k+1),hhs(laux)) -        &
                &          AMAX1(zlevel(  k),-sbc)
                IF(huWB(k,j) <= HMIN) huWB(k,j) = ZERO;
                areatot(no) = areatot(no) + huWB(k,j) * dy
              ENDDO
            ENDDO

            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            uavg = qbc/areatot(no); uhWB(:,js:je) = 0.0E0
            DO j = js, je
              kmx = kmz(ij2l(i,j))
              DO k = k1, kmx
                 uhWB(k,j) = uavg * huWB(k,j)
				 
              END DO
            ENDDO

         ! East boundary
         CASE (3)
 
            ! Get i-, j- indexes for bdry. point
            i  = isbcH(no,ide_t); 
            js = jsbcH(no,ide_t); 
            je = jebcH(no,ide_t);

            ! Get max. depth in section with bdry. values specified
            hbc = -1.e6; 
            DO j = js,je
              laux = ij2l(i,j)
              IF(hhs(laux)>hbc) hbc = hhs(laux)
            ENDDO

            ! ... Define free surface location at bdry.
            sbc = MAX(s(ij2l(i,j)), -hbc+dzmin)

            ! ... Define thickness of bdry. wet cells & total area
            areatot(no) = 0.0; huEB(:,js:je) = ZERO
            DO j = js,je
              laux = ij2l(i,j)
              kmx = kmz(laux)
              DO k = k1, kmx
                huEB (k,j)=AMIN1(zlevel(k+1),hhs(laux)) -        &
                &          AMAX1(zlevel(  k),-sbc)
                IF(huEB(k,j) <= HMIN) huEB(k,j) = ZERO;
                areatot(no) = areatot(no) + huEB(k,j) * dy
              ENDDO
            ENDDO

            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            uavg = qbc/areatot(no); uhEB(:,js:je) = 0.0E0
            DO j = js, je
              kmx = kmz(ij2l(i,j))
              DO k = k1, kmx
                 uhEB(k,j) = uavg * huEB(k,j)
              END DO
            ENDDO

         ! North boundary
         CASE (2)
 
            ! Get i-, j- indexes for bdry. point
            j  = jsbcH(no,ide_t); 
            is = isbcH(no,ide_t); 
            ie = iebcH(no,ide_t)

            ! Get max. depth in section with bdry. values specified
            hbc = -1.e6; 
            DO i = is,ie
              laux = ij2l(i,j)
              IF(hhs(laux)>hbc) hbc = hhs(laux)
            ENDDO

            ! ... Define free surface location at bdry.
            sbc = MAX(s(ij2l(i,j)), -hbc+dzmin)

            ! ... Define thickness of bdry. wet cells & total area
            !$omp critical
            if(flag(no) .EQ. 0) THEN
                areatot(no) = 0.0; 
                flag(no) = 1
            end if
                
            !$omp end critical
            !hvNB(:,is:ie) = ZERO
            DO i = is,ie
              laux = ij2l(i,j)
              kmy = kmz(laux)
              DO k = k1, kmy
                hvNB (k,i)=AMIN1(zlevel(k+1),hhs(laux)) -        &
                &          AMAX1(zlevel(  k),-sbc)
                IF(hvNB(k,i) <= HMIN) hvNB(k,i) = ZERO;
                !$omp critical
                areatot(no) = areatot(no) + hvNB(k,i) * dy
                !$omp end critical
              ENDDO
            ENDDO
                !$omp critical
                    cont(no) = cont(no) + 1
                !$omp end critical
            do while(cont(no) < nopth(no))

            end do
            flag(no) = 0
            cont(no) = 0 
            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            vavg = qbc/areatot(no); vhNB(:,is:ie) = 0.0E0
            DO i = is, ie
              kmy = kmz(ij2l(i,j))
              DO k = k1, kmy
                 vhNB(k,i) = vavg * hvNB(k,i)
               END DO
            ENDDO

         ! South boundary
         CASE (4)
 
            ! Get i-, j- indexes for bdry. point
            j  = jsbcH(no,ide_t); 
            is = isbcH(no,ide_t); 
            ie = iebcH(no,ide_t)

            ! Get max. depth in section with bdry. values specified
            hbc = -1.e6; 
            DO i = is,ie
              laux = ij2l(i,j)
              IF(hhs(laux)>hbc) hbc = hhs(laux)
            ENDDO

            ! ... Define free surface location at bdry.
            sbc = MAX(s(ij2l(i,j)), -hbc+dzmin)

            ! ... Define thickness of bdry. wet cells & total area
            !$omp critical
            if(flag(no) .EQ. 0) THEN
                areatot(no) = 0.0; 
                flag(no) = 1
            end if
                
            !$omp end critical
            !hvSB(:,is:ie) = ZERO
            DO i = is,ie
              laux = ij2l(i,j)
              kmy = kmz(laux)
              DO k = k1, kmy
                hvSB (k,i)=AMIN1(zlevel(k+1),hhs(laux)) -        &
                &          AMAX1(zlevel(  k),-sbc)
                IF(hvSB(k,i) <= HMIN) hvSB(k,i) = ZERO;
                !$omp critical
                areatot(no) = areatot(no) + hvNB(k,i) * dy
                !$omp end critical
              ENDDO
            ENDDO
            !$omp critical
                    cont(no) = cont(no) + 1
                !$omp end critical
            do while(cont(no) < nopth(no))

            end do
            flag(no) = 0
            cont(no) = 0 
            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            vavg = qbc/areatot(no); vhSB(:,is:ie) = 0.0E0
            DO i = is, ie
              kmy = kmz(ij2l(i,j))
              DO k = k1, kmy
                 vhSB(k,i) = vavg * hvSB(k,i)
              END DO
            ENDDO
         END SELECT

      !.....Case 3 -- Submerged flow specified...........................
      CASE (3)

         ! Get first boundary value of flow (in units of m**3/sec)
         ! (Assume first value applies for t=0.0. It should be consistent
         !  with the initial condition for uh, vh, u,and v)
         ! Get new boundary value of flow (in units of m**3/sec) (NewOpenBC)
         dthrs_flw = dtsecOpenBC/3600.
         qbc = parab(0.,thrs,varsOpenBC(no,1,:),dthrs_flw)

         ! Identify boundary as on the west, north, east, or south
         SELECT CASE ( iside(no) )

         ! West boundary
         CASE (1)
 
            ! Get i-, j- indexes for bdry. point
            i  = isbc(no); 
            j  = jsbc(no); 
            l  = ij2l(i,j)
            ks = iebc(no); 
            ke = jebc(no); 

            ! Make sure k- location is wett, define thickness of bdry. cell
            ! & total area for outflow-inflow section
            areatotal = 0.0
            DO k = ks, ke
              huWB (k,j)= hp(k,l) 
              IF ( huWB(k,j) <= ZERO ) THEN 
                PRINT *, "  "
                PRINT *, " ****STOPPING - Submerged West bdry. on DRY CELL"
                STOP
              ENDIF
              areatotal = areatotal + huWB(k,j) * dy
            ENDDO

            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            uavg = qbc/areatotal; 
            DO k = ks, ke
              uhWB(k,j) = uavg * huWB(k,j)       
            ENDDO

         ! East boundary
         CASE (3)
 
            ! Get i-, j- indexes for bdry. point
            i  = isbc(no); 
            j  = jsbc(no); 
            l  = ij2l(i,j)
            ks = iebc(no); 
            ke = jebc(no); 

            ! Make sure k- location is wett, define thickness of bdry. cell
            ! & total area for outflow-inflow section
            areatotal = 0.0
            DO k = ks, ke
              huEB (k,j)= hp(k,l) 
              IF ( huEB(k,j) <= ZERO ) THEN 
                PRINT *, "  "
                PRINT *, " ****STOPPING - Submerged East bdry. on DRY CELL"
                STOP
              ENDIF
              areatotal = areatotal + huEB(k,j) * dy
            ENDDO

            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            uavg = qbc/areatotal; 
            DO k = ks, ke
              uhEB(k,j) = uavg * huEB(k,j)
            ENDDO

         ! North boundary
         CASE (2)
 
            ! Get i-, j- indexes for bdry. point
            j  = jsbc(no); 
            i  = isbc(no);
            l  = ij2l(i,j) 
            ks = iebc(no);
            ke = jebc(no);

            ! Make sure k- location is wett, define thickness of bdry. cell
            ! & total area for outflow-inflow section
            areatotal = 0.0
            DO k = ks, ke
              hvNB (k,i)= hp(k,l) 
              IF ( hvNB(k,i) <= ZERO ) THEN 
                PRINT *, "  "
                PRINT *, " ****STOPPING - Submerged North bdry. on a DRY CELL"
                STOP
              ENDIF
              areatotal = areatotal + hvNB(k,i) * dx
            ENDDO

            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            vavg = qbc/areatotal; 
            DO k = ks, ke
              vhNB(k,i) = vavg * hvNB(k,i)
            ENDDO

         ! South boundary
         CASE (4)
 
            ! Get i-, j- indexes for bdry. point
            j  = jsbc(no); 
            i  = isbc(no);
            l  = ij2l(i,j) 
            ks = iebc(no);
            ke = jebc(no);

            ! Make sure k- location is wett, define thickness of bdry. cell
            ! & total area for outflow-inflow section
            areatotal = 0.0
            DO k = ks, ke
              hvSB (k,i)= hp(k,l) 
              IF ( hvSB(k,i) <= ZERO ) THEN 
                PRINT *, "  "
                PRINT *, " ****STOPPING - Submerged South bdry. on a DRY CELL"
                STOP
              ENDIF
              areatotal = areatotal + hvSB(k,i) * dx
            ENDDO

            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            vavg = qbc/areatotal; 
            DO k = ks, ke
              vhSB(k,i) = vavg * hvSB(k,i)
            ENDDO

         END SELECT    

      !.....Case 4 -- Nested grid boundaries specified ..................
      CASE (4)
      
         ! ... Update counter of nested grid boundaries
         niNGB = niNGB + 1;

         ! ... Define weighting coefficients for records 
         weight  = (thrs - thrsNGBp(no))/(thrsNGB(no)-thrsNGBp(no))
!         print *,"niNGB:",niNGB,"thrs",thrs,"thrsNGBp:",thrsNGBp
!         print *,"thrsNGB:",thrsNGB,"weight:",weight
         ! Identify boundary as on the west, north, east, or south
         SELECT CASE (iside(no))
         CASE(1)
           i   = isbcH(no,ide_t); 
           js  = jsbcH(no,ide_t); 
           je  = jebcH(no,ide_t);
           icl = 0
           DO j = js, je;
             l = ij2l(i,j) 
             DO k = k1, kmz(l)
               icl = icl + 1      
               uhWB(k,j) = uhNGB (icl,no)*    weight + &
                           uhNGBp(icl,no)*(1.-weight)
               huWB(k,j)= h(k,l)
             ENDDO 
           ENDDO 
         CASE(3)
           i   = isbcH(no,ide_t); 
           js  = jsbcH(no,ide_t); 
           je  = jebcH(no,ide_t);
           icl = 0
           DO j = js, je; 
             l = ij2l(i,j) 
             DO k = k1, kmz(l)
               icl = icl + 1      
               uhEB(k,j) = uhNGB (icl,no)*    weight + &
                           uhNGBp(icl,no)*(1.-weight)
               huEB(k,j) = h(k,l)
             ENDDO
           ENDDO
         CASE(2)
           j   = jsbcH(no,ide_t); 
           is  = isbcH(no,ide_t); 
           ie  = iebcH(no,ide_t);
           icl = siptNBI(no,ide_t) -1
           DO i = is, ie; 
             l = ij2l(i,j) 
             DO k = k1, kmz(l)
               icl = icl + 1      
               vhNB(k,i) = vhNGB (icl,no)*    weight + &
                           vhNGBp(icl,no)*(1.-weight)
               hvNB(k,i) = h(k,l)
             ENDDO
             
           ENDDO
           
         CASE(4)
           j  = jsbcH(no,ide_t); 
           is = isbcH(no,ide_t); 
           ie = iebcH(no,ide_t);
           icl = siptNBI(no,ide_t) -1
           DO i = is, ie; 
             l = ij2l(i,j) 
             DO k = k1, kmz(l)
               icl = icl + 1      
               vhSB(k,i) = vhNGB (icl,no)*    weight + &
                           vhNGBp(icl,no)*(1.-weight)
               hvSB(k,i) = h(k,l)
             ENDDO
           ENDDO
         END SELECT   
 
      END SELECT
!      print *,"uhWB:",sum(uhWB(:,:)),"huWB:",sum(huWB(:,:))
   END DO

END SUBROUTINE openbcUVH



!************************************************************************
SUBROUTINE readbcNGB(thrs)
!************************************************************************
!
!  Purpose: To read in a new frame with nested grid boundary conditions
!           if needed
!
!------------------------------------------------------------------------

   REAL, INTENT(IN) :: thrs

   REAL, ALLOCATABLE, DIMENSION (:,:) :: inputvar
   INTEGER :: icl, m1, m2, niNGB, nbiid, istat, nn,no,ide_t
   REAL :: auxthrs

   IF (nopen <= 0 .OR. ioNBTOGGLE <= 0) RETURN

   ! ... Initialize counter for nested grid boundaries
   niNGB = 0; 
   ide_t = omp_get_thread_num()+1
   ! ... Loop over open boundaries
   DO nn = 1, nopenH(ide_t) 
     no = noh2no(nn,ide_t)
     ! ... Cycle if not an embedded boundary
     IF ( itype(no) < 4 ) CYCLE

     ! ... Update counter for nested grid boundaries ...    
     niNGB = niNGB + 1; 

     ! Save and read new frame if thrs > thrsNGB .......
     
     IF (thrs > thrsNGB(no)) THEN

       ! ... Save variables from previous time .........
       thrsNGBp(no) = thrsNGB(no); 
       uhNGBp(siptNBIH(no,ide_t):eiptNBIH(no,ide_t),no) = uhNGB(siptNBIH(no,ide_t):eiptNBIH(no,ide_t),no); 
       vhNGBp(siptNBIH(no,ide_t):eiptNBIH(no,ide_t),no) = vhNGB(siptNBIH(no,ide_t):eiptNBIH(no,ide_t),no); 
       scNGBp(siptNBIH(no,ide_t):eiptNBIH(no,ide_t),no) = scNGB(siptNBIH(no,ide_t):eiptNBIH(no,ide_t),no);
       IF (ntr > 0) THEN 
         trNGBp(siptNBIH(no,ide_t):eiptNBIH(no,ide_t),:,no) = trNGB(siptNBIH(no,ide_t):eiptNBIH(no,ide_t),:,no);
       ENDIF
!       uhNGBp = uhNGB; 
!       vhNGBp = vhNGB; 
!       scNGBp = scNGB; 
!       trNGBp = trNGB;

       ! ... Set file ID ...............................            
       nbiid = nbiid0 + no

       ! ... Allocate space for temporary input variable array 
       ALLOCATE( inputvar ( iptNBI(no), 5+ntr ), STAT=istat )
       IF (istat /= 0) CALL allocate_error (istat,34)

       ! ... Read variables for NEXT FRAME variables ...
       
       !$omp critical
       print *,"*********************************"
       print *,"nthrsNGB:",thrsNGB(no),"hebra:",ide_t,"thrs:",thrs
       READ(nbiid) auxthrs, & 
           ((inputvar(m1,m2),m2=4,5+ntr),m1=1,iptNBI(no))
           contNG(no)=contNG(no)+1
           
           if(contNG(no) .EQ.  nopth(no))THEN
               contNG(no) = 0
               thrsNGB(no) = auxthrs
           else
               BACKSPACE nbiid
               
           end if
       print *,"despues:",thrsNGB(no),"input1:",sum(inputvar(:,4)),"input2:",sum(inputvar(:,5))   
       print *,"siptNBI:",siptNBI(no,ide_t),"eiptNBI:",eiptNBI(no,ide_t)
       print *,"contNG:",contNG(no),"nopth:",nopth(no)
       print *,"-----------------------------------------" 
       !$omp end critical
  
       ! ... Assign variables  
       SELECT CASE (iside(no))
       CASE(1,3)
         DO icl = siptNBI(no,ide_t), eiptNBI(no,ide_t)
           uhNGB(icl,no) = inputvar(icl,4)
           scNGB(icl,no) = inputvar(icl,5) 
           IF (ntr > 0) THEN
             trNGB(icl,1:ntr,niNGB) = inputvar(icl,6:5+ntr) 
           ENDIF
         ENDDO
       CASE(2,4)
         DO icl = siptNBI(no,ide_t), eiptNBI(no,ide_t)
           vhNGB(icl,no) = inputvar(icl,4)
           scNGB(icl,no) = inputvar(icl,5) 
           IF (ntr > 0) THEN
             trNGB(icl,1:ntr,no) = inputvar(icl,6:5+ntr) 
           ENDIF
         ENDDO
       END SELECT
          
       DEALLOCATE (inputvar)
     print *,"nthrsNGB:",thrsNGB(no),"hebra:",ide_t,"thrs:",thrs
     print *,"uhNGB:",sum(uhNGB(:,:)),"trNGB:",sum(trNGB(:,:,:)),"scNGB:",sum(scNGB(:,:))
     print *,"iside:",iside(no),"vhNGB:",sum(vhNGB(:,:))
     ENDIF
     
   ENDDO

END SUBROUTINE readbcNGB

!************************************************************************
SUBROUTINE MODqqddrr4openBC(Bstart,Bend,Bqq,Bsx,Bsy)
!************************************************************************
!
!  Purpose: To adjust the matrix coefficients used in the soln of the
!           continuity equation to account for open boundary conditions,
!           either wse or flow bdries.
!
!------------------------------------------------------------------------

   INTEGER, INTENT(IN) :: Bstart,Bend
   REAL, DIMENSION (Bstart:Bend+1), INTENT(INOUT) :: Bqq,Bsx,Bsy
   

   !.....Local variables.....
   INTEGER :: i, j, nn, is, ie, js, je, ks, ke,laux,iaux,jaux
   REAL    :: dt1, dtdx1, dtdy1,no,ide_t

   !.....Constants.....
   dtdx1 = dtdx*tz; 
   dtdy1 = dtdy*tz;
!   print *,"dtdx1:",dtdx1,"dtdy1:",dtdy1
   ide_t = omp_get_thread_num()+1
   !.....Loop over open boundaries.....
   DO nn = 1, nopenH(ide_t)
      no = noh2no(nn,ide_t)
      SELECT CASE ( itype(no) )
 
      !.....Case 1 -- wse specified.....
      CASE (1)
         
         ! Identify wse boundary as on the west, north, east, or south
         SELECT CASE ( iside(no) )
         ! West boundary
         CASE (1)
            i = isbcHH(no,ide_t); js = jsbcH(no,ide_t); je = jebcH(no,ide_t)
            ! Adjust [qq] array for column of nodes just inside boundary
            DO jaux=js,je
            	Bqq(ij2l(i+1,jaux)) = Bqq(ij2l(i+1,jaux)) + Bsx(ij2l(i,jaux))*s(ij2l(i,jaux))
            END DO
            ! Set sx(i,js:je) to zero in matrix (not really necessary)
            DO jaux=js,je
            	Bsx(ij2l(i,jaux)) = 0.0
            END DO

         ! North boundary
         CASE (2)
            j = jsbcH(no,ide_t); is = isbcHH(no,ide_t); ie = iebcHH(no,ide_t)
            ! Adjust [qq] array for row of nodes just inside boundary
            DO iaux=is,ie
            	Bqq(ij2l(iaux,j-1)) = Bqq(ij2l(iaux,j-1)) + Bsy(ij2l(iaux,j-1))*s(ij2l(iaux,j))
            END DO
            ! Set sy(is:ie,j-1) to zero in matrix (necessary)
            DO iaux=is,ie
            	Bsy(ij2l(iaux,j-1)) = 0.0
            END DO

         ! East boundary
         CASE (3)
            i = isbcHH(no,ide_t); js = jsbcH(no,ide_t); je = jebcH(no,ide_t)
            DO jaux=js,je
            	Bqq(ij2l(i-1,jaux)) = Bqq(ij2l(i-1,jaux)) + Bsx(ij2l(i-1,jaux))*s(ij2l(i,jaux))
            END DO
            ! Set sx(i,js:je) to zero in matrix (not really necessary)
            DO jaux=js,je
            	Bsx(ij2l(i-1,jaux)) = 0.0
            END DO

         ! South boundary
         CASE (4)
            j = jsbcH(no,ide_t); is = isbcHH(no,ide_t); ie = iebcHH(no,ide_t)
            ! Adjust [qq] array for row of nodes just inside boundary
            DO iaux=is,ie
            	Bqq(ij2l(iaux,j+1)) = Bqq(ij2l(iaux,j+1)) + Bsy(ij2l(iaux,j))*s(ij2l(iaux,j))
            END DO
            ! Set sy(is:ie,j-1) to zero in matrix (necessary)
            DO iaux=is,ie
            	Bsy(ij2l(iaux,j)) = 0.0
            END DO
            
         END SELECT

      !.....Case 2,4  -- Free surface flow specified.....
      CASE (2,4)

         ! Identify flow boundary as on the west, north, east, or south
         SELECT CASE ( iside(no) )

         ! West boundary
         CASE (1)
            i = isbcHH(no,ide_t); js = jsbcH(no,ide_t); je = jebcH(no,ide_t)
            ! Adjust [qq] array for flow rate into the
            ! column of nodes just inside the boundary 
            DO j = js, je
               laux = ij2l(i,j)
               Bqq(laux) = Bqq(laux) + dtdx1*SUM(uhWB  (k1:kmz(laux),j)) &
                                     + dtdx1*SUM(uhWBpp(k1:kmz(laux),j))
            END DO

         ! North boundary
         CASE (2)
            j = jsbcH(no,ide_t); is = isbcHH(no,ide_t); ie = iebcHH(no,ide_t)
            ! Adjust [qq] array for flow rate into the
            ! row of nodes just inside the boundary
            DO i = is, ie
               laux = ij2l(i,j)
               Bqq(laux) = Bqq(laux) - dtdy1*SUM(vhNB  (k1:kmz(laux),i)) &
                                     - dtdy1*SUM(vhNBpp(k1:kmz(laux),i))
            END DO

         ! East boundary
         CASE (3)
            i = isbcHH(no,ide_t); js = jsbcH(no,ide_t); je = jebcH(no,ide_t)
            ! Adjust [qq] array for flow rate into the
            ! column of nodes just inside the boundary
            DO j = js, je
               laux = ij2l(i,j)
               Bqq(laux) = Bqq(laux) - dtdx1*SUM(uhEB  (k1:kmz(laux),j)) &
                                     - dtdx1*SUM(uhEBpp(k1:kmz(laux),j)) 
            END DO

         ! South boundary
         CASE (4)
            j = jsbcH(no,ide_t); is = isbcHH(no,ide_t); ie = iebcHH(no,ide_t)          
            ! Adjust [qq] array for flow rate into the
            ! row of nodes just inside the boundary
            DO i = is, ie
               laux = ij2l(i,j)
               Bqq(laux) = Bqq(laux) + dtdy1*SUM(vhSB  (k1:kmz(laux),i)) &
                                     + dtdy1*SUM(vhSBpp(k1:kmz(laux),i))
            END DO
         END SELECT

      !.....Case 3 -- Submerged flow specified.....
      CASE (3)

         ! Identify flow boundary as on the west, north, east, or south
         SELECT CASE ( iside(no) )

         ! West boundary
         CASE (1)
            i  = isbc(no); 
            j  = jsbc(no);
            ks = iebc(no);
            ke = jebc(no);
            ! Adjust [qq] array for flow rate into the
            ! column of nodes just inside the boundary 
            Bqq(laux) = Bqq(laux) + dtdx1*SUM(uhWB  (ks:ke,j)) &
                                  + dtdx1*SUM(uhWBpp(ks:ke,j))

         ! North boundary
         CASE (2)
            j  = jsbc(no); 
            i  = isbc(no); 
            ks = iebc(no);
            ke = jebc(no);
            ! Adjust [qq] array for flow rate into the
            ! row of nodes just inside the boundary
            Bqq(laux) = Bqq(laux) - dtdy1*SUM(vhNB  (ks:ke,i)) &
                                  - dtdy1*SUM(vhNBpp(ks:ke,i))

         ! East boundary
         CASE (3)
            i  = isbc(no); 
            j  = jsbc(no);
            ks = iebc(no);
            ke = jebc(no);
            ! Adjust [qq] array for flow rate into the
            ! column of nodes just inside the boundary
            Bqq(laux) = Bqq(laux) - dtdx1*SUM(uhEB  (ks:ke,j)) &
                                  - dtdx1*SUM(uhEBpp(ks:ke,j)) 

         ! South boundary
         CASE (4)
            j  = jsbc(no); 
            i  = isbc(no); 
            ks = iebc(no);
            ke = jebc(no);
            ! Adjust [qq] array for flow rate into the
            ! row of nodes just inside the boundary
            Bqq(laux) = Bqq(laux) + dtdy1*SUM(vhSB  (ks:ke,i)) &
                                  + dtdy1*SUM(vhSBpp(ks:ke,i))
         END SELECT

      END SELECT
   END DO

END SUBROUTINE MODqqddrr4openBC

!************************************************************************
SUBROUTINE MODcoef4openBC 
!************************************************************************
!
!  Purpose: To adjust the matrix coefficients used in the soln of the
!           continuity equation to account for open boundary conditions.
!
!------------------------------------------------------------------------
 
   !.....Local variables.....
   INTEGER :: i, j, nn, is, ie, js, je, m,no,ide_t


   ide_t = omp_get_thread_num()+1
   ! ... Return if no open boundaries are specified
   

   !.....Loop over open boundaries.....
   DO nn = 1, nopenH(ide_t)
     no = noh2no(nn,ide_t)
      SELECT CASE ( itype(no) )
 
      !.....Case 1 -- wse specified.....
      CASE (1)

         ! Identify wse boundary as on the west, north, east, or south
         SELECT CASE ( iside(no) )

         ! West boundary
         CASE (1)
             i  = isbcH(no,ide_t); 
             js = jsbcH(no,ide_t); 
             je = jebcH(no,ide_t);
            DO j = js, je
!              m = (i-ifirst)*ibdwd + (j-jfirst)+1
              m  = ij2l(i,j)
              coeffA(m,1) = 1.E1
              coeffA(m,2) = 0.
              coeffA(m,3) = 0.
              coeffA(m,4) = 0.
              coeffA(m,5) = 0.
              rhs(m)    = 1.E1 * s(m) 
              zeta(m)   = s(m)
            ENDDO           

         ! North boundary
         CASE (2)
            j  = jsbcH(no,ide_t); 
            is = isbcHH(no,ide_t); 
            ie = iebcHH(no,ide_t);
            DO i = is, ie
!              m = (i-ifirst)*ibdwd + (j-jfirst)+1
              m  = ij2l(i,j)
              coeffA(m,1) = 1.E1
              coeffA(m,2) = 0.
              coeffA(m,3) = 0.
              coeffA(m,4) = 0.
              coeffA(m,5) = 0.
              rhs(m)    = 1.E1 * s(m) 
              zeta(m)   = s(m)
            ENDDO        
   
         ! East boundary
         CASE (3)
            i  = isbcH(no,ide_t); 
            js = jsbcH(no,ide_t); 
            je = jebcH(no,ide_t);
            DO j = js, je
!              m = (i-ifirst)*ibdwd + (j-jfirst)+1
              m  = ij2l(i,j)
              coeffA(m,1) = 1.E1
              coeffA(m,2) = 0.
              coeffA(m,3) = 0.
              coeffA(m,4) = 0.
              coeffA(m,5) = 0.
              rhs(m)    = 1.E1 * s(m)
              zeta(m)   = s(m) 
            ENDDO           
 
         ! South boundary
         CASE (4)
            j  = jsbcH(no,ide_t); 
            is = isbcHH(no,ide_t); 
            ie = iebcHH(no,ide_t);
            DO i = is, ie
!              m = (i-ifirst)*ibdwd + (j-jfirst)+1
              m  = ij2l(i,j)
              coeffA(m,1) = 1.E1
              coeffA(m,2) = 0.
              coeffA(m,3) = 0.
              coeffA(m,4) = 0.
              coeffA(m,5) = 0.
              rhs(m)    = 1.E1*s(m) 
              zeta(m)   = s(m)
            ENDDO           

         END SELECT

      !.....Case 2&3 -- Free surface & submerged flow specified.....
      CASE (2:)

         RETURN

      END SELECT
   END DO

END SUBROUTINE MODcoef4openBC

!************************************************************************
SUBROUTINE MODexmom4openBCX 
!************************************************************************
!
!  Purpose: To recompute ex matrix for velocity columns located along or 
!           next to open boundaries 
!
!------------------------------------------------------------------------

   !.....Local variables.....
   REAL :: twodt1
   REAL :: advx, advy, uE, uW, vN, vS, wU, wD, &
           ubdry, vbdry, uhbdry, vhbdry
   INTEGER :: i, j, k, l, istat, kmx, kmy, k1x, k1y, k1ne
   INTEGER :: nn, is, ie, js, je, ks, ke, nwlayers,no,ide_t

   !.....Constant.....
   twodt1 = twodt*tz
   ide_t = omp_get_thread_num()+1
   !.....Loop over open boundaries.....
   DO nn = 1, nopenHH(ide_t)
     no = noh2noH(nn,ide_t)
     SELECT CASE ( itype(no) )
 
     !....... Case 1 -- wse specified .....................................
     CASE (1)

       ! Identify wse boundary as on the west, north, east, or south
       SELECT CASE ( iside(no) )

       ! ..... West boundary ......
       CASE (1)

         i  = isbcH(no,ide_t); 
         js = jsbcH(no,ide_t); 
         je = jebcH(no,ide_t);
         
         DO j = js, je

           ! ... Map (i,j) into l-index
           l = ij2l(i,j);
           
           ! Compute the layer number for the bottom wet u-pt
           kmx = MIN(kmz(lEC(l)), kmz(l))
           k1x =                 k1u(l)

           ! Compute number of wet layers & skip for dry cells
           nwlayers = (kmx-k1x) + 1
           IF(nwlayers < 1) CYCLE

           ! Compute explicit term (only using advection)
           ex(:,l) = 0.0
           DO k = k1x,kmx
                                                                   
            ! Horizontal advection - Upwind differencing  
            uE = uhp(k,    lEC(l) ) +uhp(k  ,    l )
            uW = uhp(k,        l  ) +uhp(k  ,    l )
            vN = vhp(k,    lEC(l) ) +vhp(k  ,    l )
            vS = vhp(k,lSC(lEC(l))) +vhp(k  ,lSC(l))
            wU = wp (k,    lEC(l) ) +wp (k  ,    l )
            wD = wp (k+1,  lEC(l) ) +wp (k+1,    l )        
            advx = ( (uE+ABS(uE))* upp(k,    l ) +           &
                     (uE-ABS(uE))* upp(k,lEC(l)) -           &
                     (uW+ABS(uW))* upp(k,    l ) -           &
                     (uW-ABS(uW))* upp(k,    l ) ) / fourdx  &
                  +( (vN+ABS(vN))* upp(k,    l ) +           &
                     (vN-ABS(vN))* upp(k,lNC(l)) -           &
                     (vS+ABS(vS))* upp(k,lSC(l)) -           &
                     (vS-ABS(vS))* upp(k,    l ) ) / fourdy  

            ! Vertical advection - Upwind for near bdry. cells
            advx=advx+((wU+ABS(wU)) * upp(k  ,l) +           &
                       (wU-ABS(wU)) * upp(k-1,l)) / 4.       &
                     -((wD+ABS(wD)) * upp(k+1,l) +           &
                       (wD-ABS(wD)) * upp(k  ,l)) / 4.
			      		                                  
            !.....Final explicit term.....
            ex(k,l) = uhpp(k,l) - twodt1*(advx*iadv)
            
           END DO
           
         END DO

       ! ..... North boundary ......
       CASE (2)

!         print *,"no hago nada"
!         RETURN

       ! ..... East boundary ......
       CASE (3)
         i  = isbcH(no,ide_t)-1; 
         js = jsbcH(no,ide_t)  ; 
         je = jebcH(no,ide_t)  ;
         DO j = js, je

           ! ... Map (i,j) into l-index
           l = ij2l(i,j);

           ! Compute the layer number for the bottom wet u-pt
           kmx = MIN(kmz(lEC(l)), kmz(l))
           k1x =                 k1u(l)

           ! Compute number of wet layers & skip for dry cells
           nwlayers = (kmx-k1x) + 1
           IF(nwlayers < 1) CYCLE

           ! Compute explicit term (only using advection)
           ex(:,l) = 0.0
           DO k = k1x,kmx
                                                                   
            ! Horizontal advection - Upwind differencing  
            uE = uhp(k,        l  ) +uhp(k  ,    l )
            uW = uhp(k,        l  ) +uhp(k  ,lWC(l)) 
            vN = vhp(k,    lEC(l) ) +vhp(k  ,    l )
            vS = vhp(k,lSC(lEC(l))) +vhp(k  ,lSC(l))
            wU = wp (k,    lEC(l) ) +wp (k  ,    l )
            wD = wp (k+1,  lEC(l) ) +wp (k+1,    l )        
            advx = ( (uE+ABS(uE))* upp(k,    l ) +           &
                     (uE-ABS(uE))* upp(k,    l ) -           &
                     (uW+ABS(uW))* upp(k,lWC(l)) -           &
                     (uW-ABS(uW))* upp(k,    l ) ) / fourdx  &
                  +( (vN+ABS(vN))* upp(k,    l ) +           &
                     (vN-ABS(vN))* upp(k,lNC(l)) -           &
                     (vS+ABS(vS))* upp(k,lSC(l)) -           &
                     (vS-ABS(vS))* upp(k,    l ) ) / fourdy  

            ! Vertical advection - Upwind for near bdry. cells
            advx=advx+((wU+ABS(wU)) * upp(k  ,l) +           &
                       (wU-ABS(wU)) * upp(k-1,l)) / 4.       &
                     -((wD+ABS(wD)) * upp(k+1,l) +           &
                       (wD-ABS(wD)) * upp(k  ,l)) / 4.
			      		                                  
            !.....Final explicit term.....             Cola Beznar
            !ex(k,l) = uhpp(k,l) - twodt1*(advx*iadv) 
  
           END DO
         ENDDO           

       ! ..... South boundary ......
       CASE (4)
!         print *,"no hago nada"
!         RETURN

       END SELECT

     !....... Case 2,4 -- Free surface flow specified ......................
     CASE (2,4)

      ! Identify wse boundary as on the west, north, east, or south
       SELECT CASE ( iside(no) )

       ! ..... West boundary ......
       CASE (1)

         i  = isbcH(no,ide_t); 
         js = jsbcH(no,ide_t); 
         je = jebcH(no,ide_t)
         DO j = js, je

           ! ... Map (i,j) into l-index
           l = ij2l(i,j); 

           ! ... Cycle if W-column is dry
           IF (.NOT. mask2d(i+1,j)) CYCLE

           ! Compute the layer number for the bottom wet u-pt
           kmx = MIN(kmz(lEC(l)), kmz(l))
           k1x =                 k1u(l)

           ! Compute number of wet layers & skip for dry cells
           nwlayers = (kmx-k1x) + 1
           IF(nwlayers < 1) CYCLE

           ! Compute explicit term (only using advection)
           ex(:,l) = 0.0
           DO k = k1x,kmx
                                                                   
            ! Horizontal advection - Upwind differencing  
            uhbdry = (uhWB(k,j)+uhWBpp(k,j))/2.
            IF (huWBpp(k,j)>ZERO) THEN
               ubdry  = uhWBpp(k,j)/huWBpp(k,j)
            ELSE
               ubdry  = 0.0
            ENDIF
            uE = uhp(k,    lEC(l) ) +uhp(k  ,    l )
            uW = uhp(k,        l  ) +uhbdry
            vN = vhp(k,    lEC(l) ) +vhp(k  ,    l )
            vS = vhp(k,lSC(lEC(l))) +vhp(k  ,lSC(l))
            wU = wp (k,    lEC(l) ) +wp (k  ,    l )
            wD = wp (k+1,  lEC(l) ) +wp (k+1,    l )        
            advx = ( (uE+ABS(uE))* upp(k,    l ) +           &
                     (uE-ABS(uE))* upp(k,lEC(l)) -           &
                     (uW+ABS(uW))* ubdry         -           &
                     (uW-ABS(uW))* upp(k,    l ) ) / fourdx  &
                  +( (vN+ABS(vN))* upp(k,    l ) +           &
                     (vN-ABS(vN))* upp(k,lNC(l)) -           &
                     (vS+ABS(vS))* upp(k,lSC(l)) -           &
                     (vS-ABS(vS))* upp(k,    l ) ) / fourdy  

            ! Vertical advection - Upwind for near bdry. cells
            advx=advx+((wU+ABS(wU)) * upp(k  ,l) +           &
                       (wU-ABS(wU)) * upp(k-1,l)) / 4.       &
                     -((wD+ABS(wD)) * upp(k+1,l) +           &
                       (wD-ABS(wD)) * upp(k  ,l)) / 4.
			      		                                  
            !.....Final explicit term.....
            ex(k,l) = uhpp(k,l) - twodt1*(advx*iadv)
                
           END DO
         END DO

       ! ..... North boundary ......
       CASE (2)

!         RETURN

       ! ..... East boundary ......
       CASE (3)

         i  = isbcH(no,ide_t)-1; 
         js = jsbcH(no,ide_t)  ; 
         je = jebcH(no,ide_t)  ; 
         DO j = js, je

           ! ... Map (i,j) into l-index
           l = ij2l(i,j); 

           ! ... Cycle if W-column is dry
           IF (.NOT. mask2d(i+1,j)) CYCLE

           ! Compute the layer number for the bottom wet u-pt
           kmx = MIN(kmz(lEC(l)), kmz(l))
           k1x =                 k1u(l)

           ! Compute number of wet layers & skip for dry cells
           nwlayers = (kmx-k1x) + 1
           IF(nwlayers < 1) CYCLE

           ! Compute explicit term (only using advection)
           ex(:,l) = 0.0
           DO k = k1x,kmx
                                                                   
            ! Horizontal advection - Upwind differencing  
            uhbdry = (uhEB(k,j)+uhEBpp(k,j))/2.
            IF (huEBpp(k,j)>ZERO) THEN
               ubdry  = uhEBpp(k,j)/huEBpp(k,j)
            ELSE
               ubdry  = 0.0
            ENDIF
            uE = uhp(k,        l  ) +uhbdry
            uW = uhp(k,    lWC(l) ) +uhp(k  ,    l ) 
            vN = vhp(k,    lEC(l) ) +vhp(k  ,    l )
            vS = vhp(k,lSC(lEC(l))) +vhp(k  ,lSC(l))
            wU = wp (k,    lEC(l) ) +wp (k  ,    l )
            wD = wp (k+1,  lEC(l) ) +wp (k+1,    l )        
            advx = ( (uE+ABS(uE))* upp(k,    l ) +           &
                     (uE-ABS(uE))* ubdry         -           &
                     (uW+ABS(uW))* upp(k,lWC(l)) -           &
                     (uW-ABS(uW))* upp(k,    l ) ) / fourdx  &
                  +( (vN+ABS(vN))* upp(k,    l ) +           &
                     (vN-ABS(vN))* upp(k,lNC(l)) -           &
                     (vS+ABS(vS))* upp(k,lSC(l)) -           &
                     (vS-ABS(vS))* upp(k,    l ) ) / fourdy  

            ! Vertical advection - Upwind for near bdry. cells
            advx=advx+((wU+ABS(wU)) * upp(k  ,l) +           &
                       (wU-ABS(wU)) * upp(k-1,l)) / 4.       &
                     -((wD+ABS(wD)) * upp(k+1,l) +           &
                       (wD-ABS(wD)) * upp(k  ,l)) / 4.
		      		                                  
            !.....Final explicit term.....
            ex(k,l) = uhpp(k,l) - twodt1*(advx*iadv) 
  
           END DO
         ENDDO           

       ! ..... South boundary ......
       CASE (4)

!         RETURN

       END SELECT

     !....... Case 3 -- Submerged Flow specified ..........................
     CASE (3)

      ! Identify wse boundary as on the west, north, east, or south
       SELECT CASE ( iside(no) )

       ! ..... West boundary ......
       CASE (1)

         i  = isbc(no); 
         j  = jsbc(no); 
         ks = iebc(no);
         ke = jebc(no);
         l = ij2l(i,j);

         ! Compute explicit term (use only advection)
         ex(ks:ke,l) = 0.0

         DO k = ks,ke
                                                                   
            ! Horizontal advection - Upwind differencing  
            uhbdry = (uhWBpp(k,j)+uhWB  (k,j))/2.
            ubdry  =  uhWBpp(k,j)/huWBpp(k,j)
            uE = uhp(k,    lEC(l) ) +uhp(k  ,    l )
            uW = uhp(k,        l  ) +uhbdry
            vN = vhp(k,    lEC(l) ) +vhp(k  ,    l )
            vS = vhp(k,lSC(lEC(l))) +vhp(k  ,lSC(l))
            wU = wp (k,    lEC(l) ) +wp (k  ,    l )
            wD = wp (k+1,  lEC(l) ) +wp (k+1,    l )        
            advx = ( (uE+ABS(uE))* upp(k,    l ) +           &
                     (uE-ABS(uE))* upp(k,lEC(l)) -           &
                     (uW+ABS(uW))* ubdry         -           &
                     (uW-ABS(uW))* upp(k,    l ) ) / fourdx  &
                  +( (vN+ABS(vN))* upp(k,    l ) +           &
                     (vN-ABS(vN))* upp(k,lNC(l)) -           &
                     (vS+ABS(vS))* upp(k,lSC(l)) -           &
                     (vS-ABS(vS))* upp(k,    l ) ) / fourdy  

            ! Vertical advection - Upwind for near bdry. cells
            advx=advx+((wU+ABS(wU)) * upp(k  ,l) +           &
                       (wU-ABS(wU)) * upp(k-1,l)) / 4.       &
                     -((wD+ABS(wD)) * upp(k+1,l) +           &
                       (wD-ABS(wD)) * upp(k  ,l)) / 4.
			      		                                  
            !.....Final explicit term.....
            ex(k,l) = uhpp(k,l) - twodt1*(advx*iadv)
                
         END DO

       ! ..... North boundary ......
       CASE (2)

         i  = isbc(no); 
         j  = jsbc(no); 
         ks = iebc(no);
         ke = jebc(no);

         ! Compute explicit term (use only advection) for u-points
		 ! to the East of a water column with North flow specified
         IF (mask2d(i+1,j)) THEN 

           ! ... Map (i,j) into l-index
           l = ij2l(i,j);
           ex(ks:ke,l) = 0.0
           DO k = ks,ke
                                                                   
             ! Horizontal advection - Upwind differencing  
             vhbdry = (vhNBpp(k,i)+vhNB  (k,i))/2.
             ! vbdry  =  vhNBpp(k,i)/hvNBpp(k,i)
             uE = uhp(k,    lEC(l) ) +uhp(k  ,    l )
             uW = uhp(k,    lWC(l) ) +uhp(k  ,    l )
             vN = vhp(k,    lEC(l) ) +vhbdry
             vS = vhp(k,lSC(lEC(l))) +vhp(k  ,lSC(l))
             wU = wp (k,    lEC(l) ) +wp (k  ,    l )
             wD = wp (k+1,  lEC(l) ) +wp (k+1,    l )        
             advx = ( (uE+ABS(uE))* upp(k,    l ) +           &
                      (uE-ABS(uE))* upp(k,lEC(l)) -           &
                      (uW+ABS(uW))* upp(k,lWC(l)) -           & ! 30062008
                      (uW-ABS(uW))* upp(k,    l ) ) / fourdx  &
                   +( (vN+ABS(vN))* upp(k,    l ) +           &
                      (vN-ABS(vN))* upp(k,lNC(l)) -           &
                      (vS+ABS(vS))* upp(k,lSC(l)) -           &
                      (vS-ABS(vS))* upp(k,    l ) ) / fourdy  

             ! Vertical advection - Upwind for near bdry. cells
             advx=advx+((wU+ABS(wU)) * upp(k  ,l) +           &
                        (wU-ABS(wU)) * upp(k-1,l)) / 4.       &
                      -((wD+ABS(wD)) * upp(k+1,l) +           &
                        (wD-ABS(wD)) * upp(k  ,l)) / 4.
			      		                                  
             !.....Final explicit term.....
             ex(k,l) = uhpp(k,l) - twodt1*(advx*iadv)
                
           END DO 

         ENDIF

         ! Compute explicit term (use only advection) for u-points
         ! to the West of a water column with North flow specified
         IF (mask2d(i-1,j)) THEN 

           l = ij2l(i-1,j);

           ! Compute explicit term (use only advection)
           ex(ks:ke,l) = 0.0
           DO k = ks,ke
                                                                   
             ! Horizontal advection - Upwind differencing  
             vhbdry = (vhNBpp(k,i)+vhNB  (k,i))/2.
             ! vbdry  =  vhNBpp(k,i)/hvNBpp(k,i)
             uE = uhp(k,    lEC(l) ) +uhp(k  ,    l )
             uW = uhp(k,    lWC(l) ) +uhp(k  ,    l )
             vN = vhbdry             +vhp(k  ,    l )
             vS = vhp(k,lSC(lEC(l))) +vhp(k  ,lSC(l))
             wU = wp (k,    lEC(l) ) +wp (k  ,    l )
             wD = wp (k+1,  lEC(l) ) +wp (k+1,    l )        
             advx = ( (uE+ABS(uE))* upp(k,    l ) +           &
                      (uE-ABS(uE))* upp(k,lEC(l)) -           &
                      (uW+ABS(uW))* upp(k,lWC(l)) -           &  ! 30062008
                      (uW-ABS(uW))* upp(k,    l ) ) / fourdx  &
                   +( (vN+ABS(vN))* upp(k,    l ) +           &
                      (vN-ABS(vN))* upp(k,lNC(l)) -           &
                      (vS+ABS(vS))* upp(k,lSC(l)) -           &
                      (vS-ABS(vS))* upp(k,    l ) ) / fourdy  

             ! Vertical advection - Upwind for near bdry. cells
             advx=advx+((wU+ABS(wU)) * upp(k  ,l) +           &
                        (wU-ABS(wU)) * upp(k-1,l)) / 4.       &
                      -((wD+ABS(wD)) * upp(k+1,l) +           &
                        (wD-ABS(wD)) * upp(k  ,l)) / 4.
			      		                                  
             !.....Final explicit term.....
             ex(k,l) = uhpp(k,l) - twodt1*(advx*iadv)
                
           END DO

         ENDIF

       ! ..... East boundary ......
       CASE (3)

         i  = isbc(no)-1; 
         j  = jsbc(no)  ; 
         ks = iebc(no)  ;
         ke = jebc(no)  ;

         ! ... Map (i,j) into l-index
         l = ij2l(i,j);

         ! Compute explicit term (only using advection)
         ex(ks:ke,l) = 0.0
         DO k = ks,ke
                                                                   
           ! Horizontal advection - Upwind differencing  
           uhbdry = (uhEB(k,j)+uhEBpp(k,j))/2.
           ubdry  = uhEBpp(k,j)/huEBpp(k,j)
           uE = uhp(k,        l  ) +uhbdry
           uW = uhp(k,        l  ) +uhp(k  ,    l )
           vN = vhp(k,    lEC(l) ) +vhp(k  ,    l )
           vS = vhp(k,lSC(lEC(l))) +vhp(k  ,lSC(l))
           wU = wp (k,    lEC(l) ) +wp (k  ,    l )
           wD = wp (k+1,  lEC(l) ) +wp (k+1,    l )        
           advx = ( (uE+ABS(uE))* upp(k,    l ) +           &
                    (uE-ABS(uE))* ubdry         -           &
                    (uW+ABS(uW))* upp(k,lWC(l)) -           &
                    (uW-ABS(uW))* upp(k,    l ) ) / fourdx  &
                 +( (vN+ABS(vN))* upp(k,    l ) +           &
                    (vN-ABS(vN))* upp(k,lNC(l)) -           &
                    (vS+ABS(vS))* upp(k,lSC(l)) -           &
                    (vS-ABS(vS))* upp(k,    l ) ) / fourdy  
			      		                                  
           ! Vertical advection - Upwind for near bdry. cells
           advx=advx+((wU+ABS(wU)) * upp(k  ,l) +           &
                      (wU-ABS(wU)) * upp(k-1,l)) / 4.       &
                    -((wD+ABS(wD)) * upp(k+1,l) +           &
                      (wD-ABS(wD)) * upp(k  ,l)) / 4.

           !.....Final explicit term.....
           ex(k,l) = uhpp(k,l) - twodt1*(advx*iadv) 
  
         END DO

       ! ..... South boundary ......
       CASE (4)

         i  = isbc(no); 
         j  = jsbc(no); 
         ks = iebc(no);
         ke = jebc(no);

         ! Compute explicit term (use only advection) for u-points
		 ! to the East of a water column with North flow specified
         IF (mask2d(i+1,j)) THEN 

           ! ... Map (i,j) into l-index
           l = ij2l(i,j);

           ! Compute explicit term (use only advection)
           ex(ks:ke,l) = 0.0

           DO k = ks,ke
                                                                   
             ! Horizontal advection - Upwind differencing  
             vhbdry = (vhSBpp(k,i)+vhSB  (k,i))/2.
             ! vbdry  =  vhSBpp(k,i)/hvSBpp(k,i)
             uE = uhp(k,    lEC(l) ) +uhp(k  ,    l )
             uW = uhp(k,    lWC(l) ) +uhp(k  ,    l )
             vN = vhp(k,    lEC(l) ) +vhp(k  ,    l )
             vS = vhp(k,lSC(lEC(l))) +vhbdry
             wU = wp (k,    lEC(l) ) +wp (k  ,    l )
             wD = wp (k+1,  lEC(l) ) +wp (k+1,    l )        
             advx = ( (uE+ABS(uE))* upp(k,    l ) +           &
                      (uE-ABS(uE))* upp(k,lEC(l)) -           &
                      (uW+ABS(uW))* upp(k,lWC(l)) -           & !30062008
                      (uW-ABS(uW))* upp(k,    l ) ) / fourdx  &
                   +( (vN+ABS(vN))* upp(k,    l ) +           &
                      (vN-ABS(vN))* upp(k,lNC(l)) -           &
                      (vS+ABS(vS))* upp(k,lSC(l)) -           &
                      (vS-ABS(vS))* upp(k,    l ) ) / fourdy   

             ! Vertical advection - Upwind for near bdry. cells
             advx=advx+((wU+ABS(wU)) * upp(k  ,l) +           &
                        (wU-ABS(wU)) * upp(k-1,l)) / 4.       &
                      -((wD+ABS(wD)) * upp(k+1,l) +           &
                        (wD-ABS(wD)) * upp(k  ,l)) / 4.
			      		                                  
             !.....Final explicit term.....
             ex(k,l) = uhpp(k,l) - twodt1*(advx*iadv)
               
           END DO

         ENDIF


         ! Compute explicit term (use only advection) for u-points
         ! to the West of a water column with North flow specified
         IF (mask2d(i-1,j)) THEN 

           ! ... Map (i,j) into l-index
           l = ij2l(i-1,j);
 
           ! Compute explicit term (use only advection)
           ex(ks:ke,l) = 0.0

           DO k = ks,ke
                                                                   
             ! Horizontal advection - Upwind differencing  
             vhbdry = (vhSBpp(k,i)+vhSB  (k,i))/2.
             ! vbdry  =  vhSBpp(k,i)/hvSBpp(k,i)
             uE = uhp(k,    lEC(l) ) +uhp(k  ,    l )
             uW = uhp(k,    lWC(l) ) +uhp(k  ,    l )
             vN = vhp(k,    lEC(l) ) +vhp(k  ,    l )
             vS = vhbdry             +vhp(k  ,lSC(l))
             wU = wp (k,    lEC(l) ) +wp (k  ,    l )
             wD = wp (k+1,  lEC(l) ) +wp (k+1,    l )        
             advx = ( (uE+ABS(uE))* upp(k,    l ) +           &
                      (uE-ABS(uE))* upp(k,lEC(l)) -           &
                      (uW+ABS(uW))* upp(k,lWC(l)) -           &
                      (uW-ABS(uW))* upp(k,    l ) ) / fourdx  &
                   +( (vN+ABS(vN))* upp(k,    l ) +           &
                      (vN-ABS(vN))* upp(k,lNC(l)) -           &
                      (vS+ABS(vS))* upp(k,lSC(l)) -           &
                      (vS-ABS(vS))* upp(k,    l ) ) / fourdy   

             ! Vertical advection - Upwind for near bdry. cells
             advx=advx+((wU+ABS(wU)) * upp(k  ,l) +           &
                        (wU-ABS(wU)) * upp(k-1,l)) / 4.       &
                      -((wD+ABS(wD)) * upp(k+1,l) +           &
                        (wD-ABS(wD)) * upp(k  ,l)) / 4.
			      		                                  
             !.....Final explicit term.....
             ex(k,l) = uhpp(k,l) - twodt1*(advx*iadv)
                
           END DO

         ENDIF

       END SELECT ! End select from boundary locations (EWNS)

     END SELECT ! End select from boundary types (wse/flow/submerged)

   END DO ! End loop over open boundaries
!$omp barrier
END SUBROUTINE MODexmom4openBCX

!************************************************************************
SUBROUTINE MODexmom4openBCY 
!************************************************************************
!
!  Purpose: To recompute ex matrix for cells located next to cells with 
!           open boundary conditions specified 
!
!------------------------------------------------------------------------

   !.....Local variables.....
   REAL :: twodt1
   REAL :: advx, advy, uE, uW, vN, vS, wU, wD, &
           ubdry, vbdry, uhbdry, vhbdry
   INTEGER :: i, j, k, l, istat, kmx, kmy, k1x, k1y, k1ne
   INTEGER :: nn, is, ie, js, je, ks, ke, nwlayers,no,ide_t

   !.....Constant.....
   twodt1 = twodt*tz
   ide_t = omp_get_thread_num()+1
   !.....Loop over open boundaries.....
   DO nn = 1, nopenHH(ide_t)
     no = noh2noH(nn,ide_t)
     SELECT CASE ( itype(no) )
 
     !....... Case 1 -- wse specified......................................
     CASE (1)

       ! Identify wse boundary as on the west, north, east, or south
       SELECT CASE ( iside(no) )

       ! ..... West boundary ......
       CASE (1)

!         RETURN

       ! ..... North boundary ......
       CASE (2)
         
         j  = jsbcH(no,ide_t)-1; 
         is = isbcH(no,ide_t)  ; 
         ie = iebcH(no,ide_t)  ;
         DO i = is, ie

           ! ... Map (i,j) into l-index
           l = ij2l(i,j);

           ! Compute the layer number for top & bottom wet v-pt
           kmy = MIN(kmz(lNC(l)), kmz(l))
           k1y =                 k1v(l)

           ! Compute number of wet layers & skip for dry cells
           nwlayers = (kmy-k1y) + 1
           IF(nwlayers < 1) CYCLE

           ! Compute explicit term
           ex(:,l) = 0.0	     
           DO k = k1y,kmy
                                                                 
            ! Horizontal advection Upwind differencing 
            uE = uhp(k,    lNC(l) ) +uhp(k,    l )
            uW = uhp(k,lWC(lNC(l))) +uhp(k,lWC(l))
            vN = vhp(k,        l  ) +vhp(k,    l )
            vS = vhp(k,        l  ) +vhp(k,lSC(l))
            wU = wp (k,    lNC(l) ) +wp (k,    l )
            wD = wp (k+1,  lNC(l) ) +wp (k+1  ,l )        
            advy = ( (uE+ABS(uE))* vpp(k,    l ) +           &
                &    (uE-ABS(uE))* vpp(k,lEC(l)) -           &
                &    (uW+ABS(uW))* vpp(k,lWC(l)) -           &
                &    (uW-ABS(uW))* vpp(k,    l ) ) / fourdx  &
                & +( (vN+ABS(vN))* vpp(k,    l ) +           &
            	&    (vN-ABS(vN))* vpp(k,    l ) -           &
                &    (vS+ABS(vS))* vpp(k,lSC(l)) -           &
                &    (vS-ABS(vS))* vpp(k,    l ) ) / fourdy  

            ! Vertical advection - Upwind for near bdry. cells
            advy=advy+((wU+ABS(wU)) * vpp(k  ,l) +           &
                       (wU-ABS(wU)) * vpp(k-1,l)) / 4.       &
                     -((wD+ABS(wD)) * vpp(k+1,l) +           &
                       (wD-ABS(wD)) * vpp(k  ,l)) / 4.
                   
            !.....Final explicit term.....
            ex(k,l) = vhpp(k,l) - twodt1*(advy*iadv)
                
           END DO
         END DO

       ! ..... East boundary
       CASE (3)

!         RETURN

       ! ..... South boundary ......
       CASE (4)

         j  = jsbcH(no,ide_t); 
         is = isbcH(no,ide_t); 
         ie = iebcH(no,ide_t)
         DO i = is, ie

           ! ... Map (i,j) into l-index
           l = ij2l(i,j);

           ! Compute the layer number for top & bottom wet v-pt
           kmy = MIN(kmz(lNC(l)), kmz(l))
           k1y =                 k1v(l)

           ! Compute number of wet layers & skip for dry cells
           nwlayers = (kmy-k1y) + 1
           IF(nwlayers < 1) CYCLE

           ! Compute explicit term
           ex(:,l) = 0.0	     
           DO k = k1y,kmy
                                                                 
            ! Horizontal advection Upwind differencing 
            uE = uhp(k,    lNC(l) ) +uhp(k,    l )
            uW = uhp(k,lWC(lNC(l))) +uhp(k,lWC(l))
            vN = vhp(k,    lNC(l) ) +vhp(k,    l )
            vS = vhp(k,        l  ) +vhp(k,    l )
            wU = wp (k,    lNC(l) ) +wp (k,    l )
            wD = wp (k+1,  lNC(l) ) +wp (k+1  ,l )        
            advy = ( (uE+ABS(uE))* vpp(k,    l ) +           &
                &    (uE-ABS(uE))* vpp(k,lEC(l)) -           &
                &    (uW+ABS(uW))* vpp(k,lWC(l)) -           &
                &    (uW-ABS(uW))* vpp(k,    l ) ) / fourdx  &
                & +( (vN+ABS(vN))* vpp(k,    l ) +           &
            	&    (vN-ABS(vN))* vpp(k,lNC(l)) -           &
                &    (vS+ABS(vS))* vpp(k,    l ) -           &
                &    (vS-ABS(vS))* vpp(k,    l ) ) / fourdy  

            ! Vertical advection - Upwind in near bdry. cells
            advy=advy+((wU+ABS(wU)) * vpp(k  ,l) +           &
                       (wU-ABS(wU)) * vpp(k-1,l)) / 4.       &
                     -((wD+ABS(wD)) * vpp(k+1,l) +           &
                       (wD-ABS(wD)) * vpp(k  ,l)) / 4.
                   
            !.....Final explicit term.....
            ex(k,l) = vhpp(k,l) - twodt1*(advy*iadv)
                
           END DO
         END DO 

       END SELECT

     !....... Case 2,4 -- Free surface flow specified.......................
     CASE (2,4)

       ! Identify wse boundary as on the west, north, east, or south
       SELECT CASE ( iside(no) )

       ! ..... West boundary ...... *** need to include code here
       CASE (1)

       !         RETURN

       ! ..... North boundary ......
       CASE (2)
         
         j  = jsbcH(no,ide_t)-1; 
         is = isbcH(no,ide_t)  ; 
         ie = iebcH(no,ide_t)  ;
         !print *,"H:",ide_t,"j:",j
         !print *,"H:",ide_t,"is:",is,"ie:",ie
         DO i = is, ie

           ! ... Map (i,j) into l-index
           l = ij2l(i,j);

           ! ... Cycle if N-column is dry
           IF ( .NOT. mask2d(i,j+1) ) CYCLE

           ! Compute the layer number for top & bottom wet v-pt
           kmy = MIN(kmz(lNC(l)), kmz(l))
           k1y =                 k1v(l)

           ! Compute number of wet layers & skip for dry cells
           nwlayers = (kmy-k1y) + 1
           IF(nwlayers < 1) CYCLE

           ! Compute explicit term
           ex(:,l) = 0.0	     
           DO k = k1y,kmy
                                                                 
            ! Horizontal advection Upwind differencing 
            vhbdry = (vhNB(k,i)+vhNBpp(k,i))/2.
            IF (hvNBpp(k,i)>ZERO) THEN
               vbdry  = vhNBpp(k,i)/hvNBpp(k,i)
            ELSE
               vbdry  = 0.0
            ENDIF
            uE = uhp(k,    lNC(l) ) +uhp(k,    l )
            uW = uhp(k,lWC(lNC(l))) +uhp(k,lWC(l))
            vN = vhp(k,        l  ) +vhbdry
            vS = vhp(k,        l  ) +vhp(k,lSC(l))
            wU = wp (k,    lNC(l) ) +wp (k,    l )
            wD = wp (k+1,  lNC(l) ) +wp (k+1  ,l )        
            advy = ( (uE+ABS(uE))* vpp(k,    l ) +           &
                &    (uE-ABS(uE))* vpp(k,lEC(l)) -           &
                &    (uW+ABS(uW))* vpp(k,lWC(l)) -           &
                &    (uW-ABS(uW))* vpp(k,    l ) ) / fourdx  &
                & +( (vN+ABS(vN))* vpp(k,    l ) +           &
            	&    (vN-ABS(vN))* vbdry         -           &
                &    (vS+ABS(vS))* vpp(k,lSC(l)) -           &
                &    (vS-ABS(vS))* vpp(k,    l ) ) / fourdy  

            ! Vertical advection - Upwind for near bdry. cells
            advy=advy+((wU+ABS(wU)) * vpp(k  ,l) +           &
                       (wU-ABS(wU)) * vpp(k-1,l)) / 4.       &
                     -((wD+ABS(wD)) * vpp(k+1,l) +           &
                       (wD-ABS(wD)) * vpp(k  ,l)) / 4.
                   
            !.....Final explicit term.....
            ex(k,l) = vhpp(k,l) - twodt1*(advy*iadv)
                
           END DO
         END DO

       ! ..... East boundary ...... 
       CASE (3)

!         RETURN

       ! ..... South boundary ......
       CASE (4)

         j  = jsbcH(no,ide_t); 
         is = isbcH(no,ide_t); 
         ie = iebcH(no,ide_t);
         DO i = is, ie

           ! ... Map (i,j) into l-index
           l = ij2l(i,j);

           ! ... Cycle if N-column is dry
           IF ( .NOT. mask2d(i,j+1) ) CYCLE

           ! Compute the layer number for top & bottom wet v-pt
           kmy = MIN(kmz(lNC(l)), kmz(l))
           k1y =                 k1v(l)

           ! Compute number of wet layers & skip for dry cells
           nwlayers = (kmy-k1y) + 1
           IF(nwlayers < 1) CYCLE

           ! Compute explicit term
           ex(:,l) = 0.0	     
           DO k = k1y,kmy
                                                                 
            ! Horizontal advection Upwind differencing 
            vhbdry = (vhSB(k,i)+vhSBpp(k,i))/2.
            IF (hvSBpp(k,i)>ZERO) THEN
               vbdry  = vhSBpp(k,i)/hvSBpp(k,i)
            ELSE
               vbdry  = 0.0
            ENDIF
            uE = uhp(k,    lNC(l) ) +uhp(k,    l )
            uW = uhp(k,lWC(lNC(l))) +uhp(k,lWC(l))
            vN = vhp(k,    lNC(l) ) +vhp(k,    l )
            vS = vhp(k,        l  ) +vhbdry
            wU = wp (k,    lNC(l) ) +wp (k,    l )
            wD = wp (k+1,  lNC(l) ) +wp (k+1  ,l )        
            advy = ( (uE+ABS(uE))* vpp(k,    l ) +           &
                &    (uE-ABS(uE))* vpp(k,lEC(l)) -           &
                &    (uW+ABS(uW))* vpp(k,lWC(l)) -           &
                &    (uW-ABS(uW))* vpp(k,    l ) ) / fourdx  &
                & +( (vN+ABS(vN))* vpp(k,    l ) +           &
            	&    (vN-ABS(vN))* vpp(k,lNC(l)) -           &
                &    (vS+ABS(vS))* vbdry         -           &
                &    (vS-ABS(vS))* vpp(k,    l ) ) / fourdy  

            ! Vertical advection - Upwind in near bdry. cells
            advy=advy+((wU+ABS(wU)) * vpp(k  ,l) +           &
                       (wU-ABS(wU)) * vpp(k-1,l)) / 4.       &
                     -((wD+ABS(wD)) * vpp(k+1,l) +           &
                       (wD-ABS(wD)) * vpp(k  ,l)) / 4.
                   
            !.....Final explicit term.....
            ex(k,l) = vhpp(k,l) - twodt1*(advy*iadv)
                
           END DO
         END DO 

       END SELECT

     !....... Case 3 -- Submerged flow specified...........................
     CASE (3)

       ! Identify wse boundary as on the west, north, east, or south
       SELECT CASE ( iside(no) )

       ! ..... West boundary ......
       CASE (1)

         i  = isbc(no); 
         j  = jsbc(no); 
         ks = iebc(no);
         ke = jebc(no);

         ! Compute explicit term (use only advection) for v-points
		 ! to the North of a water column with EW flow specified
         IF (mask2d(i,j+1)) THEN 

           ! ... Map (i,j) into l-index
           l = ij2l(i,j);

           ! Compute explicit term (use only advection)
           ex(ks:ke,l) = 0.0

           DO k = ks,ke
                                                                   
             ! Horizontal advection - Upwind differencing  
             uhbdry = (uhWBpp(k,j)+uhWB  (k,j))/2.
             ! ubdry  =  uhWBpp(k,j)/huWBpp(k,j) 
             uE = uhp(k,    lNC(l) ) +uhp(k  ,    l )
             uW = uhp(k,lWC(lNC(l))) +uhbdry
             vN = vhp(k,    lNC(l) ) +vhp(k  ,    l )
             vS = vhp(k,    lSC(l) ) +vhp(k  ,    l )
             wU = wp (k,    lNC(l) ) +wp (k  ,    l )
             wD = wp (k+1,  lNC(l) ) +wp (k+1,    l )        
             advy = ( (uE+ABS(uE))* vpp(k,    l ) +           &
                      (uE-ABS(uE))* vpp(k,lEC(l)) -           &
                      (uW+ABS(uW))* vpp(k,lWC(l)) -           & !30062008
                      (uW-ABS(uW))* vpp(k,    l ) ) / fourdx  &
                   +( (vN+ABS(vN))* vpp(k,    l ) +           &
                      (vN-ABS(vN))* vpp(k,lNC(l)) -           &
                      (vS+ABS(vS))* vpp(k,lSC(l)) -           &
                      (vS-ABS(vS))* vpp(k,    l ) ) / fourdy   

             ! Vertical advection - Upwind for near bdry. cells
             advy=advy+((wU+ABS(wU)) * vpp(k  ,l) +           &
                        (wU-ABS(wU)) * vpp(k-1,l)) / 4.       &
                      -((wD+ABS(wD)) * vpp(k+1,l) +           &
                        (wD-ABS(wD)) * vpp(k  ,l)) / 4.
			      		                                  
             !.....Final explicit term.....
             ex(k,l) = vhpp(k,l) - twodt1*(advy*iadv)
               
           END DO

         ENDIF

         ! Compute explicit term (use only advection) for v-points
		 ! to the South of a water column with EW flow specified
         IF (mask2d(i,j-1)) THEN 

           ! ... Map (i,j) into l-index
           l = ij2l(i,j-1);
 
           ! Compute explicit term (use only advection)
           ex(ks:ke,l) = 0.0

           DO k = ks,ke
                                                                   
             ! Horizontal advection - Upwind differencing  
             uhbdry = (uhWBpp(k,j)+uhWB  (k,j))/2.
             ! ubdry  =  uhWBpp(k,j)/huWBpp(k,j) 
             uE = uhp(k,    lNC(l) ) +uhp(k  ,    l )
             uW = uhp(k,    lWC(l) ) +uhbdry
             vN = vhp(k,    lNC(l) ) +vhp(k  ,    l )
             vS = vhp(k,    lSC(l) ) +vhp(k  ,    l )
             wU = wp (k,    lNC(l) ) +wp (k  ,    l )
             wD = wp (k+1,  lNC(l) ) +wp (k+1,    l )        
             advy = ( (uE+ABS(uE))* vpp(k,    l ) +           &
                      (uE-ABS(uE))* vpp(k,lEC(l)) -           &
                      (uW+ABS(uW))* vpp(k,lWC(l)) -           & !30062008
                      (uW-ABS(uW))* vpp(k,    l ) ) / fourdx  &
                   +( (vN+ABS(vN))* vpp(k,    l ) +           &
                      (vN-ABS(vN))* vpp(k,lNC(l)) -           &
                      (vS+ABS(vS))* vpp(k,lSC(l)) -           &
                      (vS-ABS(vS))* vpp(k,    l ) ) / fourdy   

             ! Vertical advection - Upwind for near bdry. cells
             advy=advy+((wU+ABS(wU)) * vpp(k  ,l) +           &
                        (wU-ABS(wU)) * vpp(k-1,l)) / 4.       &
                      -((wD+ABS(wD)) * vpp(k+1,l) +           &
                        (wD-ABS(wD)) * vpp(k  ,l)) / 4.
			      		                                  
             !.....Final explicit term.....
             ex(k,l) = vhpp(k,l) - twodt1*(advy*iadv)
                
           END DO

         ENDIF 

       ! ..... North boundary ......
       CASE (2)
         
         j  = jsbc(no)-1; 
         i  = isbc(no); 
         ks = iebc(no);
         ke = jebc(no); 
         l = ij2l(i,j);

         ! Compute explicit term
         ex(ks:ke,l) = 0.0	     
         DO k = ks,ke
                                                                 
           ! Horizontal advection Upwind differencing 
           vhbdry = (vhNB(k,i)+vhNBpp(k,i))/2.
           vbdry  = vhNBpp(k,i)/hvNBpp(k,i)
           uE = uhp(k,    lNC(l) ) +uhp(k,    l )
           uW = uhp(k,lWC(lNC(l))) +uhp(k,lWC(l))
           vN = vhp(k,        l  ) +vhbdry
           vS = vhp(k,        l  ) +vhp(k,lSC(l))
           wU = wp (k,    lNC(l) ) +wp (k,    l )
           wD = wp (k+1,  lNC(l) ) +wp (k+1  ,l )        
           advy = ( (uE+ABS(uE))* vpp(k,    l ) +           &
               &    (uE-ABS(uE))* vpp(k,lEC(l)) -           &
               &    (uW+ABS(uW))* vpp(k,lWC(l)) -           &
               &    (uW-ABS(uW))* vpp(k,    l ) ) / fourdx  &
               & +( (vN+ABS(vN))* vpp(k,    l ) +           &
               &    (vN-ABS(vN))* vbdry         -           &
               &    (vS+ABS(vS))* vpp(k,lSC(l)) -           &
               &    (vS-ABS(vS))* vpp(k,    l ) ) / fourdy  

           ! Vertical advection - Upwind for near bdry. cells
           advy=advy+((wU+ABS(wU)) * vpp(k  ,l) +           &
                      (wU-ABS(wU)) * vpp(k-1,l)) / 4.       &
                    -((wD+ABS(wD)) * vpp(k+1,l) +           &
                      (wD-ABS(wD)) * vpp(k  ,l)) / 4.
                   
           !.....Final explicit term.....
           ex(k,l) = vhpp(k,l) - twodt1*(advy*iadv)
                
         END DO

       ! ..... East boundary ......
       CASE (3)

         i  = isbc(no); 
         j  = jsbc(no); 
         ks = iebc(no);
         ke = jebc(no);

         ! Compute explicit term (use only advection) for v-points
		 ! to the North of a water column with EW flow specified
         IF (mask2d(i,j+1)) THEN 

           ! ... Map (i,j) into l-index
           l = ij2l(i,j);

           ! Compute explicit term (use only advection)
           ex(ks:ke,l) = 0.0

           DO k = ks,ke
                                                                   
             ! Horizontal advection - Upwind differencing  
             uhbdry = (uhEBpp(k,j)+uhEB  (k,j))/2.
             ! ubdry  =  uhEBpp(k,j)/huEBpp(k,j) 
             uE = uhp(k,    lNC(l) ) +uhbdry
             uW = uhp(k,lWC(lNC(l))) +uhp(k  ,lWC(l))
             vN = vhp(k,    lNC(l) ) +vhp(k  ,    l )
             vS = vhp(k,    lSC(l) ) +vhp(k  ,    l )
             wU = wp (k,    lNC(l) ) +wp (k  ,    l )
             wD = wp (k+1,  lNC(l) ) +wp (k+1,    l )        
             advy = ( (uE+ABS(uE))* vpp(k,    l ) +           &
                      (uE-ABS(uE))* vpp(k,lEC(l)) -           &
                      (uW+ABS(uW))* vpp(k,lWC(l)) -           & !30062008
                      (uW-ABS(uW))* vpp(k,    l ) ) / fourdx  &
                   +( (vN+ABS(vN))* vpp(k,    l ) +           &
                      (vN-ABS(vN))* vpp(k,lNC(l)) -           &
                      (vS+ABS(vS))* vpp(k,lSC(l)) -           &
                      (vS-ABS(vS))* vpp(k,    l ) ) / fourdy   

             ! Vertical advection - Upwind for near bdry. cells
             advy=advy+((wU+ABS(wU)) * vpp(k  ,l) +           &
                        (wU-ABS(wU)) * vpp(k-1,l)) / 4.       &
                      -((wD+ABS(wD)) * vpp(k+1,l) +           &
                        (wD-ABS(wD)) * vpp(k  ,l)) / 4.
			      		                                  
             !.....Final explicit term.....
             ex(k,l) = vhpp(k,l) - twodt1*(advy*iadv)
               
           END DO

         ENDIF

         ! Compute explicit term (use only advection) for v-points
		 ! to the South of a water column with EW flow specified
         IF (mask2d(i,j-1)) THEN 

           ! ... Map (i,j) into l-index
           l = ij2l(i,j-1);
 
           ! Compute explicit term (use only advection)
           ex(ks:ke,l) = 0.0

           DO k = ks,ke
                                                                   
             ! Horizontal advection - Upwind differencing  
             uhbdry = (uhEBpp(k,j)+uhEB  (k,j))/2.
             ! ubdry  =  uhEBpp(k,j)/huEBpp(k,j) 
             uE = uhp(k,        l  ) +uhbdry
             uW = uhp(k,lWC(lNC(l))) +uhp(k  ,lWC(l))
             vN = vhp(k,    lNC(l) ) +vhp(k  ,    l )
             vS = vhp(k,    lSC(l) ) +vhp(k  ,    l )
             wU = wp (k,    lNC(l) ) +wp (k  ,    l )
             wD = wp (k+1,  lNC(l) ) +wp (k+1,    l )        
             advy = ( (uE+ABS(uE))* vpp(k,    l ) +           &
                      (uE-ABS(uE))* vpp(k,lEC(l)) -           &
                      (uW+ABS(uW))* vpp(k,lWC(l)) -           & !30062008
                      (uW-ABS(uW))* vpp(k,    l ) ) / fourdx  &
                   +( (vN+ABS(vN))* vpp(k,    l ) +           &
                      (vN-ABS(vN))* vpp(k,lNC(l)) -           &
                      (vS+ABS(vS))* vpp(k,lSC(l)) -           &
                      (vS-ABS(vS))* vpp(k,    l ) ) / fourdy   

             ! Vertical advection - Upwind for near bdry. cells
             advy=advy+((wU+ABS(wU)) * vpp(k  ,l) +           &
                        (wU-ABS(wU)) * vpp(k-1,l)) / 4.       &
                      -((wD+ABS(wD)) * vpp(k+1,l) +           &
                        (wD-ABS(wD)) * vpp(k  ,l)) / 4.
			      		                                  
             !.....Final explicit term.....
             ex(k,l) = vhpp(k,l) - twodt1*(advy*iadv)
                
           END DO

         ENDIF 

       ! ..... South boundary ......
       CASE (4)

         j  = jsbc(no); 
         i  = isbc(no); 
         ks = iebc(no);
         ke = jebc(no);
         l = ij2l(i,j);

         ! Compute explicit term
         ex(ks:ke,l) = 0.0	     
         DO k = ks,ke
                                                                 
           ! Horizontal advection Upwind differencing 
           vhbdry = (vhSB(k,i)+vhSBpp(k,i))/2.
           vbdry  =  vhSBpp(k,i)/hvSBpp(k,i)
           uE = uhp(k,    lNC(l) ) +uhp(k,    l )
           uW = uhp(k,lWC(lNC(l))) +uhp(k,lWC(l))
           vN = vhp(k,    lNC(l) ) +vhp(k,    l )
           vS = vhp(k,        l  ) +vhbdry
           wU = wp (k,    lNC(l) ) +wp (k,    l )
           wD = wp (k+1,  lNC(l) ) +wp (k+1  ,l )        
           advy = ( (uE+ABS(uE))* vpp(k,    l ) +           &
               &    (uE-ABS(uE))* vpp(k,lEC(l)) -           &
               &    (uW+ABS(uW))* vpp(k,lWC(l)) -           &
               &    (uW-ABS(uW))* vpp(k,    l ) ) / fourdx  &
               & +( (vN+ABS(vN))* vpp(k,    l ) +           &
               &    (vN-ABS(vN))* vpp(k,lNC(l)) -           &
               &    (vS+ABS(vS))* vbdry         -           &
               &    (vS-ABS(vS))* vpp(k,    l ) ) / fourdy  

           ! Vertical advection - Upwind in near bdry. cells
           advy=advy+((wU+ABS(wU)) * vpp(k  ,l) +           &
                      (wU-ABS(wU)) * vpp(k-1,l)) / 4.       &
                    -((wD+ABS(wD)) * vpp(k+1,l) +           &
                      (wD-ABS(wD)) * vpp(k  ,l)) / 4.
                   
           !.....Final explicit term.....
           ex(k,l) = vhpp(k,l) - twodt1*(advy*iadv)
                
         END DO 

       END SELECT
     END SELECT
   END DO
 

END SUBROUTINE MODexmom4openBCY



!************************************************************************
SUBROUTINE MODvel4openBC 
!************************************************************************
!
!  Purpose: To adjust the velocity estimates provided by the code at time
!           n+1, for the presence of open boundaries
!
!------------------------------------------------------------------------
 
   !.....Local variables.....
   INTEGER :: i, j, k, l, nn, is, ie, js, je, ks, ke, m, kms,no,ide_t
   REAL    :: deltaU,dru


   ide_t = omp_get_thread_num()+1
   !.....Loop over open boundaries.....
   DO nn = 1, nopenH(ide_t)
      no = noh2no(nn,ide_t)
      SELECT CASE ( itype(no) )
 
      !....... Case 1 -- wse specified.....................................
      CASE (1)
         
         SELECT CASE ( iside(no) )

         ! ... West boundary
         CASE (1)
            i = isbcH(no,ide_t); js = jsbcH(no,ide_t); je = jebcH(no,ide_t)
            
            DO j = js, je
              l = ij2l(i,j)
              vh (:,l) = 0.0;
              wp (:,l) = wp(:, lEC(l));   
                            
            ENDDO           

         ! ... North boundary
         CASE (2)
            j = jsbcH(no,ide_t); is = isbcHH(no,ide_t); ie = iebcHH(no,ide_t)
            
            DO i = is, ie
              l = ij2l(i,j)
              uh (:,l) = 0.0
              wp (:,l) = wp(:,lSC(l));
               
            ENDDO           

         ! ... East boundary
         CASE (3)
            i = isbcH(no,ide_t); js = jsbcH(no,ide_t); je = jebcH(no,ide_t)
            
            DO j = js, je
              l  = ij2l(i,j)
              vh (:,l) = 0.0;
              wp (:,l) = wp(:, lWC(l)); 
                
            ENDDO           
 
         ! ... South boundary
         CASE (4)
            j = jsbcH(no,ide_t); is = isbcHH(no,ide_t); ie = iebcHH(no,ide_t)
            
            DO i = is, ie
              l = ij2l(i,j)
              uh (:,l) = 0.0
              wp (:,l) = wp(:,lNC(l));
                
            ENDDO           

         END SELECT

      !....... Case 2 -- Free surface flow specified ......................
      CASE (2,4)

         SELECT CASE ( iside(no) )

         ! ... West boundary 
         CASE (1)
           i = isbcHH(no,ide_t); js = jsbcH(no,ide_t); je = jebcH(no,ide_t)
           DO j = js, je
             l   = ij2l(i,j)         
             kms = kmz(l); 
             wp(kms+1,l) = 0.0
             dru = 0.0
             DO k = kms,k1,-1
               dru = dru - (-uhWB  (k,j)         &
                            -uhWBpp(k,j))/twodx
               wp(k,l) = wp(k,l) + dru 


             END DO           
            ENDDO           

         ! ... North boundary
         CASE (2)
            j = jsbcH(no,ide_t); is = isbcHH(no,ide_t); ie = iebcHH(no,ide_t)
            DO i = is, ie
             l   = ij2l(i,j)         
             kms = kmz(l); 
             wp(kms+1,l) = 0.0
             dru = 0.0
             DO k = kms,k1,-1
               dru = dru -(vhNB  (k,i)+      &
                           vhNBpp(k,i))/twody
               wp(k,l) = wp(k,l) + dru


             END DO           
            ENDDO           

         ! ... East boundary
         CASE (3)
            i = isbcHH(no,ide_t); js = jsbcH(no,ide_t); je = jebcH(no,ide_t)
            DO j = js, je
             l   = ij2l(i,j)         
             kms = kmz(l); 
             wp(kms+1,l) = 0.0
             dru = 0.0
              DO k = kms,k1,-1
                dru = dru -(uhEB  (k,j) +       &
                            uhEBpp(k,j))/twodx
                wp(k,l) = wp(k,l) + dru    


              ENDDO                    
            ENDDO

         ! ... South boundary
         CASE (4)
            j = jsbcH(no,ide_t); is = isbcHH(no,ide_t); ie = iebcHH(no,ide_t)
            DO i = is, ie
             l   = ij2l(i,j)         
             kms = kmz(l); 
             wp(kms+1,l) = 0.0
             dru = 0.0
              DO k = kms,k1,-1
                dru = dru -(-vhSB  (k,i) +     & 
                            -vhSBpp(k,i))/twody; 
                wp(k,l) = wp(k,l) + dru
             ENDDO 
            ENDDO          

         END SELECT

      !....... Case 3 -- Submerged flow specified .........................
      CASE (3)

         SELECT CASE ( iside(no) )

         ! ... West boundary 
         CASE (1)

           i  = isbc(no); 
           j  = jsbc(no);
           ks = iebc(no); 
           ke = jebc(no); 
           l  = ij2l(i,j)         
           DO k = ks, ke
             deltaU = (uhWB(k,j)+uhWBpp(k,j))/twodx
             DO kms = k1, k
               wp(kms,l) = wp(kms,l)+deltaU                         
             END DO
           END DO

         ! ... North boundary
         CASE (2)

           i  = isbc(no); 
           j  = jsbc(no); 
           ks = iebc(no); 
           ke = jebc(no);
           l  = ij2l(i,j)        
           DO k = ks, ke
             deltaU = (vhNB(k,i)+vhNBpp(k,i))/twody
             DO kms = k1, k
               wp(kms,l) = wp(kms,l)-deltaU
             END DO    
           END DO           

         ! ... East boundary
         CASE (3)

           i  = isbc(no); 
           j  = jsbc(no); 
           ks = iebc(no); 
           ke = jebc(no); 
           l  = ij2l(i,j) 
           DO k = ks, ke
             deltaU = (uhEB(k,j)+uhEBpp(k,j))/twodx
             DO kms = k1, k
               wp(kms,l) = wp(kms,l)-deltaU   
             END DO
           END DO

         ! ... South boundary
         CASE (4)

           i  = isbc(no); 
           j  = jsbc(no); 
           ks = iebc(no); 
           ke = jebc(no);
           l  = ij2l(i,j)
           DO k = ks, ke
             deltaU = (vhSB(k,i)+vhSBpp(k,i))/twody
             DO kms = k1, k
               wp(kms,l) = wp(kms,l)+deltaU
             END DO    
           END DO           

         END SELECT

      END SELECT
   END DO

END SUBROUTINE MODvel4openBC

!************************************************************************
SUBROUTINE openbcSCA(thrs)
!************************************************************************
!
!  Purpose: To assign values of active scalars along open wse boundaries 
!           at new (n+1) time level. 
!
!------------------------------------------------------------------------

   REAL, INTENT(IN) :: thrs

   !.....Local variables.....
   REAL :: salbc, ss, delsal1, delsal2, delsal, dthrs_sal, advx
   INTEGER :: i, j, k, l, ios, nn, is, ie, js, je, kms, k1s,no,ide_t


   ide_t = omp_get_thread_num()+1
   ! ... Return if no open boundaries are specified
   IF (nopenH(ide_t) <= 0) RETURN

   !.....Loop over open boundaries.....
   DO nn = 1, nopenH(ide_t)

     no = noh2no(nn,ide_t)

     SELECT CASE ( itype(no) )      
 
     !.....Case 1 -- scalar specified on a wse BC .........................
     CASE (1)

       !.....Get new boundary value of scalar.....
       ! ... Set value of scalar at OpenBC nn - 
       dthrs_sal = dtsecOpenBC/3600.
       salbc = parab(0.,thrs,varsOpenBC(no,2,:),dthrs_sal)

       ! Identify wse boundary as on the west, north, east, or south
       SELECT CASE ( iside(no) )

       ! ..... West boundary ......
       CASE (1)
 
         i  = isbcH(no,ide_t); 
         js = jsbcH(no,ide_t); 
         je = jebcH(no,ide_t);
         DO j = js, je

           l = ij2l(i,j);

           ! Compute top & bottom layers
           kms = kmz(l)
           k1s = k1z(l)

           ! Compute explicit term (only using advection)
           DO k = k1s,kms
                                                                   
            ! Horizontal advection - Upwind differencing  
            IF ( uh(k,l) > 0.0 .AND. salbc > 0) THEN 
              sal(k,l) = salbc
            ELSE
              sal(k,l) = sal(k,lEC(l))
            ENDIF       
           ENDDO
         ENDDO

       ! ..... North boundary ......
       CASE (2)
 
         j  = jsbcH(no,ide_t); 
         is = isbcHH(no,ide_t); 
         ie = iebcHH(no,ide_t);
         DO i = is, ie

           l = ij2l(i,j);

           ! Compute top & bottom layers
           kms = kmz(l)
           k1s = k1z(l)

           ! Compute explicit term (only using advection)
           DO k = k1s,kms
                                                                   
            ! Horizontal advection - Upwind differencing  
            IF ( vh(k,lSC(l)) < 0.0 .AND. salbc > 0) THEN 
              sal(k,l) = salbc 
            ELSE
              sal(k,l) = sal(k,lSC(l))
            ENDIF

           ENDDO
         ENDDO

       ! ..... East boundary ......
       CASE (3)
 
         i  = isbcH(no,ide_t); 
         js = jsbcH(no,ide_t); 
         je = jebcH(no,ide_t);
         DO j = js, je

           l = ij2l(i,j);

           ! Compute top & bottom layers
           kms = kmz(l)
           k1s = k1z(l)

           ! Compute explicit term (only using advection)
           DO k = k1s,kms
                                                                   
            ! Horizontal advection - Upwind differencing  
            IF ( uh(k,lWC(l)) < 0.0 .AND. salbc > 0) THEN 
              sal(k,l) = salbc  
            ELSE
              sal(k,l) = sal(k,lWC(l))
            ENDIF       

           ENDDO
         ENDDO

       ! ..... South boundary ......
       CASE (4)
 
         j  = jsbcH(no,ide_t); 
         is = isbcHH(no,ide_t); 
         ie = iebcHH(no,ide_t);
         DO i = is, ie

           l = ij2l(i,j);

           ! Compute top & bottom layers
           kms = kmz(l)
           k1s = k1z(l)

           ! Compute explicit term (only using advection)
           DO k = k1s,kms
                                                                   
            ! Horizontal advection - Upwind differencing  
            IF ( vh(k,l) > 0.0 .AND. salbc > 0) THEN 
              sal(k,l) = salbc  
            ELSE
              sal(k,l) = sal(k,lNC(l))
            ENDIF

           ENDDO
         ENDDO

        ENDSELECT

      ! ... Case 2,3,4,... -- scalar specified on a flow BC ...............
      CASE (2:)

      END SELECT

   ENDDO


END SUBROUTINE openbcSCA

!************************************************************************
SUBROUTINE MODexsal4openbc(Bstart,Bend,Bex,thrs)
!************************************************************************
!
!  Purpose: Modifies the ex arrays in the scalar transport equation to 
!           account for inflows/outflows. 
!
!------------------------------------------------------------------------

   INTEGER, INTENT(IN) :: Bstart,Bend
   REAL, INTENT(IN) :: thrs
   REAL, DIMENSION (1:km1,Bstart:Bend+1), INTENT(INOUT) :: Bex

   !.....Local variables.....
   REAL    :: salbc, dthrs_sal, uE, uW, vN, vS, scE, scW, scN, scS
   INTEGER :: i, j, k, l, ios, nn, is, ie, js, je, kms, k1s
   INTEGER :: icl, m1, m2, niNGB, nbid,no,ide_t
   REAL    :: weight  

   ide_t = omp_get_thread_num()+1
   ! ... Return if no open boundaries are specified
   IF (nopenH(ide_t) <= 0) RETURN

   ! ... Initialize indexes
   niNGB = 0
!   print *,"SAthrs:",thrs
   !.....Loop over open boundaries.....
   DO nn = 1, nopenH(ide_t)

     no = noh2no(nn,ide_t)

     SELECT CASE ( itype(no) )      
 
     !.....Case 1 -- scalar specified on wse BC.....
     CASE (1)

        !  RETURN

     !.....Case 2 -- scalar specified on free surface flow BC..............
     CASE (2)

       !.....Get new boundary value of scalar.....
       ! ... Set value of scalar at OpenBC nn - 
       dthrs_sal = dtsecOpenBC/3600.
       salbc = parab(0.,thrs,varsOpenBC(no,2,:),dthrs_sal)

       ! Identify wse boundary as on the west, north, east, or south
       SELECT CASE ( iside(no) )

       ! ..... West boundary ......
       CASE (1)
 
         i  = isbcH(no,ide_t); 
         js = jsbcH(no,ide_t); 
         je = jebcH(no,ide_t);
         DO j = js, je
           l = ij2l(i,j);
           kms = kmz(l)
           k1s = k1z(l)
           DO k = k1s,kms                                                                 
             ! ... Define velocity at boundary face
             uW = uhWB(k,j) + uhWBpp(k,j);                               
             ! ... Define scalar   at boundary face  
             IF ( uW >= 0.0 .AND. salbc > 0) THEN 
               scW = salbc  
             ELSE
               scW = salpp(k,l)
             ENDIF
             ! ... Include boundary flux in ex- array 
             Bex(k,l) = Bex(k,l) +  uW * scW / twodx
           ENDDO
         ENDDO

       ! ..... North boundary ......
       CASE (2)
 
         j  = jsbcH(no,ide_t); 
         is = isbcH(no,ide_t); 
         ie = iebcH(no,ide_t);
         DO i = is, ie
           l = ij2l(i,j);
           kms = kmz(l)
           k1s = k1z(l)
           DO k = k1s,kms
             ! ... Define velocity at boundary face
             vN = vhNB(k,i) + vhNBpp(k,i);                               
             ! ... Define scalar   at boundary face  
             IF ( vN < 0.0 .AND. salbc > 0) THEN 
               scN = salbc  
             ELSE
               scN = salpp(k,l)
             ENDIF
             ! ... Include boundary flux in ex 
             Bex(k,l) = Bex(k,l) -  vN * scN / twody 
           ENDDO
         ENDDO

       ! ..... East boundary ......
       CASE (3)
 
         i  = isbcH(no,ide_t); 
         js = jsbcH(no,ide_t); 
         je = jebcH(no,ide_t);
         DO j = js, je
           l = ij2l(i,j);
           kms = kmz(l)
           k1s = k1z(l)
           DO k = k1s,kms
             ! ... Define velocity at boundary face
             uE = uhEB(k,j) + uhEBpp(k,j);                               
             ! ... Define scalar   at boundary face  
             IF ( uE < 0.0 .AND. salbc > 0) THEN 
               scE = salbc  
             ELSE
               scE = salpp(k,l)
             ENDIF
             ! ... Include boundary flux in ex 
             Bex(k,l) = Bex(k,l) -  uE * scE / twodx 
           ENDDO
         ENDDO

       ! ..... South boundary ......
       CASE (4)
 
         j  = jsbcH(no,ide_t); 
         is = isbcH(no,ide_t); 
         ie = iebcH(no,ide_t);
         DO i = is, ie
           l = ij2l(i,j);
           kms = kmz(l)
           k1s = k1z(l)
           DO k = k1s,kms
             ! ... Define velocity at boundary face
             vS = vhSB(k,i) + vhSBpp(k,i);                               
             ! ... Define scalar   at boundary face  
             IF ( vS > 0.0 .AND. salbc > 0) THEN 
               scS = salbc  
             ELSE
               scS = salpp(k,l)
             ENDIF
             ! ... Include boundary flux in ex 
             Bex(k,l) = Bex(k,l) +  vS * scS / twody
           ENDDO
         ENDDO

        ENDSELECT

     ! ... Case 3 -- scalar specified on submerged flow BC ................
     CASE (3)

       !.....Get new boundary value of scalar.....
       ! ... Set value of scalar at OpenBC nn - 
       dthrs_sal = dtsecOpenBC/3600.
       salbc = parab(0.,thrs,varsOpenBC(no,2,:),dthrs_sal)

       ! Identify wse boundary as on the west, north, east, or south
       SELECT CASE ( iside(no) )

       ! ..... West boundary ......
       CASE (1)
 
         i  = isbc(no); 
         js = jsbc(no); 
         je = jsbc(no);
         k1s= iebc(no);
         kms= jebc(no);
         DO j = js, je
           l = ij2l(i,j);
           DO k = k1s,kms               
             ! ... Define velocity at boundary face
             uW = uhWB(k,j) + uhWBpp(k,j);                               
             ! ... Define scalar   at boundary face  

             IF ( uW > 0.0 .AND. salbc > 0) THEN 
               scW = salbc  
             ELSE
               scW = salpp(k,l)
             ENDIF
             ! ... Include boundary flux in ex- array 
             Bex(k,l) = Bex(k,l) +  uW * scW / twodx 
           ENDDO
         ENDDO

       ! ..... North boundary ......
       CASE (2)

         j  = jsbc(no); 
         is = isbc(no); 
         ie = isbc(no);
         k1s= iebc(no);
         kms= jebc(no);
         DO i = is, ie
           l = ij2l(i,j);
           DO k = k1s,kms              
             ! ... Define velocity at boundary face
             vN = vhNB(k,i) + vhNBpp(k,i);                               
             ! ... Define scalar   at boundary face  
             IF ( vN < 0.0 .AND. salbc > 0) THEN 
               scN = salbc  
             ELSE
               scN = salpp(k,l)
             ENDIF
             ! ... Include boundary flux in ex 
             Bex(k,l) = Bex(k,l) -  vN * scN / twody 
           ENDDO
         ENDDO

       ! ..... East boundary ......
       CASE (3)
 
         i  = isbc(no); 
         js = jsbc(no); 
         je = jsbc(no);
         k1s= iebc(no);
         kms= jebc(no);
         DO j = js, je
           l = ij2l(i,j);
           DO k = k1s,kms               
             ! ... Define velocity at boundary face
             uE = uhEB(k,j) + uhEBpp(k,j);                               
             ! ... Define scalar   at boundary face  
             IF ( uE < 0.0 .AND. salbc > 0) THEN 
               scE = salbc  
             ELSE
               scE = salpp(k,l)
             ENDIF
             ! ... Include boundary flux in ex 
             Bex(k,l) = Bex(k,l) -  uE * scE / twodx 
           ENDDO
         ENDDO

       ! ..... South boundary ......
       CASE (4)

         i  = isbc(no); 
         j  = jsbc(no); 
         k1s= iebc(no);
         kms= jebc(no);
         l = ij2l(i,j);
         DO k = k1s, kms              
           ! ... Define velocity at boundary face
           vS = vhSB(k,i) + vhSBpp(k,i);                               
           ! ... Define scalar   at boundary face  
           IF ( vS > 0.0 .AND. salbc > 0) THEN 
             scS = salbc  
           ELSE
             scS = salpp(k,l)
           ENDIF
           ! ... Include boundary flux in ex 
           Bex(k,l) = Bex(k,l) +  vS * scS / twody 
         ENDDO
 
       END SELECT

       ! ... Case 4 - Active scalars at nested boundaries .................
       CASE (4) 

         ! ... Counter of nested grid boundaries
         niNGB = niNGB + 1;

         ! ... Define weighting coefficients for records 
         weight  = (thrs - thrsNGBp(no))/(thrsNGB(no)-thrsNGBp(no))
!         print *,"niNGB:",no,"thrs:",thrs,"thrsNGBp:",thrsNGBp
!         print *,"thrsNGB:",thrsNGB,"weight:",weight,"twodx:",twodx
         ! Identify boundary as on the west, north, east, or south
         SELECT CASE (iside(no))
         CASE(1)
           i = isbcH(no,ide_t); js = jsbcH(no,ide_t); je = jebcH(no,ide_t);
           icl = 0
           DO j = js, je; 
             l   = ij2l(i,j)
             DO k = k1, kmz(l)
             icl = icl + 1      
             ! ... Define scalar at boundary face  
             !scW = scNGB (icl,niNGB)*    weight + & ! FJR??
             !      scNGBp(icl,niNGB)*(1.-weight) ! FJR ??
             scW = scNGB (icl,no)*    weight + &
                   scNGBp(icl,no)*(1.-weight)            
             ! ... Define velocity at boundary face
             uW  = uhWB(k,j) + uhWBpp(k,j); 
             ! ... Re eefine scalar at boundary face if needed                               
             IF ( uW <= 0.0 ) scW = salpp(k,l)
             ! ... Include boundary flux in ex- array 
             Bex(k,l) = Bex(k,l) +  uW * scW / twodx
           ENDDO; ENDDO 
         CASE(3)
           i = isbcH(no,ide_t); js = jsbcH(no,ide_t); je = jebcH(no,ide_t);
           icl = 0
           DO j = js, je; 
             l   = ij2l(i,j)
             DO k = k1, kmz(l)
             icl = icl + 1 
             ! ... Define scalar at boundary face  
             !scE = scNGB (icl,niNGB)*    weight + & FJR??
             !      scNGBp(icl,niNGB)*(1.-weight) FJR??
             scE = scNGB (icl,no)*    weight + &
                   scNGBp(icl,no)*(1.-weight)            
             ! ... Define velocity at boundary face
             uE  = uhEB(k,j) + uhEBpp(k,j);                               
             ! ... Redefine scalar at boundary face if needed 
             IF ( uE >= 0.0 ) scE = salpp(k,l)
             ! ... Include boundary flux in ex 
             Bex(k,l) = Bex(k,l) -  uE * scE / twodx       
           ENDDO; ENDDO
         CASE(2)
           j = jsbcH(no,ide_t); is = isbcH(no,ide_t); ie = iebcH(no,ide_t);
           icl = siptNBI(no,ide_t) -1
           DO i = is, ie; 
             l   = ij2l(i,j)
             DO k = k1, kmz(l)
             icl = icl + 1      
             ! ... Define scalar at boundary face  
             scN = scNGB (icl,no)*    weight + &
                   scNGBp(icl,no)*(1.-weight)
             ! ... Define velocity at boundary face
             vN  = vhNB(k,i) + vhNBpp(k,i);                               
             ! ... Redefine scalar at boundary face if needed  
             IF ( vN >= 0.0 ) scN = salpp(k,l)
             ! ... Include boundary flux in ex 
             Bex(k,l) = Bex(k,l) -  vN * scN / twody 
           ENDDO; ENDDO
         CASE(4)
           j = jsbcH(no,ide_t); is = isbcH(no,ide_t); ie = iebcH(no,ide_t);
           icl = siptNBI(no,ide_t) -1
           DO i = is, ie; 
             l   = ij2l(i,j)
             DO k = k1, kmz(l)
             icl = icl + 1      
             ! ... Define scalar at boundary face  
             scS = scNGB (icl,no)*    weight + &
                   scNGBp(icl,no)*(1.-weight)
             ! ... Define velocity at boundary face
             vS  = vhSB(k,i) + vhSBpp(k,i);                               
             ! ... Define scalar   at boundary face  
             IF ( vS <= 0.0 ) scS = salpp(k,l)
             ! ... Include boundary flux in ex 
             Bex(k,l) = Bex(k,l) +  vS * scS / twody
           ENDDO; ENDDO

         END SELECT        

     END SELECT

   ENDDO

END SUBROUTINE MODexsal4openbc

!************************************************************************
SUBROUTINE openbcTracer (itr,thrs)
!************************************************************************
!
!  Purpose: To assign values of active scalars along open wse boundaries 
!           at new (n+1) time level. 
!
!------------------------------------------------------------------------

   ! ... Arguments
   INTEGER, INTENT (IN) :: itr
   REAL, INTENT(IN) :: thrs

   !.....Local variables.....
   INTEGER :: i, j, k, l, ios, nn, is, ie, js, je, kms, k1s,ide_t,no
   REAL    :: salbc, ss, delsal1, delsal2, delsal, dthrs_sal, &
              advx, uE, uW, vN, vS, wU, wD, twodt1

   ide_t = omp_get_thread_num()+1
   ! ... Return if no open boundaries are specified
   IF (nopenH(ide_t) <= 0) RETURN

   ! ... Constants used in solution
   twodt1 = twodt*tz

   !.....Loop over open boundaries.....
   DO nn = 1, nopenH(ide_t)
     no = noh2no(nn,ide_t)

     SELECT CASE ( itype(no) )      
 
     !.....Case 1 -- scalar specified on a wse BC.....
     CASE (1)

       !.....Get new boundary value of scalar.....
       ! ... Set value of scalar at OpenBC nn - 
       dthrs_sal = dtsecOpenBC/3600.
       salbc = parab(0.,thrs,varsOpenBC(no,2+itr,:),dthrs_sal)

       ! Identify wse boundary as on the west, north, east, or south
       SELECT CASE ( iside(no) )

       ! ..... West boundary ......
       CASE (1)
 
         i  = isbcH(no,ide_t); 
         js = jsbcH(no,ide_t); 
         je = jebcH(no,ide_t);
         DO j = js, je

           ! ... Map 3D into 2D array index
           l = ij2l(i,j);
           ! Compute top & bottom layers
           kms = kmz(l)
           k1s = k1z(l)
           ! Compute explicit term (only using advection)
           DO k = k1s,kms                                                                 
             ! Horizontal advection - Upwind differencing  
             IF ( uh(k,l) > 0.0 .AND. salbc > 0.0) THEN 
               tracer(k,l,itr) = salbc  
             ELSE
               tracer(k,l,itr) = tracer(k,lEC(l),itr)
             ENDIF       
           ENDDO
         ENDDO

       ! ..... North boundary ......
       CASE (2)
 
         j  = jsbcH(no,ide_t); 
         is = isbcHH(no,ide_t); 
         ie = iebcHH(no,ide_t);
         DO i = is, ie
           ! ... Map 3D into 2D array counter
           l = ij2l(i,j);
           ! Compute top & bottom layers
           kms = kmz(l)
           k1s = k1z(l)
           ! Compute explicit term (only using advection)
           DO k = k1s,kms                                                                  
             ! Horizontal advection - Upwind differencing  
             IF ( vh(k,lSC(l)) < 0.0 .AND. salbc > 0.0) THEN 
               tracer(k,l,itr) = salbc  
             ELSE
               tracer(k,l,itr) = tracer(k,lSC(l),itr)
             ENDIF
           ENDDO
         ENDDO

       ! ..... East boundary ......
       CASE (3)
 
         i  = isbcH(no,ide_t); 
         js = jsbcH(no,ide_t); 
         je = jebcH(no,ide_t);
         DO j = js, je
           ! ... Map 3D into 2D array index
           l = ij2l(i,j);
           ! Compute top & bottom layers
           kms = kmz(l)
           k1s = k1z(l)
           ! Compute explicit term (only using advection)
           DO k = k1s,kms                                                        
             ! Horizontal advection - Upwind differencing  
             IF ( uh(k,lWC(l)) < 0.0 .AND. salbc > 0.0) THEN 
               tracer(k,l,itr) = salbc  
             ELSE
               tracer(k,l,itr) = tracer(k,lWC(l),itr)
             ENDIF       
           ENDDO
         ENDDO

       ! ..... South boundary ......
       CASE (4)
 
         j  = jsbcH(no,ide_t); 
         is = isbcHH(no,ide_t); 
         ie = iebcHH(no,ide_t);
         DO i = is, ie
           ! ... Map 3D into 2D array index
           l = ij2l(i,j);
           ! Compute top & bottom layers
           kms = kmz(l)
           k1s = k1z(l)
           ! Compute explicit term (only using advection)
           DO k = k1s,kms                                                               
             ! Horizontal advection - Upwind differencing  
             IF ( vh(k,l) > 0.0 .AND. salbc > 0.0) THEN 
               tracer(k,l,itr) = salbc  
             ELSE
               tracer(k,l,itr) = tracer(k,lNC(l),itr)
             ENDIF
           ENDDO
         ENDDO

        ENDSELECT

      CASE (2:)

        CYCLE

      END SELECT

   ENDDO


END SUBROUTINE openbcTracer

!************************************************************************
SUBROUTINE MODexTracer4openbc (itr,Bstart,Bend,Bex,thrs)
!************************************************************************
!
!  Purpose: Modifies the ex arrays in the scalar transport equation to 
!           account for inflows/outflows. 
!
!------------------------------------------------------------------------

   ! ... Arguments
   INTEGER, INTENT (IN) :: itr
   REAL, INTENT(IN) :: thrs
   INTEGER, INTENT(IN) :: Bstart,Bend
   REAL, DIMENSION (1:km1,Bstart:Bend+1), INTENT(INOUT) :: Bex

   !.....Local variables.....
   REAL    :: salbc, dthrs_sal, uE, uW, vN, vS, scE, scW, scN, scS
   INTEGER :: i, j, k, l, ios, nn, is, ie, js, je, kms, k1s
   INTEGER :: icl, niNGB,ide_t,no
   REAL    :: weight  

   ide_t = omp_get_thread_num()+1
   ! ... Return if no open boundaries are specified
   IF (nopenH(ide_t) <= 0) RETURN

   niNGB = 0
!   print *,"TTthrs:",thrs
   !.....Loop over open boundaries.....
   DO nn = 1, nopenH(ide_t)
     no = noh2no(nn,ide_t)

     SELECT CASE ( itype(no) )      

     CASE (1) 

       CYCLE
 
     !.....Case 2 -- scalar specified on free surface flow BC............
     CASE (2)

       !.....Get new boundary value of scalar.....
       ! ... Set value of scalar at OpenBC nn - 
       dthrs_sal = dtsecOpenBC/3600.
       salbc = parab(0.,thrs,varsOpenBC(no,2+itr,:),dthrs_sal)

       ! Identify wse boundary as on the west, north, east, or south
       SELECT CASE ( iside(no) )   ! Changed to 'no' from 'nn' 12/2010 SWA

       ! ..... West boundary ......
       CASE (1)
         i  = isbcH(no,ide_t); 
         js = jsbcH(no,ide_t); 
         je = jebcH(no,ide_t);
         DO j = js, je
           l = ij2l(i,j);
           kms = kmz(l)
           k1s = k1z(l)
           DO k = k1s,kms                                                                 
             ! ... Define velocity at boundary face
             uW = uhWB(k,j) + uhWBpp(k,j);                               
             ! ... Define scalar   at boundary face  
             IF ( uW > 0.0 .AND. salbc > 0.0) THEN 
               scW = salbc  
             ELSE
               scW = tracerpp(k,l,itr)
             ENDIF
             ! ... Include boundary flux in ex- array 
             Bex(k,l) = Bex(k,l) +  uW * scW / twodx
           ENDDO
         ENDDO

       ! ..... North boundary ......
       CASE (2)
         j  = jsbcH(no,ide_t); 
         is = isbcH(no,ide_t); 
         ie = iebcH(no,ide_t);
         DO i = is, ie
           l = ij2l(i,j);
           kms = kmz(l)
           k1s = k1z(l)
           DO k = k1s,kms
             ! ... Define velocity at boundary face
             vN = vhNB(k,i) + vhNBpp(k,i);                               
             ! ... Define scalar   at boundary face  
             IF ( vN < 0.0 .AND. salbc > 0.0) THEN 
               scN = salbc  
             ELSE
               scN = tracerpp(k,l,itr)
             ENDIF
             ! ... Include boundary flux in ex 
             Bex(k,l) = Bex(k,l) -  vN * scN / twody 
           ENDDO
         ENDDO

       ! ..... East boundary ......
       CASE (3)
         i  = isbcH(no,ide_t); 
         js = jsbcH(no,ide_t); 
         je = jebcH(no,ide_t);
         DO j = js, je
           l = ij2l(i,j);
           kms = kmz(l)
           k1s = k1z(l)
           DO k = k1s,kms
             ! ... Define velocity at boundary face
             uE = uhEB(k,j) + uhEBpp(k,j);                               
             ! ... Define scalar   at boundary face  
             IF ( uE < 0.0 .AND. salbc > 0.0) THEN 
               scE = salbc  
             ELSE
               scE = tracerpp(k,l,itr)
             ENDIF
             ! ... Include boundary flux in ex 
             Bex(k,l) = Bex(k,l) -  uE * scE / twodx 
           ENDDO
         ENDDO

       ! ..... South boundary ......
       CASE (4)
         j  = jsbcH(no,ide_t); 
         is = isbcH(no,ide_t); 
         ie = iebcH(no,ide_t);
         DO i = is, ie
           l = ij2l(i,j);
           kms = kmz(l)
           k1s = k1z(l)
           DO k = k1s,kms
             ! ... Define velocity at boundary face
             vS = vhSB(k,i) + vhSBpp(k,i);                               
             ! ... Define scalar   at boundary face  
             IF ( vS > 0.0 .AND. salbc > 0.0) THEN 
               scS = salbc  
             ELSE
               scS = tracerpp(k,l,itr)
             ENDIF
             ! ... Include boundary flux in ex 
             Bex(k,l) = Bex(k,l) +  vS * scS / twody
           ENDDO
         ENDDO

        ENDSELECT

     ! ... Case 3 -- scalar specified on submerged flow BC ..............
     CASE (3)

       !.....Get new boundary value of scalar.....
       ! ... Set value of scalar at OpenBC nn - 
       dthrs_sal = dtsecOpenBC/3600.
       salbc = parab(0.,thrs,varsOpenBC(no,2+itr,:),dthrs_sal)

       ! Identify wse boundary as on the west, north, east, or south
       SELECT CASE ( iside(no) )

       ! ..... West boundary ......
       CASE (1)
         i   = isbc(no); 
         j   = jsbc(no); 
         k1s = iebc(no);
         kms = jebc(no);
         l   = ij2l(i,j);
         DO k = k1s,kms               
           ! ... Define velocity at boundary face
           uW = uhWB(k,j) + uhWBpp(k,j);                               
           ! ... Define scalar   at boundary face  
           IF ( uW > 0.0 .AND. salbc > 0.0) THEN 
             scW = salbc  
           ELSE
             scW = tracerpp(k,l,itr)
           ENDIF
           ! ... Include boundary flux in ex- array 
           Bex(k,l) = Bex(k,l) +  uW * scW / twodx 
         ENDDO

       ! ..... North boundary ......
       CASE (2)
         j   = jsbc(no); 
         i   = isbc(no); 
         k1s = iebc(no);
         kms = jebc(no);
         l   = ij2l(i,j);
         DO k = k1s,kms              
           ! ... Define velocity at boundary face
           vN = vhNB(k,i) + vhNBpp(k,i);                               
           ! ... Define scalar   at boundary face  
           IF ( vN < 0.0 .AND. salbc > 0.0) THEN 
             scN = salbc  
           ELSE
             scN = tracerpp(k,l,itr)
           ENDIF
           ! ... Include boundary flux in ex 
           Bex(k,l) = Bex(k,l) -  vN * scN / twody 
         ENDDO

       ! ..... East boundary ......
       CASE (3)
         i   = isbc(no); 
         j   = jsbc(no); 
         k1s = iebc(no);
         kms = jebc(no);
         l   = ij2l(i,j);
         DO k = k1s,kms               
           ! ... Define velocity at boundary face
           uE = uhEB(k,j) + uhEBpp(k,j);                               
           ! ... Define scalar   at boundary face  
           IF ( uE < 0.0 .AND. salbc > 0.0) THEN 
             scE = salbc  
           ELSE
             scE = tracerpp(k,l,itr)
           ENDIF
           ! ... Include boundary flux in ex 
           Bex(k,l) = Bex(k,l) -  uE * scE / twodx 
         ENDDO

       ! ..... South boundary ......
       CASE (4)
         i   = isbc(no); 
         j   = jsbc(no); 
         k1s = iebc(no);
         kms = jebc(no);
         l   = ij2l(i,j);
         DO k = k1s, kms              
           ! ... Define velocity at boundary face
           vS = vhSB(k,i) + vhSBpp(k,i);                               
           ! ... Define scalar   at boundary face  
           IF ( vS > 0.0 .AND. salbc > 0.0) THEN 
             scS = salbc  
           ELSE
             scS = tracerpp(k,l,itr)
           ENDIF
           ! ... Include boundary flux in ex 
           Bex(k,l) = Bex(k,l) +  vS * scS / twody 
         ENDDO
 
       END SELECT

     ! ... Case 4 - Nested grid boundaries ...............................
     CASE (4)

       ! ... Counter of nested grid boundaries
       niNGB = niNGB + 1;

       ! ... Define weighting coefficients for records 
       weight  = (thrs - thrsNGBp(no))/(thrsNGB(no)-thrsNGBp(no))
!       print *,"TTniNGB:",niNGB,"thrs:",thrs,"thrsNGBp:",thrsNGBp
!       print *,"TTthrsNGB:",thrsNGB,"weight:",weight,"twodx:",twodx
       ! Identify boundary as on the west, north, east, or south
       SELECT CASE (iside(no))

       ! ..... West boundary ......
       CASE (1)
         i = isbcH(no,ide_t); js = jsbcH(no,ide_t); je = jebcH(no,ide_t);
         icl = 0
         DO j = js, je; 
           l   = ij2l(i,j)
           DO k = k1, kmz(l)
           icl = icl + 1      
           ! ... Define scalar at boundary face  
           scW = trNGB (icl,niNGB,itr)*    weight + &
                 trNGBp(icl,niNGB,itr)*(1.-weight)
           ! ... Define velocity at boundary face
           uW  = uhWB(k,j) + uhWBpp(k,j); 
           ! ... Re eefine scalar at boundary face if needed                               
           IF ( uW <= 0.0 ) scW = tracerpp(k,l,itr)
           ! ... Include boundary flux in ex- array 
           Bex(k,l) = Bex(k,l) +  uW * scW / twodx
         ENDDO; ENDDO 

       ! ..... East boundary ......
       CASE (3)
         i = isbcH(no,ide_t); js = jsbcH(no,ide_t); je = jebcH(no,ide_t);
         icl = 0
         DO j = js, je; 
           l   = ij2l(i,j)
           DO k = k1, kmz(l)
           icl = icl + 1 
           ! ... Define scalar at boundary face  
           scE = trNGB (icl,niNGB,itr)*    weight + &
                 trNGBp(icl,niNGB,itr)*(1.-weight)
           ! ... Define velocity at boundary face
           uE  = uhEB(k,j) + uhEBpp(k,j);                               
           ! ... Redefine scalar at boundary face if needed 
           IF ( uE >= 0.0 ) scE = tracerpp(k,l,itr)
           ! ... Include boundary flux in ex 
           Bex(k,l) = Bex(k,l) -  uE * scE / twodx       
         ENDDO; ENDDO

       ! ..... North boundary ......
       CASE (2)
         j = jsbcH(no,ide_t); is = isbcH(no,ide_t); ie = iebcH(no,ide_t);
           icl = siptNBI(no,ide_t) -1
         DO i = is, ie; 
           l   = ij2l(i,j)
           DO k = k1, kmz(l)
           icl = icl + 1      
           ! ... Define scalar at boundary face  
           scN = trNGB (icl,niNGB,itr)*    weight + &
                 trNGBp(icl,niNGB,itr)*(1.-weight)
           ! ... Define velocity at boundary face
           vN  = vhNB(k,i) + vhNBpp(k,i);                               
           ! ... Redefine scalar at boundary face if needed  
           IF ( vN >= 0.0 ) scN = tracerpp(k,l,itr)
           ! ... Include boundary flux in ex 
           Bex(k,l) = Bex(k,l) -  vN * scN / twody 
         ENDDO; ENDDO

       ! ..... South boundary ......
       CASE(4)
         j = jsbcH(no,ide_t); is = isbcH(no,ide_t); ie = iebcH(no,ide_t);
           icl = siptNBI(no,ide_t) -1
         DO i = is, ie; 
           l   = ij2l(i,j)
           DO k = k1, kmz(l)
           icl = icl + 1      
           ! ... Define scalar at boundary face  
           scS = trNGB (icl,niNGB,itr)*    weight + &
                 trNGBp(icl,niNGB,itr)*(1.-weight)
           ! ... Define velocity at boundary face
           vS  = vhSB(k,i) + vhSBpp(k,i);                               
           ! ... Define scalar at boundary face if needed  
           IF ( vS <= 0.0 ) scS = tracerpp(k,l,itr)
           ! ... Include boundary flux in ex 
           Bex(k,l) = Bex(k,l) +  vS * scS / twody
         ENDDO; ENDDO
       END SELECT      

     END SELECT

   ENDDO

END SUBROUTINE MODexTracer4openbc

!************************************************************************
SUBROUTINE surfbc0
!************************************************************************
!
!  Purpose: This routine is called at the beginning of the program
!           to open files with heat boundary condition data, to
!           read the boundary condition time series data, to 
!           assign heat_sources at time t=0. It uses
!           the same scheme as used in openbc routines to read bc values. 
!
!------------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: ios, istat, nn, j, npsurfbc, nvsurfbc, imet
   CHARACTER(LEN=14) :: surfbcfmt, metxyfmt

   
   SELECT CASE (ifSurfBC) 

   ! ... Surface boundary conditions set to constant values (no heat flux)
   CASE (0) 

     RETURN

   ! .... Surface boundary conditions read from files - PRE-PROCESS mode
   CASE (1)

     !               ----- Open files with heatflux surface bc data-----   
     OPEN (UNIT=i53, FILE='surfbc.txt', STATUS="OLD", IOSTAT=ios)
     IF (ios /= 0) CALL open_error ( "Error opening surfbc.txt", ios )

     !               -----Read files with heatflux surface bc data-----
     ! Skip over first six header records in salinity boundary condition file 
     READ (UNIT=i53, FMT='(/////)', IOSTAT=ios)
     IF (ios /= 0) CALL input_error ( ios, 101 )
   
     ! Read number of points in file from seventh header record
     READ (UNIT=i53, FMT='(10X,I7)', IOSTAT=ios) npsurfbc
     IF (ios /= 0) CALL input_error ( ios, 102 )
   
     ! Allocate space for the array of data
     ALLOCATE ( surfbc1(nvSurfbcP,npSurfbc), STAT=istat )
     IF (istat /= 0) CALL allocate_error ( istat, 103 )
   
     ! Write the format of the data records into an internal file
     WRITE (UNIT=surfbcfmt, FMT='("(10X,",I3,"G11.2)")') nvSurfbcP
   
     ! Read data array and store it in memory
     DO j = 1, npSurfbc
       READ (UNIT=i53, FMT=surfbcfmt, IOSTAT=ios) &
            (surfbc1(nn,j), nn = 1, nvSurfbcP)
       IF (ios /= 0) CALL input_error ( ios, 104 )
     END DO
   
     !           ----- Assign heat flux terms at time t=0.0-----
     eta = surfbc1(1,1)
     Qsw = surfbc1(2,1)
     Qn  = surfbc1(3,1) 

     !           ----- Assign momentum flux terms at t=0.0 -----
     cdw = surfbc1(4,1)
     uair= surfbc1(5,1)
     vair= surfbc1(6,1)

     ! ... Set heat sources for each cell
     CALL DistributeQsw
     CALL DistributeMomentumHeatSources

   ! .... Surface boundary conditions RUN-TIME (I) mode
   CASE (2)

     !               ----- Open files with heatflux surface bc data-----   
     OPEN (UNIT=i53, FILE='surfbc.txt', STATUS="OLD", IOSTAT=ios)
     IF (ios /= 0) CALL open_error ( "Error opening surfbc.txt", ios )

     !               -----Read files with heatflux surface bc data-----
     ! Skip over first six header records in salinity boundary condition file 
     READ (UNIT=i53, FMT='(/////)', IOSTAT=ios)
     IF (ios /= 0) CALL input_error ( ios, 105 )
   
     ! Read number of points in file from seventh header record
     READ (UNIT=i53, FMT='(10X,I7)', IOSTAT=ios) npSurfbc
     IF (ios /= 0) CALL input_error ( ios, 106 )
   
     ! Allocate space for the array of data
     ALLOCATE ( surfbc1(nvSurfbcR,npSurfbc), STAT=istat )
     IF (istat /= 0) CALL allocate_error ( istat, 107 )
   
     ! Write the format of the data records into an internal file
     WRITE (UNIT=surfbcfmt, FMT='("(10X,",I3,"G11.2)")') nvSurfbcR
   
     ! Read data array and store it in memory
     DO j = 1, npSurfbc
       READ (UNIT=i53, FMT=surfbcfmt, IOSTAT=ios) &
            (surfbc1(nn,j), nn = 1, nvSurfbcR)
       IF (ios /= 0) CALL input_error ( ios, 108 )
     END DO
   
     !           ----- Assign heat flux terms at time t=0.0-----
     eta = surfbc1(1,1)
     Qsw = surfbc1(2,1)
     Ta  = surfbc1(3,1)        ! Data in oC
     Pa  = surfbc1(4,1)        ! Pascals
     Rh  = surfbc1(5,1)        ! fraction (i.e. < 1)
     Cc  = surfbc1(6,1)        ! fraction (i.e. < 1)

     !           ----- Assign momentum flux terms at t=0.0 -----
     cdw = surfbc1(7,1)
     uair= surfbc1(8,1)
     vair= surfbc1(9,1)

     ! ... Set heat sources for each cell
     CALL DistributeQsw
     CALL DistributeMomentumHeatSources

   ! .... Surface boundary conditions RUN-TIME (II) mode
   CASE (3)

     !               ----- Open files with heatflux surface bc data-----   
     OPEN (UNIT=i53, FILE='surfbc.txt', STATUS="OLD", IOSTAT=ios)
     IF (ios /= 0) CALL open_error ( "Error opening surfbc.txt", ios )

     !               -----Read files with heatflux surface bc data-----
     ! Skip over first six header records in salinity boundary condition file 
     READ (UNIT=i53, FMT='(/////)', IOSTAT=ios)
     IF (ios /= 0) CALL input_error ( ios, 109 )
   
     ! Read number of points in file from seventh header record
     READ (UNIT=i53, FMT='(10X,I7)', IOSTAT=ios) npSurfbc
     IF (ios /= 0) CALL input_error ( ios, 110 )
   
     ! Allocate space for the array of data
     ALLOCATE ( surfbc1(nvSurfbcR,npSurfbc), STAT=istat )
     IF (istat /= 0) CALL allocate_error ( istat, 111 )
   
     ! Write the format of the data records into an internal file
     WRITE (UNIT=surfbcfmt, FMT='("(10X,",I3,"G11.2)")') nvSurfbcR
   
     ! Read data array and store it in memory
     DO j = 1, npSurfbc
       READ (UNIT=i53, FMT=surfbcfmt, IOSTAT=ios) &
            (surfbc1(nn,j), nn = 1, nvSurfbcR)
       IF (ios /= 0) CALL input_error ( ios, 112 )
     END DO
   
     !           ----- Assign heat flux terms at time t=0.0-----
     eta = surfbc1(1,1)        ! m-1
     Qsw = surfbc1(2,1)        ! W/m2
     Ta  = surfbc1(3,1)        ! Data in oC
     Pa  = surfbc1(4,1)        ! Pascals
     Rh  = surfbc1(5,1)        ! fraction (i.e. < 1)
     Qlw = surfbc1(6,1)        ! W/m2

     !           ----- Assign momentum flux terms at t=0.0 -----
     cdw = surfbc1(7,1)
     uair= surfbc1(8,1)
     vair= surfbc1(9,1)

     ! ... Set heat sources for each cell
     CALL DistributeQsw
     CALL DistributeMomentumHeatSources

   CASE (10) ! Space & Time varying met. variables - Heat budget on run-time (I) mode

     !               ----- Open files with heatflux surface bc data-----
     OPEN (UNIT=i53, FILE='surfbc.txt', STATUS="OLD", IOSTAT=ios)
     IF (ios /= 0) CALL open_error ( "Error opening surfbc.txt", ios )

     !               -----Read files with heatflux surface bc data-----

     ! Skip over first six header records in boundary condition file
     READ (UNIT=i53, FMT='(/////)', IOSTAT=ios)
     IF (ios /= 0) CALL input_error ( ios, 113 )  

     ! Read number of met stations in file from seventh header record
     READ (UNIT=i53, FMT='(10X,I7)', IOSTAT=ios) nmetstat
     IF (ios /= 0) CALL input_error ( ios, 114 )

     ! Allocate space for metxy matrix containing x,y locations 
     ! of all met stations & variables that store individual met records
     ! for a given time step for each of the variables
     ALLOCATE ( metxy (nmetstat*2), &
                Qsw2D (nmetstat  ), Ta2D  (nmetstat),           &
                RH2D  (nmetstat  ), Cc2D  (nmetstat),           &
                uair2D(nmetstat  ), vair2D(nmetstat), STAT=istat )
     IF (istat /= 0) CALL allocate_error ( istat, 115 )

     ! Allocate space for weightst matrix containing weighting coefficients 
     ! assigned to each of the met stations for each grid point
     ALLOCATE ( weightst(im1,jm1,nmetstat), STAT=istat )
     IF (istat /= 0) CALL allocate_error ( istat, 116 )

     ! Write the format of the data records into an internal file
     WRITE (UNIT=metxyfmt, FMT='("(10X,",I3,"G11.2)")') nmetstat*2

     ! Read x,y position (in grid units) for each met station
     READ (UNIT=i53, FMT=metxyfmt, IOSTAT=ios)(metxy(nn), nn = 1, nmetstat*2)
     
     IF (ios /= 0) CALL input_error ( ios, 117 )
 
     ! ... Initialize Interpolation Schemes (weights for Barnes)
     CALL InitializeInterpolationMethods

     ! Read number of records in file from seventh header record
     READ (UNIT=i53, FMT='(10X,I7)', IOSTAT=ios) npSurfbc
     IF (ios /= 0) CALL input_error ( ios, 110 )
   
     ! Allocate space for the array of data
     nvSurfbc = 6 * nmetstat + 2 
     ALLOCATE ( surfbc1(nvSurfbc,npSurfbc), STAT=istat )
     IF (istat /= 0) CALL allocate_error ( istat, 111 )
   
     ! Write the format of the data records into an internal file
     WRITE (UNIT=surfbcfmt, FMT='("(10X,",I3,"G11.2)")') nvSurfbc
   
     ! Read data array and store it in memory
     DO j = 1, npSurfbc
         READ (UNIT=i53, FMT=surfbcfmt, IOSTAT=ios) &
		      (surfbc1(nn,j), nn = 1, nvSurfbc)
         IF (ios /= 0) CALL input_error ( ios, 112 )
     END DO
 
     !           ----- Assign surfBC at time t=0.0-----
     eta           = surfbc1(1             ,1)
     Pa            = surfbc1(2             ,1)
     DO imet = 1, nmetstat
       Qsw2D (imet)  = surfbc1((imet-1)*6 + 3,1)
       Ta2D  (imet)  = surfbc1((imet-1)*6 + 4,1) 
       RH2D  (imet)  = surfbc1((imet-1)*6 + 5,1) 
       Cc2D  (imet)  = surfbc1((imet-1)*6 + 6,1) 
       uair2D(imet)  = surfbc1((imet-1)*6 + 7,1)
       vair2D(imet)  = surfbc1((imet-1)*6 + 8,1)
     ENDDO

     ! ... Distribute heat and momentum sources entering through free surface
     CALL DistributeQsw
     CALL DistributeMomentumHeatSources

   CASE (11) ! Space & Time varying met. variables - Heat budget on run-time (II) mode

     !               ----- Open files with heatflux surface bc data-----
     OPEN (UNIT=i53, FILE='surfbc.txt', STATUS="OLD", IOSTAT=ios)
     IF (ios /= 0) CALL open_error ( "Error opening surfbc.txt", ios )

     !               -----Read files with heatflux surface bc data-----

     ! Skip over first six header records in boundary condition file
     READ (UNIT=i53, FMT='(/////)', IOSTAT=ios)
     IF (ios /= 0) CALL input_error ( ios, 113 )  

     ! Read number of met stations in file from seventh header record
     READ (UNIT=i53, FMT='(10X,I7)', IOSTAT=ios) nmetstat
     IF (ios /= 0) CALL input_error ( ios, 114 )

     ! Allocate space for metxy matrix containing x,y locations 
     ! of all met stations & variables that store individual met records
     ! for a given time step for each of the variables
     ALLOCATE ( metxy (nmetstat*2), weightn(nmetstat,nmetstat),          &  !new array created in types
                Qsw2D (nmetstat  ), Ta2D  (nmetstat),           &
                RH2D  (nmetstat  ), Qlw2D (nmetstat),           &
                uair2D(nmetstat  ), vair2D(nmetstat), STAT=istat )
     IF (istat /= 0) CALL allocate_error ( istat, 115 )

     ! Allocate space for weightst matrix containing weighting coefficients 
     ! assigned to each of the met stations for each grid point
     ALLOCATE ( weightst(im1,jm1,nmetstat), STAT=istat )
     IF (istat /= 0) CALL allocate_error ( istat, 116 )

     ! Write the format of the data records into an internal file
     WRITE (UNIT=metxyfmt, FMT='("(10X,",I3,"G11.2)")') nmetstat*2

     ! Read x,y position (in grid units) for each met station
     READ (UNIT=i53, FMT=metxyfmt, IOSTAT=ios)(metxy(nn), nn = 1, nmetstat*2)
     print *,"met:",metxy
     IF (ios /= 0) CALL input_error ( ios, 117 )

     do nn=1,nmetstat*2
           print *,"nmet",nn,": ",metxy(nn)
     end do
 
     ! ... Initialize Interpolation Schemes (weights for Barnes)
     CALL InitializeInterpolationMethods
     print *,"hola"
     ! Read number of records in file from seventh header record
     READ (UNIT=i53, FMT='(10X,I7)', IOSTAT=ios) npSurfbc
     IF (ios /= 0) CALL input_error ( ios, 110 )
     print *,"hola2"
     ! Allocate space for the array of data
     nvSurfbc = 6 * nmetstat + 2 
     ALLOCATE ( surfbc1(nvSurfbc,npSurfbc), STAT=istat )
     IF (istat /= 0) CALL allocate_error ( istat, 111 )
     print *,"hola3"
     ! Write the format of the data records into an internal file
     WRITE (UNIT=surfbcfmt, FMT='("(10X,",I3,"G11.2)")') nvSurfbc
     print *,"hola4"
     ! Read data array and store it in memory
     DO j = 1, npSurfbc
         READ (UNIT=i53, FMT=surfbcfmt, IOSTAT=ios) &
		      (surfbc1(nn,j), nn = 1, nvSurfbc)
         IF (ios /= 0) CALL input_error ( ios, 112 )
     END DO
     print *,"hola5"
     !           ----- Assign surfBC at time t=0.0-----
     eta           = surfbc1(1             ,1)
     Pa            = surfbc1(2             ,1)
     DO imet = 1, nmetstat
       Qsw2D (imet)  = surfbc1((imet-1)*6 + 3,1)
       Ta2D  (imet)  = surfbc1((imet-1)*6 + 4,1) 
       RH2D  (imet)  = surfbc1((imet-1)*6 + 5,1) 
       Qlw2D (imet)  = surfbc1((imet-1)*6 + 6,1) 
       uair2D(imet)  = surfbc1((imet-1)*6 + 7,1)
       vair2D(imet)  = surfbc1((imet-1)*6 + 8,1)
     ENDDO
     print *,"hola6"
     ! ... Distribute heat and momentum sources entering through free surface
     CALL DistributeQsw
     print *,"hola61"
     CALL DistributeMomentumHeatSources
     print *,"hola7"
   END SELECT
   
END SUBROUTINE surfbc0

!************************************************************************
SUBROUTINE surfbc(n,istep,thrs)
!************************************************************************
!
!  Purpose: To define heat & momentum fluxes through the free surface
!           at each time  
!
!------------------------------------------------------------------------

   INTEGER,INTENT(IN) :: n,istep
   REAL, INTENT(IN) :: thrs

   !.....Local variables.....
   REAL    :: dthrs_surfbc
   INTEGER :: i, j, k, ios, nn, is, ie, js, je, kb, isalin, itest, imet,liter,l

   SELECT CASE (ifSurfBC) 

   !               ----- No surface bc data----------------
   CASE (0)

      DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

       l = id_column(liter)

      uair(l) = -wa * SIN(pi*phi/180.);
      vair(l) = -wa * COS(pi*phi/180.);
      cdw(l)  =  cw
      END DO

   !               ----- Use surface bc data from file ----
   CASE(1) ! Heat budget on preprocess mode - shortwave radiative and 
           ! net heat fluxes (including longwave & sensible & latent) as input 
           ! Space uniform & time varying

     !.....Return from subroutine on trapezoidal steps (except if n=1).....
     IF (n > 1) THEN; IF (istep == 2) RETURN; END IF

     !.....Interpolate heat & momentum flux vars. to time n .....
     dthrs_surfbc = dtSurfbc/3600.
     eta = parab(0.,thrs,surfbc1(1,:),dthrs_surfbc)
     Qsw = parab(0.,thrs,surfbc1(2,:),dthrs_surfbc)
     Qn  = parab(0.,thrs,surfbc1(3,:),dthrs_surfbc) 
     cdw = parab(0.,thrs,surfbc1(4,:),dthrs_surfbc)
     uair= parab(0.,thrs,surfbc1(5,:),dthrs_surfbc)
     vair= parab(0.,thrs,surfbc1(6,:),dthrs_surfbc) 

     ! ... Calculate 3D-spatially variable sources 
     CALL DistributeQswH
     CALL DistributeMomentumHeatSourcesH(n,istep)

   CASE (2) ! Heat budget on run-time mode (I) - shortwave fluxes as input; 
            ! longwave & latent & sensible heat fluxes calculated.
            ! space uniform & time varying 

     !.....Return from subroutine on trapezoidal steps (except if n=1).....
     IF (n > 1) THEN; IF (istep == 2) RETURN; END IF

     !.....Interpolate heat & momentum flux vars. values to present time step .....
     dthrs_surfbc = dtSurfbc/3600.
     eta = parab(0.,thrs,surfbc1(1,:),dthrs_surfbc)
     Qsw = parab(0.,thrs,surfbc1(2,:),dthrs_surfbc)
     Ta  = parab(0.,thrs,surfbc1(3,:),dthrs_surfbc)  
     Pa  = parab(0.,thrs,surfbc1(4,:),dthrs_surfbc)
     Rh  = parab(0.,thrs,surfbc1(5,:),dthrs_surfbc)
     Cc  = parab(0.,thrs,surfbc1(6,:),dthrs_surfbc) 
     cdw = parab(0.,thrs,surfbc1(7,:),dthrs_surfbc)
     uair= parab(0.,thrs,surfbc1(8,:),dthrs_surfbc)
     vair= parab(0.,thrs,surfbc1(9,:),dthrs_surfbc) 

     ! ... Calculate 3D-spatially variable heat sources
     CALL DistributeQswH
     CALL DistributeMomentumHeatSourcesH(n,istep)

   CASE (3) ! Heat budget on run-time mode (II) - radiative (short & longwave)
            ! fluxes as input; latent & sensible heat fluxes calculated. 
            ! Space uniform & time varying 

     !.....Return from subroutine on trapezoidal steps (except if n=1).....
     IF (n > 1) THEN; IF (istep == 2) RETURN; END IF

     !.....Interpolate heat & momentum flux vars. values to present time step .....
     dthrs_surfbc = dtSurfbc/3600.
     eta = parab(0.,thrs,surfbc1(1,:),dthrs_surfbc)
     Qsw = parab(0.,thrs,surfbc1(2,:),dthrs_surfbc)
     Ta  = parab(0.,thrs,surfbc1(3,:),dthrs_surfbc)  
     Pa  = parab(0.,thrs,surfbc1(4,:),dthrs_surfbc)
     Rh  = parab(0.,thrs,surfbc1(5,:),dthrs_surfbc)
     Qlw = parab(0.,thrs,surfbc1(6,:),dthrs_surfbc) 
     cdw = parab(0.,thrs,surfbc1(7,:),dthrs_surfbc)
     uair= parab(0.,thrs,surfbc1(8,:),dthrs_surfbc)
     vair= parab(0.,thrs,surfbc1(9,:),dthrs_surfbc) 

     ! ... Calculate 3D-spatially variable heat sources
     CALL DistributeQswH
     CALL DistributeMomentumHeatSourcesH(n,istep)

   CASE (10) ! Spatially & time varying surface BC - Heat Budget calculated
              ! on RUN-TIME (I) mode

     !.....Return from subroutine on trapezoidal steps (except if n=1).....
     IF (n > 1) THEN; IF (istep == 2) RETURN; END IF
     
     !.....Interpolate heat & momentum flux vars. values to present time step .....
     DO imet = 1, nmetstat
       dthrs_surfbc  = dtSurfbc/3600.
       eta           = parab(0.,thrs,surfbc1(1,:),dthrs_surfbc)
       Pa            = parab(0.,thrs,surfbc1(2,:),dthrs_surfbc)
       Qsw2D (imet)  = parab(0.,thrs,surfbc1((imet-1)*6 + 3,:),dthrs_surfbc)
       Ta2D  (imet)  = parab(0.,thrs,surfbc1((imet-1)*6 + 4,:),dthrs_surfbc) 
       RH2D  (imet)  = parab(0.,thrs,surfbc1((imet-1)*6 + 5,:),dthrs_surfbc) 
       Cc2D  (imet)  = parab(0.,thrs,surfbc1((imet-1)*6 + 6,:),dthrs_surfbc) 
       uair2D(imet)  = parab(0.,thrs,surfbc1((imet-1)*6 + 7,:),dthrs_surfbc)
       vair2D(imet)  = parab(0.,thrs,surfbc1((imet-1)*6 + 8,:),dthrs_surfbc) 
     ENDDO

     ! ... Distribute heat and momentum sources entering through free surface
     CALL DistributeQswH
     CALL DistributeMomentumHeatSourcesH(n,istep)


   CASE (11) ! Spatially & time varying surface BC - Heat Budget calculated
              ! on RUN-TIME (II) mode

     !.....Return from subroutine on trapezoidal steps (except if n=1).....
     IF (n > 1) THEN; IF (istep == 2) RETURN; END IF
     
     !.....Interpolate heat & momentum flux vars. values to present time step .....
     DO imet = 1, nmetstat
       dthrs_surfbc  = dtSurfbc/3600.
       eta           = parab(0.,thrs,surfbc1(1,:),dthrs_surfbc)
       Pa            = parab(0.,thrs,surfbc1(2,:),dthrs_surfbc)
       Qsw2D (imet)  = parab(0.,thrs,surfbc1((imet-1)*6 + 3,:),dthrs_surfbc)
       Ta2D  (imet)  = parab(0.,thrs,surfbc1((imet-1)*6 + 4,:),dthrs_surfbc) 
       RH2D  (imet)  = parab(0.,thrs,surfbc1((imet-1)*6 + 5,:),dthrs_surfbc) 
       Qlw2D (imet)  = parab(0.,thrs,surfbc1((imet-1)*6 + 6,:),dthrs_surfbc) 
       uair2D(imet)  = parab(0.,thrs,surfbc1((imet-1)*6 + 7,:),dthrs_surfbc)
       vair2D(imet)  = parab(0.,thrs,surfbc1((imet-1)*6 + 8,:),dthrs_surfbc) 
     ENDDO

     ! ... Distribute heat and momentum sources entering through free surface
     CALL DistributeQswH
     CALL DistributeMomentumHeatSourcesH(n,istep)

   END SELECT

END SUBROUTINE surfbc

!************************************************************************
 SUBROUTINE DistributeQswH
!************************************************************************
!
!  Purpose: To apportion the solar irradiance penetrating the 
!           lake through the surface among the layers in each
!           water column
!
!------------------------------------------------------------------------

  ! ... Local variables
  INTEGER:: i, j, k , l, kb, nwlayers, k1s, kms,liter
  REAL   :: remFr, zfromt0
  REAL, DIMENSION (im1,jm1) :: htot
  REAL, DIMENSION (1  :km1) :: zfromt

  !.....Sweep over wet pressure points
  DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

       l = id_column(liter)

       ! ... Map 2D-l into 3D(i,j) indexes .....
       i = l2i(l); j = l2j(l);

       ! ... Define top & bottom wet layers.....
       kms = kmz(l);
       k1s = k1z(l);    
       nwlayers = (kms-k1s) + 1

       SELECT CASE (nwlayers)

       CASE (1)

         QswFr(k1s,l) = 1.

       CASE (2:)

         ! ... Compute the array of vertical distances from the 
         !     free surface to the top of each layer
         zfromt(k1s) = 0.0
         DO k = k1s+1, kms+1          
           zfromt(k)    = zfromt(k-1) + hp(k-1,l)
           QswFr(k-1,l) = SolarFr(zfromt(k-1)) - SolarFr(zfromt(k))
         END DO
         ! ... Redistribute Qsw not absorbed once it reaches the bottom 
         !     Probably it is not needed in deep lakes but could account for
         !     overheating for shallower lakes. 
         remFr = 1. - SUM ( QswFr (k1s:kms,l) ) 
         IF ( remFr > 1.e-6 ) THEN
           QswFr (k1s:kms,l) = QswFr (k1s:kms,l) + remFr / nwlayers
         END IF

       END SELECT

    END DO

END SUBROUTINE DistributeQswH

!************************************************************************
 SUBROUTINE DistributeQsw
!************************************************************************
!
!  Purpose: To apportion the solar irradiance penetrating the 
!           lake through the surface among the layers in each
!           water column
!
!------------------------------------------------------------------------

  ! ... Local variables
  INTEGER:: i, j, k , l, kb, nwlayers, k1s, kms
  REAL   :: remFr, zfromt0
  REAL, DIMENSION (im1,jm1) :: htot
  REAL, DIMENSION (1  :km1) :: zfromt

  !.....Sweep over wet pressure points
  DO l = 1, lm

       ! ... Map 2D-l into 3D(i,j) indexes .....
       i = l2i(l); j = l2j(l);

       ! ... Define top & bottom wet layers.....
       kms = kmz(l);
       k1s = k1z(l);    
       nwlayers = (kms-k1s) + 1

       SELECT CASE (nwlayers)

       CASE (1)

         QswFr(k1s,l) = 1.

       CASE (2:)

         ! ... Compute the array of vertical distances from the 
         !     free surface to the top of each layer
         zfromt(k1s) = 0.0
         DO k = k1s+1, kms+1          
           zfromt(k)    = zfromt(k-1) + hp(k-1,l)
           QswFr(k-1,l) = SolarFr(zfromt(k-1)) - SolarFr(zfromt(k))
         END DO
         ! ... Redistribute Qsw not absorbed once it reaches the bottom 
         !     Probably it is not needed in deep lakes but could account for
         !     overheating for shallower lakes. 
         remFr = 1. - SUM ( QswFr (k1s:kms,l) ) 
         IF ( remFr > 1.e-6 ) THEN
           QswFr (k1s:kms,l) = QswFr (k1s:kms,l) + remFr / nwlayers
         END IF

       END SELECT

    END DO

END SUBROUTINE DistributeQsw

!************************************************************************
REAL FUNCTION SolarFr ( depth ) 
!************************************************************************
!
!  Purpose: To calculate the attenuation of solar irradiance penetrating  
!           the water column through the free surface. It uses formulation
!           proposed in Henderson-Sellers' Engineering Limnology Eq. 2.25 
!           According to the authors this equation is only valid for 
!           eta (attenuation coefficient) larger than 0.1
!
!------------------------------------------------------------------------

  ! ... Arguments
  REAL, INTENT (IN) :: depth

  ! ... Local variables
  REAL            :: BetaSol, z
  REAL, PARAMETER :: zA = 0.60
  
  z = depth
  IF ( eta .GE. 0.1) THEN
    BetaSol = 0.265*LOG (eta)+0.614; 
    IF ( z<zA ) THEN
      SolarFr = (1.-BetaSol*z/zA)
    ELSE
      SolarFr = (1.-BetaSol)*EXP(-eta*(z-zA)); 
    ENDIF
  ELSE
    SolarFr = EXP(-eta*z)
  ENDIF

END FUNCTION SolarFr

!************************************************************************
 SUBROUTINE DistributeMomentumHeatSources
!************************************************************************
!
!  Purpose: To construct a 2D met field from discrete variables & 
!           to contruct heat sources for each computational cell 
!
!------------------------------------------------------------------------

  ! ... Local variables
  INTEGER         :: i, j, k,l, k1s, kms, irg, imet, ii, jj, nF
  REAL            :: esw, ea, Lv, EmAir, Qbri, Qbra
  REAL            :: coeffE, coeffH, ws, wsi, ws0
  REAL, DIMENSION(nmetstat) :: Qswresid,Taresid,RHresid,Qlwresid,uresid,vresid

  ! ... Local arrays containing parameters
  INTEGER, DIMENSION (6) :: range 
  REAL,    DIMENSION (5) :: a_d, b_d, p_d
  REAL,    DIMENSION (5) :: a_h, b_h, p_h, c_h
  REAL,    DIMENSION (5) :: a_e, b_e, p_e, c_e
  
  ! .... Parameter definition
  REAL, PARAMETER :: EmWater = 0.97
  REAL, PARAMETER :: StephanBoltzman = 5.6697e-8
  REAL, PARAMETER :: Al = 0.03
  REAL, PARAMETER :: Cp = 4181.6
  REAL, PARAMETER :: SpecificHeatAir = 1012.0

  ! ... Define parameters on first call
!  IF ( n == 1) THEN
    range = (/  0.000, 2.2000, 5.0000, 8.0000, 25.00, 50.0 /)
    a_d   = (/  0.000, 0.7710, 0.8670, 1.2000, 0.000 /)
    b_d   = (/  1.080, 0.0858, 0.0667, 0.0250, 0.073 /)
    p_d   = (/ -0.150, 1.0000, 1.0000, 1.0000, 1.000 /)
    a_h   = (/  0.000, 0.9270, 1.1500, 1.1700, 1.652 /)
    b_h   = (/  1.185, 0.0521, 0.0100, 0.0075,-0.017 /)
    c_h   = (/  0.000, 0.0000, 0.0000,-4.5E-4, 0.000 /)
    p_h   = (/ -0.157, 1.0000, 1.0000, 1.0000, 1.000 /)
    a_e   = (/  0.000, 0.9690, 1.1800, 1.1960, 1.680 /)
    b_e   = (/  1.230, 0.0521, 0.0100, 0.0080,-0.016 /)
    c_e   = (/  0.000, 0.0000, 0.0000,-4.0E-4, 0.000 /)
    p_e   = (/ -0.160, 1.0000, 1.0000, 1.0000, 1.000 /)
!  ENDIF

  ! ... Initialize Heat Source
  HeatSource = 0.0E0
  print *,"hello"
  SELECT CASE (ifSurfbc)

  CASE (1) ! Heat budget on PRE-PROCESS mode - spatially uniform conditions

    DO l = 1,lm ;
 
      ! ... Map 2D-l into 3D(i,j) indexes ........................
      i = l2i(l); j = l2j(l);

      ! ... Define top & bottom wet layers........................
      k1s = k1z(l);    
      kms = kmz(l);

      ! ... Add non-penetrative components to surface layer ......
      HeatSource(k1s,l)=HeatSource(k1s,l)+ (Qn-Qsw)
 
      ! ... Add penetrative components to water column & 
      !     express heat sources in temperature units ............
      DO k = k1s, kms
        HeatSource(k,l) = (HeatSource(k,l) + & 
          Qsw*QswFr(k,l))/((rhop(k,l)+1000.)*Cp)
       IF (Qsw<-100.) THEN; HeatSource(k,l) = 0.0E0; ENDIF
      END DO

    ENDDO

  CASE (2) ! Heat budget on run-time mode (I) - spatially uniform conditions

    ! ... Define variables used in cals. and equal to all surface cells
    ea    = saturated_vapor_pressure (Ta+273.) * Rh
    EmAir = 0.642 * ( ea / (Ta+273.) )**0.14285714 		
    EmAir = EmAir * ( 1. + 0.17 * Cc**2. ) 
    Qbri  = EmAir * StephanBoltzman * (Ta+273.) ** 4. 
    Qbra  = Qbri * ( 1. - Al ) 

    ! ... Loop over surface cells 
    DO l = 1, lm;

      ! ... Map 2D-l into 3D(i,j) indexes ........................
      i = l2i(l); j = l2j(l);

      ! ... Define bulk aerodynamic coefficients for neutral conditions 
      !     based on Kondo, 1975. Air-sea bulk transfer coefficients in diabatic
      !     conditions. In Boundary-Layer Meteorology 9 (1975) 91-112
      ws = SQRT(uair(l)**2. + vair(l)**2.)    
      IF      (ws>range(1).AND.ws<=range(2)) THEN; irg = 1
      ELSE IF (ws>range(2).AND.ws<=range(3)) THEN; irg = 2
      ELSE IF (ws>range(3).AND.ws<=range(4)) THEN; irg = 3
      ELSE IF (ws>range(4).AND.ws<=range(5)) THEN; irg = 4
      ELSE                                       ; irg = 5; END IF
      coeffH = (a_h(irg)+b_h(irg)*ws**p_h(irg)+c_h(irg)*(ws-8.)**2.)*1.e-3
      coeffE = (a_e(irg)+b_e(irg)*ws**p_e(irg)+c_e(irg)*(ws-8.)**2.)*1.e-3

      ! ... Define top & bottom wet layers........................
      k1s = k1z(l);    
      kms = kmz(l);

      ! ... Add non-penetrative components to surface cells
      esw = saturated_vapor_pressure(salp(k1s,l)+273.)
      Lv  = latent_heat_vaporization(salp(k1s,l)+273.)
      ! a. Long wave radiation 
      HeatSource(k1s,l)  = Qbra - EmWater * StephanBoltzman *    &
                           (salp(k1s,l)+273.) ** 4.
      ! b. Latent & Sensible heat fluxes
      HeatSource(k1s,l) = HeatSource(k1s,l) -                    &
      & rhoair*             Lv * coeffE*(esw-ea)*0.622/Pa * ws - & 
      & rhoair*SpecificHeatAir * coeffH*(salp(k1s,l)-Ta ) * ws

      ! ... Add penetrative components to water column & 
      !     express heat sources in temperature units ............
      DO k = k1s, kms
        HeatSource(k,l) = (HeatSource(k,l) +                     & 
          Qsw*QswFr(k,l))/((rhop(k,l)+1000.)*Cp)
        IF (Qsw<-100.) THEN; HeatSource(k,l) = 0.0E0; ENDIF
      END DO

    ENDDO

  CASE (3) ! Heat budget on run-time mode (II) - spatially uniform conditions

    ! ... Define variables used in cals. and equal to all surface cells
    ea    = saturated_vapor_pressure (Ta+273.) * Rh

    ! ... Loop over surface cells 
    DO l = 1, lm;

      ! ... Map 2D-l into 3D(i,j) indexes ........................
      i = l2i(l); j = l2j(l);

      ! ... Define bulk aerodynamic coefficients for neutral conditions 
      !     based on Kondo, 1975. Air-sea bulk transfer coefficients in diabatic
      !     conditions. In Boundary-Layer Meteorology 9 (1975) 91-112
      ws = SQRT(uair(l)**2. + vair(l)**2.)    
      IF      (ws>range(1).AND.ws<=range(2)) THEN; irg = 1
      ELSE IF (ws>range(2).AND.ws<=range(3)) THEN; irg = 2
      ELSE IF (ws>range(3).AND.ws<=range(4)) THEN; irg = 3
      ELSE IF (ws>range(4).AND.ws<=range(5)) THEN; irg = 4
      ELSE                                       ; irg = 5; END IF
      coeffH = (a_h(irg)+b_h(irg)*ws**p_h(irg)+c_h(irg)*(ws-8.)**2.)*1.e-3
      coeffE = (a_e(irg)+b_e(irg)*ws**p_e(irg)+c_e(irg)*(ws-8.)**2.)*1.e-3

      ! ... Define top & bottom wet layers........................
      k1s = k1z(l);    
      kms = kmz(l);

      ! ... Add non-penetrative components to surface cells
      esw = saturated_vapor_pressure(salp(k1s,l)+273.)
      Lv  = latent_heat_vaporization(salp(k1s,l)+273.)
      ! a. Long wave radiation (incoming LW as measured) 
      HeatSource(k1s,l)  = Qlw * (1.- Al) -                      &
                           EmWater * StephanBoltzman *           &
                           (salp(k1s,l)+273.) ** 4.
      ! b. Latent & Sensible heat fluxes
      HeatSource(k1s,l) = HeatSource(k1s,l) -                    &
      & rhoair*             Lv * coeffE*(esw-ea)*0.622/Pa * ws - & 
      & rhoair*SpecificHeatAir * coeffH*(salp(k1s,l)-Ta ) * ws

      ! ... Add penetrative components to water column & 
      !     express heat sources in temperature units ............
      DO k = k1s, kms
        HeatSource(k,l) = (HeatSource(k,l) +                     & 
          Qsw*QswFr(k,l))/((rhop(k,l)+1000.)*Cp)
        IF (Qsw<-100.) THEN; HeatSource(k,l) = 0.0E0; ENDIF
      END DO

    ENDDO

  CASE (10) ! Interpolation

    ! Calculate residuals from first Barnes pass at met station points (added 12/2010 SWA)
    CALL CalculateMetResiduals (Qswresid,Taresid,RHresid,Qlwresid,uresid,vresid)

    ! ... Loop over surface cells 
    DO l = 1, lm;

      ! ... Map 2D-l into 3D(i,j) indexes ........................
      i = l2i(l); j = l2j(l);

      ! ... Interpolate met records to grid points (added 12/2010 SWA)
      CALL InterpMetData (i,j,Qswresid,Taresid,RHresid,Qlwresid,uresid,vresid, &
                        & Qsw,Ta,RH,Qlw,uair(l),vair(l))

      ! ... Wind speed over the water column
      ws  = SQRT (uair(l)**2.+vair(l)**2.);
        
      ! ... Define bulk aerodynamic coefficients for heat transfer under neutral conditions 
      !     based on Kondo, 1975. Air-sea bulk transfer coefficients in diabatic
      !     conditions. In Boundary-Layer Meteorology 9 (1975) 91-112
      ! ... Define range of wind speed
      IF      (ws>range(1).AND.ws<=range(2)) THEN; irg = 1
      ELSE IF (ws>range(2).AND.ws<=range(3)) THEN; irg = 2
      ELSE IF (ws>range(3).AND.ws<=range(4)) THEN; irg = 3
      ELSE IF (ws>range(4).AND.ws<=range(5)) THEN; irg = 4
      ELSE                                       ; irg = 5; END IF
      coeffH = (a_h(irg)+b_h(irg)*ws**p_h(irg)+c_h(irg)*(ws-8.)**2.)*1.e-3
      coeffE = (a_e(irg)+b_e(irg)*ws**p_e(irg)+c_e(irg)*(ws-8.)**2.)*1.e-3

      ! ... Calculate drag coefficient (Amorocho & DeVries)  
      cdw(l) = 0.0015 * 1./(1.0+ EXP((12.5-ws)/1.56))+0.00104      

      ! ... Define variables used in cals. and equal to all surface cells
      ea    = saturated_vapor_pressure (Ta+273.) * Rh
      EmAir = 0.642 * ( ea / (Ta+273.) )**0.14285714 		
      EmAir = EmAir * ( 1. + 0.17 * Cc**2. ) 
      Qbri  = EmAir * StephanBoltzman * (Ta+273.) ** 4. 
      Qbra  = Qbri * ( 1. - Al )

      ! ... Define top & bottom wet layers........................
      k1s = k1z(l);    
      kms = kmz(l);

      ! ... Add non-penetrative components to surface cells
      esw = saturated_vapor_pressure(salp(k1s,l)+273.)
      Lv  = latent_heat_vaporization(salp(k1s,l)+273.)
      ! a. Long wave radiation 
      HeatSource(k1s,l)  = Qbra - EmWater * StephanBoltzman *    &
                           (salp(k1s,l)+273.) ** 4.
      ! b. Latent & Sensible heat fluxes
      HeatSource(k1s,l) = HeatSource(k1s,l) -                    &
      & rhoair*             Lv * coeffE*(esw-ea)*0.622/Pa * ws - & 
      & rhoair*SpecificHeatAir * coeffH*(salp(k1s,l)-Ta ) * ws

      ! ... Add penetrative components to water column & 
      !     express heat sources in temperature units ............
      DO k = k1s, kms
        HeatSource(k,l) = (HeatSource(k,l) +                     & 
          Qsw*QswFr(k,l))/((rhop(k,l)+1000.)*Cp)
        IF (Qsw<-100.) THEN; HeatSource(k,l) = 0.0E0; ENDIF
      END DO

    ENDDO
  
  CASE (11) ! Interpolation
    print *,"hello1"
    ! Calculate residuals from first Barnes pass at met station points (added 12/2010 SWA)
    CALL CalculateMetResiduals (Qswresid,Taresid,RHresid,Qlwresid,uresid,vresid)
    print *,"hello2"
    ! ... Interpolate other variables & compute heat fluxes
    DO l = 1, lm;

      ! ... Map 2D-l into 3D(i,j) indexes ........................
      i = l2i(l); j = l2j(l);

      ! ... Interpolate met records to grid points (added 12/2010 SWA)
      CALL InterpMetData (i,j,Qswresid,Taresid,RHresid,Qlwresid,uresid,vresid, &
                        & Qsw,Ta,RH,Qlw,uair(l),vair(l))

      ! ... Wind speed over the water column  
      ws  = SQRT (uair(l)**2.+vair(l)**2.);

      ! ... Define bulk aerodynamic coefficients for heat transfer under neutral conditions 
      !     based on Kondo, 1975. Air-sea bulk transfer coefficients in diabatic
      !     conditions. In Boundary-Layer Meteorology 9 (1975) 91-112
      ! ... Define range of wind speed 
      IF      (ws>range(1).AND.ws<=range(2)) THEN; irg = 1
      ELSE IF (ws>range(2).AND.ws<=range(3)) THEN; irg = 2
      ELSE IF (ws>range(3).AND.ws<=range(4)) THEN; irg = 3
      ELSE IF (ws>range(4).AND.ws<=range(5)) THEN; irg = 4
      ELSE                                       ; irg = 5; END IF
      coeffH = (a_h(irg)+b_h(irg)*ws**p_h(irg)+c_h(irg)*(ws-8.)**2.)*1.e-3
      coeffE = (a_e(irg)+b_e(irg)*ws**p_e(irg)+c_e(irg)*(ws-8.)**2.)*1.e-3

      ! ... Calculate drag coefficient (Amorocho & DeVries)  
      cdw(l) = 0.0015 * 1./(1.0+ EXP((12.5-ws)/1.56))+0.00104      

      ! ... Define variables used in cals. and equal to all surface cells
      Qbra  = Qlw * ( 1. - Al ) 

      ! ... Define top & bottom wet layers........................
      k1s = k1z(l);    
      kms = kmz(l);

      ! ... Add non-penetrative components to surface cells
      ea  = saturated_vapor_pressure(Ta         +273.) * Rh
      esw = saturated_vapor_pressure(salp(k1s,l)+273.)
      Lv  = latent_heat_vaporization(salp(k1s,l)+273.)
      ! a. Long wave radiation 
      HeatSource(k1s,l)  = Qbra - EmWater * StephanBoltzman *    &
                           (salp(k1s,l)+273.) ** 4.
      ! b. Latent & Sensible heat fluxes
      HeatSource(k1s,l) = HeatSource(k1s,l) -                    &
      & rhoair*             Lv * coeffE*(esw-ea)*0.622/Pa * ws - & 
      & rhoair*SpecificHeatAir * coeffH*(salp(k1s,l)-Ta ) * ws

      ! ... Add penetrative components to water column & 
      !     express heat sources in temperature units ............
      DO k = k1s, kms
        HeatSource(k,l) = (HeatSource(k,l) +                     & 
          Qsw*QswFr(k,l))/((rhop(k,l)+1000.)*Cp)
        IF (Qsw<-100.) THEN; HeatSource(k,l) = 0.0E0; ENDIF
      END DO

    ENDDO
    print *,"hello3"
  END SELECT

END SUBROUTINE DistributeMomentumHeatSources

!************************************************************************
 SUBROUTINE DistributeMomentumHeatSourcesH(n,istep)
!************************************************************************
!
!  Purpose: To construct a 2D met field from discrete variables & 
!           to contruct heat sources for each computational cell 
!
!------------------------------------------------------------------------

  INTEGER,INTENT(IN) :: n,istep

  ! ... Local variables
  INTEGER         :: i, j, k,l, k1s, kms, irg, imet, ii, jj, nF,liter
  REAL            :: esw, ea, Lv, EmAir, Qbri, Qbra
  REAL            :: coeffE, coeffH, ws, wsi, ws0
  REAL, DIMENSION(nmetstat) :: Qswresid,Taresid,RHresid,Qlwresid,uresid,vresid

  ! ... Local arrays containing parameters
  INTEGER, DIMENSION (6) :: range 
  REAL,    DIMENSION (5) :: a_d, b_d, p_d
  REAL,    DIMENSION (5) :: a_h, b_h, p_h, c_h
  REAL,    DIMENSION (5) :: a_e, b_e, p_e, c_e
  
  ! .... Parameter definition
  REAL, PARAMETER :: EmWater = 0.97
  REAL, PARAMETER :: StephanBoltzman = 5.6697e-8
  REAL, PARAMETER :: Al = 0.03
  REAL, PARAMETER :: Cp = 4181.6
  REAL, PARAMETER :: SpecificHeatAir = 1012.0

  ! ... Define parameters on first call
!  IF ( n == 1) THEN
    range = (/  0.000, 2.2000, 5.0000, 8.0000, 25.00, 50.0 /)
    a_d   = (/  0.000, 0.7710, 0.8670, 1.2000, 0.000 /)
    b_d   = (/  1.080, 0.0858, 0.0667, 0.0250, 0.073 /)
    p_d   = (/ -0.150, 1.0000, 1.0000, 1.0000, 1.000 /)
    a_h   = (/  0.000, 0.9270, 1.1500, 1.1700, 1.652 /)
    b_h   = (/  1.185, 0.0521, 0.0100, 0.0075,-0.017 /)
    c_h   = (/  0.000, 0.0000, 0.0000,-4.5E-4, 0.000 /)
    p_h   = (/ -0.157, 1.0000, 1.0000, 1.0000, 1.000 /)
    a_e   = (/  0.000, 0.9690, 1.1800, 1.1960, 1.680 /)
    b_e   = (/  1.230, 0.0521, 0.0100, 0.0080,-0.016 /)
    c_e   = (/  0.000, 0.0000, 0.0000,-4.0E-4, 0.000 /)
    p_e   = (/ -0.160, 1.0000, 1.0000, 1.0000, 1.000 /)

!  ENDIF


  ! ... Initialize Heat Source
DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

       l = id_column(liter)
  HeatSource(:,l) = 0.0E0
END DO

  SELECT CASE (ifSurfbc)

  CASE (1) ! Heat budget on PRE-PROCESS mode - spatially uniform conditions

    DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

       l = id_column(liter)
 
      ! ... Map 2D-l into 3D(i,j) indexes ........................
      i = l2i(l); j = l2j(l);

      ! ... Define top & bottom wet layers........................
      k1s = k1z(l);    
      kms = kmz(l);

      ! ... Add non-penetrative components to surface layer ......
      HeatSource(k1s,l)=HeatSource(k1s,l)+ (Qn-Qsw)
 
      ! ... Add penetrative components to water column & 
      !     express heat sources in temperature units ............
      DO k = k1s, kms
        HeatSource(k,l) = (HeatSource(k,l) + & 
          Qsw*QswFr(k,l))/((rhop(k,l)+1000.)*Cp)
       IF (Qsw<-100.) THEN; HeatSource(k,l) = 0.0E0; ENDIF
      END DO

    ENDDO

  CASE (2) ! Heat budget on run-time mode (I) - spatially uniform conditions

    ! ... Define variables used in cals. and equal to all surface cells
    ea    = saturated_vapor_pressure (Ta+273.) * Rh
    EmAir = 0.642 * ( ea / (Ta+273.) )**0.14285714 		
    EmAir = EmAir * ( 1. + 0.17 * Cc**2. ) 
    Qbri  = EmAir * StephanBoltzman * (Ta+273.) ** 4. 
    Qbra  = Qbri * ( 1. - Al ) 

    ! ... Loop over surface cells 
    DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

       l = id_column(liter)

      ! ... Map 2D-l into 3D(i,j) indexes ........................
      i = l2i(l); j = l2j(l);

      ! ... Define bulk aerodynamic coefficients for neutral conditions 
      !     based on Kondo, 1975. Air-sea bulk transfer coefficients in diabatic
      !     conditions. In Boundary-Layer Meteorology 9 (1975) 91-112
      ws = SQRT(uair(l)**2. + vair(l)**2.)    
      IF      (ws>range(1).AND.ws<=range(2)) THEN; irg = 1
      ELSE IF (ws>range(2).AND.ws<=range(3)) THEN; irg = 2
      ELSE IF (ws>range(3).AND.ws<=range(4)) THEN; irg = 3
      ELSE IF (ws>range(4).AND.ws<=range(5)) THEN; irg = 4
      ELSE                                       ; irg = 5; END IF
      coeffH = (a_h(irg)+b_h(irg)*ws**p_h(irg)+c_h(irg)*(ws-8.)**2.)*1.e-3
      coeffE = (a_e(irg)+b_e(irg)*ws**p_e(irg)+c_e(irg)*(ws-8.)**2.)*1.e-3

      ! ... Define top & bottom wet layers........................
      k1s = k1z(l);    
      kms = kmz(l);

      ! ... Add non-penetrative components to surface cells
      esw = saturated_vapor_pressure(salp(k1s,l)+273.)
      Lv  = latent_heat_vaporization(salp(k1s,l)+273.)
      ! a. Long wave radiation 
      HeatSource(k1s,l)  = Qbra - EmWater * StephanBoltzman *    &
                           (salp(k1s,l)+273.) ** 4.
      ! b. Latent & Sensible heat fluxes
      HeatSource(k1s,l) = HeatSource(k1s,l) -                    &
      & rhoair*             Lv * coeffE*(esw-ea)*0.622/Pa * ws - & 
      & rhoair*SpecificHeatAir * coeffH*(salp(k1s,l)-Ta ) * ws

      ! ... Add penetrative components to water column & 
      !     express heat sources in temperature units ............
      DO k = k1s, kms
        HeatSource(k,l) = (HeatSource(k,l) +                     & 
          Qsw*QswFr(k,l))/((rhop(k,l)+1000.)*Cp)
        IF (Qsw<-100.) THEN; HeatSource(k,l) = 0.0E0; ENDIF
      END DO

    ENDDO

  CASE (3) ! Heat budget on run-time mode (II) - spatially uniform conditions

    ! ... Define variables used in cals. and equal to all surface cells
    ea    = saturated_vapor_pressure (Ta+273.) * Rh

    ! ... Loop over surface cells 
    DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

       l = id_column(liter)

      ! ... Map 2D-l into 3D(i,j) indexes ........................
      i = l2i(l); j = l2j(l);

      ! ... Define bulk aerodynamic coefficients for neutral conditions 
      !     based on Kondo, 1975. Air-sea bulk transfer coefficients in diabatic
      !     conditions. In Boundary-Layer Meteorology 9 (1975) 91-112
      ws = SQRT(uair(l)**2. + vair(l)**2.)    
      IF      (ws>range(1).AND.ws<=range(2)) THEN; irg = 1
      ELSE IF (ws>range(2).AND.ws<=range(3)) THEN; irg = 2
      ELSE IF (ws>range(3).AND.ws<=range(4)) THEN; irg = 3
      ELSE IF (ws>range(4).AND.ws<=range(5)) THEN; irg = 4
      ELSE                                       ; irg = 5; END IF
      coeffH = (a_h(irg)+b_h(irg)*ws**p_h(irg)+c_h(irg)*(ws-8.)**2.)*1.e-3
      coeffE = (a_e(irg)+b_e(irg)*ws**p_e(irg)+c_e(irg)*(ws-8.)**2.)*1.e-3

      ! ... Define top & bottom wet layers........................
      k1s = k1z(l);    
      kms = kmz(l);

      ! ... Add non-penetrative components to surface cells
      esw = saturated_vapor_pressure(salp(k1s,l)+273.)
      Lv  = latent_heat_vaporization(salp(k1s,l)+273.)
      ! a. Long wave radiation (incoming LW as measured) 
      HeatSource(k1s,l)  = Qlw * (1.- Al) -                      &
                           EmWater * StephanBoltzman *           &
                           (salp(k1s,l)+273.) ** 4.
      ! b. Latent & Sensible heat fluxes
      HeatSource(k1s,l) = HeatSource(k1s,l) -                    &
      & rhoair*             Lv * coeffE*(esw-ea)*0.622/Pa * ws - & 
      & rhoair*SpecificHeatAir * coeffH*(salp(k1s,l)-Ta ) * ws

      ! ... Add penetrative components to water column & 
      !     express heat sources in temperature units ............
      DO k = k1s, kms
        HeatSource(k,l) = (HeatSource(k,l) +                     & 
          Qsw*QswFr(k,l))/((rhop(k,l)+1000.)*Cp)
        IF (Qsw<-100.) THEN; HeatSource(k,l) = 0.0E0; ENDIF
      END DO

    ENDDO

  CASE (10) ! Interpolation
  
    ! Calculate residuals from first Barnes pass at met station points (added 12/2010 SWA)
    CALL CalculateMetResiduals (Qswresid,Taresid,RHresid,Qlwresid,uresid,vresid)

    ! ... Loop over surface cells 
    DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

      l = id_column(liter)

      ! ... Map 2D-l into 3D(i,j) indexes ........................
      i = l2i(l); j = l2j(l);

      ! ... Interpolate met records to grid points (added 12/2010 SWA)
      CALL InterpMetData (i,j,Qswresid,Taresid,RHresid,Qlwresid,uresid,vresid, &
                        & Qsw,Ta,RH,Qlw,uair(l),vair(l))

      ! ... Wind speed over the water column
      ws  = SQRT (uair(l)**2.+vair(l)**2.);
        
      ! ... Define bulk aerodynamic coefficients for heat transfer under neutral conditions 
      !     based on Kondo, 1975. Air-sea bulk transfer coefficients in diabatic
      !     conditions. In Boundary-Layer Meteorology 9 (1975) 91-112
      ! ... Define range of wind speed 
      IF      (ws>range(1).AND.ws<=range(2)) THEN; irg = 1
      ELSE IF (ws>range(2).AND.ws<=range(3)) THEN; irg = 2
      ELSE IF (ws>range(3).AND.ws<=range(4)) THEN; irg = 3
      ELSE IF (ws>range(4).AND.ws<=range(5)) THEN; irg = 4
      ELSE                                       ; irg = 5; END IF
      coeffH = (a_h(irg)+b_h(irg)*ws**p_h(irg)+c_h(irg)*(ws-8.)**2.)*1.e-3
      coeffE = (a_e(irg)+b_e(irg)*ws**p_e(irg)+c_e(irg)*(ws-8.)**2.)*1.e-3

      ! ... Calculate drag coefficient (Amorocho & DeVries)  
      cdw(l) = 0.0015 * 1./(1.0+ EXP((12.5-ws)/1.56))+0.00104      

      ! ... Define variables used in cals. and equal to all surface cells
      ea    = saturated_vapor_pressure (Ta+273.) * Rh
      EmAir = 0.642 * ( ea / (Ta+273.) )**0.14285714 		
      EmAir = EmAir * ( 1. + 0.17 * Cc**2. ) 
      Qbri  = EmAir * StephanBoltzman * (Ta+273.) ** 4. 
      Qbra  = Qbri * ( 1. - Al ) 

      ! ... Define top & bottom wet layers........................
      k1s = k1z(l);    
      kms = kmz(l);

      ! ... Add non-penetrative components to surface cells
      esw = saturated_vapor_pressure(salp(k1s,l)+273.)
      Lv  = latent_heat_vaporization(salp(k1s,l)+273.)
      ! a. Long wave radiation 
      HeatSource(k1s,l)  = Qbra - EmWater * StephanBoltzman *    &
                           (salp(k1s,l)+273.) ** 4.
      ! b. Latent & Sensible heat fluxes
      HeatSource(k1s,l) = HeatSource(k1s,l) -                    &
      & rhoair*             Lv * coeffE*(esw-ea)*0.622/Pa * ws - & 
      & rhoair*SpecificHeatAir * coeffH*(salp(k1s,l)-Ta ) * ws

      ! ... Add penetrative components to water column & 
      !     express heat sources in temperature units ............
      DO k = k1s, kms
        HeatSource(k,l) = (HeatSource(k,l) +                     & 

          Qsw*QswFr(k,l))/((rhop(k,l)+1000.)*Cp)
        IF (Qsw<-100.) THEN; HeatSource(k,l) = 0.0E0; ENDIF
      END DO

    ENDDO
  
  CASE (11) ! Interpolation
  
    ! Calculate residuals from first Barnes pass at met station points (added 12/2010 SWA)
    CALL CalculateMetResiduals (Qswresid,Taresid,RHresid,Qlwresid,uresid,vresid)

    ! ... Interpolate other variables & compute heat fluxes
    DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

      l = id_column(liter)

      ! ... Map 2D-l into 3D(i,j) indexes ........................
      i = l2i(l); j = l2j(l);

      ! ... Interpolate met records to grid points (added 12/2010 SWA)
      CALL InterpMetData (i,j,Qswresid,Taresid,RHresid,Qlwresid,uresid,vresid, &
                        & Qsw,Ta,RH,Qlw,uair(l),vair(l))

      ! ... Wind speed over the water column
      ws  = SQRT (uair(l)**2.+vair(l)**2.);

      ! ... Define bulk aerodynamic coefficients for heat transfer under neutral conditions 
      !     based on Kondo, 1975. Air-sea bulk transfer coefficients in diabatic
      !     conditions. In Boundary-Layer Meteorology 9 (1975) 91-112
      ! ... Define range of wind speed 
      IF      (ws>range(1).AND.ws<=range(2)) THEN; irg = 1
      ELSE IF (ws>range(2).AND.ws<=range(3)) THEN; irg = 2
      ELSE IF (ws>range(3).AND.ws<=range(4)) THEN; irg = 3
      ELSE IF (ws>range(4).AND.ws<=range(5)) THEN; irg = 4
      ELSE                                       ; irg = 5; END IF
      coeffH = (a_h(irg)+b_h(irg)*ws**p_h(irg)+c_h(irg)*(ws-8.)**2.)*1.e-3
      coeffE = (a_e(irg)+b_e(irg)*ws**p_e(irg)+c_e(irg)*(ws-8.)**2.)*1.e-3

      ! ... Calculate drag coefficient (Amorocho & DeVries)  
      cdw(l) = 0.0015 * 1./(1.0+ EXP((12.5-ws)/1.56))+0.00104      

      ! ... Define variables used in cals. and equal to all surface cells
      Qbra  = Qlw * ( 1. - Al ) 

      ! ... Define top & bottom wet layers........................
      k1s = k1z(l);    
      kms = kmz(l);

      ! ... Add non-penetrative components to surface cells
      ea  = saturated_vapor_pressure(Ta         +273.) * Rh
      esw = saturated_vapor_pressure(salp(k1s,l)+273.)
      Lv  = latent_heat_vaporization(salp(k1s,l)+273.)
      ! a. Long wave radiation 
      HeatSource(k1s,l)  = Qbra - EmWater * StephanBoltzman *    &
                           (salp(k1s,l)+273.) ** 4.
      ! b. Latent & Sensible heat fluxes
      HeatSource(k1s,l) = HeatSource(k1s,l) -                    &
      & rhoair*             Lv * coeffE*(esw-ea)*0.622/Pa * ws - & 
      & rhoair*SpecificHeatAir * coeffH*(salp(k1s,l)-Ta ) * ws

      ! ... Add penetrative components to water column & 
      !     express heat sources in temperature units ............
      DO k = k1s, kms
        HeatSource(k,l) = (HeatSource(k,l) +                     & 
          Qsw*QswFr(k,l))/((rhop(k,l)+1000.)*Cp)
        IF (Qsw<-100.) THEN; HeatSource(k,l) = 0.0E0; ENDIF
      END DO

    ENDDO
  
  END SELECT

END SUBROUTINE DistributeMomentumHeatSourcesH

!************************************************************************
REAL FUNCTION saturated_vapor_pressure (T)
!************************************************************************

  IMPLICIT NONE
  REAL, INTENT (IN) :: T

  saturated_vapor_pressure = 2.1718e10 * EXP( -4157. / (T - 33.91) ) 

  END FUNCTION saturated_vapor_pressure 

!************************************************************************
REAL FUNCTION latent_heat_vaporization (T)
!************************************************************************

  IMPLICIT NONE
  REAL, INTENT (IN) :: T

  latent_heat_vaporization = 1.91846e6 * ( T / (T - 33.91) ) ** 2. 

  END FUNCTION latent_heat_vaporization

!***********************************************************************
SUBROUTINE InitializeInterpolationMethods
!***********************************************************************
!
!  Purpose: Solves for the weighting parameters used in the interpolation scheme
!  Modified 12/2010 by SWA to include 3 types of interpolation
!-----------------------------------------------------------------------

  ! ... Local variables
  INTEGER :: i, j, imet, istat, counter1, jmet
  REAL    :: xi, yi, di, da, sumw
  REAL    :: x1,y1,x2,y2,dist1(nmetstat-1),kappa1

  SELECT CASE (iinterp)


  CASE (1) ! Original Barnes interpolation code
    DO i = i1, im; DO j = j1, jm

      ! Find average distance from grip point to all met stations
      da = 0.0
      DO imet = 1, nmetstat
        xi = metxy((imet-1)*2+1) - FLOAT(i);  
        yi = metxy((imet-1)*2+2) - FLOAT(j);
        di = SQRT(xi**2.+yi**2.)
        da = da + di
      ENDDO
      da = da / nmetstat

      ! Estimate weight based on distance to met station imet & da
      sumw = 0.0E0
      DO imet = 1, nmetstat
        xi = metxy((imet-1)*2+1) - FLOAT(i);  
        yi = metxy((imet-1)*2+2) - FLOAT(j);
        di = SQRT(xi**2.+yi**2.)
        weightst(i,j,imet) = EXP(-4.60517018598809 * (di**2.) / (da**2.))
        sumw = sumw + weightst(i,j,imet)
      ENDDO
    
      ! Normalize the weight
      DO imet = 1, nmetstat
        weightst(i,j,imet)=weightst(i,j,imet)/sumw
      ENDDO

    ENDDO; ENDDO


  CASE (2) ! New Barnes interpolation, with adjustable parameters gammaB, delNfactor
  
    ! Find average distance of a met station to closest other met station
    ! Distance is in grid units
    da=0.0
    dist1(:)=1.0E5;
    DO imet = 1, nmetstat
      x1 = metxy((imet-1)*2+1)
      y1 = metxy((imet-1)*2+2)
      counter1=1
      DO i = 1, nmetstat
        IF (imet /= i) THEN
          x2 = metxy((i-1)*2+1)
          y2 = metxy((i-1)*2+2)
          dist1(counter1) = SQRT((x1-x2)**2.0+(y1-y2)**2.0)
          counter1=counter1+1
        ENDIF
      ENDDO
      da = da + MINVAL(dist1)
    ENDDO
    da = da / REAL(nmetstat)

    ! Calculate smoothing scale length (from Koch 1983)
    IF (gammaB<=0.2) THEN ! 5.0515 factor calculated for gammaB=0.2
       kappa1 = 5.0515*4.0/pi**2.0 * (delNfactor*da)**2.0
    ELSEIF (gammaB<=0.3) THEN ! 3.5132 factor calculated for gammaB=0.3
       kappa1 = 3.5132*4.0/pi**2.0 * (delNfactor*da)**2.0
    ELSEIF (gammaB<=0.5) THEN
       kappa1 = 2.3838*4.0/pi**2.0 * (delNfactor*da)**2.0
    ELSEIF (gammaB<=0.6) THEN
       kappa1 = 2.1153*4.0/pi**2.0 * (delNfactor*da)**2.0
    ELSEIF (gammaB>0.6) THEN ! 1.7826 factor calculated for gammaB=0.8
       kappa1 = 1.7826*4.0/pi**2.0 * (delNfactor*da)**2.0
    ENDIF
    
    DO i = i1, im; DO j = j1, jm

    ! Estimate weight based on distance to met station imet and kappa1
    ! NOTE: no weight normalization for this case
    DO imet = 1, nmetstat
      xi = metxy((imet-1)*2+1) - FLOAT(i);
      yi = metxy((imet-1)*2+2) - FLOAT(j);
      di = xi**2.+yi**2.
      weightst(i,j,imet) = EXP(-di / kappa1)
    ENDDO

    ENDDO; ENDDO
    ! Estimate weight based on distance to met station points MAC
    DO imet = 1, nmetstat
      i = NINT(metxy((imet-1)*2+1))
      j = NINT(metxy((imet-1)*2+2))
      DO jmet = 1, nmetstat
	      
	      xi = metxy((jmet-1)*2+1) - FLOAT(i);
	      yi = metxy((jmet-1)*2+2) - FLOAT(j);
	      di = xi**2.+yi**2.
	      weightn(imet,jmet) = EXP(-di / kappa1)
    ENDDO; ENDDO
    

  CASE (3) ! Inverse distance interpolation

    DO i = i1, im; DO j = j1, jm

    ! Estimate weight based on inverse distance to met station
    sumw = 0.0E0
    DO imet = 1, nmetstat
      xi = metxy((imet-1)*2+1) - FLOAT(i);
      yi = metxy((imet-1)*2+2) - FLOAT(j);
      di = xi**2.+yi**2.
      IF (di > 0.0) THEN
        weightst(i,j,imet) = 1.0/di
      ELSE
        weightst(i,j,imet) = 1.0
      ENDIF
      sumw = sumw + weightst(i,j,imet)
    ENDDO
    
    ! Normalize the weight
    DO imet = 1, nmetstat
      weightst(i,j,imet)=weightst(i,j,imet)/sumw
    ENDDO

    ENDDO; ENDDO
    
  
  END SELECT

END SUBROUTINE InitializeInterpolationMethods

!***********************************************************************
SUBROUTINE CalculateMetResiduals (Qswresid,Taresid,RHresid,Qlwresid,uresid,vresid)
!***********************************************************************
!
!  Purpose: Calculates residuals from interpolation of met parameters
!           at the met station points
!  Created 12/2010 by SWA
!-----------------------------------------------------------------------

  ! Local variables
  INTEGER :: k,i,j,imet
  REAL, DIMENSION(nmetstat), INTENT (OUT) :: Qswresid,Taresid,RHresid,Qlwresid
  REAL, DIMENSION(nmetstat), INTENT (OUT) :: uresid,vresid
  REAL :: Qswdum,Tadum,RHdum,Qlwdum,uairdum,vairdum,sumw

  ! Initialize
  Qswresid(:)=0.0;  Taresid(:)=0.0;  RHresid(:)=0.0
  Qlwresid(:)=0.0;  uresid(:)=0.0;   vresid(:)=0.0
  IF (iinterp == 2) THEN
    
    ! Calculate residuals from first Barnes pass at met station points
    DO k = 1, nmetstat
      i = NINT(metxy((k-1)*2+1))
      j = NINT(metxy((k-1)*2+2))
      sumw=SUM(Weightn(k,:))   !USE weightn MAC
      Qswdum = 0.0; Tadum = 0.0; RHdum = 0.0; Qlwdum = 0.0;
      uairdum=0.0; vairdum=0.0
      DO imet = 1, nmetstat
        Qswdum  = Qswdum  + Weightn(k,imet) / sumw * Qsw2D (imet) !USE weightn MAC
        Tadum   = Tadum   + Weightn(k,imet) / sumw * Ta2D  (imet) !USE weightn MAC
        RHdum   = RHdum   + Weightn(k,imet) / sumw * RH2D  (imet) !USE weightn MAC
        Qlwdum  = Qlwdum  + Weightn(k,imet) / sumw * Qlw2D (imet) !USE weightn MAC
        uairdum = uairdum + Weightn(k,imet) / sumw * uair2D(imet) !USE weightn MAC
        vairdum = vairdum + Weightn(k,imet) / sumw * vair2D(imet) !USE weightn MAC
      ENDDO
      Qswresid(k) = Qsw2D(k)  - Qswdum
      Taresid(k)  = Ta2D(k)   - Tadum
      RHresid(k)  = RH2D(k)   - RHdum
      Qlwresid(k) = Qlw2D(k)  - Qlwdum
      uresid(k)   = uair2D(k) - uairdum
      vresid(k)   = vair2D(k) - vairdum
    ENDDO
  ELSE ! Residuals not used in calculations
    Qswresid(:) = -9999.0
    Taresid(:)  = -9999.0
    RHresid(:)  = -9999.0
    Qlwresid(:) = -9999.0
    uresid(:)   = -9999.0
    vresid(:)   = -9999.0
  ENDIF

END SUBROUTINE CalculateMetResiduals

!***********************************************************************
SUBROUTINE InterpMetData (iin,jin,Qswresid,Taresid,RHresid,Qlwresid,uresid,vresid,Qsw,Ta,RH,Qlw,uair,vair)
!***********************************************************************
!
!  Purpose: Interpolated met data at cell i,j
!  Created 12/2010 by SWA
!-----------------------------------------------------------------------

  INTEGER, INTENT (IN) :: iin,jin
  REAL, DIMENSION(nmetstat), INTENT (IN) :: Qswresid,Taresid,RHresid,Qlwresid
  REAL, DIMENSION(nmetstat), INTENT (IN) :: uresid,vresid
  REAL, INTENT (OUT) :: Qsw,Ta,RH,Qlw,uair,vair
  ! Local variables
  INTEGER :: imet
  REAL :: sumw

  ! Initialize
  Qsw = 0.0; Ta = 0.0; RH = 0.0; Qlw = 0.0; 
  uair = 0.0E0; vair = 0.0E0;
  ! Interpolation
  DO imet = 1, nmetstat
    Qsw  = Qsw    + WeightSt(iin,jin,imet) * Qsw2D (imet)
    Ta   = Ta     + WeightSt(iin,jin,imet) * Ta2D  (imet)
    RH   = RH     + WeightSt(iin,jin,imet) * RH2D  (imet)
    Qlw  = Qlw    + WeightSt(iin,jin,imet) * Qlw2D (imet)
    uair = uair + WeightSt(iin,jin,imet) * uair2D(imet)
    vair = vair + WeightSt(iin,jin,imet) * vair2D(imet)
  ENDDO
  IF (iinterp == 2) THEN ! Modify for new Barnes interpolation
    sumw=SUM(Weightst(iin,jin,:)) ! These weights weren't initially normalized
    Qsw  = Qsw  / sumw
    RH   =  RH  / sumw
    Ta   =  Ta  / sumw
    Qlw  = Qlw  / sumw
    uair = uair / sumw
    vair = vair / sumw
    ! Second pass Barnes interpolation
    sumw=sum(WeightSt(iin,jin,:)**(1.0/gammaB))
    DO imet = 1, nmetstat
      Qsw  = Qsw    + WeightSt(iin,jin,imet)**(1.0/gammaB)/sumw * Qswresid (imet)
      Ta   = Ta     + WeightSt(iin,jin,imet)**(1.0/gammaB)/sumw * Taresid  (imet)
      RH   = RH     + WeightSt(iin,jin,imet)**(1.0/gammaB)/sumw * RHresid  (imet)
      Qlw  = Qlw    + WeightSt(iin,jin,imet)**(1.0/gammaB)/sumw * Qlwresid (imet)
      uair = uair + WeightSt(iin,jin,imet)**(1.0/gammaB)/sumw * uresid(imet)
      vair = vair + WeightSt(iin,jin,imet)**(1.0/gammaB)/sumw * vresid(imet)
    ENDDO
  ENDIF

END SUBROUTINE InterpMetData

END MODULE si3d_BoundaryConditions
