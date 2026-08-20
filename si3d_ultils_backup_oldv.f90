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

  INTEGER :: nn, inn, i, j, k, kk, l, k1s, kms, nwl, ksrc, plmdim, itr,innH
  INTEGER :: NITERMAX, kint, e1, jct, FLAG
  INTEGER, SAVE :: NITERPLUME, NLI, NLO, ktop
  REAL    :: areatot, Tsource, Rsource, sumQpss, qdenom
  REAL(8) :: wselev, dfelev, dflgth, hcell, rjulday, lambnot, diamm, linot
  REAL(8) :: elevt, qwt, tplumt, comgpt, tplum0, comgp0, qscfm, frconot
  REAL(8) :: ERRORDP, ERRORDPA, toler, depthed, salplut, qwed, QpssED
  REAL(8) :: tplumti,comgpti,salpluti,tplumed,salplued,comgped
  REAL(8) :: alphaii,alphaaa,alphaoo,gammapp,froudeii,froudeoo,lambdaa
  real    :: rho_amb, rho_source, depth
  REAL(8), DIMENSION (1:km1) :: zamb, Tamb, DOamb, UA, VA ! B.C. for plume model
  REAL(8), DIMENSION (1:km1) :: qwd, qwdi, qwdo               ! Outflow rate for plume
  REAL(8), DIMENSION (1:km1) :: bwd, bwdi, bwdo                 ! Radius of plume
  REAL(8), DIMENSION (1:km1) :: lwdi, lwdo                 ! Radius of plume
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
        ! Inicializamos Qpss, Tpss, Rpss !cintia cambio

        Qpss(:,inn)  =0.0 ;
        Tpss(:,inn)  =0.0 ;
        Rpss(:,inn,:)=0.0 ;

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
      ! ... Update forcing variables before deciding if diffuser is OFF
      CALL PointSourceSinkForcing (nn,thrs)
      IF (ABS(flpss(nn)) < qthrs(nn)) THEN  ! Diffuser OFF
        DiffON = .FALSE.  
        
        DO innH = 1, iopssH(omp_get_thread_num ( )+1)
          inn = ioph2iop(innH,omp_get_thread_num ( )+1)  
          IF (iodev(inn) .NE. nn) CYCLE 
          Qpss (inn,:) = 0.0E0
          kdetr(inn  ) = km1
        ENDDO
      ELSE                                  ! Update diffuser FLOWS          
        IF & 
        ( (DiffON == .FALSE.)       .OR.  & ! Diffuser is TURNED ON
        (n == 1)                  .OR.   &
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

            ! ... Define ambient temperatures
            Tamb(k1s:kms) = salpp(k1s:kms,l)
            Tamb(kms+1  ) = Tamb(kms);
            Tamb(1      ) = Tamb(k1s);

            ! ... Define DO concentrations
            DOamb(k1s:kms) = tracerpp(k1s:kms,l,LDO) / 1000.0 ! Concentrations are in mg/m3 and need to be g/m3 for plume model
            DOamb(kms+1  ) = DOamb(kms);
            DOamb(1      ) = DOamb(k1s);

            ! ... Define ambient U water velocity
            UA(k1s:kms) = upp(k1s:kms,l)
            UA(kms+1  ) = UA(kms);
            UA(1      ) = UA(k1s);

            ! ... Define ambient V water velocity
            VA(k1s:kms) = vpp(k1s:kms,l)
            VA(kms+1  ) = VA(kms);
            VA(1      ) = VA(k1s);

            ! ... Depths for cells in plume column from datum
            zamb(k1s  ) = hp(k1s,l)/2.
            DO k = k1s+1, kms
              zamb(k) = zamb(k-1) + (hp(k-1,l)+hp(k,l))/2.
            END DO
            zamb(kms+1) =  zamb(kms)+hp(kms,l)
            zamb(1    ) = -zamb(k1s)

            ! ... Inputs for plume model
            dfLgth  = real(dfL(nn),8)         ;       ! Length of diffuser
            rjulday = doy             ;       ! Julian day (arbitrary)
            wselev  = 0.0000          ;       ! Elevation of free surface
            ksrc    = kms-1           ;       ! Layer No. where diffuser is located
            dfelev  = -zamb(ksrc)     ;       ! Elevation of diffuser
            hcell   = real(ddz,8)             ;       ! Pressumed constant - thickess of cells
            qwd     = 0.0E0           ;       ! Initialize qwd
            bwd     = 0.0E0                   ! Initialize perimeter (FJRplume)
            bwdi    = 0.0E0           ! Initialize perimeter (JCT)
            bwdo    = 0.0E0           ! Initialize perimeter (JCT)
            lwdi    = 0.0E0           ! Initialize perimeter (JCT)
            lwdo    = 0.0E0           ! Initialize perimeter (JCT) 
            qscfm   = real(flpss(nn),8)       ;       ! Air flow rate
            frconot = 1.00            ;       ! Fraction of O2 in air (not used?)
            lambnot = lambdanot(nn)   ;       ! Half-width
            !linot   = lnot(nn)        ;       ! Diffuser length
            linot   = dfLgth        ;       ! Diffuser length
            ! PRINT*,'linot ',linot
            diamm   = diammb(nn)      ;       ! Initial bubble diameter
            alphaii  = alphai(nn)          ! Entrainment coefficient inner plume (-)
            alphaoo  = alphao(nn)          ! Entrainment coefficient outer plume (-)
            alphaaa  = alphaa(nn)          ! Entrainment coefficient from ambient (-)
            lambdaa  = lambda(nn)          ! Fraccion of plume occupied by bubble
            froudeii = froudei(nn)         ! Froude number inner plume
            froudeoo = froudeo(nn)         ! Froude number outer plume
            gammapp  = gammap(nn)          !

            IF (ptype(nn) < 2) THEN
              plmdim = 1 ! Linear Plume


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

            ! ... CIRCULAR PLUME
            ELSEIF (ptype(nn) ==2) THEN 
              plmdim = 2 ! Circular Plume


              ! ... Run plume model
              CALL CIRCULAR_PLUME(iyr,rjulday,wselev,dfelev,kms,    &
                                  lambnot,salamb,patm,diamm,        &
                                  qscfm,frconot,                    &
                                  ksrc, hcell,                      &
                                  zamb (1:kms),                     &
                                  Tamb (1:kms),                     &
                                  DOamb (1:kms),                    &
                                  UA(1:kms),                        &
                                  VA(1:kms),                        &
                                  elevt, qwt, tplumt, comgpt,       &
                                  qwd(1:ksrc), bwd(1:ksrc), ktop)

              PRINT*, '******************************************************'
              PRINT*, '----OUTPUT FROM PLUME ROUTINES (CIRCULAR) ------------'
              PRINT*, '******************************************************'
              PRINT*, 'Results of Plume Model elev',elevt,dfelev
              PRINT*, 'Results of Plume Model qwt ',qwt,ktop
              PRINT*, 'Results of Plume Model T&O ',tplumt,comgpt
              PRINT*, 'Results of Plume Model pwd ',ksrc, bwdi(ksrc),bwdi(ktop+1)
              PRINT*, '******************************************************'


            ! ... CIRCULAR DOUBLE PLUME
            ELSEIF (ptype(nn) == 3 .OR. ptype(nn)==4 ) THEN ! Double plume - circular

              ! ... Run plume model
              NITERMAX = 10
              NITERPLUME = 0
              ERRORDPA = 100.0
              !             innerplumeold = 0
              innerplume = 0
              outerplume = 0

              toler = 10         ! in            

              DO WHILE (NITERPLUME .LT. 3)
                NITERPLUME = NITERPLUME+1
                PRINT*, 'NITERPLUME', NITERPLUME

                !... Run plume model
                CALL INNER_PLUME(iyr,rjulday,wselev,dfelev,kms,       &
                                lambnot,salamb,patm,diamm,           &
                                qscfm, frconot, ksrc, hcell,         &
                                zamb (1:kms),                        &
                                Tamb (1:kms),                        &
                                DOamb(1:kms),                        &
                                UA(1:kms),                           &
                                VA(1:kms),                           &
                                elevt, qwt, tplumt, comgpt, salplut, &
                                qwdi(1:ksrc),                        &
                                bwdi(1:ksrc),                        &
                                ktop,kint,depthed, NLI, NLO,         &
                                NITERPLUME,                          &
                                innerplume,                          &
                                outerplume,                          &
                                alphaii, alphaoo, alphaaa, lambdaa,      &
                                froudeii, froudeoo, gammapp)

                PRINT*, '******************************************************'
                PRINT*, '----OUTPUT FROM PLUME ROUTINES (CIRCULAR-inner) ------'
                PRINT*, '******************************************************'
                PRINT*, 'Results of Plume Model elev',elevt,dfelev
                PRINT*, 'Results of Plume Model qwt ',qwt,ktop,kint
                PRINT*, 'Results of Plume Model T&O ',tplumt,comgpt
                PRINT*, 'Results of Plume Model pwd ',ksrc, bwdi(ksrc),bwdi(ktop+1)
                PRINT*, '******************************************************'


                ! ... Save information from the inner plume to the outer plume
                tplumti  = tplumt
                comgpti  = comgpt
                salpluti = salplut

                ! ... Define flow at cells
                Qpss(:,inn) = 0.0

                ! ... Define flow at entrainment cells
                DO k = ktop+1,ksrc
                  Qpss(k,inn) = -real(qwdi(k))*dy/real(dfLgth)
                ENDDO

                ! ... Define flow at detrainment cell to force volume conservation
                Qpss(ktop,inn) = -SUM(Qpss(ktop+1:ksrc,inn))

                CALL OUTER_PLUME(iyr,rjulday,wselev,dfelev,kms,          &
                                bwdi(ktop+1),salamb,patm,            &
                                qwt, frconot, ksrc, hcell,           &
                                zamb (1:kms),                        &
                                Tamb (1:kms),                        &
                                DOamb(1:kms),                        &
                                UA(1:kms),                           &
                                VA(1:kms),                           &
                                elevt, qwed,                         &
                                qwdo(1:ksrc),                        &
                                bwdo(1:ksrc),                        &
                                ktop,kint, tplumed,salplued,comgped, &
                                depthed, NLI, NLO,                   &
                                NITERPLUME,                          &
                                innerplume,                          &
                                outerplume,                          &
                                alphaii, alphaoo, alphaaa, lambdaa,      &
                                froudeii, froudeoo, gammapp)


                PRINT*, '******************************************************'
                PRINT*, '----OUTPUT FROM PLUME ROUTINES (CIRCULAR-outer) ------'
                PRINT*, '******************************************************'
                PRINT*, 'Results of Plume Model elev',depthed, elevt
                PRINT*, 'Results of Plume Model qwt ',qwed,ktop,kint
                PRINT*, 'Results of Plume Model T&O ',tplumed,comgped
                PRINT*, 'Results of Plume Model pwd ',ksrc, bwdo(ktop),bwdo(kint+1),bwdo(kint)
                PRINT*, '******************************************************'
                PRINT*, 'NLI,NLO', NLI,NLO
              ENDDO
              ! ... Define flow at entrainment cells
              qwd = 0.0
              qwd(ktop:kint-1)  = -qwdo(ktop:kint-1) ! JCT_2022
              qwd(kint:ksrc)    = qwdi(kint:ksrc)

              ! DO k = ktop+1, ksrc
              ! OPEN (UNIT = 55, FILE="check_outter.txt", POSITION="APPEND")
              ! WRITE(UNIT = 55, FMT = '(4F12.6)') zamb(k),qwdi(k),qwdo(k),qwd(k)
              ! CLOSE(UNIT = 55)
              ! ENDDO


              ! ... DOUBLE PLUME RECTANGULAR:
            ELSE

              ! ... Run plume model 
              NITERMAX = 10
              NITERPLUME = 0
              ERRORDPA = 100.0
              innerplume = 0
              outerplume = 0
              toler = 10         ! in

              DO WHILE (NITERPLUME .LT. 3)
                NITERPLUME = NITERPLUME+1
                ! PRINT*, 'NITERPLUME', NITERPLUME

                !... Run plume model
                ! PRINT*,'lambnot',lambnot, lambda
                CALL INNER_PLUME_RECT(iyr,rjulday,wselev,dfelev,kms,  &
                                      lambnot,salamb,patm,diamm,           &
                                      qscfm, frconot, ksrc, hcell,         &
                                      zamb (1:kms),                        &
                                      Tamb (1:kms),                        &
                                      DOamb(1:kms),                        &
                                      UA(1:kms),                           &
                                      VA(1:kms),                           &
                                      elevt, qwt, tplumt, comgpt, salplut, &
                                      qwdi(1:ksrc),                        &
                                      bwdi(1:ksrc),                        &
                                      ktop,kint,depthed, NLI, NLO,         &
                                      NITERPLUME,                          &
                                      innerplume,                          &
                                      outerplume,                          &
                                      alphaii, alphaoo, alphaaa, lambdaa,  &
                                      froudeii, froudeoo, gammapp,linot,   &
                                      lwdi(1:ksrc))

                ! PRINT*, '******************************************************'
                ! PRINT*, '----OUTPUT FROM PLUME ROUTINES (RECTANGULAR-inner) ---'
                ! PRINT*, '******************************************************'
                ! PRINT*, 'Results of Plume Model elev',elevt,dfelev
                ! PRINT*, 'Results of Plume Model qwt ',qwt,ktop,kint
                ! PRINT*, 'Results of Plume Model T&O ',tplumt,comgpt
                ! PRINT*, 'Results of Plume Model pwd ',ksrc, bwdi(ksrc),bwdi(ktop+1)
                ! PRINT*, '******************************************************'

                ! DO k = ktop+1, ksrc
                !   OPEN (UNIT = 54, FILE="check_inner.txt", POSITION="APPEND")
                !   WRITE(UNIT = 54, FMT = '(3F12.6)') zamb(k),qwdi(k),bwdi(k)
                !   CLOSE(UNIT = 54)
                ! ENDDO

                ! ... Save information from the inner plume to the outer plume
                tplumti  = tplumt
                comgpti  = comgpt
                salpluti = salplut

                ! ... Define flow at cells
                Qpss(:,inn) = 0.0

                ! ... Define flow at entrainment cells
                DO k = ktop+1,ksrc
                  Qpss(k,inn) = -(real(qwdi(k)))*dy/real(dfLgth)
                ENDDO

                ! ... Define flow at detrainment cell to force volume conservation
                Qpss(ktop,inn) = -SUM(Qpss(ktop+1:ksrc,inn))

                ! CALL OUTER_PLUME_RECT(iyr,rjulday,wselev,dfelev,kms,     &
                ! lambnot,linot, bwdi(ktop+1),         &
                ! salamb,patm,                         &
                ! qwt, frconot, ksrc, hcell,           &
                ! zamb (1:kms),                        &
                ! Tamb (1:kms),                        &
                ! DOamb(1:kms),                        &
                ! UA(1:kms),                           &
                ! VA(1:kms),                           &
                ! elevt, qwed,                         &
                ! qwdo(1:ksrc),                        &
                ! bwdo(1:ksrc),                        &
                ! ktop,kint, tplumed,salplued,comgped, &
                ! depthed, NLI, NLO,                   &
                ! NITERPLUME,                          &
                ! innerplume,                          &
                ! outerplume,                          &
                ! alphaii, alphaoo, alphaaa, lambdaa,  &
                ! froudeii, froudeoo, gammapp,         &
                ! lwdo(1:ksrc))


                CALL OUTER_PLUME_RECT2(iyr,rjulday,wselev,dfelev,kms,    &
                                      bwdi(ktop+1),salamb,patm,            &
                                      qwt, frconot, ksrc, hcell,           &
                                      zamb (1:kms),                        &
                                      Tamb (1:kms),                        &
                                      DOamb(1:kms),                        &
                                      UA(1:kms),                           &
                                      VA(1:kms),                           &
                                      elevt, qwed,                         &
                                      qwdo(1:ksrc),                        &
                                      bwdo(1:ksrc),                        &
                                      ktop,kint, tplumed,salplued,comgped, &
                                      depthed, NLI, NLO,                   &
                                      NITERPLUME,                          &
                                      innerplume,                          &
                                      outerplume,                          &
                                      alphaii, alphaoo, alphaaa, lambdaa,      &
                                      froudeii, froudeoo, gammapp,         &
                                      lwdo(1:ksrc))

                ! CALL OUTER_PLUME(iyr,rjulday,wselev,dfelev,kms,          &
                ! bwdi(ktop+1),salamb,patm,            &
                ! qwt, frconot, ksrc, hcell,           &
                ! zamb (1:kms),                        &
                ! Tamb (1:kms),                        &
                ! DOamb(1:kms),                        &
                ! UA(1:kms),                           &
                ! VA(1:kms),                           &
                ! elevt, qwed,                         &
                ! qwdo(1:ksrc),                        &
                ! bwdo(1:ksrc),                        &
                ! ktop,kint, tplumed,salplued,comgped, &
                ! depthed, NLI, NLO,                   &
                ! NITERPLUME,                          &
                ! innerplume,                          &
                ! outerplume,                          &
                ! alphaii, alphaoo, alphaaa, lambdaa,      &
                ! froudeii, froudeoo, gammapp)
                ! PRINT*, '******************************************************'
                ! PRINT*, '----OUTPUT FROM PLUME ROUTINES (RECTANGULAR-outer) ---'
                ! PRINT*, '******************************************************'
                ! PRINT*, 'Results of Plume Model elev',depthed, elevt
                ! PRINT*, 'Results of Plume Model qwt ',qwed,ktop,kint
                ! PRINT*, 'Results of Plume Model T&O ',tplumed,comgped
                ! PRINT*, 'Results of Plume Model pwd ',ksrc, bwdo(ktop),bwdo(kint+1),bwdo(kint)
                ! PRINT*, '******************************************************'

                ! DO k = ktop+1, ksrc
                !   OPEN (UNIT = 55, FILE="check_outter.txt", POSITION="APPEND")
                !   WRITE(UNIT = 55, FMT = '(3F12.6)') zamb(k),qwdo(k),bwdo(k)
                !   CLOSE(UNIT = 55)
                ! ENDDO


              ENDDO

              ! ... Define flow at entrainment cells
              qwd = 0.0
              qwd(ktop:kint-1)  = qwdo(ktop:kint-1) ! JCT_2020
              qwd(kint:ksrc)    = qwdi(kint:ksrc)


            ENDIF  ! End plume type 

            ! ... ******* (The rest is common to all CASE(1:)) *******************
            ! ... Define flow at cells
            Qpss(:,inn) = 0.0

            IF (ptype(nn)<3) THEN
              kdetr(inn) = ktop
              ! ... Define flow at entrainment cells
              DO kk = ktop+1, ksrc
                Qpss(kk,inn) = real(qwd(kk))*dy/real(dfLgth)
              ENDDO
            ELSE
              ! ... Define flow at entrainment cells
              DO kk = ktop+1,ksrc
                Qpss(kk,inn) = -real(qwd(kk))*dy/real(dfLgth)
              ENDDO
              IF (ptype(nn) ==3 .OR. ptype(nn)==5) THEN ! detrainment at end of outter plume JCT_2020
                kdetr(inn) = kint -1
                !! ... Define flow at entrainment cells
                ! DO kk = ktop +1, ksrc
                !    Qpss(kk,inn)=-qwd(kk)
                ! ENDDO
                !! ... Define flow at detrainment cell to force volume conservation
                !    Qpss(kdetr(inn),inn) = -SUM(Qpss(ktop+1:ksrc,inn)) + Qpss(kdetr(inn),inn) ! JCT_2017

              ELSEIF (ptype(nn) ==4 .OR. ptype(nn)==6) THEN  ! detrainment at equilibrium depth JCT_2020
                ! PRINT*, 'flag_plume7'
                Tsource = 0.0 ! FJRPlumes
                DO kk = ktop+1,kms
                  Tsource = Tsource + salpp(kk,l)*Qpss(kk,inn)
                ENDDO
                sumQpss = SUM(Qpss(ktop+1:kms,inn))
                IF (ABS(sumQpss) > 1.0E-12) THEN
                  Tsource = Tsource / sumQpss
                ELSE
                  Tsource = salpp(ktop,l)
                END IF
                FLAG = 0
                kdetr(inn) = MAX(ktop, MIN(kms, kint-1)) ! fallback if no equilibrium depth is found
                DO kk=1,kms
                  if (zlevel(kk) == -100) then
                    depth = 0.5 * hp(kk,l)
                  else
                    depth = zlevel(kk) + 0.5 * hp(kk,l)
                  end if
                  rho_amb = densty_s(real(Tamb(kk)), 0.00004, depth)
                  rho_source = densty_s(real(Tsource), 0.00004, depth)
                  IF ((rho_amb .ge. rho_source) .AND. (FLAG .EQ. 0)) THEN
                    kdetr(inn) = kk ! JCT_2020
                    FLAG = 1
                  ENDIF
                ENDDO
                IF (FLAG .EQ. 0) THEN
                  PRINT*, 'WARNING: no equilibrium depth found for outer plume, fallback kdetr=', kdetr(inn), 'ktop=', ktop, 'kint=', kint
                ENDIF
                ! PRINT*, 'Tsource_fin', Tsource, kdetr(inn),Tamb(MAX(1,kdetr(inn)-1)),Tamb(kdetr(inn)),Tamb(MIN(kms,kdetr(inn)+1)),Tamb(33),Tamb(34),Tamb(35)
              ENDIF

              !! ... Define flow at entrainment cells .-DEFINES ABOVE ALREADY FJR 2021 01 24
              !DO kk = ktop+1,ksrc
              !Qpss(kk,inn) = -qwd(kk)
              !ENDDO

              ! ... Define flow at detrainment cell to force volume conservation - ACC June 2026
              qdenom = 0.0
              IF (ksrc >= ktop+1) THEN
                qdenom = SUM(Qpss(ktop+1:ksrc,inn))
              END IF
              IF (kdetr(inn) < ktop) kdetr(inn) = ktop
              IF (kdetr(inn) > ksrc) kdetr(inn) = ksrc
              IF (ABS(qdenom) > 1.0E-12) THEN
                Qpss(kdetr(inn),inn) = -qdenom + Qpss(kdetr(inn),inn)
              ELSE
                Qpss(kdetr(inn),inn) = 0.0
              END IF

              ! print*, 'flag_plume8'

              ! OPEN (UNIT=56, FILE="doubleplume.txt", POSITION="APPEND")
              ! WRITE(UNIT=56, FMT = '(7F12.6)') elevt,kdetr(inn),-zamb(kdetr(inn)),Qpss(kdetr(inn),inn),Tsource,Tamb(kdetr(inn)),-zamb(kint)
              ! CLOSE(UNIT=56)

              ! PRINT*, 'kk', kk, 'ktop',ktop,'ksrc',ksrc,'inn',inn,'kdetr(inn)',kdetr(inn),'Tamb(kdetr(inn))',Tamb(kdetr(inn))

              ! print*, 'flag_plume9'

              ! OPEN (UNIT=57, FILE="check_plumes.txt", POSITION="APPEND")
              ! DO kk = 1,ksrc
                
              !   !WRITE (UNIT=56, FMT='(3I3,8F8.2)') kk,ktop,ksrc,-zamb(kk),Qpss(kk,inn),salpp(kk,l),-zamb(kdetr(inn)),Qpss(kdetr(inn),inn),Tsource,Tamb(kdetr(inn)),-zamb(kint)
              !   WRITE (UNIT=57, FMT='(4I3,14F10.2)') kk,ktop,ksrc,kdetr(inn),-zamb(kk),Qpss(kk,inn),qwdi(kk),qwdo(kk),salpp(kk,l),-zamb(kdetr(inn)),Qpss(kdetr(inn),inn),Tsource,Tamb(kdetr(inn)),-zamb(kint),bwdi(kk),bwdo(kk),lwdi(kk),lwdo(kk)
              
              ! ENDDO
              ! CLOSE (UNIT=57)

              ! print*, 'flag_plume10'

            ENDIF
          ENDDO

        ENDIF
      ENDIF


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

        IF (ptype(nn) < 3) THEN
          k = ktop;
        ELSE
          k = kdetr(inn);
        ENDIF

        Tsource = 0.0 ! FJRPlumes
        IF (k > k1s) THEN  
          DO kk = ktop+1,kms
            Tsource  = Tsource + salpp(kk,l)*Qpss(kk,inn)
          ENDDO
          Tsource  = Tsource - salpp(k,l) * Qpss(k,inn)
          qdenom  = SUM(Qpss(ktop+1:kms,inn)) - Qpss(k,inn)
          IF (ABS(qdenom) > 1.0E-12) THEN
            Tsource  = Tsource / qdenom
          ELSE
            Tsource  = salpp(k,l)
          ENDIF
        ELSE
          Tsource = salpp(k,l)
        ENDIF
        Tpss(k,inn) = Tsource  ! Tpss conection plume <-> 3D
        !PRINT *, 'FJR junk', ktop, Tsource
        ! Rpss(:, inn, :) = 0.0
        IF (ntr > 0) THEN
          DO itr = 1, ntr
            ! if (itr .ne. LDO) then
              ! print*, 'SV - itr != LDO, Rsource and Rpss not calculated for this tracer'
              ! Rpss(:, inn, itr) = tracerpp(:, l, itr)
            ! else
              Rpss(:,inn,itr) = tracerpp(:,l,itr) / 1000.0 ! Concentrations are in mg/m3 and need to be g/m3 for plume model
              IF (ptype(nn) < 3) THEN
                k = ktop;
              ELSE
                k = kdetr(inn);
              ENDIF
              ! *********** Amisk NO PLUME MIXING
              Rsource = 0.0 
              IF (k > k1s) THEN  
                DO kk = ktop+1,kms
                  Rsource  = Rsource + (tracerpp(kk,l,itr) / 1000.0) * Qpss(kk,inn)
                ENDDO
                Rsource = Rsource - (tracerpp(k,l,itr) / 1000.0) * Qpss(k,inn)
                qdenom = SUM(Qpss(ktop+1:kms,inn)) - Qpss(k,inn)
                IF (ABS(Qpss(k,inn)) > 1.0E-12) THEN
                  IF (ABS(qdenom) > 1.0E-12) THEN   
                    Rsource = Rsource / qdenom + trpss(nn,itr) * dy / dfL(nn) / Qpss(k,inn)
                  ELSE
                    Rsource = tracerpp(k,l,itr) / 1000.0
                  ENDIF
                ELSE
                  Rsource = tracerpp(k,l,itr) / 1000.0
                ENDIF
              ELSE
                IF (ABS(Qpss(k,inn)) > 1.0E-12) THEN
                  Rsource = trpss(nn,itr) * dy / dfL(nn) / Qpss(k,inn)
                ELSE
                  Rsource = tracerpp(k,l,itr) / 1000.0
                ENDIF
              ENDIF
              Rpss(k,inn,itr) = Rsource * 1000.0! is in g/m3
            ! end if
            ! print*, 'Rsource and Rpss updated for tracer', itr, 'LDO', LDO, 'device ', inn
            ! print*, 'k', k, 'ktop', ktop, 'kdetr', kdetr(inn), 'ktop+1', ktop+1, 'kms', kms
            ! print*, 'Rpss', Rpss(k, inn, itr)
            ! print*, 'trpss', trpss(nn,itr)
            ! print*, 'dy', dy
            ! print*, 'dfL', dfL(nn)
            ! print*, 'Qpss', Qpss(k,inn)
            ! print*, 'calc', trpss(nn,itr)*dy/dfL(nn)/Qpss(k,inn)
            ! *********** Amisk NO PLUME MIXING
            !DO kk  = k1,kms
            !   Qpss(kk,inn) = 0.0
            !ENDDO
            !Qpss(k-1,inn)= -0.1
            !Qpss(k  ,inn)=  0.1
            !Tpss(k  ,inn)= salpp(k,l)
            !Rpss(k,inn,itr)= Qpss(k,inn)*tracerpp(k,l,itr)+trpss(nn,itr)*dy/dfL(nn)/Qpss(k,inn)
            ! *********** Amisk NO PLUME MIXING
          ENDDO
        ENDIF
      ENDDO
    END SELECT

  ENDDO
             

END SUBROUTINE PointSourceSinkSolve