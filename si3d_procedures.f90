!************************************************************************
                          MODULE si3d_procedures
!************************************************************************
!
!  Purpose: Procedures for the semi-implicit 3-D (si3d) hydrodynamic
!           model.
!
!-------------------------------------------------------------------------

   USE omp_lib
   USE si3d_Types
   USE si3d_ecomod
   USE si3d_BoundaryConditions
   USE si3d_Mixing
   USE si3d_Utils

   IMPLICIT NONE
   SAVE

CONTAINS


!***********************************************************************
SUBROUTINE init
!***********************************************************************
!
!  Purpose: To define initial conditions for the simulation  .
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: i, j, k, l, ios,     &
              kmx  , kmy  , kms,   &
              nwlsp, nwlup, nwlvp, js, c
   REAL :: x, rhoz, salz, zeta0, amp, deltZ

   ! ... Initialize time step counter, time in seconds and hours from &
   !     start of simulations (these are global variables defined in
   !     si3d_types)
   n = 0; its = 0; thrs = 0.0E0

   !.....Define initial water surface elevations.....
   SELECT CASE (testcase)

   CASE(1)

     PRINT *, '**** Free surface initilized for SW test case ****'
     amp = 0.05;
     DO i = 1, im1
       x = REAL(i*idx) - 1.5*dx
       zeta0 = amp*COS(pi*x/xl)
       deltZ = zlevel(k1+1); ! Move initial free surface to first interface
       DO j = 1, jm1
         sp  (ij2l(i,j)) = zeta0  - deltZ;
         spp (ij2l(i,j)) = zeta0  - deltZ;
         s   (ij2l(i,j)) = zeta0  - deltZ;
       END DO
     END DO

   CASE DEFAULT

     s    = zetainit;
     sp   = zetainit;
     spp  = zetainit;

   END SELECT

   ! ... Define thickness of cells a s-,u- & v-points
   hup = ZERO;
   hvp = ZERO;
   hp  = ZERO;
   k1z = km1;
   k1u = km1;
   k1v = km1;
   DO l = 1, lm

      ! ... Map 3D-(i,j) from 2D-l indexes
      i = l2i(l); j = l2j(l);

      ! ... At zeta-points
      kms = kmz(l)
      nwlsp = 0
      DO k = k1, kms
          hp (k,l)=AMIN1(zlevel(k+1),hhs(l)) -            &
          &        AMAX1(zlevel(  k),-sp(l))
          IF(hp(k,l) > HMIN) THEN
            nwlsp = nwlsp + 1;
            IF(nwlsp==1) k1z(l) = k
          ELSE
            hp(k,l)=ZERO;
          ENDIF
      ENDDO

      ! Set zeta = hhs(i,j) for columns with mask2d = TRUE (i.e.
      ! potentially wett) but intitially dry (k1z = km1).
      IF (k1z(l) == km1) THEN
         s  (l) = -hhs(l)+HMIN;
         sp (l) = -hhs(l)+HMIN;
         spp(l) = -hhs(l)+HMIN;
      ENDIF

      ! ... At u-points
      IF (mask2d(i+1,j)) THEN
        kmx = MIN(kmz(l),kmz(lEC(l)))
        nwlup = 0
        DO k = k1, kmx
          hup(k,l)=AMIN1(zlevel(k+1),hhu(l)) -            &
          &        AMAX1(zlevel(  k),-(sp(l)+sp(lEC(l)))/2.)
          IF(hup(k,l) > HMIN) THEN
            nwlup = nwlup + 1;
            IF(nwlup==1) k1u(l) = k
          ELSE
            hup(k,l)=ZERO;
          ENDIF
        ENDDO
      ENDIF

      ! ... At v-points
      IF (mask2d(i,j+1)) THEN
        kmy = MIN(kmz(l),kmz(lNC(l)))
        nwlvp = 0
        DO k = k1, kmy
          hvp(k,l)=AMIN1(zlevel(k+1),hhv(l)) -            &
          &        AMAX1(zlevel(  k),-(sp(l)+sp(lNC(l)))/2.)
          IF(hvp(k,l) > HMIN) THEN
            nwlvp = nwlvp + 1;
            IF(nwlvp==1) k1v(l) = k
          ELSE
            hvp(k,l)=ZERO;
          ENDIF
        ENDDO
      ENDIF

   ENDDO

   hupp = hup; hu = hup;
   hvpp = hvp; hv = hvp;
   hpp  = hp ; h  = hp
   contNG=0;


   !.....Initialize 1-d arrays.....
   uout  = 0.0;
   vout  = 0.0;
   wout  = 0.0;
   uhout = 0.0
   Avout = 0.0;
   Dvout = 0.0;
   sal1 = 0.0;
   ds = 0.0;

   !.....Initialize solution arrays.....
   eagx = 0.0;
   eagy = 0.0;
   earx = 0.0;
   eary = 0.0
   sx = 0.0;
   sy = 0.0;
   dd = 0.0;
   qq = 0.0;
   rr = 1.0;

   !.....Initialize velocity arrays.....
   u=u0; up=u0; upp=u0; uh=u0h0; uhp=u0h0; uhpp=u0h0
   v=v0; vp=v0; vpp=v0; vh=v0h0; vhp=v0h0; vhpp=v0h0
   wp=w0

   ! ... Initialize eddy coefficient arrays .....
   Av=0.0;
   Dv=0.0;

   !.....Initialize arrays used in the soln of the matrix mom eq.....
   ex=0.0; agx = 0.0; arx = 0.0; agy = 0.0; ary = 0.0

   !.....Initialize 3D active scalar and density arrays.....
   CALL InitializeScalarFields

   ! ... Inialize turbulence quanties for Turbulence Model
   CALL InitializeTurbulenceModel

END SUBROUTINE init

!************************************************************************
SUBROUTINE InitializeScalarFields
!************************************************************************
!
!  Purpose: Initialize scalar fields - The initial condition field
!           is either read from an ASCII file (si3d_init.txt) or it
!           is initialized using fields which will excite specific
!           hydrodynamic responses in the lake (test case = 2).
!           For test case 1 - the initial conditions are also hardcoded
!           but I use uniform temperature for that cases.
!
!-------------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: i, j, k, l, ios, imm1, jmm1, kmm1, ncols, ncols1, nc, &
              nsets, ia, ib, nn, ntr1
   REAL    :: Vamp, rhoamp, Ts, Tb,   &  ! Used to initialize IW-problem
              NBV, meandepth, length, &
              rhohere, x, z, rhos, rhob
   CHARACTER(LEN=18)  :: initfmt
   INTEGER, PARAMETER :: InitProc = 4
   REAL, ALLOCATABLE, DIMENSION(:,:) :: Scalardepthile

   SELECT CASE (initproc)

   ! ... OPTION 1 - Surface Seiche (use uniform temperatures) ------
   CASE (1)

     salp  = 15.
     sal   = salp;
     salpp = salp;

     ! ... Initialize Non-active scalar fields
     IF (ntr > 0) THEN
       DO nn = 1, ntr;
         tracer(:,:,nn) = sal;
       ENDDO
       tracerpp = tracer;
     ENDIF

   ! ... OPTION 2 - Internal Seiche (use analytical solution) ------
   CASE (2)

     Vamp =  0.10;
     Ts   = 25.00;
     Tb   = 15.00; ! Parameters used to define the solution
     meandepth = zl; ! FLOAT(km-k1+1)*ddz
     length    = FLOAT(im-i1+1)*dx
     rhos = 1028.*(1.-1.7E-4*(Ts-10.));		! Surface density
     rhob = 1028.*(1.-1.7E-4*(Tb-10.));		! Bottom  density
     drho = rhos - rhob;		        ! Change in density from top to bottom
     NBV=SQRT(-g/rhos*drho/meandepth);          ! Brunt-Vaisala frequency
     rhoamp=rhos*Vamp*NBV/g;
     DO l = 1, lm
       i = l2i(l); j = l2j(l);
       x = FLOAT(i) * dx - 1.5 * dx
       DO k = k1, km
          z = zlevel(k+1) - 0.5 * hp(k,l); ! FLOAT(km-k1+1)*ddz
          rhohere = rhos -z*drho/meandepth+rhoamp*COS(pi*x/length)*SIN(pi*z/meandepth);
          salp(k,l)= 10.-((rhohere-1028.)/1028.)/1.7E-4;
       ENDDO
     END DO
     sal = salp;
     salpp = salp;

     ! ... Initialize Non-active scalar fields
     IF (ntr > 0) THEN
       DO nn = 1, ntr;
         tracer(:,:,nn) = sal;
       ENDDO
       tracerpp = tracer;
     ENDIF
     PRINT *, '**** Scalar field initilized for IW test case ****'

   ! ... OPTION 3 - Use analytical solution in half a closed basin to test
   !                the nesting algorithms nesting. All variables defining the basin
   !                & the IW need to be the same in the fine & coarse grid -
   !                In the fine grid we only modify the length and x.
   CASE (3)

     Vamp = 0.10; Ts = 25.00; Tb = 15.0; ! Make sure these constants are as in CASE (1)
     meandepth = zl; ! FLOAT(km-k1+1)*ddz
     length    = FLOAT(im-i1+1)*dx; length = length * 2.;
     rhos = 1028.*(1.-1.7E-4*(Ts-10.));		! Surface density
     rhob = 1028.*(1.-1.7E-4*(Tb-10.));		! Bottom  density
     drho = rhos - rhob;		! Change in density from top to bottom
     NBV=SQRT(-g/rhos*drho/meandepth);
     rhoamp=rhos*Vamp*NBV/g;
     DO l = 1, lm
       i = l2i(l); j = l2j(l);
       x = FLOAT(i) * dx - 1.5 * dx; x = x + length/2.;
       DO k = k1, km
          z = zlevel(k+1) - 0.5 * hp(k,l); ! z = FLOAT(k) * ddz - 1.5 * ddz
          rhohere = rhos -z*drho/meandepth+rhoamp*COS(pi*x/length)*SIN(pi*z/meandepth);
          salp(k,l)= 10.-((rhohere-1028.)/1028.)/1.7E-4;
       ENDDO
     END DO
     sal = salp;
     salpp = salp;

     ! ... Initialize Non-active scalar fields
     IF (ntr > 0) THEN
       DO nn = 1, ntr;
         tracer(:,:,nn) = sal;
       ENDDO
       tracerpp = tracer;
     ENDIF

   ! ... All other options - Initialize from file ---------------------
   CASE DEFAULT

     !.....Open initial condition file.....
     sal_ic_file = 'si3d_init.txt'
     OPEN (UNIT=i4, FILE='si3d_init.txt', STATUS="OLD", FORM="FORMATTED", IOSTAT=ios)
     IF(ios /= 0) CALL open_error ( "Error opening "//sal_ic_file, ios )

     !.....Allocate space for local variables used to read IC ...
     ALLOCATE ( Scalardepthile (km1, ntr+1), STAT = ios )
     IF (ios /= 0) THEN; PRINT *, 'Error alloc. init. arrays'; STOP; ENDIF

     ! Skip over first five header records in open boundary condition file
     READ (UNIT=i4, FMT='(/////)', IOSTAT=ios)
     IF (ios /= 0) CALL input_error ( ios, 13 )

     ! Write the format of the data records into an internal file
     WRITE (UNIT=initfmt, FMT='("(10X,",I3,"G11.2)")') ntr+1

     ! Read data array and store it in array Scalardepthile
     print *,"km1:",km1
     DO k = 1, km1
       print *,"leo:",k
       READ (UNIT=i4, FMT=initfmt, IOSTAT=ios) &
            (Scalardepthile(k,nn), nn = 1, ntr+1)
       IF (ios /= 0) CALL input_error ( ios, 14 )
     END DO

     ! ... Initialize the active scalar field (allways)
     salp = 0.0
     DO k = 1, km1;
        salp(k,:) = Scalardepthile(k,1)
     END DO;
     sal = salp;
     salpp = salp;

     ! ... Initialize non-active scalar fields (if requested)
     IF (ntr > 0) THEN
       tracer = 0.0;
       IF (ecomod < 0 ) THEN
         CALL InitTracerCloud
       ELSE
       DO  nn = 1, ntr
         DO k = 1, km1
           tracer(k,:,nn) = Scalardepthile(k,nn+1)
         ENDDO
       END DO ! ... End loop over tracers
       END IF
       tracerpp = tracer;
     ENDIF

     ! ... Deallocate array holding scalar concs.
     DEALLOCATE ( Scalardepthile )

     ! ... Close io file
     CLOSE (i4)



   END SELECT

   ! ... Initialize density field at time n-1 & n
   DO l = 1, lm1; DO k = k1, km1;
      rhop(k,l) = densty_s ( salp(k,l), t0 ) - 1000.
   END DO; END DO

END SUBROUTINE InitializeScalarFields

!***********************************************************************
SUBROUTINE fd(n,t_exmom2,t_matmom2,t_matcon2,Bhaxpp,Bhaypp,Bth,Bth1,Bstart,Bend, &
& lSCH,lNCH,lWCH,lECH,Bex,Bth2,Beagx,Bearx,Bagx,Barx,Beagy,Beary,Bagy,Bary,Bsx,Bsy, &
& Bdd,Bqq,Brr,Bth3,Bth4,istep,lastiter, &
& ShearProduction,BuoyancyProduction, Dissipation, TKinE, &
& Qsw,Qn,Qlw,eta,Ta,Pa,RH,Cc,Qsw2dB,Qlw2dB,Ta2dB,RH2dB,Cc2dB,uair2dB,vair2dB, &
& uairB,vairB,cdwB,heatSourceB,QswFrB,iter,bclncxB,hupdrhoB,its,thrs)
!***********************************************************************
!
!  Purpose: fd  is the supervisory subroutine for the finite
!           difference method. It advances the solution one time
!           step. Vertical diffusion is treated implicitly.
!
!-----------------------------------------------------------------------

   REAL, INTENT(INOUT) :: t_exmom2,t_matmom2,t_matcon2
   REAL, INTENT(IN) :: thrs, its
   INTEGER, INTENT(IN) :: Bstart,Bend,istep,n,lastiter,iter
   REAL, DIMENSION (1:km1,Bstart:Bend+1), INTENT(INOUT) :: Bhaxpp, Bhaypp, Bth, Bth1, Bex, heatSourceB,QswFrB
   REAL, DIMENSION (1:km1,Bstart:Bend+1), INTENT(INOUT) :: Bth2, Bagx, Barx,Bagy,Bary,Bth3,Bth4
   REAL, DIMENSION (Bstart:Bend+1), INTENT(INOUT) :: Beagx, Bearx,Beagy,Beary,Bsx,Bsy
   REAL, DIMENSION (Bstart:Bend+1), INTENT(INOUT) :: Bdd, Bqq,Brr,uairB,vairB,cdwB
   INTEGER, DIMENSION (Bstart:Bend+1), INTENT(IN) :: lSCH,lNCH,lWCH,lECH
   REAL, INTENT(INOUT) :: Qsw,Qn,Qlw,eta,Ta,Pa,RH,Cc
   REAL(real_G1), INTENT(INOUT) :: ShearProduction,BuoyancyProduction, Dissipation, TKinE
   REAL, DIMENSION (nmetstat), INTENT(INOUT) :: Qsw2dB,Qlw2dB,Ta2dB,RH2dB,Cc2dB,uair2dB,vair2dB
   REAL,DIMENSION (1:km1), INTENT(INOUT) :: bclncxB,hupdrhoB

   ! ... Local variables
   INTEGER :: itr,l,lol,liter
   REAL, EXTERNAL :: TIMER
   REAL :: tsbar,tebar

   ! ... Define tz used in all following routines ...............
   tz = 1.0/istep
   !print *,"hihihi:",omp_get_thread_num()
   !$omp barrier
   ! ... Read in nested boundary conditions, if specified .......
   CALL readbcNGB(thrs)
   !$omp barrier
   !print *,"hihihi2:",omp_get_thread_num()
   !.....Assign new values of s or u/v along open boundaries.....
   CALL openbcUVH(thrs)

   !$omp barrier
   !print *,"hihihi3:",omp_get_thread_num()
   !.... Find magnitude of sources and sinks and their scalar loads
   IF (iopss > 0) CALL PointSourceSinkSolve(n,istep,thrs)

   !$omp barrier
   ! ... Assign boundary conditions at free surface .............
   CALL surfbc(n,istep,thrs)

   !$omp barrier
   !print *,"hihihi4:",omp_get_thread_num()
   !.....Assing horizonal eddy viscosity and diffusivity at n ...........
   CALL UpdateHorizontalMixingCoefficients
   !$omp barrier
   !print *,"hihihi5:",omp_get_thread_num()
    !.....Evaluate explicit terms in x-momentum eq................
   CALL exmom(1)
   !$omp barrier
   !print *,"ex1:",sum(ex(:,:))

   !print *,"hihihi6:",omp_get_thread_num()
   !.....Solve a tridiagonal system at each horizontal node
   !     to obtain the matrices for the x-momentum equations.....
   CALL matmom(1,t_matmom2,Bstart,Bend,Bex,Beagx,Bearx,Bagx,Barx,Beagy,Beary,Bagy,Bary,uairB,vairB,cdwB,bclncxB,hupdrhoB)
   !$omp barrier
   !print *,"hihihi7:",omp_get_thread_num()
   !.....Evaluate explicit terms in y-momentum eq................
   CALL exmom(2)

   !$omp barrier
   !print *,"hihihi8:",omp_get_thread_num()
   !$omp barrier
!!   if(omp_get_thread_num() .EQ. 0) THEN
!!   print *,"ex2:",sum(ex(:,:))
!!   end if
   !eagy=0
   !Beagy=0
!!   DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
!!     l = id_column(liter)

!!     Bagy(:,l) = 0.0
!!   END DO

   !.....Solve a tridiagonal system at each horizontal node
   !     to obtain the matrices for the y-momentum equations.....
   CALL matmom(2,t_matmom2,Bstart,Bend,Bex,Beagx,Bearx,Bagx,Barx,Beagy,Beary,Bagy,Bary,uairB,vairB,cdwB,bclncxB,hupdrhoB)

   !print *,"hihihi9:",omp_get_thread_num()
   !$omp barrier

   CALL matcon(t_matcon2,Bstart,Bend,lWCH,lSCH,Beagx,Bearx,Beagy,Beary,Bsx,Bsy,Bdd,Bqq,Brr)
   !print *,"hihihi10:",omp_get_thread_num()
   !$omp barrier
   !.....Solve implicit system of equations for zeta............
   !CALL SolverBlock ! Original formulation writen by P.E. Smith
   CALL SolverSparse(n,Bstart,Bend,lWCH,lSCH,Bsx,Bsy,Bqq,Brr,iter,istep,thrs) ! Formulation by F.J. Rueda
   !print *,"hihihi11:",omp_get_thread_num()
   !.....Reassign new values of s or u/v along open boundaries.....
   CALL openbcUVH(thrs)
   !
!!   if(omp_get_thread_num() .EQ. 0) THEN
!!   print *,"vhavel:",sum(vh(:,:))
!!   vh(:,:)=0.0
!!   end if

!!   DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
!!     l = id_column(liter)
!!     ary(:,l)=Bary(:,l)
!!     agy(:,l)=Bagy(:,l)
!!   END DO
   !$omp barrier
!!   if(omp_get_thread_num() .EQ. 0) THEN
!!   print *,"ary:",sum(ary(:,:))
!!   print *,"agy:",sum(agy(:,:))
!!   end if
   !.....Solve for velocities explicitly. If DRYING occurs
   !     at this point, dry cells are removed and k1z/k1u/k1v
   !     recalculated (wetting is done after finishing the
   !     calculations for a given time step or iteration) ......
   CALL vel(Bstart,Bend,Bagx,Barx,Bagy,Bary)
   !print *,"hihihi12:",omp_get_thread_num()
   !$omp barrier
!!   if(omp_get_thread_num() .EQ. 0) THEN
!!   print *,"vhdvel:",sum(vh(:,:))
!!   end if
   ! ... Solve for active scalar transport
    IF (isal /= 0) THEN
     CALL exsal(Bstart,Bend,lSCH,lNCH,lECH,lWCH,Bhaxpp,Bhaypp,Bth3,Bth4,Bth2,Bex,thrs)
     !$omp barrier
     CALL imsal(Bstart,Bend,Bex,heatSourceB)
     !$omp barrier
     CALL openbcSCA(thrs)
   END IF

   !print *,"hihihi13:",omp_get_thread_num()


   !.....Solve for non-active scalar transport
   IF (ntr      >  0 .AND. &
       niter    >  0 .AND. &
       lastiter == 1      ) THEN
       IF      (ecomod <  0) THEN
         CALL srcsnk00
       ELSE IF (ecomod == 0) THEN
         CALL srcsnk00
       ELSE IF (ecomod == 1) THEN
         CALL srcsnkWQ
       ELSE IF (ecomod == 2) THEN
         CALL srcsnkSZ
       ELSE IF (ecomod == 3) THEN
         CALL srcsnkSD
       ENDIF
       DO itr = 1, ntr
         IF (ecomod < 0 .AND. ( trct0(itr) > n .OR. trctn(itr) < n ) ) CYCLE
         CALL exTracer (itr,Bstart,Bend,Bhaxpp,Bhaypp,Bth3,Bth4,Bth2,lSCH,lNCH,lECH,lWCH,Bex,thrs)
         CALL imTracer (itr,Bstart,Bend,Bex)
         CALL openbctracer (itr,thrs)
       ENDDO
   ENDIF
   !$omp barrier
!!   if(omp_get_thread_num() .EQ. 0) THEN
!!   print *,"vhavel2:",sum(vh(:,:))
!!   end if
   ! ... Recalculate near surface velocity to account
   !     for WETTING & recalculate k1z, k1u & k1v ..............
   CALL vel2
   !$omp barrier
!!   if(omp_get_thread_num() .EQ. 0) THEN
!!   print *,"u:",sum(u(:,:))
!!   print *,"v:",sum(v(:,:))
!!   print *,"wp:",sum(wp(:,:))
!!   print *,"vhp:",sum(vhp(:,:))
!!   print *,"hvp:",sum(hvp(:,:))
!!   print *,"vh:",sum(vh(:,:))
!!   end if
   !print *,"hihihi14:",omp_get_thread_num()
   !.....Smooth solution on leapfrog step if ismooth>=1.........
   IF (ismooth >= 1 .AND. istep == 1) CALL smooth

   !.....Assing eddy viscosity and diffusivity at n+1 ...........
   CALL UpdateMixingCoefficients(Bstart,Bend,istep,uairB,vairB,cdwB, &
   & ShearProduction,BuoyancyProduction, Dissipation, TKinE)
   !print *,"hihihi15:",omp_get_thread_num()
   !$omp barrier
END SUBROUTINE fd

!***********************************************************************
SUBROUTINE exmom ( ieq  )
!***********************************************************************
!
!  Purpose: To evaluate the explicit terms (advection, Coriolis, and
!           horizontal diffusion) in the momentum equations. The sum
!           of these terms are saved in the 3-D array  ex(i,j,k)
!           which is the primary output from this subroutine. The
!           horizontal advection terms can be evaluated by either upwind
!           or centered differencing.
!           (note: ex(i,j,k)  is a temporary workspace array that is
!           used for storing the explicit terms from both the x- and y-
!           momentum eqs and is also used later for storing the explicit
!           terms from the salinity eq in SUBR. exsal.)
!
!  Dummy argument:
!  ieq    = Parameter indicating whether the explicit terms in the
!           X  or  Y  momentum equation are to be evaluated.
!           (1=x-momentum, 2=y-momentum)
!
!-----------------------------------------------------------------------
   !.....Argument.....
   INTEGER, INTENT(IN) :: ieq


   !.....Local variables.....
   REAL :: twodt1
   REAL :: corx, cory, advx, advy, hdx, hdy, uE, uW, vN, vS, wU, wD, &
           scW, scE, scN, scS, scU, scD
   INTEGER :: i, j, k, l, istat, kmx, kmy, k1x, k1y, k1ne,liter
   INTEGER :: kmne, nwlayers, nn, is, ie, js, je,no

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)
   n_exmom = n_exmom + 1

   !.....Constant.....
   twodt1 = twodt*tz

   !.....Choose to evaluate explicit terms for either the x- or y-mom eq.....
   SELECT CASE (ieq)


   ! -----X-momentum equation-----
   CASE (1)

      !.....Calculate coefficient arrays haxpp&haypp for use
      !     in horizontal diffusion term & th1,th for use
      !     in vertical advection term in the x-momentum
      haxpp(:,lm1) = 0.0;
      haypp(:,lm1) = 0.0
!      print *, "lEC127:",lEC(127),"maskij:",mask2d(l2i(127)+1,l2j(127))

      DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
     l = id_column(liter)

        !.....Compute layer number for the bottom wet u-pt.......
        kmx = MIN(kmz(lEC(l)), kmz(l))
        k1x =                 k1u(l)
        if(k1x > 1) THEN
        haypp(1:k1x-1,l) = 0.0;
        haxpp(1:k1x-1,l) = 0.0;
        th   (1:k1x-1,l) = 0.0;
        th1  (1:k1x-1,l) = 0.0;
        end if
        if(kmx < km1) THEN
        haypp(kmx+1:km1,l) = 0.0;
        haxpp(kmx+1:km1,l) = 0.0;
        th   (kmx+1:km1,l) = 0.0;
        th1  (kmx+1:km1,l) = 0.0;
        end if
        ! ... Map 3D-(i,j) from 2D-l indexes
        i = l2i(l); j = l2j(l);

        ! ... Cycle if E-column is dry
        IF ( .NOT. mask2d(i+1,j) ) THEN
        haypp(k1x:kmx,l) = 0.0;
        haxpp(k1x:kmx,l) = 0.0;
        th   (k1x:kmx,l) = 0.0;
        th1  (k1x:kmx,l) = 0.0;
        CYCLE
        end if

        ! ... Horizontal diffusion ...........................
        IF ( ihd == 1 ) THEN ! Constant
          DO k = k1x,kmx
            haxpp(k,l)=Ax0*MIN(hupp(k,lEC(l)),hupp(k,l))
            haypp(k,l)=Ay0*MIN(hupp(k,lNC(l)),hupp(k,l))
          ENDDO
        ELSEIF ( ihd > 1) THEN ! Smagorinsky
          DO k = k1x,kmx
            haxpp(k,l)= kh(k,lEC(l))*MIN(hupp(k,lEC(l)),hupp(k,l))
            haypp(k,l)=(kh(k,lEC(    l ))+   &
                        kh(k,        l  )+   &
                        kh(k,lNC(    l ))+   &
                        kh(k,lEC(lNC(l))) )* &
                        0.25 * MIN(hupp(k,lNC(l)),hupp(k,l))
          ENDDO
        ENDIF

        !.....Calculate weighting arrays for vertical advection
        DO k = k1x, kmx
          th (k,l) = hup(k-1,l)/(hup(k-1,l)+hup(k,l))
          th1(k,l) = 1.-th(k,l)
        ENDDO

        !.....Set th=th1 at the free surface & bottom
        th (k1x  ,l) = 0.0;
        th1(k1x  ,l) = 0.0;
        th (kmx+1,l) = 0.5;
        th1(kmx+1,l) = 0.5;

      END DO

      IF(omp_get_thread_num ( )>0)THEN
      DO liter = lhiW(omp_get_thread_num ( )+1), lhfW(omp_get_thread_num ( )+1)
        l = id_column(liter)
        !.....Compute layer number for the bottom wet u-pt.......
        kmx = MIN(kmz(lEC(l)), kmz(l))
        k1x =                 k1u(l)

        if(k1x > 1) THEN
        haypp(1:k1x-1,l) = 0.0;
        haxpp(1:k1x-1,l) = 0.0;
        th   (1:k1x-1,l) = 0.0;
        th1  (1:k1x-1,l) = 0.0;
        end if
        if(kmx < km1) THEN
        haypp(kmx+1:km1,l) = 0.0;
        haxpp(kmx+1:km1,l) = 0.0;
        th   (kmx+1:km1,l) = 0.0;
        th1  (kmx+1:km1,l) = 0.0;
        end if

        ! ... Map 3D-(i,j) from 2D-l indexes
        i = l2i(l); j = l2j(l);

        ! ... Cycle if E-column is dry
        IF ( .NOT. mask2d(i+1,j) ) THEN
        haypp(k1x:kmx,l) = 0.0;
        haxpp(k1x:kmx,l) = 0.0;
        th   (k1x:kmx,l) = 0.0;
        th1  (k1x:kmx,l) = 0.0;
        CYCLE
        end if

        ! ... Horizontal diffusion ...........................
        IF ( ihd == 1 ) THEN ! Constant
          DO k = k1x,kmx
            haxpp(k,l)=Ax0*MIN(hupp(k,lEC(l)),hupp(k,l))
            haypp(k,l)=Ay0*MIN(hupp(k,lNC(l)),hupp(k,l))
          ENDDO
        ELSEIF ( ihd > 1) THEN ! Smagorinsky
         DO k = k1x,kmx
            haxpp(k,l)= kh(k,lEC(l))*MIN(hupp(k,lEC(l)),hupp(k,l))
            haypp(k,l)=(kh(k,lEC(    l ))+   &
                        kh(k,        l  )+   &
                        kh(k,lNC(    l ))+   &
                        kh(k,lEC(lNC(l))) )* &
                        0.25 * MIN(hupp(k,lNC(l)),hupp(k,l))
          ENDDO
        ENDIF

        !.....Calculate weighting arrays for vertical advection
        DO k = k1x+1, kmx
          th (k,l) = hup(k-1,l)/(hup(k-1,l)+hup(k,l))
          th1(k,l) = 1.-th(k,l)
        ENDDO

        !.....Set th=th1 at the free surface & bottom
        th (k1x  ,l) = 0.0;
        th1(k1x  ,l) = 0.0;
        th (kmx+1,l) = 0.5;
        th1(kmx+1,l) = 0.5;

        l=lWC(l)
        if(l .EQ. lm1) CYCLE
        !.....Compute layer number for the bottom wet u-pt.......
        kmx = MIN(kmz(lEC(l)), kmz(l))
        k1x =                 k1u(l)

        if(k1x > 1) THEN
        haypp(1:k1x-1,l) = 0.0;
        haxpp(1:k1x-1,l) = 0.0;
        th   (1:k1x-1,l) = 0.0;
        th1  (1:k1x-1,l) = 0.0;
        end if
        if(kmx < km1) THEN
        haypp(kmx+1:km1,l) = 0.0;
        haxpp(kmx+1:km1,l) = 0.0;
        th   (kmx+1:km1,l) = 0.0;
        th1  (kmx+1:km1,l) = 0.0;
        end if

        ! ... Map 3D-(i,j) from 2D-l indexes
        i = l2i(l); j = l2j(l);


        IF ( .NOT. mask2d(i+1,j) .OR. .NOT. mask2d(i,j) ) THEN


        haypp(k1x:kmx,l) = 0.0;
        haxpp(k1x:kmx,l) = 0.0;
        th   (k1x:kmx,l) = 0.0;
        th1  (k1x:kmx,l) = 0.0;
        CYCLE
        END IF



        ! ... Horizontal diffusion ...........................
        IF ( ihd == 1 ) THEN ! Constant
          DO k = k1x,kmx
            haxpp(k,l)=Ax0*MIN(hupp(k,lEC(l)),hupp(k,l))
            haypp(k,l)=Ay0*MIN(hupp(k,lNC(l)),hupp(k,l))
          ENDDO
        ELSEIF ( ihd > 1) THEN ! Smagorinsky
          DO k = k1x,kmx
            haxpp(k,l)= kh(k,lEC(l))*MIN(hupp(k,lEC(l)),hupp(k,l))
            haypp(k,l)=(kh(k,lEC(    l ))+   &
                        kh(k,        l  )+   &
                        kh(k,lNC(    l ))+   &
                        kh(k,lEC(lNC(l))) )* &
                        0.25 * MIN(hupp(k,lNC(l)),hupp(k,l))
          ENDDO
        ENDIF

        !.....Calculate weighting arrays for vertical advection
        DO k = k1x+1, kmx
          th (k,l) = hup(k-1,l)/(hup(k-1,l)+hup(k,l))
          th1(k,l) = 1.-th(k,l)
        ENDDO

        !.....Set th=th1 at the free surface & bottom
        th (k1x  ,l) = 0.0;
        th1(k1x  ,l) = 0.0;
        th (kmx+1,l) = 0.5;
        th1(kmx+1,l) = 0.5;

      end do
      end if


      !......Calulate the explicit terms by sweeping over interior u-pts.....

      DO liter = lhiW(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
         If(liter == 0) CYCLE
         l = id_column(liter)
          ! Compute the layer number for the bottom wet u-pt
         kmx = MIN(kmz(lEC(l)), kmz(l))
         k1x =                 k1u(l)

        if(k1x > 1) THEN
        ex(1:k1x-1,l) = 0.0;

        end if
        if(kmx < km1) THEN
        ex(kmx+1:km1,l) = 0.0;

        end if
         ! ... Map 2D-l into 3D-(i,j) indexes
         i = l2i(l); j = l2j(l);

         ! ... Cycle if E-column is dry
         IF ( .NOT. mask2d(i+1,j) ) THEN
         ex(k1x:kmx,l)=0.0
         CYCLE
         END IF

         ! Compute explicit term
         DO k = k1x,kmx

            ! ... For u-layers connecting wett & dry cells neglect
            !     contribution from advective, coriolis & diffusion
            IF ( hp(k,l) <= ZERO .OR. hp(k,lEC(l)) <= ZERO) THEN
              ex(k,l) = uhpp(k,l)
              CYCLE
            ENDIF

            !.....Coriolis.....
            corx = 0.25 * f * (vhp(k,     lEC(l) ) + vhp(k,    l )       &
                  &           +vhp(k, lSC(lEC(l))) + vhp(k,lSC(l)))

            !.....Advection
            uE = uhp(k,    lEC(l) ) +uhp(k,    l )
            uW = uhp(k,        l  ) +uhp(k,lWC(l))
            vN = vhp(k,    lEC(l) ) +vhp(k,    l )
            vS = vhp(k,lSC(lEC(l))) +vhp(k,lSC(l))
            wU = wp (k,    lEC(l) ) +wp (k,    l ); IF ( k == k1x ) wU = 0.0;
            wD = wp (k+1,  lEC(l) ) +wp (k+1  ,l ); IF ( k == kmx ) wD = 0.0;

            SELECT CASE (itrmom)

            CASE (1)   ! Centered differences with th & th1

              scE = up(k,lEC(l))+up(k,    l )
              scW = up(k,lWC(l))+up(k,    l )
              scN = up(k,lNC(l))+up(k,    l )
              scS = up(k,lSC(l))+up(k,    l )
              advx = (uE * scE - uW * scW ) / fourdx +          &
                     (vN * scN - vS * scS ) / fourdy
              advx=advx+(wU*(th (k  ,l)* up(k  ,l)  +          &
                             th1(k  ,l)* up(k-1,l)) -          &
                         wD*(th (k+1,l)* up(k+1,l)  +          &
                             th1(k+1,l)* up(k  ,l)) ) / 2.

            CASE (2) ! Upwinding all

              advx = ( (uE+ABS(uE))* upp(k,    l ) +           &
                       (uE-ABS(uE))* upp(k,lEC(l)) -           &
                       (uW+ABS(uW))* upp(k,lWC(l)) -           &
                       (uW-ABS(uW))* upp(k,    l ) ) / fourdx  &
                    +( (vN+ABS(vN))* upp(k,    l ) +           &
                       (vN-ABS(vN))* upp(k,lNC(l)) -           &
                       (vS+ABS(vS))* upp(k,lSC(l)) -           &
                       (vS-ABS(vS))* upp(k,    l ) ) / fourdy  &
                    +( (wU+ABS(wU)) * upp(k  ,l) +             &
                       (wU-ABS(wU)) * upp(k-1,l)) / 4.         &
                     -((wD+ABS(wD)) * upp(k+1,l) +             &
                       (wD-ABS(wD)) * upp(k  ,l)) / 4.

            CASE (3)  ! Centered differences - avoid computation of th1 and th

              scE = up(k,lEC(l))+up(k,    l )
              scW = up(k,lWC(l))+up(k,    l )
              scN = up(k,lNC(l))+up(k,    l )
              scS = up(k,lSC(l))+up(k,    l )
              scU = (up(k  ,l)*hup(k  ,l)+        &
                     up(k-1,l)*hup(k-1,l))/       &
                    (hup(k ,l)+hup(k-1,l))
              scD = (up(k  ,l)*hup(k  ,l)+        &
                     up(k+1,l)*hup(k+1,l))/       &
                    (hup(k ,l)+hup(k+1,l))
              advx = (uE * scE - uW * scW ) / fourdx +          &
                     (vN * scN - vS * scS ) / fourdy +          &
                     (wU * scU - wD * scD ) / 2.

            CASE (4) ! Upwinding for horizontal & centered for vertical

              advx = ( (uE+ABS(uE))* upp(k,    l ) +           &
                       (uE-ABS(uE))* upp(k,lEC(l)) -           &
                       (uW+ABS(uW))* upp(k,lWC(l)) -           &
                       (uW-ABS(uW))* upp(k,    l ) ) / fourdx  &
                    +( (vN+ABS(vN))* upp(k,    l ) +           &
                       (vN-ABS(vN))* upp(k,lNC(l)) -           &
                       (vS+ABS(vS))* upp(k,lSC(l)) -           &
                       (vS-ABS(vS))* upp(k,    l ) ) / fourdy
              advx=advx+(wU*(th (k  ,l)* up(k  ,l)  +          &
                             th1(k  ,l)* up(k-1,l)) -          &
                         wD*(th (k+1,l)* up(k+1,l)  +          &
                             th1(k+1,l)* up(k  ,l)) ) / 2.

            END SELECT

            !.....Horizontal diffusion.....
            IF ( ihd == 0) THEN
              hdx = 0.0E0
            ELSEIF ( ihd == 1) THEN
               hdx  =   (haypp(k,    l )*(upp(k,lNC(l))-upp(k,    l ))       &
                    & -  haypp(k,lSC(l))*(upp(k,    l )-upp(k,lSC(l))))/dydy &
                    & + (haxpp(k,    l )*(upp(k,lEC(l))-upp(k,    l ))       &
                    & -  haxpp(k,lWC(l))*(upp(k,    l )-upp(k,lWC(l))))/dxdx
            ELSEIF ( ihd >  1) THEN
              hdx  = 2.* (haxpp(k,   (l))*(upp(k,lEC(    l ))-upp(k,    l ))       &
                 &     -  haxpp(k,lWC(l))*(upp(k,        l  )-upp(k,lWC(l))))/dxdx &
                 &   +   (haypp(k,    l )*(upp(k,lNC(    l ))-upp(k,    l ))       &
                 &     -  haypp(k,lSC(l))*(upp(k,        l  )-upp(k,lSC(l))))/dydy &
                 &   +   (haypp(k,    l )*(vpp(k,lEC(    l ))-vpp(k,    l ))       &
                 &     -  haypp(k,lSC(l))*(vpp(k,lEC(lSC(l)))-vpp(k,lSC(l))))/dxdy
            ENDIF

            ! ... Adjust terms to account for boundary conditions - HardCoded for NortheDelta Study
            !IF((                           j>=  jm1 - 20          )   .OR. &
            !   (i<=591               .AND. j<=  60                )   .OR. &
            !   (i<=20                                             )   .OR. &
            !   (i>=990 .AND. i<=1010 .AND. j>=  74 .AND. j <=  81 )   .OR. &
            !   (                           j<=  20                )   .OR. &
            !   (i>=678 .AND. i<= 682 .AND. j>= 100 .AND. j <= 115 )) THEN
            !   corx = 0.0
            !   advx = 0.0
            !   hdx  = 2.*hdx
            !ENDIF

            !IF ((i<=86                 .AND. j>991                 )   .OR. & !CINTIA (SUTTER SLOUGH)
            !    (i>=613                .AND. j<=27                 )   .OR. & !CINTIA (SUTTER GEO)
            !    (i>679                 .AND. j<=114                ))  THEN  !CINTIA (SUTTER DCC)
            !  !corx = 0.0
            !  !advx = 0.0
            !  hdx = 2.*hdx
            !ENDIF

            ! ... Needed to keep simulations stable near the boundaries - TAHOE MAC
            !IF(i >= 2 .AND. i <= 10) THEN
            !  hdx   = 4.*hdx
            !  advx  = 0.0
            !  corx  = 0.
            !ENDIF
            !IF(j >= 2 .AND. j <= 13) THEN
            !  hdx   = 4.*hdx
            !  advx  = 0.0
            !  corx  = 0.
            !ENDIF
            !IF(i >= 2 .AND. i <= 7) THEN
            !  hdx   = 4.*hdx
            !  advx  = 0.0
            !  corx  = 0.
            !ENDIF

            !IF(i == 2) THEN
            !  hdx   = 4.*hdx
            !  advx  = 0.0
            !  corx  = 0.
            !ENDIF
            !IF(j == 2) THEN
            !  hdx   = 4.*hdx
            !  advx  = 0.0
            !  corx  = 0.
            !ENDIF
            !IF(j == 23) THEN !13
            !  hdx   = 4.*hdx
            !  advx  = 0.0
            !  corx  = 0.
            !ENDIF

            !IF (i >= (im1 - 20) ) THEN; !BeznarCOLA
            !   hdx  = 10.*hdx;
            !   corx = 0.0;
            !   advx = 0.0;
            !end IF

            !.....Final explicit term.....
            ex(k,l) = uhpp(k,l) - twodt1*(advx*iadv-corx-hdx)

        ENDDO
      END DO

      ! ... Recalculate ex for near bdry. cells
      IF (nopen > 0 ) CALL MODexmom4openBCX

   ! -----Y-momentum equation-----
   CASE (2)

      !.....Calculate coefficient arrays haxpp&haypp for use
      !     in horizontal diffusion term & th1,th for use
      !     in vertical advection term in the y-momentum

      haxpp(:,lm1) = 0.0;
      haypp(:,lm1) = 0.0

      DO liter = lhiW(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
         If(liter == 0) CYCLE
         l = id_column(liter)

         !.....Compute layer number for the bottom wet u-pt.......
        kmy = MIN(kmz(lNC(l)), kmz(l))
        k1y =                 k1v(l)

        if(k1y > 1) THEN
        haypp(1:k1y-1,l) = 0.0;
        haxpp(1:k1y-1,l) = 0.0;
        th   (1:k1y-1,l) = 0.0;
        th1  (1:k1y-1,l) = 0.0;
        end if
        if(kmy < km1) THEN
        haypp(kmy+1:km1,l) = 0.0;
        haxpp(kmy+1:km1,l) = 0.0;
        th   (kmy+1:km1,l) = 0.0;
        th1  (kmy+1:km1,l) = 0.0;
        end if

        ! ... Map 3D-(i,j) from 2D-l indexes
        i = l2i(l); j = l2j(l);

        ! ... Cycle if E-column is dry
        IF ( .NOT. mask2d(i,j+1) ) THEN
        haypp(k1y:kmy,l) = 0.0;
        haxpp(k1y:kmy,l) = 0.0;
        th   (k1y:kmy,l) = 0.0;
        th1  (k1y:kmy,l) = 0.0;
        CYCLE
        end if

         ! ... Horizontal diffusion
         IF ( ihd == 1 ) THEN ! Constant
           DO k = k1y, kmy
             haypp(k,l)=Ay0*MIN(hvpp(k,lNC(l)),hvpp(k,l))
             haxpp(k,l)=Ax0*MIN(hvpp(k,lEC(l)),hvpp(k,l))
           ENDDO
         ELSEIF ( ihd > 1) THEN ! Smagorinsky
           DO k = k1y, kmy
             haypp(k,l)=kh(k,lNC(l))*MIN(hvpp(k,lNC(l)),hvpp(k,l))
             haxpp(k,l)=(kh(k,lEC(    l )) + &
                         kh(k,        l )  + &
                         kh(k,lNC(    l )) + &
                         kh(k,lEC(lNC(l))))* &
                         0.25 *MIN(hvpp(k,lEC(l)),hvpp(k,l))
           ENDDO
         ENDIF

         !.....Calculate weighting arrays for vertical advection
         DO k = k1y, kmy
            th (k,l) = hvp(k-1,l)/(hvp(k-1,l)+hvp(k,l))
            th1(k,l) = 1.-th(k,l)
         ENDDO

         !.....Set th=th1 at the free surface & bottom
         th (k1y  ,l) = 0.0;
         th1(k1y  ,l) = 0.0;
         th (kmy+1,l) = 0.5;
         th1(kmy+1,l) = 0.5;

      END DO


      !......Calulate the explicit terms by sweeping over interior v-pts.....

      DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

         l = id_column(liter)

         ex(:,l) = 0.0

         ! ... Map 2D-l into 3D-(i,j) indexes
         i = l2i(l); j = l2j(l);

         ! ... Cycle if N-column is dry
         IF ( .NOT. mask2d(i,j+1) ) CYCLE

         ! Compute the layer number for top & bottom wet v-pt
         kmy = MIN(kmz(lNC(l)), kmz(l))
         k1y =                 k1v(l)

         ! Compute explicit term
         DO k = k1y,kmy

            ! ... For v-layers connecting wett & dry cells neglect
            !     contribution from advective, coriolis & diffusion
            IF ( hp(k,l) <= ZERO .OR. hp(k,lNC(l)) <= ZERO) THEN
               ex(k,l) = vhpp(k,l)
               CYCLE
            ENDIF

            !.....Coriolis.....
            cory = 0.25 * f * (uhp(k,     lNC(l) ) + uhp(k,    l )       &
                  &           +uhp(k, lWC(lNC(l))) + uhp(k,lWC(l)))

            !.....Advection
            uE = uhp(k,    lNC(l) ) +uhp(k,    l );
            uW = uhp(k,lWC(lNC(l))) +uhp(k,lWC(l));
            vN = vhp(k,    lNC(l) ) +vhp(k,    l );
            vS = vhp(k,        l  ) +vhp(k,lSC(l));
            wU = wp (k,    lNC(l) ) +wp (k,    l ); IF ( k == k1y ) wU = 0.0;
            wD = wp (k+1,  lNC(l) ) +wp (k+1  ,l ); IF ( k == kmy ) wD = 0.0;

            SELECT CASE ( itrmom)

            CASE (1)  ! Centered differences using th & th1 factors

              scE = vp(k,lEC(l)) + vp(k,    l )
              scW = vp(k,lWC(l)) + vp(k,    l )
              scN = vp(k,lNC(l)) + vp(k,    l )
              scS = vp(k,lSC(l)) + vp(k,    l )
              advy = (uE * scE - uW * scW ) / fourdx +       &
                     (vN * scN - vS * scS ) / fourdy
              advy = advy  +                                 &
                     (wU*(th (k  ,l)* vp(k  ,l)  +           &
                          th1(k  ,l)* vp(k-1,l)) -           &
                      wD*(th (k+1,l)* vp(k+1,l)  +           &
                          th1(k+1,l)* vp(k  ,l)) ) / 2.

            CASE(2) ! Upwinding

              advy = ( (uE+ABS(uE)) * vpp(k,    l ) +          &
                  &    (uE-ABS(uE)) * vpp(k,lEC(l)) -          &
                  &    (uW+ABS(uW)) * vpp(k,lWC(l)) -          &
                  &    (uW-ABS(uW)) * vpp(k,    l ) ) / fourdx &
                  & +( (vN+ABS(vN)) * vpp(k,    l ) +          &
                  &    (vN-ABS(vN)) * vpp(k,lNC(l)) -          &
                  &    (vS+ABS(vS)) * vpp(k,lSC(l)) -          &
                  &    (vS-ABS(vS)) * vpp(k,    l ) ) / fourdy &
                  & +( (wU+ABS(wU)) * vpp(k  ,l) +             &
                       (wU-ABS(wU)) * vpp(k-1,l)) / 4.         &
                    -( (wD+ABS(wD)) * vpp(k+1,l) +             &
                       (wD-ABS(wD)) * vpp(k  ,l)) / 4.

            CASE (3)  ! Centered differences - avoid computation of th1 and th

              scE = vp(k,lEC(l)) + vp(k,    l )
              scW = vp(k,lWC(l)) + vp(k,    l )
              scN = vp(k,lNC(l)) + vp(k,    l )
              scS = vp(k,lSC(l)) + vp(k,    l )
              scU = (vp(k  ,l)*hvp(k  ,l)+                      &
                     vp(k-1,l)*hvp(k-1,l))/                     &
                    (hvp(k ,l)+hvp(k-1,l))
              scD = (vp(k  ,l)*hvp(k  ,l)+                      &
                     vp(k+1,l)*hvp(k+1,l))/                     &
                    (hvp(k ,l)+hvp(k+1,l))
              advy = (uE * scE - uW * scW ) / fourdx +          &
                     (vN * scN - vS * scS ) / fourdy +          &
                     (wU * scU - wD * scD ) / 2.

            CASE(4) ! Upwinding only for horizontal advection

              advy = ( (uE+ABS(uE))* vpp(k,    l ) +           &
                  &    (uE-ABS(uE))* vpp(k,lEC(l)) -           &
                  &    (uW+ABS(uW))* vpp(k,lWC(l)) -           &
                  &    (uW-ABS(uW))* vpp(k,    l ) ) / fourdx  &
                  & +( (vN+ABS(vN))* vpp(k,    l ) +           &
                  &    (vN-ABS(vN))* vpp(k,lNC(l)) -           &
                  &    (vS+ABS(vS))* vpp(k,lSC(l)) -           &
                  &    (vS-ABS(vS))* vpp(k,    l ) ) / fourdy
              advy = advy +                                    &
                     (wU*(th (k  ,l)* vp(k  ,l)  +             &
                          th1(k  ,l)* vp(k-1,l)) -             &
                      wD*(th (k+1,l)* vp(k+1,l)  +             &
                          th1(k+1,l)* vp(k  ,l)) ) / 2.

            END SELECT

            !.....Horizontal diffusion.....
            IF ( ihd == 0) THEN
              hdy = 0.0E0
            ELSEIF ( ihd == 1) THEN
              hdy  =   (haypp(k,    l )*(vpp(k,lNC(l))-vpp(k,    l ))              &
                   & -  haypp(k,lSC(l))*(vpp(k,    l )-vpp(k,lSC(l))))/dydy        &
                   & + (haxpp(k,    l )*(vpp(k,lEC(l))-vpp(k,    l ))              &
                   & -  haxpp(k,lWC(l))*(vpp(k,    l )-vpp(k,lWC(l))))/dxdx
            ELSEIF ( ihd >  1) THEN

              hdy  = 2.* (haypp(k,   (l))*(vpp(k,lNC(    l ))-vpp(k,    l ))       &
                 &     -  haypp(k,lSC(l))*(vpp(k,        l  )-vpp(k,lSC(l))))/dydy &
                 &   +   (haxpp(k,    l )*(vpp(k,lEC(l     ))-vpp(k,    l ))       &
                 &     -  haxpp(k,lWC(l))*(vpp(k,        l  )-vpp(k,lWC(l))))/dxdx &
                 &   +   (haxpp(k,    l )*(upp(k,lNC(    l ))-upp(k,    l ))       &
                 &     -  haxpp(k,lWC(l))*(upp(k,lNC(lWC(l)))-upp(k,lWC(l))))/dxdy
            ENDIF

            ! ... Adjust terms to account for boundary conditions - HardCoded for NorthDelta Study
            !IF((                           j>=  jm1 - 20          )   .OR. &
            !   (i<=591               .AND. j<=  60                )   .OR. &
            !   (i<=20                                             )   .OR. &
            !   (i>=990 .AND. i<=1010 .AND. j>=  74 .AND. j <=  81 )   .OR. &
            !   (                           j<=  20                )   .OR. &
            !   (i>=678 .AND. i<= 682 .AND. j>= 100 .AND. j <= 115 )) THEN
            !   cory = 0.0
            !   advy = 0.0
            !   hdy  = 2.*hdy
            !ENDIF

            !IF ((i<=86                 .AND. j>991                 )   .OR. & !CINTIA (SUTTER SLOUGH)
            !    (i>=613                .AND. j<=27                 )   .OR. & !CINTIA (SUTTER GEO)
            !    (i>679                 .AND. j<=114                ))  THEN  !CINTIA (SUTTER DCC)

            !  !cory = 0.0
            !  !advy = 0.0
            !  hdy = 2.*hdy
            !ENDIF

            ! ... Needed to keep simulations stable near the boundaries - Cayuga
            !IF( i >=  im1 - 20) THEN;
            !  hdy   = 4.*hdy
            !  advy  = 0.0
            !  cory  = 0.
            !ENDIF
            !IF(j >= 297 .AND. j <= 301) THEN
            !  hdy   = 4.*hdy
            !  advy  = 0.0
            !  cory  = 0.
            !ENDIF
            !IF(j >= 2 .AND. j <= 6) THEN
            !  hdy   = 4.*hdy
            !  advy  = 0.0
            !  cory  = 0.
            !ENDIF
            !IF(i >= 2 .AND. i <= 6) THEN
            !  hdy   = 4.*hdy
            !  advy  = 0.0
            !  cory  = 0.
            !ENDIF


            !IF(j == 23) THEN !13
            !  hdy   = 4.*hdy
            !  advy  = 0.0
            !  cory  = 0.
            !ENDIF
            !IF(j == 2) THEN
            !  hdy   = 4.*hdy
            !  advy  = 0.0
            !  cory  = 0.
            !ENDIF
            !IF(i == 2) THEN
            !  hdy   = 4.*hdy
            !  advy  = 0.0
            !  cory  = 0.
            !ENDIF

            ! ... Needed to keep simulations stable near the boundaries - Beznar Cola
            !IF( i >=  im1 - 20) THEN;
            !  hdy   = 10.*hdy
            !  advy  = 0.0
            !  cory  = 0.0
            !ENDIF


            !.....Final explicit term.....
            ex(k,l) = vhpp(k,l) - twodt1*(advy*iadv+cory-hdy)

         END DO
      END DO


      ! ... Recalculate ex for near bdry. cells
      IF (nopen > 0 ) CALL MODexmom4openBCY

   END SELECT

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_exmom = t_exmom + (etime - btime)

END SUBROUTINE exmom

!***********************************************************************
SUBROUTINE matmom ( ieq, t_matmom2,Bstart, Bend, Bex,Beagx,Bearx,Bagx,Barx,Beagy,Beary,Bagy,Bary,uairB,vairB,cdwB,bclncxB,hupdrhoB )
!***********************************************************************
!
!  Purpose:   To define the matrices for the momentum equations.
!
!  Algorithm: The x-momentum equations at each vertical array of
!             u-pts are first expressed in the compact matrix form
!
!               [aa] [uh] = [gg] - g*dt/dx*rho*(s(i+1,j)-s(i,j))*[hh]
!
!             by defining the three matrices  [hh],  [gg], and  [aa].
!             (Because the  [aa]  matrix is tridiagonal, only the
!             diagonals are stored.) The above system of equations is
!             then rearranged into the form
!
!               [uh] = [ag] - g*dt/dx*rho*(s(i+1,j)-s(i,j))*[ar]
!
!             by indirect solution using the tridiagonal solver  trid.
!             The matrices  [ag]  and  [ar]  are the output from the  trid
!             subroutine along with their summation over the depth at
!             each horizontal node point,  [eag]  and [ear]. The matrices
!             for the x-momentum eq are stored in the fortran arrays
!             agx,  arx,  eagx, and  earx. Everything is similar for
!             the y-momentum eq.
!
!  Dummy argument:
!  ieq    = Parameter indicating whether the matrices for the
!           x  or  y  momentum equation are to be evaluated.
!           (1=x-momentum, 2=y-momentum)
!
!  23/04/2008   F.J. Rueda    Recompute baroclinic term at bottom
!  23/04/2008   F.J. Rueda    Do not set to zero baroclinic term at
!                             wett/dry cells near the surface

!-----------------------------------------------------------------------

   !.....Argument.....
   INTEGER, INTENT(IN) :: ieq, Bstart, Bend
   REAL, INTENT(INOUT) :: t_matmom2
   REAL, DIMENSION(1:km1), INTENT(INOUT) :: bclncxB,hupdrhoB

   REAL, DIMENSION (1:km1,Bstart:Bend+1), INTENT(INOUT) :: Bex,Bagx,Barx,Bagy,Bary
   REAL, DIMENSION (Bstart:Bend+1), INTENT(INOUT) :: Beagx,Bearx,Beagy,Beary,uairB,vairB,cdwB

   !.....Local variables.....
   INTEGER :: i, j, k, l, kmx, kmy, k1x, k1y, nwlayers, inn, liter, innH
   REAL    :: twodt1, wsx0, wsy0, tausx, taubx, tausy, tauby,      &
              hpmin, usurf, vsurf,Uhvalue, Usource, Vsource,       &
              Vhvalue, cwx, cwy,lol
   REAL :: aaux,aaux2,aaux3
   REAL, DIMENSION(km1) ::  vdiffx, vdiffy,         &
         &  deltaz, Avxdudz, Avydvdz, Avx, Avy,   &
         & rhopx, rhopy, gg, hh, ar, ag
   REAL,DIMENSION(km1) :: hvpdrho,bclncy,bclncx,hupdrho
   REAL, DIMENSION(3,km1) :: aa
   REAL, PARAMETER :: rhop0 = 1000.

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)
   !print *,"mattt:",omp_get_thread_num()
   !.....Constants.....
   twodt1 = twodt*tz;

   SELECT CASE ( ieq )

   !                -----X-momentum equation-----
   CASE (1)



      !.....Loop over interior u-pts .....

      DO liter = lhiWCE(omp_get_thread_num ( )+1), lhfCE(omp_get_thread_num ( )+1)
            If(liter == 0) CYCLE

            l = id_columnCE(liter)

            ! ... Map 2D-l into 3D-(i,j) indexes - FJR uncomment - needed for pss - CHECK mario!!!!
            i = l2i(l);
            j = l2j(l);

            ! ... Skip for dry u-pts
            !IF (.NOT.mask(c2l(lEC(l)))) THEN
            !Beagx(l) = 0.0;
            !Bearx(l) = 0.0;
            !CYCLE
            !END IF

            ! ... Compute the layer number for top & bottom wet u-pt
            kmx = MIN(kmz(lEC(l)), kmz(l))
            k1x =                  k1u(l)
            nwlayers = (kmx-k1x) + 1
            IF(nwlayers < 1) CYCLE


            ! ... Compute eddy viscosity at interfaces vertically between u-pts
            Avx = 0.0
            DO k = k1x+1,kmx
              Avx(k) = 0.5 * (Av(k,lEC(l))+Av(k,l))
            ENDDO

            ! ... Define average layer density at u-pt (in kg/m**3) ...........
            rhopx(k1x:kmx) = 1000. ! Neglect vertical density variations

            ! ... Compute explicit portion of water surface slope term ........
            wsx0 = rhopx(k1x) * gdtdx * (spp(lEC(l)) - spp(l))

            ! ... Compute baroclinic term .....................................
            SELECT CASE (ibc)
            CASE (0)     ! No baroclinic term
               bclncx(k1x:kmx) = 0.0
            CASE (1:)
              DO k = k1x, kmx
                hupdrho(k) = gdtdx*hup(k,l)*(rhop(k,lEC(l))-rhop(k,l))
                !IF(hp(k,l) <= ZERO .OR. hp(k,lEC(l)) <= ZERO) hupdrho(k) = 0.0
              ENDDO
              bclncx(k1x) = hupdrho(k1x)
              IF (kmx > k1x) THEN   ! Two or more wet layers
                ! Recompute bottom layer baroclinic term along horizontal plane
                CALL bclnc_km (l, kmx, 1, hupdrho(kmx) )

                DO k = k1x+1, kmx
                   aaux = hupdrho(k-1) + hupdrho(k)
                   aaux2 = bclncx(k-1) + aaux
                   bclncx(k) = aaux2
               END DO
              END IF
            END SELECT

            ! ... Compute explicit portion of vertical diffusion term ........
            SELECT CASE (nwlayers)
            CASE (1)    ! Single wet layer
              vdiffx(k1x) = 0.0
            CASE (2:)   ! Two or more wet layers (hupp->hup)
              DO k = k1x+1,kmx
                deltaz(k) = hup(k-1,l) + hup(k,l)
                Avxdudz(k)= Avx(k)*(upp(k-1,l)-upp(k,l))/deltaz(k)
              ENDDO
              Avxdudz(k1x)      = 0.0  ! Set value at free surface to zero
              Avxdudz(kmx+1)    = 0.0  ! Set value at bottom boundary to zero
              vdiffx(k1x:kmx) = twodt*(Avxdudz(k1x:kmx) - Avxdudz(k1x+1:kmx+1)) &
                                *2.*(1.-theta) ! theta accounts for semi-implicitness
            END SELECT

            !.....Form  [hh]  matrix...........................................
            hh   (k1x:kmx) = hup(k1x:kmx,l)/rhopx(k1x:kmx)

            !.....Form [gg]  matrix............................................
            gg(k1x:kmx) = ex(k1x:kmx,l)                                       &
                         -hh(k1x:kmx  ) * (bclncx(k1x:kmx)+wsx0) * tz         &
                         +vdiffx(k1x:kmx)                        * tz

            !.....Form [aa]  matrix............................................
            SELECT CASE (nwlayers)
            CASE (1)    ! Single wet layer
              aa(2,k1x) = 1.0
            CASE (2:)   ! Two or more wet layers (hu->hup)
              ! Define upper diagonal terms
              aa(3,k1x:kmx-1)= -twodt1*Avx(k1x+1:kmx)/(hup(k1x+1:kmx,l)*      &
                                (hup(k1x:kmx-1,l)+hup(k1x+1:kmx,l)))*2.*theta
              aa(3,kmx)      =  0.0
              ! Define lower diagonal terms
              aa(1,k1x+1:kmx)= -twodt1*Avx(k1x+1:kmx)/(hup(k1x:kmx-1,l)*      &
                                (hup(k1x:kmx-1,l)+hup(k1x+1:kmx,l)))*2.*theta
              aa(1,k1x)      =  0.0
              ! Define center diagonal terms
              aa(2,k1x:kmx)  =  1.0                                           &
                        -(hup(k1x-1:kmx-1,l)/hup(k1x:kmx,l))*aa(1,k1x:kmx)    &
                        -(hup(k1x+1:kmx+1,l)/hup(k1x:kmx,l))*aa(3,k1x:kmx)
            END SELECT

            ! ... Top boundary conditions......................................
            ! a. Form wind stress term
            usurf = up(k1x,l)
            vsurf =(vp(k1v(l  ),        l  ) +                            &
                    vp(k1v(lSC(l)),    lSC(l) ) +                            &
                    vp(k1v(lEC(l)),    lEC(l) ) +                            &
                    vp(k1v(lSC(lEC(l))),lSC(lEC(l))))/4.
            cwx = cdw(l)*rhoair*SQRT((vair(l)-vsurf)**2.+ &
                                       (uair(l)-usurf)**2.)
            ! b. Modify [gg] matrix
            tausx     = cwx/rhopx(k1x)*uair(l)
            gg(  k1x) = gg(  k1x) + tausx*twodt1
            ! c. Modify [aa] matrix
            tausx     = cwx/rhopx(k1x)/hup(k1x,l)
            aa(2,k1x) = aa(2,k1x) + tausx*twodt1

            ! ... Bottom boundary conditions...................................
            ! a. Form bottom stress term
            taubx = cd*SQRT((uhpp(kmx,     l)*uhpp(kmx,        l)) +          &
                  &  ((0.25*(vhpp(kmx,lEC(l))+vhpp(kmx,        l)             &
                  &         +vhpp(kmx,lSC(l))+vhpp(kmx,lSC(lEC(l)))))**2))    &
                  &                /(hup(kmx,l)*hup(kmx,l))
            ! b. Modify [aa] matrix
            aa(2,kmx)   = aa(2,kmx) + taubx*twodt1

            ! .... Point sources and sinks ....................................
            IF ( iopssH(omp_get_thread_num ( )+1) > 0) THEN

              DO innH = 1, iopssH(omp_get_thread_num ( )+1)
                inn = ioph2iop(innH,omp_get_thread_num ( )+1)
                IF (i == ipss(inn)     .AND. &
                    j == jpss(inn)     ) THEN
                  DO k = k1x,kmx
                    IF (ABS(Qpss(k,inn))<1.E-10) CYCLE
                    ! ... Strength of Source - here it is assumed that
                    !     only half of the flow shows up in the control volume
                    !     used in the momentum equation -> factor 2 below
                    Usource = ABS(Qpss(k,inn))/(dx*dy*hup(k,l))*twodt1/2.
                    IF(ptype(iodev(inn)) == -2) Usource = 1.E2
                    ! ... Velocity of the source in E direction (positive
                    !     towards east if a source; negative or towards west if
                    !     a sink) - idetr = 1 by default;
                    Uhvalue = (Qpss(k,inn)*uEpss(iodev(inn))/dy)*idetr(iodev(inn)) !mod ACC oct11
                    aa(2,k) = aa(2,k) + Usource
                    gg(  k) = gg(  k) + Usource * Uhvalue
                  ENDDO
				ENDIF
                IF (i == ipss(inn)-1    .AND. &
                    j == jpss(inn)    ) THEN
                  DO k = k1x,kmx
                    IF (ABS(Qpss(k,inn))<1.E-10) CYCLE
                    ! ... Strength of Source - here it is assumed that
                    !     only half of the flow shows up in the control volume
                    !     used in the momentum equation -> factor 2 below
                    Usource = ABS(Qpss(k,inn))/(dx*dy*hup(k,l))*twodt1/2.
                    IF(ptype(iodev(inn)) == -2) Usource = 1.E2
                    ! ... Velocity of the source in N direction (positive
                    !     towards north if a source; negative or towards south if
                    !     a sink) - idetr = 1 by default;
                    Uhvalue = -(Qpss(k,inn)*uWpss(iodev(inn))/dy)*idetr(iodev(inn)) !mod ACC oct11
                    aa(2,k) =   aa(2,k) + Usource
                    gg(  k) =   gg(  k) + Usource * Uhvalue
                  ENDDO
                ENDIF
              ENDDO
            ENDIF


            !.....Solve tridiagonal system for  [ag]  and  [ar]  arrays........
            SELECT CASE (nwlayers)
            CASE (1)    ! Single wet layer
              ag(k1x) = gg(k1x)/aa(2,k1x)
              ar(k1x) = hh(k1x)/aa(2,k1x)
            CASE (2:)   ! Two or more wet layers
              CALL trid ( aa, gg, hh, ag, ar, k1x, kmx, kmx+1, nwlayers )
            END SELECT

            !.....Save  [ag]  and  [ar]  arrays and sum them over
            !     depth for use in solution of continuity equation.............
            Bagx(k1x:kmx,l) = ag(k1x:kmx)
            Barx(k1x:kmx,l) = ar(k1x:kmx)
            Beagx(l) = SUM(ag(k1x:kmx))
            Bearx(l) = SUM(ar(k1x:kmx))

      !.....End loop over u-pts.....
      END DO



   !                -----Y-momentum equation-----
   CASE (2)



      !.....Loop over interior v-pts
      !print *,"mattt2:",omp_get_thread_num()
      DO liter = lhiCN(omp_get_thread_num ( )+1), lhfCN(omp_get_thread_num ( )+1)


            l = id_columnCN(liter)

    	    ! ... Map 2D-l into 3D-(i,j) indexes - uncomment - needed
            i = l2i(l);
            j = l2j(l);
            !print *,"mcol:",l,":",omp_get_thread_num()
            ! ... Skip if North Column is dry
            !IF (.NOT. mask2d(i,j+1)) THEN
            !Beagy(l) = 0.0;
            !Beary(l) = 0.0;
            !CYCLE
            !END IF
            ! ... Compute layer numbers of wet v-pt
            kmy = MIN(kmz(lNC(l)), kmz(l))
            k1y =                 k1v(l)
            nwlayers = (kmy-k1y) + 1
            IF(nwlayers < 1) CYCLE

            ! .... Compute eddy viscosity at interfaces between v-pts)
            Avy = 0.0
            DO k = k1y, kmy
              Avy(k) = 0.5*(Av(k,lNC(l))+Av(k,l))
            ENDDO

            ! .... Define average layer density at v-pts (in kg/m**3) .........
            rhopy(k1y:kmy) = 1000. ! Neglect vertical density variations

            !.....Compute explicit part of water surface slope term ...........
            wsy0 = rhopy(k1y) *  gdtdy  *(spp(lNC(l)) - spp(l))

            !.... Compute baroclinic term .....................................
            SELECT CASE (ibc)
            CASE (0)    ! No baroclinic term
              bclncy(k1:kmy) = 0.0
            CASE (1:)
              DO k = k1y, kmy
                hvpdrho(k) = gdtdy*hvp(k,l)*(rhop(k,lNC(l))-rhop(k,l))
                ! IF(hp(k,l)<=ZERO .OR. hp(k,lNC(l))<=ZERO) hvpdrho(k) = 0.0
              ENDDO
              bclncy(k1y) = hvpdrho(k1y)
              IF (kmy > k1y) THEN   ! Two or more wet layers
                ! Recompute bottom layer baroclinic term along horizontal plane
                CALL bclnc_km (l, kmy, 2, hvpdrho(kmy) )
                DO k = k1y+1, kmy
                  bclncy(k) = bclncy(k-1) + hvpdrho(k-1) + hvpdrho(k)
                END DO
              END IF
            END SELECT

            ! ... Compute explicit portion of vertical diffusion term ........
            SELECT CASE (nwlayers)
            CASE (1)    ! Single wet layer
              vdiffy(k1y) = 0.0
            CASE (2:)   ! Two or more wet layers (hvpp->hvp)
              DO k = k1y+1 , kmy
                deltaz(k) = hvp(k-1,l) + hvp(k,l)
                Avydvdz(k)= Avy(k)*(vpp(k-1,l)-vpp(k,l))/deltaz(k)
              ENDDO
              Avydvdz(k1y)      = 0.0  ! Set value at free surface to zero
              Avydvdz(kmy+1)    = 0.0  ! Set value at bottom boundary to zero
              vdiffy(k1y:kmy) = twodt*(Avydvdz(k1y:kmy) - Avydvdz(k1y+1:kmy+1)) &
                                *2.*(1.-theta) ! This factor accounts for semi-implicitness
            END SELECT


            !.....Form  [hh]  matrix...........................................
            hh   (k1y:kmy) = hvp(k1y:kmy,l)/rhopy(k1y:kmy)

            !.....Compute  [gg]  matrix .......................................
            gg(k1y:kmy) = ex(k1y:kmy,l)                                       &
                        - hh(k1y:kmy  ) * (bclncy(k1y:kmy)+wsy0) * tz         &
                        + vdiffy(k1y:kmy)                        * tz

            !.....Form  [aa]  matrix...........................................
            SELECT CASE (nwlayers)
            CASE (1)    ! Single wet layer
              aa(2,k1y) = 1.0
            CASE (2:)   ! Two or more wet layers (hv->hvp)
              ! Define upper diagonal terms
              aa(3,k1y:kmy-1)=-twodt1*Avy(k1y+1:kmy)/(hvp(k1y+1:kmy,l)*       &
                             & (hvp(k1y:kmy-1,l)+hvp(k1y+1:kmy,l)))*2.*theta
              aa(3,kmy)      = 0.0
              ! Define lower diagonal terms
              aa(1,k1y+1:kmy)=-twodt1*Avy(k1y+1:kmy)/(hvp(k1y:kmy-1,l)*       &
                             & (hvp(k1y:kmy-1,l)+hvp(k1y+1:kmy,l)))*2.*theta
              aa(1,k1y)      = 0.0
              ! Define center diagonal terms
              aa(2,k1y:kmy)  = 1.0                                            &
                          -(hvp(k1y-1:kmy-1,l)/hvp(k1y:kmy,l))*aa(1,k1y:kmy)  &
                          -(hvp(k1y+1:kmy+1,l)/hvp(k1y:kmy,l))*aa(3,k1y:kmy)
            END SELECT

            ! ... Top boundary conditions .....................................
            ! a. Form wind stress term
            vsurf = vp(k1y,l)
            usurf =(up(k1u(l),        l  ) +                            &
                    up(k1u(lWC(l)),    lWC(l) ) +                            &
                    up(k1u(lNC(l)),    lNC(l) ) +                            &
                    up(k1u(lWC(lNC(l))),lWC(lNC(l))))/4.
            cwy = cdw(l)*rhoair*SQRT((uair(l)-usurf)**2.+ &
                                       (vair(l)-vsurf)**2.)
            ! b. Modify [gg] matrix
            tausy = cwy/rhopy(k1y)*vair(l)
            gg(k1y)   = gg(  k1y) + tausy*twodt1
            ! c. Modify [aa] matrix
            tausy = cwy/rhopy(k1y)/hvp(k1y,l)
            aa(2,k1y) = aa(2,k1y) + tausy*twodt1

            ! ... Bottom boundary conditions ..................................
            ! a. Form bottom stress term
            tauby = cd*SQRT((vhpp(kmy,l)*vhpp(kmy,l)) +                       &
                  &  ((0.25*(uhpp(kmy,lNC(l))+uhpp(kmy,lWC(lNC(l)))           &
                  &         +uhpp(kmy,lWC(l))+uhpp(kmy,l)))**2))              &
                  &         /(hvp(kmy,l)*hvp(kmy,l))
            ! b. Modify [aa] matrix
            aa(2,kmy) = aa(2,kmy) + tauby*twodt1

            ! .... Point sources and sinks ....................................
            IF ( iopssH(omp_get_thread_num ( )+1) > 0) THEN

              DO innH = 1, iopssH(omp_get_thread_num ( )+1)
                inn = ioph2iop(innH,omp_get_thread_num ( )+1)
                IF (i == ipss(inn)     .AND. &
                    j == jpss(inn)     ) THEN
                  DO k = k1y,kmy
                    IF (ABS(Qpss(k,inn))<1.E-10) CYCLE
                    ! ... Strength of Source - here it is assumed that
                    !     only half of the flow shows up in the control volume
                    !     used in the momentum equation -> factor 2 below
                    Vsource = ABS(Qpss(k,inn))/(dx*dy*hvp(k,l))*twodt1/2.
                    IF(ptype(iodev(inn)) == -2) Vsource = 1.E2
                    ! ... Velocity of the source in N direction (positive
                    !     towards north if a source; negative or towards south if
                    !     a sink) - idetr = 1 by default;
                    Vhvalue = (Qpss(k,inn)*vNpss(iodev(inn))/dx)*idetr(iodev(inn)) !mod ACC oct11
                    aa(2,k) = aa(2,k) + Vsource
                    gg(  k) = gg(  k) + Vsource * Vhvalue
                  ENDDO
                ENDIF
                IF (i == ipss(inn)      .AND. &
                    j == jpss(inn)-1    ) THEN
                  DO k = k1y,kmy
                    IF (ABS(Qpss(k,inn))<1.E-10) CYCLE
                    ! ... Strength of Source - here it is assumed that
                    !     only half of the flow shows up in the control volume
                    !     used in the momentum equation -> factor 2 below
                    Vsource = ABS(Qpss(k,inn))/(dx*dy*hvp(k,l))*twodt1/2.
                    IF(ptype(iodev(inn)) == -2) Vsource = 1.E2
                    ! ... Velocity of the source in S direction (negative
                    !     towards south if a source; positive or towards north if
                    !     a sink) - idetr = 1 by default;
                    Vhvalue = -(Qpss(k,inn)*vSpss(iodev(inn))/dx)*idetr(iodev(inn)) !mod ACC oct11
                    aa(2,k) = aa(2,k) + Vsource
                    gg(  k) = gg(  k) + Vsource * Vhvalue
                  ENDDO
                ENDIF
              ENDDO

	ENDIF

            !.....Solve tridiagonal system for  [ag]  and  [ar]  arrays........
            SELECT CASE (nwlayers)
            CASE (1)    ! Single wet layer
              ag(k1y) = gg(k1y)/aa(2,k1y)
              ar(k1y) = hh(k1y)/aa(2,k1y)
            CASE (2:)   ! Two or more wet layers
              CALL trid ( aa, gg, hh, ag, ar, k1y, kmy, km1, nwlayers )
            END SELECT

            !.....Save  [ag]  and  [ar]  arrays and sum them over
            !     depth for use in solution of continuity equation..............
            Bagy(k1y:kmy,l) = ag(k1y:kmy)
            Bary(k1y:kmy,l) = ar(k1y:kmy)
            Beagy(l) = SUM(ag(k1y:kmy))
            Beary(l) = SUM(ar(k1y:kmy))

         !.....End loop over v-pts.....
         END DO


   END SELECT

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_matmom2 = t_matmom2 + (etime - btime)

END SUBROUTINE matmom

!***********************************************************************
SUBROUTINE bclnc_km (l, kb, ieq, huvpdrho )
!***********************************************************************
!
!  Purpose: To recompute the baroclinic term for a bottom layer of
!           variable depth using densities interpolated (or
!           extrapolated) to locations on a horizontal surface.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!  18/09/00          P.E. Smith        Original f90 code
!
!-----------------------------------------------------------------------

   !.....Arguments.....
   INTEGER, INTENT(IN   ) :: l, kb, ieq
   REAL   , INTENT(INOUT) :: huvpdrho

   !.....Local variables.....
   REAL, PARAMETER :: eps = EPSILON(0.0)
   REAL :: dz1, dz2, rhop_e, rhop_w, rhop_n, rhop_s
   REAL :: hup_plus, hup_minus, hvp_plus, hvp_minus
   LOGICAL :: condition_e, condition_w, condition_n, condition_s

   !.....Choose x- or y-momentum equation...............................
   SELECT CASE ( ieq )


   !            -----Recompute the x-momentum term-----



   CASE (1)

      !.....Check if the depth at the press-pt on the east side of the
      !     control volume is not equal to the depth at the u-pt.....
      hup_plus = hup(kb,l)+eps;  hup_minus = hup(kb,l)-eps
      IF((hp(kb,lEC(l)) > hup_plus) .OR. (hp(kb,lEC(l)) < hup_minus)) THEN
         condition_e = .TRUE.
         dz1 = 0.5*(hp(kb-1,lEC(l)) + hp (kb,lEC(l)))
         dz2 = 0.5*(hp(kb-1,lEC(l)) + hup(kb,    l ))
         ! Interpolate (or extrapolate) for the density at the
         ! horizontal level of the u-pt
         rhop_e=rhop(kb-1,lEC(l))+(dz2/dz1)*(rhop(kb,lEC(l))-rhop(kb-1,lEC(l)))
      ELSE
         condition_e = .FALSE.
         rhop_e = rhop(kb,lEC(l))
      END IF

      !.....Check if the depth at the press-pt on the west side of the
      !     control volume is not equal to the depth at the u-pt.....
      IF((hp(kb,l) > hup_plus) .OR. (hp(kb,l) < hup_minus)) THEN
         condition_w = .TRUE.
         dz1 = 0.5*(hp(kb-1,l) + hp (kb,l))
         dz2 = 0.5*(hp(kb-1,l) + hup(kb,l))
         ! Interpolate (or extrapolate) for the density at the
         ! horizontal level of the u-pt

         rhop_w=rhop(kb-1,l)+(dz2/dz1)*(rhop(kb,l)-rhop(kb-1,l))
      ELSE
         condition_w = .FALSE.
         rhop_w = rhop(kb,l)
      END IF

      !.....If necessary, recompute the x-direction
      !     baroclinic term on a horizontal plane.....
      IF ( condition_e .OR. condition_w )  THEN
         huvpdrho = gdtdx*hup(kb,l)*(rhop_e-rhop_w)
      END IF


   !            -----Recompute the y-momentum term-----


   CASE (2)

      !.....Check if the depth at the press-pt on the north side of
      !     the control volume is not equal to the depth at the v-pt.....
      hvp_plus = hvp(kb,l)+eps;  hvp_minus = hvp(kb,l)-eps
      IF((hp(kb,lNC(l)) > hvp_plus) .OR. (hp(kb,lNC(l)) < hvp_minus)) THEN
         condition_n = .TRUE.
         dz1 = 0.5*(hp(kb-1,lNC(l)) + hp (kb,lNC(l)))
         dz2 = 0.5*(hp(kb-1,lNC(l)) + hvp(kb,    l ))
         ! Interpolate (or extrapolate) for the density at the
         ! horizontal level of the v-pt
         rhop_n=rhop(kb-1,lNC(l))+(dz2/dz1)*(rhop(kb,lNC(l))-rhop(kb-1,lNC(l)))
      ELSE
         condition_n = .FALSE.
         rhop_n = rhop(kb,lNC(l))
      END IF

      !.....Check if the depth at the press-pt on the south side of
      !     the control volume is not equal to the depth at the v-pt.....
      IF((hp(kb,l) > hvp_plus) .OR. (hp(kb,l) < hvp_minus)) THEN
         condition_s = .TRUE.
         dz1 = 0.5*(hp(kb-1,l) + hp (kb,l))
         dz2 = 0.5*(hp(kb-1,l) + hvp(kb,l))
         ! Interpolate (or extrapolate) for the density at the
         ! horizontal level of the v-pt
         rhop_s=rhop(kb-1,l)+(dz2/dz1)*(rhop(kb,l)-rhop(kb-1,l))
      ELSE
         condition_s = .FALSE.
         rhop_s = rhop(kb,l)
      END IF

      !.....If necessary, recompute the y-direction
      !     baroclinic term on a horizontal plane.....
      IF ( condition_n .OR. condition_s )  THEN
         huvpdrho = gdtdy*hvp(kb,l)*(rhop_n-rhop_s)
      END IF

   END SELECT

END SUBROUTINE bclnc_km

!***********************************************************************
SUBROUTINE matcon(t_matcon2,Bstart,Bend,lWCH,lSCH,Beagx,Bearx,Beagy,Beary,Bsx,Bsy,Bdd,Bqq,Brr)
!***********************************************************************
!
!  Purpose: To calculate the matrix coefficients for solving the
!           continuity equation for zeta.
!
!-----------------------------------------------------------------------


   REAL, INTENT(INOUT) :: t_matcon2
   INTEGER, INTENT(IN) :: Bstart,Bend
   REAL, DIMENSION (Bstart:Bend+1), INTENT(INOUT) :: Beagx,Bearx,Beagy,Beary,Bsx,Bsy,Bdd,Bqq,Brr
   INTEGER, DIMENSION(Bstart:Bend+1), INTENT(IN) :: lWCH,lSCH



   !.....Local variables.....
   REAL :: cx, cy, dt1, dtdx1, dtdy1, rho4sx, rho4sy
   INTEGER :: i, j, k, l, k1s, kms, inn, liter, innH

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)

   !.....Constants.....
   cx     = gdt2dx2*tz*tz;
   cy     = gdt2dy2*tz*tz;
   dtdx1  = dtdx*tz;
   dtdy1  = dtdy*tz;
   dt1    = dt*tz;
   rho4sx = 1000. ! Neglect density variations
   rho4sy = 1000. ! Neglect density variations

   !.....Calculate  [sx] matrix at u-pts & [sy] matrix at v-pts ....
   Bsx(Bend+1) = 0.0; Bsy(Bend+1) = 0.0;
   DO liter = lhiW(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
     If(liter == 0) CYCLE
     l = id_column(liter)
     ! ... Map 2D-l into 3D-(i,j) indexes
     i = l2i(l); j = l2j(l);
     ! ... u-pts
     IF (mask2d(i+1,j)) THEN
        Bsx(l)= cx * rho4sx * Bearx(l)
     ELSE
        Bsx(l)= 0.0
     ENDIF
   ENDDO
   DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
     l = id_column(liter)
     i = l2i(l); j = l2j(l);
     ! ... v-pts
     IF (mask2d(i,j+1)) THEN
        Bsy(l)= cy * rho4sy * Beary(l)
     ELSE
        Bsy(l)= 0.0
     ENDIF

   ENDDO

   !.....Calculate  [dd], [qq], and [rr]  matrices at zeta-pts.....
!   dd(lm1) = 0.0; qq(lm1) = 0.0; rr(lm1) = 1.0
   DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
     l = id_column(liter)

     ! ... Map 2D-l into 3D-(i,j) indexes
!     i = l2i(l); j = l2j(l);

     ! ... Top & bottom cells
     k1s = k1
     kms = kmz(l)

     ! ... Form matrices
     Bdd(l) = dt1*SUM((uhpp(k1s:kms,      l) -         &
             &         uhpp(k1s:kms,lWC(l))) /dx       &
             &        +(vhpp(k1s:kms,     l) -         &
             &          vhpp(k1s:kms,lSC(l)))/dy)
     Bqq(l) = spp(l) - (dtdx1)*(Beagx(l)-Beagx(lWCH(l)))   &
             &       - (dtdy1)*(Beagy(l)-Beagy(lSCH(l)))   &
             &       - Bdd(l)
     Brr(l) = 1 + Bsx(l) + Bsx(lWCH(l)) + Bsy(l) + Bsy(lSCH(l))

   END DO;

   ! .... Modify qq matrices to incorporate sources/sinks
   IF ( iopssH(omp_get_thread_num ( )+1) > 0 ) THEN
     DO innH = 1, iopssH(omp_get_thread_num ( )+1)
       inn = ioph2iop(innH,omp_get_thread_num ( )+1)
       IF ( ptype(iodev(inn)) > 0) CYCLE
       i = ipss(inn);
       j = jpss(inn);
       Bqq(ij2l(i,j)) = Bqq(ij2l(i,j)) +  SUM(Qpss(:,inn))/(dx*dy)*twodt*tz
     ENDDO
   ENDIF

   !DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
     !l = id_column(liter)
     !qq(l)=Bqq(l)
     !dd(l)=Bdd(l)
     !eagx(l)=Beagx(l)
     !eagy(l)=Beagy(l)
   !END DO

   !print *,"qq:",sum(qq(:))
   !print *,"spp:",sum(spp(:))
   !print *,"dd:",sum(dd(:))
   !print *,"eagx:",sum(eagx(:))
   !print *,"eagy:",sum(eagy(:))

   !.....Adjust qq & sx/sy matrices for open boundary conditions.....
   IF (nopen > 0 ) CALL MODqqddrr4openBC(Bstart,Bend,Bqq,Bsx,Bsy)

   !print *,"qqmod:",sum(qq(:))

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_matcon2 = t_matcon2 + (etime - btime)

END SUBROUTINE matcon

!***********************************************************************
SUBROUTINE SolverSparse(n,Bstart,Bend,lWCH,lSCH,Bsx,Bsy,Bqq,Brr,iter,istep,thrs)
!***********************************************************************
!
!  Purpose: To solve the system matrix for zeta using the
!           preconditioned conjugate gradient method. It uses
!           Storage format 1 (i.e. ELLPACK or iparm(12) = 1;
!           in this manner a considerable amount of time is saved
!           as we do not have to store an imxjm by imxjm matrix,
!           pentadiagonal but very large; instead we only store
!           lmxlm matrix; this is extremely useful in sparse
!           bathymetries such as in rivers. Each row in the matrix
!           has at most 5 non-zero elements, which are stored in coeffA;
!           the location of the coefficients in the matrix is stored
!           in jcoefA - see instructions.
!
!-----------------------------------------------------------------------

   INTEGER, INTENT(IN) :: Bstart,Bend,n,iter,istep
   REAL, INTENT(IN) :: thrs
   REAL, DIMENSION (Bstart:Bend+1), INTENT(IN) :: Bsx,Bsy,Bqq,Brr
   INTEGER, DIMENSION(Bstart:Bend+1), INTENT(IN) :: lWCH,lSCH

   !.....Local variables.....
   EXTERNAL mic1, jac1, cg, si
   INTEGER :: nw1a, inw1a, maxnz1a, i, j, m, &
              ier, nrA, ncA, istat, nfirstA=0, ios,aux_indice,liter
   INTEGER, SAVE:: i895 = 895
   CHARACTER(LEN=25):: solverfile="SolverAMODE.txt"

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)

   !.....Set matrix solution variables and workspace
   !     parameters on first entry into subroutine.....
   !IF ( nfirstA == 0 ) THEN
   !   nfirstA = 1
   !   ! Compute the number of rows of active cells in the grid
   !   nrA = jlast-jfirst + 1;
   !   ! Compute the number of columns of active cells in the grid
   !   ncA = ilast-ifirst + 1
   !   ! Compute number of equations to be solved
   !   ndimA = nrA*ncA ! NorthDelta
   !   ! Define the matrix bandwidth
   !   ibdwdA = nrA
   !   ! Define column width of 'coef' matrix
   !   maxnzA = mdimA
   !   ! Liberally estimate workspace requirements
   !   nwA = 20*ndimA;  inwA = 5*ndimA
   !   ! Open file
   !   OPEN (UNIT=i895, FILE=solverfile, IOSTAT=ios)
   !END IF
   IF(omp_get_thread_num ( )==0) THEN
   !.....Set matrix solution variables and workspace
   !     parameters on first entry into subroutine.....
   IF ( nfirstA == 0 ) THEN
      nfirstA = 1
      ! Compute number of equations to be solved
      ndimA = lm
      ! Define column width of 'coef' matrix
      maxnzA = mdimA
      ! Liberally estimate workspace requirements
      nwA = 20*ndimA;  inwA = 5*ndimA
      ! Open file
      !OPEN (UNIT=i895, FILE=solverfile, IOSTAT=ios)
   END IF


   !.....Reset workspace parameters.....
   nw1a = nwA+(2*ndimA); inw1a = inwA; maxnz1a = maxnzA

   !.....Allocate arrays.....
   ALLOCATE (wksp  (nw1A )      , &
          &  iwksp (inw1A)      , STAT=istat)
   IF (istat /= 0) CALL allocate_error ( istat, 23 )
   END IF

   !.....Define  coef, jcoef1  and  rhs  arrays in
   !     preparation for calling .....
   DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
     m = id_column(liter)
!     i = l2i(m);
!     j = l2j(m);
     coeffA (m,1) =  Brr(m)
     coeffA (m,2) = -Bsx(lWCH(m))
     coeffA (m,3) = -Bsy(lSCH(m))
     coeffA (m,4) = -Bsy(m)
     coeffA (m,5) = -Bsx(m)
     jcoefA (m,1) =  m
     jcoefA (m,2) =  lWC(m)
     jcoefA (m,3) =  lSC(m)
     jcoefA (m,4) =  lNC(m)
     jcoefA (m,5) =  lEC(m)
     IF(jcoefA(m,2)==lm1) jcoefA(m,2)=0;
     IF(jcoefA(m,3)==lm1) jcoefA(m,3)=0;
     IF(jcoefA(m,4)==lm1) jcoefA(m,4)=0;
     IF(jcoefA(m,5)==lm1) jcoefA(m,5)=0;
     rhs(m)       =  Bqq(m)
     ! initial guess at zeta
     zeta(m)      = sp(m)
   END DO

   !if(nopen > 0) CALL MODcoef4openBC

   !$omp barrier

   IF(omp_get_thread_num ( )==0) THEN
   !print *,"coeff1",sum(coeffA(:,1))
   !print *,"coeff4",sum(coeffA(:,4))
   !print *,"coeff5",sum(coeffA(:,5))
   !print *,"rhs",sum(rhs(:))
   !print *,"zeta",sum(zeta(:))
   !.....Set parameter defaults.....


   CALL dfault ( iparm, rparm )

   !.....Reset some default parameter values.....
   iparm(1) = 2       ! Use the default for DCC runs
   iparm(2) = 200     ! Limit maximum number of iterations to 500
   iparm(3) = 1       ! Warning messages and minimum output
   iparm(4) = i6      ! Define fortran unit number for output
   !iparm(21)= 0       ! Use scalar algorithm for matrix factorization
   rparm(1) = 1.E-6   ! Try default stopping test value
   iparm(12) = 1;     ! Storage mode use (1 = Primary format)

   !.....Solve for zeta.....
   WRITE (UNIT=i6,FMT='("**Enter nspcg   n = ", I6)') n
   CALL nspcg (mic1,cg,ndimA,mdimA,ndimA,maxnz1A,coeffA,jcoefA,jp,ip,zeta,   &
             & ubar,rhs,wksp,iwksp,nw1A,inw1A,iparm,rparm,ier)
   WRITE (UNIT=i6,FMT='("**Exit nspcg    n = ", I6)') n
   IF(ier /= 0) WRITE(UNIT=i6,FMT='("**ERROR", I5, " from nspcg")') ier

   !.....STOP program execution if a fatal error is encountered in nspcg.....
   IF (ier < 0 ) THEN
      PRINT *, " "
      PRINT '(" Fatal error in matrix solution on time step = ", I7)', n
      PRINT '(" **ERROR", I5, " from nspcg")', ier
      PRINT '(" Time = ", F10.4, " hours")', thrs
      PRINT *, " "
      PRINT *, " "
      PRINT *, " ****STOPPING si3d due to fatal error in matrix solution"
      WRITE (UNIT=i6,FMT='(" ****STOPPING si3d due to fatal matrix error")' )
      WRITE (UNIT=i6,FMT='(" Time = ", F10.4, " hours")') thrs
      STOP
   END IF

!   print *,sum(zeta(:))
  ! if( n .EQ. 300) THEN
  ! do m=1,lm
  !
  !    print *,zeta(m)
  !
  ! end do
  ! end if

   END IF
   !$omp barrier
   !.....Load matrix solution into  zeta  array.....
   IF(omp_get_thread_num ( )+1 == num_threads)THEN
   aux_indice=lhf(omp_get_thread_num ( )+1)
   ELSE
   aux_indice=lhfE(omp_get_thread_num ( )+1)
   END IF
   DO liter = lhi(omp_get_thread_num ( )+1), aux_indice
     m = id_column(liter)
!     i = l2i(m);
!     j = l2j(m);
     s(m) = zeta(m)
   END DO

   !.....Load matrix solution into  zeta  array.....
   !DO m = 1, ndimA
   !   i = (m-1)/ibdwdA + ifirst
   !   j = m + (jfirst-1) - (i-ifirst)*ibdwdA
   !   s(i,j) = zeta(m)
   !END DO

   ! ... Calculate surface layer thickness at time n+1

   CALL slayer_h

   !.....Save the workspace parameters, slightly over-dimensioning it
   !     Otherwise one gets some errors. - FJR: I really do not get this
   IF(omp_get_thread_num ( )==0) THEN
   nwA = FLOOR(nw1A*1.5); inwA = FLOOR(inw1A*1.5); maxnzA = maxnz1A
   !.....Deallocate arrays.....
   DEALLOCATE (wksp, iwksp )
   END IF

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_solver = t_solver + (etime - btime)

END SUBROUTINE SolverSparse



!***********************************************************************
SUBROUTINE slayer_h
!***********************************************************************
!
!  Purpose: To recompute the new values for the surface layer thicknesses
!           (h, hu, hv) after the zeta array is redefined. Note that
!           cells may become dry at time n+1, but no new cells will appear
!           at this time (wetting). The indexes for the surface layers
!           are not updated at this time, since they are used in subr. vel
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: i, j, k, l,liter

   REAL :: haux


   ! ... Initialize layer thickness for time n+1
   h(:,lm1)=hp(:,lm1)
   hu(:,lm1)=hup(:,lm1)
   hv(:,lm1)=hvp(:,lm1)

   ! ... Redo calculations for surface cells
   !     1.- If drying occurs redo calcs. at cell k1s+1
   !     2.- If cell k1s becomes thicker than its nominal size
   !         just ignore. Wetting is not done at this time
   !     Arrays storing surface cells are not modified at this
   !     time since they are used in subr. vel

   DO liter = lhiW(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
     If(liter == 0) CYCLE

     l = id_column(liter)
     h(1:k1z(l)-1,l)=hp(1:k1z(l)-1,l)
     hu(1:k1u(l)-1,l)=hup(1:k1u(l)-1,l)
     hv(1:k1v(l)-1,l)=hvp(1:k1v(l)-1,l)
     h(k1z(l)+2:km1,l)=hp(k1z(l)+2:km1,l)
     hu(k1u(l)+1:km1,l)=hup(k1u(l)+1:km1,l)
     hv(k1v(l)+2:km1,l)=hvp(k1v(l)+2:km1,l)
      ! ... Map 3D-(i,j) from 2D-l indexes
!      i = l2i(l); j = l2j(l);

      ! ... At zeta-points
      k = k1z(l)
      haux = AMIN1(zlevel(k+1),hhs(l)) + s(l)
      IF(haux <= HMIN) THEN
        h(k,l) = ZERO
        k = k + 1
        haux = AMIN1(zlevel(k+1),hhs(l)) + s(l)
        IF(haux <= HMIN) THEN
          h(k,l) = ZERO
        ELSE
          h(k,l) = haux
        END IF
      ELSE
        h(k,l) = haux
      END IF


      ! ... At u-points
     IF (mask(c2l(lEC(l)))) THEN
        k = k1u(l)
        haux=AMIN1(zlevel(k+1),hhu(l)) +             &
        &        MAX(s(l),s(lEC(l)))
        IF(haux <= HMIN) THEN
          hu(k,l) = ZERO
          k = k + 1; k = k1u(l)
          haux =AMIN1(zlevel(k+1),hhu(l)) +         &
          &        MAX(s(l),s(lEC(l)))
          IF(haux <= HMIN) THEN
            hu(k,l) = ZERO
          ELSE
            hu(k,l) = haux
          END IF
        ELSE
          hu(k,l) = haux
        END IF
      ELSE
         hu(k1u(l):k1u(l)+1,l) = hup(k1u(l):k1u(l)+1,l)
      END IF

      ! ... At v-points
     IF (mask(c2l(lNC(l)))) THEN
        k = k1v(l)
        haux=AMIN1(zlevel(k+1),hhv(l)) +             &
        &        MAX(s(l),s(lNC(l)))
        IF(haux <= HMIN) THEN
          hv(k,l) = ZERO
          k = k + 1;
          haux =AMIN1(zlevel(k+1),hhv(l)) +         &
          &        MAX(s(l),s(lNC(l)))
          IF(haux <= HMIN) THEN
            hv(k,l) = ZERO
          ELSE
            hv(k,l) = haux
          END IF
        ELSE
          hv(k,l) = haux
        END IF
      ELSE
         hv(k1v(l):k1v(l)+1,l) = hvp(k1v(l):k1v(l)+1,l)
      END IF

   ENDDO

END SUBROUTINE slayer_h




!***********************************************************************
SUBROUTINE layer_h
!***********************************************************************
!
!  Purpose: To recompute the new values for the surface layer thicknesses
!           (h, hu, hv) after the zeta array is redefined. For a linear
!           problem (ilin=1), the surface layer thicknesses are left as
!           constants. (Note: hu and hv on closed boundaries should be
!           considered undefined.)
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: i, j, k, l, kms, kmx, kmy, liter

   DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
      l = id_column(liter)

      ! ... Map 3D-(i,j) from 2D-l indexes
!      i = l2i(l); j = l2j(l);

      ! ... At zeta-points
      kms = kmz(l)
      DO k = k1, kms
        h (k,l)=AMIN1(zlevel(k+1),hhs(l)) -            &
        &       AMAX1(zlevel(  k),-s (l))
        IF(h (k,l) <= HMIN) h (k,l) = ZERO;
      ENDDO

      ! ... At u-points
      IF (mask(c2l(lEC(l)))) THEN
        kmx = MIN(kmz(l),kmz(lEC(l)))
        DO k = k1, kmx
          hu (k,l)=AMIN1(zlevel(k+1),hhu(l)) -            &
          &        AMAX1(zlevel(  k),-MAX(s(l),s(lEC(l))))
          IF(hu(k,l) <= HMIN) hu(k,l) = ZERO;
        ENDDO
      ENDIF

      ! ... At v-points
      IF (mask(c2l(lNC(l)))) THEN
        kmy = MIN(kmz(l),kmz(lNC(l)))
        DO k = k1, kmy
          hv (k,l)=AMIN1(zlevel(k+1),hhv(l)) -            &
          &        AMAX1(zlevel(  k),-MAX(s(l),s(lNC(l))))
          IF(hv (k,l) <= HMIN) hv (k,l) = ZERO;
        ENDDO
      ENDIF

   ENDDO

END SUBROUTINE layer_h





!***********************************************************************
SUBROUTINE layer_hp2
!***********************************************************************
!
!  Purpose: To recompute the old values for the surface layer thicknesses
!           (hp, hup, hvp) after the zeta array is smoothed. For a linear
!           problem (ilin=1), the surface layer thicknesses are left as
!           constants. (Note: hup and hvp on closed boundaries should be
!           considered undefined.)
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: i, j, k, l, kmx, kmy, kms, nwlup, nwlvp, nwlsp,liter

   DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

      l = id_column(liter)

      ! ... Map 3D-(i,j) from 2D-l indexes
      i = l2i(l); j = l2j(l);

      ! ... At zeta-points
      kms = kmz(l)
      DO k = k1, kms
          hp (k,l)=AMIN1(zlevel(k+1),hhs(l)) -            &
          &        AMAX1(zlevel(  k),-sp(l))
          IF(hp(k,l) <= HMIN) hp(k,l) = ZERO;
      ENDDO

      ! ... At u-points
       IF (mask2d(i+1,j)) THEN
        kmx = MIN(kmz(l),kmz(lEC(l)))
        DO k = k1, kmx
          hup(k,l)=AMIN1(zlevel(k+1),hhu(l)) -            &
          &        AMAX1(zlevel(  k),-MAX(sp(l),sp(lEC(l))))
          IF(hup(k,l) <= HMIN) hup(k,l) = ZERO;
        ENDDO
      ENDIF

      ! ... At v-points
      IF (mask2d(i,j+1)) THEN
        kmy = MIN(kmz(l),kmz(lNC(l)))
        DO k = k1, kmy
          hvp(k,l)=AMIN1(zlevel(k+1),hhv(l)) -            &
          &        AMAX1(zlevel(  k),-MAX(sp(l),sp(lNC(l))))
          IF(hvp(k,l) <= HMIN) hvp(k,l) = ZERO;
        ENDDO
      ENDIF

   ENDDO

END SUBROUTINE layer_hp2

!***********************************************************************
SUBROUTINE layer_hp3
!***********************************************************************
!
!  Purpose: To recompute the old values for the surface layer thicknesses
!           (hp, hup, hvp) after the zeta array is smoothed. For a linear
!           problem (ilin=1), the surface layer thicknesses are left as
!           constants. (Note: hup and hvp on closed boundaries should be
!           considered undefined.)
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: i, j, k, l, kmx, kmy, kms, nwlup, nwlvp, nwlsp,liter
   REAL :: aux

   DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

      l = id_column(liter)

      ! ... Map 3D-(i,j) from 2D-l indexes
!      i = l2i(l); j = l2j(l);

      ! ... At zeta-points
      kms = kmz(l)
      DO k = k1, kms
          aux=AMIN1(zlevel(k+1),hhs(l)) -            &
          &        AMAX1(zlevel(  k),-sp(l))
          IF(aux <= HMIN) THEN
          hp(k,l) = ZERO;
          ELSE
          hp (k,l)=aux
          END IF
      ENDDO

      ! ... At u-points
      IF (mask(c2l(lEC(l)))) THEN
        kmx = MIN(kmz(l),kmz(lEC(l)))
        DO k = k1, kmx
          aux=AMIN1(zlevel(k+1),hhu(l)) -            &
          &        AMAX1(zlevel(  k),-MAX(sp(l),sp(lEC(l))))
          IF(aux <= HMIN) THEN
          hup(k,l) = ZERO;
          ELSE
          hup(k,l)=aux
          END IF
        ENDDO
      ENDIF

      ! ... At v-points
      IF (mask(c2l(lNC(l)))) THEN
        kmy = MIN(kmz(l),kmz(lNC(l)))
        DO k = k1, kmy
          aux=AMIN1(zlevel(k+1),hhv(l)) -            &
          &        AMAX1(zlevel(  k),-MAX(sp(l),sp(lNC(l))))
          IF(aux <= HMIN) THEN
          hvp(k,l) = ZERO;
          ELSE
          hvp(k,l)=aux
          END IF
        ENDDO
      ENDIF

   ENDDO

END SUBROUTINE layer_hp3


!***********************************************************************
SUBROUTINE TopLayerIndexp2
!***********************************************************************
!
!  Purpose: To determine the top layer index given that values of hp,
!           hup and hvp have been calculated previously
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: i, j, k, l, kmx, kmy, kms, nwlup, nwlvp, nwlsp,liter


   DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

       l = id_column(liter)
       k1z(l)=km1;
       k1u(l)=km1;
       k1v(l)=km1;

      ! ... Map 3D-(i,j) from 2D-l indexes
!      i = l2i(l); j = l2j(l);

      ! ... At zeta-points
      kms = kmz(l)
      DO k = k1, kms
        IF(hp (k,l) > ZERO) THEN
          k1z(l) = k
          EXIT
        ENDIF
      ENDDO

      ! ... At u-points
      IF (mask(c2l(lEC(l)))) THEN
        kmx = MIN(kmz(l),kmz(lEC(l)))
        DO k = k1, kmx
          IF(hup(k,l) > ZERO) THEN
            k1u(l) = k
            EXIT
          ENDIF
        ENDDO
      ENDIF

      ! ... At v-points
      IF (mask(c2l(lNC(l)))) THEN
        kmy = MIN(kmz(l),kmz(lNC(l)))
        DO k = k1, kmy
          IF(hvp(k,l) > ZERO) THEN
            k1v(l) = k
            EXIT
          ENDIF
        ENDDO
      ENDIF

   ENDDO

END SUBROUTINE TopLayerIndexp2

!***********************************************************************
SUBROUTINE vel(Bstart,Bend,Bagx,Barx,Bagy,Bary)
!***********************************************************************
!
!  Purpose: To solve the momentum equations explicitly for velocity.
!
!-----------------------------------------------------------------------

   INTEGER, INTENT(IN) :: Bstart,Bend
   REAL, DIMENSION (1:km1,Bstart:Bend+1), INTENT(INOUT) :: Bagx,Barx,Bagy,Bary



   !.....Local variables.....
   REAL :: gthx1, gthy1, rho4cxx, rho4cyy, vvtemp, uutemp, cxx, cyy
   INTEGER :: i, j, k, l, istat, kmx, kmy, kms, k0x, k0y, k0s, k1s, k1ss,liter

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)

   !.....Constants.....
   gthx1 = gdtdx*tz; gthy1 = gdtdy*tz

   !                -----X-momentum equation-----

   ! ... Define constants: Ignore density variations in z
   rho4cxx = 1000.
   rho4cyy = 1000.

   ! ... Loop over cells
   DO liter = lhiWCE(omp_get_thread_num ( )+1), lhfCE(omp_get_thread_num ( )+1)
      if(liter==0) CYCLE
      l = id_columnCE(liter)

      ! .... Map l into 2D-xy space
!      i = l2i(l);
!      j = l2j(l);

      ! ... Skip if East Column is dry
      !IF(.NOT.mask(c2l(lEC(l)))) CYCLE

      !.....Top & Bottom wett u-points
      kmx = MIN(kmz(lEC(l)),kmz(l))
      k0x = k1u(l)

      !.....Solve for the water surface slope portion of the
      !     x-mom eq and save the result in the cxx array.....
      cxx = gthx1 * rho4cxx * (s(lEC(l))-s(l))

      !.....Solve the x-momentum equation for uh.....
      DO k = k0x, kmx
        uh(k,l) =  Bagx(k,l)    - cxx*Barx(k,l)
      ENDDO

      ! ... Redo near surface flux calcs. if Drying occurs
      !     hu is calculated after the solution of zeta
      !     in subr. solver. The top most layer during
      !     time n (n+1/2) remains the same through the
      !     calculations to predict n+1 from n or n+1/2
      k = k0x;
      IF (hu(k  ,l) <= ZERO) THEN
        ! ... Update fluxes
        uh(k+1,l) = uh(k,l)+uh(k+1,l)
        uh(k  ,l) = 0.0
        ! ... Update surface array
        k1u(l) = k0x+1
      ENDIF

   END DO

   !                -----Y-momentum equation-----


   ! ... Loop over cells
   DO liter = lhiCN(omp_get_thread_num ( )+1), lhfCN(omp_get_thread_num ( )+1)

      l = id_columnCN(liter)

      ! .... Map l into 2D-xy space
!      i = l2i(l);
!      j = l2j(l);

      ! ... Skip if North Column is dry
      !IF(.NOT.mask(c2l(lNC(l)))) CYCLE

      !.....Top & Bottom wett v-points .....
      kmy = MIN(kmz(lNC(l)),kmz(l))
      k0y = k1v(l)

      !.....Solve for the water surface slope portion of the
      !     y-mom eq and save the result in the cyy array.....
      cyy = gthy1 * rho4cyy * (s(lNC(l))-s(l))

      !.....Solve the y-momentum equation for vh.....
      DO k = k0y, kmy
        vh(k,l) =  Bagy(k,l)    - cyy*Bary(k,l)
      ENDDO

      ! ... Redo near surface flux calcs. if Drying occurs
      !     hu is calculated after the solution of zeta
      !     in subr. solver. The top most layer during
      !     time n (n+1/2) remains the same through the
      !     calculations to predict n+1 from n or n+1/2
      k = k0y;
      IF (hv(k  ,l) <= ZERO) THEN
        ! ... Update fluxes
        vh(k+1,l) = vh(k,l)+vh(k+1,l)
        vh(k  ,l)   = 0.0
        ! ... Update surface array
        k1v(l) = k0y+1
      ENDIF

   END DO

   !.... Calculate vertical velocities from uh & uhpp values
   !     to be used in the solution of scalar transport eq.
   CALL continuity(1)

   ! ... Update surface array at s-points if drying occurs
   DO liter = lhiW(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
      IF(liter==0) CYCLE
      l = id_column(liter)
      ! .... Map l into 2D-xy space
!      i = l2i(l); j = l2j(l);
      ! ... Top wett s-points
      k0s = k1z(l);
      ! ... Recalculate surface cell if water column
      !     is wett.
      !IF (k0s<km1 .AND. h(k0s,l)<=ZERO) k1z(i,j)=k0s+1
      IF (h (k0s,l)<=ZERO) k1z(l)=MIN(k0s+1,km1)
   ENDDO

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_vel = t_vel + (etime - btime)

END SUBROUTINE vel

!***********************************************************************
SUBROUTINE vel2
!***********************************************************************
!
!  Purpose: To recompute velocities at new time layer if wetting occurs
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!-----------------------------------------------------------------------

   !.....Local variables.....
   REAL :: vvtemp, uutemp
   INTEGER :: i, j, k, l, liter

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)

   ! ... Recalculate layer thicknesses including wetting & drying
   CALL layer_h

   ! ... Loop over cells
   DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

      l = id_column(liter)

     ! .... Map l into 2D-xy space

!     i = l2i(l); j = l2j(l);

     ! ... Skip if East Column is dry
     IF(mask(c2l(lEC(l)))) THEN
       k = k1u(l);
       IF (hu(k-1,l) >  ZERO) THEN
         uutemp    = uh(k,l)/(hu(k,l)+hu(k-1,l))
         uh(k  ,l) = uutemp * hu(k  ,l);
         uh(k-1,l) = uutemp * hu(k-1,l);
       ENDIF
     ENDIF

     ! ... Skip if North Column is dry
     IF(mask(c2l(lNC(l)))) THEN
       k = k1v(l);
       IF (hv(k-1,l) >  ZERO) THEN
         vvtemp    = vh(k,l)/(hv(k,l)+ hv(k-1,l))
         vh(k  ,l) = vvtemp * hv(k  ,l);
         vh(k-1,l) = vvtemp * hv(k-1,l);
       ENDIF
     ENDIF

   END DO

   !.....No need to recalculate vertical velocities since
   !     they are calculated in either save or setrap routines for
   !     next iteration or next time setp

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_vel = t_vel + (etime - btime)

END SUBROUTINE vel2


!***********************************************************************
SUBROUTINE continuity (ist)
!***********************************************************************
!
!  Purpose: To compute wp from horizontal components of the velocity
!           field
!
!-----------------------------------------------------------------------

   ! ... Arguments
   INTEGER, INTENT(IN):: ist


   ! ... Local variables
   INTEGER :: i, j, k, l, k1s, kms, inn,liter,innH
   REAL :: uhp_aux

   SELECT CASE (ist)

   CASE (1) ! Compute wp from uh, uhpp, vh and vhpp

     DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

       l = id_column(liter)

       ! .... Map l into 2D-xy space
!       i = l2i(l); j = l2j(l);

       ! ... Bottom wett s-points
       kms = kmz(l);
       k1s = k1z(l);

       ! ... Cycle if water column is dry
       IF( k1s > kms ) CYCLE

       ! .... Loop over cells in water colum to estimate vertical
       !      velocities which are consistent with the formulation of
       !      continuity (mass conservation). These velocities are then
       !      used for scalar transport calculations
       wp(:,l) = 0.0
       DO k = kms,k1s,-1
         wp(k,l) = wp  (k+1,l)                         &
             &   -(uh  (k  ,l)-uh  (k  ,lWC(l))+       &
                   uhpp(k  ,l)-uhpp(k  ,lWC(l)))/twodx &
             &   -(vh  (k  ,l)-vh  (k  ,lSC(l))+       &
                   vhpp(k  ,l)-vhpp(k  ,lSC(l)))/twody
       ENDDO

       ! ... Correct wp estimates for surface cell, due to
       !     advective flux from neighbouring cells above it
       DO k = k1,k1s-1
         wp(k1s,l) = wp  (k1s,l)                         &
               &   -(uh  (k  ,l)-uh  (k  ,lWC(l))+       &
                     uhpp(k  ,l)-uhpp(k  ,lWC(l)))/twodx &
               &   -(vh  (k  ,l)-vh  (k  ,lSC(l))+       &
                     vhpp(k  ,l)-vhpp(k  ,lSC(l)))/twody
       ENDDO

     ENDDO

     ! .... Modify w estimates to incorporate sources/sinks (PSS)
     IF ( iopssH(omp_get_thread_num ( )+1) > 0 ) THEN
       DO innH = 1, iopssH(omp_get_thread_num ( )+1)
         inn = ioph2iop(innH,omp_get_thread_num ( )+1)
         i = ipss(inn)
         j = jpss(inn)
         l = ij2l(i,j);
         kms = kmz(l);
         k1s = k1z(l);
         DO k = kms,k1s,-1
           wp(k,l) = wp(k,l) + SUM(Qpss(k:kms,inn))/(dx*dy)
         ENDDO
       ENDDO
     ENDIF

   ! ... Calculate wp in terms of uhp and vhp values - for computations
   !     of velocities at next time step
   CASE (2)

     DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

       l = id_column(liter)

       ! .... Map l into 2D-xy space
!       i = l2i(l); j = l2j(l);

       ! ... Bottom wett s-points
       kms = kmz(l);
       k1s = k1z(l);

       ! ... Cycle if water column is dry
       IF( k1s > kms ) CYCLE

       ! .... Loop over cells in water colum to estimate vertical
       !      velocities
       wp(:,l) = 0.0
       DO k = kms,k1s,-1
         wp(k,l) = wp(k+1,l)-(uhp(k,l)-uhp(k,lWC(l)))/dx    &
                  &         -(vhp(k,l)-vhp(k,lSC(l)))/dy
       END DO

       ! ... Correct wp estimates for surface cell, due to
       !     advective flux from neighbouring cells above it
       DO k = k1,k1s-1
         wp(k1s,l) = wp(k1s,l)-(uhp(k,l)-uhp(k,lWC(l)))/dx    &
                  &           -(vhp(k,l)-vhp(k,lSC(l)))/dy
       END DO

     END DO

     ! .... Modify w estimates to incorporate sources/sinks (PSS)
     IF ( iopssH(omp_get_thread_num ( )+1) > 0 ) THEN
       DO innH = 1, iopssH(omp_get_thread_num ( )+1)
         inn = ioph2iop(innH,omp_get_thread_num ( )+1)
         i = ipss(inn)
         j = jpss(inn)
         l = ij2l(i,j);
         kms = kmz(l);
         k1s = k1z(l);
         DO k = kms,k1s,-1
           wp(k,l) = wp(k,l) + SUM(Qpss(k:kms,inn))/(dx*dy)
         ENDDO
       ENDDO
     ENDIF

     CASE (3)

     DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

       l = id_column(liter)

       ! .... Map l into 2D-xy space
!       i = l2i(l); j = l2j(l);

       ! ... Bottom wett s-points
       kms = kmz(l);
       k1s = k1z(l);

       ! ... Cycle if water column is dry
       IF( k1s > kms ) CYCLE

       ! .... Loop over cells in water colum to estimate vertical
       !      velocities
       wp(:,l) = 0.0
       DO k = kms,k1s,-1
         if(hup(k,l)>ZERO)THEN
         uhp_aux=uh(k,lWC(l))
         else
         uhp_aux=0.0
         end if
         wp(k,l) = wp(k+1,l)-(uhp(k,l)-uhp_aux)/dx    &
                  &         -(vhp(k,l)-vhp(k,lSC(l)))/dy
       END DO

       ! ... Correct wp estimates for surface cell, due to
       !     advective flux from neighbouring cells above it
       DO k = k1,k1s-1
         if(hup(k,l)>ZERO)THEN
         uhp_aux=uh(k,lWC(l))
         else
         uhp_aux=0.0
         end if
         wp(k1s,l) = wp(k1s,l)-(uhp(k,l)-uhp_aux)/dx    &
                  &           -(vhp(k,l)-vhp(k,lSC(l)))/dy
       END DO

     END DO

   ! .... Modify w estimates to incorporate sources/sinks (PSS)
     IF ( iopssH(omp_get_thread_num ( )+1) > 0 ) THEN
       DO innH = 1, iopssH(omp_get_thread_num ( )+1)
         inn = ioph2iop(innH,omp_get_thread_num ( )+1)
         i = ipss(inn)
         j = jpss(inn)
         l = ij2l(i,j);
         kms = kmz(l);
         k1s = k1z(l);
         DO k = kms,k1s,-1
           wp(k,l) = wp(k,l) + SUM(Qpss(k:kms,inn))/(dx*dy)
         ENDDO
       ENDDO
     ENDIF

   END SELECT

   ! ... Modify velocity estimates near the boundaries to
   !     account for open boundaries.
   IF (nopen > 0) CALL MODvel4openBC


 END SUBROUTINE continuity

!***********************************************************************
SUBROUTINE exsal(Bstart,Bend,lSCH,lNCH,lECH,lWCH,Bhaxpp,Bhaypp,Bth3,Bth4,Bth2,Bex,thrs)
!***********************************************************************
!
!  Purpose: To evaluate the explicit terms (advection) in the scalar
!           transport equation using flux limiter methods. The sum of these
!           terms are saved in the array  ex(k,l)  which is the
!           primary output from this subroutine.
!
!-----------------------------------------------------------------------

   INTEGER, INTENT(IN) :: Bstart,Bend
   REAL, INTENT(IN) :: thrs
   REAL, DIMENSION (1:km1,Bstart:Bend+1), INTENT(INOUT) :: Bhaxpp, Bhaypp,Bth3,Bth4,Bth2,Bex
   INTEGER, DIMENSION (Bstart:Bend+1), INTENT(IN) :: lSCH,lNCH,lWCH,lECH

   ! ... Local variables
   INTEGER :: i, j, k, l, k1s, kms, gamma1, istat,liter
   REAL    :: vel, ratio, C_f, delz, twodt1, hd
   REAL, DIMENSION (4          ) :: ss

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)

   ! ... Constants used in solution
   twodt1 = twodt*tz

   ! ... Calculate hdxpp & hdypp arrays for diffusion terms
   Bhaxpp(:,Bend+1) = 0.0;
   Bhaypp(:,Bend+1) = 0.0;
   DO liter = lhiW(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
      If(liter == 0) CYCLE

      l = id_column(liter)

      ! ... 3D-(i,j) indexes for l
!      i = l2i(l); j = l2j(l);

      ! ... Retrieve top & bottom wet sal-pts .................
      kms = kmz(l)
      k1s = k1z(l)

      Bhaxpp(1:k1-1,l) = 0.0
      Bhaxpp(kms+1:km1,l) = 0.0
      Bhaypp(1:k1-1,l) = 0.0
      Bhaypp(kms+1:km1,l) = 0.0

      ! ... Calculate hdxpp & hdypp array at u-&v- pts ........
      !     Interfaces connecting wett & dry cells will not
      !     have diffussive transport in present formulation
      DO k = k1, kms
        Bhaxpp(k,l) = Ax0*hupp(k,l)
        Bhaypp(k,l) = Ay0*hvpp(k,l)
      ENDDO

   END DO

   Bth3(:,Bend+1) = 0
!   Bth2 = 0
   ! ... Initialize ex & flux arrays to zeros
   Bex(:,Bend+1) = 0.0;  Bth4(:,Bend+1) = 0.0;! fluxXsal(:,lm1)= 0.0;
   Bth2(:,Bend+1) = 0.0;
!    Bth3 = 0.0; Bth2 = 0.0; Bex = 0.0; Bth4 = 0.0;

   DO liter = lhiW(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
      If(liter == 0) CYCLE

      l = id_column(liter)

     ! ... Map l- into (i,j)-indexes .........................
!     i = l2i(l); j = l2j(l);

     ! ... Retrieve top & bottom wet sal-pts .................
     kms = kmz(l)
     k1s = k1z(l)

     DO k = k1s, kms;

       ! ... EW fluxes .......................................
       IF (hup(k,l)> ZERO) THEN

         ! ... Velocities at time n+1/2
         vel  = (uhpp(k,l) + uh(k,l))/2.

         ! ... Define stencil for scalar transport
         ss(2)  = salpp(k,    l );
         ss(3)  = salpp(k,lEC(l));
         IF (hpp(k,    lWC(l) )<=ZERO) THEN; ss(1) = ss(2);
          ELSE; ss(1)=salpp(k,    lWC(l) ); ENDIF;
         IF (hpp(k,lEC(lEC(l)))<=ZERO) THEN; ss(4) = ss(3);
          ELSE; ss(4)=salpp(k,lEC(lEC(l))); ENDIF;

         ! ... Calculate Cf for flux computation
         C_f    = 0.0;
         gamma1 = -SIGN (1., vel)
         IF (ss(3) - ss(2) /= 0 ) THEN
           ratio =(ss(3+gamma1)-ss(2+gamma1))/(ss(3)-ss(2))
           ! MC flux limiter (VanLeer, 1977)
           !C_f   = MAX(0., MIN( 2*ratio, (1+ratio)/2., 2. ))
           ! ... Roe's Superbee Limiter
           C_f = MAX(0., MIN(1.,2.*ratio),MIN(2.,ratio))
         ENDIF

         ! ... Calculate fluxes
         Bth2(k,l) = vel/2.*(ss(3)+ss(2))- &
         & ((1.-C_f)*ABS(vel)+vel**2.*twodt1/dx*C_f)*(ss(3)-ss(2))/2.
       ELSE
         Bth2(k,l) = 0.0
       ENDIF
       ENDDO;
   ENDDO;

   DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
     l = id_column(liter)

     ! ... Map l- into (i,j)-indexes .........................
!     i = l2i(l); j = l2j(l);

     ! ... Retrieve top & bottom wet sal-pts .................
     kms = kmz(l)
     k1s = k1z(l)

     DO k = k1s, kms;

       ! ... NS fluxes .......................................
       IF (hvp(k,l)> ZERO) THEN

         ! ... Velocities at time n+1/2
         vel  = (vhpp(k,l) + vh(k,l))/2.

         ! ... Define stencil for scalar transport
         ss(2)  = salpp(k,        l  );
         ss(3)  = salpp(k,    lNC(l) );
         IF (hpp(k,    lSC(l) )<= ZERO) THEN; ss(1) = ss(2);
           ELSE; ss(1)  = salpp(k,    lSC(l) ); ENDIF;
         IF (hpp(k,lNC(lNC(l)))<= ZERO) THEN; ss(4) = ss(3);
           ELSE; ss(4)  = salpp(k,lNC(lNC(l))); ENDIF;

         ! ... Calculate Cf for flux computation
         C_f    = 0.0
         gamma1 = -SIGN (1., vel)
         IF (ss(3) - ss(2) /= 0 ) THEN
           ratio =(ss(3+gamma1)-ss(2+gamma1))/(ss(3)-ss(2))
           ! MC flux limiter (VanLeer, 1977)
           !C_f   = MAX(0., MIN( 2*ratio, (1+ratio)/2., 2. ))
           ! ... Roe's Superbee Limiter
           C_f = MAX(0., MIN(1.,2.*ratio),MIN(2.,ratio))
         ENDIF


         ! ... Calculate fluxes
         Bth4(k,l) = vel/2.*(ss(3)+ss(2))- &
         & ((1.-C_f)*ABS(vel)+vel**2.*twodt1/dx*C_f)*(ss(3)-ss(2))/2.
       ELSE
         Bth4(k,l) = 0.0
       ENDIF

       ! ... UD fluxes .......................................
       IF (hp(k-1,l) > ZERO) THEN

         ! ... Velocities at time n + 1
         vel  = wp(k,l); IF (k == k1s) vel = 0.0;

         ! ... Define stencil for scalar transport
         ss(2) = salpp(k  ,l);
         ss(3) = salpp(k-1,l);
         IF (hpp(k-2,l)<=ZERO) THEN ; ss(4)=ss(3);
            ELSE; ss(4)=salpp(k-2,l); ENDIF;
         IF (hpp(k+1,l)<=ZERO) THEN ; ss(1)=ss(2);
            ELSE; ss(1)=salpp(k+1,l); ENDIF;

         ! ... Define C_f for flux computations
         C_f   = 1.0 ! Default method is Lax-Wendroff
         gamma1 = -SIGN (1., vel)
         IF (ss(3) - ss(2) /= 0 ) THEN
           ratio =(ss(3+gamma1)-ss(2+gamma1))/(ss(3)-ss(2))
           ! MC flux limiter (VanLeer, 1977)
           !C_f   = MAX(0., MIN( 2*ratio, (1+ratio)/2., 2. ))
           ! ... Roe's Superbee Limiter
           C_f = MAX(0., MIN(1.,2.*ratio),MIN(2.,ratio))
         ENDIF

         ! ... Calculate fluxes
         delz = (hp(k,l) + hp(k-1,l))/2.
         Bth3(k,l) = vel/2.*(ss(3)+ss(2))- &
         & ((1.-C_f)*ABS(vel)+vel**2.*twodt1/delz*C_f)*(ss(3)-ss(2))/2.
       ELSE
       	 Bth3(k,l) = 0.0
       ENDIF

     ENDDO;
   ENDDO;

   ! ... Update ex array with x-flux divergence
   DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
     l = id_column(liter)

     ! ... Map l- into (i,j)-indexes .........................
!     i = l2i(l); j = l2j(l);

     ! ... Retrieve top & bottom wet sal-pts .................
     kms = kmz(l)
     k1s = k1z(l)

     Bex(1:k1s-1,l) = 0.0
     Bex(kms+1:km1,l) = 0.0

     DO k = k1s, kms;

        !.....Horizontal diffusion.....
        hd= (Bhaxpp(k,    l )*(salpp(k,lEC(l)) - salpp(k,    l ))       &
          & -Bhaxpp(k,lWCH(l))*(salpp(k,    l ) - salpp(k,lWC(l))))/dxdx &
          &+(Bhaypp(k,    l )*(salpp(k,lNC(l)) - salpp(k,    l ))       &
          & -Bhaypp(k,lSCH(l))*(salpp(k,    l ) - salpp(k,lSC(l))))/dydy

        !.....Sum all terms
        Bex(k,l) =   hpp(k,l)*salpp(k,l)/twodt1        &
                - (Bth2(k,l) - Bth2(k,lWCH(l))) / dx &
                - (Bth4(k,l) - Bth4(k,lSCH(l))) / dy &
                - (Bth3(k,l) - Bth3(k+1,l   )) !+ hd * ihd
        IF (ihd>0) Bex(k,l) = Bex(k,l) + hd   ! Changed 12/2010 SWA

     ENDDO;

   ENDDO

   CALL MODexsal4openbc(Bstart,Bend,Bex,thrs)

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_exsal = t_exsal + (etime - btime)

END SUBROUTINE exsal

!***********************************************************************
SUBROUTINE imsal(Bstart,Bend,Bex,heatSourceB)
!***********************************************************************
!
!  Purpose: To solve for active scalar concentration.
!
!  29-may-2009	(F.J.Rueda)	 Include tpload (temp. load associated to rwps)
!
!-----------------------------------------------------------------------

   INTEGER, INTENT(IN) :: Bstart,Bend
   REAL, DIMENSION (1:km1,Bstart:Bend+1), INTENT(INOUT) :: Bex,heatSourceB

   !.....Local variables.....
   REAL :: twodt1, Tsource, Qsource
   INTEGER :: i, j, k, l, k1s, kms, kt, nwlayers, inn, kk, noc,liter,lol,innH
   REAL, DIMENSION (1:km1) :: hn
   REAL, DIMENSION (3,1:km1) :: aa
   REAL, DIMENSION (1:km) :: dsal
   REAL, DIMENSION (1:ndz) :: sal1

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)

   ! ... Constants used in solution
   twodt1 = twodt*tz

   !.....Loop over interior sal-pts to solve for
   !     matrix from the active scalar equation.....
   DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
     l = id_column(liter)

      ! ... 3D-(i,j) indexes for l - FJR - uncomment
      i = l2i(l); j = l2j(l);

      !.....Compute top & bottom layer numbers & No. of layers ....
      kms = kmz(l);
      k1s = k1z(l);
      nwlayers = kms-k1s+1


!      print *,"l:",l,"kmz:",kmz(l),"k1z:",k1z(l),"Hilo:",omp_get_thread_num ( )
      ! ... Define layer thikness at time n - The corrections for
      !     surface and recently submerged cells are needed to
      !     keep mass conservation -  The test used to check
      !     mass conservation is that of a surface seiche with
      !     an equilibrium water surface level at the level of
      !     where grids cells change from level k to k+1
      hn(k1s+1:kms) = h(k1s+1:kms,l)
      hn(k1s      ) = twodt1*wp(k1s,l)+hpp(k1s,l)
      IF (hpp(k1s,l)<= ZERO) THEN
        hn(k1s+1) = hpp(k1s+1,l)
      ENDIF

      SELECT CASE (nwlayers)
      !.....Calculate active scalar for case of a single layer.....
      CASE (1)

        ! ... Use h(k1s,l) instead of hn(k1s) -
        aa( 2,k1s) = hn(k1s)/twodt1
        dsal(   k1s) = Bex(k1s,l) + HeatSource(k1s,l)
        sal(k1s,l) = dsal(k1s  )/aa(2,k1s)

        ! ... For one-layer columns that become dry - The
        !     value of the threshold 1.E-5 is completely
        !     arbitrary but small - it is used to avoid the
        !     occurrence of errors in scalar conc. etimates
        !     arising from errors in dividing ds by aa -
        !     Note that the error allowed in estimating zeta
        !     is of O(10-6) - see SUB. SOLVER
        IF (h(k1s,l) < 1.E-2) sal(k1s,l) = salpp(k1s,l)

      !.....Calculate active scalar for case of two or more layers.....
      CASE (2:)

         !.....Form coefficient matrix [aa]
         ! Define upper diagonal terms
         aa(3,k1s:kms-1) = -Dv(k1s+1:kms,l)/(hn(k1s:kms-1)+hn(k1s+1:kms))*2.
         aa(3,kms)       =  0.0
         ! Define lower diagonal terms
         aa(1,k1s+1:kms) = -Dv(k1s+1:kms,l)/(hn(k1s:kms-1)+hn(k1s+1:kms))*2.
         aa(1,k1s)       =  0.0
         ! Define center diagonal terms
         aa(2,k1s:kms)   =  hn(k1s:kms)/twodt1-aa(1,k1s:kms)-aa(3,k1s:kms)

         !.....form r.h.s. matrix [ds].....
         DO k = k1s, kms
            dsal(k) = Bex(k,l) + HeatSource(k,l)
         ENDDO

         !.....Solve tridiagonal system for the
         !     vertical distribution of active scalar.....

         ! ....Modify matrices to take into accout mixing action of plumes
         IF ( iopssH(omp_get_thread_num ( )+1) > 0) THEN
           DO innH = 1, iopssH(omp_get_thread_num ( )+1)
             inn = ioph2iop(innH,omp_get_thread_num ( )+1)

             IF ( j /= jpss(inn) .OR. i /=ipss(inn) ) CYCLE
             DO k = k1s, kms
               IF (ABS(Qpss(k,inn))<1.E-10) CYCLE
               Qsource  = Qpss(k,inn)/(dx*dy)
               Tsource  = Tpss(k,inn)
               dsal(k)    = dsal(k)+Qsource*Tsource
             ENDDO
           ENDDO

        ENDIF

         CALL trid1 (aa, dsal, sal1, k1s, kms, km1, nwlayers)

         !.....Define scalars at new time step....
         sal(k1s:kms  ,l) = sal1(1:nwlayers)
         sal(k1 :k1s-1,l) = sal1(1         )

      END SELECT

   !.....End loop over scalar-pts.....
   END DO

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_salin = t_salin + (etime - btime)

END SUBROUTINE imsal

!***********************************************************************
SUBROUTINE trid ( acoef, g, r, ag, ar, k1, km, km1, n )
!***********************************************************************
!
!  Purpose: Tridiagonal matrix solver for the momentum equation using
!           the double-sweep method
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!   7/1/98           P.E. Smith        Original f90 code
!  2/15/99           P.E. Smith        Reset kmax from 20 to 60 layers
!  6/15/99           P.E. Smith        Reset kmax to 200 layers
!
!-----------------------------------------------------------------------

   !.....Dimensioning parameter.....
   INTEGER, PARAMETER :: kmax = 500

   !.....Arguments.....
   INTEGER, INTENT(IN) :: k1, km, km1, n
   REAL, DIMENSION(km1), INTENT(IN)    :: g, r
   REAL, DIMENSION(km1), INTENT(INOUT) :: ag, ar
   REAL, DIMENSION(3,km1), INTENT(IN)  :: acoef

   !.....Local variables.....
!   INTEGER, AUTOMATIC :: k, kk
!   REAL, AUTOMATIC, DIMENSION(kmax) :: a, b, c, d, e, c1, d1, e1, d2, e2

   !.....Local variables.....
   INTEGER, AUTOMATIC :: k, kk
   REAL, AUTOMATIC, DIMENSION(kmax) :: a, b, c, d, e, c1, d1, e1, d2, e2

   !.....Timing.....
!   REAL, EXTERNAL :: TIMER
!   REAL :: btime, etime
!   btime = TIMER(0.0)
    n_trid = n_trid + 1

   !.....Load diagonals of coefficient matrix into
   !     1-d arrays and define r.h.s. vectors.....
   k = 0
   DO kk = k1, km
      k = k + 1
      a(k) = acoef(1,kk)
      b(k) = acoef(2,kk)
      c(k) = acoef(3,kk)
      d(k) = g(kk)
      e(k) = r(kk)
   END DO

   !.....Forward sweep--transform coefficient
   !     matrix into upper bidiagonal form.....
   c1(1) = c(1)/b(1)
   d1(1) = d(1)/b(1)
   e1(1) = e(1)/b(1)
   DO k = 2, n
      c1(k) = c(k)/(b(k) - a(k)*c1(k-1))
      d1(k) = (d(k) - a(k)*d1(k-1))/(b(k) - a(k)*c1(k-1))
      e1(k) = (e(k) - a(k)*e1(k-1))/(b(k) - a(k)*c1(k-1))
   END DO

   !.....Backward sweep--transform coefficient
   !     matrix into diagonal form.....
   d2(n) = d1(n)
   e2(n) = e1(n)
   DO k = n-1, 1, -1
      d2(k) = d1(k) - c1(k)*d2(k+1)
      e2(k) = e1(k) - c1(k)*e2(k+1)
   END DO

   !.....Load r.h.s. solution vectors into  [ag]  and  [ar]
   !     arrays for passing back to the calling program.....
   k = 0
   DO kk = k1, km
      k = k + 1
      ag(kk) = d2(k)
      ar(kk) = e2(k)
   END DO

END SUBROUTINE trid

!***********************************************************************
SUBROUTINE trid1 ( acoef, dsal, sal, k1, km, km1, n )
!***********************************************************************
!
!  Purpose: Tridiagonal matrix solver for the salinity equation using
!           the double-sweep method
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!   7/1/98           P.E. Smith        Original f90 code
!   5/1/99           P.E. Smith        Reset kmax from 20 to 60 layers
!  6/15/99           P.E. Smith        Reset kmax to 200 layers
!
!-----------------------------------------------------------------------

   !.....Dimensioning parameter.....
   INTEGER, PARAMETER :: kmax = 500

   !.....Arguments.....
   INTEGER, INTENT(IN) :: k1, km, km1, n
   REAL, DIMENSION(km), INTENT(INOUT)  :: dsal
   REAL, DIMENSION(n),  INTENT(INOUT)  :: sal
   REAL, DIMENSION(3,km1), INTENT(IN)  :: acoef

   !.....Local variables.....
   INTEGER :: k, kk
   REAL, DIMENSION(kmax) :: a, b, c, d, e, f

   ! Note: n = number of 3-d layers (also number of unknowns)

   !.....Load diagonals of coefficient matrix into
   !     1-d arrays and rename r.h.s. vector.....
   k = 0
   DO kk = k1, km
      k = k + 1
      a(k) = acoef(1,kk)
      b(k) = acoef(2,kk)
      c(k) = acoef(3,kk)
      d(k) = dsal(kk)
   END DO

   !.....Initialize for forward sweep.....
      e(1) = -c(1)/b(1)
      f(1) =  d(1)/b(1)

   !.....Forward sweep (solve for  e  and  f  vectors).....
   IF(n == 2) GO TO 1
   DO k = 2, n-1
      e(k) = -c(k)/(b(k)+a(k)*e(k-1))
      f(k) = (d(k)-a(k)*f(k-1))/(b(k)+a(k)*e(k-1))
   END DO

   !.....Compute salinity in bottom layer.....
 1 sal(n) = (d(n)-a(n)*f(n-1))/(b(n)+a(n)*e(n-1))

   !.....Backward sweep (solve for salinity vector).....
   DO k = n-1, 1, -1
      sal(k) = e(k)*sal(k+1) + f(k)
   END DO

END SUBROUTINE trid1


!***********************************************************************
SUBROUTINE smooth
!***********************************************************************
!
!  Purpose: To smooth the solution from the leapfrog step with the
!           Asselin time filter (Mon. Weather Rev., v. 100, 1972,
!           p. 487-490). Smoothing is only performed if ismooth>=1.
!           The degree of smoothing is determined from the parameter
!           beta. Beta=0.05 is recommended. Values as high as 1.0
!           can be used. The choices for ismooth are:
!                If ismooth = 0 --> no smoothing
!                If ismooth = 1 --> smooth zeta and velocity
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!-----------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: i, j, k, l, kmx, kmy, kms, k1s, k1x, k1y,liter
   REAL    :: wght, wghtpp, scC, scCpp

   !.....Smooth zeta (and then recalculate hp, hup, and hvp).....
   DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
     l = id_column(liter)
!     i = l2i(l); j = l2j(l);
     sp(l) = sp(l) + (beta2)*(s(l)-2.*sp(l)+spp(l))
   ENDDO
   !$omp barrier
   CALL layer_hp3

   !.....Smooth horizontal velocity components.....
   DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
     l = id_column(liter)
     ! ... Map l- into (i,j)-indexes
!     i = l2i(l); j = l2j(l);
     ! ... At u-points
     IF(mask(c2l(lEC(l)))) THEN
       kmx = MIN(kmz(l),kmz(lEC(l)))
       DO k = k1, kmx
         IF(hup(k,l)<=0.0) CYCLE
         uhp(k,l)= uhp(k,l)+beta2*(uh(k,l)-2.*uhp(k,l)+uhpp(k,l))
         up (k,l)= uhp(k,l)/hup(k,l)
       ENDDO
     ENDIF
     ! ... At v-points
     IF(mask(c2l(lNC(l)))) THEN
       kmy = MIN(kmz(l),kmz(lNC(l)))
       DO k = k1, kmy
         IF(hvp(k,l)<=0.0) CYCLE
         vhp(k,l)= vhp(k,l)+beta2*(vh(k,l)-2.*vhp(k,l)+vhpp(k,l))
         vp (k,l)= vhp(k,l)/hvp(k,l)
       ENDDO
     ENDIF
   ENDDO
   !$omp barrier
   !.....No need to recalculate vertical velocity components
   !     since these values should be stored in wpp either in save or
   !     in settrap, which are not used in any computations -

END SUBROUTINE smooth

!***********************************************************************
SUBROUTINE settrap
!***********************************************************************
!
!  Purpose: To setup the arrays at the  n  and  n+1/2  time levels
!           for use in the first iteration of the trapezoidal step.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-----------------------------------------------------------------------



   !.....Local variables.....
   INTEGER :: i, j, k, l, kmx, kmy, kms, k1x, k1y, k1s,liter,no,ide_t,is,ie,js,je,nn
   REAL    :: uutemp, vvtemp, wght, wghtpp, scC, scCpp

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)


   ide_t = omp_get_thread_num ( )+1
   spp(lm1)=sp(lm1)
   hpp(:,lm1)=hp(:,lm1)
   hupp(:,lm1)=hup(:,lm1)
   hvpp(:,lm1)=hvp(:,lm1)
   DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
     l = id_column(liter)
   !....Zeta array.....
   spp(l)  = sp(l)
   sp(l)   = 0.5*(s(l) + spp(l))

   ! ... Save layer thickness at time n
   hpp(:,l)  = hp(:,l);
   hupp(:,l) = hup(:,l);
   hvpp(:,l) = hvp(:,l);
   END DO
   !$omp barrier
   ! ... Define layer thickness at time n+1/2 &
   !     recompute top layer index

   CALL layer_hp2
   CALL TopLayerIndexp2

   ! ... Define variable values at time n+1/2.
   !     Only define values at cells that at n+1/2
   !     are wett.
   DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
     l = id_column(liter)

     ! ... Map 3D-(i,j) from 2D-l indexes
     i = l2i(l); j = l2j(l);

     ! ... At s-points
     kms = kmz(l)
     DO k = k1, kms
       salpp(k,l) = salp(k,l);
       salp (k,l)=(sal(k,l)+salpp(k,l))/2.
       rhop (k,l)=densty_s(salp(k,l),t0)-1000.
     ENDDO

     ! ... At u-points
     IF (mask2d(i+1,j)) THEN
       kmx = MIN(kmz(l),kmz(lEC(l)))
       DO k = k1, kmx
         uhpp(k,l) = uhp(k,l)
         upp (k,l) = up(k,l)
         IF (hup(k,l)>ZERO) THEN
           uhp (k,l) = 0.5*(uh(k,l) + uhpp(k,l))
           up  (k,l) = uhp(k,l)/hup(k,l)
         ELSE
           uhp (k,l) = 0.0
           up  (k,l) = 0.0
         ENDIF
       ENDDO

       ! ... Redo near surface flux calcs. at n+1/2
       k = k1u(l);

       ! a. Wetting occurs from n+1/2 to n+1
       IF (hu  (k-1,l) > ZERO) THEN
         uhp(k,l) = uhp(k,l)+uh(k-1,l)/2.
         up (k,l) = uhp(k,l) / hup(k,l)
       ENDIF

       ! b. Drying occurs from n to n+1/2
       IF (hupp(k-1,l) > ZERO) THEN
         uhp (k  ,l) = uhp (k,l)+uhpp(k-1,l)/2.
         up  (k  ,l) = uhp(k,l) / hup(k,l)
       ENDIF

     ENDIF

     ! ... At v-points
     IF (mask2d(i,j+1)) THEN
       kmy = MIN(kmz(l),kmz(lNC(l)))
       DO k = k1, kmy
         vhpp(k,l) = vhp(k,l)
         vpp (k,l) = vp(k,l)
         IF (hvp(k,l)>ZERO) THEN
           vhp (k,l) = 0.5*(vh(k,l) + vhpp(k,l))
           vp  (k,l) = vhp(k,l)/hvp(k,l)
         ELSE
           vhp (k,l) = 0.0
           vp  (k,l) = 0.0
         ENDIF
       ENDDO

       ! ... Redo near surface flux calcs. at n+1/2
       k = k1v(l);

       ! a. Wetting occurs from n+1/2 to n+1
       IF (hv  (k-1,l) > ZERO) THEN
         vhp(k,l) = vhp(k,l)+vh(k-1,l)/2.
         vp (k,l) = vhp(k,l) / hvp(k,l)
       ENDIF

       ! b. Drying occurs from n to n+1/2
       IF (hvpp(k-1,l) > ZERO) THEN
         vhp (k,l) = vhp (k,l)+vhpp(k-1,l)/2.
         vp  (k,l) = vhp (k,l) /  hvp (k,l)
       ENDIF

     ENDIF

   ENDDO
!$omp barrier
   !.....Recalculate vertical velocity at n+1/2 -  used in
   !     computing horizontal fluxes at n+1
   CALL continuity(2)

   ! ... Save bndry. variables from n-1 into n
   IF (nopenH(ide_t) > 0) THEN
    do nn=1,nopenH(ide_t)
    no = noh2no(nn,ide_t)
    SELECT CASE ( iside(no) )

         ! West boundary
         CASE (1)
            i = isbcHH(no,ide_t); js = jsbcH(no,ide_t); je = jebcH(no,ide_t)

            DO j = js, je
               uhEBpp(:,j) = uhEBp(:,j); huEBpp(:,j) = huEBp(:,j);
               uhWBpp(:,j) = uhWBp(:,j); huWBpp(:,j) = huWBp(:,j);
               uhEBp(:,j)  = (uhEBpp(:,j) + uhEB(:,j))/2.; huEBp(:,j) = (huEBp(:,j)+huEB(:,j))/2.;
               uhWBp(:,j)  = (uhWBpp(:,j) + uhWB(:,j))/2.; huWBp(:,j) = (huWBp(:,j)+huWB(:,j))/2.;
            END DO

         ! North boundary
         CASE (2)
            j = jsbcH(no,ide_t); is = isbcHH(no,ide_t); ie = iebcHH(no,ide_t)

            DO i = is, ie
               vhNBpp(:,i) = vhNBp(:,i); hvNBpp(:,i) = hvNBp(:,i);
               vhSBpp(:,i) = vhSBp(:,i); hvSBpp(:,i) = hvSBp(:,i);
               vhNBp(:,i)  = (vhNBpp(:,i) + vhNB(:,i))/2.; hvNBp(:,i) = (hvNBp(:,i)+hvNB(:,i))/2.;
               vhSBp(:,i)  = (vhSBpp(:,i) + vhSB(:,i))/2.; hvSBp(:,i) = (hvSBp(:,i)+hvSB(:,i))/2.;
            END DO

         ! East boundary
         CASE (3)
            i = isbcHH(no,ide_t); js = jsbcH(no,ide_t); je = jebcH(no,ide_t)

            DO j = js, je
               uhEBpp(:,j) = uhEBp(:,j); huEBpp(:,j) = huEBp(:,j);
               uhWBpp(:,j) = uhWBp(:,j); huWBpp(:,j) = huWBp(:,j);
               uhEBp(:,j)  = (uhEBpp(:,j) + uhEB(:,j))/2.; huEBp(:,j) = (huEBp(:,j)+huEB(:,j))/2.;
               uhWBp(:,j)  = (uhWBpp(:,j) + uhWB(:,j))/2.; huWBp(:,j) = (huWBp(:,j)+huWB(:,j))/2.;
            END DO

         ! South boundary
         CASE (4)
            j = jsbcH(no,ide_t); is = isbcHH(no,ide_t); ie = iebcHH(no,ide_t)

            DO i = is, ie
               vhNBpp(:,i) = vhNBp(:,i); hvNBpp(:,i) = hvNBp(:,i);
               vhSBpp(:,i) = vhSBp(:,i); hvSBpp(:,i) = hvSBp(:,i);
               vhNBp(:,i)  = (vhNBpp(:,i) + vhNB(:,i))/2.; hvNBp(:,i) = (hvNBp(:,i)+hvNB(:,i))/2.;
               vhSBp(:,i)  = (vhSBpp(:,i) + vhSB(:,i))/2.; hvSBp(:,i) = (hvSBp(:,i)+hvSB(:,i))/2.;
            END DO
         END SELECT
    end do
    end if

   ! ... Work with turbulence quantities (TurbModel)
   IF (iturb>0) CALL settrap_2EqTVars

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_settrap = t_settrap + (etime - btime)

END SUBROUTINE settrap

!***********************************************************************
SUBROUTINE save(istep)
!***********************************************************************
!
!  Purpose: Save solution for next time step
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-----------------------------------------------------------------------

   INTEGER,INTENT(IN) :: istep


   !.....Local variables.....
   INTEGER :: i, j, k, l, kms, k1s,liter,nn,no,js,je,is,ie,ide_t
   REAL    :: uutemp, vvtemp

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)


   ide_t = omp_get_thread_num ( )+1
   SELECT CASE (istep)

   !                   -----After a trapezoidal step (istep=2)-----
   CASE (2)

     sp(lm1)=s(lm1)
     hp(:,lm1)=h(:,lm1)
     hup(:,lm1)=hu(:,lm1)
     hvp(:,lm1)=hv(:,lm1)
     DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

       l = id_column(liter)
     !.....Save zeta.....
     sp(l) = s(l)

     ! ... Save layer thicknesses (h is calculated from s in solver)
     hp(:,l)  = h(:,l) ;
     hup(:,l) = hu(:,l);
     hvp(:,l) = hv(:,l);
     end do

     ! ... Retrieve index for surface layer for next step
     CALL TopLayerIndexp2

     ! ... Save state variables
      DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

       l = id_column(liter)

       ! ... Map 2D-l into 3D-(i,j) indexes
       i = l2i(l); j = l2j(l)

       DO k = k1, km;
         ! ... At s-points
         salp (k,l) = sal(k,l)
         rhop (k,l) = densty_s ( salp(k,l), t0 ) - 1000.
       ENDDO

       DO k = k1, km;
         ! ... At u-points
         IF (hup(k,l)>ZERO)THEN
           uhp (k,l) = uh (k,l) ! For horiz. scalar & momentum advection
           up  (k,l) = uhp(k,l)/hup(k,l) ! For horiz. momentum advection
           u   (k,l) = up (k,l) ! For output purposes
         ELSE
           uhp (k,l) = 0.0
           up  (k,l) = 0.0
           u   (k,l) = 0.0
         ENDIF
       ENDDO

       DO k = k1, km;
         ! ... At v-points
         IF (hvp(k,l)>ZERO)THEN
           vhp (k,l) = vh (k,l) ! For horiz. scalar & momentum advection
           vp  (k,l) = vhp(k,l)/hvp(k,l) ! For horiz. momentum advection
           v   (k,l) = vp (k,l) ! For output purposes
         ELSE
           vhp (k,l) = 0.0
           vp  (k,l) = 0.0
           v   (k,l) = 0.0
         ENDIF
       ENDDO

     END DO;
!$omp barrier
     !.....Recalculate vertical velocity wp to be used
     !     in calcuation of velocities at next time step
     CALL continuity(2)

     ! ... Save bndry. variables from n-1 into n
     IF (nopenH(ide_t) > 0) THEN
       do nn=1,nopenH(ide_t)
    no = noh2no(nn,ide_t)
    SELECT CASE ( iside(no) )

         ! West boundary
         CASE (1)
            i = isbcHH(no,ide_t); js = jsbcH(no,ide_t); je = jebcH(no,ide_t)

            DO j = js, je
               uhEBp(:,j)  = uhEB(:,j) ; huEBp(:,j) = huEB(:,j) ;
               uhWBp(:,j)  = uhWB(:,j) ; huWBp(:,j) = huWB(:,j) ;
            END DO

         ! North boundary
         CASE (2)
            j = jsbcH(no,ide_t); is = isbcHH(no,ide_t); ie = iebcHH(no,ide_t)

            DO i = is, ie
               vhNBp(:,i)  = vhNB(:,i) ; hvNBp(:,i) = hvNB(:,i) ;
               vhSBp(:,i)  = vhSB(:,i) ; hvSBp(:,i) = hvSB(:,i) ;
            END DO

         ! East boundary
         CASE (3)
            i = isbcHH(no,ide_t); js = jsbcH(no,ide_t); je = jebcH(no,ide_t)


            DO j = js, je
               uhEBp(:,j)  = uhEB(:,j) ; huEBp(:,j) = huEB(:,j) ;
               uhWBp(:,j)  = uhWB(:,j) ; huWBp(:,j) = huWB(:,j) ;
            END DO

         ! South boundary
         CASE (4)
            j = jsbcH(no,ide_t); is = isbcHH(no,ide_t); ie = iebcHH(no,ide_t)

            DO i = is, ie
               vhNBp(:,i)  = vhNB(:,i) ; hvNBp(:,i) = hvNB(:,i) ;
               vhSBp(:,i)  = vhSB(:,i) ; hvSBp(:,i) = hvSB(:,i) ;
            END DO
         END SELECT
    end do
end if

     ! ... Save Turbulence variables
     IF (iturb>0) CALL save_2EqTVars(istep)

     ! ... Save tracers
     IF (ntr > 0) THEN
      tracerpp = tracer;
     ENDIF


   !                   -----After a leapfrog step (istep=1)-----
   CASE (1)

     sp(lm1)=s(lm1)
     spp(lm1)=spp(lm1)
     hp(:,lm1)=h(:,lm1)
     hup(:,lm1)=hu(:,lm1)
     hvp(:,lm1)=hv(:,lm1)
     hpp(:,lm1)=hp(:,lm1)
     hupp(:,lm1)=hup(:,lm1)
     hvpp(:,lm1)=hvp(:,lm1)
     DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

       l = id_column(liter)
     !.....Save zeta.....
     spp(l)   = sp(l)
     sp(l)    = s(l)

     !.....Save layer thickness at time n
     hpp(:,l)  = hp(:,l)
     hupp(:,l) = hup(:,l)
     hvpp(:,l) = hvp(:,l)

     ! ... Save layer thickness at time n+1
     hp(:,l)  = h(:,l)
     hup(:,l) = hu(:,l)
     hvp(:,l) = hv(:,l)
     end do
!$omp barrier
     ! ... Retrieve index for surface layer for next step
     CALL TopLayerIndexp2

     ! .... Save other variables
     DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

       l = id_column(liter)

       ! ... Maks l- into (i,j)-indexes
!       i = l2i(l); j = l2j(l);

       DO k = k1, km;
         ! ... At s-points
         salpp(k,l) = salp(k,l)
         salp (k,l) = sal (k,l);
         rhop (k,l) = densty_s ( salp(k,l), t0 ) - 1000.
       ENDDO

       DO k = k1, km;
         ! ... At u- and v- points
         uhpp(k,l) = uhp(k,l)
         vhpp(k,l) = vhp(k,l)
         upp (k,l) = up (k,l)
         vpp (k,l) = vp (k,l)
       ENDDO

       DO k = k1, km;
         ! ... At u-points
         IF (hup(k,l)>ZERO)THEN
           uhp (k,l) = uh (k,l)
           up  (k,l) = uhp(k,l)/hup(k,l)
           u   (k,l) = up (k,l) ! For output purposes
         ELSE
           uhp (k,l) = 0.0
           up  (k,l) = 0.0
           u   (k,l) = 0.0
         ENDIF
       ENDDO

       DO k = k1, km;
         ! ... At v-points
         IF (hvp(k,l)>ZERO)THEN
           vhp (k,l) = vh (k,l)
           vp  (k,l) = vhp(k,l)/hvp(k,l)
           v   (k,l) = vp (k,l) ! For output purposes
         ELSE
           vhp (k,l) = 0.0
           vp  (k,l) = 0.0
           v   (k,l) = 0.0
         ENDIF
       END DO

     END DO
!$omp barrier
     !.....Recalculate vertical velocity wp to be used
     !     in calcuation of velocities at next time step
     CALL continuity(2)

     ! ... Save bndry. variables
     IF (nopenH(ide_t) > 0) THEN
       do nn=1,nopenH(ide_t)
    no = noh2no(nn,ide_t)
    SELECT CASE ( iside(no) )

         ! West boundary
         CASE (1)
            i = isbcHH(no,ide_t); js = jsbcH(no,ide_t); je = jebcH(no,ide_t)

            DO j = js, je
               uhEBpp(:,j) = uhEBp(:,j); huEBpp(:,j) = huEBp(:,j);
               uhWBpp(:,j) = uhWBp(:,j); huWBpp(:,j) = huWBp(:,j);
               uhEBp(:,j)  = uhEB(:,j) ; huEBp(:,j)  = huEB(:,j) ;
               uhWBp(:,j)  = uhWB(:,j) ; huWBp(:,j)  = huWB(:,j) ;
            END DO

         ! North boundary
         CASE (2)
            j = jsbcH(no,ide_t); is = isbcHH(no,ide_t); ie = iebcHH(no,ide_t)

            DO i = is, ie
               vhNBpp(:,i) = vhNBp(:,i); hvNBpp(:,i) = hvNBp(:,i);
               vhSBpp(:,i) = vhSBp(:,i); hvSBpp(:,i) = hvSBp(:,i);
               vhNBp(:,i)  = vhNB(:,i) ; hvNBp(:,i)  = hvNB(:,i) ;
               vhSBp(:,i)  = vhSB(:,i) ; hvSBp(:,i)  = hvSB(:,i) ;
            END DO

         ! East boundary
         CASE (3)
            i = isbcHH(no,ide_t); js = jsbcH(no,ide_t); je = jebcH(no,ide_t)

            DO j = js, je
               uhEBpp(:,j) = uhEBp(:,j); huEBpp(:,j) = huEBp(:,j);
               uhWBpp(:,j) = uhWBp(:,j); huWBpp(:,j) = huWBp(:,j);
               uhEBp(:,j)  = uhEB(:,j) ; huEBp(:,j)  = huEB(:,j) ;
               uhWBp(:,j)  = uhWB(:,j) ; huWBp(:,j)  = huWB(:,j) ;
            END DO

         ! South boundary
         CASE (4)
            j = jsbcH(no,ide_t); is = isbcHH(no,ide_t); ie = iebcHH(no,ide_t)

            DO i = is, ie
               vhNBpp(:,i) = vhNBp(:,i); hvNBpp(:,i) = hvNBp(:,i);
               vhSBpp(:,i) = vhSBp(:,i); hvSBpp(:,i) = hvSBp(:,i);
               vhNBp(:,i)  = vhNB(:,i) ; hvNBp(:,i)  = hvNB(:,i) ;
               vhSBp(:,i)  = vhSB(:,i) ; hvSBp(:,i)  = hvSB(:,i) ;
            END DO
         END SELECT
    end do
end if

     ! ... Save Turbulence variables
     IF (iturb>0) CALL save_2EqTVars(istep)

     ! ... Save tracers
     IF (ntr > 0) THEN
      tracerpp = tracer;
     ENDIF

   CASE DEFAULT

      PRINT *, "Invalid value of ISTEP in subroutine save"

   END SELECT

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_save = t_save + (etime - btime)

END SUBROUTINE save

!***********************************************************************
SUBROUTINE settrap2
!***********************************************************************
!
!  Purpose: To setup the arrays at the n+1/2 time level for use in
!           the second and subsequent iterations of the trapezoidal
!           step. Do not use smoothing as in the version of the code
!           by P.E. Smith.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-----------------------------------------------------------------------




   !.....Local variables.....
   INTEGER :: i,j,k,l,kms,kmy,kmx,lol,liter,aux_indice
   REAL    :: wght, wghtpp

   !....Zeta array.....
   sp(lm1)=0.5*(s(lm1)+spp(lm1))
   IF(omp_get_thread_num ( )+1 == num_threads)THEN
   aux_indice=lhf(omp_get_thread_num ( )+1)
   ELSE
   aux_indice=lhfE(omp_get_thread_num ( )+1)
   END IF
   DO liter = lhi(omp_get_thread_num ( )+1),aux_indice
     l = id_column(liter)
   !....Zeta array.....
   sp(l) = 0.5*(s(l) + spp(l))
   end do

   ! ...Define layers thickness at time n+1/2

   CALL layer_hp2
   CALL TopLayerIndexp2

   ! ... 3d arrays
   DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

       l = id_column(liter)

     ! ... Map 3D-(i,j) from 2D-l indexes
     i = l2i(l); j = l2j(l);

     ! ... At s-points
     kms = kmz(l)

     DO k = k1, kms
       salp (k,l)= (sal(k,l)+salpp(k,l))/2.
       rhop (k,l)=densty_s(salp(k,l),t0)-1000.
     ENDDO


     ! ... At u-points
     IF (mask2d(i+1,j)) THEN
       kmx = MIN(kmz(l),kmz(lEC(l)))
       DO k = k1, kmx
         IF (hup(k,l)>ZERO) THEN
           uhp (k,l) = 0.5*(uh(k,l) + uhpp(k,l))
           up  (k,l) = uhp(k,l)/hup(k,l)
         ELSE
           uhp (k,l) = 0.0
           up  (k,l) = 0.0
         ENDIF
       ENDDO

       ! ... Redo near surface flux calcs.
       k = k1u(l);

       ! a. Wetting occurs from n+1/2 to n+1
       IF (hu  (k-1,l) > ZERO) THEN
         uhp(k,l) = uhp(k,l)+uh(k-1,l)/2.
         up (k,l) = uhp(k,l) / hup(k,l)
       ENDIF

       ! b. Drying occurs from n to n+1/2
       IF (hupp(k-1,l) > ZERO) THEN
         uhp (k,l) = uhp (k,l)+uhpp(k-1,l)/2.
         up  (k,l) = uhp(k,l) / hup(k,l)
       ENDIF

     ENDIF

     ! ... At v-points
     IF (mask2d(i,j+1)) THEN
       kmy = MIN(kmz(l),kmz(lNC(l)))
       DO k = k1, kmy
         IF (hvp(k,l)>ZERO) THEN
           vhp (k,l) = 0.5*(vh(k,l) + vhpp(k,l))
           vp  (k,l) = vhp(k,l)/hvp(k,l)
         ELSE
           vhp (k,l) = 0.0
           vp  (k,l) = 0.0
         ENDIF
       ENDDO

       ! ... Redo near surface flux calcs.
       k = k1v(l);

       ! a. Wetting occurs from n+1/2 to n+1
       IF (hv  (k-1,l) > ZERO) THEN
         vhp(k,l) = vhp(k,l)+ vh(k-1,l)/2.
         vp (k,l) = vhp(k,l) / hvp(k,l)
       ENDIF

       ! b. Drying occurs from n to n+1/2
       IF (hvpp(k-1,l) > ZERO) THEN
         vhp (k,l) = vhp (k,l)+vhpp(k-1,l)/2.
         vp  (k,l) = vhp (k,l) /  hvp (k,l)
       ENDIF

     ENDIF

   ENDDO
!$omp barrier

   CALL continuity(2)
   !.....Recalculate vertical velocity components to be used in
   !     calculating horizontal velocity at next iteration
!   DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

!       l = id_column(liter)
!      i = l2i(l); j = l2j(l);
!      kms = kmz(l);
!      wp(kms+1,l) = 0.0
!      DO k = kms,k1,-1
!        wp(k,l) = wp(k+1,l)-(uhp(k,l)-uhp(k,lWC(l)))/dx    &
!                 &         -(vhp(k,l)-vhp(k,lSC(l)))/dy
!      END DO
!   END DO

   ! ... Work with turbulence quantities (TurbModel)
   IF (iturb>0) CALL settrap2_2EqTVars

END SUBROUTINE settrap2

!***********************************************************************
SUBROUTINE exTracer  (nt,Bstart,Bend,Bhaxpp,Bhaypp,Bth3,Bth4,Bth2,lSCH,lNCH,lECH,lWCH,Bex,thrs)
!***********************************************************************
!
!  Purpose: To solve transport equation for tracer, using Flux-limiters.
!           nt denotes tracer number to be solved
!
!-----------------------------------------------------------------------

   ! ... Arguments
   INTEGER, INTENT (IN) :: nt
   REAL, INTENT(IN) :: thrs
   INTEGER, INTENT(IN) :: Bstart,Bend
   REAL, DIMENSION (1:km1,Bstart:Bend+1), INTENT(INOUT) :: Bhaxpp, Bhaypp,Bth3,Bth4,Bth2,Bex
   INTEGER, DIMENSION (Bstart:Bend+1), INTENT(IN) :: lSCH,lNCH,lWCH,lECH

   ! ... Local variables
   INTEGER :: i, j, k, l, k1s, kms, gamma1, istat,liter
   REAL    :: vel, ratio, C_f, delz, twodt1, hd
   REAL, DIMENSION (4          ) :: ss


   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)



   ! ... Constants used in solution
   twodt1 = twodt*tz

   ! ... Calculate hdxpp & hdypp arrays for diffusion terms &
   Bhaxpp(:,Bend+1) = 0.0; Bhaypp(:,Bend+1) = 0.0;
   DO liter = lhiW(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
      If(liter == 0) CYCLE

      l = id_column(liter)

      ! ... 3D-(i,j) indexes for l
!      i = l2i(l); j = l2j(l);

      ! ... Retrieve top & bottom wet sal-pts .................
      kms = kmz(l)
      k1s = k1z(l)

      ! ... Calculate hdxpp & hdypp array at u-&v- pts ........
      DO k = k1, kms
        Bhaxpp(k,l) = Ax0*hupp(k,l)
        Bhaypp(k,l) = Ay0*hvpp(k,l)
      ENDDO

   END DO

   Bth3(:,Bend+1) = 0.0
   ! ... Initialize ex & flux arrays to zeros
   Bex(:,Bend+1) = 0.0; Bth2(:,Bend+1)= 0.0; Bth4(:,Bend+1) = 0.0;

   DO liter = lhiW(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)
     If(liter == 0) CYCLE

     l = id_column(liter)

     ! ... Map l- into (i,j)-indexes .........................
!     i = l2i(l); j = l2j(l);

     ! ... Retrieve top & bottom wet sal-pts .................
     kms = kmz(l)
     k1s = k1z(l)

     DO k = k1s, kms;

       ! ... EW fluxes .......................................
       IF (hup(k,l)> ZERO) THEN

         ! ... Velocities at time n+1/2
         vel  = (uhpp(k,l) + uh(k,l))/2.

         ! ... Define stencil for scalar transport
         ss(2)  = tracerpp(k,    l ,nt);
         ss(3)  = tracerpp(k,lEC(l),nt);
         IF (hpp(k,lWC(l))<=ZERO)THEN; ss(1)=ss(2);
           ELSE; ss(1)=tracerpp(k,lWC(l),nt); ENDIF
         IF (hpp(k,lEC(lEC(l)))<=ZERO)THEN; ss(4)=ss(3);
           ELSE; ss(4)=tracerpp(k,lEC(lEC(l)),nt); ENDIF

         ! ... Calculate Cf factor to use in flux calculations
         gamma1 = -SIGN (1., vel)
         C_f    = 0.0;
         IF (ss(3) - ss(2) /= 0 ) THEN
           ratio =(ss(3+gamma1)-ss(2+gamma1))/(ss(3)-ss(2))
           ! MC flux limiter (VanLeer, 1977)
           !C_f   = MAX(0., MIN( 2*ratio, (1+ratio)/2., 2. ))
           ! ... Roe's Superbee Limiter
           C_f = MAX(0., MIN(1.,2.*ratio),MIN(2.,ratio))
         ENDIF

         ! ... Calculate fluxes at x-faces
         Bth2(k,l) = vel/2.*(ss(3)+ss(2))- &
         & ((1.-C_f)*ABS(vel)+vel**2.*twodt1/dx*C_f)*(ss(3)-ss(2))/2.
       ELSE
         Bth2(k,l) = 0.0
       ENDIF
       ENDDO;
   ENDDO;
   DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

     l = id_column(liter)

     ! ... Map l- into (i,j)-indexes .........................
!     i = l2i(l); j = l2j(l);

     ! ... Retrieve top & bottom wet sal-pts .................
     kms = kmz(l)
     k1s = k1z(l)

     DO k = k1s, kms;

       ! ... NS fluxes .......................................
       IF (hvp(k,l)> ZERO) THEN

         ! ... Velocities at time n+1/2
         vel  = (vhpp(k,l) + vh(k,l))/2.

         ! ... Define stencil for scalar transport
         ss(2)  = tracerpp(k,        l  ,nt );
         ss(3)  = tracerpp(k,    lNC(l) ,nt );
         IF (hpp(k,lSC(l))<=ZERO)THEN; ss(1)=ss(2);
           ELSE; ss(1)=tracerpp(k,lSC(l),nt); ENDIF
         IF (hpp(k,lNC(lNC(l)))<=ZERO)THEN; ss(4)=ss(3);
           ELSE; ss(4)=tracerpp(k,lNC(lNC(l)),nt); ENDIF

         ! ... Calculate Cf factor to use in flux calculations
         C_f    = 0.0; ! Default value is for upwinding
         gamma1 = -SIGN (1., vel)
         IF (ss(3) - ss(2) /= 0 ) THEN
           ratio =(ss(3+gamma1)-ss(2+gamma1))/(ss(3)-ss(2))
           ! MC flux limiter (VanLeer, 1977)
           !C_f   = MAX(0., MIN( 2*ratio, (1+ratio)/2., 2. ))
           ! ... Roe's Superbee Limiter
           C_f = MAX(0., MIN(1.,2.*ratio),MIN(2.,ratio))
         ENDIF

         ! ... Calculate fluxes at y-faces
         Bth4(k,l) = vel/2.*(ss(3)+ss(2))- &
         & ((1.-C_f)*ABS(vel)+vel**2.*twodt1/dx*C_f)*(ss(3)-ss(2))/2.
       ELSE
         Bth4(k,l) = 0.0
       ENDIF

       ! ... UD fluxes .......................................
       IF (hp(k-1,l)> ZERO) THEN

         ! ... Velocities at time n+1/2 (include settling velocity)
         vel  = wp(k,l) ; IF (k == k1s) vel = 0.0;

         ! ... Define stencil for scalar transport
         ss(2)  = tracerpp(k,l,nt);
         ss(3)  = tracerpp(k-1,l,nt)
         IF(hpp(k-2,l)<=ZERO)THEN;ss(4)=ss(3);
           ELSE;ss(4)=tracerpp(k-2,l,nt);ENDIF
         IF(hpp(k+1,l)<=ZERO)THEN;ss(1)=ss(2);
           ELSE;ss(1)=tracerpp(k+1,l,nt);ENDIF;
         ! ... Calculate ratio of slope of solution across interfaces &
         !     estimate flux limiter
         C_f = 1.0 ! Default method is Lax-Wendroff
         gamma1 = -SIGN (1., vel)
         IF (ss(3) - ss(2) /= 0 ) THEN
           ratio =(ss(3+gamma1)-ss(2+gamma1))/(ss(3)-ss(2))
           ! MC flux limiter (VanLeer, 1977)
           ! C_f   = MAX(0., MIN( 2*ratio, (1+ratio)/2., 2. ))
           ! ... Roe's Superbee Limiter
           C_f = MAX(0., MIN(1.,2.*ratio),MIN(2.,ratio))
         ENDIF

         ! ... Calculate fluxes at z-faces
         delz = (hp(k,l) + hp(k-1,l))/2.
         Bth3(k,l) = vel/2.*(ss(3)+ss(2))- &
         & ((1.-C_f)*ABS(vel)+vel**2.*twodt1/delz*C_f)*(ss(3)-ss(2))/2.
       ELSE
         Bth3(k,l) = 0.0
       ENDIF

     ENDDO;
   ENDDO;

   ! ... Update ex array with divergence of advective fluxes & diffusion
   DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

     l = id_column(liter)

     ! ... Map l- into (i,j)-indexes .........................
!     i = l2i(l); j = l2j(l);

     ! ... Retrieve top & bottom wet sal-pts .................
     kms = kmz(l)
     k1s = k1z(l)

     DO k = k1s, kms;

       !.....Horizontal diffusion.....
       hd= (Bhaxpp(k,    l )*(tracerpp(k,lEC(l),nt) - tracerpp(k,    l ,nt))       &
         & -Bhaxpp(k,lWCH(l))*(tracerpp(k,    l ,nt) - tracerpp(k,lWC(l),nt)))/dxdx &
         &+(Bhaypp(k,    l )*(tracerpp(k,lNC(l),nt) - tracerpp(k,    l ,nt))       &
         & -Bhaypp(k,lSCH(l))*(tracerpp(k,    l ,nt) - tracerpp(k,lSC(l),nt)))/dydy

       Bex(k,l) =   hpp(k,l)*tracerpp(k,l,nt)/twodt1   &
                - (Bth2(k,l) - Bth2(k,lWCH(l))) / dx &
                - (Bth4(k,l) - Bth4(k,lSCH(l))) / dy &
                - (Bth3(k,l) - Bth3(k+1,l   )) !+ hd * ihd
       IF (ihd>0) Bex(k,l) = Bex(k,l) + hd   ! Changed 12/2010 SWA

     ENDDO;
   ENDDO

   ! ... Modify explicit term to account for flow boundary conditions
   CALL MODextracer4openbc (nt,Bstart,Bend,Bex,thrs)



   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_exsal = t_exsal + (etime - btime)

END SUBROUTINE ExTracer

!***********************************************************************
SUBROUTINE ImTracer (nt,Bstart,Bend,Bex)
!***********************************************************************
!
!  Purpose: To solve for active scalar concentration.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-----------------------------------------------------------------------


   ! ... Arguments
   INTEGER, INTENT (IN) :: nt
   INTEGER, INTENT(IN) :: Bstart,Bend
   REAL, DIMENSION (1:km1,Bstart:Bend+1), INTENT(INOUT) :: Bex

   !.....Local variables.....
   REAL :: twodt1, Osource, Qsource
   INTEGER :: i, j, k, l, k1s, kms, kt, nwlayers, inn, kk, noc,liter,innH
   REAL, DIMENSION (1:km1) :: hn
   REAL, DIMENSION (3,1:km1) :: aa
   REAL, DIMENSION (1:km) :: ds
   REAL, DIMENSION (1:ndz) :: sal1

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)

   ! ... Constants used in solution
   twodt1 = twodt*tz

   !.....Loop over interior sal-pts to solve for
   !     matrix from the active scalar equation.....
   DO liter = lhi(omp_get_thread_num ( )+1), lhf(omp_get_thread_num ( )+1)

      l = id_column(liter)

      ! ... 3D-(i,j) indexes for l - FJR - uncomment
      i = l2i(l); j = l2j(l);

      !.....Compute top & bottom layer numbers & No. of layers ....
      kms = kmz(l)
      k1s = k1z(l)
      nwlayers = (kms-k1s)+1

      ! ... Define layer thikness at time n - The corrections for
      !     surface and recently submerged cells are needed to
      !     keep mass conservation
      hn(k1s+1:kms) = h(k1s+1:kms,l)
      hn(k1s      ) = twodt1*wp(k1s,l)+hpp(k1s,l)
      IF (hpp(k1s,l)<= ZERO) THEN
        hn(k1s+1) = hpp(k1s+1,l)
      ENDIF

      !.....Calculate active scalar for case of a single layer.....
      SELECT CASE (nwlayers)
      CASE (1)

        aa( 2,k1s) = hn(k1s)/twodt1
        ds(   k1s) = Bex(k1s,l)
        tracer(k1s,l,nt) = ds(k1s  )/aa(2,k1s)

      !.....Calculate active scalar for case of two or more layers.....
      CASE (2:)

         !.....Form coefficient matrix [aa]
         ! Define upper diagonal terms
         aa(3,k1s:kms-1) = -Dv(k1s+1:kms,l)/(hn(k1s:kms-1)+hn(k1s+1:kms))*2.
         aa(3,kms)       =  0.0
         ! Define lower diagonal terms
         aa(1,k1s+1:kms) = -Dv(k1s+1:kms,l)/(hn(k1s:kms-1)+hn(k1s+1:kms))*2.
         aa(1,k1s)       =  0.0
         ! Define center diagonal terms
         aa(2,k1s:kms)   =  hn(k1s:kms)/twodt1-aa(1,k1s:kms)-aa(3,k1s:kms)

         !.....form r.h.s. matrix [ds].....
         DO k = k1s, kms
            ds(k) = Bex(k,l)
         ENDDO

         ! ... Modify transport eqs. to accont for sources & sinks.
         IF ( iopssH(omp_get_thread_num ( )+1) > 0 ) THEN
           DO innH = 1, iopssH(omp_get_thread_num ( )+1)
             inn = ioph2iop(innH,omp_get_thread_num ( )+1)
             IF ( j /= jpss(inn) .OR. i /=ipss(inn) ) CYCLE
             DO k = k1s, kms
               IF (ABS(Qpss(k,inn))<1.E-10) CYCLE

               Qsource  = Qpss(k,inn)/(dx*dy)  ! Inflow per unit area (m/s)
               Osource  = Rpss(k,inn,nt)       ! Concentration (kg/m3)
               ds(k)=ds(k)+Qsource*Osource     ! kg/m2/s = conc.* thickness / time
             ENDDO
           ENDDO
           ! ... Include SOD when modelling oxygen plumes -
           IF (nt == ntr) ds(kms) = ds(kms) - k4sod

         ENDIF

         !.....Solve tridiagonal system for the
         !     vertical distribution of active scalar.....
         CALL trid1 (aa, ds, sal1, k1s, kms, km1, nwlayers)

         !.....Define scalars at new time step....
         tracer(k1s:kms  ,l,nt) = sal1(1:nwlayers)
         tracer(k1 :k1s-1,l,nt) = sal1(1         )

      END SELECT

   !.....End loop over scalar-pts.....
   END DO

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_salin = t_salin + (etime - btime)

END SUBROUTINE ImTracer

!***********************************************************************
SUBROUTINE ConfigThreads (depth)
!***********************************************************************
!
!  Purpose: to distribute the workload across the available threads,
!  indicating the number of columns assigned to each thread
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-----------------------------------------------------------------------

    !Arguments
    INTEGER, INTENT(INOUT) :: depth

    !Local Variables
    INTEGER :: p, counter, ind, aux, aux2, index_thread, id_thread
    INTEGER :: aux3, aux4, depth_aux, aux5, i, j

!   aux=int(im1/num_threads)
!   aux2=MOD(im1,num_threads)

!   lh(1:num_threads)=aux*jm1
!   lh_aux(1:num_threads)=aux

!   index_thread = 1
!   IF(aux2>0) THEN
!   		do p=aux2, 1, -1
!			lh(index_thread)=lh(index_thread) + jm1
!			lh_aux(index_thread)=lh_aux(index_thread) + 1
!   			index_thread=index_thread + 1
!   		end do
!   	END IF
   !to distribute the workload bases on the depth
   !solve the number of rows from east to west are assigned to each thread
   depth = 0
   do p=1,lm
        depth = depth + (kmz(p) - 1)
   end do
   depth_aux = int(depth/num_threads)
   print *,"depth_aux:",depth_aux
   aux = int(lm/num_threads)
   aux2 = 0
   aux3 = 0
   aux5 = 0
   ph(1:num_threads) = 0
   lh(1:num_threads) = 0
   lh_aux(1:num_threads) = 0
   DO p=1,num_threads-1
     if(p .EQ. 1) then
         aux4 = 1
     else
         aux4 = aux4 + lh_aux(p-1)
     end if
     print *,"aux4:",aux4,p
     DO i=aux4,im1
      	lh_aux(p)=lh_aux(p)+1
      	DO j=1,jm1
          	IF(mask2d(i,j)) THEN
          			lh(p) = lh(p) + 1
          			ph(p) = ph(p) + (kmz(ij2l(i,j)) - 1)
          	END IF
      	END DO
        print *,"ph:",ph(p),p
      	IF(ph(p) >= depth_aux) THEN
         	aux2 = aux2 + lh(p)
         	aux3 = aux3 + lh_aux(p)
         	aux5 = aux5 + ph(p)
                print *,"aux2:",aux2,"aux3:",aux3,"aux5:",aux5,p
         	EXIT
      	END IF
     END DO
   END DO
   lh(num_threads) = lm - aux2
   lh_aux(num_threads) = im - aux3
   ph(num_threads) = depth - aux5

!   lh_aux(1) = 47
!   lh_aux(2) = 16
!   lh_aux(3) = 17
!   lh_aux(4) = 29
!    lh_aux(1) = 51
!    lh_aux(2) = 12
!    lh_aux(3) = 22
!    lh_aux(4) = 24
   lh = 0
   ph = 0
   IF(num_threads > 1) THEN
   DO i=1,lh_aux(1)

      	DO j=1,jm1
          	IF(mask2d(i,j)) THEN
          			lh(1) = lh(1) + 1
          			ph(1) = ph(1) + (kmz(ij2l(i,j)) - 1)
          	END IF
      	END DO
   END DO
   DO i=lh_aux(1)+1,lh_aux(1)+lh_aux(2)

      	DO j=1,jm1
          	IF(mask2d(i,j)) THEN
          			lh(2) = lh(2) + 1
          			ph(2) = ph(2) + (kmz(ij2l(i,j)) - 1)
          	END IF
      	END DO
   END DO
   IF(num_threads > 2) THEN
   DO i=lh_aux(1)+lh_aux(2)+1,lh_aux(1)+lh_aux(2)+lh_aux(3)

      	DO j=1,jm1
          	IF(mask2d(i,j)) THEN
          			lh(3) = lh(3) + 1
          			ph(3) = ph(3) + (kmz(ij2l(i,j)) - 1)
          	END IF
      	END DO
   END DO
   DO i=lh_aux(1)+lh_aux(2)+lh_aux(3)+1,lh_aux(1)+lh_aux(2)+lh_aux(3)+lh_aux(4)

      	DO j=1,jm1
          	IF(mask2d(i,j)) THEN
          			lh(4) = lh(4) + 1
          			ph(4) = ph(4) + (kmz(ij2l(i,j)) - 1)
          	END IF

      	END DO
   END DO
   END IF
   END IF


   !to create the array id_column where stores the index of all columns
   !assigned to each thread even those belonging to the border east or west
   counter = 0
   ind = 1
   lhiE = 0
   lhfE = 0
   lhiW = 0
   lhfW = 0
   do p=1,num_threads
        if(p > 1)THEN
        lhiW(p)=counter+1
        do j=1,jm1
			if(mask2d(ind-1,j) .AND. mask2d(ind,j))THEN
			counter=counter+1
			id_column(counter)= ij2l(ind-1,j)
			end if
        end do
        lhfW(p)=counter
        end if
        lhi(p)=counter+1
   		do i=ind,lh_aux(p)+ind-1
   			do j=1,jm1

   			    IF(mask2d(i,j)) THEN
   				counter = counter + 1
   				id_column(counter) = ij2l(i,j)
   				END IF
   			end do
   		end do
   	ind = lh_aux(p) + ind
   	lhf(p)=counter
   	if(p < num_threads)THEN
   		lhiE(p)=counter+1
        	do j=1,jm1
			if(mask2d(ind,j) .AND. mask2d(ind-1,j))THEN
				counter=counter+1
				id_column(counter)= ij2l(ind,j)
			end if
        	end do
        	lhfE(p)=counter
        end if
   end do

   !As above but it must also satisfy that the east column is not dry
   counter = 0
   ind = 1
   lhiECE = 0
   lhfECE = 0
   lhiWCE = 0
   lhfWCE = 0
   do p=1,num_threads
        if(p > 1)THEN
        lhiWCE(p)=counter+1
        do j=1,jm1
			if(mask2d(ind-1,j) .AND. mask2d(ind,j))THEN
			counter=counter+1
			id_columnCE(counter)= ij2l(ind-1,j)
			end if
        end do
        lhfWCE(p)=counter
        end if
        lhiCE(p)=counter+1
   		do i=ind,lh_aux(p)+ind-1
   			do j=1,jm1
   			    IF(mask2d(i,j) .AND. mask2d(i+1,j)) THEN
   				counter = counter + 1
   				id_columnCE(counter) = ij2l(i,j)
   				END IF
   			end do
   		end do
   		ind = lh_aux(p) + ind
   		lhfCE(p)=counter
   		if(p < num_threads)THEN
   		lhiECE(p)=counter+1
        do j=1,jm1
			if(mask2d(ind,j) .AND. mask2d(ind-1,j) .AND. mask2d(ind+1,j))THEN
			counter=counter+1
			id_columnCE(counter)= ij2l(ind,j)
			end if
        end do
        lhfECE(p)=counter
        end if
   end do

   !As above but it must also satisfy that the north column is not dry
   counter = 0
   ind = 1
   lhiECN = 0
   lhfECN = 0
   lhiWCN = 0
   lhfWCN = 0
   do p=1,num_threads
        if(p > 1)THEN
        lhiWCN(p)=counter+1
        do j=1,jm1
			if(mask2d(ind-1,j) .AND. mask2d(ind,j) .AND. mask2d(ind-1,j+1))THEN
			counter=counter+1
			id_columnCN(counter)= ij2l(ind-1,j)
			end if
        end do
        lhfWCN(p)=counter
        end if
        lhiCN(p)=counter+1
   		do i=ind,lh_aux(p)+ind-1
   			do j=1,jm1
   			    IF(mask2d(i,j) .AND. mask2d(i,j+1)) THEN
   				counter = counter + 1
   				id_columnCN(counter) = ij2l(i,j)
   				END IF
   			end do
   		end do
   		ind = lh_aux(p) + ind
   		lhfCN(p)=counter
   		if(p < num_threads)THEN
   		lhiECN(p)=counter+1
        do j=1,jm1
			if(mask2d(ind,j) .AND. mask2d(ind-1,j) .AND. mask2d(ind,j+1))THEN
			counter=counter+1
			id_columnCN(counter)= ij2l(ind,j)
			end if
        end do
        lhfECN(p)=counter
        end if
   end do

END SUBROUTINE ConfigThreads

!***********************************************************************
SUBROUTINE BorderThreads(Bstart,Bend,ide_thread,mincol)
!***********************************************************************
!
!  Purpose: to set the private variables used by the threads,
!  dynamically allocates memory and initializes these variables
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-----------------------------------------------------------------------

   ! Arguments
   INTEGER, INTENT(INOUT) :: Bstart,Bend,ide_thread,mincol

   ! Local Variables
   INTEGER :: i,j

    Bstart = 0; Bend = 0;

    ! Set the first column of each thread
    if(lhiW(ide_thread) .EQ. 0) THEN
    Bstart = id_column(lhi(ide_thread))
    ELSE
    i = l2i(id_column(lhiW(ide_thread))) - 1
    mincol = id_column(lhiW(ide_thread))
    do j=1,jm1
   		if(ij2l(i,j) < mincol) mincol = ij2l(i,j)
    end do
    Bstart = mincol
    end if

    ! Set the last column of each thread
    if(lhiE(ide_thread) .EQ. 0) THEN
    Bend = id_column(lhf(ide_thread))
    ELSE
    Bend = id_column(lhfE(ide_thread))
    end if

    print *,"Bstart:",Bstart,"Bend:",Bend,"h:",omp_get_thread_num()

END SUBROUTINE BorderThreads

!***********************************************************************
SUBROUTINE InitThreads(t_exmom2,t_matmom2,t_matcon2,Bhaxpp,Bhaypp,Bth,Bth1,Bstart, &
    & Bend,ide_thread,lNCH,lSCH,lECH,lWCH,Bex, Bth2,Beagx,Bearx,Bagx,Barx,Beagy, &
    & Beary,Bagy,Bary,Bsx,Bsy,Bdd,Bqq,Brr,Bth3,Bth4,uhp2,uhp3,mincol,ShearProduction, &
    & BuoyancyProduction, Dissipation, TKinE,heatSourceB,uairB,vairB,cdwB,eta, &
    & QswFrB,Qsw2dB,Qlw2dB,Ta2dB,RH2dB,Cc2dB,uair2dB,vair2dB,Qsw,Qn,Qlw,Ta,Pa,RH,Cc)
!***********************************************************************
!
!  Purpose: to set the private variables used by the threads,
!  dynamically allocates memory and initializes these variables
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-----------------------------------------------------------------------

   ! Arguments
   REAL, INTENT(INOUT) :: t_exmom2,t_matmom2,t_matcon2
   INTEGER, INTENT(IN) :: Bstart,Bend,ide_thread,mincol
   REAL(real_G1), INTENT(INOUT) :: ShearProduction,BuoyancyProduction, Dissipation, TKinE
   REAL, DIMENSION (1:km1,1:lm1), INTENT(INOUT) :: uhp3,uhp2
   REAL, DIMENSION (1:km1,Bstart:Bend+1), INTENT(INOUT) :: Bhaxpp, Bhaypp, Bth, Bth1, Bex, heatSourceB,QswFrB
   REAL, DIMENSION (1:km1,Bstart:Bend+1), INTENT(INOUT) :: Bth2, Bagx, Barx,Bagy,Bary,Bth3,Bth4
   REAL, DIMENSION (Bstart:Bend+1), INTENT(INOUT) :: Beagx, Bearx,Beagy,Beary,Bsx,Bsy
   REAL, DIMENSION (Bstart:Bend+1), INTENT(INOUT) :: Bdd, Bqq,Brr,uairB,vairB,cdwB
   INTEGER, DIMENSION (Bstart:Bend+1), INTENT(INOUT) :: lSCH,lNCH,lWCH,lECH
   REAL, INTENT(INOUT) :: Qsw,Qn,Qlw,eta,Ta,Pa,RH,Cc
   REAL, DIMENSION (nmetstat), INTENT(INOUT) :: Qsw2dB,Qlw2dB,Ta2dB,RH2dB,Cc2dB,uair2dB,vair2dB

   ! Local Variables
   INTEGER :: i,j



    ! Make a copy of lWC, lSC, lEC, lNC to each thread
    do i=Bstart,Bend
    	if(lWC(i) .EQ. lm1) THEN
    	lWCH(i) = Bend + 1
    	ELSE
    	lWCH(i) = lWC(i)
    	end if
    	if(lNC(i) .EQ. lm1) THEN
    	lNCH(i) = Bend + 1
    	ELSE
    	lNCH(i) = lNC(i)
    	end if
    	if(lSC(i) .EQ. lm1) THEN
    	lSCH(i) = Bend + 1
    	ELSE
    	lSCH(i) = lSC(i)
    	end if
    	if(lEC(i) .EQ. lm1) THEN
    	lECH(i) = Bend + 1
    	ELSE
    	lECH(i) = lEC(i)
    	end if
    end do

    ! Initialize variables
     t_exmom2 = 0; t_matmom2 = 0; t_matcon2 = 0;
    Bth = 0; Bth1 = 0; Bex = 0; Bth2 = 0; Bth3 = 0; Bth4 = 0;
    Bhaxpp = 0; Bhaypp = 0;
    Beagx = 0; Bearx = 0; Bagx(:,Bend+1) = 0; Barx(:,Bend+1) = 0;
    Beagy = 0; Beary = 0; Bagy(:,Bend+1) = 0; Bary(:,Bend+1) = 0;
    Bsx = 0.0; Bsy(Bend+1) = 0; Bdd(Bend+1) = 0; Bqq(Bend+1) = 0; Brr(Bend+1) = 1.0;
    uhp2 = 0; uhp3 = 0;
    ShearProduction = 0; BuoyancyProduction= 0; Dissipation = 0; TKinE = 0;
    uairB(Bstart:Bend) = uair(Bstart:Bend); vairB(Bstart:Bend) = vair(Bstart:Bend);
    cdwB(Bstart:Bend) = cdw(Bstart:Bend); heatSourceB(:,Bstart:Bend) = heatSource(:,Bstart:Bend);
    uairB(Bend+1) = 0; vairB(Bend+1) = 0; cdwB(Bend+1) = 0; heatSourceB(:,Bend+1) = 0;
    QswFrB(:,Bstart:Bend) = QswFr(:,Bstart:Bend);
    Qsw2dB = Qsw2d; Qlw2dB = Qlw2d; Ta2dB = Ta2d; RH2dB = RH2d; uair2dB = uair2d; vair2dB = vair2d;
   !  Cc2dB = Cc2d;

END SUBROUTINE InitThreads

END MODULE si3d_procedures
