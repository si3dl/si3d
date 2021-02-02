
!*************************************************************************
                               PROGRAM si3d
!*************************************************************************
!
!  Purpose:  Main program for the semi-implicit 3-d hydrodynamic model.
!            The model uses a leapfrog-trapezoidal finite-difference
!            numerical scheme.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-------------------------------------------------------------------------
!                         U. S. Geological Survey
!                             Sacramento, CA
!-------------------------------------------------------------------------

   USE omp_lib
   USE si3d_procedures
   IMPLICIT NONE

   !.....External function declaration.....
   REAL, EXTERNAL :: TIMER             ! Use timing routine from NSPCG

   !.....Local variables.....
   REAL :: utime, stime, utime1, stime1, utime2, stime2, utime3, stime3,  &
         & ttime1, ttime2, atime, ltime1, ltime2
   REAL :: TimeStart, TimeEnd
   INTEGER :: maxcount = 1E6,nn,is,ie,js,je,noH,isH,ieH,la,laaux,contBC,istat,iboid,lcon,lcon2,ifrontera
   INTEGER :: iter, itemp, i, j, niter1, p, Bstart, Bend, ide_thread, depth, mincol,liter,k,aux_la,ios
   CHARACTER (LEN = 12) :: Ifile

   REAL :: t_exmom2, t_matmom2, t_matcon2, tsbar, tebar, dzino,uflow
   REAL, ALLOCATABLE, DIMENSION (:,:) :: Bhaxpp, Bhaypp, Bth, Bth1, Bth2, Bex, Bagx, Barx
   REAL, ALLOCATABLE, DIMENSION (:,:) :: Bagy, Bary, Bth3, Bth4, uhp2, uhp3, heatSourceB, QswFrB
   INTEGER, ALLOCATABLE, DIMENSION (:) :: lNCH, lSCH, lWCH, lECH
   REAL, ALLOCATABLE, DIMENSION (:) :: Beagx, Bearx, Beagy, Beary, Bsx, Bsy, Bdd, Bqq, Brr
   REAL, ALLOCATABLE, DIMENSION (:) :: Qsw2dB, Qlw2dB, Ta2dB, RH2dB, Cc2dB, uair2dB, vair2dB
   REAL, ALLOCATABLE, DIMENSION (:) :: uairB, vairB, cdwB,bclncxB,hupdrhoB


   !.....Retrieve begin time of run.....
   stime = TIMER(0.0)

   !.....Read input parameters.....
   IF (idbg == 1) PRINT *, " Before entry to SUB input"
   CALL input

   ! ... Output run parameters for checking purposes
   IF (idbg == 1) PRINT *, " Before entry to SUB outr"
   CALL outr

   !.....Read bathymetry file and setup logical mask arrays.....
   IF (idbg == 1) PRINT *, " Before entry to SUB bathy"
   CALL bathy

   !.....Define initial conditions.....
   IF (idbg == 1) PRINT *, " Before entry to SUB init"
   CALL init

   !.....Open files with lateral BC data and assign values at t=0.....
   IF (idbg == 1) PRINT *, " Before entry to SUB openbc0"
   CALL openbc0

   ! ... Open files with surface BC data and assign values at t=0.....
   IF (idbg == 1) PRINT *, " Before entry to SUB surfbc0"
   CALL surfbc0

   nn=1                               !MAC BEZNAR
   Ifile = "nbofilei0    "
       
   IF ( nn < 10 ) WRITE ( Ifile(10:12), FMT='(I1,"  ")' ) nn
   IF ( nn>= 10 ) WRITE ( Ifile( 9:12), FMT='(I2,"  ")' ) nn
       
   OPEN(unit=iboid,file=Ifile,FORM='FORMATTED',IOSTAT=ios)

   !i = 0
   !   DO l=1,cm1
   !  ! Ignore dry cells
   !  IF (.NOT. mask(l)  ) CYCLE
   !    i = i + 1
   !END DO
   ALLOCATE( out_array ( lm1 * km1 , 10 ), STAT=istat )

   !.....Output initial conditions.....
   IF (idbg == 1) PRINT *, " Before entry to SUB outt"
   IF (nnodes> 0) CALL outt(n,thrs)   ! Point depthiles
   IF (idbg == 1) PRINT *, " Before entry to SUB outv"
   IF (iox   > 0) CALL outv(n)   ! Cross-sections
   IF (idbg == 1) PRINT *, " Before entry to SUB outh"
   IF (iop   > 0) CALL outh(n)   ! Horiz-planes
   IF (idbg == 1) PRINT *, " Before entry to SUB outw"
   IF (iop   > 0) CALL outw(n)   ! Wind field
   IF (idbg == 1) PRINT *, " Before entry to SUB outz"
   IF (iotr  > 0) CALL outz(n)   ! Tracers in 3D domain
   IF (idbg == 1) PRINT *, " Before entry to SUB outs"
   IF (ipxml > 0) CALL outs(n,its)   ! Three-dimensional xml outputs
   IF (idbg == 1) PRINT *, " Before entry to SUB outp"
   IF (ipxml < 0) CALL outp(n)   ! Three-dimensional outputs for ptrack
   IF (idbg == 1) PRINT *, " Before entry to SUB outnb"
   IF (ioNBO > 0) CALL outNB(n,thrs)
   if (nopen > 0) CALL outcheckMass(n,thrs)
   IF (idbg == 1) PRINT *, " Before entry to SUB outscalarbalance"
   IF (iobal > 0) CALL OutScalarEnergyBalance(n,thrs)

   ltime1 = TIMER(0.0)                ! Retrieve begin time of loop


   CALL ConfigThreads(depth)

   nnH = 0
   nopth = 0
   p = 1
   nopenHH = 0
   iauxe(p) = lh_aux(1)
   iauxs(p) = 1 
   do p=1,num_threads
         do nn=1,nopen
               if(itype(nn) .NE. 3) THEN
               is = isbc(nn);
               ie = iebc(nn); 
               do i=is ,ie
                   if(i <= iauxe(p) .AND. i >= iauxs(p)) THEN
                        nnH(p,nn) = 1
                        nnHH(p,nn) = 1
                        nopth(nn) = nopth(nn)+1
                        EXIT
                   end if
               end do
               else
                   la=ij2l(isbc(nn),jsbc(nn))
                   if(la >= id_column(lhi(p)) .AND. la <= id_column(lhf(p))) THEN
                        nnH(p,nn) = 1
                        nopth(nn) = nopth(nn)+1
                        EXIT
                   end if
                   if(lhiW(p) .EQ. 0) THEN
                       laaux=lhi(p)
                   else
                       laaux=lhiW(p)
                   end if
                   if(la >= id_column(laaux) .AND. la <= id_column(lhf(p))) THEN
                        nnHH(p,nn) = 1
                   end if
               end if               
         end do
         if(p < num_threads)THEN
                iauxs(p+1) = iauxe(p)
                iauxe(p+1) = iauxe(p) + lh_aux(p+1)
                !iauxs(p+1) = iauxs(p) + lh_aux(p) -1  + p - 1
                !iauxe(p+1) = iauxe(p) + lh_aux(p+1)
         end if        
   end do

   iopssH = 0
   ioph2iop = 0
   do p=1,num_threads
         do nn=1,iopss
                   la=ij2l(ipss(nn),jpss(nn))
                   if(id_column(lhf(p)) .EQ. lm) THEN
                        aux_la = id_column(lhf(p)) + 1
                   ELSE
                        aux_la = id_column(lhf(p))
                   END IF
                   if(la >= id_column(lhi(p)) .AND. la <= aux_la) THEN
                        iopssH(p) = iopssH(p) + 1
                        ioph2iop(iopssH(p),p) = nn
                   end if
         end do       
   end do

   do p=1,num_threads
   		PRINT *,"lh:",lh(p),"lhaux:",lh_aux(p),"lhi:",id_column(lhi(p)),"lhf:",id_column(lhf(p))
   		PRINT *,"lhiE:",lhiE(p),"lhfE:",lhfE(p),"lhiW:",lhiW(p),"lhfW:",lhfW(p)
   		PRINT *,"depth:",depth,"ph:",ph(p)
                print *,"auxs:",iauxs(p),"auxe:",iauxe(p)
   end do
   print *,"im:",im,"jm:",jm,"lm:",lm

    !Start the parallel region

    !$omp parallel NUM_THREADS(num_threads) &
    !$omp default ( shared ) &
    !$omp private (n,lastiter,niter1,iter,istep,t_exmom2,t_matmom2,t_matcon2) &
    !$omp private (Bhaxpp,Bhaypp,Bth,Bth1,Bstart,Bend,ide_thread,i,lNCH,lSCH,lECH,lWCH) &
    !$omp private (Bex, Bth2,Beagx,Bearx,Bagx,Barx,Beagy,Beary,Bagy,Bary,Bsx,Bsy) &
    !$omp private (Bdd,Bqq,Brr,Bth3,Bth4,uhp2,uhp3,mincol) &
    !$omp private (ShearProduction, BuoyancyProduction, Dissipation, TKinE) &
    !$omp private (heatSourceB,uairB,vairB,cdwB,QswFrB,ish,ieh,noH,contBC) &
    !$omp private (Qsw2dB,Qlw2dB,Ta2dB,RH2dB,Cc2dB,uair2dB,vair2dB,bclncxB,hupdrhoB) &
    !$omp firstprivate(Qsw,Qn,Qlw,eta,Ta,Pa,RH,Cc,its,thrs)



    istep = 2
    ide_thread = omp_get_thread_num ( ) +1

    nopenH(ide_thread) = 0
    nopenHH(ide_thread) = 0
    isbcH(:,ide_thread) = 0
    iebcH(:,ide_thread) = 0
    isbcHH(:,ide_thread) = 0
    iebcHH(:,ide_thread) = 0
    jsbcH(:,ide_thread) = 0
    jebcH(:,ide_thread) = 0
    noh2no(:,ide_thread) = 0
    noh2noH(:,ide_thread) = 0
    eiptNBI(:,ide_thread) = 0
    siptNBI(:,ide_thread) = 0
    eiptNBIH(:,ide_thread) = 0
    siptNBIH(:,ide_thread) = 0
    !xxNBO = 5 ! FJR
    contBC = 0
    print *,"nopen:",nopen
    do noH=1,nopen
          if(nnH(ide_thread,noH) .EQ. 1)THEN
 		nopenH(ide_thread) = nopenH(ide_thread) + 1
                noh2no(nopenH(ide_thread),ide_thread) = noH
                siptNBI(noh,ide_thread)=1
                eiptNBI(noh,ide_thread)=iptNBI(noH)
                do isH=isbc(noH),iebc(noh)
                      if(ish >= iauxs(ide_thread) .AND. ish <= iauxe(ide_thread))THEN
                             isbcH(noH,ide_thread) = isH
                             EXIT
                      end if
                      siptNBI(noh,ide_thread)=siptNBI(noh,ide_thread)+((kmz(ij2l(ish,jsbc(noH)))-k1+1))
                end do
                do ieH=iebc(noH),isbc(noh),-1
                      if(ieh >= iauxs(ide_thread) .AND. ieh <= iauxe(ide_thread))THEN
                             iebcH(noH,ide_thread) = ieH
                             EXIT
                      end if
                      
                      eiptNBI(noh,ide_thread)=eiptNBI(noh,ide_thread)-((kmz(ij2l(ieh,jsbc(noH)))-k1+1))
 
                end do
                jsbcH(noH,ide_thread) = jsbc(noH)
                jebcH(noH,ide_thread) = jebc(noH)
          end if
          if(nnHH(ide_thread,noH) .EQ. 1)THEN
 		nopenHH(ide_thread) = nopenHH(ide_thread) + 1
                noh2noH(nopenHH(ide_thread),ide_thread) = noH                
          end if  
    end do
    !$omp barrier
    if(ide_thread > 1) iauxs(ide_thread) = iauxs(ide_thread) + 1
    do noH=1,nopen
          if(nnH(ide_thread,noH) .EQ. 1)THEN
                contBC = contBC + 1
                siptNBIH(noH,ide_thread)=1
                eiptNBIH(noH,ide_thread)=iptNBI(noH)
                do isH=isbc(noH),iebc(noH)
                      if(ish >= iauxs(ide_thread) .AND. ish <= iauxe(ide_thread))THEN
                             isbcHH(noH,ide_thread) = isH
                             print *,"salgo inicio",ish,"iauxs:",iauxs(ide_thread),"iauxe:",iauxe(ide_thread),noh
                             EXIT
                      end if
                      siptNBIH(noH,ide_thread)=siptNBIH(noH,ide_thread)+((kmz(ij2l(ish,jsbc(noH)))-k1+1))
                      print *,"sumo inicio",siptNBIH(noH,ide_thread),noh,"h:",ide_thread
                end do
                do ieH=iebc(noH),isbc(noH),-1
                      if(ieh >= iauxs(ide_thread) .AND. ieh <= iauxe(ide_thread))THEN
                             iebcHH(noH,ide_thread) = ieH
                             EXIT
                             print *,"salgo final",ieh,"iauxs:",iauxs(ide_thread),"iauxe:",iauxe(ide_thread),noh
                      end if
                      eiptNBIH(noH,ide_thread)=eiptNBIH(noH,ide_thread)-((kmz(ij2l(ieh,jsbc(noH)))-k1+1))
                      print *,"sumo final",eiptNBIH(noH,ide_thread),noh,"h:",ide_thread
                end do
          end if 
    end do

        print *,"nopennn:",nopenH,":",ide_thread
        do n=1,nopenH(ide_thread)
            j=noh2no(n,ide_thread)
            print *,"h:",ide_thread,"nopen:",nopenH(ide_thread),"isbc:",isbcH(j,ide_thread)
            print *,"iebc:",iebcH(j,ide_thread),"jsbc:",jsbcH(j,ide_thread),"jebc:",jebcH(j,ide_thread)
            print *,"iebcH:",iebcHH(j,ide_thread),"isbcH:",isbcHH(j,ide_thread)
            print *,"noh2no:",noh2no(j,ide_thread),"siptNBI:",siptNBI(j,ide_thread),"eiptNBI:",eiptNBI(j,ide_thread)
            print *,"noh2no:",noh2no(j,ide_thread),"siptNBIH:",siptNBIH(j,ide_thread),"eiptNBIH:",eiptNBIH(j,ide_thread)
            print *,"iptNBI:",iptNBI(j),"nopenHH:",nopenHH(ide_thread)
        end do

        if(iopss > 0) THEN
        do n=1,iopssH(ide_thread)
            j=ioph2iop(n,ide_thread)
            print *,"h:",ide_thread,"iopss:",iopssH(ide_thread)
            print *,"ipss:",ipss(j),"jpss:",jpss(j)
        end do
        end if


    CALL BorderThreads(Bstart,Bend,ide_thread,mincol)

    ! Allocate memory to private variables
    ALLOCATE (Bhaxpp(1:km1,Bstart:Bend+1),lSCH(Bstart:Bend),lNCH(Bstart:Bend), &
    & lWCH(Bstart:Bend), lECH(Bstart:Bend), Bhaypp(1:km1,Bstart:Bend+1), &
    & Bth(1:km1,Bstart:Bend+1), Bth1(1:km1,Bstart:Bend+1), Bex(1:km1,Bstart:Bend+1), &
    & Bth2(1:km1,Bstart:Bend+1), Bagx(1:km1,Bstart:Bend+1), Barx(1:km1,Bstart:Bend+1), &
    & Bearx(Bstart:Bend+1), Beagx(Bstart:Bend+1), Bagy(1:km1,Bstart:Bend+1), &
    & Bary(1:km1,Bstart:Bend+1), Beagy(Bstart:Bend+1), Beary(Bstart:Bend+1), &
    & Bsx(Bstart:Bend+1),Bsy(Bstart:Bend+1),Bdd(Bstart:Bend+1),Bqq(Bstart:Bend+1), &
    & Brr(Bstart:Bend+1), Bth3(1:km1,Bstart:Bend+1), Bth4(1:km1,Bstart:Bend+1),&
    & uhp2(1:km1,1:lm1),uhp3(1:km1,1:lm1),heatSourceB(1:km1,Bstart:Bend+1),&
    & uairB(Bstart:Bend+1),vairB(Bstart:Bend+1),cdwB(Bstart:Bend+1),QswFrB(1:km1,Bstart:Bend+1), &
    & Qsw2dB(nmetstat),Qlw2dB(nmetstat),Ta2dB(nmetstat),RH2dB(nmetstat),Cc2dB(nmetstat), &
    & uair2dB(nmetstat),vair2dB(nmetstat),bclncxB(1:km1),hupdrhoB(1:km1))

    CALL InitThreads(t_exmom2,t_matmom2,t_matcon2,Bhaxpp,Bhaypp,Bth,Bth1,Bstart, &
    & Bend,ide_thread,lNCH,lSCH,lECH,lWCH,Bex, Bth2,Beagx,Bearx,Bagx,Barx,Beagy, &
    & Beary,Bagy,Bary,Bsx,Bsy,Bdd,Bqq,Brr,Bth3,Bth4,uhp2,uhp3,mincol,ShearProduction, &
    & BuoyancyProduction, Dissipation, TKinE,heatSourceB,uairB,vairB,cdwB,eta, &
    & QswFrB,Qsw2dB,Qlw2dB,Ta2dB,RH2dB,Cc2dB,uair2dB,vair2dB,Qsw,Qn,Qlw,Ta,Pa,RH,Cc)

    bclncxB = 0.0;
    hupdrhoB = 0.0;


   !.....Integrate over time.....
   IF (idbg == 1) PRINT *, "Before entering loop over time"
   DO n = 1, nts
      
      its = its + idt
      thrs = its/3600.
      IF(omp_get_thread_num ( )==0)THEN
      TimeStart = TIMER(0.0)
      CALL compute_date (idt)
      END IF
      
      
      niter1 = niter
      lastiter = 0
      IF ( n == 1 ) THEN
         ! Use at least two iterations to start computations
         IF ( niter <  2 ) niter1 = 2
         IF ( niter >= 2 ) niter1 = niter
         GO TO 1
      END IF

      !.....Solve leapfrog step.....
      istep = 1
      IF (idbg == 1) PRINT *, " Before entry into SUB fd in leapfrog step"
      IF((itrap == 0) .OR. (MOD(n,MAX(itrap,1)) /= 0)) lastiter = 1


!      print *,"entro en fdfrog",omp_get_thread_num ( )
      CALL fd(n,t_exmom2,t_matmom2,t_matcon2,Bhaxpp,Bhaypp,Bth,Bth1,Bstart,Bend, &
      & lSCH,lNCH,lWCH,lECH,Bex,Bth2,Beagx,Bearx,Bagx,Barx,Beagy,Beary,Bagy,Bary, &
      & Bsx,Bsy,Bdd,Bqq,Brr,Bth3,Bth4,istep,lastiter, &
      & ShearProduction,BuoyancyProduction, Dissipation, TKinE, &
      & Qsw,Qn,Qlw,eta,Ta,Pa,RH,Cc,Qsw2dB,Qlw2dB,Ta2dB,RH2dB,Cc2dB,uair2dB,vair2dB, &
      & uairB,vairB,cdwB,heatSourceB,QswFrB,iter,bclncxB,hupdrhoB,its,thrs)

      !$omp barrier

      IF((itrap == 0) .OR. (MOD(n,MAX(itrap,1)) /= 0)) GO TO 2

      !.....Solve trapezoidal step.....
      CALL settrap
      
      !$omp barrier
    1 istep = 2
      DO iter = 1, niter1
         IF ( iter == niter1 ) lastiter = 1;
         IF (idbg == 1) PRINT *, " Before entry into SUB fd in trap step"

         CALL fd(n,t_exmom2,t_matmom2,t_matcon2,Bhaxpp,Bhaypp,Bth,Bth1,Bstart,Bend, &
         & lSCH,lNCH,lWCH,lECH,Bex,Bth2,Beagx,Bearx,Bagx,Barx,Beagy,Beary,Bagy,Bary, &
         & Bsx,Bsy,Bdd,Bqq,Brr,Bth3,Bth4,istep,lastiter, &
         & ShearProduction,BuoyancyProduction, Dissipation, TKinE, &
         & Qsw,Qn,Qlw,eta,Ta,Pa,RH,Cc,Qsw2dB,Qlw2dB,Ta2dB,RH2dB,Cc2dB,uair2dB,vair2dB, &
         & uairB,vairB,cdwB,heatSourceB,QswFrB,iter,bclncxB,hupdrhoB,its,thrs)

         !$omp barrier
         IF (idbg == 1) PRINT *, " After exit from SUB fd in trapezoidal step"
         
         IF(iter < niter1) CALL settrap2
         
         !$omp barrier
         IF(iter > 20) THEN
            PRINT *, " ERROR--Too many iterations requested"
            EXIT
         END IF
      END DO
      !$omp barrier
      !.....Save information for next time step.....
    2 CALL save(istep)

      !$omp barrier
      IF(omp_get_thread_num ( )==0)THEN

      IF(itype(1) .EQ. 4) THEN
      ! ... Water Surface Elevation --> mass                                !ACC BEZNAR
         dzino = 0.0E0
         DO j = j1,jm;
           DO i = i1, im
             IF (.NOT. mask2d(i,j)) CYCLE
             ! ... Get l- from (i,j)
             lcon   = ij2l(i,j)
             ! ... Add water surface displacements
             dzino = dzino + s(lcon)
             print *, "lcon", lcon, "s", s(lcon)
           ENDDO
         ENDDO

         uflow = 0.0E0

         ifrontera  = isbc(1)-5; 
         js = jsbc(1)  ; 
         je = jebc(1)  ;
         DO j = js, je
            IF (.NOT. mask2d(ifrontera,j)) CYCLE
            lcon2 = ij2l(ifrontera,j)
            DO k= k1, kmz(lcon2)
                uflow = uflow + uh(k,lcon2)
            PRINT *, "lcon2", lcon2, "uh", uh(k,lcon2),"i",ifrontera,"j",j,"k",k
            END DO
         END DO

       WRITE (UNIT=iboid, FMT='(5E20.11)') thrs, uflow, dzino*dx*dy         !ACC BEZNAR
      END IF

      !.....Output results.....
    3 IF((nnodes > 0) .AND. (MOD(n,MAX(ipt,  1)) == 0)) CALL outt(n,thrs)
      IF((iox    > 0) .AND. (MOD(n,MAX(iox,  1)) == 0)) CALL outv(n)
      IF((iop    > 0) .AND. (MOD(n,MAX(iop,  1)) == 0)) CALL outh(n)
      IF((iop    > 0) .AND. (MOD(n,MAX(iop,  1)) == 0)) CALL outw(n)
      IF((iotr   > 0) .AND. (MOD(n,MAX(iotr, 1)) == 0)) CALL outz(n)
      IF((ipxml  > 0) .AND. (MOD(n,MAX(ipxml,1)) == 0)) CALL outs(n,its)
      IF((ipxml  < 0) .AND. (MOD(n,MAX(apxml,1)) .GE. (apxml-100))) CALL AverageOutp(n)  ! Andrea Ptrack
      IF((ipxml  < 0) .AND. (MOD(n,MAX(apxml,1)) == 0)) CALL outp(n)
      !IF((ioNBO  > 0) .AND. (MOD(n,MAX(ioNBO,1)) == 0)) CALL outNB(n,thrs)
      IF((ioNBO  > 0) ) CALL outNB(n,thrs)
      if((nopen > 0)  )  CALL outcheckMass(n,thrs)
      IF(iobal   > 0  .AND. (MOD(n,MAX(ipt  ,1)) == 0)) CALL OutScalarEnergyBalance(n,thrs)

      !.....Write to log file and check if job should be stopped.....
      CALL outlog(thrs)
      CALL check_stopfile(n)

      !.....End loop over time.....
      IF(n > maxcount) THEN
         PRINT *, " ERROR--A maximum of 1 million time steps is allowed"
         EXIT
      END IF
      TimeEnd = TIMER(0.0)

      PRINT *, 'Time spent in step ', n, ' = ', TimeEnd - TimeStart, ' seconds'
      END IF     

   END DO
   CALL cputimes(t_exmom2,t_matmom2,t_matcon2)
   !$omp end parallel
   ltime2 = TIMER(0.0)
   PRINT '(A,F10.3,A)', " Loop time =", ltime2-ltime1, " seconds"

   !.....Output CPU times.....

   utime = TIMER(0.0)

   !.....Print total time of run.....
   PRINT '(A, F10.3, A)',  &
   & " Total CPU Time =", utime-stime,         " seconds"

   !.....End program.....
   PRINT '(A)', " "
   PRINT *, " *****PROGRAM TERMINATED NORMALLY"

END PROGRAM si3d


