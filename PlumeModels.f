	SUBROUTINE RWPS(YEAR,JULDAY,WSEL,RWPSEL,LAYERS,KRWPS,HCELL,
     +LOCAT,TE,CO2M,ELEVT,QWT,TPLUMET,COMGPT,QWD,KTOP,TPLUME0,COMGP0)
C	THIS PROGRAM IS WRITTEN TO PREDICT ENTRAINMENT AND DETRAINMENT CHARACTERISTICS 
C     OF THE PUMPED INFLOW IN SPRING HOLLOW RESERVOIR, VA.    
C     THE MODEL IS MODIFIED FROM THE WUEST ET AL. (1992) CIRCULAR BUBBLE PLUME MODEL.
C	
C	VERSION 3 (to calculate Spring Hollow RWPS plume variables for manual coupling with W2 v3.2)
C     by: Vickie Singleton
C
C     This version includes the following:
C     1.  Revised dissolved flux equations to use ambient concentrations for entrainment
C     2.  Correction of ambient salinity interpolation from input data file
C     3.  Correction for salinity units in salinity flux equations (2-2-05)  
C     4.  Determination and use of average ambient density in pressure calculations (2-2-05)   
C     5.  Revision of numerical integration method from Euler to fourth-order Runge-Kutta (2-23-05)
C     6.  Corrected equations that include salinity to make units consistent (11-7-06)
C     7.  Revised entrainment and momentum equations for single-phase round jet/plume (4-23-07)
C     8.  Revised alpha from 0.11 to 0.08 (value for single-phase plumes per Fischer et al., 1979)
C     9.  Calculated local entrainment coefficient per Fischer et al. (1979) as opposed to using constant value (4-30-07)
C    10.  Revised so that initial plume temperature is either interpolated from raw data or calculated using 
C         fitted sine curves (4-15-08)
C    11.  Revised diurnal sine curve equation for July 2003 initial plume temp based on RWPS data (5-13-08)
C    12.  Output RWPS temperature and DO to si3D_procedures to calculate temperature and DO of detrainment for 
C         mass conservative approach (5-15-09)           
C
C	May 15, 2009
C
C	VARIABLES
C
C	ALPHA=ENTRAINMENT COEFFICIENT (-)
C     ALPHAJ=JET ENTRAINMENT COEFFICIENT (-)
C     ALPHAP=PLUME ENTRAINMENT COEFFICIENT (-) 
C     B=PUMP DISCHARGE PIPE RADIUS (m)
C	CO2=DISSOLVED OXYGEN CONCENTRATION (mol/m3)
C	CN2=DISSOLVED NITROGEN CONCENTRATION (mol/m3)
C	COMG=DISSOLVED OXYGEN CONCENTRATION (g/m3)
C     COMGPT=DO CONCENTRATION OF DETRAINMENT AT TOP OF PUMPED INFLOW PLUME (g/m3) 
C	CNMG=DISSOLVED NITROGEN CONCENTRATION (g/m3)
C     DENSE20=DENSITY OF WATER AT 20 C
C	DENSEA=AMBIENT WATER DENSITY (kg/m3)
C	DENSEW=WATER DENSITY IN PLUME (kg/m3)
C     DNAMB=AMBIENT DISSOLVED NITROGEN CONCENTRATION (g/m3)
C     DOAMB=AMBIENT DISSOLVED OXYGEN CONCENTRATION (g/m3)
C	E=ENTRAINMENT FACTOR (m3/s)
C     ELEV=ELEVATION (m)
C     ELEVT=TERMINAL ELEVATION OF PUMPED INFLOW PLUME (m)
C	FDO=DISSOLVED OXYGEN FLUX (mol/s)
C	FDN=DISSOLVED NITROGEN FLUX (mol/s)
C	FSAL=SALINITY FLUX (kg/s)
C	FTEMP=TEMPERATURE FLUX (C m3/s)
C     FRCNOT=FRACTION OF NITROGEN IN ATMOSPHERE (-)
C     GAMMA=SALINITY CONVERSION FACTOR [kg/m3/(uS/cm)]
C     HCELL=HEIGHT OF CELL IN GRID (m)
C     JULDAY=JULIAN DAY FOR EACH YEAR (1-365)
C     KRWPS=Si3D GRID LAYER CORRESPONDING TO DISCHARGE DEPTH FOR RAW WATER PUMPING STATION (-)
C     KTOP=Si3D GRID LAYER CORRESPONDING TO TOP OF PUMPED INFLOW PLUME (-)
C	MOMENT=MOMENTUM (m4/s)
C     PATM=ATMOSPHERIC PRESSURE AT AVERAGE WSEL (Pa)
C	QW=FLOWRATE OF WATER (m3/s)
C     QW0=WATER FLOW RATE OF PUMPED INFLOW, EQUAL TO INITIAL PLUME WATER FLOW RATE (m3/s)
C     QWD=WATER FLOW RATE OF WITHDRAWAL ENTRAINMENT INTO PLUME (m3/s)    
C     QWT=TOTAL DETRAINMENT FLOW RATE AT TOP OF PUMPED INFLOW PLUME (m3/s)  
C     RP=PLUME RICHARDSON NUMBER (-)
C     RWPSEL=ELEVATION OF RAW WATER PUMPING STATION DISCHARGE IN RESERVOIR (m)
C     SALAMB=SALINITY OF AMBIENT WATER (uS/cm)
C	SALPLU=SALINITY OF THE PLUME (uS/cm)
C	TPLUME=PLUME TEMPERATURE (C)
C	TAMB=AMBIENT WATER TEMPERATURE (C)
C     TAVG=AVERAGE AMBIENT WATER TEMPERATURE (C)
C     TDS=TOTAL DISSOLVED SOLIDS (g/m3) [TDS=SALP*0.64 conversion eqn. from Chapra book]
C	V=WATER VELOCITY (m/s)
C     WSEL=WATER SURFACE ELEVATION (m)
C	Z=DEPTH TO DIFFUSER (m)
C
      REAL*8 ALPHA,AREA,B,CO2,COMG,CN2,CNMG,DS,DZ,DMOM,DFTEM,
     +DFDO,DFDN,DFSAL,DENSEA,DENSEW,E,FDO,FDN,FSAL,FTEMP,G,GAMMA,
     +MOMENT,PI,PZ,QW,SAL0,SALAMB,SALPLU,TAMB,TPLUME,V,VG,Z,
     +LAYER,ELEV,DT,TE(1000),XLOC,DEPTH,CO2M(1000),LOCAT(1000),
     +DOAMB,COMGP,CNMGP,WSEL,PATM,RWPSEL,DNAMB,TAVG,SUMTEMP,
     +SUMSAL,H,DYDX(6),Y(6),YOUT(6),DENSE20,
     +BUOY,R,RP,ALPHAJ,ALPHAP,ELEVT,QWT,TPLUMET,COMGPT,JULDAY,QW0,HCELL,
     +QWD(1000),HWITH,SUMQWD,TPLUME0,COMGP0,X,JD2003(5000),RWPST(5000),
     +MINUTE,JD2004(5000),HN2
      INTEGER II,IJ,IK,IN,KK,NEQN,IM,MI,JK,JL,LAYERS,KS,KRWPS,KTOP,JJ,
     +YEAR,NUMT
C           
!.....Open RWPS temperature file for August 2003 and October 2004.............
      IF(YEAR.EQ.2003.AND.JULDAY.GE.226.75.AND.JULDAY.LT.239.5)THEN
          OPEN(5,FILE='Aug03RWPS.DAT',STATUS='UNKNOWN')
          NUMT=0
   41     READ(5,*,END=40)JD2003(NUMT+1),RWPST(NUMT+1)
             NUMT=NUMT+1
             GOTO 41
   40     CLOSE(UNIT=5)          
      ELSEIF(YEAR.EQ.2004)THEN
          OPEN(6,FILE='Oct04RWPS.DAT',STATUS='UNKNOWN')
          NUMT=0
   43     READ(6,*,END=42)JD2004(NUMT+1),RWPST(NUMT+1)
             NUMT=NUMT+1
             GOTO 43
   42     CLOSE(UNIT=6)                
      ENDIF
      SUMTEMP=0.0
      DO 50 KS=1,LAYERS
          SUMTEMP=SUMTEMP+TE(KS)
   50 CONTINUE        
      TAVG=SUMTEMP/LAYERS                             
      DEPTH=WSEL-RWPSEL
      Z=DEPTH
      ELEV=RWPSEL
      X=0.
C
C     Interpolate temperature and DO boundary condition profiles to obtain jet/plume initial conditions
      XLOC=DEPTH-X
      CALL LININT(LOCAT,TE,LAYERS,XLOC,TAMB)
      CALL LININT(LOCAT,CO2M,LAYERS,XLOC,COMG)
!     Initial plume temperature based on average measured RWPS wetwell temperature
      PI=ACOS(-1.0)
      IF(YEAR.EQ.2003)THEN
!          IF(JULDAY.LE.196.)TPLUME0=21.0
          IF(JULDAY.LE.196.)THEN
              MINUTE=(JULDAY-AINT(JULDAY))*24.*60.
              TPLUME0=2.0*SIN(2.*PI/1440.*MINUTE-1.57)+17.
!          IF(JULDAY.GE.226)TPLUME0=21.5
          ELSEIF(JULDAY.GE.226.0.AND.JULDAY.LT.232.365)THEN    
              MINUTE=((JULDAY-AINT(JULDAY))-0.594)*24.*60.
              TPLUME0=1.035*SIN(2.*PI/1359.*MINUTE+0.068)+21.045
          ELSEIF(JULDAY.GE.232.365)THEN
              CALL LININT(JD2003,RWPST,NUMT,JULDAY,TPLUME0)
          ENDIF       
      ELSEIF(YEAR.EQ.2004.AND.JULDAY.GE.296.)THEN        
!          TPLUME0=16.25
              CALL LININT(JD2004,RWPST,NUMT,JULDAY,TPLUME0)
      ENDIF
      TPLUME=TPLUME0                  
!     Initial RWPS plume DO based on saturation value at measured temperature in RWPS wetwell
      COMGP0=0.0035*TPLUME**2-0.3368*TPLUME+14.406
      COMGP=COMGP0
      DOAMB=COMG      
      CO2=COMGP/32.
!     Assume constant value for ambient salinity
      IF(YEAR.EQ.2003)THEN
          IF(JULDAY.LE.196.)SALAMB=193.
          IF(JULDAY.GE.226.)SALAMB=201.
      ELSEIF(YEAR.EQ.2004)THEN    
          IF(JULDAY.GE.296.)SALAMB=218.
      ENDIF    
!     Initial salinity of RWPS plume equal to TDS value measured in river by USGS gauging station           
      IF(YEAR.EQ.2003)THEN
          IF(JULDAY.LE.196.)SALPLU=164./0.64
          IF(JULDAY.GE.226.)SALPLU=202./0.64
      ELSEIF(YEAR.EQ.2004)THEN    
          IF(JULDAY.GE.296.)SALPLU=165./0.64
      ENDIF              
C
C     CONSTANTS
!     Entrainment coefficient assumed to be equal to value for one-dimensional plumes 
!     given by Fischer et al., p. 371.
      ALPHA=0.083
      G=9.80665
      GAMMA=6.9E-4
      LAYER=0.0
      DENSE20=998.2
      FRCNATM=0.79
      PATM=96261.0
      RP=0.557
      ALPHAJ=0.0535
      ALPHAP=0.0833    
      B=48./2./12./3.281     
C      
C     AMBIENT AND AVERAGE WATER DENSITIES
C
      DENSEA=(0.059385*TAMB**3-8.56272*TAMB**2+65.4891*TAMB)*0.001
     ++999.84298+(GAMMA)*SALAMB
      DENSEW=DENSEA
C
C     Assume initial ambient dissolved nitrogen conc. equals saturated conc. at surface.
      HN2=(1.042-0.02457*TPLUME+3.1714E-4*TPLUME**2)/100000.
      CN2=(PATM*FRCNATM)*HN2      
      CNMG=CN2*28.0
      CNMGP=CNMG
      DNAMB=CNMG
C
!     INITIAL WATER VELOCITY CALCULATED FROM PUMPED WATER FLOW RATE
      IF(YEAR.EQ.2003)THEN
          IF(JULDAY.LE.182.306)THEN
              QW0=0.657
          ELSEIF(JULDAY.GE.190.354.AND.JULDAY.LE.190.583)THEN
              QW0=0.657
          ELSEIF(JULDAY.GE.191.323.AND.JULDAY.LE.196.)THEN
              QW0=0.657
          ELSEIF(JULDAY.GE.226.0.AND.JULDAY.LE.226.999)THEN
              QW0=0.657
          ELSEIF(JULDAY.GE.227.0.AND.JULDAY.LE.228.292)THEN
              QW0=1.313
          ELSEIF(JULDAY.GE.231.396.AND.JULDAY.LE.231.433)THEN
              QW0=0.657
          ELSEIF(JULDAY.GE.231.434.AND.JULDAY.LE.232.465)THEN
              QW0=1.313
          ELSEIF(JULDAY.GE.232.466.AND.JULDAY.LE.232.475)THEN
              QW0=0.657
          ELSEIF(JULDAY.GE.232.476.AND.JULDAY.LE.238.556)THEN
              QW0=1.313
          ELSE
              QW0=0.0
              GOTO 100
          ENDIF        
      ELSEIF(YEAR.EQ.2004)THEN
          IF(JULDAY.GE.296.0.AND.JULDAY.LE.309.896)THEN
              QW0=0.657
          ELSEIF(JULDAY.GE.311.681.AND.JULDAY.LE.317.740)THEN
              QW0=0.657
          ELSE
              QW0=0.0
              GOTO 100
          ENDIF    
      ELSE
          QW0=0.0
          GOTO 100
      ENDIF       
      QW=QW0
      V=QW/(PI*B**2)
!
!     CALCULATION OF INITIAL LOCAL ENTRAINMENT COEFFICIENT
!      QW=PI*B**2*V
      MOMENT=PI*B**2*V**2
!      BUOY=G*ABS(DENSEA-DENSEW)/DENSEW*QW
!      R=QW*BUOY**0.5/MOMENT**(5./4.)
!      ALPHA=ALPHAJ*DEXP(DLOG(ALPHAP/ALPHAJ)*(R/RP)**2)
C
C     VARIABLE TRANSFORMATION      
!      E=2.*(L+2.*B)*ALPHA*V   
!      MOMENT=2.*L*B*V**2
      E=2*PI*B*ALPHA*V 
      FTEMP=QW*TPLUME 
      FSAL=QW*(SALPLU*GAMMA/DENSE20)*DENSEW
C     Previous equation corrected to account for salinity units conversion.        
      FDO=QW*CO2
      FDN=QW*CN2    
      PZ=PATM+(DENSEA*G*Z)
!     Initialize lateral withdrawal/entrainment flowrate for first/lowest cell.
      JJ=0
      QWD(KRWPS)=QW            
C    
C	SOLUTION PROCEEDURE
C
      DZ=0.001
      H=0.001
      HWITH=0.0
 10   Z=Z-DZ
      X=X+DZ
      LAYER=LAYER+DZ
      ELEV=ELEV+DZ
C	
C     Interpolate temperature and DO profile input to obtain jet/plume boundary conditions
      XLOC=DEPTH-X
      CALL LININT(LOCAT,CO2M,LAYERS,XLOC,COMG)
      DOAMB=COMG
      CO2=COMG/32.
      CALL LININT(LOCAT,TE,LAYERS,XLOC,TAMB)
C                  
C     Use subroutines for Runge Kutta method solution
      NEQN=6
      Y(1)=QW
      Y(2)=MOMENT
      Y(3)=FTEMP
      Y(4)=FSAL
      Y(5)=FDO
      Y(6)=FDN
      CALL DERIVS_2(E,DENSEA,DENSEW,G,B,TAMB,SALAMB,GAMMA,DENSE20,
     +     DOAMB,PI,V,COMGP,DNAMB,CNMGP,Z,Y,DYDX)
      CALL RK4_2(E,DENSEA,DENSEW,G,B,TAMB,SALAMB,GAMMA,DENSE20,DOAMB,
     +     PI,V,COMGP,DNAMB,CNMGP,Y,DYDX,NEQN,Z,H,YOUT)
C       
      QW=YOUT(1)
      MOMENT=YOUT(2)
      FTEMP=YOUT(3)
      FSAL=YOUT(4)
      FDO=YOUT(5)
      FDN=YOUT(6)
      IF(MOMENT.LT.0.0)THEN
         TPLUME=FTEMP/QW
	   SALPLU=FSAL/(QW*DENSEW)/(GAMMA/DENSE20)
C    Previous equation corrected to consistently express salinity in uS/cm	   
	   CO2=FDO/QW
	   CN2=FDN/QW
	   GOTO 20
      ENDIF
      V=MOMENT/QW
      AREA=QW/V
      B=(AREA/PI)**0.5
      DENSEA=(0.059385*TAMB**3-8.56272*TAMB**2+65.4891*TAMB)*0.001
     ++999.84298+(GAMMA)*SALAMB
      DENSEW=(0.059385*TPLUME**3-8.56272*TPLUME**2+65.4891*TPLUME)*0.001
     ++999.84298+(GAMMA)*SALPLU
C     Previous equation re-revised to account for correct salinity units (uS/cm) in density calculations.      
!      
!     CALCULATION OF LOCAL ENTRAINMENT COEFFICIENT
!      BUOY=G*ABS(DENSEA-DENSEW)/DENSEW*QW
!      R=QW*BUOY**(1./2.)/MOMENT**(5./4.)
!      ALPHA=ALPHAJ*DEXP(DLOG(ALPHAP/ALPHAJ)*(R/RP)**2)
!      IF(ALPHA.LT.0.05351)THEN
!          WRITE(*,*)BUOY,MOMENT,R,ALPHA
!          WRITE(*,*)ELEV,QW,TPLUME,COMGP
!          PAUSE
!      ENDIF   
!                   
      E=2*PI*B*ALPHA*V
!     Add incremental entrainment to total cell entrainment/withdrawal       
      QWD(KRWPS-JJ)=QWD(KRWPS-JJ)+E*DZ
      HWITH=HWITH+DZ   
      IF(HWITH.GT.HCELL)THEN
          JJ=JJ+1
          HWITH=0.0
          QWD(KRWPS-JJ)=0.0
      ENDIF                
      TPLUME=FTEMP/QW
 	SALPLU=FSAL/(QW*DENSEW)/(GAMMA/DENSE20)
C     Previous equation corrected to consistently express salinity in uS/cm
      CO2=FDO/QW
      CN2=FDN/QW
      COMGP=CO2*32.
      CNMGP=CN2*28.
      PZ=PATM+(DENSEA*G*Z)
      IF(V.GT.1.E-6)THEN
!     Revised from "Z.GT.0.0" to "XLOC.GT.0.0" by VLS on 10-3-07       
           IF(XLOC.GT.0.001)THEN
	        GOTO 10
	     ENDIF
      ENDIF
!      
   20 ELEVT=ELEV
      SUMQWD=0.0
      DO 30 I=KRWPS-JJ,KRWPS
          SUMQWD=QWD(I)+SUMQWD
   30 CONTINUE     
      QWT=SUMQWD
!     VLS: Revised QWT to represent pumped inflow flow rate to match coding in si3d_procedures.f90 (5-13-09)
      QWT=QW0            
      TPLUMET=TPLUME
      COMGPT=COMGP
      KTOP=KRWPS-JJ
!     Subtract RWPS inflow flow rate from entrainment/withdrawal flow rate in first/lowest cell.
      QWD(KRWPS)=QWD(KRWPS)-QW0
  100 IF(QW0.EQ.0.)THEN
          ELEVT=ELEV
          QWT=QW0
          TPLUMET=TPLUME0
          COMGPT=COMGP0
          KTOP=KRWPS
          QWD(KRWPS)=0.0    
      ENDIF
      RETURN
      END

C     *******************************************************************************
	SUBROUTINE LINEPLU_v0(YEAR,JULDAY,WSEL,DIFFEL,LAYERS,LDIFF,LAMBNOT,
     +SALAMB,PATM,DIAMM,LAKE,QSCFM,FRCONOT,LOCAT,TE,CO2M,ELEVT,QWT,
     +TPLUMET,COMGPT)
C     *******************************************************************************
C	THIS SUBROUTINE IS WRITTEN TO PREDICT THE PERFORMANCE OF A LINEAR BUBBLE PLUME.  
C     THE MODEL IS BASED ON THE WUEST ET AL. (1992) CIRCULAR BUBBLE PLUME MODEL.
C     By: Vickie Singleton
C	
C	VERSION 2 (to couple with Francisco Rueda's reservoir model, constant reservoir depth)
C
C     This version includes the following:
C     1.  Revised momentum flux equations re-derived on October 4, 2004
C     2.  Revised gaseous flux equations with correct plume area 
C     3.  Revised dissolved flux equations to use ambient concentrations for entrainment
C     4.  Revised dissolved flux equations to use plume concentrations for gas transfer 
C     5.  Correction of ambient salinity interpolation from input data file
C     6.  Corrected initial bubble size correlation equation for gas flow rate per unit length (11-17-04)
C     7.  Correction for salinity units in salinity flux equations (2-2-05)  
C     8.  Revised gas holdup equation in loop to account for correct plume cross-sectional area occupied by bubbles (2-2-05)
C     9.  Determination and use of average ambient density in pressure calculations (2-2-05)   
C     10.  Correction of Bnot/diffuser source radius calculation per Wuest et al. 1992, Figure 2. (2-2-05)
C     11.  Revision of numerical integration method from Euler to fourth-order Runge-Kutta (2-23-05)
C     12.  Revision of interpolation of boundary profiles from 1 m increments to 0.1 m increments (3-2-05)
C     13.  Added calculation of initial water velocity using initial Froude number of 1.6.  (Previously, initial
C          water velocity was assumed to be 0.07 m/s per Dan's original program.) (6-22-05)
C     14.  Revised entrainment and spreading coefficients from 0.08 and 0.85, respectively, to 0.11 and 0.93, 
C          respectively, to account for top-hat profiles versus Gaussian profiles (9-19-05)
C     15.  Revised initial Froude number from 1.6 to 2.0 (refer to calculations). (10-1-05)
C     16.  Revised characteristic length in Froude number calculation from equivalent radius to initial plume width. (10-1-05) 
C     17.  Corrected equations that include salinity to make units consistent (11-7-06)
C     18.  Revision of depth increments for input files from 0.1 m to 1 m to more closely match W2 output (11-14-06) 
C     19.  Enabled program to automatically update W2 input files after each line plume model run (1-5-07)
C     20.  Enabled program to read and interpolate AGPM output files directly as line plume boundary conditions (1-23-07)
C     21.  Modified to couple with Sep 1998 data set (3-5-07)
C     22.  Revised elevation of diffuser in Segment 17 for revised bathymetry (3-22-07)
C     23.  Passed LAMBNOT, SALAMB, PATM, DIAMM, and LAKE as arguments to accomodate Amisk Lake (4-19-09)           
C
C	April 14, 2009
C
C	VARIABLES
C
C	ALPHA=ENTRAINMENT COEFFICIENT (-)
C     B=1/2 DIFFUSER WIDTH (m)
C     BAVG=AVERAGE 1/2 DIFFUSER WIDTH (m)
C     BEQUIV=EQUIVALENT RADIUS FOR RECTANGULAR PLUME IN AMISK LAKE (m) 
C	CO2=DISSOLVED OXYGEN (DO) CONCENTRATION (mol/m3)
C     CO2M=DO CONCENTRATION PROFILE FOR INPUT BOUNDARY CONDITION (g/m3)
C	CN2=DISSOLVED NITROGEN CONCENTRATION (mol/m3)
C	COMG=DISSOLVED OXYGEN CONCENTRATION (g/m3)
C     COMGPT=DO CONCENTRATION OF PLUME DETRAINMENT AT TOP OF PLUME (g/m3) 
C	CNMG=DISSOLVED NITROGEN CONCENTRATION (g/m3)
C     DENSE20=DENSITY OF WATER AT 20 C
C	DENSEA=AMBIENT WATER DENSITY (kg/m3)
C	DENSEP=DENSITY OF THE PLUME (kg/m3)
C	DENSEW=WATER DENSITY IN PLUME (kg/m3)
C	DIAMM=BUBBLE DIAMETER (mm)
C     DIFFEL=DIFFUSER ELEVATION (m)
C     DMPR=DEPTH OF MAXIMUM PLUME RISE (m)
C     DNAMB=AMBIENT DISSOLVED NITROGEN CONCENTRATION (g/m3)
C     DOAMB=AMBIENT DISSOLVED OXYGEN CONCENTRATION (g/m3)
C	E=ENTRAINMENT FACTOR (m3/s)
C     ELEV=ELEVATION (m)
C     ELEVT=TERMINAL ELEVATION OF PLUME IN SEGMENT (m)
C	FDO=DISSOLVED OXYGEN FLUX (mol/s)
C	FDN=DISSOLVED NITROGEN FLUX (mol/s)
C	FRACO=MOLE FRACTION OF OXYGEN (-)
C	FRACN=MOLE FRACTION OF NITROGEN (-)
C     FRCONOT=INITIAL MOLE FRACTION OF OXYGEN IN DIFFUSER GAS SUPPLY, 0.21 OR 0.965 (-)
C	FSAL=SALINITY FLUX (kg/s)
C	FTEMP=TEMPERATURE FLUX (C m3/s)
C	FGO=GASEOUS OXYGEN FLUX (mol/s)
C     FGONOT=INITIAL GASEOUS OXYGEN FLUX (mol/s)
C	FGN=GASEOUS NITROGEN FLUX (mol/s)
C     FRCNOT=FRACTION OF NITROGEN IN ATMOSPHERE (-)
C     FRNOT=INITIAL FROUDE NUMBER (-)
C     GAMMA=SALINITY CONVERSION FACTOR [kg/m3/(uS/cm)]
C     GROSSMT=GROSS MASS TRANSFER OF OXYGEN FROM PLUME (kg/d)
C	HO=SOLUBILITY CONSTANT FOR OXYGEN (mol/m3/Pa)
C     HOD=HYPOLIMNETIC OXYGEN DEMAND FROM LITTLE AND MCGINNIS (2001) (kg/d)
C	HN=SOLUBILITY CONSTANT FOR NITROGEN (mol/m3/Pa)
C     JULDAY=JULIAN DAY IN GIVEN YEAR
C     KK=DIFFUSER SEGMENT (-)
C	KOLO=MASS TRANSFER COEFFICIENT FOR OXYGEN (m/s)
C	KOLN=MASS TRANSFER COEFFICIENT FOR NITROGEN (m/s)
C	L=DIFFUSER LENGTH (m)
C     LAKE=LAKE AND DIFFUSER TYPE (SHR AND LINEAR=1 OR AMISK AND RECTANGULAR=2) FOR SELECTION OF LAMBDA
C     LAYERS=NUMBER OF LAYERS/DATA POINTS IN BOUNDARY CONDITION PROFILES (-)
C     LDIFF=LENGTH OF DIFFUSER (m)
C     LNOT=INITIAL DIFFUSER LENGTH (m)
C	LAMBDA=FRACTION OF PLUME OCCUPIED BY BUBBLES (-)
C     LAMBNOT=LAMBDA x INITIAL PLUME RADIUS;EQUAL TO DIFFUSER RADIUS (m) 
C     LOCAT=DEPTHS FOR INPUT BOUNDARY CONDITION PROFILES (m)  
C	MOMENT=MOMENTUM (m4/s)
C	N=NUMBER OF BUBBLES PER SECOND (1/s)
C     NETMT=NET MASS TRANSFER OF OXYGEN ABOVE OXYGEN DEMAND (kg/d)
C     OTEFF=OXYGEN TRANSFER EFFICEINCY (%)
C     PATM=ATMOSPHERIC PRESSURE AT AVERAGE WSEL (Pa)
C	PSTD=STANDARD PRESSURE (Pa)
C	QSCFM=STANDARD GAS FLOW RATE (scfm), TOTAL GAS FLOW RATE TO DIFFUSER
C	QSCMS=STANDARD GAS FLOW RATE (scms)
C     QNM3HR=STANDARD GAS FLOW RATE (Nm3/hr)
C	QW=FLOWRATE OF WATER (m3/s)
C     QWT=TOTAL DETRAINMENT FLOW RATE AT THE TOP OF THE PLUME (m3/s)  
C	RB=BUBBLE RADIUS (m)
C	RGAS=IDEAL GAS CONSTANT (J/mol/K)
C     SAL=SALINITY (uS/cm)
C     SALAMB=SALINITY OF AMBIENT WATER (uS/cm)
C	SALPLU=SALINITY OF THE PLUME (uS/cm)
C	TAMB=AMBIENT WATER TEMPERATURE (C)
C     TAVG=AVERAGE AMBIENT WATER TEMPERATURE (C)
C     TDS=TOTAL DISSOLVED SOLIDS (g/m3) [0.64 conversion factor from Chapra book]
C     TE=TEMPERATURE PROFILE FOR INPUT BOUNDARY CONDITION (C) 
C	TPLUME=PLUME TEMPERATURE (C)
C     TPLUMET=DETRAINMENT PLUME TEMPERATURE AT THE TOP OF THE PLUME (C) 
C	TSTD=STANDARD TEMPERATURE (K)
C	V=WATER VELOCITY (m/s)
C     VAVG=AVERAGE WATER VELOCITY (m/s)
C	VB=BUBBLE RISE VELOCITY (m/s)
C	VBUB=BUBBLE VOLUME (m3)
C	VGUESS=GUESSED INITIAL WATER VELOCITY (m/s)
C     WSEL=WATER SURFACE ELEVATION (m)
C	YO2=GASEOUS OXYGEN CONCENTRATION (mol/m3)
C	YN2=GASEOUS NITROGEN CONCENTRATION (mol/m3)
C	Z=DEPTH TO DIFFUSER (m)
C
      REAL*8 ALPHA,AREA,B,CO2,COMG,CN2,CNMG,DS,DZ,DMOM,DFTEM,
     +DFGO,DFGN,DFDO,DFDN,DFSAL,DENSEA,DENSEP,DENSEW,DIAMM,E,FDO,
     +FDN,FRACO,FRACN,FSAL,FTEMP,FGO,FGN,G,GAMMA,HO2,HN2,KOLN,KOLO,
     +L,LAMBDA,MOMENT,N,PI,PO,PN,PSTD,PZ,QSCFM,QSCMS,QW,QGAS,RB,
     +RGAS,SAL0,SALAMB,SALPLU,TAMB,TPLUME,TSTD,V,VB,VBUB,VG,VGUESS,YO2,
     +YN2,Z,AA,BB,CC,BNOT,LNOT,LAYER,ELEV,TEST1,TEST2,DT,TE(1000),
     +XLOC,DEPTH,CO2M(1000),LOCAT(1000),DOAMB,COMGP,CNMGP,WSEL,PATM,
     +SAL(1000),DIFFEL,GROSSMT,NETMT,HOD,FGONOT,DNAMB,DMPR,
     +TAVG,SUMTEMP,SUMSAL,LAMBNOT,FRCONOT,
     +RBNOT,H,DYDX(8),Y(8),YOUT(8),DENSE20,OTEFF,FRNOT,VDIFF,FR,
     +QNM3HR,RiNOT,BUOY,DCO2,QGFRAC,LDIFF,TDS(1000),JDAY,EL(70),
     +INPUT,DELTAC,COMGNOT,ELEVT,QWT,TPLUMET,COMGPT,X,JULDAY,BEQUIV
      INTEGER II,IJ,IK,IN,JJ,KK,LL,MM,NEQN,NN,MI,JK,JL,LAKE
     +LAYERS,KM,LM,KN,KO,KP,ROWS,KQ,KR,KU,KV,KW,KX,KY,KZ,KS,YEAR  
C      
C                  
!     FJR - LDIFF passed as argument
      QGFRAC=1.0
C
      SUMTEMP=0.0
      DO 50 KS=1,LAYERS
          SUMTEMP=SUMTEMP+TE(KS)
   50 CONTINUE        
C                   
C     Assume that gas bubbles are composed of oxygen and nitrogen only.
      FRACO=FRCONOT 
      FRACN=1.0-FRACO
      DEPTH=WSEL-DIFFEL
      Z=DEPTH
      ELEV=DIFFEL
      X=0.
C
C     Interpolate raw AGPM output to obtain line plume initial conditions
      XLOC=DEPTH-X

      CALL LININT(LOCAT,TE,LAYERS,XLOC,TAMB)
      CALL LININT(LOCAT,CO2M,LAYERS,XLOC,COMG)
      COMGP=COMG
      DOAMB=COMG
      CO2=COMG/32.
      COMGNOT=COMG
!     Assume constant value for ambient salinity for initial runs
!     SALAMB passed as an argument
!      IF(YEAR.EQ.1998)SALAMB=193.
!      IF(YEAR.EQ.2003.AND.JULDAY.LE.196.)SALAMB=193.
!      IF(YEAR.EQ.2003.AND.JULDAY.GE.222.)SALAMB=201.
!      IF(YEAR.EQ.2004)SALAMB=218.     
      SALPLU=SALAMB    
C
      TPLUME=TAMB
      VGUESS=0.07
      V=VGUESS
C
C     CONSTANTS
      ALPHA=0.11
      G=9.80665
      GAMMA=6.9E-4
      IF(LAKE.EQ.1)THEN
         LAMBDA=0.93
         FRNOT=2.0
      ELSEIF(LAKE.EQ.2)THEN
         LAMBDA=0.8
         FRNOT=1.6
      ENDIF   
      PI=ACOS(-1.0)
      PSTD=101325.
      RGAS=8.314
      TSTD=293.15
      LAYER=0.0
      HOD=50.0
      DENSE20=998.2
      FRCNATM=0.79
!     VLS: PATM passed as an argument      
!      PATM=96261.0
!      IF(YEAR.EQ.1998)PATM=96379.    
C
C     Per e-mail from Paul Gantzer dated 4-30-04, total diffuser width is 6". 
!      LAMBNOT=0.0762
!     VLS: LAMBNOT passed as an argument      
      BNOT=LAMBNOT/LAMBDA
      B=BNOT
C     Per e-mail from Paul Gantzer dated 5-4-04, total diffuser length approx. 1000 ft for 2003.
C     For 2004, total diffuser length was approx. 2000 ft but only 25/60*2000=833 ft was active.
!     Revised LNOT to account for additional length due to spreading of velocity/water plume beyond bubble plume (4-16-09)  
      LNOT=LDIFF+2.0*BNOT*(1.0-LAMBDA)
      L=LNOT
      BEQUIV=0.5*(4.*LDIFF*2.*LAMBNOT/PI)**0.5
C      
C     AMBIENT AND AVERAGE WATER DENSITIES
      DENSEA=(0.059385*TAMB**3-8.56272*TAMB**2+65.4891*TAMB)*0.001
     ++999.84298+(GAMMA)*SALAMB
      DENSEW=DENSEA
C
C     BUBBLE PROPERTIES
C     Gas flow rate per segment asssumed to be proportional to fraction of total diffuser length.
      QSCMS=QGFRAC*QSCFM/3.281**3/60.0
      QGAS=PSTD*QSCMS*(TAMB+273.15)/((PATM+DENSEA*G*Z)*TSTD)
C     For diffuser in SHR, use correlation by McGinnis and Little (2000) for initial bubble size.     
!      DIAMM=1.12+0.938*(QGAS*60.0*60.0)/(L-2.0*BNOT*(1.0-LAMBDA))
!     VLS: DIAMM passed as an argument
      RB=DIAMM/2000.
      RBNOT=RB
      IF(RB.LE.(7.5E-4))THEN
            VB=1189.0*RB**1.1945
      ELSEIF(RB.GT.(7.5E-4).AND.RB.LT.(4.8E-3))THEN
            VB=0.22
      ELSE
            VB=2.995*RB**0.489
      ENDIF
C
      KOLO=0.6*RB
      IF(KOLO.GT.(4.0E-4))THEN
            KOLO=4.0E-4
      ENDIF
      KOLN=KOLO
C
      HO2=(2.125-0.05023*TPLUME+5.7714E-4*TPLUME**2)/100000.
      HN2=(1.042-0.02457*TPLUME+3.1714E-4*TPLUME**2)/100000.
C	
C     Assume initial ambient dissolved nitrogen conc. equals saturated conc. at surface.
      CN2=(PATM*FRCNATM)*HN2      
      CNMG=CN2*28.0
      CNMGP=CNMG
      DNAMB=CNMG
C
C     CALCULATION OF INITIAL WATER VELOCITY USING FROUDE NUMBER
      VBUB=4./3.*PI*RB**3
      N=QGAS/VBUB
  9   VG=QGAS/((VGUESS+VB)*(2.*LAMBDA*B)*(L-2.0*B*(1.0-LAMBDA)))
      DENSEP=(1.0-VG)*DENSEW
      IF(LAKE.EQ.1)THEN
         V=FRNOT*(2.0*LAMBDA*B*G*(DENSEA-DENSEP)/DENSEP)**0.5
      ELSEIF(LAKE.EQ.2)THEN
         V=FRNOT*(2.0*BEQUIV*G*(DENSEA-DENSEP)/DENSEP)**0.5
      ENDIF      
      VDIFF=ABS(V-VGUESS)
      IF(VDIFF.GT.1.0E-6)THEN
         VGUESS=V
         GOTO 9
      ENDIF
C
C     VARIABLE TRANSFORMATION      
      E=2.*(L+2.*B)*ALPHA*V    
      QW=2.*L*B*V
      MOMENT=2.*L*B*V**2
      FTEMP=QW*TPLUME 
      FSAL=QW*(SALPLU*GAMMA/DENSE20)*DENSEW
C     Previous equation corrected to account for salinity units conversion.        
      FDO=QW*CO2
      FDN=QW*CN2
      FGO=PSTD*QSCMS/(RGAS*TSTD)*FRACO
      FGONOT=FGO
      FGN=PSTD*QSCMS/(RGAS*TSTD)*FRACN
C     Revised gaseous flux equations.
      YO2=FGO/(LAMBDA*2.*B*(L-2.*B*(1.-LAMBDA))*(V+VB))
      YN2=FGN/(LAMBDA*2.*B*(L-2.*B*(1.-LAMBDA))*(V+VB))           
      PZ=PATM+(DENSEA*G*Z)
      PO=PZ*FRACO
      PN=PZ*FRACN
      BUOY=(G*(DENSEA-DENSEP)/DENSEP*QW)/LNOT
      RiNOT=((QW/LNOT)**2)*BUOY**(2./3.)/(MOMENT/LNOT)**2
C
      TDS=SALPLU*0.64
C      
C	SOLUTION PROCEEDURE
      DZ=0.001
      H=0.001
      MM=1
 10   Z=Z-DZ
      X=X+DZ
      LAYER=LAYER+DZ
      ELEV=ELEV+DZ
      MM=MM+1
      COUNT2=COUNT2+1
C	
C     Interpolate raw AGPM output to obtain line plume boundary conditions
      XLOC=DEPTH-X
      CALL LININT(LOCAT,CO2M,LAYERS,XLOC,COMG)
      DOAMB=COMG
      CALL LININT(LOCAT,TE,LAYERS,XLOC,TAMB)
C                  
C     Use subroutines for Runge Kutta method solution
      NEQN=8
      Y(1)=QW
      Y(2)=MOMENT
      Y(3)=FTEMP
      Y(4)=FSAL
      Y(5)=FDO
      Y(6)=FDN
      Y(7)=FGO
      Y(8)=FGN
      CALL DERIVS(E,DENSEA,DENSEW,DENSEP,G,L,B,LAMBDA,TAMB,SALAMB,
     +            GAMMA,DENSE20,DOAMB,PI,RB,N,V,VB,KOLO,HO2,PO,
     +            COMGP,DNAMB,KOLN,HN2,PN,CNMGP,Z,Y,DYDX)
      CALL RK4(E,DENSEA,DENSEW,DENSEP,G,L,B,LAMBDA,TAMB,SALAMB,
     +         GAMMA,DENSE20,DOAMB,PI,RB,N,V,VB,KOLO,HO2,PO,
     +         COMGP,DNAMB,KOLN,HN2,PN,CNMGP,Y,DYDX,NEQN,Z,H,YOUT)
C       
      QW=YOUT(1)
      MOMENT=YOUT(2)
      FTEMP=YOUT(3)
      FSAL=YOUT(4)
      FDO=YOUT(5)
      FDN=YOUT(6)
      FGO=YOUT(7)
      FGN=YOUT(8)
      IF(MOMENT.LT.0.0)THEN
         TPLUME=FTEMP/QW
	   SALPLU=FSAL/(QW*DENSEW)/(GAMMA/DENSE20)
C    Previous equation corrected to consistently express salinity in uS/cm	   
	   CO2=FDO/QW
	   CN2=FDN/QW
	   GOTO 20
      ENDIF
      V=MOMENT/QW
      AREA=QW/V
C     SOLVE FOR DIMENSIONS USING L^2+(2Bo-Lo)L-AREA=0 USING QUADRATIC EQN.
      AA=1.0
      BB=2.*BNOT-LNOT
      CC=-1.0*AREA
      L=(-1.0*BB+(BB**2-4.0*AA*CC)**(0.5))/(2.0*AA)
      IF(L.LT.0.0)THEN
           L=(-1.0*BB-(BB**2-4.0*AA*CC)**(0.5))/(2.0*AA)
      ENDIF
      B=AREA/(2.0*L)
C      
      E=2.*(L+2.*B)*ALPHA*V
      TPLUME=FTEMP/QW
 	SALPLU=FSAL/(QW*DENSEW)/(GAMMA/DENSE20)
C     Previous equation corrected to consistently express salinity in uS/cm
      CO2=FDO/QW
      CN2=FDN/QW
      COMGP=CO2*32.
      CNMGP=CN2*28.
C     Revised gaseous flux equations.
      YO2=FGO/(LAMBDA*2.*B*(L-2.*B*(1.-LAMBDA))*(V+VB))
      YN2=FGN/(LAMBDA*2.*B*(L-2.*B*(1.-LAMBDA))*(V+VB))
C  
      PZ=PATM+(DENSEA*G*Z)
      QGAS=(FGO+FGN)*RGAS*(TPLUME+273.15)/PZ
      VBUB=QGAS/N
C
      VG=VBUB*N/((V+VB)*(2.*LAMBDA*B)*(L-2.*B*(1.-LAMBDA)))
C     Previous equation revised to account for correct plume cross-sectional area occupied by bubbles.      
      RB=(3.*QGAS/(4.*PI*N))**(1./3.)
      IF(RB.LT.0.0)THEN
           RB=1.0E-8
      ENDIF
      FRACO=FGO/(FGO+FGN)
      FRACN=1.0-FRACO
C	
      PO=PZ*FRACO
      PN=PZ*FRACN
      DENSEA=(0.059385*TAMB**3-8.56272*TAMB**2+65.4891*TAMB)*0.001
     ++999.84298+(GAMMA)*SALAMB
      DENSEW=(0.059385*TPLUME**3-8.56272*TPLUME**2+65.4891*TPLUME)*0.001
     ++999.84298+(GAMMA)*SALPLU
C     Previous equation re-revised to account for correct salinity units (uS/cm) in density calculations.      
      DENSEP=(1.0-VG)*DENSEW
C
C	BUBBLE PROPERTIES
      IF(RB.LE.(7.5E-4))THEN
           VB=1189.0*RB**1.1945
      ELSEIF(RB.GT.(7.5E-4).AND.RB.LT.(4.8E-3))THEN
           VB=0.22
      ELSE
           VB=2.995*RB**0.489
      ENDIF
C
      KOLO=0.6*RB
      IF(KOLO.GT.(4.0E-4))THEN
            KOLO=4.0E-4
      ENDIF
      KOLN=KOLO
C
      HO2=(2.125-0.05023*TPLUME+5.7714E-4*TPLUME**2)/100000.
      HN2=(1.042-0.02457*TPLUME+3.1714E-4*TPLUME**2)/100000.
C
      FR=V/(2.*LAMBDA*B*G*(DENSEA-DENSEP)/DENSEP)**0.5
      DCO2=HO2*PO-CO2
C      
      IF(V.GT.1.E-6)THEN
!     Revised from "Z.GT.0.0" to "XLOC.GT.0.0" by VLS on 10-3-07       
           IF(XLOC.GT.0.001)THEN
	        GOTO 10
	    ENDIF
      ENDIF
C
C     CALCULATION OF AVERAGE NET OXYGEN MASS TRANSFER FOR DAY
   20 GROSSMT=(FGONOT-FGO)*32./1000.*86400.
      OTEFF=(FGONOT-FGO)/FGONOT*100.
      NETMT=GROSSMT-HOD
      QNM3HR=(QSCMS*PSTD/TSTD*273.15/10**5)*3600.0   
      DELTAC=COMGP-COMGNOT
C
      ELEVT=ELEV
      QWT=QW
      TPLUMET=TPLUME
      COMGPT=COMGP
                          
!  100 CONTINUE
      RETURN
      END
C
C------------------------------------------------------------------------------
C
      SUBROUTINE RK4(E,DENSEA,DENSEW,DENSEP,G,L,B,LAMBDA,TAMB,
     +              SALAMB,GAMMA,DENSE20,DOAMB,PI,RB,N,V,VB,KOLO,HO2,
     +              PO,COMGP,DNAMB,KOLN,HN2,PN,CNMGP,Y,DYDX,NN,X,H,YOUT)

      INTEGER I,NN,NMAX
      PARAMETER (NMAX=50)
      REAL*8 E,DENSEA,DENSEW,DENSEP,G,L,B,LAMBDA,TAMB,SALAMB,GAMMA,
     +DENSE20,DOAMB,PI,RB,N,V,VB,KOLO,HO2,PO,COMGP,DNAMB,KOLN,HN2,PN,
     +CNMGP,H,X,DYDX(NN),Y(NN),YOUT(NN),H6,HH,XH,DYM(NMAX),DYT(NMAX),
     +YT(NMAX)
      EXTERNAL DERIVS
      HH=H*0.5
      H6=H/6.
      XH=X+HH
      DO 11 I=1,NN
          YT(I)=Y(I)+HH*DYDX(I)
   11 CONTINUE
      CALL DERIVS(E,DENSEA,DENSEW,DENSEP,G,L,B,LAMBDA,TAMB,SALAMB,
     +           GAMMA,DENSE20,DOAMB,PI,RB,N,V,VB,KOLO,HO2,PO,
     +           COMGP,DNAMB,KOLN,HN2,PN,CNMGP,XH,YT,DYT)
      DO 12 I=1,NN
          YT(I)=Y(I)+HH*DYT(I)
   12 CONTINUE
      CALL DERIVS(E,DENSEA,DENSEW,DENSEP,G,L,B,LAMBDA,TAMB,SALAMB,
     +            GAMMA,DENSE20,DOAMB,PI,RB,N,V,VB,KOLO,HO2,PO,
     +            COMGP,DNAMB,KOLN,HN2,PN,CNMGP,XH,YT,DYM)
      DO 13 I=1,NN
          YT(I)=Y(I)+H*DYM(I)
          DYM(I)=DYT(I)+DYM(I)
   13 CONTINUE
      CALL DERIVS(E,DENSEA,DENSEW,DENSEP,G,L,B,LAMBDA,TAMB,SALAMB,
     +            GAMMA,DENSE20,DOAMB,PI,RB,N,V,VB,KOLO,HO2,PO,
     +            COMGP,DNAMB,KOLN,HN2,PN,CNMGP,X+H,YT,DYT)
      DO 14 I=1,NN
          YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
   14 CONTINUE
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE DERIVS(E,DENSEA,DENSEW,DENSEP,G,L,B,LAMBDA,TAMB,
     +              SALAMB,GAMMA,DENSE20,DOAMB,PI,RB,N,V,VB,KOLO,HO2,PO,
     +              COMGP,DNAMB,KOLN,HN2,PN,CNMGP,X,Y,DYDX)         
      REAL*8 E,DENSEA,DENSEW,DENSEP,G,L,B,LAMBDA,TAMB,SALAMB,GAMMA,
     +DENSE20,DOAMB,PI,RB,N,V,VB,KOLO,HO2,PO,COMGP,DNAMB,KOLN,HN2,PN,
     +CNMGP,X,Y(8),DYDX(8)
C     Right-hand side of differential equations for Runge-Kutta solution
      DYDX(1)=E         
      DYDX(2)=(DENSEA-DENSEW)/DENSEP*G*L*2.*B+(DENSEW-DENSEP)/
     +DENSEP*G*LAMBDA*2.*B*(L-2.*B*(1.-LAMBDA))
      DYDX(3)=E*TAMB
      DYDX(4)=E*(SALAMB*GAMMA/DENSE20)*DENSEA
      DYDX(5)=(E*DOAMB/32.+4.0*PI*RB**2*N/(V+VB)*KOLO*
     +(HO2*PO-COMGP/32.))
      DYDX(6)=(E*DNAMB/28.+4.0*PI*RB**2*N/(V+VB)*KOLN*
     +(HN2*PN-CNMGP/28.))
      DYDX(7)=-4.0*PI*RB**2*N/(V+VB)*KOLO*(HO2*PO-COMGP/32.)
      DYDX(8)=-4.0*PI*RB**2*N/(V+VB)*KOLN*(HN2*PN-CNMGP/28.)
      RETURN
      END 
C
C----------------------------------------------------------------------
C
      SUBROUTINE LININT(XTAB,YTAB,NTAB,X,Y)
      INTEGER I,NTAB
      REAL*8 X,Y,XTAB(NTAB),YTAB(NTAB)
      IF(X.LT.XTAB(1).OR.X.GT.XTAB(NTAB))THEN
          WRITE(*,*)'X = ',X,' IS OUT OF TABLE RANGE'
          PAUSE
      ENDIF
      DO 100 I=2,NTAB
          IF(X.LE.XTAB(I)) GOTO 200
 100  CONTINUE
 200  I1=I-1
      WX=(X-XTAB(I1))/(XTAB(I1+1)-XTAB(I1))
      Y=(1.-WX)*YTAB(I1)+WX*YTAB(I1+1)
      RETURN
      END                   
                          
             
C     *******************************************************************************
	SUBROUTINE LINEPLU_v1(YEAR,JULDAY,WSEL,DIFFEL,LAYERS,LDIFF,LAMBNOT,
     +SALAMB,PATM,DIAMM,LAKE,QSCFM,FRCONOT,LAYDIFF,HCELL,LOCAT,TE,CO2M,
     +ELEVT,QWT,TPLUMET,COMGPT,QWD,LAYTOP)
C     *******************************************************************************
C	THIS SUBROUTINE IS WRITTEN TO PREDICT THE PERFORMANCE OF A LINEAR BUBBLE PLUME.  
C     THE MODEL IS BASED ON THE WUEST ET AL. (1992) CIRCULAR BUBBLE PLUME MODEL.
C     By: Vickie Singleton
C	
C	VERSION 3 (to couple with Francisco Rueda's reservoir model and for Amisk Lake areal diffuser)
C
C     This version includes the following:
C     1.  Revised momentum flux equations re-derived on October 4, 2004
C     2.  Revised gaseous flux equations with correct plume area 
C     3.  Revised dissolved flux equations to use ambient concentrations for entrainment
C     4.  Revised dissolved flux equations to use plume concentrations for gas transfer 
C     5.  Correction of ambient salinity interpolation from input data file
C     6.  Corrected initial bubble size correlation equation for gas flow rate per unit length (11-17-04)
C     7.  Correction for salinity units in salinity flux equations (2-2-05)  
C     8.  Revised gas holdup equation in loop to account for correct plume cross-sectional area occupied by bubbles (2-2-05)
C     9.  Determination and use of average ambient density in pressure calculations (2-2-05)   
C     10.  Correction of Bnot/diffuser source radius calculation per Wuest et al. 1992, Figure 2. (2-2-05)
C     11.  Revision of numerical integration method from Euler to fourth-order Runge-Kutta (2-23-05)
C     12.  Added calculation of initial water velocity using initial Froude number of 1.6.  (Previously, initial
C          water velocity was assumed to be 0.07 m/s per Dan's original program.) (6-22-05)
C     13.  Revised entrainment and spreading coefficients from 0.08 and 0.85, respectively, to 0.11 and 0.93, 
C          respectively, to account for top-hat profiles versus Gaussian profiles (9-19-05)
C     14.  Revised initial Froude number from 1.6 to 2.0 (refer to calculations). (10-1-05)
C     15.  Revised characteristic length in Froude number calculation from equivalent radius to initial plume width. (10-1-05) 
C     16.  Corrected equations that include salinity to make units consistent (11-7-06)
C     17.  Output entrainment/withdrawal per cell (7-24-07)
C     18.  Passed LAMBNOT, SALAMB, PATM, DIAMM, and LAKE as arguments to accomodate Amisk Lake (4-19-09)
C
C	April 14, 2009
C
C	VARIABLES
C
C	ALPHA=ENTRAINMENT COEFFICIENT (-)
C     B=1/2 DIFFUSER WIDTH (m)
C     BAVG=AVERAGE 1/2 DIFFUSER WIDTH (m)
C     BEQUIV=EQUIVALENT RADIUS FOR RECTANGULAR PLUME IN AMISK LAKE (m)  
C	CO2=DISSOLVED OXYGEN (DO) CONCENTRATION (mol/m3)
C     CO2M=DO CONCENTRATION PROFILE FOR INPUT BOUNDARY CONDITION (g/m3)
C	CN2=DISSOLVED NITROGEN CONCENTRATION (mol/m3)
C	COMG=DISSOLVED OXYGEN CONCENTRATION (g/m3)
C     COMGPT=DO CONCENTRATION OF PLUME DETRAINMENT AT TOP OF PLUME (g/m3) 
C	CNMG=DISSOLVED NITROGEN CONCENTRATION (g/m3)
C     DENSE20=DENSITY OF WATER AT 20 C
C	DENSEA=AMBIENT WATER DENSITY (kg/m3)
C	DENSEP=DENSITY OF THE PLUME (kg/m3)
C	DENSEW=WATER DENSITY IN PLUME (kg/m3)
C	DIAMM=BUBBLE DIAMETER (mm)
C     DIFFEL=DIFFUSER ELEVATION (m)
C     DMPR=DEPTH OF MAXIMUM PLUME RISE (m)
C     DNAMB=AMBIENT DISSOLVED NITROGEN CONCENTRATION (g/m3)
C     DOAMB=AMBIENT DISSOLVED OXYGEN CONCENTRATION (g/m3)
C	E=ENTRAINMENT FACTOR (m3/s)
C     ELEV=ELEVATION (m)
C     ELEVT=TERMINAL ELEVATION OF PLUME IN SEGMENT (m)
C	FDO=DISSOLVED OXYGEN FLUX (mol/s)
C	FDN=DISSOLVED NITROGEN FLUX (mol/s)
C	FRACO=MOLE FRACTION OF OXYGEN (-)
C	FRACN=MOLE FRACTION OF NITROGEN (-)
C     FRCONOT=INITIAL MOLE FRACTION OF OXYGEN IN DIFFUSER GAS SUPPLY, 0.21 OR 0.965 (-)
C	FSAL=SALINITY FLUX (kg/s)
C	FTEMP=TEMPERATURE FLUX (C m3/s)
C	FGO=GASEOUS OXYGEN FLUX (mol/s)
C     FGONOT=INITIAL GASEOUS OXYGEN FLUX (mol/s)
C	FGN=GASEOUS NITROGEN FLUX (mol/s)
C     FRCNOT=FRACTION OF NITROGEN IN ATMOSPHERE (-)
C     FRNOT=INITIAL FROUDE NUMBER (-)
C     GAMMA=SALINITY CONVERSION FACTOR [kg/m3/(uS/cm)]
C     GROSSMT=GROSS MASS TRANSFER OF OXYGEN FROM PLUME (kg/d)
C     HCELL=HEIGHT OF CELL IN GRID (m)
C	HO=SOLUBILITY CONSTANT FOR OXYGEN (mol/m3/Pa)
C	HN=SOLUBILITY CONSTANT FOR NITROGEN (mol/m3/Pa)
C     HWITH=HEIGHT OF WITHDRAWAL/ENTRAINMENT ZONE (m)
C     JULDAY=JULIAN DAY IN GIVEN YEAR
C	KOLO=MASS TRANSFER COEFFICIENT FOR OXYGEN (m/s)
C	KOLN=MASS TRANSFER COEFFICIENT FOR NITROGEN (m/s)
C	L=DIFFUSER LENGTH (m)
C     LAKE=LAKE AND DIFFUSER TYPE (SHR AND LINEAR=1 OR AMISK AND RECTANGULAR=2) FOR SELECTION OF LAMBDA
C     LAYERS=NUMBER OF LAYERS/DATA POINTS IN BOUNDARY CONDITION PROFILES (-)
C     LAYDIFF=GRID LAYER CORRESPONDING TO DIFFUSER DEPTH (-)
C     LAYTOP=GRID LAYER CORRESPONDING TO TOP OF PLUME (-)
C     LDIFF=LENGTH OF DIFFUSER (m)
C     LNOT=INITIAL DIFFUSER LENGTH (m)
C	LAMBDA=FRACTION OF PLUME OCCUPIED BY BUBBLES (-)
C     LAMBNOT=LAMBDA x INITIAL PLUME RADIUS;EQUAL TO DIFFUSER RADIUS (m) 
C     LOCAT=DEPTHS FOR INPUT BOUNDARY CONDITION PROFILES (m)  
C	MOMENT=MOMENTUM (m4/s)
C	N=NUMBER OF BUBBLES PER SECOND (1/s)
C     OTEFF=OXYGEN TRANSFER EFFICEINCY (%)
C     PATM=ATMOSPHERIC PRESSURE AT AVERAGE WSEL (Pa)
C	PSTD=STANDARD PRESSURE (Pa)
C	QSCFM=STANDARD GAS FLOW RATE (scfm), TOTAL GAS FLOW RATE TO DIFFUSER
C	QSCMS=STANDARD GAS FLOW RATE (scms)
C	QW=FLOWRATE OF WATER (m3/s)
C     QWT=TOTAL DETRAINMENT FLOW RATE AT THE TOP OF THE PLUME (m3/s)  
C	RB=BUBBLE RADIUS (m)
C	RGAS=IDEAL GAS CONSTANT (J/mol/K)
C     SAL=SALINITY (uS/cm)
C     SALAMB=SALINITY OF AMBIENT WATER (uS/cm)
C	SALPLU=SALINITY OF THE PLUME (uS/cm)
C	TAMB=AMBIENT WATER TEMPERATURE (C)
C     TAVG=AVERAGE AMBIENT WATER TEMPERATURE (C)
C     TDS=TOTAL DISSOLVED SOLIDS (g/m3) [0.64 conversion factor from Chapra book]
C     TE=TEMPERATURE PROFILE FOR INPUT BOUNDARY CONDITION (C) 
C	TPLUME=PLUME TEMPERATURE (C)
C     TPLUMET=DETRAINMENT PLUME TEMPERATURE AT THE TOP OF THE PLUME (C) 
C	TSTD=STANDARD TEMPERATURE (K)
C	V=WATER VELOCITY (m/s)
C     VAVG=AVERAGE WATER VELOCITY (m/s)
C	VB=BUBBLE RISE VELOCITY (m/s)
C	VBUB=BUBBLE VOLUME (m3)
C	VGUESS=GUESSED INITIAL WATER VELOCITY (m/s)
C     WSEL=WATER SURFACE ELEVATION (m)
C	YO2=GASEOUS OXYGEN CONCENTRATION (mol/m3)
C	YN2=GASEOUS NITROGEN CONCENTRATION (mol/m3)
C	Z=DEPTH TO DIFFUSER (m)
C
      REAL*8 ALPHA,AREA,B,CO2,COMG,CN2,CNMG,DS,DZ,DENSEA,DENSEP,DENSEW,
     +DIAMM,E,FDO,FDN,FRACO,FRACN,FSAL,FTEMP,FGO,FGN,G,GAMMA,HO2,HN2,
     +KOLN,KOLO,L,LAMBDA,MOMENT,N,PI,PO,PN,PSTD,PZ,QSCFM,QSCMS,QW,QGAS,
     +RB,RGAS,SALAMB,SALPLU,TAMB,TPLUME,TSTD,V,VB,VBUB,VG,VGUESS,
     +YO2,YN2,Z,AA,BB,CC,BNOT,LNOT,ELEV,DT,TE(1000),XLOC,DEPTH,
     +CO2M(1000),LOCAT(1000),DOAMB,COMGP,CNMGP,WSEL,PATM,SAL(1000),
     +DIFFEL,GROSSMT,FGONOT,DNAMB,DMPR,TAVG,SUMTEMP,SUMSAL,
     +LAMBNOT,FRCONOT,RBNOT,H,DYDX(8),Y(8),YOUT(8),
     +DENSE20,OTEFF,FRNOT,VDIFF,FR,BUOY,DCO2,QGFRAC,LDIFF,TDS(1000),
     +JDAY,EL(70),DELTAC,COMGNOT,ELEVT,QWT,TPLUMET,COMGPT,QWD(500),
     +HWITH,HCELL,JULDAY,BTOP,BEQUIV
      INTEGER II,IJ,IK,IN,JJ,LL,NEQN,NN,MI,JK,JL,LAKE,LAYTOP
     +LAYERS,KM,KN,KO,KP,ROWS,KQ,KR,KU,KV,KW,KX,KY,KZ,KS,LAYDIFF,YEAR  
C 
!     FJR - LDIFF passed as argument
      QGFRAC=1.0      
      SUMTEMP=0.0
      DO 50 KS=1,LAYERS
          SUMTEMP=SUMTEMP+TE(KS)
   50 CONTINUE        
      TAVG=SUMTEMP/LAYERS             
C                   
C     Assume that gas bubbles are composed of oxygen and nitrogen only.
      FRACO=FRCONOT 
      FRACN=1.0-FRACO
      DEPTH=WSEL-DIFFEL
      Z=DEPTH
      ELEV=DIFFEL
      X=0.
C
C     Interpolate input profiles to obtain line plume initial conditions
      XLOC=DEPTH-X

      CALL LININT(LOCAT,TE,LAYERS,XLOC,TAMB)
      CALL LININT(LOCAT,CO2M,LAYERS,XLOC,COMG)
      COMGP=COMG
      DOAMB=COMG
      CO2=COMG/32.
      COMGNOT=COMG
!     Assume constant value for ambient salinity for initial runs
!     VLS: SALAMB passed as an argument
!      IF(YEAR.EQ.1998)SALAMB=193.
!      IF(YEAR.EQ.2003.AND.JULDAY.LE.196.)SALAMB=193.
!      IF(YEAR.EQ.2003.AND.JULDAY.GE.222.)SALAMB=201.
!      IF(YEAR.EQ.2004)SALAMB=218.     
      SALPLU=SALAMB    
      TPLUME=TAMB
      VGUESS=0.07
      V=VGUESS
C
C     CONSTANTS
      ALPHA=0.11
      G=9.80665
      GAMMA=6.9E-4
      IF(LAKE.EQ.1)THEN
         LAMBDA=0.93
         FRNOT=2.0
      ELSEIF(LAKE.EQ.2)THEN
         LAMBDA=0.8
         FRNOT=1.6
      ENDIF   
      PI=ACOS(-1.0)
      PSTD=101325.
      RGAS=8.314
      TSTD=293.15
      DENSE20=998.2
      FRCNATM=0.79
!     VLS: PATM passed as an argument      
!      PATM=96261.0
!      IF(YEAR.EQ.1998)PATM=96379.    
C
C     Per e-mail from Paul Gantzer dated 4-30-04, total diffuser width is 6".
!     VLS: LAMBNOT passed as an argument 
!      LAMBNOT=0.0762      
      BNOT=LAMBNOT/LAMBDA
      B=BNOT
!     Revised LNOT to account for additional length due to spreading of velocity/water plume beyond bubble plume (4-16-09)  
      LNOT=LDIFF+2.0*BNOT*(1.0-LAMBDA)
      L=LNOT
      BEQUIV=0.5*(4.*LDIFF*2.*LAMBNOT/PI)**0.5
C      
C     AMBIENT AND AVERAGE WATER DENSITIES
      DENSEA=(0.059385*TAMB**3-8.56272*TAMB**2+65.4891*TAMB)*0.001
     ++999.84298+(GAMMA)*SALAMB
      DENSEW=DENSEA
C
C     BUBBLE PROPERTIES
C     Gas flow rate per segment asssumed to be proportional to fraction of total diffuser length.
      QSCMS=QGFRAC*QSCFM/3.281**3/60.0
      QGAS=PSTD*QSCMS*(TAMB+273.15)/((PATM+DENSEA*G*Z)*TSTD)
C     For diffuser in SHR, use correlation by McGinnis and Little (2000) for initial bubble size.     
!      DIAMM=1.12+0.938*(QGAS*60.0*60.0)/(L-2.0*BNOT*(1.0-LAMBDA))
!     VLS: DIAMM passed as argument
      RB=DIAMM/2000.
      RBNOT=RB
      IF(RB.LE.(7.5E-4))THEN
            VB=1189.0*RB**1.1945
      ELSEIF(RB.GT.(7.5E-4).AND.RB.LT.(4.8E-3))THEN
            VB=0.22
      ELSE
            VB=2.995*RB**0.489
      ENDIF
C
      KOLO=0.6*RB
      IF(KOLO.GT.(4.0E-4))THEN
            KOLO=4.0E-4
      ENDIF
      KOLN=KOLO
C
      HO2=(2.125-0.05023*TPLUME+5.7714E-4*TPLUME**2)/100000.
      HN2=(1.042-0.02457*TPLUME+3.1714E-4*TPLUME**2)/100000.
C	
C     Assume initial ambient dissolved nitrogen conc. equals saturated conc. at surface.
      CN2=(PATM*FRCNATM)*HN2      
      CNMG=CN2*28.0
      CNMGP=CNMG
      DNAMB=CNMG
C
C     CALCULATION OF INITIAL WATER VELOCITY USING FROUDE NUMBER
      VBUB=4./3.*PI*RB**3
      N=QGAS/VBUB
  9   VG=QGAS/((VGUESS+VB)*(2.*LAMBDA*B)*(L-2.0*B*(1.0-LAMBDA)))
      DENSEP=(1.0-VG)*DENSEW
      IF(LAKE.EQ.1)THEN
         V=FRNOT*(2.0*LAMBDA*B*G*(DENSEA-DENSEP)/DENSEP)**0.5
      ELSEIF(LAKE.EQ.2)THEN
!         V=FRNOT*(2.0*BEQUIV*G*(DENSEA-DENSEP)/DENSEP)**0.5
!         VLS: For testing purposes, assume that characteristic length is equal to rectangle width (6-4-09)
          V=FRNOT*(2.0*LAMBDA*B*G*(DENSEA-DENSEP)/DENSEP)**0.5
      ENDIF      
      VDIFF=ABS(V-VGUESS)
      IF(VDIFF.GT.1.0E-6)THEN
         VGUESS=V
         GOTO 9
      ENDIF
C
C     VARIABLE TRANSFORMATION      
      E=2.*(L+2.*B)*ALPHA*V    
      QW=2.*L*B*V
      MOMENT=2.*L*B*V**2
      FTEMP=QW*TPLUME 
      FSAL=QW*(SALPLU*GAMMA/DENSE20)*DENSEW
C     Previous equation corrected to account for salinity units conversion.        
      FDO=QW*CO2
      FDN=QW*CN2
      FGO=PSTD*QSCMS/(RGAS*TSTD)*FRACO
      FGONOT=FGO
      FGN=PSTD*QSCMS/(RGAS*TSTD)*FRACN
C     Revised gaseous flux equations.
      YO2=FGO/(LAMBDA*2.*B*(L-2.*B*(1.-LAMBDA))*(V+VB))
      YN2=FGN/(LAMBDA*2.*B*(L-2.*B*(1.-LAMBDA))*(V+VB))           
      PZ=PATM+(DENSEA*G*Z)
      PO=PZ*FRACO
      PN=PZ*FRACN
      BUOY=(G*(DENSEA-DENSEP)/DENSEP*QW)/LNOT
      TDS=SALPLU*0.64
!     Initialize lateral withdrawal flowrate for first/lowest cell in column/segment
      JJ=0
      QWD(LAYDIFF)=QW      
C      
C	SOLUTION PROCEEDURE
      DZ=0.001
      H=0.001
      HWITH=0.0
 10   Z=Z-DZ
      X=X+DZ
      ELEV=ELEV+DZ
C	
C     Interpolate input profiles to obtain line plume boundary conditions
      XLOC=DEPTH-X
      CALL LININT(LOCAT,CO2M,LAYERS,XLOC,COMG)
      DOAMB=COMG
      CALL LININT(LOCAT,TE,LAYERS,XLOC,TAMB)
C                  
C     Use subroutines for Runge Kutta method solution
      NEQN=8
      Y(1)=QW
      Y(2)=MOMENT
      Y(3)=FTEMP
      Y(4)=FSAL
      Y(5)=FDO
      Y(6)=FDN
      Y(7)=FGO
      Y(8)=FGN
      CALL DERIVS(E,DENSEA,DENSEW,DENSEP,G,L,B,LAMBDA,TAMB,SALAMB,
     +            GAMMA,DENSE20,DOAMB,PI,RB,N,V,VB,KOLO,HO2,PO,
     +            COMGP,DNAMB,KOLN,HN2,PN,CNMGP,Z,Y,DYDX)
      CALL RK4(E,DENSEA,DENSEW,DENSEP,G,L,B,LAMBDA,TAMB,SALAMB,
     +         GAMMA,DENSE20,DOAMB,PI,RB,N,V,VB,KOLO,HO2,PO,
     +         COMGP,DNAMB,KOLN,HN2,PN,CNMGP,Y,DYDX,NEQN,Z,H,YOUT)
C       
      QW=YOUT(1)
      MOMENT=YOUT(2)
      FTEMP=YOUT(3)
      FSAL=YOUT(4)
      FDO=YOUT(5)
      FDN=YOUT(6)
      FGO=YOUT(7)
      FGN=YOUT(8)
      IF(MOMENT.LT.0.0)THEN
         TPLUME=FTEMP/QW
	   SALPLU=FSAL/(QW*DENSEW)/(GAMMA/DENSE20)
C    Previous equation corrected to consistently express salinity in uS/cm	   
	   CO2=FDO/QW
	   CN2=FDN/QW
	   GOTO 20
      ENDIF
      V=MOMENT/QW
      AREA=QW/V
C     SOLVE FOR DIMENSIONS USING L^2+(2Bo-Lo)L-AREA=0 USING QUADRATIC EQN.
      AA=1.0
      BB=2.*BNOT-LNOT
      CC=-1.0*AREA
      L=(-1.0*BB+(BB**2-4.0*AA*CC)**(0.5))/(2.0*AA)
      IF(L.LT.0.0)THEN
           L=(-1.0*BB-(BB**2-4.0*AA*CC)**(0.5))/(2.0*AA)
      ENDIF
      B=AREA/(2.0*L)
C      
      E=2.*(L+2.*B)*ALPHA*V
!     Add incremental entrainment to total cell entrainment/withdrawal       
      QWD(LAYDIFF-JJ)=QWD(LAYDIFF-JJ)+E*DZ
      HWITH=HWITH+DZ   
      IF(HWITH.GT.HCELL)THEN
          JJ=JJ+1
          HWITH=0.0
          QWD(LAYDIFF-JJ)=0.0
      ENDIF    
      TPLUME=FTEMP/QW
 	SALPLU=FSAL/(QW*DENSEW)/(GAMMA/DENSE20)
C     Previous equation corrected to consistently express salinity in uS/cm
      CO2=FDO/QW
      CN2=FDN/QW
      COMGP=CO2*32.
      CNMGP=CN2*28.
C     Revised gaseous flux equations.
      YO2=FGO/(LAMBDA*2.*B*(L-2.*B*(1.-LAMBDA))*(V+VB))
      YN2=FGN/(LAMBDA*2.*B*(L-2.*B*(1.-LAMBDA))*(V+VB))
C  
      PZ=PATM+(DENSEA*G*Z)
      QGAS=(FGO+FGN)*RGAS*(TPLUME+273.15)/PZ
      VBUB=QGAS/N
C
      VG=VBUB*N/((V+VB)*(2.*LAMBDA*B)*(L-2.*B*(1.-LAMBDA)))
C     Previous equation revised to account for correct plume cross-sectional area occupied by bubbles.      
      RB=(3.*QGAS/(4.*PI*N))**(1./3.)
      IF(RB.LT.0.0)THEN
           RB=1.0E-8
      ENDIF
      FRACO=FGO/(FGO+FGN)
      FRACN=1.0-FRACO
C	
      PO=PZ*FRACO
      PN=PZ*FRACN
      DENSEA=(0.059385*TAMB**3-8.56272*TAMB**2+65.4891*TAMB)*0.001
     ++999.84298+(GAMMA)*SALAMB
      DENSEW=(0.059385*TPLUME**3-8.56272*TPLUME**2+65.4891*TPLUME)*0.001
     ++999.84298+(GAMMA)*SALPLU
C     Previous equation re-revised to account for correct salinity units (uS/cm) in density calculations.      
      DENSEP=(1.0-VG)*DENSEW
C
C	BUBBLE PROPERTIES
      IF(RB.LE.(7.5E-4))THEN
           VB=1189.0*RB**1.1945
      ELSEIF(RB.GT.(7.5E-4).AND.RB.LT.(4.8E-3))THEN
           VB=0.22
      ELSE
           VB=2.995*RB**0.489
      ENDIF
C
      KOLO=0.6*RB
      IF(KOLO.GT.(4.0E-4))THEN
            KOLO=4.0E-4
      ENDIF
      KOLN=KOLO
C
      HO2=(2.125-0.05023*TPLUME+5.7714E-4*TPLUME**2)/100000.
      HN2=(1.042-0.02457*TPLUME+3.1714E-4*TPLUME**2)/100000.
C
      FR=V/(2.*LAMBDA*B*G*(DENSEA-DENSEP)/DENSEP)**0.5
      DCO2=HO2*PO-CO2
C      
      IF(V.GT.1.E-6)THEN
           IF(Z.GT.0.0)THEN
	        GOTO 10
	   ENDIF
      ENDIF
C
C     CALCULATION OF AVERAGE NET OXYGEN MASS TRANSFER FOR DAY
   20 GROSSMT=(FGONOT-FGO)*32./1000.*86400.
      OTEFF=(FGONOT-FGO)/FGONOT*100.
      DELTAC=COMGP-COMGNOT
C
      ELEVT=ELEV
      QWT=QW
      TPLUMET=TPLUME
      COMGPT=COMGP
      LAYTOP=LAYDIFF-JJ
      BTOP=B
      RETURN
      END

C
C------------------------------------------------------------------------------
C
      SUBROUTINE RK4_2(E,DENSEA,DENSEW,G,B,TAMB,SALAMB,GAMMA,DENSE20,
     +DOAMB,PI,V,COMGP,DNAMB,CNMGP,Y,DYDX,NN,X,H,YOUT)
      INTEGER I,NN,NMAX
      PARAMETER (NMAX=50)
      REAL*8 E,DENSEA,DENSEW,G,B,TAMB,SALAMB,GAMMA,DENSE20,DOAMB,PI,V,
     +COMGP,DNAMB,CNMGP,H,X,DYDX(NN),Y(NN),YOUT(NN),H6,HH,XH,DYM(NMAX),
     +DYT(NMAX),YT(NMAX)
      EXTERNAL DERIVS_2
      HH=H*0.5
      H6=H/6.
      XH=X+HH
      DO 11 I=1,NN
          YT(I)=Y(I)+HH*DYDX(I)
   11 CONTINUE
      CALL DERIVS_2(E,DENSEA,DENSEW,G,B,TAMB,SALAMB,GAMMA,DENSE20,DOAMB,
     +PI,V,COMGP,DNAMB,CNMGP,XH,YT,DYT)
      DO 12 I=1,NN
          YT(I)=Y(I)+HH*DYT(I)
   12 CONTINUE
      CALL DERIVS_2(E,DENSEA,DENSEW,G,B,TAMB,SALAMB,GAMMA,DENSE20,DOAMB,
     +PI,V,COMGP,DNAMB,CNMGP,XH,YT,DYM)
      DO 13 I=1,NN
          YT(I)=Y(I)+H*DYM(I)
          DYM(I)=DYT(I)+DYM(I)
   13 CONTINUE
      CALL DERIVS_2(E,DENSEA,DENSEW,G,B,TAMB,SALAMB,GAMMA,DENSE20,DOAMB,
     +     PI,V,COMGP,DNAMB,CNMGP,X+H,YT,DYT)
      DO 14 I=1,NN
          YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
   14 CONTINUE
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE DERIVS_2(E,DENSEA,DENSEW,G,B,TAMB,SALAMB,GAMMA,DENSE20,
     +           DOAMB,PI,V,COMGP,DNAMB,CNMGP,X,Y,DYDX)         
      REAL*8 E,DENSEA,DENSEW,G,B,TAMB,SALAMB,GAMMA,DENSE20,DOAMB,PI,V,
     +COMGP,DNAMB,CNMGP,X,Y(8),DYDX(8)
C     Right-hand side of differential equations for Runge-Kutta solution
      DYDX(1)=E         
      DYDX(2)=2.*PI*B**2*G*(DENSEA-DENSEW)/DENSEW
      DYDX(3)=E*TAMB
      DYDX(4)=E*(SALAMB*GAMMA/DENSE20)*DENSEA
      DYDX(5)=E*DOAMB/32.
      DYDX(6)=E*DNAMB/28.
      RETURN
      END 
