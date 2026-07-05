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
C------------------------------------------------------------------------------
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
     +ELEVT,QWT,TPLUMET,COMGPT,QWD,BWD,LAYTOP)
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
C     VAMB = Ambient Velocity (m) used to calculate VORTEX entrainment (FJR
C     PWD  = Perimeter (m) (FJR)
C     VUP  = Upward velocity (m/2)
C     AWD  = Area (m2) of plume
C     BWD  = Width 
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
     +PWD(500),BWD(500),AWD(500),
     +HWITH,HCELL,JULDAY,BTOP,BEQUIV
      INTEGER II,IJ,IK,IN,JJ,LL,NEQN,NN,MI,JK,JL,LAKE,LAYTOP
     +LAYERS,KM,KN,KO,KP,ROWS,KQ,KR,KU,KV,KW,KX,KY,KZ,KS,LAYDIFF,YEAR 
      INTEGER NELS 
C 
      PRINT*,'LINEAR_PLUME'
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

C     CALCULATION OF INITIAL WATER VELOCITY USING FROUDE NUMBER
      VBUB=4./3.*PI*RB**3
      N=QGAS/VBUB
  9   VG=QGAS/((VGUESS+VB)*(2.*LAMBDA*B)*(L-2.0*B*(1.0-LAMBDA)))
      DENSEP=(1.0-VG)*DENSEW
      IF(LAKE.EQ.1)THEN
         V=FRNOT*(2.0*LAMBDA*B*G*(DENSEA-DENSEP)/DENSEP)**0.5
      ELSEIF(LAKE.EQ.2)THEN
         V=FRNOT*(2.0*BEQUIV*G*(DENSEA-DENSEP)/DENSEP)**0.5
!         VLS: For testing purposes, assume that characteristic length is equal to rectangle width (6-4-09)
!         V=FRNOT*(2.0*LAMBDA*B*G*(DENSEA-DENSEP)/DENSEP)**0.5
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
      PWD(LAYDIFF)= 2.*(L+2.*B)  
      AWD(LAYDIFF)= L * B     
      BWD(LAYDIFF)= SQRT(AWD(LAYDIFF)/3.141592)
C      
C	SOLUTION PROCEEDURE
      DZ=0.001
      H=0.001
      NELS = 0
      HWITH=0.0
 10   Z=Z-DZ
      X=X+DZ
      ELEV=ELEV+DZ
      NELS = NELS + 1
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
      E=2.*(L+2.*B)*ALPHA*V
!     Add incremental entrainment to total cell entrainment/withdrawal       
      QWD(LAYDIFF-JJ)=QWD(LAYDIFF-JJ)+E*DZ
      PWD(LAYDIFF-JJ)=PWD(LAYDIFF-JJ)+2.*(L+2.*B)
      AWD(LAYDIFF-JJ)=AWD(LAYDIFF-JJ)+L*B*2.
      BWD(LAYDIFF-JJ)=BWD(LAYDIFF-JJ)+SQRT(L*B*2./3.141592)
      HWITH=HWITH+DZ   
      IF(HWITH.GT.HCELL)THEN
          PWD(LAYDIFF-JJ) = PWD(LAYDIFF-JJ)/NELS
          AWD(LAYDIFF-JJ) = AWD(LAYDIFF-JJ)/NELS
          BWD(LAYDIFF-JJ) = BWD(LAYDIFF-JJ)/NELS
!         PRINT *,QWD(LAYDIFF-JJ),BWD(LAYDIFF-JJ), V
          JJ=JJ+1
          HWITH=0.0
          QWD(LAYDIFF-JJ)=0.0
          PWD(LAYDIFF-JJ)=0.0
          BWD(LAYDIFF-JJ)=0.0
          AWD(LAYDIFF-JJ)=0.0
          NELS = 0
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
             
C     *******************************************************************************
	SUBROUTINE CIRCULAR_PLUME(YEAR,JULDAY,WSEL,DIFFEL,LAYERS,LAMBNOT,
     +SALAMB,PATM,DIAMM,QSCFM,FRCONOT,LAYDIFF,HCELL,LOCAT,TE,CO2M,UA,VA,
     +ELEVT,QWT,TPLUMET,COMGPT,QWD,BWD,LAYTOP)
C     *******************************************************************************
C	THIS SUBROUTINE IS WRITTEN TO PREDICT THE PERFORMANCE OF A CIRCULAR BUBBLE PLUME.  
C     THE MODEL IS BASED ON THE WUEST ET AL. (1992) CIRCULAR BUBBLE PLUME MODEL.
C     By:
C	
C	VERSION 3 (to couple with Francisco Rueda's reservoir model and for Lake Hallwil real diffuser)
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
!     BIFF=DIFFUSER RADIUS (m) 
C     BAVG=AVERAGE 1/2 DIFFUSER WIDTH (m)
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
!	FR=FROUDE NUMBER (-)
C     FRCNOT=FRACTION OF NITROGEN IN ATMOSPHERE (-)
C     FRNOT=INITIAL FROUDE NUMBER (-)
C     GAMMA=SALINITY CONVERSION FACTOR [kg/m3/(uS/cm)]
C     GROSSMT=GROSS MASS TRANSFER OF OXYGEN FROM PLUME (kg/d)
!	HO2=SOLUBILITY CONSTANT FOR OXYGEN (mol/m3/Pa)
!	HN2=SOLUBILITY CONSTANT FOR NITROGEN (mol/m3/Pa)
C     HCELL=HEIGHT OF CELL IN GRID (m)
C	HO=SOLUBILITY CONSTANT FOR OXYGEN (mol/m3/Pa)
C	HN=SOLUBILITY CONSTANT FOR NITROGEN (mol/m3/Pa)
C     HWITH=HEIGHT OF WITHDRAWAL/ENTRAINMENT ZONE (m)
C     JULDAY=JULIAN DAY IN GIVEN YEAR
C	KOLO=MASS TRANSFER COEFFICIENT FOR OXYGEN (m/s)
C	KOLN=MASS TRANSFER COEFFICIENT FOR NITROGEN (m/s)
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
!     QGAS=TOTAL GAS FLOW TO DIFFUSER
C	QSCFM=STANDARD GAS FLOW RATE (scfm), TOTAL STANDAR GAS FLOW RATE TO DIFFUSER 
C	QSCMS=STANDARD GAS FLOW RATE (scms)
C	QW=FLOWRATE OF WATER (m3/s)
C     QWT=TOTAL DETRAINMENT FLOW RATE AT THE TOP OF THE PLUME (m3/s)  
C	RB=BUBBLE RADIUS (m)
C	RGAS=IDEAL GAS CONSTANT (J/mol/K)
C     SAL=SALINITY (uS/cm)
C     SALAMB=SALINITY OF AMBIENT WATER (uS/cm)
C	SALPLU=SALINITY OF THE PLUME (uS/cm)
!     SHEARE=SHEAR ENTRAINMENT
C	TAMB=AMBIENT WATER TEMPERATURE (C)
C     TAVG=AVERAGE AMBIENT WATER TEMPERATURE (C)
C     TDS=TOTAL DISSOLVED SOLIDS (g/m3) [0.64 conversion factor from Chapra book]
C     TE=TEMPERATURE PROFILE FOR INPUT BOUNDARY CONDITION (C) 
C	TPLUME=PLUME TEMPERATURE (C)
C     TPLUMET=DETRAINMENT PLUME TEMPERATURE AT THE TOP OF THE PLUME (C) 
C	TSTD=STANDARD TEMPERATURE (K)
!     UA = U AMBIENT VELOCITY PROFILE FOR INPUT BOUNDARY CONDITION
C     UAMB = U Ambient Velocity (m) used to calculate VORTEX entrainment (FJR
C	V=WATER VELOCITY (m/s)
!     VA = V AMBIENT VELOCITY PROFILE FOR INPUT BOUNDARY CONDITION
C     VAMB = V Ambient Velocity (m) used to calculate VORTEX entrainment (FJR
C     VAVG=AVERAGE WATER VELOCITY (m/s)
C	VB=BUBBLE RISE VELOCITY (m/s)
C	VBUB=BUBBLE VOLUME (m3)
!     VG=GAS VOLUME PER TOTAL VOLUME OF THE BUBBLE-WATER MIXTURE IN THE INNER CORE OF THE PLUME
C	VGUESS=GUESSED INITIAL WATER VELOCITY (m/s)
!     VORTEXE=VORTEX ENTREAINMENT
C     WSEL=WATER SURFACE ELEVATION (m)
C	YO2=GASEOUS OXYGEN CONCENTRATION (mol/m3)
C	YN2=GASEOUS NITROGEN CONCENTRATION (mol/m3)
C	Z=DEPTH TO DIFFUSER (m)
C     PWD  = Perimeter (m) (FJR)
C     VUP  = Upward velocity (m/2)
C     AWD  = Area (m2) of plume
C     BWD  = Width 
C
      REAL*8 ALPHA,AREA,B,CO2,COMG,CN2,CNMG,DS,DZ,DENSEA,DENSEP,DENSEW,
     +DIAMM,E,FDO,FDN,FRACO,FRACN,FSAL,FTEMP,FGO,FGN,G,GAMMA,HO2,HN2,
     +KOLN,KOLO,LAMBDA,MOMENT,N,PI,PO,PN,PSTD,PZ,QSCFM,QSCMS,QW,QGAS,
     +RB,RGAS,SALAMB,SALPLU,TAMB,TPLUME,TSTD,V,VB,VBUB,VG,VGUESS,
     +YO2,YN2,Z,AA,BB,CC,BNOT,LNOT,ELEV,DT,TE(1000),XLOC,DEPTH,
     +CO2M(1000),LOCAT(1000),DOAMB,COMGP,CNMGP,WSEL,PATM,SAL(1000),
     +DIFFEL,GROSSMT,FGONOT,DNAMB,DMPR,TAVG,SUMTEMP,SUMSAL,
     +LAMBNOT,FRCONOT,H,DYDX(8),Y(8),YOUT(8),
     +DENSE20,OTEFF,FRNOT,VDIFF,FR,BUOY,DCO2,QGFRAC,LDIFF,TDS(1000),
     +JDAY,EL(70),DELTAC,COMGNOT,ELEVT,QWT,TPLUMET,COMGPT,QWD(500),
     +PWD(500),BWD(500),AWD(500),X,SHEARE,VORTEXE,
     +HWITH,HCELL,JULDAY,BTOP,UA(1000),VA(1000),UAMB,VAMB
      INTEGER II,IJ,IK,IN,JJ,LL,NEQN,NN,MI,JK,JL,LAYTOP, 
     +LAYERS,KM,KN,KO,KP,ROWS,KQ,KR,KU,KV,KW,KX,KY,KZ,KS,LAYDIFF,YEAR 
      INTEGER NELS,ierror,M,CONT
C 
C     --------------------------------------------------------------------------
C     CONSTANTS
C     --------------------------------------------------------------------------
!
       PRINT*,'CIRCULAR_PLUME'
!      PRINT*,'Yearplume', YEAR
!      PRINT*,'Julday', JULDAY
!      PRINT*,'WSEL', WSEL
!      PRINT*,'DIFFEL', DIFFEL
!      PRINT*,'LAYERS', LAYERS
!      PRINT*,'LAMBNOT', LAMBNOT
!      PRINT*,'SALAMB', SALAMB
!      PRINT*,'PATM', PATM
!      PRINT*,'DIAMM',DIAMM
!      PRINT*,'QSCFM',QSCFM
!      PRINT*,'FRCONOT',FRCONOT
!      PRINT*,'LAYDIFF',LAYDIFF
!      PRINT*,'HCELL',HCELL
!      PRINT*,'LOCAT',LOCAT
!      PRINT*,'TE',TE
!      PRINT*,'CO2M',CO2M
!      PRINT*,'UA',UA
!      PRINT*,'VA',VA
!      PRINT*,'ELEVT',ELEVT
!      PRINT*,'QWT',QWT
!      PRINT*,'TPLUMET',TPLUMET
!      PRINT*,'COMGPT',COMGPT
!      PRINT*,'QWD',QWD
!      PRINT*,'BWD',BWD
!      PRINT*,'LAYTOP',LAYTOP

      ALPHA=0.111
      G=9.80665
      GAMMA=6.9E-4
      LAMBDA=0.8 !McGinnis
      FRNOT=1.6  !McGinnis
      PI=ACOS(-1.0)
      PSTD=101325.
      RGAS=8.314
      TSTD=293.15
      DENSE20=998.2
      FRCNATM=0.79
      QGFRAC=1.0
!     
C     --------------------------------------------------------------------------
C     PARAMETERS AND INITIALIZE THE VARIABLES 
C     --------------------------------------------------------------------------
C      
!     Calculate the initial average water temperature  
      SUMTEMP=0.0
      DO KS=1,LAYERS
          SUMTEMP=SUMTEMP+TE(KS)
      END DO       
      TAVG=SUMTEMP/LAYERS     
   
!     Geometric caracteristic 
      DEPTH=WSEL-DIFFEL   
      Z=DEPTH
      ELEV=DIFFEL
!        
C     Assume that gas bubbles are composed of oxygen and nitrogen only.
      FRACO=FRCONOT 
      FRACN=1.0-FRACO
!
C     Interpolate input profiles to obtain the plume initial conditions
      X=0.
      XLOC=DEPTH-X
      CALL LININT(LOCAT,TE,LAYERS,XLOC,TAMB)
      CALL LININT(LOCAT,CO2M,LAYERS,XLOC,COMG)
      CALL LININT(LOCAT,UA,LAYERS,XLOC,UAMB)
      CALL LININT(LOCAT,VA,LAYERS,XLOC,VAMB)
      COMGP=COMG
      DOAMB=COMG
      CO2=COMG/32.
      COMGNOT=COMG
!     At first the temperature and the salinity un the plume is the same that in the ambient water
      SALPLU=SALAMB    
      TPLUME=TAMB
!     Guess the initial velocity 
      VGUESS=0.07
      V=VGUESS
!     INITIAL DIFFUSER SIZE
!     At first the plume radius BNOT (initial plume radius) in the TOP HAT models 
      BNOT=LAMBNOT/LAMBDA
      B=BNOT
!     BDIFF=LAMBNOT
C      
C     AMBIENT AND AVERAGE WATER DENSITIES
      DENSEA=(0.059385*TAMB**3-8.56272*TAMB**2+65.4891*TAMB)*0.001
     ++999.84298+(GAMMA)*SALAMB
      DENSEW=DENSEA
C
!     Solubility constant (mol/m3/Pa)
      HO2=(2.125-0.05023*TPLUME+5.7714E-4*TPLUME**2)/100000.
      HN2=(1.042-0.02457*TPLUME+3.1714E-4*TPLUME**2)/100000.
C	
C     Assume initial ambient dissolved nitrogen conc. equals saturated conc. at surface.
      CN2=(PATM*FRCNATM)*HN2      
      CNMG=CN2*28.0
      CNMGP=CNMG
      DNAMB=CNMG
!
C     --------------------------------------------------------------------------
C     BUBBLE PROPERTIES
C     --------------------------------------------------------------------------
!
      QSCMS=QGFRAC*QSCFM/3.281**3/60.0
      QGAS=PSTD*QSCMS*(TAMB+273.15)/((PATM+DENSEA*G*Z)*TSTD)
C     Initial bubble size.     
      RB=DIAMM/2000.
!     Bubble rise velocity
      IF(RB.LE.(7.5E-4))THEN
            VB=1189.0*RB**1.1945
      ELSEIF(RB.GT.(7.5E-4).AND.RB.LT.(4.8E-3))THEN
            VB=0.22
      ELSE
            VB=2.995*RB**0.489
      ENDIF
!     Mass transfer coefficient
      KOLO=0.6*RB
      IF(KOLO.GT.(4.0E-4))THEN
            KOLO=4.0E-4
      ENDIF
      KOLN=KOLO
!
C     CALCULATION OF INITIAL WATER VELOCITY USING FROUDE NUMBER
      VBUB=4./3.*PI*RB**3
      N=QGAS/VBUB
      VDIFF=1
      DO WHILE (VDIFF.GT.1.0E-6)
         VG=QGAS/((VGUESS+VB)*(PI*(LAMBDA*B)**2))
         DENSEP=(1.0-VG)*DENSEW
         V=FRNOT*(2.0*LAMBDA*B*G*(DENSEA-DENSEP)/DENSEP)**0.5
         VDIFF=ABS(V-VGUESS)
         VGUESS=V
      END DO  
C     ------------------------------------------------------------------
C     VARIABLE TRANSFORMATION      
C     ------------------------------------------------------------------
      E=(2.*PI*B)*ALPHA*V    
      QW=V*PI*B**2
      MOMENT=(PI*B**2)*V**2
      FTEMP=QW*TPLUME 
      FSAL=QW*(SALPLU*GAMMA/DENSE20)*DENSEW
C     Previous equation corrected to account for salinity units conversion.        
      FDO=QW*CO2
      FDN=QW*CN2
      FGO=PSTD*QSCMS/(RGAS*TSTD)*FRACO
      FGONOT=FGO
      FGN=PSTD*QSCMS/(RGAS*TSTD)*FRACN
C     Revised gaseous flux equations.
      YO2=FGO/((PI*(LAMBDA*B)**2)*(V+VB))
      YN2=FGN/((PI*(LAMBDA*B)**2)*(V+VB))           
      PZ=PATM+(DENSEA*G*Z)
      PO=PZ*FRACO
      PN=PZ*FRACN
      BUOY=(G*(DENSEA-DENSEP)/DENSEP*QW)/LNOT
      TDS=SALPLU*0.64
!     Initialize lateral withdrawal flowrate for first/lowest cell in column/segment
      JJ=0
      QWD(LAYDIFF)= QW      
      PWD(LAYDIFF)= (2*PI*B)
      AWD(LAYDIFF)= (PI*B**2)    
      BWD(LAYDIFF)= SQRT(AWD(LAYDIFF)/PI)       
!     
!     -------------------------------------------------------------------------- 
C	SOLUTION PROCEEDURE
!     --------------------------------------------------------------------------
!
      DZ=0.001
      H=0.001
      NELS = 0
      HWITH=0.0
!
!     OPEN (UNIT=50, FILE="salida.txt", IOSTAT=ierror)
!     WRITE (UNIT=50, FMT='(A)') "   "
!
      M=0
      DO WHILE (V.GT.1.E-6.AND.Z.GT.0.0)
         Z=Z-DZ
         X=X+DZ
         ELEV=ELEV+DZ
         NELS = NELS + 1
C	
C        Interpolate input profiles to obtain line plume boundary conditions
         XLOC=DEPTH-X
         CALL LININT(LOCAT,CO2M,LAYERS,XLOC,COMG)
         DOAMB=COMG
         CALL LININT(LOCAT,TE,LAYERS,XLOC,TAMB)
         CALL LININT(LOCAT,UA,LAYERS,XLOC,UAMB)
         CALL LININT(LOCAT,VA,LAYERS,XLOC,VAMB)

C                  
C        Use subroutines for Runge Kutta method solution
         NEQN=8
         Y(1)=QW
         Y(2)=MOMENT
         Y(3)=FTEMP
         Y(4)=FSAL
         Y(5)=FDO
         Y(6)=FDN
         Y(7)=FGO
         Y(8)=FGN
!      PRINT *, E
!      PRINT *, DENSEA
!      PRINT *, DENSEW
!      PRINT *, DENSEP
!      PRINT *, G
!      PRINT *, B
!      PRINT *, LAMBDA
!      PRINT *, TAMB
!      PRINT *, SALAMB
!      PRINT *, GAMMA
!      PRINT *, DENSE20
!      PRINT *, DOAMB
!      PRINT *, PI
!      PRINT *, RB
!      PRINT *, N
!      PRINT *, V
!      PRINT *, VB
!      PRINT *, KOLO
!      PRINT *, HO2
!      PRINT *, PO
!      PRINT *, COMGP
!      PRINT *, DNAMB
!      PRINT *, KOLN
!      PRINT *, HN2
!      PRINT *, PN
!      PRINT *, CNMGP
!      PRINT *, Z
!      PRINT *, Y
!      PRINT *, DYDX
         CALL DERIVS_3(E,DENSEA,DENSEW,DENSEP,G,B,LAMBDA,TAMB,SALAMB,
     +            GAMMA,DENSE20,DOAMB,PI,RB,N,V,VB,KOLO,HO2,PO,
     +            COMGP,DNAMB,KOLN,HN2,PN,CNMGP,Z,Y,DYDX)
         CALL RK4_3(E,DENSEA,DENSEW,DENSEP,G,B,LAMBDA,TAMB,SALAMB,
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
C        Previous equation corrected to consistently express salinity in uS/cm	   
	    CO2=FDO/QW
	    CN2=FDN/QW
!      PRINT*,"MOMENT",MOMENT
!      PRINT*,"TPLUME",TPLUME
!      PRINT*,"FSAL",FSAL
!      PRINT*,"QW",QW
!      PRINT*,"DENSEW",DENSEW
!      PRINT*,"GAMMA",GAMMA
!      PRINT*,"DENSE20",DENSE20
!      PRINT*,"SALPLU",SALPLU
!      PRINT*,"CO2",CO2
!      PRINT*,"CN2",CN2
	    GOTO 20
         ENDIF
         V=MOMENT/QW
         AREA=QW/V
         B=SQRT(AREA/PI)
!        TOTAL ENTRAIMENT: Shear and Vortex
!         E=2.*PI*B*ALPHA*V
        SHEARE=2.*PI*B*ALPHA*V
!        Vortex Entraiment: Hypothesis Projected Area Entraiment
         VORTEXE=(2*B*ABS(UAMB)+2*B*ABS(VAMB))
!	 a) Additive hypothesis 
!      	     E=SHEARE+VORTEXE
!	 b) Maximun hypothesis
	     IF (SHEARE.GT.VORTEXE) THEN
	 	E=SHEARE
             ELSE
		E=VORTEXE
	     ENDIF

         M=M+1
         IF (M.EQ.500) THEN
!             PRINT*, XLOC,B,V,UAMB,VAMB,SHEARE,VORTEXE
             WRITE (UNIT=50, FMT='(7F12.8)') XLOC,B,V,UAMB,VAMB,SHEARE,VORTEXE
             M=0
         ENDIF
!
!        Add incremental entrainment to total cell entrainment/withdrawal       
         QWD(LAYDIFF-JJ)=QWD(LAYDIFF-JJ)+E*DZ
         PWD(LAYDIFF-JJ)=PWD(LAYDIFF-JJ)+2.*PI*B
         AWD(LAYDIFF-JJ)=AWD(LAYDIFF-JJ)+PI*B**2
         BWD(LAYDIFF-JJ)=BWD(LAYDIFF-JJ)+SQRT((PI*B**2)/PI)
         HWITH=HWITH+DZ   
         IF(HWITH.GT.HCELL)THEN
            PWD(LAYDIFF-JJ) = PWD(LAYDIFF-JJ)/NELS
            AWD(LAYDIFF-JJ) = AWD(LAYDIFF-JJ)/NELS
            BWD(LAYDIFF-JJ) = BWD(LAYDIFF-JJ)/NELS
!            PRINT *,"QWD(LAYDIFF-JJ)",QWD(LAYDIFF-JJ)
!            PRINT *,"BWD(LAYDIFF-JJ)",BWD(LAYDIFF-JJ) 
!            PRINT *,"V",V
            JJ=JJ+1
            HWITH=0.0
            QWD(LAYDIFF-JJ)=0.0
            PWD(LAYDIFF-JJ)=0.0
            BWD(LAYDIFF-JJ)=0.0
            AWD(LAYDIFF-JJ)=0.0
            NELS = 0
         ENDIF
!        Temperatura and salinity in the plume    
         TPLUME=FTEMP/QW
         SALPLU=FSAL/(QW*DENSEW)/(GAMMA/DENSE20)
!      PRINT*,"2222222",2
!      PRINT*,"FSAL",FSAL
!      PRINT*,"QW",QW
!      PRINT*,"DENSEW",DENSEW
!      PRINT*,"GAMMA",GAMMA
!      PRINT*,"DENSE20",DENSE20
!      PRINT*,"SALPLU",SALPLU
C        Previous equation corrected to consistently express salinity in uS/cm
!        Dissolved oxygen and nitrogen concentration
         CO2=FDO/QW
         CN2=FDN/QW
         COMGP=CO2*32.
         CNMGP=CN2*28.
C        Revised gaseous flux equations.
         YO2=FGO/((PI*(LAMBDA*B)**2)*(V+VB))
         YN2=FGN/((PI*(LAMBDA*B)**2)*(V+VB))
C  
!        
         PZ=PATM+(DENSEA*G*Z)
         QGAS=(FGO+FGN)*RGAS*(TPLUME+273.15)/PZ
         VBUB=QGAS/N
C
!        GAS VOLUME PER TOTAL VOLUME OF THE BUBBLE-WATER MIXTURE IN THE INNER CORE OF THE PLUME
         VG=VBUB*N/((V+VB)*(PI*(LAMBDA*B)**2))
C        Previous equation revised to account for correct plume cross-sectional area occupied by bubbles. 
!        Bubbles radius     
         RB=(3.*QGAS/(4.*PI*N))**(1./3.)
         IF(RB.LT.0.0)THEN
            RB=1.0E-8
         ENDIF
         FRACO=FGO/(FGO+FGN)
         FRACN=1.0-FRACO
C	
!        Partial Pressure
         PO=PZ*FRACO
         PN=PZ*FRACN
!        Density (ambient water, water in plume and of the plume)
      DENSEA=(0.059385*TAMB**3-8.56272*TAMB**2+65.4891*TAMB)*0.001
     ++999.84298+(GAMMA)*SALAMB
      DENSEW=(0.059385*TPLUME**3-8.56272*TPLUME**2+65.4891*TPLUME)*0.001
     ++999.84298+(GAMMA)*SALPLU
C        Previous equation re-revised to account for correct salinity units (uS/cm) in density calculations.      
         DENSEP=(1.0-VG)*DENSEW
C
C        BUBBLE PROPERTIES
!
!        Bubble rise velocity
         IF(RB.LE.(7.5E-4))THEN
            VB=1189.0*RB**1.1945
         ELSEIF(RB.GT.(7.5E-4).AND.RB.LT.(4.8E-3))THEN
            VB=0.22
         ELSE
            VB=2.995*RB**0.489
         ENDIF
C
!        Mass transfer coeffitient
         KOLO=0.6*RB
         IF(KOLO.GT.(4.0E-4))THEN
            KOLO=4.0E-4
         ENDIF
         KOLN=KOLO
C
!        Solubility constant (mol/m3/Pa)
         HO2=(2.125-0.05023*TPLUME+5.7714E-4*TPLUME**2)/100000.
         HN2=(1.042-0.02457*TPLUME+3.1714E-4*TPLUME**2)/100000.
C
!        Froude number
         FR=V/(2.*LAMBDA*B*G*(DENSEA-DENSEP)/DENSEP)**0.5
!
         DCO2=HO2*PO-CO2
C      
      END DO
C
!     --------------------------------------------------------------------
C     CALCULATION OF AVERAGE NET OXYGEN MASS TRANSFER FOR DAY
!     --------------------------------------------------------------------
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
!      PRINT*,"ELEVT",ELEVT
!      PRINT*,"QWT",QWT
!      PRINT*,"TPLUMET",TPLUMET
!      PRINT*,"COMGPT",COMGPT
!      PRINT*,"LAYTOP",LAYTOP
!      PRINT*,"BTOP",BTOP
!      CLOSE (UNIT=50)
      RETURN
      END


C
C------------------------------------------------------------------------------
C                          
             
C     *******************************************************************************
	SUBROUTINE INNER_PLUME(YEAR,JULDAY,WSEL,DIFFEL,LAYERS,LAMBNOT,
     +SALAMB,PATM,DIAMM,QSCFM,FRCONOT,LAYDIFF,HCELL,LOCAT,TE,CO2M,UA,VA,
     +ELEVT,QWT,TPLUMET,COMGPT,SALPLUMET,QWDI,
     +BWDI,LAYTOP,LAYINTR,DEPTHINTR,NLI,NLO,NITERPLUME,
     +INPLUME,OUTPLUME,ALPHAI,ALPHAO,ALPHAA,LAMBDA,FRNI,FRNO,GAMMA1,
     +OTEFF)
C     *******************************************************************************

      IMPLICIT NONE

C	THIS SUBROUTINE IS WRITTEN TO PREDICT THE PERFORMANCE OF A CIRCULAR BUBBLE PLUME.  
C     THE MODEL IS BASED ON THE WUEST ET AL. (1992) CIRCULAR BUBBLE PLUME MODEL.
C     By:
C	     
C	VERSION 3 (to couple with Francisco Rueda's reservoir model and for Lake Hallwil real diffuser)
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
C	  ALPHAI=ENTRAINMENT COEFFICIENT INNER PLUME  (-)
C	  ALPHAO=ENTRAINMENT COEFFICIENT OUTER PLUME  (-)
C	  ALPHAA=ENTRAINMENT COEFFICIENT FROM AMBIENT (-)
C     BI=1/2 INNER PLUME WIDTH (m)
!     BIFF=DIFFUSER RADIUS (m) 
C     BAVG=AVERAGE 1/2 DIFFUSER WIDTH (m)
!         C1=COEFFICIENT OF DIFERENT COUNTERFLOW ENTRAINMENTS
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
!         DEPTHINTR=INTRUSION DEPTH OF THE DOUBLE PLUME
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
!	FR=FROUDE NUMBER (-)
C     FRCNOT=FRACTION OF NITROGEN IN ATMOSPHERE (-)
C     FRNO=INITIAL FROUDE NUMBER (-)
C     GAMMA=SALINITY CONVERSION FACTOR [kg/m3/(uS/cm)]
C     GROSSMT=GROSS MASS TRANSFER OF OXYGEN FROM PLUME (kg/d)
!	HO2=SOLUBILITY CONSTANT FOR OXYGEN (mol/m3/Pa)
!	HN2=SOLUBILITY CONSTANT FOR NITROGEN (mol/m3/Pa)
C     HCELL=HEIGHT OF CELL IN GRID (m)
C	HO=SOLUBILITY CONSTANT FOR OXYGEN (mol/m3/Pa)
C	HN=SOLUBILITY CONSTANT FOR NITROGEN (mol/m3/Pa)
C     HWITH=HEIGHT OF WITHDRAWAL/ENTRAINMENT ZONE (m)
!         INPLUME(:,6): INNER PLUME CHARACTERISTIC [HEIGTH; FLOWRATE OF WATER; MOMENTUM; TEMPERATURE; SALINITY; O2] JCT
C     JULDAY=JULIAN DAY IN GIVEN YEAR
C	KOLO=MASS TRANSFER COEFFICIENT FOR OXYGEN (m/s)
C	KOLN=MASS TRANSFER COEFFICIENT FOR NITROGEN (m/s)
C     LAKE=LAKE AND DIFFUSER TYPE (SHR AND LINEAR=1 OR AMISK AND RECTANGULAR=2) FOR SELECTION OF LAMBDA
C     LAYERS=NUMBER OF LAYERS/DATA POINTS IN BOUNDARY CONDITION PROFILES (-)
C     LAYDIFF=GRID LAYER CORRESPONDING TO DIFFUSER DEPTH (-)
!         LAYINTR=GRID LAYER CORRESPONDING TO INTRUSION DEPTH OF THE DOUBLE PLUME (-)
C     LAYTOP=GRID LAYER CORRESPONDING TO TOP OF PLUME (-)
C     LDIFF=LENGTH OF DIFFUSER (m)
C     LNOT=INITIAL DIFFUSER LENGTH (m)
C	LAMBDA=FRACTION OF PLUME OCCUPIED BY BUBBLES (-)
C     LAMBNOT=LAMBDA x INITIAL PLUME RADIUS;EQUAL TO DIFFUSER RADIUS (m) 
C     LOCAT=DEPTHS FOR INPUT BOUNDARY CONDITION PROFILES (m) 
C	MOMENT=MOMENTUM (m4/s)
C	N=NUMBER OF BUBBLES PER SECOND (1/s)
!         NITERPLUME= NUMBER OF ITERATION
!         NLI
!         NLO
C     OTEFF=OXYGEN TRANSFER EFFICEINCY (%)
!         OUTPLUME(:,6)= OUTER PLUME CHARACTERISTIC [HEIGTH; FLOWRATE OF WATER; MOMENTUM; TEMPERATURE; SALINITY; O2] JCT
!         OUTERMOMENT= MOMENTUM IN THE OUTER PLUME 
!         OUTERQW=FLOWRATE OF WATER IN THE OUTER PLUME
!         OUTXLOC
!         OUTQW
!         OUTXMOMENT
!         OUTXTEMP
!         OUTO2  
C     PATM=ATMOSPHERIC PRESSURE AT AVERAGE WSEL (Pa)
C	PSTD=STANDARD PRESSURE (Pa)
!     QGAS=TOTAL GAS FLOW TO DIFFUSER
C	QSCFM=STANDARD GAS FLOW RATE (scfm), TOTAL STANDAR GAS FLOW RATE TO DIFFUSER 
C	QSCMS=STANDARD GAS FLOW RATE (scms)
C	QW=FLOWRATE OF WATER (m3/s)
C     QWT=TOTAL DETRAINMENT FLOW RATE AT THE TOP OF THE PLUME (m3/s)  
C	RB=BUBBLE RADIUS (m)
C	RGAS=IDEAL GAS CONSTANT (J/mol/K)
C     SAL=SALINITY (uS/cm)
C     SALAMB=SALINITY OF AMBIENT WATER (uS/cm)
C	SALPLU=SALINITY OF THE PLUME (uS/cm)
!         SALPLUMET=DETRAINMENT PLUME SALINITY AT THE TOP OF THE PLUME (uS/cm)
!     SHEARE=SHEAR ENTRAINMENT
C	TAMB=AMBIENT WATER TEMPERATURE (C)
C	  TARO=AROUND WATER TEMPERATURE (C)
C     TAVG=AVERAGE AMBIENT WATER TEMPERATURE (C)
C     TDS=TOTAL DISSOLVED SOLIDS (g/m3) [0.64 conversion factor from Chapra book]
C     TE=TEMPERATURE PROFILE FOR INPUT BOUNDARY CONDITION (C) 
C	TPLUME=PLUME TEMPERATURE (C)
C     TPLUMET=DETRAINMENT PLUME TEMPERATURE AT THE TOP OF THE PLUME (C) 
C	TSTD=STANDARD TEMPERATURE (K)
!     UA = U AMBIENT VELOCITY PROFILE FOR INPUT BOUNDARY CONDITION
C     UAMB = U Ambient Velocity (m) used to calculate VORTEX entrainment (FJR
C	VI= INNER PLUME WATER VELOCITY (m/s)
!     VA = V AMBIENT VELOCITY PROFILE FOR INPUT BOUNDARY CONDITION
C     VAMB = V Ambient Velocity (m) used to calculate VORTEX entrainment (FJR
C     VAVG=AVERAGE WATER VELOCITY (m/s)
C	VB=BUBBLE RISE VELOCITY (m/s)
C	VBUB=BUBBLE VOLUME (m3)
!     VG=GAS VOLUME PER TOTAL VOLUME OF THE BUBBLE-WATER MIXTURE IN THE INNER CORE OF THE PLUME
C	VGUESS=GUESSED INITIAL WATER VELOCITY (m/s)
!         VO=WATER VELOCITY IN THE OUTER PLUME (m/s)
!     VORTEXE=VORTEX ENTREAINMENT
C     WSEL=WATER SURFACE ELEVATION (m)
C	YO2=GASEOUS OXYGEN CONCENTRATION (mol/m3)
C	YN2=GASEOUS NITROGEN CONCENTRATION (mol/m3)
C	Z=DEPTH TO DIFFUSER (m)
C     PWD  = Perimeter (m) (FJR)
C     VUP  = Upward velocity (m/2)
C     AWD  = Area (m2) of plume
C     BWD  = Width 

!    GAMMA1= momentum amplification factor
!  DOARO
C
      REAL*8 AREA,BI,CO2,COMG,CN2,CNMG,DS,DZ,DENSEA,DENSEP,DENSEW,
     +DIAMM,FDO,FDN,FRACO,FRACN,FSAL,FTEMP,FGO,FGN,G,GAMMA,HO2,HN2,
     +KOLN,KOLO,LAMBDA,MOMENT,N,PI,PO,PN,PSTD,PZ,QSCFM,QSCMS,QW,QGAS,
     +RB,RGAS,SALAMB,SALPLU,TARO,TPLUME,TSTD,VI,VB,VBUB,VG,VGUESS,
     +YO2,YN2,Z,AA,BB,CC,BNOT,LNOT,ELEV,DT,TE(1000),XLOC,DEPTH,
     +CO2M(1000),LOCAT(1000),DOAMB,COMGP,CNMGP,WSEL,PATM,SAL(1000),
     +DIFFEL,GROSSMT,FGONOT,DNAMB,DMPR,SUMSAL,OTEFF,
     +LAMBNOT,FRCONOT,H,DYDX(8),Y(8),YOUT(8),FRCNATM,
     +DENSE20,FRNI,FRNO,VDIFF,FR,BUOY,DCO2,QGFRAC,LDIFF,TDS(1000),
     +JDAY,EL(70),DELTAC,COMGNOT,ELEVT,QWT,TPLUMET,COMGPT,
     +PWDI(500),BWDI(500),AWDI(500),QWDI(500),X,SHEARE,VORTEXE,
     +HWITH,HCELL,JULDAY,BTOP,UA(1000),VA(1000),UAMB,VAMB,
     +INPLUME(35000,9),OUTPLUME(35000,9),SALARO,DOARO,AU,
     +DEPTHINTR,ALPHAI,ALPHAO,ALPHAA,SALPLUMET,GAMMA1,TAMB,
     +C1,VO,EI,EO,OUTERQW,OUTERMOMENT,OUTXLOC(35000),OUTQW(35000),
     +OUTMOMENT(35000),OUTTEMP(35000),OUTO2(35000),V,OUTSAL(35000)
      INTEGER II,IJ,IK,IN,JJ,LL,NEQN,NN,MI,JK,JL,LAYTOP,
     +LAYERS,KM,KN,KO,KP,ROWS,KQ,KR,KU,KV,KW,KX,KY,KZ,KS,LAYDIFF,YEAR 
      INTEGER NELS,ierror,M,CONT,
     +LAYINTR,
     +NITERPLUME,NLI,NLO,NLIO,k
C 
C     --------------------------------------------------------------------------
C     CONSTANTS
C     --------------------------------------------------------------------------
!

      PRINT*,'INNER_PLUME'
      !PRINT*, 'INNER_PLUME', ALPHAI,ALPHAO,ALPHAA,LAMBDA,FRNI,FRNO,GAMMA1

       !IF (NITERPLUME == 2) THEN
       !  OPEN (UNIT=55, FILE="radiocaudal.txt", POSITION="APPEND")
       !    DO k = 1,35000
       !      !WRITE (UNIT=55, FMT='(7F12.6)') OUTPLUME(k,1),OUTPLUME(k,2),OUTPLUME(k,3),OUTPLUME(k,4),OUTPLUME(k,5),OUTPLUME(k,6),OUTPLUME(k,7)
       !      WRITE (UNIT=55, FMT='(2F12.6)') OUTPLUME(k,1),OUTPLUME(k,3)
       !   ENDDO
       !    DO k = 1,35000
       !      !WRITE (UNIT=55, FMT='(7F12.6)') INPLUME(k,1),INPLUME(k,2),INPLUME(k,3),INPLUME(k,4),INPLUME(k,5),INPLUME(k,6),INPLUME(k,7)
       !      WRITE (UNIT=55, FMT='(2F12.6)') INPLUME(k,1),INPLUME(k,3)
       !   ENDDO
       !  CLOSE (UNIT=55)
       !ENDIF 



            ! WRITE (UNIT=55, FMT='(12F12.6)') OUTPLUME(k,1),             &
            !                 OUTPLUME(k,2),OUTPLUME(k,3),OUTPLUME(k,4),  &
            !                 OUTPLUME(k,5),OUTPLUME(k,6),OUTPLUME(k,7),  &
            !                 INPLUME(k,1),INPLUME(k,2),                  &
            !                 INPLUME(k,3),INPLUME(k,4),INPLUME(k,5),     &
            !                 INPLUME(k,6),INPLUME(k,7)
            ! WRITE (UNIT=50, FMT='(7F12.6)') XLOC,B,V,UAMB,VAMB,SHEARE,VORTEXE

    !  PRINT*,'Yearplume', YEAR
    !  PRINT*,'Julday', JULDAY
    !  PRINT*,'WSEL', WSEL
    !  PRINT*,'DIFFEL', DIFFEL
    !  PRINT*,'LAYERS', LAYERS
    !  PRINT*,'LAMBNOT', LAMBNOT
    !  PRINT*,'SALAMB', SALAMB
    !  PRINT*,'PATM', PATM
    !  PRINT*,'DIAMM', DIAMM
    !  PRINT*,'QSCFM', QSCFM
    !  PRINT*,'FRCONOT',FRCONOT
    !  PRINT*,'LAYDIFF',LAYDIFF
    !  PRINT*,'HCELL',HCELL
    !   PRINT*,'LOCAT',LOCAT
    !   PRINT*,'TE',TE
    !   PRINT*,'CO2M',CO2M
   ! PRINT*,'UA',UA
   ! PRINT*,'VA',VA
    !  PRINT*,'ELEVT',ELEVT
    !  PRINT*,'QWT',QWT
    !  PRINT*,'TPLUMET',TPLUMET
    !  PRINT*,'COMGPT',COMGPT
    !  PRINT*,'SALPLUMET',SALPLUMET



    !  PRINT*,'LAYTOP',LAYTOP
    !  PRINT*,'LAYINTR',LAYINTR
    !  PRINT*,'DEPTHINTR',DEPTHINTR
    !  PRINT*,'NLI',NLI
    !  PRINT*,'NLO',NLO
    !  PRINT*,'NITERPLUMEinINNERPLUME', NITERPLUME
    !  PRINT*,'INPLUME',INPLUME
    !  PRINT*,'OUTPLUME',OUTPLUME

!      NITERPLUME=2
!      PRINT*,'NITERPLUME2', NITERPLUME

    
      !ALPHAI=0.111 ! Crounse et al 2007
      !ALPHAI=0.055 ! Socolofsky et al. 2008
      !ALPHAO=0.11  ! Socolofsky et al. 2008
      !ALPHAA=0.11  ! Socolofsky et al. 2008
      G=9.80665
      GAMMA=6.9E-4
      !GAMMA1=1.1 ! JCT Socolofsky
!      LAMBDA=0.8 !McGinnis
      !LAMBDA=1 !Socolofsky et al. 2008
      !FRNI=1.6  !McGinnis
      PI=ACOS(-1.0)
      PSTD=101325.
      RGAS=8.314
      TSTD=293.15
      DENSE20=998.2
      FRCNATM=0.79
      QGFRAC=1.0
!      C1=-1
      C1=0 !Socolofsky et al. 2008
!     
C     --------------------------------------------------------------------------
C     PARAMETERS AND INITIALIZE THE VARIABLES 
C     --------------------------------------------------------------------------
C      
      INPLUME(:,:)=0.0
!     Geometric caracteristic 
      DEPTH=WSEL-DIFFEL  
  !    PRINT*,'DEPTHINNER_PLUME', DEPTH
      Z=DEPTH
      ELEV=DIFFEL
!        
C     Assume that gas bubbles are composed of oxygen and nitrogen only.
      FRACO=FRCONOT 
      FRACN=1.0-FRACO
!
C     Interpolate input profiles to obtain the plume initial conditions
      X=0.
      XLOC=DEPTH-X

      CALL LININT(LOCAT,TE,LAYERS,XLOC,TAMB)
      CALL LININT(LOCAT,CO2M,LAYERS,XLOC,COMG)
      CALL LININT(LOCAT,UA,LAYERS,XLOC,UAMB)
      CALL LININT(LOCAT,VA,LAYERS,XLOC,VAMB)
      COMGP=COMG
      DOARO=COMG
      CO2=COMG/32.
      COMGNOT=COMG
!     At first the temperature and the salinity un the plume is the same that in the ambient water
      SALPLU=SALAMB    
      TPLUME=TAMB
!     Guess the initial velocity 
      VGUESS=0.07
      VI=VGUESS
!     INITIAL DIFFUSER SIZE
!     At first the plume radius BNOT (initial plume radius) in the TOP HAT models 
      BNOT=LAMBNOT/LAMBDA
      BI=BNOT
!     BDIFF=LAMBNOT
C      
C     AMBIENT AND AVERAGE WATER DENSITIES
      DENSEA=(0.059385*TAMB**3-8.56272*TAMB**2+65.4891*TAMB)*0.001
     ++999.84298 !+(GAMMA)*SALAMB
      DENSEW=DENSEA
	  
	  
C
!     Solubility constant (mol/m3/Pa)
      HO2=(2.125-0.05023*TPLUME+5.7714E-4*TPLUME**2)/100000.
      HN2=(1.042-0.02457*TPLUME+3.1714E-4*TPLUME**2)/100000.
C	
C     Assume initial ambient dissolved nitrogen conc. equals saturated conc. at surface.
      CN2=(PATM*FRCNATM)*HN2      
      CNMG=CN2*28.0
      CNMGP=CNMG
      DNAMB=CNMG
!
!     Outer plume
      OUTXLOC=OUTPLUME(:,1)
      OUTQW=OUTPLUME(:,3)
      OUTMOMENT=OUTPLUME(:,4)
      OUTTEMP=OUTPLUME(:,5)
      OUTO2=OUTPLUME(:,6)
      OUTSAL=OUTPLUME(:,7)
!
C     --------------------------------------------------------------------------
C     BUBBLE PROPERTIES
C     --------------------------------------------------------------------------
!
      QSCMS=QGFRAC*QSCFM/3.281**3/60.0
      QGAS=PSTD*QSCMS*(TAMB+273.15)/((PATM+DENSEA*G*Z)*TSTD)

      PRINT*,'QGAS',QGAS,PSTD,QSCMS,TAMB,PATM,DENSEA,G,Z,TSTD


C     Initial bubble size.     
      RB=DIAMM/2000.
      !PRINT*,'RBinicial',RB
!     Bubble rise velocity
      IF(RB.LE.(7.5E-4))THEN
            VB=1189.0*RB**1.1945
      ELSEIF(RB.GT.(7.5E-4).AND.RB.LT.(4.8E-3))THEN
            VB=0.22
      ELSE
            VB=2.995*RB**0.489
      ENDIF

      PRINT*,'RBinicial',RB,VB

!     Mass transfer coefficient
      KOLO=0.6*RB
      IF(KOLO.GT.(4.0E-4))THEN
            KOLO=4.0E-4
      ENDIF
      KOLN=KOLO
      PRINT*,'Test 1'
!
C     CALCULATION OF INITIAL WATER VELOCITY USING FROUDE NUMBER
      VBUB=4./3.*PI*RB**3
      N=QGAS/VBUB
      VDIFF=1

      PRINT*,'VIinicial',FRNI,LAMBDA,BI,G,DENSEA,QGAS,VB,PI,DENSEW


      DO WHILE (VDIFF.GT.1.0E-6)
         VG=QGAS/((VGUESS+VB)*(PI*(LAMBDA*BI)**2))
         DENSEP=(1.0-VG)*DENSEW
         VI=FRNI*(2.0*LAMBDA*BI*G*(DENSEA-DENSEP)/DENSEP)**0.5
         VDIFF=ABS(VI-VGUESS)
         VGUESS=VI
      END DO
      !PRINT*,'VIinicial',VI
      !PRINT*,'VIinicialSCOTT',VI
      !PRINT*,'VGinicial',VG

C     ------------------------------------------------------------------
C     VARIABLE TRANSFORMATION      
C     ------------------------------------------------------------------
!      PRINT*,'flag1'
!     TOTAL ENTRAIMENT: Shear and Vortex
!       E=2.*PI*BI*ALPHA*VI
!      PRINT*,"Shear" 
!      SHEARE=2.*PI*BI*ALPHA*VI
!     Vortex Entraiment: Hypothesis Projected Area Entraiment
!      VORTEXE=(2*BI*(SQRT(ABS(UAMB)*ABS(UAMB)+ABS(VAMB)*ABS(VAMB)))) 
!     a) Additive hypothesis
!         PRINT*,"Additive"  
!       E=SHEARE+VORTEXE
!     b) Maximun hypothesis
!         PRINT*,"Maximun"  
!      IF (SHEARE.GT.VORTEXE) THEN
! 	E=SHEARE
!      ELSE
!	E=VORTEXE
!      ENDIF

      VO=0
      EI=2.*PI*BI*ALPHAI*(VI+C1*VO)
      EO=0
      QW=VI*PI*BI**2
      MOMENT=(PI*BI**2)*VI**2
      FTEMP=QW*TPLUME 
      !FSAL=QW*(SALPLU*GAMMA/DENSE20)*DENSEW JCT2020_Sal
      FSAL=QW*SALPLU !JCT2020_Sal
C     Previous equation corrected to account for salinity units conversion.        
      FDO=QW*CO2
      FDN=QW*CN2
      FGO=PSTD*QSCMS/(RGAS*TSTD)*FRACO
      FGONOT=FGO
      FGN=PSTD*QSCMS/(RGAS*TSTD)*FRACN
C     Revised gaseous flux equations.
      YO2=FGO/((PI*(LAMBDA*BI)**2)*(VI+VB))
      YN2=FGN/((PI*(LAMBDA*BI)**2)*(VI+VB))           
      PZ=PATM+(DENSEA*G*Z)
      PO=PZ*FRACO
      PN=PZ*FRACN
      BUOY=(G*(DENSEA-DENSEP)/DENSEP*QW)/LNOT
      TDS=SALPLU*0.64
!     Initialize lateral withdrawal flowrate for first/lowest cell in column/segment
      JJ=0
      QWDI(LAYDIFF)= QW 
      PWDI(LAYDIFF)= (2*PI*BI)
      AWDI(LAYDIFF)= (PI*BI**2)    
      BWDI(LAYDIFF)= BI 
!

      PRINT*,'Initial',BI, VI, QW
!     -------------------------------------------------------------------------- 
C	SOLUTION PROCEEDURE
!     --------------------------------------------------------------------------
!
      DZ=0.001
      H=0.001
      NELS = 0
      HWITH=0.0
      NLI=0
      NLIO=NLO
!
      M=0
      !PRINT*,'VIinicial',VI
      !PRINT*,'MOMENTinicial',MOMENT
      DO WHILE (VI.GT.1.E-6.AND.Z.GT.0.0)
         M=M+1
      !DO WHILE (M.LT.1)
       !PRINT*,'VI',VI
       !PRINT*,'MOMENT',MOMENT
       !PRINT*,'Z',Z
!      PRINT*,'EO',EO
         Z=Z-DZ
         X=X+DZ
         ELEV=ELEV+DZ
         NELS = NELS + 1
         NLI = NLI+1
!       PRINT*,'NLI',NLI
C	
C        Interpolate input profiles to obtain line plume boundary conditions
         XLOC=DEPTH-X

         !PRINT*,'if',XLOC,DEPTHINTR,ELEVT

         CALL LININT(LOCAT,UA,LAYERS,XLOC,UAMB)
         CALL LININT(LOCAT,VA,LAYERS,XLOC,VAMB)
         IF (NITERPLUME.EQ.1)THEN
            CALL LININT(LOCAT,CO2M,LAYERS,XLOC,COMG)
            DOARO=COMG
            CALL LININT(LOCAT,TE,LAYERS,XLOC,TARO)
            TAMB=TARO
            VO=0
            SALARO=SALAMB
         ELSEIF (NITERPLUME.GT.1)THEN
            IF (-XLOC.LE.DEPTHINTR.OR.-XLOC.GT.ELEVT)THEN
               CALL LININT(LOCAT,CO2M,LAYERS,XLOC,COMG)
               DOARO=COMG
               CALL LININT(LOCAT,TE,LAYERS,XLOC,TARO)
               TAMB=TARO
               !VO=0
               SALARO=SALAMB
            ELSEIF (-XLOC.GT.DEPTHINTR.AND.-XLOC.LE.ELEVT)THEN
            !   !PRINT*,'NLI',NLI, NLIO
               COMG=OUTO2(NLIO)
               DOARO=COMG
               CALL LININT(LOCAT,TE,LAYERS,XLOC,TAMB)

               TARO=OUTTEMP(NLIO)
               !TARO=TAMB
               SALARO=OUTSAL(NLIO)
               !SALARO=SALAMB

               OUTERQW=OUTQW(NLIO)
               OUTERMOMENT=OUTMOMENT(NLIO)
               VO=OUTERMOMENT/OUTERQW
               !VO=0
               !VO=-0.03

               NLIO=NLIO+1
               !PRINT*,'VO', VO 
            ENDIF
          ELSE
            PRINT*, "---------------------------------"
            PRINT*, "-------ERROR NITERPLUME----------"
            PRINT*, "---------------------------------"
         ENDIF
C                  
C        Use subroutines for Runge Kutta method solution
         NEQN=8
!      PRINT*,'QWbeforeEQ',QW
!      PRINT*,'MOMENTbeforeEQ',MOMENT

         Y(1)=QW
         Y(2)=MOMENT
         Y(3)=FTEMP
         Y(4)=FSAL
         Y(5)=FDO
         Y(6)=FDN
         Y(7)=FGO
         Y(8)=FGN



         CALL DERIVS_4(EI,EO,DENSEA,DENSEW,DENSEP,G,BI,LAMBDA,TARO,VG,
     +            SALARO,GAMMA,DENSE20,DOARO,PI,RB,N,VI,VO,VB,KOLO,HO2,
     +            PO,GAMMA1,TPLUME,SALPLU,COMGP,DNAMB,KOLN,HN2,PN,
     +            CNMGP,Z,Y,DYDX,XLOC,TAMB)

         CALL RK4_4(EI,EO,DENSEA,DENSEW,DENSEP,G,BI,LAMBDA,TARO,VG,
     +         SALARO,GAMMA,DENSE20,DOARO,PI,RB,N,VI,VO,VB,KOLO,HO2,PO,
     +         GAMMA1,TPLUME,SALPLU,COMGP,DNAMB,KOLN,HN2,PN,CNMGP,Y,
     +         DYDX,NEQN,Z,H,YOUT,XLOC,TAMB)

!         CALL DERIVS_3(EI,DENSEA,DENSEW,DENSEP,G,BI,LAMBDA,TARO,SALARO,
!     +            GAMMA,DENSE20,DOARO,PI,RB,N,VI,VB,KOLO,HO2,PO,
!     +            COMGP,DNAMB,KOLN,HN2,PN,CNMGP,Z,Y,DYDX)

!         CALL RK4_3(EI,DENSEA,DENSEW,DENSEP,G,BI,LAMBDA,TARO,SALARO,
!     +         GAMMA,DENSE20,DOARO,PI,RB,N,VI,VB,KOLO,HO2,PO,
!     +         COMGP,DNAMB,KOLN,HN2,PN,CNMGP,Y,DYDX,NEQN,Z,H,YOUT)


    !  PRINT *, 'EI',EI
    !  PRINT *, 'DENSEA',DENSEA
    !  PRINT *, 'DENSEW',DENSEW
    !  PRINT *, 'DENSEP',DENSEP
    !  PRINT *, 'G',G
    !  PRINT *, 'BI',BI
    !  PRINT *, 'LAMBDA',LAMBDA
    !   PRINT *, 'TARO',TARO
    !  PRINT *, 'SALARO',SALARO
    !  PRINT *, 'GAMMA',GAMMA
    !  PRINT *, 'DENSE',DENSE20
    !  PRINT *, 'DOARO',DOARO
    !  PRINT *, 'PI',PI
    !  PRINT *, 'RB',RB
    !  PRINT *, 'N',N
    !  PRINT *, 'VI',VI
    !  PRINT *, 'VB',VB
    !  PRINT *, 'KOLO',KOLO
    !  PRINT *, 'HO2',HO2
    !  PRINT *, 'PO',PO
    !  PRINT *, 'COMGP',COMGP
    !  PRINT *, 'DNAMB',DNAMB
    !  PRINT *, 'KOLN',KOLN
    !  PRINT *, 'HN2',HN2
    !  PRINT *, 'PN',PN
    !  PRINT *, 'CNMGP',CNMGP
    !  PRINT *, 'Z',Z
    !  PRINT *, 'Y',Y
    !  PRINT *, 'DYDX',DYDX




!    PRINT*,'DENSEA', DENSEA
!    PRINT*,'DENSEW', DENSEW
!    PRINT*,'DENSEP', DENSEP
!    PRINT*,'G', G
!    PRINT*,'PI', PI
!    PRINT*,'BI', BI
!    PRINT*,'LAMBDA', LAMBDA

         QW=YOUT(1)
         MOMENT=YOUT(2)
         FTEMP=YOUT(3)
         FSAL=YOUT(4)
         FDO=YOUT(5)
         FDN=YOUT(6)
         FGO=YOUT(7)
         FGN=YOUT(8)


!          OPEN (UNIT=58, FILE="inner_plume_eq.txt", POSITION="APPEND")
!             WRITE (UNIT=58, FMT='(9F12.6)') QW,
!      +        MOMENT,FTEMP,
!      +        FSAL,FDO,
!      +        FDN,FGO,FGN,X
!          CLOSE (UNIT=58)



         IF(MOMENT.LT.0.0)THEN
            TPLUME=FTEMP/QW
            !SALPLU=FSAL/(QW*DENSEW)/(GAMMA/DENSE20) JCT2020_Sal
            SALPLU=FSAL/QW  !JCT2020_Sal
C        Previous equation corrected to consistently express salinity in uS/cm	   
             PRINT*,'fuera momento'
	    CO2=FDO/QW
	    CN2=FDN/QW
!           Save inner plume information
            INPLUME(NLI,1)=XLOC 
            INPLUME(NLI,2)=BI          
            INPLUME(NLI,3)=QW
            INPLUME(NLI,4)=MOMENT
            INPLUME(NLI,5)=TPLUME
            INPLUME(NLI,6)=COMGP
            INPLUME(NLI,7)=EI 
            INPLUME(NLI,8)=EO
            PRINT*,'MOMENT',MOMENT
	    GOTO 20
         ENDIF


         VI=MOMENT/QW
       !PRINT*,'VI',VI,MOMENT,QW
!      PRINT*,'QW',QW
!      PRINT*,'MOMENT',MOMENT
!      PRINT*,'VI',VI
!        PRINT*,"VI",VI
         AREA=QW/VI
         BI=SQRT(AREA/PI)
!        PRINT*,"QW",QW
!        PRINT*,"BI",BI
!        TOTAL ENTRAIMENT: Shear and Vortex
!          E=2.*PI*BI*ALPHAI*VI
!         PRINT*,"Shear" 
!         SHEARE=2.*PI*BI*ALPHAI*VI
!        Vortex Entraiment: Hypothesis Projected Area Entraiment
!         VORTEXE=(2*BI*(SQRT(ABS(UAMB)*ABS(UAMB)+ABS(VAMB)*ABS(VAMB)))) 
!	 a) Additive hypothesis 
!          PRINT*,"Additive" 
!      	     E=SHEARE+VORTEXE
!	 b) Maximun hypothesis
!          PRINT*,"Maximun" 
!	    IF (SHEARE.GT.VORTEXE) THEN
!	 	E=SHEARE
!            ELSE
!		E=VORTEXE
!            ENDIF
!         

         !AU=VO
         !VO=0.0

         EI=2.*PI*BI*ALPHAI*(VI+C1*VO)
         EO=-2.*PI*BI*ALPHAO*VO  ! JCT

         !VO=AU

         !PRINT*,XLOC,EI,EO
         !PRINT*,XLOC

         !IF (M.EQ.500) THEN
!        !    PRINT*, XLOC,BI,VI,UAMB,VAMB,SHEARE,VORTEXE
         !    WRITE (UNIT=50, FMT='(7F12.6)') XLOC,BI,VI,UAMB,VAMB,SHEARE,VORTEXE
         !    M=O
         !ENDIF

!            PRINT*, XLOC,BI,VI,UAMB,VAMB,SHEARE,VORTEXE
!        WRITE (UNIT=50, FMT='(8F8.4)') NLI,INPLUME(NLI,1),INPLUME(NLI,2),INPLUME(NLI,3),INPLUME(NLI,4),INPLUME(NLI,5),INPLUME(NLI,6),INPLUME(NLI,7)
        !WRITE (UNIT=50, FMT='(1F8.4)') NLI
         

!
!        Temperatura and salinity in the plume
         TPLUME=FTEMP/QW
         !SALPLU=FSAL/(QW*DENSEW)/(GAMMA/DENSE20) JCT2020_Sal
         SALPLU=FSAL/QW  !JCT2020_Sal
!        ACC 2026: Override SALPLU=SALAMB eliminado. Era un parche para enmascarar el error en
!        DERIVS_4 donde el termino de detrainment de salinidad usaba SALARO en lugar de SALPLU.
!        Corregido DERIVS_4: DYDX(4)=EI*SALARO-EO*SALPLU. Linea erronea comentada: ACC 2026
!         SALPLU=SALAMB
!         !!!!! ******* !!!!! ******* !!!!! ******* !!!!!*******  !!!!!

         !PRINT*, XLOC, TPLUME
C        Previous equation corrected to consistently express salinity in uS/cm
!        Dissolved oxygen and nitrogen concentration
         CO2=FDO/QW
         CN2=FDN/QW 
         COMGP=CO2*32.
         CNMGP=CN2*28.
!        Save inner plume information
         INPLUME(NLI,1)=XLOC 
         INPLUME(NLI,2)=BI  
         INPLUME(NLI,3)=QW
         INPLUME(NLI,4)=MOMENT
         INPLUME(NLI,5)=TPLUME
         INPLUME(NLI,6)=COMGP
         INPLUME(NLI,7)=EI 
         INPLUME(NLI,8)=EO
         !PRINT*,'INPLUME', INPLUME(NLI,:)


 !        OPEN (UNIT=57, FILE="inner_plume.txt", POSITION="APPEND")
 !           WRITE (UNIT=57, FMT='(8F12.6)') INPLUME(NLI,1),
 !    +        INPLUME(NLI,2),INPLUME(NLI,3),
 !    +        INPLUME(NLI,4),INPLUME(NLI,5),
 !    +        INPLUME(NLI,6),INPLUME(NLI,7),INPLUME(NLI,8)
 !        CLOSE (UNIT=57)


         !PRINT*,'Temperatura', XLOC, TAMB,TARO,TPLUME

!        Add incremental entrainment to total cell entrainment/withdrawal   
         QWDI(LAYDIFF-JJ)=QWDI(LAYDIFF-JJ)+(EI-EO)*DZ
         PWDI(LAYDIFF-JJ)=PWDI(LAYDIFF-JJ)+2.*PI*BI
         AWDI(LAYDIFF-JJ)=AWDI(LAYDIFF-JJ)+PI*BI**2
         BWDI(LAYDIFF-JJ)=BWDI(LAYDIFF-JJ)+SQRT((PI*BI**2)/PI)
         HWITH=HWITH+DZ   
         IF(HWITH.GT.HCELL)THEN
            PWDI(LAYDIFF-JJ) = PWDI(LAYDIFF-JJ)/NELS
            AWDI(LAYDIFF-JJ) = AWDI(LAYDIFF-JJ)/NELS
            BWDI(LAYDIFF-JJ) = BWDI(LAYDIFF-JJ)/NELS
            JJ=JJ+1
            HWITH=0.0
            QWDI(LAYDIFF-JJ)=0.0
            PWDI(LAYDIFF-JJ)=0.0
            BWDI(LAYDIFF-JJ)=0.0
            AWDI(LAYDIFF-JJ)=0.0
            NELS = 0
         ENDIF
!
C        Revised gaseous flux equations.
         YO2=FGO/((PI*(LAMBDA*BI)**2)*(VI+VB))
         YN2=FGN/((PI*(LAMBDA*BI)**2)*(VI+VB))
C  
!        
         PZ=PATM+(DENSEA*G*Z)
         QGAS=(FGO+FGN)*RGAS*(TPLUME+273.15)/PZ
         VBUB=QGAS/N
C
!        GAS VOLUME PER TOTAL VOLUME OF THE BUBBLE-WATER MIXTURE IN THE INNER CORE OF THE PLUME
         VG=VBUB*N/((VI+VB)*(PI*(LAMBDA*BI)**2))
!         PRINT*,"VG",XLOC,VG,VBUB,N,VI,VB,PI,LAMBDA,BI
C        Previous equation revised to account for correct plume cross-sectional area occupied by bubbles. 
!        Bubbles radius     
         RB=(3.*QGAS/(4.*PI*N))**(1./3.)
         IF(RB.LT.0.0)THEN
            RB=1.0E-8
         ENDIF
         !PRINT*,'RB',RB
         FRACO=FGO/(FGO+FGN)
         FRACN=1.0-FRACO
C	
!        Partial Pressure
         PO=PZ*FRACO
         PN=PZ*FRACN
!        Density (ambient water, water in plume and of the plume)
      DENSEA=(0.059385*TAMB**3-8.56272*TAMB**2+65.4891*TAMB)*0.001
     ++999.84298 !+(GAMMA)*SALARO
      DENSEW=(0.059385*TPLUME**3-8.56272*TPLUME**2+65.4891*TPLUME)*0.001
     ++999.84298 !+(GAMMA)*SALPLU
         !PRINT*,'JCT',SALARO,SALAMB,SALPLU

C        Previous equation re-revised to account for correct salinity units (uS/cm) in density calculations.      
         DENSEP=(1.0-VG)*DENSEW
C
C        BUBBLE PROPERTIES
!
!        Bubble rise velocity
         IF(RB.LE.(7.5E-4))THEN
            VB=1189.0*RB**1.1945
         ELSEIF(RB.GT.(7.5E-4).AND.RB.LT.(4.8E-3))THEN
            VB=0.22
         ELSE
            VB=2.995*RB**0.489
         ENDIF
C
!        Mass transfer coeffitient
         KOLO=0.6*RB
         IF(KOLO.GT.(4.0E-4))THEN
            KOLO=4.0E-4
         ENDIF
         KOLN=KOLO
C
!        Solubility constant (mol/m3/Pa)
         HO2=(2.125-0.05023*TPLUME+5.7714E-4*TPLUME**2)/100000.
         HN2=(1.042-0.02457*TPLUME+3.1714E-4*TPLUME**2)/100000.
C
!        Froude number
         FR=VI/(2.*LAMBDA*BI*G*(DENSEA-DENSEP)/DENSEP)**0.5
!
         DCO2=HO2*PO-CO2
C      
      END DO
C
!     --------------------------------------------------------------------
C     CALCULATION OF AVERAGE NET OXYGEN MASS TRANSFER FOR DAY
!     --------------------------------------------------------------------
   20 GROSSMT=(FGONOT-FGO)*32./1000.*86400.
      OTEFF=(FGONOT-FGO)/FGONOT*100.
      DELTAC=COMGP-COMGNOT
C
      ELEVT=ELEV
      QWT=QW
      TPLUMET=TPLUME
      COMGPT=COMGP
      LAYTOP=LAYDIFF-JJ
      PRINT*,"LAYTOP_inner",LAYTOP
      PRINT*,"LAYDIFF_inner",LAYDIFF
      PRINT*,"JJ_inner",JJ
      BTOP=BI
      SALPLUMET=SALPLU

!      PRINT*,"Minner",M
!      PRINT*,"NLI",NLI
!      PRINT*,"XLOC",XLOC
!      PRINT*,"ELEVT",ELEVT
!      PRINT*,"QWT",QWT
!      PRINT*,"TPLUMET",TPLUMET
!      PRINT*,"COMGPT",COMGPT
!      PRINT*,"LAYTOP",LAYTOP
!      PRINT*,"BTOP",BTOP
!      PRINT*,"INPLUMEDEPTH",INPLUME(:,1)
      PRINT*,'fuera final'
      RETURN
      END
C
C------------------------------------------------------------------------------
C                          
             
C     *******************************************************************************
	SUBROUTINE OUTER_PLUME(YEAR,JULDAY,WSEL,DIFFEL,LAYERS,BITOP,
     +SALAMB,PATM,QINTOP,FRCONOT,LAYDIFF,HCELL,LOCAT,TE,CO2M,UA,VA,
     +ELEVT,QWDET,QWDO,BWDO,LAYTOP,LAYINTR,TPLUMED,SALPLUMED,COMGPD,
     +DEPTHINTR,NLI,NLO,NITERPLUME,INPLUME,OUTPLUME,ALPHAI,ALPHAO, 
     +ALPHAA,LAMBDA,FRNI,FRNO,GAMMA1)
C     *******************************************************************************

      IMPLICIT NONE

C	THIS SUBROUTINE IS WRITTEN TO PREDICT THE PERFORMANCE OF A CIRCULAR BUBBLE PLUME.  
C     THE MODEL IS BASED ON THE WUEST ET AL. (1992) CIRCULAR BUBBLE PLUME MODEL.
C     By:
C	
C	VERSION 3 (to couple with Francisco Rueda's reservoir model and for Lake Hallwil real diffuser)
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
C	  ALPHAI=ENTRAINMENT COEFFICIENT INNER PLUME  (-)
C	  ALPHAO=ENTRAINMENT COEFFICIENT OUTER PLUME  (-)
C	  ALPHAA=ENTRAINMENT COEFFICIENT FROM AMBIENT (-)
!         BI=1/2 INNER PLUME WIDTH (m)
!         BO=1/2 OUTER PLUME WIDTH (m)
!     BIFF=DIFFUSER RADIUS (m) 
C     BAVG=AVERAGE 1/2 DIFFUSER WIDTH (m)
!         BITOP=
!         BDET          	
!         C1=COEFFICIENT OF DIFERENT COUNTERFLOW ENTRAINMENTS
C	CO2=DISSOLVED OXYGEN (DO) CONCENTRATION (mol/m3)
C     CO2M=DO CONCENTRATION PROFILE FOR INPUT BOUNDARY CONDITION (g/m3)
C	CN2=DISSOLVED NITROGEN CONCENTRATION (mol/m3)
C	COMG=DISSOLVED OXYGEN CONCENTRATION (g/m3)
C     COMGPT=DO CONCENTRATION OF PLUME DETRAINMENT AT TOP OF PLUME (g/m3) 
!         COMGPD
!         COMGPTI
C	CNMG=DISSOLVED NITROGEN CONCENTRATION (g/m3)
C     DENSE20=DENSITY OF WATER AT 20 C
C	DENSEA=AMBIENT WATER DENSITY (kg/m3)
C	DENSEP=DENSITY OF THE PLUME (kg/m3)
C	DENSEW=WATER DENSITY IN PLUME (kg/m3)
!         DEPTHINTR=INTRUSION DEPTH OF THE DOUBLE PLUME
C	DIAMM=BUBBLE DIAMETER (mm)
C     DIFFEL=DIFFUSER ELEVATION (m)
C     DMPR=DEPTH OF MAXIMUM PLUME RISE (m)
C     DNAMB=AMBIENT DISSOLVED NITROGEN CONCENTRATION (g/m3)
C     DOAMB=AMBIENT DISSOLVED OXYGEN CONCENTRATION (g/m3)
C	E=ENTRAINMENT FACTOR (m3/s)
!	  EI=ENTRAINMENT  (m3/s)
!	  EO=ENTRAINMENT  (m3/s)
!         EA=ENTRAINMENT  (m3/s)
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
!	FR=FROUDE NUMBER (-)
C     FRCNOT=FRACTION OF NITROGEN IN ATMOSPHERE (-)
C     FRNO=INITIAL FROUDE NUMBER (-)
C     GAMMA=SALINITY CONVERSION FACTOR [kg/m3/(uS/cm)]
C     GROSSMT=GROSS MASS TRANSFER OF OXYGEN FROM PLUME (kg/d)
!	HO2=SOLUBILITY CONSTANT FOR OXYGEN (mol/m3/Pa)
!	HN2=SOLUBILITY CONSTANT FOR NITROGEN (mol/m3/Pa)
C     HCELL=HEIGHT OF CELL IN GRID (m)
C	HO=SOLUBILITY CONSTANT FOR OXYGEN (mol/m3/Pa)
C	HN=SOLUBILITY CONSTANT FOR NITROGEN (mol/m3/Pa)
C     HWITH=HEIGHT OF WITHDRAWAL/ENTRAINMENT ZONE (m)
!         INPLUME(:,6): INNER PLUME CHARACTERISTIC [HEIGTH; FLOWRATE OF WATER; MOMENTUM; TEMPERATURE; SALINITY; O2] JCT
C     JULDAY=JULIAN DAY IN GIVEN YEAR
C	KOLO=MASS TRANSFER COEFFICIENT FOR OXYGEN (m/s)
C	KOLN=MASS TRANSFER COEFFICIENT FOR NITROGEN (m/s)
C     LAKE=LAKE AND DIFFUSER TYPE (SHR AND LINEAR=1 OR AMISK AND RECTANGULAR=2) FOR SELECTION OF LAMBDA
C     LAYERS=NUMBER OF LAYERS/DATA POINTS IN BOUNDARY CONDITION PROFILES (-)
C     LAYDIFF=GRID LAYER CORRESPONDING TO DIFFUSER DEPTH (-)
!         LAYINTR=GRID LAYER CORRESPONDING TO INTRUSION DEPTH OF THE DOUBLE PLUME (-)
C     LAYTOP=GRID LAYER CORRESPONDING TO TOP OF PLUME (-)
C     LDIFF=LENGTH OF DIFFUSER (m)
C     LNOT=INITIAL DIFFUSER LENGTH (m)
C	LAMBDA=FRACTION OF PLUME OCCUPIED BY BUBBLES (-)
C     LAMBNOT=LAMBDA x INITIAL PLUME RADIUS;EQUAL TO DIFFUSER RADIUS (m) 
C     LOCAT=DEPTHS FOR INPUT BOUNDARY CONDITION PROFILES (m) 
C	MOMENT=MOMENTUM (m4/s)
C	N=NUMBER OF BUBBLES PER SECOND (1/s)
C     OTEFF=OXYGEN TRANSFER EFFICEINCY (%)
!         OUTPLUME(:,6)= OUTER PLUME CHARACTERISTIC [HEIGTH; FLOWRATE OF WATER; MOMENTUM; TEMPERATURE; SALINITY; O2] JCT
!         OUTERMOMENT= MOMENTUM IN THE OUTER PLUME 
!         OUTERQW=FLOWRATE OF WATER IN THE OUTER PLUME
!         OUTXLOC
!         OUTQW
!         OUTMOMENT
!         OUTTEMP
!         OUTO2  
C     PATM=ATMOSPHERIC PRESSURE AT AVERAGE WSEL (Pa)
C	PSTD=STANDARD PRESSURE (Pa)
!     QGAS=TOTAL GAS FLOW TO DIFFUSER
!       QINTOP
C!!!!	QSCFM=STANDARD GAS FLOW RATE (scfm), TOTAL STANDAR GAS FLOW RATE TO DIFFUSER 
C	QSCMS=STANDARD GAS FLOW RATE (scms)
C	QW=FLOWRATE OF WATER (m3/s)
!         QWDO
C     QWDET=TOTAL DETRAINMENT FLOW RATE AT THE INTRUSION OF THE PLUME (m3/s)  
C	RB=BUBBLE RADIUS (m)
C	RGAS=IDEAL GAS CONSTANT (J/mol/K)
C     SAL=SALINITY (uS/cm)
C     SALAMB=SALINITY OF AMBIENT WATER (uS/cm)
C	SALPLU=SALINITY OF THE PLUME (uS/cm)
!         SALPLUMED
!         SALPLUMETI
!     SHEARE=SHEAR ENTRAINMENT
C	  TARO=AROUND WATER TEMPERATURE (C)
C     TAVG=AVERAGE AMBIENT WATER TEMPERATURE (C)
C     TDS=TOTAL DISSOLVED SOLIDS (g/m3) [0.64 conversion factor from Chapra book]
C     TE=TEMPERATURE PROFILE FOR INPUT BOUNDARY CONDITION (C) 
C	TPLUME=PLUME TEMPERATURE (C)
C     TPLUMET=DETRAINMENT PLUME TEMPERATURE AT THE TOP OF THE PLUME (C) 
!         TPLUMED
!         TSTD=STANDARD TEMPERATURE (K)
!     UA = U AMBIENT VELOCITY PROFILE FOR INPUT BOUNDARY CONDITION
C     UAMB = U Ambient Velocity (m) used to calculate VORTEX entrainment (FJR
!	  VI=WATER VELOCITY (m/s)
!         VO=WATER VELOCITY (m/s)
!     VA = V AMBIENT VELOCITY PROFILE FOR INPUT BOUNDARY CONDITION
C     VAMB = V Ambient Velocity (m) used to calculate VORTEX entrainment (FJR
C     VAVG=AVERAGE WATER VELOCITY (m/s)
C	VB=BUBBLE RISE VELOCITY (m/s)
C	VBUB=BUBBLE VOLUME (m3)
!     VG=GAS VOLUME PER TOTAL VOLUME OF THE BUBBLE-WATER MIXTURE IN THE INNER CORE OF THE PLUME
C	VGUESS=GUESSED INITIAL WATER VELOCITY (m/s)
!         VI=WATER VELOCITY IN THE INNER PLUME (m/s)
!         VO=WATER VELOCITY IN THE OUTER PLUME (m/s)
!     VORTEXE=VORTEX ENTREAINMENT
C     WSEL=WATER SURFACE ELEVATION (m)
C	YO2=GASEOUS OXYGEN CONCENTRATION (mol/m3)
C	YN2=GASEOUS NITROGEN CONCENTRATION (mol/m3)
C	Z=DEPTH TO DIFFUSER (m)
C     PWD  = Perimeter (m) (FJR)
C     VUP  = Upward velocity (m/2)
C     AWD  = Area (m2) of plume
C     BWD  = Width 
C
      REAL*8 AREA,BI,BO,CO2,COMG,CN2,CNMG,DS,DZ,DENSEA,DENSEP,DENSEW,OTEFF,
     +DIAMM,FDO,FDN,FRACO,FRACN,FSAL,FTEMP,FGO,FGN,G,GAMMA,HO2,HN2,
     +KOLN,KOLO,LAMBDA,MOMENT,N,PI,PO,PN,PSTD,PZ,QSCMS,QW,QGAS,
     +RB,RGAS,SALAMB,SALPLU,TPLUME,TSTD,VB,VBUB,VG,VGUESS,
     +YO2,YN2,Z,AA,BB,CC,BNOT,LNOT,ELEV,DT,TE(1000),XLOC,DEPTH,
     +CO2M(1000),LOCAT(1000),DOAMB,COMGP,CNMGP,WSEL,PATM,SAL(1000),
     +DIFFEL,GROSSMT,FGONOT,DNAMB,DMPR,SUMSAL,FRCNATM,
     +LAMBNOT,FRCONOT,H,DYDX(5),Y(5),YOUT(5),BTOP,
     +DENSE20,FRNI,FRNO,VDIFF,FR,DCO2,QGFRAC,LDIFF,
     +JDAY,EL(70),DELTAC,COMGNOT,ELEVT,QWDET,QWDO(500),
     +PWDO(500),BWDO(500),AWDO(500),X,SHEARE,VORTEXE,
     +HWITH,HCELL,JULDAY,BDET,UA(1000),VA(1000),UAMB,VAMB,
     +INPLUME(35000,9),OUTPLUME(35000,9),DEPTHINTR,ALPHAI,ALPHAO,ALPHAA,
     +C1,EI,EO,EA,INNERQW,INNERMOMENT,INXLOC(35000),INQW(35000),
     +INMOMENT(35000),INTEMP(35000),INO2(35000),INSAL(35000),
     +INBI(35000),BITOP,VI,VO,TPLUMED,QINTOP,TARO,DEPTHDMPR,
     +SALPLUMED,COMGPD,DENSEINNER,GAMMA1,TIN,COMGIN,SALIN,EP  
      INTEGER II,IJ,IK,JJ,LL,NEQN,NN,MI,JK,JL,LAYTOP,
     +LAYERS,KM,KN,KO,KP,ROWS,KQ,KR,KU,KV,KW,KX,KY,KZ,KS,LAYDIFF,YEAR 
      INTEGER NELS,ierror,M,CONT,LAYINTR,NLI,NLO,NITERPLUME
C 
C     --------------------------------------------------------------------------
C     CONSTANTS
C     --------------------------------------------------------------------------
!

      PRINT*,'OUTER_PLUME'

      PRINT*,'NITERPLUMEinOUTERPLUME', NITERPLUME

      PRINT*, 'OUTER_PLUME',ALPHAI,ALPHAO,ALPHAA,LAMBDA,FRNI,FRNO,GAMMA1


    !  PRINT*,'Yearplume', YEAR
    !  PRINT*,'Julday', JULDAY
    !  PRINT*,'WSEL', WSEL
    !  PRINT*,'DIFFEL', DIFFEL
    !  PRINT*,'LAYERS', LAYERS
    !  PRINT*,'BITOP', BITOP
      ! PRINT*,'SALAMB', SALAMB
    !  PRINT*,'PATM', PATM
    !  PRINT*,'QINTOP', QINTOP
    !  PRINT*,'FRCONOT',FRCONOT
    !  PRINT*,'LAYDIFF',LAYDIFF
      ! PRINT*,'HCELL',HCELL
      !PRINT*,'LOCAT',LOCAT
      !PRINT*,'TE',TE
    !  !PRINT*,'CO2M',CO2M
    !  !PRINT*,'UA',UA
    !  !PRINT*,'VA',VA
    !  PRINT*,'ELEVT',ELEVT
    !  PRINT*,'QWDET',QWDET
    !  !PRINT*,'QWDO',QWDO
    !  !PRINT*,'BWDO',BWDO
    !  PRINT*,'LAYTOP',LAYTOP
    !  PRINT*,'LAYINTR',LAYINTR
    !  PRINT*,'TPLUMED',TPLUMED
    !  PRINT*,'SALPLUMED',SALPLUMED
    !  PRINT*,'COMGPD',COMGPD
    !  PRINT*,'DEPTHINTR',DEPTHINTR
    !  PRINT*,'NLI',NLI
    !  PRINT*,'NLO',NLO
    !  PRINT*,'NITERPLUME',NITERPLUME
    !  !PRINT*,'INPLUME',INPLUME
    !  !PRINT*,'OUTPLUME',OUTPLUME


      !ALPHAI=0.055 ! Crounse et al 2007
      !ALPHAO=0.11  ! Crounse et al 2007
      !ALPHAA=0.11  ! Crounse et al 2007
      G=9.80665
      GAMMA=6.9E-4
      !GAMMA1=1.1 ! JCT Socolofsky
!      LAMBDA=0.8 !McGinnis
      !LAMBDA=1 !Socolofsky et al. 2008
      !FRNO=0.1  !Socolofsky et al. 2008
      !FRNO=-FRNO  !Socolofsky et al. 2008 Para ser coherente en los signos
      PI=ACOS(-1.0)
      PSTD=101325.
      RGAS=8.314
      TSTD=293.15
      DENSE20=998.2
      FRCNATM=0.79
      QGFRAC=1.0
!      C1=-1
      C1=0
!     
C     --------------------------------------------------------------------------
C     PARAMETERS AND INITIALIZE THE VARIABLES 
C     --------------------------------------------------------------------------
C     
      OUTPLUME(:,:)=0.0
!     Information from the INNER plume

      INXLOC=INPLUME(:,1)
      INBI=INPLUME(:,2)
      INQW=INPLUME(:,3)
      INMOMENT=INPLUME(:,4)
      INTEMP=INPLUME(:,5)
      INO2=INPLUME(:,6)
      INSAL=INPLUME(:,7)
	  
	  ! PRINT*,'NLI', NLI
	  ! PRINT*,'INXLOC', INXLOC(1), INXLOC(NLI), INXLOC(26929)
	  ! PRINT*,'INBI', INBI(1), INBI(NLI), INBI(26929)
	  ! PRINT*,'INQW', INQW(1), INQW(NLI), INQW(26929)
	  ! PRINT*,'INMOMENT', INMOMENT(1), INMOMENT(NLI), INMOMENT(26929)
	  ! PRINT*,'INTEMP', INTEMP(1), INTEMP(NLI), INTEMP(26929)
	  ! PRINT*,'INO2', INO2(1), INO2(NLI), INO2(26929)
	  ! PRINT*,'INSAL', INSAL(1), INSAL(NLI), INSAL(26929)


!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!     Dejamos como zona de transicion los 0.5 metros superiores de la pluma. JCT
!     BITOP y QINTOP dejan de ser input, los sacamos de INPLUME --> ELIMINAR  JCT
      ! NLI = NLI-10
      NLI = NLI -10
      BITOP  = INPLUME(NLI,2)
      QINTOP = -INPLUME(NLI,3) !JCT2022
	  
	  
	  ! PRINT*,'INPLUME(NLI,3)',INPLUME(NLI,3)

	  
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      

!     Geometric caracteristic 
      DEPTH=WSEL-DIFFEL
!      DEPTHDMPR=WSEL-ELEVT+0.5 ! Restamos 0.5 metros para eliminar la zona de transicion
      DEPTHDMPR=WSEL-ELEVT ! Restamos 0.5 metros para eliminar la zona de transicion
      Z=DEPTHDMPR   
      ELEV=ELEVT
!        
!!C     Assume that gas bubbles are composed of oxygen and nitrogen only.
!!      FRACO=FRCONOT 
!!      FRACN=1.0-FRACO
!
C     Interpolate input profiles to obtain the plume initial conditions
      X=0.0 
      XLOC=Z-X 
!      PRINT*,'DEPTH',DEPTH
!      PRINT*,'XLOCouter',XLOC

!      PRINT*,'INXLOC',INXLOC
!      PRINT*,'XLOC',XLOC

!      PRINT*,'Que pasa aki'
!      PRINT*,'INXLOC',INXLOC
!      PRINT*,'XLOC',XLOC
!      PRINT*,'INXLOC(1)',INXLOC(1)
    !  PRINT*,'INXLOC(INNERLENGTH)',INXLOC(INNERLENGTH)
!      PRINT*,'INXLOC(23875)',INXLOC(23875)
!      PRINT*,'INXLOC(23874)',INXLOC(23874)
!      PRINT*,'INXLOC(23877)',INXLOC(23877)
    !  PRINT*,'INNERLENGTH',INNERLENGTH

 !     PRINT*,'FLAGERROR1'
 !     PRINT*,'NLI1',NLI
 !     PRINT*,"XLOC",XLOC
 !     PRINT*,"LAYERS",LAYERS
 !     PRINT*,"LOCAT",LOCAT(1:LAYERS)
 !     PRINT*,"TE",TE(1:LAYERS)


      CALL LININT(LOCAT,TE,LAYERS,XLOC,TARO)
      TIN=INTEMP(NLI)
      CALL LININT(LOCAT,CO2M,LAYERS,XLOC,COMG)
      COMGIN=INO2(NLI)
      CALL LININT(LOCAT,UA,LAYERS,XLOC,UAMB)
      CALL LININT(LOCAT,VA,LAYERS,XLOC,VAMB)
      SALIN=INSAL(NLI)
!      CALL LININT(INXLOC,INTEMP,INNERLENGTH,XLOC,TIN) 
!      CALL LININT(INXLOC,INO2,INNERLENGTH,XLOC,COMGIN)
!      CALL LININT(INXLOC,INSAL,INNERLENGTH,XLOC,SALIN)

!      PRINT*,'TIN',TIN
!      PRINT*,'COMGIN',COMGIN
!      PRINT*,'SALIN',SALIN

      DOAMB=COMG
      CO2=COMG/32.
      COMGNOT=COMGP
!     At first temperature, DO and salinity in the outer plume is the same that in the top of the inner plume
      TPLUME=TIN
      COMGP=COMGIN
      SALPLU=SALIN  
      SALPLU=224  ! Eliminar JCT SAlinidad
	  SALIN = 224

	  PRINT*,'Init_input_frominner',QINTOP,BITOP,TPLUME

!     Guess the initial velocity 
      VGUESS=0.07
      VO=VGUESS
C      
C     AMBIENT AND AVERAGE WATER DENSITIES

      DENSEA=(0.059385*TARO**3-8.56272*TARO**2+65.4891*TARO)*0.001
     ++999.84298 !+(GAMMA)*SALAMB
!!      DENSEW=DENSEA
      DENSEP=(0.059385*TPLUME**3-8.56272*TPLUME**2+65.4891*TPLUME)*0.001
     ++999.84298 !+(GAMMA)*SALPLU
      DENSEINNER=(0.059385*TIN**3-8.56272*TIN**2+65.4891*TIN)*0.001
     ++999.84298 !+(GAMMA)*SALIN
      ! ELIMI SALINIDAD JCT2022



      ! DENSEA=(0.059385*TAMB**3-8.56272*TAMB**2+65.4891*TAMB)*0.001
      ! ++999.84298+(GAMMA)*SALARO
      ! DENSEW=(0.059385*TPLUME**3-8.56272*TPLUME**2+65.4891*TPLUME)*0.001
      ! ++999.84298+(GAMMA)*SALPLU
	 
	 
	 

!      PRINT*,'DENSEA',DENSEA
!      PRINT*,'TARO',TARO
!      PRINT*,'SALAMB',SALAMB
!      PRINT*,'DENSEP',DENSEP
!      PRINT*,'TPLUME',TPLUME
!      PRINT*,'SALPLU',SALPLU
!      PRINT*,'DENSEINNER',DENSEINNER
!      PRINT*,'TIN',TIN
!      PRINT*,'SALIN',SALIN
C
!     Solubility constant (mol/m3/Pa)
!!      HO2=(2.125-0.05023*TPLUME+5.7714E-4*TPLUME**2)/100000.
!!      HN2=(1.042-0.02457*TPLUME+3.1714E-4*TPLUME**2)/100000.
C	
C     Assume initial ambient dissolved nitrogen conc. equals saturated conc. at surface.
!!      CN2=(PATM*FRCNATM)*HN2      
!!      CNMG=CN2*28.0
!!      CNMGP=CNMG
!!      DNAMB=CNMG
!
!
C     --------------------------------------------------------------------------
C     PLUME PROPERTIES
C     --------------------------------------------------------------------------
!
!      QSCMS=QGFRAC*QSCFM/3.281**3/60.0
!      QGAS=PSTD*QSCMS*(TARO+273.15)/((PATM+DENSEA*G*Z)*TSTD)
!C     Initial bubble size.     
!      RB=DIAMM/2000.
!!     Bubble rise velocity
!      IF(RB.LE.(7.5E-4))THEN
!            VB=1189.0*RB**1.1945
!      ELSEIF(RB.GT.(7.5E-4).AND.RB.LT.(4.8E-3))THEN
!            VB=0.22
!      ELSE
!            VB=2.995*RB**0.489
!      ENDIF
!!     Mass transfer coefficient
!      KOLO=0.6*RB
!      IF(KOLO.GT.(4.0E-4))THEN
!            KOLO=4.0E-4
!      ENDIF
!      KOLN=KOLO
!
C     CALCULATION OF INITIAL WATER VELOCITY USING FROUDE NUMBER
!      VBUB=4./3.*PI*RB**3
!      N=QGAS/VBUB

      VDIFF=1
      !PRINT*,"QINTOP",QINTOP
      !PRINT*,"VO",VO
      !PRINT*,"BITOP",BITOP
      !PRINT*,"BO",BO
      !PRINT*,"DENSEA",DENSEA
      !PRINT*,"DENSEP",DENSEP
      !PRINT*,"VGUESS",VGUESS
      !PRINT*,"FRNO",FRNO


      !BITOP = 31.6
      ! PRINT*,"QINTOP",QINTOP
      ! PRINT*,"PI",PI
      ! PRINT*,"VO",VO
      ! PRINT*,"BITOP",BITOP
      ! PRINT*,"G",G
      ! PRINT*,"DENSEA",DENSEA
      ! PRINT*,"DENSEP",DENSEP


      DO WHILE (VDIFF.GT.1.0E-6)
         !BO=SQRT((QINTOP**2)/(PI*QINTOP*ABS(VO))+BITOP**2)
         ! BO=SQRT((QINTOP)/(PI*QINTOP*ABS(VO))+BITOP**2)
         BO=SQRT((QINTOP)/(PI*ABS(VO))+BITOP**2)
         BO=SQRT((QINTOP)/(PI*VO)+BITOP**2) ! JCT_2022
         VO=-FRNO*(ABS((BO-BITOP)*G*(DENSEA-DENSEP)/DENSEP))**0.5
         VDIFF=ABS(VO-VGUESS)
         VGUESS=VO
      END DO
      ! VO = -0.1021

!FroudeNoCriTerion = @(x) [x + Fo*...
!    sqrt(((sqrt(qwo/(pi*x)+bi^2))-bi)*g*abs((ra-ro)/ro));]; 

      ! PRINT*,"Tras calculo"
      ! PRINT*,"QINTOP",QINTOP
      ! PRINT*,"VOinicial",VO
      ! PRINT*,"BITOP",BITOP
      ! PRINT*,"BOinicial",BO
      ! PRINT*,"DENSEA",DENSEA
      ! PRINT*,"DENSEP",DENSEP
      ! PRINT*,"VGUESS",VGUESS
      ! PRINT*,"VO",VO
!     Water velocity in the top of the inner plume
      VI=-QINTOP/(PI*BITOP**2)
      ! PRINT*,'VI',VI
C     ------------------------------------------------------------------
C     VARIABLE TRANSFORMATION      
C     ------------------------------------------------------------------
!     TOTAL ENTRAIMENT: Shear and Vortex
!       E=2.*PI*BO*ALPHA*VO
!      PRINT*,"Shear" 
!      SHEARE=2.*PI*BO*ALPHA*VO
!     Vortex Entraiment: Hypothesis Projected Area Entraiment
!      VORTEXE=(2*BO*(SQRT(ABS(UAMB)*ABS(UAMB)+ABS(VAMB)*ABS(VAMB)))) 
!     a) Additive hypothesis
!         PRINT*,"Additive"  
!       E=SHEARE+VORTEXE
!     b) Maximun hypothesis
!         PRINT*,"Maximun"  
!      IF (SHEARE.GT.VORTEXE) THEN
! 	E=SHEARE
!      ELSE
!	E=VORTEXE
!      ENDIF

      BI=BITOP
	  
	  ! BI = 31.6085
	  ! VI = 0.0038
	  ! BO = 36.5030
	  ! VO = -0.0113
	  
      EI=+2.*PI*BI*ALPHAI*(VI+C1*VO)
      EO=-2.*PI*BI*ALPHAO*VO !JCT_2022 cambio signo
      EA=-2.*PI*BO*ALPHAA*VO !JCT_2022 cambio signo
      EP=QINTOP

      ! PRINT*,"EA,EO,EI,EP,VO",EA,EO,EI,EP, VO

      QW=QINTOP
      MOMENT=QW*VO
      PRINT*,"MOMENTinicial",MOMENT,VO,QW
	  
      FTEMP=QW*TPLUME 
      !FSAL=QW*(SALPLU*GAMMA/DENSE20)*DENSEP
      FSAL=QW*SALPLU !JCT2020_Sal
C     Previous equation corrected to account for salinity units conversion.        
      FDO=QW*CO2
!!      FDN=QW*CN2
!!      FGO=PSTD*QSCMS/(RGAS*TSTD)*FRACO
!!      FGONOT=FGO
!!      FGN=PSTD*QSCMS/(RGAS*TSTD)*FRACN
C     Revised gaseous flux equations.
!!      YO2=FGO/((PI*(LAMBDA*B)**2)*(V+VB))
!!      YN2=FGN/((PI*(LAMBDA*B)**2)*(V+VB))           
!!      PZ=PATM+(DENSEA*G*Z)
!!      PO=PZ*FRACO
!!      PN=PZ*FRACN
!     Initialize lateral withdrawal flowrate for first/lowest cell in column/segment
      JJ=0
      QWDO(LAYTOP)= QW 
      PWDO(LAYTOP)= (2*PI*BO)
      AWDO(LAYTOP)= (PI*(BO**2-BI**2))    
      BWDO(LAYTOP)= BO 
      ! PRINT*,"Init_input",QWDO(LAYTOP),VO,BO,TPLUME

!
!     -------------------------------------------------------------------------- 
C	SOLUTION PROCEEDURE
!     --------------------------------------------------------------------------
!
      DZ=0.001
      ! H=0.001
      H=-0.001
      NELS = 0
      HWITH=0.0
      NLO=NLI+1
!
      M=0

!        PRINT*,"VO",VO
!        PRINT*,"Z",Z
!        PRINT*,"DEPTH",DEPTH
!        PRINT*,"DEPTHDMPR",DEPTHDMPR
!        PRINT*,"QW",QW
!        PRINT*,"EI",EI
!        PRINT*,"EA",EA
!        PRINT*,"EO",EO

!         PRINT*,"---------------------------------------"
!         PRINT*,"QW",QW
!         PRINT*,"MOMENT",MOMENT
!         PRINT*,"FTEMP",FTEMP
!         PRINT*,"FSAL",FSAL
!         PRINT*,"FDO",FDO
!         PRINT*,"VO",VO

!         PRINT*,"---------------------------------------"
!         PRINT*,"DEPTHDMPR",DEPTHDMPR
!         PRINT*,"X",X
!         PRINT*,"DZ",DZ
!         PRINT*,"---------------------------------------"


!        PRINT*,"SOLUTION PROCEEDURE"
         !PRINT*,"DO WHILE pre", VO, '-1.E-6', Z, DEPTH 
		 
		 
		         ! PRINT*,"INQW(NLO)",INQW(NLO)
		         ! PRINT*,"INMOMENT(NLO)",INMOMENT(NLO)


       DO WHILE (VO.LT.-1.E-6.AND.Z.LT.DEPTH)
         M=M+1
 !        PRINT*,"DO WHILE", VO, '-1.E-6.', Z, DEPTH 
 !       PRINT*,"Estoy dentro"
!        PRINT*,"VO",VO
         Z=Z+DZ
         X=X+DZ
         ELEV=ELEV-DZ
         NELS = NELS + 1
         NLO=NLO-1
C	
C        Interpolate input profiles to obtain line plume boundary conditions
         XLOC=DEPTHDMPR+X

         CALL LININT(LOCAT,TE,LAYERS,XLOC,TARO)
         TIN=INTEMP(NLO)
         CALL LININT(LOCAT,CO2M,LAYERS,XLOC,COMG)
         COMGIN=INO2(NLO)
         SALIN=INSAL(NLO)
		 SALIN = 224
         CALL LININT(LOCAT,UA,LAYERS,XLOC,UAMB)
         CALL LININT(LOCAT,VA,LAYERS,XLOC,VAMB)
C        Inner plume width
         BI=INBI(NLO)
C        Inner velocity plume
         INNERQW=INQW(NLO)
         INNERMOMENT=INMOMENT(NLO)
         VI=INNERMOMENT/INNERQW
		 
		 ! PRINT*,"BI INNERQW INNERMOMENT VI", BI, INNERQW, INNERMOMENT, VI
         ! PRINT*,"EA,EO,EI,EP,VO",EA,EO,EI,EP, VO

		 
!         VI=-INNERMOMENT/INNERQW


!         CALL LININT(INXLOC,INTEMP,INNERLENGTH,XLOC,TIN) 
!         CALL LININT(INXLOC,INO2,INNERLENGTH,XLOC,COMGIN)
!         CALL LININT(INXLOC,INSAL,INNERLENGTH,XLOC,SALIN) 
!         CALL LININT(INXLOC,INBI,INNERLENGTH,XLOC,BI)
!         CALL LININT(INXLOC,INQW,INNERLENGTH,XLOC,INNERQW) 
!         CALL LININT(INXLOC,INMOMENT,INNERLENGTH,XLOC,INNERMOMENT)          
C                  

C        Use subroutines for Runge Kutta method solution
         NEQN=5
         Y(1)=QW
         Y(2)=MOMENT
         Y(3)=FTEMP
         Y(4)=FSAL
         Y(5)=FDO
!!         Y(6)=FDN


         ! PRINT *, 'MOMENT_pre',QW,MOMENT
		 
         !PRINT*,"QWantes",QW
         CALL DERIVS_5(EA,EO,EI,VI,VO,DENSEA,DENSEP,G,PI,BO,BI,TARO,TIN,
     +             TPLUME,SALAMB,SALIN,SALPLU,GAMMA,DENSE20,DENSEINNER,
     +             GAMMA1,DOAMB,COMGIN,COMGP,DNAMB,CNMGP,Z,Y,DYDX)

         CALL RK4_5(EA,EO,EI,VI,VO,DENSEA,DENSEP,G,PI,BO,BI,TARO,TIN,
     +             TPLUME,SALAMB,SALIN,SALPLU,GAMMA,DENSE20,DENSEINNER,
     +             GAMMA1,DOAMB,COMGIN,COMGP,DNAMB,CNMGP,Y,DYDX,NEQN,X,
     +             H,YOUT)
	       


          ! PRINT *, 'DENSEA',DENSEA
          ! PRINT *, 'DENSE20',DENSE20
          ! PRINT *, 'DENSEINNER',DENSEINNER
		  ! PRINT *, 'EA',EA
          ! PRINT *, 'EO',EO
          ! PRINT *, 'EI',EI
		  ! PRINT *, 'VO',VO
          ! PRINT *, 'VI',VI
		  ! PRINT *, 'BO',BO
          ! PRINT *, 'BI',BI

C       
         QW=YOUT(1)
         MOMENT=YOUT(2)
         FTEMP=YOUT(3)
         FSAL=YOUT(4)
         FDO=YOUT(5)
!         FDN=YOUT(6)
		 ! PRINT *, 'MOMENT_post',QW,MOMENT

		 
		 ! PRINT*,"Analisis EQ momento"
		 ! PRINT*,"GAMMA1",GAMMA1
		 ! PRINT*,"PI",PI
		 ! PRINT*,"G",G
		 ! PRINT*,"BO",BO
		 ! PRINT*,"BI",BI
		 ! PRINT*,"DENSEP",DENSEP
		 ! PRINT*,"DENSEA",DENSEA
		 ! PRINT*,"DENSE20",DENSE20
		 ! PRINT*,"EI",EI
		 ! PRINT*,"VO",VO
		 ! PRINT*,"EO",EO
		 ! PRINT*,"VI",VI
	 
!		       DYDX(2)=((1/GAMMA1)*(-PI*G*(BO**2-BI**2)*
!     +((DENSEP-DENSEA)/DENSE20))-EI*VO+EO*VI)
	 
	 



          ! OPEN (UNIT=59, FILE="outer_plume_eq.txt", POSITION="APPEND")
             ! WRITE (UNIT=59, FMT='(6F12.6)') QW,MOMENT,FTEMP,FSAL,FDO,X
          ! CLOSE (UNIT=59)

		  ! PRINT*,"QW, MOMENT,VO",QW,MOMENT,VO


         IF(MOMENT.LT.0.0)THEN
            TPLUME=FTEMP/QW
	    !SALPLU=FSAL/(QW*DENSEP)/(GAMMA/DENSE20)
            SALPLU=224  ! Eliminar JCT SAlinidad
C           Previous equation corrected to consistently express salinity in uS/cm	   
	        CO2=FDO/QW
            COMGP=CO2*32.
            VO=MOMENT/QW
            AREA=ABS(QW/VO)
            BO=SQRT((AREA+PI*BI**2)/PI)
!           Save outer plume information
            OUTPLUME(NLO,1)=XLOC 
            OUTPLUME(NLO,2)=BO          
            OUTPLUME(NLO,3)=QW
            OUTPLUME(NLO,4)=MOMENT
            OUTPLUME(NLO,5)=TPLUME
            OUTPLUME(NLO,6)=COMGP
            OUTPLUME(NLO,7)=SALPLU
            OUTPLUME(NLO,8)=VO
            PRINT*,"me voy por momento"
			PRINT*,"QW, MOMENT,VO",QW,MOMENT,VO

	    GOTO 20
         ENDIF

         VO=MOMENT/QW
         AREA=ABS(QW/VO)
         BO=SQRT((AREA+PI*BI**2)/PI)
		 
		 ! PRINT*,"QW, MOMENT,VO",QW,MOMENT,VO
		 ! PRINT*,"TPLUME, TARO,TIN",TPLUME,TARO,TIN



        !PRINT*,"QW",QW
!        PRINT*,"B",B
!        TOTAL ENTRAIMENT: Shear and Vortex
!          E=2.*PI*B*ALPHAI*V
!         PRINT*,"Shear" 
!         SHEARE=2.*PI*B*ALPHAI*V
!        Vortex Entraiment: Hypothesis Projected Area Entraiment
!         VORTEXE=(2*B*(SQRT(ABS(UAMB)*ABS(UAMB)+ABS(VAMB)*ABS(VAMB)))) 
!	 a) Additive hypothesis 
!          PRINT*,"Additive" 
!      	     E=SHEARE+VORTEXE
!	 b) Maximun hypothesis
!          PRINT*,"Maximun" 
!	    IF (SHEARE.GT.VORTEXE) THEN
!	 	E=SHEARE
!            ELSE
!		E=VORTEXE
!            ENDIF
!         


         EI=+2.*PI*BI*ALPHAI*(VI+C1*VO)
         !EO=-2.*PI*BI*ALPHAO*VO
         !EA=-2.*PI*BO*ALPHAA*VO
         EO=-2.*PI*BI*ALPHAO*VO ! JCT_2020 VO es negativa
         EA=-2.*PI*BO*ALPHAA*VO ! JCT_2020

         

      ! PRINT*,"EA,EO,EI,VI,VO",EA,EO,EI,VI,VO
!         IF (M.EQ.500) THEN
!!            PRINT*, XLOC,B,V,UAMB,VAMB,SHEARE,VORTEXE
!             WRITE (UNIT=50, FMT='(7F12.6)') XLOC,BO,VO,UAMB,VAMB,SHEARE,VORTEXE
!             M=O
!         ENDIF
 !     PRINT*,'FLAG_6'
!
!        Temperatura and salinity in the plume    

      ! PRINT*,"FTEMP,QW,TPLUME",FTEMP,QW,TPLUME

         TPLUME=FTEMP/QW
         SALPLU=FSAL/(QW*DENSEP)/(GAMMA/DENSE20)
         SALPLU=224  ! Eliminar JCT SAlinidad
C        Previous equation corrected to consistently express salinity in uS/cm

!        Dissolved oxygen and nitrogen concentration
         CO2=FDO/QW
!         CN2=FDN/QW 
         COMGP=CO2*32.
!         CNMGP=CN2*28.

!        Save outer plume information
         OUTPLUME(NLO,1)=XLOC 
         OUTPLUME(NLO,2)=BO          
         OUTPLUME(NLO,3)=QW
         OUTPLUME(NLO,4)=MOMENT
         OUTPLUME(NLO,5)=TPLUME
         OUTPLUME(NLO,6)=COMGP
         OUTPLUME(NLO,7)=SALPLU
         OUTPLUME(NLO,8)=VO
!        Add incremental entrainment to total cell entrainment/withdrawal   
         !QWDO(LAYTOP+JJ)=QWDO(LAYTOP+JJ)+(EA+EO-EI)*DZ
         QWDO(LAYTOP+JJ)=QWDO(LAYTOP+JJ)-(EA+EO-EI)*DZ ! JCT_2022
         PWDO(LAYTOP+JJ)=PWDO(LAYTOP+JJ)+2.*PI*BO
         AWDO(LAYTOP+JJ)=AWDO(LAYTOP+JJ)+(PI*(BO**2-BI**2))
         BWDO(LAYTOP+JJ)=BWDO(LAYTOP+JJ)+BO
         HWITH=HWITH+DZ   
          ! PRINT*,"JJ_",JJ 
         IF(HWITH.GT.HCELL)THEN
            PWDO(LAYTOP+JJ) = PWDO(LAYTOP+JJ)/NELS
            AWDO(LAYTOP+JJ) = AWDO(LAYTOP+JJ)/NELS
            BWDO(LAYTOP+JJ) = BWDO(LAYTOP+JJ)/NELS
            ! PRINT*,"JJ",JJ 
            JJ=JJ+1
            HWITH=0.0
            QWDO(LAYTOP+JJ)=0.0
            PWDO(LAYTOP+JJ)=0.0
            BWDO(LAYTOP+JJ)=0.0
            AWDO(LAYTOP+JJ)=0.0
            NELS = 0
         ENDIF

!
C        Revised gaseous flux equations.
!!         YO2=FGO/((PI*(LAMBDA*B)**2)*(V+VB))
!!         YN2=FGN/((PI*(LAMBDA*B)**2)*(V+VB))
C  
!        
!!         PZ=PATM+(DENSEA*G*Z)
!!         QGAS=(FGO+FGN)*RGAS*(TPLUME+273.15)/PZ
!!         VBUB=QGAS/N
C
!        GAS VOLUME PER TOTAL VOLUME OF THE BUBBLE-WATER MIXTURE IN THE INNER CORE OF THE PLUME
!!         VG=VBUB*N/((V+VB)*(PI*(LAMBDA*B)**2))
C        Previous equation revised to account for correct plume cross-sectional area occupied by bubbles. 
!        Bubbles radius     
!!         RB=(3.*QGAS/(4.*PI*N))**(1./3.)
!!         IF(RB.LT.0.0)THEN
!!           RB=1.0E-8
!!         ENDIF
!!         FRACO=FGO/(FGO+FGN)
!!         FRACN=1.0-FRACO
C	
!        Partial Pressure
!!         PO=PZ*FRACO
!!         PN=PZ*FRACN
!        Density (ambient water, water in plume and of the plume)
      DENSEA=(0.059385*TARO**3-8.56272*TARO**2+65.4891*TARO)*0.001
     ++999.84298 !+(GAMMA)*SALAMB
      DENSEP=(0.059385*TPLUME**3-8.56272*TPLUME**2+65.4891*TPLUME)*0.001
     ++999.84298 !+(GAMMA)*SALPLU
      DENSEINNER=(0.059385*TIN**3-8.56272*TIN**2+65.4891*TIN)*0.001
     ++999.84298 !+(GAMMA)*SALIN

         ! PRINT*,"DENSEP,DENSEA,DENSEINNER",DENSEP,DENSEA,DENSEINNER
         ! PRINT*,"SALPLU,SALAMB,SALIN",SALPLU,SALAMB,SALIN

C        Previous equation re-revised to account for correct salinity units (uS/cm) in density calculations.      
C
C        BUBBLE PROPERTIES
!
!        Bubble rise velocity
!!         IF(RB.LE.(7.5E-4))THEN
!!            VB=1189.0*RB**1.1945
!!         ELSEIF(RB.GT.(7.5E-4).AND.RB.LT.(4.8E-3))THEN
!!            VB=0.22
!!         ELSE
!!            VB=2.995*RB**0.489
!!         ENDIF
C
!        Mass transfer coeffitient
!!         KOLO=0.6*RB
!!         IF(KOLO.GT.(4.0E-4))THEN
!!            KOLO=4.0E-4
!!         ENDIF
!!         KOLN=KOLO
C
!        Solubility constant (mol/m3/Pa)
!!         HO2=(2.125-0.05023*TPLUME+5.7714E-4*TPLUME**2)/100000.
!!         HN2=(1.042-0.02457*TPLUME+3.1714E-4*TPLUME**2)/100000.
C
!        Froude number
!!         FR=V/(2.*LAMBDA*BO*G*(DENSEA-DENSEP)/DENSEP)**0.5
!
!!         DCO2=HO2*PO-CO2
C      
  !       PRINT*,"QW",QW
!         PRINT*,"VI",VI 
!         PRINT*,"VO",VO 
  !       PRINT*,"Z",Z
  !       PRINT*,"MOMENT",MOMENT

  !       PRINT*,"TIN",TIN
  !       PRINT*,"COMGIN",COMGIN
  !       PRINT*,"SALIN",SALIN
!         PRINT*,"BI",BI
!         PRINT*,"BO",BO
  !       PRINT*,"INNERQW",INNERQW
  !       PRINT*,"NLO",NLO

      END DO
C
!     --------------------------------------------------------------------
C     CALCULATION OF AVERAGE NET OXYGEN MASS TRANSFER FOR DAY
!     --------------------------------------------------------------------
      !GROSSMT=(FGONOT-FGO)*32./1000.*86400.
      !OTEFF=(FGONOT-FGO)/FGONOT*100.
      !DELTAC=COMGP-COMGNOT
C
!      ELEVT=ELEV
   20 DEPTHINTR=ELEV
      QWDET=QW
      TPLUMED=TPLUME
      COMGPD=COMGP
!!      LAYTOP=LAYDIFF-JJ
      LAYINTR=LAYTOP+JJ
!!      BTOP=B
!      PRINT*,"QWDET",QWDET
!      PRINT*,"OUTER PLUME SALIDAS"
!      PRINT*,"ELEVT",ELEVT
!      PRINT*,"QWT",QWT
!      PRINT*,"TPLUMED",TPLUMED
!      PRINT*,"COMGPD",COMGPD
!      PRINT*,"LAYTOP_outer",LAYTOP
 !     PRINT*,"LAYDIFF_outer",LAYDIFF
!      PRINT*,"LAYINTR_outer",LAYINTR
!      PRINT*,"JJ_outer",JJ
!      PRINT*,"BTOP",BTOP
!      PRINT*,"INPLUMEDEPTH",INPLUME(:,1)
      RETURN
      END
	  
	  
C
C------------------------------------------------------------------------------
C                          
             
C     *******************************************************************************
	SUBROUTINE INNER_PLUME_RECT(YEAR,JULDAY,WSEL,DIFFEL,LAYERS,LAMBNOT,
     +SALAMB,PATM,DIAMM,QSCFM,FRCONOT,LAYDIFF,HCELL,LOCAT,TE,CO2M,UA,VA,
     +ELEVT,QWT,TPLUMET,COMGPT,SALPLUMET,QWDI,
     +BWDI,LAYTOP,LAYINTR,DEPTHINTR,NLI,NLO,NITERPLUME,
     +INPLUME,OUTPLUME,ALPHAI,ALPHAO,ALPHAA,LAMBDA,FRNI,
     +FRNO,GAMMA1,LDIFF,OTEFF,LWDI)
C     *******************************************************************************

      IMPLICIT NONE

C	THIS SUBROUTINE IS WRITTEN TO PREDICT THE PERFORMANCE OF A CIRCULAR BUBBLE PLUME.  
C     THE MODEL IS BASED ON THE WUEST ET AL. (1992) CIRCULAR BUBBLE PLUME MODEL.
C     By:
C	     
C	VERSION 3 (to couple with Francisco Rueda's reservoir model and for Lake Hallwil real diffuser)
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
C	  ALPHAI=ENTRAINMENT COEFFICIENT INNER PLUME  (-)
C	  ALPHAO=ENTRAINMENT COEFFICIENT OUTER PLUME  (-)
C	  ALPHAA=ENTRAINMENT COEFFICIENT FROM AMBIENT (-)
C     BI=1/2 INNER PLUME WIDTH (m)
!     BIFF=DIFFUSER RADIUS (m) 
C     BAVG=AVERAGE 1/2 DIFFUSER WIDTH (m)
!         C1=COEFFICIENT OF DIFERENT COUNTERFLOW ENTRAINMENTS
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
!         DEPTHINTR=INTRUSION DEPTH OF THE DOUBLE PLUME
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
!	FR=FROUDE NUMBER (-)
C     FRCNOT=FRACTION OF NITROGEN IN ATMOSPHERE (-)
C     FRNO=INITIAL FROUDE NUMBER (-)
C     GAMMA=SALINITY CONVERSION FACTOR [kg/m3/(uS/cm)]
C     GROSSMT=GROSS MASS TRANSFER OF OXYGEN FROM PLUME (kg/d)
!	HO2=SOLUBILITY CONSTANT FOR OXYGEN (mol/m3/Pa)
!	HN2=SOLUBILITY CONSTANT FOR NITROGEN (mol/m3/Pa)
C     HCELL=HEIGHT OF CELL IN GRID (m)
C	HO=SOLUBILITY CONSTANT FOR OXYGEN (mol/m3/Pa)
C	HN=SOLUBILITY CONSTANT FOR NITROGEN (mol/m3/Pa)
C     HWITH=HEIGHT OF WITHDRAWAL/ENTRAINMENT ZONE (m)
!         INPLUME(:,6): INNER PLUME CHARACTERISTIC [HEIGTH; FLOWRATE OF WATER; MOMENTUM; TEMPERATURE; SALINITY; O2] JCT
C     JULDAY=JULIAN DAY IN GIVEN YEAR
C	KOLO=MASS TRANSFER COEFFICIENT FOR OXYGEN (m/s)
C	KOLN=MASS TRANSFER COEFFICIENT FOR NITROGEN (m/s)
C     LAKE=LAKE AND DIFFUSER TYPE (SHR AND LINEAR=1 OR AMISK AND RECTANGULAR=2) FOR SELECTION OF LAMBDA
C     LAYERS=NUMBER OF LAYERS/DATA POINTS IN BOUNDARY CONDITION PROFILES (-)
C     LAYDIFF=GRID LAYER CORRESPONDING TO DIFFUSER DEPTH (-)
!         LAYINTR=GRID LAYER CORRESPONDING TO INTRUSION DEPTH OF THE DOUBLE PLUME (-)
C     LAYTOP=GRID LAYER CORRESPONDING TO TOP OF PLUME (-)
C     LDIFF=LENGTH OF DIFFUSER (m)
C     LNOT=INITIAL DIFFUSER LENGTH (m)
C	LAMBDA=FRACTION OF PLUME OCCUPIED BY BUBBLES (-)
C     LAMBNOT=LAMBDA x INITIAL PLUME RADIUS;EQUAL TO DIFFUSER RADIUS (m) 
C     LOCAT=DEPTHS FOR INPUT BOUNDARY CONDITION PROFILES (m) 
C	MOMENT=MOMENTUM (m4/s)
C	N=NUMBER OF BUBBLES PER SECOND (1/s)
!         NITERPLUME= NUMBER OF ITERATION
!         NLI
!         NLO
C     OTEFF=OXYGEN TRANSFER EFFICEINCY (%)
!         OUTPLUME(:,6)= OUTER PLUME CHARACTERISTIC [HEIGTH; FLOWRATE OF WATER; MOMENTUM; TEMPERATURE; SALINITY; O2] JCT
!         OUTERMOMENT= MOMENTUM IN THE OUTER PLUME 
!         OUTERQW=FLOWRATE OF WATER IN THE OUTER PLUME
!         OUTXLOC
!         OUTQW
!         OUTXMOMENT
!         OUTXTEMP
!         OUTO2  
C     PATM=ATMOSPHERIC PRESSURE AT AVERAGE WSEL (Pa)
C	PSTD=STANDARD PRESSURE (Pa)
!     QGAS=TOTAL GAS FLOW TO DIFFUSER
C	QSCFM=STANDARD GAS FLOW RATE (scfm), TOTAL STANDAR GAS FLOW RATE TO DIFFUSER 
C	QSCMS=STANDARD GAS FLOW RATE (scms)
C	QW=FLOWRATE OF WATER (m3/s)
C     QWT=TOTAL DETRAINMENT FLOW RATE AT THE TOP OF THE PLUME (m3/s)  
C	RB=BUBBLE RADIUS (m)
C	RGAS=IDEAL GAS CONSTANT (J/mol/K)
C     SAL=SALINITY (uS/cm)
C     SALAMB=SALINITY OF AMBIENT WATER (uS/cm)
C	SALPLU=SALINITY OF THE PLUME (uS/cm)
!         SALPLUMET=DETRAINMENT PLUME SALINITY AT THE TOP OF THE PLUME (uS/cm)
!     SHEARE=SHEAR ENTRAINMENT
C	TAMB=AMBIENT WATER TEMPERATURE (C)
C	  TARO=AROUND WATER TEMPERATURE (C)
C     TAVG=AVERAGE AMBIENT WATER TEMPERATURE (C)
C     TDS=TOTAL DISSOLVED SOLIDS (g/m3) [0.64 conversion factor from Chapra book]
C     TE=TEMPERATURE PROFILE FOR INPUT BOUNDARY CONDITION (C) 
C	TPLUME=PLUME TEMPERATURE (C)
C     TPLUMET=DETRAINMENT PLUME TEMPERATURE AT THE TOP OF THE PLUME (C) 
C	TSTD=STANDARD TEMPERATURE (K)
!     UA = U AMBIENT VELOCITY PROFILE FOR INPUT BOUNDARY CONDITION
C     UAMB = U Ambient Velocity (m) used to calculate VORTEX entrainment (FJR
C	VI= INNER PLUME WATER VELOCITY (m/s)
!     VA = V AMBIENT VELOCITY PROFILE FOR INPUT BOUNDARY CONDITION
C     VAMB = V Ambient Velocity (m) used to calculate VORTEX entrainment (FJR
C     VAVG=AVERAGE WATER VELOCITY (m/s)
C	VB=BUBBLE RISE VELOCITY (m/s)
C	VBUB=BUBBLE VOLUME (m3)
!     VG=GAS VOLUME PER TOTAL VOLUME OF THE BUBBLE-WATER MIXTURE IN THE INNER CORE OF THE PLUME
C	VGUESS=GUESSED INITIAL WATER VELOCITY (m/s)
!         VO=WATER VELOCITY IN THE OUTER PLUME (m/s)
!     VORTEXE=VORTEX ENTREAINMENT
C     WSEL=WATER SURFACE ELEVATION (m)
C	YO2=GASEOUS OXYGEN CONCENTRATION (mol/m3)
C	YN2=GASEOUS NITROGEN CONCENTRATION (mol/m3)
C	Z=DEPTH TO DIFFUSER (m)
C     PWD  = Perimeter (m) (FJR)
C     VUP  = Upward velocity (m/2)
C     AWD  = Area (m2) of plume
C     BWD  = Width 

!    GAMMA1= momentum amplification factor
!  DOARO
C
      REAL*8 AREA,BI,LI,CO2,COMG,CN2,CNMG,DS,DZ,DENSEA,DENSEP,DENSEW,
     +DIAMM,FDO,FDN,FRACO,FRACN,FSAL,FTEMP,FGO,FGN,G,GAMMA,HO2,HN2,
     +KOLN,KOLO,LAMBDA,MOMENT,N,PI,PO,PN,PSTD,PZ,QSCFM,QSCMS,QW,QGAS,
     +RB,RGAS,SALAMB,SALPLU,TARO,TPLUME,TSTD,VI,VB,VBUB,VG,VGUESS,
     +YO2,YN2,Z,AA,BB,CC,BNOT,LNOT,ELEV,DT,TE(1000),XLOC,DEPTH,
     +CO2M(1000),LOCAT(1000),DOAMB,COMGP,CNMGP,WSEL,PATM,SAL(1000),
     +DIFFEL,GROSSMT,FGONOT,DNAMB,DMPR,SUMSAL,
     +LAMBNOT,FRCONOT,H,DYDX(8),Y(8),YOUT(8),FRCNATM,
     +DENSE20,OTEFF,FRNI,FRNO,VDIFF,FR,BUOY,DCO2,QGFRAC,LDIFF,TDS(1000),
     +JDAY,EL(70),DELTAC,COMGNOT,ELEVT,QWT,TPLUMET,COMGPT,
     +PWDI(500),BWDI(500),AWDI(500),QWDI(500),X,SHEARE,VORTEXE,
     +HWITH,HCELL,JULDAY,BTOP,UA(1000),VA(1000),UAMB,VAMB,
     +INPLUME(35000,9),OUTPLUME(35000,9),SALARO,DOARO,AU,
     +DEPTHINTR,ALPHAI,ALPHAO,ALPHAA,SALPLUMET,GAMMA1,TAMB,
     +C1,VO,EI,EO,OUTERQW,OUTERMOMENT,OUTXLOC(35000),OUTQW(35000),
     +OUTMOMENT(35000),OUTTEMP(35000),OUTO2(35000),V,OUTSAL(35000)
	  REAL*8 LWDI(500)
      INTEGER II,IJ,IK,IN,JJ,LL,NEQN,NN,MI,JK,JL,LAYTOP,
     +LAYERS,KM,KN,KO,KP,ROWS,KQ,KR,KU,KV,KW,KX,KY,KZ,KS,LAYDIFF,YEAR 
      INTEGER NELS,ierror,M,CONT,
     +LAYINTR,
     +NITERPLUME,NLI,NLO,NLIO,k
C 
C     --------------------------------------------------------------------------
C     CONSTANTS
C     --------------------------------------------------------------------------
!

      PRINT*,'----INNER_PLUME_RECT----'

      G=9.80665
      GAMMA=6.9E-4
      PI=ACOS(-1.0)
      PSTD=101325.
      RGAS=8.314
      TSTD=293.15
      DENSE20=998.2
      FRCNATM=0.0
      QGFRAC=1.0
!      C1=-1
      C1=0 !Socolofsky et al. 2008
!     
C     --------------------------------------------------------------------------
C     PARAMETERS AND INITIALIZE THE VARIABLES 
C     --------------------------------------------------------------------------
C      
      PRINT*,'......PARAMETERS AND INITIALIZE THE VARIABLES.....'

      INPLUME(:,:)=0.0

!     ...Geometric caracteristic 
      DEPTH=WSEL-DIFFEL  
      Z=DEPTH
      ELEV=DIFFEL
!        
C     Assume that gas bubbles are composed of oxygen and nitrogen only.
      FRACO=FRCONOT 
      FRACN=1.0-FRACO
!
C     Interpolate input profiles to obtain the plume initial conditions
      X=0.
      XLOC=DEPTH-X

      CALL LININT(LOCAT,TE,LAYERS,XLOC,TAMB)
      CALL LININT(LOCAT,CO2M,LAYERS,XLOC,COMG)
      CALL LININT(LOCAT,UA,LAYERS,XLOC,UAMB)
      CALL LININT(LOCAT,VA,LAYERS,XLOC,VAMB)
      COMGP=COMG
      DOARO=COMG
      CO2=COMG/32.
      COMGNOT=COMG
!     At first the temperature and the salinity un the plume is the same that in the ambient water
      SALPLU=SALAMB    
      TPLUME=TAMB
C     ... Guess the initial velocity 
      VGUESS=0.07
      VI=VGUESS
!     INITIAL DIFFUSER SIZE
!     At first the plume radius BNOT (initial plume HALF WIDTH) in the TOP HAT models 
      BNOT=LAMBNOT/LAMBDA
      BI=BNOT
      !LNOT=LDIFF+2.0*BNOT*(1.0-LAMBDA) ! 
      LNOT=LDIFF
      LI=LNOT

C     AMBIENT AND AVERAGE WATER DENSITIES
      DENSEA=(0.059385*TAMB**3-8.56272*TAMB**2+65.4891*TAMB)*0.001
     ++999.84298 !+(GAMMA)*SALAMB !JCT_2022
      DENSEW=DENSEA

C
!     Solubility constant (mol/m3/Pa)
      HO2=(2.125-0.05023*TPLUME+5.7714E-4*TPLUME**2)/100000.
      HN2=(1.042-0.02457*TPLUME+3.1714E-4*TPLUME**2)/100000.

      PRINT*, 'DEPTHinner, Z, ELEV,XLOC, COMGP, DOARO, CO2, TAMB, UAMB, VAMB',DEPTH, Z, ELEV, XLOC, COMGP, DOARO, CO2, TAMB, UAMB, VAMB
      PRINT*, 'SALPLU, TPLUME, VI, BI, LI, DENSEA, HO2', SALPLU, TPLUME, VI, BI, LI, DENSEA,HO2

C	
C     Assume initial ambient dissolved nitrogen conc. equals saturated conc. at surface.
      CN2=(PATM*FRCNATM)*HN2      
      CNMG=CN2*28.0
      CNMGP=CNMG
      DNAMB=CNMG
!
!     Outer plume
      OUTXLOC=OUTPLUME(:,1)
      OUTQW=OUTPLUME(:,3)
      OUTMOMENT=OUTPLUME(:,4)
      OUTTEMP=OUTPLUME(:,5)
      OUTO2=OUTPLUME(:,6)
      OUTSAL=OUTPLUME(:,7)
!
C     --------------------------------------------------------------------------
C     BUBBLE PROPERTIES
C     --------------------------------------------------------------------------
!
      PRINT*,'..........BUBBLE PROPERTIES............'
      QSCMS=QGFRAC*QSCFM/3.281**3/60.0
      QGAS=PSTD*QSCMS*(TAMB+273.15)/((PATM+DENSEA*G*Z)*TSTD)

C     Initial bubble size.     
      RB=DIAMM/2000.
!     Bubble rise velocity
      IF(RB.LE.(7.5E-4))THEN
            VB=1189.0*RB**1.1945
      ELSEIF(RB.GT.(7.5E-4).AND.RB.LT.(4.8E-3))THEN
            VB=0.22
      ELSE
            VB=2.995*RB**0.489
      ENDIF

!     Mass transfer coefficient (KOLO = beta0)
      KOLO=0.6*RB
      IF(KOLO.GT.(4.0E-4))THEN
            KOLO=4.0E-4
      ENDIF
      KOLN=KOLO
      PRINT*, 'QSCMS, QGAS, RB, VB, KOLN', QSCMS, QGAS, RB, VB, KOLN
!
C     --------------------------------------------------------------------------
C     CALCULATION OF INITIAL WATER VELOCITY USING FROUDE NUMBER
C     --------------------------------------------------------------------------
      
      PRINT*,'.............CALCULATION OF INITIAL WATER VELOCITY USING FROUDE NUMBER........'

C     ... Bubbles
      VBUB=4./3.*PI*RB**3
      N=QGAS/VBUB
      VDIFF=1
      DO WHILE (VDIFF.GT.1.0E-6)
        VG=QGAS/((VGUESS+VB)*(2.*LAMBDA*BI*LI)) 
        !VG=QGAS/((VGUESS+VB)*(PI*(LAMBDA*BI)**2))
        DENSEP=(1.0-VG)*DENSEW
        VI=FRNI*(2.0*LAMBDA*BI*G*(DENSEA-DENSEP)/DENSEP)**0.5
        VDIFF=ABS(VI-VGUESS)
        VGUESS=VI
      ENDDO      

      PRINT*, 'VG, VI, N, DENSEP', VG, VI, N, DENSEP

C     ------------------------------------------------------------------
C     VARIABLE TRANSFORMATION      
C     ------------------------------------------------------------------

      PRINT*,'........VARIABLE TRANSFORMATION.......... '
      VO = 0.00
      EI=2.*LI*ALPHAI*(VI+C1*VO)          ! ACC 2026: solo lados largos (2D linea fuente)
      EO=-2.*LI*ALPHAO*VO
      
      QW=VI*(2.*LI*BI) 
      MOMENT=(2.*LI*BI)*VI**2 
      FTEMP=QW*TPLUME 
      FSAL=QW*SALPLU 
C     Previous equation corrected to account for salinity units conversion.        
      FDO=QW*CO2
      FDN=QW*CN2
      FGO=PSTD*QSCMS/(RGAS*TSTD)*FRACO
      FGONOT=FGO
      FGN=PSTD*QSCMS/(RGAS*TSTD)*FRACN

C     Revised gaseous flux equations.
      YO2=FGO/((2.*LAMBDA*BI*LI)*(VI+VB))  ! ACC 2026
      YN2=FGN/((2.*LAMBDA*BI*LI)*(VI+VB))  ! ACC 2026
      PZ=PATM+(DENSEA*G*Z)
      PO=PZ*FRACO
      PN=PZ*FRACN
      BUOY=(G*(DENSEA-DENSEP)/DENSEP*QW)/LNOT
      TDS=SALPLU*0.64

!     Initialize lateral withdrawal flowrate for first/lowest cell in column/segment
      JJ=0
      QWDI(LAYDIFF)= QW 
!     ACC 2026: perimetro = solo lados largos (Dissanayake 2021). 
      PWDI(LAYDIFF)= 2.*LI            ! ACC 2026: solo lados largos
      AWDI(LAYDIFF)= LI*BI   
      BWDI(LAYDIFF)= BI 
	  LWDI(LAYDIFF)= LI 
!

      PRINT*,'EI, QW, MOMENT, FTEMP, FDO, FGO, VO, PZ, PO', EI, QW, MOMENT, FTEMP, FDO, FGO, VO, PZ, PO
!     -------------------------------------------------------------------------- 
C	SOLUTION PROCEDURE
!     --------------------------------------------------------------------------
!
      PRINT*,'........SOLUTION PROCEDURE.......... '

      DZ=0.01
      H=0.01
      NELS = 0
      HWITH=0.0
      NLI=0
      NLIO=NLO
!
      M=0
      DO WHILE (VI.GT.1.E-6.AND.Z.GT.0.0)
         M=M+1
         Z=Z-DZ
         X=X+DZ
         ELEV=ELEV+DZ
         NELS = NELS + 1
         NLI = NLI+1
         XLOC=DEPTH-X

         PRINT*,'...NLI,NELS, XLOC, ELEV, X, Z ', NLI,NELS, XLOC, ELEV, X, Z

C	
C        Interpolate input profiles to obtain line plume boundary conditions
         
         CALL LININT(LOCAT,UA,LAYERS,XLOC,UAMB)
         CALL LININT(LOCAT,VA,LAYERS,XLOC,VAMB)

         IF (NITERPLUME.EQ.1)THEN
            CALL LININT(LOCAT,CO2M,LAYERS,XLOC,COMG)
            DOARO=COMG
            CALL LININT(LOCAT,TE,LAYERS,XLOC,TARO)
            TAMB=TARO
            VO=0
            SALARO=SALAMB
         ELSEIF (NITERPLUME.GT.1)THEN
            IF (-XLOC.LE.DEPTHINTR.OR.-XLOC.GT.ELEVT)THEN
               CALL LININT(LOCAT,CO2M,LAYERS,XLOC,COMG)
               DOARO=COMG
               CALL LININT(LOCAT,TE,LAYERS,XLOC,TARO)
               TAMB=TARO
               !VO=0
               SALARO=SALAMB
            ELSEIF (-XLOC.GT.DEPTHINTR.AND.-XLOC.LE.ELEVT)THEN
            !   !PRINT*,'NLI',NLI, NLIO
               COMG=OUTO2(NLIO)
               DOARO=COMG
               CALL LININT(LOCAT,TE,LAYERS,XLOC,TAMB)

               TARO=OUTTEMP(NLIO)
               !TARO=TAMB
               SALARO=OUTSAL(NLIO)
               !SALARO=SALAMB

               OUTERQW=OUTQW(NLIO)
               OUTERMOMENT=OUTMOMENT(NLIO)
               VO=OUTERMOMENT/OUTERQW
               !VO=0
               !VO=-0.03

               NLIO=NLIO+1
               !PRINT*,'VO', VO 
            ENDIF
          ELSE
            PRINT*, "---------------------------------"
            PRINT*, "-------ERROR NITERPLUME----------"
            PRINT*, "---------------------------------"
         ENDIF
         
            PRINT*,'Initial QW, MOMENT, FTEMP, FDO, FGO,V0',QW, MOMENT, FTEMP, FDO, FGO, VO

C                  
C        Use subroutines for Runge Kutta method solution
         NEQN=8
         Y(1)=QW
         Y(2)=MOMENT
         Y(3)=FTEMP
         Y(4)=FSAL
         Y(5)=FDO
         Y(6)=FDN
         Y(7)=FGO
         Y(8)=FGN
 
         CALL DERIVS_6(EI,EO,DENSEA,DENSEW,DENSEP,G,BI,LI,LAMBDA,TARO,
     +            VG,SALARO,GAMMA,DENSE20,DOARO,PI,RB,N,VI,VO,VB,KOLO,
     +            HO2,PO,GAMMA1,TPLUME,SALPLU,COMGP,DNAMB,KOLN,HN2,PN,
     +            CNMGP,Z,Y,DYDX,XLOC,TAMB)

       
      !PRINT*, "B1 EI,EO,DENSEA", EI,EO,DENSEA
      !PRINT*, "B2 DENSEW,DENSEP", DENSEW,DENSEP
      !PRINT*, "B3", BI,LI,LAMBDA
      !PRINT*, "B4", TARO,VG,SALARO
      !PRINT*, "B5", GAMMA,DENSE20,DOARO
      !PRINT*, "B6", PI,RB,N
      !PRINT*, "B7 VI,VO,VB", VI,VO,VB
      !PRINT*, "B8", KOLO,HO2,PO
      !PRINT*, "B9", GAMMA1,TPLUME,SALPLU
      !PRINT*, "B10", COMGP,DNAMB,KOLN
      !PRINT*, "B11", HN2,PN,CNMGP
      !PRINT*, "B12", Y,DYDX,NEQN
      !PRINT*, "B13", Z,H,YOUT
      !PRINT*, "B14", XLOC,TAMB

         CALL RK4_6(EI,EO,DENSEA,DENSEW,DENSEP,G,BI,LI,LAMBDA,
     +         TARO,VG,SALARO,GAMMA,DENSE20,DOARO,PI,RB,N,VI,
     +         VO,VB,KOLO,HO2,PO,GAMMA1,TPLUME,SALPLU,COMGP,
     +         DNAMB,KOLN,HN2,PN,CNMGP,Y,DYDX,NEQN,Z,H,YOUT,
     +         XLOC,TAMB)

         QW=YOUT(1)
         MOMENT=YOUT(2)
         FTEMP=YOUT(3)
         FSAL=YOUT(4)
         FDO=YOUT(5)
         FDN=YOUT(6)
         FGO=YOUT(7)
         FGN=YOUT(8)

        PRINT*,'Solution QW, MOMENT, FTEMP, FDO, FGO',QW, MOMENT, FTEMP, FDO, FGO

         IF(MOMENT.LT.0.0)THEN
            TPLUME=FTEMP/QW
            SALPLU=FSAL/QW  !JCT2020_Sal
            CO2=FDO/QW
	        CN2=FDN/QW

!           Save inner plume information
            INPLUME(NLI,1)=XLOC 
            INPLUME(NLI,2)=BI          
            INPLUME(NLI,3)=QW
            INPLUME(NLI,4)=MOMENT
            INPLUME(NLI,5)=TPLUME
            INPLUME(NLI,6)=COMGP
            INPLUME(NLI,7)=EI 
            INPLUME(NLI,8)=EO
            INPLUME(NLI,9)=LI
            PRINT*,'MOMENT=0: XLOC, BI, QW, MOMENT, TPLUME, COMGP, EI, EO, LI',XLOC, BI, QW, MOMENT, TPLUME, COMGP, EI, EO, LI
	      GOTO 20
         ENDIF


         VI=MOMENT/QW
         AREA=QW/VI
         LI=LNOT          ! ACC 2026: longitud fija = longitud del difusor
         BI=AREA/(2.0*LI) ! ACC 2026: solo crece el semiancho
         EI= 2.*LI*ALPHAI*(VI+C1*VO)     ! ACC 2026: solo lados largos
         EO=-2.*LI*ALPHAO*VO             ! ACC 2026: solo lados largos
         TPLUME=FTEMP/QW
         SALPLU=FSAL/QW  !JCT2020_Sal

C        Previous equation corrected to consistently express salinity in uS/cm
!        Dissolved oxygen and nitrogen concentration
         CO2=FDO/QW
         CN2=FDN/QW 
         COMGP=CO2*32.
         CNMGP=CN2*28.
!        Save inner plume information
         INPLUME(NLI,1)=XLOC 
         INPLUME(NLI,2)=BI  
         INPLUME(NLI,3)=QW
         INPLUME(NLI,4)=MOMENT
         INPLUME(NLI,5)=TPLUME
         INPLUME(NLI,6)=COMGP
         INPLUME(NLI,7)=EI 
         INPLUME(NLI,8)=EO
         INPLUME(NLI,9)=LI

         PRINT*,'Solution: VI, LI, BI,EI,EO',VI, LI, BI,EI,EO


!        Add incremental entrainment to total cell entrainment/withdrawal   
         QWDI(LAYDIFF-JJ)=QWDI(LAYDIFF-JJ)+(EI-EO)*DZ
!         PWDI(LAYDIFF-JJ)=PWDI(LAYDIFF-JJ)+(2.*(LI+2.*BI)) ! JCT_RECT (erroneo, comentado)
         PWDI(LAYDIFF-JJ)=PWDI(LAYDIFF-JJ)+(2.*LI)          ! ACC 2026
         AWDI(LAYDIFF-JJ)=AWDI(LAYDIFF-JJ)+(LI*BI*2) 
         BWDI(LAYDIFF-JJ)=BWDI(LAYDIFF-JJ)+BI 
         LWDI(LAYDIFF-JJ)=LWDI(LAYDIFF-JJ)+LI  ! ACC 2026
         HWITH=HWITH+DZ   
         IF(HWITH.GT.HCELL)THEN
            PWDI(LAYDIFF-JJ) = PWDI(LAYDIFF-JJ)/NELS
            AWDI(LAYDIFF-JJ) = AWDI(LAYDIFF-JJ)/NELS
            BWDI(LAYDIFF-JJ) = BWDI(LAYDIFF-JJ)/NELS
            LWDI(LAYDIFF-JJ) = LWDI(LAYDIFF-JJ)/NELS
            JJ=JJ+1
            HWITH=0.0
            QWDI(LAYDIFF-JJ)=0.0
            PWDI(LAYDIFF-JJ)=0.0
            BWDI(LAYDIFF-JJ)=0.0
            LWDI(LAYDIFF-JJ)=0.0
            AWDI(LAYDIFF-JJ)=0.0
            NELS = 0
         ENDIF
!
C        Revised gaseous flux equations.
         YO2=FGO/((2.*LAMBDA*BI*LI)*(VI+VB))  ! ACC 2026
         YN2=FGN/((2.*LAMBDA*BI*LI)*(VI+VB))  ! ACC 2026
C  
!        
         PZ=PATM+(DENSEA*G*Z)
         QGAS=(FGO+FGN)*RGAS*(TPLUME+273.15)/PZ
         VBUB=QGAS/N
C
!        GAS VOLUME PER TOTAL VOLUME OF THE BUBBLE-WATER MIXTURE IN THE INNER CORE OF THE PLUME
         VG=VBUB*N/((VI+VB)*(2.*LAMBDA*BI*LI))  ! ACC 2026: area burbujas = LAMBDA*2*BI*LI
C        Previous equation revised to account for correct plume cross-sectional area occupied by bubbles.
!        Bubbles radius     
         RB=(3.*QGAS/(4.*PI*N))**(1./3.)
         IF(RB.LT.0.0)THEN
            RB=1.0E-8
         ENDIF
         FRACO=FGO/(FGO+FGN)
         FRACN=1.0-FRACO
C	
!        Partial Pressure
         PO=PZ*FRACO
         PN=PZ*FRACN
!        Density (ambient water, water in plume and of the plume)
      DENSEA=(0.059385*TAMB**3-8.56272*TAMB**2+65.4891*TAMB)*0.001
     ++999.84298 !+(GAMMA)*SALARO !JCT_2022
      DENSEW=(0.059385*TPLUME**3-8.56272*TPLUME**2+65.4891*TPLUME)*0.001
     ++999.84298 !+(GAMMA)*SALPLU  !JCT_2022
         DENSEP=(1.0-VG)*DENSEW


C
C        BUBBLE PROPERTIES
!
!        Bubble rise velocity
         IF(RB.LE.(7.5E-4))THEN
            VB=1189.0*RB**1.1945
         ELSEIF(RB.GT.(7.5E-4).AND.RB.LT.(4.8E-3))THEN
            VB=0.22
         ELSE
            VB=2.995*RB**0.489
         ENDIF


C
!        Mass transfer coeffitient
         KOLO=0.6*RB
         IF(KOLO.GT.(4.0E-4))THEN
            KOLO=4.0E-4
         ENDIF
         KOLN=KOLO
C
!        Solubility constant (mol/m3/Pa)
         HO2=(2.125-0.05023*TPLUME+5.7714E-4*TPLUME**2)/100000.
         HN2=(1.042-0.02457*TPLUME+3.1714E-4*TPLUME**2)/100000.
C
!        Froude number
         FR=VI/(2.*LAMBDA*BI*G*(DENSEA-DENSEP)/DENSEP)**0.5
!
         DCO2=HO2*PO-CO2

         PRINT*,'Solution: QGAS, VG, RB,VB, DENSEA,DENSEW,DENSEP,FR',QGAS, VG, RB,VB,DENSEA,DENSEW,DENSEP,FR

C      
      END DO
C
!     --------------------------------------------------------------------
C     CALCULATION OF AVERAGE NET OXYGEN MASS TRANSFER FOR DAY
!     --------------------------------------------------------------------

      PRINT*,'.........AVERAGE NET OXYGEN MASS TRANSFER FOR DAY..........'
   20 GROSSMT=(FGONOT-FGO)*32./1000.*86400.
      OTEFF=(FGONOT-FGO)/FGONOT*100.
      DELTAC=COMGP-COMGNOT
C
      ELEVT=ELEV
      QWT=QW
      TPLUMET=TPLUME
      COMGPT=COMGP
      LAYTOP=LAYDIFF-JJ
      BTOP=BI
      SALPLUMET=SALPLU

       PRINT*,'GROSSMT,OTEFF,ELEVT,QWT,TPLUMET,',GROSSMT,OTEFF,ELEVT,QWT,TPLUMET
       PRINT*,'COMGPT,LAYTOP,BTOP,SALPLUMET', COMGPT,LAYTOP,BTOP,SALPLUMET
       PRINT*,'Inner plume END'
      RETURN
      END
C
C------------------------------------------------------------------------------
C                          
         
C
C------------------------------------------------------------------------------
C                          
              
C     *******************************************************************************
	SUBROUTINE OUTER_PLUME_RECT(YEAR,JULDAY,WSEL,DIFFEL,LAYERS,BNOT,
     +LNOT,BITOP,SALAMB,PATM,QINTOP,FRCONOT,LAYDIFF,HCELL,LOCAT,TE,
     +CO2M,UA,VA,ELEVT,QWDET,QWDO,BWDO,LAYTOP,LAYINTR,TPLUMED,
     +SALPLUMED,COMGPD,DEPTHINTR,NLI,NLO,NITERPLUME,INPLUME, 
     +OUTPLUME,ALPHAI,ALPHAO,ALPHAA,LAMBDA,FRNI,FRNO,GAMMA1,LWDO)
C     *******************************************************************************

      IMPLICIT NONE

C	THIS SUBROUTINE IS WRITTEN TO PREDICT THE PERFORMANCE OF A CIRCULAR BUBBLE PLUME.  
C     THE MODEL IS BASED ON THE WUEST ET AL. (1992) CIRCULAR BUBBLE PLUME MODEL.
C     By:
C	
C	VERSION 3 (to couple with Francisco Rueda's reservoir model and for Lake Hallwil real diffuser)
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
C	  ALPHAI=ENTRAINMENT COEFFICIENT INNER PLUME  (-)
C	  ALPHAO=ENTRAINMENT COEFFICIENT OUTER PLUME  (-)
C	  ALPHAA=ENTRAINMENT COEFFICIENT FROM AMBIENT (-)
!         BI=1/2 INNER PLUME WIDTH (m)
!         BO=1/2 OUTER PLUME WIDTH (m)
!     BIFF=DIFFUSER RADIUS (m) 
C     BAVG=AVERAGE 1/2 DIFFUSER WIDTH (m)
!         BITOP=
!         BDET          	
!         C1=COEFFICIENT OF DIFERENT COUNTERFLOW ENTRAINMENTS
C	CO2=DISSOLVED OXYGEN (DO) CONCENTRATION (mol/m3)
C     CO2M=DO CONCENTRATION PROFILE FOR INPUT BOUNDARY CONDITION (g/m3)
C	CN2=DISSOLVED NITROGEN CONCENTRATION (mol/m3)
C	COMG=DISSOLVED OXYGEN CONCENTRATION (g/m3)
C     COMGPT=DO CONCENTRATION OF PLUME DETRAINMENT AT TOP OF PLUME (g/m3) 
!         COMGPD
!         COMGPTI
C	CNMG=DISSOLVED NITROGEN CONCENTRATION (g/m3)
C     DENSE20=DENSITY OF WATER AT 20 C
C	DENSEA=AMBIENT WATER DENSITY (kg/m3)
C	DENSEP=DENSITY OF THE PLUME (kg/m3)
C	DENSEW=WATER DENSITY IN PLUME (kg/m3)
!         DEPTHINTR=INTRUSION DEPTH OF THE DOUBLE PLUME
C	DIAMM=BUBBLE DIAMETER (mm)
C     DIFFEL=DIFFUSER ELEVATION (m)
C     DMPR=DEPTH OF MAXIMUM PLUME RISE (m)
C     DNAMB=AMBIENT DISSOLVED NITROGEN CONCENTRATION (g/m3)
C     DOAMB=AMBIENT DISSOLVED OXYGEN CONCENTRATION (g/m3)
C	E=ENTRAINMENT FACTOR (m3/s)
!	  EI=ENTRAINMENT  (m3/s)
!	  EO=ENTRAINMENT  (m3/s)
!         EA=ENTRAINMENT  (m3/s)
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
!	FR=FROUDE NUMBER (-)
C     FRCNOT=FRACTION OF NITROGEN IN ATMOSPHERE (-)
C     FRNO=INITIAL FROUDE NUMBER (-)
C     GAMMA=SALINITY CONVERSION FACTOR [kg/m3/(uS/cm)]
C     GROSSMT=GROSS MASS TRANSFER OF OXYGEN FROM PLUME (kg/d)
!	HO2=SOLUBILITY CONSTANT FOR OXYGEN (mol/m3/Pa)
!	HN2=SOLUBILITY CONSTANT FOR NITROGEN (mol/m3/Pa)
C     HCELL=HEIGHT OF CELL IN GRID (m)
C	HO=SOLUBILITY CONSTANT FOR OXYGEN (mol/m3/Pa)
C	HN=SOLUBILITY CONSTANT FOR NITROGEN (mol/m3/Pa)
C     HWITH=HEIGHT OF WITHDRAWAL/ENTRAINMENT ZONE (m)
!         INPLUME(:,6): INNER PLUME CHARACTERISTIC [HEIGTH; FLOWRATE OF WATER; MOMENTUM; TEMPERATURE; SALINITY; O2] JCT
C     JULDAY=JULIAN DAY IN GIVEN YEAR
C	KOLO=MASS TRANSFER COEFFICIENT FOR OXYGEN (m/s)
C	KOLN=MASS TRANSFER COEFFICIENT FOR NITROGEN (m/s)
C     LAKE=LAKE AND DIFFUSER TYPE (SHR AND LINEAR=1 OR AMISK AND RECTANGULAR=2) FOR SELECTION OF LAMBDA
C     LAYERS=NUMBER OF LAYERS/DATA POINTS IN BOUNDARY CONDITION PROFILES (-)
C     LAYDIFF=GRID LAYER CORRESPONDING TO DIFFUSER DEPTH (-)
!         LAYINTR=GRID LAYER CORRESPONDING TO INTRUSION DEPTH OF THE DOUBLE PLUME (-)
C     LAYTOP=GRID LAYER CORRESPONDING TO TOP OF PLUME (-)
C     LDIFF=LENGTH OF DIFFUSER (m)
C     LNOT=INITIAL DIFFUSER LENGTH (m)
C	LAMBDA=FRACTION OF PLUME OCCUPIED BY BUBBLES (-)
C     LAMBNOT=LAMBDA x INITIAL PLUME RADIUS;EQUAL TO DIFFUSER RADIUS (m) 
C     LOCAT=DEPTHS FOR INPUT BOUNDARY CONDITION PROFILES (m) 
C	MOMENT=MOMENTUM (m4/s)
C	N=NUMBER OF BUBBLES PER SECOND (1/s)
C     OTEFF=OXYGEN TRANSFER EFFICEINCY (%)
!         OUTPLUME(:,6)= OUTER PLUME CHARACTERISTIC [HEIGTH; FLOWRATE OF WATER; MOMENTUM; TEMPERATURE; SALINITY; O2] JCT
!         OUTERMOMENT= MOMENTUM IN THE OUTER PLUME 
!         OUTERQW=FLOWRATE OF WATER IN THE OUTER PLUME
!         OUTXLOC
!         OUTQW
!         OUTMOMENT
!         OUTTEMP
!         OUTO2  
C     PATM=ATMOSPHERIC PRESSURE AT AVERAGE WSEL (Pa)
C	PSTD=STANDARD PRESSURE (Pa)
!     QGAS=TOTAL GAS FLOW TO DIFFUSER
!       QINTOP
C!!!!	QSCFM=STANDARD GAS FLOW RATE (scfm), TOTAL STANDAR GAS FLOW RATE TO DIFFUSER 
C	QSCMS=STANDARD GAS FLOW RATE (scms)
C	QW=FLOWRATE OF WATER (m3/s)
!         QWDO
C     QWDET=TOTAL DETRAINMENT FLOW RATE AT THE INTRUSION OF THE PLUME (m3/s)  
C	RB=BUBBLE RADIUS (m)
C	RGAS=IDEAL GAS CONSTANT (J/mol/K)
C     SAL=SALINITY (uS/cm)
C     SALAMB=SALINITY OF AMBIENT WATER (uS/cm)
C	SALPLU=SALINITY OF THE PLUME (uS/cm)
!         SALPLUMED
!         SALPLUMETI
!     SHEARE=SHEAR ENTRAINMENT
C	  TARO=AROUND WATER TEMPERATURE (C)
C     TAVG=AVERAGE AMBIENT WATER TEMPERATURE (C)
C     TDS=TOTAL DISSOLVED SOLIDS (g/m3) [0.64 conversion factor from Chapra book]
C     TE=TEMPERATURE PROFILE FOR INPUT BOUNDARY CONDITION (C) 
C	TPLUME=PLUME TEMPERATURE (C)
C     TPLUMET=DETRAINMENT PLUME TEMPERATURE AT THE TOP OF THE PLUME (C) 
!         TPLUMED
!         TSTD=STANDARD TEMPERATURE (K)
!     UA = U AMBIENT VELOCITY PROFILE FOR INPUT BOUNDARY CONDITION
C     UAMB = U Ambient Velocity (m) used to calculate VORTEX entrainment (FJR
!	  VI=WATER VELOCITY (m/s)
!         VO=WATER VELOCITY (m/s)
!     VA = V AMBIENT VELOCITY PROFILE FOR INPUT BOUNDARY CONDITION
C     VAMB = V Ambient Velocity (m) used to calculate VORTEX entrainment (FJR
C     VAVG=AVERAGE WATER VELOCITY (m/s)
C	VB=BUBBLE RISE VELOCITY (m/s)
C	VBUB=BUBBLE VOLUME (m3)
!     VG=GAS VOLUME PER TOTAL VOLUME OF THE BUBBLE-WATER MIXTURE IN THE INNER CORE OF THE PLUME
C	VGUESS=GUESSED INITIAL WATER VELOCITY (m/s)
!         VI=WATER VELOCITY IN THE INNER PLUME (m/s)
!         VO=WATER VELOCITY IN THE OUTER PLUME (m/s)
!     VORTEXE=VORTEX ENTREAINMENT
C     WSEL=WATER SURFACE ELEVATION (m)
C	YO2=GASEOUS OXYGEN CONCENTRATION (mol/m3)
C	YN2=GASEOUS NITROGEN CONCENTRATION (mol/m3)
C	Z=DEPTH TO DIFFUSER (m)
C     PWD  = Perimeter (m) (FJR)
C     VUP  = Upward velocity (m/2)
C     AWD  = Area (m2) of plume
C     BWD  = Width 
C
      REAL*8 AREA,BI,BO,LI,LO,CO2,COMG,CN2,CNMG,DS,DZ,DENSEA,DENSEP,DENSEW,
     +DIAMM,FDO,FDN,FRACO,FRACN,FSAL,FTEMP,FGO,FGN,G,GAMMA,HO2,HN2,
     +KOLN,KOLO,LAMBDA,MOMENT,N,PI,PO,PN,PSTD,PZ,QSCMS,QW,QGAS,
     +RB,RGAS,SALAMB,SALPLU,TPLUME,TSTD,VB,VBUB,VG,VGUESS,
     +YO2,YN2,Z,AA,BB,CC,BNOT,LNOT,ELEV,DT,TE(1000),XLOC,DEPTH,
     +CO2M(1000),LOCAT(1000),DOAMB,COMGP,CNMGP,WSEL,PATM,SAL(1000),
     +DIFFEL,GROSSMT,FGONOT,DNAMB,DMPR,SUMSAL,FRCNATM,
     +LAMBNOT,FRCONOT,H,DYDX(5),Y(5),YOUT(5),BTOP,
     +DENSE20,OTEFF,FRNI,FRNO,VDIFF,FR,DCO2,QGFRAC,LDIFF,
     +JDAY,EL(70),DELTAC,COMGNOT,ELEVT,QWDET,QWDO(500),
     +PWDO(500),BWDO(500),AWDO(500),X,SHEARE,VORTEXE,
     +HWITH,HCELL,JULDAY,BDET,UA(1000),VA(1000),UAMB,VAMB,
     +INPLUME(35000,9),OUTPLUME(35000,9),DEPTHINTR,ALPHAI,ALPHAO,ALPHAA,
     +C1,EI,EO,EA,INNERQW,INNERMOMENT,INXLOC(35000),INQW(35000),
     +INMOMENT(35000),INTEMP(35000),INO2(35000),INSAL(35000),
     +INBI(35000),INLI(35000),BITOP,LITOP,VI,VO,TPLUMED,QINTOP,TARO,
     +DEPTHDMPR,SALPLUMED,COMGPD,DENSEINNER,GAMMA1,TIN,COMGIN,SALIN,EP,
     +LWDO(500)	 
      INTEGER II,IJ,IK,JJ,LL,NEQN,NN,MI,JK,JL,LAYTOP,
     +LAYERS,KM,KN,KO,KP,ROWS,KQ,KR,KU,KV,KW,KX,KY,KZ,KS,LAYDIFF,YEAR 
      INTEGER NELS,ierror,M,CONT,LAYINTR,NLI,NLO,NITERPLUME
C 
C     --------------------------------------------------------------------------
C     CONSTANTS
C     --------------------------------------------------------------------------
!

      PRINT*,'-----OUTER_PLUME_RECT------'

      PRINT*,'NITERPLUMEinOUTERPLUME', NITERPLUME

      PRINT*, 'OUTER_PLUME',ALPHAI,ALPHAO,ALPHAA,LAMBDA,FRNI,FRNO,GAMMA1


    !  PRINT*,'Yearplume', YEAR
    !  PRINT*,'Julday', JULDAY
    !  PRINT*,'WSEL', WSEL
    !  PRINT*,'DIFFEL', DIFFEL
    !  PRINT*,'LAYERS', LAYERS
    !  PRINT*,'BITOP', BITOP
    !  PRINT*,'SALAMB', SALAMB
    !  PRINT*,'PATM', PATM
    !  PRINT*,'QINTOP', QINTOP
    !  PRINT*,'FRCONOT',FRCONOT
    !  PRINT*,'LAYDIFF',LAYDIFF
    !  PRINT*,'HCELL',HCELL
      !PRINT*,'LOCAT',LOCAT
      !PRINT*,'TE',TE
    !  !PRINT*,'CO2M',CO2M
    !  !PRINT*,'UA',UA
    !  !PRINT*,'VA',VA
    !  PRINT*,'ELEVT',ELEVT
    !  PRINT*,'QWDET',QWDET
    !  !PRINT*,'QWDO',QWDO
    !  !PRINT*,'BWDO',BWDO
    !  PRINT*,'LAYTOP',LAYTOP
    !  PRINT*,'LAYINTR',LAYINTR
    !  PRINT*,'TPLUMED',TPLUMED
    !  PRINT*,'SALPLUMED',SALPLUMED
    !  PRINT*,'COMGPD',COMGPD
    !  PRINT*,'DEPTHINTR',DEPTHINTR
    !  PRINT*,'NLI',NLI
    !  PRINT*,'NLO',NLO
    !  PRINT*,'NITERPLUME',NITERPLUME
    !  !PRINT*,'INPLUME',INPLUME
    !  !PRINT*,'OUTPLUME',OUTPLUME


      !ALPHAI=0.055 ! Crounse et al 2007
      !ALPHAO=0.11  ! Crounse et al 2007
      !ALPHAA=0.11  ! Crounse et al 2007
      G=9.80665
      GAMMA=6.9E-4
      !GAMMA1=1.1 ! JCT Socolofsky
!      LAMBDA=0.8 !McGinnis
      !LAMBDA=1 !Socolofsky et al. 2008
      !FRNO=0.1  !Socolofsky et al. 2008
      !FRNO=-FRNO  !Socolofsky et al. 2008 Para ser coherente en los signos
      PI=ACOS(-1.0)
      PSTD=101325.
      RGAS=8.314
      TSTD=293.15
      DENSE20=998.2
      FRCNATM=0.79
      QGFRAC=1.0
!      C1=-1
      C1=0
!     
C     --------------------------------------------------------------------------
C     PARAMETERS AND INITIALIZE THE VARIABLES 
C     --------------------------------------------------------------------------
C     
      OUTPLUME(:,:)=0.0
!     Information from the INNER plume

      INXLOC=INPLUME(:,1)
      INBI=INPLUME(:,2)
      INLI=INPLUME(:,9) ! JCT_RECT
      INQW=INPLUME(:,3)
      INMOMENT=INPLUME(:,4)
      INTEMP=INPLUME(:,5)
      INO2=INPLUME(:,6)
      INSAL=INPLUME(:,7)


!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!     Dejamos como zona de transicion los 0.5 metros superiores de la pluma. JCT
!     BITOP y QINTOP dejan de ser input, los sacamos de INPLUME --> ELIMINAR  JCT
      NLI = NLI -10
      BITOP  = INPLUME(NLI,2)
      LITOP  = INPLUME(NLI,9) ! JCT_RECT	  
      QINTOP = -INPLUME(NLI,3) ! JCT_2022
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
!     Geometric caracteristic 
      DEPTH=WSEL-DIFFEL
      DEPTHDMPR=WSEL-ELEVT ! Restamos 0.1 metros para eliminar la zona de transicion
      Z=DEPTHDMPR   
      ELEV=ELEVT
!        
!!C     Assume that gas bubbles are composed of oxygen and nitrogen only.
!!      FRACO=FRCONOT 
!!      FRACN=1.0-FRACO
!
C     Interpolate input profiles to obtain the plume initial conditions
      X=0.0 
      XLOC=Z-X 

      CALL LININT(LOCAT,TE,LAYERS,XLOC,TARO)
      TIN=INTEMP(NLI)
      CALL LININT(LOCAT,CO2M,LAYERS,XLOC,COMG)
      COMGIN=INO2(NLI)
      CALL LININT(LOCAT,UA,LAYERS,XLOC,UAMB)
      CALL LININT(LOCAT,VA,LAYERS,XLOC,VAMB)
      SALIN=INSAL(NLI)

      DOAMB=COMG
      CO2=COMG/32.
      COMGNOT=COMGP
!     At first temperature, DO and salinity in the outer plume is the same that in the top of the inner plume
      TPLUME=TIN
      COMGP=COMGIN
      SALPLU=SALIN  
      SALPLU=224  ! Eliminar JCT SAlinidad
	  SALIN = 224

!     Guess the initial velocity 
      VGUESS=-0.07
      VO=VGUESS
C      
C     AMBIENT AND AVERAGE WATER DENSITIES

      DENSEA=(0.059385*TARO**3-8.56272*TARO**2+65.4891*TARO)*0.001
     ++999.84298 !+(GAMMA)*SALAMB !JCT_2022
!!      DENSEW=DENSEA
      DENSEP=(0.059385*TPLUME**3-8.56272*TPLUME**2+65.4891*TPLUME)*0.001
     ++999.84298 !+(GAMMA)*SALPLU !JCT_2022
      DENSEINNER=(0.059385*TIN**3-8.56272*TIN**2+65.4891*TIN)*0.001
     ++999.84298 !+(GAMMA)*SALIN !JCT_2022

C	
        PRINT*,"DENSEP,DENSEA",DENSEP,DENSEA,TPLUME,TARO

C     Assume initial ambient dissolved nitrogen conc. equals saturated conc. at surface.
!!      CN2=(PATM*FRCNATM)*HN2      
!!      CNMG=CN2*28.0
!!      CNMGP=CNMG
!!      DNAMB=CNMG
!
!
C     --------------------------------------------------------------------------
C     PLUME PROPERTIES
C     --------------------------------------------------------------------------

!
C     CALCULATION OF INITIAL WATER VELOCITY USING FROUDE NUMBER

      VDIFF=1
      ! PRINT*,"QINTOP",QINTOP
      ! PRINT*,"PI",PI
      ! PRINT*,"VO",VO
      ! PRINT*,"BITOP",BITOP
      ! PRINT*,"G",G
      ! PRINT*,"DENSEA",DENSEA
      ! PRINT*,"DENSEP",DENSEP
	  
		  
		  
      DO WHILE (VDIFF.GT.1.0E-6)
         !BO=SQRT((QINTOP**2)/(PI*QINTOP*ABS(VO))+BITOP**2)
         !BO=SQRT((QINTOP)/(PI*QINTOP*ABS(VO))+BITOP**2)
		 AREA = QINTOP/VO+(LITOP*2.*BITOP) ! JCT_RECT
		 !AREA = QINTOP/VO ! JCT_RECT_2022
         !SOLVE FOR DIMENSIONS USING L^2+(2Bo-Lo)L-AREA=0 USING QUADRATIC EQN.
          AA=1.0
          ! BB=2.*BNOT-LNOT
            ! BB=2.*BNOT-LNOT
            BB=LI-2.*BI ! JCT_RECT_2022
          CC=-1.0*AREA
          LO=(-1.0*BB+(BB**2-4.0*AA*CC)**(0.5))/(2.0*AA)
          IF(LO.LT.0.0)THEN
             LO=(-1.0*BB-(BB**2-4.0*AA*CC)**(0.5))/(2.0*AA)
          ENDIF		 
		  BO=AREA/(2.0*LO)
		  ! VO=-FRNO*(ABS((BO-BITOP)*G*(DENSEA-DENSEP)/DENSEP))**0.5
		  ! VO=-FRNO*(QINTOP/VO/(LITOP + LO)*G*(ABS((DENSEA-DENSEP)/DENSEP)))**0.5 ! JCT_rect_2022
!         ACC 2026: ERROR - escala de longitud Froude usaba (BO-BITOP) en lugar de (BO-BI),
!         semiancho del anillo exterior. Ademas DENSEP es densidad de la pluma interna, no
!         de la exterior. La pluma externa es mas densa que el ambiente (DENSEP>DENSEA).
!         Referencia: Socolofsky et al. (2008) Eq.12. Linea erronea comentada:
!          VO=-FRNO*(ABS((BO-BITOP)*G*(DENSEA-DENSEP)/DENSEP))**0.5
          VO=-FRNO*(ABS((BO-BITOP)*G*(DENSEP-DENSEA)/DENSEA))**0.5  ! ACC 2026
		  VDIFF=ABS(VO-VGUESS)
		  VGUESS=VO
      END DO
      PRINT*,"QINTOP",QINTOP
      PRINT*,"VOinicial",VO
!     Water velocity in the top of the inner plume
      !VI=QINTOP/(PI*BITOP**2)
      VI=-QINTOP/(LITOP*BITOP*2) ! JCT_RECT JCT_2022
	  
	  
C     ------------------------------------------------------------------
C     VARIABLE TRANSFORMATION      
C     ------------------------------------------------------------------
!     TOTAL ENTRAIMENT: Shear and Vortex
!       E=2.*PI*BO*ALPHA*VO
!      PRINT*,"Shear" 
!      SHEARE=2.*PI*BO*ALPHA*VO
!     Vortex Entraiment: Hypothesis Projected Area Entraiment
!      VORTEXE=(2*BO*(SQRT(ABS(UAMB)*ABS(UAMB)+ABS(VAMB)*ABS(VAMB)))) 
!     a) Additive hypothesis
!         PRINT*,"Additive"  
!       E=SHEARE+VORTEXE
!     b) Maximun hypothesis
!         PRINT*,"Maximun"  
!      IF (SHEARE.GT.VORTEXE) THEN
! 	E=SHEARE
!      ELSE
!	E=VORTEXE
!      ENDIF
      BI=BITOP
      LI=LITOP ! JCT_RECT
      !EI=+2.*PI*BI*ALPHAI*(VI+C1*VO)
      !EO= 2.*PI*BI*ALPHAO*VO
      !EA= 2.*PI*BO*ALPHAA*VO
	  EI=+(2.*(LI+2.*BI))*ALPHAI*(VI+C1*VO) ! JCT_RECT
      EO=-(2.*(LI+2.*BI))*ALPHAO*VO ! JCT_RECT JCT_2022
      EA=-(2.*(LO+2.*BO))*ALPHAA*VO ! JCT_RECT JCT_2022
      EP=QINTOP
	   
      PRINT*,"LI,LO,BI,BO,VI,VO", LI,LO,BI,BO,VI,VO
      PRINT*,"ALPHAI,ALPHAO,ALPHAA", ALPHAI,ALPHAO,ALPHAA

      PRINT*,"EA,EO,EI",EA,EO,EI, VO

      QW=QINTOP
      MOMENT=QW*VO
      PRINT*,"MOMENTinicial",MOMENT,VO,QW
      PRINT*,"TPLUME QW",TPLUME,QW
      FTEMP=QW*TPLUME 
      !FSAL=QW*(SALPLU*GAMMA/DENSE20)*DENSEP
      FSAL=QW*SALPLU !JCT2020_Sal
C     Previous equation corrected to account for salinity units conversion.        
      FDO=QW*CO2
!!      FDN=QW*CN2
!!      FGO=PSTD*QSCMS/(RGAS*TSTD)*FRACO
!!      FGONOT=FGO
!!      FGN=PSTD*QSCMS/(RGAS*TSTD)*FRACN
C     Revised gaseous flux equations.
!!      YO2=FGO/((PI*(LAMBDA*B)**2)*(V+VB))
!!      YN2=FGN/((PI*(LAMBDA*B)**2)*(V+VB))           
!!      PZ=PATM+(DENSEA*G*Z)
!!      PO=PZ*FRACO
!!      PN=PZ*FRACN
!     Initialize lateral withdrawal flowrate for first/lowest cell in column/segment
      JJ=0
      QWDO(LAYTOP)= QW 
      !PWDO(LAYTOP)= (2*PI*BO)
      PWDO(LAYTOP)= (2.*(LO+2.*BO)) ! JCT_RECT
      !AWDO(LAYTOP)= (PI*(BO**2-BI**2))    
      AWDO(LAYTOP)= (LO*BO*2) - (LI*BI*2)   ! JCT_RECT 
      BWDO(LAYTOP)= BO 
      LWDO(LAYTOP)= LO 

!
!     -------------------------------------------------------------------------- 
C	SOLUTION PROCEEDURE
!     --------------------------------------------------------------------------
!
      DZ=0.001
	  ! H=0.001
      H= -0.001 ! JCT_2022
      NELS = 0
      HWITH=0.0
      NLO=NLI+1
!
      M=0

!        PRINT*,"SOLUTION PROCEEDURE"
       DO WHILE (VO.LT.-1.E-6.AND.Z.LT.DEPTH)
       ! DO WHILE (M.LT.1)
         M=M+1
         Z=Z+DZ
         X=X+DZ
         ELEV=ELEV-DZ
         NELS = NELS + 1
         NLO=NLO-1
C	
C        Interpolate input profiles to obtain line plume boundary conditions
         XLOC=DEPTHDMPR+X

         CALL LININT(LOCAT,TE,LAYERS,XLOC,TARO)
         TIN=INTEMP(NLO)
         CALL LININT(LOCAT,CO2M,LAYERS,XLOC,COMG)
         COMGIN=INO2(NLO)
         SALIN=INSAL(NLO)
		 SALIN=224
         CALL LININT(LOCAT,UA,LAYERS,XLOC,UAMB)
         CALL LININT(LOCAT,VA,LAYERS,XLOC,VAMB)
C        Inner plume width
         BI=INBI(NLO)
         LI=INLI(NLO)
C        Inner velocity plume
         INNERQW=INQW(NLO)
         INNERMOMENT=INMOMENT(NLO)
         VI=INNERMOMENT/INNERQW

C        Use subroutines for Runge Kutta method solution
         NEQN=5
         Y(1)=QW
         Y(2)=MOMENT
         Y(3)=FTEMP
         Y(4)=FSAL
         Y(5)=FDO
		 
 
	    IF(M.EQ.1)THEN
	       PRINT*," DYDX(2)",DYDX(2), GAMMA1,G,LO,BO,LI,BI,DENSEP,DENSEA
		   PRINT*,DENSE20,EI,VO,EO,VI
        ENDIF


         CALL DERIVS_7(EA,EO,EI,VI,VO,DENSEA,DENSEP,G,PI,BO,BI,LO,LI,
     +             TARO,TIN,TPLUME,SALAMB,SALIN,SALPLU,GAMMA,DENSE20,
     +             DENSEINNER,GAMMA1,DOAMB,COMGIN,COMGP,DNAMB,
     +             CNMGP,Z,Y,DYDX)

         CALL RK4_7(EA,EO,EI,VI,VO,DENSEA,DENSEP,G,PI,BO,BI,LO,LI,TARO,
     +             TIN,TPLUME,SALAMB,SALIN,SALPLU,GAMMA,DENSE20,
     +             DENSEINNER,GAMMA1,DOAMB,COMGIN,COMGP,DNAMB,CNMGP,Y,
     +             DYDX,NEQN,X,H,YOUT)


C       
         QW=YOUT(1)
         MOMENT=YOUT(2)
         FTEMP=YOUT(3)
         FSAL=YOUT(4)
         FDO=YOUT(5)
		 
         IF(MOMENT.LT.0.0)THEN
            TPLUME=FTEMP/QW
	    !SALPLU=FSAL/(QW*DENSEP)/(GAMMA/DENSE20)
            SALPLU=224  ! Eliminar JCT SAlinidad
C           Previous equation corrected to consistently express salinity in uS/cm	   
	        CO2=FDO/QW
            COMGP=CO2*32.
            VO=MOMENT/QW
            !AREA=ABS(QW/VO)
			
            AREA=ABS(QW/VO) +(LI*BI*2) ! JCT_RECT_2022
			
            !BO=SQRT((AREA+PI*BI**2)/PI)
            !BO=SQRT((AREA+(LI*BI*2))/PI)

	  		 ! JCT_RECT  ....
            !SOLVE FOR DIMENSIONS USING L^2+(2Bo-Lo)L-AREA=0 USING QUADRATIC EQN.
            AA=1.0
            ! BB=2.*BNOT-LNOT
            BB=LI-2.*BI ! JCT_RECT_2022
            CC=-1.0*AREA
            LO=(-1.0*BB+(BB**2-4.0*AA*CC)**(0.5))/(2.0*AA)
            IF(LO.LT.0.0)THEN
               LO=(-1.0*BB-(BB**2-4.0*AA*CC)**(0.5))/(2.0*AA)
            ENDIF
            BO=AREA/(2.0*LO)
		  
!           Save outer plume information
            OUTPLUME(NLO,1)=XLOC 
            OUTPLUME(NLO,2)=BO          
            OUTPLUME(NLO,3)=QW
            OUTPLUME(NLO,4)=MOMENT
            OUTPLUME(NLO,5)=TPLUME
            OUTPLUME(NLO,6)=COMGP
            OUTPLUME(NLO,7)=SALPLU
            OUTPLUME(NLO,8)=VO
            OUTPLUME(NLO,9)=LO
			
            PRINT*,"me voy por momento"
	    GOTO 20
         ENDIF

         VO=MOMENT/QW
         !AREA=ABS(QW/VO)
         AREA=ABS(QW/VO) +  (LI*BI*2)! JCT_RECT  incluir area pluma interior para obtener el perimetro real
         ! BO=SQRT((AREA+(LI*BI*2))/PI)
		 
	  		 ! JCT_RECT  ....
         !SOLVE FOR DIMENSIONS USING L^2+(2Bo-Lo)L-AREA=0 USING QUADRATIC EQN.
          AA=1.0
          ! BB=2.*BNOT-LNOT
          BB=LI-2.*BI ! JCT_RECT_2022 
		  CC=-1.0*AREA
          LO=(-1.0*BB+(BB**2-4.0*AA*CC)**(0.5))/(2.0*AA)
          IF(LO.LT.0.0)THEN
             LO=(-1.0*BB-(BB**2-4.0*AA*CC)**(0.5))/(2.0*AA)
          ENDIF
          BO=AREA/(2.0*LO)
		  

        !PRINT*,"QW",QW
!        PRINT*,"B",B
!        TOTAL ENTRAIMENT: Shear and Vortex
!          E=2.*PI*B*ALPHAI*V
!         PRINT*,"Shear" 
!         SHEARE=2.*PI*B*ALPHAI*V
!        Vortex Entraiment: Hypothesis Projected Area Entraiment
!         VORTEXE=(2*B*(SQRT(ABS(UAMB)*ABS(UAMB)+ABS(VAMB)*ABS(VAMB)))) 
!	 a) Additive hypothesis 
!          PRINT*,"Additive" 
!      	     E=SHEARE+VORTEXE
!	 b) Maximun hypothesis
!          PRINT*,"Maximun" 
!	    IF (SHEARE.GT.VORTEXE) THEN
!	 	E=SHEARE
!            ELSE
!		E=VORTEXE
!            ENDIF
!         

         !EI=+2.*PI*BI*ALPHAI*(VI+C1*VO)
         !EO=-2.*PI*BI*ALPHAO*VO ! JCT_2020 VO es negativa
         !EA=-2.*PI*BO*ALPHAA*VO ! JCT_2020
		 
		 ! JCT_RECT ....
         EI=+(2.*(LI+2.*BI))*ALPHAI*(VI+C1*VO)
         EO=-(2.*(LI+2.*BI))*ALPHAO*VO ! JCT_2020 VO es negativa
         EA=-(2.*(LO+2.*BO))*ALPHAA*VO ! JCT_2020
         ! JCT_RECT ....
         
!
!        Temperatura and salinity in the plume    
         TPLUME=FTEMP/QW
         SALPLU=FSAL/(QW*DENSEP)/(GAMMA/DENSE20)
         SALPLU=224  ! Eliminar JCT SAlinidad
C        Previous equation corrected to consistently express salinity in uS/cm

!        Dissolved oxygen and nitrogen concentration
         CO2=FDO/QW
!         CN2=FDN/QW 
         COMGP=CO2*32.
!         CNMGP=CN2*28.

!        Save outer plume information
         OUTPLUME(NLO,1)=XLOC 
         OUTPLUME(NLO,2)=BO          
         OUTPLUME(NLO,3)=QW
         OUTPLUME(NLO,4)=MOMENT
         OUTPLUME(NLO,5)=TPLUME
         OUTPLUME(NLO,6)=COMGP
         OUTPLUME(NLO,7)=SALPLU
         OUTPLUME(NLO,8)=VO
         OUTPLUME(NLO,9)=LO
		 
!        Add incremental entrainment to total cell entrainment/withdrawal   
         !QWDO(LAYTOP+JJ)=QWDO(LAYTOP+JJ)+(EA+EO-EI)*DZ 
         QWDO(LAYTOP+JJ)=QWDO(LAYTOP+JJ)-(EA+EO-EI)*DZ !JCT_2022
         PWDO(LAYTOP+JJ)=PWDO(LAYTOP+JJ)+(2.*(LO+2.*BO))
         AWDO(LAYTOP+JJ)=AWDO(LAYTOP+JJ)+(LO*BO*2) - (LI*BI*2)
         BWDO(LAYTOP+JJ)=BWDO(LAYTOP+JJ)+BO
         LWDO(LAYTOP+JJ)=LWDO(LAYTOP+JJ)+LO
         HWITH=HWITH+DZ   
         ! PRINT*,"JJ",JJ 
         IF(HWITH.GT.HCELL)THEN
            PWDO(LAYTOP+JJ) = PWDO(LAYTOP+JJ)/NELS
            AWDO(LAYTOP+JJ) = AWDO(LAYTOP+JJ)/NELS
            BWDO(LAYTOP+JJ) = BWDO(LAYTOP+JJ)/NELS
            LWDO(LAYTOP+JJ) = LWDO(LAYTOP+JJ)/NELS
            !PRINT*,"JJ",JJ 
            JJ=JJ+1
            HWITH=0.0
            QWDO(LAYTOP+JJ)=0.0
            PWDO(LAYTOP+JJ)=0.0
            BWDO(LAYTOP+JJ)=0.0
            LWDO(LAYTOP+JJ)=0.0
            AWDO(LAYTOP+JJ)=0.0
            NELS = 0
         ENDIF

!
C        Revised gaseous flux equations.
!!         YO2=FGO/((PI*(LAMBDA*B)**2)*(V+VB))
!!         YN2=FGN/((PI*(LAMBDA*B)**2)*(V+VB))
C  
!        
!!         PZ=PATM+(DENSEA*G*Z)
!!         QGAS=(FGO+FGN)*RGAS*(TPLUME+273.15)/PZ
!!         VBUB=QGAS/N
C
!        GAS VOLUME PER TOTAL VOLUME OF THE BUBBLE-WATER MIXTURE IN THE INNER CORE OF THE PLUME
!!         VG=VBUB*N/((V+VB)*(PI*(LAMBDA*B)**2))
C        Previous equation revised to account for correct plume cross-sectional area occupied by bubbles. 
!        Bubbles radius     
!!         RB=(3.*QGAS/(4.*PI*N))**(1./3.)
!!         IF(RB.LT.0.0)THEN
!!           RB=1.0E-8
!!         ENDIF
!!         FRACO=FGO/(FGO+FGN)
!!         FRACN=1.0-FRACO
C	
!        Partial Pressure
!!         PO=PZ*FRACO
!!         PN=PZ*FRACN
!        Density (ambient water, water in plume and of the plume)
      DENSEA=(0.059385*TARO**3-8.56272*TARO**2+65.4891*TARO)*0.001
     ++999.84298 !+(GAMMA)*SALAMB !JCT_2022
      DENSEP=(0.059385*TPLUME**3-8.56272*TPLUME**2+65.4891*TPLUME)*0.001
     ++999.84298 !+(GAMMA)*SALPLU !JCT_2022
      DENSEINNER=(0.059385*TIN**3-8.56272*TIN**2+65.4891*TIN)*0.001
     ++999.84298 !+(GAMMA)*SALIN  !JCT_2022

        IF(M.EQ.1)THEN
        PRINT*,"DENSEP,DENSEA",DENSEP,DENSEA, TPLUME,TARO
        ENDIF
		
		IF(M.EQ.100)THEN
        PRINT*,"DENSEP,DENSEA",DENSEP,DENSEA, TPLUME,TARO
        ENDIF
		
		IF(M.EQ.1000)THEN
        PRINT*,"DENSEP,DENSEA",DENSEP,DENSEA, TPLUME,TARO
        ENDIF
		
		IF(M.EQ.5000)THEN
        PRINT*,"DENSEP,DENSEA",DENSEP,DENSEA, TPLUME,TARO
        ENDIF
C        Previous equation re-revised to account for correct salinity units (uS/cm) in density calculations.      
C
C        BUBBLE PROPERTIES
!
!        Bubble rise velocity
!!         IF(RB.LE.(7.5E-4))THEN
!!            VB=1189.0*RB**1.1945
!!         ELSEIF(RB.GT.(7.5E-4).AND.RB.LT.(4.8E-3))THEN
!!            VB=0.22
!!         ELSE
!!            VB=2.995*RB**0.489
!!         ENDIF
C
!        Mass transfer coeffitient
!!         KOLO=0.6*RB
!!         IF(KOLO.GT.(4.0E-4))THEN
!!            KOLO=4.0E-4
!!         ENDIF
!!         KOLN=KOLO
C
!        Solubility constant (mol/m3/Pa)
!!         HO2=(2.125-0.05023*TPLUME+5.7714E-4*TPLUME**2)/100000.
!!         HN2=(1.042-0.02457*TPLUME+3.1714E-4*TPLUME**2)/100000.
C
!        Froude number
!!         FR=V/(2.*LAMBDA*BO*G*(DENSEA-DENSEP)/DENSEP)**0.5
!
!!         DCO2=HO2*PO-CO2
C      
  !       PRINT*,"QW",QW
!         PRINT*,"VI",VI 
!         PRINT*,"VO",VO 
  !       PRINT*,"Z",Z
  !       PRINT*,"MOMENT",MOMENT

  !       PRINT*,"TIN",TIN
  !       PRINT*,"COMGIN",COMGIN
  !       PRINT*,"SALIN",SALIN
!         PRINT*,"BI",BI
!         PRINT*,"BO",BO
  !       PRINT*,"INNERQW",INNERQW
  !       PRINT*,"NLO",NLO

      END DO
C
!     --------------------------------------------------------------------
C     CALCULATION OF AVERAGE NET OXYGEN MASS TRANSFER FOR DAY
!     --------------------------------------------------------------------
      !GROSSMT=(FGONOT-FGO)*32./1000.*86400.
      !OTEFF=(FGONOT-FGO)/FGONOT*100.
      !DELTAC=COMGP-COMGNOT
C
!      ELEVT=ELEV
   20 DEPTHINTR=ELEV
      QWDET=QW
      TPLUMED=TPLUME
      COMGPD=COMGP
!!      LAYTOP=LAYDIFF-JJ
      LAYINTR=LAYTOP+JJ
!!      BTOP=B
!      PRINT*,"QWDET",QWDET
!      PRINT*,"OUTER PLUME SALIDAS"
!      PRINT*,"ELEVT",ELEVT
!      PRINT*,"QWT",QWT
!      PRINT*,"TPLUMED",TPLUMED
!      PRINT*,"COMGPD",COMGPD
!      PRINT*,"LAYTOP_outer",LAYTOP
 !     PRINT*,"LAYDIFF_outer",LAYDIFF
!      PRINT*,"LAYINTR_outer",LAYINTR
!      PRINT*,"JJ_outer",JJ
!      PRINT*,"BTOP",BTOP
!      PRINT*,"INPLUMEDEPTH",INPLUME(:,1)
      RETURN
      END
C     
C----------------------------------------------------------------------
C

             
C     *******************************************************************************
	SUBROUTINE OUTER_PLUME_RECT2(YEAR,JULDAY,WSEL,DIFFEL,LAYERS,BITOP,
     +SALAMB,PATM,QINTOP,FRCONOT,LAYDIFF,HCELL,LOCAT,TE,CO2M,UA,VA,
     +ELEVT,QWDET,QWDO,BWDO,LAYTOP,LAYINTR,TPLUMED,SALPLUMED,COMGPD,
     +DEPTHINTR,NLI,NLO,NITERPLUME,INPLUME,OUTPLUME,ALPHAI,ALPHAO, 
     +ALPHAA,LAMBDA,FRNI,FRNO,GAMMA1,LWDO)
C     *******************************************************************************

      IMPLICIT NONE

C	THIS SUBROUTINE IS WRITTEN TO PREDICT THE PERFORMANCE OF A CIRCULAR BUBBLE PLUME.  
C     THE MODEL IS BASED ON THE WUEST ET AL. (1992) CIRCULAR BUBBLE PLUME MODEL.
C     By:
C	
C	VERSION 3 (to couple with Francisco Rueda's reservoir model and for Lake Hallwil real diffuser)
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
C	  ALPHAI=ENTRAINMENT COEFFICIENT INNER PLUME  (-)
C	  ALPHAO=ENTRAINMENT COEFFICIENT OUTER PLUME  (-)
C	  ALPHAA=ENTRAINMENT COEFFICIENT FROM AMBIENT (-)
!         BI=1/2 INNER PLUME WIDTH (m)
!         BO=1/2 OUTER PLUME WIDTH (m)
!     BIFF=DIFFUSER RADIUS (m) 
C     BAVG=AVERAGE 1/2 DIFFUSER WIDTH (m)
!         BITOP=
!         BDET          	
!         C1=COEFFICIENT OF DIFERENT COUNTERFLOW ENTRAINMENTS
C	CO2=DISSOLVED OXYGEN (DO) CONCENTRATION (mol/m3)
C     CO2M=DO CONCENTRATION PROFILE FOR INPUT BOUNDARY CONDITION (g/m3)
C	CN2=DISSOLVED NITROGEN CONCENTRATION (mol/m3)
C	COMG=DISSOLVED OXYGEN CONCENTRATION (g/m3)
C     COMGPT=DO CONCENTRATION OF PLUME DETRAINMENT AT TOP OF PLUME (g/m3) 
!         COMGPD
!         COMGPTI
C	CNMG=DISSOLVED NITROGEN CONCENTRATION (g/m3)
C     DENSE20=DENSITY OF WATER AT 20 C
C	DENSEA=AMBIENT WATER DENSITY (kg/m3)
C	DENSEP=DENSITY OF THE PLUME (kg/m3)
C	DENSEW=WATER DENSITY IN PLUME (kg/m3)
!         DEPTHINTR=INTRUSION DEPTH OF THE DOUBLE PLUME
C	DIAMM=BUBBLE DIAMETER (mm)
C     DIFFEL=DIFFUSER ELEVATION (m)
C     DMPR=DEPTH OF MAXIMUM PLUME RISE (m)
C     DNAMB=AMBIENT DISSOLVED NITROGEN CONCENTRATION (g/m3)
C     DOAMB=AMBIENT DISSOLVED OXYGEN CONCENTRATION (g/m3)
C	E=ENTRAINMENT FACTOR (m3/s)
!	  EI=ENTRAINMENT  (m3/s)
!	  EO=ENTRAINMENT  (m3/s)
!         EA=ENTRAINMENT  (m3/s)
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
!	FR=FROUDE NUMBER (-)
C     FRCNOT=FRACTION OF NITROGEN IN ATMOSPHERE (-)
C     FRNO=INITIAL FROUDE NUMBER (-)
C     GAMMA=SALINITY CONVERSION FACTOR [kg/m3/(uS/cm)]
C     GROSSMT=GROSS MASS TRANSFER OF OXYGEN FROM PLUME (kg/d)
!	HO2=SOLUBILITY CONSTANT FOR OXYGEN (mol/m3/Pa)
!	HN2=SOLUBILITY CONSTANT FOR NITROGEN (mol/m3/Pa)
C     HCELL=HEIGHT OF CELL IN GRID (m)
C	HO=SOLUBILITY CONSTANT FOR OXYGEN (mol/m3/Pa)
C	HN=SOLUBILITY CONSTANT FOR NITROGEN (mol/m3/Pa)
C     HWITH=HEIGHT OF WITHDRAWAL/ENTRAINMENT ZONE (m)
!         INPLUME(:,6): INNER PLUME CHARACTERISTIC [HEIGTH; FLOWRATE OF WATER; MOMENTUM; TEMPERATURE; SALINITY; O2] JCT
C     JULDAY=JULIAN DAY IN GIVEN YEAR
C	KOLO=MASS TRANSFER COEFFICIENT FOR OXYGEN (m/s)
C	KOLN=MASS TRANSFER COEFFICIENT FOR NITROGEN (m/s)
C     LAKE=LAKE AND DIFFUSER TYPE (SHR AND LINEAR=1 OR AMISK AND RECTANGULAR=2) FOR SELECTION OF LAMBDA
C     LAYERS=NUMBER OF LAYERS/DATA POINTS IN BOUNDARY CONDITION PROFILES (-)
C     LAYDIFF=GRID LAYER CORRESPONDING TO DIFFUSER DEPTH (-)
!         LAYINTR=GRID LAYER CORRESPONDING TO INTRUSION DEPTH OF THE DOUBLE PLUME (-)
C     LAYTOP=GRID LAYER CORRESPONDING TO TOP OF PLUME (-)
C     LDIFF=LENGTH OF DIFFUSER (m)
C     LNOT=INITIAL DIFFUSER LENGTH (m)
C	LAMBDA=FRACTION OF PLUME OCCUPIED BY BUBBLES (-)
C     LAMBNOT=LAMBDA x INITIAL PLUME RADIUS;EQUAL TO DIFFUSER RADIUS (m) 
C     LOCAT=DEPTHS FOR INPUT BOUNDARY CONDITION PROFILES (m) 
C	MOMENT=MOMENTUM (m4/s)
C	N=NUMBER OF BUBBLES PER SECOND (1/s)
C     OTEFF=OXYGEN TRANSFER EFFICEINCY (%)
!         OUTPLUME(:,6)= OUTER PLUME CHARACTERISTIC [HEIGTH; FLOWRATE OF WATER; MOMENTUM; TEMPERATURE; SALINITY; O2] JCT
!         OUTERMOMENT= MOMENTUM IN THE OUTER PLUME 
!         OUTERQW=FLOWRATE OF WATER IN THE OUTER PLUME
!         OUTXLOC
!         OUTQW
!         OUTMOMENT
!         OUTTEMP
!         OUTO2  
C     PATM=ATMOSPHERIC PRESSURE AT AVERAGE WSEL (Pa)
C	PSTD=STANDARD PRESSURE (Pa)
!     QGAS=TOTAL GAS FLOW TO DIFFUSER
!       QINTOP
C!!!!	QSCFM=STANDARD GAS FLOW RATE (scfm), TOTAL STANDAR GAS FLOW RATE TO DIFFUSER 
C	QSCMS=STANDARD GAS FLOW RATE (scms)
C	QW=FLOWRATE OF WATER (m3/s)
!         QWDO
C     QWDET=TOTAL DETRAINMENT FLOW RATE AT THE INTRUSION OF THE PLUME (m3/s)  
C	RB=BUBBLE RADIUS (m)
C	RGAS=IDEAL GAS CONSTANT (J/mol/K)
C     SAL=SALINITY (uS/cm)
C     SALAMB=SALINITY OF AMBIENT WATER (uS/cm)
C	SALPLU=SALINITY OF THE PLUME (uS/cm)
!         SALPLUMED
!         SALPLUMETI
!     SHEARE=SHEAR ENTRAINMENT
C	  TARO=AROUND WATER TEMPERATURE (C)
C     TAVG=AVERAGE AMBIENT WATER TEMPERATURE (C)
C     TDS=TOTAL DISSOLVED SOLIDS (g/m3) [0.64 conversion factor from Chapra book]
C     TE=TEMPERATURE PROFILE FOR INPUT BOUNDARY CONDITION (C) 
C	TPLUME=PLUME TEMPERATURE (C)
C     TPLUMET=DETRAINMENT PLUME TEMPERATURE AT THE TOP OF THE PLUME (C) 
!         TPLUMED
!         TSTD=STANDARD TEMPERATURE (K)
!     UA = U AMBIENT VELOCITY PROFILE FOR INPUT BOUNDARY CONDITION
C     UAMB = U Ambient Velocity (m) used to calculate VORTEX entrainment (FJR
!	  VI=WATER VELOCITY (m/s)
!         VO=WATER VELOCITY (m/s)
!     VA = V AMBIENT VELOCITY PROFILE FOR INPUT BOUNDARY CONDITION
C     VAMB = V Ambient Velocity (m) used to calculate VORTEX entrainment (FJR
C     VAVG=AVERAGE WATER VELOCITY (m/s)
C	VB=BUBBLE RISE VELOCITY (m/s)
C	VBUB=BUBBLE VOLUME (m3)
!     VG=GAS VOLUME PER TOTAL VOLUME OF THE BUBBLE-WATER MIXTURE IN THE INNER CORE OF THE PLUME
C	VGUESS=GUESSED INITIAL WATER VELOCITY (m/s)
!         VI=WATER VELOCITY IN THE INNER PLUME (m/s)
!         VO=WATER VELOCITY IN THE OUTER PLUME (m/s)
!     VORTEXE=VORTEX ENTREAINMENT
C     WSEL=WATER SURFACE ELEVATION (m)
C	YO2=GASEOUS OXYGEN CONCENTRATION (mol/m3)
C	YN2=GASEOUS NITROGEN CONCENTRATION (mol/m3)
C	Z=DEPTH TO DIFFUSER (m)
C     PWD  = Perimeter (m) (FJR)
C     VUP  = Upward velocity (m/2)
C     AWD  = Area (m2) of plume
C     BWD  = Width 
C
      REAL*8 AREA,BI,BO,CO2,COMG,CN2,CNMG,DS,DZ,DENSEA,DENSEP,DENSEW,
     +DIAMM,FDO,FDN,FRACO,FRACN,FSAL,FTEMP,FGO,FGN,G,GAMMA,HO2,HN2,
     +KOLN,KOLO,LAMBDA,MOMENT,N,PI,PO,PN,PSTD,PZ,QSCMS,QW,QGAS,
     +RB,RGAS,SALAMB,SALPLU,TPLUME,TSTD,VB,VBUB,VG,VGUESS,
     +YO2,YN2,Z,AA,BB,CC,BNOT,LNOT,ELEV,DT,TE(1000),XLOC,DEPTH,
     +CO2M(1000),LOCAT(1000),DOAMB,COMGP,CNMGP,WSEL,PATM,SAL(1000),
     +DIFFEL,GROSSMT,FGONOT,DNAMB,DMPR,SUMSAL,FRCNATM,
     +LAMBNOT,FRCONOT,H,DYDX(5),Y(5),YOUT(5),BTOP,
     +DENSE20,OTEFF,FRNI,FRNO,VDIFF,FR,DCO2,QGFRAC,LDIFF,
     +JDAY,EL(70),DELTAC,COMGNOT,ELEVT,QWDET,QWDO(500),
     +PWDO(500),BWDO(500),AWDO(500),X,SHEARE,VORTEXE,
     +HWITH,HCELL,JULDAY,BDET,UA(1000),VA(1000),UAMB,VAMB,
     +INPLUME(35000,9),OUTPLUME(35000,9),DEPTHINTR,ALPHAI,ALPHAO,ALPHAA,
     +C1,EI,EO,EA,INNERQW,INNERMOMENT,INXLOC(35000),INQW(35000),
     +INMOMENT(35000),INTEMP(35000),INO2(35000),INSAL(35000),
     +INBI(35000),BITOP,VI,VO,TPLUMED,QINTOP,TARO,DEPTHDMPR,
     +SALPLUMED,COMGPD,DENSEINNER,GAMMA1,TIN,COMGIN,SALIN,EP,
     +INLI(35000),LITOP,LWDO(500),LI,LO
      INTEGER II,IJ,IK,JJ,LL,NEQN,NN,MI,JK,JL,LAYTOP,
     +LAYERS,KM,KN,KO,KP,ROWS,KQ,KR,KU,KV,KW,KX,KY,KZ,KS,LAYDIFF,YEAR
      INTEGER NELS,ierror,M,CONT,LAYINTR,NLI,NLO,NITERPLUME,NLIMAX ! ACC 2026
C
C     --------------------------------------------------------------------------
C     CONSTANTS
C     --------------------------------------------------------------------------
!

      PRINT*,'------OUTER_PLUME_RECT2-------'
      PRINT*,'NLI',NLI
      PRINT*,'NLO',NLO
      G=9.80665
      GAMMA=6.9E-4
      PI=ACOS(-1.0)
      PSTD=101325.
      RGAS=8.314
      TSTD=293.15
      DENSE20=998.2
      FRCNATM=0.00
      QGFRAC=1.0
!      C1=-1
      C1=0
!     
C     --------------------------------------------------------------------------
C     PARAMETERS AND INITIALIZE THE VARIABLES 
C     --------------------------------------------------------------------------
C     
      OUTPLUME(:,:)=0.0
!     Information from the INNER plume

      INXLOC   =INPLUME(:,1)
      INBI     =INPLUME(:,2)
	  INLI     =INPLUME(:,9) ! JCT_RECT_2022
      INQW     =INPLUME(:,3)
      INMOMENT =INPLUME(:,4)
      INTEMP   =INPLUME(:,5)
      INO2     =INPLUME(:,6)
      INSAL    =INPLUME(:,7)
	  

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Dejamos como zona de transicion los 0.1 metros superiores de la pluma. JCT
!     BITOP y QINTOP dejan de ser input, los sacamos de INPLUME --> ELIMINAR  JCT
!     ACC 2026 C5: offset fijo -10 reemplazado por búsqueda basada en profundidad (zona de transición ~0.1 m)
!      NLI = NLI -10                                                ! ACC 2026
      NLIMAX = NLI                                                  ! ACC 2026
      DO WHILE (NLI > 1)                                            ! ACC 2026
        IF ((INXLOC(NLIMAX)-INXLOC(NLI)) >= 0.1D0) EXIT           ! ACC 2026
        NLI = NLI - 1                                               ! ACC 2026
      END DO                                                         ! ACC 2026
      BITOP  = INPLUME(NLI,2)
      QINTOP = -INPLUME(NLI,3) !JCT2022
	  LITOP  = INPLUME(NLI,9) ! JCT_RECT	  
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
!     Geometric caracteristic 
      DEPTH=WSEL-DIFFEL
      DEPTHDMPR=WSEL-ELEVT+0.1 ! Restamos 0.1 metros para eliminar la zona de transicion
      Z=DEPTHDMPR   
      ELEV=ELEVT
!        
!!C     Assume that gas bubbles are composed of oxygen and nitrogen only.
!!      FRACO=FRCONOT 
!!      FRACN=1.0-FRACO
!
C     Interpolate input profiles to obtain the plume initial conditions
      X=0.0 
      XLOC=Z-X 

      CALL LININT(LOCAT,TE,LAYERS,XLOC,TARO)
      TIN=INTEMP(NLI)
      CALL LININT(LOCAT,CO2M,LAYERS,XLOC,COMG)
      COMGIN=INO2(NLI)
      CALL LININT(LOCAT,UA,LAYERS,XLOC,UAMB)
      CALL LININT(LOCAT,VA,LAYERS,XLOC,VAMB)
      SALIN=INSAL(NLI)
      DOAMB=COMG
      CO2=COMG/32.
      COMGNOT=COMGP
!     At first temperature, DO and salinity in the outer plume is the same that in the top of the inner plume
      TPLUME=TIN
      COMGP=COMGIN
      SALPLU=SALIN  
      SALPLU=224  ! Eliminar JCT SAlinidad

!     Guess the initial velocity 
      VGUESS=0.07
      VO=VGUESS
C      
C     AMBIENT AND AVERAGE WATER DENSITIES

      DENSEA=(0.059385*TARO**3-8.56272*TARO**2+65.4891*TARO)*0.001
     ++999.84298 !+(GAMMA)*SALAMB
!!      DENSEW=DENSEA
      DENSEP=(0.059385*TPLUME**3-8.56272*TPLUME**2+65.4891*TPLUME)*0.001
     ++999.84298 !+(GAMMA)*SALPLU
      DENSEINNER=(0.059385*TIN**3-8.56272*TIN**2+65.4891*TIN)*0.001
     ++999.84298 !+(GAMMA)*SALIN
      ! ELIMI SALINIDAD JCT2022

C     --------------------------------------------------------------------------
C     PLUME PROPERTIES
C     --------------------------------------------------------------------------
!
C     CALCULATION OF INITIAL WATER VELOCITY USING FROUDE NUMBER

      VDIFF=1
	   !PRINT*,"QINTOP",QINTOP
       !PRINT*,"VO",VO
       !PRINT*,"LITOP",LITOP
       !PRINT*,"BITOP",BITOP
      DO WHILE (VDIFF.GT.1.0E-6)
         AREA = QINTOP/VO+(LITOP*2.*BITOP)	  
         !SOLVE FOR DIMENSIONS USING L^2+(2Bo-Lo)L-AREA=0 USING QUADRATIC EQN.
         AA=1.0
         BB=2.*BITOP-LITOP
         CC=-1.0*AREA
         LO=(-1.0*BB+(BB**2-4.0*AA*CC)**(0.5))/(2.0*AA)
         IF(LO.LT.0.0)THEN
            LO=(-1.0*BB-(BB**2-4.0*AA*CC)**(0.5))/(2.0*AA)
         ENDIF		 
		 BO=AREA/(2.0*LO)
!        ACC 2026: denominador DENSEA fisicamente correcto (pluma externa mas densa que ambiente)
!        Linea erronea comentada: VO=-FRNO*(ABS((BO-BITOP)*G*(DENSEA-DENSEP)/DENSEP))**0.5
         VO=-FRNO*(ABS((BO-BITOP)*G*(DENSEA-DENSEP)/DENSEP))**0.5  ! ACC 2026
         VDIFF=ABS(VO-VGUESS)
         VGUESS=VO
      END DO

!     Water velocity in the top of the inner plume
	  VI=-QINTOP/(LITOP*BITOP*2) ! JCT_RECT JCT_2022

       PRINT*,"VI,VO,BO,LO",VI, VO,BO,LO
	  
C     ------------------------------------------------------------------
C     VARIABLE TRANSFORMATION      
C     ------------------------------------------------------------------

      BI=BITOP
      LI=LITOP ! JCT_RECT_2022
	  
!      EA=-(2.*(LO+2.*BO))*ALPHAA*VO ! JCT_RECT JCT_2022
	  EI=+2.*LI*ALPHAI*(VI+C1*VO)  ! ACC 2026: solo lados largos (2D linea fuente)
      EO=-2.*LI*ALPHAO*VO           ! ACC 2026: solo lados largos
      EA=-2.*LO*ALPHAA*VO           ! ACC 2026: solo lados largos
      EP=QINTOP

      PRINT*,"EA,EO,EI,EP,VO",EA,EO,EI,EP, VO

      QW=QINTOP
      MOMENT=QW*VO
      PRINT*,"MOMENTinicial",MOMENT,VO,QW
	  FTEMP=QW*TPLUME 
      FSAL=QW*SALPLU !JCT2020_Sal
C     Previous equation corrected to account for salinity units conversion.        
      FDO=QW*CO2
!!    FDN=QW*CN2
!     Initialize lateral withdrawal flowrate for first/lowest cell in column/segment
      JJ=0
      QWDO(LAYTOP)= QW
	  PWDO(LAYTOP)= (2.*LO)  ! ACC 2026: solo lados largos (consistente con EA, Dissanayake 2021)
      AWDO(LAYTOP)= (LO*BO*2)-(LI*BI*2)   ! JCT_RECT_2022
      BWDO(LAYTOP)= BO
      LWDO(LAYTOP)= LO
!
!     -------------------------------------------------------------------------- 
C	SOLUTION PROCEEDURE
!     --------------------------------------------------------------------------
!
      DZ=0.01
      ! H=0.001
      H=-0.01 ! JCT_RECT_2022 
      NELS = 0
      HWITH=0.0
      NLO=NLI+1
!
      M=0
      DO WHILE (VO.LT.-1.E-6.AND.Z.LT.DEPTH)
         M=M+1
         Z=Z+DZ
         X=X+DZ
         ELEV=ELEV-DZ
         NELS = NELS + 1
         NLO=NLO-1
C	
C        Interpolate input profiles to obtain line plume boundary conditions
         XLOC=DEPTHDMPR+X

         CALL LININT(LOCAT,TE,LAYERS,XLOC,TARO)
         TIN=INTEMP(NLO)
         CALL LININT(LOCAT,CO2M,LAYERS,XLOC,COMG)
         COMGIN=INO2(NLO)
         SALIN=INSAL(NLO)
		 !SALIN = 224
         CALL LININT(LOCAT,UA,LAYERS,XLOC,UAMB)
         CALL LININT(LOCAT,VA,LAYERS,XLOC,VAMB)
C        Inner plume width
         BI=INBI(NLO)
		 LI=INLI(NLO) ! JCT_RECT_2022
C        Inner velocity plume
         INNERQW=INQW(NLO)
         INNERMOMENT=INMOMENT(NLO)
         VI=INNERMOMENT/INNERQW

C        Use subroutines for Runge Kutta method solution
         NEQN=5
         Y(1)=QW
         Y(2)=MOMENT
         Y(3)=FTEMP
         Y(4)=FSAL
         Y(5)=FDO


C         CALL DERIVS_7(EA,EO,EI,VI,VO,DENSEA,DENSEP,G,PI,BO,BI,LO,LI,
C     +             TARO,TIN,TPLUME,SALAMB,SALIN,SALPLU,GAMMA,DENSE20,
C     +             DENSEINNER,GAMMA1,DOAMB,COMGIN,COMGP,DNAMB,
C     +             CNMGP,Z,Y,DYDX)
C
C         CALL RK4_7(EA,EO,EI,VI,VO,DENSEA,DENSEP,G,PI,BO,BI,LO,LI,TARO,
C     +             TIN,TPLUME,SALAMB,SALIN,SALPLU,GAMMA,DENSE20,
C     +             DENSEINNER,GAMMA1,DOAMB,COMGIN,COMGP,DNAMB,CNMGP,Y,
C     +             DYDX,NEQN,X,H,YOUT) 

         CALL DERIVS_8(EA,EO,EI,VI,VO,DENSEA,DENSEP,G,PI,BO,BI,TARO,TIN,
     +             TPLUME,SALAMB,SALIN,SALPLU,GAMMA,DENSE20,LO,LI,
     +             GAMMA1,DOAMB,COMGIN,COMGP,DNAMB,CNMGP,Z,Y,DYDX)

         CALL RK4_8(EA,EO,EI,VI,VO,DENSEA,DENSEP,G,PI,BO,BI,TARO,TIN,
     +             TPLUME,SALAMB,SALIN,SALPLU,GAMMA,DENSE20,LO,LI,
     +             GAMMA1,DOAMB,COMGIN,COMGP,DNAMB,CNMGP,Y,DYDX,NEQN,X,
     +             H,YOUT)	 
C       
         QW=YOUT(1)
         MOMENT=YOUT(2)
         FTEMP=YOUT(3)
         FSAL=YOUT(4)
         FDO=YOUT(5) 
	 
          ! OPEN (UNIT=59, FILE="outer_plume_eq.txt", POSITION="APPEND")
             ! WRITE (UNIT=59, FMT='(6F12.6)') QW,MOMENT,FTEMP,FSAL,FDO,X
          ! CLOSE (UNIT=59)

         IF(MOMENT.LT.0.0)THEN
            TPLUME=FTEMP/QW
	    !SALPLU=FSAL/(QW*DENSEP)/(GAMMA/DENSE20)
            SALPLU=224  ! Eliminar JCT SAlinidad
C           Previous equation corrected to consistently express salinity in uS/cm	   
	        CO2=FDO/QW
            COMGP=CO2*32.
            VO=MOMENT/QW
			
			
            AREA = QW/VO+(LI*2.*BI)
            !AREA = QW/VO	  
            !SOLVE FOR DIMENSIONS USING L^2+(2Bo-Lo)L-AREA=0 USING QUADRATIC EQN.
            AA=1.0
            BB=2.*BI-LI
            CC=-1.0*AREA
            LO=(-1.0*BB+(BB**2-4.0*AA*CC)**(0.5))/(2.0*AA)
            IF(LO.LT.0.0)THEN
               LO=(-1.0*BB-(BB**2-4.0*AA*CC)**(0.5))/(2.0*AA)
            ENDIF		 
		    BO=AREA/(2.0*LO)
			
			
!           Save outer plume information
            OUTPLUME(NLO,1)=XLOC 
            OUTPLUME(NLO,2)=BO          
            OUTPLUME(NLO,3)=QW
            OUTPLUME(NLO,4)=MOMENT
            OUTPLUME(NLO,5)=TPLUME
            OUTPLUME(NLO,6)=COMGP
            OUTPLUME(NLO,7)=SALPLU
            OUTPLUME(NLO,8)=VO
            OUTPLUME(NLO,9)=LO

	    GOTO 20
        ENDIF

         VO=MOMENT/QW
            AREA = QW/VO+(LI*2.*BI)   
            !AREA = QW/VO
            !SOLVE FOR DIMENSIONS USING L^2+(2Bo-Lo)L-AREA=0 USING QUADRATIC EQN.
            AA=1.0
            BB=2.*BI-LI
            CC=-1.0*AREA
            LO=(-1.0*BB+(BB**2-4.0*AA*CC)**(0.5))/(2.0*AA)
            IF(LO.LT.0.0)THEN
               LO=(-1.0*BB-(BB**2-4.0*AA*CC)**(0.5))/(2.0*AA)
            ENDIF		 
		    BO=AREA/(2.0*LO)

         EI=+2.*LI*ALPHAI*(VI+C1*VO)  ! ACC 2026: solo lados largos
         EO=-2.*LI*ALPHAO*VO           ! ACC 2026: solo lados largos
         EA=-2.*LO*ALPHAA*VO           ! ACC 2026: solo lados largos

        !
!        Temperature and salinity in the plume
         TPLUME=FTEMP/QW
         SALPLU=FSAL/(QW*DENSEP)/(GAMMA/DENSE20)
         SALPLU=224  ! Eliminar JCT SAlinidad
C        Previous equation corrected to consistently express salinity in uS/cm

!        Dissolved oxygen and nitrogen concentration
         CO2=FDO/QW
!         CN2=FDN/QW 
         COMGP=CO2*32.
!         CNMGP=CN2*28.

!        Save outer plume information
         OUTPLUME(NLO,1)=XLOC 
         OUTPLUME(NLO,2)=BO          
         OUTPLUME(NLO,3)=QW
         OUTPLUME(NLO,4)=MOMENT
         OUTPLUME(NLO,5)=TPLUME
         OUTPLUME(NLO,6)=COMGP
         OUTPLUME(NLO,7)=SALPLU
         OUTPLUME(NLO,8)=VO
         OUTPLUME(NLO,9)=LO
		 
!        Add incremental entrainment to total cell entrainment/withdrawal
         QWDO(LAYTOP+JJ)=QWDO(LAYTOP+JJ)-(EA+EO-EI)*DZ ! JCT_2022
         PWDO(LAYTOP+JJ)=PWDO(LAYTOP+JJ)+(2.*LO)  ! ACC 2026: solo lados largos (consistente con EA)
         AWDO(LAYTOP+JJ)=AWDO(LAYTOP+JJ)+(LO*BO*2) - (LI*BI*2)
         BWDO(LAYTOP+JJ)=BWDO(LAYTOP+JJ)+BO
		 LWDO(LAYTOP+JJ)=LWDO(LAYTOP+JJ)+LO
         HWITH=HWITH+DZ   
         IF(HWITH.GT.HCELL)THEN
            PWDO(LAYTOP+JJ) = PWDO(LAYTOP+JJ)/NELS
            AWDO(LAYTOP+JJ) = AWDO(LAYTOP+JJ)/NELS
            BWDO(LAYTOP+JJ) = BWDO(LAYTOP+JJ)/NELS
            LWDO(LAYTOP+JJ) = LWDO(LAYTOP+JJ)/NELS
            JJ=JJ+1
            HWITH=0.0
            QWDO(LAYTOP+JJ)=0.0
            PWDO(LAYTOP+JJ)=0.0
            BWDO(LAYTOP+JJ)=0.0
            AWDO(LAYTOP+JJ)=0.0
            LWDO(LAYTOP+JJ)=0.0
            NELS = 0
         ENDIF

!        Density (ambient water, water in plume and of the plume)
      DENSEA=(0.059385*TARO**3-8.56272*TARO**2+65.4891*TARO)*0.001
     ++999.84298 !+(GAMMA)*SALAMB
      DENSEP=(0.059385*TPLUME**3-8.56272*TPLUME**2+65.4891*TPLUME)*0.001
     ++999.84298 !+(GAMMA)*SALPLU
      DENSEINNER=(0.059385*TIN**3-8.56272*TIN**2+65.4891*TIN)*0.001
     ++999.84298 !+(GAMMA)*SALIN

      END DO
C
!     --------------------------------------------------------------------
C     CALCULATION OF AVERAGE NET OXYGEN MASS TRANSFER FOR DAY
!     --------------------------------------------------------------------

   20 DEPTHINTR=ELEV
      QWDET=QW
      TPLUMED=TPLUME
      COMGPD=COMGP
      LAYINTR=LAYTOP+JJ
      RETURN
      END
	  
	  
	  

C     
C----------------------------------------------------------------------
C
      SUBROUTINE RK4_8(EA,EO,EI,VI,VO,DENSEA,DENSEP,G,PI,BO,BI,
     +             TARO,TIN,TPLUME,SALAMB,SALIN,SALPLU,GAMMA,
     +             DENSE20,LO,LI,GAMMA1,DOAMB,COMGIN,COMGP,
     +             DNAMB,CNMGP,Y,DYDX,NN,X,H,YOUT)

      INTEGER I,NN,NMAX
      PARAMETER (NMAX=50)
      REAL*8 EA,EO,EI,DENSEA,DENSEP,G,PI,BO,BI,TARO,TIN,TPLUME,
     +SALAMB,SALIN,SALPLU,GAMMA,DENSE20,DOAMB,COMGP,
     +DNAMB,CNMGP,H,X,DYDX(NN),Y(NN),YOUT(NN),H6,HH,XH,
     +DYM(NMAX),DYT(NMAX),YT(NMAX),GAMMA1,VI,VO,COMGIN,LO,LI

      EXTERNAL DERIVS_8
      HH=H*0.5
      H6=H/6.
      XH=X+HH
      DO 11 I=1,NN
          YT(I)=Y(I)+HH*DYDX(I)
   11 CONTINUE
      CALL DERIVS_8(EA,EO,EI,VI,VO,DENSEA,DENSEP,G,PI,BO,BI,TARO,TIN,
     +             TPLUME,SALAMB,SALIN,SALPLU,GAMMA,DENSE20,LO,LI,
     +             GAMMA1,DOAMB,COMGIN,COMGP,DNAMB,CNMGP,XH,YT,DYT)
      DO 12 I=1,NN
          YT(I)=Y(I)+HH*DYT(I)
   12 CONTINUE
      CALL DERIVS_8(EA,EO,EI,VI,VO,DENSEA,DENSEP,G,PI,BO,BI,TARO,TIN,
     +             TPLUME,SALAMB,SALIN,SALPLU,GAMMA,DENSE20,LO,LI,
     +             GAMMA1,DOAMB,COMGIN,COMGP,DNAMB,CNMGP,XH,YT,DYM)
      DO 13 I=1,NN
          YT(I)=Y(I)+H*DYM(I)
          DYM(I)=DYT(I)+DYM(I)
   13 CONTINUE
      CALL DERIVS_8(EA,EO,EI,VI,VO,DENSEA,DENSEP,G,PI,BO,BI,TARO,TIN,
     +             TPLUME,SALAMB,SALIN,SALPLU,GAMMA,DENSE20,LO,LI,
     +             GAMMA1,DOAMB,COMGIN,COMGP,DNAMB,CNMGP,X+H,YT,DYT)
      DO 14 I=1,NN
          YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
   14 CONTINUE
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE DERIVS_8(EA,EO,EI,VI,VO,DENSEA,DENSEP,G,PI,BO,BI,
     +             TARO,TIN,TPLUME,SALAMB,SALIN,SALPLU,GAMMA,DENSE20,
     +             LO,LI,GAMMA1,DOAMB,COMGIN,COMGP,DNAMB,CNMGP,
     +             X,Y,DYDX)         
      
      REAL*8 EA,EO,EI,DENSEA,DENSEP,G,PI,BO,BI,TARO,TIN,TPLUME,SALAMB,
     +SALIN,SALPLU,GAMMA,DENSE20,DOAMB,COMGP,GAMMA1,
     +DNAMB,CNMGP,X,Y(5),DYDX(5),VI,VO,COMGIN,MOM1,MOM2,MOM3,LO,LI

      DYDX(1)=EA+EO-EI
!	  DYDX(2)=((1/GAMMA1)*(-PI*G*(BO**2-BI**2)*
!     +((DENSEP-DENSEA)/DENSE20))-EI*VO+EO*VI)
!     ACC 2026 (DERIVS_8): ERROR - densidad de referencia era DENSE20 (constante 998.2 kg/m3).
!     Debe ser DENSEA (densidad ambiente variable). Lineas erroneas comentadas:
!      DYDX(2)=((1/GAMMA1)*(-G*(LO*2*BO-LI*2*BI)*
!     +((DENSEP-DENSEA)/DENSE20))-EI*VO+EO*VI)
      DYDX(2)=((1/GAMMA1)*(-G*(LO*2*BO-LI*2*BI)*
     +((DENSEP-DENSEA)/DENSE20))-EI*VO+EO*VI)  ! ACC 2026
      DYDX(3)=EA*TARO+EO*TIN-EI*TPLUME
      DYDX(4)=EA*SALAMB+EO*SALIN-EI*SALPLU
      DYDX(5)=EA*DOAMB/32.+EO*COMGIN/32.-EI*COMGP/32.
 !     DYDX(6)=EA*DNAMB/28.+EO*???/28.-EI*CNMGP/28.
      RETURN
      END


C     
C     
C----------------------------------------------------------------------
C	  
	  
	  
C     
C----------------------------------------------------------------------
C
      SUBROUTINE RK4_7(EA,EO,EI,VI,VO,DENSEA,DENSEP,G,PI,BO,BI,
     +             LO,LI,TARO,TIN,TPLUME,SALAMB,SALIN,SALPLU,GAMMA,
     +             DENSE20,DENSEINNER,GAMMA1,DOAMB,COMGIN,COMGP,
     +             DNAMB,CNMGP,Y,DYDX,NN,X,H,YOUT)

      INTEGER I,NN,NMAX
      PARAMETER (NMAX=50)
      REAL*8 EA,EO,EI,DENSEA,DENSEP,G,PI,BO,BI,LO,LI,TARO,TIN,TPLUME,
     +SALAMB,SALIN,SALPLU,GAMMA,DENSE20,DENSEINNER,DOAMB,COMGP,
     +DNAMB,CNMGP,H,X,DYDX(NN),Y(NN),YOUT(NN),H6,HH,XH,
     +DYM(NMAX),DYT(NMAX),YT(NMAX),GAMMA1,VI,VO,COMGIN

      EXTERNAL DERIVS_7
      HH=H*0.5
      H6=H/6.
      XH=X+HH
	  
	        ! PRINT *, 'DENSEA_rk',DENSEA


      DO 11 I=1,NN
          YT(I)=Y(I)+HH*DYDX(I)
   11 CONTINUE
      CALL DERIVS_7(EA,EO,EI,VI,VO,DENSEA,DENSEP,G,PI,BO,BI,LO,LI,TARO,
     +             TIN,TPLUME,SALAMB,SALIN,SALPLU,GAMMA,DENSE20,
     +             DENSEINNER,GAMMA1,DOAMB,COMGIN,COMGP,DNAMB,CNMGP,XH,
     +             YT,DYT)
      DO 12 I=1,NN
          YT(I)=Y(I)+HH*DYT(I)
   12 CONTINUE
      CALL DERIVS_7(EA,EO,EI,VI,VO,DENSEA,DENSEP,G,PI,BO,BI,LO,LI,TARO,
     +             TIN,TPLUME,SALAMB,SALIN,SALPLU,GAMMA,DENSE20,
     +             DENSEINNER,GAMMA1,DOAMB,COMGIN,COMGP,DNAMB,CNMGP,XH,
     +             YT,DYT)
      DO 13 I=1,NN
          YT(I)=Y(I)+H*DYM(I)
          DYM(I)=DYT(I)+DYM(I)
   13 CONTINUE
      CALL DERIVS_7(EA,EO,EI,VI,VO,DENSEA,DENSEP,G,PI,BO,BI,LO,LI,TARO,
     +             TIN,TPLUME,SALAMB,SALIN,SALPLU,GAMMA,DENSE20,
     +             DENSEINNER,GAMMA1,DOAMB,COMGIN,COMGP,DNAMB,CNMGP,XH,
     +             YT,DYT)
      DO 14 I=1,NN
          YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
   14 CONTINUE
      RETURN
      END
C
C----------------------------------------------------------------------
      SUBROUTINE DERIVS_7(EA,EO,EI,VI,VO,DENSEA,DENSEP,G,PI,BO,BI,LO,
     +             LI,TARO,TIN,TPLUME,SALAMB,SALIN,SALPLU,GAMMA,
     +             DENSE20,DENSEINNER,GAMMA1,DOAMB,COMGIN,COMGP,DNAMB,
     +             CNMGP,X,Y,DYDX)         
      
      REAL*8 EA,EO,EI,DENSEA,DENSEP,G,PI,BO,BI,LO,LI,TARO,TIN,TPLUME,
     +SALAMB,SALIN,SALPLU,GAMMA,DENSE20,DENSEINNER,DOAMB,COMGP,GAMMA1,
     +DNAMB,CNMGP,X,Y(5),DYDX(5),VI,VO,COMGIN,MOM1,MOM2,MOM3

    !  PRINT *, 'EA',EA
    !  PRINT *, 'EO',EO
    !  PRINT *, 'EI',EI
      ! PRINT *, 'DENSEA_derivs',DENSEA
    !  PRINT *, 'DENSEP',DENSEP
    !  PRINT *, 'G',G
    !  PRINT *, 'BI',BI
    !  PRINT *, 'BI',BI
    !  PRINT *, 'BI',BI
    !  PRINT *, 'LAMBDA',LAMBDA
    !   PRINT *, 'TARO',TARO
    !   PRINT *, 'TINderivs',TIN
    !   PRINT *, 'TPLUME',TPLUME
    !  PRINT *, 'SALAMB',SALAMB
    !  PRINT *, 'GAMMA',GAMMA
    !  PRINT *, 'DENSE20',DENSE20
    !  PRINT *, 'DOAMB',DOAMB
    !  PRINT *, 'PI',PI
    !  PRINT *, 'RB',RB
    !  PRINT *, 'N',N
    !  PRINT *, 'VI',VI
    !  PRINT *, 'VB',VB
    !  PRINT *, 'KOLO',KOLO
    !  PRINT *, 'HO2',HO2
    !  PRINT *, 'PO',PO
    !  PRINT *, 'COMGP',COMGP
    !  PRINT *, 'DNAMB',DNAMB
    !  PRINT *, 'KOLN',KOLN
    !  PRINT *, 'HN2',HN2
    !  PRINT *, 'PN',PN
    !  PRINT *, 'CNMGP',CNMGP
    !  PRINT *, 'Z',Z
    !  PRINT *, 'Y',Y
    !  PRINT *, 'DYDX',DYDX

      DYDX(1)=EA+EO-EI
!      DYDX(2)=PI*G*(BO**2-BI**2)*((DENSEP-DENSEA)/(DENSE20*GAMMA1))+
!     +EI*VO-EO*VI

 !     DYDX(2)=(1/GAMMA1)*(PI*G*(BO**2-BI**2)*((DENSEP-DENSEA)/DENSE20)+
 !    +EI*VO-EO*VI)*(-1)

!     ACC 2026 (DERIVS_7): ERROR - densidad de referencia era DENSE20 (constante 998.2 kg/m3)
!     en lugar de DENSEA (densidad ambiente variable con profundidad). Inconsistente con
!     DERIVS_6 que usa DENSEP. Para conservacion de momento de pluma externa usar DENSEA.
!     Referencia: Wuest (1992), Socolofsky (2008). Lineas erroneas comentadas:
!      DYDX(2)=((1/GAMMA1)*(-G*(LO*2*BO-LI*2*BI)*
!     +((DENSEP-DENSEA)/DENSE20))-EI*VO+EO*VI)
      DYDX(2)=((1/GAMMA1)*(-G*(LO*2*BO-LI*2*BI)*
     +((DENSEP-DENSEA)/DENSE20))-EI*VO+EO*VI)  ! ACC 2026
	 
	    ! PRINT*," DYDX(2)_drevis",DYDX(2), GAMMA1,G,LO,BO,LI,BI
		! PRINT*,DENSEP,DENSEA,DENSE20,EI,VO,EO,VI
		

!      DYDX(2)=(-1/GAMMA1)*(PI*G*(BO**2-BI**2)*((DENSEP-DENSEA)/DENSE20))
!      DYDX(2)=(1/GAMMA1)*(PI*G*(BO**2-BI**2)*((DENSEP-DENSEA)/DENSE20)+
!     +0)*(-1)
      !MOM1=-(1/GAMMA1)*(PI*G*(BO**2-BI**2)*((DENSEP-DENSEA)/(DENSE20))
      !MOM2=-(1/GAMMA1)*EI*VO
      !MOM3=+(1/GAMMA1)*EO*VI
      !PRINT *, 'M', EI, VO, EO, VI
      !PRINT *, 'M', DYDX(2), MOM1, MOM2, MOM3, EO, VI
      DYDX(3)=EA*TARO+EO*TIN-EI*TPLUME
    !  DYDX(4)=EA*(SALAMB*GAMMA/DENSE20)*DENSEA+
   !  +EO*(SALIN*GAMMA/DENSE20)*DENSEINNER-
   !  +EI*(SALPLU*GAMMA/DENSE20)*DENSEP
      DYDX(4)=EA*SALAMB+EO*SALIN-EI*SALPLU
      DYDX(5)=EA*DOAMB/32.+EO*COMGIN/32.-EI*COMGP/32.
 !     DYDX(6)=EA*DNAMB/28.+EO*???/28.-EI*CNMGP/28.
      RETURN
      END 

C     
C     
C----------------------------------------------------------------------
C

      SUBROUTINE RK4_6(EI,EO,DENSEA,DENSEW,DENSEP,G,BI,LI,LAMBDA,TARO,
     +              VG,SALARO,GAMMA,DENSE20,DOAMB,PI,RB,N,VI,VO,VB,KOLO,
     +              HO2,PO,GAMMA1,TPLUME,SALPLU,COMGP,DNAMB,KOLN,HN2,PN,
     +              CNMGP,Y,DYDX,NN,X,H,YOUT,XLOC,TAMB)

      INTEGER I,NN,NMAX
      PARAMETER (NMAX=50)
      REAL*8 EI,EO,DENSEA,DENSEW,DENSEP,G,BI,LI,LAMBDA,TARO,VG,SALARO,
     +GAMMA,DENSE20,DOAMB,PI,RB,N,VI,VO,VB,KOLO,HO2,PO,COMGP,DNAMB,
     +KOLN,HN2,PN,TPLUME,CNMGP,H,X,DYDX(NN),Y(NN),YOUT(NN),H6,HH,XH,
     +DYM(NMAX),DYT(NMAX),YT(NMAX),GAMMA1,SALPLU,XLOC,TAMB
      EXTERNAL DERIVS_6


      !PRINT*, "A1", EI,EO,DENSEA
      !PRINT*, "A2", DENSEW,DENSEP,G
      !PRINT*, "A3", BI,LI,LAMBDA
      !PRINT*, "A4", TARO,VG,SALARO
      !PRINT*, "A5", GAMMA,DENSE20,DOAMB
      !PRINT*, "A6", PI,RB,N
      !PRINT*, "A7", VI,VO,VB
      !PRINT*, "A8", KOLO,HO2,PO
      !PRINT*, "a9", GAMMA1,TPLUME,SALPLU
      !PRINT*, "A10", COMGP,DNAMB,KOLN
      !PRINT*, "a11", HN2,PN,CNMGP
      !PRINT*, "A12", Y,DYDX,NN
      !PRINT*, "A13", X,H,YOUT
      !PRINT*, "A14", XLOC,TAMB

      HH=H*0.5
      H6=H/6.
      XH=X+HH

      !NN=8

      !PRINT*, "NNa", NN
      DO 11 I=1,NN
          YT(I)=Y(I)+HH*DYDX(I)
   11 CONTINUE
      CALL DERIVS_6(EI,EO,DENSEA,DENSEW,DENSEP,G,BI,LI,LAMBDA,TARO,VG,
     +           SALARO,GAMMA,DENSE20,DOAMB,PI,RB,N,VI,VO,VB,KOLO,
     +           HO2,PO,GAMMA1,TPLUME,SALPLU,COMGP,DNAMB,KOLN,HN2,
     +           PN,CNMGP,XH,YT,DYT,XLOC,TAMB)


      DO 12 I=1,NN
          YT(I)=Y(I)+HH*DYT(I)
   12 CONTINUE
      CALL DERIVS_6(EI,EO,DENSEA,DENSEW,DENSEP,G,LI,BI,LAMBDA,TARO,VG,
     +            SALARO,GAMMA,DENSE20,DOAMB,PI,RB,N,VI,VO,VB,KOLO,
     +            HO2,PO,GAMMA1,TPLUME,SALPLU,COMGP,DNAMB,KOLN,HN2,
     +            PN,CNMGP,XH,YT,DYM,XLOC,TAMB)
      DO 13 I=1,NN
          YT(I)=Y(I)+H*DYM(I)
          DYM(I)=DYT(I)+DYM(I)
   13 CONTINUE
      CALL DERIVS_6(EI,EO,DENSEA,DENSEW,DENSEP,G,BI,LI,LAMBDA,TARO,VG,
     +            SALARO,GAMMA,DENSE20,DOAMB,PI,RB,N,VI,VO,VB,KOLO,
     +            HO2,PO,GAMMA1,TPLUME,SALPLU,COMGP,DNAMB,KOLN,HN2,
     +            PN,CNMGP,X+H,YT,DYT,XLOC,TAMB)
      DO 14 I=1,NN
          YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
   14 CONTINUE
      RETURN
      END
C
C
C----------------------------------------------------------------------
C
      SUBROUTINE DERIVS_6(EI,EO,DENSEA,DENSEW,DENSEP,G,BI,LI,LAMBDA,
     +              TARO,VG,SALARO,GAMMA,DENSE20,DOAMB,PI,RB,N,
     +              VI,VO,VB,KOLO,HO2,PO,GAMMA1,TPLUME,SALPLU,
     +              COMGP,DNAMB,KOLN,HN2,PN,CNMGP,X,Y,DYDX,
     +              XLOC,TAMB)         

      REAL*8 EI,EO,DENSEA,DENSEW,DENSEP,G,LI,BI,LAMBDA,TARO,SALARO,
     +GAMMA,DENSE20,DOAMB,PI,RB,N,VI,VO,VB,KOLO,HO2,PO,COMGP,SALPLU,
     +DNAMB,KOLN,HN2,PN,TPLUME,CNMGP,X,Y(8),DYDX(8),VG,GAMMA1,MOM1,
     +MOM2,MOM3,MOM4,MOM0,XLOC,TAMB

      !PRINT *, BI,LI

      !PRINT*, "D1", EI,EO,DENSEA
      !PRINT*, "D2", DENSEW,DENSEP,G
      !PRINT*, "D3", BI,LI,LAMBDA
      !PRINT*, "D4", TARO,VG,SALARO
      !PRINT*, "D5", GAMMA,DENSE20,DOAMB
      !PRINT*, "D6", PI,RB,N
      !PRINT*, "D7", VI,VO,VB
      !PRINT*, "D8", KOLO,HO2,PO
      !PRINT*, "D9", GAMMA1,TPLUME,SALPLU
      !PRINT*, "D10", COMGP,DNAMB,KOLN
      !PRINT*, "D11", HN2,PN,CNMGP
      !PRINT*, "D12", Y,DYDX,NN
      !PRINT*, "D13", X,H,YOUT
      !PRINT*, "D14", XLOC,TAMB


C     Right-hand side of differential equations for Runge-Kutta solution
      DYDX(1)=EI-EO      
!      PRINT *, EI,EO  
C      DYDX(2)=(DENSEA-DENSEW)/DENSEP*G*(LI*2.0*BI)+(DENSEW-DENSEP)/
C     +DENSEP*G*(LAMBDA*LI*2*LAMBDA*BI)

!      Hasta ahora he utilizado el de abajo, pero creo k esta mal 16-5-13
!      DYDX(2)=((DENSEA-DENSEW)/DENSEP)*G*PI*BI**2*(1-LAMBDA**2)+
!     +(DENSEW-DENSEP)/DENSEP*G*(PI*(LAMBDA*BI)**2)

!     Comporbar si cambian los resultados  16-5-13
!      DYDX(2)=((DENSEA-DENSEW)/DENSEP)*G*PI*BI**2*(1-LAMBDA**2)+
!     +(DENSEA-DENSEP)/DENSEP*G*(PI*(LAMBDA*BI)**2)


 !     DYDX(2)=(1/GAMMA1)*((PI*G*BI**2/DENSE20)*(LAMBDA**2*VG*(DENSEA-0)+
 !    +LAMBDA**2*(1-VG)*(DENSEA-DENSEW))+EI*VO-EO*VI)

!     ACC 2026: ERROR - area de burbujas usaba LAMBDA^2 (formulacion circular/axisimetrica)
!     en lugar de LAMBDA (formulacion rectangular/linea fuente). Para pluma 2D: solo el
!     semiancho BI escala con LAMBDA; la fraccion de area con burbujas es LAMBDA*Area,
!     no LAMBDA^2*Area. Referencia: Dissanayake (2021). Lineas erroneas comentadas:
!      DYDX(2)=(1/GAMMA1)*(((DENSEA-DENSEW)/DENSEP)*G*(LI*2.0*BI)*
!     +(1-LAMBDA**2)+((DENSEA-DENSEP)/DENSEP)*G*(LI*2.0*BI)*
!     +LAMBDA**2)+EI*VO-EO*VI

!      DYDX(2)=(1/GAMMA1)*(((DENSEA-DENSEW)/DENSEP)*G*(LI*2.0*BI)*
!     +(1-LAMBDA)+((DENSEA-DENSEP)/DENSEP)*G*(LI*2.0*BI)*
!     +LAMBDA)+EI*VO-EO*VI  ! ACC 2026

      DYDX(2)=(1/GAMMA1)*(((DENSEA-DENSEW)/DENSEP)*G*(LI*2.0*BI)*
     +(1-LAMBDA)+((DENSEA-DENSEP)/DENSEP)*G*(LI*2.0*BI)*
     +LAMBDA)+EI*VO-EO*VI  ! ACC 2026

!      DYDX(2)=(1/GAMMA1)*((PI*G*BI**2/DENSE20)*(LAMBDA**2*VG*(DENSEA-0)+
!     +LAMBDA**2*(1-VG)*(DENSEA-DENSEW)))

      MOM0=(1/GAMMA1)*(PI*G*BI**2/DENSE20)*(LAMBDA**2*VG*(DENSEA-0)+
     +LAMBDA**2*(1-VG)*(DENSEA-DENSEW))
      MOM1=(1/GAMMA1)*(PI*G*BI**2/DENSE20)*(LAMBDA**2*VG*(DENSEA-0))
      MOM2=(1/GAMMA1)*(PI*G*BI**2/DENSE20)*(LAMBDA**2*(1-VG)*
     +(DENSEA-DENSEW))
      MOM3=+(1/GAMMA1)*EI*VO
      MOM4=-(1/GAMMA1)*EO*VI

      !PRINT *, XLOC,MOM3,MOM4
      !PRINT *, GAMMA1, EI, VO, EO, VI
      !PRINT *, GAMMA1, PI, G, DENSE20, LAMBDA
!      OPEN (UNIT=50, FILE="salida.txt", POSITION="APPEND")
!        !WRITE (UNIT=50, FMT='(7F12.6)') XLOC,MOM0,MOM1,MOM2,MOM3,MOM4,VG
!         WRITE (UNIT=50, FMT='(9F12.6)') XLOC,MOM0,BI,VG,DENSEA,
!     +DENSEW,TAMB,TPLUME,TARO
!      CLOSE (UNIT=50)


!     +LAMBDA2**2*(1-VG)*(DENSEA-DENSEW))+EI*VO-EO*VI

      DYDX(3)=EI*TARO-EO*TPLUME
      !DYDX(4)=EI*(SALARO*GAMMA/DENSE20)*DENSEA-EO*(SALPLU*GAMMA/DENSE20)*DENSEP
!     ACC 2026 (DERIVS_6): ERROR - detrainment usaba SALARO (salinidad ambiente) en lugar de
!     SALPLU (salinidad pluma). Ecuacion correcta: dF_sal/dz = EI*S_amb - EO*S_plume
!     Referencia: Wuest (1992), Dissanayake (2021). Linea erronea comentada:
!      DYDX(4)=EI*SALARO-EO*SALARO
      DYDX(4)=EI*SALARO-EO*SALPLU  ! ACC 2026
      DYDX(5)=EI*DOAMB/32.-EO*COMGP/32.+4.0*PI*RB**2*N/(VI+VB)*KOLO*
     +(HO2*PO-COMGP/32.)
      DYDX(6)=EI*DNAMB/28.-EO*CNMGP/28.+4.0*PI*RB**2*N/(VI+VB)*KOLN*
     +(HN2*PN-CNMGP/28.)
      DYDX(7)=-4.0*PI*RB**2*N/(VI+VB)*KOLO*(HO2*PO-COMGP/32.)
      DYDX(8)=-4.0*PI*RB**2*N/(VI+VB)*KOLN*(HN2*PN-CNMGP/28.)

      !PRINT *, 'DYDX', DYDX(1),DYDX(2),DYDX(3),DYDX(4),DYDX(5),DYDX(6)

      RETURN
      END
C
C------------------------------------------------------------------------------
C







C     
C----------------------------------------------------------------------
C
      SUBROUTINE RK4_5(EA,EO,EI,VI,VO,DENSEA,DENSEP,G,PI,BO,BI,
     +             TARO,TIN,TPLUME,SALAMB,SALIN,SALPLU,GAMMA,
     +             DENSE20,DENSEINNER,GAMMA1,DOAMB,COMGIN,COMGP,
     +             DNAMB,CNMGP,Y,DYDX,NN,X,H,YOUT)

      INTEGER I,NN,NMAX
      PARAMETER (NMAX=50)
      REAL*8 EA,EO,EI,DENSEA,DENSEP,G,PI,BO,BI,TARO,TIN,TPLUME,
     +SALAMB,SALIN,SALPLU,GAMMA,DENSE20,DENSEINNER,DOAMB,COMGP,
     +DNAMB,CNMGP,H,X,DYDX(NN),Y(NN),YOUT(NN),H6,HH,XH,
     +DYM(NMAX),DYT(NMAX),YT(NMAX),GAMMA1,VI,VO,COMGIN

      EXTERNAL DERIVS_5
      HH=H*0.5
      H6=H/6.
      XH=X+HH
      DO 11 I=1,NN
          YT(I)=Y(I)+HH*DYDX(I)
   11 CONTINUE
      CALL DERIVS_5(EA,EO,EI,VI,VO,DENSEA,DENSEP,G,PI,BO,BI,TARO,TIN,
     +             TPLUME,SALAMB,SALIN,SALPLU,GAMMA,DENSE20,DENSEINNER,
     +             GAMMA1,DOAMB,COMGIN,COMGP,DNAMB,CNMGP,XH,YT,DYT)
      DO 12 I=1,NN
          YT(I)=Y(I)+HH*DYT(I)
   12 CONTINUE
      CALL DERIVS_5(EA,EO,EI,VI,VO,DENSEA,DENSEP,G,PI,BO,BI,TARO,TIN,
     +             TPLUME,SALAMB,SALIN,SALPLU,GAMMA,DENSE20,DENSEINNER,
     +             GAMMA1,DOAMB,COMGIN,COMGP,DNAMB,CNMGP,XH,YT,DYM)
      DO 13 I=1,NN
          YT(I)=Y(I)+H*DYM(I)
          DYM(I)=DYT(I)+DYM(I)
   13 CONTINUE
      CALL DERIVS_5(EA,EO,EI,VI,VO,DENSEA,DENSEP,G,PI,BO,BI,TARO,TIN,
     +             TPLUME,SALAMB,SALIN,SALPLU,GAMMA,DENSE20,DENSEINNER,
     +             GAMMA1,DOAMB,COMGIN,COMGP,DNAMB,CNMGP,X+H,YT,DYT)
      DO 14 I=1,NN
          YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
   14 CONTINUE
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE DERIVS_5(EA,EO,EI,VI,VO,DENSEA,DENSEP,G,PI,BO,BI,
     +             TARO,TIN,TPLUME,SALAMB,SALIN,SALPLU,GAMMA,DENSE20,
     +             DENSEINNER,GAMMA1,DOAMB,COMGIN,COMGP,DNAMB,CNMGP,
     +             X,Y,DYDX)         
      
      REAL*8 EA,EO,EI,DENSEA,DENSEP,G,PI,BO,BI,TARO,TIN,TPLUME,SALAMB,
     +SALIN,SALPLU,GAMMA,DENSE20,DENSEINNER,DOAMB,COMGP,GAMMA1,
     +DNAMB,CNMGP,X,Y(5),DYDX(5),VI,VO,COMGIN,MOM1,MOM2,MOM3

      DYDX(1)=EA+EO-EI   
	  DYDX(2)=((1/GAMMA1)*(-PI*G*(BO**2-BI**2)*
     +((DENSEP-DENSEA)/DENSE20))-EI*VO+EO*VI)
      DYDX(3)=EA*TARO+EO*TIN-EI*TPLUME
      DYDX(4)=EA*SALAMB+EO*SALIN-EI*SALPLU
      DYDX(5)=EA*DOAMB/32.+EO*COMGIN/32.-EI*COMGP/32.
 !     DYDX(6)=EA*DNAMB/28.+EO*???/28.-EI*CNMGP/28.
      RETURN
      END 


C     
C     
C----------------------------------------------------------------------
C
      SUBROUTINE RK4_4(EI,EO,DENSEA,DENSEW,DENSEP,G,BI,LAMBDA,TARO,VG,
     +              SALARO,GAMMA,DENSE20,DOAMB,PI,RB,N,VI,VO,VB,KOLO,
     +              HO2,PO,GAMMA1,TPLUME,SALPLU,COMGP,DNAMB,KOLN,HN2,PN,
     +              CNMGP,Y,DYDX,NN,X,H,YOUT,XLOC,TAMB)

      INTEGER I,NN,NMAX
      PARAMETER (NMAX=50)
      REAL*8 EI,EO,DENSEA,DENSEW,DENSEP,G,BI,LAMBDA,TARO,VG,SALARO,
     +GAMMA,DENSE20,DOAMB,PI,RB,N,VI,VO,VB,KOLO,HO2,PO,COMGP,DNAMB,
     +KOLN,HN2,PN,TPLUME,CNMGP,H,X,DYDX(NN),Y(NN),YOUT(NN),H6,HH,XH,
     +DYM(NMAX),DYT(NMAX),YT(NMAX),GAMMA1,SALPLU,XLOC,TAMB
      EXTERNAL DERIVS_4
      HH=H*0.5
      H6=H/6.
      XH=X+HH
      DO 11 I=1,NN
          YT(I)=Y(I)+HH*DYDX(I)
   11 CONTINUE
      CALL DERIVS_4(EI,EO,DENSEA,DENSEW,DENSEP,G,BI,LAMBDA,TARO,VG,
     +           SALARO,GAMMA,DENSE20,DOAMB,PI,RB,N,VI,VO,VB,KOLO,
     +           HO2,PO,GAMMA1,TPLUME,SALPLU,COMGP,DNAMB,KOLN,HN2,
     +           PN,CNMGP,XH,YT,DYT,XLOC,TAMB)

      DO 12 I=1,NN
          YT(I)=Y(I)+HH*DYT(I)
   12 CONTINUE
      CALL DERIVS_4(EI,EO,DENSEA,DENSEW,DENSEP,G,BI,LAMBDA,TARO,VG,
     +            SALARO,GAMMA,DENSE20,DOAMB,PI,RB,N,VI,VO,VB,KOLO,
     +            HO2,PO,GAMMA1,TPLUME,SALPLU,COMGP,DNAMB,KOLN,HN2,
     +            PN,CNMGP,XH,YT,DYM,XLOC,TAMB)
      DO 13 I=1,NN
          YT(I)=Y(I)+H*DYM(I)
          DYM(I)=DYT(I)+DYM(I)
   13 CONTINUE
      CALL DERIVS_4(EI,EO,DENSEA,DENSEW,DENSEP,G,BI,LAMBDA,TARO,VG,
     +            SALARO,GAMMA,DENSE20,DOAMB,PI,RB,N,VI,VO,VB,KOLO,
     +            HO2,PO,GAMMA1,TPLUME,SALPLU,COMGP,DNAMB,KOLN,HN2,
     +            PN,CNMGP,X+H,YT,DYT,XLOC,TAMB)
      DO 14 I=1,NN
          YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
   14 CONTINUE
      RETURN
      END
C
C
C----------------------------------------------------------------------
C
      SUBROUTINE DERIVS_4(EI,EO,DENSEA,DENSEW,DENSEP,G,BI,LAMBDA,TARO,
     +              VG,SALARO,GAMMA,DENSE20,DOAMB,PI,RB,N,VI,VO,VB,
     +              KOLO,HO2,PO,GAMMA1,TPLUME,SALPLU,COMGP,DNAMB,KOLN,
     +              HN2,PN,CNMGP,X,Y,DYDX,XLOC,TAMB)         
      REAL*8 EI,EO,DENSEA,DENSEW,DENSEP,G,L,BI,LAMBDA,TARO,SALARO,GAMMA,
     +DENSE20,DOAMB,PI,RB,N,VI,VO,VB,KOLO,HO2,PO,COMGP,SALPLU,DNAMB,
     +KOLN,HN2,PN,TPLUME,CNMGP,X,Y(8),DYDX(8),VG,GAMMA1,MOM1,MOM2,
     +MOM3,MOM4,MOM0,XLOC,TAMB

C     Right-hand side of differential equations for Runge-Kutta solution
      DYDX(1)=EI-EO      
!      PRINT *, EI,EO  
C      DYDX(2)=(DENSEA-DENSEW)/DENSEP*G*(PI*BI**2)+(DENSEW-DENSEP)/
C     +DENSEP*G*(PI*(LAMBDA*BI)**2)

!      Hasta ahora he utilizado el de abajo, pero creo k esta mal 16-5-13
!      DYDX(2)=((DENSEA-DENSEW)/DENSEP)*G*PI*BI**2*(1-LAMBDA**2)+
!     +(DENSEW-DENSEP)/DENSEP*G*(PI*(LAMBDA*BI)**2)

!     Comporbar si cambian los resultados  16-5-13
!      DYDX(2)=((DENSEA-DENSEW)/DENSEP)*G*PI*BI**2*(1-LAMBDA**2)+
!     +(DENSEA-DENSEP)/DENSEP*G*(PI*(LAMBDA*BI)**2)


 !     DYDX(2)=(1/GAMMA1)*((PI*G*BI**2/DENSE20)*(LAMBDA**2*VG*(DENSEA-0)+
 !    +LAMBDA**2*(1-VG)*(DENSEA-DENSEW))+EI*VO-EO*VI)

      DYDX(2)=(1/GAMMA1)*(((DENSEA-DENSEW)/DENSEP)*G*PI*BI**2*
     +(1-LAMBDA**2)+((DENSEA-DENSEP)/DENSEP)*G*PI*BI**2*
     +LAMBDA**2)+EI*VO-EO*VI


!      DYDX(2)=(1/GAMMA1)*((PI*G*BI**2/DENSE20)*(LAMBDA**2*VG*(DENSEA-0)+
!     +LAMBDA**2*(1-VG)*(DENSEA-DENSEW)))

      MOM0=(1/GAMMA1)*(PI*G*BI**2/DENSE20)*(LAMBDA**2*VG*(DENSEA-0)+
     +LAMBDA**2*(1-VG)*(DENSEA-DENSEW))
      MOM1=(1/GAMMA1)*(PI*G*BI**2/DENSE20)*(LAMBDA**2*VG*(DENSEA-0))
      MOM2=(1/GAMMA1)*(PI*G*BI**2/DENSE20)*(LAMBDA**2*(1-VG)*
     +(DENSEA-DENSEW))
      MOM3=+(1/GAMMA1)*EI*VO
      MOM4=-(1/GAMMA1)*EO*VI

      !PRINT *, XLOC,MOM3,MOM4
      !PRINT *, GAMMA1, EI, VO, EO, VI
      !PRINT *, GAMMA1, PI, G, DENSE20, LAMBDA
!      OPEN (UNIT=50, FILE="salida.txt", POSITION="APPEND")
!        !WRITE (UNIT=50, FMT='(7F12.6)') XLOC,MOM0,MOM1,MOM2,MOM3,MOM4,VG
!         WRITE (UNIT=50, FMT='(9F12.6)') XLOC,MOM0,BI,VG,DENSEA,
!     +DENSEW,TAMB,TPLUME,TARO
!      CLOSE (UNIT=50)


!     +LAMBDA2**2*(1-VG)*(DENSEA-DENSEW))+EI*VO-EO*VI

      DYDX(3)=EI*TARO-EO*TPLUME
      !DYDX(4)=EI*(SALARO*GAMMA/DENSE20)*DENSEA-EO*(SALPLU*GAMMA/DENSE20)*DENSEP
!     ACC 2026 (DERIVS_4): ERROR - detrainment usaba SALARO (salinidad ambiente) en lugar de
!     SALPLU (salinidad pluma). Ecuacion correcta: dF_sal/dz = EI*S_amb - EO*S_plume
!     Referencia: Wuest (1992). Linea erronea comentada:
!      DYDX(4)=EI*SALARO-EO*SALARO
      DYDX(4)=EI*SALARO-EO*SALPLU  ! ACC 2026
      DYDX(5)=EI*DOAMB/32.-EO*COMGP/32.+4.0*PI*RB**2*N/(VI+VB)*KOLO*
     +(HO2*PO-COMGP/32.)
      DYDX(6)=EI*DNAMB/28.-EO*CNMGP/28.+4.0*PI*RB**2*N/(VI+VB)*KOLN*
     +(HN2*PN-CNMGP/28.)
      DYDX(7)=-4.0*PI*RB**2*N/(VI+VB)*KOLO*(HO2*PO-COMGP/32.)
      DYDX(8)=-4.0*PI*RB**2*N/(VI+VB)*KOLN*(HN2*PN-CNMGP/28.)
      RETURN
      END
C
C------------------------------------------------------------------------------
C
C     
C----------------------------------------------------------------------
C
      SUBROUTINE RK4_3(E,DENSEA,DENSEW,DENSEP,G,B,LAMBDA,TAMB,
     +              SALAMB,GAMMA,DENSE20,DOAMB,PI,RB,N,V,VB,KOLO,HO2,
     +              PO,COMGP,DNAMB,KOLN,HN2,PN,CNMGP,Y,DYDX,NN,X,H,YOUT)

      INTEGER I,NN,NMAX
      PARAMETER (NMAX=50)
      REAL*8 E,DENSEA,DENSEW,DENSEP,G,B,LAMBDA,TAMB,SALAMB,GAMMA,
     +DENSE20,DOAMB,PI,RB,N,V,VB,KOLO,HO2,PO,COMGP,DNAMB,KOLN,HN2,PN,
     +CNMGP,H,X,DYDX(NN),Y(NN),YOUT(NN),H6,HH,XH,DYM(NMAX),DYT(NMAX),
     +YT(NMAX)
      EXTERNAL DERIVS_3
      HH=H*0.5
      H6=H/6.
      XH=X+HH
      DO 11 I=1,NN
          YT(I)=Y(I)+HH*DYDX(I)
   11 CONTINUE
      CALL DERIVS_3(E,DENSEA,DENSEW,DENSEP,G,B,LAMBDA,TAMB,SALAMB,
     +           GAMMA,DENSE20,DOAMB,PI,RB,N,V,VB,KOLO,HO2,PO,
     +           COMGP,DNAMB,KOLN,HN2,PN,CNMGP,XH,YT,DYT)
      DO 12 I=1,NN
          YT(I)=Y(I)+HH*DYT(I)
   12 CONTINUE
      CALL DERIVS_3(E,DENSEA,DENSEW,DENSEP,G,B,LAMBDA,TAMB,SALAMB,
     +            GAMMA,DENSE20,DOAMB,PI,RB,N,V,VB,KOLO,HO2,PO,
     +            COMGP,DNAMB,KOLN,HN2,PN,CNMGP,XH,YT,DYM)
      DO 13 I=1,NN
          YT(I)=Y(I)+H*DYM(I)
          DYM(I)=DYT(I)+DYM(I)
   13 CONTINUE
      CALL DERIVS_3(E,DENSEA,DENSEW,DENSEP,G,B,LAMBDA,TAMB,SALAMB,
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
      SUBROUTINE DERIVS_3(E,DENSEA,DENSEW,DENSEP,G,B,LAMBDA,TAMB,
     +              SALAMB,GAMMA,DENSE20,DOAMB,PI,RB,N,V,VB,KOLO,HO2,PO,
     +              COMGP,DNAMB,KOLN,HN2,PN,CNMGP,X,Y,DYDX)         
      REAL*8 E,DENSEA,DENSEW,DENSEP,G,L,B,LAMBDA,TAMB,SALAMB,GAMMA,
     +DENSE20,DOAMB,PI,RB,N,V,VB,KOLO,HO2,PO,COMGP,DNAMB,KOLN,HN2,PN,
     +CNMGP,X,Y(8),DYDX(8)
C     Right-hand side of differential equations for Runge-Kutta solution
      DYDX(1)=E         
C      DYDX(2)=(DENSEA-DENSEW)/DENSEP*G*(PI*B**2)+(DENSEW-DENSEP)/
C     +DENSEP*G*(PI*(LAMBDA*B)**2)

!      Hasta ahora he utilizado el de abajo, pero creo k esta mal 16-5-13
!      DYDX(2)=((DENSEA-DENSEW)/DENSEP)*G*PI*B**2*(1-LAMBDA**2)+
!     +(DENSEW-DENSEP)/DENSEP*G*(PI*(LAMBDA*B)**2)

!     Comporbar si cambian los resultados  16-5-13
      DYDX(2)=((DENSEA-DENSEW)/DENSEP)*G*PI*B**2*(1-LAMBDA**2)+
     +(DENSEA-DENSEP)/DENSEP*G*(PI*(LAMBDA*B)**2)
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
