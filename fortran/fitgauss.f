      SUBROUTINE FITGAUSS(X,Y,NY,PAR,NPAR,ERR,ICON,SIGMA,VERBOSE)
*=====================================================================
*
* Uses LMM fitting package to perform a fit of a gaussian function
* to a set of X,Y points. The function must be in a subroutine called FUNCT
* and require NPAR parameters to be fitted.
*--
*i X........the data array
*i Y........the data array.
*i NY.......the number of points in Y.
*b PAR......the parameters array.
*i NPAR.....the total number parameters in PAR
*o ERR......the 1 sigma errors on PAR.
*
* LN Oct 1992  - `channans' software
*--
* External Subroutines: LMM
* Internal Function   : FUNCT
*---
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      EXTERNAL FUNCT1,FUNCT2,FUNCT3
*
      PARAMETER ( NPMX = 15, NPTS = 8200, IW = NPMX*(NPMX+20)/2 )
      DIMENSION X(NY),Y(NY),RES(NPTS),RPLANE(NPTS)
      DIMENSION PAR(NPAR),GRAD(NPMX),ERR(NPMX),WK(IW)
      INTEGER   IWK(NPMX),ICON(NPAR)
      CHARACTER PARAM(NPMX)*9,BELL*1
      LOGICAL VERBOSE
*-----------------------------------
      BELL = CHAR(7)
      LOUT = 6

C-- Give names to parameters, set max number of iterations and tolerance
      PARAM(1) = 'Amplitude'
      PARAM(2) = 'Mean     '
      PARAM(3) = 'Sigma    '
      PARAM(4) = PARAM(1)
      PARAM(5) = PARAM(2)
      PARAM(6) = PARAM(3)
      PARAM(7) = PARAM(1)
      PARAM(8) = PARAM(2)
      PARAM(9) = PARAM(3)
      PARAM(10) = PARAM(1)
      PARAM(11) = PARAM(2)
      PARAM(12) = PARAM(3)
      PARAM(13) = PARAM(1)
      PARAM(14) = PARAM(2)
      PARAM(15) = PARAM(3)
      NITS     = 30
      TOL      = 0.01D0

C-- Set Abscissas
*
      IF(NPAR.EQ.3)THEN
        CALL LMM(PAR,GRAD,SUMSQ,NY,NPAR,TOL,2.0D0,0.4D0,NITS,-1,
     .      LOUT,ERR,RES,RPLANE,WK,IW,IWK,ICON,
     .      FUNCT1,X,1,Y,IERR)
      ELSEIF(NPAR.EQ.6)THEN
        CALL LMM(PAR,GRAD,SUMSQ,NY,NPAR,TOL,2.0D0,0.4D0,NITS,-1,
     .      LOUT,ERR,RES,RPLANE,WK,IW,IWK,ICON,
     .      FUNCT2,X,1,Y,IERR)
      ELSEIF(NPAR.EQ.9)THEN
        CALL LMM(PAR,GRAD,SUMSQ,NY,NPAR,TOL,2.0D0,0.4D0,NITS,-1,
     .      LOUT,ERR,RES,RPLANE,WK,IW,IWK,ICON,
     .      FUNCT3,X,1,Y,IERR)
      ELSE
        WRITE(*,'('' I was expecting npar 3 6 or 9 '')')
	STOP
      ENDIF

      IF (IERR .EQ. 0 .AND. VERBOSE) WRITE(*,'(/'' Convergence OK '')')
      IF (IERR .NE. 0) 
     +    WRITE(*,'(/'' LMM failed: IERR = '',i3)') IERR
*
      DO 150 I = 1,NPAR
  150 IF (ICON(I) .EQ. 1) ERR(I) = 0.D0
*
      IF(VERBOSE)THEN
        WRITE (6,2050) NITS,SUMSQ
 2050   FORMAT (' Number of iterations   :',I3/
     .          ' Residual sum of squares:',G13.6)
        WRITE (6,2060) (I,PARAM(I),PAR(I),ERR(I),I = 1,NPAR)
 2060   FORMAT ( /' FITTED PARAMETERS FROM "LMM" FIT PROGRAM'/
     .            ' ----------------------------------------'/
     .            '    PARAMETER       VALUE        ERROR'/
     .             10(I3,') ',A9,2X,G13.6,2X,G12.5/)/)
      ENDIF
*
      IF (IERR .NE. 0) THEN
	SIGMA=9999.
	GOTO 250
      ENDIF
      SUMWS = 0.
      SUMW  = 0.
      DO 200 I = 1,NY
	 IF(NPAR.EQ.3)THEN
           CALL FUNCT1(I,X(I),1,Y(I),WEIGHT,PAR,NPAR,RES(I),GRAD,
     .     ICON,IFL)
	 ELSEIF(NPAR.EQ.6)THEN
           CALL FUNCT2(I,X(I),1,Y(I),WEIGHT,PAR,NPAR,RES(I),GRAD,
     .     ICON,IFL)
	 ELSEIF(NPAR.EQ.9)THEN
           CALL FUNCT3(I,X(I),1,Y(I),WEIGHT,PAR,NPAR,RES(I),GRAD,
     .     ICON,IFL)
	 ENDIF
         SUMWS = SUMWS + WEIGHT*RES(I)**2
         SUMW  = SUMW  + WEIGHT
  200 CONTINUE
*
      SIGMA = SQRT(SUMWS/SUMW)
      IF(VERBOSE)WRITE (6,2112) SIGMA
 2112 FORMAT(' Weighted RMS of residuals =',G13.6/)
*
  250 CONTINUE
*
      RETURN
      END
*
*---
      SUBROUTINE FUNCT1( I,X,NX,Y,W,PAR,NPAR,RES,GRAD,ICON,IFL )
*==============================================================
*  This function is to fit a gaussian to a sequence of amplitudes
*  PAR(1) = Amplitude
*  PAR(2) = Mean
*  PAR(3) = Standard Deviation
*  W      = Weight
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PAR(NPAR),GRAD(NPAR),ICON(NPAR)
      DATA SQT2PI / 2.5066282746310D0 /
      AMP1=0D0

      X1 = -0.5*((X-PAR(2))/PAR(3))**2
      IF(X1.GT.-20D0) AMP1 = PAR(1)/SQT2PI * DEXP(X1)

      AMP=AMP1
      RES    = Y - AMP
      W=1d0

      RETURN
      END
*---
      SUBROUTINE FUNCT2( I,X,NX,Y,W,PAR,NPAR,RES,GRAD,ICON,IFL )
*==============================================================
*  This function is to fit a gaussian to a sequence of amplitudes
*  PAR(1) = Amplitude
*  PAR(2) = Mean
*  PAR(3) = Standard Deviation
*  W      = Weight
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PAR(NPAR),GRAD(NPAR),ICON(NPAR)
      DATA SQT2PI / 2.5066282746310D0 /
      AMP1=0D0
      AMP2=0D0

      X1 = -0.5*((X-PAR(2))/PAR(3))**2
      IF(X1.GT.-20D0) AMP1 = PAR(1)/SQT2PI * DEXP(X1)

      IF (PAR(6).ne.0D0)THEN
        X2 = -0.5*((X-PAR(5))/PAR(6))**2
        IF(X2.GT.-20D0) AMP2 = PAR(4)/SQT2PI * DEXP(X2)
      ENDIF

      AMP=AMP1+AMP2
      RES    = Y - AMP
      W=1d0

      RETURN
      END
*---
      SUBROUTINE FUNCT3( I,X,NX,Y,W,PAR,NPAR,RES,GRAD,ICON,IFL )
*==============================================================
*  This function is to fit a gaussian to a sequence of amplitudes
*  PAR(1) = Amplitude
*  PAR(2) = Mean
*  PAR(3) = Standard Deviation
*  W      = Weight
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PAR(NPAR),GRAD(NPAR),ICON(NPAR)
      DATA SQT2PI / 2.5066282746310D0 /
      AMP1=0D0
      AMP2=0D0
      AMP3=0D0

      X1 = -0.5*((X-PAR(2))/PAR(3))**2
      IF(X1.GT.-20D0) AMP1 = PAR(1)/SQT2PI * DEXP(X1)

      IF (PAR(6).ne.0D0)THEN
          X2 = -0.5*((X-PAR(5))/PAR(6))**2
          IF(X2.GT.-20D0) AMP2 = PAR(4)/SQT2PI * DEXP(X2)
      ENDIF	

      IF(PAR(9).ne.0D0)THEN
        X3 = -0.5*((X-PAR(8))/PAR(9))**2
        IF(X3.GT.-20D0) AMP3 = PAR(7)/SQT2PI * DEXP(X3)
      ENDIF

      AMP=AMP1+AMP2+AMP3
      RES    = Y - AMP
      W=1d0

      RETURN
      END
