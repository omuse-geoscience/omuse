C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 1999 by UCAR                   *
C     *                                                               *
C     *       UNIVERSITY CORPORATION for ATMOSPHERIC RESEARCH         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                      FISHPACK version 4.1                     *
C     *                                                               *
C     *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *
C     *                                                               *
C     *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *
C     *                                                               *
C     *                             BY                                *
C     *                                                               *
C     *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *
C     *                                                               *
C     *                             OF                                *
C     *                                                               *
C     *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *
C     *                                                               *
C     *                BOULDER, COLORADO  (80307)  U.S.A.             *
C     *                                                               *
C     *                   WHICH IS SPONSORED BY                       *
C     *                                                               *
C     *              THE NATIONAL SCIENCE FOUNDATION                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C
C     PROGRAM TO ILLUSTRATE THE USE OF HWSCSP
C
      DIMENSION       F(48,33)   ,BDTF(33)   ,W(775)     ,R(33)      ,
     1                THETA(48)
C
      PI = PIMACH(DUM)
      INTL = 0
      TS = 0.
      TF = PI/2.
      M = 36
      MBDCND = 6
      RS = 0.
      RF = 1.
      N = 32
      NBDCND = 5
      ELMBDA = 0.
      IDIMF = 48
C
C     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING THE
C     BOUNDARY DATA AND THE RIGHT SIDE OF THE EQUATION.
C
      MP1 = M+1
      DTHETA = TF/FLOAT(M)
      DO 101 I=1,MP1
         THETA(I) = FLOAT(I-1)*DTHETA
  101 CONTINUE
      NP1 = N+1
      DR = 1./FLOAT(N)
      DO 102 J=1,NP1
         R(J) = FLOAT(J-1)*DR
  102 CONTINUE
C
C     GENERATE NORMAL DERIVATIVE DATA AT EQUATOR
C
      DO 103 J=1,NP1
         BDTF(J) = 0.
  103 CONTINUE
C
C     COMPUTE BOUNDARY DATA ON THE SURFACE OF THE SPHERE
C
      DO 104 I=1,MP1
         F(I,N+1) = COS(THETA(I))**4
  104 CONTINUE
C
C     COMPUTE RIGHT SIDE OF EQUATION
C
      DO 106 I=1,MP1
         CI4 = 12.*COS(THETA(I))**2
         DO 105 J=1,N
            F(I,J) = CI4*R(J)**2
  105    CONTINUE
  106 CONTINUE
C
      CALL HWSCSP (INTL,TS,TF,M,MBDCND,BDTS,BDTF,RS,RF,N,NBDCND,BDRS,
     1             BDRF,ELMBDA,F,IDIMF,PERTRB,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR
C
      ERR = 0.
      DO 108 I=1,MP1
         CI4 = COS(THETA(I))**4
         DO 107 J=1,N
            Z = ABS(F(I,J)-CI4*R(J)**4)
            IF (Z .GT. ERR) ERR = Z
  107    CONTINUE
  108 CONTINUE
      IW = INT(W(1))
      PRINT 1001 , IERROR,ERR,IW
C
C     THE FOLLOWING PROGRAM ILLUSTRATES THE USE OF HWSCSP TO SOLVE
C     A THREE DIMENSIONAL PROBLEM WHICH HAS LONGITUDNAL DEPENDENCE
C
      MBDCND = 2
      NBDCND = 1
      DPHI = PI/72.
      ELMBDA = -2.*(1.-COS(DPHI))/DPHI**2
C
C     COMPUTE BOUNDARY DATA ON THE SURFACE OF THE SPHERE
C
      DO 109 I=1,MP1
         F(I,N+1) = SIN(THETA(I))
  109 CONTINUE
C
C     COMPUTE RIGHT SIDE OF THE EQUATION
C
      DO 111 J=1,N
         DO 110 I=1,MP1
            F(I,J) = 0.
  110    CONTINUE
  111 CONTINUE
C
      CALL HWSCSP (INTL,TS,TF,M,MBDCND,BDTS,BDTF,RS,RF,N,NBDCND,BDRS,
     1             BDRF,ELMBDA,F,IDIMF,PERTRB,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR   (FOURIER COEFFICIENTS)
C
      ERR = 0
      DO 113 I=1,MP1
         SI = SIN(THETA(I))
         DO 112 J=1,NP1
            Z = ABS(F(I,J)-R(J)*SI)
            IF (Z .GT. ERR) ERR = Z
  112    CONTINUE
  113 CONTINUE
C
      IW = INT(W(1))
      PRINT 1002 , IERROR,ERR,IW
      STOP
C
 1001 FORMAT (1H1,20X,27HSUBROUTINE HWSCSP EXAMPLE 1///
     1        10X,46HTHE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS//
     2        32X,10HIERROR = 0/
     3        18X,34HDISCRETIZATION ERROR = 7.99842E-04/
     4        12X,32HREQUIRED LENGTH OF W ARRAY = 775//
     5        10X,32HTHE OUTPUT FROM YOUR COMPUTER IS//
     6        32X,8HIERROR =,I2/18X,22HDISCRETIZATION ERROR =,E12.5/
     7        12X,28HREQUIRED LENGTH OF W ARRAY =,I4)
 1002 FORMAT (1H1,20X,27HSUBROUTINE HWSCSP EXAMPLE 2///
     1        10X,46HTHE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS//
     2        32X,10HIERROR = 0/
     3        18X,34HDISCRETIZATION ERROR = 5.86824E-05/
     4        12X,32HREQUIRED LENGTH OF W ARRAY = 775//
     5        10X,32HTHE OUTPUT FROM YOUR COMPUTER IS//
     6        32X,8HIERROR =,I2/18X,22HDISCRETIZATION ERROR =,E12.5/
     7        12X,28HREQUIRED LENGTH OF W ARRAY =,I4)
C
      END
