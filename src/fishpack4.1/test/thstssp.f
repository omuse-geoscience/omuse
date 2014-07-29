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
C     PROGRAM TO ILLUSTRATE THE USE OF HSTSSP TO SOLVE POISSON"S
C     EQUATION
C
C             (1/SIN(THETA))(D/DTHETA)(SIN(THETA)*DU/DTHETA) +
C
C    (1/SIN(THETA)**2)(D/DPHI)(DU/DPHI) = 2 - 6*(SIN(THETA)*SIN(PHI))**2
C     ON THE NORTHERN HEMISPHERE SUBJECT TO EQUATORIAL SYMMETRY, I.E.
C     THE DERIVATIVE OF THE SOLUTION AT THETA = PI/2 IS ZERO.  A
C     5-DEGREE GRID IS TO BE USED.
C
C     THE EXACT SOLUTION IS NOT UNIQUE.  ANY FUNCTION OF THE FORM
C
C           U(THETA,PHI) = (SIN(THETA)*SIN(PHI))**2 + CONSTANT
C
C     IS A SOLUTION.
C
      DIMENSION       F(18,72)   ,BDB(72)    ,SINT(18)   ,SINP(72)   ,
     1                W(630)
C
C     THE VALUE OF IDIMF IS THE FIRST DIMENSION OF F. W IS DIMENSIONED
C     (13 + INT(LOG2(N)))*M + 4*N
C
      PI = PIMACH(DUM)
      A = 0.
      B = PI/2.
      M = 18
      MBDCND = 6
      C = 0.
      D = 2.*PI
      N = 72
      NBDCND = 0
      ELMBDA = 0.
      IDIMF = 18
C
C     GENERATE SINES FOR USE IN SUBSEQUENT COMPUTATIONS
C
      DTHETA = B/FLOAT(M)
      DO 101 I=1,M
         SINT(I) = SIN((FLOAT(I)-0.5)*DTHETA)
  101 CONTINUE
      DPHI = D/FLOAT(N)
      DO 102 J=1,N
         SINP(J) = SIN((FLOAT(J)-0.5)*DPHI)
  102 CONTINUE
C
C     COMPUTE RIGHT SIDE OF EQUATION AND STORE IN F
C
      DO 104 J=1,N
         DO 103 I=1,M
            F(I,J) = 2.-6.*(SINT(I)*SINP(J))**2
  103    CONTINUE
  104 CONTINUE
C
C     STORE DERIVATIVE DATA AT THE EQUATOR
C
      DO 105 J=1,N
         BDB(J) = 0.
  105 CONTINUE
C
C     BDA, BDC, AND BDD ARE DUMMY VARIABLES.
C
      CALL HSTSSP (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,ELMBDA,F,
     1             IDIMF,PERTRB,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR. SINCE PROBLEM IS SINGULAR, THE
C     SOLUTION MUST BE NORMALIZED.
C
      ERR = 0.
      DO 107 J=1,N
         DO 106 I=1,M
            Z = ABS(F(I,J)-(SINT(I)*SINP(J))**2-F(1,1))
            IF (Z .GT. ERR) ERR = Z
  106    CONTINUE
  107 CONTINUE
C
      PRINT 1001 , IERROR,PERTRB,ERR,W(1)
      STOP
C
 1001 FORMAT (1H1,20X,25HSUBROUTINE HSTSSP EXAMPLE///
     1        10X,46HTHE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS//
     2        32X,10HIERROR = 0/32X,20HPERTRB = 6.35830E-04/
     3        18X,34HDISCRETIZATION ERROR = 3.37523E-03/
     4        12X,32HREQUIRED LENGTH OF W ARRAY = 540//
     5        10X,32HTHE OUTPUT FROM YOUR COMPUTER IS//
     6        32X,8HIERROR =,I2/32X,8HPERTRB =,E12.5/
     7        18X,22HDISCRETIZATION ERROR =,E12.5/
     8        12X,28HREQUIRED LENGTH OF W ARRAY =,F4.0)
C
      END
