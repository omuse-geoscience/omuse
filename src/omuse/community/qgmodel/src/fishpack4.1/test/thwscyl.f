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
C          PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE HWSCYL TO SOLVE
C     THE EQUATION
C
C     (1/R)(D/DR)(R*(DU/DR)) + (D/DZ)(DU/DZ)
C
C     = (2*R*Z)**2*(4*Z**2 + 3*R**2)
C
C     ON THE RECTANGLE 0 .LT. R .LT. 1, 0 .LT. Z .LT. 1 WITH THE
C     BOUNDARY CONDITIONS
C
C     U(0,Z) UNSPECIFIED
C                                            0 .LE. Z .LE. 1
C     (DU/DR)(1,Z) = 4*Z**4
C
C     AND
C
C     (DU/DZ)(R,0) = 0
C                                            0 .LE. R .LE. 1
C     (DU/DZ)(R,1) = 4*R**4 .
C
C          THE R-INTERVAL WILL BE DIVIDED INTO 50 PANELS AND THE
C     Z-INTERVAL WILL BE DIVIDED INTO 100 PANELS.
C
      DIMENSION       F(75,105)  ,BDA(101)   ,BDB(101)   ,BDC(51)    ,
     1                BDD(51)    ,W(1373)    ,R(51)      ,Z(101)
C
C     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF.  ALSO NOTE THAT W
C     IS DIMENSIONED 4*(N+1) + (13 + INT(LOG2(N+1)))*(M+1) .
C
      IDIMF = 75
      A = 0.
      B = 1.
      M = 50
      MBDCND = 6
      C = 0.
      D = 1.
      N = 100
      NBDCND = 3
      ELMBDA = 0.
C
C     AUXILIARY QUANTITIES.
C
      MP1 = M+1
      NP1 = N+1
C
C     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING
C     BOUNDARY DATA AND THE RIGHT SIDE OF THE POISSON EQUATION.
C
      DO 101 I=1,MP1
         R(I) = FLOAT(I-1)/50.
  101 CONTINUE
      DO 102 J=1,NP1
         Z(J) = FLOAT(J-1)/100.
  102 CONTINUE
C
C     GENERATE BOUNDARY DATA.
C
      DO 103 J=1,NP1
         BDB(J) = 4.*Z(J)**4
  103 CONTINUE
      DO 104 I=1,MP1
         BDC(I) = 0.
         BDD(I) = 4.*R(I)**4
  104 CONTINUE
C
C     BDA IS A DUMMY VARIABLE.
C
C
C     GENERATE RIGHT SIDE OF EQUATION.
C
      DO 106 I=1,MP1
         DO 105 J=1,NP1
            F(I,J) = 4.*R(I)**2*Z(J)**2*(4.*Z(J)**2+3.*R(I)**2)
  105    CONTINUE
  106 CONTINUE
      CALL HWSCYL (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,ELMBDA,F,
     1             IDIMF,PERTRB,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR BY MINIMIZING OVER ALL A THE FUNCTION
C     NORM(F(I,J) - A*1 - U(R(I),Z(J))).  THE EXACT SOLUTION IS
C                U(R,Z) = (R*Z)**4 + ARBITRARY CONSTANT.
C
      X = 0.
      DO 108 I=1,MP1
         DO 107 J=1,NP1
            X = X+F(I,J)-(R(I)*Z(J))**4
  107    CONTINUE
  108 CONTINUE
      X = X/FLOAT(NP1*MP1)
      DO 110 I=1,MP1
         DO 109 J=1,NP1
            F(I,J) = F(I,J)-X
  109    CONTINUE
  110 CONTINUE
      ERR = 0.
      DO 112 I=1,MP1
         DO 111 J=1,NP1
            X = ABS(F(I,J)-(R(I)*Z(J))**4)
            IF (X .GT. ERR) ERR = X
  111    CONTINUE
  112 CONTINUE
      PRINT 1001 , IERROR,PERTRB,ERR,W(1)
      STOP
C
 1001 FORMAT (1H1,20X,25HSUBROUTINE HWSCYL EXAMPLE///
     1        10X,46HTHE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS//
     2        32X,10HIERROR = 0/32X,20HPERTRB = 2.26734E-04/
     3        18X,34HDISCRETIZATION ERROR = 3.73672E-04/
     4        12X,33HREQUIRED LENGTH OF W ARRAY = 1118//
     5        10X,32HTHE OUTPUT FROM YOUR COMPUTER IS//
     6        32X,8HIERROR =,I2/32X,8HPERTRB =,E12.5/
     7        18X,22HDISCRETIZATION ERROR =,E12.5/
     8        12X,28HREQUIRED LENGTH OF W ARRAY =,F5.0)
C
      END
