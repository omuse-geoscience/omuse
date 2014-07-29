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
C          PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE HWSCRT TO SOLVE
C     THE EQUATION
C
C     (D/DX)(DU/DX) + (D/DY)(DU/DY) - 4*U
C
C     = (2 - (4 + PI**2/4)*X**2)*COS((Y+1)*PI/2)
C
C     WITH THE BOUNDARY CONDITIONS
C     ON THE RECTANGLE 0 .LT. X .LT. 2, -1 .LT. Y .LT. 3 WITH THE
C
C     U(0,Y) = 0
C                                          -1 .LE. Y .LE. 3
C     (DU/DX)(2,Y) = 4*COS((Y+1)*PI/2)
C
C     AND WITH U PERIODIC IN Y.
C          THE X-INTERVAL WILL BE DIVIDED INTO 40 PANELS AND THE
C     Y-INTERVAL WILL BE DIVIDED INTO 80 PANELS.
C
      DIMENSION       F(45,82)   ,BDB(81)    ,W(1103)    ,X(41)      ,
     1                Y(81)
C
C     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF.  ALSO NOTE THAT W
C     IS DIMENSIONED 4*(N+1) + (13 + INT(LOG2(N+1)))*(M+1) .
C
      IDIMF = 45
      A = 0.
      B = 2.
      M = 40
      MBDCND = 2
      C = -1.
      D = 3.
      N = 80
      NBDCND = 0
      ELMBDA = -4.
C
C     AUXILIARY QUANTITIES.
C
      PI = PIMACH(DUM)
      PIBY2 = PI/2.
      PISQ = PI**2
      MP1 = M+1
      NP1 = N+1
C
C     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING
C     BOUNDARY DATA AND THE RIGHT SIDE OF THE HELMHOLTZ EQUATION.
C
      DO 101 I=1,MP1
         X(I) = FLOAT(I-1)/20.
  101 CONTINUE
      DO 102 J=1,NP1
         Y(J) = -1.+FLOAT(J-1)/20.
  102 CONTINUE
C
C     GENERATE BOUNDARY DATA.
C
      DO 103 J=1,NP1
         BDB(J) = 4.*COS((Y(J)+1.)*PIBY2)
  103 CONTINUE
C
C     BDA, BDC, AND BDD ARE DUMMY VARIABLES.
C
      DO 104 J=1,NP1
         F(1,J) = 0.
  104 CONTINUE
C
C     GENERATE RIGHT SIDE OF EQUATION.
C
      DO 106 I=2,MP1
         DO 105 J=1,NP1
            F(I,J) = (2.-(4.+PISQ/4.)*X(I)**2)*COS((Y(J)+1.)*PIBY2)
  105    CONTINUE
  106 CONTINUE
      CALL HWSCRT (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,ELMBDA,F,
     1             IDIMF,PERTRB,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
C                U(X,Y) = X**2*COS((Y+1)*PIBY2)
C
      ERR = 0.
      DO 108 I=1,MP1
         DO 107 J=1,NP1
            Z = ABS(F(I,J)-X(I)**2*COS((Y(J)+1.)*PIBY2))
            IF (Z .GT. ERR) ERR = Z
  107    CONTINUE
  108 CONTINUE
      PRINT 1001 , IERROR,ERR,W(1)
      STOP
C
 1001 FORMAT (1H1,20X,25HSUBROUTINE HWSCRT EXAMPLE///
     1        10X,46HTHE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS//
     2        32X,10HIERROR = 0/
     3        18X,34HDISCRETIZATION ERROR = 5.36508E-04/
     4        12X,32HREQUIRED LENGTH OF W ARRAY = 880//
     5        10X,32HTHE OUTPUT FROM YOUR COMPUTER IS//
     6        32X,8HIERROR =,I2/18X,22HDISCRETIZATION ERROR =,E12.5/
     7        12X,28HREQUIRED LENGTH OF W ARRAY =,F4.0)
C
      END
