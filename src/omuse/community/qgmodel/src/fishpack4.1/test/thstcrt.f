C
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
C     PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE HSTCRT TO SOLVE
C     THE EQUATION
C
C    (D/DX)(DU/DX) + (D/DY)(DU/DY) - 2*U = -2(PI**2+1)SIN(PI*X)COS(PI*Y)
C
C     WHERE  1 .LE. X .LE. 3 AND -1 .LE. Y .LE. 1 AND THE BOUNDARY
C     CONDITIONS ARE
C
C     U = 0 ON X = 1,  DU/DX = -PI*COS(PI*Y) ON X = 3
C
C     AND U IS PERIODIC IN Y .
C
C     WE WANT TO HAVE 48 UNKNOWNS IN THE X-INTERVAL AND 53 UNKNOWNS
C     IN THE Y-INTERVAL.
C
      DIMENSION       F(50,53)   ,BDA(53)    ,BDB(53)    ,W(1076)    ,
     1                X(48)      ,Y(53)
C
C     FROM THE DIMENSION STATEMENT WE GET IDIMF = 50.  ALSO NOTE THAT
C     W IS DIMENSIONED (13 + INT(LOG2(N))*M + 4*N.
C
      IDIMF = 50
      A = 1.
      B = 3.
      M = 48
      DX = (B-A)/FLOAT(M)
      MBDCND = 2
      C = -1.
      D = 1.
      N = 53
      DY = (D-C)/FLOAT(N)
      NBDCND = 0
      ELMBDA = -2.
C
C     AUXILIARY QUANTITIES
C
      PI = PIMACH(DUM)
      PISQ = PI*PI
C
C     GENERATE AND STORE GRID POINTS FOR COMPUTATION OF BOUNDARY DATA
C     AND THE RIGHT SIDE OF THE HELMHOLTZ EQUATION.
C
      DO 101 I=1,M
         X(I) = A+(FLOAT(I)-0.5)*DX
  101 CONTINUE
      DO 102 J=1,N
         Y(J) = C+(FLOAT(J)-0.5)*DY
  102 CONTINUE
C
C     GENERATE BOUNDARY DATA.
C
      DO 103 J=1,N
         BDA(J) = 0.
         BDB(J) = -PI*COS(PI*Y(J))
  103 CONTINUE
C
C     BDC AND BDD ARE DUMMY ARGUMENTS IN THIS EXAMPLE.
C
C     GENERATE RIGHT SIDE OF EQUATION.
C
      T = -2.*(PISQ+1.)
      DO 105 I=1,M
         DO 104 J=1,N
            F(I,J) = T*SIN(PI*X(I))*COS(PI*Y(J))
  104    CONTINUE
  105 CONTINUE
      CALL HSTCRT (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,ELMBDA,F,
     1             IDIMF,PERTRB,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
C
C               U(X,Y) = SIN(PI*X)*COS(PI*Y) .
C
      ERR = 0.
      DO 107 I=1,M
         DO 106 J=1,N
            T = ABS(F(I,J)-SIN(PI*X(I))*COS(PI*Y(J)))
            IF (T .GT. ERR) ERR = T
  106    CONTINUE
  107 CONTINUE
      PRINT 1001 , IERROR,ERR,W(1)
      STOP
C
 1001 FORMAT (1H1,20X,25HSUBROUTINE HSTCRT EXAMPLE///
     1        10X,46HTHE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS//
     2        32X,10HIERROR = 0/
     3        18X,34HDISCRETIZATION ERROR = 1.26001E-03/
     4        12X,32HREQUIRED LENGTH OF W ARRAY = 884//
     5        10X,32HTHE OUTPUT FROM YOUR COMPUTER IS//
     6        32X,8HIERROR =,I2/18X,22HDISCRETIZATION ERROR =,E12.5/
     7        12X,28HREQUIRED LENGTH OF W ARRAY =,F4.0)
C
      END
