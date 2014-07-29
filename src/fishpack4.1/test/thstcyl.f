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
C     PROGRAM TO ILLUSTRATE THE USE OF HSTCYL TO SOLVE THE EQUATION
C
C    (1/R)(D/DR)(R*DU/DR) + (D/DZ)(DU/DZ) = (2*R*Z)**2*(4*Z**2 + 3*R**2)
C
C     ON THE RECTANGLE 0 .LT. R .LT. 1 , 0 .LT. Z .LT. 1 WITH THE
C     BOUNDARY CONDITIONS
C
C     (DU/DR)(1,Z) = 4*Z**2  FOR  0 .LE. Z .LE. 1
C
C     AND
C
C     (DU/DZ)(R,0) = 0 AND (DU/DZ)(R,1) = 4*R**2  FOR  0 .LE. R .LE. 1 .
C
C     THE SOLUTION TO THIS PROBLEM IS NOT UNIQUE.  IT IS A
C     ONE-PARAMETER FAMILY OF SOLUTIONS GIVEN BY
C
C            U(R,Z) = (R*Z)**4 + ARBITRARY CONSTANT .
C
C     THE R-INTERVAL WILL CONTAIN 50 UNKNOWNS AND THE Z-INTERVAL WILL
C     CONTAIN 52 UNKNOWNS.
C
      DIMENSION       F(51,52)   ,BDB(52)    ,BDC(50)    ,BDD(50)    ,
     1                W(1108)    ,R(50)      ,Z(52)
C
C     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF.  ALSO NOTE THAT W
C     IS DIMENSIONED (13 + INT(LOG2(N)))*M + 4*N .
C
      IDIMF = 51
      A = 0.
      B = 1.
      M = 50
      MBDCND = 6
      C = 0.
      D = 1.
      N = 52
      NBDCND = 3
      ELMBDA = 0.
C
C     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING
C     BOUNDARY DATA AND THE RIGHT SIDE OF THE POISSON EQUATION.
C
      DO 101 I=1,M
         R(I) = (FLOAT(I)-0.5)/50.
  101 CONTINUE
      DO 102 J=1,N
         Z(J) = (FLOAT(J)-0.5)/52.
  102 CONTINUE
C
C     GENERATE BOUNDARY DATA.
C
      DO 103 J=1,N
         BDB(J) = 4.*Z(J)**4
  103 CONTINUE
      DO 104 I=1,M
         BDC(I) = 0.
         BDD(I) = 4.*R(I)**4
  104 CONTINUE
C
C     BDA IS A DUMMY VARIABLE.
C
C     GENERATE RIGHT SIDE OF EQUATION.
C
      DO 106 I=1,M
         DO 105 J=1,N
            F(I,J) = 4.*R(I)**2*Z(J)**2*(4.*Z(J)**2+3.*R(I)**2)
  105    CONTINUE
  106 CONTINUE
      CALL HSTCYL (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,ELMBDA,F,
     1             IDIMF,PERTRB,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR BY MINIMIZING OVER ALL A THE FUNCTION
C     NORM(F(I,J) - A*1 - U(R(I),Z(J))).  THE EXACT SOLUTION IS
C                U(R,Z) = (R*Z)**4 + ARBITRARY CONSTANT.
C
      X = 0.
      DO 108 I=1,M
         DO 107 J=1,N
            X = X+F(I,J)-(R(I)*Z(J))**4
  107    CONTINUE
  108 CONTINUE
      X = X/FLOAT(M*N)
      DO 110 I=1,M
         DO 109 J=1,N
            F(I,J) = F(I,J)-X
  109    CONTINUE
  110 CONTINUE
      ERR = 0.
      DO 112 I=1,M
         DO 111 J=1,N
            X = ABS(F(I,J)-(R(I)*Z(J))**4)
            IF (X .GT. ERR) ERR = X
  111    CONTINUE
  112 CONTINUE
      PRINT 1001 , IERROR,PERTRB,ERR,W(1)
      STOP
C
 1001 FORMAT (1H1,20X,25HSUBROUTINE HSTCYL EXAMPLE///
     1        10X,46HTHE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS//
     2        32X,10HIERROR = 0/32X,20HPERTRB =-4.43114E-04/
     3        18X,34HDISCRETIZATION ERROR = 7.52796E-05/
     4        12X,32HREQUIRED LENGTH OF W ARRAY = 958//
     5        10X,32HTHE OUTPUT FROM YOUR COMPUTER IS//
     6        32X,8HIERROR =,I2/32X,8HPERTRB =,E12.5/
     7        18X,22HDISCRETIZATION ERROR =,E12.5/
     8        12X,28HREQUIRED LENGTH OF W ARRAY =,F4.0)
C
      END
