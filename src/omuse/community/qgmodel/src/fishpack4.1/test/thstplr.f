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
C     PROGRAM TO ILLUSTRATE THE USE OF HSTPLR TO SOLVE THE EQUATION
C
C     (1/R)(D/DR)(R*DU/DR) + (1/R**2)(D/DTHETA)(DU/DTHETA) = 16*R**2
C
C     ON THE QUARTER-DISK 0 .LT. R .LT. 1 AND 0 .LT. THETA .LT. PI/2
C     WITH THE BOUNDARY CONDITIONS
C
C     U(1,THETA) = 1 - COS(4*THETA)  FOR  0 .LE. THETA .LE. PI/2
C
C     AND
C
C     (DU/DR)(R,0) = (DU/DR)(R,PI/2) = 0  FOR  O .LT. R .LT. 1 .
C
C     NOTE THAT U AT THE ORIGIN IS UNSPECIFIED.  THE EXACT SOLUTION TO
C     THIS PROBLEM IS
C
C               U(R,THETA) = (R**4)(1-COS(4*THETA)) .
C
C     WE WILL USE 50 UNKNOWNS IN THE R-INTERVAL AND 48 UNKNOWNS IN
C     THE THETA-INTERVAL.
C
      DIMENSION       F(51,50)   ,BDB(48)    ,BDC(50)    ,BDD(50)    ,
     1                W(1092)    ,R(50)      ,THETA(48)
C
C     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF.  ALSO NOTE THAT W
C     IS DIMENSIONED (13 + LOG2(N))*M + 4*N .
C
      IDIMF = 51
      A = 0.
      B = 1.
      M = 50
      MBDCND = 5
      C = 0.
      PI = PIMACH(DUM)
      D = PI/2.
      N = 48
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
         THETA(J) = (FLOAT(J)-0.5)*PI/96.
  102 CONTINUE
C
C     GENERATE BOUNDARY DATA.
C
      DO 103 J=1,N
         BDB(J) = 1.-COS(4.*THETA(J))
  103 CONTINUE
      DO 104 I=1,M
         BDC(I) = 0.
         BDD(I) = 0.
  104 CONTINUE
C
C     BDA IS A DUMMY VARIABLE.
C
C
C     GENERATE RIGHT SIDE OF EQUATION.
C
      DO 106 I=1,M
         DO 105 J=1,N
            F(I,J) = 16.*R(I)**2
  105    CONTINUE
  106 CONTINUE
      CALL HSTPLR (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,ELMBDA,F,
     1             IDIMF,PERTRB,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
C
C                U(R,THETA) = R**4*(1 - COS(4*THETA))
C
      ERR = 0.
      DO 108 I=1,M
         DO 107 J=1,N
            Z = ABS(F(I,J)-R(I)**4*(1.-COS(4.*THETA(J))))
            IF (Z .GT. ERR) ERR = Z
  107    CONTINUE
  108 CONTINUE
      PRINT 1001 , IERROR,ERR,W(1)
      STOP
C
 1001 FORMAT (1H1,20X,25HSUBROUTINE HSTPLR EXAMPLE///
     1        10X,46HTHE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS//
     2        32X,10HIERROR = 0/
     3        18X,34HDISCRETIZATION ERROR = 1.13038E-03/
     4        12X,33HREQUIRED LENGTH OF W ARRAY = 1042//
     5        10X,32HTHE OUTPUT FROM YOUR COMPUTER IS//
     6        32X,8HIERROR =,I2/18X,22HDISCRETIZATION ERROR =,E12.5/
     7        12X,28HREQUIRED LENGTH OF W ARRAY =,F5.0)
C
      END
