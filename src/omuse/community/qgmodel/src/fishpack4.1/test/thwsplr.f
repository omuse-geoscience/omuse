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
C          PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE HWSPLR TO SOLVE
C     THE EQUATION
C
C     (1/R)(D/DR)(R*(DU/DR)) + (1/R**2)(D/DTHETA)(DU/DTHETA) = 16*R**2
C
C     ON THE QUARTER-DISK 0 .LT. R .LT. 1, 0 .LT. THETA .LT. PI/2 WITH
C     WITH THE BOUNDARY CONDITIONS
C
C     U(1,THETA) = 1 - COS(4*THETA), 0 .LE. THETA .LE. 1
C
C     AND
C
C     (DU/DTHETA)(R,0) = (DU/DTHETA)(R,PI/2) = 0,  0 .LE. R .LE. 1.
C
C     (NOTE THAT THE SOLUTION U IS UNSPECIFIED AT R = 0.)
C          THE R-INTERVAL WILL BE DIVIDED INTO 50 PANELS AND THE
C     THETA-INTERVAL WILL BE DIVIDED INTO 48 PANELS.
C
      DIMENSION       F(100,50)  ,BDC(51)    ,BDD(51)    ,W(1114)    ,
     1                R(51)      ,THETA(49)
C
C     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF.  ALSO NOTE THAT W
C     IS DIMENSIONED 4*(N+1) + (13 + INT(LOG2(N+1)))*(M+1) .
C
      IDIMF = 100
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
         THETA(J) = FLOAT(J-1)*PI/96.
  102 CONTINUE
C
C     GENERATE BOUNDARY DATA.
C
      DO 103 I=1,MP1
         BDC(I) = 0.
         BDD(I) = 0.
  103 CONTINUE
C
C     BDA AND BDB ARE DUMMY VARIABLES.
C
      DO 104 J=1,NP1
         F(MP1,J) = 1.-COS(4.*THETA(J))
  104 CONTINUE
C
C     GENERATE RIGHT SIDE OF EQUATION.
C
      DO 106 I=1,M
         DO 105 J=1,NP1
            F(I,J) = 16.*R(I)**2
  105    CONTINUE
  106 CONTINUE
      CALL HWSPLR (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,ELMBDA,F,
     1             IDIMF,PERTRB,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
C                U(R,THETA) = R**4*(1 - COS(4*THETA))
C
      ERR = 0.
      DO 108 I=1,MP1
         DO 107 J=1,NP1
            Z = ABS(F(I,J)-R(I)**4*(1.-COS(4.*THETA(J))))
            IF (Z .GT. ERR) ERR = Z
  107    CONTINUE
  108 CONTINUE
      PRINT 1001 , IERROR,ERR,W(1)
      STOP
C
 1001 FORMAT (1H1,20X,25HSUBROUTINE HWSPLR EXAMPLE///
     1        10X,46HTHE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS//
     2        32X,10HIERROR = 0/
     3        18X,34HDISCRETIZATION ERROR = 6.19134E-04/
     4        12X,32HREQUIRED LENGTH OF W ARRAY = 882//
     5        10X,32HTHE OUTPUT FROM YOUR COMPUTER IS//
     6        32X,8HIERROR =,I2/18X,22HDISCRETIZATION ERROR =,E12.5/
     7        12X,28HREQUIRED LENGTH OF W ARRAY =,F4.0)
C
      END
