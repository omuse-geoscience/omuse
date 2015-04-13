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
C     THIS PROGRAM ILLUSTRATES THE USE OF SUBROUTINE HSTCSP TO SOLVE
C     THE EQUATION
C
C               (1/R**2)(D/DR)(R**2(DU/DR)) +
C
C     (1/R**2*SIN(THETA))(D/DTHETA)(SIN(THETA)(DU/DTHETA))
C
C                     = 12*(R*COS(THETA))**2
C
C     ON THE RECTANGLE 0 .LT. THETA .LT. PI , 0 .LT. R .LT. 1
C     WITH THE BOUNDARY CONDITIONS
C
C     U(THETA,1) = COS(THETA)**4 , O .LE. THETA .LE. PI
C
C     AND THE SOLUTION UNSPECIFIED ON THE REMAINING BOUNDARIES.
C     WE WILL USE 45 UNKNOWNS IN THE THETA-INTERVAL AND 15 UNKNOWNS
C     IN THE R-INTERVAL.
C
      DIMENSION       F(47,16)   ,BDD(45)    ,W(615)     ,THETA(45)  ,
     1                R(15)      ,COST(45)
C
C     NOTE THAT FROM DIMENSION STATEMENT WE GET THAT IDIMF = 47  AND
C     THAT W IS DIMENSIONED ACCORDING TO THE STATEMENT IN THE
C     DESCRIPTION OF W.
C
      IDIMF = 47
      A = 0.
      B = PIMACH(DUM)
C
C     NOTE THAT B IS SET TO PI USING THE FUNCTION PIMACH AS REQUIRED.
C
      M = 45
      MBDCND = 9
      DT = (B-A)/FLOAT(M)
C
C     DEFINE GRID POINTS THETA(I) AND COS(THETA(I))
C
      DO 101 I=1,M
         THETA(I) = A+(FLOAT(I)-0.5)*DT
         COST(I) = COS(THETA(I))
  101 CONTINUE
      C = 0.
      D = 1.
      N = 15
      NBDCND = 5
      DR = (D-C)/FLOAT(N)
C
C     DEFINE GRID POINTS R(J)
C
      DO 102 J=1,N
         R(J) = C+(FLOAT(J)-0.5)*DR
  102 CONTINUE
C
C     DEFINE BOUNDARY ARRAY BDD.  BDA, BDB, AND BDC ARE DUMMY
C     VARIABLES IN THIS EXAMPLE.
C
      DO 103 I=1,M
         BDD(I) = COST(I)**4
  103 CONTINUE
      ELMBDA = 0.
C
C     DEFINE RIGHT SIDE F
C
      DO 105 I=1,M
         DO 104 J=1,N
            F(I,J) = 12.*(R(J)*COST(I))**2
  104    CONTINUE
  105 CONTINUE
      INTL = 0
      CALL HSTCSP (INTL,A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
     1             ELMBDA,F,IDIMF,PERTRB,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
C
C     U(THETA,R) = (R*COS(THETA))**4
C
      ERR = 0.
      DO 107 I=1,M
         DO 106 J=1,N
            Z = ABS(F(I,J)-(R(J)*COST(I))**4)
            IF (Z .GT. ERR) ERR = Z
  106    CONTINUE
  107 CONTINUE
      PRINT 1001 , IERROR,ERR,W(1)
      STOP
C
 1001 FORMAT (1H1,20X,25HSUBROUTINE HSTCSP EXAMPLE///
     1        10X,46HTHE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS//
     2        32X,10HIERROR = 0/
     3        18X,34HDISCRETIZATION ERROR = 5.58432E-03/
     4        12X,32HREQUIRED LENGTH OF W ARRAY = 583//
     5        10X,32HTHE OUTPUT FROM YOUR COMPUTER IS//
     6        32X,8HIERROR =,I2/18X,22HDISCRETIZATION ERROR =,E12.5/
     7        12X,28HREQUIRED LENGTH OF W ARRAY =,F4.0)
C
      END
