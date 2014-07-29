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
C
C     PROGRAM TO ILLUSTRATE THE USE OF HWSSSP
C
      DIMENSION       F(19,73)   ,BDTF(73)   ,SINT(19)   ,SINP(73)   ,
     1                W(600)
C
      PI = PIMACH(DUM)
      TS = 0
      TF = PI/2.
      M = 18
      MBDCND = 6
      PS = 0
      PF = PI+PI
      N = 72
      NBDCND = 0
      ELMBDA = 0.
      IDIMF = 19
C
C     GENERATE SINES FOR USE IN SUBSEQUENT COMPUTATIONS
C
      DTHETA = TF/FLOAT(M)
      MP1 = M+1
      DO 101 I=1,MP1
         SINT(I) = SIN(FLOAT(I-1)*DTHETA)
  101 CONTINUE
      DPHI = (PI+PI)/FLOAT(N)
      NP1 = N+1
      DO 102 J=1,NP1
         SINP(J) = SIN(FLOAT(J-1)*DPHI)
  102 CONTINUE
C
C     COMPUTE RIGHT SIDE OF EQUATION AND STORE IN F
C
      DO 104 J=1,NP1
         DO 103 I=1,MP1
            F(I,J) = 2.-6.*(SINT(I)*SINP(J))**2
  103    CONTINUE
  104 CONTINUE
C
C     STORE DERIVATIVE DATA AT THE EQUATOR
C
      DO 105 J=1,NP1
         BDTF(J) = 0.
  105 CONTINUE
C
      CALL HWSSSP (TS,TF,M,MBDCND,BDTS,BDTF,PS,PF,N,NBDCND,BDPS,BDPF,
     1             ELMBDA,F,IDIMF,PERTRB,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR. SINCE PROBLEM IS SINGULAR, THE
C     SOLUTION MUST BE NORMALIZED.
C
      ERR = 0
      DO 107 J=1,NP1
         DO 106 I=1,MP1
            Z = ABS(F(I,J)-(SINT(I)*SINP(J))**2-F(1,1))
            IF (Z .GT. ERR) ERR = Z
  106    CONTINUE
  107 CONTINUE
C
      IW = INT(W(1))
      PRINT 1001 , IERROR,ERR,IW
      STOP
C
C
 1001 FORMAT (1H1,20X,25HSUBROUTINE HWSSSP EXAMPLE///
     1        10X,46HTHE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS//
     2        32X,10HIERROR = 0/
     3        18X,34HDISCRETIZATION ERROR = 3.38107E-03/
     4        12X,32HREQUIRED LENGTH OF W ARRAY = 600//
     5        10X,32HTHE OUTPUT FROM YOUR COMPUTER IS//
     6        32X,8HIERROR =,I2/18X,22HDISCRETIZATION ERROR =,E12.5 /
     7        12X,28HREQUIRED LENGTH OF W ARRAY =,I4)
C
      END
