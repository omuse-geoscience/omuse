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
C        PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE HW3CRT TO
C     SOLVE THE HELMHOLTZ EQUATION
C
C     (D/DX)(DU/DX) + (D/DY)(DU/DY) + (D/DZ)(DU/DZ) - 3U
C
C       =  4X**2*(3-X**2):SIN(Y)*COS(Z)
C
C     ON THE PARALLELEPIPED 0 .LT. X .LT. 1, 0 .LT. Y .LT. 2*PI,
C     0 .LT. Z .LT. PI/2 WITH THE BOUNDARY CONDITIONS
C
C     U(0,Y,Z) = 0
C                           0 .LE. Y .LE. 2*PI , 0 .LE. Z .LE. PI/2
C     U(1,Y,Z) = SIN(Y)*COS(Z)
C
C     U PERIODIC IN Y,
C
C     U(X,Y,0) = X**4*SIN(Y)
C                           0 .LE. X .LE. 1 , 0 .LE. Y .LE. 2*PI
C     (DU/DX)(X,Y,PI/2) = -X**4*SIN(Y)
C
C     USING A FINITE DIFFERENCE GRID WITH PANEL WIDTHS
C
C     DELTAX (=DX) = 1/10 , DELTAY (=DY) = 1/40 , DELTAZ (=DZ) = 1/15.
C
C        THE EXACT SOLUTION OF THIS PROBLEM IS
C
C     U(X,Y,Z) = X**4*SIN(Y)*COS(Z) .
C
      DIMENSION       F(11,41,16),BDZF(11,41),W(370)     ,X(11)      ,
     1                Y(41)      ,Z(16)
C
C        FROM THE DESCRIPTION OF THE PROBLEM GIVEN ABOVE, WE DEFINE
C     THE FOLLOWING QUANTITIES
C
      ELMBDA = -3.
      XS = 0.
      XF = 1.
      LBDCND = 1
      YS = 0.
      PI = PIMACH(DUM)
      YF = 2.*PI
      MBDCND = 0
      ZS = 0.
      ZF = PI/2.
      NBDCND = 2
      L = 10
      M = 40
      N = 15
C
C     FROM THE DIMENSION STATEMENT ABOVE WE DEFINE
C
      LDIMF = 11
      MDIMF = 41
C
C     ALSO NOTE THAT W HAS BEEN DIMENSIONED
C        30+L+M+5*N+MAX(L,M,N)+7*(INT((L+1)/2)+INT((M+1)/2))
C      = 30+10+40+75+40+7*(5+20) = 370
C     WE DEFINE THE GRID POINTS FOR LATER USE.
C
      LP1 = L+1
      DX = (XF-XS)/FLOAT(L)
      DO 101 I=1,LP1
         X(I) = XS+FLOAT(I-1)*DX
  101 CONTINUE
      MP1 = M+1
      DY = (YF-YS)/FLOAT(M)
      DO 102 J=1,MP1
         Y(J) = YS+FLOAT(J-1)*DY
  102 CONTINUE
      NP1 = N+1
      DZ = (ZF-ZS)/FLOAT(N)
      DO 103 K=1,NP1
         Z(K) = ZS+FLOAT(K-1)*DZ
  103 CONTINUE
C
C     WE DEFINE THE ARRAY OF DERIVATIVE BOUNDARY VALUES.
C
      DO 105 I=1,LP1
         DO 104 J=1,MP1
            BDZF(I,J) = -X(I)**4*SIN(Y(J))
  104    CONTINUE
  105 CONTINUE
C
C     NOTE THAT FOR THIS EXAMPLE ALL OTHER BOUNDARY ARRAYS ARE
C     DUMMY VARIABLES.
C     WE DEFINE THE FUNCTION BOUNDARY VALUES IN THE F ARRAY.
C
      DO 107 J=1,MP1
         DO 106 K=1,NP1
            F(1,J,K) = 0.
            F(LP1,J,K) = SIN(Y(J))*COS(Z(K))
  106    CONTINUE
  107 CONTINUE
      DO 109 I=1,LP1
         DO 108 J=1,MP1
            F(I,J,1) = X(I)**4*SIN(Y(J))
  108    CONTINUE
  109 CONTINUE
C
C     WE NOW DEFINE THE VALUES OF THE RIGHT SIDE OF THE HELMHOLTZ
C     EQUATION.
C
      DO 112 I=2,L
         DO 111 J=1,MP1
            DO 110 K=2,NP1
               F(I,J,K) = 4.*X(I)**2*(3.-X(I)**2)*SIN(Y(J))*COS(Z(K))
  110       CONTINUE
  111    CONTINUE
  112 CONTINUE
C
C     CALL HW3CRT TO GENERATE AND SOLVE THE FINITE DIFFERENCE EQUATION.
C
      CALL HW3CRT (XS,XF,L,LBDCND,BDXS,BDXF,YS,YF,M,MBDCND,BDYS,BDYF,
     1             ZS,ZF,N,NBDCND,BDZS,BDZF,ELMBDA,LDIMF,MDIMF,F,
     2             PERTRB,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION TO THE
C     PROBLEM IS
C
C        U(X,Y,Z) = X**4*SIN(Y)*COS(Z)
C
      ERR = 0.
      DO 115 I=1,LP1
         DO 114 J=1,MP1
            DO 113 K=1,NP1
               T = ABS(F(I,J,K)-X(I)**4*SIN(Y(J))*COS(Z(K)))
               IF (T .GT. ERR) ERR = T
  113       CONTINUE
  114    CONTINUE
  115 CONTINUE
      PRINT 1001 , IERROR,ERR
      STOP
C
 1001 FORMAT (1H1,20X,25HSUBROUTINE HW3CRT EXAMPLE///
     1        10X,46HTHE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS//
     2        32X,10HIERROR = 0/
     3        18X,34HDISCRETIZATION ERROR = 9.64802E-03//
     4        10X,32HTHE OUTPUT FROM YOUR COMPUTER IS//
     5        32X,8HIERROR =,I2/18X,22HDISCRETIZATION ERROR =,E12.5)
C
      END
