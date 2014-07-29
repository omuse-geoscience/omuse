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
C     PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE CMGNBN TO SOLVE
C     THE EQUATION
C
C     (1+X)**2*(D/DX)(DU/DX) - 2(1+X)(DU/DX) + (D/DY)(DU/DY)
C
C             - SQRT(-1)*U = (3 - SQRT(-1))*(1+X)**4*SIN(Y)         (1)
C
C     ON THE RECTANGLE 0 .LT. X .LT. 1 AND -PI .LT. Y .LT. PI
C     WITH THE BOUNDARY CONDITIONS
C
C     (DU/DX)(0,Y) = 4SIN(Y)                               (2)
C                                -PI .LE. Y .LE. PI
C     U(1,Y) = 16SIN(Y)                                    (3)
C
C     AND WITH U PERIODIC IN Y USING FINITE DIFFERENCES ON A
C     GRID WITH DELTAX (= DX) = 1/20 AND DELTAY (= DY) = PI/20.
C        TO SET UP THE FINITE DIFFERENCE EQUATIONS WE DEFINE
C     THE GRID POINTS
C
C     X(I) = (I-1)DX            I=1,2,...,21
C
C     Y(J) = -PI + (J-1)DY      J=1,2,...,41
C
C     AND LET V(I,J) BE AN APPROXIMATION TO U(X(I),Y(J)).
C     NUMBERING THE GRID POINTS IN THIS FASHION GIVES THE SET
C     OF UNKNOWNS AS V(I,J) FOR I=1,2,...,20 AND J=1,2,...,40.
C     HENCE, IN THE PROGRAM M = 20 AND N = 40.  AT THE INTERIOR
C     GRID POINT (X(I),Y(J)), WE REPLACE ALL DERIVATIVES IN
C     EQUATION (1) BY SECOND ORDER CENTRAL FINITE DIFFERENCES,
C     MULTIPLY BY DY**2, AND COLLECT COEFFICIENTS OF V(I,J) TO
C     GET THE FINITE DIFFERENCE EQUATION
C
C     A(I)V(I-1,J) + B(I)V(I,J) + C(I)V(I+1,J)
C
C     + V(I,J-1) - 2V(I,J) + V(I,J+1) = F(I,J)            (4)
C
C     WHERE S = (DY/DX)**2, AND FOR I=2,3,...,19
C
C     A(I) = (1+X(I))**2*S + (1+X(I))*S*DX
C
C     B(I) = -2(1+X(I))**2*S - SQRT(-1)*DY**2
C
C     C(I) = (1+X(I))**2*S - (1+X(I))*S*DX
C
C     F(I,J) = (3 - SQRT(-1))*(1+X(I))**4*DY**2*SIN(Y(J))
C              FOR J=1,2,...,40.
C
C        TO OBTAIN EQUATIONS FOR I = 1, WE REPLACE THE
C     DERIVATIVE IN EQUATION (2) BY A SECOND ORDER CENTRAL
C     FINITE DIFFERENCE APPROXIMATION, USE THIS EQUATION TO
C     ELIMINATE THE VIRTUAL UNKNOWN V(0,J) IN EQUATION (4)
C     AND ARRIVE AT THE EQUATION
C
C     B(1)V(1,J) + C(1)V(2,J) + V(1,J-1) - 2V(1,J) + V(1,J+1)
C
C                       = F(1,J)
C
C     WHERE
C
C     B(1) = -2S - SQRT(-1)*DY**2 , C(1) = 2S
C
C     F(1,J) = (11-SQRT(-1)+8/DX)*DY**2*SIN(Y(J)),  J=1,2,...,40.
C
C     FOR COMPLETENESS, WE SET A(1) = 0.
C        TO OBTAIN EQUATIONS FOR I = 20, WE INCORPORATE
C     EQUATION (3) INTO EQUATION (4) BY SETTING
C
C     V(21,J) = 16SIN(Y(J))
C
C     AND ARRIVE AT THE EQUATION
C
C     A(20)V(19,J) + B(20)V(20,J)
C
C     + V(20,J-1) - 2V(20,J) + V(20,J+1) = F(20,J)
C
C     WHERE
C
C     A(20) = (1+X(20))**2*S + (1+X(20))*S*DX
C
C     B(20) = -2*(1+X(20))**2*S - SQRT(-1)*DY**2
C
C     F(20,J) = ((3-SQRT(-1))*(1+X(20))**4*DY**2 - 16(1+X(20))**2*S
C                + 16(1+X(20))*S*DX)*SIN(Y(J))
C
C                    FOR J=1,2,...,40.
C
C     FOR COMPLETENESS, WE SET C(20) = 0.  HENCE, IN THE
C     PROGRAM MPEROD = 1.
C        THE PERIODICITY CONDITION ON U GIVES THE CONDITIONS
C
C     V(I,0) = V(I,40) AND V(I,41) = V(I,1) FOR I=1,2,...,20.
C
C     HENCE, IN THE PROGRAM NPEROD = 0.
C
C          THE EXACT SOLUTION TO THIS PROBLEM IS
C
C                  U(X,Y) = (1+X)**4*SIN(Y) .
C
      COMPLEX         F          ,A          ,B          ,C          ,W
      DIMENSION       F(22,40)   ,A(20)      ,B(20)      ,C(20)      ,
     1                X(21)      ,Y(41)      ,W(380)
C
C     FROM THE DIMENSION STATEMENT WE GET THAT IDIMF = 22 AND THAT W
C     HAS BEEN DIMENSIONED
C
C     4N + (10+INT(LOG2(N)))M = 4*20 + (10+5)*20 = 380 .
C
      IDIMF = 22
      M = 20
      MP1 = M+1
      MPEROD = 1
      DX = 0.05
      N = 40
      NPEROD = 0
      PI = PIMACH(DUM)
      DY = PI/20.
C
C     GENERATE GRID POINTS FOR LATER USE.
C
      DO 101 I=1,MP1
         X(I) = FLOAT(I-1)*DX
  101 CONTINUE
      DO 102 J=1,N
         Y(J) = -PI+FLOAT(J-1)*DY
  102 CONTINUE
C
C     GENERATE COEFFICIENTS.
C
      S = (DY/DX)**2
      DO 103 I=2,19
         T = 1.+X(I)
         TSQ = T**2
         A(I) = CMPLX((TSQ+T*DX)*S,0.)
         B(I) = -2.*TSQ*S-(0.,1.)*DY**2
         C(I) = CMPLX((TSQ-T*DX)*S,0.)
  103 CONTINUE
      A(1) = (0.,0.)
      B(1) = -2.*S-(0.,1.)*DY**2
      C(1) = CMPLX(2.*S,0.)
      B(20) = -2.*S*(1.+X(20))**2-(0.,1.)*DY**2
      A(20) = CMPLX(S*(1.+X(20))**2+(1.+X(20))*DX*S,0.)
      C(20) = (0.,0.)
C
C     GENERATE RIGHT SIDE.
C
      DO 105 I=2,19
         DO 104 J=1,N
            F(I,J) = (3.,-1.)*(1.+X(I))**4*DY**2*SIN(Y(J))
  104    CONTINUE
  105 CONTINUE
      T = 1.+X(20)
      TSQ = T**2
      T4 = TSQ**2
      DO 106 J=1,N
         F(1,J) = ((11.,-1.)+8./DX)*DY**2*SIN(Y(J))
         F(20,J) = ((3.,-1.)*T4*DY**2-16.*TSQ*S+16.*T*S*DX)*SIN(Y(J))
  106 CONTINUE
      CALL CMGNBN (NPEROD,N,MPEROD,M,A,B,C,IDIMF,F,IERROR,W)
C
C     COMPUTE DISCRETIAZATION ERROR.  THE EXACT SOLUTION IS
C
C            U(X,Y) = (1+X)**4*SIN(Y) .
C
      ERR = 0.
      DO 108 I=1,M
         DO 107 J=1,N
            T = CABS(F(I,J)-(1.+X(I))**4*SIN(Y(J)))
            IF (T .GT. ERR) ERR = T
  107    CONTINUE
  108 CONTINUE
      T = REAL(W(1))
      PRINT 1001 , IERROR,ERR,T
      STOP
C
 1001 FORMAT (1H1,20X,25HSUBROUTINE CMGNBN EXAMPLE///
     1        10X,46HTHE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS//
     2        32X,10HIERROR = 0/
     3        18X,34HDISCRETIZATION ERROR = 9.16200E-03/
     4        12X,32HREQUIRED LENGTH OF W ARRAY = 380//
     5        10X,32HTHE OUTPUT FROM YOUR COMPUTER IS//
     6        32X,8HIERROR =,I2/18X,22HDISCRETIZATION ERROR =,E12.5/
     7        12X,28HREQUIRED LENGTH OF W ARRAY =,F4.0)
C
      END
