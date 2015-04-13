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
C     PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE POISTG TO
C     SOLVE THE EQUATION
C
C     (1/COS(X))(D/DX)(COS(X)(DU/DX)) + (D/DY)(DU/DY) =
C
C           2*Y**2*(6-Y**2)*SIN(X)
C
C     ON THE RECTANGLE -PI/2 .LT. X .LT. PI/2 AND
C     0 .LT. Y .LT. 1 WITH THE BOUNDARY CONDITIONS
C
C     (DU/DX) (-PI/2,Y) = (DU/DX)(PI/2,Y) = 0 , 0 .LE. Y .LE. 1  (2)
C
C     U(X,0) = 0                                           (3)
C                                 -PI/2 .LE. X .LE. PI/2
C     (DU/DY)(X,1) = 4SIN(X)                               (4)
C
C     USING FINITE DIFFERENCES ON A STAGGERED GRID WITH
C     DELTAX (= DX) = PI/40 AND DELTAY (= DY) = 1/20 .
C        TO SET UP THE FINITE DIFFERENCE EQUATIONS WE DEFINE
C     THE GRID POINTS
C
C     X(I) = -PI/2 + (I-0.5)DX            I=1,2,...,40
C
C     Y(J) = (J-O.5)DY                    J=1,2,...,20
C
C     AND LET V(I,J) BE AN APPROXIMATION TO U(X(I),Y(J)).
C     NUMBERING THE GRID POINTS IN THIS FASHION GIVES THE SET
C     OF UNKNOWNS AS V(I,J) FOR I=1,2,...,40 AND J=1,2,...,20.
C     HENCE, IN THE PROGRAM M = 40 AND N = 20.  AT THE INTERIOR
C     GRID POINT (X(I),Y(J)), WE REPLACE ALL DERIVATIVES IN
C     EQUATION (1) BY SECOND ORDER CENTRAL FINITE DIFFERENCES,
C     MULTIPLY BY DY**2, AND COLLECT COEFFICIENTS OF V(I,J) TO
C     GET THE FINITE DIFFERENCE EQUATION
C
C     A(I)V(I-1,J) + B(I)V(I,J) + C(I)V(I+1,J)
C
C     + V(I,J-1) - 2V(I,J) + V(I,J+1) = F(I,J)            (5)
C
C     WHERE S = (DY/DX)**2, AND FOR I=2,3,...,39
C
C     A(I) = S*COS(X(I)-DX/2)
C
C     B(I) = -S*(COS(X(I)-DX/2)+COS(X(I)+DX/2))
C
C     C(I) = S*COS(X(I)+DX/2)
C
C     F(I,J) = 2DY**2*Y(J)**2*(6-Y(J)**2)*SIN(X(I)) , J=1,2,...,19.
C
C        TO OBTAIN EQUATIONS FOR I = 1, WE REPLACE EQUATION (2)
C     BY THE SECOND ORDER APPROXIMATION
C
C     (V(1,J)-V(0,J))/DX = 0
C
C     AND USE THIS EQUATION TO ELIMINATE V(0,J) IN EQUATION (5)
C     TO ARRIVE AT THE EQUATION
C
C     B(1)V(1,J) + C(1)V(2,J) + V(1,J-1) - 2V(1,J) + V(1,J+1)
C
C                       = F(1,J)
C
C     WHERE
C
C     B(1) = -S*(COS(X(1)-DX/2)+COS(X(1)+DX/2))
C
C     C(1) = -B(1)
C
C     FOR COMPLETENESS, WE SET A(1) = 0.
C        TO OBTAIN EQUATIONS FOR I = 40, WE REPLACE THE DERIVATIVE
C     IN EQUATION (2) AT X=PI/2 IN A SIMILAR FASHION, USE THIS
C     EQUATION TO ELIMINATE THE VIRTUAL UNKNOWN V(41,J) IN EQUATION
C     (5) AND ARRIVE AT THE EQUATION
C
C     A(40)V(39,J) + B(40)V(40,J)
C
C     + V(40,J-1) - 2V(40,J) + V(40,J+1) = F(40,J)
C
C     WHERE
C
C     A(40) = -B(40) = -S*(COS(X(40)-DX/2)+COS(X(40)+DX/2))
C
C     FOR COMPLETENESS, WE SET C(40) = 0.  HENCE, IN THE
C     PROGRAM MPEROD = 1.
C        FOR J = 1, WE REPLACE EQUATION (3) BY THE SECOND ORDER
C     APPROXIMATION
C
C                (V(I,0) + V(I,1))/2 = 0
C
C     TO ARRIVE AT THE CONDITION
C
C                V(I,0) = -V(I,1) .
C
C     FOR J = 20, WE REPLACE EQUATION (4) BY THE SECOND ORDER
C     APPROXIMATION
C
C                (V(I,21) - V(I,20))/DY = 4*SIN(X)
C
C     AND COMBINE THIS EQUATION WITH EQUATION (5) TO ARRIVE AT
C     THE EQUATION
C
C     A(I)V(I-1,20) + B(I)V(I,20) + C(I)V(I+1,20)
C
C     + V(I,19) - 2V(I,20) + V(I,21) = F(I,20)
C
C     WHERE
C
C     V(I,21) = V(I,20)  AND
C
C     F(I,20) = 2*DY**2*Y(J)**2*(6-Y(J)**2)*SIN(X(I)) - 4*DY*SIN(X(I))
C
C     HENCE, IN THE PROGRAM NPEROD = 2 .
C        THE EXACT SOLUTION TO THIS PROBLEM IS
C
C        U(X,Y) = Y**4*COS(X) .
C
      DIMENSION       F(42,20)   ,A(40)      ,B(40)      ,C(40)      ,
     1                W(600)     ,X(40)      ,Y(20)
C
C     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF = 42.  ALSO
C     NOTE THAT W HAS BEEN DIMENSIONED
C     9M + 4N + M(INT(LOG2(N))) = 360 + 80 + 160 = 600 .
C
      IDIMF = 42
      MPEROD = 1
      M = 40
      PI = PIMACH(DUM)
      DX = PI/FLOAT(M)
      NPEROD = 2
      N = 20
      DY = 1./FLOAT(N)
C
C     GENERATE AND STORE GRID POINTS FOR COMPUTATION.
C
      DO 101 I=1,M
         X(I) = -PI/2.+(FLOAT(I)-0.5)*DX
  101 CONTINUE
      DO 102 J=1,N
         Y(J) = (FLOAT(J)-0.5)*DY
  102 CONTINUE
C
C     GENERATE COEFFICIENTS .
C
      S = (DY/DX)**2
      A(1) = 0.
      B(1) = -S*COS(-PI/2.+DX)/COS(X(1))
      C(1) = -B(1)
      DO 103 I=2,M
         A(I) = S*COS(X(I)-DX/2.)/COS(X(I))
         C(I) = S*COS(X(I)+DX/2.)/COS(X(I))
         B(I) = -(A(I)+C(I))
  103 CONTINUE
      A(40) = -B(40)
      C(40) = 0.
C
C     GENERATE RIGHT SIDE OF EQUATION.
C
      DO 105 I=1,M
         DO 104 J=1,N
            F(I,J) = 2.*DY**2*Y(J)**2*(6.-Y(J)**2)*SIN(X(I))
  104    CONTINUE
  105 CONTINUE
      DO 106 I=1,M
         F(I,N) = F(I,N)-4.*DY*SIN(X(I))
  106 CONTINUE
      CALL POISTG (NPEROD,N,MPEROD,M,A,B,C,IDIMF,F,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
C
C     U(X,Y) = Y**4*SIN(X)
C
      ERR = 0.
      DO 108 I=1,M
         DO 107 J=1,N
            T = ABS(F(I,J)-Y(J)**4*SIN(X(I)))
            IF (T .GT. ERR) ERR = T
  107    CONTINUE
  108 CONTINUE
      PRINT 1001 , IERROR,ERR,W(1)
      STOP
C
 1001 FORMAT (1H1,20X,25HSUBROUTINE POISTG EXAMPLE///
     1        10X,46HTHE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS//
     2        32X,10HIERROR = 0/
     3        18X,34HDISCRETIZATION ERROR = 5.64171E-04/
     4        12X,32HREQUIRED LENGTH OF W ARRAY = 560//
     5        10X,32HTHE OUTPUT FROM YOUR COMPUTER IS//
     6        32X,8HIERROR =,I2/18X,22HDISCRETIZATION ERROR =,E12.5/
     7        12X,28HREQUIRED LENGTH OF W ARRAY =,F4.0)
C
      END
