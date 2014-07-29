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
C     THIS PROGRAM ILLUSTRATES THE USE OF THE SUBROUTINE POIS3D TO SOLVE
C     THE EQUATION
C
C     (D/DX)(DU/DX) + (D/DY)(DU/DY) + (1+Z)**2*(D/DZ)(DU/DZ)
C
C     - 2*(1+Z)*(DU/DZ) = 2*SIN(X)*SIN(Y)*(1+Z)**4                   (1)
C
C     ON THE PARALLELEPIPED -PI .LT. X .LT. PI, -PI .LT. Y .LT. PI,
C     0 .LT. Z .LT. 1 WITH BOUNDARY CONDITIONS
C
C     U PERIODIC IN X
C
C     U PERIODIC IN Y
C
C     (DU/DZ)(X,Y,0) = 4*SIN(X)*SIN(Y)    -PI .LT. X,Y .LT. PI       (2)
C
C     U(X,Y,1) = 16*SIN(X)*SIN(Y)         -PI .LT. X,Y .LT. PI       (3)
C
C     USING A FINITE DIFFERENCE GRID WITH DELTAX (= DX) = 2*PI/30 ,
C     DELTAY (= DY) = 2*PI/30, AND DELTAZ (= DZ) = 1/10.
C         TO SET UP THE FINITE DIFFERENCE EQUATIONS WE DEFINE THE GRID
C     POINTS
C
C         X(I) = -PI + (I-1)*DX         I=1,2,...,31
C
C         Y(J) = -PI + (J-1)*DY         J=1,2,...,31
C
C         Z(K) = (K-1)*DZ               K=1,2,...,11
C
C     AND LET V(I,J,K) BE AN APPROXIMATION TO U(X(I),Y(J),Z(K)).
C     NUMBERING THE GRID POINTS IN THIS FASHION GIVES THE SET OF
C     UNKNOWNS AS V(I,J,K) FOR I=1,2,...,30, J=1,2,...,30, K=1,2,...,10.
C     HENCE, IN THE PROGRAM L=30, M = 30, AND N = 10.  AT THE INTERIOR
C     GRID POINT (X(I),Y(J),Z(K)), WE REPLACE ALL DERIVATIVES IN
C     EQUATION (1) BY SECOND ORDER CENTRAL FINITE DIFFERENCES AND
C     COLLECT COEFFICIENTS OF V(I,J,K) TO GET THE FINITE DIFFERENCE
C     EQUATION
C
C        (V(I-1,J,K) - 2V(I,J,K) + V(I+1,J,K))/DX**2
C
C      + (V(I,J-1,K) - 2V(I,J,K) + V(I,J+1,K))/DY**2
C
C      + A(K)V(I,J,K-1) + B(K)V(I,J,K) + C(K)V(I,J,K+1) = F(I,J,K)   (4)
C
C     WHERE FOR K=2,3,...,9
C
C     A(K) = (1+Z(K))**2/DZ**2 + (1+Z(K))/DZ
C
C     B(K) = -2(1+Z(K))**2/DZ**2
C
C     C(K) = (1+Z(K))**2/DZ**2 - (1+Z(K))/DZ
C
C     F(I,J,K) = 2SIN(X(I))*SIN(Y(J))*(1+Z(K))**4  FOR I,J=1,2,...,30.
C
C         TO OBTAIN EQUATIONS FOR K=1, WE REPLACE THE DERIVATIVE IN
C     EQUATION (2) BY A SECOND ORDER CENTRAL FINITE DIFFERENCE APPROX-
C     IMATION, USE THIS EQUATION TO ELIMINATE THE VIRTUAL UNKNOWN
C     V(I,J,0) IN EQUATION (4) AND ARRIVE AT THE EQUATION
C
C              (V(I-1,J,1) -2V(I,J,1) + V(I+1,J,1))/DX**2
C
C            + (V(I,J-1,1) -2V(I,J,1) + V(I,J+1,1))/DY**2
C
C            + B(1)V(I,J,1) + C(1)V(I,J,2) = F(I,J,1)
C
C     WHERE
C            B(1) = -C(1) = -2(1+Z(1))**2/DZ**2 = -2/DZ**2
C
C            F(I,J,1) = (10 + 8/DZ)SIN(X(I))*SIN(Y(J))
C
C            FOR I,J=1,2,...,30.  FOR COMPLETENESS WE SET A(1) = 0.
C
C         TO OBTAIN EQUATIONS FOR K=10, WE INCORPORATE EQUATION (3) INTO
C     EQUATION (4) BY SETTING
C
C            V(I,J,11) = U(X(I),Y(J),1) = 16SIN(X(I))*SIN(Y(J))
C
C     AND ARRIVE AT THE EQUATION
C
C              (V(I-1,J,10) - 2V(I,J,10) + V(I+1,J,10))/DX**2
C
C            + (V(I,J-1,10) - 2V(I,J,10) + V(I,J+1,10))/DY**2
C
C            + A(10)V(I,J,9) + B(10)V(I,J,10) = F(I,J,10)
C
C     WHERE
C
C            A(10) = (1+Z(10))**2/DZ**2 + (1+Z(10))/DZ
C
C            B(10) = -2(1+Z(10))**2/DZ**2
C
C            F(I,J,10) = 2SIN(X(I))*SIN(Y(J))*((1+Z(10))**4
C                        -8*((1+Z(10))**2/DZ**2 - (1+Z(10))/DZ))
C
C                        FOR I,J=1,2,...,30.
C
C     FOR COMPLETENESS, WE SET C(10) = 0.  HENCE, IN THE PROGRAM,
C     NPEROD = 1.
C         THE PERIODICITY CONDITIONS ON U GIVE THE CONDITIONS
C
C             V(0,J,K) = V(30,J,K) AND V(31,J,K) = V(1,J,K)
C                        FOR J=1,2,...,30 AND K=1,2,...,10,
C             AND
C             V(I,0,K) = V(I,30,K) AND V(I,31,K) = V(I,1,K)
C                        FOR I=1,2,...,30 AND K=1,2,...,10.
C
C     HENCE, IN THE PROGRAM LPEROD = MPEROD = 0.
C
      DIMENSION       F(32,33,10),A(10)      ,B(10)      ,C(10)      ,
     1                W(350)     ,X(30)      ,Y(30)      ,Z(10)
C
C     FROM THE DIMENSION STATEMENT WE GET THAT LDIMF = 32, MDIMF = 33,
C     AND NOTE THAT W HAS BEEN DIMENSIONED ACCORDING TO ITS DESCRIPTION.
C
      LDIMF = 32
      MDIMF = 33
      PI = PIMACH(DUM)
      LPEROD = 0
      L = 30
      DX = 2.*PI/FLOAT(L)
      C1 = 1./DX**2
      MPEROD = 0
      M = 30
      DY = 2.*PI/FLOAT(M)
      C2 = 1./DY**2
      NPEROD = 1
      N = 10
      DZ = 1./FLOAT(N)
      DZSQ = 1./DZ**2
C
C     GENERATE GRID POINTS FOR LATER USE.
C
      DO 101 I=1,L
         X(I) = -PI+FLOAT(I-1)*DX
  101 CONTINUE
      DO 102 J=1,M
         Y(J) = -PI+FLOAT(J-1)*DY
  102 CONTINUE
C
C     GENERATE COEFFICIENTS
C
      A(1) = 0.
      B(1) = -2.*DZSQ
      C(1) = -B(1)
      Z(1) = 0.
      DO 103 K=2,N
         Z(K) = FLOAT(K-1)*DZ
         T = 1.+Z(K)
         A(K) = T**2*DZSQ+T/DZ
         B(K) = -2.*T**2*DZSQ
         C(K) = T**2*DZSQ-T/DZ
  103 CONTINUE
C
C     GENERATE RIGHT SIDE OF EQUATION
C
      DO 106 I=1,L
         DO 105 J=1,M
            DO 104 K=2,N
               F(I,J,K) = 2.*SIN(X(I))*SIN(Y(J))*(1.+Z(K))**4
  104       CONTINUE
  105    CONTINUE
  106 CONTINUE
      DO 108 I=1,L
         DO 107 J=1,L
            F(I,J,1) = (10.+8./DZ)*SIN(X(I))*SIN(Y(J))
            F(I,J,N) = F(I,J,N)-C(N)*16.*SIN(X(I))*SIN(Y(J))
  107    CONTINUE
  108 CONTINUE
      C(N) = 0.
C
C     CALL POIS3D TO SOLVE EQUATIONS.
C
      CALL POIS3D (LPEROD,L,C1,MPEROD,M,C2,NPEROD,N,A,B,C,LDIMF,MDIMF,
     1             F,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
C
C              U(X,Y,Z) = SIN(X)*SIN(Y)*(1+Z)**4
C
      ERR = 0.
      DO 111 I=1,L
         DO 110 J=1,M
            DO 109 K=1,N
               T = ABS(F(I,J,K)-SIN(X(I))*SIN(Y(J))*(1.+Z(K))**4)
               IF (T .GT. ERR) ERR = T
  109       CONTINUE
  110    CONTINUE
  111 CONTINUE
      PRINT 1001 , IERROR,ERR
      STOP
C
 1001 FORMAT (1H1,20X,25HSUBROUTINE POIS3D EXAMPLE///
     1        10X,46HTHE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS//
     2        32X,10HIERROR = 0/
     3        18X,34HDISCRETIZATION ERROR = 2.93277E-02//
     4        10X,32HTHE OUTPUT FROM YOUR COMPUTER IS//
     5        32X,8HIERROR =,I2/18X,22HDISCRETIZATION ERROR =,E12.5)
C
      END
