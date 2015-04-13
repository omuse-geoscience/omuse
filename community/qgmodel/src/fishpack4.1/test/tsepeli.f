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
C     AN EXAMPLE SHOWING THE USE OF SEPELI TO SOLVE THE SEPARABLE
C     ELLIPTIC PARTIAL DIFFERENTIAL EQUATION . . .
C     (X+1)**2*UXX+2*(X+1)*UX+EXP(Y)*UYY-(X+Y)*U = G(X,Y) ON
C     0.LE.X.LE.1, 0.LE.Y.LE.1  WITH SPECIFIED BOUNDARY CONDITIONS
C     AT Y=0,1 AND MIXED BOUNDARY CONDITIONS OF THE FORM
C     UX(0,Y)+U(0,Y), UX(1,Y)+U(1,Y) AT X=0,1.
C     THE APPROXIMATION IS GENERATED ON A UNIFORM 33 BY 33 GRID.
C     THE EXACT SOLUTION U(X,Y)=(X*Y)**3+1 IS USED TO SET THE
C     RIGHT HAND SIDE, BOUNDARY CONDITIONS, AND COMPUTE  SECOND AND
C     FOURTH ORDER DISCRETIZATION ERROR
C     THE EXACT WORK SPACE LENGTH REQUIRED IS 1118 WORDS.
C     THIS WAS DETERMINED BY A PREVIOUS CALL TO SEPELI  AND PRINT
C     OUT OF W(1).
C
      DIMENSION       USOL(33,33),GRHS(33,33),BDA(33)    ,BDB(33)    ,
     1                W(1118)
C
C     DECLARE COEFFICIENT SUBROUTINES EXTERNAL
C
      EXTERNAL        COFX       ,COFY
C
C     DEFINE ARITHMETIC FUNCTIONS GIVING EXACT SOLUTION
C
      UE(S,T) = (S*T)**3+1.0
      UXE(S,T) = 3.0*S**2*T**3
      UXXE(S,T) = 6.0*S*T**3
      UYE(S,T) = 3.0*S**3*T**2
      UYYE(S,T) = 6.0*S**3*T
C
C     SET LIMITS ON REGION
C
      A = 0.0
      B = 1.0
      C = 0.0
      D = 1.0
C
C     SET GRID SIZE
C
      M = 32
      N = 32
      DLX = (B-A)/FLOAT(M)
      DLY = (D-C)/FLOAT(N)
      NX = M+1
      NY = N+1
      DO 102 I=1,NX
         X = A+FLOAT(I-1)*DLX
C
C     SET SPECIFIED BOUNDARY CONDITIONS AT Y=C,D
C
         USOL(I,1) = UE(X,C)
         USOL(I,NY) = UE(X,D)
         CALL COFX (X,AF,BF,CF)
         DO 101 J=1,NY
            Y = C+FLOAT(J-1)*DLY
            CALL COFY (Y,DF,EF,FF)
C
C     SET RIGHT HAND SIDE
C
            GRHS(I,J) = AF*UXXE(X,Y)+BF*UXE(X,Y)+CF*UE(X,Y)+
     1                  DF*UYYE(X,Y)+EF*UYE(X,Y)+FF*UE(X,Y)
  101    CONTINUE
  102 CONTINUE
C
C     SET MIXED BOUNDARY CONDITIONS AT X=A,B
C
      ALPHA = 1.0
      BETA = 1.0
      DO 103 J=1,NY
         Y = C+FLOAT(J-1)*DLY
         BDA(J) = UXE(A,Y)+ALPHA*UE(A,Y)
         BDB(J) = UXE(B,Y)+BETA*UE(B,Y)
  103 CONTINUE
C
C     SET BOUNDARY SWITHCES
C
      MBDCND = 3
      NBDCND = 1
C
C     SET FIRST DIMENSION OF USOL,GRHS AND WORK SPACE LENGTH
C
      IDMN = 33
      W(1) = 1118.
C
C     SET WORK SPACE LENGTH IN FIRST WORD
C     SET INITIAL CALL PARAMETER TO ZERO
C
      INTL = 0
C
C     OBTAIN SECOND ORDER APPROXIMATION
C
      IORDER = 2
      CALL SEPELI (INTL,IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,
     1             NBDCND,DUM,DUM,DUM,DUM,COFX,COFY,GRHS,USOL,IDMN,W,
     2             PERTRB,IERROR)
      ERR = 0.0
      DO 105 I=1,NX
         X = A+FLOAT(I-1)*DLX
         DO 104 J=1,NY
            Y = C+FLOAT(J-1)*DLY
            ERR = AMAX1(ERR,ABS((USOL(I,J)-UE(X,Y))/UE(X,Y)))
  104    CONTINUE
  105 CONTINUE
      ERR2 = ERR
C
C     OBTAIN FOURTH ORDER APPROXIMATION
C
      IORDER = 4
C
C     NON-INITIAL CALL
C
      INTL = 1
      CALL SEPELI (INTL,IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,
     1             NBDCND,DUM,DUM,DUM,DUM,COFX,COFY,GRHS,USOL,IDMN,W,
     2             PERTRB,IERROR)
C
C     COMPUTE DISCRETIZATION ERROR
C
      ERR = 0.0
      DO 107 J=1,NY
         Y = C+FLOAT(J-1)*DLY
         DO 106 I=1,NX
            X = A+FLOAT(I-1)*DLX
            ERR = AMAX1(ERR,ABS((USOL(I,J)-UE(X,Y))/UE(X,Y)))
  106    CONTINUE
  107 CONTINUE
      ERR4 = ERR
      IW = INT(W(1))
      PRINT 1001 , IERROR,ERR2,ERR4,IW
C
C
 1001 FORMAT (1H1,20X,25HSUBROUTINE SEPELI EXAMPLE///
     1        20X,46HTHE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS//
     2        20X,10HIERROR = 0/
     3        20X,47HSECOND ORDER DISCRETIZATION ERROR = 9.78910E-05/
     4        20X,47HFOURTH ORDER DISCRETIZATION ERROR = 1.47351E-06/
     5        20X,33HREQUIRED LENGTH OF W ARRAY = 1118//
     6        20X,32HTHE OUTPUT FROM YOUR COMPUTER IS//
     7        20X,8HIERROR =,I2/
     8        20X,36HSECOND ORDER DISCRETIZATION ERROR = , E12.5/
     9        20X,36HFOURTH ORDER DISCRETIZATION ERROR = , E12.5/
     +        20X,29HREQUIRED LENGTH OF W ARRAY = , I4)
C
      END
      SUBROUTINE COFX (X,AF,BF,CF)
C
C     SET COEFFICIENTS IN THE X-DIRECTION.
C
      AF = (X+1.)**2
      BF = 2.0*(X+1.)
      CF = -X
      RETURN
      END
      SUBROUTINE COFY (Y,DF,EF,FF)
C
C     SET COEFFICIENTS IN Y DIRECTION
C
      DF = EXP(Y)
      EF = 0.0
      FF = -Y
      RETURN
      END
