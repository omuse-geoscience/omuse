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
C     EXAMPLE SHOWING THE USE OF SEPX4 TO SOLVE THE ELLIPTIC PDE
C     (X+1)**2*UXX+2*(X+1)*UX+UYY-X*U=G(X,Y) ON THE REGION
C     0.LE.X.LE.1, 0.LE.Y.LE.1  WITH SPECIFIED BOUNDARY CONDITIONS
C     AT Y=0,1 AND MIXED BOUNDARY CONDITIONS OF THE FORM
C     UX(0,Y)+U(0,Y), UX(1,Y)+U(1,Y) AT X=0,1.
C     THE APPROXIMATION IS GENERATED ON A UNIFORM 33 BY 33 GRID.
C     THE EXACT SOLUTION U(X,Y)=(X*Y)**3+1 IS USED TO SET THE
C     RIGHT HAND SIDE, BOUNDARY CONDITIONS, AND COMPUTE  SECOND AND
C     FOURTH ORDER DISCRETIZATION ERROR
C     THE EXACT WORK SPACE LENGTH REQUIRED IS 1005 WORDS.  THIS
C     WAS DETERMINED BY A PREVIOUS CALL TO SEPX4 AND PRINT OUT OF
C     W(1).
C
      DIMENSION       USOL(33,33),GRHS(33,33),BDA(33)    ,BDB(33)    ,
     1                W(1024)
      EXTERNAL COFX4
C
C     DEFINE ARITHMETIC FUNCTIONS GIVING EXACT SOLUTION
C
      UE(S,T)=(S*T)**3+1.0
      UXE(S,T)=3.0*S**2*T**3
      UXXE(S,T)=6.0*S*T**3
      UYE(S,T)=3.0*S**3*T**2
      UYYE(S,T)=6.0*S**3*T
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
          CALL COFX4 (X,AF,BF,CF)
              DO 101 J=1,NY
              Y = C+FLOAT(J-1)*DLY
C
C     SET RIGHT HAND SIDE
C
              GRHS(I,J) = AF*UXXE(X,Y)+BF*UXE(X,Y)+CF*UE(X,Y)+UYYE(X,Y)
  101         CONTINUE
  102     CONTINUE
C
C     SET MIXED BOUNDARY CONDITIONS AT X=A,B
C
      ALPHA = 1.0
      BETA = 1.0
          DO 103 J=1,NY
          Y = C+FLOAT(J-1)*DLY
          BDA(J) = UXE(A,Y)+ALPHA*UE(A,Y)
          BDB(J) = UXE(B,Y)+BETA*UE(B,Y)
  103     CONTINUE
C
C     SET BOUNDARY SWITHCES
C
      MBDCND = 3
      NBDCND = 1
C
C     SET FIRST DIMENSION OF USOL,GRHS AND WORK SPACE LENGTH
C
      IDMN = 33
      W(1) = 1024.
C
C     OBTAIN SECOND ORDER APPROXIMATION
C
      IORDER = 2
      CALL SEPX4 (IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,NBDCND,
     1            DUM,DUM,COFX4,GRHS,USOL,IDMN,W,PERTRB,IERROR)
C
C     COMPUTE SECOND ORDER DISCRETIZATION ERROR (RELATIVE)
C     ALSO RESET SPECIFIED BOUNDARIES AND RIGHT HAND SIDE.
C
      ERR = 0.0
          DO 105 I=1,NX
          X = A+FLOAT(I-1)*DLX
          USOL(I,1) = UE(X,C)
          USOL(I,NY) = UE(X,D)
          CALL COFX4 (X,AF,BF,CF)
              DO 104 J=1,NY
              Y = C+FLOAT(J-1)*DLY
              ERR = AMAX1(ERR,ABS((USOL(I,J)-UE(X,Y))/UE(X,Y)))
              GRHS(I,J) = AF*UXXE(X,Y)+BF*UXE(X,Y)+CF*UE(X,Y)+UYYE(X,Y)
  104         CONTINUE
  105     CONTINUE
      ERR2=ERR
C
C     OBTAIN FOURTH ORDER APPROXIMATION
C
      IORDER = 4
      CALL SEPX4 (IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,NBDCND,
     1            DUM,DUM,COFX4,GRHS,USOL,IDMN,W,PERTRB,IERROR)
C
C     COMPUTE FOURTH ORDER DISCRETIZATION ERROR (RELATIVE)
C
      ERR = 0.0
          DO 107 J=1,NY
          Y = C+FLOAT(J-1)*DLY
              DO 106 I=1,NX
              X = A+FLOAT(I-1)*DLX
              ERR = AMAX1(ERR,ABS((USOL(I,J)-UE(X,Y))/UE(X,Y)))
  106         CONTINUE
  107     CONTINUE
      ERR4=ERR
      IW = INT(W(1)+0.5)
      PRINT 1001,IERROR,ERR2,ERR4,IW
 1001 FORMAT(1H1,20X,25HSUBROUTINE SEPX4  EXAMPLE ///
     120X,46HTHE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS //
     220X,10HIERROR = 0 /
     320X,48HSECOND ORDER DISCRETIZATION ERROR =  1.5985E-04  /
     420X,48HFOURTH ORDER DISCRETIZATION ERROR =  1.85749E-06  /
     520X,33HREQUIRED LENGTH OF W ARRAY = 1024 //
     620X, 32HTHE OUTPUT FROM YOUR COMPUTER IS //
     720X, 8HIERROR = I2 /
     820X,36HSECOND ORDER DISCRETIZATION ERROR = E12.5 /
     920X,36HFOURTH ORDER DISCRETIZATION ERROR = E12.5 /
     920X,29HREQUIRED LENGTH OF W ARRAY = I5)
C
      END
      SUBROUTINE COFX4(X,AF,BF,CF)
C
C     SET COEFFICIENTS IN THE X-DIRECTION.
C
      AF = (X+1.)**2
      BF = 2.0*(X+1.)
      CF = -X
      RETURN
      END
