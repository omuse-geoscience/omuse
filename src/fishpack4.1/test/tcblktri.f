c
c     file tblktri.f
c
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
C     THIS PROGRAM ILLUSTRATES THE USE OF SUBROUTINE CBLKTR WHICH IS
C     THE COMPLEX VERSION OF BLKTRI. THE PROGRAM SOLVES THE EQUATION
C
C     .5/S*(D/DS)(.5/S*DU/DS)+.5/T*(D/DT)(.5/T*DU/DT)-SQRT(-1)*U
C                                                          (1)
C            = 15/4*S*T*(S**4+T**4)-SQRT(-1)*(S*T)**5
C
C     ON THE RECTANGLE 0 .LT. S .LT. 1 AND 0 .LT. T .LT. 1
C     WITH THE BOUNDARY CONDITIONS
C
C     U(0,T) = 0
C                            0 .LE. T .LE. 1
C     U(1,T) = T**5
C
C     AND
C
C     U(S,0) = 0
C                            0 .LE. S .LE. 1
C     U(S,1) = S**5
C
C     THE EXACT SOLUTION OF THIS PROBLEM IS U(S,T) = (S*T)**5
C
C     DEFINE THE INTEGERS M = 50 AND N = 63. THEN DEFINE THE
C     GRID INCREMENTS DELTAS = 1/(M+1) AND DELTAT = 1/(N+1).
C
C     THE GRID IS THEN GIVEN BY S(I) = I*DELTAS FOR I = 1,...,M
C     AND T(J) = J*DELTAT FOR J = 1,...,N.
C
C     THE APPROXIMATE SOLUTION IS GIVEN AS THE SOLUTION TO
C     THE FOLLOWING FINITE DIFFERENCE APPROXIMATION OF EQUATION (1).
C
C     .5/(S(I)*DELTAS)*((U(I+1,J)-U(I,J))/(2*S(I+.5)*DELTAS)
C                     -(U(I,J)-U(I-1,J))/(2*S(I-.5)*DELTAS))
C     +.5/(T(I)*DELTAT)*((U(I,J+1)-U(I,J))/(2*T(I+.5)*DELTAT) (2)
C                     -(U(I,J)-U(I,J-1))/(2*T(I-.5)*DELTAT))
C                        -SQRT(-1)*U(I,J)
C               = 15/4*S(I)*T(J)*(S(I)**4+T(J)**4)
C                        -SQRT(-1)*(S(I)*T(J))**5
C
C             WHERE S(I+.5) = .5*(S(I+1)+S(I))
C                   S(I-.5) = .5*(S(I)+S(I-1))
C                   T(I+.5) = .5*(T(I+1)+T(I))
C                   T(I-.5) = .5*(T(I)+T(I-1))
C
C     THE APPROACH IS TO WRITE EQUATION (2) IN THE FORM
C
C     AM(I)*U(I-1,J)+BM(I,J)*U(I,J)+CM(I)*U(I+1,J)
C       +AN(J)*U(I,J-1)+BN(J)*U(I,J)+CN(J)*U(I,J+1)      (3)
C           = Y(I,J)
C
C     AND THEN CALL SUBROUTINE CBLKTR TO DETERMINE U(I,J)
C
C
C
C
      DIMENSION       Y(75,105)  ,AM(75)     ,BM(75)     ,CM(75)     ,
     1                AN(105)    ,BN(105)    ,CN(105)    ,S(75)      ,
     2                T(105)     ,W(1123)
      COMPLEX         Y          ,AM         ,BM         ,CM
C
      IFLG = 0
      NP = 1
      N = 63
      MP = 1
      M = 50
      IDIMY = 75
C
C     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING THE
C     COEFFICIENTS AND THE ARRAY Y.
C
      DELTAS = 1./FLOAT(M+1)
      DO 101 I=1,M
         S(I) = FLOAT(I)*DELTAS
  101 CONTINUE
      DELTAT = 1./FLOAT(N+1)
      DO 102 J=1,N
         T(J) = FLOAT(J)*DELTAT
  102 CONTINUE
C
C     COMPUTE THE COEFFICIENTS AM,BM,CM CORRESPONDING TO THE S DIRECTION
C
      HDS = DELTAS/2.
      TDS = DELTAS+DELTAS
      DO 103 I=1,M
         TEMP1 = 1./(S(I)*TDS)
         TEMP2 = 1./((S(I)-HDS)*TDS)
         TEMP3 = 1./((S(I)+HDS)*TDS)
         AM(I) = CMPLX(TEMP1*TEMP2,0.)
         CM(I) = CMPLX(TEMP1*TEMP3,0.)
         BM(I) = -(AM(I)+CM(I))-(0.,1.)
  103 CONTINUE
C
C     COMPUTE THE COEFFICIENTS AN,BN,CN CORRESPONDING TO THE T DIRECTION
C
      HDT = DELTAT/2.
      TDT = DELTAT+DELTAT
      DO 104 J=1,N
         TEMP1 = 1./(T(J)*TDT)
         TEMP2 = 1./((T(J)-HDT)*TDT)
         TEMP3 = 1./((T(J)+HDT)*TDT)
         AN(J) = TEMP1*TEMP2
         CN(J) = TEMP1*TEMP3
         BN(J) = -(AN(J)+CN(J))
  104 CONTINUE
C
C     COMPUTE RIGHT SIDE OF EQUATION
C
      DO 106 J=1,N
         DO 105 I=1,M
            Y(I,J) = 3.75*S(I)*T(J)*(S(I)**4+T(J)**4)-
     1               (0.,1.)*(S(I)*T(J))**5
  105    CONTINUE
  106 CONTINUE
C
C     THE NONZERO BOUNDARY CONDITIONS ENTER THE LINEAR SYSTEM VIA
C     THE RIGHT SIDE Y(I,J). IF THE EQUATIONS (3) GIVEN ABOVE
C     ARE EVALUATED AT I=M AND J=1,...,N THEN THE TERM CM(M)*U(M+1,J)
C     IS KNOWN FROM THE BOUNDARY CONDITION TO BE CM(M)*T(J)**5.
C     THEREFORE THIS TERM CAN BE INCLUDED IN THE RIGHT SIDE Y(M,J).
C     THE SAME ANALYSIS APPLIES AT J=N AND I=1,..,M. NOTE THAT THE
C     CORNER AT J=N,I=M INCLUDES CONTRIBUTIONS FROM BOTH BOUNDARIES.
C
      DO 107 J=1,N
         Y(M,J) = Y(M,J)-CM(M)*T(J)**5
  107 CONTINUE
      DO 108 I=1,M
         Y(I,N) = Y(I,N)-CN(N)*S(I)**5
  108 CONTINUE
C
  109 CALL CBLKTR (IFLG,NP,N,AN,BN,CN,MP,M,AM,BM,CM,IDIMY,Y,IERROR,W)
      IFLG = IFLG+1
      IF (IFLG-1) 109,109,110
C
C     COMPUTE DISCRETIZATION ERROR
C
  110 ERR = 0.
      DO 112 J=1,N
         DO 111 I=1,M
            Z = CABS(Y(I,J)-(S(I)*T(J))**5)
            IF (Z .GT. ERR) ERR = Z
  111    CONTINUE
  112 CONTINUE
      IW = INT(W(1))
      PRINT 1001 , IERROR,ERR,IW
      STOP
C
 1001 FORMAT (1H1,20X,25HSUBROUTINE CBLKTR EXAMPLE,//
     1        10X,"** important for TCBLKTRI example: to avoid  **",/
     2        10X,"** large discretization error, compile it    **",/
     3        10X,"** and FISHPACK routines in double precision **",///
     4        10X,46HTHE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS//
     5        32X,10HIERROR = 0/
     6        18X,34HDISCRETIZATION ERROR = 1.64572E-05/
     7        12X,33HREQUIRED LENGTH OF W ARRAY = 1123//
     8        10X,32HTHE OUTPUT FROM YOUR COMPUTER IS//
     9        32X,8HIERROR =,I2/18X,22HDISCRETIZATION ERROR =,E12.5/
     1        12X,28HREQUIRED LENGTH OF W ARRAY =,I5)
C
      END
