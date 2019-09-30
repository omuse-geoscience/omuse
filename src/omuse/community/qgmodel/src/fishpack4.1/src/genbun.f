C
C     file genbun.f
C
      SUBROUTINE GENBUN (NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,IERROR,W)
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
C
C
C DIMENSION OF           A(M),B(M),C(M),Y(IDIMY,N),
C                        W(SEE PARAMETER LIST)
C ARGUMENTS
C
C LATEST REVISION        NOVEMBER 1988
C
C PURPOSE                THE NAME OF THIS PACKAGE IS A MNEMONIC FOR THE
C                        GENERALIZED BUNEMAN ALGORITHM.
C
C                        IT SOLVES THE REAL LINEAR SYSTEM OF EQUATIONS
C
C                        A(I)*X(I-1,J) + B(I)*X(I,J) + C(I)*X(I+1,J)
C                        + X(I,J-1) - 2.*X(I,J) + X(I,J+1) = Y(I,J)
C
C                        FOR I = 1,2,...,M  AND  J = 1,2,...,N.
C
C                        INDICES I+1 AND I-1 ARE EVALUATED MODULO M,
C                        I.E., X(0,J) = X(M,J) AND X(M+1,J) = X(1,J),
C                        AND X(I,0) MAY EQUAL 0, X(I,2), OR X(I,N),
C                        AND X(I,N+1) MAY EQUAL 0, X(I,N-1), OR X(I,1)
C                        DEPENDING ON AN INPUT PARAMETER.
C
C USAGE                  CALL GENBUN (NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,
C                                     IERROR,W)
C
C ARGUMENTS
C
C ON INPUT               NPEROD
C
C                          INDICATES THE VALUES THAT X(I,0) AND
C                          X(I,N+1) ARE ASSUMED TO HAVE.
C
C                          = 0  IF X(I,0) = X(I,N) AND X(I,N+1) =
C                               X(I,1).
C                          = 1  IF X(I,0) = X(I,N+1) = 0  .
C                          = 2  IF X(I,0) = 0 AND X(I,N+1) = X(I,N-1).
C                          = 3  IF X(I,0) = X(I,2) AND X(I,N+1) =
C                               X(I,N-1).
C                          = 4  IF X(I,0) = X(I,2) AND X(I,N+1) = 0.
C
C                        N
C                          THE NUMBER OF UNKNOWNS IN THE J-DIRECTION.
C                          N MUST BE GREATER THAN 2.
C
C                        MPEROD
C                          = 0 IF A(1) AND C(M) ARE NOT ZERO
C                          = 1 IF A(1) = C(M) = 0
C
C                        M
C                          THE NUMBER OF UNKNOWNS IN THE I-DIRECTION.
C                          N MUST BE GREATER THAN 2.
C
C                        A,B,C
C                          ONE-DIMENSIONAL ARRAYS OF LENGTH M THAT
C                          SPECIFY THE COEFFICIENTS IN THE LINEAR
C                          EQUATIONS GIVEN ABOVE.  IF MPEROD = 0
C                          THE ARRAY ELEMENTS MUST NOT DEPEND UPON
C                          THE INDEX I, BUT MUST BE CONSTANT.
C                          SPECIFICALLY, THE SUBROUTINE CHECKS THE
C                          FOLLOWING CONDITION .
C
C                            A(I) = C(1)
C                            C(I) = C(1)
C                            B(I) = B(1)
C
C                          FOR I=1,2,...,M.
C
C                        IDIMY
C                          THE ROW (OR FIRST) DIMENSION OF THE
C                          TWO-DIMENSIONAL ARRAY Y AS IT APPEARS
C                          IN THE PROGRAM CALLING GENBUN.
C                          THIS PARAMETER IS USED TO SPECIFY THE
C                          VARIABLE DIMENSION OF Y.
C                          IDIMY MUST BE AT LEAST M.
C
C                        Y
C                          A TWO-DIMENSIONAL COMPLEX ARRAY THAT
C                          SPECIFIES THE VALUES OF THE RIGHT SIDE
C                          OF THE LINEAR SYSTEM OF EQUATIONS GIVEN
C                          ABOVE.
C                          Y MUST BE DIMENSIONED AT LEAST M*N.
C
C                        W
C                          A ONE-DIMENSIONAL ARRAY THAT MUST
C                          BE PROVIDED BY THE USER FOR WORK
C                          SPACE.  W MAY REQUIRE UP TO 4*N +
C                          (10 + INT(LOG2(N)))*M LOCATIONS.
C                          THE ACTUAL NUMBER OF LOCATIONS USED IS
C                          COMPUTED BY GENBUN AND IS RETURNED IN
C                          LOCATION W(1).
C
C
C  ON OUTPUT             Y
C
C                          CONTAINS THE SOLUTION X.
C
C                        IERROR
C                          AN ERROR FLAG WHICH INDICATES INVALID
C                          INPUT PARAMETERS  EXCEPT FOR NUMBER
C                          ZERO, A SOLUTION IS NOT ATTEMPTED.
C
C                          = 0  NO ERROR.
C                          = 1  M .LE. 2  .
C                          = 2  N .LE. 2
C                          = 3  IDIMY .LT. M
C                          = 4  NPEROD .LT. 0 OR NPEROD .GT. 4
C                          = 5  MPEROD .LT. 0 OR MPEROD .GT. 1
C                          = 6  A(I) .NE. C(1) OR C(I) .NE. C(1) OR
C                               B(I) .NE. B(1) FOR
C                               SOME I=1,2,...,M.
C                          = 7  A(1) .NE. 0 OR C(M) .NE. 0 AND
C                                 MPEROD = 1
C
C                        W
C                          W(1) CONTAINS THE REQUIRED LENGTH OF W.
C
C SPECIAL CONDITONS      NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       COMF AND GNBNAUX FROM FISHPACK
C FILES
C
C LANGUAGE               FORTRAN
C
C HISTORY                WRITTEN IN 1979 BY ROLAND SWEET OF NCAR'S
C                        SCIENTIFIC COMPUTING DIVISION.  MADE AVAILABLE
C                        ON NCAR'S PUBLIC LIBRARIES IN JANUARY, 1980.
C
C ALGORITHM              THE LINEAR SYSTEM IS SOLVED BY A CYCLIC
C                        REDUCTION ALGORITHM DESCRIBED IN THE
C                        REFERENCE.
C
C PORTABILITY            FORTRAN 77 --
C                        THE MACHINE DEPENDENT CONSTANT PI IS
C                        DEFINED IN FUNCTION PIMACH.
C
C REFERENCES             SWEET, R., "A CYCLIC REDUCTION ALGORITHM FOR
C                        SOLVING BLOCK TRIDIAGONAL SYSTEMS OF ARBITRARY
C                        DIMENSIONS," SIAM J. ON NUMER. ANAL., 14 (1977)
C                        PP. 706-720.
C
C ACCURACY               THIS TEST WAS PERFORMED ON A CDC 7600:
C
C                        A UNIFORM RANDOM NUMBER GENERATOR WAS USED
C                        TO CREATE A SOLUTION ARRAY X FOR THE SYSTEM
C                        GIVEN IN THE 'PURPOSE' DESCRIPTION ABOVE
C                        WITH
C                          A(I) = C(I) = -0.5*B(I) = 1, I=1,2,...,M
C
C                        AND, WHEN MPEROD = 1
C
C                          A(1) = C(M) = 0
C                          A(M) = C(1) = 2.
C
C                        THE SOLUTION X WAS SUBSTITUTED INTO THE
C                        GIVEN SYSTEM  AND, USING DOUBLE PRECISION
C                        A RIGHT SIDE Y WAS COMPUTED.
C                        USING THIS ARRAY Y, SUBROUTINE GENBUN
C                        WAS CALLED TO PRODUCE APPROXIMATE
C                        SOLUTION Z.  THEN RELATIVE ERROR
C                          E = MAX(ABS(Z(I,J)-X(I,J)))/
C                              MAX(ABS(X(I,J)))
C                        WAS COMPUTED, WHERE THE TWO MAXIMA ARE TAKEN
C                        OVER I=1,2,...,M AND J=1,...,N.
C
C                        THE VALUE OF E IS GIVEN IN THE TABLE
C                        BELOW FOR SOME TYPICAL VALUES OF M AND N.
C
C                   M (=N)    MPEROD    NPEROD    T(MSECS)    E
C                   ------    ------    ------    --------  ------
C
C                     31        0         0          36     6.E-14
C                     31        1         1          21     4.E-13
C                     31        1         3          41     3.E-13
C                     32        0         0          29     9.E-14
C                     32        1         1          32     3.E-13
C                     32        1         3          48     1.E-13
C                     33        0         0          36     9.E-14
C                     33        1         1          30     4.E-13
C                     33        1         3          34     1.E-13
C                     63        0         0         150     1.E-13
C                     63        1         1          91     1.E-12
C                     63        1         3         173     2.E-13
C                     64        0         0         122     1.E-13
C                     64        1         1         128     1.E-12
C                     64        1         3         199     6.E-13
C                     65        0         0         143     2.E-13
C                     65        1         1         120     1.E-12
C                     65        1         3         138     4.E-13
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DIMENSION       Y(IDIMY,1)
      DIMENSION       W(*)       ,B(*)       ,A(*)       ,C(*)
C
      IERROR = 0
      IF (M .LE. 2) IERROR = 1
      IF (N .LE. 2) IERROR = 2
      IF (IDIMY .LT. M) IERROR = 3
      IF (NPEROD.LT.0 .OR. NPEROD.GT.4) IERROR = 4
      IF (MPEROD.LT.0 .OR. MPEROD.GT.1) IERROR = 5
      IF (MPEROD .EQ. 1) GO TO 102
      DO 101 I=2,M
         IF (A(I) .NE. C(1)) GO TO 103
         IF (C(I) .NE. C(1)) GO TO 103
         IF (B(I) .NE. B(1)) GO TO 103
  101 CONTINUE
      GO TO 104
  102 IF (A(1).NE.0. .OR. C(M).NE.0.) IERROR = 7
      GO TO 104
  103 IERROR = 6
  104 IF (IERROR .NE. 0) RETURN
      MP1 = M+1
      IWBA = MP1
      IWBB = IWBA+M
      IWBC = IWBB+M
      IWB2 = IWBC+M
      IWB3 = IWB2+M
      IWW1 = IWB3+M
      IWW2 = IWW1+M
      IWW3 = IWW2+M
      IWD = IWW3+M
      IWTCOS = IWD+M
      IWP = IWTCOS+4*N
      DO 106 I=1,M
         K = IWBA+I-1
         W(K) = -A(I)
         K = IWBC+I-1
         W(K) = -C(I)
         K = IWBB+I-1
         W(K) = 2.-B(I)
         DO 105 J=1,N
            Y(I,J) = -Y(I,J)
  105    CONTINUE
  106 CONTINUE
      MP = MPEROD+1
      NP = NPEROD+1
      GO TO (114,107),MP
  107 GO TO (108,109,110,111,123),NP
  108 CALL POISP2 (M,N,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWB2),
     1             W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS),
     2             W(IWP))
      GO TO 112
  109 CALL POISD2 (M,N,1,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWW1),
     1             W(IWD),W(IWTCOS),W(IWP))
      GO TO 112
  110 CALL POISN2 (M,N,1,2,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWB2),
     1             W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS),
     2             W(IWP))
      GO TO 112
  111 CALL POISN2 (M,N,1,1,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWB2),
     1             W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS),
     2             W(IWP))
  112 IPSTOR = W(IWW1)
      IREV = 2
      IF (NPEROD .EQ. 4) GO TO 124
  113 GO TO (127,133),MP
  114 CONTINUE
C
C     REORDER UNKNOWNS WHEN MP =0
C
      MH = (M+1)/2
      MHM1 = MH-1
      MODD = 1
      IF (MH*2 .EQ. M) MODD = 2
      DO 119 J=1,N
         DO 115 I=1,MHM1
            MHPI = MH+I
            MHMI = MH-I
            W(I) = Y(MHMI,J)-Y(MHPI,J)
            W(MHPI) = Y(MHMI,J)+Y(MHPI,J)
  115    CONTINUE
         W(MH) = 2.*Y(MH,J)
         GO TO (117,116),MODD
  116    W(M) = 2.*Y(M,J)
  117    CONTINUE
         DO 118 I=1,M
            Y(I,J) = W(I)
  118    CONTINUE
  119 CONTINUE
      K = IWBC+MHM1-1
      I = IWBA+MHM1
      W(K) = 0.
      W(I) = 0.
      W(K+1) = 2.*W(K+1)
      GO TO (120,121),MODD
  120 CONTINUE
      K = IWBB+MHM1-1
      W(K) = W(K)-W(I-1)
      W(IWBC-1) = W(IWBC-1)+W(IWBB-1)
      GO TO 122
  121 W(IWBB-1) = W(K+1)
  122 CONTINUE
      GO TO 107
C
C     REVERSE COLUMNS WHEN NPEROD = 4.
C
  123 IREV = 1
      NBY2 = N/2
  124 DO 126 J=1,NBY2
         MSKIP = N+1-J
         DO 125 I=1,M
            A1 = Y(I,J)
            Y(I,J) = Y(I,MSKIP)
            Y(I,MSKIP) = A1
  125    CONTINUE
  126 CONTINUE
      GO TO (110,113),IREV
  127 CONTINUE
      DO 132 J=1,N
         DO 128 I=1,MHM1
            MHMI = MH-I
            MHPI = MH+I
            W(MHMI) = .5*(Y(MHPI,J)+Y(I,J))
            W(MHPI) = .5*(Y(MHPI,J)-Y(I,J))
  128    CONTINUE
         W(MH) = .5*Y(MH,J)
         GO TO (130,129),MODD
  129    W(M) = .5*Y(M,J)
  130    CONTINUE
         DO 131 I=1,M
            Y(I,J) = W(I)
  131    CONTINUE
  132 CONTINUE
  133 CONTINUE
C
C     RETURN STORAGE REQUIREMENTS FOR W ARRAY.
C
      W(1) = IPSTOR+IWP-1
      RETURN
      END
      SUBROUTINE POISD2 (MR,NR,ISTAG,BA,BB,BC,Q,IDIMQ,B,W,D,TCOS,P)
C
C     SUBROUTINE TO SOLVE POISSON'S EQUATION FOR DIRICHLET BOUNDARY
C     CONDITIONS.
C
C     ISTAG = 1 IF THE LAST DIAGONAL BLOCK IS THE MATRIX A.
C     ISTAG = 2 IF THE LAST DIAGONAL BLOCK IS THE MATRIX A+I.
C
      DIMENSION       Q(IDIMQ,1) ,BA(*)      ,BB(*)      ,BC(*)      ,
     1                TCOS(*)    ,B(*)       ,D(*)       ,W(*)       ,
     2                P(*)
      M = MR
      N = NR
      JSH = 0
      FI = 1./FLOAT(ISTAG)
      IP = -M
      IPSTOR = 0
      GO TO (101,102),ISTAG
  101 KR = 0
      IRREG = 1
      IF (N .GT. 1) GO TO 106
      TCOS(1) = 0.
      GO TO 103
  102 KR = 1
      JSTSAV = 1
      IRREG = 2
      IF (N .GT. 1) GO TO 106
      TCOS(1) = -1.
  103 DO 104 I=1,M
         B(I) = Q(I,1)
  104 CONTINUE
      CALL TRIX (1,0,M,BA,BB,BC,B,TCOS,D,W)
      DO 105 I=1,M
         Q(I,1) = B(I)
  105 CONTINUE
      GO TO 183
  106 LR = 0
      DO 107 I=1,M
         P(I) = 0.
  107 CONTINUE
      NUN = N
      JST = 1
      JSP = N
C
C     IRREG = 1 WHEN NO IRREGULARITIES HAVE OCCURRED, OTHERWISE IT IS 2.
C
  108 L = 2*JST
      NODD = 2-2*((NUN+1)/2)+NUN
C
C     NODD = 1 WHEN NUN IS ODD, OTHERWISE IT IS 2.
C
      GO TO (110,109),NODD
  109 JSP = JSP-L
      GO TO 111
  110 JSP = JSP-JST
      IF (IRREG .NE. 1) JSP = JSP-L
  111 CONTINUE
C
C     REGULAR REDUCTION
C
      CALL COSGEN (JST,1,0.5,0.0,TCOS)
      IF (L .GT. JSP) GO TO 118
      DO 117 J=L,JSP,L
         JM1 = J-JSH
         JP1 = J+JSH
         JM2 = J-JST
         JP2 = J+JST
         JM3 = JM2-JSH
         JP3 = JP2+JSH
         IF (JST .NE. 1) GO TO 113
         DO 112 I=1,M
            B(I) = 2.*Q(I,J)
            Q(I,J) = Q(I,JM2)+Q(I,JP2)
  112    CONTINUE
         GO TO 115
  113    DO 114 I=1,M
            T = Q(I,J)-Q(I,JM1)-Q(I,JP1)+Q(I,JM2)+Q(I,JP2)
            B(I) = T+Q(I,J)-Q(I,JM3)-Q(I,JP3)
            Q(I,J) = T
  114    CONTINUE
  115    CONTINUE
         CALL TRIX (JST,0,M,BA,BB,BC,B,TCOS,D,W)
         DO 116 I=1,M
            Q(I,J) = Q(I,J)+B(I)
  116    CONTINUE
  117 CONTINUE
C
C     REDUCTION FOR LAST UNKNOWN
C
  118 GO TO (119,136),NODD
  119 GO TO (152,120),IRREG
C
C     ODD NUMBER OF UNKNOWNS
C
  120 JSP = JSP+L
      J = JSP
      JM1 = J-JSH
      JP1 = J+JSH
      JM2 = J-JST
      JP2 = J+JST
      JM3 = JM2-JSH
      GO TO (123,121),ISTAG
  121 CONTINUE
      IF (JST .NE. 1) GO TO 123
      DO 122 I=1,M
         B(I) = Q(I,J)
         Q(I,J) = 0.
  122 CONTINUE
      GO TO 130
  123 GO TO (124,126),NODDPR
  124 DO 125 I=1,M
         IP1 = IP+I
         B(I) = .5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))+P(IP1)+Q(I,J)
  125 CONTINUE
      GO TO 128
  126 DO 127 I=1,M
         B(I) = .5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))+Q(I,JP2)-Q(I,JP1)+Q(I,J)
  127 CONTINUE
  128 DO 129 I=1,M
         Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
  129 CONTINUE
  130 CALL TRIX (JST,0,M,BA,BB,BC,B,TCOS,D,W)
      IP = IP+M
      IPSTOR = MAX0(IPSTOR,IP+M)
      DO 131 I=1,M
         IP1 = IP+I
         P(IP1) = Q(I,J)+B(I)
         B(I) = Q(I,JP2)+P(IP1)
  131 CONTINUE
      IF (LR .NE. 0) GO TO 133
      DO 132 I=1,JST
         KRPI = KR+I
         TCOS(KRPI) = TCOS(I)
  132 CONTINUE
      GO TO 134
  133 CONTINUE
      CALL COSGEN (LR,JSTSAV,0.,FI,TCOS(JST+1))
      CALL MERGE (TCOS,0,JST,JST,LR,KR)
  134 CONTINUE
      CALL COSGEN (KR,JSTSAV,0.0,FI,TCOS)
      CALL TRIX (KR,KR,M,BA,BB,BC,B,TCOS,D,W)
      DO 135 I=1,M
         IP1 = IP+I
         Q(I,J) = Q(I,JM2)+B(I)+P(IP1)
  135 CONTINUE
      LR = KR
      KR = KR+L
      GO TO 152
C
C     EVEN NUMBER OF UNKNOWNS
C
  136 JSP = JSP+L
      J = JSP
      JM1 = J-JSH
      JP1 = J+JSH
      JM2 = J-JST
      JP2 = J+JST
      JM3 = JM2-JSH
      GO TO (137,138),IRREG
  137 CONTINUE
      JSTSAV = JST
      IDEG = JST
      KR = L
      GO TO 139
  138 CALL COSGEN (KR,JSTSAV,0.0,FI,TCOS)
      CALL COSGEN (LR,JSTSAV,0.0,FI,TCOS(KR+1))
      IDEG = KR
      KR = KR+JST
  139 IF (JST .NE. 1) GO TO 141
      IRREG = 2
      DO 140 I=1,M
         B(I) = Q(I,J)
         Q(I,J) = Q(I,JM2)
  140 CONTINUE
      GO TO 150
  141 DO 142 I=1,M
         B(I) = Q(I,J)+.5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
  142 CONTINUE
      GO TO (143,145),IRREG
  143 DO 144 I=1,M
         Q(I,J) = Q(I,JM2)+.5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
  144 CONTINUE
      IRREG = 2
      GO TO 150
  145 CONTINUE
      GO TO (146,148),NODDPR
  146 DO 147 I=1,M
         IP1 = IP+I
         Q(I,J) = Q(I,JM2)+P(IP1)
  147 CONTINUE
      IP = IP-M
      GO TO 150
  148 DO 149 I=1,M
         Q(I,J) = Q(I,JM2)+Q(I,J)-Q(I,JM1)
  149 CONTINUE
  150 CALL TRIX (IDEG,LR,M,BA,BB,BC,B,TCOS,D,W)
      DO 151 I=1,M
         Q(I,J) = Q(I,J)+B(I)
  151 CONTINUE
  152 NUN = NUN/2
      NODDPR = NODD
      JSH = JST
      JST = 2*JST
      IF (NUN .GE. 2) GO TO 108
C
C     START SOLUTION.
C
      J = JSP
      DO 153 I=1,M
         B(I) = Q(I,J)
  153 CONTINUE
      GO TO (154,155),IRREG
  154 CONTINUE
      CALL COSGEN (JST,1,0.5,0.0,TCOS)
      IDEG = JST
      GO TO 156
  155 KR = LR+JST
      CALL COSGEN (KR,JSTSAV,0.0,FI,TCOS)
      CALL COSGEN (LR,JSTSAV,0.0,FI,TCOS(KR+1))
      IDEG = KR
  156 CONTINUE
      CALL TRIX (IDEG,LR,M,BA,BB,BC,B,TCOS,D,W)
      JM1 = J-JSH
      JP1 = J+JSH
      GO TO (157,159),IRREG
  157 DO 158 I=1,M
         Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))+B(I)
  158 CONTINUE
      GO TO 164
  159 GO TO (160,162),NODDPR
  160 DO 161 I=1,M
         IP1 = IP+I
         Q(I,J) = P(IP1)+B(I)
  161 CONTINUE
      IP = IP-M
      GO TO 164
  162 DO 163 I=1,M
         Q(I,J) = Q(I,J)-Q(I,JM1)+B(I)
  163 CONTINUE
  164 CONTINUE
C
C     START BACK SUBSTITUTION.
C
      JST = JST/2
      JSH = JST/2
      NUN = 2*NUN
      IF (NUN .GT. N) GO TO 183
      DO 182 J=JST,N,L
         JM1 = J-JSH
         JP1 = J+JSH
         JM2 = J-JST
         JP2 = J+JST
         IF (J .GT. JST) GO TO 166
         DO 165 I=1,M
            B(I) = Q(I,J)+Q(I,JP2)
  165    CONTINUE
         GO TO 170
  166    IF (JP2 .LE. N) GO TO 168
         DO 167 I=1,M
            B(I) = Q(I,J)+Q(I,JM2)
  167    CONTINUE
         IF (JST .LT. JSTSAV) IRREG = 1
         GO TO (170,171),IRREG
  168    DO 169 I=1,M
            B(I) = Q(I,J)+Q(I,JM2)+Q(I,JP2)
  169    CONTINUE
  170    CONTINUE
         CALL COSGEN (JST,1,0.5,0.0,TCOS)
         IDEG = JST
         JDEG = 0
         GO TO 172
  171    IF (J+L .GT. N) LR = LR-JST
         KR = JST+LR
         CALL COSGEN (KR,JSTSAV,0.0,FI,TCOS)
         CALL COSGEN (LR,JSTSAV,0.0,FI,TCOS(KR+1))
         IDEG = KR
         JDEG = LR
  172    CONTINUE
         CALL TRIX (IDEG,JDEG,M,BA,BB,BC,B,TCOS,D,W)
         IF (JST .GT. 1) GO TO 174
         DO 173 I=1,M
            Q(I,J) = B(I)
  173    CONTINUE
         GO TO 182
  174    IF (JP2 .GT. N) GO TO 177
  175    DO 176 I=1,M
            Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))+B(I)
  176    CONTINUE
         GO TO 182
  177    GO TO (175,178),IRREG
  178    IF (J+JSH .GT. N) GO TO 180
         DO 179 I=1,M
            IP1 = IP+I
            Q(I,J) = B(I)+P(IP1)
  179    CONTINUE
         IP = IP-M
         GO TO 182
  180    DO 181 I=1,M
            Q(I,J) = B(I)+Q(I,J)-Q(I,JM1)
  181    CONTINUE
  182 CONTINUE
      L = L/2
      GO TO 164
  183 CONTINUE
C
C     RETURN STORAGE REQUIREMENTS FOR P VECTORS.
C
      W(1) = IPSTOR
      RETURN
      END
      SUBROUTINE POISN2 (M,N,ISTAG,MIXBND,A,BB,C,Q,IDIMQ,B,B2,B3,W,W2,
     1                   W3,D,TCOS,P)
C
C     SUBROUTINE TO SOLVE POISSON'S EQUATION WITH NEUMANN BOUNDARY
C     CONDITIONS.
C
C     ISTAG = 1 IF THE LAST DIAGONAL BLOCK IS A.
C     ISTAG = 2 IF THE LAST DIAGONAL BLOCK IS A-I.
C     MIXBND = 1 IF HAVE NEUMANN BOUNDARY CONDITIONS AT BOTH BOUNDARIES.
C     MIXBND = 2 IF HAVE NEUMANN BOUNDARY CONDITIONS AT BOTTOM AND
C     DIRICHLET CONDITION AT TOP.  (FOR THIS CASE, MUST HAVE ISTAG = 1.)
C
      DIMENSION       A(*)       ,BB(*)      ,C(*)       ,Q(IDIMQ,*) ,
     1                B(*)       ,B2(*)      ,B3(*)      ,W(*)       ,
     2                W2(*)      ,W3(*)      ,D(*)       ,TCOS(*)    ,
     3                K(4)       ,P(*)
      EQUIVALENCE     (K(1),K1)  ,(K(2),K2)  ,(K(3),K3)  ,(K(4),K4)
      FISTAG = 3-ISTAG
      FNUM = 1./FLOAT(ISTAG)
      FDEN = 0.5*FLOAT(ISTAG-1)
      MR = M
      IP = -MR
      IPSTOR = 0
      I2R = 1
      JR = 2
      NR = N
      NLAST = N
      KR = 1
      LR = 0
      GO TO (101,103),ISTAG
  101 CONTINUE
      DO 102 I=1,MR
         Q(I,N) = .5*Q(I,N)
  102 CONTINUE
      GO TO (103,104),MIXBND
  103 IF (N .LE. 3) GO TO 155
  104 CONTINUE
      JR = 2*I2R
      NROD = 1
      IF ((NR/2)*2 .EQ. NR) NROD = 0
      GO TO (105,106),MIXBND
  105 JSTART = 1
      GO TO 107
  106 JSTART = JR
      NROD = 1-NROD
  107 CONTINUE
      JSTOP = NLAST-JR
      IF (NROD .EQ. 0) JSTOP = JSTOP-I2R
      CALL COSGEN (I2R,1,0.5,0.0,TCOS)
      I2RBY2 = I2R/2
      IF (JSTOP .GE. JSTART) GO TO 108
      J = JR
      GO TO 116
  108 CONTINUE
C
C     REGULAR REDUCTION.
C
      DO 115 J=JSTART,JSTOP,JR
         JP1 = J+I2RBY2
         JP2 = J+I2R
         JP3 = JP2+I2RBY2
         JM1 = J-I2RBY2
         JM2 = J-I2R
         JM3 = JM2-I2RBY2
         IF (J .NE. 1) GO TO 109
         JM1 = JP1
         JM2 = JP2
         JM3 = JP3
  109    CONTINUE
         IF (I2R .NE. 1) GO TO 111
         IF (J .EQ. 1) JM2 = JP2
         DO 110 I=1,MR
            B(I) = 2.*Q(I,J)
            Q(I,J) = Q(I,JM2)+Q(I,JP2)
  110    CONTINUE
         GO TO 113
  111    CONTINUE
         DO 112 I=1,MR
            FI = Q(I,J)
            Q(I,J) = Q(I,J)-Q(I,JM1)-Q(I,JP1)+Q(I,JM2)+Q(I,JP2)
            B(I) = FI+Q(I,J)-Q(I,JM3)-Q(I,JP3)
  112    CONTINUE
  113    CONTINUE
         CALL TRIX (I2R,0,MR,A,BB,C,B,TCOS,D,W)
         DO 114 I=1,MR
            Q(I,J) = Q(I,J)+B(I)
  114    CONTINUE
C
C     END OF REDUCTION FOR REGULAR UNKNOWNS.
C
  115 CONTINUE
C
C     BEGIN SPECIAL REDUCTION FOR LAST UNKNOWN.
C
      J = JSTOP+JR
  116 NLAST = J
      JM1 = J-I2RBY2
      JM2 = J-I2R
      JM3 = JM2-I2RBY2
      IF (NROD .EQ. 0) GO TO 128
C
C     ODD NUMBER OF UNKNOWNS
C
      IF (I2R .NE. 1) GO TO 118
      DO 117 I=1,MR
         B(I) = FISTAG*Q(I,J)
         Q(I,J) = Q(I,JM2)
  117 CONTINUE
      GO TO 126
  118 DO 119 I=1,MR
         B(I) = Q(I,J)+.5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
  119 CONTINUE
      IF (NRODPR .NE. 0) GO TO 121
      DO 120 I=1,MR
         II = IP+I
         Q(I,J) = Q(I,JM2)+P(II)
  120 CONTINUE
      IP = IP-MR
      GO TO 123
  121 CONTINUE
      DO 122 I=1,MR
         Q(I,J) = Q(I,J)-Q(I,JM1)+Q(I,JM2)
  122 CONTINUE
  123 IF (LR .EQ. 0) GO TO 124
      CALL COSGEN (LR,1,0.5,FDEN,TCOS(KR+1))
      GO TO 126
  124 CONTINUE
      DO 125 I=1,MR
         B(I) = FISTAG*B(I)
  125 CONTINUE
  126 CONTINUE
      CALL COSGEN (KR,1,0.5,FDEN,TCOS)
      CALL TRIX (KR,LR,MR,A,BB,C,B,TCOS,D,W)
      DO 127 I=1,MR
         Q(I,J) = Q(I,J)+B(I)
  127 CONTINUE
      KR = KR+I2R
      GO TO 151
  128 CONTINUE
C
C     EVEN NUMBER OF UNKNOWNS
C
      JP1 = J+I2RBY2
      JP2 = J+I2R
      IF (I2R .NE. 1) GO TO 135
      DO 129 I=1,MR
         B(I) = Q(I,J)
  129 CONTINUE
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      IP = 0
      IPSTOR = MR
      GO TO (133,130),ISTAG
  130 DO 131 I=1,MR
         P(I) = B(I)
         B(I) = B(I)+Q(I,N)
  131 CONTINUE
      TCOS(1) = 1.
      TCOS(2) = 0.
      CALL TRIX (1,1,MR,A,BB,C,B,TCOS,D,W)
      DO 132 I=1,MR
         Q(I,J) = Q(I,JM2)+P(I)+B(I)
  132 CONTINUE
      GO TO 150
  133 CONTINUE
      DO 134 I=1,MR
         P(I) = B(I)
         Q(I,J) = Q(I,JM2)+2.*Q(I,JP2)+3.*B(I)
  134 CONTINUE
      GO TO 150
  135 CONTINUE
      DO 136 I=1,MR
         B(I) = Q(I,J)+.5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
  136 CONTINUE
      IF (NRODPR .NE. 0) GO TO 138
      DO 137 I=1,MR
         II = IP+I
         B(I) = B(I)+P(II)
  137 CONTINUE
      GO TO 140
  138 CONTINUE
      DO 139 I=1,MR
         B(I) = B(I)+Q(I,JP2)-Q(I,JP1)
  139 CONTINUE
  140 CONTINUE
      CALL TRIX (I2R,0,MR,A,BB,C,B,TCOS,D,W)
      IP = IP+MR
      IPSTOR = MAX0(IPSTOR,IP+MR)
      DO 141 I=1,MR
         II = IP+I
         P(II) = B(I)+.5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
         B(I) = P(II)+Q(I,JP2)
  141 CONTINUE
      IF (LR .EQ. 0) GO TO 142
      CALL COSGEN (LR,1,0.5,FDEN,TCOS(I2R+1))
      CALL MERGE (TCOS,0,I2R,I2R,LR,KR)
      GO TO 144
  142 DO 143 I=1,I2R
         II = KR+I
         TCOS(II) = TCOS(I)
  143 CONTINUE
  144 CALL COSGEN (KR,1,0.5,FDEN,TCOS)
      IF (LR .NE. 0) GO TO 145
      GO TO (146,145),ISTAG
  145 CONTINUE
      CALL TRIX (KR,KR,MR,A,BB,C,B,TCOS,D,W)
      GO TO 148
  146 CONTINUE
      DO 147 I=1,MR
         B(I) = FISTAG*B(I)
  147 CONTINUE
  148 CONTINUE
      DO 149 I=1,MR
         II = IP+I
         Q(I,J) = Q(I,JM2)+P(II)+B(I)
  149 CONTINUE
  150 CONTINUE
      LR = KR
      KR = KR+JR
  151 CONTINUE
      GO TO (152,153),MIXBND
  152 NR = (NLAST-1)/JR+1
      IF (NR .LE. 3) GO TO 155
      GO TO 154
  153 NR = NLAST/JR
      IF (NR .LE. 1) GO TO 192
  154 I2R = JR
      NRODPR = NROD
      GO TO 104
  155 CONTINUE
C
C      BEGIN SOLUTION
C
      J = 1+JR
      JM1 = J-I2R
      JP1 = J+I2R
      JM2 = NLAST-I2R
      IF (NR .EQ. 2) GO TO 184
      IF (LR .NE. 0) GO TO 170
      IF (N .NE. 3) GO TO 161
C
C     CASE N = 3.
C
      GO TO (156,168),ISTAG
  156 CONTINUE
      DO 157 I=1,MR
         B(I) = Q(I,2)
  157 CONTINUE
      TCOS(1) = 0.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 158 I=1,MR
         Q(I,2) = B(I)
         B(I) = 4.*B(I)+Q(I,1)+2.*Q(I,3)
  158 CONTINUE
      TCOS(1) = -2.
      TCOS(2) = 2.
      I1 = 2
      I2 = 0
      CALL TRIX (I1,I2,MR,A,BB,C,B,TCOS,D,W)
      DO 159 I=1,MR
         Q(I,2) = Q(I,2)+B(I)
         B(I) = Q(I,1)+2.*Q(I,2)
  159 CONTINUE
      TCOS(1) = 0.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 160 I=1,MR
         Q(I,1) = B(I)
  160 CONTINUE
      JR = 1
      I2R = 0
      GO TO 194
C
C     CASE N = 2**P+1
C
  161 CONTINUE
      GO TO (162,170),ISTAG
  162 CONTINUE
      DO 163 I=1,MR
         B(I) = Q(I,J)+.5*Q(I,1)-Q(I,JM1)+Q(I,NLAST)-Q(I,JM2)
  163 CONTINUE
      CALL COSGEN (JR,1,0.5,0.0,TCOS)
      CALL TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
      DO 164 I=1,MR
         Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))+B(I)
         B(I) = Q(I,1)+2.*Q(I,NLAST)+4.*Q(I,J)
  164 CONTINUE
      JR2 = 2*JR
      CALL COSGEN (JR,1,0.0,0.0,TCOS)
      DO 165 I=1,JR
         I1 = JR+I
         I2 = JR+1-I
         TCOS(I1) = -TCOS(I2)
  165 CONTINUE
      CALL TRIX (JR2,0,MR,A,BB,C,B,TCOS,D,W)
      DO 166 I=1,MR
         Q(I,J) = Q(I,J)+B(I)
         B(I) = Q(I,1)+2.*Q(I,J)
  166 CONTINUE
      CALL COSGEN (JR,1,0.5,0.0,TCOS)
      CALL TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
      DO 167 I=1,MR
         Q(I,1) = .5*Q(I,1)-Q(I,JM1)+B(I)
  167 CONTINUE
      GO TO 194
C
C     CASE OF GENERAL N WITH NR = 3 .
C
  168 DO 169 I=1,MR
         B(I) = Q(I,2)
         Q(I,2) = 0.
         B2(I) = Q(I,3)
         B3(I) = Q(I,1)
  169 CONTINUE
      JR = 1
      I2R = 0
      J = 2
      GO TO 177
  170 CONTINUE
      DO 171 I=1,MR
         B(I) = .5*Q(I,1)-Q(I,JM1)+Q(I,J)
  171 CONTINUE
      IF (NROD .NE. 0) GO TO 173
      DO 172 I=1,MR
         II = IP+I
         B(I) = B(I)+P(II)
  172 CONTINUE
      GO TO 175
  173 DO 174 I=1,MR
         B(I) = B(I)+Q(I,NLAST)-Q(I,JM2)
  174 CONTINUE
  175 CONTINUE
      DO 176 I=1,MR
         T = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
         Q(I,J) = T
         B2(I) = Q(I,NLAST)+T
         B3(I) = Q(I,1)+2.*T
  176 CONTINUE
  177 CONTINUE
      K1 = KR+2*JR-1
      K2 = KR+JR
      TCOS(K1+1) = -2.
      K4 = K1+3-ISTAG
      CALL COSGEN (K2+ISTAG-2,1,0.0,FNUM,TCOS(K4))
      K4 = K1+K2+1
      CALL COSGEN (JR-1,1,0.0,1.0,TCOS(K4))
      CALL MERGE (TCOS,K1,K2,K1+K2,JR-1,0)
      K3 = K1+K2+LR
      CALL COSGEN (JR,1,0.5,0.0,TCOS(K3+1))
      K4 = K3+JR+1
      CALL COSGEN (KR,1,0.5,FDEN,TCOS(K4))
      CALL MERGE (TCOS,K3,JR,K3+JR,KR,K1)
      IF (LR .EQ. 0) GO TO 178
      CALL COSGEN (LR,1,0.5,FDEN,TCOS(K4))
      CALL MERGE (TCOS,K3,JR,K3+JR,LR,K3-LR)
      CALL COSGEN (KR,1,0.5,FDEN,TCOS(K4))
  178 K3 = KR
      K4 = KR
      CALL TRI3 (MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
      DO 179 I=1,MR
         B(I) = B(I)+B2(I)+B3(I)
  179 CONTINUE
      TCOS(1) = 2.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 180 I=1,MR
         Q(I,J) = Q(I,J)+B(I)
         B(I) = Q(I,1)+2.*Q(I,J)
  180 CONTINUE
      CALL COSGEN (JR,1,0.5,0.0,TCOS)
      CALL TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
      IF (JR .NE. 1) GO TO 182
      DO 181 I=1,MR
         Q(I,1) = B(I)
  181 CONTINUE
      GO TO 194
  182 CONTINUE
      DO 183 I=1,MR
         Q(I,1) = .5*Q(I,1)-Q(I,JM1)+B(I)
  183 CONTINUE
      GO TO 194
  184 CONTINUE
      IF (N .NE. 2) GO TO 188
C
C     CASE  N = 2
C
      DO 185 I=1,MR
         B(I) = Q(I,1)
  185 CONTINUE
      TCOS(1) = 0.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 186 I=1,MR
         Q(I,1) = B(I)
         B(I) = 2.*(Q(I,2)+B(I))*FISTAG
  186 CONTINUE
      TCOS(1) = -FISTAG
      TCOS(2) = 2.
      CALL TRIX (2,0,MR,A,BB,C,B,TCOS,D,W)
      DO 187 I=1,MR
         Q(I,1) = Q(I,1)+B(I)
  187 CONTINUE
      JR = 1
      I2R = 0
      GO TO 194
  188 CONTINUE
C
C     CASE OF GENERAL N AND NR = 2 .
C
      DO 189 I=1,MR
         II = IP+I
         B3(I) = 0.
         B(I) = Q(I,1)+2.*P(II)
         Q(I,1) = .5*Q(I,1)-Q(I,JM1)
         B2(I) = 2.*(Q(I,1)+Q(I,NLAST))
  189 CONTINUE
      K1 = KR+JR-1
      TCOS(K1+1) = -2.
      K4 = K1+3-ISTAG
      CALL COSGEN (KR+ISTAG-2,1,0.0,FNUM,TCOS(K4))
      K4 = K1+KR+1
      CALL COSGEN (JR-1,1,0.0,1.0,TCOS(K4))
      CALL MERGE (TCOS,K1,KR,K1+KR,JR-1,0)
      CALL COSGEN (KR,1,0.5,FDEN,TCOS(K1+1))
      K2 = KR
      K4 = K1+K2+1
      CALL COSGEN (LR,1,0.5,FDEN,TCOS(K4))
      K3 = LR
      K4 = 0
      CALL TRI3 (MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
      DO 190 I=1,MR
         B(I) = B(I)+B2(I)
  190 CONTINUE
      TCOS(1) = 2.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 191 I=1,MR
         Q(I,1) = Q(I,1)+B(I)
  191 CONTINUE
      GO TO 194
  192 DO 193 I=1,MR
         B(I) = Q(I,NLAST)
  193 CONTINUE
      GO TO 196
  194 CONTINUE
C
C     START BACK SUBSTITUTION.
C
      J = NLAST-JR
      DO 195 I=1,MR
         B(I) = Q(I,NLAST)+Q(I,J)
  195 CONTINUE
  196 JM2 = NLAST-I2R
      IF (JR .NE. 1) GO TO 198
      DO 197 I=1,MR
         Q(I,NLAST) = 0.
  197 CONTINUE
      GO TO 202
  198 CONTINUE
      IF (NROD .NE. 0) GO TO 200
      DO 199 I=1,MR
         II = IP+I
         Q(I,NLAST) = P(II)
  199 CONTINUE
      IP = IP-MR
      GO TO 202
  200 DO 201 I=1,MR
         Q(I,NLAST) = Q(I,NLAST)-Q(I,JM2)
  201 CONTINUE
  202 CONTINUE
      CALL COSGEN (KR,1,0.5,FDEN,TCOS)
      CALL COSGEN (LR,1,0.5,FDEN,TCOS(KR+1))
      IF (LR .NE. 0) GO TO 204
      DO 203 I=1,MR
         B(I) = FISTAG*B(I)
  203 CONTINUE
  204 CONTINUE
      CALL TRIX (KR,LR,MR,A,BB,C,B,TCOS,D,W)
      DO 205 I=1,MR
         Q(I,NLAST) = Q(I,NLAST)+B(I)
  205 CONTINUE
      NLASTP = NLAST
  206 CONTINUE
      JSTEP = JR
      JR = I2R
      I2R = I2R/2
      IF (JR .EQ. 0) GO TO 222
      GO TO (207,208),MIXBND
  207 JSTART = 1+JR
      GO TO 209
  208 JSTART = JR
  209 CONTINUE
      KR = KR-JR
      IF (NLAST+JR .GT. N) GO TO 210
      KR = KR-JR
      NLAST = NLAST+JR
      JSTOP = NLAST-JSTEP
      GO TO 211
  210 CONTINUE
      JSTOP = NLAST-JR
  211 CONTINUE
      LR = KR-JR
      CALL COSGEN (JR,1,0.5,0.0,TCOS)
      DO 221 J=JSTART,JSTOP,JSTEP
         JM2 = J-JR
         JP2 = J+JR
         IF (J .NE. JR) GO TO 213
         DO 212 I=1,MR
            B(I) = Q(I,J)+Q(I,JP2)
  212    CONTINUE
         GO TO 215
  213    CONTINUE
         DO 214 I=1,MR
            B(I) = Q(I,J)+Q(I,JM2)+Q(I,JP2)
  214    CONTINUE
  215    CONTINUE
         IF (JR .NE. 1) GO TO 217
         DO 216 I=1,MR
            Q(I,J) = 0.
  216    CONTINUE
         GO TO 219
  217    CONTINUE
         JM1 = J-I2R
         JP1 = J+I2R
         DO 218 I=1,MR
            Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
  218    CONTINUE
  219    CONTINUE
         CALL TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
         DO 220 I=1,MR
            Q(I,J) = Q(I,J)+B(I)
  220    CONTINUE
  221 CONTINUE
      NROD = 1
      IF (NLAST+I2R .LE. N) NROD = 0
      IF (NLASTP .NE. NLAST) GO TO 194
      GO TO 206
  222 CONTINUE
C
C     RETURN STORAGE REQUIREMENTS FOR P VECTORS.
C
      W(1) = IPSTOR
      RETURN
      END
      SUBROUTINE POISP2 (M,N,A,BB,C,Q,IDIMQ,B,B2,B3,W,W2,W3,D,TCOS,P)
C
C     SUBROUTINE TO SOLVE POISSON EQUATION WITH PERIODIC BOUNDARY
C     CONDITIONS.
C
      DIMENSION       A(*)       ,BB(*)      ,C(*)       ,Q(IDIMQ,1) ,
     1                B(*)       ,B2(*)      ,B3(*)      ,W(*)       ,
     2                W2(*)      ,W3(*)      ,D(*)       ,TCOS(*)    ,
     3                P(*)
      MR = M
      NR = (N+1)/2
      NRM1 = NR-1
      IF (2*NR .NE. N) GO TO 107
C
C     EVEN NUMBER OF UNKNOWNS
C
      DO 102 J=1,NRM1
         NRMJ = NR-J
         NRPJ = NR+J
         DO 101 I=1,MR
            S = Q(I,NRMJ)-Q(I,NRPJ)
            T = Q(I,NRMJ)+Q(I,NRPJ)
            Q(I,NRMJ) = S
            Q(I,NRPJ) = T
  101    CONTINUE
  102 CONTINUE
      DO 103 I=1,MR
         Q(I,NR) = 2.*Q(I,NR)
         Q(I,N) = 2.*Q(I,N)
  103 CONTINUE
      CALL POISD2 (MR,NRM1,1,A,BB,C,Q,IDIMQ,B,W,D,TCOS,P)
      IPSTOR = W(1)
      CALL POISN2 (MR,NR+1,1,1,A,BB,C,Q(1,NR),IDIMQ,B,B2,B3,W,W2,W3,D,
     1             TCOS,P)
      IPSTOR = MAX0(IPSTOR,INT(W(1)))
      DO 105 J=1,NRM1
         NRMJ = NR-J
         NRPJ = NR+J
         DO 104 I=1,MR
            S = .5*(Q(I,NRPJ)+Q(I,NRMJ))
            T = .5*(Q(I,NRPJ)-Q(I,NRMJ))
            Q(I,NRMJ) = S
            Q(I,NRPJ) = T
  104    CONTINUE
  105 CONTINUE
      DO 106 I=1,MR
         Q(I,NR) = .5*Q(I,NR)
         Q(I,N) = .5*Q(I,N)
  106 CONTINUE
      GO TO 118
  107 CONTINUE
C
C     ODD  NUMBER OF UNKNOWNS
C
      DO 109 J=1,NRM1
         NRPJ = N+1-J
         DO 108 I=1,MR
            S = Q(I,J)-Q(I,NRPJ)
            T = Q(I,J)+Q(I,NRPJ)
            Q(I,J) = S
            Q(I,NRPJ) = T
  108    CONTINUE
  109 CONTINUE
      DO 110 I=1,MR
         Q(I,NR) = 2.*Q(I,NR)
  110 CONTINUE
      LH = NRM1/2
      DO 112 J=1,LH
         NRMJ = NR-J
         DO 111 I=1,MR
            S = Q(I,J)
            Q(I,J) = Q(I,NRMJ)
            Q(I,NRMJ) = S
  111    CONTINUE
  112 CONTINUE
      CALL POISD2 (MR,NRM1,2,A,BB,C,Q,IDIMQ,B,W,D,TCOS,P)
      IPSTOR = W(1)
      CALL POISN2 (MR,NR,2,1,A,BB,C,Q(1,NR),IDIMQ,B,B2,B3,W,W2,W3,D,
     1             TCOS,P)
      IPSTOR = MAX0(IPSTOR,INT(W(1)))
      DO 114 J=1,NRM1
         NRPJ = NR+J
         DO 113 I=1,MR
            S = .5*(Q(I,NRPJ)+Q(I,J))
            T = .5*(Q(I,NRPJ)-Q(I,J))
            Q(I,NRPJ) = T
            Q(I,J) = S
  113    CONTINUE
  114 CONTINUE
      DO 115 I=1,MR
         Q(I,NR) = .5*Q(I,NR)
  115 CONTINUE
      DO 117 J=1,LH
         NRMJ = NR-J
         DO 116 I=1,MR
            S = Q(I,J)
            Q(I,J) = Q(I,NRMJ)
            Q(I,NRMJ) = S
  116    CONTINUE
  117 CONTINUE
  118 CONTINUE
C
C     RETURN STORAGE REQUIREMENTS FOR P VECTORS.
C
      W(1) = IPSTOR
      RETURN
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
C-----------------------------------------------------------------------
      END
