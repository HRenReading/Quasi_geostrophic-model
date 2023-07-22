      SUBROUTINE HFFT2A (COEFU, NX, NY, H, GH, LDXGH, BCTY, BD1, BD2,
     *                   BD3, BD4, IORDER, U, LDXU, WORK, NWORK, INFO)
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C    ---------------      4TH ORDER ACCURATE FAST DIRECT SOLUTION
C    PACKAGE :  HFFT      OF THE HELMHOLTZ EQUATION ON RECTANGULAR
C    ---------------      DOMAINS IN TWO AND THREE DIMENSIONS
C
C       HFFT2   =  2 DIMENSIONS, PROBLEM DEFINED BY FUNCTIONS
C       HFFT2A  =  2 DIMENSIONS, PROBLEM DEFINED BY ARRAYS
C       HFFT3   =  3 DIMENSIONS, PROBLEM DEFINED BY FUNCTIONS
C       HFFT3A  =  3 DIMENSIONS, PROBLEM DEFINED BY ARRAYS
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C
C   P U R P O S E
C   -------------
C
C   HFFT2A SOLVES THE HELMHOLTZ EQUATION IN CARTESIAN COORDINATES
C   ON A TWO-DIMENSIONAL RECTANGULAR DOMAIN WITH ANY COMBINATION
C   OF DIRICHLET, NEUMANN, OR AND PERIODIC BOUNDARY CONDITIONS.
C
C
C   D E S C R I P T I O N
C   ---------------------
C
C   HFFT2A SOLVES THE EQUATION
C
C                   2       2
C                  D U  +  D U
C                  ---  +  ---  +  COEFU*U   =   G
C                    2       2
C                  DX      DY
C
C   ON THE RECTANGULAR DOMAIN (AX,BX)X(AY,BY)
C
C                         TOP  (SIDE 4)
C              Y=BY ----------------------------
C                   :                          :
C                   :                          :
C              LEFT :                          : RIGHT
C                   :                          :
C          (SIDE 3) :                          : (SIDE 1)
C                   :                          :
C                   :                          :
C              Y=AY ----------------------------
C                  X=AX   BOTTOM  (SIDE 2)   X=BX
C
C   WITH SOME COMBINATION OF DIRICHLET (SOLUTION PRESCRIBED), NEUMANN
C   (FIRST DERIVATIVE PRESCRIBED), OR PERIODIC BOUNDARY CONDITIONS.
C
C   THE OUTPUT OF THIS PROGRAM IS A TWO-DIMENSIONAL ARRAY GIVING
C   ESTIMATES OF THE SOLUTION AT A SET OF GRID POINTS (X(I),Y(J)), FOR
C   I=1,..,NX AND J=1,..,NY, WHERE
C
C                      X(I) = AX + (I-1)*H
C                      Y(J) = AY + (J-1)*H
C
C   AX IS THE X-COORDINATE OF THE LEFT EDGE OF THE RECTANGLE AND AY IS
C   THE Y-COORDINATE OF THE BOTTOM EDGE OF THE RECTANGLE.
C
C   WHEN COEFU=0 AND ONLY NEUMANN OR PERIODIC BOUNDARY CONDITIONS ARE
C   PRESCRIBED, THEN ANY CONSTANT MAY BE ADDED TO THE SOLUTION TO
C   OBTAIN ANOTHER SOLUTION TO THE PROBLEM. IN THIS CASE THE SOLUTION
C   OF MINIMUM INFINITY NORM IS RETURNED.
C
C   THE SOLUTION IS COMPUTED USING EITHER A SECOND OR FOURTH ORDER
C   ACCURATE FINITE DIFFERENCE APPROXIMATION OF THE CONTINUOUS EQUATION.
C   THE RESULTING SYSTEM OF LINEAR ALGEBRAIC EQUATIONS IS SOLVED USING
C   FAST FOURIER TRANSFORM TECHNIQUES. THE ALGORITHM RELIES UPON THE
C   FACT THAT NX-1 IS HIGHLY COMPOSITE (THE PRODUCT OF SMALL PRIMES).
C
C
C   P A R A M E T E R S
C   -------------------
C
C   COEFU    REAL SCALAR (INPUT)
C            THE COEFFICIENT OF U IN THE PARTIAL DIFFERENTIAL EQUATION.
C
C   NX       INTEGER SCALAR (INPUT) .GE. 4
C            THE NUMBER OF GRID LINES IN X.
C            (THE NUMBER OF VERTICAL PANELS IS THEN NX-1.)
C            NX SHOULD BE CHOSEN SO THAT NX-1 IS HIGHLY COMPOSITE
C            (THE PRODUCT OF SMALL PRIMES) TO INCREASE THE SPEED
C            OF THE FOURIER TRANSFORM ALGORITHM.
C
C   NY       INTEGER SCALAR (INPUT) .GE. 4
C            THE NUMBER OF GRID LINES IN Y.
C            (THE NUMBER OF HORIZONTAL PANELS IS THEN NY-1.)
C
C   H        REAL SCALAR (INPUT)
C            THE SPACE BETWEEN ADJACENT GRID LINES (IN BOTH X AND Y).
C
C   GH       REAL ARRAY OF SIZE LDXGH BY NY+1 (INPUT)
C            THE RIGHT HAND SIDE OF THE PDE AT HALF GRID POINTS.
C            G(I+1,J+1) IS THE VALUE OF THE RIGHT HAND SIDE AT THE
C            (I,J)TH HALF GRID POINT, I=1,..,NX-1, J=1,..,NY-1. THAT IS,
C
C                    GH(I+1,J+1) = G(XH(I),YH(J))
C                          XH(I) = AX + (I-0.5)*H
C                          YH(J) = AY + (J-0.5)*H
C
C            THE FIRST AND LAST ROWS AND COLUMNS OF GH ARE USED AS
C            WORKING STORAGE.  GH IS NOT USED WHEN IORDER=2.
C
C   LDXGH    INTEGER SCALAR (INPUT) .GE. NX+1
C            THE LEADING DIMENSION OF THE ARRAY GH EXACTLY AS SPECIFIED
C            IN THE CALLING PROGRAM.
C
C   BCTY     INTEGER ARRAY OF SIZE 4 (INPUT)
C            INDICATES TYPE OF BOUNDARY CONDITION ON EACH SIDE OF THE
C            DOMAIN AS FOLLOWS.
C
C               BCTY(1) = TYPE ON RIGHT SIDE  (X=BX)
C               BCTY(2) = TYPE ON BOTTOM SIDE (Y=AY)
C               BCTY(3) = TYPE ON LEFT SIDE   (X=AX)
C               BCTY(4) = TYPE ON TOP SIDE    (Y=BY)
C
C            POSSIBLE VALUES ARE
C
C               1 == U PRESCRIBED
C               2 == DU/DX PRESCRIBED (FOR BCTY(1) AND BCTY(3)) OR
C                    DU/DY PRESCRIBED (FOR BCTY(2) AND BCTY(4))
C               3 == PERIODIC
C
C   BD1      REAL ARRAY OF SIZE NY (INPUT)
C            THE VALUES OF THE BOUNDARY CONDITION AT GRID POINTS ON
C            THE RIGHT SIDE OF THE DOMAIN.  THE VALUE STORED DEPENDS
C            UPON THE TYPE OF BOUNDARY CONDITION AS FOLLOWS
C
C                      BCTY(1)       VALUE STORED
C                    ------------------------------
C                         1               U
C                         2             DU/DX
C                         3             NONE
C
C            BD1(J) GIVES THE VALUE AT THE POINT (BX,Y(J)).
C
C   BD2      REAL ARRAY OF SIZE NX (INPUT)
C            THE VALUES OF THE BOUNDARY CONDITION AT GRID POINTS ON
C            THE BOTTOM SIDE OF THE DOMAIN.  THE VALUE STORED DEPENDS
C            UPON THE TYPE OF BOUNDARY CONDITION AS FOLLOWS
C
C                      BCTY(2)       VALUE STORED
C                    ------------------------------
C                         1               U
C                         2             DU/DY
C                         3             NONE
C
C            BD2(I) GIVES THE VALUE AT THE POINT (X(I),AY).
C
C   BD3      REAL ARRAY OF SIZE NY (INPUT)
C            THE VALUES OF THE BOUNDARY CONDITION AT GRID POINTS ON
C            THE LEFT SIDE OF THE DOMAIN.  THE VALUE STORED DEPENDS
C            UPON THE TYPE OF BOUNDARY CONDITION AS FOLLOWS
C
C                      BCTY(3)       VALUE STORED
C                    ------------------------------
C                         1               U
C                         2             DU/DX
C                         3             NONE
C
C            BD3(J) GIVES THE VALUE AT THE POINT (AX,Y(J)).
C
C   BD4      REAL ARRAY OF SIZE NX (INPUT)
C            THE VALUES OF THE BOUNDARY CONDITION AT GRID POINTS ON
C            THE BOTTOM SIDE OF THE DOMAIN.  THE VALUE STORED DEPENDS
C            UPON THE TYPE OF BOUNDARY CONDITION AS FOLLOWS
C
C                      BCTY(4)       VALUE STORED
C                    ------------------------------
C                         1               U
C                         2             DU/DY
C                         3             NONE
C
C            BD4(I) GIVES THE VALUE AT THE POINT (X(I),BY).
C
C   IORDER   INTEGER SCALAR (INPUT)
C            THE ORDER OF ACCURACY OF THE FINITE DIFFERENCE
C            APPROXIMATION OF THE PARTIAL DIFFERENTIAL EQUATION.
C            POSSIBLE VALUES ARE
C
C              2 == 2ND ORDER ACCURATE 9-POINT COMPACT DIFFERENCES
C              4 == 4TH ORDER ACCURATE 9-POINT COMPACT DIFFERENCES
C
C   U        REAL ARRAY OF SIZE LDXU BY NY+2 (INPUT/OUTPUT)
C            ON INPUT, U(I,J) IS THE VALUE OF THE RIGHT HAND SIDE OF
C            THE PARTIAL DIFFERENTIAL EQUATION AT THE (I,J)TH GRID
C            POINT.  THAT IS,
C
C                U(I,J) = G(X(I),Y(J)),  I=1,..,NX, J=1,..,NY
C
C            ON OUTPUT, U(I,J) IS THE VALUE OF THE COMPUTED SOLUTION
C            AT (X(I),Y(J)).  ROWS NX+1 AND NX+2 AND COLUMNS NY+1 AND
C            NY+2 ARE USED FOR WORKING STORAGE.
C
C   LDXU     INTEGER SCALAR (INPUT) .GE. NX+2
C            THE LEADING DIMENSION OF THE ARRAY U EXACTLY AS SPECIFIED
C            IN THE CALLING PROGRAM.
C
C   WORK     REAL ARRAY OF SIZE NWORK.
C            WORKING STORAGE REQUIRED BY HFFT2A.
C
C   NWORK    INTEGER SCALAR (INPUT) .GE. 5*NX + 5*NY + NX/2 + 15
C            THE LENGTH OF THE ARRAY WORK AS DECLARED IN THE CALLING
C            PROGRAM. THIS MAY BE REDUCED BY 4*NY IS COEFU.LE.0.
C
C   INFO     INTEGER SCALAR (OUTPUT)
C            INDICATES STATUS OF COMPUTED SOLUTION.
C            POSSIBLE VALUES ARE
C
C               2 == WARNING. NO SOLUTION EXISTS UNLESS A CONSISTENCY
C                    CONDITION IS SATISFIED (SEE REF. 4). THE DISCRETE
C                    PROBLEM IS ADJUSTED (BY ADDING A CONSTANT TO THE
C                    RIGHT SIDE) SO THAT THIS CONDITION IS SATISFIED.
C                    THIS CONSTANT IS RETURNED IN WORK(1). IF IT IS NOT
C                    SMALL THEN THE PROBLEM MAY NOT BE WELL-POSED. IN
C                    ADDITION, THE SOLUTION IS UNIQUE ONLY UP TO AN
C                    ADDITIVE CONSTANT.
C               1 == WARNING. COEFU.GT.0  A SOLUTION MAY NOT EXIST IF
C                    COEFU IS AN EIGENVALUE OF THE LAPLACIAN.  IF COEFU
C                    IS NEAR ONE OF THESE VALUES THEN THE COMPUTED
C                    SOLUTION MAY BE UNRELIABLE.
C               0 == SUCCESS. SUBPROGRAM RAN TO COMPLETION.
C              -1 == ERROR. NX.LT.4
C              -2 == ERROR. NY.LT.4
C              -3 == ERROR. LDXU.LT.NX+2
C              -4 == ERROR. IORDER NOT 2 OR 4.
C              -5 == ERROR. ELEMENT OF BCTY NOT 1, 2 OR 3.
C              -6 == ERROR. PERIODIC BOUNDARY CONDITIONS SPECIFIED ON
C                           ONE SIDE OF DOMAIN BUT NOT THE OPPOSITE SIDE
C              -7 == ERROR. NWORK TOO SMALL.
C              -8 == ERROR. LDXGH.LT.NX+1
C
C
C   E X T E R N A L   R E F E R E N C E S
C   -------------------------------------
C
C     FDIS2,HDIS2,FD2N,
C     STORD2,FD2D,HD2N,REFL2,
C     MDALG2,EVDISC,TRISOL,
C     TRSALL,TRSOLG,TRSOLP,
C     SGPSL,FFTI,FFTB,FFTF     -- THIS PACKAGE
C
C     RFFTI,RFFTF,RFFTB,
C     SINTI,SINT,COSTI,COST,
C     SINQI,SINQF,SINQB,
C     COSQI,COSQF,COSQB        -- FFTPACK (SEE REF. 2)
C
C     R1MACH                   -- MACHINE CONSTANTS (SEE REF. 3)
C
C
C   P O R T A B I L I T Y
C   ---------------------
C
C   THIS PACKAGE IS WRITTEN IN ANSI STANDARD FORTRAN (1977).
C   ALL MACHINE-DEPENDENT QUANTITIES ARE OBTAINED FROM THE
C   FUNCTION R1MACH (SEE REF. 3).
C
C
C   R E F E R E N C E S
C   -------------------
C
C     1) R. BOISVERT, A FOURTH ORDER ACCURATE FAST DIRECT METHOD
C        FOR THE HELMHOLTZ EQUATION, IN ELLIPTIC PROBLEM SOLVERS
C        II (G. BIRKHOFF AND A. SCHOENSTADT, EDS.), ACADEMIC PRESS,
C        ORLANDO, FLA., 1984, 35-44.
C
C     2) THE FFT PACKAGE USED IS A SLIGHTLY MODIFIED VERSION OF THE
C        PACKAGE FFTPACK WRITTEN BY PAUL SWARZTRAUBER. FFTPACK IS
C        ALSO USED BY THE SUBPROGRAM HW3CRT OF FISHPAK (VERSION 3).
C        FFTPACK IS DESCRIBED IN
C
C          P.N. SWARZTRAUBER, VECTORIZING THE FFTS, IN PARALLEL
C          COMPUTATIONS (G. RODRIGUE, ED.), ACADEMIC PRESS, 1982,
C          PP. 51-83.
C
C        FOR FURTHER INFORMATION WRITE INFORMATION SERVICES OFFICE,
C        COMPUTING FACILITY, NATIONAL CENTER FOR ATMOSPHERIC RESEARCH,
C        BOX 3000, BOULDER, CO 80303, USA.
C
C     3) P. FOX, A. HALL, AND N. SCHRYER, ALGORITHM 528: FRAMEWORK
C        FOR A PORTABLE LIBRARY, ACM TRANS. MATH. SOFT. 4 (1978),
C        PP. 177-188.
C
C     4) S. G. MIKHLIN (ED.), LINEAR EQUATIONS OF MATHEMATICAL PHYSICS,
C        HOLT, RINEHART AND WINSTON, NEW YORK, 1967.
C
C
C   A U T H O R  /  V E R S I O N
C   -----------------------------
C
C     RONALD F. BOISVERT
C     CENTER FOR APPLIED MATHEMATICS
C     NATIONAL BUREAU OF STANDARDS
C     GAITHERSBURG, MD 20899
C     USA
C
C     ORIGINAL  DECEMBER 1985
C     REVISED   APRIL 1987
C
C
C  ------------
C  DECLARATIONS
C  ------------
C
C     ... PARAMETERS
C
      INTEGER NX, NY, IORDER, LDXGH, LDXU, BCTY(4), NWORK, INFO
      REAL
     *     H, COEFU, GH(LDXGH,*), U(LDXU,*), WORK(NWORK),
     *     BD1(NY), BD2(NX), BD3(NY), BD4(NX)
C
C     ... LOCAL VARIABLES
C
      LOGICAL SNGULR, PRDX, PRDY
      INTEGER NXM1, NYM1, NM, NEEDED
      REAL
     *     EPSM, R1MACH, A, B, C, FX, FY, FACTOR, SCALE, PERTRB, UNORM
C
C     ... LOCAL CONSTANTS
C
      INTEGER DRCH, NEUM, PRDC, LEFT, RIGHT, TOP, BOTTOM
      PARAMETER (DRCH=1, NEUM =2, PRDC=3,
     *           LEFT=3, RIGHT=1, TOP =4, BOTTOM=2)
C
C
C  ---------------
C  INITIALIZATIONS
C  ---------------
C
      EPSM = R1MACH(4)
C
      NXM1 = NX - 1
      NYM1 = NY - 1
C
      SNGULR = (BCTY(RIGHT ) .NE. DRCH) .AND.
     *         (BCTY(BOTTOM) .NE. DRCH) .AND.
     *         (BCTY(LEFT  ) .NE. DRCH) .AND.
     *         (BCTY(TOP   ) .NE. DRCH) .AND.
     *         (ABS(COEFU) .LT. R1MACH(4))
      PRDX = BCTY(RIGHT ) .EQ. PRDC
      PRDY = BCTY(BOTTOM) .EQ. PRDC
C
      IL = 1
      IR = NX
      JL = 1
      JR = NY
      IF (BCTY(RIGHT ) .EQ. DRCH)  IR = NXM1
      IF (BCTY(BOTTOM) .EQ. DRCH)  JL = 2
      IF (BCTY(LEFT  ) .EQ. DRCH)  IL = 2
      IF (BCTY(TOP   ) .EQ. DRCH)  JR = NYM1
      IF (PRDX)  IR = NXM1
      IF (PRDY)  JR = NYM1
      NM = (IR-IL+1)*(JR-JL+1)
C
C
C  ---------------------------
C  CHECK VALIDITY OF ARGUMENTS
C  ---------------------------
C
      INFO = 0
      IF (NX .LT. 4)  GO TO 901
      IF (NY .LT. 4)  GO TO 902
      IF (LDXU .LT. NX+2)  GO TO 903
      IF ((IORDER .NE. 2) .AND. (IORDER .NE. 4))  GO TO 904
      DO 10 K=1,4
         IF ((BCTY(K) .LT. 1) .OR. (BCTY(K) .GT. 3))  GO TO 905
   10 CONTINUE
      IF (((BCTY(RIGHT ) .EQ. PRDC) .AND. (BCTY(LEFT  ) .NE. PRDC)) .OR.
     *    ((BCTY(LEFT  ) .EQ. PRDC) .AND. (BCTY(RIGHT ) .NE. PRDC)) .OR.
     *    ((BCTY(BOTTOM) .EQ. PRDC) .AND. (BCTY(TOP   ) .NE. PRDC)) .OR.
     *    ((BCTY(TOP   ) .EQ. PRDC) .AND. (BCTY(BOTTOM) .NE. PRDC))    )
     *   GO TO 906
      NEEDED = 5*NX + 5*NY + NX/2 + 15
      IF (NWORK .LT. NEEDED)  GO TO 907
      IF (LDXGH .LT. NX+1)  GO TO 908
C
      IF (COEFU .GT. 0.0E0)  INFO = 1
      IF (SNGULR)  INFO = 2
C
C
C  ----------
C  DISCRETIZE
C  ----------
C
      IF (IORDER .EQ. 4)  THEN
         CALL HDIS2(NX,NY,H,COEFU,GH,LDXGH-1,BCTY,BD1,BD2,BD3,BD4,
     *              A,B,C,U,LDXU-1,WORK)
      ELSE
         CALL FDIS2(NX,NY,H,COEFU,BCTY,BD1,BD2,BD3,BD4,
     *              A,B,C,U,LDXU)
      ENDIF
C
C  ---------------------------------------
C  ADJUST FOR CONSISTENCY IN SINGULAR CASE
C  ---------------------------------------
C
      PERTRB = 0.0E0
      IF (SNGULR) THEN
         SCALE = 0.0E0
         DO 100 J=JL,JR
            FY = 1.0E0
            IF (.NOT.PRDY.AND.(J.NE.1).AND.(J.NE.NY)) FY = 2.0E0
            DO 100 I=IL,IR
               FX = 1.0E0
               IF (.NOT.PRDX.AND.(I.NE.1).AND.(I.NE.NX)) FX = 2.0E0
               FACTOR = FX*FY
               PERTRB = PERTRB + FACTOR*U(I,J)
               SCALE = SCALE + FACTOR
  100    CONTINUE
         PERTRB = -PERTRB/SCALE
         DO 110 J=JL,JR
            DO 110 I=IL,IR
               U(I,J) = U(I,J) + PERTRB
  110    CONTINUE
      ENDIF
C
C
C  ------------------------------
C  MATRIX DECOMPOSITION USING FFT
C  ------------------------------
C
      CALL MDALG2(A,B,C,BCTY,U,LDXU,IL,IR,JL,JR,WORK)
C
C
C  -----------------------------------------
C  SELECT MIN NORM SOLUTION IN SINGULAR CASE
C  -----------------------------------------
C
      IF (SNGULR) THEN
         UNORM = 0.0E0
         DO 210 J=JL,JR
            DO 210 I=IL,IR
               UNORM = UNORM + U(I,J)
  210    CONTINUE
         UNORM = UNORM/REAL(NM)
         DO 220 J=JL,JR
            DO 220 I=IL,IR
               U(I,J) = U(I,J) - UNORM
  220    CONTINUE
      ENDIF
C
C
C  --------------------------------------
C  COPY IDENTICAL LINES IN PERIODIC CASES
C  --------------------------------------
C
      IF (PRDX) THEN
         DO 325 J=1,NY
            U(NX,J) = U(1,J)
  325    CONTINUE
      ENDIF
C
      IF (PRDY) THEN
         DO 375 I=1,NX
            U(I,NY) = U(I,1)
  375    CONTINUE
      ENDIF
C
C
C  -----------
C  NORMAL EXIT
C  -----------
C
      WORK(1) = PERTRB
      GO TO 999
C
C
C  -----------
C  ERROR EXITS
C  -----------
C
  901 CONTINUE
      INFO = -1
      GO TO 999
C
  902 CONTINUE
      INFO = -2
      GO TO 999
C
  903 CONTINUE
      INFO = -3
      GO TO 999
C
  904 CONTINUE
      INFO = -4
      GO TO 999
C
  905 CONTINUE
      INFO = -5
      GO TO 999
C
  906 CONTINUE
      INFO = -6
      GO TO 999
C
  907 CONTINUE
      INFO = -7
      GO TO 999
C
  908 CONTINUE
      INFO = -8
      GO TO 999
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE EVDISC (KBCL, KBCR, EIGEN, N)
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C    ---------------      4TH ORDER ACCURATE FAST DIRECT SOLUTION
C    PACKAGE :  HFFT      OF THE HELMHOLTZ EQUATION ON RECTANGULAR
C    ---------------      DOMAINS IN TWO AND THREE DIMENSIONS
C
C    INTERNAL MODULE
C
C    RONALD F. BOISVERT
C    NATIONAL BUREAU OF STANDARDS
C    DECEMBER 1985
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C
C   EVDISC COMPUTES THE EIGENVALUES OF THE FOLLOWING MATRIX
C
C                  :--                   --:
C                  :   0  R            T   :
C                  :   1  0  1             :
C                  :     1  0  1           :
C                  :       .  .  .         :
C                  :         .  .  .       :
C                  :           1  0  1     :
C                  :             1  0  1   :
C                  :   T            S  0   :
C                  :--                   --:
C
C   WHERE THE SCALARS R, S, AND T DEPEND UPON THE PARAMETERS KBCL
C   AND KBCR IN THE FOLLOWING WAY.
C
C           R = 1  (UNLESS KBCL=2 IN WHICH CASE R=2)
C           S = 1  (UNLESS KBCR=2 IN WHICH CASE S=2)
C           T = 0  (UNLESS KBCR=3 IN WHICH CASE T=1)
C
C
C   P A R A M E T E R S
C   -------------------
C
C     KBCL, KBCR    INTEGER SCALARS (INPUT)
C                   INDICATE THE TYPE OF BOUNDARY CONDITIONS AT
C                   THE LEFT AND RIGHT ENDS OF AN INTERVAL.
C                   USES SAME CONVENTIONS AS BCTY IN HFFT2A.
C
C     EIGEN         REAL ARRAY OF SIZE N (OUTPUT)
C                   THE EIGENVALUES OF THE MATRIX SHOWN ABOVE.
C
C     N             INTEGER SCALAR (INPUT)
C                   THE ORDER OF THE MATRIX.
C
C
C  ------------
C  DECLARATIONS
C  ------------
C
C     ... PARAMETERS
C
      INTEGER KBCL, KBCR, N
      REAL
     *     EIGEN(N)
C
C     ... LOCAL VARIABLES
C
      REAL
     *     PI, FACTOR
C
C     ... LOCAL CONSTANTS
C
      INTEGER DRCH, NEUM, PRDC, LEFT, RIGHT, TOP, BOTTOM
      PARAMETER (DRCH=1, NEUM =2, PRDC=3,
     *           LEFT=3, RIGHT=1, TOP =4, BOTTOM=2)
C
C
      PI = 4.0E0*ATAN(1.0E0)
C
C  -----------
C  SELECT CASE
C  -----------
C
      GO TO (10,20,30), KBCL
   10 GO TO (100,200), KBCR
   20 GO TO (200,300), KBCR
   30 GO TO 400
C
C     CASE : DIRICHLET/DIRICHLET
C
  100 CONTINUE
      FACTOR = PI/REAL(N+1)
      DO 110 I=1,N
         EIGEN(I) = 2.0E0*COS(REAL(I)*FACTOR)
  110 CONTINUE
      GO TO 500
C
C     CASE : DIRICHLET/NEUMANN
C
  200 CONTINUE
      FACTOR = PI/REAL(2*N)
      DO 210 I=1,N
         EIGEN(I) = 2.0E0*COS(REAL(2*I-1)*FACTOR)
  210 CONTINUE
      GO TO 500
C
C     CASE : NEUMANN/NEUMANN
C
  300 CONTINUE
      FACTOR = PI/REAL(N-1)
      DO 310 I=1,N
         EIGEN(I) = 2.0E0*COS(REAL(I-1)*FACTOR)
  310 CONTINUE
      GO TO 500
C
C     CASE : PERIODIC
C
  400 CONTINUE
      FACTOR = 2.0E0*PI/REAL(N)
      DO 410 I=1,N
         K = I/2
         EIGEN(I) = 2.0E0*COS(REAL(K)*FACTOR)
  410 CONTINUE
      GO TO 500
C
C  ----
C  EXIT
C  ----
C
  500 CONTINUE
      RETURN
      END
      SUBROUTINE FDIS2 (NX, NY, H, COEFU, BCTY, BD1, BD2, BD3, BD4,
     *                  A, B, C, U, LDXU)
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C    ---------------      4TH ORDER ACCURATE FAST DIRECT SOLUTION
C    PACKAGE :  HFFT      OF THE HELMHOLTZ EQUATION ON RECTANGULAR
C    ---------------      DOMAINS IN TWO AND THREE DIMENSIONS
C
C    INTERNAL MODULE
C
C    RONALD F. BOISVERT
C    NATIONAL BUREAU OF STANDARDS
C    DECEMBER 1985
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C
C   FDIS2 COMPUTES A SECOND ORDER FINITE DIFFERENCE DISCRETIZATION
C   FOR A TWO-DIMENSIONAL RECTANGULAR DOMAIN
C
C
C   P A R A M E T E R S
C   -------------------
C
C     NX, NY    INTEGER SCALARS (INPUT)
C               SEE HFFT2A.
C
C     H, COEFU  REAL SCALARS (INPUT)
C               SEE HFFT2A.
C
C     BCTY      INTEGER ARRAY OF SIZE 4 (INPUT)
C               SEE HFFT2A.
C
C     BD1, BD3  REAL ARRAYS OF SIZE NY (INPUT)
C               SEE HFFT2A.
C
C     BD2, BD4  REAL ARRAYS OF SIZE NX (INPUT)
C               SEE HFFT2A.
C
C     A, B, C   REAL SCALARS (OUTPUT)
C               GIVES VALUES IN THE BASIC FINITE DIFFERENCE STENCIL
C               (SCALED TO O(1))
C
C                        C  B  C
C                        B  A  B  U  =  RIGHT SIDE
C                        C  B  C
C
C     U         REAL ARRAY OF SIZE LDXU BY NY (INPUT/OUTPUT)
C               ON INPUT, U(I,J) IS THE RIGHT HAND SIDE OF THE PDE
C               EVALUATED AT THE (I,J)TH GRID POINT.
C               ON OUTPUT, U(I,J) IS THE RIGHT HAND SIDE OF THE
C               DISCRETE PDE AT THE (I,J)TH GRID POINT.
C
C     LDXU      INTEGER SCALAR (INPUT)
C               SEE HFFT2A.
C
C
C  ------------
C  DECLARATIONS
C  ------------
C
C     ... PARAMETERS
C
      INTEGER BCTY(4), NX, NY, LDXU
      REAL
     *     BD1(NY), BD2(NX), BD3(NY), BD4(NX), U(LDXU,NY),
     *     COEFU, H, A, B, C
C
C     ... LOCAL VARIABLES
C
      LOGICAL HAVED, HAVEN
      INTEGER I, J, NXM1, NYM1
      REAL
     *     BETA0, GAMMA0, DELTA0, DELTA1, H2, F
C
      COMMON /FD2COM/ GAMMA0, DELTA0, DELTA1
C
C     ... LOCAL CONSTANTS
C
      INTEGER DRCH, NEUM, PRDC, LEFT, RIGHT, TOP, BOTTOM
      PARAMETER (DRCH=1, NEUM =2, PRDC=3,
     *           LEFT=3, RIGHT=1, TOP =4, BOTTOM=2)
C
C
C  ---------------
C  INITIALIZATIONS
C  ---------------
C
      NXM1 = NX - 1
      NYM1 = NY - 1
C
      HAVED = (BCTY(RIGHT) .EQ. DRCH) .OR. (BCTY(BOTTOM) .EQ. DRCH) .OR.
     *        (BCTY(LEFT ) .EQ. DRCH) .OR. (BCTY(TOP   ) .EQ. DRCH)
      HAVEN = (BCTY(RIGHT) .EQ. NEUM) .OR. (BCTY(BOTTOM) .EQ. NEUM) .OR.
     *        (BCTY(LEFT ) .EQ. NEUM) .OR. (BCTY(TOP   ) .EQ. NEUM)
C
      H2 = H*H
      F = -H2*COEFU
      A = -(20.0E0 + 6.0E0*F)
      B =    4.0E0
      C =    1.0E0
      BETA0 = 6.0E0*H2
      GAMMA0 = -12.0E0*H
      DELTA0 = -10.0E0*H
      DELTA1 = -2.0E0*H
C
C
C  ---------------------------------
C  DISCRETIZE RIGHT HAND SIDE OF PDE
C  ---------------------------------
C
      ISTRT = 2
      ISTOP = NXM1
      IF (BCTY(LEFT  ) .NE. DRCH)  ISTRT = 1
      IF (BCTY(RIGH T) .EQ. NEUM)  ISTOP = NX
      JSTRT = 2
      JSTOP = NYM1
      IF (BCTY(BOTTOM) .NE. DRCH)  JSTRT = 1
      IF (BCTY(TOP   ) .EQ. NEUM)  JSTOP = NY
C
      DO 100 J=JSTRT,JSTOP
         DO 100 I=ISTRT,ISTOP
            U(I,J) =  BETA0*U(I,J)
  100 CONTINUE
C
C  ---------------------------------------
C  DISCRETIZE POINTS ON NEUMANN BOUNDARIES
C  ---------------------------------------
C
      IF (HAVEN)  CALL FD2N(NX,NY,BCTY,BD1,BD2,BD3,BD4,U,LDXU)
C
C  ----------------------------------------------
C  ADJUST POINTS ADJACENT TO DIRICHLET BOUNDARIES
C  ----------------------------------------------
C
      IF (HAVED) THEN
         CALL STORD2(NX,NY,BCTY,BD1,BD2,BD3,BD4,U,LDXU)
         CALL FD2D(NX,NY,U,LDXU,BCTY,B,C,U)
      ENDIF
C
C
C  ----
C  EXIT
C  ----
C
      RETURN
      END
      SUBROUTINE FD2D (NX, NY, UD, LDXU, BCTY, B, C, U)
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C    ---------------      4TH ORDER ACCURATE FAST DIRECT SOLUTION
C    PACKAGE :  HFFT      OF THE HELMHOLTZ EQUATION ON RECTANGULAR
C    ---------------      DOMAINS IN TWO AND THREE DIMENSIONS
C
C    INTERNAL MODULE
C
C    RONALD F. BOISVERT
C    NATIONAL BUREAU OF STANDARDS
C    DECEMBER 1985
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C
C   FD2D ELIMINATES KNOWN TERMS FROM EQUATIONS CORRESPONDING TO POINTS
C   NEAR DIRICHLET BOUNDARIES OF A TWO-DIMENSIONAL RECTANGULAR DOMAIN.
C
C
C   P A R A M E T E R S
C   -------------------
C
C     NX, NY    INTEGER SCALARS (INPUT)
C               SEE HFFT2A.
C
C     UD        REAL ARRAY OF SIZE LDXU BY NY (INPUT)
C               ENTRIES CORRESPONDING TO DIRICHLET POINTS CONTAIN
C               KNOWN VALUES OF THE SOLUTION.
C
C     LDXU      INTEGER SCALAR (INPUT)
C               THE LEADING DIMENSION OF THE ARRAYS UD AND U EXACTLY
C               AS DECLARED IN THE CALLING PROGRAM.
C
C     BCTY      INTEGER ARRAY OF SIZE 4 (INPUT)
C               SEE HFFT2A.
C
C     B, C      REAL SCALARS (INPUT)
C               FINITE DIFFERENCE STENCIL COEFFICIENTS TO USE IN THE
C               ELIMINATION
C
C                             C B C
C                             B * B
C                             C B C
C
C     U         REAL ARRAY OF SIZE LDXU BY NY (INPUT/OUTPUT)
C               ON INPUT, ENTRIES CORRESPONDING TO POINTS WHERE THE
C               SOLUTION IS TO BE DETERMINED CONTAIN THE RIGHT HAND
C               SIDE OF A FINITE DIFFERENCE DISCRETIZATION.
C               ON EXIT, THESE ENTRIES ARE UPDATED SUCH THAT KNOWN
C               TERMS (DIRICHLET POINTS) ARE ELIMINATED FROM THE
C               LEFT HAND SIDE OF THE EQUATION.
C
C
C  ------------
C  DECLARATIONS
C  ------------
C
C     ... PARAMETERS
C
      INTEGER BCTY(4), NX, NY
      REAL
     *     UD(LDXU,NY), U(LDXU,NY), B, C
C
C     ... LOCAL VARIABLES
C
      LOGICAL PRDX, PRDY
      INTEGER I, J, NXM1, NXM2, NYM1, NYM2
C
C     ... LOCAL CONSTANTS
C
      INTEGER DRCH, NEUM, PRDC, LEFT, RIGHT, TOP, BOTTOM
      PARAMETER (DRCH=1, NEUM =2, PRDC=3,
     *           LEFT=3, RIGHT=1, TOP =4, BOTTOM=2)
C
C
C  ---------------
C  INITIALIZATIONS
C  ---------------
C
      NXM1 = NX - 1
      NYM1 = NY - 1
      NXM2 = NX - 2
      NYM2 = NY - 2
      PRDX = BCTY(RIGHT ) .EQ. PRDC
      PRDY = BCTY(BOTTOM) .EQ. PRDC
C
C
C  -----------------------------------------------
C  ADJUST POINTS ADJACENT TO DIRICHLET BOUNDARIES
C  -----------------------------------------------
C
      IL = 2
      IR = NXM1
      JL = 2
      JR = NYM1
      IF (BCTY(LEFT  ) .EQ. DRCH)  IL = 3
      IF (BCTY(RIGHT ) .EQ. DRCH)  IR = NXM2
      IF (BCTY(BOTTOM) .EQ. DRCH)  JL = 3
      IF (BCTY(TOP   ) .EQ. DRCH)  JR = NYM2
C
C     ... INTERIOR POINTS NEAR BOTTOM EDGE
C
      IF (BCTY(BOTTOM) .EQ. DRCH) THEN
         DO 505 I=IL,IR
            U(I,2) = U(I,2) - (C*UD(I-1,1)+B*UD(I,1)+C*UD(I+1,1))
  505    CONTINUE
      ENDIF
C
C     ... INTERIOR POINTS NEAR TOP EDGE
C
      IF (BCTY(TOP) .EQ. DRCH) THEN
         DO 515 I=IL,IR
            U(I,NYM1)= U(I,NYM1) -(C*UD(I-1,NY)+B*UD(I,NY)+C*UD(I+1,NY))
  515    CONTINUE
      ENDIF
C
C     ... INTERIOR POINTS NEAR LEFT EDGE
C
      IF (BCTY(LEFT) .EQ. DRCH) THEN
         DO 525 J=JL,JR
            U(2,J) = U(2,J) - (C*UD(1,J-1)+B*UD(1,J)+C*UD(1,J+1))
  525    CONTINUE
      ENDIF
C
C     ... INTERIOR POINTS NEAR RIGHT EDGE
C
      IF (BCTY(RIGHT) .EQ. DRCH) THEN
         DO 535 J=JL,JR
            U(NXM1,J)= U(NXM1,J) -(C*UD(NX,J-1)+B*UD(NX,J)+C*UD(NX,J+1))
  535    CONTINUE
      ENDIF
C
C     ... INTERIOR POINTS NEAR CORNERS
C
      IF ((BCTY(RIGHT) .EQ. DRCH) .AND. (BCTY(BOTTOM) .EQ. DRCH))
     *   U(NXM1,2) = U(NXM1,2) - B*(UD(NX,2) + UD(NXM1,1))
     *                         - C*(UD(NX,3) + UD(NX,1) + UD(NXM2,1))
C
      IF ((BCTY(BOTTOM) .EQ. DRCH) .AND. (BCTY(LEFT) .EQ. DRCH))
     *   U(2,2) = U(2,2) - B*(UD(1,2) + UD(2,1))
     *                   - C*(UD(1,3) + UD(1,1) + UD(3,1))
C
      IF ((BCTY(LEFT) .EQ. DRCH) .AND. (BCTY(TOP) .EQ. DRCH))
     *   U(2,NYM1) = U(2,NYM1) - B*(UD(1,NYM1) + UD(2,NY))
     *                         - C*(UD(1,NYM2) + UD(1,NY) + UD(3,NY))
C
      IF ((BCTY(TOP) .EQ. DRCH) .AND. (BCTY(RIGHT) .EQ. DRCH))
     *   U(NXM1,NYM1) = U(NXM1,NYM1) - B*(UD(NXM1,NY) + UD(NX,NYM1))
     *                         - C*(UD(NXM2,NY)+UD(NX,NY)+UD(NX,NYM2))
C
C     ... NEUMANN POINTS NEAR CORNERS ON BOTTOM EDGE
C
      IF (BCTY(BOTTOM) .EQ. NEUM) THEN
         IF (BCTY(LEFT) .EQ. DRCH)
     *      U(2,1) = U(2,1) - (B*UD(1,1) + 2.0E0*C*UD(1,2))
         IF (BCTY(RIGHT) .EQ. DRCH)
     *      U(NXM1,1) = U(NXM1,1) - (B*UD(NX,1) + 2.0E0*C*UD(NX,2))
      ENDIF
C
C     ... NEUMANN POINTS NEAR CORNERS ON TOP EDGE
C
      IF (BCTY(TOP) .EQ. NEUM) THEN
         IF (BCTY(LEFT) .EQ. DRCH)
     *      U(2,NY) = U(2,NY) - (B*UD(1,NY) + 2.0E0*C*UD(1,NYM1))
         IF (BCTY(RIGHT) .EQ. DRCH)
     *      U(NXM1,NY)= U(NXM1,NY) - (B*UD(NX,NY) + 2.0E0*C*UD(NX,NYM1))
      ENDIF
C
C     ... NEUMANN POINTS NEAR CORNERS ON LEFT EDGE
C
      IF (BCTY(LEFT) .EQ. NEUM) THEN
         IF (BCTY(BOTTOM) .EQ. DRCH)
     *      U(1,2) = U(1,2) - (B*UD(1,1) + 2.0E0*C*UD(2,1))
         IF (BCTY(TOP) .EQ. DRCH)
     *      U(1,NYM1) = U(1,NYM1) - (B*UD(1,NY) + 2.0E0*C*UD(2,NY))
      ENDIF
C
C     ... NEUMANN POINTS NEAR CORNERS ON RIGHT EDGE
C
      IF (BCTY(RIGHT) .EQ. NEUM) THEN
         IF (BCTY(BOTTOM) .EQ. DRCH)
     *      U(NX,2) = U(NX,2) - (B*UD(NX,1) + 2.0E0*C*UD(NXM1,1))
         IF (BCTY(TOP) .EQ. DRCH)
     *      U(NX,NYM1)= U(NX,NYM1) - (B*UD(NX,NY) + 2.0E0*C*UD(NXM1,NY))
      ENDIF
C
C     ... PERIODIC POINTS ON LEFT EDGE
C
      IF (PRDX) THEN
         IF (BCTY(BOTTOM) .EQ. DRCH)
     *      U(1,2) = U(1,2) - (C*UD(NXM1,1) + B*UD(1,1) + C*UD(2,1))
         IF (BCTY(TOP) .EQ. DRCH)
     *      U(1,NYM1)= U(1,NYM1) - (C*UD(NXM1,NY)+B*UD(1,NY)+C*UD(2,NY))
      ENDIF
C
C     ... PERIODIC POINTS ON BOTTOM EDGE
C
      IF (PRDY) THEN
         IF (BCTY(LEFT) .EQ. DRCH)
     *      U(2,1) = U(2,1) - (C*UD(1,NYM1) + B*UD(1,1) + C*UD(1,2))
         IF (BCTY(RIGHT) .EQ. DRCH)
     *      U(NXM1,1)= U(NXM1,1) - (C*UD(NX,NYM1)+B*UD(NX,1)+C*UD(NX,2))
      ENDIF
C
C
C  ----
C  EXIT
C  ----
C
      RETURN
      END
      SUBROUTINE FD2N (NX, NY, BCTY, BD1, BD2, BD3, BD4, U, LDXU)
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C    ---------------      4TH ORDER ACCURATE FAST DIRECT SOLUTION
C    PACKAGE :  HFFT      OF THE HELMHOLTZ EQUATION ON RECTANGULAR
C    ---------------      DOMAINS IN TWO AND THREE DIMENSIONS
C
C    INTERNAL MODULE
C
C    RONALD F. BOISVERT
C    NATIONAL BUREAU OF STANDARDS
C    DECEMBER 1985
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C
C  FD2N COMPUTES THE SECOND ORDER FINITE DIFFERENCE DISCRETIZATION
C  AT ALL BOUNDARY POINTS OF A TWO-DIMENSIONAL RECTANGULAR DOMAIN
C  WHERE NEUMANN BOUNDARY CONDITIONS HAVE BEEN SPECIFIED.
C
C
C   P A R A M E T E R S
C   -------------------
C
C     NX, NY    INTEGER SCALARS (INPUT)
C               SEE HFFT2A.
C
C     BCTY      INTEGER ARRAY OF SIZE 4 (INPUT)
C               SEE HFFT2A.
C
C     BD1, BD3  INTEGER ARRAYS OF SIZE NY (INPUT)
C               SEE HFFT2A.
C
C     BD2, BD4  INTEGER ARRAYS OF SIZE NX (INPUT)
C               SEE HFFT2A.
C
C     U         REAL ARRAY OF SIZE LDXU BY NY (OUTPUT)
C               ON EXIT, ENTRIES OF U CORRESPONDING TO NEUMANN BOUNDARY
C               POINTS CONTAIN THE RIGHT HAND SIDE OF THE SECOND ORDER
C               FINITE DIFFERENCE DISCRETIZATION
C
C     LDXU      INTEGER SCALAR (INPUT)
C               SEE HFFT2A.
C
C
C  ------------
C  DECLARATIONS
C  ------------
C
C     ... PARAMETERS
C
      INTEGER BCTY(4), NX, NY, LDXU
      REAL
     *     BD1(NY), BD2(NX), BD3(NY), BD4(NX), U(LDXU,NY)
C
C     ... LOCAL VARIABLES
C
      LOGICAL PRDX, PRDY
      INTEGER I, J, NXM1, NYM1
      REAL
     *     GAMMA0, DELTA0, DELTA1
C
      COMMON /FD2COM/ GAMMA0, DELTA0, DELTA1
C
C     ... LOCAL CONSTANTS
C
      INTEGER DRCH, NEUM, PRDC, LEFT, RIGHT, TOP, BOTTOM
      PARAMETER (DRCH=1, NEUM =2, PRDC=3,
     *           LEFT=3, RIGHT=1, TOP =4, BOTTOM=2)
C
C
C  ---------------
C  INITIALIZATIONS
C  ---------------
C
      NXM1 = NX - 1
      NYM1 = NY - 1
      PRDX = BCTY(RIGHT ) .EQ. PRDC
      PRDY = BCTY(BOTTOM) .EQ. PRDC
C
C
C  -----------------------------
C  NEUMANN POINT DISCRETIZATIONS
C  -----------------------------
C
C     ... BOTTOM EDGE
C
      IF (BCTY(BOTTOM) .EQ. NEUM) THEN
         DO 305 I=2,NXM1
            U(I,1) = U(I,1) - GAMMA0*BD2(I)
  305    CONTINUE
      ENDIF
C
C     ... TOP EDGE
C
      IF (BCTY(TOP) .EQ. NEUM) THEN
         DO 315 I=2,NXM1
            U(I,NY) = U(I,NY) + GAMMA0*BD4(I)
  315    CONTINUE
      ENDIF
C
C     ... LEFT EDGE
C
      IF (BCTY(LEFT) .EQ. NEUM) THEN
         DO 325 J=2,NYM1
            U(1,J) = U(1,J) - GAMMA0*BD3(J)
  325    CONTINUE
      ENDIF
C
C     ... RIGHT EDGE
C
      IF (BCTY(RIGHT) .EQ. NEUM) THEN
         DO 335 J=2,NYM1
            U(NX,J) = U(NX,J) + GAMMA0*BD1(J)
  335    CONTINUE
      ENDIF
C
C     ... LOWER RIGHT CORNER
C
      IF ((BCTY(RIGHT) .EQ. NEUM) .AND. (BCTY(BOTTOM) .EQ. NEUM))
     *    U(NX,1) = U(NX,1)
     *              + DELTA0*BD1(1)
     *              + DELTA1*BD1(2)
     *              - DELTA0*BD2(NX)
     *              - DELTA1*BD2(NXM1)
C
C     ... LOWER LEFT CORNER
C
      IF ((BCTY(BOTTOM) .EQ. NEUM) .AND. (BCTY(LEFT) .EQ. NEUM))
     *    U(1,1) = U(1,1)
     *             - DELTA0*BD2(1)
     *             - DELTA1*BD2(2)
     *             - DELTA0*BD3(1)
     *             - DELTA1*BD3(2)
C
C     ... UPPER LEFT CORNER
C
      IF ((BCTY(LEFT) .EQ. NEUM) .AND. (BCTY(TOP) .EQ. NEUM))
     *    U(1,NY) = U(1,NY)
     *              - DELTA0*BD3(NY)
     *              - DELTA1*BD3(NYM1)
     *              + DELTA0*BD4(1)
     *              + DELTA1*BD4(2)
C
C     ... UPPER RIGHT CORNER
C
      IF ((BCTY(TOP) .EQ. NEUM) .AND. (BCTY(RIGHT) .EQ. NEUM))
     *    U(NX,NY) = U(NX,NY)
     *               + DELTA0*BD1(NY)
     *               + DELTA1*BD1(NYM1)
     *               + DELTA0*BD4(NX)
     *               + DELTA1*BD4(NXM1)
C
C
C  -----------------------------
C  PERIODIC POINT DISCRETIZATION
C  -----------------------------
C
C     ... LOWER LEFT CORNER (PERIODIC/NEUMANN)
C
      IF (PRDX .AND. (BCTY(BOTTOM) .EQ. NEUM))
     *   U(1,1) = U(1,1) - GAMMA0*BD2(1)
C
C     ... LOWER LEFT CORNER (NEUMANN/PERIODIC)
C
      IF (PRDY .AND. (BCTY(LEFT) .EQ. NEUM))
     *   U(1,1) = U(1,1) - GAMMA0*BD3(1)
C
C     ... UPPER LEFT CORNER (PERIODIC/NEUMANN)
C
      IF (PRDX .AND. (BCTY(TOP) .EQ. NEUM))
     *   U(1,NY) = U(1,NY) + GAMMA0*BD4(1)
C
C     ... LOWER RIGHT CORNER (NEUMANN/PERIODIC)
C
      IF (PRDY .AND. (BCTY(RIGHT) .EQ. NEUM))
     *   U(NX,1) = U(NX,1) + GAMMA0*BD1(1)
C
C
C  ----
C  EXIT
C  ----
C
      RETURN
      END
      SUBROUTINE FFTB (KBCL, KBCR, U, LDXU, IL, IR, JL, JR,  WSAVE)
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C    ---------------      4TH ORDER ACCURATE FAST DIRECT SOLUTION
C    PACKAGE :  HFFT      OF THE HELMHOLTZ EQUATION ON RECTANGULAR
C    ---------------      DOMAINS IN TWO AND THREE DIMENSIONS
C
C    INTERNAL MODULE
C
C    RONALD F. BOISVERT
C    NATIONAL BUREAU OF STANDARDS
C    DECEMBER 1985
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C
C   FFTB PERFORMS M=JR-JL+1 REAL ONE-DIMENSIONAL INVERSE DISCRETE
C   FOURIER TRANSFORMS OF LENGTH N=IR-IL+1 (FOURIER SYNTHESIS).
C   THE TYPE OF TRANSFORM SELECTED DEPENDS UPON THE BOUNDARY CONDITIONS.
C
C
C   P A R A M E T E R S
C   -------------------
C
C     KBCL, KBCR   INTEGER SCALARS (INPUT)
C                  GIVE THE TYPE OF BOUNDARY CONDITIONS AT THE LEFT AND
C                  RIGHT ENDPOINTS OF THE INTERVAL. POSSIBLE VALUES ARE
C                  THE SAME AS IN THE VECTOR BCTY OF HFFT2A.
C
C     U            REAL ARRAY OF SIZE LDXU BY JR (INPUT/OUTPUT)
C                  ON INPUT, U CONTAINS THE SEQUENCES TO BE TRANSFORMED
C                  IN COLUMNS JL TO JR.  WITHIN EACH COLUMN THE
C                  SEQUENCES ARE STORED IN POSITIONS IL TO IR.
C                  ON OUTPUT THESE VALUES ARE REPLACED BY THEIR DISCRETE
C                  FOURIER TRANSFORMS.
C
C     LDXU         INTEGER SCALAR (INPUT)
C                  THE LEADING DIMENSION OF THE ARRAY U EXACTLY AS
C                  SPECIFIED IN THE CALLING PROGRAM.
C
C     IL, IR       INTEGER SCALARS (INPUT)
C     JL, JR       GIVES THE SUBSET OF THE ARRAY U WHICH CONTAINS THE
C                  SEQUENCES TO BE TRANSFORMED, I.E., POSITIONS I=
C                  IL,..,IR OF COLUMNS J=JL,..,JR.
C
C     WSAVE        REAL ARRAY OF SIZE 3*N + N/2 + 15 (INPUT)
C                  THE WORK ARRAY EXACTLY AS RETURNED FROM THE ROUTINE
C                  FFTI.
C
C
C  ------------
C  DECLARATIONS
C  ------------
C
C     ... PARAMETERS
C
      INTEGER KBCL, KBCR, IL, IR, JL, JR
      REAL
     *     U(LDXU,*), WSAVE(*)
C
C     ... LOCAL VARIABLES
C
      INTEGER ICASE
      REAL
     *     SCALE
C
C     ... LOCAL CONSTANTS
C
      INTEGER DRCH, NEUM, PRDC, LEFT, RIGHT, TOP, BOTTOM
      PARAMETER (DRCH=1, NEUM =2, PRDC=3,
     *           LEFT=3, RIGHT=1, TOP =4, BOTTOM=2)
C
      N = IR - IL + 1
C
C  ---------------------------------------------
C  SELECT TRANSFORM BASED ON BOUNDARY CONDITIONS
C  ---------------------------------------------
C
      ICASE = 2*(KBCL-1) + KBCR
      IF (KBCL .EQ. PRDC)  ICASE = 5
      GO TO (100,200,300,400,500), ICASE
C
C     CASE :  DIRICHLET/DIRICHLET
C
  100 CONTINUE
      SCALE = 0.50E0/REAL(N+1)
      DO 150 J=JL,JR
         CALL SINT(N,U(IL,J),WSAVE)
  150 CONTINUE
      GO TO 600
C
C     CASE :  DIRICHLET/NEUMANN
C
  200 CONTINUE
      SCALE = 0.250E0/REAL(N)
      DO 250 J=JL,JR
         CALL SINQB(N,U(IL,J),WSAVE)
  250 CONTINUE
      GO TO 600
C
C     CASE :  NEUMANN/DIRICHLET
C
  300 CONTINUE
      SCALE = 0.250E0/REAL(N)
      DO 350 J=JL,JR
         CALL COSQB(N,U(IL,J),WSAVE)
  350 CONTINUE
      GO TO 600
C
C     CASE :  NEUMANN/NEUMANN
C
  400 CONTINUE
      SCALE = 0.50E0/REAL(N-1)
      DO 450 J=JL,JR
         CALL COST(N,U(IL,J),WSAVE)
  450 CONTINUE
      GO TO 600
C
C     CASE :  PERIODIC
C
  500 CONTINUE
      SCALE = 1.0E0/REAL(N)
      DO 550 J=JL,JR
         CALL RFFTB(N,U(IL,J),WSAVE)
  550 CONTINUE
C
C  -----------------------------------
C  SCALE RESULT TO GET CORRECT INVERSE
C  -----------------------------------
C
  600 CONTINUE
      DO 650 J=JL,JR
         DO 650 I=IL,IR
            U(I,J) = SCALE*U(I,J)
  650 CONTINUE
C
C
C  ----
C  EXIT
C  ----
C
      RETURN
      END
      SUBROUTINE FFTF (KBCL, KBCR, U, LDXU, IL, IR, JL, JR,  WSAVE)
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C    ---------------      4TH ORDER ACCURATE FAST DIRECT SOLUTION
C    PACKAGE :  HFFT      OF THE HELMHOLTZ EQUATION ON RECTANGULAR
C    ---------------      DOMAINS IN TWO AND THREE DIMENSIONS
C
C    INTERNAL MODULE
C
C    RONALD F. BOISVERT
C    NATIONAL BUREAU OF STANDARDS
C    DECEMBER 1985
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C
C   FFTF PERFORMS M=JR-JL+1 REAL ONE-DIMENSIONAL DISCRETE FOURIER
C   TRANSFORMS OF LENGTH N=IR-IL+1 (FOURIER ANALYSIS).
C   THE TYPE OF TRANSFORM SELECTED DEPENDS UPON THE BOUNDARY CONDITIONS.
C
C
C   P A R A M E T E R S
C   -------------------
C
C     KBCL, KBCR   INTEGER SCALARS (INPUT)
C                  GIVE THE TYPE OF BOUNDARY CONDITIONS AT THE LEFT AND
C                  RIGHT ENDPOINTS OF THE INTERVAL. POSSIBLE VALUES ARE
C                  THE SAME AS IN THE VECTOR BCTY OF HFFT2A.
C
C     U            REAL ARRAY OF SIZE LDXU BY JR (INPUT/OUTPUT)
C                  ON INPUT, U CONTAINS THE SEQUENCES TO BE TRANSFORMED
C                  IN COLUMNS JL TO JR.  WITHIN EACH COLUMN THE
C                  SEQUENCES ARE STORED IN POSITIONS IL TO IR.
C                  ON OUTPUT THESE VALUES ARE REPLACED BY THEIR DISCRETE
C                  FOURIER TRANSFORMS.
C
C     LDXU         INTEGER SCALAR (INPUT)
C                  THE LEADING DIMENSION OF THE ARRAY U EXACTLY AS
C                  SPECIFIED IN THE CALLING PROGRAM.
C
C     IL, IR       INTEGER SCALARS (INPUT)
C     JL, JR       GIVES THE SUBSET OF THE ARRAY U WHICH CONTAINS THE
C                  SEQUENCES TO BE TRANSFORMED, I.E., POSITIONS I=
C                  IL,..,IR OF COLUMNS J=JL,..,JR.
C
C     WSAVE        REAL ARRAY OF SIZE 3*N + N/2 + 15 (INPUT)
C                  THE WORK ARRAY EXACTLY AS RETURNED FROM THE ROUTINE
C                  FFTI.
C
C
C  ------------
C  DECLARATIONS
C  ------------
C
C     ... PARAMETERS
C
      INTEGER KBCL, KBCR, IL, IR, JL, JR
      REAL
     *     U(LDXU,*), WSAVE(*)
C
C     ... LOCAL VARIABLES
C
      INTEGER ICASE
C
C     ... LOCAL CONSTANTS
C
      INTEGER DRCH, NEUM, PRDC, LEFT, RIGHT, TOP, BOTTOM
      PARAMETER (DRCH=1, NEUM =2, PRDC=3,
     *           LEFT=3, RIGHT=1, TOP =4, BOTTOM=2)
C
      N = IR - IL + 1
C
C  ---------------------------------------------
C  SELECT TRANSFORM BASED ON BOUNDARY CONDITIONS
C  ---------------------------------------------
C
      ICASE = 2*(KBCL-1) + KBCR
      IF (KBCL .EQ. PRDC)  ICASE = 5
      GO TO (100,200,300,400,500), ICASE
C
C     CASE :  DIRICHLET/DIRICHLET
C
  100 CONTINUE
      DO 150 J=JL,JR
         CALL SINT(N,U(IL,J),WSAVE)
  150 CONTINUE
      GO TO 600
C
C     CASE :  DIRICHLET/NEUMANN
C
  200 CONTINUE
      DO 250 J=JL,JR
         CALL SINQF(N,U(IL,J),WSAVE)
  250 CONTINUE
      GO TO 600
C
C     CASE :  NEUMANN/DIRICHLET
C
  300 CONTINUE
      DO 350 J=JL,JR
         CALL COSQF(N,U(IL,J),WSAVE)
  350 CONTINUE
      GO TO 600
C
C     CASE :  NEUMANN/NEUMANN
C
  400 CONTINUE
      DO 450 J=JL,JR
         CALL COST(N,U(IL,J),WSAVE)
  450 CONTINUE
      GO TO 600
C
C     CASE :  PERIODIC
C
  500 CONTINUE
      DO 550 J=JL,JR
         CALL RFFTF(N,U(IL,J),WSAVE)
  550 CONTINUE
C
C
C  ----
C  EXIT
C  ----
C
  600 CONTINUE
      RETURN
      END
      SUBROUTINE FFTI (KBCL, KBCR, N, WSAVE)
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C    ---------------      4TH ORDER ACCURATE FAST DIRECT SOLUTION
C    PACKAGE :  HFFT      OF THE HELMHOLTZ EQUATION ON RECTANGULAR
C    ---------------      DOMAINS IN TWO AND THREE DIMENSIONS
C
C    INTERNAL MODULE
C
C    RONALD F. BOISVERT
C    NATIONAL BUREAU OF STANDARDS
C    DECEMBER 1985
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C
C   FFTI INITIALIZES THE ONE-DIMENSIONAL FOURIER TRANSFORM SOFTWARE.
C   THE TYPE OF TRANSFORM SELECTED DEPENDS UPON THE BOUNDARY CONDITIONS.
C
C
C   P A R A M E T E R S
C   -------------------
C
C     KBCL, KBCR   INTEGER SCALARS (INPUT)
C                  GIVE THE TYPE OF BOUNDARY CONDITIONS AT THE LEFT AND
C                  RIGHT ENDPOINTS OF THE INTERVAL. POSSIBLE VALUES ARE
C                  THE SAME AS IN THE VECTOR BCTY OF HFFT2A.
C
C     N            INTEGER SCALAR (INPUT)
C                  THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED.
C
C     WSAVE        REAL ARRAY OF SIZE 3*N + N/2 + 15 (OUTPUT)
C                  CONTAINS INFORMATION WHICH MUST BE PASSED TO THE
C                  SUBROUTINES FFTF AND FFTB WHEN ACTUALLY PERFORMING
C                  THE TRANSFORMS.
C
C
C  ------------
C  DECLARATIONS
C  ------------
C
C     ... PARAMETERS
C
      INTEGER KBCL, KBCR, N
      REAL
     *     WSAVE(*)
C
C     ... LOCAL VARIABLES
C
      INTEGER ICASE
C
C     ... LOCAL CONSTANTS
C
      INTEGER DRCH, NEUM, PRDC, LEFT, RIGHT, TOP, BOTTOM
      PARAMETER (DRCH=1, NEUM =2, PRDC=3,
     *           LEFT=3, RIGHT=1, TOP =4, BOTTOM=2)
C
C
C  ---------------------------------------------
C  SELECT TRANSFORM BASED ON BOUNDARY CONDITIONS
C  ---------------------------------------------
C
      ICASE = 2*(KBCL-1) + KBCR
      IF (KBCL .EQ. PRDC)  ICASE = 5
      GO TO (110,120,130,140,150), ICASE
C
C     CASE :  DIRICHLET/DIRICHLET
C
  110 CONTINUE
      CALL SINTI(N,WSAVE)
      GO TO 200
C
C     CASE :  DIRICHLET/NEUMANN
C
  120 CONTINUE
      CALL SINQI(N,WSAVE)
      GO TO 200
C
C     CASE :  NEUMANN/DIRICHLET
C
  130 CONTINUE
      CALL COSQI(N,WSAVE)
      GO TO 200
C
C     CASE :  NEUMANN/NEUMANN
C
  140 CONTINUE
      CALL COSTI(N,WSAVE)
      GO TO 200
C
C     CASE :  PERIODIC
C
  150 CONTINUE
      CALL RFFTI(N,WSAVE)
      GO TO 200
C
C
C  ----
C  EXIT
C  ----
C
  200 CONTINUE
      RETURN
      END
      SUBROUTINE HDIS2 (NX, NY, H, COEFU, GH, LMXGH, BCTY,
     *                  BD1, BD2, BD3, BD4, A, B, C, U, LMXU, WORK)
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C    ---------------      4TH ORDER ACCURATE FAST DIRECT SOLUTION
C    PACKAGE :  HFFT      OF THE HELMHOLTZ EQUATION ON RECTANGULAR
C    ---------------      DOMAINS IN TWO AND THREE DIMENSIONS
C
C    INTERNAL MODULE
C
C    RONALD F. BOISVERT
C    NATIONAL BUREAU OF STANDARDS
C    DECEMBER 1985
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C
C   HDIS2 COMPUTES A FOURTH ORDER FINITE DIFFERENCE DISCRETIZATION
C   FOR A TWO-DIMENSIONAL RECTANGULAR DOMAIN
C
C
C   P A R A M E T E R S
C   -------------------
C
C     NX, NY    INTEGER SCALARS (INPUT)
C               SEE HFFT2A.
C
C     H, COEFU  REAL SCALARS (INPUT)
C               SEE HFFT2A.
C
C     GH        REAL ARRAY OF SIZE LMXGH+1 BY NY+1 (INPUT)
C               SEE HFFT2A.
C
C     LMXGH     INTEGER SCALAR (INPUT)
C               UPPER LIMIT OF FIRST DIMENSION OF ARRAY GH.  MUST BE
C               LDXGH-1 (LDXGH IS DEFINED IN HFFT2A).
C
C     BCTY      INTEGER ARRAY OF SIZE 4 (INPUT)
C               SEE HFFT2A.
C
C     BD1, BD3  REAL ARRAYS OF SIZE NY (INPUT)
C               SEE HFFT2A.
C
C     BD2, BD4  REAL ARRAYS OF SIZE NX (INPUT)
C               SEE HFFT2A.
C
C     A, B, C   REAL SCALARS (OUTPUT)
C               GIVES VALUES IN THE BASIC FINITE DIFFERENCE STENCIL
C               (SCALED TO O(1))
C
C                        C  B  C
C                        B  A  B  U  =  RIGHT SIDE
C                        C  B  C
C
C     U         REAL ARRAY OF SIZE LMXU+1 BY NY+2 (INPUT/OUTPUT)
C               ON INPUT, U(I,J) IS THE RIGHT HAND SIDE OF THE PDE
C               EVALUATED AT THE (I,J)TH GRID POINT.
C               ON OUTPUT, U(I,J) IS THE RIGHT HAND SIDE OF THE
C               DISCRETE PDE AT THE (I,J)TH GRID POINT.
C
C     LMXU      INTEGER SCALAR (INPUT)
C               UPPER LIMIT OF FIRST DIMENSION OF ARRAY U.  MUST BE
C               LDXU-1 (LDXU IS DEFINED IN HFFT2A).
C
C     WORK      REAL ARRAY OF SIZE NX+1 BY 2
C               WORKING STORAGE FOR HDIS2.
C
C
C     ******************************************************************
C     *                                                                *
C     *   NOTE -- THE ARRAYS U AND GH ARE INDEXED DIFFERENTLY IN       *
C     *           THIS ROUTINE:  U(0:LMXU,0:*) AND GH(0:LMXGH,0:*)     *
C     *                                                                *
C     ******************************************************************
C
C
C  ------------
C  DECLARATIONS
C  ------------
C
C     ... PARAMETERS
C
      INTEGER BCTY(4), NX, NY, LMXGH, LMXU
      REAL
     *     BD1(NY), BD2(NX), BD3(NY), BD4(NX),
     *     GH(0:LMXGH,0:*), U(0:LMXU,0:*), WORK(0:NX,0:1),
     *     COEFU, H, A, B, C
C
C     ... LOCAL VARIABLES
C
      LOGICAL PRDX, PRDY, HELMHZ, HAVED, HAVEN
      INTEGER I, J, JM1, NXM1, NYM1, C0, CM1
      REAL
     *     BETA0, BETA1, BETA2, DELTA0, DELTA1, DELTA2,
     *     DELTA3, GAMMA0, GAMMA1, H2, F, F2
C
      COMMON /HD2COM/ GAMMA0, GAMMA1, DELTA0, DELTA1, DELTA2, DELTA3
C
C     ... LOCAL CONSTANTS
C
      INTEGER DRCH, NEUM, PRDC, LEFT, RIGHT, TOP, BOTTOM
      PARAMETER (DRCH=1, NEUM =2, PRDC=3,
     *           LEFT=3, RIGHT=1, TOP =4, BOTTOM=2)
C
C
C  ---------------
C  INITIALIZATIONS
C  ---------------
C
      LDXU = LMXU + 1
      NXM1 = NX - 1
      NYM1 = NY - 1
C
      HELMHZ = COEFU .NE. 0.0E0
      HAVED = (BCTY(RIGHT) .EQ. DRCH) .OR. (BCTY(BOTTOM) .EQ. DRCH) .OR.
     *        (BCTY(LEFT ) .EQ. DRCH) .OR. (BCTY(TOP   ) .EQ. DRCH)
      HAVEN = (BCTY(RIGHT) .EQ. NEUM) .OR. (BCTY(BOTTOM) .EQ. NEUM) .OR.
     *        (BCTY(LEFT ) .EQ. NEUM) .OR. (BCTY(TOP   ) .EQ. NEUM)
      PRDX = BCTY(RIGHT ) .EQ. PRDC
      PRDY = BCTY(BOTTOM) .EQ. PRDC
C
      H2 = H*H
      F = -H2*COEFU
      F2 = F*F
      A = -(480.0E0 + 118.0E0*F + 5.0E0*F2)/24.0E0
      B =  (192.0E0 -   8.0E0*F +       F2)/48.0E0
      C =  ( 48.0E0 -   5.0E0*F)/48.0E0
      BETA0 = (48.0E0 + 5.0E0*F)/24.0E0
      BETA1 = -F/48.0E0
      BETA2 = 1.0E0
      GAMMA0 = -(12.0E0 + 13.0E0*F/12.0E0)
      GAMMA1 = -F/12.0E0
      DELTA2 = 2.0E0 + F/48.0E0
      DELTA3 = ((6.0E0 + F)/3.0E0 - 4.0E0*DELTA2)/18.0E0
      DELTA1 = 8.0E0*DELTA2 + 33.0E0*DELTA3 - (10.0E0 + F)*BETA2
      DELTA0 = -2.0E0*(BETA0 + 2.0E0*BETA1) - (DELTA2 + 4.0E0*DELTA3)
     *            - (4.0E0 + F*0.50E0)*BETA2
C
      BETA0 = H2*BETA0
      BETA1 = H2*BETA1
      BETA2 = H2*BETA2
      GAMMA0 = H*GAMMA0
      GAMMA1 = H*GAMMA1
      DELTA0 = H*DELTA0
      DELTA1 = H*DELTA1
      DELTA2 = H*DELTA2
      DELTA3 = H*DELTA3
C
C
C  -----------------------------
C  SHIFT VALUES OF G STORED IN U
C  -----------------------------
C
      DO 50 J=NY,1,-1
         DO 50 I=NX,1,-1
            U(I,J) = U(I-1,J-1)
   50 CONTINUE
C
C
C  -----------------------------------------
C  REFLECT FUNCTIONS G AND GH OUTSIDE DOMAIN
C  -----------------------------------------
C
      IF (HELMHZ)  CALL REFL2(1,NX,NY,PRDX,PRDY,U,LMXU)
      CALL REFL2(0,NXM1,NYM1,PRDX,PRDY,GH,LMXGH)
C
C
C  ----------------------------
C  DISCRETIZE RIGHT SIDE OF PDE
C  ----------------------------
C
      ISTRT = 2
      ISTOP = NXM1
      IF (BCTY(LEFT  ) .NE. DRCH)  ISTRT = 1
      IF (BCTY(RIGHT ) .EQ. NEUM)  ISTOP = NX
      JSTRT = 2
      JSTOP = NYM1
      IF (BCTY(BOTTOM) .NE. DRCH)  JSTRT = 1
      IF (BCTY(TOP   ) .EQ. NEUM)  JSTOP = NY
C
      IF (HELMHZ) THEN
C
C        CASE :  HELMHOLTZ EQUATION
C
         C0 = 0
         JM1 = JSTRT - 1
         DO 60 I=0,NX
            WORK(I,C0) = U(I,JM1)
   60    CONTINUE
         DO 100 J=JSTRT,JSTOP
            CM1 = C0
            C0  = 1 - CM1
            DO 80 I=0,NX
               WORK(I,C0) = U(I,J)
   80       CONTINUE
            DO 100 I=ISTRT,ISTOP
               U(I,J) =   BETA0*U(I,J)
     *                  + BETA1*( U(I+1,J) + U(I,J+1)
     *                          + WORK(I,CM1) + WORK(I-1,C0) )
     *                  + BETA2*( GH(I,J) + GH(I-1,J) +
     *                            GH(I-1,J-1) + GH(I,J-1) )
  100    CONTINUE
C
      ELSE
C
C        CASE :  POISSON EQUATION
C
         DO 200 J=JSTRT,JSTOP
            DO 200 I=ISTRT,ISTOP
               U(I,J) =  BETA0*U(I,J)
     *                 + BETA2*( GH(I,J) + GH(I-1,J) +
     *                           GH(I-1,J-1) + GH(I,J-1) )
  200    CONTINUE
C
      ENDIF
C
C
C  -------------------------
C  REMOVE SHIFT FROM ARRAY U
C  -------------------------
C
      DO 250 J=0,NYM1
         DO 250 I=0,NXM1
            U(I,J) = U(I+1,J+1)
  250 CONTINUE
C
C
C  ---------------------------------------
C  DISCRETIZE POINTS ON NEUMANN BOUNDARIES
C  ---------------------------------------
C
      IF (HAVEN)  CALL HD2N(NX,NY,BCTY,BD1,BD2,BD3,BD4,U,LDXU)
C
C
C  -----------------------------------------------
C  ADJUST POINTS ADJACENT TO DIRICHLET BOUNDARIES
C  -----------------------------------------------
C
      IF (HAVED) THEN
         CALL STORD2(NX,NY,BCTY,BD1,BD2,BD3,BD4,U,LDXU)
         CALL FD2D(NX,NY,U,LDXU,BCTY,B,C,U)
      ENDIF
C
C
C  ----
C  EXIT
C  ----
C
      RETURN
      END
      SUBROUTINE HD2N (NX, NY, BCTY, BD1, BD2, BD3, BD4, U, LDXU)
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C    ---------------      4TH ORDER ACCURATE FAST DIRECT SOLUTION
C    PACKAGE :  HFFT      OF THE HELMHOLTZ EQUATION ON RECTANGULAR
C    ---------------      DOMAINS IN TWO AND THREE DIMENSIONS
C
C    INTERNAL MODULE
C
C    RONALD F. BOISVERT
C    NATIONAL BUREAU OF STANDARDS
C    DECEMBER 1985
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C
C  HD2N COMPUTES THE FOURTH ORDER FINITE DIFFERENCE DISCRETIZATION
C  AT ALL BOUNDARY POINTS OF A TWO-DIMENSIONAL RECTANGULAR DOMAIN
C  WHERE NEUMANN BOUNDARY CONDITIONS HAVE BEEN SPECIFIED.
C
C
C   P A R A M E T E R S
C   -------------------
C
C     NX, NY    INTEGER SCALARS (INPUT)
C               SEE HFFT2A.
C
C     BCTY      INTEGER ARRAY OF SIZE 4 (INPUT)
C               SEE HFFT2A.
C
C     BD1, BD3  INTEGER ARRAYS OF SIZE NY (INPUT)
C               SEE HFFT2A.
C
C     BD2, BD4  INTEGER ARRAYS OF SIZE NX (INPUT)
C               SEE HFFT2A.
C
C     U         REAL ARRAY OF SIZE LDXU BY NY (OUTPUT)
C               ON EXIT, ENTRIES CORRESPONDING TO NEUMANN BOUNDARY
C               POINTS CONTAIN THE RIGHT HAND SIDE OF THE FINITE
C               DIFFERENCE DISCRETIZTION
C
C     LDXU      INTEGER SCALAR (INPUT)
C               SEE HFFT2A.
C
C
C  ------------
C  DECLARATIONS
C  ------------
C
C     ... PARAMETERS
C
      INTEGER BCTY(4), NX, NY, LDXU
      REAL
     *     BD1(NY), BD2(NX), BD3(NY), BD4(NX),
     *     U(LDXU,NY)
C
C     ... LOCAL VARIABLES
C
      LOGICAL PRDX, PRDY
      INTEGER I, J, NXM1, NYM1
      REAL
     *       DELTA0, DELTA1, DELTA2, DELTA3, GAMMA0, GAMMA1
C
      COMMON /HD2COM/ GAMMA0, GAMMA1, DELTA0, DELTA1, DELTA2, DELTA3
C
C     ... LOCAL CONSTANTS
C
      INTEGER DRCH, NEUM, PRDC, LEFT, RIGHT, TOP, BOTTOM
      PARAMETER (DRCH=1, NEUM =2, PRDC=3,
     *           LEFT=3, RIGHT=1, TOP =4, BOTTOM=2)
C
C
C  ---------------
C  INITIALIZATIONS
C  ---------------
C
      NXM1 = NX - 1
      NYM1 = NY - 1
      NXM2 = NX - 2
      NYM2 = NY - 2
      NXM3 = NX - 3
      NYM3 = NY - 3
      PRDX = BCTY(RIGHT ) .EQ. PRDC
      PRDY = BCTY(BOTTOM) .EQ. PRDC
C
C
C  ----------------------------
C  NEUMANN POINT DISCRETIZATION
C  ----------------------------
C
C     ... BOTTOM EDGE
C
      IF (BCTY(BOTTOM) .EQ. NEUM) THEN
         DO 305 I=2,NXM1
            U(I,1) =  U(I,1)
     *                - GAMMA0*BD2(I)
     *                - GAMMA1*(BD2(I-1) + BD2(I+1))
  305    CONTINUE
      ENDIF
C
C     ... TOP EDGE
C
      IF (BCTY(TOP) .EQ. NEUM) THEN
         DO 315 I=2,NXM1
            U(I,NY) =  U(I,NY)
     *                 + GAMMA0*BD4(I)
     *                 + GAMMA1*(BD4(I-1)+BD4(I+1))
  315    CONTINUE
      ENDIF
C
C     ... LEFT EDGE
C
      IF (BCTY(LEFT) .EQ. NEUM) THEN
         DO 325 J=2,NYM1
            U(1,J) =  U(1,J)
     *                - GAMMA0*BD3(J)
     *                - GAMMA1*(BD3(J-1) + BD3(J+1))
  325    CONTINUE
      ENDIF
C
C     ... RIGHT EDGE
C
      IF (BCTY(RIGHT) .EQ. NEUM) THEN
         DO 335 J=2,NYM1
            U(NX,J) =  U(NX,J)
     *                 + GAMMA0*BD1(J)
     *                 + GAMMA1*(BD1(J-1) + BD1(J+1))
  335    CONTINUE
      ENDIF
C
C     ... LOWER RIGHT CORNER
C
      IF ((BCTY(RIGHT) .EQ. NEUM) .AND. (BCTY(BOTTOM) .EQ. NEUM))
     *   U(NX,1) =   U(NX,1)
     *               + DELTA0*BD1(1)
     *               + DELTA1*BD1(2)
     *               + DELTA2*BD1(3)
     *               + DELTA3*BD1(4)
     *               - DELTA0*BD2(NX)
     *               - DELTA1*BD2(NXM1)
     *               - DELTA2*BD2(NXM2)
     *               - DELTA3*BD2(NXM3)
C
C     ... LOWER LEFT CORNER
C
      IF ((BCTY(BOTTOM) .EQ. NEUM) .AND. (BCTY(LEFT) .EQ. NEUM))
     *    U(1,1) =   U(1,1)
     *               - DELTA0*BD2(1)
     *               - DELTA1*BD2(2)
     *               - DELTA2*BD2(3)
     *               - DELTA3*BD2(4)
     *               - DELTA0*BD3(1)
     *               - DELTA1*BD3(2)
     *               - DELTA2*BD3(3)
     *               - DELTA3*BD3(4)
C
C     ... UPPER LEFT CORNER
C
      IF ((BCTY(LEFT) .EQ. NEUM) .AND. (BCTY(TOP) .EQ. NEUM))
     *   U(1,NY) =   U(1,NY)
     *               - DELTA0*BD3(NY)
     *               - DELTA1*BD3(NYM1)
     *               - DELTA2*BD3(NYM2)
     *               - DELTA3*BD3(NYM3)
     *               + DELTA0*BD4(1)
     *               + DELTA1*BD4(2)
     *               + DELTA2*BD4(3)
     *               + DELTA3*BD4(4)
C
C     ... UPPER RIGHT CORNER
C
      IF ((BCTY(TOP) .EQ. NEUM) .AND. (BCTY(RIGHT) .EQ. NEUM))
     *   U(NX,NY) =   U(NX,NY)
     *                + DELTA0*BD1(NY)
     *                + DELTA1*BD1(NYM1)
     *                + DELTA2*BD1(NYM2)
     *                + DELTA3*BD1(NYM3)
     *                + DELTA0*BD4(NX)
     *                + DELTA1*BD4(NXM1)
     *                + DELTA2*BD4(NXM2)
     *                + DELTA3*BD4(NXM3)
C
C
      IF (.NOT. (PRDX .OR. PRDY))  GO TO 999
C
C
C  -------------------------------------
C  NEUMANN/PERIODIC POINT DISCRETIZATION
C  -------------------------------------
C
C     ... LOWER LEFT CORNER (PERIODIC/NEUMANN)
C
      IF (PRDX .AND. (BCTY(BOTTOM) .EQ. NEUM))
     *   U(1,1) =   U(1,1)
     *              - GAMMA0*BD2(1)
     *              - GAMMA1*(BD2(2) + BD2(NXM1))
C
C     ... LOWER LEFT CORNER (NEUMANN/PERIODIC)
C
      IF (PRDY .AND. (BCTY(LEFT) .EQ. NEUM))
     *   U(1,1) =   U(1,1)
     *              - GAMMA0*BD3(1)
     *              - GAMMA1*(BD3(2) + BD3(NYM1))
C
C     ... UPPER LEFT CORNER (PERIODIC/NEUMANN)
C
      IF (PRDX .AND. (BCTY(TOP) .EQ. NEUM))
     *   U(1,NY) =   U(1,NY)
     *               + GAMMA0*BD4(1)
     *               + GAMMA1*(BD4(2) + BD4(NXM1))
C
C     ... LOWER RIGHT CORNER (NEUMANN/PERIODIC)
C
      IF (PRDY .AND. (BCTY(RIGHT) .EQ. NEUM))
     *   U(NX,1) =   U(NX,1)
     *               + GAMMA0*BD1(1)
     *               + GAMMA1*(BD1(2) + BD1(NYM1))
C
C  ----
C  EXIT
C  ----
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE MDALG2 (A, B, C, BCTY, U, LDXU, IL, IR, JL, JR, WORK)
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C    ---------------      4TH ORDER ACCURATE FAST DIRECT SOLUTION
C    PACKAGE :  HFFT      OF THE HELMHOLTZ EQUATION ON RECTANGULAR
C    ---------------      DOMAINS IN TWO AND THREE DIMENSIONS
C
C    INTERNAL MODULE
C
C    RONALD F. BOISVERT
C    NATIONAL BUREAU OF STANDARDS
C
C    DECEMBER 1985  (REVISED APRIL 1987)
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C
C   MDALG2 IMPLEMENTS THE MATRIX DECOMPOSITION ALGORITHM (FOURIER
C   METHOD FOR A NINE-POINT DIFFERENCE OPERATOR ON A TWO-DIMENSIONAL
C   RECTANGULAR GRID.
C
C
C   P A R A M E T E R S
C   -------------------
C
C     A, B, C      REAL SCALARS (INPUT)
C                  GIVE THE BASIC FINITE DIFFERENCE STENCIL THAT IS
C                  USED TO APPROXIMATE THE PDE.
C
C                         C  B  C
C                         B  A  B  U  = RIGHT HAND SIDE
C                         C  B  C
C
C     BCTY         INTEGER ARRAY OF SIZE 4 (INPUT)
C                  SEE HFFT2A.
C
C     U            REAL ARRAY OF SIZE LDXU BY JR (INPUT/OUTPUT)
C                  ON INPUT, U CONTAINS THE RIGHT HAND SIDE OF THE
C                  DISCRETE APPROXIMATION TO THE PDE FOR EACH POINT
C                  POINT AT WHICH THE SOLUTION IS TO BE DETERMINED,
C                  I.E., (I,J), I=IL,..,IR, J=JL,..,JR.
C                  ON OUTPUT THESE VALUES ARE REPLACED BY THE COMPUTED
C                  SOLUTION.
C
C     LDXU         INTEGER SCALAR (INPUT)
C                  THE LEADING DIMENSION OF THE ARRAY U EXACTLY AS
C                  SPECIFIED IN THE CALLING PROGRAM.
C
C     IL, IR,      INTEGER SCALARS (INPUT)
C     JL, JR       GIVES THE SUBSET OF GRID POINTS AT WHICH THE
C                  SOLUTION IS TO BE DETERMINED, I.E., THAT SET
C                  OF INDICES (I,J) WITH I=IL,..,IR AND J=JL,..,JR.
C
C     WORK         REAL ARRAY OF SIZE 5*N + 5*M + N/2 + 15 (WORKSPACE)
C                  HERE N=IR-IL+1 AND M=JR-JL+1. THE LENGTH OF THIS
C                  ARRAY MAY BE REDUCED BY M WHEN THE SOLUTION IS
C                  NOT PERIODIC IN Y.  IT MAY BE REDUCED BY 4*M IF
C                  THE COEFFICIENT OF U IN THE PDE IS .LE. 0.
C
C
C  ------------
C  DECLARATIONS
C  ------------
C
C     ... PARAMETERS
C
      INTEGER BCTY(4), LDXU, IL, IR, JL, JR
      REAL
     *     A, B, C, U(LDXU,*), WORK(*)
C
C     ... LOCAL VARIABLES
C
      INTEGER N, M, LOCEWK, LOCFWK, LOCTWK, TOTAL
C
C     ... LOCAL CONSTANTS
C
      INTEGER DRCH, NEUM, PRDC, LEFT, RIGHT, TOP, BOTTOM
      PARAMETER (DRCH=1, NEUM =2, PRDC=3,
     *           LEFT=3, RIGHT=1, TOP =4, BOTTOM=2)
C
C
C  --------------
C  INITIALIZATION
C  --------------
C
      N = IR - IL + 1
      M = JR - JL + 1
      LOC EWK = 1
      LOC FWK = LOC EWK + N
      LOC TWK = LOC FWK + 3*N + N/2 + 15
      TOTAL   = LOC TWK + N + 5*M - 1
C
      CALL FFTI(BCTY(LEFT),BCTY(RIGHT),N,WORK(LOCFWK))
C
C
C  ------------------
C  FORWARD TRANSFORMS
C  ------------------
C
      CALL FFTF(BCTY(LEFT),BCTY(RIGHT),U,LDXU,IL,IR,JL,JR,WORK(LOCFWK))
C
C
C  ------------------
C  TRIDIAGONAL SOLVES
C  ------------------
C
      CALL TRISOL(A,B,C,BCTY,U,LDXU,IL,IR,JL,JR,WORK(LOCEWK),
     *            WORK(LOCTWK))
C
C
C  -------------------
C  BACKWARD TRANSFORMS
C  -------------------
C
      CALL FFTB(BCTY(LEFT),BCTY(RIGHT),U,LDXU,IL,IR,JL,JR,WORK(LOCFWK))
C
C
C  ----
C  EXIT
C  ----
C
      RETURN
      END
      SUBROUTINE REFL2 (KDIST, NX, NY, PRDX, PRDY, G, LMXG)
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C    ---------------      4TH ORDER ACCURATE FAST DIRECT SOLUTION
C    PACKAGE :  HFFT      OF THE HELMHOLTZ EQUATION ON RECTANGULAR
C    ---------------      DOMAINS IN TWO AND THREE DIMENSIONS
C
C    INTERNAL MODULE
C
C    RONALD F. BOISVERT
C    NATIONAL BUREAU OF STANDARDS
C    DECEMBER 1985
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C
C   REFL2 EXTENDS A TWO-DIMENSIONAL GRID FUNCTION SO THAT IT IS
C   DEFINED ONE GRID LINE OUTSIDE ITS ORGINAL DOMAIN.  THE EXTENTION
C   IS DONE BY REFLECTION THROUGH THE BOUNDARY EXCEPT WHERE
C   PERIODICITY IS SPECIFIED.
C
C
C   P A R A M E T E R S
C   -------------------
C
C     KDIST   INTEGER SCALAR (INPUT)
C             INDICATES HOW THE GRID FUNCTION IS DEFINED.
C             POSSIBLE VALUES ARE
C
C                1 == FUNCTION DEFINED AT CENTER OF GRID SQUARES
C                2 == FUNCTION DEFINED AT GRID POINTS
C
C     NX      INTEGER SCALAR (INPUT)
C             THE NUMBER OF GRID FUNCTION VALUES IN THE X DIRECTION.
C
C     NY      INTEGER SCALAR (INPUT)
C             THE NUMBER OF GRID FUNCTION VALUES IN THE Y DIRECTION.
C
C     PRDX    LOGICAL SCALAR (INPUT)
C             .TRUE. IF THE SOLUTION IS TO BE EXTENDED PERIODICALLY
C             IN X.
C
C     PRDY    LOGICAL SCALAR (INPUT)
C             .TRUE. IF THE SOLUTION IS TO BE EXTENDED PERIODICALLY
C             IN Y.
C
C     G       REAL ARRAY OF SIZE LMXG+1 BY NY+1 (INPUT/OUTPUT)
C             ON INPUT, THE GRID FUNCTION OCCUPIES G(I,J), I=1,..,NX,
C             J=1,..,NY.
C             ON OUTPUT, THE FUNCTION HAS BEEN EXTENDED TO INCLUDE THE
C             POINTS G(0,J), G(NX+1,J), G(I,0), G(I,NY+1), I=0,..,NX+1,
C             J=0,..,NY+1.
C
C     LMXG    INTEGER SCALAR (INPUT)
C             THE UPPER LIMIT OF THE FIRST DIMENSION OF THE ARRAY G.
C             MUST BE SET TO LDXU-1, WHERE LDXU IS THE ACTUAL LENGTH OF
C             THE FIRST DIMENSION OF G AS DECLARED IN THE CALLING
C             PROGRAM.
C
C
C  ------------
C  DECLARATIONS
C  ------------
C
C     ... PARAMETERS
C
      LOGICAL PRDX, PRDY
      INTEGER KDIST, NX, NY, LMXG
      REAL
     *     G(0:LMXG,0:*)
C
C
C  ---------------
C  INITIALIZATIONS
C  ---------------
C
      NXP1 = NX + 1
      NYP1 = NY + 1
C
      I0 = 1 + KDIST
      J0 = 1 + KDIST
      I1 = NX - KDIST
      J1 = NY - KDIST
C
C
C  -----------
C  REFLECTIONS
C  -----------
C
C     ... IN Y DIRECTION
C
      DO 10 I=1,NX
         G(I,0) = G(I,J0)
         G(I,NYP1) = G(I,J1)
   10 CONTINUE
C
C     ... IN X DIRECTION
C
      DO 20 J=0,NYP1
         G(0,J) = G(I0,J)
         G(NXP1,J) = G(I1,J)
   20 CONTINUE
C
C
C  -------------------
C  PERIODIC EXTENSIONS
C  -------------------
C
C     ... IN X DIRECTION
C
      IF (PRDX) THEN
         DO 30 J=0,NYP1
            G(0,J) = G(NXP1,J)
   30     CONTINUE
      ENDIF
C
C     ... IN Y DIRECTION
C
      IF (PRDY) THEN
         DO 40 I=0,NXP1
            G(I,0) = G(I,NYP1)
   40    CONTINUE
      ENDIF
C
C
C  ----
C  EXIT
C  ----
C
      RETURN
      END
      SUBROUTINE STORD2 (NX, NY, BCTY, BD1, BD2, BD3, BD4, U, LDXU)
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C    ---------------      4TH ORDER ACCURATE FAST DIRECT SOLUTION
C    PACKAGE :  HFFT      OF THE HELMHOLTZ EQUATION ON RECTANGULAR
C    ---------------      DOMAINS IN TWO AND THREE DIMENSIONS
C
C    INTERNAL MODULE
C
C    RONALD F. BOISVERT
C    NATIONAL BUREAU OF STANDARDS
C    DECEMBER 1985
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C
C   STORD2 STORES GIVEN DIRICHLET BOUNDARY DATA IN THE SOLUTION ARRAY U.
C
C
C   P A R A M E T E R S
C   -------------------
C
C     NX, NY    INTEGER SCALARS (INPUT)
C               SEE HFFT2A.
C
C     BCTY      INTEGER ARRAY OF SIZE 4 (INPUT)
C               SEE HFFT2A.
C
C     BD1, BD3  REAL ARRAYS OF SIZE NY (INPUT)
C               SEE HFFT2A.
C
C     BD2, BD4  REAL ARRAYS OF SIZE NX (INPUT)
C               SEE HFFT2A.
C
C     U        REAL ARRAY OF SIZE LDXU BY NY (OUTPUT)
C              ON OUTPUT, ENTRIES CORRESPONDING TO DIRICHLET BOUNDARY
C              POINTS CONTAIN THE KNOWN VALUES OF THE SOLUTION.
C
C     LDXU     INTEGER SCALAR (INPUT)
C              SEE HFFT2A.
C
C
C  ------------
C  DECLARATIONS
C  ------------
C
C     ... PARAMETERS
C
      INTEGER NX, NY, BCTY(4), LDXU
      REAL
     *     BD1(NY), BD2(NX), BD3(NY), BD4(NX), U(LDXU,NY)
C
C     ... LOCAL CONSTANTS
C
      INTEGER DRCH, NEUM, PRDC, LEFT, RIGHT, TOP, BOTTOM
      PARAMETER (DRCH=1, NEUM =2, PRDC=3,
     *           LEFT=3, RIGHT=1, TOP =4, BOTTOM=2)
C
C
C  ----------------------
C  HANDLE DIRICHLET SIDES
C  ----------------------
C
C     ... RIGHT SIDE
C
      IF (BCTY(RIGHT) .EQ. DRCH) THEN
         DO 210 J = 1,NY
            U(NX,J) = BD1(J)
  210    CONTINUE
      ENDIF
C
C     ... BOTTOM SIDE
C
      IF (BCTY(BOTTOM) .EQ. DRCH) THEN
         DO 220 I=1,NX
            U(I,1) = BD2(I)
  220    CONTINUE
      ENDIF
C
C     ... LEFT SIDE
C
      IF (BCTY(LEFT) .EQ. DRCH) THEN
         DO 230 J=1,NY
            U(1,J) = BD3(J)
  230    CONTINUE
      ENDIF
C
C     ... TOP SIDE
C
      IF (BCTY(TOP) .EQ. DRCH) THEN
         DO 240 I=1,NX
            U(I,NY) = BD4(I)
  240    CONTINUE
      ENDIF
C
C  ------------------------
C  HANDLE DIRICHLET CORNERS
C  ------------------------
C
      IF ((BCTY(TOP) .EQ. DRCH) .AND. (BCTY(RIGHT) .EQ. DRCH))
     *   U(NX,NY) = 0.50E0*( BD4(NX) + BD1(NY) )
      IF ((BCTY(RIGHT) .EQ. DRCH) .AND. (BCTY(BOTTOM) .EQ. DRCH))
     *   U(NX,1) = 0.50E0*( BD1(1) + BD2(NX) )
      IF ((BCTY(BOTTOM) .EQ. DRCH) .AND. (BCTY(LEFT) .EQ. DRCH))
     *   U(1,1) = 0.50E0*( BD2(1) + BD3(1) )
      IF ((BCTY(LEFT) .EQ. DRCH) .AND. (BCTY(TOP) .EQ. DRCH))
     *   U(1,NY) = 0.50E0*( BD3(NY) + BD4(1) )
C
C  ----
C  EXIT
C  ----
C
      RETURN
      END
      SUBROUTINE TRISOL (A, B, C, BCTY, U, LDXU, IL, IR, JL, JR,
     *                   EIGEN, WORK)
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C    ---------------      4TH ORDER ACCURATE FAST DIRECT SOLUTION
C    PACKAGE :  HFFT      OF THE HELMHOLTZ EQUATION ON RECTANGULAR
C    ---------------      DOMAINS IN TWO AND THREE DIMENSIONS
C
C    INTERNAL MODULE
C
C    RONALD F. BOISVERT
C    JOAN M. BAUMANN
C    NATIONAL BUREAU OF STANDARDS
C
C    DECEMBER 1985   (REVISED APRIL 1987)
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C
C   TRISOL SOLVES THE SET OF TRIDIAGONAL LINEAR SYSTEMS OF EQUATIONS
C   OBTAINED FROM THE MATRIX DECOMPOSITION ALGORITHM (FOURIER METHOD).
C
C   N TRIDIAGONAL SYSTEMS OF ORDER M ARE SOLVED, WHERE N=IR-IL+1
C   AND M=JR-JL+1.  EACH SYSTEM IS SYMMETRIC WITH CONSTANT DIAGONALS
C   (THE TOP LEFT AND BOTTOM RIGHT ELEMENTS MAY BE MULTIPLIED BY 1/2,
C   DEPENDING UPON THE BOUNDARY CONDITIONS).  THE MATRIX COEFFICIENTS
C   DEPEND UPON THE FINITE DIFFERENCE STENCIL COEFFICIENTS A, B, C, AND
C   THE BOUNDARY CONDITIONS GIVEN BY BCTY.  THE RIGHT HAND SIDES OF THE
C   SYSTEMS ARE STORED IN THE ROWS OF U.
C
C
C   P A R A M E T E R S
C   -------------------
C
C     A, B, C     REAL SCALARS (INPUT)
C                 SEE MDALG2.
C
C     BCTY        INTEGER ARRAY OF SIZE 4 (INPUT)
C                 SEE HFFT2A.
C
C     U           REAL ARRAY OF SIZE LDXU BY M (INPUT/OUTPUT)
C                 ON INPUT THE RIGHT HAND SIDES OF THE LINEAR SYSTEMS
C                 TO BE SOLVED ARE STORED IN THE ROWS OF U, I.E., THE
C                 ITH RIGHT HAND SIDE IS IN U(I,J), J=JL,..,JR, FOR
C                 I=IL,..,IR.  ON OUTPUT THESE ARE REPLACED BY THE
C                 SOLUTIONS OF THE LINEAR SYSTEMS.
C
C     LDXU        INTEGER SCALAR (INPUT)
C                 THE LEADING DIMENSION OF THE ARRAY U EXACTLY AS
C                 SPECIFIED IN THE CALLING PROGRAM.
C
C     IL, IR,     INTEGER SCALARS (INPUT)
C     JL, JR      SPECIFY THE SUBARRAY OF U IN WHICH THE RIGHT HAND
C                 SIDES OF THE SYSTEMS ARE STORED.  SEE U ABOVE.
C
C     EIGEN       REAL ARRAY OF SIZE N (WORKSPACE)
C
C     WORK        REAL ARRAY OF SIZE 5*M + N. THIS LENGTH MAY BE REDUCED
C                 BY 4*M IF THE COEFFICIENT OF U IN THE PDE IS .LE. 0.
C
C
C  ------------
C  DECLARATIONS
C  ------------
C
C     ... PARAMETERS
C
      INTEGER BCTY(4), LDXU, IL, IR, JL, JR
      REAL
     *     A, B, C, U(LDXU,*), EIGEN(*), WORK(*)
C
C     ... LOCAL VARIABLES
C
      INTEGER N, M
      REAL
     *     MU, NU
C
C     ... LOCAL CONSTANTS
C
      INTEGER DRCH, NEUM, PRDC, LEFT, RIGHT, TOP, BOTTOM
      PARAMETER (DRCH=1, NEUM =2, PRDC=3,
     *           LEFT=3, RIGHT=1, TOP =4, BOTTOM=2)
C
C
C  ---------------
C  INITIALIZATIONS
C  ---------------
C
      N = IR - IL + 1
      M = JR - JL + 1
C
      MU = 1.0E0
      NU = 1.0E0
      IF (BCTY(BOTTOM) .EQ. NEUM)  MU = 0.5E0
      IF (BCTY(TOP   ) .EQ. NEUM)  NU = 0.5E0
      CALL EVDISC(BCTY(LEFT),BCTY(RIGHT),EIGEN,N)
C
C
C  -------------------------
C  SOLVE TRIDIAGONAL SYSTEMS
C  -------------------------
C
C     ... CALCULATE THE CONSTANT DIAGONALS AND OFF-DIAGONALS OF THE
C         SYSTEMS TO BE SOLVED, STORING THEM IN THE ARRAYS EIGEN
C         AND WORK, RESPECTIVELY
C
      DO 100 K=1,N
         WORK(K)  = B + C*EIGEN(K)
         EIGEN(K) = A + B*EIGEN(K)
  100 CONTINUE
C
C     ... SOLVE SYSTEMS
C
      CALL TRSALL (EIGEN, WORK, U, LDXU, MU, NU, IL, IR, JL, JR,
     *             BCTY(BOTTOM).EQ.PRDC, WORK(N+1) )
C
C
C  ----
C  EXIT
C  ----
C
      RETURN
      END
      SUBROUTINE TRSALL (DIAG, OFFDG, U, LDXU, MU, NU, IL, IR, JL, JR,
     *                   PRDC, WORK )
C
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C    ---------------      4TH ORDER ACCURATE FAST DIRECT SOLUTION
C    PACKAGE :  HFFT      OF THE HELMHOLTZ EQUATION ON RECTANGULAR
C    ---------------      DOMAINS IN TWO AND THREE DIMENSIONS
C
C    INTERNAL MODULE
C
C    RONALD F. BOISVERT
C    JOAN M. BAUMANN
C    NATIONAL BUREAU OF STANDARDS
C
C    DECEMBER 1985   (REVISED APRIL 1987)
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C
C   TRSALL SOLVES N TRIDIAGONAL SYSTEMS OF LINEAR EQUATIONS OF SIZE
C   M, WHERE N=IR=IL+1 AND M=JR-JL+1.  ALL OF THE SYSTEMS AU=G MUST
C   HAVE MATRICES A OF ONE OF THE FOLLOWING FORMS.
C
C                      :--                      --:
C                      :  MU*A  B                 :
C                      :     B  A  B              :
C                      :        B  A  B           :
C               A  =   :          .  .  .         :
C                      :            .  .  .       :
C                      :              B  A  B     :
C                      :                 B  A*NU  :
C                      :--                      --:
C
C   AND MU AND NU HAVING THE VALUES 0.5 OR 1.0, OR
C
C                      :--                      --:
C                      :     A  B           B     :
C                      :     B  A  B              :
C                      :        B  A  B           :
C               A  =   :          .  .  .         :
C                      :            .  .  .       :
C                      :              B  A  B     :
C                      :     B           B  A     :
C                      :--                      --:
C
C
C   ALL SYSTEMS MUST BE OF THE SAME FORM, ALTHOUGH THE SCALARS A AND B
C   MAY BE DIFFERENT FOR EACH SYSTEM.
C
C   ONE OF TWO ALGORITHMS IS USED DEPENDING UPON THE VALUES OF A AND B.
C   IF ABS(A) .GE. ABS(2*B) THEN AN EXTENSION OF AN ALGORITHM DUE TO EVA
C   IS USED; OTHERWISE GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING FOR
C   GENERAL TRIDIAGONAL MATRICES IS USED.
C
C
C   P A R A M E T E R S
C   -------------------
C
C     DIAG        REAL ARRAY OF SIZE N (INPUT)
C                 CONSTANT DIAGONALS FOR THE SYSTEMS TO BE SOLVED.
C                 THE DIAGONAL ELEMENT A FOR THE I-TH SYSTEM IS STORED
C                 IN DIAG(K).
C
C     OFFDG       REAL ARRAY OF SIZE N (INPUT)
C                 CONSTANT OFF-DIAGONALS FOR THE SYSTEMS TO BE SOLVED.
C                 THE OFF-DIAGONAL ELEMENT A FOR THE I-TH SYSTEM IS
C                 STORED IN OFFDG(I).
C
C     U           REAL ARRAY OF SIZE LDXU BY M (INPUT/OUTPUT)
C                 ON INPUT THE RIGHT HAND SIDES OF THE LINEAR  SYSTEMS
C                 TO BE SOLVED ARE STORED IN THE ROWS OF U, I.E., THE
C                 ITH RIGHT HAND SIDE IS IN U(I,J), J=JL,..,JR, FOR
C                 I=IL,..,IR.
C                 ON OUTPUT THESE ROWS ARE REPLACED BY THE SOLUTIONS
C                 OF THE LINEAR SYSTEMS.
C
C     LDXU        INTEGER SCALAR (INPUT)
C                 THE LEADING DIMENSION OF THE ARRAY U EXACTLY AS
C                 SPECIFIED IN THE CALLING PROGRAM.
C
C     MU,NU       REAL  SCALARS (INPUT)
C                 MULTIPLICATIVE FACTORS FOR THE FIRST AND LAST
C                 DIAGONAL ELEMENTS, RESPECTIVELY, FOR EACH OF THE
C                 LINEAR SYSTEMS.  THESE ARE IGNORED IF PRDC=.TRUE.
C
C     IL, IR,     INTEGER SCALARS (INPUT)
C     JL, JR      SPECIFY THE SUBARRAY OF U IN WHICH THE RIGHT HAND
C                 SIDES OF THE SYSTEMS ARE STORED.  SEE U ABOVE.
C
C     PRDC        LOGICAL VARIABLE (INPUT)
C                 INDICATES THE FORM OF THE MATRICES.  IS .TRUE. FOR
C                 THE PERIODIC CASE (SECOND CASE ABOVE), AND .FALSE.
C                 OTHERWISE.
C
C     WORK        REAL WORK ARRAY OF SIZE 5*M. THIS CAN BE REDUCED TO
C                 M IF ABS(DIAG(K)) .GE. 2*ABS(OFFDG(K)) FOR ALL K.
C
C
C  ------------
C  DECLARATIONS
C  ------------
C
C     ... PARAMETERS
C
      INTEGER  LDXU, IL, IR, JL, JR
      REAL DIAG(*), OFFDG(*), U(LDXU,*), WORK(*), MU, NU
      LOGICAL PRDC
C
C     ... LOCAL VARIABLES
C
      INTEGER K, KX, IPOS, INFO
      REAL CORNER
C
C     ... LOCAL CONSTANTS
C
      INTEGER M, LOCC, LOCD, LOCE, LOCW
C
      M = JR - JL + 1
      LOCC =    1 + M
      LOCD = LOCC + M
      LOCE = LOCD + M
      LOCW = LOCE + M
C
C
C  -------------------------
C  SOLVE TRIDIAGONAL SYSTEMS
C  -------------------------
C
      K = 0
      DO 500 I=IL,IR
         K = K + 1
C
C        ... CONSTRUCT RIGHT HAND SIDE
C
         KX = 0
         DO 100 J=JL,JR
            KX = KX + 1
            WORK(KX) = U(I,J)
  100    CONTINUE
         WORK(1) = WORK(1)*MU
         WORK(M) = WORK(M)*NU
C
C        ... SOLVE
C
         IF (ABS(DIAG(K)) .GE. 2.0*ABS(OFFDG(K))) THEN
C
C           ... CASE OF COEFU .LE. 0   --   USE EVANS ALGORITHM
C
            IF (PRDC) THEN
               CALL TRSOLP(DIAG(K), OFFDG(K), WORK, M, INFO)
            ELSE
               CALL TRSOLG (DIAG(K), OFFDG(K), MU, NU, WORK, M, INFO)
            ENDIF
         ELSE
C
C           ... CASE OF COEFU .GT. 0  --   USE GAUSS WITH PIVOTING
C
            DO 200 IPOS=0,M-1
               WORK(LOCC+IPOS) = OFFDG(K)
               WORK(LOCE+IPOS) = OFFDG(K)
               WORK(LOCD+IPOS) = DIAG(K)
  200       CONTINUE
            WORK(LOCD) = MU*WORK(LOCD)
            WORK(LOCD+M-1) = NU*WORK(LOCD+M-1)
            CORNER = 0.0E0
            IF (PRDC)  CORNER = OFFDG(K)
            CALL SGPSL(M, WORK(LOCC), WORK(LOCD), WORK(LOCE),CORNER,
     +                 CORNER, WORK, WORK(LOCW), INFO)
         ENDIF
C
C        ... REPLACE ROW
C
         KX = 0
         DO 300 J=JL,JR
            KX = KX + 1
            U(I,J) = WORK(KX)
  300    CONTINUE
  500 CONTINUE
C
C
C  ----
C  EXIT
C  ----
C
      RETURN
      END
      SUBROUTINE TRSOLG (A, B, MU, NU, G, NG, INFO)
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C    ---------------      4TH ORDER ACCURATE FAST DIRECT SOLUTION
C    PACKAGE :  HFFT      OF THE HELMHOLTZ EQUATION ON RECTANGULAR
C    ---------------      DOMAINS IN TWO AND THREE DIMENSIONS
C
C    INTERNAL MODULE
C
C    RONALD F. BOISVERT
C    JOAN M. BAUMANN
C    NATIONAL BUREAU OF STANDARDS
C    DECEMBER 1985
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C
C   TRSOLG SOLVES THE TRIDIAGONAL SYSTEM AX = G, WHERE
C
C                       :--                      --:
C                       :  MU*A  B                 :
C                       :     B  A  B              :
C                       :        B  A  B           :
C                A  =   :          .  .  .         :
C                       :            .  .  .       :
C                       :              B  A  B     :
C                       :                 B  A*NU  :
C                       :--                      --:
C
C   WHERE ABS(A).GE.2*ABS(B)) AND MU AND NU HAVE ONE OF THE VALUES
C   0.5 OR 1.0.
C
C
C   ALGORITHM
C
C   THE ALGORITHM USED IS N EXTENSION OF AN INTERLOCKING FACTORIZATION
C   METHOD DUE TO EVANS.
C
C   REFERENCE:  D. J. EVANS, AN ALGORITHM FOR THE SOLUTION OF CERTAIN
C   TRIDIAGONAL SYSTEMS OF LINEAR EQUATIONS, THE COMPUTER JOURNAL,
C   VOL. 15, PP. 356-359.
C
C
C   DEGENERATE CASES
C
C   WHEN MU=NU=0.5 AND ABS(A).EQ.2*ABS(B) THE MATRIX A IS SINGULAR.
C   IN THESE CASES A CERTAIN CONSISTENCY CONDITION MUST BE SATISFIED;
C   IF IT IS, THERE ARE AN INFINITE NUMBER OF SOLUTIONS, EACH DIFFERING
C   BY AN ADDITIVE CONSTANT.  THE CONSISTENCY CONDITIONS ARE
C
C        CASE -A=2B :  SUM(I=1,..,NG) G(I) = 0
C        CASE  A=2B :  SUM(I ODD) G(I) + SUM(I EVEN) G(I) = 0
C
C   WE ASSUME THESE CONDITIONS HOLD AND SELECT THE UNIQUE SOLUTION WITH
C   G(NG) = 0.0.
C
C
C   P A R A M E T E R S
C   -------------------
C
C     A           REAL SCALAR (INPUT)
C                 CONSTANT DIAGONAL FOR THE SYSTEM TO BE SOLVED.
C
C     B           REAL SCALAR (INPUT)  ABS(A) .GE. 2*ABS(B)
C                 CONSTANT OFF-DIAGONAL FOR THE SYSTEM TO BE SOLVED.
C
C     MU,NU       REAL  SCALARS (INPUT)  .EQ. 0.5 OR 1.0
C                 MULTIPLICATIVE FACTORS FOR FIRST AND LAST DIAGONAL
C                 ELEMENTS, RESPECTIVELY.
C
C     G           REAL ARRAY OF SIZE M (INPUT/OUTPUT)
C                 ON INPUT, CONTAINS THE RIGHT HAND SIDE OF THE LINEAR
C                 SYSTEM.  ON OUTPUT THIS IS REPLACED BY THE SOLUTIONONS
C                 VECTOR.
C
C     NG          INTEGER SCALAR (INPUT)  .GE. 2
C                 THE NUMBER OF ROWS IN THE SYSTEMS TO BE SOLVED.
C
C     INFO        INTEGER SCALAR (OUTPUT)
C                 INDICATES STATUS OF COMPUTED SOLUTION.
C                 POSSIBLE VALUES ARE
C
C                    0 == SUBROUTINE RAN TO COMPLETION (SUCCESS)
C                    1 == ERROR. ABS(A).LT.2*ABS(B)
C                    2 == ERROR. N.LT.2
C                    3 == ERROR. MU OR NU .NE. 0.5 OR 1.0
C
C
C  ------------
C  DECLARATIONS
C  ------------
C
C     ... PARAMETERS
C
      REAL A, B, MU, NU, G(NG)
      INTEGER NG, INFO
C
C     ... LOCAL VARIABLES
C
      REAL AA
      INTEGER KASE, N
C
C     ... LOCAL CONSTANTS
C
      REAL ALPHA, ALPHA2, BETA, BETA2, GAMMA2, DELTA,
     *     EPS, XN, CON1, CON2, ZERO, ONE, TWO, FOUR, HALF, EPMACH
C
      SAVE  ZERO,  ONE,   TWO,   FOUR,  HALF
      DATA  ZERO,  ONE,   TWO,   FOUR,  HALF
     *    / 0.0E0, 1.0E0, 2.0E0, 4.0E0, 0.50E0 /
C
C
C  ---------------
C  INITIALIZATIONS
C  ---------------
C
      N = NG
      EPMACH = R1MACH(4)
      EPS = 50.0E0*EPMACH
      IF (ABS(MU - HALF) .LT. EPS)  MU = HALF
      IF (ABS(MU - ONE ) .LT. EPS)  MU = ONE
      IF (ABS(NU - HALF) .LT. EPS)  NU = HALF
      IF (ABS(NU - ONE ) .LT. EPS)  NU = ONE
C
C     ...  CHECK FOR ILLEGAL INPUT PARAMETERS
C
      INFO = 0
      IF (((MU .NE. HALF) .AND. (MU .NE. ONE)) .OR.
     *    ((NU .NE. HALF) .AND. (NU .NE. ONE))      )  INFO = 3
      IF (N .LT. 2)  INFO = 2
      IF ((ABS(A) - TWO*ABS(B)) .LT. ZERO)  INFO = 1
      IF (INFO .NE. 0)  GO TO 999
C
C
C  -----------------------
C  CHECK FOR B=0, CASE = 6
C  -----------------------
C
      IF (ABS(B/A) .LT. EPMACH)  GO TO 600
C
C
C  ---------------------------
C  PREPROCESSING FOR CASES 1-5
C  ---------------------------
C
      ALPHA = -TWO*B/(A + SQRT(MAX(ZERO,A*A - FOUR*B*B)))
      IF (ABS(ALPHA) .GT. ONE)  ALPHA = ONE/ALPHA
      ALPHA2 = ALPHA*ALPHA
      BETA2  = MU*(ONE + ALPHA2) - ALPHA2
      BETA   = SQRT(BETA2)
      GAMMA2 = NU*(ONE + ALPHA2) - ONE
C
C     ... RESCALE SYSTEM TO TRIDIAGONAL(-ALPHA,1+ALPHA**2,-ALPHA)
C
      CON1   = (ONE + ALPHA2)/A
      DO 10 I=1,N
         G(I) = CON1*G(I)
 10   CONTINUE
C
C
C  --------------
C  DETERMINE CASE
C  --------------
C
C     KASE = 1  ==  ABS(ALPHA) < 1.0
C     KASE = 2  ==  ABS(ALPHA) = 1.0,  MU = 1.0,  NU = 1.0
C     KASE = 3  ==  ABS(ALPHA) = 1.0,  MU =  .5,  NU =  .5
C     KASE = 4  ==  ABS(ALPHA) = 1.0,  MU =  .5,  NU = 1.0
C     KASE = 5  ==  ABS(ALPHA) = 1.0,  MU = 1.0,  NU =  .5
C
      IF ((-ONE + EPS .LE. ALPHA) .AND. (ALPHA .LE. ONE - EPS)) THEN
            KASE = 1
         ELSE IF ((MU .EQ. ONE) .AND. (NU .EQ. ONE)) THEN
            KASE = 2
         ELSE IF ((MU .EQ. HALF) .AND. (NU .EQ. HALF)) THEN
            KASE = 3
         ELSE IF ((MU .EQ. HALF) .AND. (NU .EQ. ONE)) THEN
            KASE = 4
         ELSE IF ((MU .EQ. ONE) .AND. (NU .EQ. HALF)) THEN
            KASE = 5
         ELSE
            INFO = 3
            GO TO 999
      ENDIF
      GO TO (100, 200, 300, 400, 500), KASE
C
C
C  ------
C  CASE 1
C  ------
C
C     ... PREPROCESSING
C
 100  CONTINUE
      CON2 = (ONE - ALPHA2 - BETA2)/BETA2
      AA = ALPHA
      XN = G(N)
      DO 110 J=1,N-2
         XN = XN + AA*G(N-J)
         AA = AA*ALPHA
 110  CONTINUE
      XN = XN + AA*(ONE - ALPHA2)/BETA2*G(1)
      DO 120 J=2,N
         AA = AA*ALPHA
         XN = XN + CON2*AA*G(J)
 120  CONTINUE
      DELTA = ONE - ALPHA2 + GAMMA2 + GAMMA2*CON2*AA
      XN = XN/DELTA
      GOTO 220
C
C
C  ------
C  CASE 2
C  ------
C
C     ... PREPROCESSING
C
  200 CONTINUE
      XN = ZERO
      IF (MOD(N-1,2) .EQ. 1) THEN
         AA = ALPHA
      ELSE
         AA = ONE
      ENDIF
      DO 210 K=1,N-1
         XN = XN + G(K)*REAL(K)*AA
         AA = AA*ALPHA
  210 CONTINUE
      XN = (XN + G(N)*REAL(N))/REAL(N+1)
C
C     ... BACK SUBSTITUTION (CASES 1 AND 2)
C
 220  CONTINUE
      G(N) = G(N) - GAMMA2*XN
      DO 230 I=N-1,2,-1
         G(I) = G(I) + ALPHA*G(I+1)
 230  CONTINUE
      G(1) = (G(1) + ALPHA*G(2))/BETA
C
C     ... FORWARD ELIMINATION (CASES 1 AND 2)
C
      G(1) = G(1)/BETA
      DO 240 I=2,N-1
         G(I) = G(I) + ALPHA*G(I-1)
 240   CONTINUE
      G(N)=XN
      GO TO 999
C
C
C  --------------------------
C  CASE 3 : CONVERT TO CASE 4
C  --------------------------
C
 300  CONTINUE
      G(N) = ZERO
      N = N-1
C
C  ------
C  CASE 4
C  ------
C
C     ... FORWARD ELIMINATION
C
 400  CONTINUE
      G(1) = -ALPHA*G(1)
      DO 410 I=2,N
         G(I) = -ALPHA*(G(I) - G(I-1))
 410  CONTINUE
C
C     ... BACK SUBSTITUTION
C
      G(N) = -ALPHA*G(N)
      DO 420 I=N-1,1,-1
         G(I) = -ALPHA*(G(I) - G(I+1))
 420  CONTINUE
      GO TO 999
C
C
C  ------
C  CASE 5
C  ------
C
C     ... BACK SUBSTITUTION
C
 500  CONTINUE
      DO 510 I=N-1,1,-1
         G(I) = G(I) + ALPHA*G(I+1)
 510  CONTINUE
C
C     ... FORWARD ELMINATION
C
      DO 520 I=2,N
         G(I) = G(I) + ALPHA*G(I-1)
 520  CONTINUE
      GO TO 999
C
C
C  ------------------------
C  CASE 6 : DIAGONAL SYSTEM
C  ------------------------
C
 600  CONTINUE
      DO 610 I=1,N
         G(I) = G(I)/A
 610  CONTINUE
      G(1) = G(1)/MU
      G(N) = G(N)/NU
      GO TO 999
C
C
C  ----
C  EXIT
C  ----
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE TRSOLP (A, B, G, NG, INFO)
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C    ---------------      4TH ORDER ACCURATE FAST DIRECT SOLUTION
C    PACKAGE :  HFFT      OF THE HELMHOLTZ EQUATION ON RECTANGULAR
C    ---------------      DOMAINS IN TWO AND THREE DIMENSIONS
C
C    INTERNAL MODULE
C
C    RONALD F. BOISVERT
C    JOAN M. BAUMANN
C    NATIONAL BUREAU OF STANDARDS
C
C    DECEMBER 1985   (REVISED APRIL 1987)
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C
C   TRSOLP SOLVES THE TRIDIAGONAL LINEAR SYSTEM AX = G, WHERE
C
C                        :--                  --:
C                        :   A  B           B   :
C                        :   B  A  B            :
C                        :      B  A  B         :
C                 A  =   :        .  .  .       :
C                        :          .  .  .     :
C                        :            B  A  B   :
C                        :   B           B  A   :
C                        :--                  --:
C
C   WHERE ABS(A).GE.ABS(2*B).
C
C
C   ALGORITHM
C
C    THE ALGORITHM USED IS AN INTERLOCKING FACTORIZATION METHOD DUE TO
C    EVANS.
C
C   REFERENCE: D. J. EVANS, FAST ADI METHODS FOR THE SOLUTION OF LINEAR
C   PARABOLIC PARTIAL DIFFERENTIAL EQUATIONS INVOLVING 2 SPACE
C   DIMENSIONS, BIT, VOL. 17, P.486-491.
C
C
C   DEGENERATE CASES
C
C   WHEN -A=2B OR WHEN A=2B AND NG IS EVEN THE MATRIX A IS SINGULAR.
C   IN THESE CASES THERE ARE NO SOLUTIONS UNLESS G SATISFIES A CERTAIN
C   CONSISTENCY CONDITION.  WHEN IT DOES, THERE ARE AN INFINITE NUMBER
C   OF SOLUTIONS, EACH DIFFERENING BY AN ADDITVE CONSTANT.  THE
C   CONSISTENCY CONDITIONS ARE
C
C     CASE -A=2B          :  SUM(I=1,..,NG) G(I) = 0
C     CASE  A=2B, NG EVEN :  SUM(I ODD) G(I) - SUM(I EVEN) G(I) = 0
C
C   WE ASSUME THESE CONDITIONS HOLD AND SELECT THE UNIQUE SOLUTION WITH
C   G(NG) = 0.0.
C
C
C   -------------------
C   P A R A M E T E R S
C   -------------------
C
C     A           REAL SCALAR (INPUT)
C                 CONSTANT DIAGONAL FOR THE SYSTEM TO BE SOLVED
C
C     B           REAL SCALAR (INPUT)
C                 CONSTANT OFF-DIAGONAL FOR THE SYSTEM TO BE SOLVED
C
C     G           REAL ARRAY OF SIZE M (INPUT/OUTPUT)
C                 ON INPUT THE RIGHT HAND SIDES OF THE LINEAR
C                 SYSTEMS TO BE SOLVED ARE STORED IN THE ROWS
C                 OF G, I.E., THE ITH RIGHT HAND SIDE IS IN G(I).
C                 ON OUTPUT THESE ARE REPLACED BY THE SOLUTIONS
C                 OF THE LINEAR SYSTEMS.
C
C     NG          INTEGER SCALAR (INPUT)
C                 THE NUMBER OF ROWS IN THE SYSTEMS TO BE SOLVED.
C
C     INFO        INTEGER SCALAR (OUTPUT)
C                 INDICATES STATUS OF COMPUTED SOLUTION.
C                 POSSIBLE VALUES ARE
C
C                    0 == SUBROUTINE RAN TO COMPLETION (SUCCESS)
C                    1 == ERROR. ABS(A).LT.2*ABS(B)
C                    2 == ERROR. N.LT.3
C
C
C  ------------
C  DECLARATIONS
C  ------------
C
C     ... PARAMETERS
C
      REAL A, B, G(NG)
      INTEGER NG, INFO
C
C     ... LOCAL VARIABLES
C
      REAL AA
      INTEGER N
C
C     ... LOCAL CONSTANTS
C
      REAL ALPHA, EPS, FACTOR, ZERO, ONE, TWO, FOUR, EPMACH
C
      SAVE  ZERO,  ONE,   TWO,   FOUR
      DATA  ZERO,  ONE,   TWO,   FOUR
     *    / 0.0E0, 1.0E0, 2.0E0, 4.0E0 /
C
C
C  ---------------
C  INITIALIZATIONS
C  ---------------
C
      N = NG
      EPMACH = R1MACH(4)
      EPS = 50.0E0*EPMACH
C
C     ...  CHECK FOR ILLEGAL INPUT PARAMETERS
C
      INFO = 0
C      IF ((ABS(A) - 2.0E0*ABS(B)) .LT. ZERO)  INFO = 1
      IF (N .LT. 3)  INFO = 2
      IF (INFO .NE. 0)  GO TO 999
C
C
C  ------------------
C  CHECK FOR B=0 CASE
C  ------------------
C
      IF (ABS(B/A) .LT. EPMACH)  GO TO 200
C
C
C  --------------------------------------------------------
C  RESCALE PROBLEM TO TRIDIAGONAL(-ALPHA,1+ALPHA**2,-ALPHA)
C  --------------------------------------------------------
C
      ALPHA = -TWO*B/(A + SQRT(MAX(ZERO,A*A - FOUR*B*B)))
      IF (ABS(ALPHA) .GT. ONE) ALPHA = ONE/ALPHA
      FACTOR = (ONE + ALPHA*ALPHA)/A
      DO 10 I = 1,N
         G(I) = FACTOR*G(I)
   10 CONTINUE
C
C
C  -----------------------
C  CHECK FOR SINGULAR CASE
C  -----------------------
C
      IF (((ABS(ONE + ALPHA) .LE. EPS) .AND. (MOD(N,2) .EQ. 0)) .OR.
     *    (ABS(ONE - ALPHA) .LE. EPS))  GO TO 300
C
C
C  -------------
C  STANDARD CASE
C  -------------
C
C     ... PREPROCESSING
C
      AA = ALPHA
      DO 20 I = N,2,-1
         G(1) = G(1) + AA*G(I)
         AA = AA*ALPHA
   20 CONTINUE
      G(1) = G(1)/(ONE - AA)
C
C    ... FORWARD ELIMINATION
C
      DO 40 I = 2,N
         G(I) = G(I) + ALPHA*G(I-1)
   40 CONTINUE
C
C     ... DETERMINE INTERLOCKING ELEMENT XN
C
      AA = ALPHA
      DO 45 I = 1,N-1
         G(N) = G(N) + AA*G(I)
         AA = AA*ALPHA
   45 CONTINUE
      G(N) = G(N)/(ONE - AA)
C
C     ... BACK SUBSTITUTION
C
      DO 50 I = N-1,1,-1
         G(I) = G(I) + ALPHA*G(I+1)
   50 CONTINUE
      GO TO 999
C
C
C  ---------------------
C  SPECIAL CASE :  B = 0
C  ---------------------
C
  200 CONTINUE
      DO 210 I=1,N
         G(I) = G(I)/A
  210 CONTINUE
      GO TO 999
C
C
C  --------------------------------
C  SPECIAL CASE :  SINGULAR PROBLEM
C  --------------------------------
C
  300 CONTINUE
      G(N) = ZERO
      CALL TRSOLG (ONE + ALPHA*ALPHA, -ALPHA, ONE, ONE, G, N-1, INFO)
      GO TO 999
C
C
C  ----
C  EXIT
C  ----
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE SGPSL (N,C,D,E,C0,E0,B,W,INFO)
C
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C    ---------------      4TH ORDER ACCURATE FAST DIRECT SOLUTION
C    PACKAGE :  HFFT      OF THE HELMHOLTZ EQUATION ON RECTANGULAR
C    ---------------      DOMAINS IN TWO AND THREE DIMENSIONS
C
C    INTERNAL MODULE
C
C    RONALD F. BOISVERT
C    NATIONAL BUREAU OF STANDARDS
C
C    APRIL 1987
C
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*
C
C  SGPSL SOLVES A LINEAR SYSTEM OF EQUATIONS DEFINED BY A GENERAL
C  TRIDIAGONAL MATRIX WITH ADDITIONAL NONZEROS IN THE (1,N) AND (N,1)
C  POSITIONS USING GAUSS ELIMINATION WITH PARTIAL PIVOTING.
C
C
C  P A R A M E T E R S
C  -------------------
C
C  ON ENTRY
C
C     N       INTEGER (N.GE.3)
C             IS THE ORDER OF THE TRIDIAGONAL MATRIX.
C
C     C       REAL(N)
C             IS THE SUBDIAGONAL OF THE TRIDIAGONAL MATRIX.
C             C(2) THROUGH C(N) SHOULD CONTAIN THE SUBDIAGONAL.
C             ON OUTPUT, C IS DESTROYED.
C
C     D       REAL(N)
C             IS THE DIAGONAL OF THE TRIDIAGONAL MATRIX.
C             ON OUTPUT, D IS DESTROYED.
C
C     E       REAL(N)
C             IS THE SUPERDIAGONAL OF THE TRIDIAGONAL MATRIX.
C             E(1) THROUGH E(N-1) SHOULD CONTAIN THE SUPERDIAGONAL.
C             ON OUTPUT, E IS DESTROYED.
C
C     C0      REAL
C             THE NONZERO ELEMENT A(1,N).
C
C     E0      REAL
C             THE NONZERO ELEMENT A(N,1).
C
C     B       REAL(N)
C             IS THE RIGHT HAND SIDE VECTOR.
C
C     W       REAL(N)
C             WORKSPACE.
C
C  ON RETURN
C
C     B       IS THE SOLUTION VECTOR.
C
C     INFO    INTEGER
C             = -1 IF N .LT. 3
C             = 0  NORMAL VALUE.
C             = K  IF THE K-TH ELEMENT OF THE DIAGONAL BECOMES
C                  EXACTLY ZERO.  THE SUBROUTINE RETURNS WHEN
C                  THIS IS DETECTED.
C
C  THIS IS A MODIFICATION OF THE LINPACK ROUTINE SGTSL.
C
C
C  ------------
C  DECLARATIONS
C  ------------
C
C     ... PARAMETERS
C
      INTEGER N, INFO
      REAL C(N), D(N), E(N), B(N), W(N)
C
C     ... LOCAL VARIABLES
C
      LOGICAL LASROW
      INTEGER K, KBIG, KP1, KP2, NM1, NM2, NM3
      REAL T, RK, RKP1, RKP2, A(3,3)
C
C
C  ---------------
C  INITIALIZATIONS
C  ---------------
C
      INFO = -1
      IF (N .LE. 3)  GO TO 100
C
      NM1 = N - 1
      NM2 = N - 2
      NM3 = N - 3
C
      ZERO = 0.0E0
      LASROW = E0 .NE. ZERO
      RKP1 = E0
      RKP2 = ZERO
      W(1) = C0
      DO 10 I=2,NM2
         W(I) = ZERO
   10 CONTINUE
      W(NM1) = E(NM1)
      W(N) = D(N)
C
      C(1) = D(1)
      D(1) = E(1)
      E(1) = ZERO
      E(N) = ZERO
C
C  -------------------
C  FORWARD ELIMINATION
C  -------------------
C
      DO 30 K = 1, NM3
         INFO = K
         KP1 = K + 1
         KP2 = K + 2
         RK = RKP1
         RKP1 = RKP2
         RKP2 = ZERO
         IF (K .EQ. NM3)  RKP2 = C(N)
C
C        ... FIND THE LARGEST OF THE TWO ROWS
C
         IF (ABS(C(KP1)) .GT. ABS(C(K))) THEN
C
C           ... INTERCHANGE ROW
C
            T = C(KP1)
            C(KP1) = C(K)
            C(K) = T
            T = D(KP1)
            D(KP1) = D(K)
            D(K) = T
            T = E(KP1)
            E(KP1) = E(K)
            E(K) = T
            T = W(KP1)
            W(KP1) = W(K)
            W(K) = T
            T = B(KP1)
            B(KP1) = B(K)
            B(K) = T
         ENDIF
C
C        ... CHECK FOR SINGULARITY
C
         IF (C(K) .EQ. ZERO) GO TO 100
C
C        ... ELIMINATE IN ROW K+1
C
         T = -C(KP1)/C(K)
         C(KP1) = D(KP1) + T*D(K)
         D(KP1) = E(KP1) + T*E(K)
         E(KP1) = ZERO
         W(KP1) = W(KP1) + T*W(K)
         B(KP1) = B(KP1) + T*B(K)
C
C        ... ELIMINATE IN LAST ROW
C
         IF (LASROW) THEN
            T = -RK/C(K)
            RKP1 = RKP1 + T*D(K)
            RKP2 = RKP2 + T*E(K)
            W(N) = W(N) + T*W(K)
            B(N) = B(N) + T*B(K)
         ENDIF
   30 CONTINUE
C
C     ... DO LAST 3 BY 3 BLOCK
C
      A(1,1) = C(NM2)
      A(1,2) = D(NM2)
      A(1,3) = W(NM2)
      A(2,1) = C(NM1)
      A(2,2) = D(NM1)
      A(2,3) = W(NM1)
      A(3,1) = RKP1
      A(3,2) = RKP2
      A(3,3) = W(N)
C
C     === STEP N-2 ===
C
      INFO = NM2
      KBIG = 1
      IF (ABS(A(2,1)) .GT. ABS(A(1,1)))    KBIG = 2
      IF (ABS(A(3,1)) .GT. ABS(A(KBIG,1))) KBIG = 3
      IF (KBIG .NE. 1) THEN
C
C        ... PIVOT
C
         T = A(KBIG,1)
         A(KBIG,1) = A(1,1)
         A(1,1) = T
         T = A(KBIG,2)
         A(KBIG,2) = A(1,2)
         A(1,2) = T
         T = A(KBIG,3)
         A(KBIG,3) = A(1,3)
         A(1,3) = T
         K = NM3 + KBIG
         T = B(K)
         B(K) = B(NM2)
         B(NM2) = T
      ENDIF
      IF (A(1,1) .EQ. ZERO)  GO TO 100
C
C     ... ELIMINATE
C
      T = -A(2,1)/A(1,1)
      A(2,2) = A(2,2) + T*A(1,2)
      A(2,3) = A(2,3) + T*A(1,3)
      B(NM1) = B(NM1) + T*B(NM2)
      T = -A(3,1)/A(1,1)
      A(3,2) = A(3,2) + T*A(1,2)
      A(3,3) = A(3,3) + T*A(1,3)
      B(N) = B(N) + T*B(NM2)
C
C     === STEP N-1 ===
C
      INFO = NM1
      IF (ABS(A(3,2)) .GT. ABS(A(2,2))) THEN
C
C        ... PIVOT
C
         T = A(3,2)
         A(3,2) = A(2,2)
         A(2,2) = T
         T = A(3,3)
         A(3,3) = A(2,3)
         A(2,3) = T
         T = B(N)
         B(N) = B(NM1)
         B(NM1) = T
      ENDIF
      IF (A(2,2) .EQ. ZERO)  GO TO 100
C
C     ... ELIMINATE
C
      T = -A(3,2)/A(2,2)
      A(3,3) = A(3,3) + T*A(2,3)
      B(N) = B(N) + T*B(NM1)
C
C     === STEP N ===
C
      INFO = N
      IF (A(3,3) .EQ. ZERO)  GO TO 100
C
C  ----------
C  BACK SOLVE
C  ----------
C
      B(N) = B(N)/A(3,3)
      B(NM1) = (B(NM1) - A(2,3)*B(N))/A(2,2)
      B(NM2) = (B(NM2) - A(1,2)*B(NM1) - A(1,3)*B(N))/A(1,1)
      DO 60 K = NM3, 1, -1
         B(K) = (B(K) - D(K)*B(K+1) - E(K)*B(K+2) - W(K)*B(N))/C(K)
   60 CONTINUE
      INFO = 0
C
C  ----
C  EXIT
C  ----
C
  100 CONTINUE
      RETURN
      END
*
      SUBROUTINE COSQB(N,X,WSAVE)
C***BEGIN PROLOGUE  COSQB
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  860115   (YYMMDD)
C***CATEGORY NO.  J1A3
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  UNNORMALIZED INVERSE OF COSQF.
C***DESCRIPTION
C
C  SUBROUTINE COSQB COMPUTES THE FAST FOURIER TRANSFORM OF QUARTER
C  WAVE DATA. THAT IS, COSQB COMPUTES A SEQUENCE FROM ITS
C  REPRESENTATION IN TERMS OF A COSINE SERIES WITH ODD WAVE NUMBERS.
C  THE TRANSFORM IS DEFINED BELOW AT OUTPUT PARAMETER X.
C
C  COSQB IS THE UNNORMALIZED INVERSE OF COSQF SINCE A CALL OF COSQB
C  FOLLOWED BY A CALL OF COSQF WILL MULTIPLY THE INPUT SEQUENCE X
C  BY 4*N.
C
C  THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE COSQB MUST BE
C  INITIALIZED BY CALLING SUBROUTINE COSQI(N,WSAVE).
C
C
C  INPUT PARAMETERS
C
C  N       THE LENGTH OF THE ARRAY X TO BE TRANSFORMED.  THE METHOD
C          IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
C
C  X       AN ARRAY WHICH CONTAINS THE SEQUENCE TO BE TRANSFORMED
C
C  WSAVE   A WORK ARRAY THAT MUST BE DIMENSIONED AT LEAST 3*N+15
C          IN THE PROGRAM THAT CALLS COSQB.  THE WSAVE ARRAY MUST BE
C          INITIALIZED BY CALLING SUBROUTINE COSQI(N,WSAVE), AND A
C          DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
C          VALUE OF N.  THIS INITIALIZATION DOES NOT HAVE TO BE
C          REPEATED SO LONG AS N REMAINS UNCHANGED.  THUS SUBSEQUENT
C          TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
C
C  OUTPUT PARAMETERS
C
C  X       FOR I=1,...,N
C
C               X(I)= THE SUM FROM K=1 TO K=N OF
C
C                 4*X(K)*COS((2*K-1)*(I-1)*PI/(2*N))
C
C               A CALL OF COSQB FOLLOWED BY A CALL OF
C               COSQF WILL MULTIPLY THE SEQUENCE X BY 4*N.
C               THEREFORE COSQF IS THE UNNORMALIZED INVERSE
C               OF COSQB.
C
C  WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT
C          BE DESTROYED BETWEEN CALLS OF COSQB OR COSQF.
C
C  *********************************************************************
C  *                                                                   *
C  *   SUBPROGRAM REVISION HISTORY                                     *
C  *                                                                   *
C  *   06/01/79  -  ORIGINAL VERSION BY PAUL SWARZTRAUBER.             *
C  *                DISTRIBUTED BY NCAR (REF. 1).                      *
C  *   04/01/83  -  SLATEC COMMON MATH LIBRARY SUBCOMMITTEE.           *
C  *                MODIFIED TO USE SLATEC LIBRARY SOURCE FILE FORMAT. *
C  *                DISTRIBUTED IN THE SLATEC LIBRARY (REF. 2).        *
C  *   01/15/86  -  RON BOISVERT, NATIONAL BUREAU OF STANDARDS.        *
C  *                MODIFIED TO CONVERT TO PORTABLE FORTRAN 77.        *
C  *                                                                   *
C  *   THE CHANGES INTRODUCED IN THE MOST RECENT MODIFICATION ARE      *
C  *                                                                   *
C  *   (A) DUMMY ARRAY SIZE DECLARATIONS (1) CHANGED TO (*)            *
C  *   (B) REFERENCES TO INTRINSIC FUNCTION FLOAT CHANGED TO REAL      *
C  *   (C) MATHEMATICAL CONSTANTS PREVIOUSLY CODED IN DATA STATE-      *
C  *       MENTS NOW COMPUTED AT RUNTIME USING FORTRAN INTRINSIC       *
C  *       FUNCTIONS.  THE AFFECTED VARIABLES ARE                      *
C  *                                                                   *
C  *          PI      SQRT2   SQRT3   TAUR    TR11    TR12             *
C  *          PIH     TSQRT2          TAUI    TI11    TI12             *
C  *          TPI     HSQT2                                            *
C  *                                                                   *
C  *   REFERENCES                                                      *
C  *                                                                   *
C  *   1. P.N. SWARZTRAUBER, VECTORIZING THE FFTS, IN PARALLEL         *
C  *      COMPUTATIONS (G. RODRIGUE, ED.), ACADEMIC PRESS, 1982,       *
C  *      PP. 51-83.                                                   *
C  *   2. B.L. BUZBEE, THE SLATEC COMMON MATH LIBRARY, IN SOURCES      *
C  *      AND DEVELOPMENT OF MATHEMATICAL SOFTWARE (W. COWELL, ED.),   *
C  *      PRENTICE-HALL, 1984, PP. 302-318.                            *
C  *                                                                   *
C  *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  COSQB1
C***END PROLOGUE  COSQB
      DIMENSION       X(*)       ,WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  COSQB
      TSQRT2 = 2.*SQRT(2.)
      IF (N-2) 101,102,103
  101 X(1) = 4.*X(1)
      RETURN
  102 X1 = 4.*(X(1)+X(2))
      X(2) = TSQRT2*(X(1)-X(2))
      X(1) = X1
      RETURN
  103 CALL COSQB1 (N,X,WSAVE,WSAVE(N+1))
      RETURN
      END
      SUBROUTINE COSQF(N,X,WSAVE)
C***BEGIN PROLOGUE  COSQF
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  860115   (YYMMDD)
C***CATEGORY NO.  J1A3
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  FORWARD COSINE TRANSFORM WITH ODD WAVE NUMBERS.
C***DESCRIPTION
C
C  SUBROUTINE COSQF COMPUTES THE FAST FOURIER TRANSFORM OF QUARTER
C  WAVE DATA. THAT IS, COSQF COMPUTES THE COEFFICIENTS IN A COSINE
C  SERIES REPRESENTATION WITH ONLY ODD WAVE NUMBERS.  THE TRANSFORM
C  IS DEFINED BELOW AT OUTPUT PARAMETER X
C
C  COSQF IS THE UNNORMALIZED INVERSE OF COSQB SINCE A CALL OF COSQF
C  FOLLOWED BY A CALL OF COSQB WILL MULTIPLY THE INPUT SEQUENCE X
C  BY 4*N.
C
C  THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE COSQF MUST BE
C  INITIALIZED BY CALLING SUBROUTINE COSQI(N,WSAVE).
C
C
C  INPUT PARAMETERS
C
C  N       THE LENGTH OF THE ARRAY X TO BE TRANSFORMED.  THE METHOD
C          IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
C
C  X       AN ARRAY WHICH CONTAINS THE SEQUENCE TO BE TRANSFORMED
C
C  WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15
C          IN THE PROGRAM THAT CALLS COSQF.  THE WSAVE ARRAY MUST BE
C          INITIALIZED BY CALLING SUBROUTINE COSQI(N,WSAVE), AND A
C          DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
C          VALUE OF N.  THIS INITIALIZATION DOES NOT HAVE TO BE
C          REPEATED SO LONG AS N REMAINS UNCHANGED.  THUS SUBSEQUENT
C          TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
C
C  OUTPUT PARAMETERS
C
C  X       FOR I=1,...,N
C
C               X(I) = X(1) PLUS THE SUM FROM K=2 TO K=N OF
C
C                  2*X(K)*COS((2*I-1)*(K-1)*PI/(2*N))
C
C               A CALL OF COSQF FOLLOWED BY A CALL OF
C               COSQB WILL MULTIPLY THE SEQUENCE X BY 4*N.
C               THEREFORE COSQB IS THE UNNORMALIZED INVERSE
C               OF COSQF.
C
C  WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT
C          BE DESTROYED BETWEEN CALLS OF COSQF OR COSQB.
C
C  *********************************************************************
C  *                                                                   *
C  *   SUBPROGRAM REVISION HISTORY                                     *
C  *                                                                   *
C  *   06/01/79  -  ORIGINAL VERSION BY PAUL SWARZTRAUBER.             *
C  *                DISTRIBUTED BY NCAR (REF. 1).                      *
C  *   04/01/83  -  SLATEC COMMON MATH LIBRARY SUBCOMMITTEE.           *
C  *                MODIFIED TO USE SLATEC LIBRARY SOURCE FILE FORMAT. *
C  *                DISTRIBUTED IN THE SLATEC LIBRARY (REF. 2).        *
C  *   01/15/86  -  RON BOISVERT, NATIONAL BUREAU OF STANDARDS.        *
C  *                MODIFIED TO CONVERT TO PORTABLE FORTRAN 77.        *
C  *                                                                   *
C  *   THE CHANGES INTRODUCED IN THE MOST RECENT MODIFICATION ARE      *
C  *                                                                   *
C  *   (A) DUMMY ARRAY SIZE DECLARATIONS (1) CHANGED TO (*)            *
C  *   (B) REFERENCES TO INTRINSIC FUNCTION FLOAT CHANGED TO REAL      *
C  *   (C) MATHEMATICAL CONSTANTS PREVIOUSLY CODED IN DATA STATE-      *
C  *       MENTS NOW COMPUTED AT RUNTIME USING FORTRAN INTRINSIC       *
C  *       FUNCTIONS.  THE AFFECTED VARIABLES ARE                      *
C  *                                                                   *
C  *          PI      SQRT2   SQRT3   TAUR    TR11    TR12             *
C  *          PIH     TSQRT2          TAUI    TI11    TI12             *
C  *          TPI     HSQT2                                            *
C  *                                                                   *
C  *   REFERENCES                                                      *
C  *                                                                   *
C  *   1. P.N. SWARZTRAUBER, VECTORIZING THE FFTS, IN PARALLEL         *
C  *      COMPUTATIONS (G. RODRIGUE, ED.), ACADEMIC PRESS, 1982,       *
C  *      PP. 51-83.                                                   *
C  *   2. B.L. BUZBEE, THE SLATEC COMMON MATH LIBRARY, IN SOURCES      *
C  *      AND DEVELOPMENT OF MATHEMATICAL SOFTWARE (W. COWELL, ED.),   *
C  *      PRENTICE-HALL, 1984, PP. 302-318.                            *
C  *                                                                   *
C  *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  COSQF1
C***END PROLOGUE  COSQF
      DIMENSION       X(*)       ,WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  COSQF
      SQRT2 = SQRT(2.)
      IF (N-2) 102,101,103
  101 TSQX = SQRT2*X(2)
      X(2) = X(1)-TSQX
      X(1) = X(1)+TSQX
  102 RETURN
  103 CALL COSQF1 (N,X,WSAVE,WSAVE(N+1))
      RETURN
      END
      SUBROUTINE COSQI(N,WSAVE)
C***BEGIN PROLOGUE  COSQI
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  860115   (YYMMDD)
C***CATEGORY NO.  J1A3
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  INITIALIZE FOR COSQF AND COSQB.
C***DESCRIPTION
C
C  SUBROUTINE COSQI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
C  BOTH COSQF AND COSQB.  THE PRIME FACTORIZATION OF N TOGETHER WITH
C  A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
C  STORED IN WSAVE.
C
C  INPUT PARAMETER
C
C  N       THE LENGTH OF THE ARRAY TO BE TRANSFORMED.  THE METHOD
C          IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
C
C  OUTPUT PARAMETER
C
C  WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
C          THE SAME WORK ARRAY CAN BE USED FOR BOTH COSQF AND COSQB
C          AS LONG AS N REMAINS UNCHANGED.  DIFFERENT WSAVE ARRAYS
C          ARE REQUIRED FOR DIFFERENT VALUES OF N.  THE CONTENTS OF
C          WSAVE MUST NOT BE CHANGED BETWEEN CALLS OF COSQF OR COSQB.
C
C  *********************************************************************
C  *                                                                   *
C  *   SUBPROGRAM REVISION HISTORY                                     *
C  *                                                                   *
C  *   06/01/79  -  ORIGINAL VERSION BY PAUL SWARZTRAUBER.             *
C  *                DISTRIBUTED BY NCAR (REF. 1).                      *
C  *   04/01/83  -  SLATEC COMMON MATH LIBRARY SUBCOMMITTEE.           *
C  *                MODIFIED TO USE SLATEC LIBRARY SOURCE FILE FORMAT. *
C  *                DISTRIBUTED IN THE SLATEC LIBRARY (REF. 2).        *
C  *   01/15/86  -  RON BOISVERT, NATIONAL BUREAU OF STANDARDS.        *
C  *                MODIFIED TO CONVERT TO PORTABLE FORTRAN 77.        *
C  *                                                                   *
C  *   THE CHANGES INTRODUCED IN THE MOST RECENT MODIFICATION ARE      *
C  *                                                                   *
C  *   (A) DUMMY ARRAY SIZE DECLARATIONS (1) CHANGED TO (*)            *
C  *   (B) REFERENCES TO INTRINSIC FUNCTION FLOAT CHANGED TO REAL      *
C  *   (C) MATHEMATICAL CONSTANTS PREVIOUSLY CODED IN DATA STATE-      *
C  *       MENTS NOW COMPUTED AT RUNTIME USING FORTRAN INTRINSIC       *
C  *       FUNCTIONS.  THE AFFECTED VARIABLES ARE                      *
C  *                                                                   *
C  *          PI      SQRT2   SQRT3   TAUR    TR11    TR12             *
C  *          PIH     TSQRT2          TAUI    TI11    TI12             *
C  *          TPI     HSQT2                                            *
C  *                                                                   *
C  *   REFERENCES                                                      *
C  *                                                                   *
C  *   1. P.N. SWARZTRAUBER, VECTORIZING THE FFTS, IN PARALLEL         *
C  *      COMPUTATIONS (G. RODRIGUE, ED.), ACADEMIC PRESS, 1982,       *
C  *      PP. 51-83.                                                   *
C  *   2. B.L. BUZBEE, THE SLATEC COMMON MATH LIBRARY, IN SOURCES      *
C  *      AND DEVELOPMENT OF MATHEMATICAL SOFTWARE (W. COWELL, ED.),   *
C  *      PRENTICE-HALL, 1984, PP. 302-318.                            *
C  *                                                                   *
C  *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  RFFTI
C***END PROLOGUE  COSQI
      DIMENSION       WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  COSQI
      PIH = 2.*ATAN(1.)
      DT = PIH/REAL(N)
      FK = 0.
      DO 101 K=1,N
         FK = FK+1.
         WSAVE(K) = COS(FK*DT)
  101 CONTINUE
      CALL RFFTI (N,WSAVE(N+1))
      RETURN
      END
      SUBROUTINE COSTI(N,WSAVE)
C***BEGIN PROLOGUE  COSTI
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  860115   (YYMMDD)
C***CATEGORY NO.  J1A3
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  INITIALIZE FOR COST.
C***DESCRIPTION
C
C  SUBROUTINE COSTI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
C  SUBROUTINE COST.  THE PRIME FACTORIZATION OF N TOGETHER WITH
C  A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
C  STORED IN WSAVE.
C
C  INPUT PARAMETER
C
C  N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED.  THE METHOD
C          IS MOST EFFICIENT WHEN N-1 IS A PRODUCT OF SMALL PRIMES.
C
C  OUTPUT PARAMETER
C
C  WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
C          DIFFERENT WSAVE ARRAYS ARE REQUIRED FOR DIFFERENT VALUES
C          OF N.  THE CONTENTS OF WSAVE MUST NOT BE CHANGED BETWEEN
C          CALLS OF COST.
C
C  *********************************************************************
C  *                                                                   *
C  *   SUBPROGRAM REVISION HISTORY                                     *
C  *                                                                   *
C  *   06/01/79  -  ORIGINAL VERSION BY PAUL SWARZTRAUBER.             *
C  *                DISTRIBUTED BY NCAR (REF. 1).                      *
C  *   04/01/83  -  SLATEC COMMON MATH LIBRARY SUBCOMMITTEE.           *
C  *                MODIFIED TO USE SLATEC LIBRARY SOURCE FILE FORMAT. *
C  *                DISTRIBUTED IN THE SLATEC LIBRARY (REF. 2).        *
C  *   01/15/86  -  RON BOISVERT, NATIONAL BUREAU OF STANDARDS.        *
C  *                MODIFIED TO CONVERT TO PORTABLE FORTRAN 77.        *
C  *                                                                   *
C  *   THE CHANGES INTRODUCED IN THE MOST RECENT MODIFICATION ARE      *
C  *                                                                   *
C  *   (A) DUMMY ARRAY SIZE DECLARATIONS (1) CHANGED TO (*)            *
C  *   (B) REFERENCES TO INTRINSIC FUNCTION FLOAT CHANGED TO REAL      *
C  *   (C) MATHEMATICAL CONSTANTS PREVIOUSLY CODED IN DATA STATE-      *
C  *       MENTS NOW COMPUTED AT RUNTIME USING FORTRAN INTRINSIC       *
C  *       FUNCTIONS.  THE AFFECTED VARIABLES ARE                      *
C  *                                                                   *
C  *          PI      SQRT2   SQRT3   TAUR    TR11    TR12             *
C  *          PIH     TSQRT2          TAUI    TI11    TI12             *
C  *          TPI     HSQT2                                            *
C  *                                                                   *
C  *   REFERENCES                                                      *
C  *                                                                   *
C  *   1. P.N. SWARZTRAUBER, VECTORIZING THE FFTS, IN PARALLEL         *
C  *      COMPUTATIONS (G. RODRIGUE, ED.), ACADEMIC PRESS, 1982,       *
C  *      PP. 51-83.                                                   *
C  *   2. B.L. BUZBEE, THE SLATEC COMMON MATH LIBRARY, IN SOURCES      *
C  *      AND DEVELOPMENT OF MATHEMATICAL SOFTWARE (W. COWELL, ED.),   *
C  *      PRENTICE-HALL, 1984, PP. 302-318.                            *
C  *                                                                   *
C  *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  RFFTI
C***END PROLOGUE  COSTI
      DIMENSION       WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  COSTI
      PI = 4.*ATAN(1.)
      IF (N .LE. 3) RETURN
      NM1 = N-1
      NP1 = N+1
      NS2 = N/2
      DT = PI/REAL(NM1)
      FK = 0.
      DO 101 K=2,NS2
         KC = NP1-K
         FK = FK+1.
         WSAVE(K) = 2.*SIN(FK*DT)
         WSAVE(KC) = 2.*COS(FK*DT)
  101 CONTINUE
      CALL RFFTI (NM1,WSAVE(N+1))
      RETURN
      END
      SUBROUTINE COST(N,X,WSAVE)
C***BEGIN PROLOGUE  COST
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  860115   (YYMMDD)
C***CATEGORY NO.  J1A3
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  COSINE TRANSFORM OF A REAL, EVEN SEQUENCE.
C***DESCRIPTION
C
C  SUBROUTINE COST COMPUTES THE DISCRETE FOURIER COSINE TRANSFORM
C  OF AN EVEN SEQUENCE X(I).  THE TRANSFORM IS DEFINED BELOW AT OUTPUT
C  PARAMETER X.
C
C  COST IS THE UNNORMALIZED INVERSE OF ITSELF SINCE A CALL OF COST
C  FOLLOWED BY ANOTHER CALL OF COST WILL MULTIPLY THE INPUT SEQUENCE
C  X BY 2*(N-1).  THE TRANSFORM IS DEFINED BELOW AT OUTPUT PARAMETER X.
C
C  THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE COST MUST BE
C  INITIALIZED BY CALLING SUBROUTINE COSTI(N,WSAVE).
C
C  INPUT PARAMETERS
C
C  N       THE LENGTH OF THE SEQUENCE X.  N MUST BE GREATER THAN 1.
C          THE METHOD IS MOST EFFICIENT WHEN N-1 IS A PRODUCT OF
C          SMALL PRIMES.
C
C  X       AN ARRAY WHICH CONTAINS THE SEQUENCE TO BE TRANSFORMED
C
C  WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15
C          IN THE PROGRAM THAT CALLS COST.  THE WSAVE ARRAY MUST BE
C          INITIALIZED BY CALLING SUBROUTINE COSTI(N,WSAVE), AND A
C          DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
C          VALUE OF N.  THIS INITIALIZATION DOES NOT HAVE TO BE
C          REPEATED SO LONG AS N REMAINS UNCHANGED.  THUS SUBSEQUENT
C          TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
C
C  OUTPUT PARAMETERS
C
C  X       FOR I=1,...,N
C
C             X(I) = X(1)+(-1)**(I-1)*X(N)
C
C               + THE SUM FROM K=2 TO K=N-1
C
C                   X(K)*COS((K-1)*(I-1)*PI/(N-1))
C
C               A CALL OF COST FOLLOWED BY ANOTHER CALL OF
C               COST WILL MULTIPLY THE SEQUENCE X BY 2*(N-1).
C               HENCE COST IS THE UNNORMALIZED INVERSE
C               OF ITSELF.
C
C  WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT BE
C          DESTROYED BETWEEN CALLS OF COST.
C
C  *********************************************************************
C  *                                                                   *
C  *   SUBPROGRAM REVISION HISTORY                                     *
C  *                                                                   *
C  *   06/01/79  -  ORIGINAL VERSION BY PAUL SWARZTRAUBER.             *
C  *                DISTRIBUTED BY NCAR (REF. 1).                      *
C  *   04/01/83  -  SLATEC COMMON MATH LIBRARY SUBCOMMITTEE.           *
C  *                MODIFIED TO USE SLATEC LIBRARY SOURCE FILE FORMAT. *
C  *                DISTRIBUTED IN THE SLATEC LIBRARY (REF. 2).        *
C  *   01/15/86  -  RON BOISVERT, NATIONAL BUREAU OF STANDARDS.        *
C  *                MODIFIED TO CONVERT TO PORTABLE FORTRAN 77.        *
C  *                                                                   *
C  *   THE CHANGES INTRODUCED IN THE MOST RECENT MODIFICATION ARE      *
C  *                                                                   *
C  *   (A) DUMMY ARRAY SIZE DECLARATIONS (1) CHANGED TO (*)            *
C  *   (B) REFERENCES TO INTRINSIC FUNCTION FLOAT CHANGED TO REAL      *
C  *   (C) MATHEMATICAL CONSTANTS PREVIOUSLY CODED IN DATA STATE-      *
C  *       MENTS NOW COMPUTED AT RUNTIME USING FORTRAN INTRINSIC       *
C  *       FUNCTIONS.  THE AFFECTED VARIABLES ARE                      *
C  *                                                                   *
C  *          PI      SQRT2   SQRT3   TAUR    TR11    TR12             *
C  *          PIH     TSQRT2          TAUI    TI11    TI12             *
C  *          TPI     HSQT2                                            *
C  *                                                                   *
C  *   REFERENCES                                                      *
C  *                                                                   *
C  *   1. P.N. SWARZTRAUBER, VECTORIZING THE FFTS, IN PARALLEL         *
C  *      COMPUTATIONS (G. RODRIGUE, ED.), ACADEMIC PRESS, 1982,       *
C  *      PP. 51-83.                                                   *
C  *   2. B.L. BUZBEE, THE SLATEC COMMON MATH LIBRARY, IN SOURCES      *
C  *      AND DEVELOPMENT OF MATHEMATICAL SOFTWARE (W. COWELL, ED.),   *
C  *      PRENTICE-HALL, 1984, PP. 302-318.                            *
C  *                                                                   *
C  *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  RFFTF
C***END PROLOGUE  COST
      DIMENSION       X(*)       ,WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  COST
      NM1 = N-1
      NP1 = N+1
      NS2 = N/2
      IF (N-2) 106,101,102
  101 X1H = X(1)+X(2)
      X(2) = X(1)-X(2)
      X(1) = X1H
      RETURN
  102 IF (N .GT. 3) GO TO 103
      X1P3 = X(1)+X(3)
      TX2 = X(2)+X(2)
      X(2) = X(1)-X(3)
      X(1) = X1P3+TX2
      X(3) = X1P3-TX2
      RETURN
  103 C1 = X(1)-X(N)
      X(1) = X(1)+X(N)
      DO 104 K=2,NS2
         KC = NP1-K
         T1 = X(K)+X(KC)
         T2 = X(K)-X(KC)
         C1 = C1+WSAVE(KC)*T2
         T2 = WSAVE(K)*T2
         X(K) = T1-T2
         X(KC) = T1+T2
  104 CONTINUE
      MODN = MOD(N,2)
      IF (MODN .NE. 0) X(NS2+1) = X(NS2+1)+X(NS2+1)
      CALL RFFTF (NM1,X,WSAVE(N+1))
      XIM2 = X(2)
      X(2) = C1
      DO 105 I=4,N,2
         XI = X(I)
         X(I) = X(I-2)-X(I-1)
         X(I-1) = XIM2
         XIM2 = XI
  105 CONTINUE
      IF (MODN .NE. 0) X(N) = XIM2
  106 RETURN
      END
      SUBROUTINE RFFTB(N,R,WSAVE)
C***BEGIN PROLOGUE  RFFTB
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  860115   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  BACKWARD TRANSFORM OF A REAL COEFFICIENT ARRAY.
C***DESCRIPTION
C
C  SUBROUTINE RFFTB COMPUTES THE REAL PERODIC SEQUENCE FROM ITS
C  FOURIER COEFFICIENTS (FOURIER SYNTHESIS).  THE TRANSFORM IS DEFINED
C  BELOW AT OUTPUT PARAMETER R.
C
C  INPUT PARAMETERS
C
C  N       THE LENGTH OF THE ARRAY R TO BE TRANSFORMED.  THE METHOD
C          IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
C          N MAY CHANGE SO LONG AS DIFFERENT WORK ARRAYS ARE PROVIDED.
C
C  R       A REAL ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE
C          TO BE TRANSFORMED
C
C  WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 2*N+15
C          IN THE PROGRAM THAT CALLS RFFTB.  THE WSAVE ARRAY MUST BE
C          INITIALIZED BY CALLING SUBROUTINE RFFTI(N,WSAVE), AND A
C          DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
C          VALUE OF N.  THIS INITIALIZATION DOES NOT HAVE TO BE
C          REPEATED SO LONG AS N REMAINS UNCHANGED.  THUS SUBSEQUENT
C          TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
C          THE SAME WSAVE ARRAY CAN BE USED BY RFFTF AND RFFTB.
C
C
C  OUTPUT PARAMETERS
C
C  R       FOR N EVEN AND FOR I = 1,...,N
C
C               R(I) = R(1)+(-1)**(I-1)*R(N)
C
C                    PLUS THE SUM FROM K=2 TO K=N/2 OF
C
C                     2.*R(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
C
C                    -2.*R(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
C
C          FOR N ODD AND FOR I = 1,...,N
C
C               R(I) = R(1) PLUS THE SUM FROM K=2 TO K=(N+1)/2 OF
C
C                    2.*R(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
C
C                   -2.*R(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
C
C   *****  NOTE:
C               THIS TRANSFORM IS UNNORMALIZED SINCE A CALL OF RFFTF
C               FOLLOWED BY A CALL OF RFFTB WILL MULTIPLY THE INPUT
C               SEQUENCE BY N.
C
C  WSAVE   CONTAINS RESULTS WHICH MUST NOT BE DESTROYED BETWEEN
C          CALLS OF RFFTB OR RFFTF.
C
C  *********************************************************************
C  *                                                                   *
C  *   SUBPROGRAM REVISION HISTORY                                     *
C  *                                                                   *
C  *   06/01/79  -  ORIGINAL VERSION BY PAUL SWARZTRAUBER.             *
C  *                DISTRIBUTED BY NCAR (REF. 1).                      *
C  *   04/01/83  -  SLATEC COMMON MATH LIBRARY SUBCOMMITTEE.           *
C  *                MODIFIED TO USE SLATEC LIBRARY SOURCE FILE FORMAT. *
C  *                DISTRIBUTED IN THE SLATEC LIBRARY (REF. 2).        *
C  *   01/15/86  -  RON BOISVERT, NATIONAL BUREAU OF STANDARDS.        *
C  *                MODIFIED TO CONVERT TO PORTABLE FORTRAN 77.        *
C  *                                                                   *
C  *   THE CHANGES INTRODUCED IN THE MOST RECENT MODIFICATION ARE      *
C  *                                                                   *
C  *   (A) DUMMY ARRAY SIZE DECLARATIONS (1) CHANGED TO (*)            *
C  *   (B) REFERENCES TO INTRINSIC FUNCTION FLOAT CHANGED TO REAL      *
C  *   (C) MATHEMATICAL CONSTANTS PREVIOUSLY CODED IN DATA STATE-      *
C  *       MENTS NOW COMPUTED AT RUNTIME USING FORTRAN INTRINSIC       *
C  *       FUNCTIONS.  THE AFFECTED VARIABLES ARE                      *
C  *                                                                   *
C  *          PI      SQRT2   SQRT3   TAUR    TR11    TR12             *
C  *          PIH     TSQRT2          TAUI    TI11    TI12             *
C  *          TPI     HSQT2                                            *
C  *                                                                   *
C  *   REFERENCES                                                      *
C  *                                                                   *
C  *   1. P.N. SWARZTRAUBER, VECTORIZING THE FFTS, IN PARALLEL         *
C  *      COMPUTATIONS (G. RODRIGUE, ED.), ACADEMIC PRESS, 1982,       *
C  *      PP. 51-83.                                                   *
C  *   2. B.L. BUZBEE, THE SLATEC COMMON MATH LIBRARY, IN SOURCES      *
C  *      AND DEVELOPMENT OF MATHEMATICAL SOFTWARE (W. COWELL, ED.),   *
C  *      PRENTICE-HALL, 1984, PP. 302-318.                            *
C  *                                                                   *
C  *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  RFFTB1
C***END PROLOGUE  RFFTB
      DIMENSION       R(*)       ,WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  RFFTB
      IF (N .EQ. 1) RETURN
      CALL RFFTB1 (N,R,WSAVE,WSAVE(N+1),WSAVE(2*N+1))
      RETURN
      END
      SUBROUTINE RFFTF(N,R,WSAVE)
C***BEGIN PROLOGUE  RFFTF
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  860115   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  FORWARD TRANSFORM OF A REAL, PERIODIC SEQUENCE.
C***DESCRIPTION
C
C  SUBROUTINE RFFTF COMPUTES THE FOURIER COEFFICIENTS OF A REAL
C  PERODIC SEQUENCE (FOURIER ANALYSIS).  THE TRANSFORM IS DEFINED
C  BELOW AT OUTPUT PARAMETER R.
C
C  INPUT PARAMETERS
C
C  N       THE LENGTH OF THE ARRAY R TO BE TRANSFORMED.  THE METHOD
C          IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
C          N MAY CHANGE SO LONG AS DIFFERENT WORK ARRAYS ARE PROVIDED
C
C  R       A REAL ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE
C          TO BE TRANSFORMED
C
C  WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 2*N+15
C          IN THE PROGRAM THAT CALLS RFFTF.  THE WSAVE ARRAY MUST BE
C          INITIALIZED BY CALLING SUBROUTINE RFFTI(N,WSAVE), AND A
C          DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
C          VALUE OF N.  THIS INITIALIZATION DOES NOT HAVE TO BE
C          REPEATED SO LONG AS N REMAINS UNCHANGED.  THUS SUBSEQUENT
C          TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
C          THE SAME WSAVE ARRAY CAN BE USED BY RFFTF AND RFFTB.
C
C
C  OUTPUT PARAMETERS
C
C  R       R(1) = THE SUM FROM I=1 TO I=N OF R(I)
C
C          IF N IS EVEN SET L = N/2; IF N IS ODD SET L = (N+1)/2
C
C            THEN FOR K = 2,...,L
C
C               R(2*K-2) = THE SUM FROM I = 1 TO I = N OF
C
C                    R(I)*COS((K-1)*(I-1)*2*PI/N)
C
C               R(2*K-1) = THE SUM FROM I = 1 TO I = N OF
C
C                   -R(I)*SIN((K-1)*(I-1)*2*PI/N)
C
C          IF N IS EVEN
C
C               R(N) = THE SUM FROM I = 1 TO I = N OF
C
C                    (-1)**(I-1)*R(I)
C
C   *****  NOTE:
C               THIS TRANSFORM IS UNNORMALIZED SINCE A CALL OF RFFTF
C               FOLLOWED BY A CALL OF RFFTB WILL MULTIPLY THE INPUT
C               SEQUENCE BY N.
C
C  WSAVE   CONTAINS RESULTS WHICH MUST NOT BE DESTROYED BETWEEN
C          CALLS OF RFFTF OR RFFTB.
C
C  *********************************************************************
C  *                                                                   *
C  *   SUBPROGRAM REVISION HISTORY                                     *
C  *                                                                   *
C  *   06/01/79  -  ORIGINAL VERSION BY PAUL SWARZTRAUBER.             *
C  *                DISTRIBUTED BY NCAR (REF. 1).                      *
C  *   04/01/83  -  SLATEC COMMON MATH LIBRARY SUBCOMMITTEE.           *
C  *                MODIFIED TO USE SLATEC LIBRARY SOURCE FILE FORMAT. *
C  *                DISTRIBUTED IN THE SLATEC LIBRARY (REF. 2).        *
C  *   01/15/86  -  RON BOISVERT, NATIONAL BUREAU OF STANDARDS.        *
C  *                MODIFIED TO CONVERT TO PORTABLE FORTRAN 77.        *
C  *                                                                   *
C  *   THE CHANGES INTRODUCED IN THE MOST RECENT MODIFICATION ARE      *
C  *                                                                   *
C  *   (A) DUMMY ARRAY SIZE DECLARATIONS (1) CHANGED TO (*)            *
C  *   (B) REFERENCES TO INTRINSIC FUNCTION FLOAT CHANGED TO REAL      *
C  *   (C) MATHEMATICAL CONSTANTS PREVIOUSLY CODED IN DATA STATE-      *
C  *       MENTS NOW COMPUTED AT RUNTIME USING FORTRAN INTRINSIC       *
C  *       FUNCTIONS.  THE AFFECTED VARIABLES ARE                      *
C  *                                                                   *
C  *          PI      SQRT2   SQRT3   TAUR    TR11    TR12             *
C  *          PIH     TSQRT2          TAUI    TI11    TI12             *
C  *          TPI     HSQT2                                            *
C  *                                                                   *
C  *   REFERENCES                                                      *
C  *                                                                   *
C  *   1. P.N. SWARZTRAUBER, VECTORIZING THE FFTS, IN PARALLEL         *
C  *      COMPUTATIONS (G. RODRIGUE, ED.), ACADEMIC PRESS, 1982,       *
C  *      PP. 51-83.                                                   *
C  *   2. B.L. BUZBEE, THE SLATEC COMMON MATH LIBRARY, IN SOURCES      *
C  *      AND DEVELOPMENT OF MATHEMATICAL SOFTWARE (W. COWELL, ED.),   *
C  *      PRENTICE-HALL, 1984, PP. 302-318.                            *
C  *                                                                   *
C  *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  RFFTF1
C***END PROLOGUE  RFFTF
      DIMENSION       R(*)       ,WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  RFFTF
      IF (N .EQ. 1) RETURN
      CALL RFFTF1 (N,R,WSAVE,WSAVE(N+1),WSAVE(2*N+1))
      RETURN
      END
      SUBROUTINE RFFTI(N,WSAVE)
C***BEGIN PROLOGUE  RFFTI
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  860115   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  INITIALIZE FOR RFFTF AND RFFTB.
C***DESCRIPTION
C
C  SUBROUTINE RFFTI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
C  BOTH RFFTF AND RFFTB.  THE PRIME FACTORIZATION OF N TOGETHER WITH
C  A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
C  STORED IN WSAVE.
C
C  INPUT PARAMETER
C
C  N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED.
C
C  OUTPUT PARAMETER
C
C  WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 2*N+15.
C          THE SAME WORK ARRAY CAN BE USED FOR BOTH RFFTF AND RFFTB
C          AS LONG AS N REMAINS UNCHANGED.  DIFFERENT WSAVE ARRAYS
C          ARE REQUIRED FOR DIFFERENT VALUES OF N.  THE CONTENTS OF
C          WSAVE MUST NOT BE CHANGED BETWEEN CALLS OF RFFTF OR RFFTB.
C
C  *********************************************************************
C  *                                                                   *
C  *   SUBPROGRAM REVISION HISTORY                                     *
C  *                                                                   *
C  *   06/01/79  -  ORIGINAL VERSION BY PAUL SWARZTRAUBER.             *
C  *                DISTRIBUTED BY NCAR (REF. 1).                      *
C  *   04/01/83  -  SLATEC COMMON MATH LIBRARY SUBCOMMITTEE.           *
C  *                MODIFIED TO USE SLATEC LIBRARY SOURCE FILE FORMAT. *
C  *                DISTRIBUTED IN THE SLATEC LIBRARY (REF. 2).        *
C  *   01/15/86  -  RON BOISVERT, NATIONAL BUREAU OF STANDARDS.        *
C  *                MODIFIED TO CONVERT TO PORTABLE FORTRAN 77.        *
C  *                                                                   *
C  *   THE CHANGES INTRODUCED IN THE MOST RECENT MODIFICATION ARE      *
C  *                                                                   *
C  *   (A) DUMMY ARRAY SIZE DECLARATIONS (1) CHANGED TO (*)            *
C  *   (B) REFERENCES TO INTRINSIC FUNCTION FLOAT CHANGED TO REAL      *
C  *   (C) MATHEMATICAL CONSTANTS PREVIOUSLY CODED IN DATA STATE-      *
C  *       MENTS NOW COMPUTED AT RUNTIME USING FORTRAN INTRINSIC       *
C  *       FUNCTIONS.  THE AFFECTED VARIABLES ARE                      *
C  *                                                                   *
C  *          PI      SQRT2   SQRT3   TAUR    TR11    TR12             *
C  *          PIH     TSQRT2          TAUI    TI11    TI12             *
C  *          TPI     HSQT2                                            *
C  *                                                                   *
C  *   REFERENCES                                                      *
C  *                                                                   *
C  *   1. P.N. SWARZTRAUBER, VECTORIZING THE FFTS, IN PARALLEL         *
C  *      COMPUTATIONS (G. RODRIGUE, ED.), ACADEMIC PRESS, 1982,       *
C  *      PP. 51-83.                                                   *
C  *   2. B.L. BUZBEE, THE SLATEC COMMON MATH LIBRARY, IN SOURCES      *
C  *      AND DEVELOPMENT OF MATHEMATICAL SOFTWARE (W. COWELL, ED.),   *
C  *      PRENTICE-HALL, 1984, PP. 302-318.                            *
C  *                                                                   *
C  *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  RFFTI1
C***END PROLOGUE  RFFTI
      DIMENSION       WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  RFFTI
      IF (N .EQ. 1) RETURN
      CALL RFFTI1 (N,WSAVE(N+1),WSAVE(2*N+1))
      RETURN
      END
      SUBROUTINE SINQB(N,X,WSAVE)
C***BEGIN PROLOGUE  SINQB
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  860115   (YYMMDD)
C***CATEGORY NO.  J1A3
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  UNNORMALIZED INVERSE OF SINQF.
C***DESCRIPTION
C
C  SUBROUTINE SINQB COMPUTES THE FAST FOURIER TRANSFORM OF QUARTER
C  WAVE DATA.  THAT IS, SINQB COMPUTES A SEQUENCE FROM ITS
C  REPRESENTATION IN TERMS OF A SINE SERIES WITH ODD WAVE NUMBERS.
C  THE TRANSFORM IS DEFINED BELOW AT OUTPUT PARAMETER X.
C
C  SINQF IS THE UNNORMALIZED INVERSE OF SINQB SINCE A CALL OF SINQB
C  FOLLOWED BY A CALL OF SINQF WILL MULTIPLY THE INPUT SEQUENCE X
C  BY 4*N.
C
C  THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE SINQB MUST BE
C  INITIALIZED BY CALLING SUBROUTINE SINQI(N,WSAVE).
C
C
C  INPUT PARAMETERS
C
C  N       THE LENGTH OF THE ARRAY X TO BE TRANSFORMED.  THE METHOD
C          IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
C
C  X       AN ARRAY WHICH CONTAINS THE SEQUENCE TO BE TRANSFORMED
C
C  WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15
C          IN THE PROGRAM THAT CALLS SINQB.  THE WSAVE ARRAY MUST BE
C          INITIALIZED BY CALLING SUBROUTINE SINQI(N,WSAVE), AND A
C          DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
C          VALUE OF N.  THIS INITIALIZATION DOES NOT HAVE TO BE
C          REPEATED SO LONG AS N REMAINS UNCHANGED.  THUS SUBSEQUENT
C          TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
C
C  OUTPUT PARAMETERS
C
C  X       FOR I=1,...,N
C
C               X(I)= THE SUM FROM K=1 TO K=N OF
C
C                 4*X(K)*SIN((2K-1)*I*PI/(2*N))
C
C               A CALL OF SINQB FOLLOWED BY A CALL OF
C               SINQF WILL MULTIPLY THE SEQUENCE X BY 4*N.
C               THEREFORE SINQF IS THE UNNORMALIZED INVERSE
C               OF SINQB.
C
C  WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT
C          BE DESTROYED BETWEEN CALLS OF SINQB OR SINQF.
C
C  *********************************************************************
C  *                                                                   *
C  *   SUBPROGRAM REVISION HISTORY                                     *
C  *                                                                   *
C  *   06/01/79  -  ORIGINAL VERSION BY PAUL SWARZTRAUBER.             *
C  *                DISTRIBUTED BY NCAR (REF. 1).                      *
C  *   04/01/83  -  SLATEC COMMON MATH LIBRARY SUBCOMMITTEE.           *
C  *                MODIFIED TO USE SLATEC LIBRARY SOURCE FILE FORMAT. *
C  *                DISTRIBUTED IN THE SLATEC LIBRARY (REF. 2).        *
C  *   01/15/86  -  RON BOISVERT, NATIONAL BUREAU OF STANDARDS.        *
C  *                MODIFIED TO CONVERT TO PORTABLE FORTRAN 77.        *
C  *                                                                   *
C  *   THE CHANGES INTRODUCED IN THE MOST RECENT MODIFICATION ARE      *
C  *                                                                   *
C  *   (A) DUMMY ARRAY SIZE DECLARATIONS (1) CHANGED TO (*)            *
C  *   (B) REFERENCES TO INTRINSIC FUNCTION FLOAT CHANGED TO REAL      *
C  *   (C) MATHEMATICAL CONSTANTS PREVIOUSLY CODED IN DATA STATE-      *
C  *       MENTS NOW COMPUTED AT RUNTIME USING FORTRAN INTRINSIC       *
C  *       FUNCTIONS.  THE AFFECTED VARIABLES ARE                      *
C  *                                                                   *
C  *          PI      SQRT2   SQRT3   TAUR    TR11    TR12             *
C  *          PIH     TSQRT2          TAUI    TI11    TI12             *
C  *          TPI     HSQT2                                            *
C  *                                                                   *
C  *   REFERENCES                                                      *
C  *                                                                   *
C  *   1. P.N. SWARZTRAUBER, VECTORIZING THE FFTS, IN PARALLEL         *
C  *      COMPUTATIONS (G. RODRIGUE, ED.), ACADEMIC PRESS, 1982,       *
C  *      PP. 51-83.                                                   *
C  *   2. B.L. BUZBEE, THE SLATEC COMMON MATH LIBRARY, IN SOURCES      *
C  *      AND DEVELOPMENT OF MATHEMATICAL SOFTWARE (W. COWELL, ED.),   *
C  *      PRENTICE-HALL, 1984, PP. 302-318.                            *
C  *                                                                   *
C  *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  COSQB
C***END PROLOGUE  SINQB
      DIMENSION       X(*)       ,WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  SINQB
      IF (N .GT. 1) GO TO 101
      X(1) = 4.*X(1)
      RETURN
  101 NS2 = N/2
      DO 102 K=2,N,2
         X(K) = -X(K)
  102 CONTINUE
      CALL COSQB (N,X,WSAVE)
      DO 103 K=1,NS2
         KC = N-K
         XHOLD = X(K)
         X(K) = X(KC+1)
         X(KC+1) = XHOLD
  103 CONTINUE
      RETURN
      END
      SUBROUTINE SINQF(N,X,WSAVE)
C***BEGIN PROLOGUE  SINQF
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  860115   (YYMMDD)
C***CATEGORY NO.  J1A3
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  FORWARD SINE TRANSFORM WITH ODD WAVE NUMBERS.
C***DESCRIPTION
C
C  SUBROUTINE SINQF COMPUTES THE FAST FOURIER TRANSFORM OF QUARTER
C  WAVE DATA.  THAT IS, SINQF COMPUTES THE COEFFICIENTS IN A SINE
C  SERIES REPRESENTATION WITH ONLY ODD WAVE NUMBERS.  THE TRANSFORM
C  IS DEFINED BELOW AT OUTPUT PARAMETER X.
C
C  SINQB IS THE UNNORMALIZED INVERSE OF SINQF SINCE A CALL OF SINQF
C  FOLLOWED BY A CALL OF SINQB WILL MULTIPLY THE INPUT SEQUENCE X
C  BY 4*N.
C
C  THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE SINQF MUST BE
C  INITIALIZED BY CALLING SUBROUTINE SINQI(N,WSAVE).
C
C
C  INPUT PARAMETERS
C
C  N       THE LENGTH OF THE ARRAY X TO BE TRANSFORMED.  THE METHOD
C          IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
C
C  X       AN ARRAY WHICH CONTAINS THE SEQUENCE TO BE TRANSFORMED
C
C  WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15
C          IN THE PROGRAM THAT CALLS SINQF.  THE WSAVE ARRAY MUST BE
C          INITIALIZED BY CALLING SUBROUTINE SINQI(N,WSAVE), AND A
C          DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
C          VALUE OF N.  THIS INITIALIZATION DOES NOT HAVE TO BE
C          REPEATED SO LONG AS N REMAINS UNCHANGED.  THUS SUBSEQUENT
C          TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
C
C  OUTPUT PARAMETERS
C
C  X       FOR I=1,...,N
C
C               X(I) = (-1)**(I-1)*X(N)
C
C                  + THE SUM FROM K=1 TO K=N-1 OF
C
C                  2*X(K)*SIN((2*I-1)*K*PI/(2*N))
C
C               A CALL OF SINQF FOLLOWED BY A CALL OF
C               SINQB WILL MULTIPLY THE SEQUENCE X BY 4*N.
C               THEREFORE SINQB IS THE UNNORMALIZED INVERSE
C               OF SINQF.
C
C  WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT
C          BE DESTROYED BETWEEN CALLS OF SINQF OR SINQB.
C
C  *********************************************************************
C  *                                                                   *
C  *   SUBPROGRAM REVISION HISTORY                                     *
C  *                                                                   *
C  *   06/01/79  -  ORIGINAL VERSION BY PAUL SWARZTRAUBER.             *
C  *                DISTRIBUTED BY NCAR (REF. 1).                      *
C  *   04/01/83  -  SLATEC COMMON MATH LIBRARY SUBCOMMITTEE.           *
C  *                MODIFIED TO USE SLATEC LIBRARY SOURCE FILE FORMAT. *
C  *                DISTRIBUTED IN THE SLATEC LIBRARY (REF. 2).        *
C  *   01/15/86  -  RON BOISVERT, NATIONAL BUREAU OF STANDARDS.        *
C  *                MODIFIED TO CONVERT TO PORTABLE FORTRAN 77.        *
C  *                                                                   *
C  *   THE CHANGES INTRODUCED IN THE MOST RECENT MODIFICATION ARE      *
C  *                                                                   *
C  *   (A) DUMMY ARRAY SIZE DECLARATIONS (1) CHANGED TO (*)            *
C  *   (B) REFERENCES TO INTRINSIC FUNCTION FLOAT CHANGED TO REAL      *
C  *   (C) MATHEMATICAL CONSTANTS PREVIOUSLY CODED IN DATA STATE-      *
C  *       MENTS NOW COMPUTED AT RUNTIME USING FORTRAN INTRINSIC       *
C  *       FUNCTIONS.  THE AFFECTED VARIABLES ARE                      *
C  *                                                                   *
C  *          PI      SQRT2   SQRT3   TAUR    TR11    TR12             *
C  *          PIH     TSQRT2          TAUI    TI11    TI12             *
C  *          TPI     HSQT2                                            *
C  *                                                                   *
C  *   REFERENCES                                                      *
C  *                                                                   *
C  *   1. P.N. SWARZTRAUBER, VECTORIZING THE FFTS, IN PARALLEL         *
C  *      COMPUTATIONS (G. RODRIGUE, ED.), ACADEMIC PRESS, 1982,       *
C  *      PP. 51-83.                                                   *
C  *   2. B.L. BUZBEE, THE SLATEC COMMON MATH LIBRARY, IN SOURCES      *
C  *      AND DEVELOPMENT OF MATHEMATICAL SOFTWARE (W. COWELL, ED.),   *
C  *      PRENTICE-HALL, 1984, PP. 302-318.                            *
C  *                                                                   *
C  *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  COSQF
C***END PROLOGUE  SINQF
      DIMENSION       X(*)       ,WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  SINQF
      IF (N .EQ. 1) RETURN
      NS2 = N/2
      DO 101 K=1,NS2
         KC = N-K
         XHOLD = X(K)
         X(K) = X(KC+1)
         X(KC+1) = XHOLD
  101 CONTINUE
      CALL COSQF (N,X,WSAVE)
      DO 102 K=2,N,2
         X(K) = -X(K)
  102 CONTINUE
      RETURN
      END
      SUBROUTINE SINQI(N,WSAVE)
C***BEGIN PROLOGUE  SINQI
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  860115   (YYMMDD)
C***CATEGORY NO.  J1A3
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  INITIALIZE FOR SINQF AND SINQB.
C***DESCRIPTION
C
C  SUBROUTINE SINQI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
C  BOTH SINQF AND SINQB.  THE PRIME FACTORIZATION OF N TOGETHER WITH
C  A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
C  STORED IN WSAVE.
C
C  INPUT PARAMETER
C
C  N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED.  THE METHOD
C          IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
C
C  OUTPUT PARAMETER
C
C  WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
C          THE SAME WORK ARRAY CAN BE USED FOR BOTH SINQF AND SINQB
C          AS LONG AS N REMAINS UNCHANGED.  DIFFERENT WSAVE ARRAYS
C          ARE REQUIRED FOR DIFFERENT VALUES OF N.  THE CONTENTS OF
C          WSAVE MUST NOT BE CHANGED BETWEEN CALLS OF SINQF OR SINQB.
C
C  *********************************************************************
C  *                                                                   *
C  *   SUBPROGRAM REVISION HISTORY                                     *
C  *                                                                   *
C  *   06/01/79  -  ORIGINAL VERSION BY PAUL SWARZTRAUBER.             *
C  *                DISTRIBUTED BY NCAR (REF. 1).                      *
C  *   04/01/83  -  SLATEC COMMON MATH LIBRARY SUBCOMMITTEE.           *
C  *                MODIFIED TO USE SLATEC LIBRARY SOURCE FILE FORMAT. *
C  *                DISTRIBUTED IN THE SLATEC LIBRARY (REF. 2).        *
C  *   01/15/86  -  RON BOISVERT, NATIONAL BUREAU OF STANDARDS.        *
C  *                MODIFIED TO CONVERT TO PORTABLE FORTRAN 77.        *
C  *                                                                   *
C  *   THE CHANGES INTRODUCED IN THE MOST RECENT MODIFICATION ARE      *
C  *                                                                   *
C  *   (A) DUMMY ARRAY SIZE DECLARATIONS (1) CHANGED TO (*)            *
C  *   (B) REFERENCES TO INTRINSIC FUNCTION FLOAT CHANGED TO REAL      *
C  *   (C) MATHEMATICAL CONSTANTS PREVIOUSLY CODED IN DATA STATE-      *
C  *       MENTS NOW COMPUTED AT RUNTIME USING FORTRAN INTRINSIC       *
C  *       FUNCTIONS.  THE AFFECTED VARIABLES ARE                      *
C  *                                                                   *
C  *          PI      SQRT2   SQRT3   TAUR    TR11    TR12             *
C  *          PIH     TSQRT2          TAUI    TI11    TI12             *
C  *          TPI     HSQT2                                            *
C  *                                                                   *
C  *   REFERENCES                                                      *
C  *                                                                   *
C  *   1. P.N. SWARZTRAUBER, VECTORIZING THE FFTS, IN PARALLEL         *
C  *      COMPUTATIONS (G. RODRIGUE, ED.), ACADEMIC PRESS, 1982,       *
C  *      PP. 51-83.                                                   *
C  *   2. B.L. BUZBEE, THE SLATEC COMMON MATH LIBRARY, IN SOURCES      *
C  *      AND DEVELOPMENT OF MATHEMATICAL SOFTWARE (W. COWELL, ED.),   *
C  *      PRENTICE-HALL, 1984, PP. 302-318.                            *
C  *                                                                   *
C  *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  COSQI
C***END PROLOGUE  SINQI
      DIMENSION       WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  SINQI
      CALL COSQI (N,WSAVE)
      RETURN
      END
      SUBROUTINE SINTI(N,WSAVE)
C***BEGIN PROLOGUE  SINTI
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  860115   (YYMMDD)
C***CATEGORY NO.  J1A3
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  INITIALIZE FOR SINT.
C***DESCRIPTION
C
C  SUBROUTINE SINTI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
C  SUBROUTINE SINT.  THE PRIME FACTORIZATION OF N TOGETHER WITH
C  A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
C  STORED IN WSAVE.
C
C  INPUT PARAMETER
C
C  N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED.  THE METHOD
C          IS MOST EFFICIENT WHEN N+1 IS A PRODUCT OF SMALL PRIMES.
C
C  OUTPUT PARAMETER
C
C  WSAVE   A WORK ARRAY WITH AT LEAST INT(3.5*N+16) LOCATIONS.
C          DIFFERENT WSAVE ARRAYS ARE REQUIRED FOR DIFFERENT VALUES
C          OF N.  THE CONTENTS OF WSAVE MUST NOT BE CHANGED BETWEEN
C          CALLS OF SINT.
C
C  *********************************************************************
C  *                                                                   *
C  *   SUBPROGRAM REVISION HISTORY                                     *
C  *                                                                   *
C  *   06/01/79  -  ORIGINAL VERSION BY PAUL SWARZTRAUBER.             *
C  *                DISTRIBUTED BY NCAR (REF. 1).                      *
C  *   04/01/83  -  SLATEC COMMON MATH LIBRARY SUBCOMMITTEE.           *
C  *                MODIFIED TO USE SLATEC LIBRARY SOURCE FILE FORMAT. *
C  *                DISTRIBUTED IN THE SLATEC LIBRARY (REF. 2).        *
C  *   01/15/86  -  RON BOISVERT, NATIONAL BUREAU OF STANDARDS.        *
C  *                MODIFIED TO CONVERT TO PORTABLE FORTRAN 77.        *
C  *                                                                   *
C  *   THE CHANGES INTRODUCED IN THE MOST RECENT MODIFICATION ARE      *
C  *                                                                   *
C  *   (A) DUMMY ARRAY SIZE DECLARATIONS (1) CHANGED TO (*)            *
C  *   (B) REFERENCES TO INTRINSIC FUNCTION FLOAT CHANGED TO REAL      *
C  *   (C) MATHEMATICAL CONSTANTS PREVIOUSLY CODED IN DATA STATE-      *
C  *       MENTS NOW COMPUTED AT RUNTIME USING FORTRAN INTRINSIC       *
C  *       FUNCTIONS.  THE AFFECTED VARIABLES ARE                      *
C  *                                                                   *
C  *          PI      SQRT2   SQRT3   TAUR    TR11    TR12             *
C  *          PIH     TSQRT2          TAUI    TI11    TI12             *
C  *          TPI     HSQT2                                            *
C  *                                                                   *
C  *   REFERENCES                                                      *
C  *                                                                   *
C  *   1. P.N. SWARZTRAUBER, VECTORIZING THE FFTS, IN PARALLEL         *
C  *      COMPUTATIONS (G. RODRIGUE, ED.), ACADEMIC PRESS, 1982,       *
C  *      PP. 51-83.                                                   *
C  *   2. B.L. BUZBEE, THE SLATEC COMMON MATH LIBRARY, IN SOURCES      *
C  *      AND DEVELOPMENT OF MATHEMATICAL SOFTWARE (W. COWELL, ED.),   *
C  *      PRENTICE-HALL, 1984, PP. 302-318.                            *
C  *                                                                   *
C  *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  RFFTI
C***END PROLOGUE  SINTI
      DIMENSION       WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  SINTI
      PI = 4.*ATAN(1.)
      IF (N .LE. 1) RETURN
      NP1 = N+1
      NS2 = N/2
      DT = PI/REAL(NP1)
      KS = N+2
      KF = KS+NS2-1
      FK = 0.
      DO 101 K=KS,KF
         FK = FK+1.
         WSAVE(K) = 2.*SIN(FK*DT)
  101 CONTINUE
      CALL RFFTI (NP1,WSAVE(KF+1))
      RETURN
      END
      SUBROUTINE SINT(N,X,WSAVE)
C***BEGIN PROLOGUE  SINT
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  860115   (YYMMDD)
C***CATEGORY NO.  J1A3
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  SINE TRANSFORM OF A REAL, ODD SEQUENCE.
C***DESCRIPTION
C
C  SUBROUTINE SINT COMPUTES THE DISCRETE FOURIER SINE TRANSFORM
C  OF AN ODD SEQUENCE X(I).  THE TRANSFORM IS DEFINED BELOW AT
C  OUTPUT PARAMETER X.
C
C  SINT IS THE UNNORMALIZED INVERSE OF ITSELF SINCE A CALL OF SINT
C  FOLLOWED BY ANOTHER CALL OF SINT WILL MULTIPLY THE INPUT SEQUENCE
C  X BY 2*(N+1).
C
C  THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE SINT MUST BE
C  INITIALIZED BY CALLING SUBROUTINE SINTI(N,WSAVE).
C
C  INPUT PARAMETERS
C
C  N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED.  THE METHOD
C          IS MOST EFFICIENT WHEN N+1 IS THE PRODUCT OF SMALL PRIMES.
C
C  X       AN ARRAY WHICH CONTAINS THE SEQUENCE TO BE TRANSFORMED
C
C
C  WSAVE   A WORK ARRAY WITH DIMENSION AT LEAST INT(3.5*N+16)
C          IN THE PROGRAM THAT CALLS SINT.  THE WSAVE ARRAY MUST BE
C          INITIALIZED BY CALLING SUBROUTINE SINTI(N,WSAVE), AND A
C          DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
C          VALUE OF N.  THIS INITIALIZATION DOES NOT HAVE TO BE
C          REPEATED SO LONG AS N REMAINS UNCHANGED.  THUS SUBSEQUENT
C          TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
C
C  OUTPUT PARAMETERS
C
C  X       FOR I=1,...,N
C
C               X(I)= THE SUM FROM K=1 TO K=N
C
C                    2*X(K)*SIN(K*I*PI/(N+1))
C
C               A CALL OF SINT FOLLOWED BY ANOTHER CALL OF
C               SINT WILL MULTIPLY THE SEQUENCE X BY 2*(N+1).
C               HENCE SINT IS THE UNNORMALIZED INVERSE
C               OF ITSELF.
C
C  WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT BE
C          DESTROYED BETWEEN CALLS OF SINT.
C
C  *********************************************************************
C  *                                                                   *
C  *   SUBPROGRAM REVISION HISTORY                                     *
C  *                                                                   *
C  *   06/01/79  -  ORIGINAL VERSION BY PAUL SWARZTRAUBER.             *
C  *                DISTRIBUTED BY NCAR (REF. 1).                      *
C  *   04/01/83  -  SLATEC COMMON MATH LIBRARY SUBCOMMITTEE.           *
C  *                MODIFIED TO USE SLATEC LIBRARY SOURCE FILE FORMAT. *
C  *                DISTRIBUTED IN THE SLATEC LIBRARY (REF. 2).        *
C  *   01/15/86  -  RON BOISVERT, NATIONAL BUREAU OF STANDARDS.        *
C  *                MODIFIED TO CONVERT TO PORTABLE FORTRAN 77.        *
C  *                                                                   *
C  *   THE CHANGES INTRODUCED IN THE MOST RECENT MODIFICATION ARE      *
C  *                                                                   *
C  *   (A) DUMMY ARRAY SIZE DECLARATIONS (1) CHANGED TO (*)            *
C  *   (B) REFERENCES TO INTRINSIC FUNCTION FLOAT CHANGED TO REAL      *
C  *   (C) MATHEMATICAL CONSTANTS PREVIOUSLY CODED IN DATA STATE-      *
C  *       MENTS NOW COMPUTED AT RUNTIME USING FORTRAN INTRINSIC       *
C  *       FUNCTIONS.  THE AFFECTED VARIABLES ARE                      *
C  *                                                                   *
C  *          PI      SQRT2   SQRT3   TAUR    TR11    TR12             *
C  *          PIH     TSQRT2          TAUI    TI11    TI12             *
C  *          TPI     HSQT2                                            *
C  *                                                                   *
C  *   REFERENCES                                                      *
C  *                                                                   *
C  *   1. P.N. SWARZTRAUBER, VECTORIZING THE FFTS, IN PARALLEL         *
C  *      COMPUTATIONS (G. RODRIGUE, ED.), ACADEMIC PRESS, 1982,       *
C  *      PP. 51-83.                                                   *
C  *   2. B.L. BUZBEE, THE SLATEC COMMON MATH LIBRARY, IN SOURCES      *
C  *      AND DEVELOPMENT OF MATHEMATICAL SOFTWARE (W. COWELL, ED.),   *
C  *      PRENTICE-HALL, 1984, PP. 302-318.                            *
C  *                                                                   *
C  *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  RFFTF
C***END PROLOGUE  SINT
      DIMENSION       X(*)       ,WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  SINT
      SQRT3 = SQRT(3.)
      IF (N-2) 101,102,103
  101 X(1) = X(1)+X(1)
      RETURN
  102 XH = SQRT3*(X(1)+X(2))
      X(2) = SQRT3*(X(1)-X(2))
      X(1) = XH
      RETURN
  103 NP1 = N+1
      NS2 = N/2
      WSAVE(1) = 0.
      KW = NP1
      DO 104 K=1,NS2
1        KW = KW+1
         KC = NP1-K
         T1 = X(K)-X(KC)
         T2 = WSAVE(KW)*(X(K)+X(KC))
         WSAVE(K+1) = T1+T2
         WSAVE(KC+1) = T2-T1
  104 CONTINUE
      MODN = MOD(N,2)
      IF (MODN .NE. 0) WSAVE(NS2+2) = 4.*X(NS2+1)
      NF = NP1+NS2+1
      CALL RFFTF (NP1,WSAVE,WSAVE(NF))
      X(1) = .5*WSAVE(1)
      DO 105 I=3,N,2
         X(I-1) = -WSAVE(I)
         X(I) = X(I-2)+WSAVE(I-1)
  105 CONTINUE
      IF (MODN .NE. 0) RETURN
      X(N) = -WSAVE(N+1)
      RETURN
      END
      SUBROUTINE COSQB1(N,X,W,XH)
C***BEGIN PROLOGUE  COSQB1
C***REFER TO  COSQB
C***ROUTINES CALLED  RFFTB
C***END PROLOGUE  COSQB1
      DIMENSION       X(*)       ,W(*)       ,XH(*)
C***FIRST EXECUTABLE STATEMENT  COSQB1
      NS2 = (N+1)/2
      NP2 = N+2
      DO 101 I=3,N,2
         XIM1 = X(I-1)+X(I)
         X(I) = X(I)-X(I-1)
         X(I-1) = XIM1
  101 CONTINUE
      X(1) = X(1)+X(1)
      MODN = MOD(N,2)
      IF (MODN .EQ. 0) X(N) = X(N)+X(N)
      CALL RFFTB (N,X,XH)
      DO 102 K=2,NS2
         KC = NP2-K
         XH(K) = W(K-1)*X(KC)+W(KC-1)*X(K)
         XH(KC) = W(K-1)*X(K)-W(KC-1)*X(KC)
  102 CONTINUE
      IF (MODN .EQ. 0) X(NS2+1) = W(NS2)*(X(NS2+1)+X(NS2+1))
      DO 103 K=2,NS2
         KC = NP2-K
         X(K) = XH(K)+XH(KC)
         X(KC) = XH(K)-XH(KC)
  103 CONTINUE
      X(1) = X(1)+X(1)
      RETURN
      END
      SUBROUTINE COSQF1(N,X,W,XH)
C***BEGIN PROLOGUE  COSQF1
C***REFER TO  COSQF
C***ROUTINES CALLED  RFFTF
C***END PROLOGUE  COSQF1
      DIMENSION       X(*)       ,W(*)       ,XH(*)
C***FIRST EXECUTABLE STATEMENT  COSQF1
      NS2 = (N+1)/2
      NP2 = N+2
      DO 101 K=2,NS2
         KC = NP2-K
         XH(K) = X(K)+X(KC)
         XH(KC) = X(K)-X(KC)
  101 CONTINUE
      MODN = MOD(N,2)
      IF (MODN .EQ. 0) XH(NS2+1) = X(NS2+1)+X(NS2+1)
      DO 102 K=2,NS2
         KC = NP2-K
         X(K) = W(K-1)*XH(KC)+W(KC-1)*XH(K)
         X(KC) = W(K-1)*XH(K)-W(KC-1)*XH(KC)
  102 CONTINUE
      IF (MODN .EQ. 0) X(NS2+1) = W(NS2)*XH(NS2+1)
      CALL RFFTF (N,X,XH)
      DO 103 I=3,N,2
         XIM1 = X(I-1)-X(I)
         X(I) = X(I-1)+X(I)
         X(I-1) = XIM1
  103 CONTINUE
      RETURN
      END
      SUBROUTINE RADBG(IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
C***BEGIN PROLOGUE  RADBG
C***REFER TO  RFFTB
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  RADBG
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,
     1                C1(IDO,L1,IP)          ,C2(IDL1,IP),
     2                CH2(IDL1,IP)           ,WA(*)
C***FIRST EXECUTABLE STATEMENT  RADBG
      TPI = 8.*ATAN(1.)
      ARG = TPI/REAL(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IF (IDO .LT. L1) GO TO 103
      DO 102 K=1,L1
         DO 101 I=1,IDO
            CH(I,K,1) = CC(I,1,K)
  101    CONTINUE
  102 CONTINUE
      GO TO 106
  103 DO 105 I=1,IDO
         DO 104 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
  106 DO 108 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 107 K=1,L1
            CH(1,K,J) = CC(IDO,J2-2,K)+CC(IDO,J2-2,K)
            CH(1,K,JC) = CC(1,J2-1,K)+CC(1,J2-1,K)
  107    CONTINUE
  108 CONTINUE
      IF (IDO .EQ. 1) GO TO 116
      IF (NBD .LT. L1) GO TO 112
      DO 111 J=2,IPPH
         JC = IPP2-J
         DO 110 K=1,L1
CDIR$ IVDEP
            DO 109 I=3,IDO,2
               IC = IDP2-I
               CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)
               CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)
               CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)
               CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)
  109       CONTINUE
  110    CONTINUE
  111 CONTINUE
      GO TO 116
  112 DO 115 J=2,IPPH
         JC = IPP2-J
CDIR$ IVDEP
         DO 114 I=3,IDO,2
            IC = IDP2-I
            DO 113 K=1,L1
               CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)
               CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)
               CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)
               CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)
  113       CONTINUE
  114    CONTINUE
  115 CONTINUE
  116 AR1 = 1.
      AI1 = 0.
      DO 120 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 117 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+AR1*CH2(IK,2)
            C2(IK,LC) = AI1*CH2(IK,IP)
  117    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 119 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 118 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+AR2*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)+AI2*CH2(IK,JC)
  118       CONTINUE
  119    CONTINUE
  120 CONTINUE
      DO 122 J=2,IPPH
         DO 121 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  121    CONTINUE
  122 CONTINUE
      DO 124 J=2,IPPH
         JC = IPP2-J
         DO 123 K=1,L1
            CH(1,K,J) = C1(1,K,J)-C1(1,K,JC)
            CH(1,K,JC) = C1(1,K,J)+C1(1,K,JC)
  123    CONTINUE
  124 CONTINUE
      IF (IDO .EQ. 1) GO TO 132
      IF (NBD .LT. L1) GO TO 128
      DO 127 J=2,IPPH
         JC = IPP2-J
         DO 126 K=1,L1
CDIR$ IVDEP
            DO 125 I=3,IDO,2
               CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)
               CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)
               CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)
               CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      GO TO 132
  128 DO 131 J=2,IPPH
         JC = IPP2-J
         DO 130 I=3,IDO,2
            DO 129 K=1,L1
               CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)
               CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)
               CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)
               CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)
  129       CONTINUE
  130    CONTINUE
  131 CONTINUE
  132 CONTINUE
      IF (IDO .EQ. 1) RETURN
      DO 133 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  133 CONTINUE
      DO 135 J=2,IP
         DO 134 K=1,L1
            C1(1,K,J) = CH(1,K,J)
  134    CONTINUE
  135 CONTINUE
      IF (NBD .GT. L1) GO TO 139
      IS = -IDO
      DO 138 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 137 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 136 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  136       CONTINUE
  137    CONTINUE
  138 CONTINUE
      GO TO 143
  139 IS = -IDO
      DO 142 J=2,IP
         IS = IS+IDO
         DO 141 K=1,L1
            IDIJ = IS
CDIR$ IVDEP
            DO 140 I=3,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  140       CONTINUE
  141    CONTINUE
  142 CONTINUE
  143 RETURN
      END
      SUBROUTINE RADB2(IDO,L1,CC,CH,WA1)
C***BEGIN PROLOGUE  RADB2
C***REFER TO  RFFTB
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  RADB2
      DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           ,
     1                WA1(*)
C***FIRST EXECUTABLE STATEMENT  RADB2
      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(IDO,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(IDO,2,K)
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      IF((IDO-1)/2.LT.L1) GO TO 108
      DO 104 K=1,L1
CDIR$ IVDEP
         DO 103 I=3,IDO,2
            IC = IDP2-I
            CH(I-1,K,1) = CC(I-1,1,K)+CC(IC-1,2,K)
            TR2 = CC(I-1,1,K)-CC(IC-1,2,K)
            CH(I,K,1) = CC(I,1,K)-CC(IC,2,K)
            TI2 = CC(I,1,K)+CC(IC,2,K)
            CH(I-1,K,2) = WA1(I-2)*TR2-WA1(I-1)*TI2
            CH(I,K,2) = WA1(I-2)*TI2+WA1(I-1)*TR2
  103    CONTINUE
  104 CONTINUE
      GO TO 111
  108 DO 110 I=3,IDO,2
         IC = IDP2-I
CDIR$ IVDEP
         DO 109 K=1,L1
            CH(I-1,K,1) = CC(I-1,1,K)+CC(IC-1,2,K)
            TR2 = CC(I-1,1,K)-CC(IC-1,2,K)
            CH(I,K,1) = CC(I,1,K)-CC(IC,2,K)
            TI2 = CC(I,1,K)+CC(IC,2,K)
            CH(I-1,K,2) = WA1(I-2)*TR2-WA1(I-1)*TI2
            CH(I,K,2) = WA1(I-2)*TI2+WA1(I-1)*TR2
  109    CONTINUE
  110 CONTINUE
  111 IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         CH(IDO,K,1) = CC(IDO,1,K)+CC(IDO,1,K)
         CH(IDO,K,2) = -(CC(1,2,K)+CC(1,2,K))
  106 CONTINUE
  107 RETURN
      END
      SUBROUTINE RADB3(IDO,L1,CC,CH,WA1,WA2)
C***BEGIN PROLOGUE  RADB3
C***REFER TO  RFFTB
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  RADB3
      DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           ,
     1                WA1(*)     ,WA2(*)
C***FIRST EXECUTABLE STATEMENT  RADB3
      TAUR = -.5
      TAUI = .5*SQRT(3.)
      DO 101 K=1,L1
         TR2 = CC(IDO,2,K)+CC(IDO,2,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         CI3 = TAUI*(CC(1,3,K)+CC(1,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      IF((IDO-1)/2.LT.L1) GO TO 104
      DO 103 K=1,L1
CDIR$ IVDEP
         DO 102 I=3,IDO,2
            IC = IDP2-I
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,3,K)-CC(IC,2,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,3,K)-CC(IC-1,2,K))
            CI3 = TAUI*(CC(I,3,K)+CC(IC,2,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
  102    CONTINUE
  103 CONTINUE
      RETURN
  104 DO 106 I=3,IDO,2
         IC = IDP2-I
CDIR$ IVDEP
         DO 105 K=1,L1
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,3,K)-CC(IC,2,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,3,K)-CC(IC-1,2,K))
            CI3 = TAUI*(CC(I,3,K)+CC(IC,2,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
  105    CONTINUE
  106 CONTINUE
      RETURN
      END
      SUBROUTINE RADB4(IDO,L1,CC,CH,WA1,WA2,WA3)
C***BEGIN PROLOGUE  RADB4
C***REFER TO  RFFTB
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  RADB4
      DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           ,
     1                WA1(*)     ,WA2(*)     ,WA3(*)
C***FIRST EXECUTABLE STATEMENT  RADB4
      SQRT2 = SQRT(2.)
      DO 101 K=1,L1
         TR1 = CC(1,1,K)-CC(IDO,4,K)
         TR2 = CC(1,1,K)+CC(IDO,4,K)
         TR3 = CC(IDO,2,K)+CC(IDO,2,K)
         TR4 = CC(1,3,K)+CC(1,3,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,2) = TR1-TR4
         CH(1,K,3) = TR2-TR3
         CH(1,K,4) = TR1+TR4
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      IF((IDO-1)/2.LT.L1) GO TO 108
      DO 104 K=1,L1
CDIR$ IVDEP
         DO 103 I=3,IDO,2
            IC = IDP2-I
            TI1 = CC(I,1,K)+CC(IC,4,K)
            TI2 = CC(I,1,K)-CC(IC,4,K)
            TI3 = CC(I,3,K)-CC(IC,2,K)
            TR4 = CC(I,3,K)+CC(IC,2,K)
            TR1 = CC(I-1,1,K)-CC(IC-1,4,K)
            TR2 = CC(I-1,1,K)+CC(IC-1,4,K)
            TI4 = CC(I-1,3,K)-CC(IC-1,2,K)
            TR3 = CC(I-1,3,K)+CC(IC-1,2,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1-TR4
            CR4 = TR1+TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-2)*CR2-WA1(I-1)*CI2
            CH(I,K,2) = WA1(I-2)*CI2+WA1(I-1)*CR2
            CH(I-1,K,3) = WA2(I-2)*CR3-WA2(I-1)*CI3
            CH(I,K,3) = WA2(I-2)*CI3+WA2(I-1)*CR3
            CH(I-1,K,4) = WA3(I-2)*CR4-WA3(I-1)*CI4
            CH(I,K,4) = WA3(I-2)*CI4+WA3(I-1)*CR4
  103    CONTINUE
  104 CONTINUE
      GO TO 111
  108 DO 110 I=3,IDO,2
         IC = IDP2-I
CDIR$ IVDEP
         DO 109 K=1,L1
            TI1 = CC(I,1,K)+CC(IC,4,K)
            TI2 = CC(I,1,K)-CC(IC,4,K)
            TI3 = CC(I,3,K)-CC(IC,2,K)
            TR4 = CC(I,3,K)+CC(IC,2,K)
            TR1 = CC(I-1,1,K)-CC(IC-1,4,K)
            TR2 = CC(I-1,1,K)+CC(IC-1,4,K)
            TI4 = CC(I-1,3,K)-CC(IC-1,2,K)
            TR3 = CC(I-1,3,K)+CC(IC-1,2,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1-TR4
            CR4 = TR1+TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-2)*CR2-WA1(I-1)*CI2
            CH(I,K,2) = WA1(I-2)*CI2+WA1(I-1)*CR2
            CH(I-1,K,3) = WA2(I-2)*CR3-WA2(I-1)*CI3
            CH(I,K,3) = WA2(I-2)*CI3+WA2(I-1)*CR3
            CH(I-1,K,4) = WA3(I-2)*CR4-WA3(I-1)*CI4
            CH(I,K,4) = WA3(I-2)*CI4+WA3(I-1)*CR4
  109    CONTINUE
  110 CONTINUE
  111 IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         TI1 = CC(1,2,K)+CC(1,4,K)
         TI2 = CC(1,4,K)-CC(1,2,K)
         TR1 = CC(IDO,1,K)-CC(IDO,3,K)
         TR2 = CC(IDO,1,K)+CC(IDO,3,K)
         CH(IDO,K,1) = TR2+TR2
         CH(IDO,K,2) = SQRT2*(TR1-TI1)
         CH(IDO,K,3) = TI2+TI2
         CH(IDO,K,4) = -SQRT2*(TR1+TI1)
  106 CONTINUE
  107 RETURN
      END
      SUBROUTINE RADB5(IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
C***BEGIN PROLOGUE  RADB5
C***REFER TO  RFFTB
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  RADB5
      DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           ,
     1                WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)
C***FIRST EXECUTABLE STATEMENT  RADB5
      PI = 4.*ATAN(1.)
      TR11 = SIN(.1*PI)
      TI11 = SIN(.4*PI)
      TR12 = -SIN(.3*PI)
      TI12 = SIN(.2*PI)
      DO 101 K=1,L1
         TI5 = CC(1,3,K)+CC(1,3,K)
         TI4 = CC(1,5,K)+CC(1,5,K)
         TR2 = CC(IDO,2,K)+CC(IDO,2,K)
         TR3 = CC(IDO,4,K)+CC(IDO,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI5 = TI11*TI5+TI12*TI4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(1,K,5) = CR2+CI5
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      IF((IDO-1)/2.LT.L1) GO TO 104
      DO 103 K=1,L1
CDIR$ IVDEP
         DO 102 I=3,IDO,2
            IC = IDP2-I
            TI5 = CC(I,3,K)+CC(IC,2,K)
            TI2 = CC(I,3,K)-CC(IC,2,K)
            TI4 = CC(I,5,K)+CC(IC,4,K)
            TI3 = CC(I,5,K)-CC(IC,4,K)
            TR5 = CC(I-1,3,K)-CC(IC-1,2,K)
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
            TR4 = CC(I-1,5,K)-CC(IC-1,4,K)
            TR3 = CC(I-1,5,K)+CC(IC-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
            CH(I-1,K,4) = WA3(I-2)*DR4-WA3(I-1)*DI4
            CH(I,K,4) = WA3(I-2)*DI4+WA3(I-1)*DR4
            CH(I-1,K,5) = WA4(I-2)*DR5-WA4(I-1)*DI5
            CH(I,K,5) = WA4(I-2)*DI5+WA4(I-1)*DR5
  102    CONTINUE
  103 CONTINUE
      RETURN
  104 DO 106 I=3,IDO,2
         IC = IDP2-I
CDIR$ IVDEP
         DO 105 K=1,L1
            TI5 = CC(I,3,K)+CC(IC,2,K)
            TI2 = CC(I,3,K)-CC(IC,2,K)
            TI4 = CC(I,5,K)+CC(IC,4,K)
            TI3 = CC(I,5,K)-CC(IC,4,K)
            TR5 = CC(I-1,3,K)-CC(IC-1,2,K)
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
            TR4 = CC(I-1,5,K)-CC(IC-1,4,K)
            TR3 = CC(I-1,5,K)+CC(IC-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
            CH(I-1,K,4) = WA3(I-2)*DR4-WA3(I-1)*DI4
            CH(I,K,4) = WA3(I-2)*DI4+WA3(I-1)*DR4
            CH(I-1,K,5) = WA4(I-2)*DR5-WA4(I-1)*DI5
            CH(I,K,5) = WA4(I-2)*DI5+WA4(I-1)*DR5
  105    CONTINUE
  106 CONTINUE
      RETURN
      END
      SUBROUTINE RADFG(IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
C***BEGIN PROLOGUE  RADFG
C***REFER TO  RFFTF
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  RADFG
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,
     1                C1(IDO,L1,IP)          ,C2(IDL1,IP),
     2                CH2(IDL1,IP)           ,WA(*)
C***FIRST EXECUTABLE STATEMENT  RADFG
      TPI = 8.*ATAN(1.)
      ARG = TPI/REAL(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IPPH = (IP+1)/2
      IPP2 = IP+2
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IF (IDO .EQ. 1) GO TO 119
      DO 101 IK=1,IDL1
         CH2(IK,1) = C2(IK,1)
  101 CONTINUE
      DO 103 J=2,IP
         DO 102 K=1,L1
            CH(1,K,J) = C1(1,K,J)
  102    CONTINUE
  103 CONTINUE
      IF (NBD .GT. L1) GO TO 107
      IS = -IDO
      DO 106 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 105 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 104 K=1,L1
               CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J)+WA(IDIJ)*C1(I,K,J)
               CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J)-WA(IDIJ)*C1(I-1,K,J)
  104       CONTINUE
  105    CONTINUE
  106 CONTINUE
      GO TO 111
  107 IS = -IDO
      DO 110 J=2,IP
         IS = IS+IDO
         DO 109 K=1,L1
            IDIJ = IS
CDIR$ IVDEP
            DO 108 I=3,IDO,2
               IDIJ = IDIJ+2
               CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J)+WA(IDIJ)*C1(I,K,J)
               CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J)-WA(IDIJ)*C1(I-1,K,J)
  108       CONTINUE
  109    CONTINUE
  110 CONTINUE
  111 IF (NBD .LT. L1) GO TO 115
      DO 114 J=2,IPPH
         JC = IPP2-J
         DO 113 K=1,L1
CDIR$ IVDEP
            DO 112 I=3,IDO,2
               C1(I-1,K,J) = CH(I-1,K,J)+CH(I-1,K,JC)
               C1(I-1,K,JC) = CH(I,K,J)-CH(I,K,JC)
               C1(I,K,J) = CH(I,K,J)+CH(I,K,JC)
               C1(I,K,JC) = CH(I-1,K,JC)-CH(I-1,K,J)
  112       CONTINUE
  113    CONTINUE
  114 CONTINUE
      GO TO 121
  115 DO 118 J=2,IPPH
         JC = IPP2-J
         DO 117 I=3,IDO,2
            DO 116 K=1,L1
               C1(I-1,K,J) = CH(I-1,K,J)+CH(I-1,K,JC)
               C1(I-1,K,JC) = CH(I,K,J)-CH(I,K,JC)
               C1(I,K,J) = CH(I,K,J)+CH(I,K,JC)
               C1(I,K,JC) = CH(I-1,K,JC)-CH(I-1,K,J)
  116       CONTINUE
  117    CONTINUE
  118 CONTINUE
      GO TO 121
  119 DO 120 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  120 CONTINUE
  121 DO 123 J=2,IPPH
         JC = IPP2-J
         DO 122 K=1,L1
            C1(1,K,J) = CH(1,K,J)+CH(1,K,JC)
            C1(1,K,JC) = CH(1,K,JC)-CH(1,K,J)
  122    CONTINUE
  123 CONTINUE
C
      AR1 = 1.
      AI1 = 0.
      DO 127 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 124 IK=1,IDL1
            CH2(IK,L) = C2(IK,1)+AR1*C2(IK,2)
            CH2(IK,LC) = AI1*C2(IK,IP)
  124    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 126 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 125 IK=1,IDL1
               CH2(IK,L) = CH2(IK,L)+AR2*C2(IK,J)
               CH2(IK,LC) = CH2(IK,LC)+AI2*C2(IK,JC)
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      DO 129 J=2,IPPH
         DO 128 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+C2(IK,J)
  128    CONTINUE
  129 CONTINUE
C
      IF (IDO .LT. L1) GO TO 132
      DO 131 K=1,L1
         DO 130 I=1,IDO
            CC(I,1,K) = CH(I,K,1)
  130    CONTINUE
  131 CONTINUE
      GO TO 135
  132 DO 134 I=1,IDO
         DO 133 K=1,L1
            CC(I,1,K) = CH(I,K,1)
  133    CONTINUE
  134 CONTINUE
  135 DO 137 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 136 K=1,L1
            CC(IDO,J2-2,K) = CH(1,K,J)
            CC(1,J2-1,K) = CH(1,K,JC)
  136    CONTINUE
  137 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IF (NBD .LT. L1) GO TO 141
      DO 140 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 139 K=1,L1
CDIR$ IVDEP
            DO 138 I=3,IDO,2
               IC = IDP2-I
               CC(I-1,J2-1,K) = CH(I-1,K,J)+CH(I-1,K,JC)
               CC(IC-1,J2-2,K) = CH(I-1,K,J)-CH(I-1,K,JC)
               CC(I,J2-1,K) = CH(I,K,J)+CH(I,K,JC)
               CC(IC,J2-2,K) = CH(I,K,JC)-CH(I,K,J)
  138       CONTINUE
  139    CONTINUE
  140 CONTINUE
      RETURN
  141 DO 144 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 143 I=3,IDO,2
            IC = IDP2-I
            DO 142 K=1,L1
               CC(I-1,J2-1,K) = CH(I-1,K,J)+CH(I-1,K,JC)
               CC(IC-1,J2-2,K) = CH(I-1,K,J)-CH(I-1,K,JC)
               CC(I,J2-1,K) = CH(I,K,J)+CH(I,K,JC)
               CC(IC,J2-2,K) = CH(I,K,JC)-CH(I,K,J)
  142       CONTINUE
  143    CONTINUE
  144 CONTINUE
      RETURN
      END
      SUBROUTINE RADF2(IDO,L1,CC,CH,WA1)
C***BEGIN PROLOGUE  RADF2
C***REFER TO  RFFTF
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  RADF2
      DIMENSION       CH(IDO,2,L1)           ,CC(IDO,L1,2)           ,
     1                WA1(*)
C***FIRST EXECUTABLE STATEMENT  RADF2
      DO 101 K=1,L1
         CH(1,1,K) = CC(1,K,1)+CC(1,K,2)
         CH(IDO,2,K) = CC(1,K,1)-CC(1,K,2)
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      IF((IDO-1)/2.LT.L1) GO TO 108
      DO 104 K=1,L1
CDIR$ IVDEP
         DO 103 I=3,IDO,2
            IC = IDP2-I
            TR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            TI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            CH(I,1,K) = CC(I,K,1)+TI2
            CH(IC,2,K) = TI2-CC(I,K,1)
            CH(I-1,1,K) = CC(I-1,K,1)+TR2
            CH(IC-1,2,K) = CC(I-1,K,1)-TR2
  103    CONTINUE
  104 CONTINUE
      GO TO 111
  108 DO 110 I=3,IDO,2
         IC = IDP2-I
CDIR$ IVDEP
         DO 109 K=1,L1
            TR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            TI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            CH(I,1,K) = CC(I,K,1)+TI2
            CH(IC,2,K) = TI2-CC(I,K,1)
            CH(I-1,1,K) = CC(I-1,K,1)+TR2
            CH(IC-1,2,K) = CC(I-1,K,1)-TR2
  109    CONTINUE
  110 CONTINUE
  111 IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         CH(1,2,K) = -CC(IDO,K,2)
         CH(IDO,1,K) = CC(IDO,K,1)
  106 CONTINUE
  107 RETURN
      END
      SUBROUTINE RADF3(IDO,L1,CC,CH,WA1,WA2)
C***BEGIN PROLOGUE  RADF3
C***REFER TO  RFFTF
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  RADF3
      DIMENSION       CH(IDO,3,L1)           ,CC(IDO,L1,3)           ,
     1                WA1(*)     ,WA2(*)
C***FIRST EXECUTABLE STATEMENT  RADF3
      TAUR = -.5
      TAUI = .5*SQRT(3.)
      DO 101 K=1,L1
         CR2 = CC(1,K,2)+CC(1,K,3)
         CH(1,1,K) = CC(1,K,1)+CR2
         CH(1,3,K) = TAUI*(CC(1,K,3)-CC(1,K,2))
         CH(IDO,2,K) = CC(1,K,1)+TAUR*CR2
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      IF((IDO-1)/2.LT.L1) GO TO 104
      DO 103 K=1,L1
CDIR$ IVDEP
         DO 102 I=3,IDO,2
            IC = IDP2-I
            DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            CR2 = DR2+DR3
            CI2 = DI2+DI3
            CH(I-1,1,K) = CC(I-1,K,1)+CR2
            CH(I,1,K) = CC(I,K,1)+CI2
            TR2 = CC(I-1,K,1)+TAUR*CR2
            TI2 = CC(I,K,1)+TAUR*CI2
            TR3 = TAUI*(DI2-DI3)
            TI3 = TAUI*(DR3-DR2)
            CH(I-1,3,K) = TR2+TR3
            CH(IC-1,2,K) = TR2-TR3
            CH(I,3,K) = TI2+TI3
            CH(IC,2,K) = TI3-TI2
  102    CONTINUE
  103 CONTINUE
      RETURN
  104 DO 106 I=3,IDO,2
         IC = IDP2-I
CDIR$ IVDEP
         DO 105 K=1,L1
            DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            CR2 = DR2+DR3
            CI2 = DI2+DI3
            CH(I-1,1,K) = CC(I-1,K,1)+CR2
            CH(I,1,K) = CC(I,K,1)+CI2
            TR2 = CC(I-1,K,1)+TAUR*CR2
            TI2 = CC(I,K,1)+TAUR*CI2
            TR3 = TAUI*(DI2-DI3)
            TI3 = TAUI*(DR3-DR2)
            CH(I-1,3,K) = TR2+TR3
            CH(IC-1,2,K) = TR2-TR3
            CH(I,3,K) = TI2+TI3
            CH(IC,2,K) = TI3-TI2
  105    CONTINUE
  106 CONTINUE
      RETURN
      END
      SUBROUTINE RADF4(IDO,L1,CC,CH,WA1,WA2,WA3)
C***BEGIN PROLOGUE  RADF4
C***REFER TO  RFFTF
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  RADF4
      DIMENSION       CC(IDO,L1,4)           ,CH(IDO,4,L1)           ,
     1                WA1(*)     ,WA2(*)     ,WA3(*)
C***FIRST EXECUTABLE STATEMENT  RADF4
      HSQT2 = .5*SQRT(2.)
      DO 101 K=1,L1
         TR1 = CC(1,K,2)+CC(1,K,4)
         TR2 = CC(1,K,1)+CC(1,K,3)
         CH(1,1,K) = TR1+TR2
         CH(IDO,4,K) = TR2-TR1
         CH(IDO,2,K) = CC(1,K,1)-CC(1,K,3)
         CH(1,3,K) = CC(1,K,4)-CC(1,K,2)
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      IF((IDO-1)/2.LT.L1) GO TO 111
      DO 104 K=1,L1
CDIR$ IVDEP
         DO 103 I=3,IDO,2
            IC = IDP2-I
            CR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            CI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            CR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            CI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            CR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
            CI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
            TR1 = CR2+CR4
            TR4 = CR4-CR2
            TI1 = CI2+CI4
            TI4 = CI2-CI4
            TI2 = CC(I,K,1)+CI3
            TI3 = CC(I,K,1)-CI3
            TR2 = CC(I-1,K,1)+CR3
            TR3 = CC(I-1,K,1)-CR3
            CH(I-1,1,K) = TR1+TR2
            CH(IC-1,4,K) = TR2-TR1
            CH(I,1,K) = TI1+TI2
            CH(IC,4,K) = TI1-TI2
            CH(I-1,3,K) = TI4+TR3
            CH(IC-1,2,K) = TR3-TI4
            CH(I,3,K) = TR4+TI3
            CH(IC,2,K) = TR4-TI3
  103    CONTINUE
  104 CONTINUE
      GO TO 110
  111 DO 109 I=3,IDO,2
         IC = IDP2-I
CDIR$ IVDEP
         DO 108 K=1,L1
            CR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            CI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            CR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            CI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            CR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
            CI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
            TR1 = CR2+CR4
            TR4 = CR4-CR2
            TI1 = CI2+CI4
            TI4 = CI2-CI4
            TI2 = CC(I,K,1)+CI3
            TI3 = CC(I,K,1)-CI3
            TR2 = CC(I-1,K,1)+CR3
            TR3 = CC(I-1,K,1)-CR3
            CH(I-1,1,K) = TR1+TR2
            CH(IC-1,4,K) = TR2-TR1
            CH(I,1,K) = TI1+TI2
            CH(IC,4,K) = TI1-TI2
            CH(I-1,3,K) = TI4+TR3
            CH(IC-1,2,K) = TR3-TI4
            CH(I,3,K) = TR4+TI3
            CH(IC,2,K) = TR4-TI3
  108    CONTINUE
  109 CONTINUE
  110 IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         TI1 = -HSQT2*(CC(IDO,K,2)+CC(IDO,K,4))
         TR1 = HSQT2*(CC(IDO,K,2)-CC(IDO,K,4))
         CH(IDO,1,K) = TR1+CC(IDO,K,1)
         CH(IDO,3,K) = CC(IDO,K,1)-TR1
         CH(1,2,K) = TI1-CC(IDO,K,3)
         CH(1,4,K) = TI1+CC(IDO,K,3)
  106 CONTINUE
  107 RETURN
      END
      SUBROUTINE RADF5(IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
C***BEGIN PROLOGUE  RADF5
C***REFER TO  RFFTF
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  RADF5
      DIMENSION       CC(IDO,L1,5)           ,CH(IDO,5,L1)           ,
     1                WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)
C***FIRST EXECUTABLE STATEMENT  RADF5
      PI = 4.*ATAN(1.)
      TR11 = SIN(.1*PI)
      TI11 = SIN(.4*PI)
      TR12 = -SIN(.3*PI)
      TI12 = SIN(.2*PI)
      DO 101 K=1,L1
         CR2 = CC(1,K,5)+CC(1,K,2)
         CI5 = CC(1,K,5)-CC(1,K,2)
         CR3 = CC(1,K,4)+CC(1,K,3)
         CI4 = CC(1,K,4)-CC(1,K,3)
         CH(1,1,K) = CC(1,K,1)+CR2+CR3
         CH(IDO,2,K) = CC(1,K,1)+TR11*CR2+TR12*CR3
         CH(1,3,K) = TI11*CI5+TI12*CI4
         CH(IDO,4,K) = CC(1,K,1)+TR12*CR2+TR11*CR3
         CH(1,5,K) = TI12*CI5-TI11*CI4
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      IF((IDO-1)/2.LT.L1) GO TO 104
      DO 103 K=1,L1
CDIR$ IVDEP
         DO 102 I=3,IDO,2
            IC = IDP2-I
            DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            DR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
            DI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
            DR5 = WA4(I-2)*CC(I-1,K,5)+WA4(I-1)*CC(I,K,5)
            DI5 = WA4(I-2)*CC(I,K,5)-WA4(I-1)*CC(I-1,K,5)
            CR2 = DR2+DR5
            CI5 = DR5-DR2
            CR5 = DI2-DI5
            CI2 = DI2+DI5
            CR3 = DR3+DR4
            CI4 = DR4-DR3
            CR4 = DI3-DI4
            CI3 = DI3+DI4
            CH(I-1,1,K) = CC(I-1,K,1)+CR2+CR3
            CH(I,1,K) = CC(I,K,1)+CI2+CI3
            TR2 = CC(I-1,K,1)+TR11*CR2+TR12*CR3
            TI2 = CC(I,K,1)+TR11*CI2+TR12*CI3
            TR3 = CC(I-1,K,1)+TR12*CR2+TR11*CR3
            TI3 = CC(I,K,1)+TR12*CI2+TR11*CI3
            TR5 = TI11*CR5+TI12*CR4
            TI5 = TI11*CI5+TI12*CI4
            TR4 = TI12*CR5-TI11*CR4
            TI4 = TI12*CI5-TI11*CI4
            CH(I-1,3,K) = TR2+TR5
            CH(IC-1,2,K) = TR2-TR5
            CH(I,3,K) = TI2+TI5
            CH(IC,2,K) = TI5-TI2
            CH(I-1,5,K) = TR3+TR4
            CH(IC-1,4,K) = TR3-TR4
            CH(I,5,K) = TI3+TI4
            CH(IC,4,K) = TI4-TI3
  102    CONTINUE
  103 CONTINUE
      RETURN
  104 DO 106 I=3,IDO,2
         IC = IDP2-I
CDIR$ IVDEP
         DO 105 K=1,L1
            DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            DR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
            DI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
            DR5 = WA4(I-2)*CC(I-1,K,5)+WA4(I-1)*CC(I,K,5)
            DI5 = WA4(I-2)*CC(I,K,5)-WA4(I-1)*CC(I-1,K,5)
            CR2 = DR2+DR5
            CI5 = DR5-DR2
            CR5 = DI2-DI5
            CI2 = DI2+DI5
            CR3 = DR3+DR4
            CI4 = DR4-DR3
            CR4 = DI3-DI4
            CI3 = DI3+DI4
            CH(I-1,1,K) = CC(I-1,K,1)+CR2+CR3
            CH(I,1,K) = CC(I,K,1)+CI2+CI3
            TR2 = CC(I-1,K,1)+TR11*CR2+TR12*CR3
            TI2 = CC(I,K,1)+TR11*CI2+TR12*CI3
            TR3 = CC(I-1,K,1)+TR12*CR2+TR11*CR3
            TI3 = CC(I,K,1)+TR12*CI2+TR11*CI3
            TR5 = TI11*CR5+TI12*CR4
            TI5 = TI11*CI5+TI12*CI4
            TR4 = TI12*CR5-TI11*CR4
            TI4 = TI12*CI5-TI11*CI4
            CH(I-1,3,K) = TR2+TR5
            CH(IC-1,2,K) = TR2-TR5
            CH(I,3,K) = TI2+TI5
            CH(IC,2,K) = TI5-TI2
            CH(I-1,5,K) = TR3+TR4
            CH(IC-1,4,K) = TR3-TR4
            CH(I,5,K) = TI3+TI4
            CH(IC,4,K) = TI4-TI3
  105    CONTINUE
  106 CONTINUE
      RETURN
      END
      SUBROUTINE RFFTB1(N,C,CH,WA,FAC)
C***BEGIN PROLOGUE  RFFTB1
C***REFER TO  RFFTB
C***ROUTINES CALLED  RADB2,RADB3,RADB4,RADB5,RADBG
C***END PROLOGUE  RFFTB1
      DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,FAC(*)
C***FIRST EXECUTABLE STATEMENT  RFFTB1
      NF = int(FAC(2))
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = int(FAC(K1+2))
         L2 = IP*L1
         IDO = N/L2
         IDL1 = IDO*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL RADB4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL RADB4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL RADB2 (IDO,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL RADB2 (IDO,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 107
         CALL RADB3 (IDO,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL RADB3 (IDO,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 110
         CALL RADB5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL RADB5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL RADBG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL RADBG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (IDO .EQ. 1) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDO
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      DO 117 I=1,N
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END
      SUBROUTINE RFFTF1(N,C,CH,WA,FAC)
C***BEGIN PROLOGUE  RFFTF1
C***REFER TO  RFFTF
C***ROUTINES CALLED  RADF2,RADF3,RADF4,RADF5,RADFG
C***END PROLOGUE  RFFTF1
      DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,FAC(*)
C***FIRST EXECUTABLE STATEMENT  RFFTF1
      NF = int(fac(2))
      NA = 1
      L2 = N
      IW = N
      DO 111 K1=1,NF
         KH = NF-K1
         IP = int(FAC(KH+3))
         L1 = L2/IP
         IDO = N/L2
         IDL1 = IDO*L1
         IW = IW-(IP-1)*IDO
         NA = 1-NA
         IF (IP .NE. 4) GO TO 102
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL RADF4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  101    CALL RADF4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  102    IF (IP .NE. 2) GO TO 104
         IF (NA .NE. 0) GO TO 103
         CALL RADF2 (IDO,L1,C,CH,WA(IW))
         GO TO 110
  103    CALL RADF2 (IDO,L1,CH,C,WA(IW))
         GO TO 110
  104    IF (IP .NE. 3) GO TO 106
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 105
         CALL RADF3 (IDO,L1,C,CH,WA(IW),WA(IX2))
         GO TO 110
  105    CALL RADF3 (IDO,L1,CH,C,WA(IW),WA(IX2))
         GO TO 110
  106    IF (IP .NE. 5) GO TO 108
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 107
         CALL RADF5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  107    CALL RADF5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  108    IF (IDO .EQ. 1) NA = 1-NA
         IF (NA .NE. 0) GO TO 109
         CALL RADFG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         NA = 1
         GO TO 110
  109    CALL RADFG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
         NA = 0
  110    L2 = L1
  111 CONTINUE
      IF (NA .EQ. 1) RETURN
      DO 112 I=1,N
         C(I) = CH(I)
  112 CONTINUE
      RETURN
      END
      SUBROUTINE RFFTI1(N,WA,FAC)
C***BEGIN PROLOGUE  RFFTI1
C***REFER TO  RFFTI
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  RFFTI1
      DIMENSION       WA(*)      ,FAC(*)    ,NTRYH(4)
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
C***FIRST EXECUTABLE STATEMENT  RFFTI1
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      FAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         FAC(IB+2) = FAC(IB+1)
  106 CONTINUE
      FAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      FAC(1) = N
      FAC(2) = NF
      TPI = 8.*ATAN(1.)
      ARGH = TPI/REAL(N)
      IS = 0
      NFM1 = NF-1
      L1 = 1
      IF (NFM1 .EQ. 0) RETURN
      DO 110 K1=1,NFM1
         IP = int(FAC(K1+2))
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IPM = IP-1
         DO 109 J=1,IPM
            LD = LD+L1
            I = IS
            ARGLD = REAL(LD)*ARGH
            FI = 0.
            DO 108 II=3,IDO,2
               I = I+2
               FI = FI+1.
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IS = IS+IDO
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END
      REAL FUNCTION R1MACH(I)
C***BEGIN PROLOGUE  R1MACH
C***DATE WRITTEN   790101   (YYMMDD)
C***REVISION DATE  831014   (YYMMDD)
C***CATEGORY NO.  R1
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  FOX, P. A., (BELL LABS)
C           HALL, A. D., (BELL LABS)
C           SCHRYER, N. L., (BELL LABS)
C***PURPOSE  RETURNS SINGLE PRECISION MACHINE DEPENDENT CONSTANTS
C***DESCRIPTION
C
C     R1MACH CAN BE USED TO OBTAIN MACHINE-DEPENDENT PARAMETERS
C     FOR THE LOCAL MACHINE ENVIRONMENT.  IT IS A FUNCTION
C     SUBROUTINE WITH ONE (INPUT) ARGUMENT, AND CAN BE CALLED
C     AS FOLLOWS, FOR EXAMPLE
C
C          A = R1MACH(I)
C
C     WHERE I=1,...,5.  THE (OUTPUT) VALUE OF A ABOVE IS
C     DETERMINED BY THE (INPUT) VALUE OF I.  THE RESULTS FOR
C     VARIOUS VALUES OF I ARE DISCUSSED BELOW.
C
C  SINGLE-PRECISION MACHINE CONSTANTS
C  R1MACH(1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C  R1MACH(2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C  R1MACH(3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C  R1MACH(4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C  R1MACH(5) = LOG10(B)
C***REFERENCES  FOX, P.A., HALL, A.D., SCHRYER, N.L, *FRAMEWORK FOR
C                 A PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHE-
C                 MATICAL SOFTWARE, VOL. 4, NO. 2, JUNE 1978,
C                 PP. 177-188.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  R1MACH
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
C
      REAL RMACH(5)
C
      EQUIVALENCE (RMACH(1),SMALL(1))
      EQUIVALENCE (RMACH(2),LARGE(1))
      EQUIVALENCE (RMACH(3),RIGHT(1))
      EQUIVALENCE (RMACH(4),DIVER(1))
      EQUIVALENCE (RMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C     DATA RMACH(1) / Z400800000 /
C     DATA RMACH(2) / Z5FFFFFFFF /
C     DATA RMACH(3) / Z4E9800000 /
C     DATA RMACH(4) / Z4EA800000 /
C     DATA RMACH(5) / Z500E730E8 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS.
C
C     DATA RMACH(1) / O1771000000000000 /
C     DATA RMACH(2) / O0777777777777777 /
C     DATA RMACH(3) / O1311000000000000 /
C     DATA RMACH(4) / O1301000000000000 /
C     DATA RMACH(5) / O1157163034761675 /
C
C
C     MACHINE CONSTANTS FOR THE CDC CYBER 170 SERIES (FTN5).
C
      DATA RMACH(1) / O"00014000000000000000" /
      DATA RMACH(2) / O"37767777777777777777" /
      DATA RMACH(3) / O"16404000000000000000" /
      DATA RMACH(4) / O"16414000000000000000" /
      DATA RMACH(5) / O"17164642023241175720" /
C
C     MACHINE CONSTANTS FOR THE CDC CYBER 200 SERIES
C
C     DATA RMACH(1) / X'9000400000000000' /
C     DATA RMACH(2) / X'6FFF7FFFFFFFFFFF' /
C     DATA RMACH(3) / X'FFA3400000000000' /
C     DATA RMACH(4) / X'FFA4400000000000' /
C     DATA RMACH(5) / X'FFD04D104D427DE8' /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
C
C     DATA RMACH(1) / 00564000000000000000B /
C     DATA RMACH(2) / 37767777777777777776B /
C     DATA RMACH(3) / 16414000000000000000B /
C     DATA RMACH(4) / 16424000000000000000B /
C     DATA RMACH(5) / 17164642023241175720B /
C
C     MACHINE CONSTANTS FOR THE CRAY 1
C
C     DATA RMACH(1) / 200034000000000000000B /
C     DATA RMACH(2) / 577767777777777777776B /
C     DATA RMACH(3) / 377224000000000000000B /
C     DATA RMACH(4) / 377234000000000000000B /
C     DATA RMACH(5) / 377774642023241175720B /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
C     STATIC RMACH(5)
C
C     DATA SMALL/20K,0/,LARGE/77777K,177777K/
C     DATA RIGHT/35420K,0/,DIVER/36020K,0/
C     DATA LOG10/40423K,42023K/
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA SMALL(1),SMALL(2) / '20000000, '00000201 /
C     DATA LARGE(1),LARGE(2) / '37777777, '00000177 /
C     DATA RIGHT(1),RIGHT(2) / '20000000, '00000352 /
C     DATA DIVER(1),DIVER(2) / '20000000, '00000353 /
C     DATA LOG10(1),LOG10(2) / '23210115, '00000377 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.
C
C     DATA RMACH(1) / O402400000000 /
C     DATA RMACH(2) / O376777777777 /
C     DATA RMACH(3) / O714400000000 /
C     DATA RMACH(4) / O716400000000 /
C     DATA RMACH(5) / O776464202324 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C
C     3 WORD DOUBLE PRECISION WITH FTN4
C
C     DATA SMALL(1), SMALL(2) / 40000B,       1 /
C     DATA LARGE(1), LARGE(2) / 77777B, 177776B /
C     DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
C     DATA DIVER(1), DIVER(2) / 40000B,    327B /
C     DATA LOG10(1), LOG10(2) / 46420B,  46777B /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     4 WORD DOUBLE PRECISION WITH FTN4
C
C     DATA SMALL(1), SMALL(2) / 40000B,       1 /
C     DATA LARGE91), LARGE(2) / 77777B, 177776B /
C     DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
C     DATA DIVER(1), DIVER(2) / 40000B,    327B /
C     DATA LOG10(1), LOG10(2) / 46420B,  46777B /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86  AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA RMACH(1) / Z00100000 /
C     DATA RMACH(2) / Z7FFFFFFF /
C     DATA RMACH(3) / Z3B100000 /
C     DATA RMACH(4) / Z3C100000 /
C     DATA RMACH(5) / Z41134413 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR).
C
C     DATA RMACH(1) / "000400000000 /
C     DATA RMACH(2) / "377777777777 /
C     DATA RMACH(3) / "146400000000 /
C     DATA RMACH(4) / "147400000000 /
C     DATA RMACH(5) / "177464202324 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1) /    8388608 /
C     DATA LARGE(1) / 2147483647 /
C     DATA RIGHT(1) /  880803840 /
C     DATA DIVER(1) /  889192448 /
C     DATA LOG10(1) / 1067065499 /
C
C     DATA RMACH(1) / O00040000000 /
C     DATA RMACH(2) / O17777777777 /
C     DATA RMACH(3) / O06440000000 /
C     DATA RMACH(4) / O06500000000 /
C     DATA RMACH(5) / O07746420233 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGERS  (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1),SMALL(2) /   128,     0 /
C     DATA LARGE(1),LARGE(2) / 32767,    -1 /
C     DATA RIGHT(1),RIGHT(2) / 13440,     0 /
C     DATA DIVER(1),DIVER(2) / 13568,     0 /
C     DATA LOG10(1),LOG10(2) / 16282,  8347 /
C
C     DATA SMALL(1),SMALL(2) / O000200, O000000 /
C     DATA LARGE(1),LARGE(2) / O077777, O177777 /
C     DATA RIGHT(1),RIGHT(2) / O032200, O000000 /
C     DATA DIVER(1),DIVER(2) / O032400, O000000 /
C     DATA LOG10(1),LOG10(2) / O037632, O020233 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C     DATA RMACH(1) / O000400000000 /
C     DATA RMACH(2) / O377777777777 /
C     DATA RMACH(3) / O146400000000 /
C     DATA RMACH(4) / O147400000000 /
C     DATA RMACH(5) / O177464202324 /
C
C     MACHINE CONSTANTS FOR THE VAX 11/780
C    (EXPRESSED IN INTEGER AND HEXADECIMAL)
C  ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS***
C  *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS***
C
C     DATA SMALL(1) /       128 /
C     DATA LARGE(1) /    -32769 /
C     DATA RIGHT(1) /     13440 /
C     DATA DIVER(1) /     13568 /
C     DATA LOG10(1) / 547045274 /
C
C     DATA SMALL(1) / Z00000080 /
C     DATA LARGE(1) / ZFFFF7FFF /
C     DATA RIGHT(1) / Z00003480 /
C     DATA DIVER(1) / Z00003500 /
C     DATA LOG10(1) / Z209B3F9A /
C
C     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
C
C     DATA SMALL(1),SMALL(2) /     0,    256/
C     DATA LARGE(1),LARGE(2) /    -1,   -129/
C     DATA RIGHT(1),RIGHT(2) /     0,  26880/
C     DATA DIVER(1),DIVER(2) /     0,  27136/
C     DATA LOG10(1),LOG10(2) /  8347,  32538/
C
C
C***FIRST EXECUTABLE STATEMENT  R1MACH
C
      R1MACH = RMACH(I)
      RETURN
C
      END
C
C  EXAMPLE OF HFFT3 USAGE
C
C-----------------------------------------------------------------------
C
C  PROBLEM   UXX + UYY + UZZ - PI**2*U = G(X,Y)
C
C            U = 1                    ON X=0
C            UX= Z*EXP(Z)+COS(2*PI*Y) ON X=1
C            UZ= X                    ON Z=0
C            U = EXP(X)+X*COS(2*PI*Y) ON Z=1
C
C            U(X,0,Z) = U(X,1,Z)
C
C            G = (X**2+Z**2-PI**2)*EXP(X*Z) - 5*PI**2*X*COS(2*PI*Y)
C
C  SOLUTION  U = EXP(X*Z) + X*COS(2*PI*Y)
C
C-----------------------------------------------------------------------
C
C  ------------
C  DECLARATIONS
C  ------------
C
C     ... CONSTANTS
C
C      PARAMETER (NMAX=29)
C      PARAMETER (NWORK = 3+NMAX*(7+NMAX*(11+NMAX)))
C      PARAMETER (LDXU = NMAX+2)
C
C     ... GLOBAL VARIABLES
C
C      COMMON /GLOBAL/ TWOPI, PISQR
C      EXTERNAL PRHS, BRHS
C
C     ... LOCAL VARIABLES
C
C      INTEGER BCTY(6)
C      REAL U(LDXU,LDXU,LDXU), WORK(NWORK)
C
C---------------------------------------------------------------------
C
C      PI = 4.0*ATAN(1.0)
C      TWOPI = 2.0*PI
C      PISQR = PI*PI
C
C  -------------
C  SETUP PROBLEM
C  -------------
C
C      AX = 0.0
C      BX = 1.0
C      AY = 0.0
C      BY = 1.0
C      AZ = 0.0
C      BZ = 1.0
C      COEFU = -PISQR
C      BCTY(1) = 2
C      BCTY(2) = 3
C      BCTY(3) = 1
C      BCTY(4) = 3
C      BCTY(5) = 1
C      BCTY(6) = 2
C      IORDER = 4
C
C  --------------------------
C  SOLVE ON SEQUENCE OF GRIDS
C  --------------------------
C
C      PRINT 2000
C      DO 500 N=5,29,4
C         NX = N
C         NY = N
C         NZ = N
C
C        ... SOLVE PDE
C
C         CALL HFFT3(COEFU,PRHS,BRHS,AX,BX,AY,BY,AZ,BZ,NX,NY,NZ,
C     *              BCTY,IORDER,U,LDXU,LDXU,WORK,NWORK,INFO)
C         IF (INFO .LT. 0)  GO TO 900
C
C        ... EVALUATE ERROR
C
C         H = (BX-AX)/REAL(NX-1)
C         ERRMAX = 0.0E0
C         DO 100 K=1,NZ
C            Z = AZ + REAL(K-1)*H
C            DO 100 J=1,NY
C               Y = AY + REAL(J-1)*H
C               DO 100 I=1,NX
C                  X = AX + REAL(I-1)*H
C                  TRUSOL = TRUE(X,Y,Z)
C                  ERROR = ABS(TRUSOL-U(I,J,K))
C                  ERRMAX = MAX(ERROR,ERRMAX)
C  100    CONTINUE
C         PRINT 2001,N,ERRMAX
C  500 CONTINUE
C      STOP
CC
CC     ... ERROR EXIT
CC
C  900 CONTINUE
C      PRINT 2002,INFO
C      STOP
C
C 2000 FORMAT(/'   GRID    MAX-ERROR' / '  -------------------' /)
C 2001 FORMAT(4X,I2,4X,1PE10.3)
C 2002 FORMAT(/'   HFFT3 RETURNED INFO = ',I2)
CC
C      END
C      FUNCTION PRHS (X,Y,Z)
C
C     ... RIGHT HAND SIDE OF PDE (USER-SUPPLIED)
C
C      COMMON /GLOBAL/ TWOPI, PISQR
C      PRHS = (X*X+Z*Z-PISQR)*EXP(X*Z) - 5.0*PISQR*X*COS(TWOPI*Y)
C      RETURN
C      END
C      FUNCTION BRHS(K,X,Y,Z)
C
C     ... RIGHT HAND SIDE OF BOUNDARY CONDITIONS (USER-SUPPLIED)
C
C      COMMON /GLOBAL/ TWOPI, PISQR
C      GO TO (1,2,3,4,5,6),K
C      GO TO 999
C
C    1 CONTINUE
C      BRHS = Z*EXP(Z) + COS(TWOPI*Y)
C      GO TO 999
C
C    2 CONTINUE
C      GO TO 999
C
C    3 CONTINUE
C      BRHS = 1.0
C      GO TO 999
C
C    4 CONTINUE
C      GO TO 999
C
C    5 CONTINUE
C      BRHS = EXP(X) + X*COS(TWOPI*Y)
C      GO TO 999
C
C    6 CONTINUE
C      BRHS = X
C      GO TO 999
C
C  999 CONTINUE
C      RETURN
C      END
C      FUNCTION TRUE (X,Y,Z)
C      COMMON /GLOBAL/ TWOPI, PISQR
C      TRUE = EXP(X*Z) + X*COS(TWOPI*Y)
C      RETURN
C      END
C     PROGRAM TESTH3 (INPUT,OUTPUT,TAPE5=INPUT,TAPE6=OUTPUT)
C
C                       ----------------------
C                       TEST DRIVER FOR HFFT3
C                       ----------------------
C
C   SOLVE SAMPLE PROBLEMS ON A SEQUENCE OF GRIDS AND PRINT STATISTICS.
C   USER SELECTS ORDER OF ACCURACY (2 OR 4) AND PROBLEM NUMBER (1-11).
C
C   NOTE -- THIS PROGRAM CALLS THE SUBROUTINE TIMER (INCLUDED) TO OBTAIN
C           THE ELAPSED CPU TIME IN SECONDS SINCE THE START OF THE RUN.
C           ROUTINE MUST BE REPLACED WHEN IMPLEMENTING THIS PROGRAM ON
C           A NEW COMPUTER.
C
C   THE PROBLEMS ARE
C
C     1.  SMOOTH HOMOGENEOUS DIRICHLET PROBLEM FOR THE POISSON EQUATION.
C         (PROBLEM A3 OF REFERENCE)
C
C         COEFU = 0
C         U GIVEN ON AX=0, BX=1, AY=0, BY=1, AZ=0, BZ=1
C
C         U = XYZ(1-X)(1-Y)(1-Z) EXP(X+Y+Z)
C
C     2.  SMOOTH HELMHOLTZ PROBLEM WITH SOME DERIVATIVE BOUNDARY
C         CONDITIONS; SOLUTION SAME AS IN PROBLEM 1.
C         (PROBLEM B3 OF REFERENCE)
C
C         COEFU = -5
C         U GIVEN ON BY=1, BZ=1
C         UX GIVEN ON AX=0, BX=1
C         UY GIVEN ON AY=0
C         UZ GIVEN ON AZ=0
C
C         U = XYZ(1-X)(1-Y)(1-Z) EXP(X+Y+Z)
C
C     3.  POISSON EQUATION; SOLUTION IS A WAVE FRONT ALONG A RIGHT
C         ANGLE JOINING TWO REGIONS WHERE IT IS A CONSTANT; SOLUTION
C         HAS DISCONTINUOUS 3RD DERIVATIVES.
C
C         COEFU = 0
C         U GIVEN ON AX=0, AY=0, AZ=0
C         UX GIVEN ON BX=1
C         UY GIVEN ON BY=1
C         UZ GIVEN ON BZ=1
C
C         U = P(X)*P(Y)*P(Z)
C
C         P(X) IS 1 FOR X.LT.0.15 AND 0 FOR X.GT.0.85.  BETWEEN THESE
C         IT IS DEFINED AS THE QUINTIC POLYNOMIAL WHICH JOINS THESE
C         TWO FLAT REGIONS AND GIVES P TWO CONTINUOUS DERIVATIVES.
C
C     4.  SMOOTH HELMHOLTZ PROBLEM WITH SOME DERIVATIVE BOUNDARY
C         CONDITIONS.
C
C         COEFU = -100
C         U GIVEN ON BX=1, BY=1, BZ=1
C         UX GIVEN ON AX=0
C         UY GIVEN ON AY=0
C         UZ GIVEN ON AZ=0
C
C         U = (P(X;10)+P(Y;20)+P(Z;5))/3
C
C         P(X;A) = COSH(A*X)/COSH(A)
C
C     5.  DIRICHLET PROBLEM FOR POISSON EQUATION; SOLUTION HAS
C         SINGULAR SECOND DERIVATIVES ALONG X=0, Y=0, Z=0.
C         (PROBLEM F3 OF THE REFERENCE)
C
C         COEFU = 0
C         U GIVEN ON AX=0, BX=1, AY=0, AY=1, AZ=0, BZ=1
C
C         U = (X**1.5-X)*(Y**1.5-Y)*(Z**1.5-Z)
C
C     6.  SMOOTH, PERIODIC HELMHOLTZ PROBLEM.
C         (PROBLEM C3 OF THE REFERENCE)
C
C         COEFU = -20
C         U PERIODIC ON AX=0, BX=PI, AY=0, BY=PI, AZ=0, BZ=PI
C
C         U = COS(4Y) + SIN(4(X-Y)) + COS(4Z)
C
C     7.  SMOOTH POISSON PROBLEM WITH SOME DERIVATIVE BOUNDARY
C         CONDITIONS.
C
C         COEFU = 0
C         U GIVEN ON AX=0, BX=1, AZ=-1, BZ=1
C         UY GIVEN ON AY=-1, BY=1
C
C         U = T(Y)*(A(X) + T(Y)*B(X))/2 + S(Z)*(A(X) + S(Z)*B(X))/2
C
C         T(Y) = 1-Y**2,  S(Z) = 1-Z,  A(X) = R*C(X) + EXP(SQRT(Q*X)),
C         B(X) = (7-P)*R/16/C(X),  P = 14+SQRT(133),  Q = 14-SQRT(133),
C         R = (7-Q)/R/SQRT(133),  C(X) = EXP(SQRT(P*X))-EXP(SQRT(Q*X))
C
C     8.  POISSON PROBLEM WITH SOME DERIVATIVE BOUNDARY CONDITIONS;
C         SOLUTION HAS SINGULAR THIRD DERIVATIVES ALONG X=0, Y=0, Z=0.
C         (PROBLEM E2 OF THE REFERNCE)
C
C         COEFU = 0
C         U GIVEN ON BX=1, BY=1
C         UX GIVEN ON AX=0
C         UY GIVEN ON AY=0
C
C         U = (X*Y)**2.5
C
C     9.  SMOOTH HELMHOLTZ PROBLEM WITH DIRICHLET, NEUMANN AND PERIODIC
C         BOUNDARY CONDITIONS; SAME SOLUTION AS PROBLEM 6.
C         (PROBLEM D2 OF THE REFERENCE)
C
C         COEFU = -20
C         U PERIODIC ON AX=0, BX=PI
C         U GIVEN ON AY=0, AZ=0
C         UY GIVEN ON BY=PI
C         UZ GIVEN ON BZ=PI
C
C         U = COS(4Y) + SIN(4(X-Y))
C
C    10.  HELMHOLTZ EQUATION WITH SOME DERIVATIVE BOUNDARY CONDITIONS;
C         SOLUTION IS CUBIC MONOMIAL; ALL METHODS SHOULD BE EXACT.
C
C         COEFU = -5
C         U GIVEN ON AY=0
C         UX GIVEN ON AX=0, BX=1
C         UY GIVEN ON BY=1
C         UZ GIVEN ON AZ=0, BZ=1
C
C         U = X*Y*Z
C
C    11.  HOMOGENEOUS NEUMANNT PROBLEM; INCLUDED TO COLLECT TIMING
C         DATA WHEN COST OF FUNCTION EVALUATIONS IS MINIMUM.
C
C         COEFU = 0
C         UX=0 ON AX=0, BX=1
C         UY=0 ON AY=0, BY=1
C         UZ=0 ON AZ=0, AZ=1
C
C
C-----------------------------------------------------------------------
C
C
C     ... CONSTANTS
C
C         NMAX = MAX GRID SIZE FOR THIS PROGRAM
C         NWORK = WORKSPACE REQUIRED FOR NMAX BY NMAX GRID
C         LDXU,LDYU,LDZU = REQUIRED DIMENSIONS FOR U
C         LUIN = INPUT UNIT FOR READS
C         LUOUT = OUTPUT UNIT FOR WRITES
C
C      PARAMETER (NMAX=29)
C      PARAMETER (NWORK = 33 + NMAX*(21 + NMAX*(12 + NMAX)))
C      PARAMETER (LDXU=NMAX+2, LDYU=NMAX+2, LDZU=NMAX+2)
C      PARAMETER (LUIN=5, LUOUT=6)
C
C     ... VARIABLES
C
C      INTEGER BCTY(6)
C      REAL
C     *     COEFU, AX, BX, AY, BY, AZ, BZ, U(LDXU,LDYU,LDZU),
C     *     WORK(NWORK), H, ABSERR, RELERR, TRUMAX, HOLD, EOLD, RATE
C
C      EXTERNAL PRHS, BRHS
C
C
C-----------------------------------------------------------------------
C
C
C  SELECT ORDER OF ACCURACY
C
C    5 CONTINUE
C      WRITE(LUOUT,*) 'ENTER ORDER OF ACCURACY (2 OR 4) '
C      READ(LUIN,*) IORDER
C      IF ((IORDER .NE. 2) .AND. (IORDER .NE. 4))  GO TO 5
C
C  SELECT PROBLEM
C
C   10 CONTINUE
C      WRITE(LUOUT,*) ' ENTER PROBLEM NUMBER (1-11, 0 TO STOP) '
C      READ(LUIN,*) KPROB
C      IF (KPROB .EQ. 0)  STOP
C      IF ((KPROB .LT. 0) .OR. (KPROB .GT. 11))  GO TO 10
C
C  SOLVE ON SEQUENCE OF GRIDS
C
C      WRITE(LUOUT,2000) KPROB,IORDER
C      EOLD = 0.0E0
C      NRUNS = (NMAX-1)/4
C      NX = 1
C      DO 500 K=1,NRUNS
C
C        ... SETUP PROBLEM
C
C         NX = NX + 4
C         CALL SETUP(KPROB,NX,NY,NZ,COEFU,AX,BX,AY,BY,AZ,BZ,BCTY)
C         IF (NY .GT. NMAX)  GO TO 10
C         IF (NZ .GT. NMAX)  GO TO 10
C
C        ... SOLVE PROBLEM
C
C         CALL TIMER(T0)
C         CALL HFFT3(COEFU,PRHS,BRHS,AX,BX,AY,BY,AZ,BZ,NX,NY,NZ,BCTY,
CC     *              IORDER,U,LDXU,LDYU,WORK,NWORK,INFO)
C         CALL TIMER(T1)
C         CPTIME = T1 - T0
C
C         IF (INFO .LT. 0) THEN
C            WRITE(6,2002) INFO
C            GO TO 10
C         ENDIF
C
C        ... COMPUTE ERROR
C
C         H = (BX-AX)/REAL(NX-1)
C         CALL GETERR(AX,AY,AZ,H,NX,NY,NZ,U,LDXU,LDYU,ABSERR,RELERR,
C     *               TRUMAX)
C
C        ... COMPUTE CONVERGENCE RATE
C
C         IF (RELERR .LT. EOLD) THEN
C            RATE = ABS(ALOG(RELERR/EOLD)/ALOG(H/HOLD))
C         ELSE
C            RATE = 0.0E0
C         ENDIF
C
C        ... PRINT SUMMARY
C
C         WRITE(LUOUT,2001) NX,NY,NZ,H,INFO,ABSERR,RELERR,TRUMAX,RATE,
C     *                     CPTIME
C         HOLD = H
C         EOLD = RELERR
C  500 CONTINUE
C      GO TO 10
C
C 2000 FORMAT(///'  PROBLEM ',I2,7X,'MODULE HFFT3',7X,'IORDER = ',I1
C     *       //'   NX  NY  NZ',7X,'H',5X,'INFO',3X,'ABSERR',6X,'RELERR',
C     *       6X,'TRUMAX',4X,'RATE',3X,'SECS' / 2X,39('--') /)
C 2001 FORMAT(2X,I3,1X,I3,1X,I3,2X,1P,E11.4,1X,I2,3E12.4,
C     *       2X,0P,F4.1,1X,F6.2)
C 2002 FORMAT(/' HFFT3 RETURNS INFO = ',I3 /)
C
C      END
      SUBROUTINE GETERR (AX,AY,AZ,H,NX,NY,NZ,U,LDXU,LDYU,ABSERR,
     *                   RELERR,TRUMAX)
C
C  ----------------------------------------------------------
C  COMPUTE MAX(ABSOLUTE ERROR) AND MAX(TRUE SOLUTION) ON GRID
C  ----------------------------------------------------------
C
      REAL AX,AY,AZ,H,U(LDXU,LDYU,*),ABSERR,RELERR,TRUMAX
      REAL X,Y,Z,TRUSOL,DIFF
C
      ABSERR = 0.0E0
      TRUMAX = 0.0E0
      DO 100 K=1,NZ
         Z = AZ + REAL(K-1)*H
         DO 100 J=1,NY
            Y = AY + H*REAL(J-1)
            DO 100 I=1,NX
               X = AX + H*REAL(I-1)
               TRUSOL = TRUE(X,Y,Z)
               ABSERR = MAX(ABS(TRUSOL-U(I,J,K)),ABSERR)
               TRUMAX = MAX(ABS(TRUSOL),TRUMAX)
  100 CONTINUE
C
      IF (TRUMAX .GT. 0.0E0) THEN
         RELERR = ABSERR/TRUMAX
      ELSE
         RELERR = ABSERR
      ENDIF
      RETURN
      END
      SUBROUTINE SETUP (KPROB,NX,NY,NZ,COEFU,AX,BX,AY,BY,AZ,BZ,BCTY)
C
C  -----------------------------------------------------------
C  SETUP PROBLEM KPROB (NX,NY,NZ,COEFU,AX,BX,AY,BY,AZ,BZ,BCTY)
C  -----------------------------------------------------------
C
      INTEGER KPROB,NX,NY,NZ,BCTY(6)
      REAL
     *     COEFU,AX,BX,AY,BY,AZ,BZ
C
      COMMON /SELECT/ IPROB
C
      IPROB = KPROB
C
C  DEFAULT VALUES
C
      NY = NX
      NZ = NX
      AX = 0.0E0
      BX = 1.0E0
      AY = 0.0E0
      BY = 1.0E0
      AZ = 0.0E0
      BZ = 1.0E0
      COEFU = 0.0E0
      BCTY(1) = 1
      BCTY(2) = 1
      BCTY(3) = 1
      BCTY(4) = 1
      BCTY(6) = 1
      BCTY(5) = 1
C
C  SELECT OPTIONS FOR EACH CASE
C
      GO TO (100,200,300,400,500,600,700,800,900,1000,1100), IPROB
      GO TO 9999
C
  100 CONTINUE
      GO TO 9999
C
  200 CONTINUE
      BCTY(1) = 2
      BCTY(2) = 2
      BCTY(3) = 2
      BCTY(6) = 2
      COEFU = -5.0E0
      GO TO 9999
C
  300 CONTINUE
      BCTY(1) = 2
      BCTY(4) = 2
      BCTY(5) = 2
      GO TO 9999
C
  400 CONTINUE
      BCTY(2) = 2
      BCTY(3) = 2
      BCTY(6) = 2
      COEFU = -100.0E0
      GO TO 9999
C
  500 CONTINUE
      GO TO 9999
C
  600 CONTINUE
      BX = 4.0E0*ATAN(1.0E0)
      BY = BX
      BZ = BX
      COEFU = -20.0E0
      DO 605 I=1,6
         BCTY(I) = 3
  605 CONTINUE
      GO TO 9999
C
  700 CONTINUE
      NY = 2*NX - 1
      NZ = NY
      AY = -1.0E0
      AZ = -1.0E0
      BCTY(2) = 2
      BCTY(4) = 2
      GO TO 9999
C
  800 CONTINUE
      BCTY(2) = 2
      BCTY(3) = 2
      GO TO 9999
C
  900 CONTINUE
      BX = 4.0E0*ATAN(1.0E0)
      BY = BX
      BZ = BX
      COEFU = -20.0E0
      BCTY(1) = 3
      BCTY(2) = 1
      BCTY(3) = 3
      BCTY(4) = 2
      BCTY(6) = 1
      BCTY(5) = 2
      GO TO 9999
C
 1000 CONTINUE
      BCTY(1) = 2
      BCTY(3) = 2
      BCTY(4) = 2
      BCTY(6) = 2
      BCTY(5) = 2
      COEFU = -5.0E0
      GO TO 9999
C
 1100 CONTINUE
      BCTY(1) = 2
      BCTY(2) = 2
      BCTY(3) = 2
      BCTY(4) = 2
      BCTY(6) = 2
      BCTY(5) = 2
      GO TO 9999
C
 9999 CONTINUE
      RETURN
      END
      REAL FUNCTION PRHS(X,Y,Z)
C
C  ----------------------------------------
C  RIGHT-HAND SIDE OF DIFFERENTIAL EQUATION
C  ----------------------------------------
C
      REAL X,Y,Z
C
      COMMON /SELECT/ IPROB
C
      GO TO (100,200,300,400,500,600,700,800,600,1000,1100), IPROB
C
         PRHS = 0.0E0
         GO TO 9999
C
  100 CONTINUE
      XMX2 = X*(1.0E0 - X)
      YMY2 = Y*(1.0E0 - Y)
      ZMZ2 = Z*(1.0E0 - Z)
      PRHS = -EXP(X+Y+Z)*( X*(X+3.0E0)*YMY2*ZMZ2 +
     *                     XMX2*Y*(Y+3.0E0)*ZMZ2 +
     *                     XMX2*YMY2*Z*(Z+3.0E0)   )
      GO TO 9999
C
  200 CONTINUE
      XMX2 = X*(1.0E0 - X)
      YMY2 = Y*(1.0E0 - Y)
      ZMZ2 = Z*(1.0E0 - Z)
      PRHS = -EXP(X+Y+Z)*( X*(X+3.0E0)*YMY2*ZMZ2 +
     *                     XMX2*Y*(Y+3.0E0)*ZMZ2 +
     *                     XMX2*YMY2*Z*(Z+3.0E0)   )
      PRHS = PRHS - 5.0E0*TRUE(X,Y,Z)
      GO TO 9999
C
  300 CONTINUE
      PRHS = D2P(X)*P(Y)*P(Z) + P(X)*D2P(Y)*P(Z) + P(X)*P(Y)*D2P(Z)
      GO TO 9999
C
  400 CONTINUE
      PRHS = 100.0E0*Q(Y,20.0E0) - 25.0E0*Q(Z,5.0E0)
      GO TO 9999
C
  500 CONTINUE
      C = -0.18750E0
      A = 0.750E0
      AM2 = A - 2.0E0
      IF ((X .NE. 0.0E0) .AND. (Y .NE. 0.0E0) .AND. (Z .NE. 0.0E0)) THEN
         PRHS = C*( X**AM2*(Y**A-Y)*(Z**A-Z)
     *            + Y**AM2*(X**A-X)*(Z**A-Z)
     *            + Z**AM2*(X**A-X)*(Y**A-Y) )
      ELSE
         PRHS = 0.0E0
      ENDIF
      GO TO 9999
C
  600 CONTINUE
      B = 4.0E0
      B2 = B*B
      PRHS = - 2.0E0*(B2+10.0E0)*SIN(B*(X-Y))
     *       -       (B2+20.0E0)*(COS(B*Y)+COS(B*Z))
      GO TO 9999
C
  700 CONTINUE
      S133 = SQRT(133.0E0)
      RK1 = SQRT(14.0E0+S133)
      RK2 = SQRT(14.0E0-S133)
      A = (-7.0E0+S133)/(2.0E0*S133)
      B = (-7.0E0-S133)*A/16.0E0
      EK1X = EXP(RK1*X)
      EK2X = EXP(RK2*X)
      EDIFF = EK1X-EK2X
      F1 = A*EDIFF + EK2X
      F2 = B*EDIFF
      Y21 = 1.0E0 - Y*Y
      Z1 = 1.0E0 - Z
      DEK1X = RK1*RK1*EK1X
      DEK2X = RK2*RK2*EK2X
      DEDIFF = DEK1X-DEK2X
      DDF1 = A*DEDIFF + DEK2X
      DDF2 = B*DEDIFF
      UXX = 0.50E0*( Y21*(DDF1+Y21*DDF2) + Z1*(DDF1+Z1*DDF2) )
      UYY = -(F1 + 2.0E0*(1.0E0-3.0E0*Y*Y)*F2)
      UZZ = F2
      PRHS = UXX + UYY + UZZ
      GO TO 9999
C
  800 CONTINUE
      PRHS= 3.750E0*(SQRT(X*(Y*Z)**5)+SQRT(Y*(X*Z)**5)+SQRT(Z*(X*Y)**5))
      GO TO 9999
C
 1000 CONTINUE
      PRHS = -5.0E0*X*Y*Z
      GO TO 9999
C
 1100 CONTINUE
      PRHS = 0.0E0
      GO TO 9999
C
 9999 CONTINUE
      RETURN
      END
      REAL FUNCTION BRHS (K,X,Y,Z)
C
C  --------------------------------------
C  RIGHT-HAND SIDE OF BOUNDARY CONDITIONS
C  --------------------------------------
C
      REAL X,Y,Z
C
      COMMON / SELECT / IPROB
C
      GO TO (100,200,300,400,100,100,700,800,900,1000,1100), IPROB
C
  100 CONTINUE
      BRHS = TRUE(X,Y,Z)
      GO TO 9999
C
  200 CONTINUE
      IF ((K .EQ. 1) .OR. (K .EQ. 3)) THEN
         BRHS = (1.0-X-X**2)*Y*(1.0-Y)*Z*(1.0-Z)*EXP(X+Y+Z)
      ELSE IF (K .EQ. 2) THEN
              BRHS = (1.0-Y-Y**2)*X*(1.0-X)*Z*(1.0-Z)*EXP(X+Y+Z)
           ELSE IF (K .EQ. 6) THEN
                   BRHS = (1.0-Z-Z**2)*X*(1.0-X)*Y*(1.0-Y)*EXP(X+Y+Z)
      ELSE
         BRHS = TRUE(X,Y,Z)
      ENDIF
      GO TO 9999
C
  300 CONTINUE
      IF ((K .EQ. 1) .OR. (K .EQ. 4) .OR. (K .EQ. 5)) THEN
         BRHS = 0.0E0
      ELSE
         BRHS = TRUE(X,Y,Z)
      ENDIF
      GO TO 9999
C
  400 CONTINUE
      IF ((K .EQ. 2) .OR. (K .EQ. 3) .OR. (K .EQ. 6)) THEN
         BRHS = 0.0E0
      ELSE
         BRHS = TRUE(X,Y,Z)
      ENDIF
      GO TO 9999
C
  700 CONTINUE
      IF ((K .EQ. 2) .OR. (K .EQ. 4)) THEN
         S133 = SQRT(133.0E0)
         RK1 = SQRT(14.0E0+S133)
         RK2 = SQRT(14.0E0-S133)
         A = (-7.0E0+S133)/(2.0E0*S133)
         EK1X = EXP(RK1*X)
         EK2X = EXP(RK2*X)
         EDIFF = EK1X-EK2X
         F1 = A*EDIFF + EK2X
         BRHS = -Y*F1
      ELSE
         BRHS = TRUE(X,Y,Z)
      ENDIF
      GO TO 9999
C
  800 CONTINUE
      IF ((K .EQ. 2) .OR. (K .EQ. 3)) THEN
         BRHS = 0.0E0
      ELSE
         BRHS = TRUE(X,Y,Z)
      ENDIF
      GO TO 9999
C
  900 CONTINUE
      B = 4.0E0
      IF ((K .EQ. 2) .OR. (K .EQ. 6)) THEN
         BRHS = TRUE(X,Y,Z)
      ELSE IF (K .EQ. 4) THEN
         BRHS = -B*( SIN(B*Y) + COS(B*(X-Y)) )
      ELSE
         BRHS = -B*SIN(B*Z)
      ENDIF
      GO TO 9999
C
 1000 CONTINUE
      IF ((K .EQ. 1) .OR. (K .EQ. 3)) THEN
         BRHS = Y*Z
      ELSE IF (K .EQ. 4) THEN
         BRHS = X*Z
      ELSE IF ((K .EQ. 5) .OR. (K .EQ. 6)) THEN
         BRHS = X*Y
      ELSE
         BRHS = TRUE(X,Y,Z)
      ENDIF
      GO TO 9999
C
 1100 CONTINUE
      BRHS = 0.0E0
      GO TO 9999
C
 9999 CONTINUE
      RETURN
      END
      REAL FUNCTION TRUE(X,Y,Z)
C
C  -------------
C  TRUE SOLUTION
C  -------------
C
      REAL X,Y,Z
C
      COMMON /SELECT/ IPROB
C
      GO TO (100,200,300,400,500,600,700,800,600,1000,1100), IPROB
C
         TRUE = 0.0E0
         GO TO 9999
C
  100 CONTINUE
      TRUE = EXP(X+Y+Z)*X*Y*Z*(1.0E0-X)*(1.0E0-Y)*(1.0E0-Z)
      GO TO 9999
C
  200 CONTINUE
      TRUE = EXP(X+Y+Z)*X*Y*Z*(1.0E0-X)*(1.0E0-Y)*(1.0E0-Z)
      GO TO 9999
C
  300 CONTINUE
      TRUE = P(X)*P(Y)*P(Z)
      GO TO 9999
C
  400 CONTINUE
      TRUE = ( Q(X,10.0E0) + Q(Y,20.0E0) + Q(Z,5.0E0) )/3.0E0
      GO TO 9999
C
  500 CONTINUE
      A = 0.750E0
      TRUE = (X**A-X)*(Y**A-Y)*(Z**A-Z)
      GO TO 9999
C
  600 CONTINUE
      B = 4.0E0
      TRUE = COS(B*Y) + SIN(B*(X-Y)) + COS(B*Z)
      GO TO 9999
C
  700 CONTINUE
      S133 = SQRT(133.0E0)
      RK1 = SQRT(14.0E0+S133)
      RK2 = SQRT(14.0E0-S133)
      A = (-7.0E0+S133)/(2.0E0*S133)
      B = (-7.0E0-S133)*A/16.0E0
      EK1X = EXP(RK1*X)
      EK2X = EXP(RK2*X)
      EDIFF = EK1X-EK2X
      F1 = A*EDIFF + EK2X
      F2 = B*EDIFF
      Y21 = 1.0E0 - Y*Y
      Z1 = 1.0E0 - Z
      TRUE = 0.50E0*( Y21*(F1 + Y21*F2) + Z1*(F1 + Z1*F2) )
      GO TO 9999
C
  800 CONTINUE
      TRUE = SQRT((X*Y*Z)**5)
      GO TO 9999
C
 1000 CONTINUE
      TRUE = X*Y*Z
      GO TO 9999
C
 1100 CONTINUE
      TRUE = 0.0E0
      GO TO 9999
C
 9999 CONTINUE
      RETURN
      END
      REAL FUNCTION P (Z)
      Z1 = 0.150E0
      Z2 = 0.850E0
      IF (Z .LE. Z1)  GO TO 10
      IF (Z .GE. Z2)  GO TO 20
      DZ = Z2 - Z1
      P = 1.0E0 - (Z-Z1)**3*(1.0E0 - 3.0E0*(Z-Z2)/DZ
     *              + 6.0E0*(Z-Z2)**2/DZ**2)/DZ**3
      RETURN
   10 CONTINUE
      P = 1.0E0
      RETURN
   20 CONTINUE
      P = 0.0E0
      RETURN
      END
      REAL FUNCTION D2P (Z)
      Z1 = 0.150E0
      Z2 = 0.850E0
      IF (Z .LE. Z1)  GO TO 10
      IF (Z .GE. Z2)  GO TO 10
      DZ = Z2 - Z1
      C3 = -1.0E0/DZ**3
      C4 = 3.0E0/DZ**4
      C5 = -6.0E0/DZ**5
      ZMZ1 = Z - Z1
      ZMZ2 = Z - Z2
      ZMZ12 = ZMZ1*ZMZ1
      ZMZ22 = ZMZ2*ZMZ2
      D2P = 6.0E0*(C3*ZMZ1 + C4*ZMZ1*ZMZ2 + C4*ZMZ12 + C5*ZMZ1*ZMZ22
     *               + 2.0E0*C5*ZMZ12*ZMZ2) + 2.0E0*C5*ZMZ12*ZMZ1
      RETURN
   10 CONTINUE
      D2P = 0.0E0
      RETURN
      END
      REAL FUNCTION Q (Z,A)
      Q = COSH(A*Z)/COSH(A)
      RETURN
      END
      SUBROUTINE TIMER (T)
C
C  ------------------------------------------------
C  RETURNS ELAPSED CP TIME SINCE START OF JOB (SEC)
C  ------------------------------------------------
C
      REAL T
      T = SECOND()
      RETURN
      END
