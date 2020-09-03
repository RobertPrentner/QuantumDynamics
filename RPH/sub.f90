       SUBROUTINE NORMEQ ( nx, x, f, nfit, alpha, beta, y, chisq, funcs, &
                         vv, b )

!  30.5.1996
! SET UP THE 'NORMAL EQUATIONS' and solve them.

! Input:
!   NX     : TOTAL NUMBER OF POINTS
!   X      : coordinate values
!   F      : function values
!   NFIT   : number of parameters
!   FUNCS  : model function (external)
! Output:
!   Y      : covariance matrix
!   BETA   : parameter values
! Working arrays:
!   ALPHA, VV, B

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION x(nx), f(nx), alpha(nfit,nfit), beta(nfit), &
               y(nfit,nfit), vv(nfit), b(nfit)
      external funcs

      DO 250 K = 1 , nfit
         BETA(K) = 0.d0
         DO 250 L = 1 , K
  250       ALPHA(K,L)=0.d0

      DO 270 I=1,nx
         call funcs ( x(i), vv, b, nfit )
         DO 270 K=1,nfit
            BETA(K)=BETA(K)+f(i)*vv(k)
            DO 270 L=1,K
  270          ALPHA(K,L)=ALPHA(K,L)+vv(k)*vv(L)
      DO 280 K = 1 , nfit
         DO 280 L = K+1 , nfit
  280       ALPHA(K,L) = ALPHA(L,K)

! RATIONAL CHOLESKY-DECOMPOSITION AND BACKSUBSTITUTION FOR MATRIX INVERSION
      CALL RCDCMP (Alpha,Nfit,Nfit,VV,IERR)
      CALL RCBKSB (Alpha,Nfit,Nfit,Beta,VV)
      DO 314 J=1,Nfit
         DO 313 I=1,Nfit
  313       B(I) = 0.0d0
         b(j) = 1.0d0
         CALL RCBKSB (Alpha,Nfit,Nfit,B,VV)
         DO 314 I=1,Nfit
  314       Y(I,J)=B(I)

      CHISQ=0.d0
      DO 370 I=1,nx
         call funcs ( x(i), vv, b, nfit )
         df = f(i)
         DO 360 K=1,nfit
  360       df = df - beta(k)*vv(k)
         if (i.eq.1) then
            mind = dabs(df)
            maxd = dabs(df)
         else
            if(mind.lt.dabs(df))then
             mind = mind
            else
             mind = dabs(df)
            endif
            if(maxd.gt.dabs(df))then
             maxd = maxd
            else
             maxd = dabs(df)
            endif
!            mind = min(mind,dabs(df))
!            maxd = max(maxd,dabs(df))
         end if
  370    chisq = chisq + df**2

      write (*,*) '     min. deviation = ',mind
      write (*,*) '     max. deviation = ',maxd
      write (*,*) '     rms  deviation = ',dsqrt(chisq/nx)

      END

!************************************************************************

      SUBROUTINE RCDCMP(A, NDIM, N, VV, IERR)
!
! Rational Cholesky decomposition
! from GWDG Marquardt routine (U. Dierk)

! Ulrich Schmitt, 16 May 89

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      DIMENSION A(NDIM,N), VV(N)

      IERR = 0

      DO 14 K = 1, N

        DO 12 L = 1, K
          L1 = L - 1
          S = A(K,L)
          IF(L .LT. K) THEN
            IF(L1 .GT. 0) THEN
              DO 11 I = 1, L1
11              S = S - A(I,L) * A(I,K)
            END IF
            A(L,K) = S
          END IF
12      CONTINUE

        IF(S .EQ. 0.0D0) THEN
          VV(K) = 0.0D0
        ELSE
          IF(L1 .GT. 0) THEN
            DO 13 I = 1, L1
              T = A(I,K)
              A(I,K) = T * VV(I)
13            S = S - T * A(I,K)
          END IF
          IF(S .LE. 0.0D0) THEN
            IERR = 3
            write (*,*) ' RCDCMP: A not positive definite.'
            RETURN
          END IF
          VV(K) = 1.0D0 / S
        END IF

14    CONTINUE

      END

!************************************************************************

      SUBROUTINE RCBKSB(A, NDIM, N, B, VV)

! Rational Cholesky backsubstitution
! from GWDG Marquardt routine (U. Dierk)

! Ulrich Schmitt, 16 May 89

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      DIMENSION A(NDIM,N), B(N), VV(N)

      IF(N .GT. 1) THEN
        DO 12 L = 2, N
          DO 12 K = 1, L-1
12          B(L) = B(L) - A(K,L) * B(K)
      END IF

      B(N) = B(N) * VV(N)

      IF(N .GT. 1) THEN
        DO 14 L = N-1, 1, -1
          B(L) = B(L) * VV(L)
          DO 14 K = L+1, N
14          B(L) = B(L) - A(L,K) * B(K)
      END IF

      END
      SUBROUTINE TRED2(NM,N,A,D,E,Z)
!
      INTEGER I,J,K,L,N,II,NM,JP1
      DOUBLE PRECISION A(NM,N),D(N),E(N),Z(NM,N)
      DOUBLE PRECISION F,G,H,HH,SCALE
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2,
!     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A
!     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING
!     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
!
!     ON INPUT
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT.
!
!        N IS THE ORDER OF THE MATRIX.
!
!        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE
!          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
!
!     ON OUTPUT
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.
!
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
!          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
!
!        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX
!          PRODUCED IN THE REDUCTION.
!
!        A AND Z MAY COINCIDE.  IF DISTINCT, A IS UNALTERED.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
!     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
!
!     THIS VERSION DATED AUGUST 1983.
!
!     ------------------------------------------------------------------
!
      DO 100 I = 1, N
!
         DO 80 J = I, N
   80    Z(J,I) = A(J,I)
!
         D(I) = A(N,I)
  100 CONTINUE
!
      IF (N .EQ. 1) GO TO 510
!     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0D0
         SCALE = 0.0D0
         IF (L .LT. 2) GO TO 130
!     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + DABS(D(K))
!
         IF (SCALE .NE. 0.0D0) GO TO 140
  130    E(I) = D(L)
!
         !WRITE (*,'(4x,a,f9.4)') 'scale zero',e(i)

         DO 135 J = 1, L
            D(J) = Z(L,J)
            Z(I,J) = 0.0D0
            Z(J,I) = 0.0D0
  135    CONTINUE
!
         GO TO 290
!
  140    DO 150 K = 1, L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
  150    CONTINUE
!
         F = D(L)
         G = -DSIGN(DSQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
!     .......... FORM A*U ..........
         DO 170 J = 1, L
  170    E(J) = 0.0D0
!
         DO 240 J = 1, L
            F = D(J)
            Z(J,I) = F
            G = E(J) + Z(J,J) * F
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
!
            DO 200 K = JP1, L
               G = G + Z(K,J) * D(K)
               E(K) = E(K) + Z(K,J) * F
  200       CONTINUE
!
  220       E(J) = G
  240    CONTINUE
!     .......... FORM P ..........
         F = 0.0D0
!
         DO 245 J = 1, L
            E(J) = E(J) / H
            F = F + E(J) * D(J)
  245    CONTINUE
!
         HH = F / (H + H)
!     .......... FORM Q ..........
         DO 250 J = 1, L
  250    E(J) = E(J) - HH * D(J)
!     .......... FORM REDUCED A ..........
         DO 280 J = 1, L
            F = D(J)
            G = E(J)
!
            DO 260 K = J, L
  260       Z(K,J) = Z(K,J) - F * E(K) - G * D(K)
!
            D(J) = Z(L,J)
            Z(I,J) = 0.0D0
  280    CONTINUE
!
  290    D(I) = H
  300 CONTINUE
!     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........
      DO 500 I = 2, N
         L = I - 1
         Z(N,L) = Z(L,L)
         Z(L,L) = 1.0D0
         H = D(I)
         IF (H .EQ. 0.0D0) GO TO 380
!
         DO 330 K = 1, L
  330    D(K) = Z(K,I) / H
!
         DO 360 J = 1, L
            G = 0.0D0
!
            DO 340 K = 1, L
  340       G = G + Z(K,I) * Z(K,J)
!
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * D(K)
  360    CONTINUE
!
  380    DO 400 K = 1, L
  400    Z(K,I) = 0.0D0
!
  500 CONTINUE
!
  510 DO 520 I = 1, N
         D(I) = Z(N,I)
         Z(N,I) = 0.0D0
  520 CONTINUE
!
      Z(N,N) = 1.0D0
      E(1) = 0.0D0
      RETURN
      END
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
!
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      DOUBLE PRECISION D(N),E(N),Z(NM,N)
      DOUBLE PRECISION C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2,TST1,TST2,PYTHAG
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,
!     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
!     WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
!
!     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
!     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.
!     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
!     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
!     FULL MATRIX TO TRIDIAGONAL FORM.
!
!     ON INPUT
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT.
!
!        N IS THE ORDER OF THE MATRIX.
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
!
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
!          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
!
!        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
!          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
!          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
!          THE IDENTITY MATRIX.
!
!      ON OUTPUT
!
!        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
!          ERROR EXIT IS MADE, THE EIGENVALUES ARE !ORRECT BUT
!          UNORDERED FOR INDICES 1,2,...,IERR-1.
!
!        E HAS BEEN DESTROYED.
!
!        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
!          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
!          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
!          EIGENVALUES.
!
!        IERR IS SET TO
!          ZERO       FOR NORMAL RETURN,
!          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
!                     DETERMINED AFTER 30 ITERATIONS.
!
!     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
!     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
!
!     THIS VERSION DATED AUGUST 1983.
!
!     ------------------------------------------------------------------
!
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
!
      DO 100 I = 2, N
  100 E(I-1) = E(I)
!
      F = 0.0D0
      TST1 = 0.0D0
      E(N) = 0.0D0
!
      DO 240 L = 1, N
         J = 0
         H = DABS(D(L)) + DABS(E(L))
         IF (TST1 .LT. H) TST1 = H
!     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         DO 110 M = L, N
            TST2 = TST1 + DABS(E(M))
            IF (TST2 .EQ. TST1) GO TO 120
!     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
!
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
!     .......... FORM SHIFT ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * E(L))
         R = PYTHAG(P,1.0D0)
         D(L) = E(L) / (P + DSIGN(R,P))
         D(L1) = E(L) * (P + DSIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145
!
         DO 140 I = L2, N
  140    D(I) = D(I) - H
!
  145    F = F + H
!     .......... QL TRANSFORMATION ..........
         P = D(M)
         C = 1.0D0
         C2 = C
         EL1 = E(L1)
         S = 0.0D0
         MML = M - L
!     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            R = PYTHAG(P,E(I))
            E(I+1) = S * R
            S = E(I) / R
            C = P / R
            P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
!     .......... FORM VECTOR ..........
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
!
  200    CONTINUE
!
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         TST2 = TST1 + DABS(E(L))
         IF (TST2 .GT. TST1) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
!     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
!
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
!
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
!
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
!
  300 CONTINUE
!
      GO TO 1001
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L

	write(*,*)'AN ERROR OCCURED IN THE TRED2',J
 1001 RETURN
      END
      DOUBLE PRECISION FUNCTION PYTHAG(A,B)
      DOUBLE PRECISION A,B
!
!     FINDS DSQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW
!
      DOUBLE PRECISION P,R,S,T,U
      P = DMAX1(DABS(A),DABS(B))
      IF (P .EQ. 0.0D0) GO TO 20
      R = (DMIN1(DABS(A),DABS(B))/P)**2
   10 CONTINUE
         T = 4.0D0 + R
         IF (T .EQ. 4.0D0) GO TO 20
         S = R/T
         U = 1.0D0 + 2.0D0*S
         P = U*P
         R = (S/U)**2 * R
      GO TO 10
   20 PYTHAG = P
      RETURN
      END
!      SUBROUTINE BAKVEC(NM,N,T,E,M,Z,IERR)
!
!      INTEGER I,J,M,N,NM,IERR
!      DOUBLE PRECISION T(NM,3),E(N),Z(NM,M)
!
!     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A NONSYMMETRIC
!     TRIDIAGONAL MATRIX BY BACK TRANSFORMING THOSE OF THE
!     CORRESPONDING SYMMETRIC MATRIX DETERMINED BY  FIGI.
!
!     ON INPUT
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT.
!
!        N IS THE ORDER OF THE MATRIX.
!
!        T CONTAINS THE NONSYMMETRIC MATRIX.  ITS SUBDIAGONAL IS
!          STORED IN THE LAST N-1 POSITIONS OF THE FIRST COLUMN,
!          ITS DIAGONAL IN THE N POSITIONS OF THE SECOND COLUMN,
!          AND ITS SUPERDIAGONAL IN THE FIRST N-1 POSITIONS OF
!          THE THIRD COLUMN.  T(1,1) AND T(N,3) ARE ARBITRARY.
!
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE SYMMETRIC
!          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
!
!        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED.
!
!        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
!          IN ITS FIRST M COLUMNS.
!
!-------------------------------------------------------------     ************
!                                                                      MINV
!                                                                  ************
      SUBROUTINE MINV(AB,N,ND,SCRATCH,DET,EPS,M,MODE)
!
!     A subroutine that calculates the determinant and inverse of
!          a matrix, as well as solving systems of linear equations.
!     Martin J. McBride.  11/25/85.
!     General Electric CRD, Information System Operation.
!
      INTEGER N,ND,M,MODE,OUTER,ROW,COL,I,SCOL,SROW,PIVCNT
      DOUBLE PRECISION AB(ND,1),SCRATCH(1),DET,EPS,MULT,COLNUM,TEMP

!  Initialize scratch space, with 1 to N holding the diagonal of the identity
!  matrix used to compute the inverse and N+1 to 2N holding the positions of
!  the first N columns of the matrix (for use when pivot occurs).
      DO 5 I = 1,N
    5    SCRATCH(I) = 1.0
      COLNUM = 1.0
      DO 6 I = N+1,2*N
         SCRATCH(I) = COLNUM
         COLNUM = COLNUM + 1.0
    6 CONTINUE

!  Make left, square matrix an upper triangular matrix.
      DET = 0.0
      PIVCNT = 0
      DO 10 OUTER = 1,N-1
         IF (ABS(AB(OUTER,OUTER)) .LE. EPS) THEN
            CALL PIVOT(AB,N,ND,OUTER,SCRATCH,EPS)
            IF (AB(OUTER,OUTER) .EQ. 0.0d0) THEN
               PRINT*
               PRINT*,'*************************************'
               PRINT*,'  MINV called with singular matrix.'
               PRINT*,'*************************************'
               PRINT*
               STOP
            ENDIF
            PIVCNT = PIVCNT + 1
         ENDIF
         DO 20 ROW = OUTER+1,N
            MULT = AB(ROW,OUTER)/AB(OUTER,OUTER)
            DO 30 COL = OUTER,N+M
   30          AB(ROW,COL) = AB(ROW,COL) - AB(OUTER,COL)*MULT
            DO 25 SCOL = 1,OUTER-1
   25          AB(ROW,SCOL) = AB(ROW,SCOL) - AB(OUTER,SCOL)*MULT
            AB(ROW,OUTER) = AB(ROW,OUTER) - SCRATCH(OUTER)*MULT
   20    CONTINUE
   10 CONTINUE

!  Compute determinant.
      DET = AB(1,1)
      DO 40 I = 2,N
   40    DET = DET*AB(I,I)
      DET = (-1.0d0)**PIVCNT * DET

!  Return if inverse is not to be found and there are no systems of equations
!  to solve.
      IF (MODE .EQ. 0 .AND. M .EQ. 0) RETURN

!  Place ones in diagonal of square matrix A.
      DO 80 ROW = 1,N
         DIV = AB(ROW,ROW)
         DO 90 COL = 1,N+M
            AB(ROW,COL) = AB(ROW,COL)/DIV
   90    CONTINUE
         SCRATCH(ROW) = SCRATCH(ROW)/DIV
   80 CONTINUE

!  Reduce upper triangle to zeros to give matrix A = I.
      DO 50 OUTER = 2,N
         DO 60 ROW = OUTER-1,1,-1
            MULT = AB(ROW,OUTER)/AB(OUTER,OUTER)
            DO 70 COL = OUTER,N+M
   70          AB(ROW,COL) = AB(ROW,COL) - AB(OUTER,COL)*MULT
            DO 65 SCOL = 1,ROW-1
   65          AB(ROW,SCOL) = AB(ROW,SCOL) - AB(OUTER,SCOL)*MULT
            SCRATCH(ROW) = SCRATCH(ROW) - AB(OUTER,ROW)*MULT
            DO 63 SCOL = ROW+1,OUTER-1
   63          AB(ROW,SCOL) = AB(ROW,SCOL) - AB(OUTER,SCOL)*MULT
            AB(ROW,OUTER) = AB(ROW,OUTER) - SCRATCH(OUTER)*MULT
   60    CONTINUE
   50 CONTINUE

!  Move diagonals of inverse to matrix AB.
      DO 85 I = 1,N
   85    AB(I,I) = SCRATCH(I)

!  If pivot was made, switch rows corresponding to the columns that were
!  pivoted.
      IF (PIVCNT .EQ. 0) RETURN
      ROW = 1
      DO 95 I = 1,N-1
         SROW = INT(SCRATCH(ROW+N))
         IF (SROW .NE. ROW) THEN
            DO 92 COL = 1,N+M
               TEMP = AB(ROW,COL)
               AB(ROW,COL) = AB(SROW,COL)
               AB(SROW,COL) = TEMP
   92       CONTINUE
            TEMP = SCRATCH(ROW+N)
            SCRATCH(ROW+N) = SCRATCH(SROW+N)
            SCRATCH(SROW+N) = TEMP
         ELSE
            ROW = ROW + 1
         ENDIF
   95 CONTINUE

      RETURN
      END

      SUBROUTINE PIVOT(AB,N,ND,OUTER,SCRATCH,EPS)
!
!     This subroutine switches two columns of a matrix to get
!         a nonzero entry in the diagonal.
!     Martinrow OUTER.                                                    
!     General Electric CRD, Information System Operation.
!
      INTEGER N,ND,COL,OUTER,I
      DOUBLE PRECISION AB(ND,1),SCRATCH(1),TEMP,EPS       
!  Get first column with non-zero element in row OUTER.
      COL = OUTER + 1
   10 IF (COL .GT. N) GO TO 90
      IF (ABS(AB(OUTER,COL)) .GT. EPS) GO TO 20 
         COL = COL + 1
         GO TO 10         

!  Switch column OUTER with column COL, which has non-zero element in 
!  row OUTER.
   20 DO 30 I = 1,N
         TEMP = AB(I,OUTER)
         AB(I,OUTER) = AB(I,COL)
         AB(I,COL) = TEMP
   30 CONTINUE
      TEMP = SCRATCH(N+OUTER)
      SCRATCH(N+OUTER) = SCRATCH(N+COL)
      SCRATCH(N+COL) = TEMP

   90 CONTINUE
      RETURN
      END
      SUBROUTINE  SSCAL(N,DA,DX,INCX)
!
!     SCALES A VECTOR BY A CONSTANT.
!     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
!     JACK DONGARRA, 

      DOUBLE PRECISION DA,DX(*)
      INTEGER I,INCX,M,MP1,N,IX

      IF (N .LT. 0) STOP     
      IF (N .EQ. 0) RETURN
      IF (INCX .EQ. 1) GO TO 20
!
!        CODE FOR INCREMENT NOT EQUAL TO 1
!
      IX = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        DX(IX) = DA*DX(IX)
        IX = IX + INCX
   10 CONTINUE                             
      return
!
!        CODE FOR INCREMENT EQUAL TO 1
!
!
!        CLEAN-UP LOOP
!
   20 M =MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN       
   40 MP1 = M + 1 
      DO 50 I = MP1,N,5       
        DX(I) = DA*DX(I)
        DX(I+1) = DA*DX(I+1)
        DX(I+2) = DA*DX(I+2)
        DX(I+3) = DA*DX(I+3)
        DX(I+4) = DA*DX(I+4)
   50 CONTINUE
      RETURN
      END 

























