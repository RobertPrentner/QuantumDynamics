	MODULE marqfit
	
	CONTAINS
	
	SUBROUTINE SetupFit(tau,fseries,v0,nq,order,flag,ntol,nopt,rmsr,sigmav0)
	
!	Parameters
	DOUBLE PRECISION, PARAMETER ::  tol = 1.d32
	INTEGER, PARAMETER :: maxopt = 1140, maxtry=14,maxnfit=2048      
!	Input				!	
	INTEGER :: nq,order		!  nq = number of points on rp
	INTEGER :: flag
	DOUBLE PRECISION :: tau(nq+1),fseries(order)
	DOUBLE PRECISION :: v0(nq)
	DOUBLE PRECISION, OPTIONAL :: sigmav0(nq)
!	Output
!	Marquardt fitted Fourie series in fseries, deviations from initial function in sigmav0

!	working variables		
	INTEGER :: q,i,j,ifit,ntol,nopt,newq=0
	DOUBLE PRECISION :: f(maxnfit),p(order),g(order,nq),sigma(nq), &
					 alpha(order,order), beta(order),x(maxnfit)
	DOUBLE PRECISION :: pi,dtau,chisq,sumw,rms,rmsr,dev
	CHARACTER*8 func

        pi=DACOS(-1.d0)
	DO q=1,nq	!initialize, set uncertainties to default value
		f(q)=0.d0
		Do ifit=1,order
		g(ifit,q)=0.d0
		ENDDO
		IF(PRESENT(sigmav0)) THEN; sigma(q)=sigmav0(q)
		ELSE; sigma(q) = 1.d0
		ENDIF
	ENDDO
	
	SELECT CASE (flag)
		CASE (1:2)
		tau(nq+1)=360.d0
			DO ifit=1,order			 !Initialize parameter with discrete cos coefficients
				p(ifit)=0.d0
				DO q=1,nq
					dtau=ABS(tau(q+1)-tau(q))
					IF(ifit.eq.1) THEN
						p(ifit) = p(ifit) + v0(q)/180.d0*dtau
					ELSEIF(ifit.eq.2) THEN
						p(ifit) = + v0(q)/180.d0*&
						    (DCOS(0.5*pi*tau(q)/180.d0)+pi/2)*dtau 
					ELSE 
						p(ifit) = p(ifit) + v0(q)/180.d0*&
						    (DCOS((ifit-2)*pi*tau(q)/180.d0))*dtau 	     
					ENDIF
				ENDDO
			ENDDO
			p(1) =p(1)/2.d0

		CASE (0)
		tau(nq+1)=360.d0
			DO ifit=1,order			 !Initialize parameter with discrete fourier coefficients
				p(ifit)=0.d0
				DO q=1,nq
					dtau=ABS(tau(q+1)-tau(q))
					IF(ifit.eq.1) p(ifit) = p(ifit) + v0(q)/180.d0*dtau
					IF(MOD(ifit,2).EQ.0) THEN
						p(ifit) = p(ifit) + v0(q)/180.d0*&
							(DCOS((ifit-1)*pi*tau(q)/180.d0))*dtau 	     
					ELSE
						p(ifit) = p(ifit) + v0(q)/180.d0*&
							(DSIN((ifit-1)*pi*tau(q)/180.d0))*dtau 	  
					ENDIF				
				ENDDO
			ENDDO

		CASE (-2:-1)
		tau(nq+1)=360.d0
!			DO ifit=1,order			 !Initialize parameter with discrete sin coefficients + constant offset
!				p(ifit)=0.d0
!				DO q=1,nq
!					dtau=ABS(tau(q+1)-tau(q))
!					p(ifit) = p(ifit) + v0(q)/180.d0*&
!						    (DSIN((ifit)*pi*tau(q)/180.d0))*dtau 	     
!				ENDDO
!			ENDDO
			DO ifit=1,order			 !Initialize parameter with discrete sin coefficients + constant offset
				p(ifit)=0.d0
				DO q=1,nq
					dtau=ABS(tau(q+1)-tau(q))
					IF(ifit.eq.1) THEN 
						p(ifit) = p(ifit) + v0(q)/180.d0*dtau
					ELSEIF(ifit.eq.2) THEN
						p(ifit) = p(ifit)+ v0(q)/180.d0*&
						    (DSIN(0.5*pi*tau(q)/180.d0)+pi/2)*dtau 	
					ELSE
						p(ifit) = p(ifit) + v0(q)/180.d0*&
						    (DSIN((ifit-2)*pi*tau(q)/180.d0))*dtau 	     
					ENDIF	    
				ENDDO
			ENDDO
			p(1) =0.d0 !no constant for sine
		CASE DEFAULT 
			STOP 'marqfit: error with symmetry'
	END SELECT	

!	WRITE (*,*) 'Initial parameters (derived from cos series):'
!	WRITE (*,'(5f12.5)') (p(ifit),ifit=1,order)
!	WRITE (*,*) 'initial guess for function'
!	WRITE (*,*) 'coordinate	function	estimate'
!	WRITE (*,*)
	!p(3)=0.0
	DO q=1,nq
		f(q)= expansion(flag,tau(q),p(1),g(1,q),order)
!				f(q) = fCos(tau(q),p(1),g(1,q),order)
!		WRITE (*,'(f12.3,2f12.5)') tau(q),v0(q),f(q)
	ENDDO
        
	CALL NORMEQ (tol,chisq,sumw,rms,order,nq,ntol, &
                               v0,f(1),g(1,1),sigma,alpha,beta )

	!beta(3) = 0.0
	CALL SolveNEQ (TOL,CHISQ,SUMW,RMS,order,nq,NTOL, &
                          v0,f(1),g(1,1),sigma,alpha,beta,p,tau,flag,nopt,rmsr)

!	IF(PRESENT(nfit)) THEN
		
!		WRITE (*,*) 'overwrite function values with fitted ones'
!		DO q=1,nfit
!		     x(q)=q*10.d0!dble(q-1)*360.d0/dble(nfit)+360.d0/(2.d0*dble(nfit))
!		     f(q)=fCos(x(q),p(1),g(1,1),order)
!		     
!		ENDDO
!		DO q=1,nq; WRITE (*,'(f12.3,2g14.6)') tau(q),v0(q),f(q); ENDDO
!	ENDIF
	DO ifit=1,order; fseries(ifit)=p(ifit); ENDDO
	END SUBROUTINE SetupFit
	
!****************************************************************************
!  fNAME is a double precision function that provides functional values
!  with respect to parameters at given argument values.
!  fNAME is called by marqfit and has to have the syntax:
!
!               fNAME ( tau, p, g, order )
!
! Input:
!    tau         : reaction path coordinate
!    p   	        : parameters
!    order     : order of a series expansion fit
!
! Output:
!    fNAME     : functional value
!   g		:gradients iwth respect to paramters (dfNAME/dp(i) )
!***************************************************************************
	DOUBLE PRECISION FUNCTION Expansion (flag,tau, p, g, order)
	
	INTEGER :: order,flag
	double precision :: tau
	DOUBLE PRECISION, DIMENSION(order) :: g, p
		
	SELECT Case (flag)
	CASE(2); expansion=fCos2(tau,p,g,order)
	CASE(1); expansion=fCos(tau,p,g,order)
	CASE(-1); expansion=fSin(tau,p,g,order)
	CASE(-2); expansion=fSin2(tau,p,g,order)
	CASE(0); expansion=fSeries(tau,p,g,order)
	CASE DEFAULT; STOP 'EXPANSION: Specify parity wrt reaction path'
	END SELECT
	
	END FUNCTION Expansion

!	SUBROUTINE Expansion (f,flag,tau, p, g, order)
!	
!	INTEGER :: order,flag
!	double precision :: tau,f
!	DOUBLE PRECISION, DIMENSION(order) :: g, p
!		
!	SELECT Case (flag)
!	CASE(2); f=fCos2(tau,p,g,order)
!	CASE(1); f=fCos(tau,p,g,order)
!	CASE(-1); f=fSin(tau,p,g,order)
!	CASE(0); f=fSeries(tau,p,g,order)
!	CASE DEFAULT; STOP 'EXPANSION: Specify parity wrt reaction path'
!	END SELECT
!	
!	END SUBROUTINE Expansion
!*******************************
!  Fourier cosin series:  fCos  *
!*******************************

	double precision function fCos (tau, p, g, order )
	INTEGER :: order,i
	double precision :: tau,a, s
	DOUBLE PRECISION, DIMENSION(order) :: g, p
		
	s=0.d0
	a = tau * DACOS(-1.d0)/180.d0 ! q*pi/180
	g(1) = 1.d0
	g(2) = DCOS(0.5*a+DACOS(-1.d0)/2.) !cos(0.5*x + pi/2)
	s = p(1)*g(1)+p(2)*g(2)
	do i = 3 , order
        g(i) = DCOS((i-2)*a) ! cos(i*q*pi/180)
        s = s + g(i)*p(i)
	end do
	fCos = s
! change fx_fourier in fx
	end FUNCTION fCos

!*******************************
!  Fourier cosin series:  fCos2, symmetric around pi/2 (f(pi/2-x) = f(pi/2+x) *
!*******************************

	double precision function fCos2 (tau, p, g, order )
	INTEGER :: order,i
	double precision :: tau,a, s
	DOUBLE PRECISION, DIMENSION(order) :: g, p
		
	s=0.d0
	a = tau * DACOS(-1.d0)/180.d0 ! q*pi/180
	g(1) = 1.d0
	g(2) = 0.
	s = p(1)*g(1)
	do i = 3 , order,2
        g(i) = DCOS((i-2)*a)
	g(i+1) = DCOS((i-1)*a)
	!g(i+1) = 0.0
	p(i) = 0.0
        s = s + g(i)*p(i)+g(i+1)*p(i+1)
	end do
	fCos2 = s
! change fx_fourier in fx
	end FUNCTION fCos2
	
!*******************************
!  Fourier sin series:  fSin2  *
!*******************************

	double precision function fSin2 (tau, p, g, order )
	INTEGER :: order,i
	double precision :: tau,a, s
	DOUBLE PRECISION, DIMENSION(order) :: g, p
	
	s=0.d0
	a = tau * DACOS(-1.d0)/180.d0 ! q*pi/180
	g(1) = 0.0
	g(2) = 0.0
	s = p(1)*g(1)
	do i = 3, order,2
        g(i) = DSIN((i-2)*a) ! sin(i*q*pi/180)
	g(i+1) = DSIN((i-1)*a)
	p(i)=0.
        s = s + g(i)*p(i)+g(i+1)*p(i+1)
	end do
	fSin2 = s
! change fx_fourier in fx
	end FUNCTION fSin2

!*******************************
!  Fourier sin series:  fSin  *
!*******************************

	double precision function fSin (tau, p, g, order )
	INTEGER :: order,i
	double precision :: tau,a, s
	DOUBLE PRECISION, DIMENSION(order) :: g, p
	
	s=0.d0
	a = tau * DACOS(-1.d0)/180.d0 ! q*pi/180
	g(1) = 0.0
	g(2) = DSIN(0.5*a+DACOS(-1.d0)/2.) !sin(0.5*x + pi/2)
	s = p(1)*g(1)+p(2)*g(2)
	do i = 3 , order
        g(i) = DSIN((i-2)*a) ! cos(i*q*pi/180)
        s = s + g(i)*p(i)
	end do
	fSin = s
! change fx_fourier in fx
	end FUNCTION fSin

!*******************************
!  Fourier series: fSeries  *
!*******************************

	double precision function fSeries ( x, g, p, np )
	implicit double precision ( a-h,o-z )
	dimension g(*), p(*)
	
	s=0.d0
	a = x * dacos(-1.d0)/180.d0
	g(1) = 1.d0
	s = p(1)
	do i = 2 , np, 2
        g(i) = dcos(i/2*a)
        s = s + g(i) * p(i)
        if (i.lt.np) then
          g(i+1) = dsin(i/2*a)
          s = s + g(i+1) * p(i+1)
        end if
	end do
	fSeries = s
! change fx_fourier in fx
	end FUNCTION fSeries
	
!************************************************************************
! DOUBLE PRECISION VERSION OF  28.2.1989, by David Luckhaus

! SET UP THE 'NORMAL EQUATIONS' FOR THE MARQUARDT-FIT :
!              ALPHA * X = BETA

! TOL    : MAXIMUM ALLOWED DEVIATION OF OBSERVED AND CALCULATED POSITIONS
! nx     : TOTAL NUMBER OF POINTS
! NTOL   : NUMBER OF ASSIGNMENTS CONSIDERED : ABS( NU(OBS) - NU(CALC) ) < TOL
! f      : function values
! SIGMA  : UNCERTAINTY OF function values
! fcal   : 
! grad   : gradient of fcal with resp. to parameters
! SUMW   : SUM OF STATISTICAL WEIGHTS
! CHISQ  : (WEIGHTED) SUM OF SQUARED DEVIATIONS
! RMS    : ROOT MEAN SQUARE DEVIATION (NON-WEIGHTED)

	SUBROUTINE NORMEQ (TOL,CHISQ,SUMW,RMS,nfit,nx,NTOL, &
                          f,fcal,grad,sigma,alpha,beta )

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION f(nx), fcal(nx), grad(nfit,nx), sigma(nx), &
		          alpha(nfit,nfit), beta(nfit)

	NTOL=0
	CHISQ=0.d0
	RMS=0.d0
	SUMW=0.d0
	DO 250 K=1,nfit
         BETA(K)=0.d0
         DO 250 L=1,K
  250       ALPHA(K,L)=0.d0

	DO 270 I=1,nx
         SIG2I=1./(SIGMA(I)**2)
         DY=f(i)-fcal(i)

! SKIP LINES WHERE DEVIATIONS ARE TOO LARGE
         IF (ABS(DY).LE.TOL) THEN
            RMS=RMS+DY**2
            NTOL=NTOL+1
         ELSE
            SIG2I=0.
         END IF
         SUMW=SUMW+SIG2I
         CHISQ=CHISQ+SIG2I*DY**2

         DO 270 K=1,nfit
            GK=grad(k,I)*SIG2I
            BETA(K)=BETA(K)+DY*GK
            DO 270 L=1,K
  270          ALPHA(K,L)=ALPHA(K,L)+GK*grad(L,I)

! FILL UPPER TRIANGLE :
	DO 280 K=1,nfit-1
         DO 280 L=K+1,nfit
  280       ALPHA(K,L)=ALPHA(L,K)

	RMS=SQRT(RMS/NTOL)
	IF (SUMW.GT.0.) CHISQ=CHISQ/SUMW


	RETURN
	END SUBROUTINE NORMEQ

!************************************************************************
	SUBROUTINE SolveNEQ (TOL,CHISQ,SUMW,RMS,order,nq,NTOL, &
                          v0,f,grad,sigma,alpha,beta,p,tau,flag,nopt,rmsr)
			  
!	Solve the normal equations an return optimized paramter
!	adapted from mrqfit_sub written by David Luckhaus	
	
	INTEGER, PARAMETER :: maxopt = 1140, maxtry=21
	DOUBLE PRECISION, PARAMETER :: dchbnd=5.d-10
	INTEGER :: i,ii,ij,q,ntry,nopt
	INTEGER, INTENT(IN) :: order,nq,ntol, flag
	DOUBLE PRECISION, INTENT(IN) :: tau(nq), v0(nq)
!	
	DOUBLE PRECISION, DIMENSION(nq) :: f,sigma
	DOUBLE PRECISION, DIMENSION(order) :: beta, beta1,p,pnew,pold,vv
	DOUBLE PRECISION :: alamda,tol,chisq,chisq1,sumw,rms,rmsr,dchi
	DOUBLE PRECISION :: grad(order,nq),alpha(order**2),alpha1(order**2)

	NOPT=1
        NTRY=1
	alamda=1.d-3
  290 CONTINUE	
	DO 300 i=1,order**2
  300    alpha1(i) = alpha(i)
        do 301 i = 1 , order
         ii = (i-1)*order+i
         alpha1(ii) = alpha1(ii) * (1.+alamda)
  301    beta1(i) = beta(i)
	

  305 CALL RCDCMP1 (ALPHA1,order,order,VV,IERR)
	CALL RCBKSB1 (ALPHA1,order,order,BETA1,VV)

! make a predicition with new parameters and calculate CHISQ
        DO i=1,order
		PNEW(i)=P(i)+BETA1(I)
!		write(*,*) PNEW(I), P(I), BETA1(I)
        ENDDO
	DO q=1,nq	!calculate fit for new parameters pnew
		f(q) = expansion(flag,tau(q),pnew(1),grad(1,q),order)
	ENDDO

	CALL NORMEQ (tol,chisq1,sumw,rms,order,nq,ntol,  &
                               v0,f(1),grad(1,1),sigma,alpha1,beta1)

!	convergence test
        IF (chisq1.GT.chisq.AND.ntry.LT.maxtry) then
         ALAMDA=10.*ALAMDA
         NTRY=NTRY+1
         go to 290
        else if (chisq1.GT.chisq) then
          WRITE (*,*) ' no convergence after', ntry, 'cycles! '
          stop 
        ELSE
         NOPT=NOPT+1
         dchi = abs(chisq1-chisq)
         CHISQ=CHISQ1
         DO 400 I=1,order
            P(i)=PNEW(i)
  400       BETA(I)=BETA1(I)
         DO 401 ij=1,order**2
  401       ALPHA(ij)=ALPHA1(IJ)
         alamda = max( 0.1d0 * alamda , 0.001d0 )
         ntry = 1
        END IF
        IF ( dchi.GT.dchbnd*chisq.AND.nopt.LE.maxopt) go to 305
	
!	IF (dchi.LE.dchbnd*chisq) WRITE(*,*) dchi,dchbnd*chisq,dchbnd,chisq

	rmsr = 0.
	do 420 q = 1 , nq
            dev = DABS(v0(q)-f(q))
  420    if ( DABS(f(q)).gt.1.d-10 .and. dev.lt.tol .and. v0(q).gt.0. ) &
		rmsr = rmsr + ( (v0(q)-f(q))/v0(q) )**2
        rmsr = SQRT( rmsr/ntol )
	
	!IF(rmsr.gt.0.01) WRITE (*,'(i4,17x,i4,10x,g11.3)') NOPT,ntol,rmsr
	END SUBROUTINE SolveNEQ
	
	
	SUBROUTINE RCDCMP1(A, NDIM, N, VV, IERR)
!
! Rational Cholesky decomposition
! from GWDG Marquardt routine (U. Dierk)
!
! Ulrich Schmitt, 16 May 89
!
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      DIMENSION A(NDIM,N), VV(N)
!
      IERR = 0
!
      DO 14 K = 1, N
!
        DO 12 L = 1, K
          L1 = L - 1
          S = A(K,L)
          IF(L .LT. K) THEN
            IF(L1 .GT. 0) THEN
              DO 11 I = 1, L1
                S = S - A(I,L) * A(I,K)
11            CONTINUE
            END IF
            A(L,K) = S
          END IF
12      CONTINUE
!
        IF(S .EQ. 0.0D0) THEN
          VV(K) = 0.0D0
        ELSE
          IF(L1 .GT. 0) THEN
            DO 13 I = 1, L1
              T = A(I,K)
              A(I,K) = T * VV(I)
              S = S - T * A(I,K)
13          CONTINUE
          END IF
          IF(S .LE. 0.0D0) THEN
! Not positive definite
            IERR = 3
            RETURN
          END IF
          VV(K) = 1.0D0 / S
        END IF
!
14    CONTINUE
      
      RETURN
      END SUBROUTINE RCDCMP1

!*******************************************************************

      SUBROUTINE RCBKSB1(A, NDIM, N, B, VV)
!
! Rational Cholesky backsubstitution
! from GWDG Marquardt routine (U. Dierk)
!
! Ulrich Schmitt, 16 May 89
!
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      DIMENSION A(NDIM,N), B(N), VV(N)
!
      IF(N .GT. 1) THEN
        DO 12 L = 2, N
          DO 11 K = 1, L-1
            B(L) = B(L) - A(K,L) * B(K)
11        CONTINUE
12      CONTINUE
      END IF
!
      B(N) = B(N) * VV(N)
!
      IF(N .GT. 1) THEN
        DO 14 L = N-1, 1, -1
          B(L) = B(L) * VV(L)
          DO 13 K = L+1, N
            B(L) = B(L) - A(L,K) * B(K)
13        CONTINUE
14      CONTINUE
      END IF
!
      RETURN
      END SUBROUTINE RCBKSB1

!************************************************************************

	END MODULE marqfit
	