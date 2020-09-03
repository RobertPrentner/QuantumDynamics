     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!This program reads in fitted geometries, tensor of inertia and calculates the G-Matrix, 
!!!!!!!!!normal coordinates and frequencies after projection of the force filed matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!input is given as a series of cosine/sine or polynomial (not yet included) functions!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!multiplication and derivatives can therefore be calculated analytically!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PROGRAM RPH
      
      USE STRUCTS
      USE MARQFIT
      INTEGER :: nx,nq,natoms,nx2,nfit,nfitout, nfitout2,flag1,flag2,flagout,nq2      
      INTEGER :: nsch, nvib(3*MAXATOMS),ntime(3),ntime0(3)
      DOUBLE PRECISION ::  g0(maxpoints), h(maxpoints**2*4), h0(maxpoints**2*4),  &
				        u(maxpoints),p(2*maxpoints-1), v(maxpoints), wat(maxpoints), &
				        w(maxpoints*4), v0(maxpoints), &
					omega(3*maxatoms,maxpoints),  g(4,4,maxpoints), &
					b(12,12,4,MAXPOINTS), jsym(maxpoints), &
					pg(maxpoints),pw(maxpoints),pu(maxpoints),pv(maxpoints), &
				        gw(maxpoints), tau(maxpoints),e(MAXPOINTS),wfit,dnsum 
      DOUBLE PRECISION :: qmin,qmax, par1(4*MAXORDER),dq,sg,su,sv,sw,zpe,&
					svib(3*MAXATOMS),sh(3*MAXATOMS,10),elapsed(2),et,f
      DOUBLE PRECISION, PARAMETER :: hbar8=16.85763128d0/4.d0                      ! hbar^2/8hc
    
      
      TYPE(fitgeom) :: s_fit 
      TYPE(gcontra) :: s_Gcontra(4,4) 
      TYPE(vpseudo) :: s_Vpseudo 
      TYPE(energies) :: s_energies(3*maxatoms-7) 
      TYPE(bkl) :: s_bkl(3*maxatoms-7,3*maxatoms-7) 

      CALL itime(ntime0)
      WRITE (*,'(x,a)') 'User input:'
      WRITE (*,'(x,a)') 'number of atoms'
      READ (*,*) natoms
      WRITE (*,'(i4)') natoms
      WRITE (*,'(x,a)') 'number of grid points'
      READ (*,*) nq
      WRITE (*,'(i4)') nq
      WRITE (*,'(x,a)') 'values for qmin, qmax'
      READ (*,*) qmin, qmax
      WRITE (*,'(2F9.3)') qmin, qmax
      WRITE (*,'(x,a)') 'scaling factor of potential'
      READ (*,*) f
      
      
      nx=3*natoms
      !nsel = nx-7
      WRITE (*,*) 'read in data'
      CALL ReadData(nx,nfit,s_fit,s_energies,s_Vpseudo,s_Gcontra,s_bkl(1:nx-7,1:nx-7))
      WRITE (*,*) nfit
!!!!!!!!!!!!!!expand on the grid      
     
      dq = dble((qmax-qmin)/(nq))
      DO q=1,nq
        tau(q)=dble(q-1)*360.d0/dble(nq)+360.d0/(2.d0*dble(nq))!qmin + (qmax-qmin)/(nq)*(q-0.5)
	v(q) = expansion(s_fit%flagv0,tau(q),s_fit%fitv0(1:nfit),par1(1:nfit),nfit)
	DO i=1,4
	    g(i,i,q) = expansion(s_Gcontra(i,i)%flag,tau(q),s_Gcontra(i,i)%series(1:nfit),par1,nfit) 
	  DO j=1,i-1
	    g(i,j,q) = expansion(s_Gcontra(i,j)%flag,tau(q),s_Gcontra(i,j)%series(1:nfit),par1,nfit) 
	    g(j,i,q) = g(i,j,q)
	  ENDDO
	ENDDO
	DO i=1,nx-7
	   omega(i,q) = expansion(s_energies(i)%flag,tau(q),s_energies(i)%series(1:nfit),par1,nfit)	!symmetrize
	   b(i,i,1,q) = 0.0 !must be due to symmetry
	   b(i,i,2,q) = 0.0
	   b(i,i,3,q) = 0.0
	   b(i,i,4,q) = expansion(s_bkl(i,i)%flag4,tau(q),s_bkl(i,i)%series4(1:nfit),par1,nfit)
	   DO j=1,i-1
	    b(i,j,1,q) = expansion(s_bkl(i,j)%flag1,tau(q),s_bkl(i,j)%series1(1:nfit),par1,nfit)
	    b(i,j,2,q) = expansion(s_bkl(i,j)%flag2,tau(q),s_bkl(i,j)%series2(1:nfit),par1,nfit)
	    b(i,j,3,q) = expansion(s_bkl(i,j)%flag3,tau(q),s_bkl(i,j)%series3(1:nfit),par1,nfit)
	    b(i,j,4,q) = expansion(s_bkl(i,j)%flag4,tau(q),s_bkl(i,j)%series4(1:nfit),par1,nfit)
	    b(j,i,1,q)=b(i,j,1,q)
	    b(j,i,2,q)=b(i,j,2,q)
	    b(j,i,3,q)=b(i,j,3,q)
	    b(j,i,4,q)=b(i,j,4,q)
	   ENDDO
	ENDDO
	u(q) = expansion(s_Vpseudo%flag,tau(q),s_Vpseudo%series(1:nfit),par1,nfit) 
	g0(q) = g(4,4,q)
	wat(q) = -hbar8 * ( g(1,1,q) + g(2,2,q) + g(3,3,q) )   ! Watson's term
	v0(q) = v(q) + wat(q) + u(q)
      ENDDO
! tests: 
!      DO i=1,3 
!      WRITE (*,'(4f15.9)') (g(i,j,100),j=1,i    )
!      ENDDO
!      WRITE (*,'(4f15.9)') (g(4,j,100)*(nq/360.)**2,j=1,4)
        !nfit =nq
        !CALL SetupFit(tau,pg,g0,nq,nfit,1,ntol,nopt,rmsr)
	!CALL SetupFit(tau,pv,v0,nq,nfit,1,ntol,nopt,rmsr)
	
	!Do i=1,nq
	!    v0(i) = expansion(1,tau(i),pv,par1,nfit)
	!    v0(nq+1-i) = v0(i)
	!    g0(i) = expansion(1,tau(i),pg,par1,nfit)
	!    g0(nq+1-i) = g0(i)
	!    DO j=1,nx-7
	!      omega(j,nq+1-i) = omega(j,i)
	!    ENDDO
	!ENDDO
	
!  test symmetry w.r.t. inversion at the center of the grid, i.e. (NDAT+1)/2.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      sg = 0.
      su = 0.
      sv = 0.
      sw = 0.
      svib= 0.
      write (*,'(/2x,82("="))')
      write (*,'(13x, "tau", 14x, "G",15x,"V",15x,"U",13x,"Watson"/2x,82("-"))')
      do q = 1 , nq
        sv = sv + (v(q)-v(nq+1-q))**2
        su = su + (u(q)-u(nq+1-q))**2
        sg = sg + (g0(q)-g0(nq+1-q))**2
        sw = sw + (wat(q)-wat(nq+1-q))**2
	DO i=1,nx-7
	   svib(i) = svib(i) + (omega(i,q)-omega(i,nq+1-q))**2
	ENDDO
        write (*,'(i3,5f16.6)') q, tau(q), g0(q)/dq**2, v(q), u(q), wat(q)
        write (99,'(3f16.6)') g0(q), u(q), wat(q)
	WRITE (66,'(5f16.6)') (omega(j,q),j=1,nx-7)
      end do
!      WRITE (99,*)
      write (*,'(2x,82("-"))')
      write (*,*) ' Inversion Symmetry of G : ',sqrt(sg)
      write (*,*) '                       V : ',sqrt(sv)
      write (*,*) '                       U : ',sqrt(su)
      write (*,*) '                       W : ',sqrt(sw)
      write (*,*) ' Inversion Symmetry of vibrational energies: '
      WRITE (*,'(2x,8E10.2)') (Sqrt(svib(i)),i=1,nx-7)
      write (*,'(2x,82("-")/)')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

      call Initp ( nq, p )

      write (*,'(/2x,82("="))')
      write (*,'(66x,"VIBANG")')
      write (*,'(29x,"Tunnelling levels",16x,"min   max   sym")')
      write (*,'(16x,43("-"),2x,18("-"))')


      do nsel = 0, 0!nx-7
        do i = 1 , nx-7
          nvib(i) = 0
        end do
        if ( nsel.gt.0 ) nvib(nsel) = 1

! add harmonic vibrational energy to the potential energy
        call AddVib ( nq, nx-7, nvib, omega(1:(nx-7),1:nq), v0, v, tau,nsel )
        
! Add adiabatic Coriolis ( vibrational angular momentum ) contributions
        call VibAng ( nq, nx-7, nvib, omega(1:(nx-7),1:nq), b(1:nx-7,1:nx-7,1:4,1:nq), g(1:4,1:4,1:nq), w )
        
	!CALL SetupFit(tau,pw,w,nq,nq,1,ntol,nopt,rmsr)
	!CALL SetupFit(tau,pg,g0,nq,nq,1,ntol,nopt,rmsr)
	sc = 0.0
        wmin = w(1)
        wmax = w(1)
	DO i=1,nq/2
	!w(i) = w(nq+1-i)
	ENDDO
        do i = 1 , nq
	  !w(i) = expansion(1,tau(i),pw,par1,nq)
	  !g0(i) = expansion(1,tau(i),pg,par1,nq)
	  write (99,'(i3,f16.6)')i, w(i)
	  v(i) = f*(v(i) + w(i))
          wmin = min(wmin,w(i))
          wmax = max(wmax,w(i))
        end do
	!WRITE (99,*) 
	Do i=1,nq
          sc = sc + 0.5*(w(i)-w(nq+1-i))**2
	ENDDO
        write (*,*) 
        write (*,*) ' VIBANG : Min/Max = ', wmin, wmax
        write (*,*) '          Symmetry w.r.t. inversion : ',sqrt(sc)
       
	! overlap between n.c. wavefunctions along the path in H
        call Overlap ( h, nq, nx-7, nvib, omega(1:(nx-7),1:nq ) )
	! store hamiltonian in h0
        call Ham1d ( nq, e, h, h0,jsym, p, g0, v, w, dq )

        if ( nsel.eq.0 ) then 
          zpe = e(1)                  ! zero point energy
          write (*,'(/6x," Zero point energy : ",g16.8/)') zpe
        end if

        write(*,'(1x,i3,f10.2,4f11.5,2f6.1,2x,ES)') &
                nsel, e(1)-zpe, ( e(i)-e(1), i=2,5), wmin, wmax, sc
        !write(*,'(1x,i3,5g16.8,2f6.1,2x,f4.2)') &
        !        nsel, ( e(2*i)-e(2*i-1), i=1,5), wmin, wmax, sc
        
	
	write(*,"(A,i3,E17.8)")('state',i,e(i)-zpe, i=1,40,2)
        DO ie = 1,20
          WRITE(*,"(A,2i3,E15.6)") 'tunnel splitting: ',ie*2,ie*2-1, e(ie*2)-e(ie*2-1) 
        ENDDO


        !CALL SymmetrizeFunctions(tau,nq,h,nfit)
	!CALL TestSym(nq,h,h0,e)
	
	    write (45,*) nsel, nq
        write (45,*) (v(i),i=1,nq)
        write (45,*) (e(i),i=1,nq)
        write (45,*) (h(i),i=1,nq**2)
      	WRITE (445,*) nsel,0
	write (445,*) ((e(i)-zpe,1),i=1,40)
	
	
	!!!!!!!!!!!!!!Determine Prity of the wave function!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DO j=1,20,2
	sh(nsel,j) = 0.
	sh(nsel,j+1) = 0.
	 DO i=1,INT(nq/2)
	     sh(nsel,j) = sh(nsel,j) + 0.5*( h((j-1)*nq+i) - h(j*nq+1-i) )**2 
	     sh(nsel,j+1) = sh(nsel,j+1) + 0.5*( h((j)*nq+i) + h((j+1)*nq+1-i) )**2 
	 ENDDO
	
	ENDDO
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do
      write (*,'(2x,82("-")//)')

      
!      do nsel = 0, 0!nx-7
!	WRITE (*,'(a,2i4,2ES15.3)') ('E*', nsel,j, sh(nsel,j),sh(nsel,j+1) ,j=1,20,2)
!     ENDDO
     
! test for othonormality
        WRITE(*,*) 'test for orthonormality'
	do i=1,20
	     DO j = 1,i
		    dnsum = 0.0
            do q=1,nq
			   dnsum = dnsum + h((j-1)*nq+q)*h((i-1)*nq+q)
		    end do
		    IF(DABS(dnsum).gt.1D-32.AND.i.ne.j) &
     	   WRITe(*,'(a,2i4,a,ES12.3)') 'state ',i,j,' not normal: scalarproduct = ',dnsum
        END DO 
        ENDDO   
	
!    test values of eigenfunctions at tau=180
!     Do i=1,20
!      WRITE (*,'(i4,3E15.3)') i, h((i-1)*nq+(nq+1)/2-1), &
!      h((i-1)*nq+(nq+1)/2),h((i-1)*nq+(nq+1)/2+1)
!     END DO
     
     write (*,'(2x,82("-")//)')
     CALL itime(ntime)
     WRITE (*,*) "Rechenzeit:"
     WRITE (*,*) ntime-ntime0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
      
      
      END PROGRAM RPH
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE ReadData(nx,nfit,s_fit,s_energies,s_Vpseudo,s_Gcontra,s_bkl)
      
      USE structs
      TYPE(fitgeom) :: s_fit 
      TYPE(gcontra) :: s_Gcontra(4,4)
      TYPE(vpseudo) :: s_Vpseudo
      TYPE(energies) :: s_energies(nx-7)
      TYPE(bkl) :: s_bkl(nx-7,nx-7)
      INTEGER :: nx,nfit,nfitout,nfitout2,nfitout3     
      Character :: xa*3
      
      nfit = 0
      READ(1,'(i3)')  s_fit%flagv0			        ! Read in flag for potential
      IF(s_fit%flagv0.NE.1) STOP 'error reading in the potential' 
      DO j=1,maxorder+1,3					!Read in fitted potential and determine order of fit
                READ(1,*) xa
		IF(xa.NE."@") THEN; BACKSPACE (1)
		ELSE;  nfit=j-1; EXIT; ENDIF
		READ (1,'(3ES)') s_fit%fitv0(j:j+2)
      ENDDO
      WRITE (*,'(x,a)') 'potential successfully read'
      !nfitout = 2*nfit
!!!!!!!!!!!!!!!!!!!!!!!!!!Read energies!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      DO j=1,nx-7
	READ(1,'(i3)')  s_energies(j)%flag
	READ(1,'(3ES)') s_energies(j)%series(1:nfit)	
	READ(1,*)
      ENDDO
      WRITE (*,'(x,a)') 'energies successfully read'
      	
!!!!!!!!!!!!!!!!!!!!!!!!!!contra variant metric tensor!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
     Do i=1,4
	Do j=1,i
	  READ(1,'(i3)')  s_Gcontra(i,j)%flag
	  READ(1,'(3ES)') s_Gcontra(i,j)%series(1:nfit)	
	  READ(1,*)
	ENDDO
     ENDDO
      WRITE (*,'(x,a)') 'contravariant tensor G succesfully read'	
!!!!!!!!!!!!!!!!!!!!!!!!!!contra variant metric tensor!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
     Do k=1,nx-7
	Do l=1,k
	  READ(1,'(i3)')  s_bkl(k,l)%flag1
	  READ(1,'(3ES)') s_bkl(k,l)%series1(1:nfit)
	  READ(1,*)
          READ(1,'(i3)')  s_bkl(k,l)%flag2
	  READ(1,'(3ES)') s_bkl(k,l)%series2(1:nfit)
	  READ(1,*)
          READ(1,'(i3)')  s_bkl(k,l)%flag3
	  READ(1,'(3ES)') s_bkl(k,l)%series3(1:nfit)
	  READ(1,*)
	  READ(1,'(i3)')  s_bkl(k,l)%flag4
	  READ(1,'(3ES)') s_bkl(k,l)%series4(1:nfit)	  
	  READ(1,*)
	ENDDO
     ENDDO
     WRITE (*,'(x,a)') 'Bkl matrix succesfully read'	
!!!!!!!!!!!!!!!!!!!!!!!!!!Write out pseudio potential!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      READ(1,'(i3)')  s_Vpseudo%flag			        ! Read in flag for potential
      IF(s_fit%flagv0.NE.1) STOP 'error writing in the potential' 
      READ(1,'(3ES)') s_Vpseudo%series(1:nfit)	
      READ(1,*)
      WRITE (*,'(x,a)') 'pseudo potential U successfully written'

      
      WRITE(*,2001)
      WRITE(*,*)
      WRITE(*,*) 'Read fitted input '
      WRITE(*,*) 'Number of series parameters: ', nfit      
      WRITE (*,*)  
2001  format (80('!'))	  
      END SUBROUTINE ReadData
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SymmetrizeFunctions(tau,nq,h,nfit)
      
      USE STRUCTS
      USE MARQFIT
      INTEGER :: nfit,nq,flagf(maxorder)
      DOUBLE PRECISION :: tau(nq),h(nq,nq),sh(maxpoints),fitf(maxorder),par1(maxorder)

        nfit=maxorder
      	!!!!!!!!!!!!!!Determine Prity of the wave function!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DO j=1,30,2
	sh(j) = 0.
	sh(j+1) = -1.
	flagf(j)=1
	flagf(j+1)=-1
	!!!!!!!!!!!!!!!!!!determine symmetry befor fit!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 DO q=1,INT(nq/2)
	     sh(j) = sh(j) + 0.5*( h(q,j) - h(nq+1-q,j) )**2 !=0 for E*=1, =2 for E*=-1
	     sh(j+1) = sh(j+1) + 0.5*( h(q,j+1) - h(nq+1-q,j+1) )**2
	     !sh(nsch,j+1) = sh(nsch,j+1) +  4*h((j-1)*nq+i)**2
	    ! IF(svib(j+1).GT.0.0) svib(j) = svib(j)/svib(j+1)
	 ENDDO
!	 WRITE (*,'(a,i4,2ES15.3)') 'E*', j, sh(j),sh(j+1) 
	ENDDO
	
	!!!!!!!fit and expand!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	WRITE (*,*) 'fit and expand eigenfunctions'
	DO j=1,20
	  CALL SetupFit(tau,fitf(1:nfit),h(1:nq,j),nq,nfit,flagf(j),ntol,nopt,rmsr)
	  IF(rmsr.gt.0.001) WRITE (*,'(i4,5x,a,i4,17x,i4,10x,g11.3)') j, 'warning: big rmsr',NOPT,ntol,rmsr
	  DO q=1,nq
	   h(q,j)= expansion(flagf(j),tau(q),fitf(1:nfit),par1,nfit)
	  ENDDO
	ENDDO
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE TestSym(nq,l,h,e)

        !this SR checks the diagonalization E = L'HL
	
	INTEGER :: nq,j,k,m
	DOUBLE PRECISION :: l(nq,nq),h(nq,nq),e(nq),teste(nq)
	
        WRITE (*,*) 'test diagonalisation'
	DO j = 1,nq
	teste(j)=0.
	 DO k=1,nq	
	  DO m=1,nq
		teste(j) = teste(j) + l(k,j)*h(k,m)*l(m,j)	
	  ENDDO	
	 ENDDO
	ENDDO
	WRITE (*,'(2f12.3)') e(1),teste(1)
	DO j=1,20,2
		WRITE (*,'(2ES12.3)') e(2*j)-e(2*j-1),teste(2*j)-teste(2*j-1)
	ENDDO	

        
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Initp ( n, p )

!    06.06.2010
! Called by MAIN

!  Initialize momentum vector P/SQRT(-2*hc) , i.e. P(N+I-J) is the matrix
!  element (I,J) of the operator { -SQRT[h/(8*pi**2*c)] * d/dQ } where
!  h is Planck's constant, c the speed of light, and d/dQ the derivative
!  with respect to the grid index.
!  h/(8*pi**2*c) = 16.85763128 [ u * Angstrom**2 * cm**(-1) ]
!     ( u : atomic mass unit )

! Input :
!   N : number of grid points

! Output:
!   P : linear momentum operator

      implicit NONE
      INTEGER, INTENT(IN) :: n
      INTEGER :: imj
      DOUBLE PRECISION :: p(2*n-1),w,scale      
      !dimension p(1)

      if ( mod(n,2).ne.1 ) stop ' INITP: Number of points must be odd !!'

      w = dacos(-1.d0)/n
      scale =  -sqrt( 16.85763128d0 ) * w

      p(n) = 0.0
      do imj = 1 , n-1
         p(n+imj) = scale*(-1.)**imj/sin(w*imj)
         p(n-imj) = -p(n+imj)
      ENDDO	

      end SUBROUTINE Initp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine AddVib ( n, nsch, nvib, omega, v0, v,tau,nsel )

!    06.06.2010
! Called by MAIN
! Add vibrational energy in the orthogonal coordinates to the potential

! Input:
!  N     :   number of grid points
!  NSCH  :   number of orthogonal coordinates (normal modes)
!  NVIB  :   number of quanta excited in each normal mode
!  OMEGA :   harmonic wavenumbers
!  V0    :   potential energy

! Output:
!  V     :   zero point energy corrected potential
      
      !implicit real*8 (a-h,o-z)
      INTEGER, INTENT(IN) :: n,nsch,nvib(nsch)
      DOUBLE PRECISION :: v0(n),v(n),omega(nsch,n),pw(n),gw(n),vmin
	
       !v0min = v0(1)	
       do i = 1 , n
        v(i) = v0(i)
        do j = 1 , nsch
	     v(i) = v(i) + omega(j,i) * ( dble(nvib(j)) +0.5d0 )
        end do
        write (99,'(3f16.6)') v0(i),v(i),v(i)-v0(i)
      end do
      !WRITE (99,*)

      end SUBROUTINE AddVib

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine VibAng ( np, nw, nv, w, b, g, v )

! 07.06.2010
! Called by MAIN
! adiabatic Coriolis ( vibrational angular momentum ) contributions
! to the effective adiabatic potential

! Input: 
!    NP : number of grid points
!    NW : number of normal modes = 3 * (number of atoms) - 7
!    NV : number of quanta in each normal mode
!    W  : harmonic wavenumbers
!    B  : generalized Coriolis coupling matrix
!    G  : effective G-tensor (after inversion)

! Output:
!    V  : vibrational angular momentum (Coriolis) terms

      INTEGER :: i,j,k,l,ip,np,nw,nv(nw)
      DOUBLE PRECISION :: w(nw,np),b(nw,nw,4,np),g(4,4,np), v(np),sum,s,fac 
      DOUBLE PRECISION, PARAMETER :: hbar2=16.85763128d0 

!!!!!!!!!!!!!!!!!!!!test!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ip =101 
     DO i=1,10
	 DO j=1,i-1
	 !WRITE (*,'(3i3,4f9.3)') ip, i, j, b(i,j,1:4,ip)
	ENDDO
	ENDDO
      DO i=1,4
       !WRITE (*,'(i3,4f9.3)') ip,(g(i,j,ip),j=1,i)
      ENDDO	
      DO j=1,4
      ! g(4,j,ip) = g(4,j,ip)/(360./ np)
      ENDDO
      !g(4,4,ip) = g(4,4,ip)/ (360./ np)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
      do ip = 1 , np
	sum = 0.d0
        do k = 1, nw
          do l = 1, k-1
            fac = ( ( w(l,ip)/w(k,ip) + w(k,ip)/w(l,ip) ) &
                 * ( dble(nv(k))+0.5d0 )*(dble(nv(l))+0.5d0 ) - 0.5d0 )
            s = 0.0d0
            do i = 1, 4
              do j = 1 , i-1
              s = s + g(i,j,ip) * b(k,l,i,ip) * b(k,l,j,ip) 
             !PRINT *,'g zero', g(i,j,ip),g(i,j,np+1-ip),i,j
             !PRINT *,'b zero', b(k,l,j,ip)- b(k,l,j,np+1-ip),i,j    	
              end do
               s = 2.d0*s + g(i,i,ip) * b(k,l,i,ip)**dble(2.)
            end do
            sum = sum + s * fac
          end do
        end do
        v(ip) = sum * hbar2
        
      end do
      


      end SUBROUTINE VibAng
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Overlap ( h, n, nw, nvib, w )

!     06.06.2010
! Called by MAIN
! Calculated the overlap between consecutive (a)diabatic wavefunctions

! Input:
!   N     : number of grid points
!   NW    : number of orthogonal normal modes
!   NVIB  : vibr. quantum numbers
!   W     : harmonic wavenumbers

! Output:
!   H     : overlap matrix

      INTEGER, INTENT(IN) :: n,nw,nvib(nw)
      INTEGER :: k,j
      DOUBLE PRECISION :: w(nw,n),h(n,n)

      do i = 1, n
        do j = 1, i
          h(i,j) = 1.d0
          do k = 1, nw
            if ( nvib(k).eq.1 ) then   ! first excited state
              h(i,j) = h(i,j) * 2.d0 &
                    * dsqrt( dsqrt(w(k,i)) * dsqrt(w(k,j)) ) &
                    / dsqrt( 2.d0 * w(k,i) + 2.d0 * w(k,j) )
            else                       ! overlap of zero point wave function
              h(i,j) = h(i,j) &
                    * ( 2.d0 / ( w(k,i) + w(k,j) ) )**(3.d0/2.d0) &
                    * ( w(k,i) * w(k,j) )**(3.d0/4.d0)
            end if
          end do
        end do
      end do

!TEST      call matout ( n, -1, h )

      end SUBROUTINE Overlap
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      subroutine Ham1d ( n, e, h, h0, jsym, p, g, v, w, dq )

!    06.06.2010
!    Called by MAIN
!   Set up and diagonalize the effective 1D Hamiltonian

! Input:
!  N      : number of grid points
!  H      : overlap matrix (see subroutine OVERLAP)
!  P      : linear momentum operator devided by 
!           SQRT(-2*hc) : P(I,J) --> P(I-J+N)
!  V      : potential energy (pseudopotential included)
!  G      : vibrational G-matrix element

! Output:
! E      : eigenvalues
! H      : eigenvectors
! H0     : hamilton Matrix
! JSYM   : symmetry index for inversion at the origin
!                  (0,1,2=none,symmetric,antisymmetric)

! Local   :
!  W      : working array for TRED2 and TQL2

      !implicit real*8 ( a-h, o-z )
      INTEGER :: n,jsym(n),j,k
      DOUBLE PRECISION :: e(n),h(n,n),g(n),p(2*n-1),v(n),w(n),dq 
      DOUBLE PRECISION :: sum,s,h0(n,n)
      
      !!! hack output
      Do i=1,n
      WRITE (666,*) v(i),g(i)/dq**2
      ENDDO
   
      !!!!!!!!!!!!
      
      do i = 1 , n
       
	do j = 1 , i
          sum = 0.0
          do k = 1 , n
            sum = sum - p(i-k+n) * g(k)/dq**2 * p(k-j+n)
          end do
          h(i,j) = h(i,j) * sum
          h(j,i) = h(i,j)
        end do
!	 WRITE (*,*) "diag T", h(i,i), i
	h(i,i) = h(i,i) + v(i)
      end do
      
	do i = 1 , n
	 do j = 1 , i
          h0(i,j) = h(i,j)
	  h0(j,i) = h(j,i)
	 ENDDO	
	ENDDO
       
      call tred2 ( n, n, h, e, w, h )
      call  tql2 ( n, n, e, w, h, ierr )

!  find symmetry with respect to inversion at the centre
      do j = 1 , n
        s = 0.0
        do i = 1 , n
          s = s + h(i,j)*h(n+1-i,j)
        end do
        if (s.gt.0.99) then
          jsym(j) = 1
        else if (s.lt.-0.99) then
          jsym(j) = 2
        else
          jsym(j) = 0
        end if
      end do

      end SUBROUTINE Ham1d
      
