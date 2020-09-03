!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!This program reads in fitted geometries, tensor of inertia and calculates the G-Matrix, 
!!!!!!!!!normal coordinates and frequencies after projection of the force filed matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!input is given as a series of cosine/sine or polynomial (not yet included) functions!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!multiplication and derivatives can therefore be calculated analytically!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PROGRAM RPH_HAM
      
      USE STRUCTS
      USE MARQFIT
      INTEGER :: nx,nq,natoms,nx2,nfit,nfitout, flag1,flag2,flagout,nfreq,ng,nb,nf,nl      
      DOUBLE PRECISION :: amass(maxatoms), tau(maxpoints)
      DOUBLE PRECISION :: par1(maxorder),par2(maxorder),parout (2*maxorder)

      
      TYPE(fitgeom) :: s_fit 
      TYPE(gcontra) :: s_Gcontra(4,4)
      TYPE(gco) :: s_gco(4,4)
      TYPE(vpseudo) :: s_Vpseudo
      TYPE(energies) :: s_energies(3*maxatoms)
      TYPE(nc) :: s_nc(3*maxatoms)
      TYPE(bkl) :: s_bkl(3*maxatoms-7,3*maxatoms-7)

      WRITE (*,*) 'User input:'
      READ (*,*) natoms
      READ (*,*) (amass(i),i=1,natoms)
      WRITE (*,*) 'number of atoms'; WRITE (*,'(i4)') natoms
      WRITE(*,*) 'atomic masses :'; WRITE (*,'(f9.3)') (amass(i),i=1,natoms)
      IF ( natoms.lt.1 ) stop 'error: number of atoms 0'
      IF (natoms.GT.maxatoms) STOP 'error: to many atoms for calculation'
      DO i=1,natoms; IF(amass(i).LE.0.0) STOP 'error: check atomic masses'; ENDDO
      
      
      
      nx=3*natoms
      nx2=nx**2
      
      CALL ReadGeomFit(nx,nfit,s_fit)      !!!!!!!!!read in points for fitting
      
      !fitting order of i) frequencies (s_energies) ii) G-Matrix iii) B Matrix
      READ (*,*) nfreq, ng ,nl, nb
      IF (nfreq.eq.0) nfreq = Int(nfit/2)
      IF (ng.eq.0) ng = Int(nfit/2)
      IF (nl.eq.0) ng = Int(nfit)
      IF (nb.eq.0) nb =  Int(nfit/4)

      
      READ(*,*) nq
      nq = MAX(nq,nfit) 			!never use less poins than in the ab initio calculation
      
      WRITE (*,*) 'grid points used for numerical calculations: '; WRITE(*,'(i4)') nq
      WRITE(*,*)

      

      !Calculate normal coodinates, energies, diabt rotation, fit
     
      CALL CalcFreq(nx,natoms,amass,nfit,nfreq,nl,s_fit,s_energies,s_nc,nq)

      !Calculate Gmatrix, mu tensor, correction potential U
      
      CALL SetupGmat(nx,nfit,ng,s_fit,natoms,amass,s_Gcontra,s_Vpseudo,nq)
      
      !Calculate B-Matrix
      
      CALL CalcBmat(natoms,nfit,nb,s_bkl(1:nx-7,1:nx-7),s_nc,nq)
      
      nf = MAX(nfit,nfreq,ng,nb)
      CALL WriteOut(nx,nf,s_fit,s_energies,s_Vpseudo,s_Gcontra,s_bkl(1:nx-7,1:nx-7))
      CALL WriteHam(nq,nx,nf,s_fit,s_energies,s_Vpseudo,s_Gcontra,s_bkl(1:nx-7,1:nx-7))
      

      
      END PROGRAM RPH_HAM
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE ReadGeomFit (nx,nfit,s_fit)
      
      USE structs
      TYPE (fitgeom) :: s_fit
      INTEGER :: nx,nfit    
      character  :: xa*3
      
      
      READ(1,'(i3)')  s_fit%flagv0			        ! Read in flag for potential
      IF(s_fit%flagv0.NE.1) STOP 'error reading in the potential' 
      DO j=1,maxorder	+1,3		
                READ(1,*) xa
		IF(xa.NE."@") THEN; BACKSPACE (1)
		ELSE;  nfit=j-1; EXIT; ENDIF
		READ (1,'(3ES)') s_fit%fitv0(j:j+2)
      ENDDO
      
      WRITE (*,'(x,a,i4)') 'number of fitting parameters: ', nfit
      WRITE (*,'(x,a)') 'potential successfully read'
      s_fit%fitv0(1:nfit) = hartree * s_fit%fitv0(1:nfit)

      
      DO j=1,nx
	READ(1,'(i3)')  s_fit%flagx0(j)
	DO i=1,nfit,3
		READ (1,'(3ES)') s_fit%fitx0(j,i:i+2)
	ENDDO
	!s_fit%fitx0(j,1:nfit) = s_fit%fitx0(j,1:nfit)*a0
	READ (1,*)
!	WRITE (*,*) s_fit%fitx0(j,1:nfit)
      ENDDO
      WRITE (*,'(x,a)') 'geometries successfully read'	
      
      DO j=1,nx
	READ(1,'(i3)')  s_fit%flaggrad(j)
	DO i=1,nfit,3
		READ (1,'(3ES)') s_fit%fitgrad(j,i:i+2)
	ENDDO
	READ (1,*)
      ENDDO
      WRITE (*,'(x,a)') 'gradients successfully read'	
      
      DO j=1,nx
	DO k=1,j
		READ(1,'(i3)')  s_fit%flagff(j,k)
		DO i=1,nfit,3
			READ (1,'(3ES)') s_fit%fitff(j,k,i:i+2)
		ENDDO
		READ (1,*)
	ENDDO	
      ENDDO
      WRITE (*,'(x,a)') 'force field successfully read'	
      
      	DO j=1,3
	   READ(1,'(i3)') s_fit%flagrot(j)
	   DO i=1,nfit,3
		 READ (1,'(3ES)') s_fit%fitrot(j,i:i+2)
	    ENDDO
	    READ (1,*)
	ENDDO	
	DO j=1,3
	 Do i=1,j
	      READ(1,'(i3)') s_fit%flagxi(j,i)
	      DO k=1,nfit,3
		  READ (1,'(3ES)') s_fit%fitxi(j,i,k:k+2)
	      ENDDO
	      READ (1,*)
	      s_fit%flagxi(i,j)  = s_fit%flagxi(j,i)
	      s_fit%fitxi(i,j,1:nfit) = s_fit%fitxi(j,i,1:nfit)
	 ENDDO	
	ENDDO
	WRITE (*,'(x,a)') 'inertia data successfully read'	
      
      WRITE(*,2001)
      WRITE(*,*)
      WRITE(*,*) 'Read fitted input '
      WRITE(*,*) 'Number of series parameters: ', nfit      
      WRITE (*,*)  
2001  format (80('!'))	  
      END SUBROUTINE ReadGeomFit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
	
      SUBROUTINE WriteOut(nx,nfit,s_fit,s_energies,s_Vpseudo,s_Gcontra,s_bkl)
      
      USE structs
      TYPE(fitgeom) :: s_fit 
      TYPE(gcontra) :: s_Gcontra(4,4)
      TYPE(vpseudo) :: s_Vpseudo
      TYPE(energies) :: s_energies(nx)
      TYPE(bkl) :: s_bkl(nx-7,nx-7)
      INTEGER :: nx,nfit
      

!!!!!!!!!!!all energiers in cm-1, geometire in aB!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!Write out potential!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      WRITE(2,'(i3)')  s_fit%flagv0			        ! Read in flag for potential
      IF(s_fit%flagv0.NE.1) STOP 'error writing in the potential' 
      WRITE(2,'(3ES)') s_fit%fitv0(1:nfit)	
      WRITE(2,*) '@'
      WRITE (*,'(x,a)') 'potential successfully written'

!!!!!!!!!!!!!!!!!!!!!!!!!!Write out energies!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      DO j=8,nx
	WRITE(2,'(i3)')  s_energies(j)%flag
	WRITE(2,'(3ES)') s_energies(j)%series(1:nfit)	
	WRITE(2,*) '@'
      ENDDO
      WRITE (*,'(x,a)') 'energies successfully written'	
!!!!!!!!!!!!!!!!!!!!!!!!!!contra variant metric tensor!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
     Do i=1,4
	Do j=1,i
	  WRITE(2,'(i3)')  s_Gcontra(i,j)%flag
	  WRITE(2,'(3ES)') s_Gcontra(i,j)%series(1:nfit)	
	  WRITE(2,*) '@'
	ENDDO
     ENDDO
      WRITE (*,'(x,a)') 'contravariant tensor G written'	
!!!!!!!!!!!!!!!!!!!!!!!!!!b matrix elements!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
     Do k=1,nx-7
	Do l=1,k
	  WRITE(2,'(i3)')  s_bkl(k,l)%flag1
	  WRITE(2,'(3ES)') s_bkl(k,l)%series1(1:nfit)
	  WRITE(2,*) '@'
          WRITE(2,'(i3)')  s_bkl(k,l)%flag2
	  WRITE(2,'(3ES)') s_bkl(k,l)%series2(1:nfit)
	  WRITE(2,*) '@'
          WRITE(2,'(i3)')  s_bkl(k,l)%flag3
	  WRITE(2,'(3ES)') s_bkl(k,l)%series3(1:nfit)
	  WRITE(2,*) '@'
	  WRITE(2,'(i3)')  s_bkl(k,l)%flag4
	  WRITE(2,'(3ES)') s_bkl(k,l)%series4(1:nfit)	  
	  WRITE(2,*) '@'
	ENDDO
     ENDDO
     WRITE (*,'(x,a)') 'Bkl matrix written written'	
!!!!!!!!!!!!!!!!!!!!!!!!!!Write out pseudio potential!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      WRITE(2,'(i3)')  s_Vpseudo%flag			        ! Read in flag for potential
      IF(s_fit%flagv0.NE.1) STOP 'error writing in the pseudo potential' 
      WRITE(2,'(3ES)') s_Vpseudo%series(1:nfit)	
      WRITE(2,*) '@'
      
      WRITE (*,'(x,a)') 'pseudo potential U successfully written'
!!!!!!!!!!!!!!!!!!!!!!!!!!Write out rotational constants !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      DO j=1,3
	WRITE(22,'(i3)')  s_fit%flagrot(j)			        ! Read in flag for potential
	IF(s_fit%flagv0.NE.1) STOP 'error writing rotational constants' 
	WRITE(22,'(3ES)') s_fit%fitrot(j,1:nfit)	
	WRITE(22,*) '@'
      ENDDO	
      
      WRITE (*,'(x,a)') 'rotational constantes successfully written'
      END SUBROUTINE WriteOut

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE WriteHam(nq,nx,nf,s_fit,s_energies,s_Vpseudo,s_Gcontra,s_bkl)
      
      USE structs
      USE marqfit
      TYPE(fitgeom) :: s_fit 
      TYPE(gcontra) :: s_Gcontra(4,4)
      TYPE(vpseudo) :: s_Vpseudo
      TYPE(energies) :: s_energies(nx)
      TYPE(bkl) :: s_bkl(nx-7,nx-7)
      INTEGER :: nx,nf,nq
      DOUBLE PRECISION :: tau(5001),v,par1(nf),g(4,4)
      open (32,form='formatted',status='unknown')
     
      DO q=1,nq
      tau(q)=dble(q-1)*360.d0/dble(nq)+360.d0/(2.d0*dble(nq))
      !WRITE(*,*) q,tau(q)
      ENDDO
     
      write (32,*) nq, nx-7
      write (32,*)

      write (32,"(ES)") (   expansion(s_fit%flagv0,tau(q),s_fit%fitv0,par1,nf), q=1, nq )
      write (32,*)
      DO q=1,nq
        
	DO i=1,4
          Do j=1,i
	  g(i,j) =  expansion(s_Gcontra(i,j)%flag,tau(q),s_Gcontra(i,j)%series,par1,nf)!*0.52917720859**2
	  g(j,i)=g(i,j)
	  ENDDO  
        ENDDO
        
	DO i=1,3 	
          write (32,"(4ES)") ( g(i,j), j=1,4) 
	ENDDO  
	 write (32,"(4ES)") ( g(4,j), j=1,3),g(4,4)*(dble(nq)/360.)**2
      ENDDO
      
      WRITE (*,*) 'v,g ok' 
      write (32,*)
      write (32,"(ES)") (   expansion(s_Vpseudo%flag,tau(q),s_Vpseudo%series,par1,nf), q=1, nq )
      write (32,*)
      write (32,"(11ES)") ( ( expansion(s_energies(i)%flag,tau(q),s_energies(i)%series,par1,nf), i=8,nx ), q=1,nq )
      write (32,*)
      WRITE (*,*) 'vps,freq ok'
      DO i=1,nx-7
       DO j= 1,nx-7
       !WRITE (*,*) s_bkl(i,j)%flag1,s_bkl(i,j)%flag2,s_bkl(i,j)%flag3,s_bkl(i,j)%flag4 
       ENDDO
      ENDDO
      DO q=1,nq
      write (32,"(ES)") (  ( expansion(s_bkl(i,j)%flag1,tau(q),s_bkl(i,j)%series1,par1,nf) ,i=1,nx-7) ,j=1,nx-7) ! i=1, ndat*nsch*nsch*4 )
      write (32,"(ES)") (  ( expansion(s_bkl(i,j)%flag2,tau(q),s_bkl(i,j)%series2,par1,nf) ,i=1,nx-7) ,j=1,nx-7) ! i=1, ndat*nsch*nsch*4 )
      write (32,"(ES)") (  ( expansion(s_bkl(i,j)%flag3,tau(q),s_bkl(i,j)%series3,par1,nf) ,i=1,nx-7) ,j=1,nx-7) ! i=1, ndat*nsch*nsch*4 )
      write (32,"(ES)") (  ( expansion(s_bkl(i,j)%flag4,tau(q),s_bkl(i,j)%series4,par1,nf) ,i=1,nx-7) ,j=1,nx-7) ! i=1, ndat*nsch*nsch*4 )
      ENDDO
      WRITE (*,*) 'bkl ok' 
      close (32)
      
      ENd SUBROUTINE WriteHam

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SetupGmat(nx,nfit,ng,s_fit,natoms,amass,s_Gcontra,s_detg,nq)
      
      USE marqfit
      USE structs
      TYPE (fitgeom), INTENT(IN) :: s_fit
      TYPE (gco) :: s_gco(4,4)
      TYPE (gcontra) :: s_Gcontra(4,4)
      TYPE (vpseudo) :: s_detg
      INTEGER, INTENT(IN) :: nx,nfit, natoms
      INTEGER :: e0(3,3), e1(3,3),e2(3,3),e3(3,3),nfitout,nfitout2,nfitout3,flag(3)
      INTEGER :: x0flag(natoms,3),dr1flag(3),rqflag(3),rr1flag(3),ggflag(3),rr2flag(3),rr3flag(3),dr2flag(3),dr3flag(3), &
	               par1flag,par2flag,par3flag,nq
      DOUBLE PRECISION :: amass(natoms),eps,dg(nq),dgfit(nq)
      DOUBLE PRECISION :: x0(natoms,3,nfit),G(4,4),g0(4,4,nq),gg(3,2*nfit),rq(3,nfit), &
				       dr1(3,nfit),rr1(3,2*nfit),dr2(3,nfit), &
				       rr2(3,2*nfit),dr3(3,nfit),rr3(3,2*nfit),scratch(4*2)
     
      DOUBLE PRECISION :: tau(maxpoints),par1(maxorder),par2(maxorder),par3(maxorder)
      DOUBLE PRECISION :: DR(4,4,maxpoints),S,symg(4,4),normsymg(4,4),flagf
      DOUBLE PRECISION :: d(nq),u(nq),fser(MAXORDER),dser(MAXORDER),user(MAXORDER)
      
      DO i=1,natoms														!mass weight cartesian coordinates
         ij=3*(i-1)+1
      	 x0(i,1,1:nfit)   =   s_fit%fitx0(ij,1:nfit)*DSqrt(amass(i))
      	 x0flag(i,1) 	  =   s_fit%flagx0(ij)
      	 x0(i,2,1:nfit)   =   s_fit%fitx0(ij+1,1:nfit)*DSqrt(amass(i))
      	 x0flag(i,2) 	  =   s_fit%flagx0(ij+1)
      	 x0(i,3,1:nfit)   =   s_fit%fitx0(ij+2,1:nfit)*DSqrt(amass(i))
      	 x0flag(i,3) 	  =   s_fit%flagx0(ij+2)
      ENDDO
      
      
      
      DO q=1,nq
         tau(q)=(q-0.5)*360./nq
	 DO K=1,natoms
          DR(1,2,K)=  expansion(x0flag(K,3),tau(q),x0(k,3,1:nfit),par1,nfit)
          DR(1,3,K)= -expansion(x0flag(K,2),tau(q),x0(k,2,1:nfit),par1,nfit)
          DR(2,3,K)=  expansion(x0flag(K,1),tau(q),x0(k,1,1:nfit),par1,nfit)
          DR(2,1,K)= -DR(1,2,K)
          DR(3,1,K)= -DR(1,3,K)
          DR(3,2,K)= -DR(2,3,K)
          Do j=1,3
		rq(j,:) = DfDx(nfit,x0(K,j,:),x0flag(K,j),rqflag(j))  !rq = (Dx/dQ,Dy/Dq,Dz/Dq) 
		DR(j,j,K) = 0.0
                DR(j,4,K) = expansion(rqflag(j),tau(q),rq(j,1:nfit),par1,nfit)
          end do
	 ENDDO

         DO I = 1, 4
          DO J = I, 4
            S=0.
            DO K = 1, natoms
              S = S + DR(1,I,K) * DR(1,J,K)  &
                       + DR(2,I,K) * DR(2,J,K)   &
                       + DR(3,I,K) * DR(3,J,K)
             end do
            g0(i,j,q) = S
            g0(j,i,q) = S
          end do
	 end do
          DO i=1,3
	  WRITE (51,'(4F18.10)')     (g0(i,j,q)/a0**2,j=1,i)
          ENDDO
	  WRITE (51,'(4F18.10)')     (g0(4,j,q)/a0**2,j=1,4)
	ENDDO
	
	!!!!!!!!!!1Check flags for covaraint g matrix!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Do i=1,4
	 Do j=1,i
	  symg(i,j) =0.
	  normsymg(i,j)=0.0
	  Do q=1,(nq-1)/2
           IF(g0(j,i,q).NE.0.0) THEN	
	   symg(j,i) = symg(j,i) + (g0(j,i,q) - g0(j,i,nq+1-q))**2
	   normsymg(j,i) = normsymg(j,i) 	+ 4*g0(j,i,q)**2
	  ENDIF
	 ENDDO
	flagf = symg(j,i)/normsymg(j,i)
        IF (flagf.LT.0.0.OR.flagf.GT.1.001) THEN; WRITE (*,*) 'ChkSym: Error',flagf; s_gco(j,i)%flag=2
	ELSEIF (flagf.LE.0.05) THEN; s_gco(j,i)%flag=1
	ELSEIF (flagf.GE.0.95) THEN; s_gco(j,i)%flag=-1
	ELSE; flag=0
	END IF
	s_gco(i,j)%flag = s_gco(j,i)%flag
       ENDDO
      ENDDO
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      WRITE (*,'(x,a)') 'symmetry flags of the covariant tensor g'
      Do i=1,4
      WRITE (*,'(4i4)') (s_gco(i,j)%flag,j=1,i)
      ENDDO
      
      !!!!!!!!!!!!!!!!!!!!!!!fit and expand co matrix!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Do q=1,nq
      DO i=1,3
        !g0(i,i,q) = expansion(s_Gcontra(i,i)%flag,tau(q),s_Gcontra(i,i)%series(:),par2,nfitout)
      	Do j=1,i-1
      	!g0(i,j,q) = 0.! expansion(s_Gcontra(i,j)%flag,tau(q),s_Gcontra(i,j)%series(:),par2,nfitout)
      	!g0(j,i,q) = g0(i,j,q)
	ENDDO
      ENDDO
        !g0(4,1,q)=0.
	!g0(4,3,q)=0.
	!g0(1,4,q)=0.
	!g0(3,4,q)=0.
      !Do i=1,4
      !WRITE (*,'(4F16.8)')     (g0(i,j,q),j=1,i)
      !ENDDO
      ENDDO    


            
!!!!!!!!!!!!!!!!!!!!!!!Test multiplication!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   WRITE (*,*) 'test multiplication' 
!   j=1;K=1
!   rq(j,:) = DfDx(nfit,x0(K,j,:),x0flag(K,j),rqflag(j))  !rq = (Dx/dQ,Dy/Dq,Dz/Dq) 
!  
!   rr1(j,:) = mult(nfit,rq(j,:),rqflag(j),rq(j,:),rqflag(j),rr1flag(j),nfitout) 
!   DO q=1,nq
!     WRITE (*,'(3f12.6)') tau(q),expansion(rqflag(j),tau(q),rq(j,1:nfit),par1,nfit)**2 ,&
!			expansion(rr1flag(j),tau(q),rr1(j,1:nfitout),par1,nfitout)
!   ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!alternative: calculate covariant tensor analytically!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     e0 = reshape( (/1,0,0,		0,1,0,		0,0,1/), (/3,3/) )
!     e1 = reshape( (/0,0,0,		0,0,1,		0,-1,0/), (/3,3/) )
!     e2 = reshape( (/0,0,-1,	0,0,0,	        1,0,0/), (/3,3/) )
!     e3 = reshape( (/0,1,0,        -1,0,0,		0,0,0/), (/3,3/) )
!     
!     ! Calculate G tensor, Inversion-> mu tensor, fit mu tensor
!     
!      DO i=1,3
!        s_gco(i,i)%flag=s_fit%flagrot(i)
!     	s_gco(i,i)%series(1:nfit) = s_fit%fitrot(i,1:nfit)
!     	DO j=1,i-1												 	!Molecule in Main Axes frame (diagonal tensor of inertia)
!     	   s_gco(i,j)%flag = 1
!     	   s_gco(j,i)%flag = 1
!     	   s_gco(i,j)%series(:) = 0
!     	   s_gco(j,i)%series(:) = 0
!     	ENDDO
!      ENDDO
!      
!      Do i=1,natoms
!      	Do j=1,3
!		rq(j,:) = DfDx(nfit,x0(i,j,:),x0flag(i,j),rqflag(j))  !rq = (Dx/dQ,Dy/Dq,Dz/Dq) 
!        end do	  
        !end do
	!ENDDO
	!WRITE (*,*) rqflag
!	dr1 =      MATMUL (e1,x0(i,1:3,:))	      !dr1 = e1 x (x,y,z)
!	dr1flag = MATMUL (e1,(/1,2,3/))	      !ordering
!	dr2 =      MATMUL (e2,x0(i,1:3,:))	      !dr2 = e2 x (x,y,z)
!	dr2flag = MATMUL (e2,(/1,2,3/))	      !ordering
!	dr3 =      MATMUL (e3,x0(i,1:3,:))	      !dr3 = e3 x (x,y,z)
	!dr3flag = MATMUL (e3,(/1,2,3/))	      !ordering
	!WRITE (*,*) drflag
	!WRITE (*,*) dr
!	DO k=1,3
!		IF(dr1flag(k).NE.0) THEN
!			dr1flag(k) = x0flag(i,ABS(dr1flag(k)))
!		ELSE 
!			dr1flag(k)=-1
!	        ENDIF		
!	        IF(dr2flag(k).NE.0) THEN
!			dr2flag(k) = x0flag(i,ABS(dr2flag(k)))
!		ELSE 
!			dr2flag(k)=-1
!	        ENDIF	
!		IF(dr3flag(k).NE.0) THEN
!			dr3flag(k) = x0flag(i,ABS(dr3flag(k)))
!		ELSE 
!			dr3flag(k)=-1
!	        ENDIF				
!	ENDDO
!	WRITE (*,*) drflag

!	Do j=1,3
!       WRITE (*,*) j, drflag(j),rqflag(j)
!		rr1(j,:) = mult(nfit,dr1(j,:),dr1flag(j),rq(j,:),rqflag(j),rr1flag(j),nfitout) !rr1 = dr1(j)*rq(j) j=1,3
!		rr2(j,:) = mult(nfit,dr2(j,:),dr2flag(j),rq(j,:),rqflag(j),rr2flag(j),nfitout) !rr2 = dr2(j)*rq(j) j=1,3
!		rr3(j,:) = mult(nfit,dr3(j,:),dr3flag(j),rq(j,:),rqflag(j),rr3flag(j),nfitout) !rr3 = dr3(j)*rq(j) j=1,3
!		gg(j,:) = mult(nfit,rq(j,:),rqflag(j),rq(j,:),rqflag(j),ggflag(j),nfitout)     !gg = rq(j)*rq(j) j=1,3
!	ENDDO
	!WRITE (*,*) rr2(:,2)
!	IF(rr1flag(1).NE.rr1flag(2).OR.rr1flag(2).NE.rr1flag(3)) STOP 'construct g: erorr line 278'
!	IF(rr2flag(1).NE.rr2flag(2).OR.rr2flag(2).NE.rr2flag(3)) STOP 'construct g: erorr line 279'
!	IF(rr3flag(1).NE.rr3flag(2).OR.rr3flag(2).NE.rr3flag(3)) STOP 'construct g: erorr line 280'
!	IF(ggflag(1).NE.ggflag(2).OR.ggflag(2).NE.ggflag(3)) STOP 'construct g: erorr line 281'
!	s_gco(4,1)%series(:) = s_gco(4,1)%series(:) + rr1(1,:) +  rr1(2,:) + rr1(3,:) 		!sum over all atoms
!	s_gco(4,2)%series(:) = s_gco(4,2)%series(:) + rr2(1,:) +  rr2(2,:) + rr2(3,:) 		!g(4,2) = rr(1) + rr(2) + rr(3)
!	s_gco(4,3)%series(:) = s_gco(4,3)%series(:) + rr3(1,:) +  rr3(2,:) + rr3(3,:) 
!	s_gco(4,4)%series(:) = s_gco(4,4)%series(:) + gg(1,:) + gg(2,:) + gg(3,:) 		        !g(4,1) = rr(1) + rr(2) + rr(3)
!	s_gco(4,1)%flag = rr1flag(1)
!	s_gco(4,2)%flag = rr2flag(1)
!	s_gco(4,3)%flag = rr3flag(1)
!	s_gco(4,4)%flag = ggflag(1)
!	s_gco(1,4)%series(:)= s_gco(4,1)%series(:)
!	s_gco(2,4)%series(:)=s_gco(4,2)%series(:)
!	s_gco(3,4)%series(:)=s_gco(4,3)%series(:)
!	s_gco(1,4)%flag=s_gco(4,1)%flag
!	s_gco(2,4)%flag= s_gco(4,2)%flag
!	s_gco(3,4)%flag=s_gco(4,3)%flag
        !WRITE (*,*) s_gco(4,1)%series(:)
!      ENDDO
!      DO q=1,nq
!            DO i=1,4
!            G(i,i) = expansion(s_gco(i,i)%flag,tau(q),s_gco(i,i)%series(:),par2,nfitout)
!	    Do j=1,i-1
!		G(i,j) = expansion(s_gco(i,j)%flag,tau(q),s_gco(i,j)%series(:),par2,nfitout)
!		G(j,i) = G(i,j)
!		!g0(i,j,q) = G(i,j)
!		!g0(j,i,q) = G(i,j)
!	    ENDDO
!	    !g0(i,i,q) = G(i,i)
!	    WRITE (52,'(4F18.10)')     (G(i,j)/a0**2,j=1,i)
!          ENDDO
!      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	
        DO q=1,nq	
         !!!Calculate the contravariant G tensor by matrix inversion and calc of determinante of the covariant g tensor (see external SR minv)
         eps = 1.e-12
         call minv ( g0(1,1,q), 4, 4, scratch, dg(q), eps, 0, 1 )
         !
	 WRITE (481,'(f12.6)') dg(q)
      
      !CALL SetupFit(tau,s_fit%fitv0,dg(:),nx,nfit,1,ntol,nopt,rmsr)
       
       DO i=1,3
	   WRITE (48,'(4F16.8)')     (g0(i,j,q),j=1,i)
       ENDDO
       WRITE (48,'(4F16.8)') (g0(4,j,q),j=1,4)
       
      ENDDO
      
      !!!!!!!!!!1Check flags for covaraint g matrix!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
        Do i=1,4
	 Do j=1,i
	  symg(i,j) =0.
	  normsymg(i,j)=0.0
	  Do q=1,(nq-1)/2
           IF(g0(j,i,q).NE.0.0) THEN	
	   symg(j,i) = symg(j,i) + (g0(j,i,q) - g0(j,i,nq+1-q))**2
	   normsymg(j,i) = normsymg(j,i) 	+ 4*g0(j,i,q)**2
	  ENDIF
	 ENDDO
	
	flagf = symg(j,i)/normsymg(j,i)
        IF (flagf.LT.0.0.OR.flagf.GT.1.001) THEN; WRITE (*,*) 'ChkSym: Error',flagf; s_Gcontra(j,i)%flag=2
	ELSEIF (flagf.LE.0.05) THEN; s_Gcontra(j,i)%flag=1
	ELSEIF (flagf.GE.0.95) THEN; s_Gcontra(j,i)%flag=-1
	ELSE; flag=0
	END IF
	s_Gcontra(i,j)%flag = s_Gcontra(j,i)%flag
       ENDDO
      ENDDO
      
      
      WRITE (*,'(x,a)') 'symmetry flags of the contravariant tensor G'
      Do i=1,4
      WRITE (*,'(4i4)') (s_Gcontra(i,j)%flag,j=1,i)
      ENDDO
            
      !!!!!!!!!Fit the contravariant G tensor of inertia!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      nfitout=2*nfit-1
      nfitout2 = 2*nfitout-1
      Do i=1,4
	CALL SetupFit(tau,s_Gcontra(i,i)%series,g0(i,i,1:nq),nq,ng,s_Gcontra(i,i)%flag,ntol,nopt,rmsr)
	IF(rmsr.gt.0.001) WRITE (*,'(2i4,5x,a,i4,17x,i4,10x,g11.3)') i,i, 'warning: big rmsr',NOPT,ntol,rmsr
        DO j=1,i-1
         CALL SetupFit(tau,s_Gcontra(i,j)%series,g0(i,j,1:nq),nq,ng,s_Gcontra(i,j)%flag,ntol,nopt,rmsr)
	 IF(rmsr.gt.0.001) WRITE (*,'(2i4,5x,a,i4,17x,i4,10x,g11.3)') i,j,'warning: big rmsr',NOPT,ntol,rmsr
      ENDDO
      ENDDO
      
      !Test g Matrix
      !Do q=1,nq,10
      !DO i=1,4
      !  g0(i,i,q) = expansion(s_Gcontra(i,i)%flag,tau(q),s_Gcontra(i,i)%series(:),par2,nfitout)
      !	Do j=1,i-1
      !	g0(i,j,q) = expansion(s_Gcontra(i,j)%flag,tau(q),s_Gcontra(i,j)%series(:),par2,nfitout)
      !	 g0(j,i,q) = g0(i,j,q)
      !ENDDO
      !ENDDO 
      !Do i=1,4
      !WRITE (*,'(4F16.8)')     (g0(i,j,q),j=1,i)
      !ENDDO
      !ENDDO      
      
      s_detg%flag=1
      CALL SetupFit(tau,s_detg%series,dg(1:nq),nq,ng,s_detg%flag,ntol,nopt,rmsr)
      WRITE (481,'(i4)') s_detg%flag
      WRITE(481,'(3ES)') s_detg%series(1:ng)
      WRITE (*,'(a,i4)') 'symmetry flag of det(g):', 1       
      IF(rmsr.gt.0.01) WRITE (*,'(a,i4,17x,i4,10x,g11.3)') 'warning: big rmsr',NOPT,ntol,rmsr
      !Do q=1,nq
       !dgfit(q) = expansion(s_detg%flag,tau(q),s_detg%series(:),par1,ng)
       !WRITE (*,'(3F12.6)') q,  dg(q), dgfit(q)
      !ENDDO
      
      
      !Calculate the pseudo potential
      ! U=1/4 * (d/dq (G(4,4)*d(ln(det(g)/dq))+1/4*G(D ln(det(g))/dq)^2) ) = 1/4 (G*ln(detg)')' + 1/16 ln(detg)'G ln(detg)'
      
      do q = 1 , nq
        d(q) = 0.25 * log(dg(q))		!d = 0.25 * log(detg)
      end do
     
     
      fser = DfDx(ng,s_detg%series(1:ng),s_detg%flag,flag(1))      			       			! fser =(detg)' (fit)
      fser(ng+1:MAXORDER) = 0.
      
      do q = 1 , nq
        d(q) = 0.25 * expansion(flag(1),tau(q),fser(1:ng),par1,ng) / dg(q) *g0(4,4,q)		!d = 0.25 * (log(detg))'G = 0.25 (detg)'/detg * G 
	WRITE(482,'(ES)') d(q)
      end do
          
      flag(2) = flag(1)
      CALL SetupFit(tau,dser(:),d(1:nq),nq,ng,flag(2),ntol,nopt,rmsr)  						!dser = d (fit)
      WRITE (482,'(i4)') flag(2)
      WRITE(482,'(3ES)') dser(1:ng)
      
      IF(rmsr.gt.0.01) WRITE (*,'(a,i4,17x,i4,10x,g11.3)') 'd/dq (1/4*log(d))	warning: big rmsr',NOPT,ntol,rmsr
      
     DO q=1,nq
        d(q) = 0.25*d(q)/dg(q) * expansion(flag(1),tau(q),fser(1:ng),par1,ng)  		        !d =d*fser= 1/16 log(detg)' * G * log(detg)' = 1/4 (detg)'/detg * G * 1/4 (detg)'/detg 
	WRITE(483,'(ES)') d(q)
      ENDDO
      
      fser=DfDx(ng,dser,flag(2),flag(3)) 			                  				                !fser = dser'= (G fser)' = 1/4 * (log(detg)'*G)' =  1/4 * (detg'/detg * G )'
     
      
      IF (flag(1)*flag(2).NE.flag(3)) STOP 'pseudo potential mismatch' 
      Do q=1,nq
         d(q) = d(q) + expansion(flag(3),tau(q),fser(1:ng),par1,ng) 						!d = 1/16 log(detg)' * G * log(detg)' +  1/4 * (G *log(detg)')'
         WRITE (480,'(f9.3,f12.6)') d(q)
      ENDDO

      
      IF (flag(3).NE.s_detg%flag) STOP 'pseudo potential mismatch'  
      !_detg%flag=s_Gcontra(4,4)%flag
      WRITE (*,'(a,i4)') 'symmetry flag of pseudo potential:', s_detg%flag
      CALL SetupFit(tau,s_detg%series(:),d(1:nq),nq,ng,s_detg%flag,ntol,nopt,rmsr)  			       !fit pseudo potential d and overwrite determinante of g  
      IF(rmsr.gt.0.01) WRITE (*,'(a,i4,17x,i4,10x,g11.3)') 'pseudo potential warning: big rmsr',NOPT,ntol,rmsr
     
      END SUBROUTINE SetupGmat
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE CalcFreq(nx,natoms,amass,nfit,nfreq,nl,s_fit,s_energies,s_nc,nq)
      
      
      USE marqfit
      USE structs
      INTEGER, INTENT(IN) :: nx,natoms,nfit,nfreq,nl
      DOUBLE PRECISION, INTENT(IN) :: amass(natoms)
      TYPE (fitgeom) :: s_fit
      TYPE (energies) :: s_energies(nx)
      TYPE (nc) :: s_nc(nx)
      DOUBLE PREcISION :: gxseries(nx,nfit),par1(2*nfit),p(nx,nx),angle,fp(nx,nx),sum,symff(nx,nx),normff(nx,nx),flagf
      DOUBLE PRECISION :: x0(nx,nq),gx(nx,nq),ff(nx,nx,nq),xi(9,nq),tau(nq),c0(nx,7),w0(nx,nq),rot(3,nq),rotsym(3)
      DOUBLE PRECISION :: gnc(nx**2*nq),v0(nq),s(nx**2),cn(nx**2),b(3*natoms-7,3*natoms-7,4,nq),symv !for diabatic rotations
      INTEGER :: ij,gxflag(nx),mi,mj,ic,jj,nfitout,nq
      INTEGER :: i0, iref, nsel, isel(nx) !for diabatic rotations
      DOUBLE PRECISION, PARAMETER :: pi   = acos(-1.d0),bohr = 5.29177249d-1, &
							   h = 6.6260755d0, amu  = 1.6605402d0, cvac = 2.99792458d0
      
      !mass weight gradient, forces and forcefield
      nx2 = nx**2
      rotsym(1:3)=0.
      
      DO i=1,nx
        gxseries(i,1:nfit) = DfDx(nfit,s_fit%fitx0(i,1:nfit),s_fit%flagx0(i),gxflag(i))
      ENDDO
      
      Do q=1,nq
      !create coordinates etc on the grid and mass weight!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      tau(q)=dble(q-1)*360.d0/dble(nq)+360.d0/(2.d0*dble(nq))
      v0(q) = expansion(s_fit%flagv0,tau(q),s_fit%fitv0,par1,nfit)
      
      
      Do i=1,nx
         mi=INT( (i-1)/3)+1
         x0(i,q)= expansion(s_fit%flagx0(i),tau(q),s_fit%fitx0(i,1:nfit),par1,nfit)!*DSqrt(amass(mi))	already mass weighted!
      ENDDO
      
      !WRITE (*,'(3f9.6)')   (x0(i,q),i=1,nx) 
      !WRITE (*,'(3f9.6)')   (xi(i,ip)/0.529177249,i=1,9)

      Do i=1,nx
         mi=INT( (i-1)/3)+1
         gx(i,q)= expansion(gxflag(i),tau(q),gxseries(i,1:nfit),par1,nfit)*DSqrt(amass(mi))
      ENDDO

      DO i=1,nx
      DO j=1,i
         mi = INT((i-1)/3)+1
         mj = INT((j-1)/3)+1	 
         ff(i,j,q) = expansion(s_fit%flagff(i,j),tau(q),s_fit%fitff(i,j,1:nfit),par1,nfit)/DSqrt(amass(mi)*amass(mj))
	 ff(j,i,q)=ff(i,j,q)
      ENDDO
      ENDDO
      
      DO i=1,3
       xi(3*(i-1)+1,q) = expansion(s_fit%flagxi(i,1),tau(q),s_fit%fitxi(i,1,1:nfit),par1,nfit)
       xi(3*(i-1)+2,q) = expansion(s_fit%flagxi(i,2),tau(q),s_fit%fitxi(i,2,1:nfit),par1,nfit)
       xi(3*(i-1)+3,q) = expansion(s_fit%flagxi(i,3),tau(q),s_fit%fitxi(i,3,1:nfit),par1,nfit)
       rot(i,q) = expansion(s_fit%flagrot(i),tau(q),s_fit%fitrot(i,1:nfit),par1,nfit)
      ENDDO
      
      
      !generate mass weighted displacement and rotation vectors to project out (trarot) in C0(*,1:6)
      
      !WRITE (*,*) 'x0'
      !WRITE(*,'(3f12.6)') x0(1:12,q)/0.529

      call trarot ( natoms, amass, x0(1,q), c0, xi(1,q) )
      
      !WRITE (*,*) 'cp'
      !WRITE (*,'(7f9.3)') c0
      
      !project out from mass weighted gradient wrt RP (nrmvec) and normalzie
      
      call nrmvec ( nx, c0, gx(1,q), c0(1,7), p, angle )
      
      !Project forecefield  
     
      call project ( nx, 7, c0, ff(1,1,q), fp, p )
       
       !Diagonalize forcefield
      
      call freq ( nx, fp, ff(1,1,q), w0(1,q), p )
        do ic = 1 , 6
          do i = 1 , nx
            ff(i,ic,q) = c0(i,ic)        
          end do
        end do
      
      !!!!!!!!!!!!!!! normalised mass-weighted Cartesian displacements along the path in FF(*,7,ip)
        sum = 0.0
        do i = 1 , nx
          sum = sum + gx(i,q)**2        
        end do
        sum = 1.d0/dsqrt(sum)
        do i = 1 , nx
          ff(i,7,q) = gx(i,q)*sum
        end do
	
	
	do i = 1 , 3
          w0(i,q)   = 0.0
          w0(3+i,q) =  1.d3*h/(8.*pi**2.*cvac*amu)/rot(i,q)!** !!= expansion(s_fit%flagxi(i,1),tau(q),s_fit%fitxi(i,1,1:nfit),par1,nfit)
        end do
        w0(7,q) = 0.0
        
      ENDDO
      
      !!!check symmetries
      DO q=1,nq 
	symv = symv + 0.5*(v0(q) - v0(nq+1-q))**2
	DO i=1,3
	 rotsym(i) = rotsym(i) + 0.5*(rot(i,q)-rot(i,nq+1-q))**2
	ENDDO
      ENDDO
      !WRITE (*,*) symv
      IF(DABS(rotsym(i)).GE.0.01) THEN
	WRITE (*,*) 'error with symmetry of rot constants!'
	WRITE (*,'(3e15.3)') (rotsym(j),j=1,3)
      ENDIF	
      
      
      Do q=1,nq
      write (*,'(i4,2x,f7.2,2x,3f7.3,2x,6f10.1)') &
	         q,angle,(w0(i,q),i=4,6),(w0(i,q),i=8,min(13,nx))
        if ( nx.gt.12 ) write (*,'(38x,6f10.1)') (w0(i,q),i=14,nx)
       ENDDO	
     
      WRITE (*,'(x,a)') 'fit rotational constants to pi periodic cosine series'
      Do i=1,3
        s_fit%flagrot(i)=1
        s_fit%fitrot(i,1:nfreq) = 0.
	CALL SetupFit(tau,s_fit%fitrot(i,1:nfreq),w0(3+i,1:nq),nq,nfreq,s_fit%flagrot(i),ntol,nopt,rmsr)	 ! fit energies to a pi/2 periodic cosine series
	IF(rmsr.gt.0.01) WRITE (*,'(i4,5x,a,i4,17x,i4,10x,g11.3)') 3+i,'warning: big rmsr',NOPT,ntol,rmsr
      ENDDO
            
      Do q=1,nq	!write out energies befor diabatic rotations
              write (7,'(i4,8f12.3)') q, ( w0(i,q), i=8,min(nx,15) )     
              if (nx.gt.15) write (7,'(4x,8f12.3)') ( w0(i,q), i=16,nx )
	      WRITE (40,'(3ES)') (w0(i,q),i=4,6)
      ENDDO
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Check phase & phase match of normal coordinates!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Do q=2,nq
	   Do k=8,nx
	   sum=0.
		DO i=1,3*natoms
			sum=sum+ff(i,k,q)*ff(i,k,q-1) !check overlap between geometries sum_i=1,3N l_i_k (q) * l_i_k (q-1) (k-th normal coordinate vector)
		ENDDO
		IF (sum.LT.0.0) THEN !phase change 1 -> -1
		 Do i=1,3*natoms
		     ff(i,k,q) = -ff(i,k,q)
		 ENDDO
		ENDIF 
	   ENDDO 
      ENDDO
     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!diabtatic roation routines!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
	
! potential minimum as reference point for diabatic rotations if not specified
! otherwise in the input
! GNC temporary used for harmonic wavenumbers (non-diagonal !)

      do q = 1 ,nq
        i0 = (q-1)*nx2
        do i = 1 , nx2
          gnc(i0+i) = 0.0
        end do
        do i = 1 , nx
          gnc(i0+i+nx*(i-1)) = w0(i,q)
        end do
      end do
 11 read (*,*,END=19) iref, nsel, (isel(i),i=1,nsel)
    ! IF (nsel.GT.0) READ (*,*) (isel(i),i=1,nsel)
     if ( iref.lt.1 .or. iref.gt.nq ) then
        iref = 1
        do q = 2 , nq
         if ( v0(q).lt.v0(iref) ) iref = q
        end do
      end if
      call diabat ( nx, nq, iref, nsel, isel, ff(1:nx,1:nx,1:nq), gnc, s, p, cn )!!!!!!!!!!!!!!!!!!!!!!!call external function diabat see diabat.f
      go to 11
   19 continue
      do q = 1 , nq
        i0 = (q-1)*nx2
        do i = 1 , nx
          w0(i,q) = gnc(i0+i+nx*(i-1))
        end do
      end do
      
      Do q=1,nq						!WRITE out energies AFTER diabat rot																															
              write (8,'(i4,8f12.3)') q, ( w0(i,q), i=8,min(nx,15) )     
              if (nx.gt.15) write (8,'(4x,8f12.3)') ( w0(i,q), i=16,nx )
      ENDDO
      
  
     
 
      WRITE (*,'(x,a)') 'Chk symmetry flags for normal coordinates:'
      WRITE (*,'(x,a)') 'row index: i-th atom (x,y,z), col: k-kth normal coordinate (k=8,3N)'
      WRITe (*,'(2x,a,40i3)') 'i', (k,k=8,nx) 
      
      Do i=1,3*natoms									!nuclear coordinate
	 Do j=8,nx							        !k-th normal coordinate (1-3 translation,3-6 rotational,  7 reaction path, 8-3N vibrational)
	  symff(i,j) =0.0
	  normff(i,j)=0.0
	  Do q=1,(nq-1)/2
           IF( ff(i,j,q).NE.0.0 ) THEN	
	   symff(i,j) = symff(i,j) + (ff(i,j,q) - ff(i,j,nq+1-q))**2
	   normff(i,j) = normff(i,j) 	+ 4*ff(i,j,q)**2
	  ENDIF
	 ENDDO
	 flagf = symff(i,j)/normff(i,j)
         IF (flagf.LT.0.0.OR.flagf.GT.1.001) THEN; WRITE (*,*) 'ChkSym: Error',flagf; s_nc(i)%flag(j)=2
	ELSEIF (flagf.LE.0.05) THEN; s_nc(i)%flag(j)=1
	ELSEIF (flagf.GE.0.95) THEN; s_nc(i)%flag(j)=-1
	ELSE; flag=0;WRITE (*,*) flagf
	END IF
	
     ENDDO
       WRITE (*,'(40i3)') i,(s_nc(i)%flag(j),j=8,nx)
      ENDDO
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!Fit normal coordinates!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      WRITE (*,'(x,a)') 'fit normal coordinates to pi periodic (co)sine series'
      nfitout = nl
      Do k=8,nx
      Do i=1,3*natoms
        !WRITE (*,'(3i8)') i,k,s_nc(i)%flag(k)
        CALL SetupFit(tau,s_nc(i)%series(k,1:nfitout), ff(i,k,1:nq),nq,nfitout,s_nc(i)%flag(k),ntol,nopt,rmsr) !l_i_j(q) with i index of i-th coordinate and j-th normal coordinate 
        IF(rmsr.gt.0.01) WRITE (*,'(2i4,5x,a,i4,17x,i4,10x,g11.3)') i,k,'warning: big rmsr',NOPT,ntol,rmsr

!!!!!!!!!!!!!!!!!!!!!!!1 calculate the derivatives of normal coordinates
	s_nc(i)%gncseries(k,1:nfitout) = DfDx(nfitout,s_nc(i)%series(k,1:nfitout),s_nc(i)%flag(k),s_nc(i)%flaggnc(k))
       ENDDO	
      ENDDO
      
      !!!Test!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      k=8
      Do q=1,nq
       WRITE (78,'(20F15.8)') (ff(i,k,q),expansion(s_nc(i)%flag(k),tau(q),s_nc(i)%series(k,1:nfitout),par1,nfitout), i=1,10)
      ENDDO
!      DO k=8,nx
!        DO j=1,8
!	     WRITE (78,'(i4)') s_nc(
!      ENDDO
!      ENDDO      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Do q=1,nq
       WRITE (49,*) q
       Do i=1,3*natoms
	WRITE (49,'(40F15.8)') (ff(i,j,q),j=8,nx)!, ff(2,10,q),ff(3,10,q) !write out unfitted normal coordinates
       ENDDO	
       WRITE (49,*)	
      ENDDO
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      Do q=1,nq
	WRITE (77,'(40f15.7)') ((ff(i,j,q)-expansion(s_nc(i)%flag(j),tau(q),s_nc(i)%series(j,1:nfitout),par1,nfitout),j=8,nx),i=1,10)!, ff(2,10,q),ff(3,10,q) !write out unfitted normal coordinates
      ENDDO
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Fit energies!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      WRITE (*,'(x,a)') 'fit energies to pi periodic cosine series'
      Do i=1,nx
        s_energies(i)%flag = 1
        s_energies(i)%series(:) = 0.
	CALL SetupFit(tau,s_energies(i)%series(1:nfreq),w0(i,1:nq),nq,nfreq,s_energies(i)%flag,ntol,nopt,rmsr)	 ! fit energies to a pi/2 periodic cosine series
	IF(rmsr.gt.0.01) WRITE (*,'(i4,5x,a,i4,17x,i4,10x,g11.3)') i-7,'warning: big rmsr',NOPT,ntol,rmsr
	s_energies(i)%flag=1
	Do q=1,nq
		w0(i,q) = expansion(s_energies(i)%flag,tau(q),s_energies(i)%series(1:nfit),par1,nfit)
	ENDDO
      ENDDO
      
      END SUBROUTINE CalcFreq
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

     SUBROUTINE CalcBmat(natoms,nfit,nb,s_bkl,s_nc,nq)
     
     USE Structs
     USE marqfit
     INTEGER, INTENT(IN) :: natoms,nfit,nb
     INTEGER :: nfitout,nfitout2,jj,kk,testflag(4,3*natoms),flag(4)
     DOUBLE PRECISION :: par1(MAXORDER),par2(MAXORDER),par3(MAXORDER),par4(MAXORDER),tau(nq)
     TYPE (nc), INTENT(IN) :: s_nc(3*natoms)
     TYPE (bkl) :: s_bkl(3*natoms-7,3*natoms-7)
     DOUBLE PRECISION :: b(3*natoms-7,3*natoms-7,4,nq)
     
     nfitout = 2*nfit-1
     nfitout2 = 2*nfitout-1
     nx=3*natoms
     WRITE(*,'(a)') 'calculate and fit Bmat'
     
          b = 0.
          Do j=8,nx
	    jj = j-7
            DO k=8,j			!B_kk = 0 (Lik x Lik = 0)
	       kk=k-7
	       testflag(1:4,k) = 0
	       Do i=1,natoms  
	          DO q=1,nq
		  tau(q) = 360./nq*(q-0.5)
	          ij =3*(i-1)
		  b(jj,kk,1,q) = b(jj,kk,1,q) + expansion(s_nc(ij+2)%flag(j),tau(q),s_nc(ij+2)%series(j,1:nfitout),par1,nfitout) & 
			                                    *expansion(s_nc(ij+3)%flag(k),tau(q),s_nc(ij+3)%series(k,:),par1,nfitout) &
		     				          - expansion(s_nc(ij+3)%flag(j),tau(q),s_nc(ij+3)%series(j,:),par1,nfitout)  &
							  *expansion(s_nc(ij+2)%flag(k),tau(q),s_nc(ij+2)%series(k,:),par1,nfitout)
                   b(jj,kk,2,q) = b(jj,kk,2,q) - expansion(s_nc(ij+1)%flag(j),tau(q),s_nc(ij+1)%series(j,:),par1,nfitout) &
			                                   *expansion(s_nc(ij+3)%flag(k),tau(q),s_nc(ij+3)%series(k,:),par1,nfitout) &
		     				        + expansion(s_nc(ij+3)%flag(j),tau(q),s_nc(ij+3)%series(j,:),par1,nfitout)  &
							  *expansion(s_nc(ij+1)%flag(k),tau(q),s_nc(ij+1)%series(k,:),par1,nfitout)
		   b(jj,kk,3,q) = b(jj,kk,3,q) + expansion(s_nc(ij+1)%flag(j),tau(q),s_nc(ij+1)%series(j,:),par1,nfitout) &
			                                   *expansion(s_nc(ij+2)%flag(k),tau(q),s_nc(ij+2)%series(k,:),par1,nfitout) &
		     				        - expansion(s_nc(ij+2)%flag(j),tau(q),s_nc(ij+2)%series(j,:),par1,nfitout)  &
							  *expansion(s_nc(ij+1)%flag(k),tau(q),s_nc(ij+1)%series(k,:),par1,nfitout)
!		   bkl = sum_i Lik'. Lil								  
		   b(jj,kk,4,q) = b(jj,kk,4,q) + expansion(s_nc(ij+1)%flaggnc(j),tau(q),s_nc(ij+1)%gncseries(j,:),par1,nfitout)  &
							   *expansion(s_nc(ij+1)%flag(k),tau(q),s_nc(ij+1)%series(k,:),par1,nfitout) 		&
		     				          +expansion(s_nc(ij+2)%flaggnc(j),tau(q),s_nc(ij+2)%gncseries(j,:),par1,nfitout)   &
			                                   *expansion(s_nc(ij+2)%flag(k),tau(q),s_nc(ij+2)%series(k,:),par1,nfitout)		&
		     				          +expansion(s_nc(ij+3)%flaggnc(j),tau(q),s_nc(ij+3)%gncseries(j,:),par1,nfitout) 	&
			                                   *expansion(s_nc(ij+3)%flag(k),tau(q),s_nc(ij+3)%series(k,:),par1,nfitout) 							   
		 ENDDO
		 !!!!!flags are tested by summing them up and dividing by the number of atoms. If flags != +1 or -1 something is wrong
		 testflag(1,k) = testflag(1,k) + (s_nc(ij+2)%flag(j)* s_nc(ij+3)%flag(k)+ &
							        s_nc(ij+2)%flag(k)* s_nc(ij+3)%flag(j))/2
		 testflag(2,k) = testflag(2,k) + (s_nc(ij+1)%flag(j)* s_nc(ij+3)%flag(k)+ &
								s_nc(ij+1)%flag(k)* s_nc(ij+3)%flag(j))/2
		 testflag(3,k) = testflag(3,k) + (s_nc(ij+1)%flag(j)* s_nc(ij+2)%flag(k)+ &
								s_nc(ij+1)%flag(k)* s_nc(ij+2)%flag(j))/2
		 testflag(4,k) = testflag(4,k) + (s_nc(ij+1)%flaggnc(j)*s_nc(ij+1)%flag(k) + &
								s_nc(ij+2)%flaggnc(j)*s_nc(ij+2)%flag(k) + &
								s_nc(ij+3)%flaggnc(j)*s_nc(ij+3)%flag(k))/3
	       ENDDO 	 
	       testflag(1:4,k)=testflag(1:4,k)/natoms
	       Do i=1,4; IF(ABS(testflag(i,k)).NE.1) STOP 'Bmat: symmetry error';ENDDO
	       s_bkl(jj,kk)%flag1 = testflag(1,k)
	       s_bkl(jj,kk)%flag2 = testflag(2,k)
	       s_bkl(jj,kk)%flag3 = testflag(3,k)
	       s_bkl(jj,kk)%flag4 = testflag(4,k)
	       s_bkl(kk,jj)%flag1 = s_bkl(jj,kk) %flag1
	       s_bkl(kk,jj)%flag2 = s_bkl(jj,kk)%flag2 
	       s_bkl(kk,jj)%flag3 = s_bkl(jj,kk)%flag3
	       s_bkl(kk,jj)%flag4 = s_bkl(jj,kk)%flag4
	       CALL SetupFit(tau,s_bkl(jj,kk)%series1(:),b(jj,kk,1,1:nq),nq,nb,s_bkl(jj,kk)%flag1,ntol,nopt,rmsr)	 ! fit energies to a pi/2 periodic cosine series
	       IF(rmsr.gt.0.01) WRITE (*,'(i4,5x,2i4,5x,a,i4,17x,i4,10x,g11.3)') 1, jj, kk, 'warning: big rmsr',NOPT,ntol,rmsr
	       CALL SetupFit(tau,s_bkl(jj,kk)%series2(:),b(jj,kk,2,1:nq),nq,nb,s_bkl(jj,kk)%flag2,ntol,nopt,rmsr)	 ! fit energies to a pi/2 periodic cosine series
	       IF(rmsr.gt.0.01) WRITE (*,'(i4,5x,2i4,5x,a,i4,17x,i4,10x,g11.3)') 2, jj, kk, 'warning: big rmsr',NOPT,ntol,rmsr
	       CALL SetupFit(tau,s_bkl(jj,kk)%series3(:),b(jj,kk,3,1:nq),nq,nb,s_bkl(jj,kk)%flag3,ntol,nopt,rmsr)	 ! fit energies to a pi/2 periodic cosine series
	       IF(rmsr.gt.0.01) WRITE (*,'(i4,5x,2i4,5x,a,i4,17x,i4,10x,g11.3)') 3, jj, kk, 'warning: big rmsr',NOPT,ntol,rmsr
	       CALL SetupFit(tau,s_bkl(jj,kk)%series4(:),b(jj,kk,4,1:nq),nq,nb,s_bkl(jj,kk)%flag4,ntol,nopt,rmsr)	 ! fit energies to a pi/2 periodic cosine series
	       IF(rmsr.gt.0.01) WRITE (*,'(i4,5x,2i4,5x,a,i4,17x,i4,10x,g11.3)') 4, jj, kk, 'warning: big rmsr',NOPT,ntol,rmsr
	       s_bkl(kk,jj)%series1(:) = s_bkl(jj,kk)%series1(:)	       
	       s_bkl(kk,jj)%series2(:) = s_bkl(jj,kk)%series2(:)
	       s_bkl(kk,jj)%series3(:) = s_bkl(jj,kk)%series3(:)
	       s_bkl(kk,jj)%series4(:) = s_bkl(jj,kk)%series4(:)
	    ENDDO
	    !WRITe (79,'(i4,x,200i4)') j,(testflag(1:4,k),k=8,j)
          ENDDO
     DO q=1,nq
       WRITE (50,'(i4)') q
       DO i=1,4
         Do k=1,nx-7
           WRITE  (50,'(12f24.18)') (b(k,l,i,q),l=1,k)	
         ENDDO
       ENDDO
     ENDDO

     
     END SUBROUTINE CalcBmat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      subroutine trarot ( nat, amass, x0, c0, xi )

!   31.5.2010 adapted from version of 14.06.2000 by D. Luckhaus
! Generate rotational and translational mass-weighted displacement vectors
! in the  cartesian reference coordinate system, which is assumed to be a
! principle axis system.

! Input:
!   NAT : number of atoms
!   AMASS : atomic masses
!   X0 : cartesian reference geometry
!   XI : inverse square root of the tensor of inertia

! Output:
!   C0 : mass-weighted translational and rotational displacement vectors

         
      INTEGER :: iat,i,j,k,nat
      DOUBLE PRECISION, INTENT(IN) :: amass(nat), xi(3,3), x0(3,nat) 
      DOUBLE PRECISION :: s, c0(3,nat,6),sm,x1,x2,x3  
      
      s=amass(1)
      do iat = 2 , nat
        s = s + amass(iat)
      end do
      do j = 1 , 6
        do iat = 1 , nat
          do i = 1 , 3
            c0(i,iat,j) = 0.0
          end do
        end do
      end do

! translational vectors
      do j = 1 , 3
        do iat = 1 , nat
          c0(j,iat,j) = Dsqrt(amass(iat)/s)
        end do
      end do

! rotational vectors
      do iat = 1 , nat
        sm = Dsqrt(amass(iat))
        x1 = sm * x0(1,iat) 
        x2 = sm * x0(2,iat)
        x3 = sm * x0(3,iat) 
	!WRITE (*,*) x1, x2, x3
        do l = 1 , 3
          c0(1,iat,3+l) =  xi(l,2) * x3 - xi(l,3) * x2
          c0(2,iat,3+l) =  xi(l,3) * x1 - xi(l,1) * x3
          c0(3,iat,3+l) =  xi(l,1) * x2 - xi(l,2) * x1
        end do
      end do

!TEST
!      WRITE (*,*) 'TRAROT: XI'
!      Do i=1,3
!      WRITE (*,'(3f16.6)') (xi(i,j),j=1,3)
!      ENDDO
!      write (*,*) ' TRAROT: X0'
!      write (*,'(3f16.6)') x0
!      write (*,*) ' TRAROT: C0'
!      do iat = 1 , nat
!        do i = 1 , 3
!          write (*,'(2i3,6f12.6)') iat,i,(c0(i,iat,k),k=1,6)
!        end do
!      end do

      end SUBROUTINE trarot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine nrmvec ( nx, c, f0, f, p, angle )

!   31.5.2010 adapted from nrmvec 06.06.2000
! Project rotational and translational components C out of a vector F0
! and normalize it to produce F.

!* Input:
!*  NX : number of cartesian coordinates
!*  C  : vectors to project out of F0
!*  F0 : input vector

! Output:
!*  F     : projected vector
!*  ANGLE : projection angle in degrees

!* Local:
!*  P  : projection matrix

      INTEGER :: nx,i,j,k
      INTEGER, PARAMETER :: nc=6
      DOUBLE PRECISION :: c(nx,nc), f0(nx), f(nx), p(nx,nx)
      DOUBLE PRECISION :: sf0,sf,angle

!TEST 
!*      do k = 1 , nc
!*        s = 0
!*        s1 = 0
!*        s2 = 0
!*        do i = 1 , nx
!*          s1 = s1 + c(i,k)**2
!*          s2 = s2 + f0(i)**2
!*          s = s + c(i,k)*f0(i)
!*        end do
!*        write (*,*) k,s/sqrt(s1*s2)
!*      end do
!      WRITE (*,*) 'nrmvec c:'
!      Do i=1,nx
!      WRITE (*,'(6f11.5)') c(i,1:6)
!      ENDDO
      do i = 1 , nx
        do j = 1 , i
          p(i,j) = 0.0
          if (j.eq.i) p(i,j) = 1.0
          do k = 1 , nc
            p(i,j) = p(i,j) - c(i,k)*c(j,k)
          end do
          p(j,i) = p(i,j)
        end do
      end do

      sf0 = 0.0        ! norm of original vector
      sf  = 0.0        ! norm of projected vector
      do i = 1 , nx
        f(i) = 0.0
        do k = 1 , nx
          f(i) = f(i) + p(i,k)*f0(k)
        end do
        sf0 = sf0 + f0(i)**2
        sf = sf + f(i)**2
      end do
!      WRITE (*,*) 'nrmvec f0:'
!      WRITE (*,'(3F11.5)') f0
!      WRITE (*,*) 'nrmvec f:'
!      WRITE (*,'(3F11.5)') f
!* normalize projected vector
      angle = 0.0               ! projection angle
      if ( sf.gt.0.0 ) then
        sf0 = 1.d0 / sqrt( sf0 )
        sf  = 1.d0 / sqrt( sf )
        do i = 1 , nx
          f(i) = f(i) * sf
          angle = angle + f(i)*f0(i)*sf0
        end do
        if ( angle.lt.-1.d0 ) then
          angle = -180.0
        else if ( angle.gt.1.d0 ) then
          angle = 180.0
        else
          angle = acos( angle ) / acos( -1.d0 ) * 180.d0
        end if
      end if
      
      end SUBROUTINE nrmvec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine project ( nx, nc, c, f0, f, p )
   
! Project the NC coordinates C out of F0 :  F = [1-C*C(t)]*F0*[1-C*C(t)]

! Input:
!*  NX : number of cartesian coordinates
!*  NC : number of vectors C
!*  C  : vectors to project out of F0
!*  F0 : cartesian force constant matrix

!* Output:
!*  F  : projected force constant matrix

!* Local:
!*  P  : projection matrix

      INTEGER :: nx,i,j,k,l
      DOUBLE PRECISION ::  c(nx,nc), f0(nx,nx), f(nx,nx), p(nx,nx)

      do i = 1 , nx
        do j = 1 , i
          p(i,j) = 0.0
          if (j.eq.i) p(i,j) = 1.0
          do k = 1 , nc
            p(i,j) = p(i,j) - c(i,k)*c(j,k)
          end do
          p(j,i) = p(i,j)
        end do
      end do

      do i = 1 , nx
        do j = 1 , i
          f(i,j) = 0.0
          do k = 1 , nx
            do l = 1 , nx
              f(i,j) = f(i,j) + p(i,k)*f0(k,l)*p(l,j)
            end do
          end do
          f(j,i) = f(i,j)
        end do
      end do

      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11


      subroutine freq ( nx, f, c, e, w )

!* diagonalize the force constant matrix

!* Input:
!*  NX : dimension of F
!*  F  : mass-weighted cartesian force constant matrix [a.u. = Hartree/(u Bohr)]

!* Output:
!*  E : harmanic wavenumbers in [cm^-1]
!*  C : corresponding eigenvectors

!* Local:
!*  W : used by TRED2/TQL2
      INTEGER :: nx,j
      DOUBLE PRECISION ::  c(nx,*), e(nx), w(nx), f(nx,nx),ierr,scale
      DOUBLE PRECISION, PARAMETER :: pi = acos(-1.d0), bohr= 5.29177249d-1,rydberg = 1.0973731534d0, &
							   h  = 6.6260755d0, 									 &! Planck's constant (mantissa only)
							   hartree = 2.d0*rydberg*h,amu = 1.6605402d0,          			 &! atomic mass unit (mantissa only)
							   cvac = 2.99792458d0         								   ! vacuum speed of light (mantissa only)


      call tred2 ( nx, nx, f, e, w, c )
      call  tql2 ( nx, nx, e, w, c, ierr )
      
      
      
      scale = sqrt(hartree/(amu*bohr**2*cvac))*1.d4/pi*.5d0
      do j = 1 , nx
        e(j) = sqrt(abs(e(j))) * scale

!TEST
!        write (*,*) j, e(j)

      end do

      end SUBROUTINE freq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!