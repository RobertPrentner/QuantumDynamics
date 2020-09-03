      MODULE structs
      
      IMPLICIT NONE
      INTEGER, PARAMETER :: maxorder=84,maxpoints=2048, maxatoms=64, &
			                maxdimensions=3*maxatoms
      INTEGER, PARAMETER :: inp=1,iout=2,std=6,sym=11,pout=7,intout=8 ! input/output units
      INTEGER :: i,j,k,l,m,n,q					       ! counting vaiables, q for recation path coordinate
      INTEGER :: nopt,ntol						       !later used for fitting
      DOUBLE PRECISION :: rmsr	
      
      DOUBLE PRECISION, PARAMETER :: eh2cm=21947.46313705,hartree=2.d0*1.0973731534d5, &
      								 a0=0.529177249d0
      
      TYPE geom
           DOUBLE PRECISION :: x0(maxdimensions)
           DOUBLE PRECISION :: v0
           DOUBLE PRECISION :: grad(maxdimensions)
           DOUBLE PRECISION :: ff(maxdimensions,maxdimensions)
      END TYPE
      
      TYPE fitflag										!determine flag for fits: 
          INTEGER :: flagx0(maxdimensions)				        !flag = 0 <> cos series
          INTEGER :: flagv0								!flag = 1 <> sin series
          INTEGER :: flaggrad(maxdimensions)				!0 < flag < 1 fourier series
          INTEGER :: flagff(maxdimensions,maxdimensions)
      END TYPE
      
      TYPE fitgeom
        INTEGER :: flagx0(maxdimensions)				!flag = 0 <> cos series
        INTEGER :: flagv0								!flag = 1 <> sin series
        INTEGER :: flaggrad(maxdimensions)				!0 < flag < 1 fourier series
        INTEGER :: flagff(maxdimensions,maxdimensions)
		INTEGER :: flagrot(3)
		INTEGER :: flagxi(3,3)
	  	DOUBLE PRECISION :: fitx0(maxdimensions,maxorder)
	  	DOUBLE PRECISION :: fitv0(maxorder)   	  
	  	DOUBLE PRECISION :: fitgrad(maxdimensions,maxorder)   
	  	DOUBLE PRECISION :: fitff(maxdimensions,maxdimensions,maxorder)   
		DOUBLE PRECISION :: fitrot(3,maxorder)
		DOUBLE PRECISION :: fitxi(3,3,maxorder)   
      END TYPE
       
      TYPE fitdip
	INTEGER :: flagdip(3)
	INTEGER :: flagdipderiv(3,maxdimensions)
	DOUBLE PRECISION :: dip(3,maxorder)
	DOUBLE PRECISION :: dipderiv(3,maxorder,maxdimensions)
      END TYPE
      
      TYPE Gcontra
		INTEGER :: flag
		DOUBLE PRECISION :: series(maxorder)
      END TYPE
      
      TYPE gco
		INTEGER :: flag
		DOUBLE PRECISION :: series(maxorder)
      END TYPE
      
      TYPE vpseudo
		INTEGER :: flag
		DOUBLE PRECISION :: series(maxorder)
      END TYPe
      
      TYPE nc
		INTEGER :: flag(3*maxatoms-7)
		INTEGER :: flaggnc(3*maxatoms-7)
		DOUBLE PRECISION :: series(3*maxatoms-7,maxorder)
		DOUBLE PRECISION :: gncseries(3*maxatoms-7,maxorder)
      END TYPE
      
      TYPE energies
		INTEGER :: flag
		DOUBLE PRECISION :: series(maxorder)
      END TYPE
      
      TYPE  bkl
	       INTEGER :: flag1
	       INTEGER :: flag2
	       INTEGER :: flag3
	       INTEGER :: flag4
	       DOUBLE PRECISION :: series1(maxorder)
	       DOUBLE PRECISION :: series2(maxorder)
	       DOUBLE PRECISION :: series3(maxorder)
	       DOUBLE PRECISION :: series4(maxorder)
      ENDTYPE
      
      CONTAINS	

        FUNCTION DfDx(nfit,p1,flag1,flagout)

	! This function calculates the derivative of a fitted function given as parameter = flag

	  INTEGER :: nfit, flag1, flagout, i, j, k, nout
	  DOUBLE PRECISION, INTENT(IN)  :: p1(nfit)
	  DOUBLE PRECISION, DIMENSION(nfit) :: pout,DfDx
	  pout= 0.0
	  flagout = -1.*flag1 !new flag: Cos -> Sin, Sin-> -Cos
	  select case (flag1)
		case(1:4)  !either DCos/Dx -> DSin/Dx
			!WRITE (*,*) 'derive cos series' 	!p_i Cos(i*pi/180*x) -> -p_i*i*pi/180 Sin(ix)
		    pout(1) = 0.0						! constant derivative = 0
                    pout(2) =-p1(2)*0.5*DACOS(-1.d0)/180.d0		!DCos(0.5*x+pi/2) = -0.5  Cos(0.5*x) = -0.5*Sin(0.5*x+pi/2)    
		    DO i=3,nfit
		    	pout(i)=-p1(i)*(i-2)*DACOS(-1.d0)/180.d0
		    ENDDO
!		CASE (0)
!			WRITE (*,*) 'sin/cos series'
		CASE (-4:-1)
		    !WRITE (*,*) 'derive sin series' !p_i Sin(i*pi/180*x) -> p_i*i*pi/180 Cos(ix)
			pout(1) = 0.0
			pout(2)=p1(2)*0.5*DACOS(-1.d0)/180.d0
			DO i=3,nfit
				pout(i)=p1(i)*(i-2)*DACOS(-1.d0)/180.d0
			ENDDO
		CASE DEFAULT 
			STOP 'error in function derivation'			
        END SELECT
        DfDx=pout
	END FUNCTION DfDx
	
	FUNCTION mult(nfit,p1,flag1,p2,flag2,flagout,nfitout)
	
!!!!!!!!!!!This SR multiplies fitted functions by using sin/cos addition theorem
!!!!!!!!!!!!e.g. Cos(n*x)*Cos(m*x) = Cos*(m-n)*x) + Cos*(m+n)*x)
	
	  INTEGER, INTENT(IN):: nfit, flag1, flag2
	  INTEGER :: flagout, i, j, k, nout, flag1new,flag2new,nfitout
	  DOUBLE PRECISION, INTENT(IN)  :: p1(nfit),p2(nfit)
	  DOUBLE PRECISION :: pout(2*nfit),p1new(nfit),p2new(nfit),mult(2*nfit)
	
	  flagout = flag1*flag2 !new flag, compare rules for cos/sin multiplication
	  nfitout = 2*nfit-1
	  select case (flagout)
		case(1:4)  !either cos*cos or sin*sin
		!WRITE (*,*) 'cos/cos series'
		        DO i=1,nfitout;pout(i)=0.0;ENDDO
			IF (flag1.EQ.1) THEN							!cos*cos
			    pout(1) =  p1(1)*p2(1)					!i=1.j=1
			    DO i=2, nfit
				pout(i) = p1(1)*p2(i) + p1(i)*p2(1)		        !i=1,j>1 and i>1, j=1
			    ENDDO
			    DO i=2,nfit
                               pout(1) = pout(1) + 0.5*p1(i)*p2(i)				 					!p1*cos(i*x)*p2*cos(i*x)= 1/2*(p1*p2 + p1*p2 *cos(2i*x))
			       pout(2*i-1) = pout(2*i-1) + 0.5*p1(i)*p2(i)	
			       DO j=2,i-1
			         pout(i-j+1) = pout(i-j+1)+ 0.5*(p1(i)*p2(j)+p1(j)*p2(i)  )         !p1*cos(i*x)*p2*cos(j*x)= 1/2*(p1*p2*cos((i-j) *x)+ p1*p2 *cos((i+j)*x)) (+double count)
			    	 pout(i+j-1) = pout(i+j-1)+ 0.5*(p1(i)*p2(j)+p1(j)*p2(i))
			       ENDDO
			    ENDDO
			ELSEIF (flag1.EQ.-1) THEN 						!sin*sin
			    pout(1) =  p1(1)*p2(1)					!i=1.j=1 should be zero
			    DO i=2, nfit
				pout(i) = p1(1)*p2(i) + p1(i)*p2(1)		        !i=1,j>1 and i>1, j=1 should be zero
			    ENDDO
			    DO i=2,nfit
                               pout(1) = pout(1) + 0.5*p1(i)*p2(i)			    
			       pout(2*i-1) = pout(2*i-1) -0.5*p1(i)*p2(i)				!p1*sin(i*x)*p2*sin(i*x)= 1/2*(p1*p2* - p1*p2 *cos(2i*x))
			       DO j=2,i-1
			         pout(i-j+1) = pout(i-j+1)+0.5*(p1(i)*p2(j)+p1(j)*p2(i))            !p1*sin(i*x)*p2*sin(j*x)= 1/2*(p1*p2*cos((i-j) *x)- p1*p2 *cos((i+j)*x)) (+double count)
				 pout(i+j-1) = pout(i+j-1) -0.5*(p1(i)*p2(j)+p1(j)*p2(i))
			       ENDDO
			    ENDDO
			ENDIF
!		CASE (0)
!			WRITE (*,*) 'sin/cos series'
		CASE (-4:-1)
		        !WRITE (*,*) 'sin'
			IF(flag1.EQ.-1) THEN				!reorder such that p1 = cos, p2 = sin series	
				flag2new=flag1
				flag1new= flag2
				p2new(1:nfit)=p1(1:nfit)
				p1new(1:nfit)=p2(1:nfit)
			ELSEIF(flag1.EQ.1) THEN			!reorder such that p1 = cos, p2 = sin series	
				flag1new=flag1
				flag2new= flag2
				p1new(1:nfit)=p1(1:nfit)
				p2new(1:nfit)=p2(1:nfit)
			ENDIF
			DO i=1,nfitout;pout(i)=0.0;ENDDO
			pout(1) =  p1new(1)*p2new(1)					!i=1.j=1 should be zero
			DO i=2, nfit
				pout(i) = p1new(1)*p2new(i) + p1new(i)*p2new(1)		        !i=1,j>1 and i>1, j=1 should be zero
			ENDDO
			DO i=2,nfit
			       pout(2*i-1) = pout(2*i-1) +0.5*p1new(i)*p2new(i)				!p1*cos(i*x)*p2*sin(i*x)= 1/2*(p1*p2 *sin(2i*x))
			       DO j=2,i-1
			         pout(i-j+1) = pout(i-j+1)+0.5*(p1new(j)*p2new(i)-p1new(i)*p2new(j))            !p1*cos(i*x)*p2*sin(j*x)= 1/2*(-p1*p2*sin((i-j) *x)+ p1*p2 *sin((i+j)*x)) (+double count and (sin-x)=-sin(x))
				 pout(i+j-1) = pout(i+j-1) +0.5*(p1new(i)*p2new(j)+p1new(j)*p2new(i))
			       ENDDO
			ENDDO
		CASE DEFAULT 
			STOP 'error in function multiplication'			
	  END SELECT
          mult = pout
	  END FUNCTION mult
	  
      END MODULE structs