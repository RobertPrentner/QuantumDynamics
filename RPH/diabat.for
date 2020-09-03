************************************************************************

      subroutine diabat ( nx, ndat, istart, nsel, isel,
     #                    ff, w0, s, p, cn )

* 26.10.2000
* Called by MAIN
* deperturb avoided crossings of harmonic wavenumbers by diabatic rotations
* in the normal coordinate space ( of displacement vectors)

* Input:
*   NX     : number of Cartesian coordinates
*   NDAT   : number of points along the reaction path
*   ISTART : reference point for diabatic rotations
*   NSEL   : number of channels selected for diabatic rotation
*   ISEL   : channel indices selected for diabatic rotations
*   FF     : normal coordinates in terms of mass-weighted displacement vectors
*   W0     : harmonic wavenumbers (non-diagonal !)

* Output:
*   FF     : rotated normal coordinates
*   W0     : rotated harmonic wavenumbers

* Local:
*   S      : overlap of displacement vectors of consecutive points
*   P, CN  : used in ASNADC

      implicit double precision (a-h,o-z)
      parameter ( m_rot=100000 )   ! max. number of rotations
      dimension ff(nx,nx,ndat), w0(nx,nx,ndat), s(*), p(*), cn(*)
      dimension isel(*)
      common /cross_rot/ demax, rotmin, maxrot

      nx2 = nx**2
      nbas = nx-7


!
!    DEMAX  : max. energy difference of channels to be rotated
!    ROTMIN : min. required rotation matrix element
!    MAXROT : max. number of rotations per ray

      read (*,*) demax, rotmin, maxrot
      WRITE (*,2001)
      write (*,*) 'Diabatic rotation parameters: DEMAX  =',demax
      write (*,*) '                              ROTMIN ='  ,rotmin
      write (*,*) '                              MAXROT =' ,maxrot
      WRITE(*,2001)

* starting loop from first minimum
      do istep = -1 , 1 , 2
        ir_end = 1
        if ( istep.gt.0 ) ir_end = ndat

        do ir = istart+istep , ir_end , istep
          jr = ir - istep                 ! reference ray

* ...overlap in S
          call overlap ( nx-7, nx, ff(1,8,ir), ff(1,8,jr), s )
*TEST
*          call matout ( nx-7, nx-7, s )

* ...rotated (diabatic) normal coordinates in FF, diabatic wavenumbers in W0
          call asnadc ( nx-7, nx, w0(1,1,ir), ff(1,8,ir), s, p, cn,
     #                  nrot, m_rot, nsel, isel )
*
          if (nrot.gt.0) write (*,*) ' DIABAT: NROT(',ir,') = ', nrot 

        end do                           ! end loop over grid points

      end do                             ! end loop over direction ISTEP
2001  format (80('!'))
      end

************************************************************************

      subroutine asnadc ( ne, nbas, w0, c0, s, e, c, nrot, m_rot,
     #                    nsel, isel )

*  27.10.2000
*  Called by INIT_ADC
*  Find phase factors and channel numbers for a set of row ray vectors and
*  determine 2x2 rotations within the set of row vectors to deperturb crossings
*  The column ray vectors serve as reference.

*  Input:
*   NE     : Number of ray eigenvectors 
*   NBAS   : dimension of ray eigenvectors
*   W0     : adiabatic row ray eigenvalues
*   C0     : adiabatic row ray eigenvectors
*   S      : overlap matrix between row and column ray vectors
*   M_ROT  : max. total number of rotations
*   NSEL   : number of channels selected for diabatic rotation
*   ISEL   : indices of channels selected for diabatic rotation

*  Output:
*   W0     : rotated diabatic ray eigenvectors
*   C0     : rotated diabatic ray eigenvectors
*   NROT   : number of row rotations

* Working arrays:
*   E      : ray energy matrix
*   C      : ray eigenvector matrix

*  Local arrays:
*   ICH1   : row channel numbers
*   ICH2   : used to mark column ray vectors already assigned
*   INDROT : ray eigenvector indices for each row rotation
*   CROT   : first row of each row rotation matrix


      implicit double precision (a-h,o-z)
      logical flag
      parameter ( me=200 )
      dimension ich1(me), ich2(me), crot(2), indrot(2)
      dimension s(ne,ne), e(ne,ne), c(nbas,ne), c0(nbas,ne),
     #          w0(nbas,nbas), isel(*), isel1(me)

*  The following parameters determine when a crossing is to be deperturbed
*   DEMAX  : max. energy difference of channels
*   ROTMIN : min. required rotation matrix element
*   MAXROT : max. number of rotations per ray

      common /cross_rot/ demax, rotmin, maxrot

      nrot = 0
      if (ne.lt.1)  stop ' ASNADC: No ray eigenvectors ?!'
      if (ne.gt.me) stop ' ASNADC: Too many eigenvalues !'
      if (nsel.gt.me) stop ' ASNADC: Too many selected channels !'

      eps = dsqrt(1.d0/(ne+1.d0))       ! min. overlap required for assignment

      do i = 1 , ne
        do j = 1 , ne
          e(i,j) = w0(7+i,7+j)
        end do
        do j = 1 , nbas
          c(j,i) = c0(j,i)
        end do
      end do

      do ie = 1 , ne
        ich1(ie) = 0
        ich2(ie) = 1
      end do

      do 50 ie = 1 , ne
        mrot = 0                         ! rotations per ray

*  find the largest remaining overlap
        iasn1 = 0
        iasn2 = 0
        sm = -1.0
        do ie1 = 1 , ne
          if (ich1(ie1).eq.0) then           ! not yet assigned
            do ie2 = 1 , ne
              if (ich2(ie2).gt.0.and.abs(s(ie1,ie2)).gt.sm) then
                iasn1 = ie1
                iasn2 = ie2
                sm = abs(s(ie1,ie2))
              end if
            end do
          end if
        end do

*TEST
*      write (*,*) ' ASNADC: IASN1,IASN2,S = ',
*     #      iasn1, iasn2, s(iasn1,iasn2)

        if (iasn1.eq.0.or.sm.lt.eps) go to 51     ! no column ray vectors left
        ie1 = iasn1
        ie2 = iasn2
        ich1(ie1) =  ie2
        ich2(ie2) = -ich2(ie2)                    ! mark assigned column vector

        flag = .false.                            ! selected channels only
        do i = 1 , nsel
          flag = flag .or. ( ie2.eq.isel(i) )
        end do
        if (nsel.gt.0 .and. .not.flag) go to 50
 
* find second largest overlap of the selected column vector with an
* unassignes row vector
   35   if (mrot.ge.maxrot) then
          write (*,*) ' ASNADC: MAXROT reached !'
          go to 50              ! max. number of rot. reached
        end if 
        iasn = 0
        sm = -1.0
        do i = 1 , ne
          if ( ich1(i).eq.0 ) then          ! unassiged ray vector
            flag = ( nsel.le.0 )            ! only selected channels
            do l = 1 , nsel
              k = isel(l)
              do j = 1 , ne
                if ( abs(s(j,isel(l))) .gt. abs(s(k,isel(l))) ) k = j
              end do
              flag = flag .or. ( i.eq.k )
            end do
            if ( flag .and. abs(s(i,ie2)).gt.sm .and. 
     #           abs(e(ie1,ie1)-e(i,i)).le.demax ) then
              iasn = i
              sm = abs(s(i,ie2))
            end if
          end if
        end do
*TEST
*      write (*,'(" ASNADC: IE1,IASN,IE2,S,E = ",3i4,2f8.3,2f9.2)')
*     #      ie1, iasn, ie2, s(ie1,ie2), s(iasn,ie2),
*     #      e(ie1,ie1), e(iasn,iasn)

        if (iasn.eq.0) go to 50                    ! no row vectors left

*  deperturb the crossing
        ss = dsqrt( s(iasn,ie2)**2 + s(ie1,ie2)**2 )
        indrot(1) = ie1
        indrot(2) = iasn
        crot(1) = s(ie1,ie2)/ss
        crot(2) = s(iasn,ie2)/ss
        if ( abs(crot(2)).lt.rotmin ) go to 50         ! insignificant rotation

* limit rotation angle
        de = abs(e(ie1,ie1)-e(iasn,iasn))
        if ( de .gt. 0. ) then
*          call rotlim ( crot(1), crot(2), e(ie1,ie1),
*     #                  e(iasn,iasn), demax, rotmin )
          if ( abs(crot(2)).lt.rotmin ) go to 50       ! insignificant rotation
        end if
        nrot = nrot + 1

        if (nrot.gt.m_rot) stop ' ASNADC: Too many rotations.'

*TEST
*      write (*,*) ' ASNADC: CROT ', crot
*      write (*,*) ' Overlap : ',s(ie1,ie2), s(iasn,ie2)
*      write (*,'(20x,2f12.6)') ((e(indrot(i),indrot(j)),j=1,2),i=1,2)

*  rotate adiabatic overlap matrix
        call rotmat ( 1, indrot, crot, 0, indrot, crot, ne, ne, s )

*  rotate adiabatic channel potential matrix
        call rotmat ( 1, indrot, crot, 1, indrot, crot, ne, ne, e )

*TEST
*      write (*,*)
*      write (*,*) ' Overlap : ',s(ie1,ie2), s(iasn,ie2)
*      write (*,'(20x,2f12.6)') ((e(indrot(i),indrot(j)),j=1,2),i=1,2)
*      write (*,*)

*  9.5.2000: rotate row ray eigenvectors (CHECK !!!)
        call rotmat ( 0, indrot, crot, 1, indrot, crot, nbas, ne, c )

        mrot = mrot + 1                              ! rotations per row vector
        go to 35                      
   50 continue
   51 continue

*  assign the lowest possible channel number to remaining row ray vectors
      ie2 = 1
      do i = 1 , ne
        if ( ich1(i).eq.0 ) then             ! not yet assigned
   60     if ( ich2(ie2).gt.0 ) then         ! lowest possible assignment
            ich1(i) = ie2
          else if ( ie2.lt.ne ) then        ! find next possible assignment
            ie2 = ie2 + 1
            go to 60
          else                               ! no further assignments available 
            stop ' ASNADC: assignment incomplete !!'
          end if
        end if
      end do

* correct for phase jumps and rearrange row ray expectation values and
* eigenvectors
      do ie1 = 1 , ne
        ie2 = ich1(ie1)
        do je1 = 1 , ne
          je2 = ich1(je1)
          w0(7+ie2,7+je2) = e(ie1,je1)
        end do
        if ( s(ie1,ie2).lt.0.0 ) then
          phase = -1.d0
        else
          phase = 1.d0
        end if
        do i = 1 , nbas
          c0(i,ie2) = phase * c(i,ie1)
        end do

*TEST
*        write (*,*) '          ', ich1(ie1), w0(ie2,ie2), phase

      end do

      end

************************************************************************

      subroutine rotmat ( nr, indr, cr, nc, indc, cc, ni, nj, s )

*  26.11.1997
*  Called by ASNADC
*  Transform a ray block S of a matrix by a succession of 2x2 rotations:
*  CR(transpose).S.CC --> S(transformed) 

* Input:
*  NR/C : number of row/column rotations
*  INDR/C : indices for each row/column rotation
*  CR/C : first column of each row/column rotation
*  NI/J : number of rows/columns in S
*  S : ray block to be transformed

* Output:
*  S : transformed ray block  

      implicit double precision (a-h,o-z)
      dimension indr(2,*), indc(2,*), cr(2,*), cc(2,*), s(ni,*)

*  row rotations
      do ir = 1 , nr
        i1 = indr(1,ir)
        i2 = indr(2,ir)
        do j = 1 , nj
          a1 =  cr(1,ir)*s(i1,j) + cr(2,ir)*s(i2,j)
          a2 = -cr(2,ir)*s(i1,j) + cr(1,ir)*s(i2,j)
          s(i1,j) = a1
          s(i2,j) = a2
        end do
      end do

*  column rotations
      do ic = 1 , nc
        j1 = indc(1,ic)
        j2 = indc(2,ic)
        do i = 1 , ni
          a1 =  cc(1,ic)*s(i,j1) + cc(2,ic)*s(i,j2)
          a2 = -cc(2,ic)*s(i,j1) + cc(1,ic)*s(i,j2)
          s(i,j1) = a1
          s(i,j2) = a2
        end do
      end do

      end

************************************************************************


      subroutine overlap ( nc, nx, c, c0, s )

* 16.05.2000
* Calculate overlap matrix between between two sets of vectors: S = C^t.C0

* Input:
*  NC : number of vectors in each set
*  NX : basis size
*  C  : 1st set of vectors
*  C0 : 2nd set of vectors

* Output:
*  S  : overlap matrix

      implicit double precision (a-h,o-z)
      dimension c(nx,nc), c0(nx,nc), s(nc,nc)

      do ic = 1 , nc
        do jc = 1 , nc
          sum = 0.0
          do i = 1, nx
            sum = sum + c(i,ic)*c0(i,jc)
          end do
          s(ic,jc) = sum
        end do
      end do

      end

************************************************************************

      subroutine gradf ( n1, n, n3, y , y1, idif )

* 11.10.2000
* Numerical first derivative with respect second grid index
* Called by MAIN, GMAT, MAKE_NC, VPSEUD

* Input:
*   N1 , N3 : inner/outer dimension
*   N       : number of grid points
*   Y       : function value
*   IDIF    : differentiation scheme:
*                =1 : 2nd order finite differences
*                =2 : periodic (Fourier) grid derivative

* Output:
*   Y1      : first derivative

      implicit double precision  (a-h,o-z)
      logical flag
      parameter ( m1=1001 )
      dimension y(n1,n3,n), y1(n1,n,n3)
      dimension ind(m1), f(m1), df(m1)

      if (n.lt.1) stop ' FINDIF: Programming error.'
      if (n.gt.m1) stop ' FINDIF: Too many points.'

      do i3 = 1 , n3
        do i1 = 1 , n1

          do i = 1 , n
            ind(i) = i
            f(i) = y(i1,i3,i)
          end do
          if ( idif.eq.0 ) then
            call findif ( n, ind, f, df )
          else
            flag = .false.
            call fourdif ( n, f, df, flag )
          end if
          do i = 1 , n
            y1(i1,i,i3) = df(i)
          end do

        end do
      end do

      call smooth ( n1, n, n3, y1 )

      end

************************************************************************

      subroutine findif ( n, ind, f, df )

* 11.10.2000
* First derivative by second order finite differences

* Input:
*  N   : number of grid points
*  IND : grid indices
*  F   : function values

* Output:
*  DF  : derivative of F w.r.t. the grid index

      implicit double precision  (a-h,o-z)
      dimension ind(*), f(*), df(*)

      if (n.eq.1) then
        df(1) = 0.0
      else if (n.eq.2) then
        df(1) = (f(2)-f(1))/(ind(2)-ind(1))
        if ( n.gt.1 ) df(2) = df(1)
      else if (n.eq.3) then
        call df1 ( 1, 2, 3, ind, f, df(1) )
        if (n.gt.1) call df1 ( 2, 1, 3, ind, f, df(2) )
        if (n.gt.2) call df1 ( 3, 1, 2, ind, f, df(3) )
      else
        call df2 ( 1, 2, 3, 2, 4, ind, f, df(1) )
        call df2 ( 2, 1, 3, 1, 4, ind, f, df(2) )
        do i = 3 , n-2
          call df2 ( i, i-1, i+1, i-2, i+2, ind, f, df(i) )
        end do
        call df2 ( n  , n-1, n-2, n-2, n-3, ind, f, df(n  ) )
        call df2 ( n-1, n-2, n  , n-3, n  , ind, f, df(n-1) )
      end if

      end

************************************************************************

      subroutine df1( i0, i1, j1, ind, Y, Y1 )

*  30.6.1999
*  first order finite difference (FD) formular for an incomplete equidistant
*  grid with unit mesh size

*  Input:
*    I0          : reference index for center point
*    I1,J1       : reference indices for finite differences ( I1 must be
*                  different from J1 )
*    IND         : grid indices
*    Y           : function values

*  Output:
*    Y1          : 1st order finite difference approximation to the
*                  first derivative of Y

      implicit double precision (a-h,o-z)
      dimension ind(*), y(*)
      
      i = ind(i1) - ind(i0)
      j = ind(j1) - ind(i0)
      fi = (y(i1)-y(i0))/i
      fj = (y(j1)-y(i0))/j
      y1 = ( j*fi - i*fj ) / (j-i)

      end

************************************************************************

      subroutine df2( i0, i1, j1, i2, j2, ind, Y, Y1 )

*  30.6.1999
*  second order finite difference (FD) formular for an incomplete equidistant
*  grid with unit mesh size

*  Input:
*    I0          : reference index for center point
*    I1,J1,I2,J2 : reference indices for finite differences ( I1 must be
*                  different from J1 and I2 must be different from J2 )
*    IND         : grid indices
*    Y           : function values

*  Output:
*    Y1          : 2nd order finite difference approximation to the
*                  first derivative of Y

      implicit double precision (a-h,o-z)
      dimension ind(*), y(*)
      
      i = ind(i1) - ind(i0)
      j = ind(j1) - ind(i0)
      k = ind(i2) - ind(i0)
      l = ind(j2) - ind(i0)
      ij = i*j
      kl = k*l
      fi = (y(i1)-y(i0))/i
      fj = (y(j1)-y(i0))/j
      fk = (y(i2)-y(i0))/k
      fl = (y(j2)-y(i0))/l
      fij = ( j*fi - i*fj ) / (j-i)
      fkl = ( l*fk - k*fl ) / (l-k)
      y1 = ( kl*fij - ij*fkl ) / (kl-ij)

      end

************************************************************************

      subroutine fourdif ( n, f, df, init )

*  11.10.2000
*  First derivative on a grid with respect to the grid index for periodic
*  functions.
*  The linear correction has been switched off for ISYM=0 and -1.
*  In the present form it causes artifacts near the boundaries since 
*  non-symmetric functions need not be the same at the two boundaries.

*  INPUT:
*    ID   : grid index
*    N    : number of points
*    F    : function values
*    INIT : true --> initialise derivative operator for N points in [-pi,pi]
*                    ( Return without performing any differentiation.)

*  OUTPUT:
*    DF : 1st derivative

      implicit double precision (a-h,o-z)
      logical init
      parameter ( m1=1001 )
      dimension f(*), df(*)
      common /fourier/ dd(0:m1), sk(m1)

      if ( init ) then
        if (n.gt.m1)   stop ' FOURDIF: Too many points.'
        pi = dacos(-1.d0)
        w = pi/n
        dd(0) = 0.0
        do imj = 1 , n-1
          dd(imj) = (-1)**imj*w/dsin(w*imj)
        end do
* correction for linear term
        do j = 1 , n
          sk(j) = 1.0
          do k = 1 , n
            sk(j) = sk(j) - k*sign(1,j-k)*dd(iabs(j-k))
          end do
        end do
        return
      end if

      if (n.lt.1) then
        stop ' FOURDIF: Programming error.'
      else if (n.eq.1) then
        df(1) = 0.0
        return
      end if

*      a = ( f(n) - f(1) ) / ( n - 1 )     ! correction for linear term
      do  i = 1 , n
        df(i) = 0.0
*        df(i) = a*sk(i)                   ! correction for linear term
        do  k = 1 , n
          df(i) = df(i) + sign(1,i-k)*dd(iabs(i-k))*f(k)
        end do
      end do

      end

************************************************************************

      subroutine smooth ( n1, n, n3, f )

* 3.11.2000
* Smooth a function F given on a Fourier-DVR grid with N points by
* interpolation on to a coarser grid with FAC*N points and back on
* to the original grid (should be equivalent to a Wien filter).

* Input:
*   FAC  : smoothing factor ( 0.le.FAC.le.1 )
*   N    : number of grid points
*   F    : function to smooth
*   INIT : flag for initialization

* Output:
*   F    : smoothed function

      implicit double precision ( a-h,o-z)
      logical init
c change (m1=maxpt??)
      parameter ( m1=1001, m2=m1**2 )
      dimension f(n1,n,n3), q(m1), q0(m1), d(m2,2)
      common /c_smooth/ fac

      n0 = nint( fac * n )
      if ( mod(n0,2).eq.0 ) n0 = n0 + 1
      if ( n0.lt.1 .or. n0.gt.n ) stop ' SMOOTH: FAC out of range.'
      write (*,*) ' SMOOTH: FAC = ', fac, '      N0 = ', n0
      do i = 1 , n
        q(i) = i-0.5d0
      end do
      do i = 1 , n0
        q0(i) = ( n * q(i) ) / n0
      end do
      call intpol ( n, q, n0, q0, d(1,1) )
      call intpol ( n0, q0, n, q, d(1,2) )
      do i3 = 1 , n3
        do i1 = 1 , n1
          do i = 1 , n
            q(i) = f(i1,i,i3)
          end do
          call vec_x_mat ( n, q, d(1,1), n0, q0 )
          call vec_x_mat ( n0, q0, d(1,2), n, q )
          do i = 1 , n
            f(i1,i,i3) = q(i)
          end do
        end do
      end do
      end

************************************************************************

      subroutine vec_x_mat ( n, f, d, n0, f0 )

* 3.11.2000
* Row vector F times matrix D

* Input:
*  N : dimension of F
*  F : row vector
*  D : matrix
*  N0 : dimension of F0

* Output:
*  F0 : transformed row vector

      implicit double precision ( a-h,o-z )
      dimension f(n), f0(n0), d(n,n0)
      do i0 = 1 , n0
        f0(i0) = 0.d0
        do i = 1 , n
          f0(i0) = f0(i0) + f(i) * d(i,i0)
        end do
      end do
      end

************************************************************************

      subroutine intpol ( ns, qs, nq, q, d )

C Matrix for interpolation from a grid of NS points to one of NQ points using
C Meyer's equidistant DVR: If S is a row vector on the original grid with NS
C points, the interpolated row vector F is given by F=S*D.

C input:
C  NS : number of points on the original grid
C  NQ : number of points on the new grid
C  QS : original grid coordinates
C  Q  : new grid coordinates

C output:
C  D  : interpolation matrix. D(I,J) is the value of the I-th DVR basis
C       function at Q(J)

      implicit double precision (a-h,o-z)
      dimension  d(ns,nq), q(nq), qs(ns)

      scale = 2.*dacos(-1.d0)/(ns*(qs(2)-qs(1)))
      eps = 2.d-6*dacos(-1.d0)/ns
      do 5 is = 1 , ns
         do 5 iq = 1 , nq
            xmxj = (q(iq)-qs(is))*scale
            if ( abs(xmxj).ge.eps) then
               d(is,iq) = dsin(ns*0.5d0*xmxj)/(ns*sin(0.5d0*xmxj))
            else
               d(is,iq) = 1.d0
            end if
   5  continue

      end


      subroutine b_mat ( ndat, nat, c, g, b )

* 07.06.2000
* Called by MAIN
* Calculate the generalized Coriolis interaction parameters B(K,L,IQ) between
* normal modes K and L. IQ=1,2,3 are the Cartesian axes and IQ=4 the reaction
* path coordinate. The gradient of the Cartesian coordinates w.r.t. the
* reaction path has to be orthogonal to the displacement vectors C

* Input:
*   NDAT    : number of gridpoints
*   NAT     : number of atoms
*   C       : Cartesian displacement vectors for the normal coordinates at each
*             grid point
*   G       : gradient of normalcoordinates along the grid s

* Output:
*   B   : generalized Coriolis interaction parameters

      implicit double precision (a-h,o-z)
      dimension c(3*nat,3*nat,ndat), g(3,nat,3*nat,ndat),
     #          b(3*nat-7,3*nat-7,4,ndat)

      nx = 3 * nat           ! number of Cartesian coordinates

      do i = 1, ndat
        do k = 8, nx
          kk = k-7
          do ib = 1 , 4
            b(kk,kk,ib,i) = 0.d0
          end do
          do l = 8, k-1
            ll = l-7
            do ib = 1 , 4
              b(kk,ll,ib,i) = 0.d0
            end do
            do j = 1, nat
	      c1k = c(3*(j-1)+1,k,i)
	      c2k = c(3*(j-1)+2,k,i)
	      c3k = c(3*(j-1)+3,k,i)
	      c1l = c(3*(j-1)+1,l,i)
	      c2l = c(3*(j-1)+2,l,i)
	      c3l = c(3*(j-1)+3,l,i)
              b(kk,ll,1,i) = b(kk,ll,1,i) + c2k * c3l
     #                                    - c3k * c2l
              b(kk,ll,2,i) = b(kk,ll,2,i) + c3k * c1l
     #                                    - c1k * c3l
              b(kk,ll,3,i) = b(kk,ll,3,i) + c1k * c2l
     #                                    - c2k * c1l

* exploit the fact that the derivative of the norm vanishes to reduce
* inaccuracies due to the differentiation of the displacement vectors
              b(kk,ll,4,i) = b(kk,ll,4,i) + g(1,j,k,i) * c1l
     #                                    + g(2,j,k,i) * c2l
     #                                    + g(3,j,k,i) * c3l
     #                                    - c1k * g(1,j,l,i)
     #                                    - c2k * g(2,j,l,i)
     #                                    - c3k * g(3,j,l,i)
            end do

            b(kk,ll,4,i) = 0.5d0 * b(kk,ll,4,i)
            do ib = 1 , 4
              b(ll,kk,ib,i) = -b(kk,ll,ib,i)
            end do

          end do
        end do
      end do

      end
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      