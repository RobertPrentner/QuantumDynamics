      program MKINRPH

*  last update: 21.1.2009, Robert Prentner 
*  Read Gaussian archive blocks and generate an appropriate input file for
*  RPH_CALC. Written by David Luckhaus

      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      INTEGER bdip
      parameter ( ma=50, mvar=100, mx=3*ma, mf=mx*(mx+1)/2,
     #            inp=1, iout=2, izfile=3 )

      character fname*79, rec*160, qname*8
      dimension x(mx), ff(mf), g(mx),dip(5)

* used by READ_ZMAT

      dimension am(ma), inda(ma), indz(3*ma), izvar(3,mvar),
     #          z(3,ma), q(mvar), qname(mvar), fac(mvar)

      read (*,*) eref
      write (*,*) 'Reference energy (EREF) = ', eref
      READ (*,*) bdip
      WRITE (*,*) 'Print Dipole:', bdip
      READ(*,*) potflag

*  read the Z-matrix
      call read_zmat ( izfile, na, am, inda, naz, indz, z, ma,
     #                 nvar, izvar, mvar, qname )
      do i = 1 , nvar
        iq = izvar(3,i)
        fac(i) = 1.d0
        if (iq.lt.0) then
          fac(i) = -fac(i)
          iq = -iq
          izvar(3,i) = iq
        end if
      end do

    1 read (*,*,end=999) fname, icoord
      open (inp,file=fname,form='formatted',status='unknown')
      read (inp,'(1x,a70)',end=99) rec
      read (inp,'(1x,a70)',end=99) rec(71:140)
      k = 0
      nx = 0
      nf = 0
      ng = 0
      nq = -1

    5 k = k + 1
      if ( rec(k:k).eq.'@' ) go to 99

      if ( k.gt.70 ) then
        k = k - 70
        rec(1:70) = rec(71:140)                 
        read (inp,'(1x,a70)',end=99) rec(71:140)
      end if

*-----------------------------------------------------------------------

      if ( rec(k:k+5).eq.'\\0,1\' ) then           ! Cartesian coordinates next
        k = k + 6
        if ( icoord.eq.0 ) call get_cart ( inp, k, rec, nx, x )
        nq = 0
        go to 5

*-----------------------------------------------------------------------

      else if ( rec(k:k+1).eq.'\\' .and. nq.eq.0 .and.     ! Z-matrix variables
     #          icoord.ne.0 ) then
        k = k + 2
        call get_q ( inp, k, rec, nq, q, nvar, qname )
        WRITe(*,*) 'internal coordinate values'
        do i = 1 , nq
          write (*,'(1x,a8,f16.10)') qname(i), q(i)
        end do
        do i = 1 , nvar
          z(izvar(2,i),izvar(1,i)) = q(izvar(3,i)) * fac(i)
        end do
        call corgen_z ( na, am, x, inda, naz, indz, z, tau )
        call x_gaussian ( na, x )
        nx = 3*na
        go to 5

*-----------------------------------------------------------------------

      else if ( rec(k:k+4).eq.'\MP2=' ) then      ! electronic energy
        k = k + 5
        l = k
   25   if ( rec(l:l).ne.'\' ) then    ! end of data set
          l = l + 1
          go to 25
        end if
        read (rec(k:l-1),*) e
        k = l
	goto 5
*-----------------------------------------------------------------------    
      else if ( rec(k:k+7).eq.'\Dipole='.AND.bdip.ne.0 ) then ! dipoles
        WRITE (*,*) 'dipole found'
        k = k + 8
        l = k
        do idip = 1,3
   27   if ( rec(l:l).ne.'\' .and. 
     #       rec(l:l).ne.','    ) then    ! end of data set 
          l = l + 1
          go to 27 
        end if
         read (rec(k:l-1),*)ee
          dip(idip) = ee
        WRITE (*,*) 'value',ee,idip,rec(l:l)
         if ( rec(l:l).eq.',') l = l + 1
        k=l
        end do
        go to 5 

*-----------------------------------------------------------------------

      else if ( rec(k:k+6).eq.'\NImag=' ) then      ! force field
        k = k + 10
   30   if ( k.gt.70 ) then
          k = k - 70
          rec(1:70) = rec(71:140)                 
          read (inp,'(1x,a70)',end=99) rec(71:140)
        end if
        l = k
   35   if ( rec(l:l).ne.',' .and. rec(l:l).ne.'\' ) then        ! next number
          l = l + 1
          go to 35
        end if
        nf = nf + 1
        read (rec(k:l-1),*) ff(nf)
        k = l + 1
        if ( nf.lt.nx*(nx+1)/2) go to 30
        go to 5      

*-----------------------------------------------------------------------

      else if ( nf.gt.0 .and. ng.eq.0 ) then                      !  gradient
        k = k + 1
   40   if ( k.gt.70 ) then
          k = k - 70
          rec(1:70) = rec(71:140)                 
          read (inp,'(1x,a70)',end=40) rec(71:140) !end=99
        end if
        l = k
   45   if ( rec(l:l).ne.',' .and. rec(l:l).ne.'\' ) then        ! next number
          l = l + 1
          go to 45
        end if
        ng = ng + 1
        read (rec(k:l-1),*) g(ng)
        k = l + 1
        if ( ng.lt.nx ) go to 40
      else

        go to 5
      end if

   99 IF(potflag.NE.0) THEN 
		READ (11,*) e
		e = -e
		eref = eref
	ENDIF
      write (iout,'(3f18.12)') (x(i),i=1,nx)
      !write (iout,'(/f18.12/)') e - eref
      write (iout,'(/f18.12/)') e - eref
      write (iout,'(3f18.12)') (g(i),i=1,ng)
      write (iout,*)
      write (iout,'(5f18.12)') (ff(i),i=1,nf)
      write (iout,*) 
      write (*,'(1x,a,f18.12)') 'V(tau) - EREF', e-eref
      WRITE(*,*)
      
*------------If dipole is processed (bdip = true)------------------------- 
C     Change Dipole from a.u. to Debye
C     Rotate Dipole around Tau
      IF(bdip.ne.0) THEN 
 	     WRITE (*,'(A,4f16.10)') 'dip before rot',dip(1),dip(2),dip(3),tau

    	  DIP1=COS(TAU/2.)*DIP(2)
     # 	      -SIN(TAU/2.)*DIP(1)


      	  DIP(1)=SIN(TAU/2.)*DIP(2)
     #    		  +COS(TAU/2.)*DIP(1)
     
	 DIP(2)=DIP1
	
       	WRITE (*,'(A,4f16.10)') 'dip after rot',dip(1),dip(2),dip(3),tau


        IF(abs(DIP(1)).LT.10**(-5))DIP(1)=0.00000000000
        IF(abs(DIP(2)).LT.10**(-5))DIP(2)=0.00000000000


      
      	Write (iout,'(3f18.12)') (2.541765*dip(i),i=1,3)
      	WRITE (iout,*) 
      ENDIF		
      close (inp)

      go to 1

  999 continue

      end

************************************************************************

      subroutine x_gaussian ( na, x )

*  rearrange Cartesian axes to conform with the Z-matrix convention of
*  GAUSSIAN: (x,y,z) --> (y,-z,x)

      implicit double precision (a-h,o-z)
      dimension x(3,na), xg(3)
      do ia = 1 , na
        xg(1) =  x(2,ia)
        xg(2) = -x(3,ia)
        xg(3) =  x(1,ia)
        x(1,ia) = xg(1)
        x(2,ia) = xg(2)
        x(3,ia) = xg(3)
      end do
      end

************************************************************************

      subroutine get_cart ( inp, k, rec, nx, x )

      implicit double precision (a-h,o-z)
      character rec*160
      dimension x(*) 

   10   if ( k.gt.70 ) then
          k = k - 70
          rec(1:70) = rec(71:140)                 
          read (inp,'(1x,a70)',end=99) rec(71:140)
        end if
        l = k
   15   if ( rec(l:l).ne.'\' ) then    ! end of data set
          l = l + 1
          go to 15
        end if
        if ( k.lt.l ) then             ! read one set of Cartesian coordinates
   16     if ( rec(k:k).ne.',' ) then
            k = k + 1
            go to 16
          end if
          read (rec(k+1:l-1),*) (x(i),i=nx+1,nx+3)
          nx = nx + 3
          k = l + 1
          go to 10
        end if
        k = l

   99 continue
      end

************************************************************************

      subroutine get_q ( inp, k, rec, nq, q, nvar, qname )

      implicit double precision (a-h,o-z)
      character rec*160, qname*8, test*8
      dimension q(*) , qname(*)

   10 if ( k.gt.70 ) then
        k = k - 70
        rec(1:70) = rec(71:140)                 
        read (inp,'(1x,a70)',end=99) rec(71:140)
      end if
      l = k
   15 if ( rec(l:l).ne.'\' ) then    ! end of data set
        l = l + 1
        go to 15
      end if
      i1 = k
      if ( k.lt.l ) then             ! read one set of Cartesian coordinates
   16   if ( rec(k:k).ne.'=' ) then
          k = k + 1
          go to 16
        end if
        do iq = 1 , nvar
          test = qname(iq)
          if ( rec(i1:k-1).eq.test(1:k-i1) )
     #      read (rec(k+1:l-1),*) q(iq)
        end do
        nq = nq + 1
        k = l + 1
        if ( nq.lt.nvar ) go to 10
      end if
      k = l

   99 continue
      end

************************************************************************

      subroutine READ_ZMAT ( input, na0, am, inda, na, indz, z, maxa,
     #                       nzvar, izvar, maxvar, qname )

*  24.8.1999
*  READ a Z-matrix and generate a list of variables

*  input :
*    MAXA   : max. allowed number of rows in the Z-matrix
*    MAXVAR : max. allowed number of variables

*  output :
*    NA0    : number of atoms (without pseudo-atoms)
*    AM     : atomic masses (a.m.u.)
*    INDA   : row number in the Z-matrix for each atom
*    NA     : number of rows in the Z-matrix (i.e. incl. pseudo-atoms)
*    INDZ   : Z-matrix definition ( reference atoms )
*    Z      : Z-Matrix coordinates ( values in Angstrom,Degrees)
*    NZVAR  : number of variable Z-matrix elements
*    IZVAR  : For each variable Z-matrix element: (1) row and (2) column
*             in the Z-matrix, and (3) index of variable

*  input read from unit INPUT:
*    The name of each variable occuring in the Z-matrix.
*    Floating variables must precede fixed ones.
*    If a variable is to be held at fixed value, the name must be followed
*    by the corresponding value in the same record.
*    The List of variables is finished by a blank line.
*    This is followed by the symbolic Z-matrix. Each atom is defined by its
*    name (upper and lower case are distinguished) followed by its atomic mass.
*    Pseudo atoms are defined through negative mass.
*    The following Z-matrix for CH2O uses the umbrella angle TAU to describe
*    the out of plane motion and BETA for (the planar projection of) the
*    symmetric OCH bend:
*------------------------------------------------------------------------------
*      BETA 
*      TAU
*      RCH  1.1033
*      RCO  1.2096
*
*      C  12.00000000d0
*      X  -1            , C , 1.0
*      O  15.99491463d0 , C , RCO ,  X , TAU
*      H1 1.007825035d0 , C , RCH ,  X , TAU ,  O ,  BETA
*      H2 1.007825035d0 , C , RCH ,  X , TAU ,  O , -BETA
*------------------------------------------------------------------------------

      implicit real*8 (a-h,o-z)
      logical ifnum, eor, fix
      character empty*8, qname*8, record*79, ref*8, sym*8, var*8
      parameter ( len=79, ma=50, mq=3*ma )
      dimension indz(3,maxa), z(3,maxa), am(maxa), inda(maxa),
     #          izvar(3,maxvar)
      dimension sym(ma), ref(3,ma), qname(mq), qval(mq), fix(mq)

      empty = '        '
      nq = 0
      nzvar = 0
      na0 = 0
      na = 0

      if (maxa.gt.ma) then
        write (*,'(A,i4)') 'READ_ZMAT: increase MA to', maxa
        stop
      end if

    1 read (input,'(a79)') record
      write (*,*) record
      i1 = 1
      var = empty
      call getrec ( var, dummy, record, i1, len, ifnum, eor )
      if (eor) go to 9                                        ! empty record
      if (nq.lt.mq) then
        nq = nq + 1
        qname(nq) = var
        PRINT *,'VAR NAME',VAR
        fix(iq) = .false.
      else
        write (*,*) ' READ_ZMAT : Too many variables.'
        stop
      end if
      call getrec ( var, qval(nq), record, i1, len, ifnum, eor )
      if (ifnum) fix(nq) = .true.
      go to 1
    9 continue

   10 read (input,'(a79)',end=19) record
      write (*,*) record
      na = na + 1
      i1 = 1
      sym(na) = empty
      call getrec ( sym(na), dummy, record, i1, len, ifnum, eor )
      call getrec ( sym(na), amass, record, i1, len, ifnum, eor )
      do i = 1 , min(3,na-1)
        ref(i,na+1) = empty
        call getrec ( ref(i,na), dummy, record, i1, len, ifnum, eor )
        var = empty
        call getrec ( var, z(i,na), record, i1, len, ifnum, eor )
        if (.not.ifnum) then
          iq = 0
          il = 0
          do l = 1 , nq
            if (.not.fix(l)) il = il + 1     ! count floating variables
            if ( qname(l).eq.var ) then
              iq = l
              iqq = il
            end if
          end do
          if ( iq.eq.0 ) then
            write (*,*) ' READ_ZMAT : undefined variable ',var
            stop
          end if
          if (fix(iq)) then
            z(i,na) = z(i,na) * qval(iq)
          else
            nzvar = nzvar + 1           
            if (nzvar.gt.maxvar) then
              write (*,*) ' READ_ZMAT : Too many variables.'
              stop
            end if
            izvar(1,nzvar) = na
            izvar(2,nzvar) = i
            izvar(3,nzvar) = iqq * nint(z(i,na))
           end if
        end if
      end do
      do ia = 1 , na-1
        if ( sym(ia).eq.sym(na) ) then
          write (*,*) ' READ_ZMAT : Ambiguous atom symbol.'
          stop
        end if
      end do
      if (amass.gt.0.0) then
         na0 = na0 + 1
         inda(na0) = na
         am(na0) = amass
      end if
      if (na.lt.maxa) then
         go to 10 
      else
         stop ' READ_ZMAT : Too many rows in the Z-Matrix.'
      end if
   19 write (*,*) na,' atoms.'

*  find indices of reference atoms
      do ia = 1 , na
        do i = 1 , min(3,ia-1)
          iref = 0
          do l = 1 , na
            if ( ref(i,ia).eq.sym(l) ) iref = l
          end do
          if ( iref.eq.0 ) then
            write (*,*) ' READ_ZMAT : ',ref(i,ia),' not defined.'
            stop
          end if
          indz(i,ia) = iref
        end do
      end do

      do ia = 1 , na
        amass = 0.
        do i0 = 1 , na0
          if (inda(i0).eq.ia) amass = am(i0)
        end do
        write (*,1001) ia,amass,(indz(i,ia),z(i,ia),i=1,min(3,ia-1))
      end do

      do i = 1 , nzvar
        write (*,*) i,(izvar(k,i),k=1,3)
      end do

 1001 format (i4,f12.6,3(i5,f15.8))

      end

************************************************************************

      subroutine getrec ( sym, x, record,  i1, len, ifnum, eor )

*  7.6.1999
*  Read a word from RECORD starting at I1. Blank and comma are treated as
*  separators. 

*  Input:
*    RECORD : Record to read from
*    I1     : first character to read
*    LEN    : total length of RECORD

*  Output:
*    SYM    : string ("symbolic input", max. 6 characters)
*    X      : numerical value (or sign factor)
*    I1     : next character to read
*    IFNUM  : .TRUE. if input was numeric
*    EOR    : .TRUE. if the end of record was found

      implicit real*8 (a-h,o-z)
      character sym*8, record*79
      logical ifnum, eor, ifd1, ifd2
      character s1*1, s2*1
      parameter ( len_sym=8 )                 ! number of characters in SYM

      eor = .false.

*  find start of the next word
   11 s1 = record(i1:i1)
      if (s1.eq.' '.or.s1.eq.',') then
        i1 = i1 + 1
        if (i1.lt.len) go to 11
      end if

*  record ended without data
      if (s1.eq.' '.or.s1.eq.',') then
        eor = .true.
        return
      end if

* find end of word
      i2 = i1
   21 if (i2.lt.len.and.record(i2+1:i2+1).ne.' '.and.
     #                  record(i2+1:i2+1).ne.','      ) then
        i2 = i2 + 1
        if (i2.lt.len) go to 21
      end if
      
      call numdig ( s1, ifd1 )
      if (i2.gt.i1) then
        s2 = record(i1+1:i1+1)
        call numdig ( s2, ifd2 )
      else
        s2 = ' '
        ifd2 = .false.
      end if

      if ( ifd1 .or. s1.eq.'.' .or.
     #     ifd2 .and. (s1.eq.'+' .or. s1.eq.'-') ) then
        ifnum = .true.                                   ! numerical input
        read (record(i1:i2),*) x
      else
        ifnum = .false.                                  ! symbolic input
        x = 1.d0
        if (s1.eq.'+'.and.i2.gt.i1) then                 ! signed variables (+)
          x = 1.d0
          i1 = i1 + 1
        else if (s1.eq.'-'.and.i2.gt.i1) then            ! signed variables (-)
          x = -1.d0
          i1 = i1 + 1
        end if
        sym = record(i1:min(i2,i1+len_sym-1))
      end if

      i1 = i2 + 1

      end

************************************************************************

      subroutine numdig ( s1, ifdig )

*  7.6.1999
*  check for numerical input

*  Input:
*   S1 : character to check

*  Output:
*   IFDIG : .TRUE. if S1 is a digit

      logical ifdig
      character*1 s1

      ifdig =  s1.eq.'0' .or. s1.eq.'1' .or. s1.eq.'2' .or.
     #         s1.eq.'3' .or. s1.eq.'4' .or. s1.eq.'5' .or.
     #         s1.eq.'6' .or. s1.eq.'7' .or. s1.eq.'8'. or .s1.eq.'9'

      end

************************************************************************

      subroutine CORGEN_Z ( na0, am, x0, inda, na, indz, z, tau )

*  14.12.1997
*  Convert a structure given in Z-matrix coordinates to one in cartesian
*  coordinates. To start with the first atom is put at the origin, the
*  second on the X-axis (X>0), and the third in the (X,Y)-plane (Y>0). 
*  The remaining atoms are then positioned in this reference system.
*  Finally the coordinate system is transformed to the centre of mass 
*  system.

*  input :
*    NA0  : number of atoms (without pseudo-atoms)
*    AM   : atomic masses [u] ( for mass-weighting, currently not active )
*    INDA : row number in the Z-matrix for each atom
*    NA   : number of rows in the Z-matrix (i.e. incl. pseudo-atoms)
*    INDZ : Z-matrix definition ( reference atoms )
*    Z    : Z-Matrix coordinates ( values in Angstrom,Degrees)
*  output :
*    X0   : centre of mass coordinates

      implicit real*8 (a-h,o-z)
      parameter (maxa=100)
      dimension am(*), x0(3,na0), inda(*), indz(3,*), z(3,na)
      dimension x(3,maxa), xcm(3)

      if (na.gt.maxa) stop ' CORGEN : Too many rows in the Z-Matrix.'
      pi = acos(-1.d0)

*  cartesian coordinates in Z-Matrix orientation (Angstrom)
      x(1,1) = 0.0
      x(2,1) = 0.0
      x(3,1) = 0.0
      x(1,2) = z(1,2)
      x(2,2) = 0.0
      x(3,2) = 0.0
      x(1,3) = x(1,indz(1,3))
      if ( indz(1,3).eq.2 ) then
         x(1,3) = x(1,3) - z(1,3)*cos(z(2,3)/180.d0*pi)
      else if ( indz(1,3).eq.1 ) then
* Vorzeichen am 15.1.1996 korrigiert
         x(1,3) = x(1,3) + z(1,3)*cos(z(2,3)/180.d0*pi)
      else
         stop ' CORGEN : First bond angle undefined.'
      end if
      x(2,3) = z(1,3)*sin(z(2,3)/180.d0*pi)
      x(3,3) = 0.0
      do 10 ia = 4 , na
         alpha = z(2,ia)/180.d0*pi
         tau   = z(3,ia)/180.d0*pi
         rp    = z(1,ia)
         call zmat( x(1,ia), x(1,indz(1,ia)), x(1,indz(2,ia)),
     #                       x(1,indz(3,ia)), rp, alpha, tau )
   10 continue

*  centre of mass coordinates
      do 15 k = 1 , 3
   15    xcm(k) = 0.0
      s = 0.
      do 20 i0 = 1 , na0
         s = s + am(i0)
         do 20 k = 1 , 3
            x0(k,i0) = x(k,inda(i0))
   20       xcm(k) = xcm(k) + x0(k,i0)*am(i0)
      do 25 k = 1 , 3
   25    xcm(k) = xcm(k)/s
c      do 30 ia = 1 , na0
c         do 30 i = 1 , 3
c   30       x0(i,ia) = x0(i,ia)-xcm(i)

      end

************************************************************************

      subroutine ZMAT ( p, a, b, c, r, alpha, tau )

*  28.4.1999
*  Calculate cartesian coordinates of a point P from Z-matrix coordinates
*  and the cartesian coordinates of the corresponding reference points.
*  Right handed coordinate systems are assumed.

*  input :
*   A,B,C : cartesian coordinates of reference points defining the local
*           coordinate system
*   R     : distance between P and A
*   ALPHA : polar angle PAB
*   TAU   : dihedral angle between planes PAB and ABC correponding to a
*           right-handed rotation around A-->B of P out of the ABC plane
*            (sense of rotation corrected on 17-jan-1996)

*  ouput :
*   P : cartesian coordinates of point P

      implicit real*8 (a-h,o-z)
      parameter ( eps=1.d-12 )
      dimension p(3), a(3), b(3), c(3), e(3,3)

*  vanishing bond length
      if (r.lt.eps) then
        p(1) = a(1)
        p(2) = a(2)
        p(3) = a(3)
        return
      end if

*  local cartesian coordinates
      x =  r * cos(alpha)
      sa = sin(alpha)
      y =  r * sa * cos(tau)
      z =  r * sa * sin(tau)

*  local axis system (origin at A):

*   X =: A-->B
      rab = sqrt( (b(1)-a(1))**2 + (b(2)-a(2))**2 + (b(3)-a(3))**2 )
      if (rab.eq.0.0) stop ' ZMAT: vanishing RAB.'
      e(1,1) = (b(1)-a(1))/rab
      e(2,1) = (b(2)-a(2))/rab
      e(3,1) = (b(3)-a(3))/rab

*  linear connection (corrected: 28.4.1999)
      if ( abs(sa).lt.eps ) then
        p(1) = x*e(1,1) + a(1)
        p(2) = x*e(2,1) + a(2)
        p(3) = x*e(3,1) + a(3)
        return
      end if

*   Z =: (A-->B).x.(A-->C)
      e(1,2) = c(1)-a(1)
      e(2,2) = c(2)-a(2)
      e(3,2) = c(3)-a(3)
      e(1,3) =  e(2,1)*e(3,2) - e(3,1)*e(2,2)
      e(2,3) = -e(1,1)*e(3,2) + e(3,1)*e(1,2)
      e(3,3) =  e(1,1)*e(2,2) - e(2,1)*e(1,2)
      rz = sqrt( e(1,3)**2 + e(2,3)**2 + e(3,3)**2 )
      if ( rz.eq.0.0 ) stop ' ZMAT: undefined TAU.'
      e(1,3) = e(1,3)/rz
      e(2,3) = e(2,3)/rz
      e(3,3) = e(3,3)/rz

*   Y =: Z.x.X 
      e(1,2) =  e(2,3)*e(3,1) - e(3,3)*e(2,1)
      e(2,2) = -e(1,3)*e(3,1) + e(3,3)*e(1,1)
      e(3,2) =  e(1,3)*e(2,1) - e(2,3)*e(1,1)


      p(1) = x*e(1,1) + y*e(1,2) + z*e(1,3) + a(1)
      p(2) = x*e(2,1) + y*e(2,2) + z*e(2,3) + a(2)
      p(3) = x*e(3,1) + y*e(3,2) + z*e(3,3) + a(3)

      end     
