**************************************************************
c     This file is part of the StarCrash code
c     Version 1.0
c
c     Copyright 2003, by Joshua Faber
c     This code is publicly available under the terms of the
c     GNU General Public License.  See our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

      SUBROUTINE QGRAV
      IMPLICIT NONE
*****************************************************
c     Release 1.0
C     Calculates gravitational forces by FFTs on a grid.
C     Cloud in cell is used to set up 
C     the density on the grid and interpolate the forces back to
C     particles in this version.
C     Forces on particles outside the grid are calculated by simple
C     monopole approximation. 
C     This is the parallel version of the code, using MPI to 
c     beat the performance of serial routines.  As such, the FFT has
c     been replaced by the FFTW (Fastest Fourier Transform in the West),
c     developed by Steven Johnson and Matteo Frigo at MIT. For more info,
c     please check http://www.fftw.org
c     called by BALVDOTS,CLAVDOTS,NEWVDOTS
**********************************************************
      include 'spha.h'
      include 'mpif.h'
      INTEGER NAMX
      REAL TGR
      PARAMETER (NAMX=2000,TGR=5.0)
      REAL AMX(NAMX)

      LOGICAL OUT(N),CHANGENGR,xfix,yfix,zfix 
      INTEGER plan1(4),plan2(4),plan3(4),plan4(4),plan5(4)
      DOUBLE PRECISION  RHOG(NNgrav,NNgrav,nnp),
     $     GREEN(NNgrav,NNgrav,nnp),
     $     PHI(NNgrav,NNgrav,nnp)
      DOUBLE COMPLEX RHOGT(n1grav,NNgrav,nnp),
     $     GreenT(n1grav,NNgrav,nnp),
     $     PHIT(n1grav,NNgrav,nnp)
      REAL FGRX(Ngrav,Ngrav,nnp),FGRY(Ngrav,Ngrav,nnp)
      REAL mygrpot(N),myGx(N),myGy(N),myGz(N)

      REAL pm,amtot,amxtot,dxmin,dxmax,xgrmax2,xgrmin2
      REAL ygrmax2,ygrmin2,zgrmax2,zgrmin2,fnp,hx,hy,hz,rip3
      REAL xmin,xmax,ymin,ymax,zmin,zmax,xc,yc,zc,rc2,rip
      REAL xp,yp,zp,ain,ajn,akn,scalefactor,hx2,hy2,hz2,rc
      REAL amin,xcmin,ycmin,zcmin,grpottot,gxtot,gytot,gztot
      INTEGER INIT,ip,ix,ixp,i,j,k,k1,i1,j1,ierr,npout
      INTEGER in,jn,kn,knproc,kn1proc,knmod,kn1mod,in1,jn1,kn1
      INTEGER knc,kn1c,kn2proc,kn2c,knm1proc,knm1c

      DATA INIT /0/
      SAVE INIT,GreenT,plan1,plan2,plan3,plan4,plan5

      amtot=0.
      CHANGENGR=.FALSE.
      xfix=.false.
      yfix=.false.
      zfix=.false.
      
c      IF (NPART.GT.NPMAX) STOP 'QGRAV: NPART>NPMAX ???'
      
C     Calculate boundaries of mesh according to NGR:
C     NGR=0  ->  no gravity
C     NGR=1  ->  fixed mesh (boundaries provided by calling unit)
C     NGR=2  ->  mesh boundaries recalculated to cover all particles
C     NGR=21 ->  behave initially like NGR=2, then switch to NGR=1
c     after time TGR has passed
C     NGR=P with 90 <= P <= 99  -> Like NGR=1 at first, then 
c     mesh boundaries recalculated
C     to cover P% of the mass in all 3 directions once a particle
c     leaves the box defined by XGRLIM,YGRLIM,ZGRLIM
C     NGR=P with 290 <= P <= 299  -> mesh boundaries recalculated here
C     initially to cover all particles and then to cover (P-200)% of the
C     mass once a single particle leaves the box defined by XGRLIM,
C     YGRLIM and ZGRLIM
c     NGR=398 -> Like NGR=98 at first, but the box size is limited such
c     to size 2*XGRLIM,2*YGRLIM, and 2*ZGRLIM.

      IF (NGR.GT.1) THEN
         XGRMAX=-1.E30
         YGRMAX=-1.E30
         ZGRMAX=-1.E30
         XGRMIN=1.E30
         YGRMIN=1.E30
         ZGRMIN=1.E30
         DO IP=1,N
            XGRMAX=MAX(X(IP),XGRMAX)
            YGRMAX=MAX(Y(IP),YGRMAX)
            ZGRMAX=MAX(Z(IP),ZGRMAX)
            XGRMIN=MIN(X(IP),XGRMIN)
            YGRMIN=MIN(Y(IP),YGRMIN)
            ZGRMIN=MIN(Z(IP),ZGRMIN)
         ENDDO
      ENDIF
      
      IF (NGR.EQ.21 .AND. T.GT.TGR) THEN
         XGRMAX=XGRLIM
         YGRMAX=YGRLIM
         ZGRMAX=ZGRLIM
         XGRMIN=-XGRLIM
         YGRMIN=-YGRLIM
         ZGRMIN=-ZGRLIM
         CHANGENGR=.TRUE.
      ENDIF

      IF (NGR.GT.89) THEN
         IF(XGRMAX.GT.XGRLIM .OR. XGRMIN.LT.-XGRLIM)
     $        xfix=.true.
         if(YGRMAX.GT.YGRLIM .OR. YGRMIN.LT.-YGRLIM)
     $        yfix=.true.
         if(ZGRMAX.GT.ZGRLIM .OR. ZGRMIN.LT.-ZGRLIM)
     $        zfix=.true.
         if((xfix.or.yfix.or.zfix).and.ngr.gt.200.and.
     $        ngr.lt.300) then
            NGR=NGR-200
            if(myrank.eq.0)WRITE(6,*)
     $           ' QGRAV: NGR HAS CHANGED TO',NGR
         ENDIF
      ENDIF
      
      IF ((NGR.GT.89 .AND. NGR.LT.100).or.ngr.eq.398)THEN
         PM=0.5*(1.-FLOAT(mod(NGR,100))/1.E2)
C     Adjust x-boundaries:
C     Calculate mass distribution in x:
         if (xfix) then
            DO IX=1,NAMX
               AMX(IX)=0.
            ENDDO
            AMTOT=0.
            DO IP=1,N
               IXP=INT((X(IP)-XGRMIN)*FLOAT(NAMX)/(XGRMAX-XGRMIN))
               IXP=MIN(IXP+1,NAMX)
               AMX(IXP)=AMX(IXP)+AM(IP)
               AMTOT=AMTOT+AM(IP)
            ENDDO
C     Scan from left:
            AMXTOT=0.
            DO IX=1,NAMX
               AMXTOT=AMXTOT+AMX(IX)
               IF (AMXTOT/AMTOT.GT.PM) GOTO 11
            ENDDO
 11         DXMIN=FLOAT(IX-1)*(XGRMAX-XGRMIN)/FLOAT(NAMX)
C     Scan from right:
            AMXTOT=0.
            DO IX=NAMX,1,-1
               AMXTOT=AMXTOT+AMX(IX)
               IF (AMXTOT/AMTOT.GT.PM) GOTO 12
            ENDDO
 12         DXMAX=FLOAT(IX)*(XGRMAX-XGRMIN)/FLOAT(NAMX)
C     Adjust boundaries:
            XGRMAX2=XGRMIN+DXMAX
            XGRMIN2=XGRMIN+DXMIN
            if(XGRMAX.gt.xgrlim) then
               if(xgrmax2.gt.xgrlim) then
                  xgrmax=xgrmax2
               else
                  xgrmax=xgrlim
               endif
            endif
            if(xgrmin.lt.-xgrlim) then
               if(xgrmin2.lt.-xgrlim) then
                  xgrmin=xgrmin2
               else
                  xgrmin=-xgrlim
               endif
            endif
         endif
C     Adjust y-boundaries:
C     Calculate mass distribution in y:
         if (yfix) then
            DO IX=1,NAMX
               AMX(IX)=0.
            ENDDO
            DO IP=1,N
               IXP=INT((Y(IP)-YGRMIN)*FLOAT(NAMX)/(YGRMAX-YGRMIN))
               IXP=MIN(IXP+1,NAMX)
               AMX(IXP)=AMX(IXP)+AM(IP)
            ENDDO
C     Scan from left:
            AMXTOT=0.
            DO IX=1,NAMX
               AMXTOT=AMXTOT+AMX(IX)
               IF (AMXTOT/AMTOT.GT.PM) GOTO 21
            ENDDO
 21         DXMIN=FLOAT(IX-1)*(YGRMAX-YGRMIN)/FLOAT(NAMX)
C     Scan from right:
            AMXTOT=0.
            DO IX=NAMX,1,-1
               AMXTOT=AMXTOT+AMX(IX)
               IF (AMXTOT/AMTOT.GT.PM) GOTO 22
            ENDDO
 22         DXMAX=FLOAT(IX)*(YGRMAX-YGRMIN)/FLOAT(NAMX)
C     Adjust boundaries:
            YGRMAX2=YGRMIN+DXMAX
            YGRMIN2=YGRMIN+DXMIN
            if(YGRMAX.gt.ygrlim) then
               if(ygrmax2.gt.ygrlim) then
                  ygrmax=ygrmax2
               else
                  ygrmax=ygrlim
               endif
            endif
            if(ygrmin.lt.-ygrlim) then
               if(ygrmin2.lt.-ygrlim) then
                  ygrmin=ygrmin2
               else
                  ygrmin=-ygrlim
               endif
            endif
         endif
C     Adjust z-boundaries:
C     Calculate mass distribution in z:
         if (zfix) then
            DO IX=1,NAMX
               AMX(IX)=0.
            ENDDO
            DO IP=1,N
               IXP=INT((Z(IP)-ZGRMIN)*FLOAT(NAMX)/(ZGRMAX-ZGRMIN))
               IXP=MIN(IXP+1,NAMX)
               AMX(IXP)=AMX(IXP)+AM(IP)
            ENDDO
C     Scan from left:
            AMXTOT=0.
            DO IX=1,NAMX
               AMXTOT=AMXTOT+AMX(IX)
               IF (AMXTOT/AMTOT.GT.PM) GOTO 31
            ENDDO
 31         DXMIN=FLOAT(IX-1)*(ZGRMAX-ZGRMIN)/FLOAT(NAMX)
C     Scan from right:
            AMXTOT=0.
            DO IX=NAMX,1,-1
               AMXTOT=AMXTOT+AMX(IX)
               IF (AMXTOT/AMTOT.GT.PM) GOTO 32
            ENDDO
 32         DXMAX=FLOAT(IX)*(ZGRMAX-ZGRMIN)/FLOAT(NAMX)
C     Adjust boundaries:
            ZGRMAX2=ZGRMIN+DXMAX
            ZGRMIN2=ZGRMIN+DXMIN
            if(ZGRMAX.gt.zgrlim) then
               if(zgrmax2.gt.zgrlim) then
                  zgrmax=zgrmax2
               else
                  zgrmax=zgrlim
               endif
            endif
            if(zgrmin.lt.-zgrlim) then
               if(zgrmin2.lt.-zgrlim) then
                  zgrmin=zgrmin2
               else
                  zgrmin=-zgrlim
               endif
            endif
         endif
      ENDIF

      if(ngr.eq.398) then
         xgrmin=max(xgrmin,-2.0*xgrlim)
         xgrmax=min(xgrmax,2.0*xgrlim)
         ygrmin=max(ygrmin,-2.0*ygrlim)
         ygrmax=min(ygrmax,2.0*ygrlim)
         zgrmin=max(zgrmin,-2.0*zgrlim)
         zgrmax=min(zgrmax,2.0*zgrlim)
      endif

C     Correct size and calculate grid cell widths so that one layer of
C     empty cells surrounds the mesh on all sides (this is necessary to
C     do the force assignment):
      FNP=1./( FLOAT(NGRAV)-2.2)
      HX=(XGRMAX-XGRMIN)*FNP
      HY=(YGRMAX-YGRMIN)*FNP
      HZ=(ZGRMAX-ZGRMIN)*FNP
      XMIN=XGRMIN-HX*1.1
      XMAX=XGRMAX+HX*1.1
      YMIN=YGRMIN-HY*1.1
      YMAX=YGRMAX+HY*1.1
      ZMIN=ZGRMIN-HZ*1.1
      ZMAX=ZGRMAX+HZ*1.1

      IF (NGR.GT.1.and.myrank.eq.0) WRITE (6,101) XGRMIN,XGRMAX,
     $     YGRMIN,YGRMAX,ZGRMIN,ZGRMAX
 101  FORMAT (' QGRAV: MESH HAS BEEN ADJUSTED:',/,
     $     '   XGRMIN   ','   XGRMAX   ','   YGRMIN   ',
     $     '   YGRMAX   ','   ZGRMIN   ','   ZGRMAX   ',/,
     $     6E12.4)   
      
c     initialize the FFTW

      if(init.eq.0) then
         if(myrank.eq.0)WRITE (6,*) 'QGRAV: INITIALIZING GT'
         call makeplans(plan1,plan2,plan3,plan4,plan5,nngrav)
      endif

C     If this is the first call, or if mesh has been adjusted,
C     calculate transform of Green's function:

      IF ((INIT.EQ.0).OR.(NGR.GT.1)) THEN
         if(myrank.eq.0)WRITE (6,*) 'QGRAV: Calculating GREENT'
         DO K=1,NNp
            DO J=1,nngrav
               DO I=1,nngrav
                  GREEN(I,J,K)=0.
               ENDDO
            ENDDO
         ENDDO
         
         INIT=1
C     Calculate values of Green's function on the mesh
C     ...first on inner mesh (where G=1/r):
         DO K=1,nnp
            k1=k-1+kstart
            ZC=FLOAT(K1-1)*HZ
            DO J=1,NGRAV+1
               YC=FLOAT(J-1)*HY
               DO I=1,NGRAV+1
                  XC=FLOAT(I-1)*HX
                  RC=SQRT(XC**2+YC**2+ZC**2)+1.E-10
                  GREEN(I,J,K)=-1./(4.*PI*RC)
               ENDDO
            ENDDO
         ENDDO
************************************************
C     this is only for the first process!!!!
         if(myrank.eq.0) then
            rc2=sqrt((hx**2+hy**2+hz**2)/3.0)
            GREEN(1,1,1)=-1.0/(4.*pi*rc2)*4.0
         endif
C
************************************************
C     ... then on extended mesh (by symmetries):
         DO K=1,nnp
            DO J=1,nngrav
               DO I=1,nngrav

C***************************************************
C     The reflection in the z-axis is no longer necessary

                  IF (I.GT.NGRAV+1) THEN
                     I1=NNGRAV-I+2
                  ELSE
                     I1=I
                  ENDIF
                  IF (J.GT.NGRAV+1) THEN
                     J1=NNGRAV-J+2
                  ELSE
                     J1=J
                  ENDIF
                  GREEN(I,J,K)=GREEN(I1,J1,K)
               ENDDO
            ENDDO
         ENDDO

         call rcfft(plan1,plan3,plan4,nnp,nngrav,green,greent)
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ENDIF

************************************************************************
C     Initialize density values:
      DO K=1,nnp
         DO J=1,nngrav
            DO I=1,nngrav
               RHOG(I,J,K)=0.
            ENDDO
         ENDDO
      ENDDO
      
C     Calculate new density values by cloud in cell scheme:
      NPOUT=0
      DO IP=1,N
         XP=X(IP)-XMIN
         YP=Y(IP)-YMIN
         ZP=Z(IP)-ZMIN
         AIN=XP/HX+1.
         AJN=YP/HY+1.
         AKN=ZP/HZ+1.
         IN=INT(AIN)
         JN=INT(AJN)
         KN=INT(AKN)
         knproc=(kn-1)/nnp
         kn1proc=kn/nnp
         knmod=mod(kn-1,nnp)+1
         kn1mod=mod(kn,nnp)+1
         IN1=IN+1
         JN1=JN+1
         KN1=KN+1
         IF ((IN.GT.0).AND.(IN1.LT.N1grav).AND.
     $        (JN.GT.0).AND.(JN1.LT.N1grav).AND.
     $        (KN.GT.0).AND.(KN1.LT.N1grav)) THEN
C     (particle is inside the mesh)
            OUT(IP)=.FALSE.
            if(myrank.eq.knproc) then
               RHOG(IN,JN,knmod)=RHOG(IN,JN,knmod)+
     $              (IN1-AIN)*(JN1-AJN)*(KN1-AKN)*AM(IP)
               RHOG(IN,JN1,knmod)=RHOG(IN,JN1,knmod)+
     $              (IN1-AIN)*(AJN-JN)*(KN1-AKN)*AM(IP)
               RHOG(IN1,JN,knmod)=RHOG(IN1,JN,knmod)+
     $              (AIN-IN)*(JN1-AJN)*(KN1-AKN)*AM(IP)
               RHOG(IN1,JN1,knmod)=RHOG(IN1,JN1,knmod)+
     $              (AIN-IN)*(AJN-JN)*(KN1-AKN)*AM(IP)
            endif
            if(myrank.eq.kn1proc) then
               RHOG(IN,JN,kn1mod)=RHOG(IN,JN,kn1mod)+
     $              (IN1-AIN)*(JN1-AJN)*(AKN-KN)*AM(IP)
               RHOG(IN,JN1,kn1mod)=RHOG(IN,JN1,kn1mod)+
     $              (IN1-AIN)*(AJN-JN)*(AKN-KN)*AM(IP)
               RHOG(IN1,JN,kn1mod)=RHOG(IN1,JN,kn1mod)+
     $              (AIN-IN)*(JN1-AJN)*(AKN-KN)*AM(IP)
               RHOG(IN1,JN1,kn1mod)=RHOG(IN1,JN1,kn1mod)+
     $              (AIN-IN)*(AJN-JN)*(AKN-KN)*AM(IP)
            endif
         ELSE
            OUT(IP)=.TRUE.
            NPOUT=NPOUT+1
         ENDIF
      ENDDO
      
C     Renormalize density values for FFT:
      DO K=1,nnp
         DO J=1,Ngrav
            DO I=1,Ngrav
               RHOG(I,J,K)=RHOG(I,J,K)*4.*PI
            ENDDO
         ENDDO
      ENDDO

c     forward transform
      call rcfft(plan1,plan3,plan4,nnp,nngrav,rhog,rhogt)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

C     Do the convolution:
      DO K=1,nnp
         DO J=1,nngrav
            DO I=1,nngrav
               PHIt(I,J,K)=RHOGt(I,J,K)*Greent(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      
c     inverse transform
      call crfft(plan2,plan3,plan5,nnp,nngrav,phit,phi)
      SCALEFACTOR=1.0/FLOAT(NNgrav)**3
      do k=1,nnp
         do j=1,nngrav
            do i=1,nngrav
               phi(i,j,k)=phi(i,j,k)*scalefactor
            enddo
         enddo
      enddo

C     Calculate force values on the mesh:

      do k=1,nnp
         do j=1,ngrav
            do i=1,ngrav
               fgrx(i,j,k)=0.0
               fgry(i,j,k)=0.0
            enddo
         enddo
      enddo
      HX2=2.*HX
      HY2=2.*HY
      HZ2=2.*HZ
      DO K=1,nnp
         DO J=2,Ngrav
            DO I=2,Ngrav
               FGRX(I,J,k)=(PHI(I-1,J,k)-PHI(I+1,J,k))/HX2
               FGRY(I,J,k)=(PHI(I,J-1,K)-PHI(I,J+1,K))/HY2
            ENDDO
         ENDDO
      ENDDO
C     Calculate forces and potential at particle positions:
C     Calculate forces for particles inside grid by finite diff:
      DO IP=1,N
         mygrpot(ip)=0.0
         mygx(ip)=0.0
         mygy(ip)=0.0
         mygz(ip)=0.0
      enddo
      DO IP=1,N
         XP=X(IP)-XMIN
         YP=Y(IP)-YMIN
         ZP=Z(IP)-ZMIN
         AIN=XP/HX+1.
         AJN=YP/HY+1.
         AKN=ZP/HZ+1.
         IN=INT(AIN)
         JN=INT(AJN)
         KN=INT(AKN)
         knproc=(kn-1)/nnp
         knc=mod((kn-1),nnp)+1
         IN1=IN+1
         JN1=JN+1
         KN1=KN+1
         kn1proc=(kn/nnp)
         kn1c=mod(kn,nnp)+1
         kn2proc=(kn+1)/nnp
         kn2c=mod((kn+1),nnp)+1
         knm1proc=(kn-2)/nnp
         knm1c=mod((kn-2),nnp)+1
         IF (.NOT.OUT(IP)) THEN
c     The data is distributed over processors in the z-direction, we
c     must keep careful track of which data values are where.  Finite
c     differencing requires a slightly different treatment in the 
c     z-direction
            IF (myrank.eq.knproc) then
               mygrpot(IP)= (KN1-AKN)*(
     $              (JN1-AJN)*((IN1-AIN)*PHI(IN,JN,knc)+
     $              (AIN-IN)*PHI(IN1,JN,knc))+
     $              (AJN-JN)*((IN1-AIN)*PHI(IN,JN1,knc)+
     $              (AIN-IN)*PHI(IN1,JN1,knc)))
               myGx(IP)= (KN1-AKN)*(
     $              (JN1-AJN)*((IN1-AIN)*FGRX(IN,JN,knc)+
     $              (AIN-IN)*FGRX(IN1,JN,knc))+
     $              (AJN-JN)*((IN1-AIN)*FGRX(IN,JN1,knc)+
     $              (AIN-IN)*FGRX(IN1,JN1,knc)))
               myGy(IP)= (KN1-AKN)*(
     $              (JN1-AJN)*((IN1-AIN)*FGRY(IN,JN,KNc)+
     $              (AIN-IN)*FGRY(IN1,JN,KNc))+
     $              (AJN-JN)*((IN1-AIN)*FGRY(IN,JN1,KNc)+
     $              (AIN-IN)*FGRY(IN1,JN1,KNc)))
               myGz(IP)=(aKN-KN)/hz2*(
     $              (JN1-AJN)*((IN1-AIN)*phi(IN,JN,KNc)+ 
     $              (AIN-IN)*phi(IN1,JN,KNc))+
     $              (AJN-JN)*((IN1-AIN)*phi(IN,JN1,KNc)+
     $              (AIN-IN)*phi(IN1,JN1,KNc)))
            else
               mygrpot(IP)=0.0
               myGx(IP)=0.0
               myGy(IP)=0.0
               myGz(IP)=0.0
            endif
            if(myrank.eq.kn1proc) then
               mygrpot(IP)=mygrpot(IP)+(AKN-KN)*(
     $              (JN1-AJN)*((IN1-AIN)*PHI(IN,JN,kn1c)+
     $              (AIN-IN)*PHI(IN1,JN,kn1c))+
     $              (AJN-JN)*((IN1-AIN)*PHI(IN,JN1,kn1c)+
     $              (AIN-IN)*PHI(IN1,JN1,kn1c)))
               myGx(IP)=myGx(IP)+(AKN-KN)*(
     $              (JN1-AJN)*((IN1-AIN)*FGRX(IN,JN,KN1c)+
     $              (AIN-IN)*FGRX(IN1,JN,KN1c))+
     $              (AJN-JN)*((IN1-AIN)*FGRX(IN,JN1,KN1c)+
     $              (AIN-IN)*FGRX(IN1,JN1,KN1c)))
               myGy(IP)=myGy(IP)+(AKN-kn)*(
     $              (JN1-AJN)*((IN1-AIN)*FGRY(IN,JN,KN1c)+
     $              (AIN-IN)*FGRY(IN1,JN,KN1c))+
     $              (AJN-JN)*((IN1-AIN)*FGRY(IN,JN1,KN1c)+
     $              (AIN-IN)*FGRY(IN1,JN1,KN1c)))
               myGz(IP)=myGz(IP)+(AKN-KN1)/hz2*(
     $              (JN1-AJN)*((IN1-AIN)*phi(IN,JN,KN1c)+
     $              (AIN-IN)*phi(IN1,JN,KN1c))+
     $              (AJN-JN)*((IN1-AIN)*phi(IN,JN1,KN1c)+
     $              (AIN-IN)*phi(IN1,JN1,KN1c)))
            endif
            if(myrank.eq.kn2proc) then
               myGz(IP)=myGz(IP)+(KN-AKN)/hz2*(
     $              (JN1-AJN)*((IN1-AIN)*phi(IN,JN,KN2c)+
     $              (AIN-IN)*phi(IN1,JN,KN2c))+
     $              (AJN-JN)*((IN1-AIN)*phi(IN,JN1,KN2c)+
     $              (AIN-IN)*phi(IN1,JN1,KN2c)))
            endif
            if(myrank.eq.knm1proc) then
               myGz(IP)=myGz(IP)+(kn1-akn)/hz2*(
     $              (JN1-AJN)*((IN1-AIN)*phi(IN,JN,KNm1c)+
     $              (AIN-IN)*phi(IN1,JN,KNm1c))+
     $              (AJN-JN)*((IN1-AIN)*phi(IN,JN1,KNm1c)+
     $              (AIN-IN)*phi(IN1,JN1,KNm1c)))
            endif
         endif
      enddo

      CALL MPI_ALLREDUCE(mygrpot,grpot,N,MPI_REAL,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(myGx,Gx,N,MPI_REAL,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(myGy,Gy,N,MPI_REAL,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(myGz,Gz,N,MPI_REAL,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)

      IF(CHANGENGR .AND.NGR.EQ.21) THEN
         NGR=1
         if(myrank.eq.0)WRITE(6,*) ' QGRAV: NGR HAS CHANGED TO',NGR
      ENDIF

C     Exit if no particle is outside the grid:
      IF (NPOUT.EQ.0) THEN
         RETURN
      ENDIF

C     Calculate forces for particles outside grid by monopole approx:
      AMIN=0.
      XCMIN=0.
      YCMIN=0.
      ZCMIN=0.

      DO IP=1,N
         IF (.NOT.OUT(IP)) THEN
            AMIN=AMIN+AM(IP)
            XCMIN=XCMIN+AM(IP)*X(IP)
            YCMIN=YCMIN+AM(IP)*Y(IP)
            ZCMIN=ZCMIN+AM(IP)*Z(IP)
         ENDIF
      ENDDO
      XCMIN=XCMIN/AMIN
      YCMIN=YCMIN/AMIN
      ZCMIN=ZCMIN/AMIN
      if(myrank.eq.0)WRITE (6,*) 
     $     'QGRAV: NPOUT=',NPOUT,'  AMIN=',AMIN

      IF(NPOUT.GT.N/5) THEN
         if(myrank.eq.0)WRITE(6,*) 
     $        'WARNING!!!!! In qgrav:NPOUT.GT.N/5  ???'
C         STOP
      ENDIF

      grpotTOT=0.
      GxTOT=0.
      GyTOT=0.
      GzTOT=0.
      DO IP=1,N
         IF (OUT(IP)) THEN
            RIP=((X(IP)-XCMIN)**2+(Y(IP)-YCMIN)**2
     $           +(Z(IP)-ZCMIN)**2)**0.5
            RIP3=RIP**3
            grpot(IP)=-AMIN/RIP
            grpotTOT=grpotTOT-AM(IP)/RIP
            Gx(IP)=AMIN*(XCMIN-X(IP))/RIP3
            Gy(IP)=AMIN*(YCMIN-Y(IP))/RIP3
            Gz(IP)=AMIN*(ZCMIN-Z(IP))/RIP3
            GxTOT=GxTOT-AM(IP)*(XCMIN-X(IP))/RIP3
            GyTOT=GyTOT-AM(IP)*(YCMIN-Y(IP))/RIP3
            GzTOT=GzTOT-AM(IP)*(ZCMIN-Z(IP))/RIP3
         ENDIF
      ENDDO
      DO IP=1,N
         IF(.NOT.OUT(IP)) THEN
            grpot(IP)=grpot(IP)+grpotTOT
            Gx(IP)=Gx(IP)+GxTOT
            Gy(IP)=Gy(IP)+GyTOT
            Gz(IP)=Gz(IP)+GzTOT
         ENDIF
      ENDDO

      RETURN
      END
***********************************************************************








