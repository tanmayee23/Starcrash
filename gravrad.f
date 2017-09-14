**************************************************************
c     This file is part of the StarCrash code
c     Version 1.0
c
c     Copyright 2003, by Joshua Faber
c     This code is publicly available under the terms of the
c     GNU General Public License.  See our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

      SUBROUTINE QGRAVRAD
      IMPLICIT NONE
*****************************************************
c     Release 1.0
C     Calculates radiation reaction potential by FFTs on a grid.
c     For more info on the equation we are solving, see
c     Eq. 5.14 of Blanchet, Damour, and Schaefer MNRAS 242,289 (1990).
c
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
      LOGICAL OUT(N),CHANGENGR,xfix,yfix,zfix 
      INTEGER plan1(4),plan2(4),plan3(4),plan4(4),plan5(4)
      DOUBLE PRECISION  RHOG(NNgrav,NNgrav,nnp),
     $     GREEN(NNgrav,NNgrav,nnp),
     $     PHI(NNgrav,NNgrav,nnp)
      DOUBLE COMPLEX RHOGT(n1grav,NNgrav,nnp),
     $     GreenT(n1grav,NNgrav,nnp),
     $     PHIT(n1grav,NNgrav,nnp)
      REAL FGRX(Ngrav,Ngrav,nnp),FGRY(Ngrav,Ngrav,nnp)
      REAL RHO1(Ngrav,Ngrav,nnp),rzup(Ngrav,Ngrav,nnp)
      REAL RZDN(Ngrav,Ngrav,nnp)
      REAL myGx(N),myGy(N),myGz(N)

      REAL pm,amtot,amxtot,dxmin,dxmax,xgrmax2,xgrmin2
      REAL ygrmax2,ygrmin2,zgrmax2,zgrmin2,fnp,hx,hy,hz
      REAL xmin,xmax,ymin,ymax,zmin,zmax,xc,yc,zc,rc2
      REAL xp,yp,zp,ain,ajn,akn,scalefactor,hx2,hy2,hz2,rc
      REAL yj,yjp1,yjm1,xi,xip1,xim1,zk
      INTEGER kn2mod,knm1mod
      INTEGER INIT,ip,ix,ixp,i,j,k,k1,i1,j1,ierr,npout
      INTEGER in,jn,kn,knproc,kn1proc,knmod,kn1mod,in1,jn1,kn1
      INTEGER knc,kn1c,kn2proc,kn2c,knm1proc,knm1c

      INTEGER NAMX
      REAL TGR
      PARAMETER (NAMX=2000,TGR=5.0)
      REAL AMX(NAMX)
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
      
C     If this is the first call, or if mesh has been adjusted,
C     initialize FFT and calculate transform of Green's function:

      if(init.eq.0) then
         if(myrank.eq.0)WRITE (6,*) 'QGRAV: INITIALIZING GT'
         call makeplans(plan1,plan2,plan3,plan4,plan5,nngrav)
      endif

      IF ((INIT.EQ.0).OR.(NGR.GT.1)) THEN
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
         do j=1,ngrav
            do i=1,ngrav
               rho1(i,j,k)=0.
               rzup(i,j,k)=0.
               rzdn(i,j,k)=0.
            enddo
         enddo
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
         kn2proc=(kn+1)/nnp
         kn2mod=mod(kn+1,nnp)+1
         knm1proc=(kn-2)/nnp
         knm1mod=mod(kn-2,nnp)+1
         IN1=IN+1
         JN1=JN+1
         KN1=KN+1

c     This long and complicated method allows us to calculate
c     d/dx rho in a way that conserves the total integrated source as 0.

         IF ((IN.GT.1).AND.(IN1.LT.N1grav).AND.
     $        (JN.GT.1).AND.(JN1.LT.N1grav).AND.
     $        (KN.GT.1).AND.(KN1.LT.N1grav)) THEN
C     (particle is inside the mesh)
            OUT(IP)=.FALSE.
            if(myrank.eq.knproc) then
               RHO1(IN,JN,knmod)=RHO1(IN,JN,knmod)+
     $              (IN1-AIN)*(JN1-AJN)*(KN1-AKN)*AM(IP)
               RHO1(IN,JN1,knmod)=RHO1(IN,JN1,knmod)+
     $              (IN1-AIN)*(AJN-JN)*(KN1-AKN)*AM(IP)
               RHO1(IN1,JN,knmod)=RHO1(IN1,JN,knmod)+
     $              (AIN-IN)*(JN1-AJN)*(KN1-AKN)*AM(IP)
               RHO1(IN1,JN1,knmod)=RHO1(IN1,JN1,knmod)+
     $              (AIN-IN)*(AJN-JN)*(KN1-AKN)*AM(IP)

               RZDN(IN,JN,knmod)=RZDN(IN,JN,knmod)+
     $              (IN1-AIN)*(JN1-AJN)*(AKN-KN)*AM(IP)
               RZDN(IN,JN1,knmod)=RZDN(IN,JN1,knmod)+
     $              (IN1-AIN)*(AJN-JN)*(AKN-KN)*AM(IP)
               RZDN(IN1,JN,knmod)=RZDN(IN1,JN,knmod)+
     $              (AIN-IN)*(JN1-AJN)*(AKN-KN)*AM(IP)
               RZDN(IN1,JN1,knmod)=RZDN(IN1,JN1,knmod)+
     $              (AIN-IN)*(AJN-JN)*(AKN-KN)*AM(IP)

            endif
            if(myrank.eq.kn1proc) then
               RHO1(IN,JN,kn1mod)=RHO1(IN,JN,kn1mod)+
     $              (IN1-AIN)*(JN1-AJN)*(AKN-KN)*AM(IP)
               RHO1(IN,JN1,kn1mod)=RHO1(IN,JN1,kn1mod)+
     $              (IN1-AIN)*(AJN-JN)*(AKN-KN)*AM(IP)
               RHO1(IN1,JN,kn1mod)=RHO1(IN1,JN,kn1mod)+
     $              (AIN-IN)*(JN1-AJN)*(AKN-KN)*AM(IP)
               RHO1(IN1,JN1,kn1mod)=RHO1(IN1,JN1,kn1mod)+
     $              (AIN-IN)*(AJN-JN)*(AKN-KN)*AM(IP)

               RZUP(IN,JN,kn1mod)=RZUP(IN,JN,kn1mod)+
     $              (IN1-AIN)*(JN1-AJN)*(KN1-AKN)*AM(IP)
               RZUP(IN,JN1,kn1mod)=RZUP(IN,JN1,kn1mod)+
     $              (IN1-AIN)*(AJN-JN)*(KN1-AKN)*AM(IP)
               RZUP(IN1,JN,kn1mod)=RZUP(IN1,JN,kn1mod)+
     $              (AIN-IN)*(JN1-AJN)*(KN1-AKN)*AM(IP)
               RZUP(IN1,JN1,kn1mod)=RZUP(IN1,JN1,kn1mod)+
     $              (AIN-IN)*(AJN-JN)*(KN1-AKN)*AM(IP)

            endif
            if(myrank.eq.kn2proc) then
               RZUP(IN,JN,kn2mod)=RZUP(IN,JN,kn2mod)+
     $              (IN1-AIN)*(JN1-AJN)*(AKN-KN)*AM(IP)
               RZUP(IN,JN1,kn2mod)=RZUP(IN,JN1,kn2mod)+
     $              (IN1-AIN)*(AJN-JN)*(AKN-KN)*AM(IP)
               RZUP(IN1,JN,kn2mod)=RZUP(IN1,JN,kn2mod)+
     $              (AIN-IN)*(JN1-AJN)*(AKN-KN)*AM(IP)
               RZUP(IN1,JN1,kn2mod)=RZUP(IN1,JN1,kn2mod)+
     $              (AIN-IN)*(AJN-JN)*(AKN-KN)*AM(IP)
            endif
            if(myrank.eq.knm1proc) then
               RZDN(IN,JN,knm1mod)=RZDN(IN,JN,knm1mod)+
     $              (IN1-AIN)*(JN1-AJN)*(KN1-AKN)*AM(IP)
               RZDN(IN,JN1,knm1mod)=RZDN(IN,JN1,knm1mod)+
     $              (IN1-AIN)*(AJN-JN)*(KN1-AKN)*AM(IP)
               RZDN(IN1,JN,knm1mod)=RZDN(IN1,JN,knm1mod)+
     $              (AIN-IN)*(JN1-AJN)*(KN1-AKN)*AM(IP)
               RZDN(IN1,JN1,knm1mod)=RZDN(IN1,JN1,knm1mod)+
     $              (AIN-IN)*(AJN-JN)*(KN1-AKN)*AM(IP)
            endif
         ELSE
            OUT(IP)=.TRUE.
            NPOUT=NPOUT+1
         ENDIF
      ENDDO

      do k=1,nnp
         zk=(k-1+nnp*myrank)*hz+zmin
         do j=2,ngrav
            yj=(j-1)*hy+ymin
            yjp1=j*hy+ymin
            yjm1=(j-2)*hy+ymin
            do i=2,ngrav
               xi=(i-1)*hx+xmin
               xip1=i*hx+xmin
               xim1=(i-2)*hx+xmin

               rhog(i-1,j,k)=rhog(i-1,j,k)+rho1(i,j,k)/(2.0*hx)*
     $              (xim1*q3xx+yj*q3xy+zk*q3xz)
               rhog(i+1,j,k)=rhog(i+1,j,k)-rho1(i,j,k)/(2.0*hx)*
     $              (xip1*q3xx+yj*q3xy+zk*q3xz)
               rhog(i,j-1,k)=rhog(i,j-1,k)+rho1(i,j,k)/(2.0*hy)*
     $              (xi*q3xy+yjm1*q3yy+zk*q3yz)
               rhog(i,j+1,k)=rhog(i,j+1,k)-rho1(i,j,k)/(2.0*hy)*
     $              (xi*q3xy+yjp1*q3yy+zk*q3yz)
               rhog(i,j,k)=rhog(i,j,k)+(rzdn(i,j,k)-rzup(i,j,k))
     $              /(2.0*hz)*(xi*q3xz+yj*q3yz+zk*q3zz)
            enddo
         enddo
      enddo

C     Renormalize density values for FFT:
      DO K=1,nnp
         DO J=1,Ngrav
            DO I=1,Ngrav
               RHOG(I,J,K)=RHOG(I,J,K)*4.*PI
            ENDDO
         ENDDO
      ENDDO

      call rcfft(plan1,plan3,plan4,nnp,nngrav,rhog,rhogt)

C     Do the convolution:
      DO K=1,nnp
         DO J=1,nngrav
            DO I=1,nngrav
               PHIt(I,J,K)=RHOGt(I,J,K)*Greent(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      
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

      do ip=1,n
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
            IF (myrank.eq.knproc) then
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
               myGx(IP)=0.0
               myGy(IP)=0.0
               myGz(IP)=0.0
            endif
            if(myrank.eq.kn1proc) then
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

      CALL MPI_ALLREDUCE(myGx,dxpnr,N,MPI_REAL,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(myGy,dypnr,N,MPI_REAL,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(myGz,dzpnr,N,MPI_REAL,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)

      IF(CHANGENGR .AND.NGR.EQ.21) THEN
         NGR=1
         if(myrank.eq.0)WRITE(6,*) ' QGRAV: NGR HAS CHANGED TO',NGR
      ENDIF

      RETURN
      END
***********************************************************************








