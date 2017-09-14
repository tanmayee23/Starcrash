**************************************************************
c     This file is part of the StarCrash code
c     Version 1.0
c
c     Copyright 2003, by Joshua Faber
c     This code is publicly available under the terms of the
c     GNU General Public License.  See our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

      SUBROUTINE POLY(AN,AM,R,NR,RHO,AK)
*****************************************************
c     Release 1.0
C     Calculates specific entropy A for polytrope with index N,
c     mass AM, radius R, and gives density RHO at NR points in radius
c     This routine does not "include spha.h", thus rho is unambiguous
c     Called by SETUP1ES,SETUP1EM,SETUP2CS,SETUP2CM
************************************************** 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (XMAX=7.D0,MAXSTP=15,PI=3.14159265359D0,NRM=20000)
      REAL*4 AN,AM,R,RHO(NR),AK
      DIMENSION XTAB(NRM),YTAB(NRM),YPTAB(NRM)
      COMMON/POLIDX/ANDBL

      ANDBL=AN

      XTAB(1)=0.
      YTAB(1)=1.
      YPTAB(1)=0.
      XS=XMAX
      DO 20 NSTP=1,MAXSTP
C advance solution to r=dr using series expansion near origin:
       DR=XS/DBLE(NR-1)
       Y1=1.D0-DR**2/6.D0+ANDBL*DR**4/120.D0
       YP1=-DR/3.D0+ANDBL*DR**3/30.D0
       X1=DR
       CALL RKTAB(Y1,YP1,X1,XS,NR-1,XTAB(2),YTAB(2),YPTAB(2))
       DO 10 IR=2,NR
10      IF (YTAB(IR).LT.0.D0) GOTO 11
11     XS=XTAB(IR)
       YS=YTAB(IR)
       YPS=YPTAB(IR)
       IF (IR.EQ.NR) GOTO 21
20    CONTINUE
      STOP 'POLY: NO CONVERGENCE ???'
21    CONTINUE
      YTAB(NR)=0.D0

      RHOC=XS*DBLE(AM)/(4.D0*PI*DABS(YPS)*DBLE(R)**3)
      AK=SNGL(4.D0*PI*DBLE(R)**2*RHOC**(1.D0-1.D0/ANDBL)/((ANDBL+1.D0)
     +                                                        *XS**2))
      DO 30 IR=1,NR
30     RHO(IR)=SNGL(RHOC*YTAB(IR)**ANDBL)

      RETURN
      END
************************************************************************
      SUBROUTINE DERIVS(X,V,DV)
      DOUBLE PRECISION X,V(2),DV(2),AN
      COMMON/POLIDX/AN
      DV(1)=V(2)
      DV(2)=-2.D0*V(2)/X-DABS(V(1))**AN
      RETURN
      END
************************************************************************
      SUBROUTINE RKTAB(Y,YP,X1,X2,NTAB,XTAB,YTAB,YPTAB)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION XTAB(NTAB),YTAB(NTAB),YPTAB(NTAB)
      DIMENSION V(2),DV(2)

      XTAB(1)=X1
      YTAB(1)=Y
      YPTAB(1)=YP
      X=X1
      H=(X2-X1)/DBLE(NTAB-1)
      V(1)=Y
      V(2)=YP
      DO 13 K=1,NTAB-1
        CALL DERIVS(X,V,DV)
        CALL RK4(V,DV,2,X,H,V)
        IF (X+H.EQ.X) STOP 'RKTAB: STEPSIZE NOT SIGNIFICANT ???'
        X=X+H
        XTAB(K+1)=X
        YTAB(K+1)=V(1)
        YPTAB(K+1)=V(2)
13    CONTINUE
      RETURN
      END
************************************************************************
C This is the standard 4th order Runge-Kutta algorithm:
      SUBROUTINE RK4(Y,DYDX,N,X,H,YOUT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NMAX=10)
      DIMENSION Y(N),DYDX(N),YOUT(N),YT(NMAX),DYT(NMAX),DYM(NMAX)

      IF (N.GT.NMAX) STOP 'RK4: N>NMAX ???'

      HH=H*0.5D0
      H6=H/6.D0
      XH=X+HH
      DO 11 I=1,N
        YT(I)=Y(I)+HH*DYDX(I)
11    CONTINUE
      CALL DERIVS(XH,YT,DYT)
      DO 12 I=1,N
        YT(I)=Y(I)+HH*DYT(I)
12    CONTINUE
      CALL DERIVS(XH,YT,DYM)
      DO 13 I=1,N
        YT(I)=Y(I)+H*DYM(I)
        DYM(I)=DYT(I)+DYM(I)
13    CONTINUE
      CALL DERIVS(X+H,YT,DYT)
      DO 14 I=1,N
        YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.D0*DYM(I))
14    CONTINUE
      RETURN
      END
