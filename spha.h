**************************************************************
c     This file is part of the StarCrash code
c     Version 1.0
c
c     Copyright 2003, by Joshua Faber
c     This code is publicly available under the terms of the
c     GNU General Public License.  See our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

C     NUMBER OF PROCESSORS
      INTEGER NPROCS
      PARAMETER (NPROCS=2)
c     Particle Parameters
      INTEGER NMAX,NMAXP,NNMAX
      PARAMETER (NMAX=210000,NMAXP=(NMAX-1)/NPROCS+1,NNMAX=200)
c     Interpolation Parameters for kernel summation
      INTEGER NTAB
      PARAMETER (NTAB=100000)                                           
c     Real constants
      REAL CN,PI
      PARAMETER (CN=0.15,PI=3.14159265359)
c     Gravity grid parameters
      INTEGER NGRAV,NNGRAV,N1GRAV,NNP
      PARAMETER (NGRAV=64,NNGRAV=NGRAV*2)
      PARAMETER (N1GRAV=NGRAV+1,NNP=NNGRAV/NPROCS)
c     Energy evolution eqn. choice
      INTEGER NENERGY
      PARAMETER (NENERGY=1)

c     ARTVIS: artificial viscosity quantities
      INTEGER NAV
      REAL ALPHA,BETA,ETA2 
C     GRAV: Gravity grid values
      INTEGER NGR
      REAL XGRMIN,XGRMAX,YGRMIN,YGRMAX,ZGRMIN,ZGRMAX,
     $             XGRLIM,YGRLIM,ZGRLIM
C     RELAXP: Relaxation values and related quantities
      INTEGER NRELAX
      REAL TRELAX,TRELOFF,OMEGA2
C     BINARYP: Binary parameters
      INTEGER NLEFT
      REAL SEP0,QDAR,AMNS,RNS,RP,VPEAK
C     INTPAR: Input parameters
      INTEGER N,NNOPT,NTIMESTEPPER
      REAL GAM,HMIN,HMAX
C     OUTP: Output parameters
      INTEGER NOUT,NIT
      REAL DT,T,VXS(NMAX),VYS(NMAX),VZS(NMAX),TF,DTOUT,
     $             XM2(NMAX),YM2(NMAX),ZM2(NMAX)
c     WTABUL: tabulated kernel functions
      REAL WTAB(NTAB),DWTAB(NTAB),CTAB
c     PART: particle-based quantities
      REAL X(NMAX),Y(NMAX),Z(NMAX),VX(NMAX),VY(NMAX),VZ(NMAX),
     $     AM(NMAX),HP(NMAX),A(NMAX),RHO(NMAX),POR2(NMAX),GRPOT(NMAX),
     $     EN(NMAX),U(NMAX)
      INTEGER NREF(NMAX)
c     DYN: dynamical quantities 
      REAL VXDOT(NMAX),VYDOT(NMAX),VZDOT(NMAX),ADOT(NMAX),
     $      GX(NMAX),GY(NMAX),GZ(NMAX),ux(NMAX),uy(NMAX),uz(NMAX),
     $      ENDOT(NMAX),UDOT(NMAX)
C     NEIGH: Neighbor list quantities 
      REAL XIJ(NNMAX,NMAXP),YIJ(NNMAX,NMAXP),ZIJ(NNMAX,NMAXP)
      INTEGER NNI(NNMAX,NMAXP),NN(NMAX)
C     PARALLEL: parallelization parameters
      INTEGER myrank,n_lower,n_upper,kstart,koffset
C     q3: 3rd time derivative of quadrupole 
      REAL q3xx,q3xy,q3xz,q3yy,q3yz,q3zz
C     q2: 2nd time derivative of quadrupole 
      REAL q2xx,q2xy,q2xz,q2yy,q2yz,q2zz
C     gravrad: quantities for computing gravitational radiation losses 
      INTEGER ngravrad,initgr
      REAL sol,freacx(NMAX),freacy(NMAX),freacz(NMAX),gxx(NMAX),
     $      gxy(NMAX),gyy(NMAX),gxz(NMAX),gyz(NMAX),gzz(NMAX),
     $      dxpnr(NMAX),dypnr(NMAX),dzpnr(NMAX)

      COMMON/ARTVIS/ NAV,ALPHA,BETA,ETA2
      COMMON/GRAV/ NGR,XGRMIN,XGRMAX,YGRMIN,YGRMAX,ZGRMIN,
     $             ZGRMAX,XGRLIM,YGRLIM,ZGRLIM
      COMMON/RELAXP/ NRELAX,TRELAX,TRELOFF,OMEGA2
      COMMON/BINARYP/ NLEFT,SEP0,QDAR,AMNS,RNS,RP,VPEAK
      COMMON/INTPAR/ GAM,N,NNOPT,NTIMESTEPPER,HMIN,HMAX
      COMMON/OUTP/ DT,T,VXS,VYS,VZS,TF,DTOUT,NOUT,NIT,xm2,ym2,zm2     
      COMMON/WTABUL/ WTAB,DWTAB,CTAB                        
      COMMON/PART/X,Y,Z,VX,VY,VZ,AM,HP,A,RHO,POR2,GRPOT,EN,NREF,U  
      COMMON/DYN/ VXDOT,VYDOT,VZDOT,ADOT,GX,GY,GZ,ux,uy,uz,EnDOT,UDOT
      COMMON/NEIGH/ XIJ,YIJ,ZIJ,NNI,NN
      COMMON/PARALLEL/ myrank,n_lower,n_upper,kstart,koffset
      COMMON/q3/ q3xx,q3xy,q3xz,q3yy,q3yz,q3zz
      COMMON/q2/ q2xx,q2xy,q2xz,q2yy,q2yz,q2zz
      COMMON/gravrad/ ngravrad,initgr,sol,
     $      freacx,freacy,freacz,gxx,gxy,gyy,gxz,gyz,gzz,
     $      dxpnr,dypnr,dzpnr
      NAMELIST/INPUT/ TF,DTOUT,GAM,N,NNOPT,NAV,ALPHA,BETA,ETA2,
     $      NGR,XGRMIN,XGRMAX,YGRMIN,YGRMAX,ZGRMIN,
     $      ZGRMAX,XGRLIM,YGRLIM,ZGRLIM,HMIN,HMAX,NRELAX,TRELAX,SEP0,
     $      QDAR,AMNS,RNS,RP,VPEAK,TRELOFF,NGRAVRAD,SOL,NTIMESTEPPER
cccccccc     $      NENERGY
