**************************************************************
c     This file is part of the StarCrash code
c     Version 1.0
c
c     Copyright 2003, by Joshua Faber
c     This code is publicly available under the terms of the
c     GNU General Public License.  See our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

      SUBROUTINE setup1es                                                
      IMPLICIT NONE
********************************************************
c     Release 1.0
C     Sets up a single NS of index 1/(GAM-1), thus P=a*rho^gam
c     represented by a uniform distribution of ~N particles with 
c     varying masses (HP(I)=const) in a hexagonal close-packed lattice.
c     Note that the actual number of particles will differ from N, generally
c     slightly lower.
c     Called by INIT
c     Calls POLY,GRAVQUANT,LFSTART
**********************************************************
      INCLUDE 'spha.h'

      INTEGER nrgrid
      PARAMETER(NRGRID=10000)                                         
      REAL RHOPOL(NRGRID),avec(3),bvec(3),cvec(3)

      REAL ak,space,xrand,xp,yp,zp,r,cden,ammin,ammax
      REAL amtot,ri,rhoi,hc,amtot2,amxt,amyt,amzt
      INTEGER maxn,i,idum,k,l,m

      REAL ran1
 
C Get parameters from input file:
      OPEN(12,FILE='sph.input',ERR=100)
      READ(12,INPUT)
      CLOSE(12)

      if(n.gt.nmax) then
         if(myrank.eq.0)write(6,*)'ERROR: N>NMAX'
         stop
      endif

C Get density and pressure profiles:                                  
      CALL POLY(1./(GAM-1.),amns,rns,NRGRID,RHOPOL,AK)                 

C Lay down particles uniformly within first of two stars:             
C    (use hexagonal lattice)                             

      avec(1)=1.0
      avec(2)=0.0
      avec(3)=0.0
      bvec(1)=0.5
      bvec(2)=0.5*sqrt(3.0)
      bvec(3)=0.0
      cvec(1)=0.5
      cvec(2)=1.0/sqrt(12.0)
      cvec(3)=sqrt(2.0/3.0)

      space=rns/(1.0*n/2.0/6.0)**0.333
c The factor 1.23 makes sure to catch all points

      maxn=int(1.23*rns/space)
      if(myrank.eq.0)write (6,*)'1es: maxn: ',maxn,' spacing: ',space
      i=0
      idum=-2391
      xrand=0.001*space

      do k=-maxn,maxn
         do l=-maxn,maxn
            do m=-maxn,maxn
               xp=space*(k*avec(1)+l*bvec(1)+m*cvec(1))
               yp=space*(k*avec(2)+l*bvec(2)+m*cvec(2))
               zp=space*(k*avec(3)+l*bvec(3)+m*cvec(3))
               xp=xp+xrand*(-1.+2.*RAN1(IDUM))
               yp=yp+xrand*(-1.+2.*RAN1(IDUM))
               zp=zp+xrand*(-1.+2.*RAN1(IDUM))
               r=sqrt(xp**2+yp**2+zp**2)
               if (r.lt.rns) then
                  i=i+1
                  if(i.le.nmax) then
                     x(i)=xp
                     y(i)=yp
                     z(i)=zp
                  else
                     if(myrank.eq.0)write (6,*)'1es: Reached nmax!'
                     stop
                  endif
               endif
            enddo
         enddo
      enddo
      n=i
      if(myrank.eq.0)write (6,*)'number of particles within rns: ',i

C Assign particle masses (to represent density):                      
      CDEN=(4.*PI*rns**3)/(3.*FLOAT(N))                                      
      AMMIN=1.E10                                                     
      AMMAX=0.                                                        
      AMTOT=0.                                                        
      DO I=1,N                                                        
        RI=SQRT(X(I)**2+Y(I)**2+Z(I)**2)/rns                              
        RHOI=RHOPOL(1+INT(RI*FLOAT(NRGRID-1)))                        
        AM(I)=CDEN*RHOI                                               
        AMMIN=MIN(AMMIN,AM(I))                                        
        AMMAX=MAX(AMMAX,AM(I))                                        
        AMTOT=AMTOT+AM(I)                                             
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        A(I)=AK
      ENDDO                                                           
      if(myrank.eq.0)
     $     WRITE (6,*) '1es: TOTAL MASS WAS',AMTOT,' changing to ',amns
                      
C Renormalize total mass and adjust CoM:

      amtot2=0.
      amxt=0.
      amyt=0.
      amzt=0.
      DO I=1,N                                                        
        AM(I)=AM(I)*amns/AMTOT                                             
        amtot2=amtot2+am(i)
        amxt=amxt+am(i)*x(i)
        amyt=amyt+am(i)*y(i)
        amzt=amzt+am(i)*z(i)
      ENDDO                                                           
      do i=1,n
         x(i)=x(i)-amxt/amtot2
         y(i)=y(i)-amyt/amtot2
         z(i)=z(i)-amzt/amtot2
      enddo
      AMMIN=AMMIN*amns/AMTOT                                               
      AMMAX=AMMAX*amns/AMTOT                                               
      if(myrank.eq.0) then
         WRITE (6,*) 
     $        '1es: MINIMUM PARTICLE MASS=',AMMIN              
         WRITE (6,*) 
     $        '1es: MAXIMUM PARTICLE MASS=',AMMAX              
      endif                                                  
      DO I=1,N                                                        
        VX(I)=0.                                                      
        VY(I)=0.                                                      
        VZ(I)=0.                                                      
      ENDDO                                                      
      NLEFT=N
                                                                      
C Set kernel support (constant kernel width)                     
      HC=rns*(FLOAT(NNOPT)/FLOAT(4*N))**0.3333                            
      if(myrank.eq.0)WRITE (6,*) '1es: H=',HC
C     (should give a number of nearest neighbors close to NNOPT)      
      DO I=1,N                                                        
        HP(I)=HC                                                      
      ENDDO                                                           

      if(myrank.eq.0)write(6,*)'1es:nrelax=',nrelax,' trelax=',trelax
C prepare leap-frog scheme for first iteration                       
      call gravquant
      CALL LFSTART                                                    
                                                                      
      RETURN                                                          
C     Error condition:
 100  STOP '1es: ERROR OPENING INPUT FILE ???'
      END              

************************************************************************
      SUBROUTINE setup1em                                                
      IMPLICIT NONE
********************************************************
c     Release 1.0
C     Sets up a single NS of index 1/(GAM-1), thus P=a*rho^gam,
c     represented by equal-mass particles laid down Monte-Carlo
c     style weighted by the polytropic density profile.
c     In general, exactly N particles should be used.
c     Called by INIT
c     Calls POLY,GRAVQUANT,OPTHP2,LFSTART
**********************************************************
      INCLUDE 'spha.h'

      INTEGER nrgrid,maxtry
      PARAMETER(MAXTRY=10000000,NRGRID=10000)         
      REAL RHOPOL(NRGRID),ak,rhomax,xtry,ytry,ztry,rtry,rhoex
      REAL rhotry,cden,amtot,amtot2,amxt,amyt,amzt
      integer idum,ip,ntry,i

      REAL ran1

C Get parameters from input file:
      OPEN(12,FILE='sph.input',ERR=100)
      READ(12,INPUT)
      CLOSE(12)

      if(n.gt.nmax) then
         if(myrank.eq.0)write(6,*)'ERROR: N>NMAX'
         stop
      endif

C Get density and pressure profiles:                                  
      CALL POLY(1./(GAM-1.),amns,rns,NRGRID,RHOPOL,AK)                 

      IDUM=-2391         
      RHOMAX=RHOPOL(1)   
      IP=0               
      DO NTRY=1,MAXTRY   
         XTRY=-1.+2.*RAN1(IDUM)  
         YTRY=-1.+2.*RAN1(IDUM)  
         ZTRY=-1.+2.*RAN1(IDUM)
         RTRY=SQRT(XTRY**2+YTRY**2+ZTRY**2) 
         IF (RTRY.LT.1.) THEN          
            RHOEX=RHOPOL(1+INT(RTRY*FLOAT(NRGRID-1)))
            RHOTRY=RHOMAX*RAN1(IDUM)      
            IF (RHOTRY.LT.RHOEX) THEN     
C     (particle is accepted)    
               IP=IP+1                 
               X(IP)=XTRY*rns       
               Y(IP)=YTRY*rns        
               Z(IP)=ZTRY*rns            
               A(IP)=AK
            ENDIF                   
         ENDIF                     
         IF (IP.EQ.N) goto 1
      ENDDO
      if(myrank.eq.0)WRITE(6,*) 'ONLY GOT TO IP=',IP
      STOP '1em: NOT ENOUGH PARTICLES ???'
 1    CONTINUE                
      
C     Assign all equal particle masses:
      CDEN=amns/FLOAT(N) 
      amtot=0
      DO I=1,N                     
         AM(I)=CDEN                 
         amtot=amtot+am(i)
      ENDDO                        
      if(myrank.eq.0)
     $     WRITE (6,*) '1em: TOTAL MASS WAS',AMTOT,' changing to ',amns
      
      amtot2=0.
      amxt=0.
      amyt=0.
      amzt=0.
      do i=1,n
         am(i)=am(i)*amns/amtot
         amtot2=amtot2+am(i)
         amxt=amxt+am(i)*x(i)
         amyt=amyt+am(i)*y(i)
         amzt=amzt+am(i)*z(i)
      ENDDO                                                           
      do i=1,n
         x(i)=x(i)-amxt/amtot2
         y(i)=y(i)-amyt/amtot2
         z(i)=z(i)-amzt/amtot2
      enddo
      
      NLEFT=N

      DO I=1,N                                                        
        VX(I)=0.                                                      
        VY(I)=0.                                                      
        VZ(I)=0.                                                      
      ENDDO                                                      
C There are now 2N particles (N for star on left) 
      NLEFT=N

c     set kernel support:
      call gravquant
      CALL OPTHP2

      if(myrank.eq.0)write(6,*)'1es:nrelax=',nrelax,' trelax=',trelax
C prepare leap-frog scheme for first iteration                       
      CALL LFSTART                                                    
                                                                      
      RETURN                                                          
C     Error condition:
 100  STOP '1em: ERROR OPENING INPUT FILE ???'
      END              

******************************************************************
      SUBROUTINE OPTHP2
      PARAMETER(NRGRID=10000)
      REAL RHOPOL(NRGRID)   
************************************************************
c     Release 1.0
C     For a spherical distribution of matter, this routine determines 
C     the kernel widths HP(I) of all particles, in such a way that   
C     the numbers of nearest neighbors NN(I) are all approximately NNOPT. 
c     Called by SETUP1EM,SETUP2EM
*********************************************************************
      INCLUDE 'spha.h' 

C get density and pressure profiles: 
      CALL POLY(1./(gam-1.),amns,rns,NRGRID,RHOPOL,AK) 
C     (constructs polytrope of total mass=1 and radius=1)
      if(myrank.eq.0)WRITE (6,*) 'OPTHP2: INITIALIZING HP...'   
                                                
      CONST=(3./32./3.14159265358979)**(1./3.)
      DO IP=1,N
         RIP=SQRT(X(IP)**2+Y(IP)**2+Z(IP)**2)/rns
         ANUMDEN=RHOPOL(1+INT(RIP*FLOAT(NRGRID-1)))*N/amns
         HP(IP)=0.9*CONST*(NNOPT/ANUMDEN)**(1./3.)*rns
      ENDDO
      do nloop=1,10
         HPMIN=1.E30                             
         HPMAX=0.                                
         CALL NENE
         CALL ADJUST
         DO I=1,N                                
            HPMIN=MIN(HPMIN,HP(I))                
            HPMAX=MAX(HPMAX,HP(I))                
         ENDDO                                  
         if(myrank.eq.0)WRITE (6,*) 'OPTHP2: HPMIN=',HPMIN       
         if(myrank.eq.0)WRITE (6,*) '       HPMAX=',HPMAX       
      enddo
                                              
      HPMIN=1.E30                             
      HPMAX=0.                                
      DO I=1,N                                
        HPMIN=MIN(HPMIN,HP(I))                
        HPMAX=MAX(HPMAX,HP(I))                
      ENDDO                                  
      if(myrank.eq.0)WRITE (6,*) 'OPTHP2: HPMIN=',HPMIN       
      if(myrank.eq.0)WRITE (6,*) '       HPMAX=',HPMAX       
                                       
      RETURN                           
      END
