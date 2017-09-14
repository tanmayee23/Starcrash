**************************************************************
c     This file is part of the StarCrash code
c     Version 1.0
c
c     Copyright 2003, by Joshua Faber
c     This code is publicly available under the terms of the
c     GNU General Public License.  See our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

      SUBROUTINE setup2qs                                                
      IMPLICIT NONE
********************************************************
c     Release 1.0
C     Sets up a synchronized binary with unequal-mass NS of 
c     mass ratio q<1.0, both w/index 1/(GAM-1), thus P=a*rho^gam
c     represented by a uniform distribution of ~N particles per star with 
c     varying masses (HP(I)=const) in a hexagonal close-packed lattice.
c     Note that the actual number of particles will differ from 2*N, 
c     generally slightly lower.
c     Called by INIT
c     Calls POLY,GRAVQUANT,LFSTART
**********************************************************
      INCLUDE 'spha.h'

      INTEGER nrgrid
      PARAMETER(NRGRID=10000)                                         
      REAL RHOPOL(NRGRID),avec(3),bvec(3),cvec(3)

      REAL ak,space,xrand,xp,yp,zp,r,cden,ammin,ammax
      REAL amtot,ri,rhoi,hc,al,amns2,rns2,dum
      INTEGER maxn,i,idum,k,l,m

      REAL ran1
 
C Get parameters from input file:
      OPEN(12,FILE='sph.input',ERR=100)
      READ(12,INPUT)
      CLOSE(12)

      amns2=qdar*amns
      if(gam.ge.4.0/3.0) then
         rns2=rns*qdar**((gam-2.0)/(3.0*gam-4.0))
      else
         write(6,*)'illegal gamma of',gam
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
      if(myrank.eq.0)write (6,*)'2qs: maxn: ',maxn,' spacing: ',space
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
                     if(myrank.eq.0)write (6,*)'2qs: Reached nmax!'
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
     $     WRITE (6,*) '2qs: Star 1 MASS WAS',AMTOT,' changing to ',amns
                      
C Renormalize total mass:                                             
      DO I=1,N                                                        
        AM(I)=AM(I)*amns/AMTOT                                             
      ENDDO                                                           
      AMMIN=AMMIN*amns/AMTOT                                               
      AMMAX=AMMAX*amns/AMTOT                                               
      if(myrank.eq.0) then
         WRITE (6,*) 
     $        '2qs: star 1 MINIMUM PARTICLE MASS=',AMMIN              
         WRITE (6,*) 
     $        '2qs: star 1 MAXIMUM PARTICLE MASS=',AMMAX              
      endif                                                  
      DO I=1,N                                                        
        VX(I)=0.                                                      
        VY(I)=0.                                                      
        VZ(I)=0.                                                      
      ENDDO                                                      
C There are now 2N particles (N for star on left) 
      NLEFT=N

      maxn=int(1.23*rns2/space)
      if(myrank.eq.0)write (6,*)'2qs: maxn: ',maxn,' spacing: ',space
      i=nleft

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
               if (r.lt.rns2) then
                  i=i+1
                  if(i.le.nmax) then
                     x(i)=xp
                     y(i)=yp
                     z(i)=zp
                  else
                     if(myrank.eq.0)write (6,*)'2qs: Reached nmax!'
                     stop
                  endif
               endif
            enddo
         enddo
      enddo
      n=i
      if(myrank.eq.0)write (6,*)'total number of particles: ',i

      CALL POLY(1./(GAM-1.),amns2,rns2,NRGRID,RHOPOL,dum)                 

C Assign particle masses (to represent density):                      
      CDEN=(4.*PI*rns2**3)/(3.*FLOAT(N-nleft))
      AMMIN=1.E10                                                     
      AMMAX=0.                                                        
      AMTOT=0.                                                        
      DO I=nleft+1,N                                                        
        RI=SQRT(X(I)**2+Y(I)**2+Z(I)**2)/rns2                              
        RHOI=RHOPOL(1+INT(RI*FLOAT(NRGRID-1)))                        
        AM(I)=CDEN*RHOI                                               
        AMMIN=MIN(AMMIN,AM(I))                                        
        AMMAX=MAX(AMMAX,AM(I))                                        
        AMTOT=AMTOT+AM(I)                                             
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        A(I)=AK
      ENDDO                                                           
      if(myrank.eq.0)
     $     WRITE (6,*) '2qs:star 2 MASS WAS',AMTOT,' changing to ',amns2
                      
C Renormalize total mass:                                             
      DO I=nleft+1,N                                                        
        AM(I)=AM(I)*amns2/AMTOT                                             
      ENDDO                                                           
      AMMIN=AMMIN*amns2/AMTOT                                               
      AMMAX=AMMAX*amns2/AMTOT                                               
      if(myrank.eq.0) then
         WRITE (6,*) 
     $        '2qs: star 2 MINIMUM PARTICLE MASS=',AMMIN              
         WRITE (6,*) 
     $        '2qs: star 2 MAXIMUM PARTICLE MASS=',AMMAX              
      endif                                                  

      DO I=1,N                                                        
        VX(I)=0.                                                      
        VY(I)=0.                                                      
        VZ(I)=0.                                                      
      ENDDO                                                      
                                                                      
C Set kernel support (constant kernel width)                     
      HC=rns*(FLOAT(NNOPT)/FLOAT(4*Nleft))**0.3333                            
      if(myrank.eq.0)WRITE (6,*) '2qs: H=',HC
C     (should give a number of nearest neighbors close to NNOPT)      
      DO I=1,N                                                        
        HP(I)=HC                                                      
      ENDDO                                                           

      if(n.gt.NMAX) STOP 'SETUPDAR: 2N>NMAX ???'  
      if(myrank.eq.0) write(6,*)'NLEFT=',nleft,'  N=',n
      call gravquant
      
C Shift stars along x-axis so that CM is at the origin
C  and separation = SEP0
      AL=qdar*sep0/(1.0+qdar)
      if(myrank.eq.0) write(6,*)'SEP0=',sep0,'  AL=',al
      DO I=1,Nleft                                                    
         X(I)=X(I)-AL
      enddo
      al=sep0/(1.0+qdar)
      do i=nleft+1,n
         X(i)=x(i)+AL
      ENDDO                                                           

      if(myrank.eq.0)write(6,*)'2qs:nrelax=',nrelax,' trelax=',trelax
C prepare leap-frog scheme for first iteration                       
      CALL LFSTART                                                    
                                                                      
      RETURN                                                          
C     Error condition:
 100  STOP '2qs: ERROR OPENING INPUT FILE ???'
      END              

************************************************************************
      SUBROUTINE setup2qm                                                
      IMPLICIT NONE
********************************************************
c     Release 1.0
C     Sets up a synchronized binary with unequal-mass NS of mass ratio
c     q<1.0, both with index 1/(GAM-1), thus P=a*rho^gam,
c     represented by equal-mass particles laid down Monte-Carlo
c     style weighted by the polytropic density profile.
c     Called by INIT
c     Calls POLY,GRAVQUANT,OPTHP2,LFSTART
**********************************************************

      INCLUDE 'spha.h'

      INTEGER nrgrid,maxtry
      PARAMETER(MAXTRY=10000000,NRGRID=10000)         
      REAL RHOPOL(NRGRID),ak,rhomax,xtry,ytry,ztry,rtry,rhoex
      REAL rhotry,cden,amtot,al,amns2,rns2,dum
      integer idum,ip,ntry,i,n2

      REAL ran1

C Get parameters from input file:
      OPEN(12,FILE='sph.input',ERR=100)
      READ(12,INPUT)
      CLOSE(12)

      n2=int(qdar*n)
      amns2=qdar*amns
      if(gam.ge.4.0/3.0) then
         rns2=rns*qdar**((gam-2.0)/(3.0*gam-4.0))
      else
         write(6,*)'illegal gamma of',gam
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
C         (particle is accepted)    
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
      STOP '2cm: NOT ENOUGH PARTICLES ???'
1     CONTINUE                

C Assign all equal particle masses:
      CDEN=amns/FLOAT(N) 
      amtot=0
      DO I=1,N                     
        AM(I)=CDEN                 
        amtot=amtot+am(i)
      ENDDO                        
      if(myrank.eq.0)
     $     WRITE (6,*) '2cm: TOTAL MASS WAS',AMTOT,' changing to ',amns
      do i=1,n
         am(i)=am(i)*amns/amtot
      enddo

      NLEFT=N
      n=nleft+n2

c     set kernel support:
      call gravquant
      CALL OPTHP2

      CALL POLY(1./(GAM-1.),amns2,rns2,NRGRID,RHOPOL,dum)                 

      RHOMAX=RHOPOL(1)   
      IP=nleft               
      DO NTRY=1,MAXTRY   
        XTRY=-1.+2.*RAN1(IDUM)  
        YTRY=-1.+2.*RAN1(IDUM)  
        ZTRY=-1.+2.*RAN1(IDUM)
        RTRY=SQRT(XTRY**2+YTRY**2+ZTRY**2) 
        IF (RTRY.LT.1.) THEN          
          RHOEX=RHOPOL(1+INT(RTRY*FLOAT(NRGRID-1)))
          RHOTRY=RHOMAX*RAN1(IDUM)      
          IF (RHOTRY.LT.RHOEX) THEN     
C         (particle is accepted)    
            IP=IP+1                 
            X(IP)=XTRY*rns2       
            Y(IP)=YTRY*rns2        
            Z(IP)=ZTRY*rns2            
            A(IP)=AK
          ENDIF                   
        ENDIF                     
        IF (IP.EQ.N) goto 2
      ENDDO
      if(myrank.eq.0)WRITE(6,*) 'ONLY GOT TO IP=',IP
      STOP '2cm: NOT ENOUGH PARTICLES ???'
 2    CONTINUE                

C Assign all equal particle masses:
      CDEN=amns2/FLOAT(N2) 
      amtot=0
      DO I=nleft+1,N                     
        AM(I)=CDEN                 
        amtot=amtot+am(i)
      ENDDO                        
      if(myrank.eq.0)
     $     WRITE (6,*)'2cm:TOTAL MASS WAS',AMTOT,' changing to ',amns2
      do i=nleft+1,n
         am(i)=am(i)*amns2/amtot
      enddo

c     set kernel support:
      call gravquant
      CALL OPTHP3

C There are now 2N particles (N for star on left) 
      if(n.gt.NMAX) STOP '2qm: N>NMAX ???'  
      if(myrank.eq.0) write(6,*)'NLEFT=',nleft,'  N=',n
      call gravquant
      
C Shift stars along x-axis so that CM is at the origin
C  and separation = SEP0
      AL=qdar*sep0/(1.0+qdar)
      if(myrank.eq.0) write(6,*)'SEP0=',sep0,'  AL=',al
      DO I=1,Nleft                                                    
         X(I)=X(I)-AL
      enddo
      AL=sep0/(1.0+qdar)
      do i=nleft+1,n
         X(i)=x(i)+AL
      ENDDO                                                           

C Set all velocities equal to zero
      DO IP=1,N
         VX(IP)=0.                                                      
         VY(IP)=0.                                                      
         VZ(IP)=0.                                                      
      ENDDO                                                           

      if(myrank.eq.0)write(6,*)'2cm:nrelax=',nrelax,' trelax=',trelax
C prepare leap-frog scheme for first iteration                       
      CALL LFSTART                                                    
                                                                      
      RETURN                                                          
C     Error condition:
 100  STOP '2cm: ERROR OPENING INPUT FILE ???'
      END              

*********************************************************************8
      SUBROUTINE OPTHP3
      PARAMETER(NRGRID=10000)
      REAL RHOPOL(NRGRID)
************************************************************
c     Release 1.0
C     For a spherical distribution of matter, this routine determines
C     the kernel widths HP(I) of all particles in the second star, 
c     as OPTHP2 does for the primary, in such a way that
C     the numbers of nearest neighbors NN(I) are all approximately NNOPT.
c     Called by SETUP1EM,SETUP2EM
*********************************************************************
      INCLUDE 'spha.h'

      n2=int(qdar*n)
      amns2=qdar*amns
      if(gam.ge.4.0/3.0) then
         rns2=rns*qdar**((gam-2.0)/(3.0*gam-4.0))
      else
         write(6,*)'illegal gamma of',gam
         stop
      endif

C get density and pressure profiles:
      CALL POLY(1./(gam-1.),amns2,rns2,NRGRID,RHOPOL,AK)
C     (constructs polytrope of total mass=1 and radius=1)
      if(myrank.eq.0)WRITE (6,*) 'OPTHP3: INITIALIZING HP...'

      CONST=(3./32./3.14159265358979)**(1./3.)
      DO IP=nleft+1,N
         RIP=SQRT(X(IP)**2+Y(IP)**2+Z(IP)**2)/rns2
         ANUMDEN=RHOPOL(1+INT(RIP*FLOAT(NRGRID-1)))*N/amns2
         HP(IP)=0.9*CONST*(NNOPT/ANUMDEN)**(1./3.)
      ENDDO
      do i=1,5
         CALL NENE
         CALL ADJUST
      enddo

      HPMIN=1.E30
      HPMAX=0.
      DO I=1,N
        HPMIN=MIN(HPMIN,HP(I))
        HPMAX=MAX(HPMAX,HP(I))
      ENDDO
      if(myrank.eq.0)WRITE (6,*) 'OPTHP3: HPMIN=',HPMIN
      if(myrank.eq.0)WRITE (6,*) '       HPMAX=',HPMAX

      RETURN
      END

********************************************************************
      SUBROUTINE setup2qr                                                
      IMPLICIT NONE
********************************************************
c     Release 1.0
C     Sets up a synchronized binary with unequal-mass NS,
c     using previously calculated relaxed single-star models,
c     processed by b2pos and b2pos2 into files called "pos.sph"
c     and "pos2.sph", respectively.
c     Called by INIT
c     Calls GRAVQUANT,LFSTART
**********************************************************
      INCLUDE 'spha.h'

      INTEGER nold,nnoptold,i,ip,n2old
      REAL hminold,hmaxold,gamold,trelaxold,al,amtot1,amtot2

C Get parameters from input file:
      OPEN(12,FILE='sph.input',ERR=100)
      READ(12,INPUT)
      CLOSE(12)

      amtot1=0.
      amtot2=0.

      open(13,file='pos.sph', form='formatted')
      read(13,'(2i6,4e15.7)')nold,nnoptold,
     $     hminold,hmaxold,gamold,trelaxold

      if(2*nold.gt.nmax) then
         if(myrank.eq.0)write(6,*)'ERROR: 2N>NMAX'
         stop
      endif
      if(nnoptold.ne.nnopt)write(6,*)'star1:nnopt changed from ',
     $     nnoptold,' to ',nnopt
      if(hminold.ne.hmin)write(6,*)'star1:hmin changed from ',
     $     hminold,' to ',hmin
      if(hmaxold.ne.hmax)write(6,*)'star1:hmax changed from ',
     $     hmaxold,' to ',hmax
      if(gamold.ne.gam)write(6,*)'star1:gam changed from ',
     $     gamold,' to ',gam
      if(trelaxold.ne.trelax)write(6,*)'star1:trelax changed from ',
     $     trelaxold,' to ',trelax
      
      DO I=1,Nold
         READ (13,'(7e15.7)') X(i),Y(i),z(i),am(i),hp(i),rho(i),a(i)
         amtot1=amtot1+am(i)
      ENDDO
      close(13)

C There are now 2N particles (N for star on left) 
      NLEFT=N
                                                                      
      open(13,file='pos2.sph', form='formatted')
      read(13,'(2i6,4e15.7)')n2old,nnoptold,
     $     hminold,hmaxold,gamold,trelaxold

      n=nleft+n2old

      if(nold+n2old.gt.nmax) then
         if(myrank.eq.0)write(6,*)'ERROR: 2N>NMAX'
         stop
      endif
      if(nnoptold.ne.nnopt)write(6,*)'star2:nnopt changed from ',
     $     nnoptold,' to ',nnopt
      if(hminold.ne.hmin)write(6,*)'star2:hmin changed from ',
     $     hminold,' to ',hmin
      if(hmaxold.ne.hmax)write(6,*)'star2:hmax changed from ',
     $     hmaxold,' to ',hmax
      if(gamold.ne.gam)write(6,*)'star2:gam changed from ',
     $     gamold,' to ',gam
      if(trelaxold.ne.trelax)write(6,*)'star2:trelax changed from ',
     $     trelaxold,' to ',trelax
      
      DO I=nleft+1,N
         READ (13,'(7e15.7)') X(i),Y(i),z(i),am(i),hp(i),rho(i),a(i)
         amtot2=amtot2+am(i)
      ENDDO
      close(13)

      if(myrank.eq.0) write(6,*)'qdar=',qdar,' from sph.input'
      if(myrank.eq.0) write(6,*)
     $     'm1=',amtot1,' m2=',amtot2,'m2/m1=',amtot2/amtot1
      
      if(myrank.eq.0) write(6,*)'NLEFT=',nleft,'  N=',n
      call gravquant
      
C Shift stars along x-axis so that CM is at the origin
C  and separation = SEP0
      AL=qdar*sep0/(1.0+qdar)
      if(myrank.eq.0) write(6,*)'SEP0=',sep0,'  AL=',al
      DO I=1,Nleft                                                    
         X(I)=X(I)-AL
      enddo
      AL=sep0/(2.0*(1.0+qdar))
      do i=nleft+1,n
         X(i)=x(i)+AL
      ENDDO                                                           

C Set all velocities equal to zero
      DO IP=1,N
         VX(IP)=0.                                                      
         VY(IP)=0.                                                      
         VZ(IP)=0.                                                      
      ENDDO                                                           

      if(myrank.eq.0)write(6,*)'2cr:nrelax=',nrelax,' trelax=',trelax
C prepare leap-frog scheme for first iteration                       
      CALL LFSTART                                                    
                                                                      
      RETURN                                                          
C     Error condition:
 100  STOP '2cr: ERROR OPENING INPUT FILE ???'
      END              

