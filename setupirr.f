**************************************************************
c     This file is part of the StarCrash code
c     Version 1.0
c
c     Copyright 2003, by Joshua Faber
c     This code is publicly available under the terms of the
c     GNU General Public License.  See our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

      SUBROUTINE setup2i1                                                
      IMPLICIT NONE
********************************************************
c     Release 1.0
C     Sets up a irrotational binary with equal-mass NS,
c     using previously calculated relaxed single-star models,
c     processed by b2pos into a file called "pos.sph".
c     We also need the file "axes1.sph", which specifies the
c     axis ratios of the stars.
c     Called by INIT
c     Calls GRAVQUANT,LFSTART
**********************************************************
      INCLUDE 'spha.h'

      INTEGER nold,nnoptold,i,ip,navold,i2,ip2
      REAL hminold,hmaxold,gamold,trelaxold,ar,a1,a2,a3,x0(nmax)
      REAL afact,amxt,amyt,amzt,amtot

C Get parameters from input file:
      OPEN(12,FILE='sph.input',ERR=100)
      READ(12,INPUT)
      CLOSE(12)

      open(13,file='pos.sph', form='formatted')
      read(13,'(2i6,4e15.7)')
     $     nold,nnoptold,hminold,hmaxold,gamold,trelaxold
      if(myrank.eq.0)write(6,*)
     $     'n_old=',nold,'nnopt_old=',nnoptold,'hmin_old=',hminold,
     $     'hmax_old=',hmaxold,'gam_old=',gamold,'trelax_old=',trelaxold


      if(2*nold.gt.nmax) then
         if(myrank.eq.0)write(6,*)'ERROR: 2N>NMAX'
         stop
      endif
      if(nnoptold.ne.nnopt)write(6,*)'nnopt changed from ',
     $     nnoptold,' to ',nnopt
      if(hminold.ne.hmin)write(6,*)'hmin changed from ',
     $     hminold,' to ',hmin
      if(hmaxold.ne.hmax)write(6,*)'hmax changed from ',
     $     hmaxold,' to ',hmax
      if(gamold.ne.gam)write(6,*)'gam changed from ',
     $     gamold,' to ',gam
      if(trelaxold.ne.trelax)write(6,*)'trelax changed from ',
     $     trelaxold,' to ',trelax

      amtot=0.
      amxt=0.
      amyt=0.
      amzt=0.
      DO I=1,Nold
         READ (13,'(7e15.7)') X(i),Y(i),z(i),am(i),hp(i),rho(i),a(i)
         amtot=amtot+am(i)
         amxt=amxt+am(i)*x(i)
         amyt=amyt+am(i)*y(i)
         amzt=amzt+am(i)*z(i)
      ENDDO
      close(13)
      n=nold

c     adjust CoM of initial configuration
      amxt=amxt/amtot
      amyt=amyt/amtot
      amzt=amzt/amtot
      write(6,*)'adjusting center of mass of pos.sph from'
      write(6,*)amxt,amyt,amzt,' to origin'
      do i=1,n
         x(i)=x(i)-amxt
         y(i)=y(i)-amyt
         z(i)=z(i)-amzt
      enddo

C There will be 2N particles (N for star on left) 
      NLEFT=N
                                                                      
      do i=1,n
         i2=i+n
         x(i2)=x(i)
         y(i2)=y(i)
         z(i2)=z(i)
         a(i2)=a(i)
         am(i2)=am(i)
         hp(i2)=hp(i)
         rho(i2)=rho(i)
      enddo

      nleft=n
      n=2*n
      if(myrank.eq.0) write(6,*)'NLEFT=',nleft,'  N=',n
      call gravquant
      
c     Deform stars into new axis ratios
      
      open(13,file='axes1.sph',form='formatted')
      read(13,'(3e15.7)')a1,a2,a3
      close(13)

      afact=2.0/(a1**2+a2**2)

      do ip=1,n
         x(ip)=x(ip)*a1
         y(ip)=y(ip)*a2
         z(ip)=z(ip)*a3
      enddo

C Shift stars along x-axis so that CM is at the origin
C  and separation = SEP0
      AR=sep0/2.0
      if(myrank.eq.0) write(6,*)'SEP0=',sep0,'  AR=',ar
      DO IP=1,Nleft                                                    
         ip2=ip+Nleft
         x0(ip)=x(ip)
         x0(ip2)=x(ip)
         X(IP)=X(IP)+AR
         X(ip2)=x(ip2)-AR
      ENDDO                                                           

c     adjust smoothing lengths toward optimal configuration
      do i=1,10
         call nene
         call adjust
      enddo
      call rhos

c     set velocity profile to irrotational condition
c     See Lombardi, Rasio, and Shapiro for details

c     Initialize velocities
      DO IP=1,N
         VX(IP)=0.                                                      
         VY(IP)=0.                                                      
         VZ(IP)=0.                                                      
      ENDDO                                                           

c     we need to calculate omega, finding inward accelerations by
c     pretending temporarily to use corotating relaxation techniques
      navold=nav
      nav=0
      nrelax=2
      trelax=1000.0
      call vdots
      
      nrelax=0
      nav=navold

      if(myrank.eq.0)write(6,*)'omega=',sqrt(omega2)

      DO IP=1,N
         VX(IP)=sqrt(omega2)*y(ip)*(afact*a1**2-1.0)
         VY(IP)=sqrt(omega2)*(x(ip)-afact*a2**2*x0(ip))
      ENDDO                                                           

      omega2=0

C prepare leap-frog scheme for first iteration                       
      CALL LFSTART                                                    
                                                                      
      RETURN                                                          
C     Error condition:
 100  STOP '2i1: ERROR OPENING INPUT FILE ???'
      END              

***********************************************************************

********************************************************************
      SUBROUTINE setup2iq                                                
      IMPLICIT NONE
********************************************************
c     Release 1.0
C     Sets up a irrotational binary with unequal-mass NS,
c     using previously calculated relaxed single-star models,
c     processed by b2pos and b2pos2 into files called "pos.sph"
c     and "pos2.sph".
c     We also need the files "axes1.sph" and 'axes2.sph", 
c     which specifies the axis ratios of the stars.
c     Called by INIT
c     Calls GRAVQUANT,LFSTART
**********************************************************
      INCLUDE 'spha.h'

      INTEGER nold,nnoptold,i,ip,navold
      REAL hminold,hmaxold,gamold,trelaxold,ar,a1,a2,a3,x0(nmax)
      REAL b1,b2,b3,bfact
      REAL afact,amxt,amyt,amzt,amtot

C Get parameters from input file:
      OPEN(12,FILE='sph.input',ERR=100)
      READ(12,INPUT)
      CLOSE(12)

      open(13,file='pos.sph', form='formatted')
      read(13,'(2i6,4e15.7)')
     $     nold,nnoptold,hminold,hmaxold,gamold,trelaxold
      if(myrank.eq.0)write(6,*)
     $     'n_old=',nold,'nnopt_old=',nnoptold,'hmin_old=',hminold,
     $     'hmax_old=',hmaxold,'gam_old=',gamold,'trelax_old=',trelaxold

      if(nnoptold.ne.nnopt)write(6,*)'nnopt changed from ',
     $     nnoptold,' to ',nnopt
      if(hminold.ne.hmin)write(6,*)'hmin changed from ',
     $     hminold,' to ',hmin
      if(hmaxold.ne.hmax)write(6,*)'hmax changed from ',
     $     hmaxold,' to ',hmax
      if(gamold.ne.gam)write(6,*)'gam changed from ',
     $     gamold,' to ',gam
      if(trelaxold.ne.trelax)write(6,*)'trelax changed from ',
     $     trelaxold,' to ',trelax

      amtot=0.
      amxt=0.
      amyt=0.
      amzt=0.
      DO I=1,Nold
         READ (13,'(7e15.7)') X(i),Y(i),z(i),am(i),hp(i),rho(i),a(i)
         amtot=amtot+am(i)
         amxt=amxt+am(i)*x(i)
         amyt=amyt+am(i)*y(i)
         amzt=amzt+am(i)*z(i)
      ENDDO
      close(13)
      nleft=nold

c     adjust CoM of initial configuration
      amxt=amxt/amtot
      amyt=amyt/amtot
      amzt=amzt/amtot
      write(6,*)'adjusting center of mass of pos.sph from'
      write(6,*)amxt,amyt,amzt,' to origin'
      do i=1,n
         x(i)=x(i)-amxt
         y(i)=y(i)-amyt
         z(i)=z(i)-amzt
      enddo

C There will be 2N particles (N for star on left) 
      NLEFT=N

      open(13,file='pos2.sph', form='formatted')
      read(13,'(2i6,4e15.7)')
     $     nold,nnoptold,hminold,hmaxold,gamold,trelaxold
      if(myrank.eq.0)write(6,*)
     $     'n_old=',nold,'nnopt_old=',nnoptold,'hmin_old=',hminold,
     $     'hmax_old=',hmaxold,'gam_old=',gamold,'trelax_old=',trelaxold


      if(nleft+nold.gt.nmax) then
         if(myrank.eq.0)write(6,*)'ERROR: 2N>NMAX'
         stop
      endif
      if(nnoptold.ne.nnopt)write(6,*)'nnopt changed from ',
     $     nnoptold,' to ',nnopt
      if(hminold.ne.hmin)write(6,*)'hmin changed from ',
     $     hminold,' to ',hmin
      if(hmaxold.ne.hmax)write(6,*)'hmax changed from ',
     $     hmaxold,' to ',hmax
      if(gamold.ne.gam)write(6,*)'gam changed from ',
     $     gamold,' to ',gam
      if(trelaxold.ne.trelax)write(6,*)'trelax changed from ',
     $     trelaxold,' to ',trelax

      amtot=0.
      amxt=0.
      amyt=0.
      amzt=0.
      DO I=nleft+1,nleft+Nold
         READ (13,'(7e15.7)') X(i),Y(i),z(i),am(i),hp(i),rho(i),a(i)
         amtot=amtot+am(i)
         amxt=amxt+am(i)*x(i)
         amyt=amyt+am(i)*y(i)
         amzt=amzt+am(i)*z(i)
      ENDDO
      close(13)

c     adjust CoM of initial configuration
      amxt=amxt/amtot
      amyt=amyt/amtot
      amzt=amzt/amtot
      write(6,*)'adjusting center of mass of pos2.sph from'
      write(6,*)amxt,amyt,amzt,' to origin'
      do i=nleft+1,nleft+nold         
         x(i)=x(i)-amxt
         y(i)=y(i)-amyt
         z(i)=z(i)-amzt
      enddo

C There will be 2N particles (N for star on left) 
      N=nleft+nold
                                                                      
      if(myrank.eq.0) write(6,*)'NLEFT=',nleft,'  N=',n
      call gravquant
      
c     Deform stars into new axis ratios
      
      open(13,file='axes1.sph',form='formatted')
      read(13,'(3e15.7)')a1,a2,a3
      close(13)

      afact=2.0/(a1**2+a2**2)

      do ip=1,nleft
         x(ip)=x(ip)*a1
         y(ip)=y(ip)*a2
         z(ip)=z(ip)*a3
      enddo

      open(13,file='axes2.sph',form='formatted')
      read(13,'(3e15.7)')b1,b2,b3
      close(13)

      bfact=2.0/(b1**2+b2**2)

      do ip=nleft+1,n
         x(ip)=x(ip)*b1
         y(ip)=y(ip)*b2
         z(ip)=z(ip)*b3
      enddo

C Shift stars along x-axis so that CM is at the origin
C  and separation = SEP0


      AR=qdar*sep0/(1.0+qdar)
      if(myrank.eq.0) write(6,*)'SEP0=',sep0,'  AR=',ar
      DO IP=1,Nleft                                                    
         x0(ip)=x(ip)
         X(IP)=X(IP)-AR
      enddo
      ar=sep0/(1.0+qdar)
      do ip=nleft+1,n
         x0(ip)=x(ip)
         X(ip)=x(ip)+AR
      ENDDO                                                           

c     adjust smoothing lengths toward optimal configuration
      do i=1,10
         call nene
         call adjust
      enddo
      call rhos

c     set velocity profile to irrotational condition
c     See Lombardi, Rasio, and Shapiro for details

c     Initialize velocities
      DO IP=1,N
         VX(IP)=0.                                                      
         VY(IP)=0.                                                      
         VZ(IP)=0.                                                      
      ENDDO                                                           

c     we need to calculate omega, finding inward accelerations by
c     pretending temporarily to use corotating relaxation techniques
      navold=nav
      nav=0
      nrelax=2
      trelax=1000.0
      call vdots
      
      nrelax=0
      nav=navold

      if(myrank.eq.0)write(6,*)'omega=',sqrt(omega2)

      DO IP=1,Nleft
         VX(IP)=sqrt(omega2)*y(ip)*(afact*a1**2-1.0)
         VY(IP)=sqrt(omega2)*(x(ip)-afact*a2**2*x0(ip))
      ENDDO                                                           
      do ip=nleft+1,n
         VX(IP)=sqrt(omega2)*y(ip)*(bfact*b1**2-1.0)
         VY(IP)=sqrt(omega2)*(x(ip)-bfact*b2**2*x0(ip))
      ENDDO                                                           

      omega2=0

C prepare leap-frog scheme for first iteration                       
      CALL LFSTART                                                    
                                                                      
      RETURN                                                          
C     Error condition:
 100  STOP '2iq: ERROR OPENING INPUT FILE ???'
      END              

