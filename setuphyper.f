**************************************************************
c     This file is part of the StarCrash code
c     Version 1.0
c
c     Copyright 2003, by Joshua Faber
c     This code is publicly available under the terms of the
c     GNU General Public License.  See our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

      subroutine setuphyp1
****************************************************8888
c     Release 1.0
C     Sets up a hyperbolic binary orbit with equal-mass NS,
c     using a previously calculated relaxed single-star model,
c     processed by b2pos into a file called "pos.sph".
c     Called by INIT
c     Calls GRAVQUANT,LFSTART
************************************************************
c      implicit none

      INCLUDE 'spha.h'

      REAL amtot,amu,rimp,vphi,etot,altot,vmin,rmin,ke,pe,v0,vt,vr
      REAL theta,vxip,vyip,vzip,vmag,vrip,rip,hminold,hmaxold
      REAL gamold,trelaxold
      INTEGER i,i2,nold,nnoptold

      OPEN(12,FILE='sph.input')
      READ(12,INPUT)
      CLOSE(12)

      open(13,file='pos.sph', form='formatted')
      read(13,'(2i6,4e15.7)')
     $     nold,nnoptold,hminold,hmaxold,gamold,trelaxold

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

      amtot=0.0
      DO I=1,Nold
         READ (13,'(7e15.7)') X(i),Y(i),z(i),am(i),hp(i),rho(i),a(i)
         amtot=amtot+am(i)
      ENDDO
      close(13)

C There are now 2N particles (N for star on left)
      NLEFT=Nold

      do i=1,nleft
         i2=i+nleft
         x(i2)=x(i)
         y(i2)=y(i)
         z(i2)=z(i)
         a(i2)=a(i)
         am(i2)=am(i)
         hp(i2)=hp(i)
         rho(i2)=rho(i)
      enddo

      amu=amtot/2.0
      amtot=2.0*amtot

      rimp=sqrt(rp**2+2*amtot/vpeak**2*rp)
      
      vphi=sqrt(amtot/rimp)
      etot=0.5*amu*vpeak**2
      altot=amu*vpeak*rimp
      
      vmin=vphi**2/vpeak+sqrt(vpeak**2+vphi**4/vpeak**2)
      rmin=vpeak*rimp/vmin
      
      ke=0.5*amu*vmin**2
      pe=amtot*amu/rmin
      
      v0=sqrt(vpeak**2+2*amtot/sep0)
      vt=vpeak*rimp/sep0
      vr=sqrt(v0**2-vt**2)
      theta=180.0/3.141592654*atan(vt/vr)
      
      if(myrank.eq.0)write(6,*)'Impact Parameter:',rimp
      if(myrank.eq.0)write(6,*)'Rmin and check:',rmin,rp
      if(myrank.eq.0)write(6,*)'v(rmin):',vmin,'v(sep0)',v0
      if(myrank.eq.0)write(6,*)'Initial velocity-transverse,radial'
      if(myrank.eq.0)write(6,*)'vt=:',vt,' vr=:',vr,' theta=',theta 

c     vt in direction (-1/sqrt(6),-1/sqrt(6),2/sqrt(6)) from neg corner
c     vr in direction (1/sqrt(3) ,1/sqrt(3) ,1/sqrt(3)) from neg corner

      vxip=vr/sqrt(3.0)-vt/sqrt(6.0)
      vyip=vxip
      vzip=vr/sqrt(3.0)+2.0*vt/sqrt(6.0)
      vmag=sqrt(vxip**2+vyip**2+vzip**2)
      vrip=(vxip+vyip+vzip)/sqrt(3.0)
      rip=sep0/sqrt(12.0)

      if(myrank.eq.0)write(6,*)
     $     'v0=',vxip,vyip,vzip,'vr=',vrip,'vmag=',vmag
      
      do i=1,nleft
         x(i+nleft)=x(i)+rip
         x(i)=x(i)-rip
         y(i+nleft)=y(i)+rip
         y(i)=y(i)-rip
         z(i+nleft)=z(i)+rip
         z(i)=z(i)-rip
         vx(i)=vxip/2.0
         vx(i+nleft)=-1.0*vxip/2.0
         vy(i)=vyip/2.0
         vy(i+nleft)=-1.0*vyip/2.0
         vz(i)=vzip/2.0
         vz(i+nleft)=-1.0*vzip/2.0
         am(i+nleft)=am(i)
         hp(i+nleft)=hp(i)
         rho(i+nleft)=rho(i)
         a(i+nleft)=a(i)
         vxdot(i)=0.
         vydot(i)=0.
         vxdot(i)=0.
         vzdot(i+nleft)=0.
         vydot(i+nleft)=0.
         vzdot(i+nleft)=0.
         adot(i)=0.
         adot(i+nleft)=0.
      enddo
      n=2*nleft
      t=0

      call gravquant
      call lfstart

      return
      end

*************************************************************************

      subroutine setuphypq
****************************************************8888
C     Sets up a hyperbolic binary orbit with unequal-mass NS,
c     using previously calculated relaxed single-star models,
c     processed by b2pos and b2pos2 into files called "pos.sph"
c     and "pos2.sph".
c     Called by INIT
c     Calls GRAVQUANT,LFSTART
************************************************************
      implicit none

      INCLUDE 'spha.h'

      REAL amtot,amu,rimp,vphi,etot,altot,vmin,rmin,ke,pe,vt,vr
      REAL theta,vxip,vyip,vzip,vmag,vrip,rip,hminold,hmaxold
      REAL gamold,trelaxold,amtot1,amtot2,v0
      INTEGER i,i2,nold,nnoptold

      OPEN(12,FILE='sph.input')
      READ(12,INPUT)
      CLOSE(12)

      open(13,file='pos.sph', form='formatted')
      read(13,'(2i6,4e15.7)')
     $     nold,nnoptold,hminold,hmaxold,gamold,trelaxold

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

      amtot1=0.0
      DO I=1,Nold
         READ (13,'(7e15.7)') X(i),Y(i),z(i),am(i),hp(i),rho(i),a(i)
         amtot1=amtot1+am(i)
      ENDDO
      close(13)

      nleft=nold

      open(13,file='pos2.sph', form='formatted')
      read(13,'(2i6,4e15.7)')
     $     nold,nnoptold,hminold,hmaxold,gamold,trelaxold

      if(nold+nleft.gt.nmax) then
         if(myrank.eq.0)write(6,*)'ERROR: N>NMAX'
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

      amtot2=0.0
      DO I2=1,Nold
         i=nleft+i2
         READ (13,'(7e15.7)') X(i),Y(i),z(i),am(i),hp(i),rho(i),a(i)
         amtot2=amtot2+am(i)
      ENDDO
      close(13)

      qdar=amtot2/amtot1
      if(myrank.eq.0)write(6,*)'Mass ratio q=',qdar

      amu=amtot1*amtot2/(amtot1+amtot2)
      amtot=amtot1+amtot2

      rimp=sqrt(rp**2+2*amtot/vpeak**2*rp)
      
      vphi=sqrt(amtot/rimp)
      etot=0.5*amu*vpeak**2
      altot=amu*vpeak*rimp
      
      vmin=vphi**2/vpeak+sqrt(vpeak**2+vphi**4/vpeak**2)
      rmin=vpeak*rimp/vmin
      
      ke=0.5*amu*vmin**2
      pe=amtot*amu/rmin
      
      v0=sqrt(vpeak**2+2*amtot/sep0)
      vt=vpeak*rimp/sep0
      vr=sqrt(v0**2-vt**2)
      theta=180.0/3.141592654*atan(vt/vr)
      
      if(myrank.eq.0)write(6,*)'Impact Parameter:',rimp
      if(myrank.eq.0)write(6,*)'Rmin and check:',rmin,rp
      if(myrank.eq.0)write(6,*)'v(rmin):',vmin,'v(sep0)',v0
      if(myrank.eq.0)write(6,*)'Initial velocity-transverse,radial'
      if(myrank.eq.0)write(6,*)'vt=:',vt,' vr=:',vr,' theta=',theta 

c     vt in direction (-1/sqrt(6),-1/sqrt(6),2/sqrt(6)) from neg corner
c     vr in direction (1/sqrt(3) ,1/sqrt(3) ,1/sqrt(3)) from neg corner

      vxip=vr/sqrt(3.0)-vt/sqrt(6.0)
      vyip=vxip
      vzip=vr/sqrt(3.0)+2.0*vt/sqrt(6.0)
      vmag=sqrt(vxip**2+vyip**2+vzip**2)
      vrip=(vxip+vyip+vzip)/sqrt(3.0)
      rip=sep0/sqrt(3.0)

      if(myrank.eq.0)write(6,*)
     $     'v0=',vxip,vyip,vzip,'vr=',vrip,'vmag=',vmag
      
      do i=1,nleft
         x(i)=x(i)-rip*qdar/(1.0+qdar)
         y(i)=y(i)-rip*qdar/(1.0+qdar)
         z(i)=z(i)-rip*qdar/(1.0+qdar)
         vx(i)=vxip*qdar/(1.0+qdar)
         vy(i)=vyip*qdar/(1.0+qdar)
         vz(i)=vzip*qdar/(1.0+qdar)
         vxdot(i)=0.
         vydot(i)=0.
         vxdot(i)=0.
         adot(i)=0.
      enddo
      do i=nleft+1,nleft+nold
         x(i)=x(i)+rip/(1.0+qdar)
         y(i)=y(i)+rip/(1.0+qdar)
         z(i)=z(i)+rip/(1.0+qdar)
         vx(i)=-1.0*vxip/(1.0+qdar)
         vy(i)=-1.0*vyip/(1.0+qdar)
         vz(i)=-1.0*vzip/(1.0+qdar)
         vxdot(i)=0.
         vydot(i)=0.
         vzdot(i)=0.
         adot(i)=0.
      enddo

      n=nleft+nold
      t=0
      
      call gravquant
      call lfstart

      return
      end

      
      
