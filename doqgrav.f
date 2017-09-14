**************************************************************
c     This file is part of the StarCrash code
c     Version 1.0
c
c     Copyright 2003, by Joshua Faber
c     This code is publicly available under the terms of the
c     GNU General Public License.  See our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

      SUBROUTINE DOQGRAV
      IMPLICIT NONE
******************************************************
c     Release 1.0
c     Calculate accelerations including radiation reaction
c     according to BDS formalism
c     called by BALVDOTS,CLAVDOTS,NEWVDOTS
c     calls QGRAV2,QGRAVRAD
********************************************************
      INCLUDE 'spha.h'
      INCLUDE 'mpif.h'

      REAL q2xxold,q2xyold,q2xzold,q2yyold,q2yzold,q2zzold
      REAL pxx,pxy,pxz,pyy,pyz,pzz,pyx,pzx,pzy
      REAL frxl,fryl,frzl,frxr,fryr,frzr
      real ajdot,dxu5(n),dyu5(n),dzu5(n)
      INTEGER i,ip

c     WE need to wait an iteration after starting the dynamical run to 
c     calculate the quadrupole derivatives accurately.  INITGR counts up
c     to 3 over the course of two iterations, and then starts including
c     radiation reaction damping
      if(initgr.eq.0) then
c     initialize 2nd quadrupole derivatives
         initgr=1
         q2xx=0.0
         q2xy=0.0
         q2xz=0.0
         q2yy=0.0
         q2yz=0.0
         q2zz=0.0
      endif
      q2xxold=q2xx
      q2xyold=q2xy
      q2xzold=q2xz
      q2yyold=q2yy
      q2yzold=q2yz
      q2zzold=q2zz

c     QGRAV2 does everything that QGRAV does, as well as calculating 2nd
c     spatial derivatives of the gravitational potential
      call qgrav2
      
c     pxx et al. are the quadrupole moments, 
c     qxx et al. are the traceless quadrupole moments
      pxx=0.0
      pxy=0.0
      pxz=0.0
      pyx=0.0
      pyy=0.0
      pyz=0.0
      pzx=0.0
      pzy=0.0
      pzz=0.0
      do ip=1,n
         pxx=pxx+2*(am(ip)*(vx(ip)*vx(ip)+x(ip)*gx(ip)))
         pxy=pxy+2*(am(ip)*(vx(ip)*vy(ip)+x(ip)*gy(ip)))
         pxz=pxz+2*(am(ip)*(vx(ip)*vz(ip)+x(ip)*gz(ip)))
         pyx=pyx+2*(am(ip)*(vy(ip)*vx(ip)+y(ip)*gx(ip)))
         pyy=pyy+2*(am(ip)*(vy(ip)*vy(ip)+y(ip)*gy(ip)))
         pyz=pyz+2*(am(ip)*(vy(ip)*vz(ip)+y(ip)*gz(ip)))
         pzx=pzx+2*(am(ip)*(vz(ip)*vx(ip)+z(ip)*gx(ip)))
         pzy=pzy+2*(am(ip)*(vz(ip)*vy(ip)+z(ip)*gy(ip)))
         pzz=pzz+2*(am(ip)*(vz(ip)*vz(ip)+z(ip)*gz(ip)))
      enddo
      q2xx=(2*pxx-pyy-pzz)/3.0
      q2xy=(pxy+pyx)/2.0
      q2xz=(pxz+pzx)/2.0
      q2yy=(2*pyy-pxx-pzz)/3.0
      q2yz=(pyz+pzy)/2.0
      q2zz=(2*pzz-pxx-pyy)/3.0
      q3xx=(q2xx-q2xxold)/dt
      q3xy=(q2xy-q2xyold)/dt
      q3xz=(q2xz-q2xzold)/dt
      q3yy=(q2yy-q2yyold)/dt
      q3yz=(q2yz-q2yzold)/dt
      q3zz=(q2zz-q2zzold)/dt
      if(myrank.eq.0) write(6,*)'Q2:',q2xx,q2xy,q2xz,q2yy,q2yz,q2zz
      if(myrank.eq.0) write(6,*)'Q3:',q3xx,q3xy,q3xz,q3yy,q3yz,q3zz
      

c     QGRAVRAD solves for the field R, defined by Blanchet, Damour, 
c     and Schaefer, Eq. (5.14)
      CALL QGRAVRAD

c     Now we calculate BDS, Eqs. 5.15, 5.19, and 5.20c for the radiative vdot
      do ip=1,n
         dxu5(ip)=q3xx*gx(ip)+q3xy*gy(ip)+q3xz*gz(ip)+
     $        gxx(ip)*(q3xx*x(ip)+q3xy*y(ip)+q3xz*z(ip))+
     $        gxy(ip)*(q3xy*x(ip)+q3yy*y(ip)+q3yz*z(ip))+
     $        gxz(ip)*(q3xz*x(ip)+q3yz*y(ip)+q3zz*z(ip))
         dyu5(ip)=q3xy*gx(ip)+q3yy*gy(ip)+q3yz*gz(ip)+
     $        gxy(ip)*(q3xx*x(ip)+q3xy*y(ip)+q3xz*z(ip))+
     $        gyy(ip)*(q3xy*x(ip)+q3yy*y(ip)+q3yz*z(ip))+
     $        gyz(ip)*(q3xz*x(ip)+q3yz*y(ip)+q3zz*z(ip))
         dzu5(ip)=q3xz*gx(ip)+q3yz*gy(ip)+q3zz*gz(ip)+
     $        gxz(ip)*(q3xx*x(ip)+q3xy*y(ip)+q3xz*z(ip))+
     $        gyz(ip)*(q3xy*x(ip)+q3yy*y(ip)+q3yz*z(ip))+
     $        gzz(ip)*(q3xz*x(ip)+q3yz*y(ip)+q3zz*z(ip))
         
         if(initgr.ge.3) then
            freacx(ip)=.4*(dxpnr(ip)-dxu5(ip))/sol**5
            freacy(ip)=.4*(dypnr(ip)-dyu5(ip))/sol**5
            freacz(ip)=.4*(dzpnr(ip)-dzu5(ip))/sol**5
         else
            freacx(ip)=0.
            freacy(ip)=0.
            freacz(ip)=0.
         endif
      enddo

c     track the resulting force on each star
      frxl=0.
      fryl=0.
      frzl=0.
      frxr=0.
      fryr=0.
      frzr=0.
      ajdot=0.
      do i=1,nleft
         frxl=frxl+am(i)*freacx(i)
         fryl=fryl+am(i)*freacy(i)
         frzl=frzl+am(i)*freacz(i)
         ajdot=ajdot+am(i)*(x(i)*freacy(i)-y(i)*freacx(i))
      enddo
      do i=nleft+1,n
         frxr=frxr+am(i)*freacx(i)
         fryr=fryr+am(i)*freacy(i)
         frzr=frzr+am(i)*freacz(i)
         ajdot=ajdot+am(i)*(x(i)*freacy(i)-y(i)*freacx(i))
      enddo

      if(myrank.eq.0) write(6,*)'RAD. REAC. Force, left:',
     $     frxl,fryl,frzl
      if(myrank.eq.0) write(6,*)'RAD. REAC. Force, right:',
     $     frxr,fryr,frzr
      if(myrank.eq.0) write(6,*)'ANG. MOM LOSS:',ajdot

      if(initgr.lt.3)initgr=initgr+1
      
c     if(myrank.eq.0)write(6,*) 'doqgrav: done'
      return
      end















