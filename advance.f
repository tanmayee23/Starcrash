**************************************************************
c     This file is part of the StarCrash code
c     Version 1.0
c
c     Copyright 2003, by Joshua Faber
c     This code is publicly available under the terms of the
c     GNU General Public License.  See our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

      SUBROUTINE TSTEP
      IMPLICIT NONE
***********************************************************
c     Release 1.0
c     Calculate the timestep
c     Called by INIT, MAINIT
***************************************************************
      INCLUDE 'spha.h'                                          
      INCLUDE 'mpif.h'

      REAL UIJMAX(NMAX),DIVV(NMAX),TINY
      COMMON/UIJMAX/ UIJMAX
      COMMON/COMMDIVV/DIVV
      PARAMETER (TINY=1.E-10)

      INTEGER i
      REAL dtmin,ci2,dtivel,ai,dtiacc,dti

      DTIVEL=0.
      DTMIN=1500000.
      DO I=1,N
         IF (NN(I).NE.0) THEN
c     calculate sound-crossing time over a smoothing length, including AV
            CI2=GAM*A(I)*RHO(I)**(GAM-1.)
            IF(NAV.EQ.0) THEN
               DTIVEL=HP(I)/SQRT(CI2)
            ELSE IF (NAV.EQ.2) THEN
               DTIVEL=HP(I)/(SQRT(CI2)*(1.+1.2*ALPHA)+1.2*BETA*
     $              HP(I)*MAX(0.,-DIVV(I)))
            ELSE IF (NAV.EQ.1 .OR. NAV.GE.3) THEN
               DTIVEL=HP(I)/(SQRT(CI2)*(1.+1.2*ALPHA)+1.2*
     $              BETA*UIJMAX(I))
            ENDIF
            AI=SQRT(VXDOT(I)**2+VYDOT(I)**2+VZDOT(I)**2)
c     calculate timescale on which velocity changes
            DTIACC=SQRT(HP(I)/(AI+TINY))
            DTI=MIN(DTIVEL,DTIACC)
            DTMIN=MIN(DTMIN,DTI)
         ENDIF
      ENDDO
      DT=CN*DTMIN
      IF (DT+T .GT. TF) THEN
         DT=TF-T
      ENDIF
      RETURN
      END
************************************************************************
      SUBROUTINE ADVANCE
      IMPLICIT NONE
*****************************************************************
c     Release 1.0
c     Advance hydro quantities by a timestep, using a leapfrog algorithm,
c     i.e. time derivatives are calculated at mid-timestep for second
c     order accuracy
c     Called by MAINIT
c     Calls NENE,RHOS,ADOTS,CMADJ,VDOTS
***************************************************************
      INCLUDE 'spha.h'                                          
      REAL UTINY
      PARAMETER (UTINY=1.E-10)

      REAL AO(NMAX),VXO(NMAX),VYO(NMAX),VZO(NMAX)
      REAL DTH
      INTEGER I

      DTH=0.5*DT

C     Remember old values:
      IF (NAV.GE.1.and.nrelax.eq.0) THEN
         DO I=1,N
            AO(I)=A(I)
         ENDDO  
      ENDIF
C     Advance entropies and positions to half-timestep:
      IF (NAV.GE.1.and.nrelax.eq.0) THEN
         DO I=1,N
            A(I)=A(I)+DTH*ADOT(I)
         ENDDO
      ENDIF
      DO I=1,N   
c     If radiation reaction is on
         if(initgr.eq.3) then
            X(I)=X(I)+DTH*UX(I)
            Y(I)=Y(I)+DTH*UY(I)
            Z(I)=Z(I)+DTH*UZ(I)
         else 
            X(I)=X(I)+DTH*VX(I)
            Y(I)=Y(I)+DTH*VY(I)
            Z(I)=Z(I)+DTH*VZ(I)
         endif
      ENDDO

      IF (NAV.GE.1.and.nrelax.eq.0) THEN      
C     Calculate interaction lists and densities at half-timestep:
         CALL NENE
         CALL RHOS
C     Calculate Adot (only if AV is on, and relaxation is off):
         CALL ADOTS
C     Advance entropies from original values:
         DO I=1,N
            A(I)=AO(I)+DT*ADOT(I)
         ENDDO
      ENDIF
      
C     Advance positions to next timestep:
      DO I=1,N
         if(initgr.eq.3) then
c     ux is defined below, it is used only when Radiation reaction is on
            X(I)=X(I)+DTH*UX(I)
            Y(I)=Y(I)+DTH*UY(I)
            Z(I)=Z(I)+DTH*UZ(I)
         else
            X(I)=X(I)+DTH*VX(I)
            Y(I)=Y(I)+DTH*VY(I)
            Z(I)=Z(I)+DTH*VZ(I)
         endif
      ENDDO
      
C     If binary relaxation problem, adjust center of mass positions: 
      IF (NRELAX.GE.1) CALL CMADJ                                    
      
C     Recalculate interaction lists and densities:
      CALL NENE
      CALL RHOS

      DO I=1,N  
         VXO(I)=VX(I)
         VYO(I)=VY(I)
         VZO(I)=VZ(I)
      ENDDO
      
C     Advance velocities to half-timestep:
      DO I=1,N
         VX(I)=VX(I)+DTH*VXDOT(I)
         VY(I)=VY(I)+DTH*VYDOT(I)
         VZ(I)=VZ(I)+DTH*VZDOT(I)
      ENDDO
      
C     Calculate accelerations (pressure gradients + gravity):
      CALL VDOTS
     
C     Advance velocities:
      DO I=1,N
         VX(I)=VXO(I)+DT*VXDOT(I)
         VY(I)=VYO(I)+DT*VYDOT(I)
         VZ(I)=VZO(I)+DT*VZDOT(I)
      ENDDO

      do i=1,n
         if(initgr.eq.3) then
c     this is taken from Blanchet, Damour, and Schaefer
c     MNRAS 242, 289 (1990), Eq. (5.16),  with beta=A_i=0 in the Newtonian 
c     limit. Note that our code's "v" is written there as "w", 
c     and our code's "u" is written as "v"
            ux(i)=vx(i)+0.8/sol**5*
     $           (q3xx*vx(i)+q3xy*vy(i)+q3xz*vz(i))
            uy(i)=vy(i)+0.8/sol**5*
     $           (q3xy*vx(i)+q3yy*vy(i)+q3yz*vz(i))
            uz(i)=vz(i)+0.8/sol**5*
     $           (q3xz*vx(i)+q3yz*vy(i)+q3zz*vz(i))
         else
            ux(i)=vx(i)
            uy(i)=vy(i)
            uz(i)=vz(i)
         endif
      enddo

      T=T+DT

      RETURN
      END
************************************************************************
      SUBROUTINE CMADJ                                                 
      IMPLICIT NONE
******************************************************************
c     Release 1.0
c     Adjust center-of-mass for binary relaxation
c     Called by ADVANCE
*******************************************************************
      INCLUDE 'spha.h'                                          

      REAL xcm1,ycm1,zcm1,am1,xcm2,ycm2,zcm2,am2
      REAL delx1,dely1,delz1,delx2,dely2,delz2
      INTEGER i

      XCM1=0.
      YCM1=0.
      ZCM1=0.
      AM1=0.                                                           
      xcm2=0.
      ycm2=0.
      zcm2=0.
      am2=0.
      DO I=1,NLEFT   
         AM1=AM1+AM(I)
         XCM1=XCM1+AM(I)*X(I)                                           
         YCM1=YCM1+AM(I)*Y(I)                                           
         ZCM1=ZCM1+AM(I)*Z(I)                                           
      ENDDO                                                            
      XCM1=XCM1/AM1                                                    
      YCM1=YCM1/AM1                                                    
      ZCM1=ZCM1/AM1  

      if(nrelax.eq.2) then
         DO I=NLEFT+1,N   
            AM2=AM2+AM(I)
            XCM2=XCM2+AM(I)*X(I)                                           
            YCM2=YCM2+AM(I)*Y(I)                                           
            ZCM2=ZCM2+AM(I)*Z(I)                                           
         ENDDO                                                            
         XCM2=XCM2/AM2                                                    
         YCM2=YCM2/AM2                                                    
         ZCM2=ZCM2/AM2 
      endif

      if(nrelax.eq.1) then
         DELX1=-XCM1
      else if(nrelax.eq.2) then
         DELX1=-XCM1-AM2*SEP0/(AM1+AM2)                                   
      endif
      DELY1=-YCM1
      DELZ1=-ZCM1
      if(nrelax.eq.2) then
         DELX2=-XCM2+AM1*SEP0/(AM1+AM2)                                   
         DELY2=-YCM2
         DELZ2=-ZCM2
      endif

      if(myrank.eq.0)write(6,*)'CMADJ: Star 1:',delx1,dely1,delz1
      if(nrelax.eq.2.and.myrank.eq.0)
     $     write(6,*)'CMADJ: Star 2:',delx2,dely2,delz2
                                                      
      DO I=1,NLEFT
        X(I)=X(I)+DELX1
        Y(I)=Y(I)+DELY1
        Z(I)=Z(I)+DELZ1
      ENDDO                                                            
      if(nrelax.eq.2) then
         DO I=NLEFT+1,N
            X(I)=X(I)+DELX2
            Y(I)=Y(I)+DELY2
            Z(I)=Z(I)+DELZ2                                                
         ENDDO                                                            
      endif                                  
                                     
      RETURN                                                           
      END
************************************************************************
      SUBROUTINE RELAX                                                 
      IMPLICIT NONE
***************************************************************
c     Release 1.0
c     Calculate drag forces for runs with relaxation, as well as
c     centrifugal acceleration for synchronized binaries.
c     Called by ADVANCE
************************************************************
      INCLUDE 'spha.h'                                          

      REAL xcm1,am1,gcm1,fcm1,xcm2,am2,gcm2,fcm2
      INTEGER i

C     Note: TRELAX and other parameters must be set-up previously.
C     NRELAX=0 not a relaxation calculation                    
C     NRELAX=1 for simple drag force
C     NRELAX=2 for drag + centrifugal at constant binary
C     separation (set-up binary configuration)
      
      IF (NRELAX.EQ.1) THEN
         DO I=1,N
            VXDOT(I)=VXDOT(I)-VX(I)/TRELAX
            VYDOT(I)=VYDOT(I)-VY(I)/TRELAX
            VZDOT(I)=VZDOT(I)-VZ(I)/TRELAX
         ENDDO
      ELSE IF (NRELAX.EQ.2) THEN             
C     Calculate Omega so that the total net force (pressure+
C     gravity+centrifugal) on the centers of mass is zero)
         XCM1=0.     
         AM1=0.                                                         
         GCM1=0.                                                        
         FCM1=0.
         DO I=1,NLEFT 
            AM1=AM1+AM(I)                                                
            XCM1=XCM1+AM(I)*X(I)                                         
            GCM1=GCM1+AM(I)*GX(I) 
            FCM1=FCM1+AM(I)*VXDOT(I)
         ENDDO                                                          
         XCM1=XCM1/AM1                                                  
         GCM1=GCM1/AM1                                                  
         FCM1=FCM1/AM1
         XCM2=0.                                                        
         AM2=0.                                                         
         GCM2=0.                                                        
         FCM2=0.
         DO I=NLEFT+1,N                                                   
            AM2=AM2+AM(I)                                                
            XCM2=XCM2+AM(I)*X(I)                                         
            GCM2=GCM2+AM(I)*GX(I)
            FCM2=FCM2+AM(I)*VXDOT(I)
         ENDDO                                                          
         XCM2=XCM2/AM2                                                  
         GCM2=GCM2/AM2             
         FCM2=FCM2/AM2
         if(myrank.eq.0)WRITE (6,*) 
     $        'CDRAG: GCM1/XCM1=',ABS(GCM1/XCM1),           
     $        ' GCM2/XCM2=',ABS(GCM2/XCM2)             
         if(myrank.eq.0)WRITE (6,*) 
     $        'CDRAG: FCM1/XCM1=',ABS(FCM1/XCM1),           
     $        ' FCM2/XCM2=',ABS(FCM2/XCM2)             
C     (Note that the two terms should be equal by Newton's 3rd law)
         OMEGA2=0.5*(ABS(FCM1/XCM1)+ABS(FCM2/XCM2))                     
      ELSE                                                             
         STOP 'RELAX: NRELAX UNKNOWN ???'                               
      ENDIF                                                            

      IF (NRELAX.EQ.2) THEN
c     add in centrifugal acceleration, since VX is measured in 
c     corotating frame
         Do i=1,N
            VXDOT(I)=VXDOT(I)-VX(I)/TRELAX+OMEGA2*X(I)                     
            VYDOT(I)=VYDOT(I)-VY(I)/TRELAX+OMEGA2*Y(I)                     
            VZDOT(I)=VZDOT(I)-VZ(I)/TRELAX
         ENDDO
      ENDIF
      
      RETURN                                                           
      END
************************************************************************
      SUBROUTINE RHOS
************************************************************
c     Release 1.0
c     Calculate SPH density at each particle position
c     Called by ADVANCE
***************************************************************

      IMPLICIT NONE

C Compute densities and P/rho^2:
      INCLUDE 'spha.h'                                         
      INCLUDE 'mpif.h'

      REAL rhotiny
      PARAMETER (RHOTINY=1.E-10)
      real myrho(n)

      REAL ami,hpi,h2,h3,rhoi,r2,wijini
      INTEGER i,ii,in,j,itab,ierr

C Initialize:
      DO I=1,N
        MYRHO(I)=0.
        rho(i)=0.
      ENDDO

C  For each particle:
      DO I=n_lower,n_upper
        II=I-n_lower+1
        AMI=AM(I)
        HPI=HP(I)
        H2=HPI**2
        H3=HPI*H2
C     Accumulate contributions:
C     "Gather" part of the sum (i-j pair contributes to rho_i):
        RHOI=0.
        DO IN=1,NN(I)
          J=NNI(IN,II)
          R2=XIJ(IN,II)**2.+YIJ(IN,II)**2.+ZIJ(IN,II)**2.
          ITAB=INT(CTAB*R2/H2)+1
          WIJINI=WTAB(ITAB)/H3
          RHOI=RHOI+AM(J)*WIJINI
        ENDDO
        MYRHO(I)=MYRHO(I)+0.5*RHOI
C     "Scatter" part of the sum (i-j pair contributes to rho_j):
        DO IN=1,NN(I)
           J=NNI(IN,II)
           R2=XIJ(IN,II)**2.+YIJ(IN,II)**2.+ZIJ(IN,II)**2.
           ITAB=INT(CTAB*R2/H2)+1
           WIJINI=WTAB(ITAB)/H3
           MYRHO(J)=MYRHO(J)+0.5*AMI*WIJINI
        ENDDO
      ENDDO

c     We have the density subtotals, calculated from i=n_lower to n_upper
c     Now, we sum them all, and distribute the result to all processes
      CALL MPI_ALLREDUCE(myrho,rho,n,MPI_REAL,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)

C Do something crude for particles with NN(I)=0 (to avoid
C having RHO(I)=0):
      DO I=1,n
         RHO(I)=MAX(RHO(I), RHOTINY )
      ENDDO

C Calculate p/rho**2 for all particles:
      DO I=1,N
         POR2(I)=A(I)*RHO(I)**(GAM-2.)
      ENDDO

      RETURN
      END
************************************************************************











