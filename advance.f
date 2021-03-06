
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
      REAL dtmin,ci2,dtivel,ai,dtiacc,dti,dten

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
c            dten = en(i)/(abs(endot(i))+TINY)
cy            dtmin = min(dtmin,dten)
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
          WIJINI=WTAB(ITAB)*2*PI/3.0/HPI!WTAB(ITAB)/H3
          RHOI=RHOI+AM(J)*WIJINI
        ENDDO
        MYRHO(I)=MYRHO(I)+0.5*RHOI
C     "Scatter" part of the sum (i-j pair contributes to rho_j):
        DO IN=1,NN(I)
           J=NNI(IN,II)
           R2=XIJ(IN,II)**2.+YIJ(IN,II)**2.+ZIJ(IN,II)**2.
           ITAB=INT(CTAB*R2/H2)+1
           WIJINI=WTAB(ITAB)*2*PI/3.0/HPI!WTAB(ITAB)/H3
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
         if(rho(i).le.0)write(6,*)'rhos: Nonpositive rho:',i,rho(i)


         RHO(I)=MAX(RHO(I), RHOTINY )
      ENDDO

C Calculate p/rho**2 for all particles:
      DO I=1,N
         if(a(i).le.0)write(6,*)'rhos: Nonpositive a:',i,a(i)

         POR2(I)=A(I)*RHO(I)**(GAM-2.)
      ENDDO

      RETURN
      END
************************************************************************
      SUBROUTINE RHOS_1
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
          WIJINI=WTAB(ITAB)*2*PI/3.0/HPI !WTAB(ITAB)/H3
          RHOI=RHOI+AM(J)*WIJINI
        ENDDO
        MYRHO(I)=MYRHO(I)+0.5*RHOI
C     "Scatter" part of the sum (i-j pair contributes to rho_j):
        DO IN=1,NN(I)
           J=NNI(IN,II)
           R2=XIJ(IN,II)**2.+YIJ(IN,II)**2.+ZIJ(IN,II)**2.
           ITAB=INT(CTAB*R2/H2)+1
           WIJINI=WTAB(ITAB)*2*PI/3.0/HPI ! WTAB(ITAB)/H3
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
      SUBROUTINE ADVANCE_EULER
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
      CALL NENE
      CALL RHOS
      CALL VDOTS

      IF (NAV.GE.1.and.nrelax.eq.0) THEN
C     Calculate Adot (only if AV is on, and relaxation is off):
         CALL ADOTS
C     Advance entropies from original values:
         DO I=1,N
            A(I)=A(I)+DT*ADOT(I)
         ENDDO
      ENDIF
                                                                                  
C     Advance positions to full timestep:
      DO I=1,N
         if(initgr.eq.3) then
c     ux is defined below, it is used only when Radiation reaction is on
            X(I)=X(I)+DT*UX(I)
            Y(I)=Y(I)+DT*UY(I)
            Z(I)=Z(I)+DT*UZ(I)
         else
            X(I)=X(I)+DT*VX(I)
            Y(I)=Y(I)+DT*VY(I)
            Z(I)=Z(I)+DT*VZ(I)
         endif
      ENDDO
      
C     If binary relaxation problem, adjust center of mass positions: 
      IF (NRELAX.GE.1) CALL CMADJ                                    
      
      
C     Advance velocities to full-timestep:
      DO I=1,N
         VX(I)=VX(I)+DT*VXDOT(I)
         VY(I)=VY(I)+DT*VYDOT(I)
         VZ(I)=VZ(I)+DT*VZDOT(I)
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
***********************************************************************                                                                                      
      SUBROUTINE ADVANCE_RK2
      IMPLICIT NONE
*****************************************************************                                                                                             
c     Release 1.0                                                                                                   
c     Advance hydro quantities by a timestep, using second order Runge Kutta method
c     Called by MAINIT                                                                                                                                       
c     Calls NENE,RHOS,ADOTS,CMADJ,VDOTS                                                                                                                      
***************************************************************                                                                                               
      INCLUDE 'spha.h'
      REAL UTINY
      PARAMETER (UTINY=1.E-10)

      REAL AOLD(NMAX), K1_X(NMAX),K1_Y(NMAX),K1_Z(NMAX)
      REAL DTH,K2_X(NMAX),K2_Y(NMAX),K2_Z(NMAX)
      REAL K1_VX(NMAX),K1_VY(NMAX),K1_VZ(NMAX)
      REAL K2_VX(NMAX),K2_VY(NMAX),K2_VZ(NMAX)
      REAL Xold(NMAX),Yold(NMAX),Zold(NMAX),A2(NMAX)
      REAL VXold(NMAX),VYold(NMAX),VZold(NMAX),A1(NMAX)
      INTEGER I

      DTH=0.5*DT
      CALL NENE
      CALL RHOS
      CALL VDOTS

      IF (NAV.GE.1.and.nrelax.eq.0) THEN
C     Calculate Adot (only if AV is on, and relaxation is off):                                                                                               
         CALL ADOTS
C     Advance entropies from original values:                                                                                                                 
!         DO I=1,N
!            A(I)=A(I)+DT*ADOT(I)
!         ENDDO
      ENDIF
C     Advance positions to full timestep:

      DO I=1,N
            Xold(I)=X(I)
            Yold(I)=Y(I)
            Zold(I)=Z(I)

            VXold(I)=VX(I)
            VYold(I)=VY(I)
            VZold(I)=VZ(I)
            
            Aold(I)=A(I)


       END DO
       DO I=1,N 
            K1_X(I)=VX(I)
            K1_Y(I)=VY(I)
            K1_Z(I)=VZ(I)

            K1_VX(I)=VXDOT(I)
            K1_VY(I)=VYDOT(I)
            K1_VZ(I)=VZDOT(I)

            A1(I)=ADOT(I)

       END DO     

 
       DO I=1,N
      
            X(I)=XOLD(I)+DTH*K1_X(I)
            Y(I)=YOLD(I)+DTH*K1_Y(I)
            Z(I)=ZOLD(I)+DTH*K1_Z(I)
            
            VX(I)=VXold(I)+DTh*K1_VX(I)
            VY(I)=VYold(I)+DTH*K1_VY(I)
            VZ(I)=VZold(I)+DTH*K1_VZ(I)
            
            A(I)=Aold(I)+DTH*A1(I)
      ENDDO
      
      CALL NENE
      CALL RHOS
      CALL VDOTS
      CALL ADOTS
      
       DO I=1,N

            K2_X(I)=VX(I)
            K2_Y(I)=VY(I)
            K2_Z(I)=VZ(I)
          
            K2_VX(I)=VXDOT(I)
            K2_VY(I)=VYDOT(I)
            K2_VZ(I)=VZDOT(I)
       
            A2(I)=ADOT(I)
         ENDDO
      
        DO I=1,N
            X(I)=Xold(I)+DT*K2_X(I)
            Y(I)=Yold(I)+DT*K2_Y(I)
            Z(I)=Zold(I)+DT*K2_Z(I)

      ENDDO


C     If binary relaxation problem, adjust center of mass positions:                                                                                          
      IF (NRELAX.GE.1) CALL CMADJ
C     Advance velocities to full-timestep:                                                                                                                    
      DO I=1,N
         VX(I)=VXold(I)+DT*K2_VX(I)
         VY(I)=VYold(I)+DT*K2_VY(I)
         VZ(I)=VZold(I)+DT*K2_VZ(I)
      
         A(I)=Aold(I)+A2(I)*DT
      ENDDO

      do i=1,n
         if(initgr.eq.3) then
c     this is taken from Blanchet, Damour, and Schaefer                                                                                                      
c     MNRAS 242, 289 (1990), Eq. (5.16),  with beta=A_i=0 in the Newtonian                                                                                   
c     limit. Note that our code's "v" is written there as "w",                                                                                               c     and our code's "u" is written as "v"                                                                                                                    
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
****************************************************************
      SUBROUTINE ADVANCE_RK4
      IMPLICIT NONE
*****************************************************************                                                                
c     Release 1.0                                                                                                                                         
c     Advance hydro quantities by a timestep, using second order Runge Kutta method                                                                  
c     Called by MAINIT                                                                                                                                       
c     Calls NENE,RHOS,ADOTS,CMADJ,VDOTS                                                                                                                      
***************************************************************                                                                                         
      INCLUDE 'spha.h'
      REAL UTINY
      PARAMETER (UTINY=1.E-10)

      REAL AO(NMAX), K1_X(NMAX),K1_Y(NMAX),K1_Z(NMAX)
      REAL DTH,K2_X(NMAX),K2_Y(NMAX),K2_Z(NMAX)
      REAL K1_VX(NMAX),K1_VY(NMAX),K1_VZ(NMAX)
      REAL K2_VX(NMAX),K2_VY(NMAX),K2_VZ(NMAX)
      REAL K3_X(NMAX),K3_Y(NMAX),K3_Z(NMAX)
      REAL K3_VX(NMAX),K3_VY(NMAX),K3_VZ(NMAX)
      REAL K4_X(NMAX),K4_Y(NMAX),K4_Z(NMAX)
      REAL K4_VX(NMAX),K4_VY(NMAX),K4_VZ(NMAX)
      REAL Xold(NMAX),Yold(NMAX),Zold(NMAX)
      REAL VXold(NMAX),VYold(NMAX),VZold(NMAX)
      REAL AOLD(NMAX),A1(NMAX),A2(NMAX)
      REAL A3(NMAX),A4(NMAX)
      INTEGER I

      DTH=0.5*DT
      CALL NENE
      CALL RHOS
      CALL VDOTS

      IF (NAV.GE.1.and.nrelax.eq.0) THEN
C     Calculate Adot (only if AV is on, and relaxation is off):                                                                                           
         CALL ADOTS
C     Advance entropies from original values:                                                                                                           
!         DO I=1,N
 !           A(I)=A(I)+DT*ADOT(I)
 !        ENDDO
      ENDIF
      DO I=1,N
            Xold(I)=X(I)
            Yold(I)=Y(I)
            Zold(I)=Z(I)

            VXold(I)=VX(I)
            VYold(I)=VY(I)
            VZold(I)=VZ(I)
 
           Aold(I)=A(I) 
      END do
      
      DO I=1,N
            K1_X(I)=VX(I)
            K1_Y(I)=VY(I)
            K1_Z(I)=VZ(I)

            K1_VX(I)=VXDOT(I)
            K1_VY(I)=VYDOT(I)
            K1_VZ(I)=VZDOT(I)

            A1(I)=ADOT(I)


            END DO

            DO I=1,N

            X(I)=XOLD(I)+DTH*K1_X(I)
            Y(I)=YOLD(I)+DTH*K1_Y(I)
            Z(I)=ZOLD(I)+DTH*K1_Z(I)

            VX(I)=VXOLD(I)+DTH*K1_VX(I)
            VY(I)=VYOLD(I)+DTH*K1_VY(I)
            VZ(I)=VZOLD(I)+DTH*K1_VZ(I)

            A(I)=AOLD(I)+DTH*A1(I)
      END DO
       
      CALL NENE
      CALL RHOS
      CALL VDOTS
      CALL ADOTS
       
      DO I=1,N
         K2_X(I)=VX(I)
         K2_Y(I)=VY(I)
         K2_Z(I)=VZ(I)

         K2_VX(I)=VXDOT(I)
         K2_VY(I)=VYDOT(I)
         K2_VZ(I)=VZDOT(I)
         
         A2(I)=ADOT(I)
      END DO

      DO I=1,N
          
          X(I)=Xold(I)+K2_X(I)*DTH
          Y(I)=Yold(I)+K2_Y(I)*DTH
          Z(I)=Zold(I)+K2_Z(I)*DTH
          
          VX(I)=VXOLD(I)+K2_VX(I)*DTH
          VY(I)=VYOLD(I)+K2_VY(I)*DTH
          VZ(I)=VZOLD(I)+K2_VZ(I)*DTH

          A(I)=AOLD(I)+A2(I)*DTH
      END DO

      CALL NENE
      CALL RHOS
      CALL VDOTS
      CALL ADOTS

      DO I=1,N

         K3_X(I)=VX(I)
         K3_Y(I)=VY(I)
         K3_Z(I)=VZ(I)


         K3_VX(I)=VXDOT(I)
         K3_VY(I)=VYDOT(I)
         K3_VZ(I)=VZDOT(I)

         A3(I)=ADOT(I)
      
         END DO

         DO I=1,N

          X(I)=Xold(I)+K3_X(I)*DT
          Y(I)=Yold(I)+K3_Y(I)*DT
          Z(I)=Zold(I)+K3_Z(I)*DT

          VX(I)=VXOLD(I)+K3_VX(I)*DT
          VY(I)=VYOLD(I)+K3_VY(I)*DT
          VZ(I)=VZOLD(I)+K3_VZ(I)*DT

          A(I)=AOLD(I)+A3(I)*DT

      END DO

      CALL NENE
      CALL RHOS
      CALL VDOTS
      CALL ADOTS
      
      DO I=1,N 

         K4_X(I)=VX(I)
         K4_Y(I)=VY(I)
         K4_Z(I)=VZ(I)


         K4_VX(I)=VXDOT(I)
         K4_VY(I)=VYDOT(I)
         K4_VZ(I)=VZDOT(I)

         A4(I)=ADOT(I)
      END DO

!       FINAL POSITIONS and VELOCITIES
      DO I=1,N
          
         X(I)=Xold(I)+DT*(K1_X(I)+2*K2_X(I)+2*K3_X(I)+K4_X(I))/6
         Y(I)=Yold(I)+DT*(K1_Y(I)+2*K2_Y(I)+2*K3_Y(I)+K4_Y(I))/6
         Z(I)=ZolD(I)+DT*(K1_Z(I)+2*K2_Z(I)+2*K3_Z(I)+K4_Z(I))/6

         VX(I)=VXold(I)+DT*(K1_VX(I)+2*K2_VX(I)+2*K3_VX(I)+K4_VX(I))/6
         VY(I)=VYold(I)+DT*(K1_VY(I)+2*K2_VY(I)+2*K3_VY(I)+K4_VY(I))/6
         VZ(I)=VZold(I)+DT*(K1_VZ(I)+2*K2_VZ(I)+2*K3_VZ(I)+K4_VZ(I))/6

         A(I)=AOLD(I)+DT*(A1(I)+2*A2(I)+2*A3(I)+A4(I))/6
      END DO
C     If binary relaxation problem, adjust center of mass positions:                                                                                          
      IF (NRELAX.GE.1) CALL CMADJ

      do i=1,n
         if(initgr.eq.3) then
c     this is taken from Blanchet, Damour, and Schaefer                                                                                                    
c     MNRAS 242, 289 (1990), Eq. (5.16),  with beta=A_i=0 in the Newtonian                                                                                 
C     limit. Note that our code's "v" is written there as "w",                                                                                              
C     and our code's "u" is written as "v"                                                                                                                   
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
!----------------------------------------------------------------------------------!
****************************************************************
      SUBROUTINE ADVANCE_SHOCK_NEW
      IMPLICIT NONE
*****************************************************************                                                                
c     Release 1.0                                                                                                                                         c     Advance hydro quantities by a timestep, using second order Runge Kutta method                                                                  
c     Called by MAINIT                                                                                                                                       
c     Calls NENE,RHOS,ADOTS,CMADJ,VDOTS                                                                                                                      
***************************************************************                                                                                         
      INCLUDE 'shock.h'
      INCLUDE 'spha.h'
!      INCLUDE 'sph.input'
!      OPEN(124,FILE='splash_ghost.txt')
      
      REAL UTINY
      PARAMETER (UTINY=1.E-10)

      REAL K1_X(NNMAX),K1_Y(NMAX),K1_Z(NMAX)
      REAL DTH,K2_X(NMAX),K2_Y(NMAX),K2_Z(NMAX)
      REAL K1_VX(NMAX),K1_VY(NMAX),K1_VZ(NMAX)
      REAL K2_VX(NMAX),K2_VY(NMAX),K2_VZ(NMAX)
      REAL K3_X(NMAX),K3_Y(NMAX),K3_Z(NMAX)
      REAL K3_VX(NMAX),K3_VY(NMAX),K3_VZ(NMAX)
      REAL K4_X(NMAX),K4_Y(NMAX),K4_Z(NMAX)
      REAL K4_VX(NMAX),K4_VY(NMAX),K4_VZ(NMAX)
      REAL Xold(NMAX),Yold(NMAX),Zold(NMAX)
      REAL VXold(NMAX),VYold(NMAX),VZold(NMAX)
      REAL AOLD(NMAX),A1(NMAX),A2(NMAX)
      REAL A3(NMAX),A4(NMAX),ENold(NMAX),u1
      REAL EN1(NMAX),EN2(NMAX),EN3(NMAX),EN4(NMAX)
      INTEGER I,J,k

      OPEN(124,FILE='splash_ghost.txt')
 !     OPEN(12,FILE='sph.input')
 !     READ(12,INPUT)
 !     CLOSE(12)


c     Number of particles in y and z directions 
      nwidth = 100 !2*nnopt
c     Number of particles in x direction 
      nx = 1000
c     Number of real particles
      np = nx
c     Number of particles+ghost particles 
      ntube=np
c     Number of particles on one end                                                                                                                        
      nend = nwidth
c     Total number of particles                                                                                                                              
      ntot = ntube+2*nwidth
      ntimestepper = -1
      space=1.0/(nx)

c      N= 10000
c      NENERGY=  2
 !     GAM=2.599999


      DTH=0.5*DT
      CALL NENE
      CALL RHOS
      CALL VDOTS
      IF (NAV.GE.1.and.nrelax.eq.0) THEN
         
         CALL EDOTS
         
      ENDIF 

!      do i=1,np
! 1    write(124,*),t,i          !,endot(I)
!      end do

      DO I=1,np
c Make sure we never reference 'xold' for indices > np in advance_shock
            Xold(I)=X(I)
            Yold(I)=Y(I)
            Zold(I)=Z(I)

            VXold(I)=VX(I)
            VYold(I)=VY(I)
            VZold(I)=VZ(I)

            if (nenergy.eq.1)then
               Aold(I)=A(I)
            end if

!            if (nenergy.eq.2)then
!               ENold(I)=EN(I)

!            end if
            
         END do

      DO I=1,np
            K1_X(I)=VX(I)
            K1_Y(I)=VY(I)
            K1_Z(I)=VZ(I)

            K1_VX(I)=VXDOT(I)
            K1_VY(I)=VYDOT(I)
            K1_VZ(I)=VZDOT(I)

            if (nenergy.eq.1)then
               A1(I)=ADOT(I)
            endif

!            if (nenergy.eq.2)then
!               EN1(I)=ENDOT(I)
!            endif

            END DO

      DO I=1,np
            X(I)=XOLD(I)+DTH*K1_X(I)
            VX(I)=VXOLD(I)+DTH*K1_VX(I)
            
            if (nenergy.eq.1)then
               A(I)=AOLD(I)+DTH*A1(I)
               call e_from_a
            end if
            
 !           if (nenergy.eq.2)then
 !           EN(I)=ENOLD(I)+DTH*EN1(I)
  !          if(en(i).le.0)write(6,*)'as1:bad en',i,en(i),enold(i),en1(i)
   !         call a_from_e
   
 !       end if

      END DO

      CALL NENE
      CALL RHOS
      CALL VDOTS
      CALL EDOTS
       
      DO I=1,np
         K2_X(I)=VX(I)
         K2_Y(I)=VY(I)
         K2_Z(I)=VZ(I)

         K2_VX(I)=VXDOT(I)
         K2_VY(I)=VYDOT(I)
         K2_VZ(I)=VZDOT(I)
         
         if (nenergy.eq.1)then
            A2(I)=ADOT(I)
         end if
         
  !       if (nenergy.eq.2)then
  !        EN2(I)=ENDOT(I)
   !   end if
      
      END DO

    

      DO I=1,np
          X(I)=Xold(I)+K2_X(I)*DTH
          VX(I)=VXOLD(I)+K2_VX(I)*DTH
          
          if (nenergy.eq.1)then          
             A(I)=AOLD(I)+A2(I)*DTH
             call e_from_a
          end if
          
  !        if (nenergy.eq.2)then
  !           EN(I)=ENOLD(I)+EN2(I)*DTH
  !          if(en(i).le.0)write(6,*)'as2:bad en',i,en(i),enold(i),en2(i)
 !        call a_from_e
 !        end if 
          
       END DO

      CALL NENE
      CALL RHOS
      CALL VDOTS
      CALL EDOTS

      DO I=1,np

         K3_X(I)=VX(I)
         K3_Y(I)=VY(I)
         K3_Z(I)=VZ(I)

         K3_VX(I)=VXDOT(I)
         K3_VY(I)=VYDOT(I)
         K3_VZ(I)=VZDOT(I)

         
         if (nenergy.eq.1)then
            A3(I)=ADOT(I)
         end if
         
  !       if (nenergy.eq.2)then
  !          EN3(I)=ENDOT(I)
  !       end if
      END DO

      DO I=1,np
         
         X(I)=Xold(I)+K3_X(I)*DT
          VX(I)=VXOLD(I)+K3_VX(I)*DT
          
          if (nenergy.eq.1)then
             A(I)=AOLD(I)+A3(I)*DT
             call e_from_a
          end if
          
   !       if (nenergy.eq.2)then
   !          EN(I)=ENOLD(I)+EN3(I)*DT
   !         if(en(i).le.0)write(6,*)'as3:bad en',i,en(i),enold(i),en3(i)
   !          call a_from_e
  !    end if
       
       END DO

      CALL NENE
      CALL RHOS
      CALL VDOTS
      CALL EDOTS
      
      DO I=1,np

         K4_X(I)=VX(I)
         K4_Y(I)=VY(I)
         K4_Z(I)=VZ(I)


         K4_VX(I)=VXDOT(I)
         K4_VY(I)=VYDOT(I)
         K4_VZ(I)=VZDOT(I)

         if (nenergy.eq.1)then
            A4(I)=ADOT(I)
         end if
         
   !      if (nenergy.eq.2)then
   !         EN4(I)=ENDOT(I)
    !     end if
      END DO

!       FINAL POSITIONS and VELOCITIES
      DO I=1,np
      
         X(I)=Xold(I)+DT*(K1_X(I)+2*K2_X(I)+2*K3_X(I)+K4_X(I))/6
         Y(I)=Yold(I)
         Z(I)=ZolD(I)

         VX(I)=VXold(I)+DT*(K1_VX(I)+2*K2_VX(I)+2*K3_VX(I)+K4_VX(I))/6
         VY(I)=0
         VZ(I)=0

         if (nenergy.eq.1)then
            A(I)=AOLD(I)+DT*(A1(I)+2*A2(I)+2*A3(I)+A4(I))/6 
!        
            call e_From_a
         end if
         
!         if (nenergy.eq.2)then
!            EN(I)=ENOLD(I)+DT*(EN1(I)+2*EN2(I)+2*EN3(I)+EN4(I))/6
!            if(en(i).le.0)write(6,*)'as4:bad en',i,en(i),enold(i),en4(i)

!            u = EN(I)-0.5*(vx(I)**2+vy(I)**2+vz(I)**2)

c  Move this block to SUBROUTINE a_from_e
 !           if (u1.lt.0)then
 !              write(6,*)"as:Bad u value!", i,u,EN(I),vx(I),vy(I),vz(I)
  !          end if

   !         call a_from_e
  !       end if
      
      END DO

      
!      write(*,*),"Gamma",GAM

      do i=ntube,ntot
         x(i)=x(i)+dt*vx(I)
         y(i)=y(i)
         z(i)=z(i)
      end do


C     If binary relaxation problem, adjust center of mass positions:                                                                                          
      IF (NRELAX.GE.1) CALL CMADJ

      do i=1,n
         if(initgr.eq.3) then
c     this is taken from Blanchet, Damour, and Schaefer                                                                                                    
c     MNRAS 242, 289 (1990), Eq. (5.16),  with beta=A_i=0 in the Newtonian                                                                                 
C     limit. Note that our code's "v" is written there as "w",                                                                                              
C     and our code's "u" is written as "v"                                                                                                                   
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

****************************************************************
      SUBROUTINE ADVANCE_SHOCK
      IMPLICIT NONE
*****************************************************************                                                                
c     Release 1.0                                                                                                                                         c     Advance hydro quantities by a timestep, using second order Runge Kutta method                                                                  
c     Called by MAINIT                                                                                                                                       
c     Calls NENE,RHOS,ADOTS,CMADJ,VDOTS                                                                                                                      
***************************************************************                                                                                         
      INCLUDE 'shock.h'
      INCLUDE 'spha.h'
      !OPEN(124,FILE='splash_ghost.txt')

      REAL UTINY
      PARAMETER (UTINY=1.E-10)

      REAL AO(NMAX), K1_X(NMAX),K1_Y(NMAX),K1_Z(NMAX)
      REAL DTH,K2_X(NMAX),K2_Y(NMAX),K2_Z(NMAX)
      REAL K1_VX(NMAX),K1_VY(NMAX),K1_VZ(NMAX)
      REAL K2_VX(NMAX),K2_VY(NMAX),K2_VZ(NMAX)
      REAL K3_X(NMAX),K3_Y(NMAX),K3_Z(NMAX)
      REAL K3_VX(NMAX),K3_VY(NMAX),K3_VZ(NMAX)
      REAL K4_X(NMAX),K4_Y(NMAX),K4_Z(NMAX)
      REAL K4_VX(NMAX),K4_VY(NMAX),K4_VZ(NMAX)
      REAL Xold(NMAX),Yold(NMAX),Zold(NMAX)
      REAL VXold(NMAX),VYold(NMAX),VZold(NMAX)
      REAL AOLD(NMAX),A1(NMAX),A2(NMAX)
      REAL A3(NMAX),A4(NMAX),ENold(NMAX),u1
      REAL EN1(NMAX),EN2(NMAX),EN3(NMAX),EN4(NMAX)
      INTEGER I,J,k

c     Number of particles in y and z directions                                                                                                               
      nwidth = 100 !2*nnopt                                                                                                                                   
c     Number of particles in x direction                                                                                                                      
      nx = 1000
c     Number of real particles                                                                                                                                
      np = nx
c     Number of particles+ghost particles                                                                                                                     
      ntube=np
c     Number of particles on one end                                                                                                                          
      nend = nwidth
c     Total number of particles                                                                                                                              
      ntot = ntube+2*nwidth
      ntimestepper = -1
      space=1.0/(nx)

      DTH=0.5*DT
      CALL NENE
      CALL RHOS
      CALL VDOTS

      IF (NAV.GE.1.and.nrelax.eq.0) THEN

         CALL ADOTS

      ENDIF 
 
       DO I=1,np
c Make sure we never reference 'xold' for indices > np in advance_shock                                                                                       
            Xold(I)=X(I)
            Yold(I)=Y(I)
            Zold(I)=Z(I)

            VXold(I)=VX(I)
            VYold(I)=VY(I)
            VZold(I)=VZ(I)

            if (nenergy.eq.1)then
               Aold(I)=A(I)
            end if
             if (nenergy.eq.2)then
               ENold(I)=EN(I)
            end if

         END do

      DO I=1,np
            K1_X(I)=VX(I)
            K1_Y(I)=VY(I)
            K1_Z(I)=VZ(I)

            K1_VX(I)=VXDOT(I)
            K1_VY(I)=VYDOT(I)
            K1_VZ(I)=VZDOT(I)

            
         END DO

         do i=1,np
            X(I)=XOLD(I)+DTH*K1_X(I)
            VX(I)=VXOLD(I)+DTH*K1_VX(I)
            if (nenergy.eq.1)then
               A(I)=AOLD(I)+DTH*A1(I)
               call e_from_a
            end if
            if (nenergy.eq.2)then                                                                                                                            
               EN(I)=ENOLD(I)+DTH*EN1(I)                                                                                                                      
               call a_from_e                                                                                                                                  
            end if

         end do
      CALL NENE
      CALL RHOS
      CALL VDOTS
      CALL EDOTS

      DO I=1,np
         K2_X(I)=VX(I)
         K2_Y(I)=VY(I)
         K2_Z(I)=VZ(I)

         K2_VX(I)=VXDOT(I)
         K2_VY(I)=VYDOT(I)
         K2_VZ(I)=VZDOT(I)

         if (nenergy.eq.1)then
            A2(I)=ADOT(I)
         end if

         if (nenergy.eq.2)then                                                                                                                               
          EN2(I)=ENDOT(I)                                                                                                                                    
        end if                                                                                                                                                 

      END DO


      DO I=1,np
          X(I)=Xold(I)+K2_X(I)*DTH
          VX(I)=VXOLD(I)+K2_VX(I)*DTH

          if (nenergy.eq.1)then
             A(I)=AOLD(I)+A2(I)*DTH
             call e_from_a
          end if

          if (nenergy.eq.2)then
             EN(I)=ENOLD(I)+EN2(I)*DTH
             call a_from_e
          end if
          
       END DO
!-----------------------------------------------------------------!
      CALL NENE
      CALL RHOS
      CALL VDOTS
      CALL EDOTS

      DO I=1,np

         K3_X(I)=VX(I)
         K3_Y(I)=VY(I)
         K3_Z(I)=VZ(I)

         K3_VX(I)=VXDOT(I)
         K3_VY(I)=VYDOT(I)
         K3_VZ(I)=VZDOT(I)


         if (nenergy.eq.1)then
            A3(I)=ADOT(I)
         end if

        if (nenergy.eq.2)then                                                                                                                               
            EN3(I)=ENDOT(I)                                                                                                                                  
         end if                                                                                                                                              
      END DO

      DO I=1,np

         X(I)=Xold(I)+K3_X(I)*DT
          VX(I)=VXOLD(I)+K3_VX(I)*DT

          if (nenergy.eq.1)then
             A(I)=AOLD(I)+A3(I)*DT
             call e_from_a
          end if

          if (nenergy.eq.2)then
             EN(I)=ENOLD(I)+EN3(I)*DT
             call a_from_e
          end if

      END DO

      CALL NENE
      CALL RHOS
      CALL VDOTS
      CALL EDOTS
      
      DO I=1,np

         K4_X(I)=VX(I)
         K4_Y(I)=VY(I)
         K4_Z(I)=VZ(I)


         K4_VX(I)=VXDOT(I)
         K4_VY(I)=VYDOT(I)
         K4_VZ(I)=VZDOT(I)

         if (nenergy.eq.1)then
            A4(I)=ADOT(I)
         end if

        if (nenergy.eq.2)then                                                                                                                               
            EN4(I)=ENDOT(I)                                                                                                                                  
         end if                                                                                                                                              
      END DO   
 
      DO I=1,np

         X(I)=Xold(I)+DT*(K1_X(I)+2*K2_X(I)+2*K3_X(I)+K4_X(I))/6
         Y(I)=Yold(I)
         Z(I)=ZolD(I)

         VX(I)=VXold(I)+DT*(K1_VX(I)+2*K2_VX(I)+2*K3_VX(I)+K4_VX(I))/6
         VY(I)=0
         VZ(I)=0

         if (nenergy.eq.1)then
            A(I)=AOLD(I)+DT*(A1(I)+2*A2(I)+2*A3(I)+A4(I))/6
!                                                                                                                                                             
            call e_From_a
         end if

         if (nenergy.eq.2)then
            EN(I)=ENOLD(I)+DT*(EN1(I)+2*EN2(I)+2*EN3(I)+EN4(I))/6
            
            call a_From_e
         end if

!       FINAL POSITIONS and VELOCITIES

!         if(I.gt.ntube)then!

!         X(I)=Xold(I)
!         Y(I)=Yold(I)       
!         Z(I)=ZolD(I)   

!         VX(I)=VXold(I)
!         VY(I)=0
!         VZ(I)=0 

!         A(I)=AOLD(I)
!         write(124,*),I,t,y(I),vy(I),z(I),vz(I)
!      end if
      END DO
      do i=ntube,ntot
         x(i)=x(i)+dt*vx(I)
         y(i)=y(i)
         z(i)=z(i)
      end do


C     If binary relaxation problem, adjust center of mass positions:                                                                                          
      IF (NRELAX.GE.1) CALL CMADJ

      do i=1,n
         if(initgr.eq.3) then
c     this is taken from Blanchet, Damour, and Schaefer                                                                                                    
c     MNRAS 242, 289 (1990), Eq. (5.16),  with beta=A_i=0 in the Newtonian                                                                                 
C     limit. Note that our code's "v" is written there as "w",                                                                                              
C     and our code's "u" is written as "v"                                                                                                                   
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
