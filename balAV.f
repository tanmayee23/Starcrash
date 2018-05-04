
**************************************************************
c     This file is part of the StarCrash code
c     Version 1.0
c
c     Copyright 2003, by Joshua Faber
c     This code is publicly available under the terms of the
c     GNU General Public License.  See our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

      SUBROUTINE BALADOTS
      IMPLICIT NONE
*****************************************************
c     Release 1.0
c     This routine calculates ADOT using the Balsara AV formulation.
c     See Balsara J.Comp.Phys. 121, 357 (1995) and
c     Lombardi, Rasio, and Shapiro J.Comp.Phys 152, 687 (1999)
c     Called by ADOTS
c     Calls CALCDIVV,GETCURLV
************************************************** 

C Evaluate right-hand sides of dA/dt equations:
      INCLUDE 'spha.h'
      INCLUDE 'mpif.h' 
      REAL DIVV(NMAX),CURLVX(NMAX),CURLVY(NMAX),CURLVZ(NMAX)
      COMMON/COMMDIVV/DIVV
      COMMON/COMCURLV/CURLVX,CURLVY,CURLVZ

      REAL DWIJ(NNMAX),PIJ(NNMAX),SIJ(NNMAX),DIV(NNMAX),myadot(N)
      REAL hpi,h2,h5,ami,ci2,rhoi,por2i,divvi,abscurli
      REAL vxi,vyi,vzi,r2,hij,ci,cj,fi,fj,udbij,adoti
      INTEGER i,ii,in,itab,j,ierr

      divvi=0.
      abscurli=0.

C  Initialize dA/dt before accumulating:
      DO I=1,N
        myADOT(I)=0.
      ENDDO

C  dA/dt=0 if GAM=1. (isothermal calculation):
      IF (GAM.EQ.1.) RETURN

      CALL CALCDIVV
      CALL GETCURLV

C  For each particle:
      DO I=n_lower,n_upper
         II=I-n_lower+1
         HPI=HP(I)
         H2=HPI**2
         H5=HPI*H2**2
         AMI=AM(I)
         CI2=GAM*A(I)*RHO(I)**(GAM-1.)
         RHOI=RHO(I)
         POR2I=POR2(I)
         DIVVI=DIVV(I)
         ABSCURLI=SQRT(CURLVX(I)**2+CURLVY(I)**2+CURLVZ(I)**2)

         VXI=VX(I)
         VYI=VY(I)
         VZI=VZ(I)
         
C     Calculate gradWij's to all neighbors:
         DO IN=1,NN(I)
            R2=XIJ(IN,II)**2.+YIJ(IN,II)**2.+ZIJ(IN,II)**2.
            ITAB=INT(CTAB*R2/H2)+1
            DWIJ(IN)=DWTAB(ITAB)*2.0*PI/3.0/HPI**3 !H5
         ENDDO
         
C     Calculate (Xi-Xj)*(VXi-VXj) to all neighbors:
         DO IN=1,NN(I)
            J=NNI(IN,II)
            DIV(IN)=XIJ(IN,II)*(vx(i)-vx(j))+
     $           YIJ(IN,II)*(vy(i)-vy(j))+
     $           ZIJ(IN,II)*(vz(i)-vz(j))
         ENDDO
         
         IF (NAV.NE.0) THEN
C     Calculate artificial viscosity contribution to pressure gradients:
            DO IN=1,NN(I)
               IF (DIV(IN).LT.0.) THEN
                  J=NNI(IN,II)
                  R2=XIJ(IN,II)**2.+YIJ(IN,II)**2.+ZIJ(IN,II)**2.
                  HIJ=0.5*(HPI+HP(J))
                  CI=SQRT(CI2)
                  CJ=SQRT(GAM*A(J)*RHO(J)**(GAM-1.))
C     Calculate Dinshaw Balsara's (DB) \mu_{ij}
                  FI=ABS(DIVVI)/(ABS(DIVVI)+ABSCURLI+
     $                 0.00001*CI/HPI)
                  FJ=ABS(DIVV(J))/(ABS(DIVV(J))+
     $                 SQRT(CURLVX(J)**2+CURLVY(J)**2+
     $                 CURLVZ(J)**2)+0.00001*CJ/HP(J))
                  UDBIJ=HIJ*DIV(IN)/((CI+CJ)*(R2+ETA2*HIJ**2))
     $                 *(FI+FJ)
                  PIJ(IN)=0.5*GAM*(POR2I+POR2(J))*
     $                 (-ALPHA*UDBIJ+BETA*UDBIJ**2)
               ELSE
                  PIJ(IN)=0.
               ENDIF
            ENDDO
            DO IN=1,NN(I)
               SIJ(IN)=0.5*(GAM-1.)*PIJ(IN)*DWIJ(IN)*DIV(IN)
            ENDDO
         ELSE
            DO IN=1,NN(I)
               SIJ(IN)=0.
            ENDDO
         ENDIF
         
C     Finally compute ADOT:
C     "Gather" part of the sum (i-j pair contributes to adot_i):
         ADOTI=0.
         DO IN=1,NN(I)
            J=NNI(IN,II)
            ADOTI=ADOTI+AM(J)*SIJ(IN)/RHOI**(GAM-1.)
         ENDDO
         myADOT(I)=myADOT(I)+0.5*ADOTI
         
C     "Scatter" part of the sum (i-j pair contributes to adot_j):
         DO IN=1,NN(I)
            J=NNI(IN,II)
            myADOT(J)=myADOT(J)+0.5*AMI*SIJ(IN)/RHO(J)**(GAM-1.)
         ENDDO
         
      ENDDO

c     Use MPI_ALLREDUCE to collect subtotals, sum them, and redistribute
c     to all processes

      CALL MPI_ALLREDUCE(myadot,adot,n,MPI_REAL,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)
      
      RETURN
      END
************************************************************************
      SUBROUTINE BALVDOTS
      IMPLICIT NONE
****************************************************************
c     Release 1.0
c     This routine calculates VDOT using the Balsara AV formulation.
c     See Balsara J.Comp.Phys. 121, 357 (1995) and
c     Lombardi, Rasio, and Shapiro J.Comp.Phys 152, 687 (1999)
c     Called by VDOTS
c     Calls CALCDIVV,GETCURLV,DOQGRAV,QGRAV
************************************************** 
      INCLUDE 'spha.h'                                         
      INCLUDE 'mpif.h'

      REAL UIJMAX(NMAX),DIVV(NMAX)
      REAL CURLVX(NMAX),CURLVY(NMAX),CURLVZ(NMAX)
      COMMON/UIJMAX/ UIJMAX
      COMMON/COMMDIVV/DIVV
      COMMON/COMCURLV/CURLVX,CURLVY,CURLVZ
      REAL DWIJ(NNMAX),PIJ(NNMAX),SXIJ(NNMAX),SYIJ(NNMAX),SZIJ(NNMAX),
     $     DIV(NNMAX),SXIJ_DISS(NNMAX),SYIJ_DISS(NNMAX),
     $     SZIJ_DISS(NNMAX),DIV_V
      REAL myuijmax(N),myvxdot(N),myvydot(N),myvzdot(N)

      REAL hpi,h2,h5,ami,ci2,rhoi,por2i,vxi,vyi,vzi,divvi
      REAL abscurli,r2,hij,ci,cj,fi,fj,udbij,csij,csij_diss
      real vxdoti,vydoti,vzdoti,al,rhoavg,vsig,r1
      INTEGER i,ii,in,itab,j,ierr

      divvi=0.
      abscurli=0.
      al=0.651
C  Initialize dv/dt before accumulating:
      DO I=1,N
        myVXDOT(I)=0.
        myVYDOT(I)=0.
        myVZDOT(I)=0.
      ENDDO

c     this routine is called when AV is off, don't do extra work
      if(nav.ge.1.and.nrelax.eq.0) then
         CALL CALCDIVV
         CALL GETCURLV
      endif

      do i=1,n
         uijmax(i)=0.
         myuijmax(i)=0.
      enddo

C  For each particle:
      DO I=n_lower,n_upper
         II=I-n_lower+1
         HPI=HP(I)
         H2=HPI**2
         H5=HPI*H2**2
         AMI=AM(I)
         CI2=GAM*A(I)*RHO(I)**(GAM-1.)

         RHOI=RHO(I)
         POR2I=POR2(I)
         VXI=VX(I)
         VYI=VY(I)
         VZI=VZ(I)
         if(nav.ne.0)then
            DIVVI=DIVV(I)
            ABSCURLI=SQRT(CURLVX(I)**2+CURLVY(I)**2+CURLVZ(I)**2)
         endif
C     Calculate gradWij's to all neighbors:
         DO IN=1,NN(I)
            R2=XIJ(IN,II)**2.+YIJ(IN,II)**2.+ZIJ(IN,II)**2.
            ITAB=INT(CTAB*R2/H2)+1
            DWIJ(IN)=DWTAB(ITAB)*2.0*PI/3.0/HPI**3   !/H5
         ENDDO
         
         IF (NAV.NE.0.and.nrelax.eq.0) THEN
C     Calculate artificial viscosity contribution to pressure gradients:
            DO IN=1,NN(I)
               J=NNI(IN,II)
               DIV(IN)=XIJ(IN,II)*(vx(i)-vx(j))+
     $              YIJ(IN,II)*(vy(i)-vy(j))+
     $              ZIJ(IN,II)*(vz(i)-vz(j))     
            ENDDO
            DO IN=1,NN(I)
               IF (DIV(IN).LT.0.) THEN
                  J=NNI(IN,II)
                  
                  R2=XIJ(IN,II)**2.+YIJ(IN,II)**2.+ZIJ(IN,II)**2.
                  HIJ=0.5*(HPI+HP(J))
                  CI=SQRT(CI2)
                  CJ=SQRT(GAM*A(J)*RHO(J)**(GAM-1.))
C     Calculate Dinshaw Balsara's (DB) \mu_{ij}
                  FI=ABS(DIVVI)/(ABS(DIVVI)+ABSCURLI+
     $                 0.00001*CI/HPI)
                  FJ=ABS(DIVV(J))/(ABS(DIVV(J))+
     $                 SQRT(CURLVX(J)**2+CURLVY(J)**2+
     $                 CURLVZ(J)**2)+0.00001*CJ/HP(J))
                  UDBIJ=HIJ*DIV(IN)/((CI+CJ)*(R2+ETA2*HIJ**2))
     $                 *(FI+FJ)
                  myUIJMAX(I)=MAX(myUIJMAX(I),ABS(UDBIJ))
                  PIJ(IN)=0.5*GAM*(POR2I+POR2(J))*
     $                 (-ALPHA*UDBIJ+BETA*UDBIJ**2)
               ELSE
                  PIJ(IN)=0.
               ENDIF
            ENDDO
            DO IN=1,NN(I)
               J=NNI(IN,II)
               if (i.ne.j)then 
               CSIJ=(POR2I+POR2(J))*DWIJ(IN)
!              CSIJ=(POR2I+POR2(J)+PIJ(IN))*DWIJ(IN
               SXIJ(IN)=CSIJ*XIJ(IN,II)
               SYIJ(IN)=CSIJ*YIJ(IN,II)
               SZIJ(IN)=CSIJ*ZIJ(IN,II)
               
C            With dissapation term
               CI=SQRT(CI2)
               CJ=SQRT(GAM*A(J)*RHO(J)**(GAM-1.))
             !  R2=XIJ(IN,II)**2.+YIJ(IN,II)**2.+ZIJ(IN,II)**2.
               R1=SQRT((x(i)-x(j))**2+(y(I)-Y(J))**2+(Z(I)-Z(J))**2)
               DIV_V=(X(I)-x(J))*(vx(i)-vx(j))+
     $              (Y(I)-y(J))*(vy(i)-vy(j))+
     $              (Z(I)-z(J))*(vz(i)-vz(j))
                VSIG=CI+CJ-(DIv_V/R1)
                RHOAVG=(RHO(I)+RHO(J))/2
                CSIJ_DISS=(AL*VSIG*(DIV_V/R1)*DWIJ(IN))/RHOAVG
!                write(6,*),IN,vsig,r1,i,j,r2
                SXIJ_DISS(IN)=CSIJ_DISS*XIJ(IN,II)
                SYIJ_DISS(IN)=CSIJ_DISS*YIJ(IN,II)
                SZIJ_DISS(IN)=CSIJ_DISS*ZIJ(IN,II)
             end if
             ENDDO
         ELSE
            DO IN=1,NN(I)
               J=NNI(IN,II)
               CSIJ=(POR2I+POR2(J))*DWIJ(IN)
               SXIJ(IN)=CSIJ*XIJ(IN,II)
               SYIJ(IN)=CSIJ*YIJ(IN,II)
               SZIJ(IN)=CSIJ*ZIJ(IN,II)

C            With dissapation term                                                                                                                            
               CI=SQRT(CI2)
               CJ=SQRT(GAM*A(J)*RHO(J)**(GAM-1.))
             !  R2=XIJ(IN,II)**2.+YIJ(IN,II)**2.+ZIJ(IN,II)**2.                                                                                               
               R1=SQRT((x(i)-x(j))**2+(y(I)-Y(J))**2+(Z(I)-Z(J))**2)
               DIV-V=(X(I)-x(J))*(vx(i)-vx(j))+
     $              (Y(I)-y(J))*(vy(i)-vy(j))+
     $              (Z(I)-z(J))*(vz(i)-vz(j))
                VSIG=CI+CJ-(DIV_V/R1)
                RHOAVG=(RHO(I)+RHO(J))/2
                CSIJ_DISS=(AL*VSIG*(DIV_V/R1)*DWIJ(IN))/RHOAVG
!                write(6,*),IN,vsig,r1,i,j,r2                                                                                                                 
                SXIJ_DISS(IN)=CSIJ_DISS*XIJ(IN,II)
                SYIJ_DISS(IN)=CSIJ_DISS*YIJ(IN,II)
                SZIJ_DISS(IN)=CSIJ_DISS*ZIJ(IN,II)
            ENDDO
         ENDIF
         
C     Finally compute each component of dv/dt:
C     "Gather" part of the sum (i-j pair contributes to vdot_i):
         VXDOTI=0.
         VYDOTI=0.
         VZDOTI=0.
         DO IN=1,NN(I)
            J=NNI(IN,II)
            VXDOTI=VXDOTI-AM(J)*SXIJ(IN)+AM(J)*SXIJ_DISS(IN)
            VYDOTI=VYDOTI-AM(J)*SYIJ(IN)+AM(J)*SYIJ_DISS(IN)
            VZDOTI=VZDOTI-AM(J)*SZIJ(IN)+AM(J)*SZIJ_DISS(IN)
         ENDDO
         myVXDOT(I)=myVXDOT(I)+0.5*VXDOTI
         myVYDOT(I)=myVYDOT(I)+0.5*VYDOTI
         myVZDOT(I)=myVZDOT(I)+0.5*VZDOTI
         
C     "Scatter" part of the sum (i-j pair contributes to vdot_j):
         DO IN=1,NN(I)
            J=NNI(IN,II)
            myVXDOT(J)=myVXDOT(J)+0.5*AMI*SXIJ(IN)-0.5*AMI*SXIJ_DISS(IN)
            myVYDOT(J)=myVYDOT(J)+0.5*AMI*SYIJ(IN)-0.5*AMI*SYIJ_DISS(IN)
            myVZDOT(J)=myVZDOT(J)+0.5*AMI*SZIJ(IN)-0.5*AMI*SZIJ_DISS(IN)
         ENDDO
         
      ENDDO

c     Use MPI_ALLREDUCE to collect subtotals, sum them, and redistribute
c     to all processes

      CALL MPI_ALLREDUCE(myvxdot,vxdot,n,MPI_REAL,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(myvydot,vydot,n,MPI_REAL,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(myvzdot,vzdot,n,MPI_REAL,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(myuijmax,uijmax,n,MPI_REAL,MPI_MAX,
     $     MPI_COMM_WORLD,ierr)

C Calculate self gravity of fluid:
      IF (NGR.GT.0) THEN 
c     doqgrav is the version with radiation reaction on, qgrav with it off
         if(ngravrad.eq.1.and.nrelax.eq.0) then            
            call doqgrav
         else
           CALL QGRAV
         endif
         DO I=1,N
            VXDOT(I)=VXDOT(I)+GX(I)
            VYDOT(I)=VYDOT(I)+GY(I)
            VZDOT(I)=VZDOT(I)+GZ(I)
            if(ngravrad.eq.1) then
               vxdot(i)=vxdot(i)+freacx(i)
               vydot(i)=vydot(i)+freacy(i)
               vzdot(i)=vzdot(i)+freacz(i)
            endif
         ENDDO
      ENDIF

C     If this is a relaxation calculation, add centrifugal and drag forces:
      IF (NRELAX.NE.0) CALL RELAX                                      

      RETURN
      END

************************************************************************
      SUBROUTINE ENDOTS
      IMPLICIT NONE
****************************************************************
c     Release 1.0
c     This routine calculates VDOT using the Balsara AV formulation.
c     See Balsara J.Comp.Phys. 121, 357 (1995) and
c     Lombardi, Rasio, and Shapiro J.Comp.Phys 152, 687 (1999)
c     Called by VDOTS
c     Calls CALCDIVV,GETCURLV,DOQGRAV,QGRAV
************************************************** 
      INCLUDE 'spha.h'                                         
      INCLUDE 'mpif.h'

      REAL UIJMAX(NMAX),DIVV(NMAX)
      REAL CURLVX(NMAX),CURLVY(NMAX),CURLVZ(NMAX)
      COMMON/UIJMAX/ UIJMAX
      COMMON/COMMDIVV/DIVV
      COMMON/COMCURLV/CURLVX,CURLVY,CURLVZ
      REAL DWIJ(NNMAX),PIJ(NNMAX),DIV(NNMAX),SEIJ(NNMAX)
      REAL myuijmax(N),myendot(N)

      REAL hpi,h2,h5,ami,ci2,rhoi,por2i,vxi,vyi,vzi,divvi,e_star_j
      REAL abscurli,r2,hij,ci,cj,fi,fj,udbij,csxij,csyij,cszij
      real endoti,vxijavg,vyijavg,vzijavg,csxij_diss,csyij_diss
      real cszij_diss,u_I,u_J,div_I,div_J,vsig,vsig_mu,e_star_i
      real al,rhoavg,r1 
      INTEGER i,ii,in,itab,j,ierr

      write(6,*)'In ENDOTS'

      divvi=0.
      abscurli=0.
      al=0.651
      
C  Initialize dv/dt before accumulating:
      DO I=1,N
        myENDOT(I)=0.
      ENDDO

c     this routine is called when AV is off, don't do extra work
      if(nav.ge.1.and.nrelax.eq.0) then
         write(6,*)'calling curl and div'
         CALL CALCDIVV
         CALL GETCURLV
         write(6,*)'callied curl and div'
      endif

      do i=1,n
         uijmax(i)=0.
         myuijmax(i)=0.
      enddo

C  For each particle:
      DO I=n_lower,n_upper
         II=I-n_lower+1
         HPI=HP(I)
         H2=HPI**2
         H5=HPI*H2**2
         AMI=AM(I)

         if(rho(i).le.0.0)write(6,*)'nonpositive rho!',i,rho(i)

         CI2=GAM*A(I)*RHO(I)**(GAM-1.)
         RHOI=RHO(I)

         if(por2(i).le.0.0)write(6,*)'nonpositive por2!',i,por2(i)

         POR2I=POR2(I)
         VXI=VX(I)
         VYI=VY(I)
         VZI=VZ(I)
         if(nav.ne.0)then
            DIVVI=DIVV(I)
            ABSCURLI=SQRT(CURLVX(I)**2+CURLVY(I)**2+CURLVZ(I)**2)
         endif
C     Calculate gradWij's to all neighbors:
         DO IN=1,NN(I)
            R2=XIJ(IN,II)**2.+YIJ(IN,II)**2.+ZIJ(IN,II)**2.
            ITAB=INT(CTAB*R2/H2)+1
            DWIJ(IN)=DWTAB(ITAB)*2.0*PI/3.0/HPI**3 !/H5
         ENDDO
         
         IF (NAV.NE.0.and.nrelax.eq.0) THEN
C     Calculate artificial viscosity contribution to pressure gradients:
            DO IN=1,NN(I)
               J=NNI(IN,II)
               DIV(IN)=XIJ(IN,II)*(vx(i)-vx(j))+
     $              YIJ(IN,II)*(vy(i)-vy(j))+
     $              ZIJ(IN,II)*(vz(i)-vz(j))
            ENDDO
            DO IN=1,NN(I)
               IF (DIV(IN).LT.0.) THEN
                  J=NNI(IN,II)
                  
                  R2=XIJ(IN,II)**2.+YIJ(IN,II)**2.+ZIJ(IN,II)**2.
                  HIJ=0.5*(HPI+HP(J))
                  CI=SQRT(CI2)
                  CJ=SQRT(GAM*A(J)*RHO(J)**(GAM-1.))
C     Calculate Dinshaw Balsara's (DB) \mu_{ij}
                  FI=ABS(DIVVI)/(ABS(DIVVI)+ABSCURLI+
     $                 0.00001*CI/HPI)
                  FJ=ABS(DIVV(J))/(ABS(DIVV(J))+
     $                 SQRT(CURLVX(J)**2+CURLVY(J)**2+
     $                 CURLVZ(J)**2)+0.00001*CJ/HP(J))
                  UDBIJ=HIJ*DIV(IN)/((CI+CJ)*(R2+ETA2*HIJ**2))
     $                 *(FI+FJ)
                  myUIJMAX(I)=MAX(myUIJMAX(I),ABS(UDBIJ))
                  PIJ(IN)=0.5*GAM*(POR2I+POR2(J))*
     $                 (-ALPHA*UDBIJ+BETA*UDBIJ**2)
               ELSE
                  PIJ(IN)=0.
               ENDIF
            ENDDO
            DO IN=1,NN(I)
               J=NNI(IN,II)
               IF (I.NE.J)THEN
                  VXIJAVG=(VX(I)+VX(J))/2
                  VYIJAVG=(VY(I)+VY(J))/2
                  VZIJAVG=(VZ(I)+VZ(J))/2
                  CSXIJ=(POR2I*VX(J)+POR2(J)*VX(I))!+PIJ(IN)*VXIJAVG)
                  CSYIJ=(POR2I*VY(J)+POR2(J)*VY(I))!+PIJ(IN)*VYIJAVG)
                  CSZIJ=(POR2I*VZ(J)+POR2(J)*VZ(I))!+PIJ(IN)*VZIJAVG)
!                  SEIJ(IN)=(CSXIJ*XIJ(IN,II)+CSYIJ*YIJ(IN,II)
!     $                 +CSZIJ*ZIJ(IN,II))*DWIJ(IN)

                                !Dissipation term 
                  CI=SQRT(CI2)
                  CJ=SQRT(GAM*A(J)*RHO(J)**(GAM-1.))
                  R1=SQRT((x(i)-x(j))**2+(y(I)-Y(J))**2+(Z(I)-Z(J))**2)
                  DIV(IN)=(X(I)-x(J))*(vx(i)-vx(j))+
     $                 (Y(I)-y(J))*(vy(i)-vy(j))+
     $                 (Z(I)-z(J))*(vz(i)-vz(j))
                  VSIG=CI+CJ-(DIV(IN)/R1)
                  RHOAVG=(RHO(I)+RHO(J))/2
                  VSIG_MU=SQRT((ABS(A(I)*RHO(I)-A(J)*RHO(J)))/RHOAVG)
                  U_I=EN(I)-(VX(I)**2+VY(I)**2+VZ(I)**2)/2
                  U_J=EN(J)-(VX(J)**2+VY(J)**2+VZ(J)**2)/2
                  DIV_I=(X(I)-x(J))*(vx(i))+
     $                 (Y(I)-y(J))*(vy(i))+
     $                 (Z(I)-z(J))*(vz(i))
                  DIV_J=(X(J)-x(I))*(vx(j))+
     $                 (Y(J)-y(I))*(vy(j))+
     $                 (Z(J)-z(I))*(vz(j))
                  
                  E_STAR_I=(0.5*AL*VSIG*((DIV_I)/R1)**2)+
     $                 (AL*VSIG_MU*U_I)
                  E_STAR_J=(0.5*AL*VSIG*((DIV_J)/R1)**2)+
     $                 (AL*VSIG_MU*U_J)
                  CSXIJ_DISS=((E_STAR_I-E_STAR_J)/RHOAVG)*(X(I)-X(J))
     $                 /R1
                  CSYIJ_DISS=((E_STAR_I-E_STAR_J)/RHOAVG)*(Y(I)-Y(J))
     $                 /R1
                 CSZIJ_DISS=((E_STAR_I-E_STAR_J)/RHOAVG)*(Z(I)-Z(J))
     $                 /R1
                  SEIJ(IN)=((CSXIJ-CSXIJ_DISS)*XIJ(IN,II)+
     $                 (CSYIJ-CSYIJ_DISS)*YIJ(IN,II)
     $                 +(CSZIJ-CSZIJ_DISS)*ZIJ(IN,II))*DWIJ(IN)
               end if 
              ENDDO
          ELSE
            DO IN=1,NN(I)
               J=NNI(IN,II)
               CSXIJ=(POR2I*VX(J)+POR2(J)*VX(I))
               CSYIJ=(POR2I*VY(J)+POR2(J)*VY(I))
               CSZIJ=(POR2I*VZ(J)+POR2(J)*VZ(I))
               SEIJ(IN)=(CSXIJ*XIJ(IN,II)+CSYIJ*YIJ(IN,II)
     $              +CSZIJ*ZIJ(IN,II))*DWIJ(IN)
            ENDDO
         ENDIF
         
C     Finally compute each component of dv/dt:
C     "Gather" part of the sum (i-j pair contributes to vdot_i):
         ENDOTI=0
         DO IN=1,NN(I)
            J=NNI(IN,II)
            ENDOTI=ENDOTI-AM(J)*SEIJ(IN)
         ENDDO
         myENDOT(I)=myENDOT(I)+0.5*ENDOTI
         
C     "Scatter" part of the sum (i-j pair contributes to vdot_j):
         DO IN=1,NN(I)
            J=NNI(IN,II)
            myENDOT(J)=myENDOT(J)+0.5*AMI*SEIJ(IN)
         ENDDO
         
      ENDDO

c     Use MPI_ALLREDUCE to collect subtotals, sum them, and redistribute
c     to all processes

      CALL MPI_ALLREDUCE(myendot,endot,n,MPI_REAL,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)
C      CALL MPI_ALLREDUCE(myvydot,vydot,n,MPI_REAL,MPI_SUM,
C     $     MPI_COMM_WORLD,ierr)
C      CALL MPI_ALLREDUCE(myvzdot,vzdot,n,MPI_REAL,MPI_SUM,
C     $     MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(myuijmax,uijmax,n,MPI_REAL,MPI_MAX,
     $     MPI_COMM_WORLD,ierr)

C Calculate self gravity of fluid:
      IF (NGR.GT.0) THEN 
c     doqgrav is the version with radiation reaction on, qgrav with it off
         if(ngravrad.eq.1.and.nrelax.eq.0) then            
            call doqgrav
         else
           CALL QGRAV
         endif
         DO I=1,N
            VXDOT(I)=VXDOT(I)+GX(I)
            VYDOT(I)=VYDOT(I)+GY(I)
            VZDOT(I)=VZDOT(I)+GZ(I)
            if(ngravrad.eq.1) then
               vxdot(i)=vxdot(i)+freacx(i)
               vydot(i)=vydot(i)+freacy(i)
               vzdot(i)=vzdot(i)+freacz(i)
            endif
         ENDDO
      ENDIF

C     If this is a relaxation calculation, add centrifugal and drag forces:
      IF (NRELAX.NE.0) CALL RELAX                                      

      RETURN
      END

************************************************************************                                                                                      
      SUBROUTINE EDOTS
      IMPLICIT NONE
************************************************************************                                                                                              
c     Release 1.0 
c     See Balsara J.Comp.Phys. 121, 357 (1995) and 
c     Lombardi, Rasio, and Shapiro J.Comp.Phys 152, 687 (1999)                                                                                             
c     This routine is used to decide whoch method is used in advance.f 
c     Either adots or endots
************************************************************************                                                                                                           
      INCLUDE 'spha.h'
      INCLUDE 'mpif.h'

      if(nenergy.eq.1) then
         CALL ADOTS
      else if(nenergy.eq.2) then 
         CALL ENDOTS
      endif
      END

***********************************************************************                                                                                  
      SUBROUTINE A_FROM_E
      IMPLICIT NONE
************************************************************************ 
c     Release 1.0
c     This routine calculates a(i) from e(i)
************************************************************************
      INCLUDE 'spha.h'
      INCLUDE 'mpif.h'

      INTEGER I
      
      do i=1,N
         a(i)=(en(i)-((vx(i))**2+(vy(i))**2+(vz(I)**2))/2)/
     $        ((gam-1)*((rho(i))**(gam-1)))
         
         if(a(i).lt.0)write(6,*)'afrome: nonpostive a!',i,a(i),
     $        en(i),vx(i),vy(i),vz(i),rho(i)
         end do 

      END
************************************************************************                                                                                     
      SUBROUTINE E_FROM_A
      IMPLICIT NONE
************************************************************************                                                                                     
c     Release 1.0                                                                                                                                            
c     This routine calculates e(i) from a(i)                                                                                                                
************************************************************************                                                                                    
      INCLUDE 'spha.h'
      INCLUDE 'mpif.h'
      
      INTEGER i
      REAL u1
       
!      OPEN(135,FILE='splash_ghost.txt')
      
      do i=1,N
         en(i)=(gam-1)*a(i)*((rho(i))**(gam-1))+
     $        ((vx(i))**2+(vy(i))**2+(vz(i))**2)/2
         
         u1 = EN(I)-0.5*(vx(I)**2+vy(I)**2+vz(I)**2)
         if (u1.lt.0)then                                                                                                                                  
            write(6,*),"Bad u value!", i,u,EN(I),vx(I),vy(I),vz(I)                                                                                        
         end if 
 !        write(135,*),u
      end do
      END
************************************************************************
