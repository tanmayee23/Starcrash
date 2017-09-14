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
            DWIJ(IN)=DWTAB(ITAB)/H5
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
     $       DIV(NNMAX)
      REAL myuijmax(N),myvxdot(N),myvydot(N),myvzdot(N)

      REAL hpi,h2,h5,ami,ci2,rhoi,por2i,vxi,vyi,vzi,divvi
      REAL abscurli,r2,hij,ci,cj,fi,fj,udbij,csij
      real vxdoti,vydoti,vzdoti
      INTEGER i,ii,in,itab,j,ierr

      divvi=0.
      abscurli=0.

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
            DWIJ(IN)=DWTAB(ITAB)/H5
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
               CSIJ=(POR2I+POR2(J)+PIJ(IN))*DWIJ(IN)
               SXIJ(IN)=CSIJ*XIJ(IN,II)
               SYIJ(IN)=CSIJ*YIJ(IN,II)
               SZIJ(IN)=CSIJ*ZIJ(IN,II)
            ENDDO
         ELSE
            DO IN=1,NN(I)
               J=NNI(IN,II)
               CSIJ=(POR2I+POR2(J))*DWIJ(IN)
               SXIJ(IN)=CSIJ*XIJ(IN,II)
               SYIJ(IN)=CSIJ*YIJ(IN,II)
               SZIJ(IN)=CSIJ*ZIJ(IN,II)
            ENDDO
         ENDIF
         
C     Finally compute each component of dv/dt:
C     "Gather" part of the sum (i-j pair contributes to vdot_i):
         VXDOTI=0.
         VYDOTI=0.
         VZDOTI=0.
         DO IN=1,NN(I)
            J=NNI(IN,II)
            VXDOTI=VXDOTI-AM(J)*SXIJ(IN)
            VYDOTI=VYDOTI-AM(J)*SYIJ(IN)
            VZDOTI=VZDOTI-AM(J)*SZIJ(IN)
         ENDDO
         myVXDOT(I)=myVXDOT(I)+0.5*VXDOTI
         myVYDOT(I)=myVYDOT(I)+0.5*VYDOTI
         myVZDOT(I)=myVZDOT(I)+0.5*VZDOTI
         
C     "Scatter" part of the sum (i-j pair contributes to vdot_j):
         DO IN=1,NN(I)
            J=NNI(IN,II)
            myVXDOT(J)=myVXDOT(J)+0.5*AMI*SXIJ(IN)
            myVYDOT(J)=myVYDOT(J)+0.5*AMI*SYIJ(IN)
            myVZDOT(J)=myVZDOT(J)+0.5*AMI*SZIJ(IN)
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





