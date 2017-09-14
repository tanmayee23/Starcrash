**************************************************************
c     This file is part of the StarCrash code
c     Version 1.0
c
c     Copyright 2003, by Joshua Faber
c     This code is publicly available under the terms of the
c     GNU General Public License.  See our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

      SUBROUTINE CLAADOTS
      IMPLICIT NONE
*****************************************************
c     Release 1.0
c     This routine calculates ADOT using the Monaghan AV formulation.
c     See Monaghan J.Comp.Phys. 82, 1 (1989) and
c     Lombardi, Rasio, and Shapiro J.Comp.Phys 152, 687 (1999)
c     Called by ADOTS
************************************************** 

C Evaluate right-hand sides of dA/dt equations:
      INCLUDE 'spha.h'                                          
      INCLUDE 'mpif.h'
      REAL DWIJ(NNMAX),PIJ(NNMAX),SIJ(NNMAX),DIV(NNMAX),myadot(N)

      REAL hpi,h2,h5,ami,ci2,rhoi,por2i,vxi,vyi,vzi
      REAL hij,cij,uij,rij,adoti,r2
      INTEGER i,ii,in,itab,j,ierr

C  Initialize dA/dt before accumulating:
      DO I=1,N
        myadot(i)=0.
      ENDDO

C  dA/dt=0 if GAM=1. (isothermal calculation):
      IF (GAM.EQ.1.) RETURN

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

C Calculate gradWij's to all neighbors:
         DO IN=1,NN(I)
            R2=XIJ(IN,II)**2.+YIJ(IN,II)**2.+ZIJ(IN,II)**2.
            ITAB=INT(CTAB*R2/H2)+1
            DWIJ(IN)=DWTAB(ITAB)/H5
         ENDDO
         
C     Calculate (Xi-Xj)*(VXi-VXj) to all neighbors:
         DO IN=1,NN(I)
            J=NNI(IN,II)
            DIV(IN)=XIJ(IN,II)*(VXI-VX(J))+
     $           YIJ(IN,II)*(VYI-VY(J))+
     $           ZIJ(IN,II)*(VZI-VZ(J))
         ENDDO
         
         IF (NAV.NE.0) THEN
C     Calculate artificial viscosity contribution to pressure gradients:
            DO IN=1,NN(I)
               IF (DIV(IN).LT.0.) THEN
                  J=NNI(IN,II)
                  
                  R2=XIJ(IN,II)**2.+YIJ(IN,II)**2.+ZIJ(IN,II)**2.
                  HIJ=0.5*(HPI+HP(J))
                  UIJ=DIV(IN)/(HIJ*(R2/HIJ**2+ETA2))
                  CIJ=(0.5*(CI2+GAM*A(J)*RHO(J)**(GAM-1.)))**0.5
                  RIJ=0.5*(RHOI+RHO(J))
                  PIJ(IN)=(-ALPHA*UIJ*CIJ+BETA*UIJ**2)/RIJ
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
      END DO

c     Use MPI_ALLREDUCE to collect subtotals, sum them, and redistribute
c     to all processes
      CALL MPI_ALLREDUCE(myadot,adot,n,MPI_REAL,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)

      RETURN
      END
************************************************************************
      SUBROUTINE CLAVDOTS
      IMPLICIT NONE
*****************************************************
c     Release 1.0
c     This routine calculates VDOT using the Monaghan AV formulation.
c     See Monaghan J.Comp.Phys. 82, 1 (1989) and
c     Lombardi, Rasio, and Shapiro J.Comp.Phys 152, 687 (1999)
c     Called by VDOTS
c     Calls DOQGRAV,QGRAV
************************************************** 
C Evaluate right-hand sides of equations of motion:
      INCLUDE 'spha.h'
      INCLUDE 'mpif.h'

      REAL UIJMAX(NMAX)
      COMMON/UIJMAX/ UIJMAX
      REAL DWIJ(NNMAX),PIJ(NNMAX),SXIJ(NNMAX),SYIJ(NNMAX),SZIJ(NNMAX),
     $       DIV(NNMAX),myuijmax(N),myvxdot(N),myvydot(N),
     $       myvzdot(N)

      REAL hpi,h2,h5,ami,ci2,rhoi,por2i,vxi,vyi,vzi
      REAL hij,cij,uij,rij,r2,csij,vxdoti,vydoti,vzdoti
      INTEGER i,ii,in,itab,j,ierr

C  Initialize dv/dt before accumulating:
      DO I=1,N
        myVXDOT(I)=0.
        myVYDOT(I)=0.
        myVZDOT(I)=0.
        myuijmax(I)=0.
      ENDDO

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

C Calculate gradWij's to all neighbors:
        DO IN=1,NN(I)
          R2=XIJ(IN,II)**2.+YIJ(IN,II)**2.+ZIJ(IN,II)**2.
          ITAB=INT(CTAB*R2/H2)+1
          DWIJ(IN)=DWTAB(ITAB)/H5
        ENDDO

        IF (NAV.NE.0.and.NRELAX.eq.0) THEN
C    Calculate artificial viscosity contribution to pressure gradients:
          DO IN=1,NN(I)
            J=NNI(IN,II)
            DIV(IN)=XIJ(IN,II)*(VXI-VX(J))+
     $              YIJ(IN,II)*(VYI-VY(J))+
     $              ZIJ(IN,II)*(VZI-VZ(J))
          ENDDO
          UIJMAX(I)=0.
          DO IN=1,NN(I)
            IF (DIV(IN).LT.0.) THEN
              J=NNI(IN,II)

              R2=XIJ(IN,II)**2.+YIJ(IN,II)**2.+ZIJ(IN,II)**2.
              HIJ=0.5*(HPI+HP(J))
              UIJ=DIV(IN)/(HIJ*(R2/HIJ**2+ETA2))
              myUIJMAX(I)=MAX(myUIJMAX(I),ABS(UIJ))
              CIJ=(0.5*(CI2+GAM*A(J)*RHO(J)**(GAM-1.)))**0.5
              RIJ=0.5*(RHOI+RHO(J))
              PIJ(IN)=(-ALPHA*UIJ*CIJ+BETA*UIJ**2)/RIJ
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

C Finally compute each component of dv/dt:
C   "Gather" part of the sum (i-j pair contributes to vdot_i):
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

C    "Scatter" part of the sum (i-j pair contributes to vdot_j):
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
c     call doqgrav is radiation reaction is on, qgrav if it's off
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

C If this is a relaxation calculation, add centrifugal and drag forces:
      IF (NRELAX.NE.0) CALL RELAX                                      

      RETURN
      END



