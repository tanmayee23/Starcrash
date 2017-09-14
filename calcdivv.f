**************************************************************
c     This file is part of the StarCrash code
c     Version 1.0
c
c     Copyright 2003, by Joshua Faber
c     This code is publicly available under the terms of the
c     GNU General Public License.  See our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

      SUBROUTINE CALCDIVV
      IMPLICIT NONE
*****************************************************
c     Release 1.0
c     This routine calculates the divergence of the velocity field
c     See Lombardi, Rasio, and Shapiro J.Comp.Phys 152, 687 (1999) Eq. 15
c     Called by BALADOTS,BALVDOTS,NEWADOTS,NEWVDOTS
************************************************** 

      INCLUDE 'spha.h'                                          
      include 'mpif.h'
      REAL DIVV(NMAX)
      COMMON/COMMDIVV/DIVV
      REAL DWIJ(NNMAX),TIJ(NNMAX),DIV(NNMAX),mydivv(n)

      REAL hpi,h2,h5,ami,ci2,rhoi,vxi,vyi,vzi,r2,divvi
      INTEGER i,ii,in,itab,j,ierr

C  Initialize before accumulating:
      DO I=1,N
        myDIVV(I)=0.
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
        VXI=VX(I)
        VYI=VY(I)
        VZI=VZ(I)

C Calculate gradWij's to all neighbors:
        DO IN=1,NN(I)
          R2=XIJ(IN,II)**2.+YIJ(IN,II)**2.+ZIJ(IN,II)**2.
          ITAB=INT(CTAB*R2/H2)+1
          DWIJ(IN)=DWTAB(ITAB)/H5
        ENDDO

C Calculate (Xi-Xj)*(VXi-VXj) to all neighbors:
        DO IN=1,NN(I)
           J=NNI(IN,II)

           DIV(IN)=XIJ(IN,II)*(vx(i)-vx(j))+
     $          YIJ(IN,II)*(vy(i)-vy(j))+
     $          ZIJ(IN,II)*(vz(i)-vz(j))
           TIJ(IN)=DWIJ(IN)*DIV(IN)
        ENDDO

C Compute divergence of velocity, DIVV(I)
C   "Gather" part of the sum (i-j pair contrinutes to div v_i):
         DIVVI=0.
         DO IN=1,NN(I)
            J=NNI(IN,II)

            DIVVI=DIVVI-AM(J)*TIJ(IN)/RHOI
         ENDDO
         myDIVV(I)=myDIVV(I)+0.5*DIVVI

C    "Scatter" part of the sum (i-j pair contributes to div v_j):
         DO IN=1,NN(I)
            J=NNI(IN,II)

            myDIVV(J)=myDIVV(J)-0.5*AMI*TIJ(IN)/RHO(J)
         ENDDO
      ENDDO

c     Use MPI_ALLREDUCE to sum over subtotals and distribute back out

      CALL MPI_ALLREDUCE(mydivv,divv,n,MPI_REAL,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)

      RETURN
      END
************************************************************************







