**************************************************************
c     This file is part of the StarCrash code
c     Version 1.0
c
c     Copyright 2003, by Joshua Faber
c     This code is publicly available under the terms of the
c     GNU General Public License.  See our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

      SUBROUTINE GETCURLV
      IMPLICIT NONE
*****************************************************
c     Release 1.0
c     This routine calculates the curl of the velocity field
c     See Lombardi, Rasio, and Shapiro J.Comp.Phys 152, 687 (1999) Eq. 19
c     Called by BALVDOTS
************************************************** 

C Calculate curl of velocity
      INCLUDE 'spha.h'                                          
      INCLUDE 'mpif.h'
      REAL CURLVX(NMAX),CURLVY(NMAX),CURLVZ(NMAX) 
      COMMON/COMCURLV/CURLVX,CURLVY,CURLVZ
      REAL DWIJ(NNMAX),CURLX(NNMAX),CURLY(NNMAX),CURLZ(NNMAX)
      REAL mycurlvx(n),mycurlvy(n),mycurlvz(n)

      REAL hpi,h2,h5,ami,ci2,rhoi,vxi,vyi,vzi,r2
      real curlvxi,curlvyi,curlvzi
      INTEGER i,ii,in,itab,ierr,j

C  Initialize before accumulating:
      DO I=1,N
        myCURLVX(I)=0.
        myCURLVY(I)=0.
        myCURLVZ(I)=0.
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
C     Calculate gradWij's to all neighbors:
         DO IN=1,NN(I)
            R2=XIJ(IN,II)**2.+YIJ(IN,II)**2.+ZIJ(IN,II)**2.
            ITAB=INT(CTAB*R2/H2)+1
            DWIJ(IN)=DWTAB(ITAB)/H5
         ENDDO
C     Calculate (Xi-Xj)*(VXi-VXj) to all neighbors:
         DO IN=1,NN(I)
            J=NNI(IN,II)
            CURLX(IN)=ZIJ(IN,II)*(vy(i)-vy(j))
     $           -YIJ(IN,II)*(vz(i)-vz(j))
            CURLY(IN)=XIJ(IN,II)*(vz(i)-vz(j))
     $           -ZIJ(IN,II)*(vx(i)-vx(j))
            CURLZ(IN)=YIJ(IN,II)*(vx(i)-vx(j))
     $           -XIJ(IN,II)*(vy(i)-vy(j))
         ENDDO
C     Compute the 3 components of the curl of the velocity:
C     CURLVX(I), CURLVY(I) and CURLVZ(I)
C     "Gather" part of the sum (i-j pair contributes to curl v_i):
         CURLVXI=0.
         CURLVYI=0.
         CURLVZI=0.
         DO IN=1,NN(I)
            J=NNI(IN,II)
            CURLVXI=CURLVXI+AM(J)*DWIJ(IN)*CURLX(IN)
            CURLVYI=CURLVYI+AM(J)*DWIJ(IN)*CURLY(IN)
            CURLVZI=CURLVZI+AM(J)*DWIJ(IN)*CURLZ(IN)
         ENDDO
         myCURLVX(I)=myCURLVX(I)+0.5*CURLVXI/RHOI
         myCURLVY(I)=myCURLVY(I)+0.5*CURLVYI/RHOI
         myCURLVZ(I)=myCURLVZ(I)+0.5*CURLVZI/RHOI
C    "Scatter" part of the sum (i-j pair contributes to curl v_j):
         DO IN=1,NN(I)
            J=NNI(IN,II)
            myCURLVX(J)=myCURLVX(J)+0.5*AMI*DWIJ(IN)*CURLX(IN)/RHO(J)
            myCURLVY(J)=myCURLVY(J)+0.5*AMI*DWIJ(IN)*CURLY(IN)/RHO(J)
            myCURLVZ(J)=myCURLVZ(J)+0.5*AMI*DWIJ(IN)*CURLZ(IN)/RHO(J)
         ENDDO
      ENDDO

c     use MPI_ALLREDUCE to sum subtotals and distribute to processes

      CALL MPI_ALLREDUCE(mycurlvx,curlvx,n,MPI_REAL,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(mycurlvy,curlvy,n,MPI_REAL,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(mycurlvz,curlvz,n,MPI_REAL,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)

      RETURN
      END








