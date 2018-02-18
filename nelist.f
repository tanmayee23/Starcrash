**************************************************************
c     This file is part of the StarCrash code
c     Version 1.0
c
c     Copyright 2003, by Joshua Faber
c     This code is publicly available under the terms of the
c     GNU General Public License.  See our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

      SUBROUTINE NENE
      IMPLICIT NONE
***********************************************************
c     Release 1.0
C     Neighbor searching routine based on the linked-list technique
C     on a single 3-D grid.
c     Called by ADVANCE,INIT,LFSTART,OPTHP2
******************************************************
      INCLUDE 'spha.h'      
      INCLUDE 'mpif.h'
C Select a maximum size for the grid:
      INTEGER ncmx,ncmy,ncmz
C      PARAMETER (NCMX=128,NCMY=128,NCMZ=128)
      PARAMETER (NCMX=512,NCMY=64,NCMZ=64)

C Linked-List and Head-Of-Cell arrays (cf. Hockney and Eastwood):
      INTEGER*4 HOC,LL(NMAX)
      COMMON/HEAD/ HOC(NCMX,NCMY,NCMZ)
      LOGICAL EXC,OUT(NMAX)
      integer ierr,mynn(N)

      REAL xmin,ymin,zmin,xmax,ymax,zmax,hpav,hh,r2max,dx,dy,dz,r2
      INTEGER i,ncx,ncy,ncz,k,l,m,npout,nex,ii,k1,k2,l1,l2,m1,m2
      INTEGER kne,lne,mne,ine

      integer nnemax,nnemin,nnesig,nneavr

C Determine system box and mean HP:
      XMIN=1.E30
      YMIN=1.E30
      ZMIN=1.E30
      XMAX=-1.E30
      YMAX=-1.E30
      ZMAX=-1.E30
      HPAV=0.
      DO I=1,N
        XMIN=MIN(XMIN,X(I))
        YMIN=MIN(YMIN,Y(I))
        ZMIN=MIN(ZMIN,Z(I))
        XMAX=MAX(XMAX,X(I))
        YMAX=MAX(YMAX,Y(I))
        ZMAX=MAX(ZMAX,Z(I))
        HPAV=HPAV+HP(I)
        mynn(I)=0
      ENDDO
      HPAV=HPAV/FLOAT(N)
      HH=HPAV*1.3
C     (This should give near-optimal search times)
C     6/29/98: The factor of 1.3 seems to give search times for typical
C     cases which are about 10-20% shorter than if HH=HPAV.
 
C Calculate number of cells in each dimension:
      NCX=INT((XMAX-XMIN)/HH)+1
      NCY=INT((YMAX-YMIN)/HH)+1
      NCZ=INT((ZMAX-ZMIN)/HH)+1

C Readjust if maximum number allowed is exceeded:
C (Note that this version assumes centering on the origin!)
      EXC=.FALSE.
      IF (NCX.GT.NCMX) THEN
         if(myrank.eq.0)write(6,*)
     $        'Warning! NELIST: EXC TRUE in x-direction',NCX,NCMX
         EXC=.TRUE.
         XMIN = -FLOAT(NCMX/2)*HH
         NCX=NCMX
      ENDIF
      IF (NCY.GT.NCMY) THEN
         if(myrank.eq.0)write(6,*)
     $        'Warning! NELIST: EXC TRUE in y-direction',NCY,NCMY
         EXC=.TRUE.
         YMIN = -FLOAT(NCMY/2)*HH
         NCY=NCMY
      ENDIF
      IF (NCZ.GT.NCMZ) THEN
         if(myrank.eq.0)write(6,*)
     $        'Warning! NELIST: EXC TRUE in z-direction',NCZ,NCMZ
         EXC=.TRUE.
         ZMIN = -FLOAT(NCMZ/2)*HH
         NCZ=NCMZ
      ENDIF

C Initialize Head-Of-Cell array:
      DO M=1,NCZ
        DO L=1,NCY
          DO K=1,NCX
            HOC(K,L,M)=0
          ENDDO
        ENDDO
      ENDDO

C Build linked lists (cf. Hockney and Eastwood):
      NPOUT=0
      IF (EXC) THEN
C    (test for particles outside grid)
         DO I=1,N
            K=INT((X(I)-XMIN)/HH)+1
            L=INT((Y(I)-YMIN)/HH)+1
            M=INT((Z(I)-ZMIN)/HH)+1
            IF ((K.LT.1).OR.(K.GT.NCMX).OR.
     $           (L.LT.1).OR.(L.GT.NCMY).OR.
     $           (M.LT.1).OR.(M.GT.NCMZ)) THEN
               OUT(I)=.TRUE.
               NPOUT=NPOUT+1
            ELSE
               OUT(I)=.FALSE.
               LL(I)=HOC(K,L,M)
               HOC(K,L,M)=I
            ENDIF
         ENDDO
      ELSE
C    (no need to check for out of bounds indices)
         DO I=1,N
            K=INT((X(I)-XMIN)/HH)+1
            L=INT((Y(I)-YMIN)/HH)+1
            M=INT((Z(I)-ZMIN)/HH)+1
            OUT(I)=.FALSE.
            LL(I)=HOC(K,L,M)
            HOC(K,L,M)=I
         ENDDO
      ENDIF
      IF (NPOUT.NE.0.and.myrank.eq.0) WRITE (6,*) 
     $     'NELIST: WARNING !!! NPOUT=',NPOUT

C Calculate neighbor lists
C The GOTO 123 skips to the next particle, either because the
C current particle is outside the grid (then NN(I)=0) or because
C the maximum number of neighbors has been reached.

      NEX=0
      DO I=n_lower,n_upper

         II=I-n_lower+1

        IF (OUT(I)) GOTO 23
        R2MAX=4.*HP(I)**2
C Find indices of cells containing potential neighbors:
        K1=MAX(1,INT((X(I)-2.*HP(I)-XMIN)/HH)+1)
        K2=MIN(NCX,INT((X(I)+2.*HP(I)-XMIN)/HH)+1)
        L1=MAX(1,INT((Y(I)-2.*HP(I)-YMIN)/HH)+1)
        L2=MIN(NCY,INT((Y(I)+2.*HP(I)-YMIN)/HH)+1)
        M1=MAX(1,INT((Z(I)-2.*HP(I)-ZMIN)/HH)+1)
        M2=MIN(NCZ,INT((Z(I)+2.*HP(I)-ZMIN)/HH)+1)
C Look for neighbors in all these cells:
        DO KNE=K1,K2
          DO LNE=L1,L2
            DO MNE=M1,M2
              INE=HOC(KNE,LNE,MNE)
              DO WHILE (INE.NE.0)
                DX=X(I)-X(INE)
                DY=Y(I)-Y(INE)
                DZ=Z(I)-Z(INE)
                R2=DX**2+DY**2+DZ**2
                IF (R2.LT.R2MAX) THEN
C           (a new neighbor has been found)
                  myNN(I)=myNN(I)+1
                  NNI(myNN(I),II)=INE
                  XIJ(myNN(I),II)=DX
                  YIJ(myNN(I),II)=DY
                  ZIJ(myNN(I),II)=DZ
                  IF (myNN(I).GE.NNMAX) THEN
                    NEX=NEX+1
                    GOTO 23
                  ENDIF
                ENDIF
                INE=LL(INE)
              ENDDO
C            (done with that cell)
            ENDDO
          ENDDO
        ENDDO
 23    CONTINUE
C (go to next particle)
      ENDDO

c     Use MPI_ALLREDUCE to sum subtotals and redistribute

      CALL MPI_ALLREDUCE(mynn,nn,n,MPI_INTEGER,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)
      IF (NEX.NE.0.and.myrank.eq.0) WRITE (6,*) 
     $     'NELIST: WARNING !!! NEX=',NEX

      NNEMIN=10000
      NNEMAX=0
      NNEAVR=0
      NNESIG=0
      DO I=1,N
         NNEMIN=MIN(NNEMIN,NN(I))
         NNEMAX=MAX(NNEMAX,NN(I))
         NNEAVR=NNEAVR+NN(I)
         NNESIG=NNESIG+NN(I)**2
      ENDDO
      NNEAVR=INT(FLOAT(NNEAVR)/FLOAT(N))
      NNESIG=INT(SQRT(FLOAT(NNESIG)/FLOAT(N)-FLOAT(NNEAVR)**2))
      if(myrank.eq.0)
     $     write (6,*)'NNMIN:',nnemin,' NNMAX:',nnemax
      if(myrank.eq.0)
     $     write(6,*)' AVG:',nneavr,' SIG:',nnesig


      RETURN
      END
****************************************************************************
      SUBROUTINE ADJUST
      IMPLICIT NONE
***********************************************************
c     Release 1.0
C     Adjusts kernel widths to maintain reasonable neighbor lists:
c     Called MAINIT
******************************************************
C     Adjust kernel widths to maintain reasonable neighbor lists:
      INCLUDE 'spha.h'                                         

      INTEGER i,nhtin,nhexc

C     Try to maintain NN as close as possible to NNOPT (cf. Katz 
C     and Hernquist):
      DO I=1,N
         IF (NN(I).NE.0)HP(I)=HP(I)*0.5*
     $        (1.+(FLOAT(NNOPT)/FLOAT(NN(I)))**(1./3.))
      ENDDO

C     Impose lower limit on h (skip if HMIN<=0):
      IF (HMIN.GT.0.) THEN
         NHTIN=0
         DO I=1,N
            IF (HP(I).LT.HMIN) THEN
               HP(I)=HMIN
               NHTIN=NHTIN+1
            ENDIF
         ENDDO
         IF (NHTIN.GT.0) WRITE (6,*) 
     $        'ADJUST: ', NHTIN,' particles had h < HMIN'
      ENDIF

C     Impose upper limit on h (skip if HMAX<=0):
      IF (HMAX.GT.0.) THEN
         NHEXC=0
         DO I=1,N
            IF (HP(I).GT.HMAX) THEN
               HP(I)=HMAX
               NHEXC=NHEXC+1
            ENDIF
         ENDDO
         IF (NHEXC.GT.0.and.myrank.eq.0) WRITE (6,*) 'ADJUST: ', 
     $        NHEXC,' particles had h > HMAX'
      ENDIF

      RETURN
      END
************************************************************************



