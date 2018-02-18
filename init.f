**************************************************************



c     This file is part of the StarCrash code
c     Version 1.0
c
c     Copyright 2003, by Joshua Faber
c     This code is publicly available under the terms of the
c     GNU General Public License.  See our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

      SUBROUTINE INIT
      IMPLICIT NONE
****************************************************
c     Release 1.0
c     The initialization subroutine
C     Called by MAIN
C     calls TABULINIT,GRAVQUANT,SETUP1EM,SETUP1ES,SETUP2CM,SETUP2CS,NENE
*********************************************************

C     Initialization of run (new or restart)
      INCLUDE 'spha.h'                                          
      LOGICAL RESTART

      REAL hminold,hmaxold,gamold,sep0old,tfold,dtoutold
      REAL alphaold,betaold,eta2old,xgrminold,xgrmaxold
      REAL xgrlimold,ygrlimold,zgrlimold
      REAL ygrminold,ygrmaxold,zgrminold,zgrmaxold,trelaxold
      INTEGER nnoptold,navold,ngrold,nrelaxold,i,nchk
      CHARACTER*3 INAME
      NAMELIST /INITT/ INAME

      
C     Compute look-up tables:
      CALL TABULINIT
      
      
C     Initialize data:
C     Check for existence of dump file:
      INQUIRE (FILE='restart.sph', EXIST=RESTART)
C     If dump file exists, read it and restart from it:

      IF (RESTART) THEN
         
C     Get parameters from input file:
         OPEN(12,FILE='sph.input',ERR=100)
         READ(12,INPUT)
         CLOSE(12)
         
C     WRITE (6,*) 'INIT: continuing run'
C     WRITE (6,*) 'INIT: reading dump file ...'
         OPEN(12,FILE='restart.sph',FORM='UNFORMATTED')
C     (The following READ sequence must match exactly the WRITE sequence
C     used in subroutine DUMP)
         READ(12) N,NNOPTOLD,HMINOLD,HMAXOLD,GAMOLD,SEP0OLD,
     $        TFOLD,DTOUTOLD,NOUT,NLEFT,NIT,T,
     $        NAVOLD,ALPHAOLD,BETAOLD,ETA2OLD,
     $        NGROLD,XGRMINOLD,XGRMAXOLD,YGRMINOLD,YGRMAXOLD,
     $        ZGRMINOLD,ZGRMAXOLD,XGRLIMOLD,YGRLIMOLD,ZGRLIMOLD,
     $        NRELAXOLD,TRELAXOLD,initgr,ngravrad,sol,ntimestepper,
     $        q2xx,q2xy,q2xz,q2yy,q2yz,q2zz
         DO I=1,N
            READ (12) X(I),Y(I),Z(I),AM(I),HP(I),RHO(I),VX(I),VY(I),
     $           VZ(I),VXDOT(I),VYDOT(I),VZDOT(I),A(I),ADOT(I),
     $           GX(I),GY(I),GZ(I),GRPOT(I),ux(i),uy(i),uz(i)
         ENDDO
         READ (12) NCHK,nleft
         CLOSE (12)
         IF (NCHK.NE.N) STOP 'INIT: PROBLEM WITH DUMP FILE ???'

         call gravquant
         if(myrank.eq.0)WRITE (6,*) 'INIT:            ... done'
C     Otherwise create new data set:
      ELSE
C     Get type of initial condition from init file (must be
C     a 3-letter code):
         OPEN(12,FILE='sph.init',ERR=100)
         READ(12,INITT)
         CLOSE(12)
         IF (INAME.EQ.'1em') THEN
            CALL setup1em
         ELSE IF (INAME.EQ.'1es') THEN
            CALL setup1es
         ELSE IF (INAME.EQ.'2cs') THEN
            CALL setup2cs
         ELSE IF (INAME.EQ.'2cm') THEN
            CALL setup2cm
         else if (iname.eq.'2cr') THEN
            CALL setup2cr
         ELSE IF (INAME.EQ.'2qs') THEN
            CALL setup2qs
         ELSE IF (INAME.EQ.'2qm') THEN
            CALL setup2qm
         else if (iname.eq.'2qr') THEN
            CALL setup2qr
         else if (iname.eq.'2i1') THEN
            CALL setup2i1
         else if (iname.eq.'2iq') THEN
            CALL setup2iq
         else if (iname.eq.'hy1') THEN
            CALL setuphyp1
         else if (iname.eq.'hyq') THEN
            CALL setuphypq
         else if (iname.eq.'sti') THEN
            CALL shock_tube 
         ELSE
            STOP 'INIT: UNKNOWN INAME ???'
         ENDIF
C     Initialize output parameters:
         NOUT=0
         NIT=0
         T=0.
      ENDIF
C     Initialize neighbor lists:
      CALL NENE

C     Write run parameters
      if(myrank.eq.0)WRITE (6,*) 'INIT: T=',T,' NIT=',NIT
      if(myrank.eq.0)WRITE (6,101) N,NNOPT,HMIN,HMAX,DTOUT,NOUT,TF,
     $     GAM, SEP0,NAV,ALPHA,BETA,ETA2,
     $     NGR,XGRMIN,XGRMAX,YGRMIN,YGRMAX,ZGRMIN,ZGRMAX,
     $     XGRLIM,YGRLIM,ZGRLIM,NRELAX,TRELAX,NGRAVRAD,SOL,NTIMESTEPPER
 101  FORMAT (' INIT: PARAMETERS FOR THIS RUN:',/,
     $     ' N=',I7,' NNOPT=',I4,' HMIN=',E12.4,' HMAX=',E12.4,/,
     $     ' DTOUT=',F8.4,' NOUT=',I4,' TF=',F9.3,/,
     $     ' GAM=',F8.4,' SEP0=',E12.4,/,
     $     ' NAV=',I2,' ALPHA=',F6.2,' BETA=',F6.2,' ETA2=',F7.4,/,
     $     ' NGR=',I3,' XGRMIN=',E12.4,' XGRMAX=',E12.4,/,
     $     '        YGRMIN=',E12.4,' YGRMAX=',E12.4,/,
     $     '        ZGRMIN=',E12.4,' ZGRMAX=',E12.4,/,
     $     ' XGRLIM=',E12.4,' YGRLIM=',E12.4,' ZGRLIM=',E12.4,/,
     $     ' NRELAX=',I2,' TRELAX=',E12.4,' NGRAVRAD=',I2,' SOL=',E12.4, 
     $     'NTIMESTEPPER=',I2,
     $       /)       


      
      RETURN
C     Error condition:
 100  STOP 'INIT:  ERROR READING INPUT FILE ???'
      END
************************************************************************
      SUBROUTINE LFSTART
****************************************************
c     Release 1.0
c     Set up the leapfrog timestep algorithm by advancing velocities
c     half a timestep
C     Called by SETUP1ES,SETUP1EM,SETUP2CS,SETUP2CM
C     CALLS GRAVQUANT,NENE,RHOS,VDOTS,TSTEP
*********************************************************
      IMPLICIT NONE
C    (Prepare for first leap-frog iteration)
      INCLUDE 'spha.h'                                          

      REAL dth
      INTEGER i
       if(myrank.eq.0)write (6,*),ntimestepper

      !ntimestepper=1
 2    if(myrank.eq.0)write (6,*) 'LFSTART: starting'
      
C Set ADOT=0 initially (assumes no shock in initial data!):
      DO I=1,N
        ADOT(I)=0.
      ENDDO
C     This call to gravquant is redundant, but safe
      call gravquant

C Advance velocities to half-timestep:
      CALL NENE
      CALL RHOS
      CALL VDOTS
      CALL TSTEP
         DTH=0.5*DT
         if(myrank.eq.0) write(6,*)'Initial dth=',dth
         if(ntimestepper.eq.0)then
         DO I=1,N
            VX(I)=VX(I)+VXDOT(I)*DTH
            VY(I)=VY(I)+VYDOT(I)*DTH
            VZ(I)=VZ(I)+VZDOT(I)*DTH
         ENDDO
         end if
         if(myrank.eq.0) write(6,*) 'LFSTART: closing'
      RETURN
      END














