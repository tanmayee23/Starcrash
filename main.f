**************************************************************
c     This file is part of the StarCrash code
c     Version 1.0
c
c     Copyright 2003, by Joshua Faber
c     This code is publicly available under the terms of the
c     GNU General Public License.  See our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.

c This program is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 2 of the License, or
c (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this program; if not, write to the Free Software
c Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
c
***************************************************************

***********************************************************
c     Release 1.0
c     The main body of the standard MPI-based code.
c     calls INIT,MAINIT
*****************************************************************
      IMPLICIT NONE

      INCLUDE 'spha.h'
      INCLUDE 'mpif.h'

      INTEGER ierr,nprocs2,i
      REAL omeg

C Set up execution parameters:

C Initialization:
      CALL MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs2,ierr)
      if(nprocs2.ne.nprocs) then
         write(6,*)'MAJOR ERROR: NPROCS wrong'
         stop
      endif
      initgr=0
      omega2=0

      CALL INIT

      if(myrank.eq.0) then
         OPEN(72,FILE='energy.sph', status='unknown')
         OPEN(73,FILE='biout.sph', status='unknown')
         OPEN(74,FILE='gwdata.sph', status='unknown')
         open(77,file='spin.sph', status='unknown')
         open(78,file='com.sph', status='unknown')
      endif

C Dump initial conditions
C      CALL DUMPINIT

C Main program loop:
      DO WHILE (.TRUE.)
C   Do one iteration:


        CALL MAINIT

C   Check for end of integration:
        IF(T.ge.treloff.and.nrelax.ge.1)then
           if(myrank.eq.0)write(6,*)'MAIN: relaxation off'

           omeg=sqrt(omega2)
           do i=1,n
              vx(i)=-1.0*omeg*y(i)
              vy(i)=omeg*x(i)
              vz(i)=0.
           enddo
           nrelax=0
           trelax=0.0
           omega2=0.0
        endif

        IF (T.GE.TF) THEN 
           if(myrank.eq.0)WRITE (6,*) 'MAIN: end of integration, t=', T
           if(myrank.eq.0) then
              close(72)
              close(73)
              close(74)
              close(77)
              close(78)
           endif
           call MPI_FINALIZE(MPI_COMM_WORLD,ierr)
           STOP 
        ENDIF
      ENDDO

      END
***********************************************************************
      SUBROUTINE MAINIT
***********************************************************
c     Release 1.0
c     One full iteration of the hydro code
c     called by MAIN
c     calls TSTEP,ADVANCE,ADJUST,CHECKPT,OUTPUT
*****************************************************************

      INCLUDE 'spha.h'                                         

C Recalculate time step:
      CALL TSTEP

C Advance particle positions and velocities:
      if(ntimestepper.eq.0)CALL ADVANCE
      if(ntimestepper.eq.1)CALL ADVANCE_EULER
      if(ntimestepper.eq.2)CALL ADVANCE_RK2   
      if(ntimestepper.eq.3)CALL ADVANCE_RK4
      if(ntimestepper.eq.-1)CALL ADVANCE_SHOCK
C Recalculate smoothing lengths:
      CALL ADJUST

      NIT=NIT+1

C Write checkpointing file:
      if(myrank.eq.0)CALL CHECKPT

C Write results:
      if(myrank.eq.0)CALL OUTPUT

      RETURN
      END
************************************************************************
      SUBROUTINE gravquant
      IMPLICIT NONE
***********************************************************
c     Release 1.0
c     Calculate quantities for the parallel code, and check numbers
c     to make sure they work properly
c     called by INIT,LFSTART,SETUP1EM,SETUP1ES,SETUP2CM,SETUP2CS
*****************************************************************

      include 'spha.h'
      include 'mpif.h'

      INTEGER stride,np2,nmaxtest,nngravtest

      if(myrank.eq.0)write(6,*)'CALLED gravquant'
      stride=(n-1)/nprocs+1
      if(stride.gt.NMAXP)then
         write(6,*)'MAJOR ERROR: STRIDE>NMAXP',stride,n,nprocs
         stop
      endif
      n_lower=stride*myrank+1
      n_upper=min(n_lower+stride-1,n)
      np2=nprocs/2
      if (myrank.lt.np2) then
         kstart=nnp*myrank+1
         koffset=1
      else
         kstart=nnp*myrank-nngrav+1
         koffset=0
      endif
      nmaxtest=nmaxp*nprocs
      if(nmaxtest.lt.nmax)then
         write(6,*)'MAJOR ERROR: NMAXP wrong',
     $        nmax,nmaxtest,nmaxp,nprocs
         stop
      endif
      if((n1grav-ngrav).ne.1)then
         write(6,*)'ERROR: NGRAV, N1GRAV'
         stop
      endif
      nngravtest=2*ngrav
      if(nngravtest.ne.nngrav)then
         write(6,*)'ERROR: NNGRAV, NGRAV'
         stop
      endif
      end







