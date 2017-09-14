**************************************************************
c     This file is part of the StarCrash code
c     Version 1.0
c
c     Copyright 2003, by Joshua Faber
c     This code is publicly available under the terms of the
c     GNU General Public License.  See our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

      program bin2pos
**************************************************************
c     Use this program to convert a binary output file into
c     an ASCII input file named pos.sph
***********************************************************
      CHARACTER*12 filename
      LOGICAL CONVERT
      real dum
      integer ndum

      write(6,*)'1 for restart.sph, 2 for out***.sph'
      READ *,i1
      If(i1.eq.1) then
         filename='restart.sph'
      else if (i1.eq.2) then
         write(6,*)'enter number as integer for outfile'
         read *,i2
         write(filename,101) i2
 101     FORMAT('out',I3.3,'.sph')
      endif
      INQUIRE (FILE=filename, EXIST=CONVERT)
      IF (.not.convert) then
         write(6,*)filename,' does not exist'
         stop
      endif
      OPEN(12,FILE=filename,FORM='UNFORMATTED')
      OPEN(14,FILE='pos.sph', status='unknown')
C    (The following READ sequence must match exactly the WRITE sequence
C     used in subroutine DUMP)
      READ(12) N,nnopt,hmin,hmax,gam,dum,
     $           dum,dum,ndum,ndum,ndum,dum,Ndum,dum,dum,dum,
     $           Ndum,dum,dum,dum,dum,dum,dum,
     $           Ndum,TRELAX
      write(14,'(2i6,4e15.7)')n,nnopt,hmin,hmax,gam,trelax
      DO I=1,N
         READ (12) X,Y,z,am,hp,rho,dum,dum,dum,dum,dum,dum,a
         write(14,'(7e15.7)') X,Y,z,am,hp,rho,a
      ENDDO
      READ (12) NCHK
      CLOSE (12)
      IF (NCHK.NE.N) STOP 'INIT: PROBLEM WITH DUMP FILE ???'
      END



























