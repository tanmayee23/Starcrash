**************************************************************
c     This file is part of the StarCrash code
c     Version 1.0
c
c     Copyright 2003, by Joshua Faber
c     This code is publicly available under the terms of the
c     GNU General Public License.  See our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

      SUBROUTINE ADOTS
      IMPLICIT NONE
********************************************
c     Release 1.0
c     Selects which AV routine to use to calculate dA/dt
c     Called by ADVANCE
c     Calls BALADOTS,CLAADOTS,NEWADOTS
******************************************
      INCLUDE 'spha.h'
C     write (6,*) 'ADOTS: starting'
      IF(nrelax.eq.0) then
         IF (NAV.EQ.0.or.NAV.EQ.1) THEN
            CALL BALADOTS
         ELSE IF (NAV.EQ.2) THEN
            CALL CLAADOTS
         ELSE IF (NAV.EQ.3) THEN
            CALL NEWADOTS
         ELSE
            write(6,*)'ERROR! ILLEGAL VALUE OF NAV=',nav
            stop
         ENDIF
      ENDIF
      RETURN
      END
************************************************************************
      SUBROUTINE VDOTS
      IMPLICIT NONE
***********************************************************
c     Release 1.0
c     Selects which AV routine to use to calculate dV/dt
c     Called by ADVANCE
c     Calls BALVDOTS,CLAVDOTS,NEWVDOTS
***************************************************************

      INCLUDE 'spha.h'
C     write (6,*) 'VDOTS: starting'
      IF (NAV.EQ.0 .OR. NAV.EQ.1) THEN
         CALL BALVDOTS
      ELSE IF (NAV.EQ.2) THEN
         CALL CLAVDOTS
      ELSE IF (NAV.EQ.3) THEN
         CALL NEWVDOTS
      ELSE
         write(6,*)'ERROR! ILLEGAL VALUE OF NAV=',nav
         stop
      ENDIF
      RETURN
      END
************************************************************************

