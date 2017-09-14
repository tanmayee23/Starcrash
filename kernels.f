**************************************************************
c     This file is part of the StarCrash code
c     Version 1.0
c
c     Copyright 2003, by Joshua Faber
c     This code is publicly available under the terms of the
c     GNU General Public License.  See our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

      SUBROUTINE TABULINIT
*****************************************************************
c     Release 1.0
c     Calculate tabulated values of smoothing kernel for SPH summation
c     Called by INIT
c     Calls functions W,DW
**********************************************************************
      INCLUDE 'spha.h'                                         

C Compute tabulations of W(u) and dW(u)/du:
      DO I=1,NTAB
        UU2=4.*FLOAT(I-1)/FLOAT(NTAB-1)+1.E-15
        WTAB(I)=W(SQRT(UU2))
        DWTAB(I)=DW(SQRT(UU2))
      ENDDO

C Compute normalization constant:
      CTAB=FLOAT(NTAB-1)/4.

      RETURN
      END
***********************************************************************
      FUNCTION W(U)
c     Release 1.0
c     Calculate smoothing kernel function W 
c     See Rasio and Shapiro, ApJ 401, 226 (1992), Eq. 4
c     Called by TABULINIT

      IF (U.LT.1.) THEN
       W=1.-1.5*U**2+0.75*U**3
      ELSE IF (U.LT.2.) THEN
       W=0.25*(2.-U)**3
      ELSE
       W=0.
      ENDIF
      W=W/3.14159265359

      RETURN
      END
************************************************************************
      FUNCTION DW(U)
c     Release 1.0
c     Calculate derivative of the smoothing kernel function W 
c     Called by TABULINIT

      IF (U.LT.1.) THEN
       DW=-3.*U+2.25*U**2
      ELSE IF (U.LT.2.) THEN
       DW=-0.75*(2.-U)**2
      ELSE
       DW=0.
      ENDIF
      DW=DW/(U*3.14159265359)

      RETURN
      END
************************************************************************
