**************************************************************
c     This file is part of the StarCrash code
c     Version 1.0
c
C     Random Number generator, returns x s.t. 0<x<1
c
c     derived from the GNU Scientific Library (GSL)
c     Copyright (C) 1996, 1997, 1998, 1999, 2000 James Theiler, Brian Gough
c
c     Called by SETUP1ES,SETUP1EM,SETUP2CS,SETUP2CM 
c
c     This code is publicly available under the terms of the
c     GNU General Public License.  See our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

      FUNCTION ran1(seed)

      integer seed,m,a,q,r,n_shuffle,n_div

      PARAMETER(m=2147483647, a = 16807, q = 127773, r = 2836)
      PARAMETER(n_shuffle=32,n_div=1 + 2147483646/N_SHUFFLE)

      REAL ran1,x_max,am
      
      PARAMETER(am=1.0/m,x_max=1.0-1.2e-7)


      INTEGER j,h,shuffle(n_shuffle),y
      SAVE shuffle,y
      DATA shuffle /N_SHUFFLE*0/, y /0/

      if (seed.le.0.or.y.eq.0) then
        seed=max(-seed,1)
        do 11 j=N_SHUFFLE+8,1,-1
          h=seed/q
          seed=a*(seed-h*q)-h*r
          if (seed.lt.0) seed=seed+m
          if (j.le.N_SHUFFLE) shuffle(j)=seed
11      continue
        y=shuffle(1)
      endif
      h=seed/q
      seed=a*(seed-h*q)-h*r
      if (seed.lt.0) seed=seed+m
      j=1+y/N_DIV
      y=shuffle(j)
      shuffle(j)=seed
      ran1=min(am*y,x_max)
      return
      END
