**************************************************************
c     This file is part of the StarCrash code
c     Version 1.0
c
c     Copyright 2003, by Joshua Faber
c     This code is publicly available under the terms of the
c     GNU General Public License.  See our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

      program main
      LOGICAL RESTART,CONVERT,convert2
      character*3 inittype

      sep0=0.
      qdar=0.0

      open(12,file='sph.input',form='formatted')
      open(13,file='sph.init',form='formatted')
      write(12,*)'&INPUT '
      write(13,*)'&INITT '

      write(6,*)'Welclome to the SPH code!  This routine will set up'
      write(6,*)'some input files for you, and tell you what others'
      write(6,*)'need to be supplied for everything to work properly.'

      write(6,*)'First, are you starting a new calculation,' 
      write(6,*)'or continuing'
      write(6,*)'from a restart file from a previous run?'
      write(6,*)'Note that duplicating a relaxed single-star model'
      write(6,*)'to use as the initial state of a binary is'
      write(6,*)'considered a new run for our purposes.'
      write(6,*)'Please type 1 for a new run, 2 to continue an old one'
      read *,nstart
      
      if(nstart.eq.1) then
         write(6,*)'How many NS in the simulation? Please type 1 or 2.'
         read *,nstar
         if(nstar.eq.1) then
            write(6,*)'Do you want equal-mass SPH particles,' 
            write(6,*)'randomly located, or do you want equally spaced'
            write(6,*)'particles with varying mass?'
            write(6,*)'Type 1 for equal-mass, 2 for equal-spacing'
            read *,nequal
            if(nequal.eq.1) then
               inittype='1em'
               call do1em
                !write(*,*), nequal
            else if (nequal.eq.2) then
               inittype='1es'
               call do1es
            else
               write(6,*)'illegal parameter'
               stop
            endif
         else if(nstar.eq.2)then
            write(6,*)'Do you want equal-mass or unequal-mass NS?'
            write(6,*)'Type 1 for equal mass NS, type 2 for unequal.'
            read *,nq
            if(nq.ne.1.and.nq.ne.2) then
               write(6,*)'Illegal Parameter!'
               stop
            endif
            write(6,*)'We have two routines for setting up binaries'
            write(6,*)'in circular orbits:'
            write(6,*)'Type 1 if the stars are initially irrotational'
            write(6,*)'Type 2 if the stars are init. synchronized'
            write(6,*)'We also have a routine that sets up a'
            write(6,*)'hyperbolic collision.  Type 3 for this routine.'
            read *,nspin
            if(nspin.eq.2) then
               if(nq.eq.1) then
                  write(6,*)'Have you calculated a relaxed'
                  write(6,*)'single-star model? Type 1 for yes,'
                  write(6,*)'2 for no.'
                  read *,n1s
                  if(n1s.eq.1) then
                     INQUIRE (FILE='pos.sph', EXIST=CONVERT)
                     if(.not.convert) then
                        write(6,*)'Please run the routine b2pos on the'
                        write(6,*)'output file containing the single'
                        write(6,*)'configuration.  It will produce a'
                        write(6,*)'file called pos.sph.  Place it in'
                        write(6,*)'this directory and rerun testinput'
                        stop
                     else
                        inittype='2cr'
                        call do2cr
                     endif
                  else if (n1s.eq.2) then
                     write(6,*)
     $            'Do you want equal-mass SPH particles,' 
                     write(6,*)
     $             'randomly located, or do you want equally spaced'
                     write(6,*)
     $             'particles with varying mass?'
                     write(6,*)
     $             'Type 1 for equal-mass, 2 for equal-spacing'
                     read *,nequal
                     if(nequal.eq.1) then
                        inittype='2cm'
                        call do2cm
                     else if (nequal.eq.2) then
                        inittype='2cs'
                        call do2cs
                     else
                        write(6,*)'illegal parameter'
                        stop
                     endif
                  else
                     write(6,*)'illegal parameter'
                     stop
                  endif
               else
                  write(6,*)'Have you calculated relaxed'
                  write(6,*)'single-star models? Type 1 for yes,'
                  write(6,*)'2 for no.'
                  read *,n2s
                  if(n2s.eq.1) then
                     INQUIRE (FILE='pos.sph', EXIST=CONVERT)
                     INQUIRE (FILE='pos2.sph', EXIST=CONVERT2)
                     if(.not.convert) then
                        write(6,*)'Please run the routine'
                        write(6,*)'b2pos on the output file containing'
                        write(6,*)'the configuration of the first star.'
                        write(6,*)'It will produce a file called'
                        write(6,*)'pos.sph. Place it in this directory,'
                        write(6,*)'and the rerun testinput.'
                     endif
                     if(.not.convert2) then
                        if(convert) then
                           write(6,*)'Please run the routine'
                        else
                           write(6,*)'Please also run the routine'
                        endif
                        write(6,*)'b2pos2 on the output file containing'
                        write(6,*)'the configuration of the 2nd star.'
                        write(6,*)'It will produce a file called'
                        write(6,*)'pos2.sph. Put it in this directory,'
                        write(6,*)'and then rerun testinput.'
                     endif
                     if(convert.and.convert2) then
                        inittype='2qr'
                        call do2qr
                     else
                        stop
                     endif
                  else if (n2s.eq.2) then
                     write(6,*)
     $           'Do you want equal-mass SPH particles,' 
                     write(6,*)
     $           'randomly located, or do you want equally spaced'
                     write(6,*)
     $           'particles with varying mass?'
                     write(6,*)
     $           'Type 1 for equal-mass, 2 for equal-spacing'
                     read *,nequal
                     if(nequal.eq.1) then
                        inittype='2qm'
                        call do2qm
                     else if (nequal.eq.2) then
                        inittype='2qs'
                        call do2qs
                     else
                        write(6,*)'illegal parameter'
                        stop
                     endif
                  else
                     write(6,*)'illegal parameter'
                     stop
                  endif
               endif
            else if (nspin.eq.1) then
               if(nq.eq.1) then
                  write(6,*)'Have you calculated a relaxed'
                  write(6,*)'single-star model? Type 1 for yes,'
                  write(6,*)'2 for no.'
                  read *,n1s
                  if(n1s.eq.1) then
                     INQUIRE (FILE='pos.sph', EXIST=CONVERT)
                     if(.not.convert) then
                        write(6,*)'Please run the routine b2pos on the'
                        write(6,*)'output file containing the single'
                        write(6,*)'configuration.  It will produce a'
                        write(6,*)'file called pos.sph.  Place it in'
                        write(6,*)'this directory and rerun testinput'
                        stop
                     else
                        inittype='2i1'
                        call do2i1
                     endif
                  else
                     write(6,*)'You will need to calculate such a model'
                     write(6,*)'first.  Please rerun testinput to'
                     write(6,*)'set up such a run.'
                     stop
                  endif
               else
                  write(6,*)'Have you calculated relaxed'
                  write(6,*)'single-star models? Type 1 for yes,'
                  write(6,*)'2 for no.'
                  read *,n1s
                  if(n1s.eq.1) then
                     INQUIRE (FILE='pos.sph', EXIST=CONVERT)
                     INQUIRE (FILE='pos2.sph', EXIST=CONVERT2)
                     if(.not.convert) then
                        write(6,*)'Please run the routine'
                        write(6,*)'b2pos on the output file containing'
                        write(6,*)'the configuration of the first star.'
                        write(6,*)'It will produce a file called'
                        write(6,*)'pos.sph. Place it in this directory,'
                        write(6,*)'and the rerun testinput.'
                     endif
                     if(.not.convert2) then
                        if(convert) then
                           write(6,*)'Please run the routine'
                        else
                           write(6,*)'Please also run the routine'
                        endif
                        write(6,*)'b2pos2 on the output file containing'
                        write(6,*)'the configuration of the 2nd star.'
                        write(6,*)'It will produce a file called'
                        write(6,*)'pos2.sph. Put it in this directory,'
                        write(6,*)'and then rerun testinput.'
                     endif
                     if(convert.and.convert2) then
                        inittype='2qr'
                        call do2qr
                     else
                        stop
                     endif
                  endif
               endif
            else if (nspin.eq.3) then
               if(nq.eq.1) then
                  write(6,*)'Have you calculated a relaxed'
                  write(6,*)'single-star model? Type 1 for yes,'
                  write(6,*)'2 for no.'
                  read *,n1s
                  if(n1s.eq.1) then
                     INQUIRE (FILE='pos.sph', EXIST=CONVERT)
                     if(.not.convert) then
                        write(6,*)'Please run the routine b2pos on the'
                        write(6,*)'output file containing the single'
                        write(6,*)'configuration.  It will produce a'
                        write(6,*)'file called pos.sph.  Place it in'
                        write(6,*)'this directory and rerun testinput'
                        stop
                     else
                        inittype='hy1'
                        call dohy1
                     endif
                  else
                     write(6,*)'You will need to calculate such a model'
                     write(6,*)'first.  Please rerun testinput to'
                     write(6,*)'set up such a run.'
                     stop
                  endif
               else
                  write(6,*)'Have you calculated relaxed'
                  write(6,*)'single-star models? Type 1 for yes,'
                  write(6,*)'2 for no.'
                  read *,n1s
                  if(n1s.eq.1) then
                     INQUIRE (FILE='pos.sph', EXIST=CONVERT)
                     INQUIRE (FILE='pos2.sph', EXIST=CONVERT2)
                     if(.not.convert) then
                        write(6,*)'Please run the routine'
                        write(6,*)'b2pos on the output file containing'
                        write(6,*)'the configuration of the first star.'
                        write(6,*)'It will produce a file called'
                        write(6,*)'pos.sph. Place it in this directory,'
                        write(6,*)'and the rerun testinput.'
                     endif
                     if(.not.convert2) then
                        if(convert) then
                           write(6,*)'Please run the routine'
                        else
                           write(6,*)'Please also run the routine'
                        endif
                        write(6,*)'b2pos2 on the output file containing'
                        write(6,*)'the configuration of the 2nd star.'
                        write(6,*)'It will produce a file called'
                        write(6,*)'pos2.sph. Put it in this directory,'
                        write(6,*)'and then rerun testinput.'
                     endif
                     if(convert.and.convert2) then
                        inittype='2qr'
                        call do2qr
                     else
                        stop
                     endif
                  endif
               endif
            else
               write(6,*)'illegal choice of spin'
               stop
            endif
         else 
            write(6,*)'illegal number of stars'
            stop
         endif
      else if(nstart.eq.2) then
         write(6,*)'If you have not done so already, take the'
         write(6,*)'output file from the previous run, named'
         write(6,*)'restart.sph or outxxx.sph for some 3-digit'
         write(6,*)'integer "xxx", rename it restart.sph,'
         write(6,*)'and place it in this directory.'
         INQUIRE (FILE='restart.sph', EXIST=CONVERT)
         if(.not.convert) then
            write(6,*)'No file named restart.sph in this directory'
            stop
         else
            inittype='res'
            call dorestart
         endif
      else
         write(6,*)'illegal choice!'
         stop
      endif

      write(13,*)'INAME=''',inittype,''''

      write(12,*)'&END '
      close(12)
      write(13,*)'&END '
      close(13)

C      return
      end

**********************************************************************
      subroutine do1em
      include "testinput.h"

      sep0=0.0
      rp=0.
      vpeak=0.0
      qdar=0.0

      call doamrns(1,amns,rns)
      
      call donnnopt(1,n,nnopt)
 
      call dohminmax(hmin,hmax)

      call dorelax(1,nrelax,trelax)

      call dotf1(tf,nrelax,treloff)

      call dodtout(dtout)
 
      call dongr(ngr,xgrmax,ygrmax,zgrmax,xgrmin,ygrmin,zgrmin,
     $     xgrlim,ygrlim,zgrlim)

      call dogam(gam)

      if(treloff.lt.tf) then
         call donav(nav,alpha,beta,eta2)
         call dongravrad(ngravrad,sol)
      else
         nav=0
         alpha=0.
         beta=0.
         eta2=0.
         ngravrad=0
         sol=1000.0
      endif

      call dofiles

      return
      end
******************************************************************************
      subroutine do1es
      include "testinput.h"

      sep0=0.0
      rp=0.
      vpeak=0.0
      qdar=0.0

      call doamrns(1,amns,rns)
      
      call donnnopt(1,n,nnopt)
 
      call dohminmax(hmin,hmax)

      call dorelax(1,nrelax,trelax)

      call dotf1(tf,nrelax,treloff)

      call dodtout(dtout)
 
      call dongr(ngr,xgrmax,ygrmax,zgrmax,xgrmin,ygrmin,zgrmin,
     $     xgrlim,ygrlim,zgrlim)

      call dogam(gam)

      if(treloff.lt.tf) then
         call donav(nav,alpha,beta,eta2)
         call dongravrad(ngravrad,sol)
      else
         nav=0
         alpha=0.
         beta=0.
         eta2=0.
         ngravrad=0
         sol=1000.0
      endif

      call dofiles

      return
      end
******************************************************************************
      subroutine do2cm
      include "testinput.h"

      rp=0.
      vpeak=0.0
      qdar=1.0

      call doamrns(2,amns,rns)
      
      call dosep0(sep0)

      call donnnopt(2,n,nnopt)
 
      call dohminmax(hmin,hmax)

      call dorelax(2,nrelax,trelax)

      call dotf1(tf,nrelax,treloff)

      call dodtout(dtout)
 
      call dongr(ngr,xgrmax,ygrmax,zgrmax,xgrmin,ygrmin,zgrmin,
     $     xgrlim,ygrlim,zgrlim)

      call dogam(gam)

      call donav(nav,alpha,beta,eta2)

      call dongravrad(ngravrad,sol)

      call dofiles

      return
      end

******************************************************************************
      subroutine do2cs
      include "testinput.h"

      rp=0.
      vpeak=0.0
      qdar=1.0

      call doamrns(2,amns,rns)
      
      call dosep0(sep0)

      call donnnopt(2,n,nnopt)
 
      call dohminmax(hmin,hmax)

      call dorelax(2,nrelax,trelax)

      call dotf1(tf,nrelax,treloff)

      call dodtout(dtout)
 
      call dongr(ngr,xgrmax,ygrmax,zgrmax,xgrmin,ygrmin,zgrmin,
     $     xgrlim,ygrlim,zgrlim)

      call dogam(gam)

      call donav(nav,alpha,beta,eta2)

      call dongravrad(ngravrad,sol)

      call dofiles
         
      return
      end

******************************************************************************
      subroutine do2cr
      include"testinput.h"

      rp=0.
      vpeak=0.0
      qdar=1.0
      amns=0.0
      rns=0.0

      amtot=0.
      rmax=0.
      open(17,file='pos.sph', form='formatted')
      read(17,*)n,nnopt,
     $     hminold,hmaxold,gam,trelaxold
      do i=1,n
         read(17,*)x,y,z,am
         amtot=amtot+am
         rmax=max(rmax,sqrt(x**2+y**2+z**2))
      enddo
      close(17)


      write(6,*)'There will be',n
      write(6,*)'particles in each star.'
      write(6,*)'Make sure NMAX is set to at least ',2*n

      write(6,*)'The mass of each particle is contained in pos.sph.' 
      write(6,*)'The total mass of each star will be ',amtot

      write(6,*)'The postion of each particle is also in pos.sph'
      write(6,*)'The radius of each star is approximately',rmax
      
      write(6,*)'The code currently assumes that the stars have a'
      write(6,*)'polytropic EOS, i.e., P=k*rho**gamma.'
      write(6,*)'The value of gamma from pos.sph is',gam

      call dosep0(sep0)

      write(6,*)'NNOPT is the optimal number of neighbors for'
      write(6,*)'each particle.  NNOPT has been set to ',nnopt

      write(6,*)'HMIN is the minimum smoothing length'
      write(6,*)'A negative value means there is no minimum.'
      write(6,*)'Do not worry, this is perfectly valid.'
      write(6,*)'HMIN was set to ',hminold
      write(6,*)'We recommend you do not change this number.'
      write(6,*)'If you feel you must, type 2,'
      write(6,*)'to leave it alone, type 1 now.'
      read *,nchoice
      if(nchoice.eq.2) then
         write(6,*)'Type in the new value of HMIN'
         read *,hmin
      else
         hmin=hminold
      endif

      write(6,*)'HMAX is the maximum smoothing length'
      write(6,*)'HMAX was set to ',hmaxold
      write(6,*)'We recommend you do not change this number.'
      write(6,*)'If you feel you must, type 2.'
      write(6,*)'To leave it alone, type 1 now.'
      read *,nchoice
      if(nchoice.eq.2) then
         write(6,*)'Type in the new value of HMAX'
         read *,hmax
         if(hmax.le.0) then
            write(6,*)'Must have positive smoothing length'
            stop
         endif
      else
         hmax=hmaxold
      endif

      write(6,*)'Since this is a new run, you will want to'
      write(6,*)'relax the particles in the star.'
      nrelax=2
      write(6,*)'You will need to set the relaxation time, so that'
      write(6,*)'the drag acceleration=particle velocity/t_relax.'
      write(6,*)'It is generally given by t_relax=???'
      write(6,*)'TRELAX was set to ',trelaxold
      write(6,*)'You probably do not need to change this number.'
      write(6,*)'If you feel you must, type 2, to leave it alone'
      write(6,*)'type 1 now.'
      read *,nchoice
      if(nchoice.eq.2) then
         write(6,*)'Type in the new value of TRELAX'
         read *,trelax
         if(trelax.le.0) then
            write(6,*)'Must have positive relaxation time'
            stop
         endif
      else
         trelax=trelaxold
      endif

      call dotf1(tf,nrelax,treloff)

      call dodtout(dtout)
 
      call dongr(ngr,xgrmax,ygrmax,zgrmax,xgrmin,ygrmin,zgrmin,
     $     xgrlim,ygrlim,zgrlim)

      call donav(nav,alpha,beta,eta2)

      call dongravrad(ngravrad,sol)

      call dofiles

      return
      end

******************************************************************************
      subroutine do2qm
      include "testinput.h"

      rp=0.
      vpeak=0.0

      call doamrns(3,amns,rns)
      
      call dogam(gam)

      write(6,*)'What is the binary mass ratio,i.e. M_2/M_1?'
      write(6,*)'This should be a number <= 1.0'
      write(6,*)'Type the value now.'
      read *,qdar
      if(qdar.lt.0.0) then
         write(6,*)'Mass ratio must be positive!'
         stop
      endif
      if(qdar.gt.1.0) then
         write(6,*)'This value must be <= 1.0!!'
         stop
      endif

      write(6,*)'The mass of the secondary will be ',amns*qdar
      write(6,*)'The radius of the secondary will be',
     $     rns*qdar**((gam-2.0)/(3.0*gam-4.0))

      call dosep0(sep0)

      write(6,*)'The center of mass of the primary and secondary will'
      write(6,*)'be at ',-qdar*sep0/(1.0+qdar),' and ',sep0/(1.0+qdar)
      write(6,*)'respectively.'

      call donnnopt(3,n,nnopt)
      write(6,*)'The secondary will have ',int(qdar*n),' particles.'
      write(6,*)'Make sure NMAX is set larger than ',n+int(qdar*n)

 
      call dohminmax(hmin,hmax)

      call dorelax(2,nrelax,trelax)

      call dotf1(tf,nrelax,treloff)

      call dodtout(dtout)
 
      call dongr(ngr,xgrmax,ygrmax,zgrmax,xgrmin,ygrmin,zgrmin,
     $     xgrlim,ygrlim,zgrlim)

      call donav(nav,alpha,beta,eta2)

      call dongravrad(ngravrad,sol)

      call dofiles

      return
      end

******************************************************************************
      subroutine do2qs
      include "testinput.h"

      rp=0.
      vpeak=0.0

      call doamrns(3,amns,rns)

      call dogam(gam)

      write(6,*)'What is the binary mass ratio,i.e. M_2/M_1?'
      write(6,*)'This should be a number <= 1.0'
      write(6,*)'Type the value now.'
      read *,qdar
      if(qdar.lt.0.0) then
         write(6,*)'Mass ratio must be positive!'
         stop
      endif
      if(qdar.gt.1.0) then
         write(6,*)'This value must be <= 1.0!!'
         stop
      endif

      write(6,*)'The mass of the secondary will be ',amns*qdar
      write(6,*)'The radius of the secondary will be',
     $     rns*qdar**((gam-2.0)/(3.0*gam-4.0))

      call dosep0(sep0)

      write(6,*)'The center of mass of the primary and secondary will'
      write(6,*)'be at ',-qdar*sep0/(1.0+qdar),' and ',sep0/(1.0+qdar)
      write(6,*)'respectively.'

      call donnnopt(3,n,nnopt)
      write(6,*)'It is hard to give a guess as to the exact number of'
      write(6,*)'particles which will be used in the calculation, but'
      write(6,*)'to be safe, set NMAX larger than ',2*n

 
      call dohminmax(hmin,hmax)

      call dorelax(2,nrelax,trelax)

      call dotf1(tf,nrelax,treloff)

      call dodtout(dtout)
 
      call dongr(ngr,xgrmax,ygrmax,zgrmax,xgrmin,ygrmin,zgrmin,
     $     xgrlim,ygrlim,zgrlim)

      call donav(nav,alpha,beta,eta2)

      call dongravrad(ngravrad,sol)

      call dofiles

      return
      end

******************************************************************************
      subroutine do2qr
      include "testinput.h"

      rp=0.
      vpeak=0.0
      amns=0.0
      rns=0.0

      amtot=0.
      rmax=0.
      open(17,file='pos.sph', form='formatted')
      read(17,*)n,nnopt,
     $     hminold,hmaxold,gam,trelaxold
      do i=1,n
         read(17,*)x,y,z,am
         amtot=amtot+am
         rmax=max(rmax,sqrt(x**2+y**2+z**2))
      enddo
      close(17)

      write(6,*)'There will be ',n
      write(6,*)'particles in the primary, i.e. the more massive star.'

      write(6,*)'The mass of each particle in the primary is contained'
      write(6,*)'in pos.sph. The total mass of the star will be ',amtot

      write(6,*)'The postion of each particle is also in pos.sph'
      write(6,*)'The radius of the star is approximately',rmax
      
      amtot2=0.
      rmax=0.
      open(17,file='pos2.sph', form='formatted')
      read(17,*)n2,nnopt2,
     $     hmin2old,hmax2old,gam2,trelax2old
      do i=1,n2
         read(17,*)x,y,z,am
         amtot2=amtot2+am
         rmax=max(rmax,sqrt(x**2+y**2+z**2))
      enddo
      close(17)

      write(6,*)'There will be ',n2
      write(6,*)'particles in the secondary,i.e. the less massive star.'

      write(6,*)'The mass of each particle in the secondary is '
      write(6,*)'contained in pos2.sph. The total mass of the star'
      write(6,*)'will be ',amtot2
      write(6,*)'This gives a mass ratio q=',amtot2/amtot
      qdar=amtot2/amtot

      write(6,*)'The postion of each particle is also in pos2.sph'
      write(6,*)'The radius of the star is approximately',rmax

      write(6,*)'The code currently assumes that the stars have a'
      write(6,*)'polytropic EOS, i.e., P=k*rho**gamma.'
      if(gam.eq.gam2) then
         write(6,*)'The value of gamma from the input files is',gam
      else
         write(6,*)'Gamma is different for the two stars, with values'
         write(6,*)gam,' and ',gam2,' respectively.  This can cause'
         write(6,*)'huge problems if the values are truly different'
         write(6,*)'and not just off at machine precision.  We suggest'
         write(6,*)'rerunning your initial conditions.  If you really'
         write(6,*)'want to do this, we are setting gam=',gam
         write(6,*)'Edit sph.input to change this value.  Good luck!'
      endif 

      call dosep0(sep0)

      write(6,*)'The center of mass of the primary and secondary will'
      write(6,*)'be at ',-qdar*sep0/(1.0+qdar),' and ',sep0/(1.0+qdar)
      write(6,*)'respectively.'

      write(6,*)'NNOPT is the optimal number of neighbors for each'
      write(6,*)'particle.'
      if(nnopt.eq.nnopt2) then
         write(6,*)'The value of NNOPT from the input files is',nnopt
      else
         write(6,*)'NNOPT is different for the two stars, with values'
         write(6,*)nnopt,' and ',nnopt2,' respectively.  This can cause'
         write(6,*)'problems if the values are very different.'
         write(6,*)'You may want to rerun your stellar models.  If'
         write(6,*)'you want to do the run anyway, we are setting'
         write(6,*)'nnopt=',nnopt
         write(6,*)'Edit sph.input to change this value.  Good luck!'
      endif 

      write(6,*)'HMIN is the minimum smoothing length'
      write(6,*)'A negative value means there is no minimum.'
      write(6,*)'Do not worry, this is perfectly valid.'
      if((hminold.lt.0.and.hmin2old.lt.0).or.hminold.eq.hmin2old) then
         write(6,*)'HMIN was set to ',hminold
         write(6,*)'We recommend you do not change this number.'
         write(6,*)'If you feel you must, type 2,'
         write(6,*)'to leave it alone, type 1 now.'
         read *,nchoice
         if(nchoice.eq.2) then
            write(6,*)'Type in the new value of HMIN'
            read *,hmin
         else
            hmin=hminold
         endif
      else
         write(6,*)'HMIN is different for the two stars, with values'
         write(6,*)hminold,' and ',hminold2,' respectively.'  
         write(6,*)'This might cause'
         write(6,*)'problems if the values are very different.'
         write(6,*)'You may want to rerun your stellar models.  If'
         write(6,*)'you want to do the run anyway, we are setting'
         write(6,*)'hmin=',min(hminold,hmin2old)
         write(6,*)'Edit sph.input to change this value.  Good luck!'
         hmin=min(hminold,hmin2old)
      endif 

      write(6,*)'HMAX is the maximum smoothing length'
      if(hmaxold.eq.hmax2old) then
         write(6,*)'HMAX was set to ',hmaxold
         write(6,*)'We recommend you do not change this number.'
         write(6,*)'If you feel you must, type 2,'
         write(6,*)'to leave it alone, type 1 now.'
         read *,nchoice
         if(nchoice.eq.2) then
            write(6,*)'Type in the new value of HMAX'
            read *,hmax
            if(hmax.lt.0) then
               write(6,*)'Smoothing length must be positive!'
               stop
            endif
         else
            hmax=hmaxold
         endif
      else
         write(6,*)'HMAX is different for the two stars, with values'
         write(6,*)hmaxold,' and ',hmaxold2,' respectively.'  
         write(6,*)'This might cause'
         write(6,*)'problems if the values are very different.'
         write(6,*)'You may want to rerun your stellar models.  If'
         write(6,*)'you want to do the run anyway, we are setting'
         write(6,*)'hmax=',max(hmaxold,hmax2old)
         write(6,*)'Edit sph.input to change this value.  Good luck!'
         hmax=max(hmaxold,hmax2old)
      endif 

      write(6,*)'Since this is a new run, you will want to'
      write(6,*)'relax the particles in the star.'
      nrelax=2
      write(6,*)'You will need to set the relaxation time, so that'
      write(6,*)'the drag acceleration=particle velocity/t_relax.'
      write(6,*)'It is generally given by t_relax=???'
      if(trelaxold.eq.trelax2old) then
         write(6,*)'TRELAX was set to ',trelaxold
         write(6,*)'We recommend you do not change this number.'
         write(6,*)'If you feel you must, type 2,'
         write(6,*)'to leave it alone, type 1 now.'
         read *,nchoice
         if(nchoice.eq.2) then
            write(6,*)'Type in the new value of TRELAX'
            read *,trelax
            if(trelax.lt.0) then
               write(6,*)'Relaxation time must be positive!'
               stop
            endif
         else
            trelax=trelaxold
         endif
      else
         write(6,*)'TRELAX is different for the two stars, with values'
         write(6,*)trelaxold,' and ',trelax2old,' respectively.'  
         write(6,*)'Type 1 to use the value ',trelaxold
         write(6,*)'Type 2 to use the value ',trelax2old
         write(6,*)'Type 3 to use a new value '
         read *,nchoice
         if(nchoice.eq.1) then
            trelax=trelaxold
         else  if(nchoice.eq.2) then
            trelax=trelax2old
         else 
            write(6,*)'Type in the new value of TRELAX'
            read *,trelax
            if(trelax.lt.0) then
               write(6,*)'Relaxation time must be positive!'
               stop
            endif
         endif
      endif 

      call dotf1(tf,nrelax,treloff)

      call dodtout(dtout)
 
      call dongr(ngr,xgrmax,ygrmax,zgrmax,xgrmin,ygrmin,zgrmin,
     $     xgrlim,ygrlim,zgrlim)

      call donav(nav,alpha,beta,eta2)

      call dongravrad(ngravrad,sol)

      call dofiles

      return
      end

******************************************************************************
      subroutine do2i1
      include "testinput.h"

      rp=0.
      vpeak=0.
      qdar=1.0
      amns=0.0
      rns=0.0

      amtot=0.
      rmax=0.
      open(17,file='pos.sph', form='formatted')
      read(17,*)n,nnopt,
     $     hminold,hmaxold,gam,trelaxold
      do i=1,n
         read(17,*)x,y,z,am
         amtot=amtot+am
         rmax=max(rmax,sqrt(x**2+y**2+z**2))
      enddo
      close(17)


      write(6,*)'There will be',n
      write(6,*)'particles in each star.'
      write(6,*)'Make sure NMAX is set to at least ',2*n

      write(6,*)'The mass of each particle is contained in pos.sph.' 
      write(6,*)'The total mass of each star will be ',amtot

      write(6,*)'The postion of each particle is also in pos.sph'
      write(6,*)'The radius of each star is approximately',rmax
      
      write(6,*)'Since there is no easy way to relax an irrotational'
      write(6,*)'configuration, we recommend that you take a relaxed'
      write(6,*)'single-star configuration and deform it into a'
      write(6,*)'triaxial ellipsoid.  Appropriate values of the axis'
      write(6,*)'can be found in the Lai, Rasio, and Shapiro papers.'

      write(6,*)'Please enter the ratio of the new semi-major axis'
      write(6,*)'length to the radius of the single star model:',rmax
      write(6,*)'in the x-direction, i.e. the line connecting the'
      write(6,*)'centers-of-mass of the two stars.'
      read *,a1
      write(6,*)'Please enter the ratio of the new semi-major axis'
      write(6,*)'length to the radius of the single star model:',rmax
      write(6,*)'in the y-direction, i.e. the direction of motion'
      write(6,*)'of each two stars.'
      read *,a2
      write(6,*)'Please enter the ratio of the new semi-major axis'
      write(6,*)'length to the radius of the single star model:',rmax
      write(6,*)'in the z-direction, i.e. the rotation axis.'
      read *,a3
      open(17,file='axes1.sph',form='formatted')
      write(17,'(3e15.7)')a1,a2,a3
      close(17)

      write(6,*)'The code currently assumes that the stars have a'
      write(6,*)'polytropic EOS, i.e., P=k*rho**gamma.'
      write(6,*)'The value of gamma from pos.sph is',gam

      call dosep0(sep0)

      write(6,*)'NNOPT is the optimal number of neighbors for'
      write(6,*)'each particle.  NNOPT has been set to ',nnopt

      write(6,*)'HMIN is the minimum smoothing length'
      write(6,*)'A negative value means there is no minimum.'
      write(6,*)'Do not worry, this is perfectly valid.'
      write(6,*)'HMIN was set to ',hminold
      write(6,*)'We recommend you do not change this number.'
      write(6,*)'If you feel you must, type 2,'
      write(6,*)'to leave it alone, type 1 now.'
      read *,nchoice
      if(nchoice.eq.2) then
         write(6,*)'Type in the new value of HMIN'
         read *,hmin
      else
         hmin=hminold
      endif

      write(6,*)'HMAX is the maximum smoothing length'
      write(6,*)'HMAX was set to ',hmaxold
      write(6,*)'We recommend you do not change this number.'
      write(6,*)'If you feel you must, type 2.'
      write(6,*)'To leave it alone, type 1 now.'
      read *,nchoice
      if(nchoice.eq.2) then
         write(6,*)'Type in the new value of HMAX'
         read *,hmax
         if(hmax.le.0) then
            write(6,*)'Must have positive smoothing length'
            stop
         endif
      else
         hmax=hmaxold
      endif

      write(6,*)'There is currently no way to relax irrotational'
      write(6,*)'binary models.'
      nrelax=0
      trelax=0.0
      treloff=0.0

      call dotf2(tf)

      call dodtout(dtout)
 
      call dongr(ngr,xgrmax,ygrmax,zgrmax,xgrmin,ygrmin,zgrmin,
     $     xgrlim,ygrlim,zgrlim)

      call donav(nav,alpha,beta,eta2)

      call dongravrad(ngravrad,sol)

      call dofiles

      return
      end

******************************************************************************
      subroutine do2iq
      include "testinput.h"

      rp=0.
      vpeak=0.0
      amns=0.0
      rns=0.0

      amtot=0.
      rmax=0.
      open(17,file='pos.sph', form='formatted')
      read(17,*)n,nnopt,
     $     hminold,hmaxold,gam,trelaxold
      do i=1,n
         read(17,*)x,y,z,am
         amtot=amtot+am
         rmax=max(rmax,sqrt(x**2+y**2+z**2))
      enddo
      close(17)

      amtot2=0.
      rmax2=0.
      open(17,file='pos2.sph', form='formatted')
      read(17,*)n2,nnopt2,
     $     hmin2old,hmax2old,gam2,trelax2old
      do i=1,n2
         read(17,*)x,y,z,am
         amtot2=amtot2+am
         rmax2=max(rmax2,sqrt(x**2+y**2+z**2))
      enddo
      close(17)

      write(6,*)'There will be ',n
      write(6,*)'particles in the primary, i.e. the more massive star.'

      write(6,*)'The mass of each particle in the primary is contained'
      write(6,*)'in pos.sph. The total mass of the star will be ',amtot

      write(6,*)'The postion of each particle is also in pos.sph'
      write(6,*)'The radius of the star is approximately',rmax

      write(6,*)'There will be ',n2
      write(6,*)'particles in the secondary,i.e. the less massive star.'

      write(6,*)'The mass of each particle in the secondary is '
      write(6,*)'contained in pos2.sph. The total mass of the star'
      write(6,*)'will be ',amtot2
      write(6,*)'This gives a mass ratio q=',amtot2/amtot
      qdar=amtot2/amtot

      write(6,*)'The postion of each particle is also in pos2.sph'
      write(6,*)'The radius of the star is approximately',rmax2
      
      write(6,*)'Since there is no easy way to relax an irrotational'
      write(6,*)'configuration, we recommend that you take relaxed'
      write(6,*)'single-star configurations and deform them into'
      write(6,*)'triaxial ellipsoids.  Appropriate values of the axis'
      write(6,*)'lengths can be found in the'
      write(6,*)'Lai, Rasio, and Shapiro papers.'

      write(6,*)'First we do the primary.'
      write(6,*)'Please enter the ratio of the new semi-major axis'
      write(6,*)'length to the radius of the single star model:',rmax
      write(6,*)'in the x-direction, i.e. the line connecting the'
      write(6,*)'centers-of-mass of the two stars, for the PRIMARY.'
      read *,a1
      write(6,*)'Please enter the ratio of the new semi-major axis'
      write(6,*)'length to the radius of the single star model:',rmax
      write(6,*)'in the y-direction, i.e. the direction of motion'
      write(6,*)'of each two stars, for the PRIMARY.'
      read *,a2
      write(6,*)'Please enter the ratio of the new semi-major axis'
      write(6,*)'length to the radius of the single star model:',rmax
      write(6,*)'in the z-direction, i.e. the rotation axis'
      write(6,*)'for the primary star.'
      read *,a3
      open(17,file='axes1.sph',form='formatted')
      write(17,'(3e15.7)')a1,a2,a3
      close(17)


      write(6,*)'Now we do the secondary.'
      write(6,*)'Please enter the ratio of the new semi-major axis'
      write(6,*)'length to the radius of the single star model:',rmax2
      write(6,*)'in the x-direction, i.e. the line connecting the'
      write(6,*)'centers-of-mass of the two stars, for the secondary.'
      read *,b1
      write(6,*)'Please enter the ratio of the new semi-major axis'
      write(6,*)'length to the radius of the single star model:',rmax2
      write(6,*)'in the y-direction, i.e. the direction of motion'
      write(6,*)'of each two stars, for the secondary.'
      read *,b2
      write(6,*)'Please enter the ratio of the new semi-major axis'
      write(6,*)'length to the radius of the single star model:',rmax2
      write(6,*)'in the z-direction, i.e. the rotation axis'
      write(6,*)'for the secondary star.'
      read *,b3
      open(17,file='axes2.sph',form='formatted')
      write(17,'(3e15.7)')b1,b2,b3
      close(17)

      write(6,*)'The code currently assumes that the stars have a'
      write(6,*)'polytropic EOS, i.e., P=k*rho**gamma.'
      if(gam.eq.gam2) then
         write(6,*)'The value of gamma from the input files is',gam
      else
         write(6,*)'Gamma is different for the two stars, with values'
         write(6,*)gam,' and ',gam2,' respectively.  This can cause'
         write(6,*)'huge problems if the values are truly different'
         write(6,*)'and not just off at machine precision.  We suggest'
         write(6,*)'rerunning your initial conditions.  If you really'
         write(6,*)'want to do this, we are setting gam=',gam
         write(6,*)'Edit sph.input to change this value.  Good luck!'
      endif 

      call dosep0(sep0)

      write(6,*)'The center of mass of the primary and secondary will'
      write(6,*)'be at ',-qdar*sep0/(1.0+qdar),' and ',sep0/(1.0+qdar)
      write(6,*)'respectively.'

      write(6,*)'NNOPT is the optimal number of neighbors for each'
      write(6,*)'particle.'
      if(nnopt.eq.nnopt2) then
         write(6,*)'The value of NNOPT from the input files is',nnopt
      else
         write(6,*)'NNOPT is different for the two stars, with values'
         write(6,*)nnopt,' and ',nnopt2,' respectively.  This can cause'
         write(6,*)'problems if the values are very different.'
         write(6,*)'You may want to rerun your stellar models.  If'
         write(6,*)'you want to do the run anyway, we are setting'
         write(6,*)'nnopt=',nnopt
         write(6,*)'Edit sph.input to change this value.  Good luck!'
      endif 

      write(6,*)'HMIN is the minimum smoothing length'
      write(6,*)'A negative value means there is no minimum.'
      write(6,*)'Do not worry, this is perfectly valid.'
      if((hminold.lt.0.and.hmin2old.lt.0).or.hminold.eq.hmin2old) then
         write(6,*)'HMIN was set to ',hminold
         write(6,*)'We recommend you do not change this number.'
         write(6,*)'If you feel you must, type 2,'
         write(6,*)'to leave it alone, type 1 now.'
         read *,nchoice
         if(nchoice.eq.2) then
            write(6,*)'Type in the new value of HMIN'
            read *,hmin
         else
            hmin=hminold
         endif
      else
         write(6,*)'HMIN is different for the two stars, with values'
         write(6,*)hminold,' and ',hminold2,' respectively.'  
         write(6,*)'This might cause'
         write(6,*)'problems if the values are very different.'
         write(6,*)'You may want to rerun your stellar models.  If'
         write(6,*)'you want to do the run anyway, we are setting'
         write(6,*)'hmin=',min(hminold,hmin2old)
         write(6,*)'Edit sph.input to change this value.  Good luck!'
         hmin=min(hminold,hmin2old)
      endif 

      write(6,*)'HMAX is the maximum smoothing length'
      if(hmaxold.eq.hmax2old) then
         write(6,*)'HMAX was set to ',hmaxold
         write(6,*)'We recommend you do not change this number.'
         write(6,*)'If you feel you must, type 2,'
         write(6,*)'to leave it alone, type 1 now.'
         read *,nchoice
         if(nchoice.eq.2) then
            write(6,*)'Type in the new value of HMAX'
            read *,hmax
            if(hmax.lt.0) then
               write(6,*)'Smoothing length must be positive!'
               stop
            endif
         else
            hmax=hmaxold
         endif
      else
         write(6,*)'HMAX is different for the two stars, with values'
         write(6,*)hmaxold,' and ',hmaxold2,' respectively.'  
         write(6,*)'This might cause'
         write(6,*)'problems if the values are very different.'
         write(6,*)'You may want to rerun your stellar models.  If'
         write(6,*)'you want to do the run anyway, we are setting'
         write(6,*)'hmax=',max(hmaxold,hmax2old)
         write(6,*)'Edit sph.input to change this value.  Good luck!'
         hmax=max(hmaxold,hmax2old)
      endif 

      write(6,*)'There is currently no way to relax irrotational'
      write(6,*)'binary models.'
      nrelax=0
      trelax=0.0
      treloff=0.0

      call dotf2(tf)

      call dodtout(dtout)
 
      call dongr(ngr,xgrmax,ygrmax,zgrmax,xgrmin,ygrmin,zgrmin,
     $     xgrlim,ygrlim,zgrlim)

      call donav(nav,alpha,beta,eta2)

      call dongravrad(ngravrad,sol)

      call dofiles

      return
      end

******************************************************************************
      subroutine dohy1

      qdar=1.0
      amns=0.0
      rns=0.0
      nrelax=0.0
      trelax=0.0
      treloff=0.0

      amtot=0.
      rmax=0.
      open(17,file='pos.sph', form='formatted')
      read(17,*)n,nnopt,
     $     hminold,hmaxold,gam,trelaxold
      do i=1,n
         read(17,*)x,y,z,am
         amtot=amtot+am
         rmax=max(rmax,sqrt(x**2+y**2+z**2))
      enddo
      close(17)


      write(6,*)'There will be',n
      write(6,*)'particles in each star.'
      write(6,*)'Make sure NMAX is set to at least ',2*n

      write(6,*)'The mass of each particle is contained in pos.sph.' 
      write(6,*)'The total mass of each star will be ',amtot

      write(6,*)'The postion of each particle is also in pos.sph'
      write(6,*)'The radius of each star is approximately',rmax
      
      write(6,*)'NNOPT is the optimal number of neighbors for'
      write(6,*)'each particle.  NNOPT has been set to ',nnopt

      write(6,*)'HMIN is the minimum smoothing length'
      write(6,*)'A negative value means there is no minimum.'
      write(6,*)'Do not worry, this is perfectly valid.'
      write(6,*)'HMIN was set to ',hminold
      write(6,*)'We recommend you do not change this number.'
      write(6,*)'If you feel you must, type 2,'
      write(6,*)'to leave it alone, type 1 now.'
      read *,nchoice
      if(nchoice.eq.2) then
         write(6,*)'Type in the new value of HMIN'
         read *,hmin
      else
         hmin=hminold
      endif

      write(6,*)'HMAX is the maximum smoothing length'
      write(6,*)'HMAX was set to ',hmaxold
      write(6,*)'We recommend you do not change this number.'
      write(6,*)'If you feel you must, type 2.'
      write(6,*)'To leave it alone, type 1 now.'
      read *,nchoice
      if(nchoice.eq.2) then
         write(6,*)'Type in the new value of HMAX'
         read *,hmax
         if(hmax.le.0) then
            write(6,*)'Must have positive smoothing length'
            stop
         endif
      else
         hmax=hmaxold
      endif

      write(6,*)'Hyperbolic collision trajectories can be defined'
      write(6,*)'by a two-parameter family of initial states once'
      write(6,*)'the stellar model is fixed.  Here, we use the'
      write(6,*)'periastron separation and '
      write(6,*)'the center of mass velocity.'
      write(6,*)'What is the periastron separation of the system,'
      write(6,*)'which we will use to extrapolate the full orbit.'
      write(6,*)'Write the value now, in code units.'
      read *,rp
      write(6,*)'The periastron separation is ',rp,' code units,'
      write(6,*)'or in stellar radii:',rp/rmax

      write(6,*)'Please enter the CoM velocity of the system.'
      write(6,*)'i.e., the time derivative of |r_1-r_2|, at'
      write(6,*)'large separations.'
      read *,vpeak      

      write(6,*)'Using these values, we can calculate the proper'
      write(6,*)'point-mass trajectory for the stars.  However,'
      write(6,*)'we lack the resources to start the run at infinite'
      write(6,*)'separation.  Thus, you must enter the binary'
      write(6,*)'separation, in code units, at which to start'
      write(6,*)'the run.  Enter the value now.'
      read *,sep0
      write(6,*)'You have chosen an initial separation of',sep0
      write(6,*)'code units, or in stellar radii, r_0=',sep0/rmax

      amtot2=2*amtot
      amu=amtot/2.0
      rimp=sqrt(rp**2+2*amtot2/vpeak**2*rp)

      vphi=sqrt(amtot2/rimp)
      etot=0.5*amu*vpeak**2
      altot=amu*vpeak*rimp

      vmin=vphi**2/vpeak+sqrt(vpeak**2+vphi**4/vpeak**2)

      v0=sqrt(vpeak**2+2*amtot2/sep0)

      if(myrank.eq.0)write(6,*)'The Impact Parameter will be:',rimp
      if(myrank.eq.0)write(6,*)'The periastron velcity is:',vmin
      if(myrank.eq.0)write(6,*)'The initial velocity is',v0

      call dotf2(tf)

      call dodtout(dtout)
 
      call dongr(ngr,xgrmax,ygrmax,zgrmax,xgrmin,ygrmin,zgrmin,
     $     xgrlim,ygrlim,zgrlim)

      call donav(nav,alpha,beta,eta2)

      call dongravrad(ngravrad,sol)

      call dofiles
         
      return
      end

******************************************************************************
      subroutine dohyq
      include "testinput.h"

      amns=0.0
      rns=0.0
      nrelax=0.0
      trelax=0.0
      treloff=0.0

      amtot1=0.
      rmax1=0.
      open(17,file='pos.sph', form='formatted')
      read(17,*)n,nnopt,
     $     hminold,hmaxold,gam,trelaxold
      do i=1,n
         read(17,*)x,y,z,am
         amtot1=amtot1+am
         rmax1=max(rmax1,sqrt(x**2+y**2+z**2))
      enddo
      close(17)

      amtot2=0.
      rmax2=0.
      open(17,file='pos2.sph', form='formatted')
      read(17,*)n2,nnopt2,
     $     hmin2old,hmax2old,gam2,trelax2old
      do i=1,n2
         read(17,*)x,y,z,am
         amtot2=amtot2+am
         rmax2=max(rmax2,sqrt(x**2+y**2+z**2))
      enddo
      close(17)

      write(6,*)'There will be ',n2
      write(6,*)'particles in the secondary,i.e. the less massive star.'

      write(6,*)'The mass of each particle in the secondary is '
      write(6,*)'contained in pos2.sph. The total mass of the star'
      write(6,*)'will be ',amtot2
      write(6,*)'This gives a mass ratio q=',amtot2/amtot1
      qdar=amtot2/amtot1

      write(6,*)'The postion of each particle is also in pos2.sph'
      write(6,*)'The radius of the star is approximately',rmax2

      write(6,*)'There will be ',n
      write(6,*)'particles in the primary, i.e. the more massive star.'

      write(6,*)'The mass of each particle in the primary is contained'
      write(6,*)'in pos.sph. The total mass of the star will be ',amtot1

      write(6,*)'The postion of each particle is also in pos.sph'
      write(6,*)'The radius of the star is approximately',rmax1

      write(6,*)'NNOPT is the optimal number of neighbors for each'
      write(6,*)'particle.'
      if(nnopt.eq.nnopt2) then
         write(6,*)'The value of NNOPT from the input files is',nnopt
      else
         write(6,*)'NNOPT is different for the two stars, with values'
         write(6,*)nnopt,' and ',nnopt2,' respectively.  This can cause'
         write(6,*)'problems if the values are very different.'
         write(6,*)'You may want to rerun your stellar models.  If'
         write(6,*)'you want to do the run anyway, we are setting'
         write(6,*)'nnopt=',nnopt
         write(6,*)'Edit sph.input to change this value.  Good luck!'
      endif 

      write(6,*)'HMIN is the minimum smoothing length'
      write(6,*)'A negative value means there is no minimum.'
      write(6,*)'Do not worry, this is perfectly valid.'
      if((hminold.lt.0.and.hmin2old.lt.0).or.hminold.eq.hmin2old) then
         write(6,*)'HMIN was set to ',hminold
         write(6,*)'We recommend you do not change this number.'
         write(6,*)'If you feel you must, type 2,'
         write(6,*)'to leave it alone, type 1 now.'
         read *,nchoice
         if(nchoice.eq.2) then
            write(6,*)'Type in the new value of HMIN'
            read *,hmin
         else
            hmin=hminold
         endif
      else
         write(6,*)'HMIN is different for the two stars, with values'
         write(6,*)hminold,' and ',hminold2,' respectively.'  
         write(6,*)'This might cause'
         write(6,*)'problems if the values are very different.'
         write(6,*)'You may want to rerun your stellar models.  If'
         write(6,*)'you want to do the run anyway, we are setting'
         write(6,*)'hmin=',min(hminold,hmin2old)
         write(6,*)'Edit sph.input to change this value.  Good luck!'
         hmin=min(hminold,hmin2old)
      endif 

      write(6,*)'HMAX is the maximum smoothing length'
      if(hmaxold.eq.hmax2old) then
         write(6,*)'HMAX was set to ',hmaxold
         write(6,*)'We recommend you do not change this number.'
         write(6,*)'If you feel you must, type 2,'
         write(6,*)'to leave it alone, type 1 now.'
         read *,nchoice
         if(nchoice.eq.2) then
            write(6,*)'Type in the new value of HMAX'
            read *,hmax
            if(hmax.lt.0) then
               write(6,*)'Smoothing length must be positive!'
               stop
            endif
         else
            hmax=hmaxold
         endif
      else
         write(6,*)'HMAX is different for the two stars, with values'
         write(6,*)hmaxold,' and ',hmaxold2,' respectively.'  
         write(6,*)'This might cause'
         write(6,*)'problems if the values are very different.'
         write(6,*)'You may want to rerun your stellar models.  If'
         write(6,*)'you want to do the run anyway, we are setting'
         write(6,*)'hmax=',max(hmaxold,hmax2old)
         write(6,*)'Edit sph.input to change this value.  Good luck!'
         hmax=max(hmaxold,hmax2old)
      endif 


      write(6,*)'The code currently assumes that the stars have a'
      write(6,*)'polytropic EOS, i.e., P=k*rho**gamma.'
      if(gam.eq.gam2) then
         write(6,*)'The value of gamma from the input files is',gam
      else
         write(6,*)'Gamma is different for the two stars, with values'
         write(6,*)gam,' and ',gam2,' respectively.  This can cause'
         write(6,*)'huge problems if the values are truly different'
         write(6,*)'and not just off at machine precision.  We suggest'
         write(6,*)'rerunning your initial conditions.  If you really'
         write(6,*)'want to do this, we are setting gam=',gam
         write(6,*)'Edit sph.input to change this value.  Good luck!'
      endif 

      write(6,*)'Hyperbolic collision trajectories can be defined'
      write(6,*)'by a two-parameter family of initial states once'
      write(6,*)'the stellar model is fixed.  Here, we use the'
      write(6,*)'periastron separation and'
      write(6,*)'the center of mass velocity.'
      write(6,*)'What is the periastron separation of the system,'
      write(6,*)'which we will use to extrapolate the full orbit.'
      write(6,*)'Write the value now, in code units.'
      read *,rp
      write(6,*)'The periastron separation is ',rp,' code units.'

      write(6,*)'Please enter the CoM velocity of the system.'
      write(6,*)'i.e., the time derivative of |r_1-r_2|, at'
      write(6,*)'large separations.'
      read *,vpeak      

      write(6,*)'Using these values, we can calculate the proper'
      write(6,*)'point-mass trajectory for the stars.  However,'
      write(6,*)'we lack the resources to start the run at infinite'
      write(6,*)'separation.  Thus, you must enter the binary'
      write(6,*)'separation, in code units, at which to start'
      write(6,*)'the run.  Enter the value now.'
      read *,sep0
      write(6,*)'You have chosen an initial separation of',sep0

      amtot=amtot+amtot2
      amu=amtot1*amtot2/amtot
      rimp=sqrt(rp**2+2*amtot/vpeak**2*rp)

      vphi=sqrt(amtot/rimp)
      etot=0.5*amu*vpeak**2
      altot=amu*vpeak*rimp

      vmin=vphi**2/vpeak+sqrt(vpeak**2+vphi**4/vpeak**2)

      v0=sqrt(vpeak**2+2*amtot/sep0)

      if(myrank.eq.0)write(6,*)'The Impact Parameter will be:',rimp
      if(myrank.eq.0)write(6,*)'The periastron velcity is:',vmin
      if(myrank.eq.0)write(6,*)'The initial velocity is',v0

      call dotf2(tf)

      call dodtout(dtout)
 
      call dongr(ngr,xgrmax,ygrmax,zgrmax,xgrmin,ygrmin,zgrmin,
     $     xgrlim,ygrlim,zgrlim)

      call donav(nav,alpha,beta,eta2)

      call dongravrad(ngravrad,sol)

      call dofiles

      return
      end

******************************************************************************
      subroutine dorestart
      include"testinput.h"

      OPEN(17,FILE='restart.sph',FORM='UNFORMATTED')
      READ(17) N,NNOPTOLD,HMINOLD,HMAXOLD,GAMOLD,SEP0OLD,
     $     TFOLD,DTOUTOLD,NOUT,NLEFT,NIT,Told,
     $     NAVOLD,ALPHAOLD,BETAOLD,ETA2OLD,
     $     NGROLD,XGRMINOLD,XGRMAXOLD,YGRMINOLD,YGRMAXOLD,
     $     ZGRMINOLD,ZGRMAXOLD,
     $     NRELAXOLD,TRELAXOLD,initgr,ngravradold,solold,
     $     q2xx,q2xy,q2xz,q2yy,q2yz,q2zz

      rp=0.
      vpeak=0.0
      qdar=0.0
      amns=0.0
      rns=0.0

      write(6,*)'The previous run ended at time T=',told
      write(6,*)'Please type the time at which this run will end.'
      read *,tf
      if(tf.le.told) then
         write(6,*)'Illegal choice, run must end after it starts!'
         stop
      endif

      nrelax=nrelaxold
      if(nrelaxold.eq.1) then
         sep0=0.
         write(6,*)'The previous run ended during a single star'
         write(6,*)'relaxation.  Should the current one do so'
         write(6,*)'as well?  Type 1 for yes, 2 for no.'
         read(6,*)nreloff
         if(nreloff.eq.1) then
            treloff=tf+0.5
         else if (nreloff.eq.2) then
            write(6,*)'Enter a time between ',told,' and ',tf
            write(6,*)'at which relaxation will be turned off.'
            read *,treloff
            if(treloff.le.t.or.treloff.gt.tf) then
               write(6,*)'Illegal value!'
               stop
            endif
         else
            write(6,*)'Illegal value!'
            stop
         endif
         write(6,*)'TRELAX is the relaxation timescale.'
         write(6,*)'TRELAX has been set to ',trelaxold
         write(6,*)'We will assume you wish to keep TRELAX the same'
         trelax=trelaxold
      else if(nrelaxold.eq.2) then
         write(6,*)'When the previous run finished, you were relaxing'
         write(6,*)'a corotating binary with separation ',sep0old
         write(6,*)'We will assume you wish to keep SEP0 the same'
         sep0=sep0old
         write(6,*)'Should the current one keep relaxation on'
         write(6,*)'until the end?  Type 1 for yes, 2 for no.'
         read(6,*)nreloff
         if(nreloff.eq.1) then
            treloff=tf+0.5
         else if (nreloff.eq.2) then
            write(6,*)'Enter a time between ',told,' and ',tf
            write(6,*)'at which relaxation will be turned off.'
            read *,treloff
            if(treloff.le.t.or.treloff.gt.tf) then
               write(6,*)'Illegal value!'
               stop
            endif
         else
            write(6,*)'Illegal value!'
            stop
         endif
         write(6,*)'TRELAX is the relaxation timescale.'
         write(6,*)'TRELAX has been set to ',trelaxold
         write(6,*)'We will assume you wish to keep TRELAX the same'
         trelax=trelaxold
      else
         sep0=0.0
         nrelax=0
         trelax=0.0
         treloff=0.0
      endif
      
      write(6,*)'There will be',n
      write(6,*)'total particles in the calculation.'
      write(6,*)'Make sure NMAX is set to at least ',n

      write(6,*)'The code currently assumes that the stars have a'
      write(6,*)'polytropic EOS, i.e., P=k*rho**gamma.'
      write(6,*)'The value of gamma from restart.sph is',gamold
      write(6,*)'We will assume that you want to keep this value'
      write(6,*)' of the parameter GAM'
      gam=gamold

      write(6,*)'Type any number to continue'
      read *,dum

      write(6,*)'NNOPT is the optimal number of neighbors for'
      write(6,*)'each particle.  NNOPT has been set to ',nnoptold
      write(6,*)'We will assume you wish to keep NNOPT the same'
      nnopt=nnoptold

      write(6,*)'HMIN is the minimum smoothing length'
      write(6,*)'A negative value means there is no minimum.'
      write(6,*)'Do not worry, this is perfectly valid.'
      write(6,*)'HMIN was set to ',hminold
      write(6,*)'We recommend you do not change this number.'
      write(6,*)'If you feel you must, type 2,'
      write(6,*)'to leave it alone, type 1 now.'
      read *,nchoice
      if(nchoice.eq.2) then
         write(6,*)'Type in the new value of HMIN'
         read *,hmin
      else
         hmin=hminold
      endif

      write(6,*)'HMAX is the maximum smoothing length'
      write(6,*)'HMAX was set to ',hmaxold
      write(6,*)'We recommend you do not change this number.'
      write(6,*)'If you feel you must, type 2.'
      write(6,*)'To leave it alone, type 1 now.'
      read *,nchoice
      if(nchoice.eq.2) then
         write(6,*)'Type in the new value of HMAX'
         read *,hmax
         if(hmax.le.0) then
            write(6,*)'Must have positive smoothing length'
            stop
         endif
      else
         hmax=hmaxold
      endif

      write(6,*)'DTOUT is the time between writing output files.'
      write(6,*)'DTOUT has been set to ',dtoutold
      write(6,*)'We will assume you wish to keep DTOUT the same'
      dtout=dtoutold

      write(6,*)'You have previously set the following values for'
      write(6,*)'the gravity grid, which we will copy for you.'
      write(6,*)'NGR=',ngrold
      ngr=ngrold
      write(6,*)'XGRMAX=',xgrmaxold
      xgrmax=xgrmaxold
      write(6,*)'XGRMIN=',xgrminold
      xgrmin=xgrminold
      write(6,*)'YGRMAX=',ygrmaxold
      ygrmax=ygrmaxold
      write(6,*)'YGRMIN=',ygrminold
      ygrmin=ygrminold
      write(6,*)'ZGRMAX=',zgrmaxold
      zgrmax=zgrmaxold
      write(6,*)'ZGRMIN=',zgrminold
      zgrmin=zgrminold
      write(6,*)'XGRLIM=',xgrlimold
      xgrlim=xgrlimold
      write(6,*)'YGRLIM=',ygrlimold
      ygrlim=ygrlimold
      write(6,*)'ZGRLIM=',zgrlimold
      zgrlim=zgrlimold

      write(6,*)'Type any number to continue'
      read *,dum

      nav=navold
      alpha=alphaold
      beta=betaold
      eta2=eta2old
      write(6,*)'NAV sets the choice of the artificial viscosity'
      write(6,*)'routine: 0 for no AV, 1 for Balsara,'
      write(6,*)'2 for Hernquist and Katz,'
      write(6,*)'3 for Monaghan'
      write(6,*)'You have previously used NAV=',nav
      if(nav.eq.0)then
         write(6,*)'AV will be off'
      else if (nav.eq.1) then
         write(6,*)'ALPHA will be set to ',alpha
         write(6,*)'We recommend alpha=0.5'
         write(6,*)'BETA will be set to ',beta
         write(6,*)'We recommend beta=0.5'
         write(6,*)'ETA2 will be set to ',eta2
         write(6,*)'We recommend eta2=0.01'
      else if (nav.eq.2) then
         write(6,*)'ALPHA will be set to ',alpha
         write(6,*)'We recommend alpha=0.5'
         write(6,*)'BETA will be set to ',beta
         write(6,*)'We recommend beta=1.0'
         write(6,*)'ETA2 will be set to ',eta2
         write(6,*)'We recommend eta2=0.01'
      else if (nav.eq.3) then
         write(6,*)'ALPHA will be set to ',alpha
         write(6,*)'We recommend alpha=0.5'
         write(6,*)'BETA will be set to ',beta
         write(6,*)'We recommend beta=0.5'
         write(6,*)'ETA2 will be set to ',eta2
         write(6,*)'We recommend eta2=0.01'
      endif

      ngravrad=ngravradold
      sol=solold
      write(6,*)'NGRAVRAD determines whether you will use radiation'
      write(6,*)'reaction: 0 for no, 1 for yes'
      write(6,*)'You have previously used NGRAVRAD=',ngravrad
      if(ngravrad.eq.1) then
         write(6,*)'SOL is the speed of light in code units.'
         write(6,*)'You have previously used SOL=',sol
      endif

      call dofiles

      return
      end

*************************************************************************

      subroutine dotf1(tf,nrelax,treloff)
      implicit none
      integer nrelax,nchoice
      real tf,treloff

      write(6,*)'Times are measured in terms of the'
      write(6,*)'Dynamical timescale, t_d=???'
      write(6,*)'At what time should this run terminate?'
      write(6,*)'Please type the value now.'
      read *,tf
      if(tf.le.0) then
         write(6,*)'Must have positive ending time'
         stop
      endif
      if(nrelax.eq.0) then
         treloff=0.0
      else
         write(6,*)'Do you want relaxation on for the entire run?'
         write(6,*)'Type one for yes, type 2 for no.'
         read *,nchoice
         if(nchoice.eq.1) then
            treloff=tf+0.05
         else if (nchoice.eq.2) then
            write(6,*)'Enter the time at which relaxation is shut off.'
            read *,treloff
            if(treloff.le.0)then
               write(6,*)'Relaxation will never start w/negative value!'
            else if (treloff.gt.tf) then
               write(6,*)'Relaxation turns off after run ends!'
            endif
         else
            write(6,*)'Illegal choice!'
            stop
         endif
      endif

      return
      end

**************************************************************************8

      subroutine dotf2(tf)
      implicit none
      real tf

      write(6,*)'Times are measured in terms of the'
      write(6,*)'Dynamical timescale, t_d=???'
      write(6,*)'At what time should this run terminate?'
      write(6,*)'Please type the value now.'
      read *,tf
      if(tf.le.0) then
         write(6,*)'Must have positive ending time'
         stop
      endif

      return
      end

******************************************************************

      subroutine dodtout(dtout)
      implicit none
      real dtout

      write(6,*)'This program is equipped to write a checkpoint'
      write(6,*)'file at periodic intervals.  How often would you'
      write(6,*)'like such a file to be written' 
      write(6,*)'(size=90 bytes/particle)'
      write(6,*)'Type the value of the time interval now.'
      read *,dtout
      if(dtout.le.0) then
         write(6,*)'Must have positive time interval'
         stop
      endif
      return
      end

******************************************************************

      subroutine dogam(gam)
      implicit none
      real gam

      write(6,*)'The code currently assumes that the stars have a'
      write(6,*)'polytropic EOS, i.e., P=k*rho**gamma.'
      write(6,*)'Please type the value of the adiabatic index gamma'
      read *,gam

      return
      end

************************************************************************

      subroutine donnnopt(nflag,n,nnopt)
      implicit none
      integer nflag,n,nnopt

      write(6,*)'How many particles will be used in'
      if(nflag.eq.1) then
         write(6,*)'the run?'
      else if (nflag.eq.2) then
         write(6,*)'each star?'
      else if (nflag.eq.3) then
         write(6,*)'the primary?'
      endif
      write(6,*)'Please type in the number as an integer'
      read *,n
      if(n.le.0) then
         write(6,*)'Must have positive number of particles!'
         stop
      endif

      write(6,*)'Please choose the optimal number of neighbors'
      write(6,*)'each particle will have.  Make sure that NNMAX'
      write(6,*)'in the spha.h file is approximately twice as big.'
      write(6,*)'Generally, a value of NNOPT=50-100 is good for a'
      write(6,*)'calculation using 10^4-10^5 particles.'
      write(6,*)'Type the number of neighbors desired.'
      read *,nnopt
      if(nnopt.le.0) then
         write(6,*)'Must have positive number of neighbors!'
         stop
      endif

      return
      end

*************************************************************************

      subroutine donav(nav,alpha,beta,eta2)
      implicit none
      integer nav
      real alpha,beta,eta2

      write(6,*)'Now we do the artificial viscosity.'
      write(6,*)'Type 0 for no AV, 1 for Balsara''s form,'
      write(6,*)'2 for Monaghan, 3 for Hernquist and Katz.'
      read *,nav
      if(nav.lt.0.or.nav.gt.3) then
         write(6,*)'illegal choice'
         stop
      else if (nav.gt.0.and.nav.le.3) then
         write(6,*)'We will need to set the AV parameters'
         write(6,*)'See astro-ph/9807290 for more info'
      endif
      if(nav.eq.0) then
         alpha=0.0
         beta=0.0
         eta2=0.0
      else if (nav.eq.1) then
         write(6,*)'Type in alpha (we recommend alpha=0.5)'
         read *,alpha
         write(6,*)'Type in beta (we recommend beta=0.5)'
         read *,beta
         write(6,*)'Type in eta^2 (we recommend eta^2=0.01)'
         read *,eta2
      else if (nav.eq.2) then
         write(6,*)'Type in alpha (we recommend alpha=0.5)'
         read *,alpha
         write(6,*)'Type in beta (we recommend beta=1.0)'
         read *,beta
         write(6,*)'Type in eta^2 (we recommend eta^2=0.01)'
         read *,eta2
      else if (nav.eq.3) then
         write(6,*)'Type in alpha (we recommend alpha=0.5)'
         read *,alpha
         write(6,*)'Type in beta (we recommend beta=0.5)'
         read *,beta
         write(6,*)'Type in eta^2 (we recommend eta^2=0.01)'
         read *,eta2
      endif

      return
      end

**************************************************************************

      subroutine dongr(ngr,xgrmax,ygrmax,zgrmax,xgrmin,ygrmin,zgrmin,
     $     xgrlim,ygrlim,zgrlim)

      implicit none
      integer ngr,ncube
      real xgrmax,ygrmax,zgrmax,xgrmin,ygrmin,zgrmin
      real xgrlim,ygrlim,zgrlim,xmax

      write(6,*)'Now we set up the gravity grid'
      write(6,*)'Please choose among the following options'
      write(6,*)'Type 1 for a fixed gravity grid'
      write(6,*)'  '
      write(6,*)'Type 2 to have the grid evolve to just contain'
      write(6,*)'     all SPH particles'
      write(6,*)'  '
      write(6,*)'Type 21 top start off lile NGR=2, but switch to'
      write(6,*)'NGR=1 after some time you will set.'
      write(6,*)'  '
      write(6,*)'Type a number between 90 and 99 to start off like'
      write(6,*)'NGR=1, but resize the grid to cover that percentage'
      write(6,*)'of the mass once a particle leaves the box'
      write(6,*)'  '
      write(6,*)'Type a number between 290 and 299 to start off like'
      write(6,*)'NGR=2, but resize the grid to cover that number,'
      write(6,*)'minus 200, percent of the mass once a particle'
      write(6,*)'leaves a box whose dimensions you set.'
      write(6,*)'  '
      write(6,*)'Type 398 to have the grid always cover 98% of the'
      write(6,*)'mass once the particle leaves a box defined by you,'
      write(6,*)'expanding until it reaches a limiting value of twice'
      write(6,*)'the original size.'
      write(6,*)'  '
      write(6,*)'Please note that option 1 is the most efficient.'
      read *,ngr

      if(ngr.eq.1)then
         write(6,*)'Do you want a cubic grid centered at the origin?'
         write(6,*)'Type 1 for yes, 2 for no.'
         read *,ncube
         if(ncube.eq.1) then
            write(6,*)'Please type the distance between the'
            write(6,*)'origin and the boundary, i.e. half the'
            write(6,*)'length of one side of the cube.'
            read *,xmax
            xgrmax=xmax
            ygrmax=xmax
            zgrmax=xmax
            xgrlim=xmax
            ygrlim=xmax
            zgrlim=xmax
            xgrmin=-1.0*xmax
            ygrmin=-1.0*xmax
            zgrmin=-1.0*xmax
         else
            write(6,*)'Please type the positive x-direction extent'
            read *,xgrmax
            write(6,*)'Please type the negative x-direction extent'
            write(6,*)'Please type it as a negative number if that is'
            write(6,*)'what you want.'
            read *,xgrmin
            if(xgrmin.gt.0.0) write(6,*)'The origin is not in the'
            if(xgrmin.gt.0.0) write(6,*)'gravity grid.  Good luck.'
            write(6,*)'Please type the positive y-direction extent'
            read *,ygrmax
            write(6,*)'Please type the negative y-direction extent'
            read *,ygrmin
            if(ygrmin.gt.0.0) write(6,*)'The origin is not in the'
            if(ygrmin.gt.0.0) write(6,*)'gravity grid.  Good luck.'
            write(6,*)'Please type the positive z-direction extent'
            read *,zgrmax
            write(6,*)'Please type the negative z-direction extent'
            read *,zgrmin
            if(zgrmin.gt.0.0) write(6,*)'The origin is not in the'
            if(zgrmin.gt.0.0) write(6,*)'gravity grid.  Good luck.'
         endif
      else if (ngr.eq.2) then
         write(6,*)'You have chosen to have the gravity grid resize'
         write(6,*)'itself every iteration to cover all SPH particles.'
         xgrmax=0.
         ygrmax=0.
         zgrmax=0.
         xgrmin=0.
         ygrmin=0.
         zgrmin=0.
         xgrlim=0.
         ygrlim=0.
         zgrlim=0.
      else if (ngr.eq.21) then
         write(6,*)'You have chosen for the grid to initially cover all'
         write(6,*)'all particles, but switch after a time to fixed'
         write(6,*)'boundaries.  This time is set with the parameter'
         write(6,*)'TGR in grav.f, grav2.f, and gravrad.f'
         write(6,*)'We assume the grid will be a cube centered on the'
         write(6,*)'origin, please edit the files if you want'
         write(6,*)'something different.'
         write(6,*)'Please type the distance between the'
         write(6,*)'origin and the boundary, i.e. half the'
         write(6,*)'length of one side of the cube.'
         read *,xmax
         xgrmax=xmax
         ygrmax=xmax
         zgrmax=xmax
         xgrlim=xmax
         ygrlim=xmax
         zgrlim=xmax
         xgrmin=-1.0*xmax
         ygrmin=-1.0*xmax
         zgrmin=-1.0*xmax
      else if (ngr.ge.90.and.ngr.le.99) then
         write(6,*)'The grid will cover ',ngr,'% of the matter after'
         write(6,*)'a particle leaves the grid defined by XGRLIM,'
         write(6,*)'YGRLIM, and ZGRLIM.  We assume here that they are'
         write(6,*)'all equal, please edit the sph.input file if they'
         write(6,*)'are not. Please type the distance between the'
         write(6,*)'origin and the boundary, i.e. half the'
         write(6,*)'length of one side of the cube.'
         read *,xmax
         xgrmax=xmax
         ygrmax=xmax
         zgrmax=xmax
         xgrlim=xmax
         ygrlim=xmax
         zgrlim=xmax
         xgrmin=-1.0*xmax
         ygrmin=-1.0*xmax
         zgrmin=-1.0*xmax
      else if (ngr.ge.290.and.ngr.le.299) then
         write(6,*)'Initially, the grid will be sized to cover all'
         write(6,*)'particles.'
         write(6,*)'The grid will cover ',ngr-200,'% of the matter once'
         write(6,*)'a particle leaves the grid defined by XGRLIM,'
         write(6,*)'YGRLIM, and ZGRLIM.  We assume here that they are'
         write(6,*)'all equal, please edit the sph.input file if they'
         write(6,*)'are not. Please type the distance between the'
         write(6,*)'origin and the boundary, i.e. half the'
         write(6,*)'length of one side of the cube.'
         read *,xmax
         xgrmax=xmax
         ygrmax=xmax
         zgrmax=xmax
         xgrlim=xmax
         ygrlim=xmax
         zgrlim=xmax
         xgrmin=-1.0*xmax
         ygrmin=-1.0*xmax
         zgrmin=-1.0*xmax
      else if (ngr.eq.398)then
         write(6,*)'The grid will cover ',ngr,'% of the matter after'
         write(6,*)'a particle leaves the grid defined by XGRLIM,'
         write(6,*)'YGRLIM, and ZGRLIM.  It will then expand until'
         write(6,*)'it reaches a size twice as large in any dimension.'
         write(6,*),'We assume the original box is cubic,'
         write(6,*)'please edit the sph.input file if it'
         write(6,*)'is not. Please type the original distance between'
         write(6,*)'the origin and the boundary, i.e. half the'
         write(6,*)'length of one side of the cube.'
         write(6,*)'The maximum box size will be twice as large.'
         read *,xmax
         xgrmax=xmax
         ygrmax=xmax
         zgrmax=xmax
         xgrlim=xmax
         ygrlim=xmax
         zgrlim=xmax
         xgrmin=-1.0*xmax
         ygrmin=-1.0*xmax
         zgrmin=-1.0*xmax
       else
         write(6,*)'please see the documentation, and set the values'
         write(6,*)'of the following parameters in sph.input by hand:'
         write(6,*)'xgrmin,xgrmax,ygrmin,ygrmax,zgrmin,zgrmax'
         write(6,*)'xgrlim,ygrlim,zgrlim'
         xmax=0.
         xgrmax=xmax
         ygrmax=xmax
         zgrmax=xmax
         xgrlim=xmax
         ygrlim=xmax
         zgrlim=xmax
         xgrmin=-1.0*xmax
         ygrmin=-1.0*xmax
         zgrmin=-1.0*xmax
      endif

      return
      end

******************************************************************************

      subroutine dohminmax(hmin,hmax)
      implicit none
      real hmin,hmax

      write(6,*)'Please set the minimum smoothing length for each'
      write(6,*)'particle.  For no fixed minimum, please enter a'
      write(6,*)'negative number.  This is perfectly acceptable.'
      write(6,*)'Type the value of the min. smoothing length.'
      read *,hmin

      write(6,*)'Please set the maximum smoothing length for each'
      write(6,*)'particle.  A value equal to the star''s radius is a'
      write(6,*)'safe choice.'
      write(6,*)'Type the value of the max. smoothing length.'
      read *,hmax
      if(hmax.le.0) then
         write(6,*)'Must have positive smoothing length'
         stop
      endif

      return
      end

**************************************************************

      subroutine dorelax(nflag,nrelax,trelax)
      implicit none
      integer nflag,nrelax
      real trelax
      
      if(nflag.eq.1) then
         write(6,*)'Since this is a new run, you probably want to'
      else if (nflag.eq.2) then
         write(6,*)'Since this is a new run, you will want to'
      endif
      write(6,*)'relax the particles in the star.'
      if(nflag.eq.1) then
         write(6,*)'Type 1 for a single-star relaxation,'
         write(6,*)'type 0 for no relaxation'
         read *,nrelax
         if(nrelax.ne.0.and.nrelax.ne.1) then
            write(6,*)'illegal parameter!'
            stop
         endif
      else if (nflag.eq.2) then
         nrelax=2
      endif   
      if(nrelax.eq.0) then
         trelax=0
      else 
         write(6,*)'You will need to set the relaxation time, so that'
         write(6,*)'the drag acceleration=particle velocity/t_relax.'
         write(6,*)'It is generally given by t_relax=???'
         write(6,*)'Please type the value of the relaxation time'
         read *,trelax
         if(trelax.le.0) then
            write(6,*)'Must have positive relaxation time'
            stop
         endif
      endif

      return
      end

************************************************************

      subroutine dosep0(sep0)
      implicit none
      real sep0
      
      write(6,*)'What is the initial separation of the stars'
      write(6,*)'with respect to each other (not with respect'
      write(6,*)'to the origin!).  Please type the value'
      write(6,*)'of the initial separation.'
      read *,sep0
      if(sep0.le.0.0) then
         write(6,*)'Initial separation must be positive!'
         stop
      endif

      return
      end

**************************************************************

      subroutine doamrns(nflag,amns,rns)
      implicit none
      integer nflag
      real amns,rns

      if(nflag.eq.1) then
         write(6,*)'What is the mass of your star, in the units'
      else if (nflag.eq.2) then
         write(6,*)'What is the mass of each star, in the units'
      else if (nflag.eq.3) then
         write(6,*)'What is the mass of the primary, i.e. the more'
         write(6,*)'massive star in the binary, in the units'
      endif
      write(6,*)'you use?  Type the value now.'
      read *,amns
      if(amns.le.0) then
         write(6,*)'Must have positive mass'
         stop
      endif

      write(6,*)'What is the radius of that star, in the units'
      write(6,*)'you use?  Type the value now.'
      read *,rns
      if(rns.le.0) then
         write(6,*)'Must have positive radius'
         stop
      endif

      return
      end

****************************************************************

      subroutine dongravrad(ngravrad,sol)
      implicit none
      integer ngravrad
      real sol

      write(6,*)'Is radiation reaction going to be included in'
      write(6,*)'this run?  It will only be used when relaxation is'
      write(6,*)'turned off.  If you want to use radiation reaction,'
      write(6,*)'Type 1, if not type 0.'
      read *,ngravrad
      if(ngravrad.eq.1)then
         write(6,*)'You will need to set the value of the speed'
         write(6,*)'of light being used in your calculation.'
         write(6,*)'Calculate the physical value of GM/rc^2 for'
         write(6,*)'your stars, where M and r are the physical'
         write(6,*)'values of the unit mass and length you use.'
         write(6,*)'Take this value to the -0.5 power.  That is the'
         write(6,*)'speed of light for the calculation.'
         write(6,*)'Type the value now.'
         read *,sol
      else if(ngravrad.eq.0)then
         sol=1000.0
      else
         write(6,*)'invalid parameter'
      endif

      return
      end

****************************************************************************

      subroutine dofiles
      implicit none
      include "testinput.h"

      write(12,*)'TF=',tf,','
      write(12,*)'DTOUT=',dtout,','
      write(12,*)'GAM=',gam,','
      write(12,*)'N=',n,','
      write(12,*)'NNOPT=',nnopt,','
      write(12,*)'NAV=',nav,','
      write(12,*)'ALPHA=',alpha,','
      write(12,*)'BETA=',beta,','
      write(12,*)'ETA2=',eta2,','
      write(12,*)'NGR=',ngr,','
      write(12,*)'XGRMIN=',xgrmin,','
      write(12,*)'XGRMAX=',xgrmax,','
      write(12,*)'YGRMIN=',ygrmin,','
      write(12,*)'YGRMAX=',ygrmax,','
      write(12,*)'ZGRMIN=',zgrmin,','
      write(12,*)'ZGRMAX=',zgrmax,','
      write(12,*)'XGRLIM=',xgrlim,','
      write(12,*)'YGRLIM=',ygrlim,','
      write(12,*)'ZGRLIM=',zgrlim,','
      write(12,*)'HMIN=',hmin,','
      write(12,*)'HMAX=',hmax,','
      write(12,*)'NRELAX=',nrelax,','
      write(12,*)'TRELAX=',trelax,','
      write(12,*)'SEP0=',sep0,','
      write(12,*)'QDAR=',qdar,','
      write(12,*)'AMNS=',amns,','
      write(12,*)'RNS=',rns,','
      write(12,*)'RP=',rp,','
      write(12,*)'VPEAK=',vpeak,','
      write(12,*)'TRELOFF=',treloff,','
      write(12,*)'NGRAVRAD=',ngravrad,','
      write(12,*)'SOL=',sol,','

      return
      end
