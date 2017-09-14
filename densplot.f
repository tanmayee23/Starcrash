**************************************************************
c     This file is part of the StarCrash code
c     Version 1.0
c
c     Copyright 2003, by Joshua Faber
c     This code is publicly available under the terms of the
c     GNU General Public License.  See our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

      subroutine densplot
      IMPLICIT NONE
*******************************************************
c     RELEASE 1.0
c     plots a 3-d SPH density grid in an integer array
c     called by OUTPUT
****************************************************
      INTEGER NGRIDX,NGRIDZ
      REAL xmin,xmax,zmin,zmax
c     The output file is an NGRIDX by NGRIDX by NGRIDZ array, with physical 
c     dimensions running from xmin < x <xmax, xmin < y < xmax, zmin < z < zmax
      PARAMETER(NGRIDX=96,XMIN=-4.0,XMAX=4.0)
      PARAMETER(NGRIDZ=48,zmin=-2.0,zmax=2.0)

      include 'spha.h'

      CHARACTER*12 outfile
      real maxrho
c     PARAMETER(MAXRHO= a value you set by hand)
      real rhogrid(ngridx,ngridx,ngridz)
      integer rhobyte(ngridx,ngridx,ngridz)

      REAL hp2,hp3,xtry,ytry,ztry,r2,uu2,wtab2,wijini
      INTEGER i,j,k,nxgr,nxrange,nygr,nzgr,ii,nxi,jj,nyj,kk,nzk

c     initialize arrays
      do i=1,ngridx
         do j=1,ngridx
            do k=1,ngridz
               rhogrid(i,j,k)=0.0
               rhobyte(i,j,k)=0
            enddo
         enddo
      enddo
      write(outfile,102) nit
 102  FORMAT('dens',I4.4,'.dat')
     
      OPEN(118,FILE=outfile, form='unformatted')

      DO I=1,N
         nxgr=int((x(i)-XMIN)/(XMAX-XMIN)*NGRIDx)
         nxrange=int(2.0*hp(i)*ngridx/(xmax-xmin))+1
         nygr=int((y(i)-XMIN)/(XMAX-XMIN)*NGRIDx)
         nzgr=int((z(i)-ZMIN)/(ZMAX-ZMIN)*NGRIDz)
         hp2=hp(i)**2
         hp3=hp(i)**3
         do ii=-nxrange,nxrange
            nxi=nxgr+ii
            if(nxi.ge.1.and.nxi.le.ngridx) then
               xtry=(xmax-xmin)*(nxi)/(1.0*ngridx)+xmin            
               do jj=-nxrange,nxrange
                  nyj=nygr+jj
                  if(nyj.ge.1.and.nyj.le.ngridx) then
                     ytry=(xmax-xmin)*(nyj)/(1.0*ngridx)+xmin
                     do kk=-nxrange,nxrange
                        nzk=nzgr+kk
                        if(nzk.ge.1.and.nzk.le.ngridz) then
                           ztry=(zmax-zmin)*(nzgr+kk)/
     $                          (1.0*ngridz)+zmin
                           r2=(x(i)-xtry)**2+(y(i)-ytry)**2+
     $                          (z(i)-ztry)**2
                           uu2=sqrt(r2/hp2)
                           if(Uu2.LT.1.) THEN
                              Wtab2=1.-1.5*Uu2**2+0.75*Uu2**3
                           ELSE IF (Uu2.LT.2.) THEN
                              Wtab2=0.25*(2.-Uu2)**3
                           ELSE
                              Wtab2=0.
                           ENDIF
                           Wtab2=Wtab2/3.14159265359
                           wijini=wtab2/hp3
c     WE calculate the density contribution of all particles to all
c     the grid points which their kernel function overlaps
                           rhogrid(nxi,nyj,nzk)=
     $                          rhogrid(nxi,nyj,nzk)+am(i)*wijini
                        endif
                     enddo
                  endif
               enddo
            endif
         ENDDO
      enddo

c     If you wist to set MAXRHO as a parameter, comment out the following 
c     block
c*******************************************************************
      maxrho = 0.0
      do i=1,ngridx
         do j=1,ngridx
            do k=1,ngridz
               maxrho = max(maxrho,rhogrid(i,j,k))
            enddo
         enddo
      enddo
c*****************************************************************

      do i=1,ngridx
         do j=1,ngridx
            do k=1,ngridz
               rhogrid(i,j,k) = (254.9/maxrho)*rhogrid(i,j,k)
               rhobyte(i,j,k) = int(rhogrid(i,j,k)) 
            enddo
         enddo
      enddo

      write (118)rhobyte
      return
      end
