      SUBROUTINE SHOCK_TUBE
      IMPLICIT NONE

c     ********************************************************
c     This program initializes the configuration for a shock
c     tube test of an SPH code. It places particles in a tube,
c     the axis of which is the x-axis. The particles are 
c     evenly spaced, with the tube extending from x=0.0 to
c     x=1.0. These initial conditions are meant for a shock 
c     along the x-direction. Every particle in the tube also 
c     has 8 "ghost" particles associated with it. In addition,
c     there are end caps that will remain stationary. 
c     The particles are ordered as follows:
c     1. The value of x is set, and the value of y is set.
c     2. For these values of x and y, particles are placed 
c        along the z-direction.
c     3. Then y is incremented, and particles are again
c        placed along the z-direction.
c     4. After that is done, x is incremented and the process
c        is repeated.
c     5. Then, in the order of the real particles already
c        placed, the ghost particles are added. The 8 ghost
c        particles are added to the first real particle, then
c        the next real particle, and so on.
c     6. Then the left end cap is added in the region to the 
c        left of x=0.0, in the same order as step 1 through 
c        step 4. It is added so that the x values increase, 
c        so they go from more negative to less negative.
c     7. Last, the right end cap is added to the right of 
c        x=1.0, again following the same procedure as step 1
c        through step 4. It is also ordered so that x 
c        increases from 1.0 onward.
c     ********************************************************
      INCLUDE 'shock.h'
      INCLUDE 'spha.h'

C      INTEGER nwidth, nx, nghost, np, ntube, nend, ntot
      INTEGER i, j, k,l         !,nstmax!  index, index2, index3, index4
C      PARAMETER (nstmax=1000000)
!      REAL*8 space
      REAL*8 z_1(nstmax)

      OPEN(129,FILE='shock.input',ERR=120)
      READ(129,SHOCKINPUT)
      CLOSE(129)
      
      OPEN(12,FILE='sph.input',ERR=100)
      READ(12,INPUT)
      CLOSE(12)

      
      OPEN(123,FILE='splash.txt')
      OPEN(125,FILE='splash_left.txt')
      OPEN(126,FILE='splash_right.txt')
      OPEN(124,FILE='test.txt')

      write(6,*)'Shock pars -- left:',vx_l,press_l,rho_l
      write(6,*)'Shock pars -- right:',vx_r,press_r,rho_r

!      vx_l= 0.0
!      vx_r= 0.0
!      press_l= 1.0
!      press_r= 0.1 
!      rho_l= 1.0
!      rho_r= 0.125

c     Number of particles in endcaps
      nwidth = 100!2*nnopt
c     Number of particles in x direction
      nx = 1000
c     Number of real particles
      np = nx
c     Number of particles+ghost particles
      ntube=np
c     Number of particles on one end
      nend = nwidth
c     Total number of particles
      ntot = nx+2*nend
      ntimestepper = -1
      
!      write(124,*),nx,nwidth,nnopt

c     Spacing for particles placed from 0 to 1 on x-axis
      space = 1.0/(nx)

      DO i=1,nx
         index=i               
         x(index)=(i-0.5)*space-0.5
         y(index)=0            
         z(index)=0            
         if (x(index).le.0) THEN
            am(index)=rho_l*space
            a(index)=press_l/rho_l**gam
            vx(index)=vx_l
            rho(index)=rho_l
         ELSE
            am(index)=rho_r*space
            a(index)=press_r/rho_r**gam
            vx(index)=vx_r
            rho(index)=rho_r
         end if
         hp(index)=space*nnopt/4.0
         vy(index)=0
         vz(index)=0
         write(123,'(i6,8E15.7)'),index,a(index),en(index),
     $        rho(index),x(index),y(index),z(index)
      ENDDO
    
      
c     Now we need the particles at the ends of the 
c     tube that do not move
      
c     First, the left end
      
      DO i=nwidth,1,-1
         index3=nx+nwidth-i+1
         x(index3)=-0.5-space*(i-0.5)
         y(index3)=0
         z(index3)=0
         am(index3)=rho_l*space
         a(index3)=press_l/rho_l**gam
         vx(index3)=vx_l
         vy(index3)=0
         vz(index3)=0
         rho(index3)=rho_l
         hp(index3)=space*nnopt/4.0
         write(125,'(i6,9E15.7)'),index3,x(index3),y(index3),
     $        z(index3),vx(index3),am(index3),a(index3),hp(index3),
     $        en(index3),rho(index3)
      ENDDO

c     Now the right end

      DO i=1,nwidth
         index4=nx+nwidth+i
         x(index4)=0.5+space*(i-0.5)
         y(index4)=0
         z(index4)=0
         am(index4)=rho_r*space
         a(index4)=press_r/rho_r**gam
         vx(index4)=vx_r
         vy(index4)=0
         vz(index4)=0 
         rho(index4)=rho_r
         hp(index4)=space*nnopt/4.0
         write(126,'(i6,9E15.7)'),index4,a(index4),en(index4),
     $        rho(index4)
      ENDDO

c     Tests 

      call e_from_a

      do i=1,ntot
      write(124,'(i6,9E15.7)'),i,x(i),y(i),z(i),vx(i),a(i),en(I)
      
      end do

      n=np+2*nwidth
      nleft=n
      call gravquant
      CALL LFSTART


      
c     Tests

 !     PRINT *, "Control Variable Values"
 !     PRINT *,

 !     PRINT *, "Number of particles along x= ", nx
 !     PRINT *, "Number of particles along y and z= ", nwidth
 !     PRINT *, "Total number of particles= ", ntot
 !     PRINT *, "spacing= ", space

 !     PRINT *
 !     PRINT *, "***************************************"
 !     PRINT *
 !     PRINT *, "All Values"
 !     PRINT *

c     I wanted to do this test because for a while I was
c     not filling the arrays properly; I had a bunch of
c     zeros at the end, before I reached ntot. Of course
c     the arrays will have zeros at the end because the
c     arrays are bigger than they need to be, but I had
c     too many zeros.

  !    DO i=1,ntot
  !       PRINT *, i, x(i), y(i), z(i)
  !    ENDDO

 !     PRINT *
 !     PRINT *, "***************************************"
 !     PRINT *
 !     PRINT *, "Values Left of Tube"
 !     PRINT *

c     I did this test specifically to make sure that
c     the stationary particles to the left of the tube
c     are placed properly, from left to right along the
c     x-axis

!      DO i=ntube+1,ntube+nend
!         PRINT *, i, x(i), y(i), z(i)
!      ENDDO
      RETURN 
 100  STOP '2cs: ERROR OPENING INPUT FILE ???'
      CLOSE(123)
 120  STOP ' ERROR OPENING SHOCK INPUT FILE ???'
      END 
