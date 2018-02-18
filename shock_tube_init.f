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
      
      OPEN(123,FILE='splash.txt')
      OPEN(124,FILE='splash_ghost.txt')
      OPEN(125,FILE='splash_left.txt')
      OPEN(126,FILE='splash_right.txt')

      write(6,*)'Shock pars -- left:',vx_l,press_l,rho_l
      write(6,*)'Shock pars -- right:',vx_r,press_r,rho_r

!      vx_l= 0.0
!      vx_r= 0.0
!      press_l= 1.0
!      press_r= 0.1 
!      rho_l= 1.0
!      rho_r= 0.125

c     Number of particles in y and z directions
      nwidth = 10
c     Number of particles in x direction
      nx = 100
c     Number of real particles
      np = nx*nwidth*nwidth
c     Number of sets of ghost particles
      nghost = 8
c     Number of particles+ghost particles
      ntube=(nghost+1)*np
c     Number of particles on one end
      nend = 9*nwidth*nwidth*nwidth
c     Total number of particles
      ntot = ntube+2*nend
      ntimestepper = -1
      
c     Spacing for particles placed from 0 to 1 on x-axis
      space = 1.0/(nx)
      OPEN(12,FILE='sph.input',ERR=100)
      READ(12,INPUT)
      CLOSE(12)

      DO i=1,nx
         DO j=1,nwidth
            DO k=1,nwidth
               index=1+(i-1)*nx+(j-1)*nwidth+(k-1)
               x(index)=(i+0.5)*space-0.5
               y(index)=space*(-(nwidth+1)/2.0+j)
               z(index)=space*(-(nwidth+1)/2.0+k)
               if (x(index).le.0) THEN
                  am(index)=rho_l*space**3
                  a(index)=press_l/rho_l**gam
                  vx(index)=vx_l
               ELSE
                  am(index)=rho_r*space**3
                  a(index)=press_r/rho_r**gam
                  vx(index)=vx_r
               end if
               hp(index)=space*(nnopt/32.0)**(1.0/3.0)
               vy(index)=0
               vz(index)=0
               write(123,'(i6,8E15.7)'),index,x(index),y(index),
     $              z(index),vx(index),am(index),a(index),hp(index)
            ENDDO
         ENDDO
      ENDDO

c     Now we need the ghost particles that correspond to
c     the real particles.
      
      do i=1,ntot
         if (i.le.np) THEN
         else if (i.gt.ntube) then
            nref(i)=0
         end if
      end do
      
      DO j=1,np
         DO k=1,nghost
            index2=np+j*nghost-8+k
            nref(index2)=j
            x(index2)=x(j)
            SELECT CASE (k)
            CASE (1)
               y(index2)=y(j)+nwidth*space
               z(index2)=z(j)
            CASE (2)
               y(index2)=y(j)+nwidth*space
               z(index2)=z(j)+nwidth*space
            CASE (3)
               y(index2)=y(j)
               z(index2)=z(j)+nwidth*space
            CASE (4)
               y(index2)=y(j)-nwidth*space
               z(index2)=z(j)+nwidth*space
            CASE (5)
               y(index2)=y(j)-nwidth*space
               z(index2)=z(j)
            CASE (6)
               y(index2)=y(j)-nwidth*space
               z(index2)=z(j)-nwidth*space
            CASE (7)
               y(index2)=y(j)
               z(index2)=z(j)-nwidth*space
            CASE (8)
               y(index2)=y(j)+nwidth*space
               z(index2)=z(j)-nwidth*space
            END SELECT
            am(index2)=am(j)
            a(index2)=a(j)
            vx(index2)=vx(j)
            hp(index2)=space*(nnopt/32.0)**(1.0/3.0)
            vy(index2)=0
            vz(index2)=0

            write(124,'(i6,8E15.7)'),index2,x(index2),y(index2),     
     $           z(index2),vx(index2),am(index2),a(index2),hp(index2)
         ENDDO
      ENDDO
      
c     Now we need the particles at the ends of the 
c     tube that do not move
      
c     First, the left end
      
      DO i=nwidth,1,-1
         DO j=1,3*nwidth
            DO k=1,3*nwidth
               index3=ntube+9*(nwidth-i)*nwidth*nwidth+3*(j-1)*nwidth
     $              +k
               x(index3)=-0.5-space*(i+0.5)
               y(index3)=space*(-(nwidth+1)/2.0+j-nwidth)
               z(index3)=space*(-(nwidth+1)/2.0+k-nwidth)
               am(index3)=rho_l*space**3
               a(index3)=press_l/rho_l**gam
               vx(index3)=vx_l
               vy(index3)=0
               vz(index3)=0
               hp(index3)=space*(nnopt/32.0)**(1.0/3.0)
               write(125,'(i6,8E15.7)'),index3,x(index3),y(index3),
     $              z(index3),vx(index3),am(index3),a(index3),hp(index3)
              ENDDO
         ENDDO
      ENDDO

c     Now the right end

      DO i=1,nwidth
         DO j=1,3*nwidth
            DO k=1,3*nwidth
              index4=ntube+nend+9*(i-1)*nwidth*nwidth+3*(j-1)*nwidth+k             
               x(index4)=0.5+space*(i+0.5)
               y(index4)=space*(-(nwidth+1)/2.0+j-nwidth)
               z(index4)=space*(-(nwidth+1)/2.0+k-nwidth)
               am(index4)=rho_r*space**3
               a(index4)=press_r/rho_r**gam
               vx(index4)=vx_r
               vy(index4)=0
               vz(index4)=0
               hp(index4)=space*(nnopt/32.0)**(1.0/3.0)
               write(126,'(i6,8E15.7)'),index4,x(index4),y(index4),
     $              z(index4),vx(index4),am(index4),a(index4),hp(index4)
            ENDDO
         ENDDO
      ENDDO


c     Tests  

      n=ntube+2*nend
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
