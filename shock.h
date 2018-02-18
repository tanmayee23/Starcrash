
C      PARAMETERS 
       INTEGER nstmax
       PARAMETER (nstmax=1000000)

C     Shock: Shock tube  
      INTEGER nwidth, nx, nghost, np, ntube, nend, ntot,
     $          index, index1, index2, index3, index4
      REAL space
C     Shock_init: Init
      REAL rho_l, rho_r, vx_l, vx_r, press_l, press_r

      COMMON/Shock/ nwidth,nx,np,ntube,nend,ntot,index,
     $      index1,index2,index3,index4,
     $      space
      COMMON/Shock_init/ rho_l,rho_r,vx_l,vx_r,press_l,press_r
      NAMELIST/SHOCKINPUT/ vx_l,vx_r,press_l,press_r,rho_l,rho_r
