**************************************************************
c     This file is part of the StarCrash code
c     Version 1.0
c
c     Copyright 2003, by Joshua Faber
c     This code is publicly available under the terms of the
c     GNU General Public License.  See our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

      SUBROUTINE CMADJ                                                 
      IMPLICIT NONE
******************************************************************
c     Release 1.0
c     Adjusts system center-of-mass for relaxation runs
c_______________________________________________________________
c     This version scans in the distance separating a binary
c     After time-TSCANON passes, the binary separation decreases linearly
c     reaching a final value of SEPFINAL at the end of the run, or 
c     the time at which the relaxation is turned off, whichever comes first
c_______________________________________________________________
c     Called by ADVANCE
*******************************************************************
      INCLUDE 'spha.h'                                          
      PARAMETER(TSCANON=5.0,SEPFINAL=3.0)
      REAL xcm1,ycm1,zcm1,am1,xcm2,ycm2,zcm2,am2
      REAL delx1,dely1,delz1,delx2,dely2,delz2
      REAL TSCANON,SEPFINAL,SEP1,tscanoff
      INTEGER i    

      tscanoff=min(TF,TRELOFF)

      if(t.gt.TSCANON) then
         sep1=sep0+(SEPFINAL-sep0)*(t-TSCANON)/(TSCANOFF-TSCANON)
         if(myrank.eq.0)write(6,*)'Binary separation ',sep1,' at time=',t
      else
         sep1=sep0
      endif

      XCM1=0.
      YCM1=0.
      ZCM1=0.
      AM1=0.                                                           
      xcm2=0.
      ycm2=0.
      zcm2=0.
      am2=0.
      DO I=1,NLEFT   
         AM1=AM1+AM(I)
         XCM1=XCM1+AM(I)*X(I)                                           
         YCM1=YCM1+AM(I)*Y(I)                                           
         ZCM1=ZCM1+AM(I)*Z(I)                                           
      ENDDO                                                            
      XCM1=XCM1/AM1                                                    
      YCM1=YCM1/AM1                                                    
      ZCM1=ZCM1/AM1  

      if(nrelax.eq.2) then
         DO I=NLEFT+1,N   
            AM2=AM2+AM(I)
            XCM2=XCM2+AM(I)*X(I)                                           
            YCM2=YCM2+AM(I)*Y(I)                                           
            ZCM2=ZCM2+AM(I)*Z(I)                                           
         ENDDO                                                            
         XCM2=XCM2/AM2                                                    
         YCM2=YCM2/AM2                                                    
         ZCM2=ZCM2/AM2 
      endif

      if(nrelax.eq.1) then
         DELX1=-XCM1
      else if(nrelax.eq.2) then
         DELX1=-XCM1-AM2*SEP1/(AM1+AM2)                                   
      endif
      DELY1=-YCM1
      DELZ1=-ZCM1
      if(nrelax.eq.2) then
         DELX2=-XCM2+AM1*SEP1/(AM1+AM2)                                   
         DELY2=-YCM2
         DELZ2=-ZCM2
      endif

      if(myrank.eq.0)write(6,*)'CMADJ: Star 1:',delx1,dely1,delz1
      if(nrelax.eq.2.and.myrank.eq.0)
     $     write(6,*)'CMADJ: Star 2:',delx2,dely2,delz2
                                                      
      DO I=1,NLEFT
        X(I)=X(I)+DELX1
        Y(I)=Y(I)+DELY1
        Z(I)=Z(I)+DELZ1
      ENDDO                                                            
      if(nrelax.eq.2) then
         DO I=NLEFT+1,N
            X(I)=X(I)+DELX2
            Y(I)=Y(I)+DELY2
            Z(I)=Z(I)+DELZ2                                                
         ENDDO                                                            
      endif                                  
                                     
      RETURN                                                           
      END
