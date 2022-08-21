*$ CREATE SOURCE.FOR
*COPY SOURCE
*
*=== source ===========================================================*
*
      SUBROUTINE SOURCE ( NOMORE )

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1990-2006      by    Alfredo Ferrari & Paola Sala  *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     New source for FLUKA9x-FLUKA200x:                                *
*                                                                      *
*     Created on 07 january 1990   by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 03-mar-06     by    Alfredo Ferrari               *
*                                                                      *
*  This is just an example of a possible user written source routine.  *
*  note that the beam card still has some meaning - in the scoring the *
*  maximum momentum used in deciding the binning is taken from the     *
*  beam momentum.  Other beam card parameters are obsolete.            *
*                                                                      *
*----------------------------------------------------------------------*
*
      INCLUDE '(BEAMCM)'
      INCLUDE '(FHEAVY)'
      INCLUDE '(FLKSTK)'
      INCLUDE '(IOIOCM)'
      INCLUDE '(LTCLCM)'
      INCLUDE '(PAPROP)'
      INCLUDE '(SOURCM)'
      INCLUDE '(SUMCOU)'
*
      LOGICAL LFIRST, trigger
*
      integer ctr,evt,pid(49999)
      real*8 ek(49999),vec(10),delta, xcirc,  extra(4,7)
      real*4 age(49999),
     .       x(49999),y(49999),z(49999),cx(49999),cy(49999),cz(49999)
      save pid,ek,age,x,y,z,cx,cy,cz

      SAVE LFIRST,ctr,delta,MM,m, c45,s45,pi
      DATA LFIRST/.true./,ctr/0/,delta/1d-6/,MM/0/,m/0/, nctr/0/

      common /aeSave/ R2
      common /src2mgd/ vec,extra,evt,iRun,multi

*======================================================================*
*                                                                      *
*                 BASIC VERSION                                        *
*                                                                      *
*======================================================================*
      NOMORE = 0
*  +-------------------------------------------------------------------*
*  |  First call initializations:
      IF ( LFIRST ) THEN
*  |  *** The following 3 cards are mandatory ***
         TKESUM = ZERZER
         LFIRST = .FALSE.
         LUSSRC = .TRUE.
*  |  *** User initialization ***
cc  -------------------------------------------------
        call OAUXFI('evts',81,'unformatted old',ierr)
        if (ierr.ne.0) goto 91
cc  -------------------------------------------------

cc
cc --- new feature  mar'12 05-10:45 ----------
cc
c       nn= 0
c       if (WHASOU(1).gt.0d0) nn= WHASOU(1)
c       call skipAbunch(81,nn)
cc -------------------------------------------

      pi= acos(-1d0)
      c45= cos(45*pi/180)
      s45= sin(45*pi/180)
      END IF
*  |
*  +-------------------------------------------------------------------*
      WEIPRI = ONEONE

    1 continue
      nctr= nctr+1
      nMuon= 0

      read(81,end=99,err=99) iRun,evt,iNcase, n,multi

      if (n.gt.50000) then
        write(*,*) ' ---L--- too many particles ',n
        stop
        endif
     
       write(*,*) ' -o-', evt,nctr,iNcase, n,multi
      call fflush(0)

      do i=1,multi
        read(81) (extra(i,j),j=1,7)
c        write(86,*) (extra(i,j),j=1,7)
         write(*,5) (extra(i,j),j=1,7)
        enddo
    5 format('ex1x7',7f13.6)
      nm= 0
      do j=1,n
        read(81) pid(j),ek(j),age(j),x(j),y(j),z(j),cx(j),cy(j),cz(j)
        
         

cc        xx= x(j)*c45 -y(j)*s45
cc        yy= x(j)*s45 +y(j)*c45
cc     ------------------------- rotate original evts by 45 degree in XY plane
cc        x(j)= xx
cc        y(j)= yy

cc        xx= cx(j)*c45 -cy(j)*s45
cc        yy= cx(j)*s45 +cy(j)*c45
cc     --------------------------- dont forget directions
cc        cx(j)= xx
cc        cy(j)= yy

          nm= nm+1
c        write(*,2) 'muwi',x(j),y(j),z(j), ek(j),age(j),1

          write(78,*) pid(j), x(j),y(j),z(j), ek(j),age(j),1
         
                    

        if (pid(j).eq.10 .or. pid(j).eq.11) then
          nm= nm+1
          write(*,2) 'muwi',x(j),y(j),z(j), ek(j),age(j),1
c Here you print out the muon information in 77 file
          write(77,*) x(j),y(j),z(j), ek(j),age(j),1
          
          if (ek(j).lt.0) ek(j)= -ek(j)
          endif
        enddo

    2 format(a7,3f12.3,2f16.9,i5)

cc ------------------
cc                 ---  muon path length in tank dl
c     trigger= .false.
c     do j=1,n
c       if (pid(j).eq.10 .or. pid(j).eq.11) then
c         call path(x(j),y(j),z(j),cx(j),cy(j),cz(j), dl)
cc       ------
c         trigger= trigger .or. (dl.gt.200 .and. ek(j).gt. 4)
c         if (dl.gt.200 .and. ek(j).gt. 4) then
c           ddl= dl
c           eek= ek(j)
c           endif
c         endif
c       enddo

c     if (trigger) then
c       write(*,*) ' -skid-  dl>200 ',evt,trigger,ddl,eek
c       goto 1
c       endif
cc ------------------

      do i=1,n
        IJBEAM= pid(i)
        write(*,*) ' ---L---  ',i,pid(i)
        if (pid(i).lt.-6) then
          ia= -pid(i)/100
          iz= mod(-pid(i),100)
          IJBEAM= -2
cc       --------------------------------- decode that ion ---
cc                                 NPFLKA == i
          IJHION= iz*1000+ia
          IJHION= IJHION*100+KXHEAV
          IONID= IJHION
          CALL DCDION(IONID)
          CALL SETION(IONID)
          ILOFLK(i)= IJHION
          LRADDC(i)= .FALSE.
          endif

*  +-------------------------------------------------------------------*
*  Push one source particle to the stack. Note that you could as well
*  push many but this way we reserve a maximum amount of space in the
*  stack for the secondaries to be generated

      NPFLKA = NPFLKA + 1
      WTFLK  (NPFLKA) = ONEONE
*  +-------------------------------------------------------------------*
*  |  Heavy ion:
      IF ( IJBEAM .EQ. -2 ) THEN
         IJHION = iz*1000+ia
         IJHION = IJHION * 100 + KXHEAV
         IONID  = IJHION
         CALL DCDION ( IONID )
         CALL SETION ( IONID )
         ILOFLK (NPFLKA) = IJHION
         LRADDC (NPFLKA) = .FALSE.
*  |
*  +-------------------------------------------------------------------*
*  |  Normal hadron:
      ELSE
         IONID = IJBEAM
         ILOFLK (NPFLKA) = IJBEAM
         LRADDC (NPFLKA) = .FALSE.
      END IF
*  |
*  +-------------------------------------------------------------------*
* From this point .....

      LOFLK  (NPFLKA) = 1
      LOUSE  (NPFLKA) = 0
      DO 100 ISPR = 1, MKBMX1
         SPAREK (ISPR,NPFLKA) = ZERZER
 100  CONTINUE
      DO 200 ISPR = 1, MKBMX2
         ISPARK (ISPR,NPFLKA) = 0
 200  CONTINUE
      ISPARK (MKBMX2,NPFLKA) = NPFLKA
      NPARMA = NPARMA + 1
      NUMPAR (NPFLKA) = NPARMA
      NEVENT (NPFLKA) = 0
      DFNEAR (NPFLKA) = +ZERZER
* ... to this point: don't change anything

* Particle age (s)
      agestk(NPFLKA)= age(i)
      AKNSHR (NPFLKA) = -TWOTWO
* Group number for "low" energy neutrons, set to 0 anyway
      IGROUP (NPFLKA) = 0

* Kinetic energy of the particle (GeV) and momentum
      TKEFLK(NPFLKA)= ek(i)
      PMOFLK(NPFLKA)= SQRT(ek(i)*(ek(i)+TWOTWO*AM(IONID)))

cc --- nov'17 11-21:50
cc     so the 14m diameter super CTF was placed such that the top is
cc      at 950.6cm and bottom at approx -450cm  BUT now we have the floor
cc       located at -109.2cm further down
cc     also change 'pull-back' to 78cm rather than the 50cm

* Cosines (tx,ty,tz)
      ccx= cx(i)
      ccy= cy(i)
      ccz= cz(i)
      r= sqrt(ccx**2+ccy**2+ccz**2)
      ccx= ccx/r
      ccy= ccy/r
      ccz= ccz/r
      txflk(NPFLKA)= ccx
      tyflk(NPFLKA)= ccy
      tzflk(NPFLKA)= ccz
* Polarization cosines:
      TXPOL  (NPFLKA) = -TWOTWO
      TYPOL  (NPFLKA) = +ZERZER
      TZPOL  (NPFLKA) = +ZERZER
* Particle coordinates
      xflk(NPFLKA)= x(i) -ccx*78d0
      yflk(NPFLKA)= y(i) -ccy*78d0
      zflk(NPFLKA)= z(i) -ccz*78d0
cc
*  Calculate the total kinetic energy of the primaries: don't change
      IF ( ILOFLK (NPFLKA) .EQ. -2 .OR. ILOFLK (NPFLKA) .GT. 100000 )
     &   THEN
         TKESUM = TKESUM + TKEFLK (NPFLKA) * WTFLK (NPFLKA)
      ELSE IF ( ILOFLK (NPFLKA) .NE. 0 ) THEN
         TKESUM = TKESUM + ( TKEFLK (NPFLKA) + AMDISC (ILOFLK(NPFLKA)) )
     &          * WTFLK (NPFLKA)
      ELSE
         TKESUM = TKESUM + TKEFLK (NPFLKA) * WTFLK (NPFLKA)
      END IF
      RADDLY (NPFLKA) = ZERZER
*  Here we ask for the region number of the hitting point.
*     NREG (NPFLKA) = ...
*  The following line makes the starting region search much more
*  robust if particles are starting very close to a boundary:
      CALL GEOCRS ( TXFLK (NPFLKA), TYFLK (NPFLKA), TZFLK (NPFLKA) )
      CALL GEOREG ( XFLK  (NPFLKA), YFLK  (NPFLKA), ZFLK  (NPFLKA),
     &              NRGFLK(NPFLKA), IDISC )
*  Do not change these cards:
      CALL GEOHSM ( NHSPNT (NPFLKA), 1, -11, MLATTC )
      NLATTC (NPFLKA) = MLATTC
      CMPATH (NPFLKA) = ZERZER

cc     -------------->>>
        enddo

      CALL SOEVSV
      RETURN

   99 continue
      NOMORE= 1
      return

*=== End of subroutine Source =========================================*

   91 continue
      write(*,*) ' --> could not open file!, exit'
      stop


cc   --------------------------------
      entry in2file(iEvt,jNcase,iN)

      write(80) multi,iEvt,jNcase,n,iN
      do i=1,multi
        write(80) (extra(i,j),j=1,7)
        enddo
      do i=1,n
        write(80) pid(i),ek(i),age(i),x(i),y(i),z(i),cx(i),cy(i),cz(i)
        enddo
      return
cc   ------------------

      end

