      subroutine MGDRAW (ICODE,MREG)

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
      INCLUDE '(CASLIM)'
      INCLUDE '(COMPUT)'
      INCLUDE '(SOURCM)'
      INCLUDE '(FHEAVY)'
      INCLUDE '(FLKSTK)'
      INCLUDE '(GENSTK)'
      INCLUDE '(MGDDCM)'
      INCLUDE '(PAPROP)'
      INCLUDE '(QUEMGD)'
      INCLUDE '(SUMCOU)'
      INCLUDE '(TRACKR)'


      LOGICAL Lfirst,keepIt,firstOne, skip
      SAVE Lfirst,keepIt,firstOne
      DATA Lfirst /.true./

      real*8 vec(10),extra(4,7)
      common /aeStop/ skip
      common /src2mgd/ vec,extra,evt,iRun,multi
      integer n,evt,pid(92:94,19999),nD(92:94),reg(92:94,19999)
      real*8 ek(92:94,19999)
      real*4 age(92:94,19999)
      real*4 x(92:94,19999),y(92:94,19999),z(92:94,19999)
      real*4 cx(92:94,19999),cy(92:94,19999),cz(92:94,19999)
      integer arr(0:999),brr(19999), krr(99)
      real*4 atime(19999)
      logical first,pmts
      save nrCTF,nrLSV,first,pmts,arr,brr,atime, Kctr,krr
      data first/.true./
      
      character*8 regnam(3)
      integer nregs(3)
cc --------------------------------------------------------------------
      data regnam/'LAr_Bath','Inner_Ar','SensVol'/

      save pid,reg,ek,age,x,y,z,cx,cy,cz, nD,Ictr,Jctr,
     .     nLAr_Bath,nInner_Ar,nSensVol
cc ---------------------------------------------------------------------

     
      DIMENSION DTQUEN ( MXTRCK, MAXQMG ), AMHELP ( MFSTCK )
*
      CHARACTER*20 FILNAM
      LOGICAL LFCOPE
      SAVE LFCOPE
      DATA LFCOPE / .FALSE. /
*
*----------------------------------------------------------------------*
*                                                                      *
*     Icode = 1: call from Kaskad                                      *
*     Icode = 2: call from Emfsco                                      *
*     Icode = 3: call from Kasneu                                      *
*     Icode = 4: call from Kashea                                      *
*     Icode = 5: call from Kasoph                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
      IF ( .NOT. LFCOPE ) THEN
         LFCOPE = .TRUE.
         IF ( KOMPUT .EQ. 2 ) THEN
            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
         ELSE
            FILNAM = CFDRAW
         END IF
         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
     &          'UNFORMATTED' )
      END IF
      WRITE (IODRAW) NTRACK, MTRACK, JTRACK, SNGL (ETRACK),
     &               SNGL (WTRACK)
      WRITE (IODRAW) ( SNGL (XTRACK (I)), SNGL (YTRACK (I)),
     &                SNGL (ZTRACK (I)), I = 0, NTRACK ),
     &               ( SNGL (DTRACK (I)), I = 1, MTRACK ),
     &                 SNGL (CTRACK)

      RETURN



cc    Boundary-(X)crossing DRAWing
cc   ------------------------------
      ENTRY BXDRAW ( ICODE, MREG, NEWREG, XSCO, YSCO, ZSCO )
cc --- nov'11 14:14 (munich)
cc     try to quickly adapt to the latest Darkside deign idea
cc     make use of (atmospheric) liquid argon as (active) veto shield

      if (NEWREG.eq.nLAr_Bath .and. NEWREG.gt.MREG) then
         nUnit= 92
      else if (NEWREG.eq.nInner_Ar .and. NEWREG.gt.MREG) then
         nUnit= 93
      else if (NEWREG.eq.nSensVol .and. NEWREG.gt.MREG) then
         nUnit= 94
         e= etrack-am(JTRACK)
c         write(*,10) MREG,JTRACK,e,atrack,XSCO,YSCO,ZSCO,Ltrack
c             write(80,*) mreg,jtrack,e,atrack,XSCO,YSCO,ZSCO,Ltrack
      else
         return
      endif
 10   format('  -c- ',2i6,2f16.9,3f11.4,i5)

cc   if we get here we got a trigger
      
      if (JTRACK.ge.-6 .and. JTRACK.lt.68) then
         myTrack= JTRACK
         if (JTRACK.eq.-2) then
            myTrack= -ibarch(-2)*100-ichrge(-2)
            write(0,*) ' -->>>-- -2 ',myTrack
         endif



         if (nD(nUnit).ge.19999) then
            write(*,*) ' -x-> nD(nUnit) >19999 ',NCASE,nUnit
         else
            nD(nUnit)= nD(nUnit)+1
            n= nD(nUnit)
            reg(nUnit,n)= mReg
            pid(nUnit,n)= myTrac
            ek(nUnit,n)= ETRACK-AM(jTrack)
            age(nUnit,n)= aTrack
            x(nUnit,n)= xsco
            y(nUnit,n)= ysco
            z(nUnit,n)= zsco
            cx(nUnit,n)= cXtrck
            cy(nUnit,n)= cYtrck
            cz(nUnit,n)= cZtrck
         endif
      else
         write(*,*) ' ===> trouble jTrack ',NCASE
      endif

      RETURN

cc    Event End DRAWing
cc   -------------------
      ENTRY EEDRAW ( ICODE )

      write(*,*) ' -xx- ',NCASE,evt
      do k=92,94
         if (nD(k).gt.0) then
            write(*,*) ' -on- ',k,NCASE,nD(k),evt,nOrig,Ictr,Jctr,skip
c            do i=1,nD(k)
c               write(*,11) pid(k,i),ek(k,i),age(k,i),x(k,i),y(k,i),z(k,i)
c            enddo
         endif
cc For geant4 run, track events as entering LAr veto, otherwise backscattering particles make evtbin energy deposition in TPC less reliable to count the number of scatterings
         if (k.eq.94) then
            write(k,*) NCASE,nD(k),evt, nOrig,nSelec,Ictr
            do i=1,nD(k)
               write(k,*) reg(k,i),pid(k,i),ek(k,i),age(k,i),
     .                     x(k,i),y(k,i),z(k,i),cx(k,i),cy(k,i),cz(k,i)
               write(k,*)NCASE,evt, reg(k,i), pid(k,i),ek(k,i),age(k,i),
     .                     x(k,i),y(k,i),z(k,i),cx(k,i),cy(k,i),cz(k,i)
            enddo
         endif
         
         if (k.eq.93) then
            write(k,*) NCASE,nD(k),evt, nOrig,nSelec,Ictr
            do i=1,nD(k)
c               write(k,*) reg(k,i),pid(k,i),ek(k,i),age(k,i),
c     .                     x(k,i),y(k,i),z(k,i),cx(k,i),cy(k,i),cz(k,i)
               write(k,*)NCASE,evt, reg(k,i), pid(k,i),ek(k,i),age(k,i),
     .                     x(k,i),y(k,i),z(k,i),cx(k,i),cy(k,i),cz(k,i)
            enddo

          endif
        enddo
   11 format('  -in- ',i6,2e16.9,3f11.4)
   
c     if (nD(94).gt.0) call record(ncase)
c      write(*,*) '--ncase--',NCASE
      write(*,*) ' -end- ',NCASE
      RETURN



cc    ENergy deposition DRAWing
cc   ------------------------------------------------------
      ENTRY ENDRAW ( ICODE, MREG, RULL, XSCO, YSCO, ZSCO)
*
      IF ( .NOT. LFCOPE ) THEN
         LFCOPE = .TRUE.
         IF ( KOMPUT .EQ. 2 ) THEN
            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
         ELSE
            FILNAM = CFDRAW
         END IF
         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
     &          'UNFORMATTED' )
      END IF
      WRITE (IODRAW)  0, ICODE, JTRACK, SNGL (ETRACK), SNGL (WTRACK)
cc    CuWall_updated_to_Ti_wall_(Titanium vessel)
      IF (MREG ==92) THEN
         WRITE (53,*)  MREG, NEVENT(NPFLKA), NCASE, ATRACK,
     .               ILOFLK (NPFLKA), SNGL (XSCO),
     .               SNGL (YSCO), SNGL (ZSCO), SNGL (RULL)
      ENDIF
cc    LAr_bath
      IF (MREG ==93) THEN
         WRITE (54,*)  MREG, NEVENT(NPFLKA), NCASE, ATRACK,
     .               ILOFLK (NPFLKA), SNGL (XSCO),
     .               SNGL (YSCO), SNGL (ZSCO), SNGL (RULL)
      ENDIF
cc   Inner Argon
      IF (MREG ==94) THEN
         WRITE (55,*)  MREG, NEVENT(NPFLKA), NCASE, ATRACK,
     .               ILOFLK (NPFLKA), SNGL (XSCO),
     .               SNGL (YSCO), SNGL (ZSCO), SNGL (RULL)
      ENDIF
cc    TPCWall
      IF (MREG ==95) THEN
         WRITE (56,*)  MREG, NEVENT(NPFLKA), NCASE, ATRACK,
     .               ILOFLK (NPFLKA), SNGL (XSCO),
     .               SNGL (YSCO), SNGL (ZSCO), SNGL (RULL)
      ENDIF
cc    PMMA top_TPC
      IF (MREG ==96) THEN
         WRITE (57,*)  MREG, NEVENT(NPFLKA), NCASE, ATRACK,
     .               ILOFLK (NPFLKA), SNGL (XSCO),
     .               SNGL (YSCO), SNGL (ZSCO), SNGL (RULL)
      ENDIF
cc   Sensitive volume inside TPC
      IF (MREG ==97) THEN
         WRITE (58,*)  MREG, NEVENT(NPFLKA), NCASE, ATRACK,
     .               ILOFLK (NPFLKA), SNGL (XSCO),
     .               SNGL (YSCO), SNGL (ZSCO), SNGL (RULL)
      ENDIF
cc   Gas_pocket
      IF (MREG ==98) THEN
         WRITE (59,*)  MREG, NEVENT(NPFLKA), NCASE, ATRACK,
     .               ILOFLK (NPFLKA), SNGL (XSCO),
     .               SNGL (YSCO), SNGL (ZSCO), SNGL (RULL)
      ENDIF
cc   PMMA bottom_TPC
      IF (MREG ==99) THEN
         WRITE (60,*)  MREG, NEVENT(NPFLKA), NCASE, ATRACK,
     .               ILOFLK (NPFLKA), SNGL (XSCO),
     .               SNGL (YSCO), SNGL (ZSCO), SNGL (RULL)
      ENDIF




      RETURN



cc    SOurce particle DRAWing
cc   -------------------------
      ENTRY SODRAW
      skip= .false.
cc -------------------- pick up regions by name here ---
      if (first) then
         first= .false.
         write(*,*) ' -r-   --selected regions--'
         do i=1,3
            call GEON2R(regnam(i),nregs(i),ierr)
            write(*,*) ' -r- ',i,' ',regnam(i),nregs(i),ierr
         enddo
         nLAr_Bath= nregs(1)
         nInner_Ar= nregs(2)
         nSensVol= nregs(3)
      endif
      do i=92,94
         nD(i)= 0
      enddo
      RETURN


cc    USer dependent DRAWing
cc   -------------------------
      ENTRY USDRAW ( ICODE, MREG, XSCO, YSCO, ZSCO)
      RETURN

      END
