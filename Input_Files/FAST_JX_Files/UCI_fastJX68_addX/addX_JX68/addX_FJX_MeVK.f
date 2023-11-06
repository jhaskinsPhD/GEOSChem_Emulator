c----addX_FJX.f  -- with user supplied subroutine that supplies X-section x Q-yield
C----     generate fast-JX 18-bin X-sections
c---------revised and updated for (mprather,2/2013)
      implicit none
      integer, parameter :: NB_ = 100
      integer, parameter :: NS_ = 40000
      integer, parameter :: NZ_ = 13550
      integer, parameter :: NJ_ = 18

      real*8 SRB(15,NJ_)
      real*8, dimension(NB_+1) :: WBIN
      real*8, dimension(NB_) :: FBIN, ABIN
      real*8, dimension(NJ_) :: FFBIN,AABIN
      integer IJX(NB_), ITT
      integer NB,I,J,J1,J2,K,K1,K2
      integer INIT
      real*8 W(NS_),F(NS_)
      integer IBINJ(NS_)
      real*8 WZ(NZ_),X(NZ_,3)
      real*8 W1,W2, TT,XP,XM, WW,XNEW
      character*6 TITLNEW
      character*4 TITLET

      open (1, file='wavel-bins.dat', status='OLD')
        SRB(:,:) = 0.d0
        read(1,'(i5)') NB
        if (NB .gt. NB_) stop
        read(1,'(5x,f8.3)') (WBIN(I), I=1,NB+1)
        read(1,*)
        read(1,*)
        read(1,'(2x,15f5.1)') ((SRB(I,J),I=1,15),J=1,8)
        read(1,*)
        read(1,'(5x,i5)') (IJX(I),I=16,NB)
      close (1)

      open (2, file='solar-p05nm-UCI.dat', status='OLD')
        read(2,*)
        read(2,*)
        read(2,'(f10.4,e10.3)') (W(J),F(J), J=1,NS_)
      close (2)

      open (3, file='XO3-p05nm-UCI.dat', status='OLD')
        read(3,*)
        read(3,*)
        read(3,'(f10.4,3e10.3)') (WZ(J),X(J,1),X(J,2),X(J,3), J=1,NZ_)
      close (3)

!!!!!!!!!!!!!!!!!!!!initialization call to user subroutine!!!!!!!!!!!!!!
        INIT = 0
      call X_MeVK (WW,TT,XP, XM, X,INIT,TITLNEW)



c---synchronize with the O3 cross sections (whether done or not)
c---will loop K=1:NZ_  (J = K - 1 + K1), wavel = WZ(K)
        do J=1,NS_
          if (WZ(1) .eq. W(J)) goto 10
        enddo
          write(6,*) ' cannot synch the solar & Xsections'
          stop
   10   K1 = J
        K2 = min( NS_, NZ_+K1-1)
c          write(6,'(a,2f10.3,i10)') ' synch:',WZ(1),W(K1),K1
c          write(6,'(a,2f10.3,i10)') ' synch:',WZ(NZ_),W(K2),K2


c---now assign bin #(I=1:77) to each p05nm microbin J (1:40000)
        IBINJ(:) = 0
      do I=1,NB
         W1 = WBIN(I)
         W2 = WBIN(I+1)
        do J=1,NS_
          if (W(J) .gt. W1) goto 11
        enddo
          J = NS_ + 1
   11     J1 = J
       do J=J1,NS_
          if (W(J) .gt. W2) goto 12
        enddo
          J = NS_ + 1
   12     J2 = J-1
        do J=J1,J2
          IBINJ(J) = I
        enddo
c          write(6,'(i5,2f9.3,2i9,2f9.3)') I, W1,W2, J1,J2,W(J1),W(J2)
      enddo


c>>>>be aware that this binning does not interpolate and is OK for large bins
c          it has 7% error in the very short wavel S-R bins of pratmo.
c>>>>it should be fine for weighting cross sections!

!  Total Q-yld from Stern-Volmer converted (as in acetone) to 3 different trop levels:
!         alt=   0 km      5 km      13 km
!           M=  2.46E+19   1.50E+19   5.8E+18
!           T=     295K      272K       220K
!           P =    999 hPa   566 hPa    177 hPa
!   Only apply a single overall quantum yield PHI
!        PHI = exp[ -0.055 * (w-308) ] /( 5.5 + 9.2e-19*[M])
!   Thus there are 3 tables for MeVK, each designated by T(K)


c---this looping is set for tropospheric VOCs - may need to change for other Js
c---nominally the interpolation is only over Temperature, user sets here based on
c           ranges of measured Xsections.  new fast-JX can use 1, 2 or 3 sets.

!!!!!!!!!!!!!!!!!!!major temperature-density loop !!!!!!!!!!!!!!!!!!!!!!
      do ITT =1,3

        if (ITT.eq.1) then
          TITLET = 'p999'
c--      NB - for some VOCs do pressure interpolation of X-section
c--          if doing T interp, then TITLET = ' 295'
c            XM or XP may be used in the X_xxx call for Stern-Volmer
          TT = 295.d0
          XP = 999.d0
          XM = 2.46d19
        else if (ITT.eq.2) then
          TITLET = 'p566'
          TT = 272.d0
          XP = 566.d0
          XM = 1.506d19
        else
          TITLET = 'p177'
          TT = 220.d0
          XP = 177.d0
          XM = 0.58d19
        endif
c---now ready to do any flux-weighted means over the 77-pratmo bins
         FBIN(:) = 0.d0
         ABIN(:) = 0.d0

!!!!!!!!!!!!!!!!!!!primary high-resolution wavelength loop!!!!!!!!!!!!!!
       do J = K1,K2
        K = J - K1 + 1
        I = IBINJ(J)
        if (I .gt. 0) then

        call X_MeVK (W(J), TT,XP,XM, XNEW, INIT, TITLNEW)

          FBIN(I) = FBIN(I) + F(J)
          ABIN(I) = ABIN(I) + F(J)*XNEW
        endif
       enddo
       do I=1,NB
        if (FBIN(I) .gt. 0.d0) ABIN(I) = ABIN(I)/FBIN(I)
       enddo
c---write out UCI std 77-bin data
c       write(6,'(a6,a4/(1p,8e10.3))') TITLNEW,TITLET, ABIN

!!!!!!!!!!!!!!!!!!!secondary sum 77-bin pratmo ==> 18-bin fast-JX!!!!!!!
c---combine fast-JX bins:
c---    non-SR bands (16:NB) are assigned a single JX bin
c---    SR bands are split (by Opacity Distrib Fn) into a range of JX bins
        FFBIN(:) = 0.d0
        AABIN(:) = 0.d0
      do I=16,NB
        J = IJX(I)
        FFBIN(J) = FFBIN(J) + FBIN(I)
        AABIN(J) = AABIN(J) + FBIN(I)*ABIN(I)
      enddo
      do I=1,15
        do J=1,NJ_
          FFBIN(J) = FFBIN(J) + FBIN(I)*SRB(I,J)
          AABIN(J) = AABIN(J) + FBIN(I)*ABIN(I)*SRB(I,J)
        enddo
      enddo
      do J=1,NJ_
        if (FFBIN(J) .gt. 0.d0) AABIN(J) = AABIN(J)/FFBIN(J)
      enddo


!!!!!!!!!!!!!!!!!!!save UCI fast-JX data bins!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(6,'(a6,a4,1p,6e10.3/10x,6e10.3/10x,6e10.3)')
     &     TITLNEW, TITLET, AABIN

      enddo
      stop
      end




c-------------sample subroutine for fast-JX Xsection generation---------
c-----------------------------------------------------------------------
      subroutine X_MeVK (WW,TT,PP,MM, XXWT,INIT,TITLNEW)
c-----------------------------------------------------------------------
c   WW = wavelength (nm)
c   TT = temerature (K) for interpolation
c   PP = pressure (hPa) can be used in Stern-Volmer formula if need be
c   MM = air density (#/cm3), ditto
c   XXWT = cross section (cm2) as a function of WW and TT (and PP, MM)
c   INIT = initialization:
c     if INIT.eq.0 .then reads in any tables and sets Xsect name to TITLNEW

      implicit none
      real*8, intent(in) :: WW, TT,PP,MM
      integer, intent(inout) :: INIT
      real*8, intent(out) :: XXWT
      character*6, intent(out) :: TITLNEW
      character*80 FTBL,TABLE,FORMW
      real*8, dimension(999) :: W,XW
      real*8  TTL,XXW,FW, XWI,XWIP1, TFACT, T3(3),XT3(999,3)
      integer NW,N, I,IW,IT
      save  W,XW,T3,XT3,NW
      TITLNEW = 'MeVK  '
      FTBL = 'XMeVK_JPL10V.dat'

!  Total Q-yld from Stern-Volmer converted (as in acetone) to 3 different trop levels:
!         alt=   0 km      5 km      13 km
!           M=  2.46E+19   1.50E+19   5.8E+18
!           T=     295K      272K       220K
!           P =    999 hPa   566 hPa    177 hPa
!   Only apply a single overall quantum yield PHI
!        PHI = exp[ -0.055 * (w-308) ] /( 5.5 + 9.2e-19*[M])
!   Thus there are 3 tables for MeVK, each designated by T(K)
      if (INIT .eq. 0) then
        open (3, file=FTBL, status='OLD')
          read(3,'(a80)') TABLE
           write(6,'(2a/a)') ' openfile=',FTBL, TABLE
          read(3,*)
          read(3,*)
          read(3,'(12x,3f8.0)') (T3(I),I=3,1,-1)
c         write(6,'(A,3F8.0)') ' Temperature:',(T3(I),I=3,1,-1)
          read(3,'(i4,1x,a)') NW,FORMW
        do N = 1,NW
          read(3,FORMW) W(N),XW(N),(XT3(N,I),I=3,1,-1)
c         write(6,'(I3,5F12.5)') int(W(N)),XW(N),(XT3(N,I),I=3,1,-1)
        enddo
        close(3)
        INIT = 1
      else
        TTL = min(T3(3), max(T3(1), TT))
        if (TT .gt. T3(2)) then
          IT = 2
        else
          IT = 1
        endif
        TFACT = (TTL - T3(IT))/(T3(IT+1) - T3(IT))
c---interpolate X-section vs. Wavelength, T=298
        IW = 1
        do I = 2,NW-1
          if (WW .gt. W(I)) IW = I
        enddo
        FW = (WW - W(IW))/(W(IW+1) - W(IW))
        FW = min(1.d0, max(0.d0, FW))
c---T dependence: set for MeVK
        XWI   = XW(IW) * (XT3(IW,IT)+TFACT*(XT3(IW,IT+1)-XT3(IW,IT)))
        XWIP1 = XW(IW+1) *
     &          (XT3(IW+1,IT)+TFACT*(XT3(IW+1,IT+1)-XT3(IW+1,IT)))
        XXW = XWI + FW*(XWIP1-XWI)

        XXWT = XXW  * 1.d-20
      endif
      return
      end
