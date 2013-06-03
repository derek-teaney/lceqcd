!   LUXURY LEVELS.
!  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
!           and Zaman, very long period, but fails many tests.
!  level 1  (p=48): considerable improvement in quality over level 0,
!           now passes the gap test, but still fails spectral test.
!  level 2  (p=97): passes all known tests, but theoretically still
!           defective.
!  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
!           correlations have very small chance of being observed.
!  level 4  (p=389): highest possible luxury, all 24 bits chaotic.
c---------------------------------------------
      SUBROUTINE ranlux(rvec, lenv)
c                    returns a vector RVEC of LEN   
c                   32-bit random floating point numbers between 
c                   zero (not included) and one (also not incl.).
c           Next call to ranlux gives next LEN of random numbers
c---------------------------------------------
      INTEGER  lenv
      REAL        rvec(lenv)
      INTEGER  iseeds(24)
      INCLUDE 'luxury.h'
      DATA ndskip/ 0, 24, 73, 199, 365 /
      DATA i24,j24,luxlev/24,10,lxdflt/
      DATA notyet/.true./
      DATA in24,kount,mkount,carry/0,0,0,0./


!  NOTYET is .TRUE. if no initialization has been performed yet.
!              Default Initialization by Multiplicative Congruential

      IF (notyet) THEN
        notyet = .false.
        jseed = jsdflt
        inseed = jseed
        WRITE (6,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ', jseed
        luxlev = lxdflt
        nskip = ndskip(luxlev)
        lp = nskip + 24
        in24 = 0
        kount = 0
        mkount = 0
      WRITE (6,'(A,I2,A,I4)') ' RANLUX DEFAULT LUXURY LEVEL= ', luxlev,   
     &                            '    p =', lp
      twom24 = 1.
        DO i = 1, 24
          twom24 = twom24 * 0.5
          k = jseed / 53668
          jseed = 40014 * (jseed-k*53668) - k * 12211
          IF (jseed.LT.0) jseed = jseed + icons
          iseeds(i) = MOD(jseed,itwo24)
        END DO
        twom12 = twom24 * 4096.
        DO i = 1, 24
          seeds(i) = REAL(iseeds(i)) * twom24
          next(i) = i - 1
        END DO
        next(1) = 24
        i24 = 24
        j24 = 10
        carry = 0.
        IF (seeds(24).EQ.0.) carry = twom24
      END IF

!          The Generator proper: "Subtract-with-borrow",
!          as proposed by Marsaglia and Zaman,
!          Florida State University, March, 1989

      DO ivec = 1, lenv
        uni = seeds(j24) - seeds(i24) - carry
        IF (uni.LT.0.) THEN
          uni = uni + 1.0
          carry = twom24
        ELSE
          carry = 0.
        END IF
        seeds(i24) = uni
        i24 = next(i24)
        j24 = next(j24)
        rvec(ivec) = uni
      !  small numbers (with less than 12 "significant" bits) are "padded".
        IF (uni.LT.twom12) THEN
          rvec(ivec) = rvec(ivec) + twom24 * seeds(j24)
      !        and zero is forbidden in case someone takes a logarithm
          IF (rvec(ivec).EQ.0.) rvec(ivec) = twom24 * twom24
        END IF
      !        Skipping to luxury.  As proposed by Martin Luscher.
        in24 = in24 + 1
        IF (in24.EQ.24) THEN
          in24 = 0
          kount = kount + nskip
          DO isk = 1, nskip
            uni = seeds(j24) - seeds(i24) - carry
            IF (uni.LT.0.) THEN
              uni = uni + 1.0
              carry = twom24
            ELSE
              carry = 0.
            END IF
            seeds(i24) = uni
            i24 = next(i24)
            j24 = next(j24)
          END DO
        END IF
      END DO
      kount = kount + lenv
      IF (kount.GE.igiga) THEN
        mkount = mkount + 1
        kount = kount - igiga
      END IF
      RETURN

      END ! SUBROUTINE ranlux

c---------------------------------------------
c           Subroutine to input and float integer seeds from previous run
      SUBROUTINE rluxin(isdext)
c                   restarts the generator from vector  
c                   ISDEXT of 25 32-bit integers 
c---------------------------------------------
      INCLUDE 'luxury.h'
      INTEGER isdext(25)
      IF (notyet) THEN
        WRITE (6,'(A)')' Proper results ONLY with ',
     &        'initialisation from 25 integers obtained with RLUXUT'
        notyet = .false.
      END IF
      
      twom24 = 1.
      DO i = 1, 24
        next(i) = i - 1
        twom24 = twom24 * 0.5
      END DO
      next(1) = 24
      twom12 = twom24 * 4096.
c      WRITE (6,'(A)') ' FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:'
c      WRITE (6,'(5X,5I12)') isdext
      DO i = 1, 24
        seeds(i) = REAL(isdext(i)) * twom24
      END DO
      carry = 0.
      IF (isdext(25).LT.0) carry = twom24
      isd = ABS(isdext(25))
      i24 = MOD(isd,100)
      isd = isd / 100
      j24 = MOD(isd,100)
      isd = isd / 100
      in24 = MOD(isd,100)
      isd = isd / 100
      luxlev = isd
      IF (luxlev.LE.maxlev) THEN
        nskip = ndskip(luxlev)
c        WRITE (6,'(A,I2)') ' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ',
c     &                       luxlev
      ELSE IF (luxlev.GE.24) THEN
        nskip = luxlev - 24
c        WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXIN TO:', luxlev
      ELSE
        nskip = ndskip(maxlev)
c        WRITE (6,'(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ', luxlev
        luxlev = maxlev
      END IF
      inseed = -1
      RETURN
      
      END ! SUBROUTINE rluxin

c---------------------------------------------
c                    Subroutine to ouput seeds as integers
      SUBROUTINE rluxut(isdext)
c      outputs the current values of the 25
c                 32-bit integer seeds, to be used for restarting
c      ISDEXT must be dimensioned 25 in the calling program 
c---------------------------------------------
      INCLUDE 'luxury.h' 
      INTEGER isdext(25)
      DO i = 1, 24
        isdext(i) = INT(seeds(i)*twop12*twop12)
      END DO
      isdext(25) = i24 + 100 * j24 + 10000 * in24 + 1000000 * luxlev
      IF (carry.GT.0.) isdext(25) = -isdext(25)
      RETURN
      
      END ! SUBROUTINE rluxut

c---------------------------------------------
!                    Subroutine to output the "convenient" restart point
      SUBROUTINE rluxat(lout, inout, k1, k2)
c---------------------------------------------
      INCLUDE 'luxury.h'
      
      lout = luxlev
      inout = inseed
      k1 = kount
      k2 = mkount
      RETURN
      
      END ! SUBROUTINE rluxat

c---------------------------------------------
c                    Subroutine to initialize from one or three integers
      SUBROUTINE rluxgo(lux, ins, k1, k2)
c                initializes the generator from  
c               one 32-bit integer INT and sets Luxury Level LUX
c               which is integer between zero and MAXLEV, or if
c               LUX .GT. 24, it sets p=LUX directly.  K1 and K2
c               should be set to zero unless restarting at a break
c               point given by output of RLUXAT
c---------------------------------------------
      INCLUDE 'luxury.h'
      INTEGER  iseeds(24)
      
      IF (lux.LT.0) THEN
        luxlev = lxdflt
      ELSE IF (lux.LE.maxlev) THEN
        luxlev = lux
      ELSE IF (lux.LT.24.OR.lux.GT.2000) THEN
        luxlev = maxlev
        WRITE (6,'(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ', lux
      ELSE
        luxlev = lux
        DO ilx = 0, maxlev
          IF (lux.EQ.ndskip(ilx)+24) luxlev = ilx
        END DO
      END IF
      IF (luxlev.LE.maxlev) THEN
        nskip = ndskip(luxlev)
        WRITE (6,'(A,I2,A,I4)') ' RANLUX LUXURY LEVEL SET BY RLUXGO :',
     &              luxlev,'     P=', nskip + 24
      ELSE
        nskip = luxlev - 24
        WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXGO TO:', luxlev
      END IF
      in24 = 0
      IF (ins.LT.0) WRITE (6,'(A)')
     &        ' Illegal initialization by RLUXGO, negative input seed'
      IF (ins.GT.0) THEN
        jseed = ins
        WRITE (6,'(A,3I12)')' RANLUX INITIALIZED BY RLUXGO FROM SEEDS',
     &                                    jseed, k1, k2
      ELSE
        jseed = jsdflt
       WRITE (6,'(A)')' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED'
      END IF
      inseed = jseed
      notyet = .false.
      twom24 = 1.
      DO i = 1, 24
        twom24 = twom24 * 0.5
        k = jseed / 53668
        jseed = 40014 * (jseed-k*53668) - k * 12211
        IF (jseed.LT.0) jseed = jseed + icons
        iseeds(i) = MOD(jseed,itwo24)
      END DO
      twom12 = twom24 * 4096.
      DO i = 1, 24
        seeds(i) = REAL(iseeds(i)) * twom24
        next(i) = i - 1
      END DO
      next(1) = 24
      i24 = 24
      j24 = 10
      carry = 0.
      IF (seeds(24).EQ.0.) carry = twom24
      !        If restarting at a break point, skip K1 + IGIGA*K2
      !        Note that this is the number of numbers delivered to
      !        the user PLUS the number skipped (if luxury .GT. 0).
      kount = k1
      mkount = k2
      IF (k1+k2.NE.0) THEN
        DO iouter = 1, k2 + 1
          inner = igiga
          IF (iouter.EQ.k2+1) inner = k1
          DO isk = 1, inner
            uni = seeds(j24) - seeds(i24) - carry
            IF (uni.LT.0.) THEN
              uni = uni + 1.0
              carry = twom24
            ELSE
              carry = 0.
            END IF
            seeds(i24) = uni
            i24 = next(i24)
            j24 = next(j24)
          END DO
        END DO
      !         Get the right value of IN24 by direct calculation
        in24 = MOD(kount,nskip+24)
        IF (mkount.GT.0) THEN
          izip = MOD(igiga, nskip+24)
          izip2 = mkount * izip + in24
          in24 = MOD(izip2, nskip+24)
        END IF
      !       Now IN24 had better be between zero and 23 inclusive
        IF (in24.GT.23) THEN
          WRITE (6,'(A/A,3I11,A,I5)') 
     &       '  Error in RESTARTING with RLUXGO:', '  The values', ins, 
     &                k1, k2, ' cannot occur at luxury level', luxlev
          in24 = 0
        END IF
      END IF
      RETURN
      
      END ! SUBROUTINE rluxgo

