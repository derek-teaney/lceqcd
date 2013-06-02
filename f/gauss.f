**************************************************
*                                                *
*     GAUSSIAN RANDOM NUMBER GENERATOR           *
*                                                *
*     X IS A VECTOR WITH LENGTH NSITE2 WHICH     *
*     CONTAINS GAUSSIAN RANDOM NUMBERS WITH      *
*     DISTRIBUTION OF EXP( -R*X**2 )             *
*                                                *
**************************************************
      SUBROUTINE GAUSS(X,R,NH)
      IMPLICIT REAL (A-H,O-Z)
      IMPLICIT INTEGER  (I-N)
c      PARAMETER(NS=8,NT=8,NDIM=3,NGRP=4)
c      PARAMETER(NVOL=NT*(NS**(NDIM-1)),NH=NVOL/2)
      DIMENSION W1(NH),W2(NH)
      DIMENSION X(NH),R(NH),RR(NH)
      TPI=2.0E+00*ACOS(-1.)
      DO 100 I=1,NH
 100  RR(I)=1.0E+00/R(I)
      DO 102 I=1,NH
      CALL RANDOM_NUMBER(RX)
 102  W1(I)=RX
      DO 104 IS=1,NH
 104  W1(IS)=-1.0E+00*RR(IS)*LOG(W1(IS))
      DO 106 I=1,NH
      CALL RANDOM_NUMBER(RX)
 106  W2(I)=RX
      DO 108 IS=1,NH
 108  X(IS)=SQRT(W1(IS)) * COS(TPI*W2(IS))
      RETURN
      END
