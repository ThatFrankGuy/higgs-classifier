      SUBROUTINE ML5_0_0_1_HELAS_CALLS_AMPB_1(P,NHEL,H,IC)
C     
C     Modules
C     
      USE ML5_0_0_1_POLYNOMIAL_CONSTANTS
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=5)
      INTEGER    NCOMB
      PARAMETER (NCOMB=16)
      INTEGER    NLOOPS, NLOOPGROUPS, NCTAMPS
      PARAMETER (NLOOPS=228, NLOOPGROUPS=228, NCTAMPS=54)
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=282)
      INTEGER    NWAVEFUNCS,NLOOPWAVEFUNCS
      PARAMETER (NWAVEFUNCS=35,NLOOPWAVEFUNCS=420)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      REAL*16     MP__ZERO
      PARAMETER (MP__ZERO=0.0E0_16)
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=1, NSQUAREDSO=1, NAMPSO=1)
C     
C     ARGUMENTS
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      INTEGER H
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J,K
      COMPLEX*16 COEFS(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)

      LOGICAL DUMMYFALSE
      DATA DUMMYFALSE/.FALSE./
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
      INCLUDE 'mp_coupl.inc'

      INTEGER HELOFFSET
      INTEGER GOODHEL(NCOMB)
      LOGICAL GOODAMP(NSQUAREDSO,NLOOPGROUPS)
      COMMON/ML5_0_0_1_FILTERS/GOODAMP,GOODHEL,HELOFFSET

      LOGICAL CHECKPHASE
      LOGICAL HELDOUBLECHECKED
      COMMON/ML5_0_0_1_INIT/CHECKPHASE, HELDOUBLECHECKED

      INTEGER SQSO_TARGET
      COMMON/ML5_0_0_1_SOCHOICE/SQSO_TARGET

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/ML5_0_0_1_SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE
     $ ,CT_REQ_SO_DONE,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE
     $ ,MP_LOOP_REQ_SO_DONE,CTCALL_REQ_SO_DONE,FILTER_SO

      INTEGER I_SO
      COMMON/ML5_0_0_1_I_SO/I_SO
      INTEGER I_LIB
      COMMON/ML5_0_0_1_I_LIB/I_LIB

      COMPLEX*16 W(20,NWAVEFUNCS)
      COMMON/ML5_0_0_1_W/W

      COMPLEX*16 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE,
     $ -1:NLOOPWAVEFUNCS)
      COMPLEX*16 PL(0:3,-1:NLOOPWAVEFUNCS)
      COMMON/ML5_0_0_1_WL/WL,PL

      COMPLEX*16 AMPL(3,NLOOPAMPS)
      COMMON/ML5_0_0_1_AMPL/AMPL

C     
C     ----------
C     BEGIN CODE
C     ----------

C     The target squared split order contribution is already reached
C      if true.
      IF (FILTER_SO.AND.CT_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

      CALL VXXXXX(P(0,1),ZERO,NHEL(1),-1*IC(1),W(1,1))
      CALL VXXXXX(P(0,2),ZERO,NHEL(2),-1*IC(2),W(1,2))
      CALL SXXXXX(P(0,3),+1*IC(3),W(1,3))
      CALL VXXXXX(P(0,4),ZERO,NHEL(4),+1*IC(4),W(1,4))
      CALL VXXXXX(P(0,5),ZERO,NHEL(5),+1*IC(5),W(1,5))
      CALL VVV1P0_1(W(1,1),W(1,2),GC_4,ZERO,ZERO,W(1,6))
      CALL VVV1P0_1(W(1,4),W(1,6),GC_4,ZERO,ZERO,W(1,7))
C     Counter-term amplitude(s) for loop diagram number 1
      CALL VVS1_0(W(1,5),W(1,7),W(1,3),R2_GGHB,AMPL(1,1))
      CALL VVV1P0_1(W(1,5),W(1,6),GC_4,ZERO,ZERO,W(1,8))
C     Counter-term amplitude(s) for loop diagram number 2
      CALL VVS1_0(W(1,4),W(1,8),W(1,3),R2_GGHB,AMPL(1,2))
      CALL VVV1P0_1(W(1,4),W(1,5),GC_4,ZERO,ZERO,W(1,9))
C     Counter-term amplitude(s) for loop diagram number 5
      CALL VVS1_0(W(1,6),W(1,9),W(1,3),R2_GGHB,AMPL(1,3))
      CALL VVV1P0_1(W(1,1),W(1,4),GC_4,ZERO,ZERO,W(1,10))
      CALL VVV1P0_1(W(1,2),W(1,5),GC_4,ZERO,ZERO,W(1,11))
C     Counter-term amplitude(s) for loop diagram number 13
      CALL VVS1_0(W(1,10),W(1,11),W(1,3),R2_GGHB,AMPL(1,4))
      CALL VVV1P0_1(W(1,5),W(1,10),GC_4,ZERO,ZERO,W(1,12))
C     Counter-term amplitude(s) for loop diagram number 15
      CALL VVS1_0(W(1,2),W(1,12),W(1,3),R2_GGHB,AMPL(1,5))
      CALL VVV1P0_1(W(1,2),W(1,10),GC_4,ZERO,ZERO,W(1,13))
C     Counter-term amplitude(s) for loop diagram number 23
      CALL VVS1_0(W(1,5),W(1,13),W(1,3),R2_GGHB,AMPL(1,6))
      CALL VVV1P0_1(W(1,2),W(1,4),GC_4,ZERO,ZERO,W(1,14))
      CALL VVV1P0_1(W(1,1),W(1,5),GC_4,ZERO,ZERO,W(1,15))
C     Counter-term amplitude(s) for loop diagram number 25
      CALL VVS1_0(W(1,14),W(1,15),W(1,3),R2_GGHB,AMPL(1,7))
      CALL VVV1P0_1(W(1,4),W(1,15),GC_4,ZERO,ZERO,W(1,16))
C     Counter-term amplitude(s) for loop diagram number 27
      CALL VVS1_0(W(1,2),W(1,16),W(1,3),R2_GGHB,AMPL(1,8))
      CALL VVV1P0_1(W(1,2),W(1,15),GC_4,ZERO,ZERO,W(1,17))
C     Counter-term amplitude(s) for loop diagram number 35
      CALL VVS1_0(W(1,4),W(1,17),W(1,3),R2_GGHB,AMPL(1,9))
      CALL VVV1P0_1(W(1,5),W(1,14),GC_4,ZERO,ZERO,W(1,18))
C     Counter-term amplitude(s) for loop diagram number 37
      CALL VVS1_0(W(1,1),W(1,18),W(1,3),R2_GGHB,AMPL(1,10))
      CALL VVV1P0_1(W(1,4),W(1,11),GC_4,ZERO,ZERO,W(1,19))
C     Counter-term amplitude(s) for loop diagram number 45
      CALL VVS1_0(W(1,1),W(1,19),W(1,3),R2_GGHB,AMPL(1,11))
      CALL VVV1P0_1(W(1,2),W(1,9),GC_4,ZERO,ZERO,W(1,20))
C     Counter-term amplitude(s) for loop diagram number 68
      CALL VVS1_0(W(1,1),W(1,20),W(1,3),R2_GGHB,AMPL(1,12))
      CALL VVVV1P0_1(W(1,2),W(1,4),W(1,5),GC_6,ZERO,ZERO,W(1,21))
      CALL VVVV3P0_1(W(1,2),W(1,4),W(1,5),GC_6,ZERO,ZERO,W(1,22))
      CALL VVVV4P0_1(W(1,2),W(1,4),W(1,5),GC_6,ZERO,ZERO,W(1,23))
C     Counter-term amplitude(s) for loop diagram number 85
      CALL VVS1_0(W(1,1),W(1,21),W(1,3),R2_GGHB,AMPL(1,13))
      CALL VVS1_0(W(1,1),W(1,22),W(1,3),R2_GGHB,AMPL(1,14))
      CALL VVS1_0(W(1,1),W(1,23),W(1,3),R2_GGHB,AMPL(1,15))
      CALL VVV1P0_1(W(1,1),W(1,14),GC_4,ZERO,ZERO,W(1,24))
C     Counter-term amplitude(s) for loop diagram number 87
      CALL VVS1_0(W(1,5),W(1,24),W(1,3),R2_GGHB,AMPL(1,16))
      CALL VVV1P0_1(W(1,1),W(1,11),GC_4,ZERO,ZERO,W(1,25))
C     Counter-term amplitude(s) for loop diagram number 89
      CALL VVS1_0(W(1,4),W(1,25),W(1,3),R2_GGHB,AMPL(1,17))
      CALL VVV1P0_1(W(1,1),W(1,9),GC_4,ZERO,ZERO,W(1,26))
C     Counter-term amplitude(s) for loop diagram number 91
      CALL VVS1_0(W(1,2),W(1,26),W(1,3),R2_GGHB,AMPL(1,18))
      CALL VVVV1P0_1(W(1,1),W(1,2),W(1,4),GC_6,ZERO,ZERO,W(1,27))
      CALL VVVV3P0_1(W(1,1),W(1,2),W(1,4),GC_6,ZERO,ZERO,W(1,28))
      CALL VVVV4P0_1(W(1,1),W(1,2),W(1,4),GC_6,ZERO,ZERO,W(1,29))
C     Counter-term amplitude(s) for loop diagram number 93
      CALL VVS1_0(W(1,5),W(1,27),W(1,3),R2_GGHB,AMPL(1,19))
      CALL VVS1_0(W(1,5),W(1,28),W(1,3),R2_GGHB,AMPL(1,20))
      CALL VVS1_0(W(1,5),W(1,29),W(1,3),R2_GGHB,AMPL(1,21))
      CALL VVVV1P0_1(W(1,1),W(1,2),W(1,5),GC_6,ZERO,ZERO,W(1,30))
      CALL VVVV3P0_1(W(1,1),W(1,2),W(1,5),GC_6,ZERO,ZERO,W(1,31))
      CALL VVVV4P0_1(W(1,1),W(1,2),W(1,5),GC_6,ZERO,ZERO,W(1,32))
C     Counter-term amplitude(s) for loop diagram number 95
      CALL VVS1_0(W(1,4),W(1,30),W(1,3),R2_GGHB,AMPL(1,22))
      CALL VVS1_0(W(1,4),W(1,31),W(1,3),R2_GGHB,AMPL(1,23))
      CALL VVS1_0(W(1,4),W(1,32),W(1,3),R2_GGHB,AMPL(1,24))
      CALL VVVV1P0_1(W(1,1),W(1,4),W(1,5),GC_6,ZERO,ZERO,W(1,33))
      CALL VVVV3P0_1(W(1,1),W(1,4),W(1,5),GC_6,ZERO,ZERO,W(1,34))
      CALL VVVV4P0_1(W(1,1),W(1,4),W(1,5),GC_6,ZERO,ZERO,W(1,35))
C     Counter-term amplitude(s) for loop diagram number 97
      CALL VVS1_0(W(1,2),W(1,33),W(1,3),R2_GGHB,AMPL(1,25))
      CALL VVS1_0(W(1,2),W(1,34),W(1,3),R2_GGHB,AMPL(1,26))
      CALL VVS1_0(W(1,2),W(1,35),W(1,3),R2_GGHB,AMPL(1,27))
C     Counter-term amplitude(s) for loop diagram number 99
      CALL VVS1_0(W(1,5),W(1,7),W(1,3),R2_GGHT,AMPL(1,28))
C     Counter-term amplitude(s) for loop diagram number 100
      CALL VVS1_0(W(1,4),W(1,8),W(1,3),R2_GGHT,AMPL(1,29))
C     Counter-term amplitude(s) for loop diagram number 103
      CALL VVS1_0(W(1,6),W(1,9),W(1,3),R2_GGHT,AMPL(1,30))
C     Counter-term amplitude(s) for loop diagram number 111
      CALL VVS1_0(W(1,10),W(1,11),W(1,3),R2_GGHT,AMPL(1,31))
C     Counter-term amplitude(s) for loop diagram number 113
      CALL VVS1_0(W(1,2),W(1,12),W(1,3),R2_GGHT,AMPL(1,32))
C     Counter-term amplitude(s) for loop diagram number 121
      CALL VVS1_0(W(1,5),W(1,13),W(1,3),R2_GGHT,AMPL(1,33))
C     Counter-term amplitude(s) for loop diagram number 123
      CALL VVS1_0(W(1,14),W(1,15),W(1,3),R2_GGHT,AMPL(1,34))
C     Counter-term amplitude(s) for loop diagram number 125
      CALL VVS1_0(W(1,2),W(1,16),W(1,3),R2_GGHT,AMPL(1,35))
C     Counter-term amplitude(s) for loop diagram number 133
      CALL VVS1_0(W(1,4),W(1,17),W(1,3),R2_GGHT,AMPL(1,36))
C     Counter-term amplitude(s) for loop diagram number 135
      CALL VVS1_0(W(1,1),W(1,18),W(1,3),R2_GGHT,AMPL(1,37))
C     Counter-term amplitude(s) for loop diagram number 143
      CALL VVS1_0(W(1,1),W(1,19),W(1,3),R2_GGHT,AMPL(1,38))
C     Counter-term amplitude(s) for loop diagram number 166
      CALL VVS1_0(W(1,1),W(1,20),W(1,3),R2_GGHT,AMPL(1,39))
C     Counter-term amplitude(s) for loop diagram number 183
      CALL VVS1_0(W(1,1),W(1,21),W(1,3),R2_GGHT,AMPL(1,40))
      CALL VVS1_0(W(1,1),W(1,22),W(1,3),R2_GGHT,AMPL(1,41))
      CALL VVS1_0(W(1,1),W(1,23),W(1,3),R2_GGHT,AMPL(1,42))
C     Counter-term amplitude(s) for loop diagram number 185
      CALL VVS1_0(W(1,5),W(1,24),W(1,3),R2_GGHT,AMPL(1,43))
C     Counter-term amplitude(s) for loop diagram number 187
      CALL VVS1_0(W(1,4),W(1,25),W(1,3),R2_GGHT,AMPL(1,44))
C     Counter-term amplitude(s) for loop diagram number 189
      CALL VVS1_0(W(1,2),W(1,26),W(1,3),R2_GGHT,AMPL(1,45))
C     Counter-term amplitude(s) for loop diagram number 191
      CALL VVS1_0(W(1,5),W(1,27),W(1,3),R2_GGHT,AMPL(1,46))
      CALL VVS1_0(W(1,5),W(1,28),W(1,3),R2_GGHT,AMPL(1,47))
      CALL VVS1_0(W(1,5),W(1,29),W(1,3),R2_GGHT,AMPL(1,48))
C     Counter-term amplitude(s) for loop diagram number 193
      CALL VVS1_0(W(1,4),W(1,30),W(1,3),R2_GGHT,AMPL(1,49))
      CALL VVS1_0(W(1,4),W(1,31),W(1,3),R2_GGHT,AMPL(1,50))
      CALL VVS1_0(W(1,4),W(1,32),W(1,3),R2_GGHT,AMPL(1,51))
C     Counter-term amplitude(s) for loop diagram number 195
      CALL VVS1_0(W(1,2),W(1,33),W(1,3),R2_GGHT,AMPL(1,52))
      CALL VVS1_0(W(1,2),W(1,34),W(1,3),R2_GGHT,AMPL(1,53))
      CALL VVS1_0(W(1,2),W(1,35),W(1,3),R2_GGHT,AMPL(1,54))
C     At this point, all CT amps needed for (QCD=8), i.e. of split
C      order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 2000

      GOTO 1001
 2000 CONTINUE
      CT_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

