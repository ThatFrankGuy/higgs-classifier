      SUBROUTINE MP_HELAS_CALLS_AMPB_1(P,NHEL,H,IC)
C     
      USE POLYNOMIAL_CONSTANTS
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=5)
      INTEGER    NCOMB
      PARAMETER (NCOMB=24)

      INTEGER NBORNAMPS
      PARAMETER (NBORNAMPS=2)
      INTEGER    NLOOPS, NLOOPGROUPS, NCTAMPS
      PARAMETER (NLOOPS=19, NLOOPGROUPS=11, NCTAMPS=32)
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=51)
      INTEGER    NWAVEFUNCS,NLOOPWAVEFUNCS
      PARAMETER (NWAVEFUNCS=14,NLOOPWAVEFUNCS=43)
      REAL*16     ZERO
      PARAMETER (ZERO=0.0E0_16)
      COMPLEX*32     IZERO
      PARAMETER (IZERO=CMPLX(0.0E0_16,0.0E0_16,KIND=16))
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=1, NSQUAREDSO=1, NAMPSO=2)
C     
C     ARGUMENTS
C     
      REAL*16 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      INTEGER H
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J,K
      COMPLEX*32 COEFS(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'mp_coupl_same_name.inc'

      INTEGER GOODHEL(NCOMB)
      LOGICAL GOODAMP(NSQUAREDSO,NLOOPGROUPS)
      COMMON/FILTERS/GOODAMP,GOODHEL

      INTEGER SQSO_TARGET
      COMMON/SOCHOICE/SQSO_TARGET

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE
     $ ,CT_REQ_SO_DONE,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE
     $ ,MP_LOOP_REQ_SO_DONE,CTCALL_REQ_SO_DONE,FILTER_SO

      COMPLEX*32 AMP(NBORNAMPS)
      COMMON/MP_AMPS/AMP
      COMPLEX*32 W(20,NWAVEFUNCS)
      COMMON/MP_W/W

      COMPLEX*32 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE,
     $ -1:NLOOPWAVEFUNCS)
      COMPLEX*32 PL(0:3,-1:NLOOPWAVEFUNCS)
      COMMON/MP_WL/WL,PL

      COMPLEX*32 AMPL(3,NCTAMPS)
      COMMON/MP_AMPL/AMPL

C     
C     ----------
C     BEGIN CODE
C     ----------

C     The target squared split order contribution is already reached
C      if true.
      IF (FILTER_SO.AND.MP_CT_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

      CALL MP_VXXXXX(P(0,1),ZERO,NHEL(1),-1*IC(1),W(1,1))
      CALL MP_OXXXXX(P(0,2),ZERO,NHEL(2),-1*IC(2),W(1,2))
      CALL MP_SXXXXX(P(0,3),+1*IC(3),W(1,3))
      CALL MP_VXXXXX(P(0,4),MDL_MW,NHEL(4),+1*IC(4),W(1,4))
      CALL MP_IXXXXX(P(0,5),ZERO,NHEL(5),-1*IC(5),W(1,5))
      CALL MP_FFV1_1(W(1,2),W(1,1),GC_5,ZERO,ZERO,W(1,6))
      CALL MP_VVS1_1(W(1,4),W(1,3),GC_31,MDL_MW,MDL_WW,W(1,7))
C     Amplitude(s) for born diagram with ID 1
      CALL MP_FFV2_0(W(1,5),W(1,6),W(1,7),GC_47,AMP(1))
      CALL MP_FFV1_2(W(1,5),W(1,1),GC_5,ZERO,ZERO,W(1,8))
C     Amplitude(s) for born diagram with ID 2
      CALL MP_FFV2_0(W(1,8),W(1,2),W(1,7),GC_47,AMP(2))
      CALL MP_FFV2_2(W(1,5),W(1,7),GC_47,ZERO,ZERO,W(1,9))
C     Counter-term amplitude(s) for loop diagram number 3
      CALL MP_R2_QQ_1_0(W(1,9),W(1,6),R2_QQQ,AMPL(1,1))
C     Counter-term amplitude(s) for loop diagram number 4
      CALL MP_FFV2_0(W(1,5),W(1,6),W(1,7),R2_SXCW,AMPL(1,2))
C     Counter-term amplitude(s) for loop diagram number 5
      CALL MP_FFV2_0(W(1,8),W(1,2),W(1,7),R2_SXCW,AMPL(1,3))
C     Counter-term amplitude(s) for loop diagram number 6
      CALL MP_FFV1_0(W(1,9),W(1,2),W(1,1),UV_GQQQ_1EPS,AMPL(2,4))
      CALL MP_FFV1_0(W(1,9),W(1,2),W(1,1),UV_GQQQ_1EPS,AMPL(2,5))
      CALL MP_FFV1_0(W(1,9),W(1,2),W(1,1),UV_GQQQ_1EPS,AMPL(2,6))
      CALL MP_FFV1_0(W(1,9),W(1,2),W(1,1),UV_GQQQ_1EPS,AMPL(2,7))
      CALL MP_FFV1_0(W(1,9),W(1,2),W(1,1),UV_GQQB,AMPL(1,8))
      CALL MP_FFV1_0(W(1,9),W(1,2),W(1,1),UV_GQQQ_1EPS,AMPL(2,9))
      CALL MP_FFV1_0(W(1,9),W(1,2),W(1,1),UV_GQQT,AMPL(1,10))
      CALL MP_FFV1_0(W(1,9),W(1,2),W(1,1),UV_GQQQ_1EPS,AMPL(2,11))
      CALL MP_FFV1_0(W(1,9),W(1,2),W(1,1),UV_GQQG_1EPS,AMPL(2,12))
      CALL MP_FFV1_0(W(1,9),W(1,2),W(1,1),R2_GQQ,AMPL(1,13))
      CALL MP_FFV2_1(W(1,2),W(1,7),GC_47,ZERO,ZERO,W(1,10))
C     Counter-term amplitude(s) for loop diagram number 11
      CALL MP_R2_QQ_1_0(W(1,8),W(1,10),R2_QQQ,AMPL(1,14))
C     Counter-term amplitude(s) for loop diagram number 12
      CALL MP_FFV1_0(W(1,5),W(1,10),W(1,1),UV_GQQQ_1EPS,AMPL(2,15))
      CALL MP_FFV1_0(W(1,5),W(1,10),W(1,1),UV_GQQQ_1EPS,AMPL(2,16))
      CALL MP_FFV1_0(W(1,5),W(1,10),W(1,1),UV_GQQQ_1EPS,AMPL(2,17))
      CALL MP_FFV1_0(W(1,5),W(1,10),W(1,1),UV_GQQQ_1EPS,AMPL(2,18))
      CALL MP_FFV1_0(W(1,5),W(1,10),W(1,1),UV_GQQB,AMPL(1,19))
      CALL MP_FFV1_0(W(1,5),W(1,10),W(1,1),UV_GQQQ_1EPS,AMPL(2,20))
      CALL MP_FFV1_0(W(1,5),W(1,10),W(1,1),UV_GQQT,AMPL(1,21))
      CALL MP_FFV1_0(W(1,5),W(1,10),W(1,1),UV_GQQQ_1EPS,AMPL(2,22))
      CALL MP_FFV1_0(W(1,5),W(1,10),W(1,1),UV_GQQG_1EPS,AMPL(2,23))
      CALL MP_FFV1_0(W(1,5),W(1,10),W(1,1),R2_GQQ,AMPL(1,24))
      CALL MP_FFV2_1(W(1,2),W(1,4),GC_47,ZERO,ZERO,W(1,11))
      CALL MP_FFV1P0_3(W(1,5),W(1,11),GC_5,ZERO,ZERO,W(1,12))
C     Counter-term amplitude(s) for loop diagram number 14
      CALL MP_VVS1_0(W(1,1),W(1,12),W(1,3),R2_GGHB,AMPL(1,25))
      CALL MP_FFV2_2(W(1,5),W(1,4),GC_47,ZERO,ZERO,W(1,13))
      CALL MP_FFV1P0_3(W(1,13),W(1,2),GC_5,ZERO,ZERO,W(1,14))
C     Counter-term amplitude(s) for loop diagram number 16
      CALL MP_VVS1_0(W(1,1),W(1,14),W(1,3),R2_GGHB,AMPL(1,26))
C     Counter-term amplitude(s) for loop diagram number 18
      CALL MP_VVS1_0(W(1,1),W(1,12),W(1,3),R2_GGHT,AMPL(1,27))
C     Counter-term amplitude(s) for loop diagram number 20
      CALL MP_VVS1_0(W(1,1),W(1,14),W(1,3),R2_GGHT,AMPL(1,28))
C     At this point, all CT amps needed for (QCD=4), i.e. of split
C      order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 2000

      GOTO 1001
 2000 CONTINUE
      MP_CT_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

