      SUBROUTINE ML5_0_1_1_HELAS_CALLS_AMPB_1(P,NHEL,H,IC)
C     
C     Modules
C     
      USE ML5_0_1_1_POLYNOMIAL_CONSTANTS
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
      PARAMETER (NLOOPS=40, NLOOPGROUPS=40, NCTAMPS=14)
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=54)
      INTEGER    NWAVEFUNCS,NLOOPWAVEFUNCS
      PARAMETER (NWAVEFUNCS=17,NLOOPWAVEFUNCS=80)
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
      COMMON/ML5_0_1_1_FILTERS/GOODAMP,GOODHEL,HELOFFSET

      LOGICAL CHECKPHASE
      LOGICAL HELDOUBLECHECKED
      COMMON/ML5_0_1_1_INIT/CHECKPHASE, HELDOUBLECHECKED

      INTEGER SQSO_TARGET
      COMMON/ML5_0_1_1_SOCHOICE/SQSO_TARGET

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/ML5_0_1_1_SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE
     $ ,CT_REQ_SO_DONE,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE
     $ ,MP_LOOP_REQ_SO_DONE,CTCALL_REQ_SO_DONE,FILTER_SO

      INTEGER I_SO
      COMMON/ML5_0_1_1_I_SO/I_SO
      INTEGER I_LIB
      COMMON/ML5_0_1_1_I_LIB/I_LIB

      COMPLEX*16 W(20,NWAVEFUNCS)
      COMMON/ML5_0_1_1_W/W

      COMPLEX*16 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE,
     $ -1:NLOOPWAVEFUNCS)
      COMPLEX*16 PL(0:3,-1:NLOOPWAVEFUNCS)
      COMMON/ML5_0_1_1_WL/WL,PL

      COMPLEX*16 AMPL(3,NLOOPAMPS)
      COMMON/ML5_0_1_1_AMPL/AMPL

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
      CALL OXXXXX(P(0,4),ZERO,NHEL(4),+1*IC(4),W(1,4))
      CALL IXXXXX(P(0,5),ZERO,NHEL(5),-1*IC(5),W(1,5))
      CALL VVV1P0_1(W(1,1),W(1,2),GC_4,ZERO,ZERO,W(1,6))
      CALL FFV1P0_3(W(1,5),W(1,4),GC_5,ZERO,ZERO,W(1,7))
C     Counter-term amplitude(s) for loop diagram number 1
      CALL VVS1_0(W(1,6),W(1,7),W(1,3),R2_GGHB,AMPL(1,1))
      CALL FFV1_1(W(1,4),W(1,1),GC_5,ZERO,ZERO,W(1,8))
      CALL FFV1P0_3(W(1,5),W(1,8),GC_5,ZERO,ZERO,W(1,9))
C     Counter-term amplitude(s) for loop diagram number 3
      CALL VVS1_0(W(1,2),W(1,9),W(1,3),R2_GGHB,AMPL(1,2))
      CALL FFV1_2(W(1,5),W(1,1),GC_5,ZERO,ZERO,W(1,10))
      CALL FFV1P0_3(W(1,10),W(1,4),GC_5,ZERO,ZERO,W(1,11))
C     Counter-term amplitude(s) for loop diagram number 5
      CALL VVS1_0(W(1,2),W(1,11),W(1,3),R2_GGHB,AMPL(1,3))
      CALL FFV1_1(W(1,4),W(1,2),GC_5,ZERO,ZERO,W(1,12))
      CALL FFV1P0_3(W(1,5),W(1,12),GC_5,ZERO,ZERO,W(1,13))
C     Counter-term amplitude(s) for loop diagram number 7
      CALL VVS1_0(W(1,1),W(1,13),W(1,3),R2_GGHB,AMPL(1,4))
      CALL FFV1_2(W(1,5),W(1,2),GC_5,ZERO,ZERO,W(1,14))
      CALL FFV1P0_3(W(1,14),W(1,4),GC_5,ZERO,ZERO,W(1,15))
C     Counter-term amplitude(s) for loop diagram number 9
      CALL VVS1_0(W(1,1),W(1,15),W(1,3),R2_GGHB,AMPL(1,5))
      CALL VVV1P0_1(W(1,2),W(1,7),GC_4,ZERO,ZERO,W(1,16))
C     Counter-term amplitude(s) for loop diagram number 14
      CALL VVS1_0(W(1,1),W(1,16),W(1,3),R2_GGHB,AMPL(1,6))
      CALL VVV1P0_1(W(1,1),W(1,7),GC_4,ZERO,ZERO,W(1,17))
C     Counter-term amplitude(s) for loop diagram number 19
      CALL VVS1_0(W(1,2),W(1,17),W(1,3),R2_GGHB,AMPL(1,7))
C     Counter-term amplitude(s) for loop diagram number 21
      CALL VVS1_0(W(1,6),W(1,7),W(1,3),R2_GGHT,AMPL(1,8))
C     Counter-term amplitude(s) for loop diagram number 23
      CALL VVS1_0(W(1,2),W(1,9),W(1,3),R2_GGHT,AMPL(1,9))
C     Counter-term amplitude(s) for loop diagram number 25
      CALL VVS1_0(W(1,2),W(1,11),W(1,3),R2_GGHT,AMPL(1,10))
C     Counter-term amplitude(s) for loop diagram number 27
      CALL VVS1_0(W(1,1),W(1,13),W(1,3),R2_GGHT,AMPL(1,11))
C     Counter-term amplitude(s) for loop diagram number 29
      CALL VVS1_0(W(1,1),W(1,15),W(1,3),R2_GGHT,AMPL(1,12))
C     Counter-term amplitude(s) for loop diagram number 34
      CALL VVS1_0(W(1,1),W(1,16),W(1,3),R2_GGHT,AMPL(1,13))
C     Counter-term amplitude(s) for loop diagram number 39
      CALL VVS1_0(W(1,2),W(1,17),W(1,3),R2_GGHT,AMPL(1,14))
C     At this point, all CT amps needed for (QCD=8), i.e. of split
C      order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 2000

      GOTO 1001
 2000 CONTINUE
      CT_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

