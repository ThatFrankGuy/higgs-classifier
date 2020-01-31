ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP1()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      REAL*16 MP__PI, MP__ZERO
      PARAMETER (MP__PI=3.1415926535897932384626433832795E0_16)
      PARAMETER (MP__ZERO=0E0_16)
      INCLUDE 'mp_input.inc'
      INCLUDE 'mp_coupl.inc'

      GC_22 = (MDL_CW*MDL_EE*MDL_COMPLEXI)/(2.000000D+00*MDL_SW)
      MP__GC_22 = (MP__MDL_CW*MP__MDL_EE*MP__MDL_COMPLEXI)/(2.000000E
     $ +00_16*MP__MDL_SW)
      GC_23 = -(MDL_EE*MDL_COMPLEXI*MDL_SW)/(6.000000D+00*MDL_CW)
      MP__GC_23 = -(MP__MDL_EE*MP__MDL_COMPLEXI*MP__MDL_SW)/(6.000000E
     $ +00_16*MP__MDL_CW)
      GC_31 = (MDL_EE__EXP__2*MDL_COMPLEXI*MDL_V)/(2.000000D+00
     $ *MDL_SW__EXP__2)
      MP__GC_31 = (MP__MDL_EE__EXP__2*MP__MDL_COMPLEXI*MP__MDL_V)
     $ /(2.000000E+00_16*MP__MDL_SW__EXP__2)
      GC_32 = MDL_EE__EXP__2*MDL_COMPLEXI*MDL_V+(MDL_CW__EXP__2
     $ *MDL_EE__EXP__2*MDL_COMPLEXI*MDL_V)/(2.000000D+00
     $ *MDL_SW__EXP__2)+(MDL_EE__EXP__2*MDL_COMPLEXI*MDL_SW__EXP__2
     $ *MDL_V)/(2.000000D+00*MDL_CW__EXP__2)
      MP__GC_32 = MP__MDL_EE__EXP__2*MP__MDL_COMPLEXI*MP__MDL_V
     $ +(MP__MDL_CW__EXP__2*MP__MDL_EE__EXP__2*MP__MDL_COMPLEXI
     $ *MP__MDL_V)/(2.000000E+00_16*MP__MDL_SW__EXP__2)
     $ +(MP__MDL_EE__EXP__2*MP__MDL_COMPLEXI*MP__MDL_SW__EXP__2
     $ *MP__MDL_V)/(2.000000E+00_16*MP__MDL_CW__EXP__2)
      GC_47 = (MDL_EE*MDL_COMPLEXI*MDL_CONJG__CKM33)/(MDL_SW
     $ *MDL_SQRT__2)
      MP__GC_47 = (MP__MDL_EE*MP__MDL_COMPLEXI*MP__MDL_CONJG__CKM33)
     $ /(MP__MDL_SW*MP__MDL_SQRT__2)
      END
