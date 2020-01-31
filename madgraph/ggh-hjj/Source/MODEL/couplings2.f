ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP2()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      R2_GGHB = 4.000000D+00*(-((MDL_COMPLEXI*MDL_YB)/MDL_SQRT__2))
     $ *(1.000000D+00/2.000000D+00)*(MDL_G__EXP__2/(8.000000D+00*PI**2)
     $ )*MDL_MB
      R2_GGHT = 4.000000D+00*(-((MDL_COMPLEXI*MDL_YT)/MDL_SQRT__2))
     $ *(1.000000D+00/2.000000D+00)*(MDL_G__EXP__2/(8.000000D+00*PI**2)
     $ )*MDL_MT
      GC_4 = -G
      GC_5 = MDL_COMPLEXI*G
      GC_6 = MDL_COMPLEXI*MDL_G__EXP__2
      END
