REAL*8 FUNCTION BBFN(X)
C     ********************************************************************
 C    *  This subroutine calculates the fractional blackbody             *
C     *  emissive power f(n*lambda*T), where X=n*lambda*T in (micro-m*K) *
C     ********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PI=4.D0*DATAN(1.D0)
      CC=1.5D1/PI**4
      C2=1.4388D4
      EPS=1.D-16

      V=C2/X
      EX=DEXP(V)

      M=0
      BBFN=0.D0
      EM=1.D0
      M=M+1
      VM=M*V
      BM=(6.D0+VM*(6.D0+VM*(3.D0+VM)))/FLOAT(M)**4
      EM=EM/EX
      BBFN=BBFN+BM*EM
      IF(VM**3*EM.GT.EPS) GOTO 5
      BBFN=CC*BBFN
      RETURN
      END

