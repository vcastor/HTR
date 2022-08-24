!***********************************************************************
!********************   HONI SOIT QUI MAL Y PENSE   ********************
!********************        NOLI ME TANGERE        ********************
!***********************************************************************
    FUNCTION f_distance2(xyz, uvw) RESULT (distance2)
    IMPLICIT NONE
    REAL (KIND=8)               :: distance2, dx, dy, dz
    REAL (KIND=8), DIMENSION(3) :: xyz, uvw

    dx = xyz(1) - uvw(1)
    dy = xyz(2) - uvw(2)
    dz = xyz(3) - uvw(3)

    distance2 = dx*dx + dy*dy + dz*dz

    ENDFUNCTION f_distance2
!***********************************************************************
    FUNCTION f_distance(xyz, uvw) RESULT (distance)
    IMPLICIT NONE
    REAL (KIND=8)               :: distance
    REAL (KIND=8), EXTERNAL     :: f_distance2
    REAL (KIND=8), DIMENSION(3) :: xyz, uvw

    distance = DSQRT(f_distance2(xyz, uvw))

    ENDFUNCTION
!***********************************************************************
    FUNCTION f_coulomb(q,Z,distance) RESULT (coulomb)
    INTEGER       :: q, Z
    REAL (KIND=8) :: coulomb, distance

    coulomb = q*Z/distance

    ENDFUNCTION
!***********************************************************************
    FUNCTION f_trace_mat(n, A) RESULT (tr)
    IMPLICIT NONE
    INTEGER                      :: i, n
    REAL(KIND=8)                 :: tr
    REAL(KIND=8), DIMENSION(n,n) :: A

    tr = 0.d0
    DO i = 1, n
      tr = tr + A(i,i)
    ENDDO

    ENDFUNCTION
!***********************************************************************
    FUNCTION f_cero(x) RESULT (y)
    IMPLICIT NONE
    INCLUDE 'parameters.h'
    REAL (KIND=8), INTENT(IN) :: x
    REAL (KIND=8)             :: y

    IF (x .GT. 1.D-10) THEN                                 !Aproximation
      y = 0.5d0*DSQRT(pi/x)*DERF(DSQRT(x))
    ELSE                                                    !Analytically
      y = 1.d0
    ENDIF

    ENDFUNCTION f_cero
!***********************************************************************
!    FUNCTION f_n(n, x) RESULT (y)
!    IMPLICIT NONE
!    INCLUDE 'parameters.h'
!    REAL (KIND=8), INTENT(IN) :: x, n
!    REAL (KIND=8)             :: y
!    REAL (KIND=8), EXTERNAL   :: DGAMI
!
!    y = GAMMA(n+0.5) - DGAMI(n+0.5,x)
!    y = y/(2.d0*x**(n+0.5))
!
!    ENDFUNCTION f_n
!***********************************************************************
!    FUNCTION DGAMI(a, x) RESULT (y)
!    IMPLICIT NONE
!    INCLUDE 'paramters.h'
!    REAL (KIND=8), INTENT(IN) :: a, x
!    REAL (KIND=8)             :: y
!    RAEL (KIND=8), EXTERNAL   :: DGAMIT
!
!    IF (a .LE. 0.d0).OR.(x .LE. 0.d0) THEN
!      y = 0.d0
!    ELSE
!      IF x .GT. 0.d0 THEN
!        y = DEXP(a*DLOG(x) + alngam(x)) * gamit(a, x)
!      ENDIF
!    ENDIF
!!***********************************************************************
!    FUNCTION DGAMIT (a, x) RESULT (y)
!    IMPLICIT NONE
!    REAL, INTENT(IN) :: a, x
!    REAL             :: aeps, ainta, algap1, alng, alx, h, sga, sgngam
!    REAL             :: t, alneps, sqeps, bot
!
!    alneps = -DLOG(3.4E+38)   !the largest magnitude
!    sqeps  =  DSQRT(1.2E-07)  !the largest relative spacing
!    bot    =  DLOG(1.2E-38)   !the smallest positive magnitude
!
!    alx = DLOG(x)
!    sga = DSIGN (1.0, a)
!
!    ainta = AINT (a+0.5*sga)
!    aeps  = a - ainta
!
!    IF ( x .LE. 1.0 ) THEN
!      IF ( a  .GE. -0.5 .OR. aeps .NE. 0.0) THEN
!        CALL algams (a+1.0, algap1, sgngam)
!      ENDIF
!      gamit = r9gmit (a, x, algap1, sgngam, alx)
!    ELSE
!      IF ( a .GE. x ) THEN
!        t = r9lgit (a, x, alngam(a+1.0))
!        gamit = DEXP(t)
!      ELSE
!        alng = r9lgic (a, x, alx)
!        ! evaluate gamit in terms of log(gamic(a,x))
!        h = 1.0
!        IF ( aeps .NE. 0.0 .or. ainta .GT. 0.0 ) THEN
!          CALL algams (a+1.0, algap1, sgngam)
!          t = DLOG(a) + alng - algap1
!          IF ( t .LE. alneps ) THEN
!            IF ( t .GT. -alneps ) THEN
!              h = 1.0 - sga*sgngam*DEXP(t)
!            ENDIF
!            t     = -a*alx + log(abs(h))
!            gamit = sign (exp(t), h)
!          ELSE
!          t     = t - a*alx
!          gamit = -sga*sgngam*exp(t)
!          ENDIF
!        ENDIF
!      ENDIF
!    ENDIF
!    ENDFUNCTION gamit
!!***********************************************************************
