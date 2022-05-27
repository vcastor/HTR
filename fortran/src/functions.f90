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

    ENDFUNCTION f_distance
!***********************************************************************
    FUNCTION f_distance(xyz, uvw) RESULT (distance)
    IMPLICIT NONE
    REAL (KIND=8)               :: distance
    REAL (KIND=8), DIMENSION(3) :: xyz, uwv

    distance = DSQRT(f_distance(xyz, uvw))

    ENDFUNCTION
!***********************************************************************
    FUNCTION f_coulomb(q,Z,distance) RESULT (coulomb)
    INTEGER       :: q, Z
    REAL (KIND=8) :: coulomb, distance

    coulomb = q*Z/distance

    ENDFUNCTION f_coulomb
!***********************************************************************
    FUNCTION f_cero(x) RESULT (y)
    IMPLICIT NONE
    INCLUDE 'parameters.h'
    REAL (KIND=8), INTENT(IN)  :: x
    REAL (KIND=8), INTENT(OUT) :: y

    y  = 0.5d0*DSQRT(pi/x)*DERF(DSQRT(x))

    ENDFUNCTION f_cero
!***********************************************************************
