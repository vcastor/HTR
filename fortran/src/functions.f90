!***********************************************************************
!********************   HONI SOIT QUI MAL Y PENSE   ********************
!********************        NOLI ME TANGERE        ********************
!***********************************************************************
    FUNCTION f_distance(xyz, uvw) RESULT (distance)
    IMPLICIT NONE
    REAL (KIND=8) :: distance, dx, dy, dz
    REAL (KIND=8), DIMENSION(3) :: xyz, uvw

    dx = xyz(1) - uvw(1)
    dy = xyz(2) - uvw(2)
    dz = xyz(3) - uvw(3)

    distance = DSQRT(dx*dx + dy*dy + dz*dz)
    ENDFUNCTION f_distance
!***********************************************************************
    FUNCTION f_coulomb(q,Z,distance) RESULT (coulomb)
    INTEGER :: q, Z
    REAL (KIND=8) :: coulomb, distance

    coulomb = q*Z/distance
    ENDFUNCTION f_coulomb
!***********************************************************************
    FUNCTION f_cero(x) RESULT (y)
    IMPLICIT NONE
    REAL (KIND=8) :: x, y, pi

    pi = DACOS(-1.d0)
    y  = 0.5d0*DSQRT(pi/x)*DERF(DSQRT(x))
    ENDFUNCTION f_cero

