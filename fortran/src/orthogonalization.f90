!-----------------------------------------------------------------------
! This SUBROUTINE computes an orthogonalization for the problem:       !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~  B = A^(-1/2)  ~~~~~~~~~~~~~~~~~~~~~~~~~!
! n    :: size of matrices                                             !
! A    :: input matrix                                                 !
! B    :: output matrix                                                !
! info :: info flag                                                    !
!-----------------------------------------------------------------------
SUBROUTINE ORTHOGONALIZATION_M(n, A, B)
IMPLICIT NONE
INTEGER                       :: i, n, info
REAL (KIND=8), DIMENSION(n)   :: vecaux
REAL (KIND=8), DIMENSION(n,n) :: A, B, mataux

mataux = A; info = 0

CALL DIAG_M(n, mataux, vecaux)

DO i = 1, n
  IF (vecaux(i) .GT. 1.d-6) THEN
    vecaux(i) = 1.d0/DSQRT(vecaux(i))
  ELSE
    info = 301                                          !too small value
  ENDIF
ENDDO

CALL WIZARD(info)

CALL AdiAtrans(n, mataux, vecaux, B)

ENDSUBROUTINE ORTHOGONALIZATION_M
