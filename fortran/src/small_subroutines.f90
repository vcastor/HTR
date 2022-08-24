!-----------------------------------------------------------------------
! This SUBROUTINE computes a matrix B as the result of the product     !
! between a matrix, a diagonal matrix (saved as a vector) and the      !
! transpose of the first one 🏳️‍⚧️ .                               !
!                              B = A di A^T                            !
! |     |    |     | |.d...| |     |^T                                 !
! |  B  |  = |  A  | |..d..| |  A  |                                   !
! |     |    |     | |...d.| |     |                                   !
!-----------------------------------------------------------------------
SUBROUTINE AdiAtrans(n, A, di, B)
IMPLICIT NONE
INTEGER                      :: n, i, j, k
REAL(KIND=8), DIMENSION(n)   :: di
REAL(KIND=8), DIMENSION(n,n) :: A, B

B(:,:) = 0.d0
DO i = 1, n
  DO j = 1, n
    DO k = 1, n
      B(i,k) = B(i,k) + A(i,j)*di(j)*A(k,j)
    ENDDO
  ENDDO
ENDDO

ENDSUBROUTINE AdiAtrans
