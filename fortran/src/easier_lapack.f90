!-----------------------------------------------------------------------
!********************   HONI SOIT QUI MAL Y PENSE   ********************
!********************        NOLI ME TANGERE        ********************
!-----------------------------------------------------------------------
!***********************************************************************
! The next SUBROUTINES are no more than the LAPACK SUBROUTINES that    !
! the RHFR program calls. We are just avoiling compute explicitly every!
! single time values like: 'work', 'iwork', 'lwork', ... in the main   !
! program.                                                             !
!***********************************************************************
    SUBROUTINE DIAG_M(n, A, B)
    ! Diagonalize the matrix A                                         !
    !               !!!The matrix A will be overwrite!!!               !
    ! n :: size of the matrix                                          !
    ! A :: input matrix, it will have the eigenvectors                 !
    ! B :: eigenvalues                                                 !
    IMPLICIT NONE
    INTEGER                                 :: n, lwork, info
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: work
    REAL(KIND=8), DIMENSION(n)              :: B
    REAL(KIND=8), DIMENSION(n,n)            :: A
    
    lwork = 3*n; ALLOCATE(work(lwork))
    
    CALL DSYEV('V', 'U', n, A, n, B, work, lwork, info)
    
    CALL WIZARD(info)

    ENDSUBROUTINE
!***********************************************************************
