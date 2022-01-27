PROGRAM try2

IMPLICIT NONE
INTEGER :: i, j, info, lwork, liwork
REAL (KIND=8), DIMENSION(2,2) :: F, S, C
!REAL (KIND=8), DIMENSION(2) :: E, w
INTEGER, DIMENSION(:), ALLOCATABLE :: iwork
REAL (KIND=8), DIMENSION(:), ALLOCATABLE :: work, w

F(1,1) = -2.61253; F(1,2) = -1.43358
F(2,1) = -1.43358; F(2,2) = -1.73439

S(1,1) = 1.000000; S(1,2) = 0.53707
S(2,1) = 0.537070; S(2,2) = 1.00000

WRITE(*,*) '****'
WRITE(*,*) F(1,1),F(1,2)
WRITE(*,*) F(1,1),F(2,2)
WRITE(*,*) '****'
WRITE(*,*) S(1,1),S(1,2)
WRITE(*,*) S(1,1),S(2,2)
WRITE(*,*) '****'

lwork  = 1 + 6*2 + 2*2**2
liwork = 3 + 5*2
ALLOCATE(w(lwork),work(lwork),iwork(liwork))

!E(:) = 0.
!E,C = eigh(F, S, eigvals_only=False)
!call sygvd(a, b, w [,itype] [,jobz] [,uplo] [,info])
!CALL SYGVD(F, S, E, 1, 'V', 'L')!, info)! [,itype] [,jobz] [,uplo] [,info])
!call ssygvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, iwork, liwork, info)
!CALL ssygvd(1, 'V', 'L', 2, F, 2, S, 2, w, work, lwork, iwork, liwork, info)
CALL dsygvd(1, 'V', 'L', 2, F, 2, S, 2, w, work, lwork, iwork, liwork, info)

WRITE(*,*) info

WRITE(*,*) w(1), w(2)
WRITE(*,*) F(1,1), F(1,2)
WRITE(*,*) F(1,1), F(2,2)
WRITE(*,*) S(1,1), S(1,2)
WRITE(*,*) S(1,1), S(2,2)

ENDPROGRAM try2
