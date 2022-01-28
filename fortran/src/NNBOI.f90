SUBROUTINE NNBOI(n_atoms, distance_mat, atomic_number, Ecoulomb, E_cou)
!-----------------------------------------------------------------------
!   This SUBROUTINE fill a matrix with $V_{n,n}$ at BO approximation
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER ::                                   i, j, n_atoms
INTEGER, DIMENSION(n_atoms) ::               atomic_number
REAL(KIND=8) ::                              E_cou
REAL(KIND=8), DIMENSION(n_atoms,n_atoms) ::  distance_mat, Ecoulomb

E_cou = 0.d0
DO i = 1, n_atoms-1
    DO j = i+1, n_atoms
        Ecoulomb(i,j) = f_coulomb(atomic_number(i), atomic_number(j), distance_mat(i,j))
        Ecoulomb(j,i) = Ecoulomb(i,j)
        E_cou = E_cou + Ecoulomb(i,j)
    ENDDO
ENDDO

ENDSUBROUTINE NNBOI
