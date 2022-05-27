!-----------------------------------------------------------------------
!   This SUBROUTINE fill a matrix with the atom and basis f distances
!-----------------------------------------------------------------------
SUBROUTINE DISTANCES(n_atoms, b_functions, n_basisf_pa, xyz, xyz_bf, distance_atoms, distance_bf)

IMPLICIT NONE
INTEGER ::                                           i, j, k, l, m, n, n_atoms, b_functions
INTEGER, DIMENSION(n_atoms) ::                       n_basisf_pa
REAL(KIND=8), DIMENSION(n_atoms, 3) ::               xyz
REAL(KIND=8), DIMENSION(b_functions, 3) ::           xyz_bf
REAL(KIND=8), DIMENSION(n_atoms, n_atoms) ::         distance_atoms
REAL(KIND=8), DIMENSION(b_functions, b_functions) :: distance_bf

DO i = 1, n_atoms                                            !distance between atoms
    DO j = i, n_atoms
        distance_atoms(i,j) = f_distance(xyz(i,:), xyz(j,:))
        distance_atoms(j,i) = distance_atoms(i,j)            !same distance between A & B than B & A
    ENDDO
ENDDO

k = 1                                                        !counter
DO i = 1, n_atoms                                            !distance between basis functions
  DO j = 1, n_basisf_pa(i)
    xyz_bf(k,:) = xyz(i,:)
    k = k + 1
  ENDDO
ENDDO

n = 1; m = 1                                                 !counters
DO i = 1, n_atoms
  DO j = 1, n_atoms
    DO k = 1, n_basisf_pa(i)
      DO l = 1, n_basisf_pa(j)
        distance_bf(m,n) = distance_atoms(i,j)
        n = n + 1                                            !only n_basisf_pa times
      ENDDO
      n = n - n_basisf_pa(j)                                 !restart every time that we jump to the next row
      m = m + 1
    ENDDO
    n = n + n_basisf_pa(j)                                   !start where we stop the last time with the new vale
    m = m - n_basisf_pa(i)                                   !restart every time that we jump to the next column
  ENDDO
  n = 1                                                      !every row
  m = m + n_basisf_pa(i)
ENDDO

ENDSUBROUTINE DISTANCES
