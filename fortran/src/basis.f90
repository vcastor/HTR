SUBROUTINE AT_BASIS(basis_set, n_atoms, atom, b_functions, atomic_number, n_basisf_pa, n_pri_bf, max_prim, n_max_zeta, zeta, d, &
                    info)
!-----------------------------------------------------------------------
!   This SUBROUTINE gets the data about the basis set
!-----------------------------------------------------------------------
IMPLICIT NONE
INTEGER ::                                          i, j, k, n_atoms, max_bf, max_prim, info, b_functions, l, n_max_zeta
INTEGER, DIMENSION(:), ALLOCATABLE ::               n_pri_bf
INTEGER, DIMENSION(n_atoms) ::                      atomic_number, n_basisf_pa
INTEGER, DIMENSION(:,:), ALLOCATABLE ::             n_prim
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE ::        zeta, d
REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE ::      zeta_atom, d_atom
CHARACTER(LEN=10) ::                                basis_set
CHARACTER(LEN=40) ::                                general_data
CHARACTER(LEN=2), DIMENSION(n_atoms) ::             atom
CHARACTER(LEN=30), DIMENSION(n_atoms) ::            file_basis_set

general_data="./basis/"//TRIM(basis_set)//"/data_of_basis.dat"

OPEN(11,FILE=general_data)
    READ(11,*) max_bf               !maximum of basis functions per atom
    READ(11,*) max_prim        !maximum of primitives per basis function
CLOSE(11)

!-----------------------------------------------------------------------
! Knowing the maximum number of basis functions and the atoms that the
! system has, we can know the maximum "i's" of z_{i,j} and d_{i,j}.
! Where "j's" is the maximum number of Gaussian primitives functions
! per basis function. Therefore at (i^{th}, j^{th}) value we have the
! value of i^{th} basis function at its j^{th} Gaussian Primitive
! function value. Also, we save as: at atom^{th}, at basis^{th}, at
! Gaussian primitive^{th}.
!-----------------------------------------------------------------------

n_max_zeta = n_atoms * max_bf

ALLOCATE(n_prim(n_atoms, max_bf), zeta(n_max_zeta,max_prim), d(n_max_zeta,max_prim), n_pri_bf(n_max_zeta), &
         zeta_atom(n_atoms, max_bf, max_prim), d_atom(n_atoms, max_bf, max_prim))        !(at atom, at basis, at Gaussian primitive)
DO i = 1, n_atoms
    file_basis_set(i) = "./basis/"//TRIM(basis_set)//"/"//TRIM(atom(i))//".dat"
ENDDO

zeta(:,:) = 0.d0; d(:,:) = 0.d0
b_functions = 0; l = 1
DO i = 1, n_atoms
  OPEN(12,FILE=file_basis_set(i))
    READ(12,*) atomic_number(i), n_basisf_pa(i)
!   (at the i'th atom, at the j'th basis of the i'th atom)
    DO j = 1, n_basisf_pa(i)
      READ (12,*) n_prim(i,j)
      n_pri_bf(l)      = n_prim(i,j)
!     (at the i'th atom, at the j'th basis of the i'th atom, at th k'th zeta/d
      DO k = 1, n_prim(i,j)
        READ(12,*) zeta_atom(i,j,k), d_atom(i,j,k)
        zeta(l,k)   = zeta_atom(i,j,k)
        d(l,k)      = d_atom(i,j,k)
      ENDDO
      l = l + 1
    ENDDO
  CLOSE(12)
  b_functions = b_functions + n_basisf_pa(i)
ENDDO
l = l - 1

IF (l .NE. b_functions) THEN
  info = 201
ENDIF

ENDSUBROUTINE AT_BASIS
