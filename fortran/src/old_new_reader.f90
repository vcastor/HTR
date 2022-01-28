SUBROUTINE NEW_READER(input_name, q, n_atoms, atom, xyz, b_functions, atomic_number, n_basisf_pa, n_pri_bf, max_prim, n_max_zeta, &
                      zeta, d, info)
!-----------------------------------------------------------------------
! This SUBROUTINE read the input file for Hartree-Fock-Roothan program
!      *COMBINTATION OF TWO SUBROUTINES, CAN HAVE ISSUES/BUGS*
!                                                                       Temporary Code, I hope I hope I hope
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER ::                                        i, j, k, q, n_atoms, max_bf, max_prim, b_functions, l, n_max_zeta, info, m, n, p
INTEGER, DIMENSION(:), ALLOCATABLE ::             n_pri_bf, atomic_number, n_basisf_pa, n_prim
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE ::      xyz, zeta, d
REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE ::    zeta_atom, d_atom
CHARACTER(LEN=100) ::                             input_name
CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE ::    atom
CHARACTER(LEN=2) ::                               dummy

OPEN(11,FILE=input_name)
    READ(11,*)                                               !dummy line
    READ(11,*) n_atoms

    ALLOCATE(atom(n_atoms), atomic_number(n_atoms), n_basisf_pa(n_atoms), xyz(n_atoms,3))     !setting how many atoms

    READ(11,*)                                               !dummy line
    DO i = 1, n_atoms                      !chemical symbol, coordinates
        READ(11,*) atom(i), atomic_number(i), (xyz(i,j), j=1, 3) 
    ENDDO
    READ(11,*)                                               !dummy line
    READ(11,*) q                                            !over charge
    READ(11,*)                                               !dummy line
    READ(11,*) b_functions                     !number of basis function
    READ(11,*)                                               !dummy line
    READ(11,*) max_prim                               !max of primitives
    READ(11,FMT='(A)')                                       !dummy line

    n_max_zeta = max_prim * b_functions
    ALLOCATE(zeta(n_max_zeta,max_prim), d(n_max_zeta,max_prim), n_pri_bf(b_functions)); zeta(:,:) = 0.d0; d(:,:) = 0.d0

    m = 1; n = 0; p = 1
    DO i = 1, b_functions
        READ(11,*) k, dummy, k, l              !dummy line, (except 'l')
        READ(11,*) n_pri_bf(i)
        DO j = 1, n_pri_bf(i)
            READ(11,*) zeta(i,j), d(i,j)
        ENDDO
        IF ( l == m ) THEN
            n = n + 1
        ELSE
            n_basisf_pa(p) = n                             !SOULD BE
            n = 1; m = l; p = p + 1                        !REFACTORIZED
        ENDIF
    ENDDO
    n_basisf_pa(n_atoms) = n                               
CLOSE(11)

xyz(:,:) = 1.8897259886*xyz(:,:)                       !angstrom to bohr

info = 0                         !We are not seeing sysntaxis issues now
ENDSUBROUTINE NEW_READER
