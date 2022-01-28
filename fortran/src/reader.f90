SUBROUTINE READER(input_name, basis_set, q, n_atoms, atom, xyz, info)
!-----------------------------------------------------------------------
!   This SUBROUTINE read the input file for Hartree-Fock-Roothan.f90
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER ::                                     i, j, q, n_atoms, info
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE ::   xyz
CHARACTER(LEN=15) ::                           basis_set, units
CHARACTER(LEN=100) ::                          input_name
CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE :: atom

OPEN(11,FILE=input_name)
    READ(11,*) basis_set
    READ(11,*) q                                            !over charge
    READ(11,*) units
    READ(11,*) n_atoms
    READ(11,*)                       !the comment as xyz standar/jobname

    ALLOCATE(atom(n_atoms), xyz(n_atoms,3))      !setting how many atoms

    DO i = 1, n_atoms                      !chemical symbol, coordinates
        READ(11,*) atom(i), (xyz(i,j), j=1, 3)
    ENDDO
CLOSE(11)

SELECTCASE (units)
    CASE ('angstrom')
        xyz(:,:) = 1.8897259886*xyz(:,:)
        info = 0
    CASE ('bohr')
        info = 0
    CASE DEFAULT
        info = 101
ENDSELECT

ENDSUBROUTINE READER
