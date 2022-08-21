!-----------------------------------------------------------------------
! This SUBROUTINE reads the coordiantes in:
!             input file for Hartree-Fock-Roothan program
! input_name :: file input name
! basis_set  :: basis set name
! na         :: number of atoms
! atom       :: atom chemical characer
! xyz        :: coordinates 
! info       :: info flag
!-----------------------------------------------------------------------
SUBROUTINE GEO_READER(input_name, na, atom, xyz, units, info)

IMPLICIT NONE
INCLUDE 'parameters.h'
INTEGER                         :: i, j, na, info
REAL(KIND=8), DIMENSION(na,3)   :: xyz
CHARACTER(LEN=15)               :: units
CHARACTER(LEN=100)              :: input_name
CHARACTER(LEN=2), DIMENSION(na) :: atom

OPEN(11,FILE=input_name)

  DO i = 1, 5; READ(11,*); ENDDO         !dummy lines

  DO i = 1, na                           !chemical symbol, coordinates
    READ(11,*) atom(i), (xyz(i,j), j=1, 3)
  ENDDO

CLOSE(11)

SELECTCASE (units)
  CASE ('angstrom')
    xyz(:,:) = atb*xyz(:,:)
    info = 0
  CASE ('bohr')
    info = 0
  CASE DEFAULT
    info = 101
ENDSELECT

ENDSUBROUTINE GEO_READER
