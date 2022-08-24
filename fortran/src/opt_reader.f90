!-----------------------------------------------------------------------
! This SUBROUTINE reads the options in the:
!             input file for Hartree-Fock-Roothan program
! input_name :: file input name
! basis_set  :: basis set name
! q          :: over charge
! na         :: number of atoms
! info       :: info flag
!-----------------------------------------------------------------------
SUBROUTINE OPT_READER(input_name, basis_set, q, na, units, info)

IMPLICIT NONE
INTEGER            :: q, na, info
CHARACTER(LEN=15)  :: basis_set, units
CHARACTER(LEN=100) :: input_name

OPEN(11,FILE=input_name)
    READ(11,*) basis_set
    READ(11,*) q 
    READ(11,*) units
    READ(11,*) na
CLOSE(11)

SELECTCASE (units)
    CASE ('angstrom')
        info = 0
    CASE ('bohr')
        info = 0
    CASE('Angstrom')
        units = 'angstrom'
        info  = 0
    CASE('Bohr')
        units = 'bohr'
        info  = 0
    CASE('ANGSTROM')
        units = 'angstrom'
        info  = 0
    CASE('BOHR')
        units = 'bohr'
        info  = 0
    CASE DEFAULT
        info = 101
ENDSELECT

IF (na .LE. 0) info = 102

ENDSUBROUTINE OPT_READER
