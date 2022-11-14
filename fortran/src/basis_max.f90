!-----------------------------------------------------------------------
!   This SUBROUTINE gets the maxima size values for the memory         !
!   about the basis set that we'll use and if we have that basis set   !
! basis_set :: basis set name                                          !
! na        :: number of atoms                                         !
! max_bf    :: maximum number of basis functions                       !
! max_prim  :: maximum number of primitives                            !
! max_zeta  :: maximum number of zeta                                  !
! info      :: info flag                                               !
!-----------------------------------------------------------------------
SUBROUTINE BASIS_MAX(basis_set, na, max_bf, max_prim, max_zeta)
IMPLICIT NONE
INTEGER           :: na, max_bf, max_prim, info, max_zeta
CHARACTER(LEN=15) :: basis_set
CHARACTER(LEN=50) :: general_data

SELECTCASE (basis_set)
  CASE('STO-3G')
    info = 0
  CASE('sto-3G')
    basis_set = 'STO-3G'
    info      = 0
  CASE('6-31G')
    info = 0
  CASE('6-31g')
    basis_set = '6-31G'
    info      = 0
  CASE DEFAULT
    info = 201
ENDSELECT

CALL WIZARD(info)                                         !Everthing ok?

general_data="../basis/"//TRIM(basis_set)//"/data_of_basis.dat"

OPEN(11,FILE=general_data)
  READ(11,*) max_bf                 !maximum of basis functions per atom
  READ(11,*) max_prim          !maximum of primitives per basis function
CLOSE(11)

!-----------------------------------------------------------------------
! Knowing the maximum number of basis functions and the atoms that the !
! system has, we can know the maximum "i's" of z_{i,j} and d_{i,j}.    !
! Where "j's" is the maximum number of Gaussian primitives functions   !
! per basis function. Therefore at (i^{th}, j^{th}) value we have the  !
! value of i^{th} basis function at its j^{th} Gaussian Primitive      !
! function value. Also, we save as: at atom^{th}, at basis^{th}, at    !
! Gaussian primitive^{th}.                                             !
!-----------------------------------------------------------------------

max_zeta = na * max_bf

ENDSUBROUTINE BASIS_MAX
