!***********************************************************************
! This SUBROUTINE write the result output values if everythin was
! computed ok
!***********************************************************************

SUBROUTINE WRITER(na, bf, atom, xyz, steps, threshold, Enn, E_ele)

IMPLICIT NONE
INTEGER                         :: na, bf, steps, i, j
REAL (KIND=8)                   :: threshold, Enn, E_ele
REAL (KIND=8), DIMENSION(3,na)  :: xyz
!REAL (KIND=8), DIMENSION(bf,bf) :: Hcore, F, P
CHARACTER(LEN=2), DIMENSION(na) :: atom

OPEN(616,FILE='../tmp/out.out')
    WRITE(616,FMT='(A28,I3)') "Number of atoms: ", na
    WRITE(616,FMT='(A28,I3)') "Number of basis functions: ", bf
    WRITE(616,*) "Coordinates of the system in bohr:"
    DO i = 1, na
      WRITE(616,FMT='(A4,3F12.8)') atom(i), (xyz(j,i), j=1,3)
    ENDDO
    WRITE(616,FMT='(A33,I12,A22,E9.3,A22)') "Calculation at the step: ", steps, "; with a threshold: ", threshold, &
                                            " for dencity matrix"
    WRITE(616,FMT='(A33,A20)')   '                    ', "ATOMIC UNITS"
    WRITE(616,FMT='(A33,F12.8)') 'N-N Coulomb Energy: ', Enn
    WRITE(616,FMT='(A33,F12.8)') 'Electronic Energy:  ', E_ele
    WRITE(616,FMT='(A33,F12.8)') 'TOTAL ENERGY:       ', E_ele + Enn
!    WRITE(616,*) "The system has: "
!    WRITE(616,*) 'Core Hamiltonian'
!    DO j = 1, bf
!        DO i = j, bf
!            WRITE(616,FMT='(2I5,F16.8)') i, j, Hcore(i,j)
!        ENDDO
!    ENDDO
!    WRITE(616,*) 'Fock Matrix'
!    DO j = 1, bf
!        DO i = j, bf
!            WRITE(616,FMT='(2I5,F16.8)') i, j, F(i,j)
!        ENDDO
!    ENDDO
!    WRITE(616,*) 'Density Matrix'
!    DO j = 1, bf
!        DO i = j, bf
!            WRITE(616,FMT='(2I5,F16.8)') i, j, P(i,j)
!        ENDDO
!    ENDDO
CLOSE(616)

ENDSUBROUTINE WRITER
