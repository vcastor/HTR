!***********************************************************************
! This SUBROUTINE will write an output with the information of the --- !
! error if it is known, if not, it will just write the number of the - !
! info flag, please report any trouble that it is not known yet to:    !
!           Victoria Castor                                            !
!               vcastorv@gmail.com                                     !
!               @vcastorv           Instagram and Twitter              !
!***********************************************************************
SUBROUTINE ERROR_WRITER(info)
IMPLICIT NONE
INTEGER            :: info
CHARACTER (LEN=10) :: fin

fin = "../tmp/"//"001"

OPEN(666,FILE=fin)
  WRITE(666,*) info
CLOSE(666)

ENDSUBROUTINE 
