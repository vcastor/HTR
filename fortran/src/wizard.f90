!***********************************************************************
! This SUBROUTINE will know if something isn't going well, in that case
! call the error writer and stop the program.
!***********************************************************************
SUBROUTINE WIZARD(info)
IMPLICIT NONE
INTEGER, INTENT(IN) :: info

IF (info .NE. 0) THEN
!    CALL EWRITER(info)
    STOP                                          !OFF TO SEE THE WIZARD
ENDIF
ENDSUBROUTINE
!                                        Fuímonooos, que akì espantan xC
