!***********************************************************************
! This SUBROUTINE will be called if the info flag in the program has   !
! information about that something is not working correctly, for no-   !
! reported issues please contact to:                                   !
!       Victoria Castor                                                !
!                vcastorv@gmail.com                                    !
!                @vcastorv    Instagram and Twitter                    !
!***********************************************************************
SUBROUTINE WIZARD(info)
IMPLICIT NONE
INTEGER, INTENT(IN) :: info

IF (info .NE. 0) THEN
    CALL ERROR_WRITER(info)
    STOP                                          !OFF TO SEE THE WIZARD
ENDIF
!                                        Fuímonooos, que akì espantan xC
ENDSUBROUTINE
