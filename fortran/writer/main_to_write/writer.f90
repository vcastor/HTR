PROGRAM writer_error

!-----------------------------------------------------------------------
!   This SUBROUTINE write a file if we have troubles 
!***********************************************************************
!       with hate,
!               Victoria Castor 2021
!-----------------------------------------------------------------------
INTEGER :: info

OPEN(666,FILE="error_HF.out")

SELECTCASE (info)

!-----------------------------------------------------------------------------------------------------------------------------------
    CASE(101)
        WRITE(666,*) 'The only 2 units that the program can read are:'
        WRITE(666,*) 'i)  angstorom'
        WRITE(666,*) 'ii) bohr'
        WRITE(666,*) 'working in more, like Imperial Units... of course not'

    CASE(201)
        WRITE(666,*) 'There are someting wrong in the basis set dat, PLEASE:'
        WRITE(666,*) 'contact to Victoria Castor [at] the next:'
        WRITE(666,*) 'e-mail:      vcastorv@gmail.com'
        WRITE(666,*) 'Instagram:   @vcastorv'
        WRITE(666,*) 'Twitter:     @vcastor'

ENDSELECT

ENDPROGRAM writer_error
