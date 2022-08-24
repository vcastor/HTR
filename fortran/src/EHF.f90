!-----------------------------------------------------------------------
! This SUBROUTINE computes the electronic energy at converget system   !
! for the program RHFR                                                 !
! bf    :: basis function, size of the matrices                        !
! P     :: density matrix                                              !
! T     :: kinetic matrix                                              !
! V     :: potential matrix                                            !
! Jm    :: coulomb integral matrix                                     !
! Km    :: exchange integral matrix                                    !
! E_ele :: electronic energy at RHF                                    !
!-----------------------------------------------------------------------
SUBROUTINE EHF(bf, P, T, V, Jm, Km, E_ele)
IMPLICIT NONE
INTEGER                        :: bf
REAL(KIND=8)                   :: E_ele, ET, EV, EJ, EK
REAL(KIND=8), EXTERNAL         :: f_trace_mat
REAL(KIND=8), DIMENSION(bf,bf) :: P, T, V, Jm, Km

ET = f_trace_mat(bf, MATMUL(P,T))
EV = f_trace_mat(bf, MATMUL(P,V))
EJ = f_trace_mat(bf, MATMUL(P,Jm))
EK = f_trace_mat(bf, MATMUL(P,Km))

E_ele = ET + EV + EJ - 0.5d0*EK
E_ele = 2.d0*E_ele

ENDSUBROUTINE EHF
