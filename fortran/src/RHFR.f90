                     PROGRAM HartreeFockRoothaan
!***********************************************************************
!             Restricted Hartree-Fock-Roothaan calculation             !
!                          (closed shell)                              !
!***********************************************************************
! PLEASE DO NOT ATTEMP TO SIMPLIFY THIS CODE.                          !
! KEEP THE SPACE SHUTTLE FLYING.                                       !
!                     ~Just kidding, all contributions are wellcome <3 !
!***********************************************************************
! We call this style 'space shuttle style'. Space shuttle style is     !
! meant to ensure that every branch and condition is considered and    !
! accounted for the same way code is written at NASA for apllications  !
! like the space shuttle.                                              !
!***********************************************************************
! PLEASE: FIRST, READ THE README FILE                                  !
!***********************************************************************
! Full program was developed thinking in Hartree Atomic Units:         !
! Reduced Planck constant: ℏ = 1 (\hbar),        atomic unit of action !
! Elementary charge:       e = 1,                atomic unit of charge !
! Bohr radius:           a_0 = 1,                atomic unit of length !
! Electron mass:         m_e = 1                 atomic unit of mass   !
!***********************************************************************
!        with love,                                                    !
!            Victoria Castor 2021                                      !
!***********************************************************************

!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER                                          :: i, j, k, l, m, n, step, info
INTEGER                                          :: q, na, bf, n_occ, n_ele, max_prim, max_zeta, max_bf
INTEGER, DIMENSION(:), ALLOCATABLE               :: a_num, n_pri_bf, n_pri_bf_aux,  n_bf_pa

REAL (KIND=8)                                    :: Enn, E_ele, threshold, aux1, aux2
REAL (KIND=8), DIMENSION (:), ALLOCATABLE        :: e
REAL (KIND=8), DIMENSION (:,:), ALLOCATABLE      :: S, T, V, P, F, C, X, Fp, Cp, Hcore, xyz, xyz_bf, zeta, d, Ennmat, Jm, Km, Pb
REAL (KIND=8), DIMENSION (:,:), ALLOCATABLE      :: dis_a, dis_bf, dis2_a, dis2_bf
REAL (KIND=8), DIMENSION (:,:,:,:), ALLOCATABLE  :: TwoEleInt

CHARACTER (LEN=2), DIMENSION(:), ALLOCATABLE     :: atom
CHARACTER (LEN=15)                               :: basis_set, units
CHARACTER (LEN=100)                              :: input_name
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------
!********                      INPUT DATA                       ********
!-----------------------------------------------------------------------
! We will get the basis set, molecule over charge, how many atoms, units
! of the coordinates and a flag if something is going wrong. The we will
! read the atoms and its coordiantes.
! Basis sets and atomic numbers are given in files with all program.
!-----------------------------------------------------------------------

CALL GETARG(1, input_name)

CALL OPT_READER(input_name, basis_set, q, na, units, info)

CALL WIZARD(info)                                         !Everthing ok?

ALLOCATE(atom(na), xyz(3,na))

CALL GEO_READER(input_name, na, atom, xyz, units, info)

CALL WIZARD(info)                                         !Everthing ok?

!-----------------------------------------------------------------------
!********                      BASIS SET                        ********
!-----------------------------------------------------------------------
! NOTE: the $d_{i,j}$ are already normalized, IMPORTANT to know how the
! program will compute the integrals.
!-----------------------------------------------------------------------


!               know the maximum values to allocate the memory propertly
CALL BASIS_MAX(basis_set, na, max_bf, max_prim, max_zeta, info)

ALLOCATE(a_num(na), n_bf_pa(na), &
      zeta(max_zeta,max_prim), d(max_zeta,max_prim), n_pri_bf(max_zeta))

CALL BASIS_AT(basis_set, na, atom, bf, a_num, n_bf_pa, n_pri_bf, & 
                              max_bf, max_prim, max_zeta, zeta, d, info)

CALL WIZARD(info)                                         !Everthing ok?

ALLOCATE(n_pri_bf_aux(bf))                              !Cleaning memory

DO i = 1, bf
  n_pri_bf_aux(i) = n_pri_bf(i)
ENDDO

DEALLOCATE(n_pri_bf); ALLOCATE(n_pri_bf(bf))

n_pri_bf(:) = n_pri_bf_aux(:); DEALLOCATE(n_pri_bf_aux)      !Cleaned it

!-----------------------------------------------------------------------
!********        DISTANCE BETWEEN ATOMS AND BASIS SETS          ********
!-----------------------------------------------------------------------

ALLOCATE(xyz_bf(3,bf), dis_a(na,na), dis2_a(na,na), dis_bf(bf,bf), &
                                                         dis2_bf(bf,bf))

CALL DISTANCES(na, bf, n_bf_pa, xyz, xyz_bf, dis_a, dis2_a, dis_bf, &
                                                                dis2_bf)

!-----------------------------------------------------------------------
!********            BORN-OPPENHEIMER APPROXIMATION             ********
!-----------------------------------------------------------------------

ALLOCATE(Ennmat(na,na))

CALL NNBOI(na, dis_a, a_num, Ennmat, Enn)

!-----------------------------------------------------------------------
!********                ELECTRONS AND ORBITALS                 ********
!-----------------------------------------------------------------------

n_ele = 0
DO i = 1, na                        !electrons as nuclei charge per atom
  n_ele = n_ele + a_num(i)  
ENDDO
n_ele = n_ele - q                         !do not forget the over charge

n_occ = n_ele / 2                     !how many occuped orbitals we have

!-----------------------------------------------------------------------
!********                       INTEGRALS                       ********
!-----------------------------------------------------------------------

ALLOCATE(S(bf,bf), T(bf,bf), V(bf,bf), TwoEleInt(bf,bf,bf,bf))

CALL INTEGRALS(bf, na, a_num, n_pri_bf, xyz, xyz_bf, dis2_bf, &
                  max_prim, max_zeta, zeta, d, S, T, V, TwoEleInt, info)

CALL WIZARD(info)                                         !Everthing ok?

!-----------------------------------------------------------------------------------------------------------------------------------
!********                                             HARTREE-FOCK-ROOTHAN ALRGORITHM                                       ********
!********                                                MAIN POINT OF THE PROGRAM                                          ********
!-----------------------------------------------------------------------------------------------------------------------------------

ALLOCATE(F(bf,bf), C(bf,bf), X(bf,bf), Fp(bf,bf), Cp(bf,bf), Hcore(bf,bf), P(bf,bf), Pb(bf,bf), Jm(bf,bf), Km(bf,bf), e(bf))

Hcore(:,:) = T(:,:) + V(:,:)                              !core Hamiltonian: H^{core} = KineticIntegral + Nuclear Atraction Integral
C(:,:) = 0.d0; Pb(:,:)= 1.d0                                                           !just in case the RAM was with something else

!-----------------------------------------------------------------------------------------------------------------------------------
!---- Here start the loop until the threshold smaller than $10^{-10}$

threshold = 1.d0; step = 0
DO WHILE (threshold .GT. 1d-10)

    P(:,:) = 0.d0                            !Restart the sum every loop
    DO i = 1, bf                             !Density matrix: P(i,j)
      DO j = 1, bf
        DO k = 1, n_occ
          P(j,i) = P(j,i) + C(j,k)*C(i,k)
        ENDDO
      ENDDO
    ENDDO

    Jm(:,:) = 0.d0; Km(:,:) = 0.d0                   !restart every loop
    DO n = 1, bf                         !J = P * (\mu\nu|\lambda\sigma)
      DO m = 1, bf                       !K = P * (\mu\lambda|\nu\sigma)
        DO l = 1, bf
          DO k = 1, bf                   !k :: \sigma
            Jm(m,n) = Jm(m,n) + P(k,l)*TwoEleInt(m,n,l,k)
            Km(m,n) = Km(m,n) + P(k,l)*TwoEleInt(m,l,n,k)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

!   F = T + V + P(2(\mu\nu|\lambda\sigma)-(\mu\lambda|\nu\sigma))
    F = Hcore + 2.d0*Jm - Km                                            !F = H^{core} + 2J - K 

!***********************************************************************
!********                       FC = SCE                        ********
!********                         ⬇️                             ********
!                                                                      !
    CALL ORTHOGONALIZATION_M(bf, S, X)             !     X = S^(-1/2)  !
    Fp = MATMUL(TRANSPOSE(X),MATMUL(F,X))          !    F' = X^T F X   !
    Cp = Fp; CALL DIAG_M(bf, Cp, e)                ! compute C' and ε  !
    C  = MATMUL(X,Cp)                              !      C' ➡️   C     !
!                                                                      !
!********                         ⬇️                             ********
!********                     F'C' = C'ε                        ********
!***********************************************************************

    aux1 = 0.d0; aux2 = 0.d0                           !      Converged?
    DO j = 1, bf                                 
      DO i = 1, bf
        aux1 = P(i,j) - Pb(i,j)
        aux2 = aux2 + aux1*aux1
      ENDDO
    ENDDO
    threshold = DSQRT(aux2/REAL(bf*bf))

    Pb = P                         !Save the Density Matrix of the steps

    step = step + 1                                        !Step counter
ENDDO                                                                   !Close the DO WHILE of the Hartree-Fock iterative part 

!-----------------------------------------------------------------------------------------------------------------------------------
!---- Finally. Compute the energy at the converged system

CALL EHF(bf, P, T, V, Jm, Km, E_ele)

!                                                                               ya mero cerramos el chiringuito/changarro xddxdxxdxd
!-----------------------------------------------------------------------------------------------------------------------------------
!********                                       END OF THE HARTREE-FOCK-ROOTHAN ALGORITHM                                   ********
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------
!********                 WRITE THE OUTPUT FILE                 ********
!-----------------------------------------------------------------------
! If all above was done correctly the work program was done, all about
! writing files in fancy style is on the bash/perl scripts, DANKE SCHÖN

CALL WRITER(na, bf, atom, xyz, step, threshold, Enn, E_ele, Hcore, F, P)

ENDPROGRAM
!                                                                 'ora sí ya cerramos el changarro/chiringuito bien /eichidi/ xdxddx
!                                                    ah nu ma' grax por leer el código si has llegado hasta acá mejor vete de fiesta
!                                                                    es más io te invito un Gin and Tonic y si no tomas pues chtm :b
!                                                                                                                   ah nocierto xdxd
!                                                                                               este programa no tiene Easter Eggs Ü
!-----------------------------------------------------------------------------------------------------------------------------------
