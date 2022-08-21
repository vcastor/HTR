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

INTEGER                                          :: i, j, k, l, lwork, liwork, steps, info
INTEGER                                          :: q, na, bf, n_occ, n_ele, max_prim, max_zeta, max_bf
INTEGER, DIMENSION(:), ALLOCATABLE               :: a_num, n_pri_bf, n_bf_pa, iwork

REAL (KIND=8)                                    :: Enn, E_ele, threshold, aux1, aux2
REAL (KIND=8), DIMENSION (:), ALLOCATABLE        :: E, work
REAL (KIND=8), DIMENSION (:,:), ALLOCATABLE      :: S, T, V, P, F, C, S_aux, Hcore, xyz, xyz_bf, zeta, d, Ennmat, Jm, Km, Pb
REAL (KIND=8), DIMENSION (:,:), ALLOCATABLE      :: dis_a, dis_bf, dis2_a, dis2_bf
REAL (KIND=8), DIMENSION (:,:,:,:), ALLOCATABLE  :: TwoEleInt

CHARACTER (LEN=2), DIMENSION(:), ALLOCATABLE     :: atom
CHARACTER (LEN=15)                               :: basis_set, units
CHARACTER (LEN=100)                              :: input_name
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------
!********                      INPUT DATA                       ********
!-----------------------------------------------------------------------
! We will get the basis set name, the molecule over charge, how many
! atoms, the chemical symbols of its atoms, the coordinates of its, and
! an information flag (to know if the input has something wrong).
! Basis sets and atomic numbers are given in files with all program.
!-----------------------------------------------------------------------

READ(UNIT=5,FMT='(100A)') input_name

CALL OPT_READER(input_name, basis_set, q, na, units, info)

CALL WIZARD(info)                                         !Everthing ok?

ALLOCATE(atom(na), xyz(na,3))

CALL GEO_READER(input_name, na, atom, xyz, units, info)

WRITE(*,*) units
DO i = 1, na
  WRITE(*,*) (xyz(i, j), j=1,3)
ENDDO
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

!ALLOCATE(n_pri_bf_aux(bf))                             !Cleaning memory
!
!DO i = 1, bf
!  n_pri_bf_aux(i) = n_pri_bf(i)
!ENDDO
!
!DEALLOCATE(n_pri_bf); ALLOCATE(n_pri_bf(bf))
!
!n_pri_bf(:) = n_pri_bf_aux(:); DEALLOCATE(n_pri_bf_aux)

!-----------------------------------------------------------------------
!********        DISTANCE BETWEEN ATOMS AND BASIS SETS          ********
!-----------------------------------------------------------------------

ALLOCATE(xyz_bf(bf,3), dis_a(na,na), dis2_a(na,na), dis_bf(bf,bf), &
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

ALLOCATE(S(bf,bf), T(bf,bf), V(bf,bf), TwoEleInt(bf,bf,bf,bf), &
                                                           S_aux(bf,bf))

CALL INTEGRALS(bf, na, a_num, n_pri_bf, xyz, xyz_bf, dis2_bf, &
                  max_prim, max_zeta, zeta, d, S, T, V, TwoEleInt, info)

CALL WIZARD(info)                                         !Everthing ok?

!-----------------------------------------------------------------------------------------------------------------------------------
!********                                             HARTREE-FOCK-ROOTHAN ALRGORITHM                                       ********
!********                                                MAIN POINT OF THE PROGRAM                                          ********
!-----------------------------------------------------------------------------------------------------------------------------------

ALLOCATE(F(bf,bf), C(bf,bf), Hcore(bf,bf), P(bf,bf), Pb(bf,bf), Jm(bf,bf), Km(bf,bf))
!                                if we rotate the matrix system, we need ALLOCATE some another matrices, descomenting the next line:
!        Fp(bf,bf), Cp(bf,bf), G(bf,bf),   X(bf,bf))

Hcore(:,:) = T(:,:) + V(:,:)                              !core Hamiltonian: H^{core} = KineticIntegral + Nuclear Atraction Integral
C(:,:) = 0.d0; Pb(:,:)= 1.d0                                                           !just in case the RAM was with something else

!-----------------------------------------------------------------------
! We will use the LIBRARY LAPACK to solve the Eigenvalue problem,
! particularlly we will CALL DSYGVD. Because of that we need compute
! before some values to ALLOCATE some arrays that DSYGVD use in the
! computing, the values does not have physical meaning.
!-----------------------------------------------------------------------

    lwork  = 1 + 6*bf + 2*bf**2                
    liwork = 3 + 5*bf
    ALLOCATE(E(bf), work(lwork), iwork(liwork))

!-----------------------------------------------------------------------------------------------------------------------------------
!---- Here start the loop until the threshold smaller than $10^-10$

threshold = 1.d0; steps = 0
DO WHILE (threshold .GT. 1d-10)

    P(:,:) = 0.d0                            !Restart the sum every loop
    DO i = 1, bf                             !Density matrix: P(i,j)
      DO j = 1, bf
        DO k = 1, n_occ
          P(i,j) = P(i,j) + C(i,k)*C(j,k)
        ENDDO
      ENDDO
    ENDDO

    Jm(:,:) = 0.d0; Km(:,:) = 0.d0       !restart every loop
    DO i = 1, bf                !J = P * (\mu\nu|\lambda\sigma)
      DO j = 1, bf              !K = P * (\mu\lambda|\nu\sigma)
        DO k = 1, bf
          DO l = 1, bf
            Jm(i,j) = Jm(i,j) + P(k,l)*TwoEleInt(i,j,l,k)
            Km(i,j) = Km(i,j) + P(k,l)*TwoEleInt(i,l,k,j)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
!   F      = T + V      + P(2(\mu\nu|\lambda\sigma)-(\mu\lambda|\nu\sigma))
    F(:,:) = Hcore(:,:) + 2.d0*Jm(:,:) - Km(:,:)                           !F = H^{core} + 2J - K 

!********                       FC = SCE                        ********
! Using LAPACK LIBRARY we will get as output the E vector, and C matrix
! in our F matrix, therefore we copy the data to C before DSYGVD calling 
! The values at S also will change, therefore we use an auxiliar S, to
! always have the correct values.

    S_aux = S; C = F
    
    CALL DSYGVD(1, 'V', 'L', bf, C, bf, S_aux, bf, E, work, lwork, iwork, liwork, info)
!                                                            If we rotate the matrix F, this is the correct LAPACK SUBROUTINE
!                                                            CALL DSYEV('V', 'L', bf, F, bf, E, work, lwork, info)
    CALL WIZARD(info)                                     !Everthing ok?

    aux1 = 0.d0; aux2 = 0.d0                                 !Converged?
    DO j = 1, bf                                 
      DO i = 1, bf
        aux1 = (P(i,j) - Pb(i,j))
        aux2 = aux2 + aux1*aux1
      ENDDO
    ENDDO

    threshold = DSQRT(aux2*1.d0/REAL(bf*bf))
    Pb = P                         !Save the Density Matrix of the steps

    steps = steps + 1                                      !Step counter
ENDDO                                                                   !Close the DO WHILE of the Hartree-Fock iterative part 

!-----------------------------------------------------------------------------------------------------------------------------------
!---- Finally. Compute the energy at the converged system

E_ele = 0.d0              
DO i = 1, bf
  DO j = 1, bf
    E_ele = E_ele + P(j,i)*(Hcore(i,j) + F(i,j))
  ENDDO
ENDDO

!                                                                               ya mero cerramos el chiringuito/changarro xddxdxxdxd
!-----------------------------------------------------------------------------------------------------------------------------------
!********                                       END OF THE HARTREE-FOCK-ROOTHAN ALGORITHM                                   ********
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------
!********                 WRITE THE OUTPUT FILE                 ********
!-----------------------------------------------------------------------
! If all above was done correctly the work program was done, all about
! writing files in fancy style is on the bash/perl scripts, DANKE SCHÖN

CALL WRITER(na, bf, atom, xyz, steps, threshold, Enn, E_ele, Hcore, F, P)

ENDPROGRAM
!                                                                 'ora sí ya cerramos el changarro/chiringuito bien /eichidi/ xdxddx
!                                                   ah nu ma' grax por leer el código si has llegado hasta acá, mejor vete de fiesta
!                                                                                               este programa no tiene Easter Eggs Ü
!-----------------------------------------------------------------------------------------------------------------------------------
