PROGRAM HartreeFockRoothaan

!***********************************************************************
!       This program do Hartree-Fock-Roothaan calculation method 
! Restricted Hartree-Fock-Roothaan calculation
! (closed shell)
!***********************************************************************
! PLEASE DO NOT ATTEMP TO SIMPLIFY THIS CODE.    
! KEEP THE SPACE SHUTTLE FLYING.     
!                      ~Just kidding, all contributions are wellcome <3
!***********************************************************************
! We call this style 'space shuttle style'. Space shuttle style is meant
! to ensure that every branch and condition is considered and accounted
! for the same way code is written at NASA for apllications like the --
! space shuttle.
!***********************************************************************
! PLEASE: FIRST, READ THE README FILE
!***********************************************************************
! Full program was developed thinking for Hartree Atomic Units:
! Reduced Planck constant: ℏ = 1 (\hbar),          atomic unit of action
! Elementary charge:       e = 1,                  atomic unit of charge
! Bohr radius:           a_0 = 1,                  atomic unit of length
! Electron mass:         m_e = 1                   atomic unit of mass
!***********************************************************************
!    with love,                                              
!        Victoria Castor 2021                                           
!***********************************************************************

!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER                                          :: i, j, k, l, q, n_atoms, b_functions, n_occ, n_ele, lwork, liwork, steps
INTEGER                                          :: max_prim, n_max_zeta, info
INTEGER, DIMENSION(:), ALLOCATABLE               :: atomic_number, iwork, n_pri_bf, n_basisf_pa

REAL (KIND=8)                                    :: E_cou, E_ele, threshold, aux1, aux2
REAL (KIND=8), DIMENSION (:), ALLOCATABLE        :: E, work
REAL (KIND=8), DIMENSION (:,:), ALLOCATABLE      :: S, T, V, P, F, C, S_aux, Hcore, xyz, xyz_bf, zeta, d, Ecoulomb, Jm, Km, Pb
REAL (KIND=8), DIMENSION (:,:), ALLOCATABLE      :: distance_atoms, distance_bf
REAL (KIND=8), DIMENSION (:,:,:,:), ALLOCATABLE  :: TwoEleInt

CHARACTER (LEN=2), DIMENSION (:), ALLOCATABLE    :: atom
CHARACTER (LEN=15)                               :: basis_set
CHARACTER (LEN=100)                              :: input_name
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------
!********                      INPUT DATA                       ********
!-----------------------------------------------------------------------
! The program will CALL a FORTRAN SUBROUTINE to read the input file.
! We will get the basis set name, the molecule over charge, how many
! atoms, the chemical symbols of its atoms, the coordinates of its, and
! an information flag (to know if the input has something wrong). The
! really first argument is the input file name.
!-----------------------------------------------------------------------

READ(UNIT=5,FMT='(100A)') input_name

CALL READER(input_name, basis_set, q, n_atoms, atom, xyz, info)

IF (info .NE. 0) THEN
    WRITE(*,*) info
    STOP                                          !OFF TO SEE THE WIZARD
ENDIF

!-----------------------------------------------------------------------
!********                      BASIS SET                        ********
!-----------------------------------------------------------------------
! The program will CALL a FORTRAN SUBROUTINE to get the necessary
! information about the basis set.
! We will give the basis set name, how many atoms and which atoms we
! have, getting the number of basis functions, how many basis functions
! per atom, how many primitive Gaussians functions per basis function,
! the $\zeta$ and $d_{i,j}$ values for each basis function and a info
! flag if we have any error.
! NOTE: the $d_{i,j}$ are already normalized, IMPORTANT to know how the
! program will compute the integrals.
!-----------------------------------------------------------------------

ALLOCATE(atomic_number(n_atoms), n_basisf_pa(n_atoms))

CALL AT_BASIS(basis_set, n_atoms, atom, b_functions, atomic_number, n_basisf_pa, n_pri_bf, max_prim, n_max_zeta, zeta, d, info)

IF (info .NE. 0) THEN
    WRITE(*,*) info
    STOP                                          !OFF TO SEE THE WIZARD
ENDIF

!-----------------------------------------------------------------------
!********        DISTANCE BETWEEN ATOMS AND BASIS SETS          ********
!-----------------------------------------------------------------------
! The protram will CALL a FORTRAN SUBROUTINE to fill a matrix where will
! be the distances between every atom with each other in the system.
!-----------------------------------------------------------------------

ALLOCATE(xyz_bf(b_functions,3), distance_atoms(n_atoms,n_atoms), distance_bf(b_functions,b_functions))

CALL DISTANCES(n_atoms, b_functions, n_basisf_pa, xyz, xyz_bf, distance_atoms, distance_bf)

!-----------------------------------------------------------------------
!********            BORN-OPPENHEIMER APPROXIMATION             ********
!-----------------------------------------------------------------------
! The program will CALL a FORTRAN SUBROUTINE to fill a matrix where will
! be the coulombian interactions for nuclei: $V_{n,n}$
!-----------------------------------------------------------------------

ALLOCATE(Ecoulomb(n_atoms,n_atoms))

CALL NNBOI(n_atoms, distance_atoms, atomic_number, Ecoulomb, E_cou)

!-----------------------------------------------------------------------
!********                ELECTRONS AND ORBITALS                 ********
!-----------------------------------------------------------------------

n_ele = 0
DO i = 1, n_atoms                   !electrons as nuclei charge per atom
    n_ele = n_ele + atomic_number(i)  
ENDDO
n_ele = n_ele - q                         !do not forget the over charge

n_occ = n_ele / 2                     !how many occuped orbitals we have

!-----------------------------------------------------------------------
!********                       INTEGRALS                       ********
!-----------------------------------------------------------------------
! The program will CALL a FORTRAN SUBROUTINE to compute all integrals
! requiered to the Hartree-Fock-Roothan calculation.
! Giving the basis functions, the atomic number, how many basis
! functions; taking back the integrals on the respective matrices.
!-----------------------------------------------------------------------

ALLOCATE(S(b_functions,b_functions), T(b_functions,b_functions), V(b_functions,b_functions), &
         TwoEleInt(b_functions,b_functions,b_functions,b_functions), S_aux(b_functions,b_functions))

CALL INTEGRALS_HF(b_functions, n_atoms, atomic_number, n_pri_bf, xyz, xyz_bf, distance_bf, max_prim, n_max_zeta, zeta, &
                  d, S, T, V, TwoEleInt, info)

IF (info .NE. 0) THEN
    WRITE(*,*) info
    STOP                                          !OFF TO SEE THE WIZARD
ENDIF

!-----------------------------------------------------------------------------------------------------------------------------------
!********                                             HARTREE-FOCK-ROOTHAN ALRGORITHM                                       ********
!********                                                MAIN POINT OF THE PROGRAM                                          ********
!-----------------------------------------------------------------------------------------------------------------------------------

ALLOCATE(F(b_functions,b_functions), C(b_functions,b_functions), Hcore(b_functions,b_functions), P(b_functions,b_functions), &
!                                if we rotate the matrix system, we need ALLOCATE some another matrices, descomenting the next line:
!        Fp(b_functions,b_functions), Cp(b_functions,b_functions), G(b_functions,b_functions),   X(b_functions,b_functions), &
         Pb(b_functions,b_functions), Jm(b_functions,b_functions), Km(b_functions,b_functions))

Hcore(:,:) = T(:,:) + V(:,:)                              !core Hamiltonian: H^{core} = KineticIntegral + Nuclear Atraction Integral
C(:,:) = 0.d0; Pb(:,:)= 1.d0                                                           !just in case the RAM was with something else

!-----------------------------------------------------------------------
! We will use the LIBRARY LAPACK to solve the Eigenvalue problem,
! particularlly we will CALL DSYGVD. Because of that we need compute
! before some values to ALLOCATE some arrays that DSYGVD use in the
! computing, the values does not have physical meaning.
! Code out of the loop to not recompute the value.
!-----------------------------------------------------------------------

    lwork  = 1 + 6*b_functions + 2*b_functions**2                
    liwork = 3 + 5*b_functions
    ALLOCATE(E(b_functions), work(lwork), iwork(liwork))

!-----------------------------------------------------------------------------------------------------------------------------------
!---- Here start the loop until the threshold smaller than $10^-10$

threshold = 1.d0; steps = 0
DO WHILE (threshold .GT. 1d-10)

    P(:,:) = 0.d0                            !Restart the sum every loop
    DO i = 1, b_functions                        !Density matrix: P(i,j)
        DO j = 1, b_functions
            DO k = 1, n_occ
                P(i,j) = P(i,j) + C(i,k)*C(j,k)
            ENDDO
        ENDDO
    ENDDO

    Jm(:,:) = 0.d0; Km(:,:) = 0.d0       !restart every loop
    DO i = 1, b_functions                !J = P * (\mu\nu|\lambda\sigma)
        DO j = 1, b_functions            !K = P * (\mu\lambda|\nu\sigma)
            DO k = 1, b_functions
                DO l = 1, b_functions
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
    
    CALL DSYGVD(1, 'V', 'L', b_functions, C, b_functions, S_aux, b_functions, E, work, lwork, iwork, liwork, info)
!                                                            If we rotate the matrix F, this is the correct LAPACK SUBROUTINE
!                                                            CALL DSYEV('V', 'L', b_functions, F, b_functions, E, work, lwork, info)
    IF (info .NE. 0) THEN
        WRITE(*,*) info
        STOP                                      !OFF TO SEE THE WIZARD
    ENDIF

    aux1 = 0.d0; aux2 = 0.d0                              !Is converged?
    DO i = 1, b_functions                                 
        DO j = 1, b_functions
            aux1 = (P(i,j) - Pb(i,j))
            aux2 = aux2 + aux1*aux1
        ENDDO
    ENDDO

    threshold = DSQRT(aux2*1.d0/REAL(b_functions*b_functions))
    Pb = P                         !Save the Density Matrix of the steps

    steps = steps + 1                                      !Step counter
ENDDO                                                                   !Close the DO WHILE of the Hartree-Fock iterative part 

!-----------------------------------------------------------------------------------------------------------------------------------
!---- Finally. Compute the energy at the converged system

E_ele = 0.d0              
DO i = 1, b_functions
    DO j = 1, b_functions
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

OPEN(616,FILE='./tmp/out.out')
    WRITE(616,FMT='(A28,I3)') "Number of atoms: ", n_atoms
    WRITE(616,FMT='(A28,I3)') "Number of basis functions: ", b_functions
    WRITE(616,*) "Coordinates of the system in bohr:"
    DO i = 1, n_atoms
        WRITE(616,FMT='(A4,3F12.8)') atom(i), (xyz(i,j), j=1,3)
    ENDDO
    WRITE(616,FMT='(A33,I12,A22,E9.3,A22)') "Calculation at the step: ", steps, "; with a threshold: ", threshold, &
                                            " for dencity matrix"
    WRITE(616,FMT='(A33,A20)')   '                    ', "ATOMIC UNITS"
    WRITE(616,FMT='(A33,F12.8)') 'N-N Coulomb Energy: ', E_cou
    WRITE(616,FMT='(A33,F12.8)') 'Electronic Energy:  ', E_ele
    WRITE(616,FMT='(A33,F12.8)') 'TOTAL ENERGY:       ', E_ele + E_cou
    WRITE(616,*) "The system has: "
    WRITE(616,*) 'Core Hamiltonian'
    DO i = 1, b_functions
        DO j = i, b_functions
            WRITE(616,FMT='(2I5,F16.8)') i, j, Hcore(i,j)
        ENDDO
    ENDDO
    WRITE(616,*) 'Fock Matrix'
    DO i = 1, b_functions
        DO j = i, b_functions
            WRITE(616,FMT='(2I5,F16.8)') i, j, F(i,j)
        ENDDO
    ENDDO
    WRITE(616,*) 'Density Matrix'
    DO i = 1, b_functions
        DO j = i, b_functions
            WRITE(616,FMT='(2I5,F16.8)') i, j, P(i,j)
        ENDDO
    ENDDO
CLOSE(616)

!                                                                 'ora sí ya cerramos el changarro/chiringuito bien /eichidi/ xdxddx
!                                                   ah nu ma' grax por leer el código si has llegado hasta acá, vete de fiesta mejor
!                                                                                               este programa no tiene Easter Eggs Ü
!-----------------------------------------------------------------------------------------------------------------------------------
ENDPROGRAM HartreeFockRoothaan
