PROGRAM HartreeFockRoothaan

!***********************************************************************
!       This program do Hartree-Fock-Roothaan calculation method 
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
INTEGER ::                                          i, j, k, l, q, n_atoms, b_functions, n_occ, n_ele, lwork, liwork, steps, &
                                                    max_prim, n_max_zeta, info
INTEGER, DIMENSION(:), ALLOCATABLE ::               atomic_number, iwork, n_pri_bf, n_basisf_pa
REAL (KIND=8) ::                                    E_cou, E_ele, threshold, aux1, aux2
REAL (KIND=8), DIMENSION (:), ALLOCATABLE ::        E, work
REAL (KIND=8), DIMENSION (:,:), ALLOCATABLE ::      S, T, V, P, F, C, S_aux, Hcore, xyz, xyz_bf, zeta, d, Ecoulomb, Jm, Km, Pb, &
                                                    distance_atoms, distance_bf
REAL (KIND=8), DIMENSION (:,:,:,:), ALLOCATABLE ::  TwoEleInt
CHARACTER (LEN=2), DIMENSION (:), ALLOCATABLE ::    atom
CHARACTER (LEN=15) ::                               basis_set
CHARACTER (LEN=100) ::                              input_name
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------
!********                INPUT DATA & BASIS SET                 ********
!-----------------------------------------------------------------------
! The program will CALL a FORTRAN SUBROUTINE to read the input file.  --
! Examples are given in a directory 'expamples_inputs'. The program,  --
! 'ask' about the input file where the system is. Actually the laucher
! stript does that.
! Giving the input name the program will get: over charge, number of  --
! atoms and which atoms are in the system, its coordinates, and also, --
! about the basis function, how many primitive gaussian are at every  --
! basis function, and $\zeta$ and $d_{i,j}$ values for each basis     --
! function and an info flag to know if the input file has any syntaxis
! error.
! NOTE: the $d_{i,j}$ are already normalized, IMPORTANT to know how the
! program will compute the integrals.
!-----------------------------------------------------------------------

READ(UNIT=5,FMT='(100A)') input_name

CALL NEW_READER(input_name, q, n_atoms, atom, xyz, b_functions, atomic_number, n_basisf_pa, n_pri_bf, max_prim, n_max_zeta, zeta, &
                d, info)

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

CALL BOCInteraction(n_atoms, distance_atoms, atomic_number, Ecoulomb, E_cou)

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
! The program was developed using all SUBROUTINES and FUNCTIONS in
! separated files. However, trying to avoild any issue that could be
! about that, all is on the same file to compile.

CONTAINS

!***********************************************************************
!********************   HONI SOIT QUI MAL Y PENSE   ********************
!********************        NOLI ME TANGERE        ********************
!***********************************************************************
SUBROUTINE READER(input_name, basis_set, q, n_atoms, atom, xyz, info)
!-----------------------------------------------------------------------
!   This SUBROUTINE read the input file for Hartree-Fock-Roothan.f90
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER ::                                     i, j, q, n_atoms, info
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE ::   xyz
CHARACTER(LEN=15) ::                           basis_set, units
CHARACTER(LEN=100) ::                          input_name
CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE :: atom

OPEN(11,FILE=input_name)
    READ(11,*) basis_set
    READ(11,*) q                                            !over charge
    READ(11,*) units
    READ(11,*) n_atoms
    READ(11,*)                       !the comment as xyz standar/jobname

    ALLOCATE(atom(n_atoms), xyz(n_atoms,3))      !setting how many atoms

    DO i = 1, n_atoms                      !chemical symbol, coordinates
        READ(11,*) atom(i), (xyz(i,j), j=1, 3)
    ENDDO
CLOSE(11)

SELECTCASE (units)
    CASE ('angstrom')
        xyz(:,:) = 1.8897259886*xyz(:,:)
        info = 0
    CASE ('bohr')
        info = 0
    CASE DEFAULT
        info = 101
ENDSELECT

ENDSUBROUTINE READER
!***********************************************************************************************************************************
SUBROUTINE AT_BASIS(basis_set, n_atoms, atom, b_functions, atomic_number, n_basisf_pa, n_pri_bf, max_prim, n_max_zeta, zeta, d, &
                    info)
!-----------------------------------------------------------------------
!   This SUBROUTINE gets the data about the basis set
!-----------------------------------------------------------------------
IMPLICIT NONE
INTEGER ::                                          i, j, k, n_atoms, max_bf, max_prim, info, b_functions, l, n_max_zeta
INTEGER, DIMENSION(:), ALLOCATABLE ::               n_pri_bf
INTEGER, DIMENSION(n_atoms) ::                      atomic_number, n_basisf_pa
INTEGER, DIMENSION(:,:), ALLOCATABLE ::             n_prim
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE ::        zeta, d
REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE ::      zeta_atom, d_atom
CHARACTER(LEN=10) ::                                basis_set
CHARACTER(LEN=40) ::                                general_data
CHARACTER(LEN=2), DIMENSION(n_atoms) ::             atom
CHARACTER(LEN=30), DIMENSION(n_atoms) ::            file_basis_set

general_data="./basis/"//TRIM(basis_set)//"/data_of_basis.dat"

OPEN(11,FILE=general_data)
    READ(11,*) max_bf               !maximum of basis functions per atom
    READ(11,*) max_prim        !maximum of primitives per basis function
CLOSE(11)

!-----------------------------------------------------------------------
! Knowing the maximum number of basis functions and the atoms that the
! system has, we can know the maximum "i's" of z_{i,j} and d_{i,j}.
! Where "j's" is the maximum number of Gaussian primitives functions
! per basis function. Therefore at (i^{th}, j^{th}) value we have the
! value of i^{th} basis function at its j^{th} Gaussian Primitive
! function value. Also, we save as: at atom^{th}, at basis^{th}, at
! Gaussian primitive^{th}.
!-----------------------------------------------------------------------

n_max_zeta = n_atoms * max_bf

ALLOCATE(n_prim(n_atoms, max_bf), zeta(n_max_zeta,max_prim), d(n_max_zeta,max_prim), n_pri_bf(n_max_zeta), &
         zeta_atom(n_atoms, max_bf, max_prim), d_atom(n_atoms, max_bf, max_prim))        !(at atom, at basis, at Gaussian primitive)
DO i = 1, n_atoms
    file_basis_set(i) = "./basis/"//TRIM(basis_set)//"/"//TRIM(atom(i))//".dat"
ENDDO

zeta(:,:) = 0.d0; d(:,:) = 0.d0
b_functions = 0; l = 1
DO i = 1, n_atoms
    OPEN(12,FILE=file_basis_set(i))
        READ(12,*) atomic_number(i), n_basisf_pa(i)
!      (at the i'th atom, at the j'th basis of the i'th atom)
        DO j = 1, n_basisf_pa(i)
            READ (12,*) n_prim(i,j)
            n_pri_bf(l)      = n_prim(i,j)
!      (at the i'th atom, at the j'th basis of the i'th atom, at th k'th zeta/d
            DO k = 1, n_prim(i,j)
                READ(12,*) zeta_atom(i,j,k), d_atom(i,j,k)
                zeta(l,k)   = zeta_atom(i,j,k)
                d(l,k)      = d_atom(i,j,k)
            ENDDO
        l = l + 1
        ENDDO
    CLOSE(12)
    b_functions = b_functions + n_basisf_pa(i)
ENDDO
l = l - 1

IF (l .NE. b_functions) THEN
    info = 201
ENDIF

ENDSUBROUTINE AT_BASIS
!***********************************************************************************************************************************
SUBROUTINE NEW_READER(input_name, q, n_atoms, atom, xyz, b_functions, atomic_number, n_basisf_pa, n_pri_bf, max_prim, n_max_zeta, &
                      zeta, d, info)
!-----------------------------------------------------------------------
! This SUBROUTINE read the input file for Hartree-Fock-Roothan program
!      *COMBINTATION OF TWO SUBROUTINES, CAN HAVE ISSUES/BUGS*
!                                                                       Temporary Code, I hope I hope I hope
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER ::                                        i, j, k, q, n_atoms, max_bf, max_prim, b_functions, l, n_max_zeta, info, m, n, p
INTEGER, DIMENSION(:), ALLOCATABLE ::             n_pri_bf, atomic_number, n_basisf_pa, n_prim
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE ::      xyz, zeta, d
REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE ::    zeta_atom, d_atom
CHARACTER(LEN=100) ::                             input_name
CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE ::    atom
CHARACTER(LEN=2) ::                               dummy

OPEN(11,FILE=input_name)
    READ(11,*)                                               !dummy line
    READ(11,*) n_atoms

    ALLOCATE(atom(n_atoms), atomic_number(n_atoms), n_basisf_pa(n_atoms), xyz(n_atoms,3))     !setting how many atoms

    READ(11,*)                                               !dummy line
    DO i = 1, n_atoms                      !chemical symbol, coordinates
        READ(11,*) atom(i), atomic_number(i), (xyz(i,j), j=1, 3) 
    ENDDO
    READ(11,*)                                               !dummy line
    READ(11,*) q                                            !over charge
    READ(11,*)                                               !dummy line
    READ(11,*) b_functions                     !number of basis function
    READ(11,*)                                               !dummy line
    READ(11,*) max_prim                               !max of primitives
    READ(11,FMT='(A)')                                       !dummy line

    n_max_zeta = max_prim * b_functions
    ALLOCATE(zeta(n_max_zeta,max_prim), d(n_max_zeta,max_prim), n_pri_bf(b_functions)); zeta(:,:) = 0.d0; d(:,:) = 0.d0

    m = 1; n = 0; p = 1
    DO i = 1, b_functions
        READ(11,*) k, dummy, k, l              !dummy line, (except 'l')
        READ(11,*) n_pri_bf(i)
        DO j = 1, n_pri_bf(i)
            READ(11,*) zeta(i,j), d(i,j)
        ENDDO
        IF ( l == m ) THEN
            n = n + 1
        ELSE
            n_basisf_pa(p) = n                             !SOULD BE
            n = 0; m = l; p = p + 1                        !REFACTORIZED
        ENDIF
    ENDDO
    n_basisf_pa(n_atoms) = n                               !SOULD 
    n_basisf_pa(:) = n_basisf_pa(:) + 1                    !BE
    n_basisf_pa(1) = n_basisf_pa(1) - 1                    !REFACTORIZED
CLOSE(11)

xyz(:,:) = 1.8897259886*xyz(:,:)                       !angstrom to bohr

info = 0                         !We are not seeing sysntaxis issues now
ENDSUBROUTINE NEW_READER
!***********************************************************************************************************************************
SUBROUTINE DISTANCES(n_atoms, b_functions, n_basisf_pa, xyz, xyz_bf, distance_atoms, distance_bf)
!-----------------------------------------------------------------------
!   This SUBROUTINE fill a matrix with the atom and basis f distances
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER ::                                           i, j, k, l, m, n, n_atoms, b_functions
INTEGER, DIMENSION(n_atoms) ::                       n_basisf_pa
REAL(KIND=8), DIMENSION(n_atoms, 3) ::               xyz
REAL(KIND=8), DIMENSION(b_functions, 3) ::           xyz_bf
REAL(KIND=8), DIMENSION(n_atoms, n_atoms) ::         distance_atoms
REAL(KIND=8), DIMENSION(b_functions, b_functions) :: distance_bf

DO i = 1, n_atoms                                            !distance between atoms
    DO j = i, n_atoms
        distance_atoms(i,j) = f_distance(xyz(i,:), xyz(j,:))
        distance_atoms(j,i) = distance_atoms(i,j)            !same distance between A & B than B & A
    ENDDO
ENDDO

k = 1                                                        !counter
DO i = 1, n_atoms                                            !distance between basis functions
   DO j = 1, n_basisf_pa(i)
       xyz_bf(k,:) = xyz(i,:)
       k = k + 1
   ENDDO
ENDDO

n = 1; m = 1                                                 !counters
DO i = 1, n_atoms
    DO j = 1, n_atoms
        DO k = 1, n_basisf_pa(i)
            DO l = 1, n_basisf_pa(j)
                distance_bf(m,n) = distance_atoms(i,j)
                n = n + 1                                    !only n_basisf_pa times
            ENDDO
            n = n - n_basisf_pa(j)                           !restart every time that we jump to the next row
            m = m + 1
        ENDDO
        n = n + n_basisf_pa(j)                               !start where we stop the last time with the new vale
        m = m - n_basisf_pa(i)                               !restart every time that we jump to the next column
    ENDDO
    n = 1                                                    !every row
    m = m + n_basisf_pa(i)
ENDDO

ENDSUBROUTINE DISTANCES
!***********************************************************************************************************************************
SUBROUTINE BOCInteraction(n_atoms, distance_mat, atomic_number, Ecoulomb, E_cou)
!-----------------------------------------------------------------------
!   This SUBROUTINE fill a matrix with $V_{n,n}$ at BO approximation
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER ::                                   i, j, n_atoms
INTEGER, DIMENSION(n_atoms) ::               atomic_number
REAL(KIND=8) ::                              E_cou
REAL(KIND=8), DIMENSION(n_atoms,n_atoms) ::  distance_mat, Ecoulomb

E_cou = 0.d0
DO i = 1, n_atoms-1
    DO j = i+1, n_atoms
        Ecoulomb(i,j) = f_coulomb(atomic_number(i), atomic_number(j), distance_mat(i,j))
        Ecoulomb(j,i) = Ecoulomb(i,j)
        E_cou = E_cou + Ecoulomb(i,j)
    ENDDO
ENDDO

ENDSUBROUTINE BOCInteraction
!***********************************************************************************************************************************
SUBROUTINE INTEGRALS_HF(b_functions, n_atoms, atomic_number, n_pri_bf, xyz, xyz_bf, distance_bf, max_prim, n_max_zeta, zeta, &
                        d, S, T, V, TwoEleInt, info)
!-----------------------------------------------------------------------
!   This SUBROUTINE compute integrals for Roothaan Hartree-Fock
!-----------------------------------------------------------------------
IMPLICIT NONE
INTEGER ::                                                i, j, k, l, m, n, o, p, n_atoms, max_prim, n_max_zeta, b_functions, &
                                                          one_int, two_int, info
INTEGER, DIMENSION(n_atoms) ::                            atomic_number, n_pri_bf
REAL (KIND=8) ::                                          zeta_escalar, xi_escalar, pi, OverLap, r_pq, Kkl, Kmn, phi_four, rho, &
                                                          zeta_prima, xi_prima, Kinetic, aux1, aux2, NucEleC, r_pm
REAL (KIND=8), DIMENSION(3) ::                            r_p, r_q
REAL (KIND=8), DIMENSION(n_atoms,3) ::                    xyz
REAL (KIND=8), DIMENSION(b_functions,3) ::                xyz_bf
REAL (KIND=8), DIMENSION(n_max_zeta,max_prim) ::          zeta, d
REAL (KIND=8), DIMENSION(b_functions,b_functions) ::      S, T, V, distance_bf
REAL (KIND=8), DIMENSION(b_functions,b_functions,b_functions,b_functions) :: TwoEleInt

pi   = DACOS(-1.d0)                                       !The best \pi value with your computer and fortran compilator ;)

!-----------------------------------------------------------------------
! How many one- & two-electron unique integrals the system has.
! Counting as Gaussian 10 years old.

one_int = b_functions*(b_functions+1)/2                   !how many overlap, kinetic and nuclear atraction integrals we have
two_int = one_int*(one_int+1)/2                           !how many two integrals we have

!-----------------------------------------------------------------------------------------------------------------------------------
!---- Overlap Integrals (S), Kinetic Integrals (T), Nucleus-Electron Coulomb Integrals (V) and Two Electron Integrals 

S(:,:) = 0.d0; T(:,:) = 0.d0; V(:,:) = 0.d0; TwoEleInt(:,:,:,:) = 0.d0      !Starting the sum

DO i = 1, b_functions
  DO j = i, b_functions
    DO k = 1, n_pri_bf(i)
      DO l = 1, n_pri_bf(j)
        zeta_escalar     = zeta(i,k) + zeta(j,l)
        xi_escalar       = zeta(i,k) * zeta(j,l)/zeta_escalar
        r_p(:)           = (zeta(i,k)*xyz_bf(i,:) + zeta(j,l)*xyz_bf(j,:))/zeta_escalar
        OverLap          = DEXP(-xi_escalar*distance_bf(i,j)**2) * (pi/zeta_escalar)**(3.d0/2.d0)
        S(i,j)           = S(i,j) + d(i,k)*d(j,l)*OverLap
        Kinetic          = xi_escalar*(3.d0-(2.d0*xi_escalar*distance_bf(i,j)**2))*OverLap
        T(i,j)           = T(i,j) + d(i,k)*d(j,l)*Kinetic
        NucEleC          = 0.d0
        DO m = 1, n_atoms
          r_pm      = f_distance(r_p(:), xyz(m,:))
          IF ( r_pm .LE. 1.D-5 ) THEN                                                ! distance of basis i, j and the nuclei is zero
            NucEleC = NucEleC +2.d0*atomic_number(m)*DSQRT(zeta_escalar/pi)*OverLap
          ELSE
            aux1    = zeta_escalar*r_pm**2                                             !zeta_escalar*f_distance(r_p(:), xyz(m,:))**2
            NucEleC = NucEleC +2.d0*atomic_number(m)*DSQRT(zeta_escalar/pi)*OverLap*f_cero(aux1)
          ENDIF
        ENDDO
        V(i,j)      = V(i,j) - d(i,k)*d(j,l)*NucEleC
        Kkl         = (DSQRT(2.d0)*(pi**(5.d0/4.d0))/zeta_escalar)*DEXP(-xi_escalar*distance_bf(i,j)**2)
        DO m = 1, b_functions
          DO n = m, b_functions
            IF ((j*(j+1)/2 +i).GE.(n*(n+1)/2 +m)) THEN
              DO o = 1, n_pri_bf(m)
                DO p = 1, n_pri_bf(n)
                  zeta_prima  = zeta(m,o) + zeta(n,p)
                  xi_prima    = zeta(m,o) * zeta(n,p)/zeta_prima
                  r_q(:)      = (zeta(m,o)*xyz_bf(m,:) + zeta(n,p)*xyz_bf(n,:))/zeta_prima
                  Kmn         = (DSQRT(2.d0)*(pi**(5.d0/4.d0))/zeta_prima)*DEXP(-xi_prima*distance_bf(m,n)**2)
                  rho         = (zeta_escalar*zeta_prima)/(zeta_escalar + zeta_prima)
                  r_pq        = f_distance(r_p(:), r_q(:))
                  IF (r_pq .LE. 1.D-5) THEN 
                    phi_four  = (Kkl*Kmn/DSQRT(zeta_escalar+zeta_prima))
                  ELSE
                    aux2      = rho*r_pq**2
                    phi_four  = (Kkl*Kmn/DSQRT(zeta_escalar+zeta_prima))*f_cero(aux2)
                  ENDIF
                  TwoEleInt(i,j,m,n) = TwoEleInt(i,j,m,n) + d(i,k)*d(j,l)*d(m,o)*d(n,p)*phi_four
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO
info = 0
!-----------------------------------------------------------------------
! Is the same interaction with the particle 1 and particle 2 as
! the particle 2 with 1 for S, V and T.
!         S_{ij} = S_{ji};   V_{ij} = V_{ji};   T_{ij} = T_{ji}
! For the two electron integrals the interaction symmetry, since our
! $\phi$ are real, is:
! (ij|kl) = (ij|lk) = (ji|kl) = (ji|lk) = (kl|ij) = (kl|ji) = (lk|ij)
!   ...   = (lk|ji)
! We just full fill the tensors with these redundancy.
! Also write on a temporal files.
!-----------------------------------------------------------------------
OPEN (501, FILE='./tmp/Overlap.int'); OPEN (502, FILE='./tmp/Potential.int')
OPEN (503, FILE='./tmp/Kinetic.int'); OPEN (504, FILE='./tmp/TwoElectron.int')  

WRITE(501,*) one_int; WRITE(502,*) one_int
WRITE(503,*) one_int; WRITE(504,*) two_int

DO i = 1, b_functions
  DO j = i, b_functions

    WRITE(501,FMT='(2I5,F24.16)') i, j, S(i,j); IF (ISNAN(S(i,j))) info = 501
    WRITE(502,FMT='(2I5,F24.16)') i, j, V(i,j); IF (ISNAN(V(i,j))) info = 502
    WRITE(503,FMT='(2I5,F24.16)') i, j, T(i,j); IF (ISNAN(T(i,j))) info = 503
    IF ((i.NE.j)) THEN
      S(j,i) = S(i,j);      V(j,i) = V(i,j);      T(j,i) = T(i,j)
    ENDIF

    DO k = 1, b_functions
      DO l = k, b_functions
        IF ((j*(j+1)/2 +i).GE.(l*(l+1)/2 +k)) THEN
          WRITE(504,FMT='(4I5,F24.16)') i, j, k, l, TwoEleInt(i,j,k,l)
          IF (ISNAN(TwoEleInt(i,j,l,k))) info = 504
          IF (.NOT.((i.EQ.j).AND.(j.EQ.k).AND.(k.EQ.l))) THEN
            TwoEleInt(i,j,l,k) = TwoEleInt(i,j,k,l)
            TwoEleInt(j,i,k,l) = TwoEleInt(i,j,k,l)
            TwoEleInt(j,i,l,k) = TwoEleInt(i,j,k,l)
            TwoEleInt(k,l,i,j) = TwoEleInt(i,j,k,l)
            TwoEleInt(k,l,j,i) = TwoEleInt(i,j,k,l)
            TwoEleInt(l,k,i,j) = TwoEleInt(i,j,k,l)
            TwoEleInt(l,k,j,i) = TwoEleInt(i,j,k,l)
          ENDIF
        ENDIF
      ENDDO
    ENDDO

  ENDDO 
  S(i,i) = 1.d0                               !S_{i,i} analitically is 1
ENDDO

CLOSE(501); CLOSE(502); CLOSE(503); CLOSE(504)  

ENDSUBROUTINE INTEGRALS_HF
!***********************************************************************
    FUNCTION f_distance(xyz, uvw) RESULT (distance)
    IMPLICIT NONE
    REAL (KIND=8) :: distance, dx, dy, dz
    REAL (KIND=8), DIMENSION(3) :: xyz, uvw 

    dx = xyz(1) - uvw(1)
    dy = xyz(2) - uvw(2)
    dz = xyz(3) - uvw(3)

    distance = DSQRT(dx*dx + dy*dy + dz*dz)
    ENDFUNCTION f_distance
!***********************************************************************
    FUNCTION f_coulomb(q,Z,distance) RESULT (coulomb)
    INTEGER :: q, Z
    REAL (KIND=8) :: coulomb, distance

    coulomb = q*Z/distance
    ENDFUNCTION f_coulomb
!***********************************************************************
    FUNCTION f_cero(x) RESULT (y)
    IMPLICIT NONE
    REAL (KIND=8) :: x, y, pi

    pi = DACOS(-1.d0)
    y  = DSQRT(pi/x)*DERF(DSQRT(x))/2.d0
    ENDFUNCTION f_cero
ENDPROGRAM HartreeFockRoothaan
