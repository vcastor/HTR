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
