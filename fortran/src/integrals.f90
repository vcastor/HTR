!-----------------------------------------------------------------------
!   This SUBROUTINE compute integrals for Roothaan Hartree-Fock
!-----------------------------------------------------------------------
! bf        :: number of basis functions
! na        :: number of atoms
! a_num     :: atomic numbers
! n_pri_bf  :: number of primitive per basis functions
! xyz       :: coordinates of atoms
! xyz_bf    :: coordinates of basis functions
! dis2_bf   :: squared distance between basis functions
! max_prim  :: maximum number of primitives
! max_zeta  :: maximum number of zeta
! zeta      :: zeta values
! d         :: d values NORMALIZED
! S         :: Overlap integral
! T         :: Kinetic integral
! V         :: Potential integral
! TwoEleInt :: Two electron integral
! info      :: info flag
!-----------------------------------------------------------------------
SUBROUTINE INTEGRALS(bf, na, a_num, n_pri_bf, xyz, xyz_bf, dis2_bf, max_prim, max_zeta, zeta, d, S, T, V, TwoEleInt, info)
IMPLICIT NONE
INCLUDE 'parameters.h'
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                                           :: i, j, k, l, m, n, o, p, na, max_prim, max_zeta, bf, one_int, two_int, info
INTEGER, DIMENSION(na)                            :: a_num
INTEGER, DIMENSION(bf)                            :: n_pri_bf
REAL (KIND=8)                                     :: OverLap, r_pq2, Kkl, Kmn, phi_four, rho
REAL (KIND=8)                                     :: Kinetic, NucEleC, r_pm2
REAL (KIND=8), EXTERNAL                           :: f_distance2, f_cero
REAL (KIND=8), DIMENSION(3)                       :: r_p, r_q
REAL (KIND=8), DIMENSION(na,3)                    :: xyz
REAL (KIND=8), DIMENSION(bf,3)                    :: xyz_bf
REAL (KIND=8), DIMENSION(max_zeta,max_prim)       :: zeta, d
REAL (KIND=8), DIMENSION(bf,bf)                   :: S, T, V, dis2_bf
REAL (KIND=8), DIMENSION(bf,bf,bf,bf)             :: TwoEleInt
REAL (KIND=8), DIMENSION(bf,max_prim,bf,max_prim) :: zeta_escalar, xi_escalar
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------
! To not recompute the \xi and \zeta values in the two electron
! integrals, firstly we compute all values and then just read it from
! memory.

DO i = 1, bf
  DO j = i, bf
    DO k = 1, n_pri_bf(i)
      DO l = 1, n_pri_bf(j)
        zeta_escalar(i,k,j,l) = zeta(i,k) + zeta(j,l)
        xi_escalar(i,k,j,l)   = zeta(i,k) * zeta(j,l)/zeta_escalar(i,k,j,l)
      ENDDO
    ENDDO
  ENDDO
ENDDO

!-----------------------------------------------------------------------
! How many one- & two-electron unique integrals the system has.
! Counting as Gaussian 10 years old.

one_int = bf*(bf+1)/2                    !how many overlap, kinetic and nuclear atraction integrals we have
two_int = one_int*(one_int+1)/2          !how many two integrals we have

!-----------------------------------------------------------------------------------------------------------------------------------
!---- Overlap Integrals (S), Kinetic Integrals (T), Nucleus-Electron Coulomb Integrals (V) and Two Electron Integrals 

S(:,:) = 0.d0; T(:,:) = 0.d0; V(:,:) = 0.d0; TwoEleInt(:,:,:,:) = 0.d0                                             !Starting the sum

DO i = 1, bf
  DO j = i, bf
    DO k = 1, n_pri_bf(i)
      DO l = 1, n_pri_bf(j)
        r_p(:)   = (zeta(i,k)*xyz_bf(i,:) + zeta(j,l)*xyz_bf(j,:))/zeta_escalar(i,k,j,l)
        OverLap  = DEXP(-xi_escalar(i,k,j,l)*dis2_bf(i,j)) * (pi/zeta_escalar(i,k,j,l))**(3.d0/2.d0)
        S(i,j)   = S(i,j) + d(i,k)*d(j,l)*OverLap
        Kinetic  = xi_escalar(i,k,j,l)*(3.d0-(2.d0*xi_escalar(i,k,j,l)*dis2_bf(i,j)))*OverLap
        T(i,j)   = T(i,j) + d(i,k)*d(j,l)*Kinetic
        NucEleC  = 0.d0
        DO m = 1, na
          r_pm2 = f_distance2(r_p(:), xyz(m,:))
          !IF ( r_pm .LE. 1.D-5 )                                                 ! distance of basis i, j and the nuclei is zero
          NucEleC = NucEleC +2.d0*a_num(m)*DSQRT(zeta_escalar(i,k,j,l)/pi)*OverLap*f_cero(zeta_escalar(i,k,j,l)*r_pm2)
        ENDDO
        V(i,j) = V(i,j) - d(i,k)*d(j,l)*NucEleC
        Kkl    = (DSQRT(2.d0)*(pi**(5.d0/4.d0))/zeta_escalar(i,k,j,l))*DEXP(-xi_escalar(i,k,j,l)*dis2_bf(i,j))
        DO m = 1, bf
          DO n = m, bf
            IF ((j*(j+1)/2 +i).GE.(n*(n+1)/2 +m)) THEN
              DO o = 1, n_pri_bf(m)
                DO p = 1, n_pri_bf(n)
                  r_q(:)    = (zeta(m,o)*xyz_bf(m,:) + zeta(n,p)*xyz_bf(n,:))/zeta_escalar(m,o,n,p)
                  Kmn       = (DSQRT(2.d0)*(pi**(5.d0/4.d0))/zeta_escalar(m,o,n,p))*DEXP(-xi_escalar(m,o,n,p)*dis2_bf(m,n))
                  rho       = (zeta_escalar(i,k,j,l)*zeta_escalar(m,o,n,p))/(zeta_escalar(i,k,j,l) + zeta_escalar(m,o,n,p))
                  r_pq2     = f_distance2(r_p(:), r_q(:))
                  phi_four  = (Kkl*Kmn/DSQRT(zeta_escalar(i,k,j,l)+zeta_escalar(m,o,n,p)))*f_cero(rho*r_pq2)
                  !IF (r_pq .LE. 1.D-5) THEN 
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
!OPEN (501, FILE='./tmp/Overlap.int'); OPEN (502, FILE='./tmp/Potential.int')
!OPEN (503, FILE='./tmp/Kinetic.int'); OPEN (504, FILE='./tmp/TwoElectron.int')  

!WRITE(501,*) one_int; WRITE(502,*) one_int
!WRITE(503,*) one_int; WRITE(504,*) two_int

DO i = 1, bf
  DO j = i, bf

!    WRITE(501,FMT='(2I5,F24.16)') i, j, S(i,j); IF (ISNAN(S(i,j))) info = 501
!    WRITE(502,FMT='(2I5,F24.16)') i, j, V(i,j); IF (ISNAN(V(i,j))) info = 502
!    WRITE(503,FMT='(2I5,F24.16)') i, j, T(i,j); IF (ISNAN(T(i,j))) info = 503
    IF ((i.NE.j)) THEN
      S(j,i) = S(i,j);      V(j,i) = V(i,j);      T(j,i) = T(i,j)
    ENDIF

    DO k = 1, bf
      DO l = k, bf
        IF ((j*(j+1)/2 +i).GE.(l*(l+1)/2 +k)) THEN
!          WRITE(504,FMT='(4I5,F24.16)') i, j, k, l, TwoEleInt(i,j,k,l)
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

!CLOSE(501); CLOSE(502); CLOSE(503); CLOSE(504)  

ENDSUBROUTINE INTEGRALS
