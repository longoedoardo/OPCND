MODULE ReferenceFunctions4D
USE TypesDef

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE cub_gausscheb_tens4D(deg, xyzw)
    !*******************************************************************************
    ! Calcola i punti e i pesi per l'integrazione numerica in 4D.
    ! Utilizza una griglia a prodotto tensoriale basata sui nodi di Gauss-Chebyshev.
    !*******************************************************************************
    INTEGER, INTENT(IN) :: deg
    REAL(dp), ALLOCATABLE, INTENT(OUT) :: xyzw(:,:)
    INTEGER :: n, n4, i, j, k, l, idx
    REAL(dp), ALLOCATABLE :: x(:), w(:)

    n = CEILING(REAL(deg + 1, dp) / 2.0_dp)
    n4 = n**4

    IF (ALLOCATED(xyzw)) DEALLOCATE(xyzw)
    ALLOCATE(xyzw(n4, 5))

    IF (ALLOCATED(x)) DEALLOCATE(x)
    ALLOCATE(x(n))

    IF (ALLOCATED(w)) DEALLOCATE(w)
    ALLOCATE(w(n))

    DO i = 1, n
        x(i) = COS( ((2.0 * REAL(i) - 1.0) * PI) / (2.0 * REAL(n)) )
        w(i) = PI / REAL(n)
    END DO

    idx = 1
    DO l = 1, n ! TAU
        DO k = 1, n ! Z
            DO j = 1, n ! Y
                DO i = 1, n ! X
                    xyzw(idx, 1) = x(i)
                    xyzw(idx, 2) = x(j)
                    xyzw(idx, 3) = x(k)
                    xyzw(idx, 4) = x(l)
                    xyzw(idx, 5) = w(i) * w(j) * w(k) * w(l)
                    idx = idx + 1
                END DO
            END DO
        END DO
    END DO

    END SUBROUTINE cub_gausscheb_tens4D







    SUBROUTINE mono_next_grlex(m, x)

    !**************************************************************************
    ! Calcola il monomio successivo in ordinamento grlex.
    !**************************************************************************
    
        INTEGER, INTENT(IN)    :: m
        INTEGER, INTENT(INOUT) :: x(m)

        INTEGER                :: i, j, t, im1

        ! Mi assicuro che M >= 1
        IF (m < 1) THEN
            WRITE(*,*) ' '
            WRITE(*,*) 'MONO_NEXT_GRLEX - Fatal error!'
            WRITE(*,*) '  M < 1, m = ', m
            ERROR STOP 'MONO_NEXT_GRLEX - Fatal error!'
        END IF

        ! Mi assicuro che X(I) => 0
        DO i = 1, m
            IF (x(i) < 0) THEN
                WRITE(*,*) ' '
                WRITE(*,*) 'MONO_NEXT_GRLEX - Fatal error!'
                WRITE(*,*) '  X(I) < 0 presso indice ', i
                ERROR STOP 'MONO_NEXT_GRLEX - Fatal error!'
            END IF
        END DO

        ! Trova I, l'indice dell'elemento più a destra diverso da zero
        i = 0
        DO j = m, 1, -1
            IF (x(j) > 0) THEN
                i = j
                EXIT ! Sostituisce il break di Matlab
            END IF
        END DO

        ! Core logico dell'avanzamento grlex
        IF (i == 0) THEN
            x(m) = 1
            RETURN
        ELSE IF (i == 1) THEN
            t = x(1) + 1
            im1 = m
        ELSE ! Equivalente a i>1
            t = x(i)
            im1 = i - 1
        END IF

        x(i) = 0
        x(im1) = x(im1) + 1
        x(m) = x(m) + t - 1

    END SUBROUTINE mono_next_grlex

   SUBROUTINE dCHEBVAND(deg, X, duples, V)
      INTEGER, INTENT(IN)                         :: deg
      REAL(dp), INTENT(IN)                        :: X(:,:)
      INTEGER, INTENT(IN)                         :: duples(:,:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)          :: V(:,:)

      INTEGER                                     :: m, d, ncols, k, col, g
      REAL(dp), ALLOCATABLE                       :: T_dim(:,:)

      m = SIZE(X, 1)
      d = SIZE(X, 2)
      ncols = SIZE(duples, 1)

      ALLOCATE(V(m, ncols))
      ALLOCATE(T_dim(m, 0:deg))
      V = 1.0_dp

      DO k = 1, d
         T_dim(:,0) = 1.0_dp
         IF (deg >= 1) THEN
            T_dim(:,1) = X(:,k)
            DO g = 2, deg
               T_dim(:,g) = 2.0_dp * X(:,k) * T_dim(:,g-1) - T_dim(:,g-2)
            END DO
         END IF
         DO col = 1, ncols
            V(:,col) = V(:,col) * T_dim(:, duples(col,k))
         END DO
      END DO
   END SUBROUTINE dCHEBVAND

   SUBROUTINE tenscheb_norm2sq(chebyshev_indices, coeffs)
      INTEGER, INTENT(IN)                         :: chebyshev_indices(:,:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)          :: coeffs(:)

      INTEGER                                     :: i, j, n, d, non_zero_counter

      n = SIZE(chebyshev_indices, 1)
      d = SIZE(chebyshev_indices, 2)
      ALLOCATE(coeffs(n))

      DO i = 1, n
         non_zero_counter = 0
         DO j = 1, d
            IF (chebyshev_indices(i,j) /= 0) non_zero_counter = non_zero_counter + 1
         END DO
         coeffs(i) = PI**REAL(d, dp) / 2.0_dp**REAL(non_zero_counter, dp)
      END DO
   END SUBROUTINE tenscheb_norm2sq

   SUBROUTINE scale_rule(XYZTW_tens_ref, dbox, XYZTW_tens)
      REAL(dp), INTENT(IN)                        :: XYZTW_tens_ref(:,:)
      REAL(dp), INTENT(IN)                        :: dbox(2,4)
      REAL(dp), ALLOCATABLE, INTENT(OUT)          :: XYZTW_tens(:,:)

      INTEGER                                     :: j
      REAL(dp)                                    :: center, half_width

      ALLOCATE(XYZTW_tens(SIZE(XYZTW_tens_ref,1), SIZE(XYZTW_tens_ref,2)))
      DO j = 1, 4
         center = 0.5_dp * (dbox(1,j) + dbox(2,j))
         half_width = 0.5_dp * (dbox(2,j) - dbox(1,j))
         XYZTW_tens(:,j) = center + half_width * XYZTW_tens_ref(:,j)
      END DO
      XYZTW_tens(:,5) = XYZTW_tens_ref(:,5)
   END SUBROUTINE scale_rule

END MODULE ReferenceFunctions4D
