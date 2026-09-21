MODULE ReferenceFunctions
USE TypesDef

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE cub_gausscheb_tens4D(deg, xyzw)
    !*******************************************************************************
    ! Calcola i punti e i pesi per l'integrazione numerica in 4D.
    ! Utilizza una griglia a prodotto tensoriale basata sui nodi di Gauss-Chebyshev.
    !*******************************************************************************

    IMPLICIT NONE
    !*******************************************************************************
    ! Argomenti
    !*******************************************************************************
    INTEGER, INTENT(IN)                                         :: deg
    REAL(dp), ALLOCATABLE, INTENT(OUT)                          :: xyzw(:,:)
    !*******************************************************************************
    ! Variabili locali
    !*******************************************************************************
    INTEGER                                                     :: n, n4, i, j, k, l, idx
    REAL(dp), ALLOCATABLE :: x(:), w(:)
    !*******************************************************************************
    
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

        !*******************************************************************************
        ! Calcola il monomio successivo in ordinamento grlex.
        !*******************************************************************************

        IMPLICIT NONE

        !*******************************************************************************
        ! Argomenti
        !*******************************************************************************
        INTEGER, INTENT(IN)                                 :: m
        INTEGER, INTENT(INOUT)                              :: x(m)
        !*******************************************************************************
        ! Variabili locali
        !*******************************************************************************
        INTEGER                                             :: i, j, t, im1
        !*******************************************************************************

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
                EXIT
            END IF
        END DO

        IF (i == 0) THEN
            x(m) = 1
            RETURN
        ELSE IF (i == 1) THEN
            t = x(1) + 1
            im1 = m
        ELSE 
            t = x(i)
            im1 = i - 1
        END IF

        x(i) = 0
        x(im1) = x(im1) + 1
        x(m) = x(m) + t - 1

    END SUBROUTINE mono_next_grlex


   SUBROUTINE dCHEBVAND(deg, X, duples, V)

        !*******************************************************************************
        ! Calcola la matrice di Chebyshev-Vandermonde per il grado totale "deg"
        ! sui punti X (m x d), usando la base tensoriale indicizzata da "duples"
        ! (N x d), shiftata sul cubo unitario [-1,1]^d.
        !*******************************************************************************

        IMPLICIT NONE

        !*******************************************************************************
        ! Argomenti
        !*******************************************************************************
        INTEGER, INTENT(IN)                                   :: deg
        REAL(dp), INTENT(IN)                                  :: X(:,:)
        INTEGER, INTENT(IN)                                   :: duples(:,:) 
        REAL(dp), ALLOCATABLE, INTENT(OUT)                    :: V(:,:)
        !*******************************************************************************
        ! Variabili locali
        !*******************************************************************************
        INTEGER                                               :: m, d, N, i, k, col, g
        REAL(dp)                                              :: min_val, max_val
        REAL(dp), ALLOCATABLE                                 :: dbox(:,:)
        REAL(dp), ALLOCATABLE                                 :: map(:,:)
        REAL(dp), ALLOCATABLE                                 :: T_dim(:,:,:)
        !*******************************************************************************

        m = SIZE(X,1)
        d = SIZE(X,2)
        N = SIZE(duples,1)

        ! dbox: cubo unitario [-1,1] per ogni dimensione 
        ALLOCATE(dbox(2,d))
        dbox(1,:) = -1.0_dp
        dbox(2,:) =  1.0_dp

        ! Porto tutto in [-1,1]
        ALLOCATE(map(m,d))
        DO i = 1, d
            min_val = dbox(1,i)
            max_val = dbox(2,i)
            map(:,i) = (2.0_dp * X(:,i) - max_val - min_val) / (max_val - min_val)
        END DO

        ALLOCATE(V(m,N)) ! Inizializzo la matrice di Vandermonde
        V = 1.0_dp

        ! Costruisco i polinomi di Chebyshev 1D per ciascuna dimensione tramite ricorrenza
        ALLOCATE(T_dim(m, 0:deg, d))
        DO k = 1, d
            T_dim(:,0,k) = 1.0_dp ! T0(x) = 1

            IF (deg >= 1) THEN
                T_dim(:,1,k) = map(:,k) ! T1(x) = x

                DO g = 2, deg
                T_dim(:,g,k) = 2.0_dp * map(:,k) * T_dim(:,g-1,k) - T_dim(:,g-2,k)
                END DO
            END IF

            ! Prodotto tensoriale, aggiorna ogni colonna di V
            DO col = 1, N
                V(:,col) = V(:,col) * T_dim(:, duples(col,k), k)
            END DO

        END DO
        DEALLOCATE(dbox, map, T_dim)
  
    END SUBROUTINE dCHEBVAND


   SUBROUTINE tenscheb_norm2sq(chebyshev_indices, coeffs)

        !*******************************************************************************
        ! Calcola la norma quadra (L^2 pesata) di ciascun polinomio ortogonale 
        ! di Chebyshev tensoriale in 4D, utilizzato come fattore di normalizzazione.
        !*******************************************************************************

        IMPLICIT NONE

        !*******************************************************************************
        ! Argomenti
        !*******************************************************************************
        INTEGER, INTENT(IN)             :: chebyshev_indices(:, :)
        REAL(dp), ALLOCATABLE, INTENT(OUT) :: coeffs(:)
        !*******************************************************************************
        ! Variabili locali
        !*******************************************************************************
        INTEGER                         :: i, n
        INTEGER                         :: non_zero_counter
        !*******************************************************************************

        n = size(chebyshev_indices, 1)

        IF(ALLOCATED(coeffs))  DEALLOCATE(coeffs)
        ALLOCATE(coeffs(n))

        DO i = 1, n

            non_zero_counter = 0

            IF (chebyshev_indices(i, 1) /= 0) &
                non_zero_counter = non_zero_counter + 1

            IF (chebyshev_indices(i, 2) /= 0) &
                non_zero_counter = non_zero_counter + 1

            IF (chebyshev_indices(i, 3) /= 0) &
                non_zero_counter = non_zero_counter + 1

            IF (chebyshev_indices(i, 4) /= 0) &
                non_zero_counter = non_zero_counter + 1

            coeffs(i) = (PI**4) / (2.0_dp**non_zero_counter)

        END DO

    END SUBROUTINE tenscheb_norm2sq

END MODULE ReferenceFunctions
