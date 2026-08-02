MODULE CubaCheap

   USE PrismQuadratureGJ
   USE PolyhedronMesh
   USE TypesDef

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE dCHEBVAND(deg, X, duples, V)

   !*******************************************************************************
   ! Calcola la matrice di Chebyshev-Vandermonde per il grado totale "deg"
   ! sui punti X (m x d), usando la base tensoriale indicizzata da "duples"
   ! (N x d), shiftata sul cubo unitario [-1,1]^d.
   !*******************************************************************************

      IMPLICIT NONE

      ! Input variables
      INTEGER, INTENT(IN)                                   :: deg
      REAL(dp), INTENT(IN)                                  :: X(:,:)
      INTEGER, INTENT(IN)                                   :: duples(:,:) 

      ! Output variables
      REAL(dp), ALLOCATABLE, INTENT(OUT)                    :: V(:,:)

      ! Local variables
      INTEGER                                               :: m, d, N, i, k, col, g
      REAL(dp)                                              :: min_val, max_val
      REAL(dp), ALLOCATABLE                                 :: dbox(:,:)
      REAL(dp), ALLOCATABLE                                 :: map(:,:)
      REAL(dp), ALLOCATABLE                                 :: T_dim(:,:)
      
      m = SIZE(X,1)
      d = SIZE(X,2)
      N = SIZE(duples,1)

      ! dbox: cubo unitario [-1,1] per ogni dimensione (4D nel caso d=4)
      ALLOCATE(dbox(2,d))
      dbox(1,:) = -1.0_dp
      dbox(2,:) =  1.0_dp

      ! Porto tutto in [-1,1] in modo vettorizzato per colonne
      ALLOCATE(map(m,d))
      DO i = 1, d
         min_val = dbox(1,i)
         max_val = dbox(2,i)
         map(:,i) = (2.0_dp * X(:,i) - max_val - min_val) / (max_val - min_val)
      END DO

      ALLOCATE(V(m,N)) ! Inizializzo la matrice di Vandermonde
      V = 1.0_dp

      ! Pre-allochiamo la matrice temporanea per i polinomi 1D (m x deg+1)
      ALLOCATE(T_dim(m, 0:deg))

      ! Ciclo sulle dimensioni (x, y, z, tau) ispirato alla versione MATLAB 4D
      DO k = 1, d
         T_dim(:,0) = 1.0_dp ! T0(x) = 1

         IF (deg >= 1) THEN
            T_dim(:,1) = map(:,k) ! T1(x) = x

            DO g = 2, deg
               T_dim(:,g) = 2.0_dp * map(:,k) * T_dim(:,g-1) - T_dim(:,g-2)
            END DO
         END IF

         ! Prodotto tensoriale, aggiorna ogni colonna di V in base ai multi-indici duples
         DO col = 1, N
            V(:,col) = V(:,col) * T_dim(:, duples(col,k))
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

      ! Input variables 
      INTEGER, INTENT(in)             :: chebyshev_indices(:, :)

      ! Output variables
      REAL(dp), ALLOCATABLE, INTENT(out) :: coeffs(:)

      ! Local variables
      INTEGER                         :: i, n, d, non_zero_counter

      n = size(chebyshev_indices, 1)
      d = size(chebyshev_indices, 2) ! Sarà 4 per il caso 4D

      IF(ALLOCATED(coeffs))  DEALLOCATE(coeffs)
      ALLOCATE(coeffs(n))

      DO i = 1, n
         non_zero_counter = 0
         IF (d >= 1 .AND. chebyshev_indices(i, 1) /= 0) non_zero_counter = non_zero_counter + 1
         IF (d >= 2 .AND. chebyshev_indices(i, 2) /= 0) non_zero_counter = non_zero_counter + 1
         IF (d >= 3 .AND. chebyshev_indices(i, 3) /= 0) non_zero_counter = non_zero_counter + 1
         IF (d >= 4 .AND. chebyshev_indices(i, 4) /= 0) non_zero_counter = non_zero_counter + 1
        
         coeffs(i) = (PI**REAL(d, dp)) / (2.0_dp**REAL(non_zero_counter, dp))
      END DO

   END SUBROUTINE tenscheb_norm2sq

   
   SUBROUTINE chebyshev_moments_polyhedron_4D(vertici_4D, Hyperfacets, ade, chebyshev_indices, dbox, moments)

   !*******************************************************************************
   ! Calcola i momenti di Chebyshev su un iperpoliedro 4D utilizzando il teorema 
   ! della divergenza esteso al caso 4D. Le facce sono prismi spaziotemporali.
   !*******************************************************************************

      IMPLICIT NONE

      ! Input variables
      REAL(dp), INTENT(IN)                  :: vertici_4D(:, :)
      TYPE(t_hyperfacet), INTENT(IN)        :: Hyperfacets(:)
      INTEGER, INTENT(IN)                   :: ade
      INTEGER, INTENT(IN)                   :: chebyshev_indices(:, :)
      REAL(dp), INTENT(IN)                  :: dbox(:, :)

      ! Output variables
      REAL(dp), INTENT(OUT), ALLOCATABLE    :: moments(:)

      ! Local variables
      INTEGER                               :: num_indici, n_facce, k, ii, num_pts_dummy
      REAL(dp)                              :: scale_x_dummy
      
      ! Variabili per la quadratura e gestione del prisma 4D
      REAL(dp), ALLOCATABLE                 :: XI_ref(:), ETA_ref(:), T_ref(:), W_ref(:)
      REAL(dp), ALLOCATABLE                 :: XYZTW(:, :), WV_CUB(:)
      REAL(dp), ALLOCATABLE                 :: V_prism(:, :)
      REAL(dp), ALLOCATABLE                 :: chebyshev_moms(:, :)
      REAL(dp), ALLOCATABLE                 :: moms_facet_raw(:)
      INTEGER, ALLOCATABLE                  :: nodes_id(:)

      num_indici = SIZE(chebyshev_indices, 1)
      n_facce = SIZE(Hyperfacets)

      ALLOCATE(chebyshev_moms(num_indici, n_facce))
      chebyshev_moms = 0.0_dp
      
      DO k = 1, n_facce
         nodes_id = Hyperfacets(k)%Vertices_ID
         ALLOCATE(V_prism(SIZE(nodes_id), 4))
         DO ii = 1, SIZE(nodes_id)
            V_prism(ii, :) = vertici_4D(nodes_id(ii), :)
         END DO

         CALL reference_prism_quadrature(ade, dbox, XI_ref, ETA_ref, T_ref, W_ref, scale_x_dummy, num_pts_dummy)
         CALL PrismQuad4D(V_prism, XI_ref, ETA_ref, T_ref, W_ref, XYZTW, WV_CUB)

         CALL cubature_tens_chebyshev_facet_4D(XYZTW, WV_CUB, ade, chebyshev_indices, dbox, moms_facet_raw)

         chebyshev_moms(:, k) = moms_facet_raw

         DEALLOCATE(V_prism, XI_ref, ETA_ref, T_ref, W_ref, XYZTW, WV_CUB, moms_facet_raw, nodes_id)
      END DO

      ! Somma i contributi di tutte le facce
      ALLOCATE(moments(num_indici))
      moments = SUM(chebyshev_moms, DIM=2)

      DEALLOCATE(chebyshev_moms)

   END SUBROUTINE chebyshev_moments_polyhedron_4D

   SUBROUTINE cubature_tens_chebyshev_facet_4D(XYZTW, WV_CUB, ade, chebyshev_indices, dbox, chebyshev_moms)

   !***************************************************************
   ! Calcola i momenti di Chebyshev su un'iperfaccia 4D
   ! Adattato alla versione MATLAB 4D
   !***************************************************************

      IMPLICIT NONE

      ! Input variables
      REAL(dp), INTENT(IN)    :: XYZTW(:, :)        
      REAL(dp), INTENT(IN)    :: WV_CUB(:)         
      INTEGER, INTENT(IN)    :: chebyshev_indices(:, :) 
      REAL(dp), INTENT(IN)    :: dbox(:, :) ! Matrice 2x4 per il 4D

      ! Output variables
      REAL(dp), ALLOCATABLE, INTENT(OUT)   :: chebyshev_moms(:)  

      ! Local variables
      INTEGER                :: n, m, ade, c, n_idx, ii
      REAL(dp)                :: scale_x
      REAL(dp), ALLOCATABLE   :: X_ref(:, :)
      REAL(dp), ALLOCATABLE   :: Tx(:, :), Ty(:, :), Tz(:, :), Tt(:, :)
      REAL(dp), ALLOCATABLE   :: IntTx(:, :)
      REAL(dp), ALLOCATABLE   :: V_poly(:, :)
      REAL(dp), ALLOCATABLE   :: w(:)

      n = SIZE(XYZTW, 1)
      m = SIZE(chebyshev_indices, 1)

      ! Allokazione e scalatura dei punti sul dominio [-1, 1]^4
      ALLOCATE(X_ref(n, 4))
      DO c = 1, 4
         X_ref(:, c) = (2.0_dp * XYZTW(:, c) - (dbox(1, c) + dbox(2, c))) / (dbox(2, c) - dbox(1, c))
      END DO

      w = WV_CUB(:)

      ! Valutazione delle basi di Chebyshev 1D
      ALLOCATE(Tx(n, ade + 2))
      ALLOCATE(Ty(n, ade + 1))
      ALLOCATE(Tz(n, ade + 1))
      ALLOCATE(Tt(n, ade + 1))

      CALL chebpolys(ade + 1, X_ref(:, 1), Tx)
      CALL chebpolys(ade,     X_ref(:, 2), Ty)
      CALL chebpolys(ade,     X_ref(:, 3), Tz)
      CALL chebpolys(ade,     X_ref(:, 4), Tt)

      ! Integrazione analitica (Primitiva) lungo la direzione X
      scale_x = (dbox(2, 1) - dbox(1, 1)) / 2.0_dp
      ALLOCATE(IntTx(n, SIZE(Tx, 2)))
      IntTx = 0.0_dp

      IntTx(:, 1) = Tx(:, 2)
      IF (ade >= 1) THEN
         IntTx(:, 2) = 0.25_dp * Tx(:, 3) + 1.0_dp
      END IF

      DO n_idx = 2, ade
         IntTx(:, n_idx + 1) = 0.5_dp * (Tx(:, n_idx + 2) / REAL(n_idx + 1, dp) - Tx(:, n_idx) / REAL(n_idx - 1, dp))
      END DO

      IntTx = IntTx * scale_x

      ! Combinazione tensoriale guidata dagli indici
      ALLOCATE(V_poly(n, m))
      DO ii = 1, m
         V_poly(:, ii) = IntTx(:, chebyshev_indices(ii, 1) + 1) * &
                         Ty(:,     chebyshev_indices(ii, 2) + 1) * &
                         Tz(:,     chebyshev_indices(ii, 3) + 1) * &
                         Tt(:,     chebyshev_indices(ii, 4) + 1)
      END DO

      IF (ALLOCATED(chebyshev_moms)) DEALLOCATE(chebyshev_moms)
      ALLOCATE(chebyshev_moms(m))

      DO ii = 1, m
         chebyshev_moms(ii) = SUM(w * V_poly(:, ii))
      END DO

      DEALLOCATE(X_ref, Tx, Ty, Tz, Tt, IntTx, V_poly)

   END SUBROUTINE cubature_tens_chebyshev_facet_4D

   SUBROUTINE chebpolys(deg, x, T)

   !*******************************************************************************
   ! Calcola i polinomi di Chebyshev di prima specie T_j(x) fino al grado "deg" 
   ! valutati su un vettore di punti "x", utilizzando la formula di ricorrenza 
   ! a tre termini e memorizzando i risultati in una matrice T.
   !*******************************************************************************

      IMPLICIT NONE

      ! Input variables
      INTEGER, INTENT(IN)    :: deg
      REAL(dp), INTENT(IN)    :: x(:)

      ! Output variables
      REAL(dp), INTENT(OUT)   :: T(:, :)

      ! Local variables
      INTEGER                :: n, j
      REAL(dp), ALLOCATABLE   :: t0(:), t1(:), t2(:)

      n = SIZE(x)

      ALLOCATE(t0(n), t1(n), t2(n))

      t0 = 1.0_dp
      T(:, 1) = t0

      IF (deg >= 1) THEN
         t1 = x
         T(:, 2) = t1
      END IF

      DO j = 2, deg
         t2 = 2.0_dp * x * t1 - t0
         T(:, j + 1) = t2
         t0 = t1
         t1 = t2
      END DO

      DEALLOCATE(t0, t1, t2)

   END SUBROUTINE chebpolys

   SUBROUTINE scale_rule(XYZTW_tens_ref, dbox, XYZW_tens)

   !*******************************************************************************
   ! Scala le coordinate dei punti di quadratura di una regola tensoriale 
   ! dal dominio di riferimento standard [-1, 1]^4 al bounding box geometrico 
   ! effettivo in 4D, lasciando inalterati i pesi di cubatura.
   !*******************************************************************************

      IMPLICIT NONE

      ! Input variables
      REAL(dp), INTENT(IN)  :: XYZTW_tens_ref(:,:)
      REAL(dp), INTENT(IN)  :: dbox(:,:) ! Matrice 2x4 per il 4D [min; max]

      ! Output variables
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: XYZW_tens(:,:)
      
      ! Local variables
      INTEGER :: n_rows, n_cols, j
      REAL(dp), ALLOCATABLE :: centro(:), semi_ampiezza(:)

      n_rows = SIZE(XYZTW_tens_ref, 1)
      n_cols = SIZE(XYZTW_tens_ref, 2)

      IF (ALLOCATED(XYZW_tens)) DEALLOCATE(XYZW_tens)
      ALLOCATE(XYZW_tens(n_rows, n_cols))

      ALLOCATE(centro(4), semi_ampiezza(4))

      ! Calcolo vettoriale del centro e della semi-ampiezza per i 4 assi
      DO j = 1, 4
         centro(j) = (dbox(1, j) + dbox(2, j)) / 2.0_dp
         semi_ampiezza(j) = (dbox(2, j) - dbox(1, j)) / 2.0_dp
      END DO

      ! Scalatura delle prime 4 colonne (coordinate spaziotemporali)
      DO j = 1, 4
         XYZW_tens(:, j) = centro(j) + semi_ampiezza(j) * XYZTW_tens_ref(:, j)
      END DO

      ! Se la matrice di riferimento ha 5 colonne, copiamo i pesi nella quinta
      IF (n_cols >= 5) THEN
         XYZW_tens(:, 5) = XYZTW_tens_ref(:, 5)
      END IF

      DEALLOCATE(centro, semi_ampiezza)

   END SUBROUTINE scale_rule

END MODULE CubaCheap