MODULE CubatureTensor

   USE TypesDef
   USE ReferenceFunctions
   USE PrismQuadrature

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE chebyshev_moments_polyhedron_4D(vertices_4d, facets, ade, &
                                           chebyshev_indices, dbox, method, moments)

      !*******************************************************************************
      !
      ! Calcola i momenti dei polinomi di Chebyshev sul poliedro quadridimensionale
      ! delimitato dalla mesh triangolare mobile definita da vertices_4d e facets.
      !
      ! I momenti volumetrici quadridimensionali sono trasformati, mediante il teorema
      ! della divergenza, in integrali di superficie sulle ipersuperfici 3D laterali
      ! del prisma triangolare spazio-temporale.
      !
      ! La quadratura sul prisma puo' essere effettuata mediante uno dei seguenti
      ! metodi:
      !
      !   'DCC'  : quadratura di Dunavant sul triangolo e Clenshaw-Curtis nel tempo;
      !   'DGL'  : quadratura di Dunavant sul triangolo e Gauss-Legendre nel tempo;
      !   'GJCC' : quadratura di Gauss-Jacobi sul triangolo e Clenshaw-Curtis nel tempo;
      !   'GJL'  : quadratura di Gauss-Jacobi sul triangolo e Gauss-Legendre nel tempo.
      !
      ! Per ciascuna faccia triangolare viene costruito il prisma 3D spazio-temporale
      ! mediante interpolazione lineare tra il triangolo iniziale e quello finale.
      !
      ! Il contributo ai momenti viene quindi calcolato sull'ipersuperficie laterale
      ! del prisma. Le facce vengono elaborate indipendentemente e parallelamente
      ! mediante OpenMP.
      !
      !*******************************************************************************

      USE OMP_LIB, ONLY: omp_get_thread_num

      IMPLICIT NONE

      !*******************************************************************************
      ! Argomenti
      !*******************************************************************************

      REAL(dp), INTENT(IN)                               :: vertices_4d(:,:)
      INTEGER, INTENT(IN)                                :: facets(:,:)
      INTEGER, INTENT(IN)                                :: ade
      INTEGER, INTENT(IN)                                :: chebyshev_indices(:,:)
      REAL(dp), INTENT(IN)                               :: dbox(2,4)
      CHARACTER(LEN=*), INTENT(IN)                       :: method
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: moments(:)

      !*******************************************************************************
      ! Variabili locali
      !*******************************************************************************

      INTEGER :: num_indici
      INTEGER :: n_facce
      INTEGER :: n_vertici
      INTEGER :: f
      INTEGER :: n1
      INTEGER :: n2
      INTEGER :: n3

      REAL(dp) :: V_prism(6,4)

      REAL(dp), ALLOCATABLE :: XI_ref(:)
      REAL(dp), ALLOCATABLE :: ETA_ref(:)
      REAL(dp), ALLOCATABLE :: T_ref(:)
      REAL(dp), ALLOCATABLE :: W_ref(:)

      REAL(dp), ALLOCATABLE :: moments_facets(:,:)

      REAL(dp), ALLOCATABLE :: XYZTW(:,:)
      REAL(dp), ALLOCATABLE :: WV_X(:)
      REAL(dp), ALLOCATABLE :: facet_moments(:)

      !*******************************************************************************

      ! Dimensioni
      num_indici = SIZE(chebyshev_indices, 1)
      n_facce = SIZE(facets, 1)
      n_vertici = SIZE(vertices_4d, 1) / 2

      ! Costruzione della regola di quadratura sul prisma di riferimento

      CALL PrismQuadratureRule(ade, method, XI_ref, ETA_ref, T_ref, W_ref)

      ALLOCATE(moments_facets(num_indici, n_facce))
      moments_facets = 0.0_dp

      ! Parallelizzazione sulle facce
      !
      ! Ogni faccia viene elaborata indipendentemente dalle altre.
      ! Le variabili temporanee relative alla singola faccia sono private
      ! a ciascun thread.

      !$OMP PARALLEL DEFAULT(SHARED) &
      !$OMP PRIVATE(f, n1, n2, n3, V_prism, XYZTW, WV_X, facet_moments)

      !$OMP DO SCHEDULE(STATIC)

      DO f = 1, n_facce

         ALLOCATE(XYZTW(SIZE(W_ref), 4))
         ALLOCATE(WV_X(SIZE(W_ref)))

         n1 = facets(f,1)
         n2 = facets(f,2)
         n3 = facets(f,3)

         V_prism(1,:) = vertices_4d(n1, :)
         V_prism(2,:) = vertices_4d(n2, :)
         V_prism(3,:) = vertices_4d(n3, :)

         V_prism(4,:) = vertices_4d(n1 + n_vertici, :)
         V_prism(5,:) = vertices_4d(n2 + n_vertici, :)
         V_prism(6,:) = vertices_4d(n3 + n_vertici, :)

         CALL ShiftingPrismQuadrature(V_prism, XI_ref, ETA_ref, T_ref, W_ref, &
                                       XYZTW, WV_X)

         IF (MAXVAL(ABS(WV_X)) > 1.0e-14_dp) THEN

            CALL cubature_tens_chebyshev_facet_4D( &
                  XYZTW, WV_X, chebyshev_indices, dbox, facet_moments)

            moments_facets(:,f) = facet_moments

            DEALLOCATE(facet_moments)

         END IF

         DEALLOCATE(XYZTW, WV_X)

      END DO

      !$OMP END DO
      !$OMP END PARALLEL

      ALLOCATE(moments(num_indici))
      moments = SUM(moments_facets, DIM=2)

      DEALLOCATE(moments_facets)

      DEALLOCATE(XI_ref)
      DEALLOCATE(ETA_ref)
      DEALLOCATE(T_ref)
      DEALLOCATE(W_ref)

   END SUBROUTINE chebyshev_moments_polyhedron_4D

   SUBROUTINE cubature_tens_chebyshev_facet_4D(XYZTW, WV_X, chebyshev_indices, dbox, chebyshev_moms)

      !*******************************************************************************
      ! Calcola i momenti di Chebyshev su una faccia 4D
      !*******************************************************************************

      IMPLICIT NONE

      !*******************************************************************************
      ! Argomenti
      !*******************************************************************************
      REAL(dp), INTENT(IN)                  :: XYZTW(:, :)
      REAL(dp), INTENT(IN)                  :: WV_X(:)
      INTEGER, INTENT(IN)                   :: chebyshev_indices(:, :)
      REAL(dp), INTENT(IN)                  :: dbox(2,4)
      REAL(dp), ALLOCATABLE, INTENT(OUT)    :: chebyshev_moms(:)
      !*******************************************************************************
      ! Variabili locali
      !*******************************************************************************
      INTEGER                              :: n, m, max_i, max_j, max_k, max_l, c, iv
      REAL(dp)                             :: B1
      REAL(dp), ALLOCATABLE                :: XN(:), YN(:), ZN(:), TN(:)
      REAL(dp), ALLOCATABLE                :: TX(:, :), TY(:, :), TZ(:, :), TT(:, :)
      REAL(dp), ALLOCATABLE                :: IntX(:, :)
      INTEGER, ALLOCATABLE                 :: idx_i(:), idx_j(:), idx_k(:), idx_l(:)
      !*******************************************************************************

      n = SIZE(XYZTW, 1)
      m = SIZE(chebyshev_indices, 1)

      ALLOCATE(XN(n), YN(n), ZN(n), TN(n))

      XN = (2.0_dp * XYZTW(:, 1) - dbox(1,1) - dbox(2,1)) / &
            (dbox(2,1) - dbox(1,1))

      YN = (2.0_dp * XYZTW(:, 2) - dbox(1,2) - dbox(2,2)) / &
            (dbox(2,2) - dbox(1,2))

      ZN = (2.0_dp * XYZTW(:, 3) - dbox(1,3) - dbox(2,3)) / &
            (dbox(2,3) - dbox(1,3))

      TN = (2.0_dp * XYZTW(:, 4) - dbox(1,4) - dbox(2,4)) / &
            (dbox(2,4) - dbox(1,4))

      B1 = (dbox(2,1) - dbox(1,1)) / 2.0_dp

      max_i = MAXVAL(chebyshev_indices(:, 1))
      max_j = MAXVAL(chebyshev_indices(:, 2))
      max_k = MAXVAL(chebyshev_indices(:, 3))
      max_l = MAXVAL(chebyshev_indices(:, 4))

      ALLOCATE(TX(n, max_i + 2))
      ALLOCATE(TY(n, max_j + 1))
      ALLOCATE(TZ(n, max_k + 1))
      ALLOCATE(TT(n, max_l + 1))

      CALL chebpolys(max_i + 1, XN, TX)
      CALL chebpolys(max_j,     YN, TY)
      CALL chebpolys(max_k,     ZN, TZ)
      CALL chebpolys(max_l,     TN, TT)

      ALLOCATE(IntX(n, max_i + 1))

      IntX = 0.0_dp

      IntX(:, 1) = XN

      IF (max_i >= 1) THEN
         IntX(:, 2) = (XN**2) / 2.0_dp
      END IF

      IF (max_i >= 2) THEN
         DO iv = 2, max_i
            IntX(:, iv + 1) = TX(:, iv+2) / (2.0_dp * (iv+1)) - TX(:, iv)   / (2.0_dp * (iv-1))
         END DO
      END IF

      ALLOCATE(idx_i(m), idx_j(m), idx_k(m), idx_l(m))

      idx_i = chebyshev_indices(:, 1) + 1
      idx_j = chebyshev_indices(:, 2) + 1
      idx_k = chebyshev_indices(:, 3) + 1
      idx_l = chebyshev_indices(:, 4) + 1

      IF (ALLOCATED(chebyshev_moms)) DEALLOCATE(chebyshev_moms)

      ALLOCATE(chebyshev_moms(m))

      DO c = 1, m
         chebyshev_moms(c) = B1 * SUM(WV_X(:) * IntX(:, idx_i(c)) * TY(:, idx_j(c)) * TZ(:, idx_k(c)) * TT(:, idx_l(c)) )
      END DO

      DEALLOCATE(XN, YN, ZN, TN, TX, TY, TZ, TT, IntX, idx_i, idx_j, idx_k, idx_l)

   END SUBROUTINE cubature_tens_chebyshev_facet_4D

   SUBROUTINE chebpolys(deg, x, T)

      !*******************************************************************************
      ! Calcola i polinomi di Chebyshev di prima specie T_j(x) fino al grado "deg" 
      ! valutati su un vettore di punti "x", utilizzando la formula di ricorrenza 
      ! a tre termini e memorizzando i risultati in una matrice T.
      !*******************************************************************************

      IMPLICIT NONE

      !*******************************************************************************
      ! Argomenti
      !*******************************************************************************
      INTEGER, INTENT(IN)                          :: deg
      REAL(dp), INTENT(IN)                         :: x(:)
      REAL(dp), INTENT(OUT)                        :: T(:, :)
      !*******************************************************************************
      ! Variabili locali
      !*******************************************************************************
      INTEGER                                      :: n, j
      REAL(dp), ALLOCATABLE                        :: t0(:), t1(:), t2(:)
      !*******************************************************************************

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

   SUBROUTINE scale_rule(XYZTW_tens_ref, dbox, XYZTW_tens)

      !*******************************************************************************
      ! Scala le coordinate dei punti di quadratura di una regola tensoriale 4D
      ! dal dominio di riferimento [-1,1]^4 al bounding box geometrico del dominio.
      ! I pesi (colonna 5) restano inalterati.
      !*******************************************************************************

      IMPLICIT NONE

      !*******************************************************************************
      ! Argomenti
      !*******************************************************************************
      REAL(dp), INTENT(IN)                            :: XYZTW_tens_ref(:,:)
      REAL(dp), INTENT(IN)                            :: dbox(2,4)             
      REAL(dp), ALLOCATABLE, INTENT(OUT)              :: XYZTW_tens(:,:)
      !*******************************************************************************
      ! Variabili locali
      !*******************************************************************************
      INTEGER  :: k
      REAL(dp) :: a, b
      !*******************************************************************************

      ALLOCATE(XYZTW_tens(SIZE(XYZTW_tens_ref, 1), SIZE(XYZTW_tens_ref, 2)))

      DO k = 1, 4
         a = dbox(1,k)
         b = dbox(2,k)
         XYZTW_tens(:,k) = (a + b) / 2.0_dp + ((b - a) / 2.0_dp) * XYZTW_tens_ref(:,k)
      END DO

      XYZTW_tens(:,5) = XYZTW_tens_ref(:,5)

   END SUBROUTINE scale_rule

END MODULE CubatureTensor