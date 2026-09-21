MODULE CubatureFunctions

   USE TypesDef

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE triangle_quadrature(method, ade, nodes, weights)

      USE TriangleQuadratureGJ
      USE TriangleQuadratureDunavant

      IMPLICIT NONE

      !********************************************************************************
      ! Argomenti
      !********************************************************************************
      CHARACTER(LEN=*), INTENT(IN) :: method
      INTEGER, INTENT(IN) :: ade
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: nodes(:,:)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: weights(:)
      !********************************************************************************
      ! Variabili locali
      !********************************************************************************
      INTEGER :: nGP
      INTEGER :: n_points
      INTEGER :: rule
      !********************************************************************************

      ! METODO GJ: Quadratura di Gauss-Jacobi
      IF (TRIM(method) == 'GJ') THEN

         ! Numero di punti di Gauss-Jacobi 1D necessari per integrare
         ! esattamente sul triangolo un polinomio di grado (ade+1).
         nGP = CEILING((ade + 2.0_dp) / 2.0_dp)

         n_points = nGP * nGP

         ALLOCATE(nodes(n_points,2))
         ALLOCATE(weights(n_points))

         CALL TriangleQuadratureGJPoints(nodes, weights, n_points, nGP)

         ! I pesi della quadratura GJ standard sono riferiti a un
         ! triangolo di area 1/2.
         weights = 2.0_dp * weights

      ! METODO D: Quadratura di Dunavant
      ELSE IF (TRIM(method) == 'D') THEN

         rule = ade + 1

         ! Calcolo dei nodi e dei pesi sul triangolo di riferimento
         CALL TriangleQuadratureDunavantPoints(rule, nodes, weights)

      ! Caso in cui il metodo inserito non sia valido
      ELSE

         WRITE(*,'(A)') 'ERROR: Metodo di quadratura non valido. Usare ''GJ'' oppure ''D''.'
         RETURN

      END IF

   END SUBROUTINE triangle_quadrature

   
   SUBROUTINE chebyshev_moments_polyhedron(vertices, facets, ade, chebyshev_indices, dbox, method, moments)

      !*******************************************************************************
      !
      ! Calcola i momenti dei polinomi di Chebyshev sul poliedro tridimensionale
      ! delimitato dalla mesh triangolare definita da vertices e facets.
      ! I momenti volumetrici sono trasformati, mediante il teorema della divergenza,
      ! in integrali di superficie sulle facce triangolari del poliedro. Gli integrali
      ! sulle singole facce vengono quindi valutati numericamente mediante una formula
      ! di quadratura sul triangolo.
      ! La quadratura sulle facce può essere effettuata mediante uno dei seguenti
      ! metodi:
      !
      !   'GJ' : quadratura di Gauss-Jacobi sul triangolo di riferimento;
      !   'D'  : quadratura di Dunavant sul triangolo.
      !
      ! Per ciascuna faccia viene calcolata la normale esterna unitaria e il relativo
      ! contributo ai momenti. I contributi ottenuti da tutte le facce vengono infine
      ! sommati per ottenere i momenti del poliedro.
      !
      !*******************************************************************************

      USE TriangleQuadratureGJ
      USE TriangleQuadratureDunavant

      IMPLICIT NONE

      !*******************************************************************************
      ! Argomenti
      !*******************************************************************************
      REAL(dp), INTENT(IN)               :: vertices(:,:)
      INTEGER, INTENT(IN)                :: facets(:,:)
      INTEGER, INTENT(IN)                :: ade
      INTEGER, INTENT(IN)                :: chebyshev_indices(:,:)
      REAL(dp), INTENT(IN)               :: dbox(6)
      CHARACTER(LEN=*), INTENT(IN)       :: method
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: moments(:)
      !*******************************************************************************
      ! Variabili locali
      !*******************************************************************************
      INTEGER :: num_indici
      INTEGER :: n_facce
      INTEGER :: n_quad
      INTEGER :: k
      INTEGER :: q

      REAL(dp) :: A(3)
      REAL(dp) :: B(3)

      REAL(dp) :: V_face(3,3)

      REAL(dp) :: norm_ext(3)

      REAL(dp), ALLOCATABLE :: XYZW(:,:)
      REAL(dp), ALLOCATABLE :: WV_CUB(:)

      REAL(dp), ALLOCATABLE :: moms_facet_raw(:)

      REAL(dp), ALLOCATABLE :: nodes(:,:)
      REAL(dp), ALLOCATABLE :: weights(:)

      REAL(dp) :: cp(3)
      REAL(dp) :: area2
      !*******************************************************************************

      ! Dimensioni
      num_indici = SIZE(chebyshev_indices, 1)
      n_facce = SIZE(facets, 1)

      CALL triangle_quadrature(method, ade, nodes, weights)
      n_quad = SIZE(weights)

      ALLOCATE(XYZW(n_quad, 3))
      ALLOCATE(WV_CUB(n_quad))
      ALLOCATE(moments(num_indici))
      moments = 0.0_dp

      DO k = 1, n_facce

         V_face(1,:) = vertices(facets(k,1),:)
         V_face(2,:) = vertices(facets(k,2),:)
         V_face(3,:) = vertices(facets(k,3),:)

         ! Calcolo della normale alla faccia
         A = V_face(2,:) - V_face(1,:)
         B = V_face(3,:) - V_face(1,:)

         cp(1) = A(2)*B(3) - A(3)*B(2)
         cp(2) = A(3)*B(1) - A(1)*B(3)
         cp(3) = A(1)*B(2) - A(2)*B(1)

         area2 = SQRT(SUM(cp**2))

         IF (area2 <= 1.0d-14) THEN
            CYCLE
         END IF

         ! Trasformazione affine dei punti dal triangolo di riferimento
         ! al triangolo fisico x = V1 + xi * (V2-V1) + eta * (V3-V1).
         DO q = 1, n_quad
            XYZW(q,:) = V_face(1,:) + nodes(q,1) * A + nodes(q,2) * B
         END DO

         WV_CUB = 0.5_dp * area2 * weights
         norm_ext = cp / area2

         ! Calcolo dei momenti sulla faccia
         CALL cubature_tens_chebyshev_facet_V(XYZW, WV_CUB, chebyshev_indices, dbox, moms_facet_raw)

         ! Contributo della faccia
         moments = moments + norm_ext(1) * moms_facet_raw

      END DO

   END SUBROUTINE chebyshev_moments_polyhedron


   SUBROUTINE cubature_tens_chebyshev_facet_V(nodes, weights, chebyshev_indices, dbox, chebyshev_moms)

      !*******************************************************************************
      ! Calcola i momenti di Chebyshev su una faccia triangolare 
      !*******************************************************************************

      IMPLICIT NONE

      !*******************************************************************************
      ! Argomenti
      !*******************************************************************************
      REAL(dp), INTENT(IN)    :: nodes(:, :)        
      REAL(dp), INTENT(IN)    :: weights(:)         
      INTEGER, INTENT(IN)    :: chebyshev_indices(:, :) 
      REAL(dp), INTENT(IN)    :: dbox(6)
      REAL(dp), ALLOCATABLE, INTENT(OUT)   :: chebyshev_moms(:)  
      !*******************************************************************************
      ! Variabili locali
      !*******************************************************************************
      INTEGER                :: n, m, max_i, max_j, max_k, c, iv
      REAL(dp)                :: B1
      REAL(dp), ALLOCATABLE   :: XN(:), YN(:), ZN(:)
      REAL(dp), ALLOCATABLE   :: TX(:, :), TY(:, :), TZ(:, :)
      REAL(dp), ALLOCATABLE   :: IntX(:, :)
      INTEGER, ALLOCATABLE   :: idx_i(:), idx_j(:), idx_k(:)
      REAL(dp), ALLOCATABLE   :: w(:)
      !*******************************************************************************

      n = SIZE(nodes, 1)
      m = SIZE(chebyshev_indices, 1)

      ALLOCATE(XN(n), YN(n), ZN(n))
      ALLOCATE(w(n))

      XN = (nodes(:, 1) - (dbox(1) + dbox(2)) / 2.0_dp) / ((dbox(2) - dbox(1)) / 2.0_dp)
      YN = (nodes(:, 2) - (dbox(3) + dbox(4)) / 2.0_dp) / ((dbox(4) - dbox(3)) / 2.0_dp)
      ZN = (nodes(:, 3) - (dbox(5) + dbox(6)) / 2.0_dp) / ((dbox(6) - dbox(5)) / 2.0_dp)
    
      B1 = (dbox(2) - dbox(1)) / 2.0_dp
      w = weights(:)

      max_i = MAXVAL(chebyshev_indices(:, 1))
      max_j = MAXVAL(chebyshev_indices(:, 2))
      max_k = MAXVAL(chebyshev_indices(:, 3))

      ALLOCATE(TX(n, max_i + 2))
      ALLOCATE(TY(n, max_j + 1))
      ALLOCATE(TZ(n, max_k + 1))

      CALL chebpolys(max_i + 1, XN, TX)
      CALL chebpolys(max_j, YN, TY)
      CALL chebpolys(max_k, ZN, TZ)

      ALLOCATE(IntX(n, max_i + 1))
      IntX = 0.0_dp
    
      IntX(:, 1) = XN 

      IF (max_i >= 1) THEN
         IntX(:, 2) = (XN**2) / 2.0_dp
      END IF

      IF (max_i >= 2) THEN
         DO iv = 2, max_i
            IntX(:, iv + 1) = TX(:, iv+2) / (2.0_dp*(iv+1)) - TX(:, iv) / (2.0_dp*(iv-1))
         END DO
      END IF

      ALLOCATE(idx_i(m), idx_j(m), idx_k(m))
      idx_i = chebyshev_indices(:, 1) + 1
      idx_j = chebyshev_indices(:, 2) + 1
      idx_k = chebyshev_indices(:, 3) + 1

      IF (ALLOCATED(chebyshev_moms)) DEALLOCATE(chebyshev_moms)
      ALLOCATE(chebyshev_moms(m))

      DO c = 1, m
         chebyshev_moms(c) = B1 * SUM(w * IntX(:, idx_i(c)) * TY(:, idx_j(c)) * TZ(:, idx_k(c)))
      END DO

      DEALLOCATE(XN, YN, ZN, w, TX, TY, TZ, IntX, idx_i, idx_j, idx_k)

   END SUBROUTINE cubature_tens_chebyshev_facet_V

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
      INTEGER, INTENT(IN)    :: deg
      REAL(dp), INTENT(IN)    :: x(:)
      REAL(dp), INTENT(OUT)   :: T(:, :)
      !*******************************************************************************
      ! Variabili locali
      !*******************************************************************************
      INTEGER                :: n, j
      REAL(dp), ALLOCATABLE   :: t0(:), t1(:), t2(:)
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

   SUBROUTINE scale_rule(XYZW_tens_ref, dbox, XYZW_tens)

      !*******************************************************************************
      ! Scala le coordinate dei punti di quadratura di una regola tensoriale 
      ! dal dominio di riferimento standard $[-1, 1]^3$ al bounding box geometrico 
      ! effettivo del dominio, lasciando inalterati i pesi di cubatura.
      !*******************************************************************************

      IMPLICIT NONE

      !*******************************************************************************
      ! Argomenti
      !*******************************************************************************
      REAL(dp), INTENT(IN)  :: XYZW_tens_ref(:,:)
      REAL(dp), INTENT(IN)  :: dbox(6)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: XYZW_tens(:,:)
      !*******************************************************************************
      ! Variabili locali
      !*******************************************************************************
      REAL(dp) :: a, b
      !*******************************************************************************
      
      IF (ALLOCATED(XYZW_tens)) DEALLOCATE(XYZW_tens)
      ALLOCATE(XYZW_tens(SIZE(XYZW_tens_ref, 1), SIZE(XYZW_tens_ref, 2)))

      ! X dimension scaling
      a = dbox(1); b = dbox(2)
      XYZW_tens(:,1) = (a + b) / 2.0_dp + ((b - a) / 2.0_dp) * XYZW_tens_ref(:,1)

      ! Y dimension scaling
      a = dbox(3); b = dbox(4)
      XYZW_tens(:,2) = (a + b) / 2.0_dp + ((b - a) / 2.0_dp) * XYZW_tens_ref(:,2)

      ! Z dimension scaling
      a = dbox(5); b = dbox(6)
      XYZW_tens(:,3) = (a + b) / 2.0_dp + ((b - a) / 2.0_dp) * XYZW_tens_ref(:,3)

      ! I pesi rimangono inalterati
      IF (SIZE(XYZW_tens_ref, 2) >= 4) THEN
         XYZW_tens(:,4) = XYZW_tens_ref(:,4)
      END IF

   END SUBROUTINE scale_rule

END MODULE CubatureFunctions
