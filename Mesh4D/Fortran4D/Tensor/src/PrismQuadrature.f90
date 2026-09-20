MODULE PrismQuadrature
  USE TypesDef
  USE triangleQuadratureGJ
  USE TriangleQuadratureDunavant
  
  IMPLICIT NONE

CONTAINS

  SUBROUTINE PrismQuadratureRule(ade, method, XI_ref, ETA_ref, T_ref, W_ref)
      INTEGER, INTENT(IN)                   :: ade
      CHARACTER(LEN=*), INTENT(IN)          :: method
      REAL(dp), ALLOCATABLE, INTENT(OUT)    :: XI_ref(:), ETA_ref(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)    :: T_ref(:), W_ref(:)

      REAL(dp), ALLOCATABLE                 :: tri_nodes(:,:), tri_weights(:)
      REAL(dp), ALLOCATABLE                 :: tau_nodes(:), tau_weights(:)
      INTEGER                               :: nGP_tri, n1D, n_tri, n_tau, i, k, idx

      SELECT CASE (TRIM(method))
      CASE ('DCC', 'DGL')
         CALL TriangleQuadratureDunavantPoints(ade + 1, tri_nodes, tri_weights)
         tri_weights = 0.5_dp * tri_weights
      CASE ('GJCC', 'GJL')
         nGP_tri = CEILING(REAL(ade + 2, dp) / 2.0_dp)
         CALL TriangleGJPoints(nGP_tri, tri_nodes, tri_weights)
      CASE DEFAULT
         ERROR STOP 'Metodo non valido. Usare DCC, DGL, GJCC oppure GJL.'
      END SELECT

      SELECT CASE (TRIM(method))
      CASE ('DCC', 'GJCC')
         CALL ClenshawCurtisInterval(ade + 4, tau_nodes, tau_weights)
      CASE ('DGL', 'GJL')
         n1D = CEILING(REAL(ade + 4, dp) / 2.0_dp)
         ALLOCATE(tau_nodes(n1D), tau_weights(n1D))
         CALL gauleg(0.0_dp, 1.0_dp, tau_nodes, tau_weights, n1D)
      END SELECT

      n_tri = SIZE(tri_weights)
      n_tau = SIZE(tau_weights)
      ALLOCATE(XI_ref(n_tri * n_tau), ETA_ref(n_tri * n_tau), T_ref(n_tri * n_tau), W_ref(n_tri * n_tau))

      idx = 1
      DO k = 1, n_tau
         DO i = 1, n_tri
            XI_ref(idx) = tri_nodes(i,1)
            ETA_ref(idx) = tri_nodes(i,2)
            T_ref(idx) = tau_nodes(k)
            W_ref(idx) = tau_weights(k) * tri_weights(i)
            idx = idx + 1
         END DO
      END DO
   END SUBROUTINE PrismQuadratureRule

  SUBROUTINE reference_prism_quadrature(ade, XI_ref, ETA_ref, T_ref, W_ref)

   !*******************************************************************************
   ! Genera i nodi e i pesi di quadratura di riferimento per il prisma 4D
   ! (prodotto cartesiano tra il triangolo 2D di riferimento e l'intervallo
   ! temporale [0,1]).
   !
   ! Grado di precisione differenziato:
   !   - nGP_tri per la base triangolare (xi, eta), via Gauss-Jacobi conico
   !   - n1D     per l'asse temporale tau, via Gauss-Legendre
   !*******************************************************************************

      IMPLICIT NONE

      ! Input variables
      INTEGER, INTENT(IN)                   :: ade

      ! Output variables
      REAL(dp), ALLOCATABLE, INTENT(OUT)    :: XI_ref(:), ETA_ref(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)    :: T_ref(:), W_ref(:)

      ! Local variables
      INTEGER                               :: nGP_tri, n1D
      INTEGER                               :: n_tri, n_tau_loc, n_tot
      INTEGER                               :: i, k, idx
      REAL(dp), ALLOCATABLE                 :: P_std(:,:), W_std(:)
      REAL(dp), ALLOCATABLE                 :: r(:), s(:)
      REAL(dp), ALLOCATABLE                 :: tau_nodes_loc(:), tau_weights_loc(:)

      ! Grado di precisione per la base triangolare
      nGP_tri = CEILING(REAL(ade + 2, dp) / 2.0_dp) + 1

      n_tri = nGP_tri * nGP_tri
      ALLOCATE(P_std(2, n_tri))
      ALLOCATE(W_std(n_tri))

      CALL TriangleQuadraturePoints(P_std, W_std, n_tri, nGP_tri)

      ALLOCATE(r(n_tri), s(n_tri))
      r = P_std(1, :)
      s = P_std(2, :)

      ! Grado di precisione per la dimensione temporale tau su [0,1]
      n1D = CEILING(REAL(ade + 4, dp) / 2.0_dp)
      n_tau_loc = n1D

      ALLOCATE(tau_nodes_loc(n_tau_loc), tau_weights_loc(n_tau_loc))
      CALL gauleg(0.0_dp, 1.0_dp, tau_nodes_loc, tau_weights_loc, n_tau_loc)

      ! Costruzione della griglia come prodotto tensoriale triangolo x tempo
      n_tot = n_tri * n_tau_loc

      IF (ALLOCATED(XI_ref))  DEALLOCATE(XI_ref)
      IF (ALLOCATED(ETA_ref)) DEALLOCATE(ETA_ref)
      IF (ALLOCATED(T_ref))   DEALLOCATE(T_ref)
      IF (ALLOCATED(W_ref))   DEALLOCATE(W_ref)

      ALLOCATE(XI_ref(n_tot))
      ALLOCATE(ETA_ref(n_tot))
      ALLOCATE(T_ref(n_tot))
      ALLOCATE(W_ref(n_tot))

      idx = 1
      DO i = 1, n_tri
         DO k = 1, n_tau_loc
            XI_ref(idx)  = r(i)
            ETA_ref(idx) = s(i)
            T_ref(idx)   = tau_nodes_loc(k)
            W_ref(idx)   = W_std(i) * tau_weights_loc(k)
            idx = idx + 1
         END DO
      END DO

      DEALLOCATE(P_std, W_std, r, s, tau_nodes_loc, tau_weights_loc)

   END SUBROUTINE reference_prism_quadrature

   SUBROUTINE TriangleGJPoints(nGP, nodes, weights)
      INTEGER, INTENT(IN)                    :: nGP
      REAL(dp), ALLOCATABLE, INTENT(OUT)     :: nodes(:,:), weights(:)

      REAL(dp), ALLOCATABLE                  :: pstd(:,:), wstd(:)
      INTEGER                                :: n_points

      n_points = nGP * nGP
      ALLOCATE(pstd(2, n_points), wstd(n_points))
      CALL TriangleQuadraturePoints(pstd, wstd, n_points, nGP)
      ALLOCATE(nodes(n_points, 2), weights(n_points))
      nodes(:,1) = pstd(1,:)
      nodes(:,2) = pstd(2,:)
      weights = wstd
   END SUBROUTINE TriangleGJPoints

   SUBROUTINE ClenshawCurtisInterval(n, nodes, weights)
      INTEGER, INTENT(IN)                    :: n
      REAL(dp), ALLOCATABLE, INTENT(OUT)     :: nodes(:), weights(:)

      INTEGER                                :: N_cc, i, j
      REAL(dp)                               :: theta_i, sum_w, fact_j, g_cc

      N_cc = n - 1
      ALLOCATE(nodes(n), weights(n))

      IF (n == 1) THEN
         nodes(1) = 0.5_dp
         weights(1) = 1.0_dp
         RETURN
      END IF

      DO i = 1, n
         theta_i = PI * REAL(i - 1, dp) / REAL(N_cc, dp)
         sum_w = 0.0_dp
         DO j = 1, N_cc / 2
            fact_j = 2.0_dp
            IF (MOD(N_cc, 2) == 0 .AND. j == N_cc / 2) fact_j = 1.0_dp
            sum_w = sum_w + fact_j * COS(2.0_dp * REAL(j, dp) * theta_i) / (1.0_dp - 4.0_dp * REAL(j, dp)**2)
         END DO
         g_cc = 1.0_dp
         IF (i == 1 .OR. i == n) g_cc = 0.5_dp
         nodes(n - i + 1) = 0.5_dp * (COS(theta_i) + 1.0_dp)
         weights(n - i + 1) = (1.0_dp / REAL(N_cc, dp)) * g_cc * (1.0_dp + sum_w)
      END DO
   END SUBROUTINE ClenshawCurtisInterval

   SUBROUTINE ShiftingPrismQuadrature(V, XI_ref, ETA_ref, T_ref, W_ref, XYZTW, WV_X)
      REAL(dp), INTENT(IN)     :: V(6,4)
      REAL(dp), INTENT(IN)     :: XI_ref(:), ETA_ref(:), T_ref(:), W_ref(:)
      REAL(dp), INTENT(OUT)    :: XYZTW(:,:), WV_X(:)

      INTEGER                  :: j
      REAL(dp)                 :: A0(4), B0(4), C0(4), A1(4), B1(4), C1(4)
      REAL(dp)                 :: lato1_base(4), lato2_base(4), lato1_top(4), lato2_top(4)
      REAL(dp), ALLOCATABLE    :: base(:,:), top(:,:)
      REAL(dp), ALLOCATABLE    :: dy_dxi(:), dz_dxi(:), dy_deta(:), dz_deta(:), Nx(:)

      A0 = V(1,:); B0 = V(2,:); C0 = V(3,:)
      A1 = V(4,:); B1 = V(5,:); C1 = V(6,:)

      lato1_base = B0 - A0
      lato2_base = C0 - A0
      lato1_top = B1 - A1
      lato2_top = C1 - A1

      ALLOCATE(base(SIZE(W_ref),4), top(SIZE(W_ref),4))
      DO j = 1, 4
         base(:,j) = A0(j) + XI_ref * lato1_base(j) + ETA_ref * lato2_base(j)
         top(:,j) = A1(j) + XI_ref * lato1_top(j) + ETA_ref * lato2_top(j)
         XYZTW(:,j) = (1.0_dp - T_ref) * base(:,j) + T_ref * top(:,j)
      END DO

      ALLOCATE(dy_dxi(SIZE(W_ref)), dz_dxi(SIZE(W_ref)), dy_deta(SIZE(W_ref)), dz_deta(SIZE(W_ref)), Nx(SIZE(W_ref)))
      dy_dxi = (1.0_dp - T_ref) * lato1_base(2) + T_ref * lato1_top(2)
      dz_dxi = (1.0_dp - T_ref) * lato1_base(3) + T_ref * lato1_top(3)
      dy_deta = (1.0_dp - T_ref) * lato2_base(2) + T_ref * lato2_top(2)
      dz_deta = (1.0_dp - T_ref) * lato2_base(3) + T_ref * lato2_top(3)

      Nx = dy_dxi * dz_deta - dz_dxi * dy_deta
      WV_X = W_ref * Nx
   END SUBROUTINE ShiftingPrismQuadrature

   SUBROUTINE PrismQuad4D(V, XI_ref, ETA_ref, T_ref, W_ref, baricentro, XYZTW, WV_X)

      !*******************************************************************************
      ! Costruisce i nodi fisici e i pesi orientati per integrare su una
      ! iperfaccia laterale spazio-temporale 3D immersa in R^4.
      !
      ! Il prisma di riferimento e':
      !
      !   K_ref = {(xi,eta,t) : xi >= 0, eta >= 0, xi + eta <= 1, 0 <= t <= 1}.
      !
      ! La mappa verso R^4 e' X(xi,eta,t) = (1-t) X0(xi,eta) + t X1(xi,eta), con:
      !
      !   X0(xi,eta) = A0 + xi (B0-A0) + eta (C0-A0),
      !   X1(xi,eta) = A1 + xi (B1-A1) + eta (C1-A1).
      !
      ! Il vettore superficie orientato locale e' N(xi,eta,t) = n(xi,eta,t) dS,
      ! ortogonale ai tre vettori tangenti dX/dxi, dX/deta, dX/dt. Per il teorema
      ! della divergenza applicato al campo F = (Phi_x, 0, 0, 0), serve integrare
      ! Phi_x * n_x dS. Quindi la subroutine restituisce direttamente:
      !
      !   WV_X(q) = W_ref(q) * N_x(q)
      !
      ! INPUT:
      !   V           matrice 6 x 4 dei vertici del prisma spazio-temporale:
      !                 V(1,:) = A0, vertice 1 della faccia a tau = 0
      !                 V(2,:) = B0, vertice 2 della faccia a tau = 0
      !                 V(3,:) = C0, vertice 3 della faccia a tau = 0
      !                 V(4,:) = A1, vertice 1 della faccia a tau = 1
      !                 V(5,:) = B1, vertice 2 della faccia a tau = 1
      !                 V(6,:) = C1, vertice 3 della faccia a tau = 1
      !               Ogni riga e' [x, y, z, tau].
      !
      !   XI_ref      vettore Nq delle coordinate xi dei nodi di quadratura
      !   ETA_ref     vettore Nq delle coordinate eta dei nodi di quadratura
      !   T_ref       vettore Nq delle coordinate temporali in [0,1]
      !   W_ref       vettore Nq dei pesi sul prisma di riferimento
      !   baricentro  vettore di 4 componenti interno al dominio 4D, usato per
      !               orientare N verso l'esterno
      !
      ! OUTPUT:
      !   XYZTW       matrice Nq x 4 dei nodi fisici 4D [x, y, z, tau]
      !   WV_X        vettore Nq dei pesi orientati: WV_X = W_ref .* N_x
      !*******************************************************************************

      IMPLICIT NONE

      ! Input variables
      REAL(dp), INTENT(IN)                  :: V(:, :)
      REAL(dp), INTENT(IN)                  :: XI_ref(:)
      REAL(dp), INTENT(IN)                  :: ETA_ref(:)
      REAL(dp), INTENT(IN)                  :: T_ref(:)
      REAL(dp), INTENT(IN)                  :: W_ref(:)
      REAL(dp), INTENT(IN)                  :: baricentro(4)

      ! Output variables
      REAL(dp), ALLOCATABLE, INTENT(OUT)    :: XYZTW(:, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT)    :: WV_X(:)

      ! Local variables
      INTEGER                               :: num_pts, i
      REAL(dp), ALLOCATABLE                 :: lato1_base(:), lato2_base(:)
      REAL(dp), ALLOCATABLE                 :: lato1_top(:), lato2_top(:)
      REAL(dp), ALLOCATABLE                 :: nodi_base_4D(:, :), nodi_top_4D(:, :)
      REAL(dp), ALLOCATABLE                 :: c(:, :), r(:, :), z(:, :)
      REAL(dp), ALLOCATABLE                 :: M1(:), M2(:), M3(:), M4(:)
      REAL(dp), ALLOCATABLE                 :: Nvec(:, :)
      REAL(dp), ALLOCATABLE                 :: direzione_uscente(:, :)
      REAL(dp), ALLOCATABLE                 :: dot_flip(:)
      LOGICAL, ALLOCATABLE                  :: flip(:)

      num_pts = SIZE(W_ref)

      IF (ALLOCATED(XYZTW)) DEALLOCATE(XYZTW)
      IF (ALLOCATED(WV_X))  DEALLOCATE(WV_X)

      ALLOCATE(XYZTW(num_pts, 4))
      ALLOCATE(WV_X(num_pts))
      ALLOCATE(lato1_base(4), lato2_base(4), lato1_top(4), lato2_top(4))
      ALLOCATE(nodi_base_4D(num_pts, 4), nodi_top_4D(num_pts, 4))
      ALLOCATE(c(num_pts, 4), r(num_pts, 4), z(num_pts, 4))
      ALLOCATE(M1(num_pts), M2(num_pts), M3(num_pts), M4(num_pts))
      ALLOCATE(Nvec(num_pts, 4))
      ALLOCATE(direzione_uscente(num_pts, 4))
      ALLOCATE(dot_flip(num_pts))
      ALLOCATE(flip(num_pts))

      ! Vettori lato della faccia iniziale (tau = 0) e finale (tau = 1)
      lato1_base = V(2, :) - V(1, :)
      lato2_base = V(3, :) - V(1, :)
      lato1_top  = V(5, :) - V(4, :)
      lato2_top  = V(6, :) - V(4, :)

      ! Nodi sulla faccia iniziale X0(xi,eta) e finale X1(xi,eta)
      DO i = 1, 4
         nodi_base_4D(:, i) = V(1, i) + XI_ref * lato1_base(i) + ETA_ref * lato2_base(i)
         nodi_top_4D(:, i)  = V(4, i) + XI_ref * lato1_top(i)  + ETA_ref * lato2_top(i)
      END DO

      ! Nodi fisici sulla faccia laterale spazio-temporale:
      ! X(xi,eta,t) = (1-t) X0(xi,eta) + t X1(xi,eta)
      DO i = 1, 4
         XYZTW(:, i) = (1.0_dp - T_ref) * nodi_base_4D(:, i) + T_ref * nodi_top_4D(:, i)
      END DO

      ! Vettori tangenti rispetto a xi ed eta
      DO i = 1, 4
         c(:, i) = (1.0_dp - T_ref) * lato1_base(i) + T_ref * lato1_top(i)
         r(:, i) = (1.0_dp - T_ref) * lato2_base(i) + T_ref * lato2_top(i)
      END DO

      ! Vettore tangente rispetto al tempo (dipende da xi ed eta:
      ! in generale la normale non e' costante se la mesh si deforma)
      z = nodi_top_4D - nodi_base_4D

      ! Minori 3 x 3 ottenuti eliminando, rispettivamente, la colonna
      ! x, y, z, tau dalla matrice jacobiana 3 x 4:
      !
      !   J = [c; r; z]
      !
      M1 = c(:,2) * (r(:,3)*z(:,4) - r(:,4)*z(:,3)) &
         - c(:,3) * (r(:,2)*z(:,4) - r(:,4)*z(:,2)) &
         + c(:,4) * (r(:,2)*z(:,3) - r(:,3)*z(:,2))

      M2 = c(:,1) * (r(:,3)*z(:,4) - r(:,4)*z(:,3)) &
         - c(:,3) * (r(:,1)*z(:,4) - r(:,4)*z(:,1)) &
         + c(:,4) * (r(:,1)*z(:,3) - r(:,3)*z(:,1))

      M3 = c(:,1) * (r(:,2)*z(:,4) - r(:,4)*z(:,2)) &
         - c(:,2) * (r(:,1)*z(:,4) - r(:,4)*z(:,1)) &
         + c(:,4) * (r(:,1)*z(:,2) - r(:,2)*z(:,1))

      M4 = c(:,1) * (r(:,2)*z(:,3) - r(:,3)*z(:,2)) &
         - c(:,2) * (r(:,1)*z(:,3) - r(:,3)*z(:,1)) &
         + c(:,3) * (r(:,1)*z(:,2) - r(:,2)*z(:,1))

      ! Vettore superficie 4D con i segni alternati:
      !
      !   N = [M1, -M2, M3, -M4] = n dS
      !
      ! Ogni riga di Nvec e' ortogonale a dX/dxi, dX/deta e dX/dt
      Nvec(:, 1) =  M1
      Nvec(:, 2) = -M2
      Nvec(:, 3) =  M3
      Nvec(:, 4) = -M4

      ! Orientamento uscente, se N punta verso il baricentro, viene ribaltato
      DO i = 1, 4
         direzione_uscente(:, i) = XYZTW(:, i) - baricentro(i)
      END DO

      dot_flip = Nvec(:,1) * direzione_uscente(:,1) + Nvec(:,2) * direzione_uscente(:,2) &
               + Nvec(:,3) * direzione_uscente(:,3) + Nvec(:,4) * direzione_uscente(:,4)

      flip = (dot_flip < 0.0_dp)

      DO i = 1, num_pts
         IF (flip(i)) THEN
            Nvec(i, :) = -Nvec(i, :)
         END IF
      END DO

      ! Peso orientato richiesto dalla formula di divergenza con
      ! F = (Phi_x, 0, 0, 0): F dot n dS = Phi_x * N_x
      WV_X = W_ref * Nvec(:, 1)

      DEALLOCATE(lato1_base, lato2_base, lato1_top, lato2_top)
      DEALLOCATE(nodi_base_4D, nodi_top_4D, c, r, z)
      DEALLOCATE(M1, M2, M3, M4, Nvec, direzione_uscente, dot_flip, flip)

   END SUBROUTINE PrismQuad4D

END MODULE PrismQuadrature
