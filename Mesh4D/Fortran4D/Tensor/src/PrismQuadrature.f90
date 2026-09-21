MODULE PrismQuadrature

   USE TypesDef
   USE TriangleQuadratureGJ
   USE TriangleQuadratureDunavant

   IMPLICIT NONE

CONTAINS

   SUBROUTINE PrismQuadratureRule(ade, method, XI_ref, ETA_ref, T_ref, W_ref)

      IMPLICIT NONE

      !*************************************************************************
      !
      ! Costruisce la regola di quadratura di riferimento sul prisma
      ! triangolare spazio-temporale, ottenuta come prodotto tensoriale
      ! tra una regola sul triangolo di riferimento e una regola
      ! sull'intervallo temporale [0,1].
      !
      ! Il metodo di quadratura puo' essere uno dei seguenti:
      !
      !   'DCC'  : Dunavant--Clenshaw--Curtis
      !   'DGL'  : Dunavant--Gauss--Legendre
      !   'GJCC' : Gauss--Jacobi--Clenshaw--Curtis
      !   'GJL'  : Gauss--Jacobi--Gauss--Legendre
      !
      !*************************************************************************
      ! Argomenti
      !*************************************************************************
      INTEGER, INTENT(IN)                :: ade
      CHARACTER(LEN=*), INTENT(IN)       :: method
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: XI_ref(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: ETA_ref(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: T_ref(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: W_ref(:)
      !*************************************************************************

      ! METODO DCC: Dunavant + Clenshaw--Curtis
      IF (TRIM(method) == 'DCC') THEN

         CALL PrismDCC(ade, XI_ref, ETA_ref, T_ref, W_ref)

      ! METODO DGL: Dunavant + Gauss--Legendre
      ELSE IF (TRIM(method) == 'DGL') THEN

         CALL PrismDGL(ade, XI_ref, ETA_ref, T_ref, W_ref)

      ! METODO GJCC: Gauss--Jacobi + Clenshaw--Curtis
      ELSE IF (TRIM(method) == 'GJCC') THEN

         CALL PrismGJCC(ade, XI_ref, ETA_ref, T_ref, W_ref)

      ! METODO GJL: Gauss--Jacobi + Gauss--Legendre
      ELSE IF (TRIM(method) == 'GJL') THEN

         CALL PrismGJL(ade, XI_ref, ETA_ref, T_ref, W_ref)

      ! Caso in cui il metodo inserito non sia valido
      ELSE

         WRITE(*,'(A)') 'ERROR: Metodo scelto non valido. Usare ''GJL'', ''GJCC'', ''DGL'' o ''DCC''.'
         ERROR STOP 'PrismQuadratureRule: metodo non valido'

      END IF

   END SUBROUTINE PrismQuadratureRule


   SUBROUTINE PrismDCC(ade, XI_ref, ETA_ref, T_ref, W_ref)

      IMPLICIT NONE

      !*************************************************************************
      !
      ! Costruisce la regola di quadratura tensoriale Dunavant--Clenshaw--Curtis
      ! (DCC) sul prisma triangolare di riferimento
      !
      !    P_ref = { (xi,eta,tau) : xi>=0, eta>=0, xi+eta<=1, 0<=tau<=1 },
      !
      ! ottenuta come prodotto tensoriale tra:
      !
      !   - la regola simmetrica di Dunavant sul triangolo di riferimento;
      !
      !   - la regola di Clenshaw--Curtis sull'intervallo temporale [0,1].
      !
      !*************************************************************************
      ! Argomenti
      !*************************************************************************
      INTEGER, INTENT(IN)                :: ade
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: XI_ref(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: ETA_ref(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: T_ref(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: W_ref(:)
      !*************************************************************************
      ! Variabili locali
      !*************************************************************************
      INTEGER               :: rule_tri
      INTEGER               :: n_tri, n_tau
      INTEGER               :: i, k, idx

      REAL(dp), ALLOCATABLE :: tri_nodes(:,:)
      REAL(dp), ALLOCATABLE :: tri_weights(:)
      REAL(dp), ALLOCATABLE :: tau_nodes(:)
      REAL(dp), ALLOCATABLE :: tau_weights(:)
      !*************************************************************************

      !*************************************************************************
      ! Quadratura spaziale
      !*************************************************************************

      rule_tri = ade + 1

      CALL TriangleQuadratureDunavantPoints(rule_tri, tri_nodes, tri_weights)

      ! I pesi di Dunavant sono riferiti a un triangolo di area 1, mentre
      ! il triangolo di riferimento ha area 1/2.
      tri_weights = 0.5_dp * tri_weights

      !*************************************************************************
      ! Quadratura temporale
      !*************************************************************************

      CALL clenshaw_curtis(ade + 4, tau_nodes, tau_weights)

      !*************************************************************************
      ! Prodotto tensoriale
      !*************************************************************************

      n_tri = SIZE(tri_weights)
      n_tau = SIZE(tau_weights)

      ALLOCATE(XI_ref(n_tri * n_tau))
      ALLOCATE(ETA_ref(n_tri * n_tau))
      ALLOCATE(T_ref(n_tri * n_tau))
      ALLOCATE(W_ref(n_tri * n_tau))

      idx = 1
      DO k = 1, n_tau
         DO i = 1, n_tri
            XI_ref(idx)  = tri_nodes(i,1)
            ETA_ref(idx) = tri_nodes(i,2)
            T_ref(idx)   = tau_nodes(k)
            W_ref(idx)   = tau_weights(k) * tri_weights(i)
            idx = idx + 1
         END DO
      END DO

   END SUBROUTINE PrismDCC


   SUBROUTINE PrismDGL(ade, XI_ref, ETA_ref, T_ref, W_ref)

      IMPLICIT NONE

      !*************************************************************************
      !
      ! Costruisce la regola di quadratura tensoriale Dunavant--Gauss--Legendre
      ! (DGL) sul prisma triangolare di riferimento
      !
      !    P_ref = { (xi,eta,tau) : xi>=0, eta>=0, xi+eta<=1, 0<=tau<=1 },
      !
      ! ottenuta come prodotto tensoriale tra:
      !
      !   - la regola simmetrica di Dunavant sul triangolo di riferimento;
      !
      !   - la regola di Gauss--Legendre sull'intervallo temporale [0,1].
      !
      !*************************************************************************
      ! Argomenti
      !*************************************************************************
      INTEGER, INTENT(IN)                :: ade
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: XI_ref(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: ETA_ref(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: T_ref(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: W_ref(:)
      !*************************************************************************
      ! Variabili locali
      !*************************************************************************
      INTEGER               :: rule_tri, n1D
      INTEGER               :: n_tri, n_tau
      INTEGER               :: i, k, idx

      REAL(dp), ALLOCATABLE :: tri_nodes(:,:)
      REAL(dp), ALLOCATABLE :: tri_weights(:)
      REAL(dp), ALLOCATABLE :: tau_nodes(:)
      REAL(dp), ALLOCATABLE :: tau_weights(:)
      !*************************************************************************

      !*************************************************************************
      ! Quadratura spaziale
      !*************************************************************************

      rule_tri = ade + 1

      CALL TriangleQuadratureDunavantPoints(rule_tri, tri_nodes, tri_weights)

      ! I pesi di Dunavant sono riferiti a un triangolo di area 1, mentre
      ! il triangolo di riferimento ha area 1/2.
      tri_weights = 0.5_dp * tri_weights

      !*************************************************************************
      ! Quadratura temporale
      !*************************************************************************

      n1D = CEILING((ade + 4.0_dp) / 2.0_dp)

      ALLOCATE(tau_nodes(n1D))
      ALLOCATE(tau_weights(n1D))

      CALL gauleg(0.0_dp, 1.0_dp, tau_nodes, tau_weights, n1D)

      !*************************************************************************
      ! Prodotto tensoriale
      !*************************************************************************

      n_tri = SIZE(tri_weights)
      n_tau = SIZE(tau_weights)

      ALLOCATE(XI_ref(n_tri * n_tau))
      ALLOCATE(ETA_ref(n_tri * n_tau))
      ALLOCATE(T_ref(n_tri * n_tau))
      ALLOCATE(W_ref(n_tri * n_tau))

      idx = 1
      DO k = 1, n_tau
         DO i = 1, n_tri
            XI_ref(idx)  = tri_nodes(i,1)
            ETA_ref(idx) = tri_nodes(i,2)
            T_ref(idx)   = tau_nodes(k)
            W_ref(idx)   = tau_weights(k) * tri_weights(i)
            idx = idx + 1
         END DO
      END DO

   END SUBROUTINE PrismDGL


   SUBROUTINE PrismGJCC(ade, XI_ref, ETA_ref, T_ref, W_ref)

      IMPLICIT NONE

      !*************************************************************************
      !
      ! Costruisce la regola di quadratura tensoriale Gauss--Jacobi--Clenshaw--
      ! Curtis (GJCC) sul prisma triangolare di riferimento
      !
      !    P_ref = { (xi,eta,tau) : xi>=0, eta>=0, xi+eta<=1, 0<=tau<=1 },
      !
      ! ottenuta come prodotto tensoriale tra:
      !
      !   - la regola di Gauss--Jacobi sul triangolo di riferimento;
      !
      !   - la regola di Clenshaw--Curtis sull'intervallo temporale [0,1].
      !
      !*************************************************************************
      ! Argomenti
      !*************************************************************************
      INTEGER, INTENT(IN)                :: ade
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: XI_ref(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: ETA_ref(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: T_ref(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: W_ref(:)
      !*************************************************************************
      ! Variabili locali
      !*************************************************************************
      INTEGER               :: nGP_tri
      INTEGER               :: n_tri, n_tau
      INTEGER               :: i, k, idx

      REAL(dp), ALLOCATABLE :: tri_nodes(:,:)
      REAL(dp), ALLOCATABLE :: tri_weights(:)
      REAL(dp), ALLOCATABLE :: tau_nodes(:)
      REAL(dp), ALLOCATABLE :: tau_weights(:)
      !*************************************************************************

      !*************************************************************************
      ! Quadratura spaziale
      !*************************************************************************

      nGP_tri = CEILING((ade + 2.0_dp) / 2.0_dp)

      n_tri = nGP_tri * nGP_tri

      ALLOCATE(tri_nodes(n_tri,2))
      ALLOCATE(tri_weights(n_tri))

      CALL TriangleQuadratureGJPoints(tri_nodes, tri_weights, n_tri, nGP_tri)

      ! I pesi della quadratura GJ sono gia' riferiti al triangolo di
      ! area 1/2: qui non serve alcun fattore di scala (a differenza di
      ! triangle_quadrature, dove i pesi vengono normalizzati ad area 1).

      !*************************************************************************
      ! Quadratura temporale
      !*************************************************************************

      CALL clenshaw_curtis(ade + 4, tau_nodes, tau_weights)

      !*************************************************************************
      ! Prodotto tensoriale
      !*************************************************************************

      n_tau = SIZE(tau_weights)

      ALLOCATE(XI_ref(n_tri * n_tau))
      ALLOCATE(ETA_ref(n_tri * n_tau))
      ALLOCATE(T_ref(n_tri * n_tau))
      ALLOCATE(W_ref(n_tri * n_tau))

      idx = 1
      DO k = 1, n_tau
         DO i = 1, n_tri
            XI_ref(idx)  = tri_nodes(i,1)
            ETA_ref(idx) = tri_nodes(i,2)
            T_ref(idx)   = tau_nodes(k)
            W_ref(idx)   = tau_weights(k) * tri_weights(i)
            idx = idx + 1
         END DO
      END DO

   END SUBROUTINE PrismGJCC


   SUBROUTINE PrismGJL(ade, XI_ref, ETA_ref, T_ref, W_ref)

      IMPLICIT NONE

      !*************************************************************************
      !
      ! Costruisce la regola di quadratura tensoriale Gauss--Jacobi--Legendre
      ! (GJL) sul prisma triangolare di riferimento
      !
      !    P_ref = { (xi,eta,tau) : xi>=0, eta>=0, xi+eta<=1, 0<=tau<=1 },
      !
      ! ottenuta come prodotto tensoriale tra:
      !
      !   - la regola di Gauss--Jacobi sul triangolo di riferimento
      !     (trasformazione di Duffy);
      !
      !   - la regola di Gauss--Legendre sull'intervallo temporale [0,1].
      !
      !*************************************************************************
      ! Argomenti
      !*************************************************************************
      INTEGER, INTENT(IN)                :: ade
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: XI_ref(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: ETA_ref(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: T_ref(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: W_ref(:)
      !*************************************************************************
      ! Variabili locali
      !*************************************************************************
      INTEGER               :: nGP_tri, n1D
      INTEGER               :: n_tri, n_tau
      INTEGER               :: i, k, idx

      REAL(dp), ALLOCATABLE :: tri_nodes(:,:)
      REAL(dp), ALLOCATABLE :: tri_weights(:)
      REAL(dp), ALLOCATABLE :: tau_nodes(:)
      REAL(dp), ALLOCATABLE :: tau_weights(:)
      !*************************************************************************

      !*************************************************************************
      ! Quadratura spaziale
      !*************************************************************************

      nGP_tri = CEILING((ade + 2.0_dp) / 2.0_dp)

      n_tri = nGP_tri * nGP_tri

      ALLOCATE(tri_nodes(n_tri,2))
      ALLOCATE(tri_weights(n_tri))

      CALL TriangleQuadratureGJPoints(tri_nodes, tri_weights, n_tri, nGP_tri)

      ! I pesi della quadratura GJ sono gia' riferiti al triangolo di
      ! area 1/2: qui non serve alcun fattore di scala.

      !*************************************************************************
      ! Quadratura temporale
      !*************************************************************************

      n1D = CEILING((ade + 4.0_dp) / 2.0_dp)

      ALLOCATE(tau_nodes(n1D))
      ALLOCATE(tau_weights(n1D))

      CALL gauleg(0.0_dp, 1.0_dp, tau_nodes, tau_weights, n1D)

      !*************************************************************************
      ! Prodotto tensoriale
      !*************************************************************************

      n_tau = SIZE(tau_weights)

      ALLOCATE(XI_ref(n_tri * n_tau))
      ALLOCATE(ETA_ref(n_tri * n_tau))
      ALLOCATE(T_ref(n_tri * n_tau))
      ALLOCATE(W_ref(n_tri * n_tau))

      idx = 1
      DO k = 1, n_tau
         DO i = 1, n_tri
            XI_ref(idx)  = tri_nodes(i,1)
            ETA_ref(idx) = tri_nodes(i,2)
            T_ref(idx)   = tau_nodes(k)
            W_ref(idx)   = tau_weights(k) * tri_weights(i)
            idx = idx + 1
         END DO
      END DO

   END SUBROUTINE PrismGJL

   SUBROUTINE clenshaw_curtis(n_tau, tau_nodes, tau_weights)
      
      IMPLICIT NONE

      !***********************************************************************
      ! Argomenti
      !***********************************************************************
      INTEGER, INTENT(IN)                :: n_tau
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: tau_nodes(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: tau_weights(:)
      !***********************************************************************
      ! Variabili locali
      !***********************************************************************
      INTEGER  :: N_cc, i, jj
      REAL(dp) :: theta_i, sum_w, fact_j, g_cc
      !***********************************************************************

      N_cc = n_tau - 1

      ALLOCATE(tau_nodes(n_tau))
      ALLOCATE(tau_weights(n_tau))

      DO i = 1, n_tau

         theta_i = PI * REAL(i - 1, dp) / REAL(N_cc, dp)

         sum_w = 0.0_dp

         DO jj = 1, N_cc / 2

            fact_j = 2.0_dp
            IF (MOD(N_cc, 2) == 0 .AND. jj == N_cc / 2) THEN
               fact_j = 1.0_dp
            END IF

            sum_w = sum_w + fact_j * COS(2.0_dp * REAL(jj, dp) * theta_i) / &
                        (1.0_dp - 4.0_dp * REAL(jj, dp)**2)

         END DO

         g_cc = 1.0_dp
         IF (i == 1 .OR. i == n_tau) g_cc = 0.5_dp

         ! Flip rispetto all'ordine naturale dei nodi Chebyshev:
         ! tau cresce da 0 a 1.
         tau_nodes(n_tau - i + 1) = 0.5_dp * (COS(theta_i) + 1.0_dp)
         tau_weights(n_tau - i + 1) = (1.0_dp / REAL(N_cc, dp)) * g_cc * (1.0_dp + sum_w)

      END DO

    END SUBROUTINE clenshaw_curtis

   SUBROUTINE ShiftingPrismQuadrature(V, XI_ref, ETA_ref, T_ref, W_ref, XYZTW, WV_X)

      IMPLICIT NONE

      !*************************************************************************
      !
      ! Costruisce i nodi fisici e i pesi orientati per integrare su una
      ! iperfaccia laterale spazio-temporale 3D immersa in R^4.
      !
      ! Il prisma di riferimento e':
      !
      !    P = {(xi,eta,t) : xi >= 0, eta >= 0, xi + eta <= 1, 0 <= t <= 1}.
      !
      ! L'iperfaccia e' ottenuta collegando linearmente due configurazioni
      ! triangolari corrispondenti mediante la parametrizzazione
      !
      !    X(xi,eta,tau) = (1-tau) X_0(xi,eta) + tau X_1(xi,eta),
      !
      ! con X_0 e X_1 parametrizzazioni affini delle facce triangolari a
      ! tau = 0 e tau = 1:
      !
      !    X_0(xi,eta) = A_0 + xi (B_0-A_0) + eta (C_0-A_0),
      !    X_1(xi,eta) = A_1 + xi (B_1-A_1) + eta (C_1-A_1).
      !
      ! Poiche' l'ultima coordinata della mappa coincide con tau, i vettori
      ! tangenti hanno la forma
      !
      !    X_xi = (x_xi, y_xi, z_xi, 0),  X_eta = (x_eta, y_eta, z_eta, 0),
      !    X_tau = (x_tau, y_tau, z_tau, 1),
      !
      ! e la componente x della normale generalizzata (minore 3x3 ottenuto
      ! eliminando la prima coordinata) si riduce a
      !
      !    N_x = y_xi z_eta - z_xi y_eta.
      !
      ! I pesi restituiti includono il fattore geometrico orientato:
      !
      !    WV_X = W_ref .* N_x
      !
      ! Implementazione: tutte le grandezze costanti (lati, differenze tra
      ! configurazione finale e iniziale) sono precalcolate come scalari; le
      ! interpolazioni lineari sono scritte nella forma
      !
      !    (1-tau) a + tau b = a + tau (b - a),
      !
      ! e i cicli sui Nq nodi sono contigui in memoria e privi di array
      ! temporanei, cosi' da poter essere vettorizzati dal compilatore.
      !
      !*************************************************************************
      ! Argomenti
      !*************************************************************************
      REAL(dp), INTENT(IN)               :: V(6,4)
      REAL(dp), INTENT(IN)               :: XI_ref(:)
      REAL(dp), INTENT(IN)               :: ETA_ref(:)
      REAL(dp), INTENT(IN)               :: T_ref(:)
      REAL(dp), INTENT(IN)               :: W_ref(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: XYZTW(:,:)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: WV_X(:)
      !*************************************************************************
      ! Variabili locali
      !*************************************************************************
      INTEGER  :: num_pts, q, j
      REAL(dp) :: t
      ! Vettori lato delle facce a tau = 0 (base) e tau = 1 (top)
      REAL(dp) :: lato1_base(4), lato2_base(4)
      REAL(dp) :: lato1_top(4),  lato2_top(4)

      ! Differenze top - base (spostamento dei lati e del vertice A)
      REAL(dp) :: dA(4), dL1(4), dL2(4)

      ! Coefficienti delle derivate tangenti in y, z:  d = a + tau * b
      REAL(dp) :: ay_xi, by_xi, az_xi, bz_xi
      REAL(dp) :: ay_eta, by_eta, az_eta, bz_eta

      REAL(dp) :: dy_dxi, dz_dxi, dy_deta, dz_deta
      !*************************************************************************

      num_pts = SIZE(W_ref)

      ALLOCATE(XYZTW(num_pts,4))
      ALLOCATE(WV_X(num_pts))

      !*************************************************************************
      ! Geometria costante
      !*************************************************************************

      ! Vettori lato della faccia iniziale (tau = 0) e finale (tau = 1)
      lato1_base = V(2,:) - V(1,:)
      lato2_base = V(3,:) - V(1,:)
      lato1_top  = V(5,:) - V(4,:)
      lato2_top  = V(6,:) - V(4,:)

      ! Poiche' X_1 - X_0 = dA + xi * dL1 + eta * dL2, il vettore tangente
      ! rispetto al tempo e' calcolabile direttamente senza costruire
      ! esplicitamente X_0 e X_1
      dA  = V(4,:) - V(1,:)
      dL1 = lato1_top - lato1_base
      dL2 = lato2_top - lato2_base

      ! Derivate tangenti (solo y,z, sufficienti per N_x):
      !    dX/dxi  = lato1_base + tau (lato1_top - lato1_base)
      !    dX/deta = lato2_base + tau (lato2_top - lato2_base)
      ay_xi  = lato1_base(2);  by_xi  = dL1(2)
      az_xi  = lato1_base(3);  bz_xi  = dL1(3)
      ay_eta = lato2_base(2);  by_eta = dL2(2)
      az_eta = lato2_base(3);  bz_eta = dL2(3)

      !*************************************************************************
      ! Nodi fisici
      !
      !    X(xi,eta,tau) = X_0(xi,eta) + tau * (X_1(xi,eta) - X_0(xi,eta))
      !
      ! Ciclo esterno sulle 4 coordinate, ciclo interno sui nodi: accesso
      ! contiguo alle colonne di XYZTW (ordinamento column-major).
      !*************************************************************************

      DO j = 1, 4
         DO q = 1, num_pts
            XYZTW(q,j) = V(1,j) + XI_ref(q) * lato1_base(j) + ETA_ref(q) * lato2_base(j) &
                       + T_ref(q) * (dA(j) + XI_ref(q) * dL1(j) + ETA_ref(q) * dL2(j))
         END DO
      END DO

      !*************************************************************************
      ! Pesi orientati
      !
      !    N_x = y_xi z_eta - z_xi y_eta,     WV_X = W_ref * N_x
      !*************************************************************************

      DO q = 1, num_pts
         t = T_ref(q)

         dy_dxi  = ay_xi  + t * by_xi
         dz_dxi  = az_xi  + t * bz_xi
         dy_deta = ay_eta + t * by_eta
         dz_deta = az_eta + t * bz_eta

         WV_X(q) = W_ref(q) * (dy_dxi * dz_deta - dz_dxi * dy_deta)
      END DO

   END SUBROUTINE ShiftingPrismQuadrature

END MODULE PrismQuadrature