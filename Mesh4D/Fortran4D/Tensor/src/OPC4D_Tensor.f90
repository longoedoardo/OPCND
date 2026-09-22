MODULE OPC4D_Tensor_Module

   USE TypesDef
   USE ReferenceFunctions
   USE CubatureTensor

   IMPLICIT NONE

CONTAINS

   SUBROUTINE OPC4D_Tensor(ade, vertices_initial, vertices_final, facets, method, XYZT, W)

      !***********************************************************************
      !
      ! SUBROUTINE OPC4D_Tensor(ade, vertices_initial, vertices_final, facets, method, XYZT, W)
      !
      ! Calcola i nodi e i pesi di una regola di cubatura quadridimensionale
      ! su un dominio poliedrico in movimento.
      ! Il dominio spazio-temporale e' il politopo 4D descritto da un poliedro
      ! 3D che si deforma linearmente, per tempo adimensionale tau in [0,1],
      ! dalla configurazione iniziale a quella finale, mantenendo invariata
      ! la connettivita' della mesh superficiale triangolare.
      !
      ! Il metodo non richiede una decomposizione volumetrica del dominio 4D,
      ! i momenti polinomiali vengono ricondotti, mediante il teorema della
      ! divergenza, a integrali sulle ipersuperfici laterali 3D (prismi
      ! triangolari spazio-temporali generati da ciascuna faccia).
      !
      !***********************************************************************
      !
      ! INPUTS:
      !
      !   - ade:
      !       Grado polinomiale totale massimo della regola di cubatura.
      !       Il numero di momenti utilizzati e'
      !       N_mom = (ade+1)(ade+2)(ade+3)(ade+4)/24.
      !
      !   - vertices_initial:
      !       Matrice N x 3 con le coordinate [x_i, y_i, z_i] dei vertici
      !       della mesh superficiale all'istante iniziale (tau = 0).
      !
      !   - vertices_final:
      !       Matrice N x 3 con le coordinate dei medesimi vertici
      !       all'istante finale (tau = 1). Deve avere lo stesso numero di
      !       righe e lo stesso ordinamento di vertices_initial.
      !
      !   - facets:
      !       Matrice M x 3 con la connettivita' della mesh triangolare.
      !       Le facce devono formare una superficie chiusa (ogni spigolo
      !       condiviso da esattamente due facce) e orientata coerentemente
      !       verso l'esterno: l'orientazione determina il segno dei
      !       contributi al teorema della divergenza.
      !
      !   - method:
      !       Metodo di quadratura sul prisma spazio-temporale di ogni faccia:
      !       'DCC'  : Dunavant sul triangolo + Clenshaw-Curtis nel tempo;
      !       'DGL'  : Dunavant sul triangolo + Gauss-Legendre nel tempo;
      !       'GJCC' : Gauss-Jacobi sul triangolo + Clenshaw-Curtis nel tempo;
      !       'GJL'  : Gauss-Jacobi sul triangolo + Gauss-Legendre nel tempo.
      !
      !***********************************************************************
      !
      ! OUTPUTS:
      !
      !   - XYZT:
      !       Matrice n_nodes x 4 con le coordinate [x_k, y_k, z_k, tau_k]
      !       dei nodi di cubatura nel dominio fisico spazio-temporale.
      !
      !   - W:
      !       Vettore n_nodes x 1 con i pesi associati ai nodi.
      !
      !***********************************************************************

      IMPLICIT NONE

      !**********************************************************************
      ! Argomenti
      !**********************************************************************
      INTEGER, INTENT(IN)                         :: ade
      REAL(dp), INTENT(IN)                        :: vertices_initial(:,:)
      REAL(dp), INTENT(IN)                        :: vertices_final(:,:)
      INTEGER, INTENT(IN)                         :: facets(:,:)
      CHARACTER(LEN=*), INTENT(IN)                :: method
      REAL(dp), ALLOCATABLE, INTENT(OUT)          :: XYZT(:,:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)          :: W(:)
      !**********************************************************************
      ! Variabili locali
      !**********************************************************************
      INTEGER                                     :: n_nodes, n_vertices, n_mom, k
      INTEGER                                     :: ind_curr(4)
      REAL(dp), ALLOCATABLE                       :: XYZTW_ref(:,:), XYZTW(:,:), X_ref(:,:)
      REAL(dp), ALLOCATABLE                       :: V_ref(:,:), coeffs(:), moments_ch(:), alpha(:)
      REAL(dp), ALLOCATABLE                       :: vertices_4d(:,:)
      INTEGER, ALLOCATABLE                        :: chebyshev_indices(:,:)
      REAL(dp)                                    :: bbox(2,4)
      !**********************************************************************

      !**********************************************************************
      ! INIZIO PARTE SHAPE-INDEPENDENT
      !**********************************************************************

      ! Griglia tensoriale di Gauss-Chebyshev sull'ipercubo di riferimento [-1,1]^4
      CALL cub_gausscheb_tens4D(2 * ade, XYZTW_ref)
      n_nodes = SIZE(XYZTW_ref, 1)

      ! Dimensione dello spazio polinomiale 4D di grado totale <= ade
      n_mom = ((ade + 1) * (ade + 2) * (ade + 3) * (ade + 4)) / 24

      ! Indici (i,j,k,l) dei polinomi tensoriali di Chebyshev, ordinamento GRLEX
      ALLOCATE(chebyshev_indices(n_mom, 4))
      chebyshev_indices(1,:) = 0
      DO k = 2, n_mom
         ind_curr = chebyshev_indices(k - 1, :)
         CALL mono_next_grlex(4, ind_curr)
         chebyshev_indices(k,:) = ind_curr
      END DO

      ! Matrice di Vandermonde-Chebyshev 4D valutata sui nodi di riferimento
      ALLOCATE(X_ref(n_nodes, 4))
      X_ref = XYZTW_ref(:,1:4)
      CALL dCHEBVAND(ade, X_ref, chebyshev_indices, V_ref)

      ! Norme quadre (pesate) dei polinomi di Chebyshev tensoriali
      CALL tenscheb_norm2sq(chebyshev_indices, coeffs)

      !**********************************************************************
      ! INIZIO PARTE SHAPE-DEPENDENT
      !**********************************************************************

      ! Assemblaggio dei vertici 4D: le prime n_vertices righe sono la
      ! configurazione iniziale (tau = 0), le successive quella finale (tau = 1)
      n_vertices = SIZE(vertices_initial, 1)
      ALLOCATE(vertices_4d(2 * n_vertices, 4))
      vertices_4d(1:n_vertices, 1:3) = vertices_initial
      vertices_4d(1:n_vertices, 4) = 0.0_dp
      vertices_4d(n_vertices+1:2*n_vertices, 1:3) = vertices_final
      vertices_4d(n_vertices+1:2*n_vertices, 4) = 1.0_dp

      ! Bounding box 4D: [xmin,ymin,zmin,tmin ; xmax,ymax,zmax,tmax]
      DO k = 1, 4
         bbox(1,k) = MINVAL(vertices_4d(:,k))
         bbox(2,k) = MAXVAL(vertices_4d(:,k))
      END DO

      ! Momenti della base di Chebyshev sul politopo 4D
      CALL chebyshev_moments_polyhedron_4D(vertices_4d, facets, ade, chebyshev_indices, bbox, method, moments_ch)

      ! Riscalamento dei nodi di riferimento sulla bounding box reale
      CALL scale_rule(XYZTW_ref, bbox, XYZTW)

      !**********************************************************************
      ! Costruzione dei pesi di cubatura
      !**********************************************************************

      ALLOCATE(alpha(n_mom), W(n_nodes), XYZT(n_nodes, 4))

      alpha = moments_ch / coeffs

      W = XYZTW(:,5) * MATMUL(V_ref, alpha)
      
      XYZT = XYZTW(:,1:4)

   END SUBROUTINE OPC4D_Tensor

END MODULE OPC4D_Tensor_Module