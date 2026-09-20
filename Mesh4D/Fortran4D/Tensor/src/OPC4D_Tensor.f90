MODULE OPC4D_Tensor_Module

   USE TypesDef
   USE ReferenceFunctions4D
   USE CubatureTensor4D

   IMPLICIT NONE

CONTAINS

   SUBROUTINE OPC4D_Tensor(ade, vertices_initial, vertices_final, facets, method, XYZT, W)
      INTEGER, INTENT(IN)                         :: ade
      REAL(dp), INTENT(IN)                        :: vertices_initial(:,:)
      REAL(dp), INTENT(IN)                        :: vertices_final(:,:)
      INTEGER, INTENT(IN)                         :: facets(:,:)
      CHARACTER(LEN=*), INTENT(IN)                :: method
      REAL(dp), ALLOCATABLE, INTENT(OUT)          :: XYZT(:,:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)          :: W(:)

      INTEGER                                     :: n_nodes, n_vertices, n_mom, k
      INTEGER                                     :: ind_curr(4)
      REAL(dp), ALLOCATABLE                       :: XYZTW_ref(:,:), XYZTW(:,:), X_ref(:,:)
      REAL(dp), ALLOCATABLE                       :: V_ref(:,:), coeffs(:), moments_ch(:), alpha(:)
      REAL(dp), ALLOCATABLE                       :: vertices_4d(:,:)
      INTEGER, ALLOCATABLE                        :: chebyshev_indices(:,:)
      REAL(dp)                                    :: bbox(2,4)

      CALL cub_gausscheb_tens4D(2 * ade, XYZTW_ref)
      n_nodes = SIZE(XYZTW_ref, 1)

      n_mom = ((ade + 1) * (ade + 2) * (ade + 3) * (ade + 4)) / 24
      ALLOCATE(chebyshev_indices(n_mom, 4))
      chebyshev_indices(1,:) = 0
      DO k = 2, n_mom
         ind_curr = chebyshev_indices(k - 1, :)
         CALL mono_next_grlex(4, ind_curr)
         chebyshev_indices(k,:) = ind_curr
      END DO

      ALLOCATE(X_ref(n_nodes, 4))
      X_ref = XYZTW_ref(:,1:4)
      CALL dCHEBVAND(ade, X_ref, chebyshev_indices, V_ref)
      CALL tenscheb_norm2sq(chebyshev_indices, coeffs)

      n_vertices = SIZE(vertices_initial, 1)
      ALLOCATE(vertices_4d(2 * n_vertices, 4))
      vertices_4d(1:n_vertices, 1:3) = vertices_initial
      vertices_4d(1:n_vertices, 4) = 0.0_dp
      vertices_4d(n_vertices+1:2*n_vertices, 1:3) = vertices_final
      vertices_4d(n_vertices+1:2*n_vertices, 4) = 1.0_dp

      DO k = 1, 4
         bbox(1,k) = MINVAL(vertices_4d(:,k))
         bbox(2,k) = MAXVAL(vertices_4d(:,k))
      END DO

      CALL chebyshev_moments_polyhedron_4D(vertices_4d, facets, ade, chebyshev_indices, bbox, method, moments_ch)
      CALL scale_rule(XYZTW_ref, bbox, XYZTW)

      ALLOCATE(alpha(n_mom), W(n_nodes), XYZT(n_nodes, 4))
      alpha = moments_ch / coeffs
      W = XYZTW(:,5) * MATMUL(V_ref, alpha)
      XYZT = XYZTW(:,1:4)
   END SUBROUTINE OPC4D_Tensor

END MODULE OPC4D_Tensor_Module
