MODULE OPC4D_MoL_Module

   USE TypesDef
   USE ReferenceFunctions
   USE CubatureFunctions
   USE TimeDiscretization

   IMPLICIT NONE

CONTAINS

   SUBROUTINE OPC4D_MoL(ade, n_tau, vertices_initial, vertices_final, facets, method, XYZtau, W)

      IMPLICIT NONE

      INTEGER, INTENT(IN)                         :: ade
      INTEGER, INTENT(IN)                         :: n_tau
      REAL(dp), INTENT(IN)                        :: vertices_initial(:,:)
      REAL(dp), INTENT(IN)                        :: vertices_final(:,:)
      INTEGER, INTENT(IN)                         :: facets(:,:)
      CHARACTER(LEN=*), INTENT(IN)                :: method
      REAL(dp), ALLOCATABLE, INTENT(OUT)          :: XYZtau(:,:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)          :: W(:)

      INTEGER                                     :: n_spaz
      INTEGER                                     :: n_mom
      INTEGER                                     :: k
      INTEGER                                     :: first
      INTEGER                                     :: last
      INTEGER                                     :: ind_curr(3)

      REAL(dp), ALLOCATABLE                       :: tau_nodes(:)
      REAL(dp), ALLOCATABLE                       :: tau_weights(:)
      REAL(dp), ALLOCATABLE                       :: xyzw_tens_ref(:,:)
      REAL(dp), ALLOCATABLE                       :: xyzw_tens(:,:)
      REAL(dp), ALLOCATABLE                       :: x_ref(:,:)
      REAL(dp), ALLOCATABLE                       :: v_ref(:,:)
      REAL(dp), ALLOCATABLE                       :: coeffs(:)
      REAL(dp), ALLOCATABLE                       :: moments_ch(:)
      REAL(dp), ALLOCATABLE                       :: alpha(:)
      REAL(dp), ALLOCATABLE                       :: w_spatial(:)
      REAL(dp), ALLOCATABLE                       :: vertices_tau(:,:)
      INTEGER, ALLOCATABLE                        :: chebyshev_indices(:,:)
      REAL(dp)                                    :: bbox_tau(6)
      REAL(dp)                                    :: tau

      IF (TRIM(method) == 'D' .AND. ade > 20) THEN
         WRITE(*,'(A)') 'WARNING: Le regole di Dunavant sono disponibili solo fino al grado 20.'
         RETURN
      END IF

      CALL ComputeClenshawCurtis(n_tau, tau_nodes, tau_weights)
      CALL cub_gausscheb_tens3D(2 * ade, xyzw_tens_ref)

      n_spaz = SIZE(xyzw_tens_ref, 1)
      n_mom = ((ade + 1) * (ade + 2) * (ade + 3)) / 6

      ALLOCATE(chebyshev_indices(n_mom, 3))
      chebyshev_indices(1,:) = 0
      DO k = 2, n_mom
         ind_curr = chebyshev_indices(k - 1, :)
         CALL mono_next_grlex(3, ind_curr)
         chebyshev_indices(k, :) = ind_curr
      END DO

      ALLOCATE(x_ref(n_spaz, 3))
      x_ref(:,1) = xyzw_tens_ref(:,1)
      x_ref(:,2) = xyzw_tens_ref(:,2)
      x_ref(:,3) = xyzw_tens_ref(:,3)

      CALL dCHEBVAND(ade, x_ref, chebyshev_indices, v_ref)
      CALL tenscheb_norm2sq(chebyshev_indices, coeffs)

      ALLOCATE(vertices_tau(SIZE(vertices_initial, 1), 3))
      ALLOCATE(alpha(n_mom))
      ALLOCATE(w_spatial(n_spaz))

      IF (ALLOCATED(XYZtau)) DEALLOCATE(XYZtau)
      IF (ALLOCATED(W)) DEALLOCATE(W)
      ALLOCATE(XYZtau(n_spaz * n_tau, 4))
      ALLOCATE(W(n_spaz * n_tau))

      DO k = 1, n_tau
         tau = tau_nodes(k)
         vertices_tau = (1.0_dp - tau) * vertices_initial + tau * vertices_final

         bbox_tau(1) = MINVAL(vertices_tau(:,1))
         bbox_tau(2) = MAXVAL(vertices_tau(:,1))
         bbox_tau(3) = MINVAL(vertices_tau(:,2))
         bbox_tau(4) = MAXVAL(vertices_tau(:,2))
         bbox_tau(5) = MINVAL(vertices_tau(:,3))
         bbox_tau(6) = MAXVAL(vertices_tau(:,3))

         CALL chebyshev_moments_polyhedron(vertices_tau, facets, ade, chebyshev_indices, bbox_tau, method, moments_ch)
         CALL scale_rule(xyzw_tens_ref, bbox_tau, xyzw_tens)

         alpha = moments_ch / coeffs
         w_spatial = xyzw_tens(:,4) * MATMUL(v_ref, alpha)

         first = (k - 1) * n_spaz + 1
         last = k * n_spaz

         XYZtau(first:last, 1:3) = xyzw_tens(:, 1:3)
         XYZtau(first:last, 4) = tau
         W(first:last) = tau_weights(k) * w_spatial
      END DO

   END SUBROUTINE OPC4D_MoL

END MODULE OPC4D_MoL_Module
