MODULE CubatureTensor4D

   USE TypesDef
   USE ReferenceFunctions4D, ONLY: dCHEBVAND, tenscheb_norm2sq, scale_rule
   USE PrismQuadrature

   IMPLICIT NONE

CONTAINS

   SUBROUTINE chebyshev_moments_polyhedron_4D(vertices_4d, facets, ade, chebyshev_indices, dbox, method, moments)
      REAL(dp), INTENT(IN)                        :: vertices_4d(:,:)
      INTEGER, INTENT(IN)                         :: facets(:,:)
      INTEGER, INTENT(IN)                         :: ade
      INTEGER, INTENT(IN)                         :: chebyshev_indices(:,:)
      REAL(dp), INTENT(IN)                        :: dbox(2,4)
      CHARACTER(LEN=*), INTENT(IN)                :: method
      REAL(dp), ALLOCATABLE, INTENT(OUT)          :: moments(:)

      INTEGER                                     :: n_mom, n_facets, n_vertices
      INTEGER                                     :: f, n1, n2, n3
      REAL(dp)                                    :: V_prism(6,4)
      REAL(dp), ALLOCATABLE                       :: XI_ref(:), ETA_ref(:), T_ref(:), W_ref(:)
      REAL(dp), ALLOCATABLE                       :: XYZTW(:,:), WV_X(:), facet_moments(:)

      n_mom = SIZE(chebyshev_indices, 1)
      n_facets = SIZE(facets, 1)
      n_vertices = SIZE(vertices_4d, 1) / 2

      ALLOCATE(moments(n_mom))
      moments = 0.0_dp

      CALL PrismQuadratureRule(ade, method, XI_ref, ETA_ref, T_ref, W_ref)

      DO f = 1, n_facets
         ALLOCATE(XYZTW(SIZE(W_ref), 4), WV_X(SIZE(W_ref)))
         n1 = facets(f,1)
         n2 = facets(f,2)
         n3 = facets(f,3)

         V_prism(1,:) = vertices_4d(n1, :)
         V_prism(2,:) = vertices_4d(n2, :)
         V_prism(3,:) = vertices_4d(n3, :)
         V_prism(4,:) = vertices_4d(n1 + n_vertices, :)
         V_prism(5,:) = vertices_4d(n2 + n_vertices, :)
         V_prism(6,:) = vertices_4d(n3 + n_vertices, :)

         CALL ShiftingPrismQuadrature(V_prism, XI_ref, ETA_ref, T_ref, W_ref, XYZTW, WV_X)

         IF (MAXVAL(ABS(WV_X)) > 1.0e-14_dp) THEN

            CALL cubature_tens_chebyshev_facet_4D(XYZTW, WV_X, chebyshev_indices, dbox, facet_moments)

            moments = moments + facet_moments

            DEALLOCATE(facet_moments)

         END IF

         DEALLOCATE(XYZTW)
         DEALLOCATE(WV_X)

      END DO


   END SUBROUTINE chebyshev_moments_polyhedron_4D

   SUBROUTINE cubature_tens_chebyshev_facet_4D(XYZTW, WV_X, chebyshev_indices, dbox, chebyshev_moms)
      REAL(dp), INTENT(IN)                        :: XYZTW(:,:)
      REAL(dp), INTENT(IN)                        :: WV_X(:)
      INTEGER, INTENT(IN)                         :: chebyshev_indices(:,:)
      REAL(dp), INTENT(IN)                        :: dbox(2,4)
      REAL(dp), ALLOCATABLE, INTENT(OUT)          :: chebyshev_moms(:)

      INTEGER                                     :: n, m, c, iv
      INTEGER                                     :: max_i, max_j, max_k, max_l
      INTEGER, ALLOCATABLE                        :: idx_i(:), idx_j(:), idx_k(:), idx_l(:)
      REAL(dp)                                    :: B1
      REAL(dp), ALLOCATABLE                       :: XN(:), YN(:), ZN(:), TN(:)
      REAL(dp), ALLOCATABLE                       :: TX(:,:), TY(:,:), TZ(:,:), TT(:,:), IntX(:,:)

      n = SIZE(XYZTW, 1)
      m = SIZE(chebyshev_indices, 1)

      ALLOCATE(XN(n), YN(n), ZN(n), TN(n))
      XN = (2.0_dp * XYZTW(:,1) - dbox(1,1) - dbox(2,1)) / (dbox(2,1) - dbox(1,1))
      YN = (2.0_dp * XYZTW(:,2) - dbox(1,2) - dbox(2,2)) / (dbox(2,2) - dbox(1,2))
      ZN = (2.0_dp * XYZTW(:,3) - dbox(1,3) - dbox(2,3)) / (dbox(2,3) - dbox(1,3))
      TN = (2.0_dp * XYZTW(:,4) - dbox(1,4) - dbox(2,4)) / (dbox(2,4) - dbox(1,4))

      max_i = MAXVAL(chebyshev_indices(:,1))
      max_j = MAXVAL(chebyshev_indices(:,2))
      max_k = MAXVAL(chebyshev_indices(:,3))
      max_l = MAXVAL(chebyshev_indices(:,4))

      ALLOCATE(TX(n, max_i + 2), TY(n, max_j + 1), TZ(n, max_k + 1), TT(n, max_l + 1))
      CALL chebpolys(max_i + 1, XN, TX)
      CALL chebpolys(max_j, YN, TY)
      CALL chebpolys(max_k, ZN, TZ)
      CALL chebpolys(max_l, TN, TT)

      ALLOCATE(IntX(n, max_i + 1))
      IntX(:,1) = XN
      IF (max_i >= 1) IntX(:,2) = 0.5_dp * XN**2
      DO iv = 2, max_i
         IntX(:, iv + 1) = TX(:,iv+2) / (2.0_dp * REAL(iv+1, dp)) - TX(:,iv) / (2.0_dp * REAL(iv-1, dp))
      END DO

      ALLOCATE(idx_i(m), idx_j(m), idx_k(m), idx_l(m))
      idx_i = chebyshev_indices(:,1) + 1
      idx_j = chebyshev_indices(:,2) + 1
      idx_k = chebyshev_indices(:,3) + 1
      idx_l = chebyshev_indices(:,4) + 1

      B1 = 0.5_dp * (dbox(2,1) - dbox(1,1))
      ALLOCATE(chebyshev_moms(m))

      DO c = 1, m
         chebyshev_moms(c) = B1 * SUM(WV_X * IntX(:,idx_i(c)) * TY(:,idx_j(c)) * TZ(:,idx_k(c)) * TT(:,idx_l(c)))
      END DO
   END SUBROUTINE cubature_tens_chebyshev_facet_4D

   SUBROUTINE chebpolys(deg, x, T)
      INTEGER, INTENT(IN)                         :: deg
      REAL(dp), INTENT(IN)                        :: x(:)
      REAL(dp), INTENT(OUT)                       :: T(:,:)

      INTEGER                                     :: j

      T(:,1) = 1.0_dp
      IF (deg >= 1) T(:,2) = x
      DO j = 2, deg
         T(:,j+1) = 2.0_dp * x * T(:,j) - T(:,j-1)
      END DO
   END SUBROUTINE chebpolys

END MODULE CubatureTensor4D
