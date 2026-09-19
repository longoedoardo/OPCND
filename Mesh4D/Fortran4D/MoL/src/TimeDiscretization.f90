MODULE TimeDiscretization
    USE TypesDef
    IMPLICIT NONE

CONTAINS

    SUBROUTINE ComputeClenshawCurtis(n_tau, tau_nodes, tau_weights)
        IMPLICIT NONE
        INTEGER, INTENT(IN)                :: n_tau
        REAL(dp), ALLOCATABLE, INTENT(OUT) :: tau_nodes(:)
        REAL(dp), ALLOCATABLE, INTENT(OUT) :: tau_weights(:)

        INTEGER  :: N_cc, i, jj
        REAL(dp) :: theta_i, sum_w, fact_j, g_cc

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
            tau_weights(n_tau - i + 1) = &
                (1.0_dp / REAL(N_cc, dp)) * g_cc * (1.0_dp + sum_w)

        END DO

    END SUBROUTINE ComputeClenshawCurtis

END MODULE TimeDiscretization