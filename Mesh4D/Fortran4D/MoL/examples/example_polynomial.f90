PROGRAM example_polynomial

   USE TypesDef
   USE OPC4D_MoL_Module

   IMPLICIT NONE

   REAL(dp)                            :: vertices_initial(8,3)
   REAL(dp)                            :: vertices_final(8,3)
   INTEGER                             :: facets(12,3)
   REAL(dp), ALLOCATABLE               :: XYZtau(:,:)
   REAL(dp), ALLOCATABLE               :: W(:)
   INTEGER                             :: n_tau

   n_tau = 10

   vertices_initial = RESHAPE([ &
      -1.0_dp, -1.0_dp, -1.0_dp, &
       1.0_dp, -1.0_dp, -1.0_dp, &
       1.0_dp,  1.0_dp, -1.0_dp, &
      -1.0_dp,  1.0_dp, -1.0_dp, &
      -1.0_dp, -1.0_dp,  1.0_dp, &
       1.0_dp, -1.0_dp,  1.0_dp, &
       1.0_dp,  1.0_dp,  1.0_dp, &
      -1.0_dp,  1.0_dp,  1.0_dp], [8,3], ORDER=[2,1])

   vertices_final = RESHAPE([ &
      -0.90_dp, -0.95_dp, -0.90_dp, &
       1.10_dp, -0.95_dp, -0.90_dp, &
       1.10_dp,  1.05_dp, -0.90_dp, &
      -0.90_dp,  1.05_dp, -0.90_dp, &
      -0.90_dp, -0.95_dp,  1.10_dp, &
       1.10_dp, -0.95_dp,  1.10_dp, &
       1.10_dp,  1.05_dp,  1.10_dp, &
      -0.90_dp,  1.05_dp,  1.10_dp], [8,3], ORDER=[2,1])

   facets = RESHAPE([ &
      1, 3, 2,  1, 4, 3, &
      5, 6, 7,  5, 7, 8, &
      1, 2, 6,  1, 6, 5, &
      4, 8, 7,  4, 7, 3, &
      1, 5, 8,  1, 8, 4, &
      2, 3, 7,  2, 7, 6], [12,3], ORDER=[2,1])

   WRITE(*,'(/,A)') '**************************************************************'
   WRITE(*,'(A)')   '                       OPC4D_MoL'
   WRITE(*,'(A)')   '                 Metodo delle Linee'
   WRITE(*,'(A)')   '            Cubatura su Cubo in Movimento'
   WRITE(*,'(A,/)') '**************************************************************'
   WRITE(*,'(A,I0)') 'Numero di vertici:     ', SIZE(vertices_initial, 1)
   WRITE(*,'(A,I0)') 'Numero di facce:       ', SIZE(facets, 1)
   WRITE(*,'(A,I0)') 'Numero di nodi tau:    ', n_tau

   CALL run_case(1, 'f_1(x,y,z,tau) = 1', 8.0_dp, 1)
   CALL run_case(2, 'f_2(x,y,z,tau) = x^2 + y^2 + z^2 + tau^2', 8.0_dp + 8.18_dp / 3.0_dp, 2)
   CALL run_case(6, 'f_3(x,y,z,tau) = x^2*y^2*z^2*tau^2', 29150263.0_dp / 283500000.0_dp, 3)

CONTAINS

   SUBROUTINE run_case(ade, label, exact, which_fun)
      INTEGER, INTENT(IN)             :: ade
      CHARACTER(LEN=*), INTENT(IN)    :: label
      REAL(dp), INTENT(IN)            :: exact
      INTEGER, INTENT(IN)             :: which_fun

      INTEGER                         :: c0
      INTEGER                         :: c1
      INTEGER                         :: rate
      REAL(dp)                        :: elapsed
      REAL(dp)                        :: integral
      REAL(dp)                        :: err_abs

      WRITE(*,'(/,A)') '--------------------------------------------------------------'
      WRITE(*,'(A,I0)') 'Funzione integranda ', which_fun
      WRITE(*,'(A)') '--------------------------------------------------------------'
      WRITE(*,'(A,I0)')      'Ade:                   ', ade
      WRITE(*,'(A,I0)')      'Numero di nodi tau:    ', n_tau
      WRITE(*,'(A,A)')       'Funzione integranda:   ', label
      WRITE(*,'(A,ES22.15)') 'Integrale esatto:      ', exact
      WRITE(*,'(/,A)')       'Inizio Cubatura...'

      CALL SYSTEM_CLOCK(c0, rate)
      CALL OPC4D_MoL(ade, n_tau, vertices_initial, vertices_final, facets, 'GJ', XYZtau, W)
      CALL SYSTEM_CLOCK(c1, rate)
      elapsed = REAL(c1 - c0, dp) / REAL(rate, dp)

      SELECT CASE (which_fun)
      CASE (1)
         integral = SUM(W)
      CASE (2)
         integral = SUM(W * (XYZtau(:,1)**2 + XYZtau(:,2)**2 + XYZtau(:,3)**2 + XYZtau(:,4)**2))
      CASE (3)
         integral = SUM(W * (XYZtau(:,1)**2 * XYZtau(:,2)**2 * XYZtau(:,3)**2 * XYZtau(:,4)**2))
      END SELECT

      err_abs = ABS(integral - exact)

      WRITE(*,'(A)')            'Fine Cubatura...'
      WRITE(*,'(/,A,I0)')       'Numero di nodi 4D:    ', SIZE(XYZtau, 1)
      WRITE(*,'(A,ES22.15)')    'Integrale numerico:   ', integral
      WRITE(*,'(A,ES12.6)')     'Errore assoluto:      ', err_abs
      WRITE(*,'(A,ES12.6,A)')   'Tempo di calcolo:     ', elapsed, ' s'
   END SUBROUTINE run_case

END PROGRAM example_polynomial
