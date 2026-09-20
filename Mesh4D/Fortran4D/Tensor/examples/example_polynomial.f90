PROGRAM example_polynomial
   USE TypesDef, ONLY: dp
   USE OPC4D_Tensor_Module, ONLY: OPC4D_Tensor
   IMPLICIT NONE

   REAL(dp) :: vertices_initial(8,3), vertices_final(8,3)
   INTEGER :: facets(12,3)

   vertices_initial = RESHAPE([ &
      -1.0_dp,-1.0_dp,-1.0_dp,  1.0_dp,-1.0_dp,-1.0_dp,  1.0_dp,1.0_dp,-1.0_dp, -1.0_dp,1.0_dp,-1.0_dp, &
      -1.0_dp,-1.0_dp, 1.0_dp,  1.0_dp,-1.0_dp, 1.0_dp,  1.0_dp,1.0_dp, 1.0_dp, -1.0_dp,1.0_dp, 1.0_dp], [8,3], ORDER=[2,1])

   vertices_final = RESHAPE([ &
      -0.90_dp,-0.95_dp,-0.90_dp,  1.10_dp,-0.95_dp,-0.90_dp,  1.10_dp,1.05_dp,-0.90_dp, -0.90_dp,1.05_dp,-0.90_dp, &
      -0.90_dp,-0.95_dp, 1.10_dp,  1.10_dp,-0.95_dp, 1.10_dp,  1.10_dp,1.05_dp, 1.10_dp, -0.90_dp,1.05_dp, 1.10_dp], [8,3], ORDER=[2,1])

   facets = RESHAPE([ &
      1,3,2, 1,4,3, 5,6,7, 5,7,8, 1,2,6, 1,6,5, &
      4,8,7, 4,7,3, 1,5,8, 1,8,4, 2,3,7, 2,7,6], [12,3], ORDER=[2,1])

   WRITE(*,'(/,A)') '**************************************************************'
   WRITE(*,'(A)')   '                 OPTIMALPOLYCUBA4D Tensor'
   WRITE(*,'(A)')   '            Cubatura Tensor su Cubo in Movimento'
   WRITE(*,'(A,/)') '**************************************************************'

   CALL run_case(1, 'f_1(x,y,z,tau) = 1', 8.0_dp, 1)
   CALL run_case(2, 'f_2(x,y,z,tau) = x^2 + y^2 + z^2 + tau^2', 8.0_dp + 8.18_dp / 3.0_dp, 2)
   CALL run_case(8, 'f_3(x,y,z,tau) = x^2*y^2*z^2*tau^2', 29150263.0_dp / 283500000.0_dp, 3)

CONTAINS

   SUBROUTINE run_case(degree, label, exact, which_fun)
      INTEGER, INTENT(IN) :: degree, which_fun
      CHARACTER(LEN=*), INTENT(IN) :: label
      REAL(dp), INTENT(IN) :: exact
      REAL(dp), ALLOCATABLE :: XYZT(:,:), weights(:)
      INTEGER :: c0, c1, rate
      REAL(dp) :: elapsed, integral

      WRITE(*,'(/,A)') '--------------------------------------------------------------'
      WRITE(*,'(A,A)') 'Funzione integranda:   ', label
      WRITE(*,'(A,I0)') 'Ade:                   ', degree
      WRITE(*,'(A,ES22.15)') 'Integrale esatto:      ', exact

      CALL SYSTEM_CLOCK(c0, rate)
      CALL OPC4D_Tensor(degree, vertices_initial, vertices_final, facets, 'GJCC', XYZT, weights)
      CALL SYSTEM_CLOCK(c1, rate)
      elapsed = REAL(c1 - c0, dp) / REAL(rate, dp)

      SELECT CASE (which_fun)
      CASE (1)
         integral = SUM(weights)
      CASE (2)
         integral = SUM(weights * (XYZT(:,1)**2 + XYZT(:,2)**2 + XYZT(:,3)**2 + XYZT(:,4)**2))
      CASE (3)
         integral = SUM(weights * XYZT(:,1)**2 * XYZT(:,2)**2 * XYZT(:,3)**2 * XYZT(:,4)**2)
      END SELECT

      WRITE(*,'(A,I0)')       'Numero di nodi 4D:    ', SIZE(XYZT,1)
      WRITE(*,'(A,ES22.15)')  'Integrale numerico:   ', integral
      WRITE(*,'(A,ES12.6)')   'Errore assoluto:      ', ABS(integral - exact)
      WRITE(*,'(A,ES12.6,A)') 'Tempo di calcolo:     ', elapsed, ' s'
   END SUBROUTINE run_case
END PROGRAM example_polynomial
