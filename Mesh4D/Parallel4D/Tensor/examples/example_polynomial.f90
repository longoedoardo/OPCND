PROGRAM example_polynomial

USE TypesDef, ONLY: dp
USE OMP_LIB, ONLY: omp_get_wtime, omp_set_num_threads, omp_set_dynamic
USE OPC4D_Parallel_Tensor_Module, ONLY: OPC4D_Parallel_Tensor

IMPLICIT NONE

!*******************************************************************************
!
!   Esempio:
!       Cubatura 4D su dominio cubico in movimento
!
!   Descrizione:
!       Questo esempio dimostra l'utilizzo del metodo
!       OptimalPolyCuba4D Tensor per l'integrazione di funzioni
!       polinomiali su un dominio spazio-temporale (x, y, z, tau),
!       ottenuto dalla deformazione lineare, per tau in [0,1], di un
!       cubo rappresentato mediante una mesh superficiale triangolare
!       chiusa e orientata secondo le normali esterne.
!
!       Il cubo iniziale [-1,1]^3 viene traslato e allungato fino alla
!       configurazione finale [-0.90,1.10] x [-0.95,1.05] x [-0.90,1.10].
!
!       Vengono considerate tre funzioni integrande di diverso grado
!       polinomiale e i risultati numerici vengono confrontati con
!       i corrispondenti valori analitici.
!
!       Il calcolo viene eseguito con 1 e 8 thread OpenMP per confrontare
!       il tempo di esecuzione della versione parallelizzata sulle facce.
!
!*******************************************************************************

!*******************************************************************************
! Dichiarazione delle variabili
!*******************************************************************************

INTEGER                         :: ade
INTEGER                         :: n_vertici
INTEGER                         :: n_facce
INTEGER                         :: i

REAL(dp)                        :: Integrale
REAL(dp)                        :: I_exact
REAL(dp)                        :: error_abs
REAL(dp)                        :: elapsedTime
REAL(dp)                        :: t_start
REAL(dp)                        :: t_end

REAL(dp), ALLOCATABLE           :: vertices_initial(:,:)
REAL(dp), ALLOCATABLE           :: vertices_final(:,:)
INTEGER, ALLOCATABLE            :: facets(:,:)

REAL(dp), ALLOCATABLE           :: XYZT(:,:)
REAL(dp), ALLOCATABLE           :: W(:)
REAL(dp), ALLOCATABLE           :: fXYZT(:)

!*******************************************************************************
! Parametri e caricamento mesh
!*******************************************************************************

WRITE(*,'(A)')
WRITE(*,'(A)') '**************************************************************'
WRITE(*,'(A)') '                       OPC4D Tensor'
WRITE(*,'(A)') '            Cubatura sul Cubo in Movimento'
WRITE(*,'(A)') '**************************************************************'
WRITE(*,'(A)')

!*******************************************************************************
! Vertici e facce del cubo in movimento
!*******************************************************************************

! Configurazione iniziale (tau = 0): cubo [-1,1]^3
ALLOCATE(vertices_initial(8,3))

vertices_initial = TRANSPOSE(RESHAPE([ &
   -1.0_dp, -1.0_dp, -1.0_dp, &
    1.0_dp, -1.0_dp, -1.0_dp, &
    1.0_dp,  1.0_dp, -1.0_dp, &
   -1.0_dp,  1.0_dp, -1.0_dp, &
   -1.0_dp, -1.0_dp,  1.0_dp, &
    1.0_dp, -1.0_dp,  1.0_dp, &
    1.0_dp,  1.0_dp,  1.0_dp, &
   -1.0_dp,  1.0_dp,  1.0_dp  &
], [3,8]))

! Configurazione finale (tau = 1): stessi vertici, stesso ordinamento
ALLOCATE(vertices_final(8,3))

vertices_final = TRANSPOSE(RESHAPE([ &
   -0.90_dp, -0.95_dp, -0.90_dp, &
    1.10_dp, -0.95_dp, -0.90_dp, &
    1.10_dp,  1.05_dp, -0.90_dp, &
   -0.90_dp,  1.05_dp, -0.90_dp, &
   -0.90_dp, -0.95_dp,  1.10_dp, &
    1.10_dp, -0.95_dp,  1.10_dp, &
    1.10_dp,  1.05_dp,  1.10_dp, &
   -0.90_dp,  1.05_dp,  1.10_dp  &
], [3,8]))

! Connettivita' delle facce (invariata nel tempo), normali esterne
ALLOCATE(facets(12,3))

facets = TRANSPOSE(RESHAPE([ &
   1, 3, 2, &
   1, 4, 3, &
   5, 6, 7, &
   5, 7, 8, &
   1, 2, 6, &
   1, 6, 5, &
   4, 8, 7, &
   4, 7, 3, &
   1, 5, 8, &
   1, 8, 4, &
   2, 3, 7, &
   2, 7, 6  &
], [3,12]))

n_vertici = SIZE(vertices_initial,1)
n_facce   = SIZE(facets,1)

WRITE(*,'(A,I0)') 'Numero di vertici:     ', n_vertici
WRITE(*,'(A,I0)') 'Numero di facce:       ', n_facce

!*******************************************************************************
! Impostazioni OpenMP
!*******************************************************************************

CALL omp_set_dynamic(.FALSE.)

!*******************************************************************************
! Definizione funzione integranda f1
!*******************************************************************************

ade = 1
I_exact = 8.0_dp

WRITE(*,'(A)')
WRITE(*,'(A)') '--------------------------------------------------------------'
WRITE(*,'(A)') 'Funzione integranda 1'
WRITE(*,'(A)') '--------------------------------------------------------------'
WRITE(*,'(A,I0)')    'Ade:                   ', ade
WRITE(*,'(A)')       'Funzione integranda:   f_1(x,y,z,tau) = 1'
WRITE(*,'(A,ES22.15)') 'Integrale esatto:      ', I_exact

!*******************************************************************************
! Calcolo con 1 thread
!*******************************************************************************

WRITE(*,'(A)')
WRITE(*,'(A)') 'Calcolo con 1 thread...'

CALL omp_set_num_threads(1)

t_start = omp_get_wtime()

CALL OPC4D_Parallel_Tensor(ade, poly%vertici, poly%vertici_new, poly%facce, 'GJCC', XYZT, W)

t_end = omp_get_wtime()

elapsedTime = t_end - t_start

ALLOCATE(fXYZT(SIZE(W)))

DO i = 1, SIZE(W)
   fXYZT(i) = f1(XYZT(i,1), XYZT(i,2), XYZT(i,3), XYZT(i,4))
END DO

Integrale = DOT_PRODUCT(W, fXYZT)
error_abs = ABS(Integrale - I_exact)

DEALLOCATE(fXYZT)

WRITE(*,'(A)')
WRITE(*,'(A,I0)')        'Thread:                ', 1
WRITE(*,'(A,I0)')        'Numero di nodi 4D:    ', SIZE(W)
WRITE(*,'(A,ES22.15)')   'Integrale numerico:    ', Integrale
WRITE(*,'(A,ES14.6)')    'Errore assoluto:       ', error_abs
WRITE(*,'(A,ES14.6,A)')  'Tempo di calcolo:      ', elapsedTime, ' s'

!*******************************************************************************
! Calcolo con 8 thread
!*******************************************************************************

WRITE(*,'(A)')
WRITE(*,'(A)') 'Calcolo con 8 thread...'

CALL omp_set_num_threads(8)

t_start = omp_get_wtime()

CALL OPC4D_Parallel_Tensor(ade, poly%vertici, poly%vertici_new, poly%facce, 'GJCC', XYZT, W)

t_end = omp_get_wtime()

elapsedTime = t_end - t_start

IF (ALLOCATED(fXYZT)) DEALLOCATE(fXYZT)
ALLOCATE(fXYZT(SIZE(W)))

DO i = 1, SIZE(W)
   fXYZT(i) = f1(XYZT(i,1), XYZT(i,2), XYZT(i,3), XYZT(i,4))
END DO

Integrale = DOT_PRODUCT(W, fXYZT)
error_abs = ABS(Integrale - I_exact)

DEALLOCATE(fXYZT)

WRITE(*,'(A)')
WRITE(*,'(A,I0)')        'Thread:                ', 8
WRITE(*,'(A,I0)')        'Numero di nodi 4D:    ', SIZE(W)
WRITE(*,'(A,ES22.15)')   'Integrale numerico:    ', Integrale
WRITE(*,'(A,ES14.6)')    'Errore assoluto:       ', error_abs
WRITE(*,'(A,ES14.6,A)')  'Tempo di calcolo:      ', elapsedTime, ' s'

!*******************************************************************************
! Definizione funzione integranda f2
!*******************************************************************************

ade = 2
I_exact = 8.0_dp + 8.18_dp / 3.0_dp

WRITE(*,'(A)')
WRITE(*,'(A)') '--------------------------------------------------------------'
WRITE(*,'(A)') 'Funzione integranda 2'
WRITE(*,'(A)') '--------------------------------------------------------------'
WRITE(*,'(A,I0)')    'Ade:                   ', ade
WRITE(*,'(A)')       'Funzione integranda:   f_2(x,y,z,tau) = x^2 + y^2 + z^2 + tau^2'
WRITE(*,'(A,ES22.15)') 'Integrale esatto:      ', I_exact

!*******************************************************************************
! Calcolo con 1 thread
!*******************************************************************************

WRITE(*,'(A)')
WRITE(*,'(A)') 'Calcolo con 1 thread...'

CALL omp_set_num_threads(1)

t_start = omp_get_wtime()

CALL OPC4D_Parallel_Tensor(ade, poly%vertici, poly%vertici_new, poly%facce, 'GJCC', XYZT, W)

t_end = omp_get_wtime()

elapsedTime = t_end - t_start

IF (ALLOCATED(fXYZT)) DEALLOCATE(fXYZT)
ALLOCATE(fXYZT(SIZE(W)))

DO i = 1, SIZE(W)
   fXYZT(i) = f2(XYZT(i,1), XYZT(i,2), XYZT(i,3), XYZT(i,4))
END DO

Integrale = DOT_PRODUCT(W, fXYZT)
error_abs = ABS(Integrale - I_exact)

DEALLOCATE(fXYZT)

WRITE(*,'(A)')
WRITE(*,'(A,I0)')        'Thread:                ', 1
WRITE(*,'(A,I0)')        'Numero di nodi 4D:    ', SIZE(W)
WRITE(*,'(A,ES22.15)')   'Integrale numerico:    ', Integrale
WRITE(*,'(A,ES14.6)')    'Errore assoluto:       ', error_abs
WRITE(*,'(A,ES14.6,A)')  'Tempo di calcolo:      ', elapsedTime, ' s'

!*******************************************************************************
! Calcolo con 8 thread
!*******************************************************************************

WRITE(*,'(A)')
WRITE(*,'(A)') 'Calcolo con 8 thread...'

CALL omp_set_num_threads(8)

t_start = omp_get_wtime()

CALL OPC4D_Parallel_Tensor(ade, poly%vertici, poly%vertici_new, poly%facce, 'GJCC', XYZT, W)

t_end = omp_get_wtime()

elapsedTime = t_end - t_start

IF (ALLOCATED(fXYZT)) DEALLOCATE(fXYZT)
ALLOCATE(fXYZT(SIZE(W)))

DO i = 1, SIZE(W)
   fXYZT(i) = f2(XYZT(i,1), XYZT(i,2), XYZT(i,3), XYZT(i,4))
END DO

Integrale = DOT_PRODUCT(W, fXYZT)
error_abs = ABS(Integrale - I_exact)

DEALLOCATE(fXYZT)

WRITE(*,'(A)')
WRITE(*,'(A,I0)')        'Thread:                ', 8
WRITE(*,'(A,I0)')        'Numero di nodi 4D:    ', SIZE(W)
WRITE(*,'(A,ES22.15)')   'Integrale numerico:    ', Integrale
WRITE(*,'(A,ES14.6)')    'Errore assoluto:       ', error_abs
WRITE(*,'(A,ES22.15)')   'Integrale numerico:    ', Integrale
WRITE(*,'(A,ES14.6)')    'Errore assoluto:       ', error_abs
WRITE(*,'(A,ES14.6,A)')  'Tempo di calcolo:      ', elapsedTime, ' s'

!*******************************************************************************
! Definizione funzione integranda f3
!*******************************************************************************

! Il grado totale di f3 e' 8: serve ade = 8 per avere una regola esatta
ade = 8

I_exact = 29150263.0_dp / 283500000.0_dp

WRITE(*,'(A)')
WRITE(*,'(A)') '--------------------------------------------------------------'
WRITE(*,'(A)') 'Funzione integranda 3'
WRITE(*,'(A)') '--------------------------------------------------------------'
WRITE(*,'(A,I0)')    'Ade:                   ', ade
WRITE(*,'(A)')       'Funzione integranda:   f_3(x,y,z,tau) = x^2*y^2*z^2*tau^2'
WRITE(*,'(A,ES22.15)') 'Integrale esatto:      ', I_exact

!*******************************************************************************
! Calcolo con 1 thread
!*******************************************************************************

WRITE(*,'(A)')
WRITE(*,'(A)') 'Calcolo con 1 thread...'

CALL omp_set_num_threads(1)

t_start = omp_get_wtime()

CALL OPC4D_Parallel_Tensor(ade, poly%vertici, poly%vertici_new, poly%facce, 'GJCC', XYZT, W)

t_end = omp_get_wtime()

elapsedTime = t_end - t_start

IF (ALLOCATED(fXYZT)) DEALLOCATE(fXYZT)
ALLOCATE(fXYZT(SIZE(W)))

DO i = 1, SIZE(W)
   fXYZT(i) = f3(XYZT(i,1), XYZT(i,2), XYZT(i,3), XYZT(i,4))
END DO

Integrale = DOT_PRODUCT(W, fXYZT)
error_abs = ABS(Integrale - I_exact)

DEALLOCATE(fXYZT)

WRITE(*,'(A)')
WRITE(*,'(A,I0)')        'Thread:                ', 1
WRITE(*,'(A,I0)')        'Numero di nodi 4D:    ', SIZE(W)
WRITE(*,'(A,ES22.15)')   'Integrale numerico:    ', Integrale
WRITE(*,'(A,ES14.6)')    'Errore assoluto:       ', error_abs
WRITE(*,'(A,ES14.6,A)')  'Tempo di calcolo:      ', elapsedTime, ' s'

!*******************************************************************************
! Calcolo con 8 thread
!*******************************************************************************

WRITE(*,'(A)')
WRITE(*,'(A)') 'Calcolo con 8 thread...'

CALL omp_set_num_threads(8)

t_start = omp_get_wtime()

CALL OPC4D_Parallel_Tensor(ade, poly%vertici, poly%vertici_new, poly%facce, 'GJCC', XYZT, W)

t_end = omp_get_wtime()

elapsedTime = t_end - t_start

IF (ALLOCATED(fXYZT)) DEALLOCATE(fXYZT)
ALLOCATE(fXYZT(SIZE(W)))

DO i = 1, SIZE(W)
   fXYZT(i) = f3(XYZT(i,1), XYZT(i,2), XYZT(i,3), XYZT(i,4))
END DO

Integrale = DOT_PRODUCT(W, fXYZT)
error_abs = ABS(Integrale - I_exact)

DEALLOCATE(fXYZT)

WRITE(*,'(A)')
WRITE(*,'(A,I0)')        'Thread:                ', 8
WRITE(*,'(A,I0)')        'Numero di nodi 4D:    ', SIZE(W)
WRITE(*,'(A,ES22.15)')   'Integrale numerico:    ', Integrale
WRITE(*,'(A,ES14.6)')    'Errore assoluto:       ', error_abs
WRITE(*,'(A,ES14.6,A)')  'Tempo di calcolo:      ', elapsedTime, ' s'

!*******************************************************************************
! Deallocazione memoria
!*******************************************************************************

IF (ALLOCATED(vertices_initial)) DEALLOCATE(vertices_initial)
IF (ALLOCATED(vertices_final))   DEALLOCATE(vertices_final)
IF (ALLOCATED(facets))            DEALLOCATE(facets)
IF (ALLOCATED(XYZT))              DEALLOCATE(XYZT)
IF (ALLOCATED(W))                DEALLOCATE(W)

CONTAINS

!*******************************************************************************
! Funzione integranda f1
!*******************************************************************************

REAL(dp) FUNCTION f1(x,y,z,tau)

  IMPLICIT NONE

  REAL(dp), INTENT(IN) :: x
  REAL(dp), INTENT(IN) :: y
  REAL(dp), INTENT(IN) :: z
  REAL(dp), INTENT(IN) :: tau

  f1 = 1.0_dp

END FUNCTION f1

!*******************************************************************************
! Funzione integranda f2
!*******************************************************************************

REAL(dp) FUNCTION f2(x,y,z,tau)

  IMPLICIT NONE

  REAL(dp), INTENT(IN) :: x
  REAL(dp), INTENT(IN) :: y
  REAL(dp), INTENT(IN) :: z
  REAL(dp), INTENT(IN) :: tau

  f2 = x**2 + y**2 + z**2 + tau**2

END FUNCTION f2

!*******************************************************************************
! Funzione integranda f3
!*******************************************************************************

REAL(dp) FUNCTION f3(x,y,z,tau)

  IMPLICIT NONE

  REAL(dp), INTENT(IN) :: x
  REAL(dp), INTENT(IN) :: y
  REAL(dp), INTENT(IN) :: z
  REAL(dp), INTENT(IN) :: tau

  f3 = x**2 * y**2 * z**2 * tau**2

END FUNCTION f3

END PROGRAM example_polynomial