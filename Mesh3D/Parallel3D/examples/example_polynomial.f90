PROGRAM example_polynomial

USE TypesDef, ONLY: dp
USE PolyhedronMesh, ONLY: t_polyhedron, MeshReader
USE OPC3D_Parallel_Module, ONLY: OPC3D_Parallel
USE OMP_LIB, ONLY: omp_get_wtime, omp_set_num_threads, omp_set_dynamic

IMPLICIT NONE

!***********************************************************************
!
!   Esempio:
!       Cubatura su dominio cubico.
!       Confronto fra versione seriale e parallela.
!
!   Descrizione:
!       Questo esempio dimostra l'utilizzo del metodo OptimalPolyCuba3D
!       per l'integrazione di funzioni polinomiali su un cubo, rappresentato
!       mediante una mesh superficiale triangolare chiusa e orientata
!       secondo le normali esterne.
!
!       Vengono considerate tre funzioni integrande di diverso grado
!       polinomiale e i risultati numerici vengono confrontati con
!       i corrispondenti valori analitici.
!
!***********************************************************************

!***********************************************************************
!   Dichiarazione delle variabili
!***********************************************************************

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

REAL(dp), ALLOCATABLE           :: vertices(:,:)
INTEGER, ALLOCATABLE            :: facets(:,:)

REAL(dp), ALLOCATABLE           :: XYZ(:,:)
REAL(dp), ALLOCATABLE           :: W(:)
REAL(dp), ALLOCATABLE           :: fXYZ(:)

!***********************************************************************
!   Parametri e caricamento mesh
!***********************************************************************

CALL omp_set_dynamic(.FALSE.)

WRITE(*,'(A)')
WRITE(*,'(A)') '**************************************************************'
WRITE(*,'(A)') '                 OPTIMALPOLYCUBA3D'
WRITE(*,'(A)') '                 Cubatura sul Cubo'
WRITE(*,'(A)') '**************************************************************'
WRITE(*,'(A)')

!***********************************************************************
!   Vertici e facce del cubo [-1,1]^3
!***********************************************************************

ALLOCATE(vertices(8,3))

vertices = TRANSPOSE(RESHAPE([ &
-1.0_dp, -1.0_dp, -1.0_dp, &
 1.0_dp, -1.0_dp, -1.0_dp, &
 1.0_dp,  1.0_dp, -1.0_dp, &
-1.0_dp,  1.0_dp, -1.0_dp, &
-1.0_dp, -1.0_dp,  1.0_dp, &
 1.0_dp, -1.0_dp,  1.0_dp, &
 1.0_dp,  1.0_dp,  1.0_dp, &
-1.0_dp,  1.0_dp,  1.0_dp  &
], [3,8]))

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

n_vertici = SIZE(vertices,1)
n_facce   = SIZE(facets,1)

WRITE(*,'(A,I0)') 'Numero di vertici:     ', n_vertici
WRITE(*,'(A,I0)') 'Numero di facce:       ', n_facce

!***********************************************************************
!   Funzione integranda f1
!***********************************************************************

ade = 20
I_exact = 8.0_dp

WRITE(*,'(A)')
WRITE(*,'(A)') '--------------------------------------------------------------'
WRITE(*,'(A)') 'Funzione integranda 1'
WRITE(*,'(A)') '--------------------------------------------------------------'
WRITE(*,'(A,I0)')       'Ade:                   ', ade
WRITE(*,'(A)')          'Funzione integranda:   f_1(x,y,z) = 1'
WRITE(*,'(A,ES22.15)')  'Integrale esatto:      ', I_exact

!***********************************************************************
!   Esecuzione seriale
!***********************************************************************

WRITE(*,'(A)')
WRITE(*,'(A)') '                     ESECUZIONE SERIALE'
WRITE(*,'(A)')

CALL omp_set_num_threads(1)

WRITE(*,'(A)') 'Inizio Cubatura...'

t_start = omp_get_wtime()

CALL OPC3D_Parallel(ade, vertices, facets, 'GJ', XYZ, W)
! CALL OPC3D_Parallel(ade, vertices, facets, 'D', XYZ, W)

t_end = omp_get_wtime()

elapsedTime = t_end - t_start

ALLOCATE(fXYZ(SIZE(W)))

DO i = 1, SIZE(W)
    fXYZ(i) = f1(XYZ(i,1), XYZ(i,2), XYZ(i,3))
END DO

Integrale = DOT_PRODUCT(W, fXYZ)

DEALLOCATE(fXYZ)

WRITE(*,'(A)') 'Fine Cubatura...'

error_abs = ABS(Integrale - I_exact)

WRITE(*,'(A,ES22.15)') 'Integrale numerico:    ', Integrale
WRITE(*,'(A,ES14.6)')  'Errore assoluto:       ', error_abs
WRITE(*,'(A,F12.6,A)') 'Tempo di calcolo:      ', elapsedTime, ' s'

!***********************************************************************
!   Esecuzione parallela con 8 thread
!***********************************************************************

WRITE(*,'(A)')
WRITE(*,'(A)') '                     ESECUZIONE PARALLELA'
WRITE(*,'(A)') '                         8 THREADS'
WRITE(*,'(A)')

CALL omp_set_num_threads(8)

WRITE(*,'(A)') 'Inizio Cubatura...'

t_start = omp_get_wtime()

CALL OPC3D_Parallel(ade, vertices, facets, 'GJ', XYZ, W)
! CALL OPC3D_Parallel(ade, vertices, facets, 'D', XYZ, W)

t_end = omp_get_wtime()

elapsedTime = t_end - t_start

ALLOCATE(fXYZ(SIZE(W)))

DO i = 1, SIZE(W)
    fXYZ(i) = f1(XYZ(i,1), XYZ(i,2), XYZ(i,3))
END DO

Integrale = DOT_PRODUCT(W, fXYZ)

DEALLOCATE(fXYZ)

WRITE(*,'(A)') 'Fine Cubatura...'

error_abs = ABS(Integrale - I_exact)

WRITE(*,'(A,ES22.15)') 'Integrale numerico:    ', Integrale
WRITE(*,'(A,ES14.6)')  'Errore assoluto:       ', error_abs
WRITE(*,'(A,F12.6,A)') 'Tempo di calcolo:      ', elapsedTime, ' s'

!***********************************************************************
!   Funzione integranda f2
!***********************************************************************

ade = 2
I_exact = 8.0_dp

WRITE(*,'(A)')
WRITE(*,'(A)') '--------------------------------------------------------------'
WRITE(*,'(A)') 'Funzione integranda 2'
WRITE(*,'(A)') '--------------------------------------------------------------'
WRITE(*,'(A,I0)')       'Ade:                   ', ade
WRITE(*,'(A)')          'Funzione integranda:   f_2(x,y,z) = x^2 + y^2 + z^2'
WRITE(*,'(A,ES22.15)')  'Integrale esatto:      ', I_exact

!***********************************************************************
!   Esecuzione seriale
!***********************************************************************

WRITE(*,'(A)')
WRITE(*,'(A)') '                     ESECUZIONE SERIALE'
WRITE(*,'(A)')

CALL omp_set_num_threads(1)

WRITE(*,'(A)') 'Inizio Cubatura...'

t_start = omp_get_wtime()

CALL OPC3D_Parallel(ade, vertices, facets, 'GJ', XYZ, W)
! CALL OPC3D_Parallel(ade, vertices, facets, 'D', XYZ, W)

t_end = omp_get_wtime()

elapsedTime = t_end - t_start

ALLOCATE(fXYZ(SIZE(W)))

DO i = 1, SIZE(W)
    fXYZ(i) = f2(XYZ(i,1), XYZ(i,2), XYZ(i,3))
END DO

Integrale = DOT_PRODUCT(W, fXYZ)

DEALLOCATE(fXYZ)

WRITE(*,'(A)') 'Fine Cubatura...'

error_abs = ABS(Integrale - I_exact)

WRITE(*,'(A,ES22.15)') 'Integrale numerico:    ', Integrale
WRITE(*,'(A,ES14.6)')  'Errore assoluto:       ', error_abs
WRITE(*,'(A,F12.6,A)') 'Tempo di calcolo:      ', elapsedTime, ' s'

!***********************************************************************
!   Esecuzione parallela con 8 thread
!***********************************************************************

WRITE(*,'(A)')
WRITE(*,'(A)') '                     ESECUZIONE PARALLELA'
WRITE(*,'(A)') '                         8 THREADS'
WRITE(*,'(A)')

CALL omp_set_num_threads(8)

WRITE(*,'(A)') 'Inizio Cubatura...'

t_start = omp_get_wtime()

CALL OPC3D_Parallel(ade, vertices, facets, 'GJ', XYZ, W)
! CALL OPC3D_Parallel(ade, vertices, facets, 'D', XYZ, W)

t_end = omp_get_wtime()

elapsedTime = t_end - t_start

ALLOCATE(fXYZ(SIZE(W)))

DO i = 1, SIZE(W)
    fXYZ(i) = f2(XYZ(i,1), XYZ(i,2), XYZ(i,3))
END DO

Integrale = DOT_PRODUCT(W, fXYZ)

DEALLOCATE(fXYZ)

WRITE(*,'(A)') 'Fine Cubatura...'

error_abs = ABS(Integrale - I_exact)

WRITE(*,'(A,ES22.15)') 'Integrale numerico:    ', Integrale
WRITE(*,'(A,ES14.6)')  'Errore assoluto:       ', error_abs
WRITE(*,'(A,F12.6,A)') 'Tempo di calcolo:      ', elapsedTime, ' s'

!***********************************************************************
!   Funzione integranda f3
!***********************************************************************

ade = 6
I_exact = 8.0_dp / 27.0_dp

WRITE(*,'(A)')
WRITE(*,'(A)') '--------------------------------------------------------------'
WRITE(*,'(A)') 'Funzione integranda 3'
WRITE(*,'(A)') '--------------------------------------------------------------'
WRITE(*,'(A,I0)')       'Ade:                   ', ade
WRITE(*,'(A)')          'Funzione integranda:   f_3(x,y,z) = x^2*y^2*z^2'
WRITE(*,'(A,ES22.15)')  'Integrale esatto:      ', I_exact

!***********************************************************************
!   Esecuzione seriale
!***********************************************************************

WRITE(*,'(A)')
WRITE(*,'(A)') '                     ESECUZIONE SERIALE'
WRITE(*,'(A)')

CALL omp_set_num_threads(1)

WRITE(*,'(A)') 'Inizio Cubatura...'

t_start = omp_get_wtime()

CALL OPC3D_Parallel(ade, vertices, facets, 'D', XYZ, W)

t_end = omp_get_wtime()

elapsedTime = t_end - t_start

ALLOCATE(fXYZ(SIZE(W)))

DO i = 1, SIZE(W)
    fXYZ(i) = f3(XYZ(i,1), XYZ(i,2), XYZ(i,3))
END DO

Integrale = DOT_PRODUCT(W, fXYZ)

DEALLOCATE(fXYZ)

WRITE(*,'(A)') 'Fine Cubatura...'

error_abs = ABS(Integrale - I_exact)

WRITE(*,'(A,ES22.15)') 'Integrale numerico:    ', Integrale
WRITE(*,'(A,ES14.6)')  'Errore assoluto:       ', error_abs
WRITE(*,'(A,F12.6,A)') 'Tempo di calcolo:      ', elapsedTime, ' s'

!***********************************************************************
!   Esecuzione parallela con 8 thread
!***********************************************************************

WRITE(*,'(A)')
WRITE(*,'(A)') '                     ESECUZIONE PARALLELA'
WRITE(*,'(A)') '                         8 THREADS'
WRITE(*,'(A)')

CALL omp_set_num_threads(8)

WRITE(*,'(A)') 'Inizio Cubatura...'

t_start = omp_get_wtime()

CALL OPC3D_Parallel(ade, vertices, facets, 'D', XYZ, W)

t_end = omp_get_wtime()

elapsedTime = t_end - t_start

ALLOCATE(fXYZ(SIZE(W)))

DO i = 1, SIZE(W)
    fXYZ(i) = f3(XYZ(i,1), XYZ(i,2), XYZ(i,3))
END DO

Integrale = DOT_PRODUCT(W, fXYZ)

DEALLOCATE(fXYZ)

WRITE(*,'(A)') 'Fine Cubatura...'

error_abs = ABS(Integrale - I_exact)

WRITE(*,'(A,ES22.15)') 'Integrale numerico:    ', Integrale
WRITE(*,'(A,ES14.6)')  'Errore assoluto:       ', error_abs
WRITE(*,'(A,F12.6,A)') 'Tempo di calcolo:      ', elapsedTime, ' s'

!***********************************************************************
!   Deallocazione memoria
!***********************************************************************

IF (ALLOCATED(vertices)) DEALLOCATE(vertices)
IF (ALLOCATED(facets))   DEALLOCATE(facets)
IF (ALLOCATED(XYZ))      DEALLOCATE(XYZ)
IF (ALLOCATED(W))        DEALLOCATE(W)

CONTAINS

!***********************************************************************
!   Funzione integranda f1
!***********************************************************************

REAL(dp) FUNCTION f1(x,y,z)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: x
    REAL(dp), INTENT(IN) :: y
    REAL(dp), INTENT(IN) :: z

    f1 = 1.0_dp

END FUNCTION f1

!***********************************************************************
!   Funzione integranda f2
!***********************************************************************

REAL(dp) FUNCTION f2(x,y,z)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: x
    REAL(dp), INTENT(IN) :: y
    REAL(dp), INTENT(IN) :: z

    f2 = x**2 + y**2 + z**2

END FUNCTION f2

!***********************************************************************
!   Funzione integranda f3
!***********************************************************************

REAL(dp) FUNCTION f3(x,y,z)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: x
    REAL(dp), INTENT(IN) :: y
    REAL(dp), INTENT(IN) :: z

    f3 = x**2 * y**2 * z**2

END FUNCTION f3

END PROGRAM example_polynomial