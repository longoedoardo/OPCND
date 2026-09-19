PROGRAM example_convex

USE TypesDef, ONLY: dp
USE PolyhedronMesh, ONLY: t_polyhedron, MeshReader
USE OPC3D_Parallel_Module, ONLY: OPC3D_Parallel
USE OMP_LIB, ONLY: omp_get_wtime, omp_set_num_threads, omp_set_dynamic

IMPLICIT NONE

!***********************************************************************
!
!   Esempio:
!       Cubatura su dominio poliedrale convesso
!       Confronto fra versione seriale e parallela.
!
!   Descrizione:
!       Questo esempio dimostra l'utilizzo del metodo
!       OptimalPolyCuba3D per l'integrazione di una funzione su un
!       dominio poliedrale convesso, rappresentato da una mesh
!       superficiale triangolare.
!
!***********************************************************************

!***********************************************************************
!   Dichiarazione delle variabili
!***********************************************************************

INTEGER                         :: ade
INTEGER                         :: n_vertici, n_facce
REAL(dp)                        :: Integrale
INTEGER                         :: i
REAL(dp), ALLOCATABLE           :: XYZ(:,:)
REAL(dp), ALLOCATABLE           :: W(:)
REAL(dp)                        :: elapsedTime
REAL(dp)                        :: t_start, t_end

TYPE(t_polyhedron)              :: poly

CHARACTER(LEN=256)              :: vertici_file
CHARACTER(LEN=256)              :: tri_file

!***********************************************************************
!   Parametri e caricamento mesh
!***********************************************************************

ade = 30

CALL omp_set_dynamic(.FALSE.)

WRITE(*,'(A)')
WRITE(*,'(A)') '**************************************************************'
WRITE(*,'(A)') '                 OPTIMALPOLYCUBA3D'
WRITE(*,'(A)') '        Cubatura su Dominio Poliedrale Convesso'
WRITE(*,'(A)') '**************************************************************'
WRITE(*,'(A)')

WRITE(*,'(A,I0)') 'Ade:                   ', ade

! Percorsi relativi alla cartella examples.
vertici_file = 'convex_vertex.dat'
tri_file     = 'convex_tri.dat'

! Lettura della mesh superficiale triangolare.
CALL MeshReader(vertici_file, tri_file, poly)

n_vertici = SIZE(poly%vertici,1)
n_facce   = SIZE(poly%facce,1)

WRITE(*,'(A,I0)') 'Numero di vertici:     ', n_vertici
WRITE(*,'(A,I0)') 'Numero di facce:       ', n_facce

!***********************************************************************
!   Definizione funzione integranda
!***********************************************************************

WRITE(*,'(A)') 'Funzione integranda:   f(x,y,z) = 1'

!***********************************************************************
!   Esecuzione seriale
!***********************************************************************

WRITE(*,'(A)')
WRITE(*,'(A)') '**************************************************************'
WRITE(*,'(A)') '                     ESECUZIONE SERIALE'
WRITE(*,'(A)') '**************************************************************'
WRITE(*,'(A)')

CALL omp_set_num_threads(1)

WRITE(*,'(A)') 'Inizio Cubatura...'

t_start = omp_get_wtime()

CALL OPC3D_Parallel(ade, poly%vertici, poly%facce, 'GJ', XYZ, W)
! CALL OPC3D_Parallel(ade, poly%vertici, poly%facce, 'D', XYZ, W)

t_end = omp_get_wtime()

elapsedTime = t_end - t_start

Integrale = 0.0_dp

DO i = 1, SIZE(W)
    Integrale = Integrale + W(i) * f(XYZ(i,1), XYZ(i,2), XYZ(i,3))
END DO

WRITE(*,'(A)') 'Fine Cubatura...'

WRITE(*,'(A)')
WRITE(*,'(A,F0.15)') 'Integrale numerico:     ', Integrale
WRITE(*,'(A,F12.6,A)') 'Tempo di calcolo:       ', elapsedTime, ' s'

!***********************************************************************
!   Esecuzione parallela con 8 thread
!***********************************************************************

WRITE(*,'(A)')
WRITE(*,'(A)') '**************************************************************'
WRITE(*,'(A)') '                  ESECUZIONE PARALLELA'
WRITE(*,'(A)') '                     8 THREADS'
WRITE(*,'(A)') '**************************************************************'
WRITE(*,'(A)')

CALL omp_set_num_threads(8)

WRITE(*,'(A)') 'Inizio Cubatura...'

t_start = omp_get_wtime()

CALL OPC3D_Parallel(ade, poly%vertici, poly%facce, 'GJ', XYZ, W)
! CALL OPC3D_Parallel(ade, poly%vertici, poly%facce, 'D', XYZ, W)

t_end = omp_get_wtime()

elapsedTime = t_end - t_start

Integrale = 0.0_dp

DO i = 1, SIZE(W)
    Integrale = Integrale + W(i) * f(XYZ(i,1), XYZ(i,2), XYZ(i,3))
END DO

WRITE(*,'(A)') 'Fine Cubatura...'

WRITE(*,'(A)')
WRITE(*,'(A,F0.15)') 'Integrale numerico:     ', Integrale
WRITE(*,'(A,F12.6,A)') 'Tempo di calcolo:       ', elapsedTime, ' s'

CONTAINS

!***********************************************************************
!   Funzione integranda
!
!   Modificare esclusivamente questa funzione per cambiare
!   l'integranda dell'esempio.
!***********************************************************************

REAL(dp) FUNCTION f(x,y,z)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: x
    REAL(dp), INTENT(IN) :: y
    REAL(dp), INTENT(IN) :: z

    f = 1.0_dp

END FUNCTION f

END PROGRAM example_convex