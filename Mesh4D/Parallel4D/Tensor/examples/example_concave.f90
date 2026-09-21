PROGRAM example_concave

USE TypesDef, ONLY: dp
USE PolyhedronMesh, ONLY: t_polyhedron, MeshReader
USE OPC4D_Parallel_Tensor_Module, ONLY: OPC4D_Parallel_Tensor
USE OMP_LIB, ONLY: omp_get_wtime, omp_get_max_threads

IMPLICIT NONE

!***********************************************************************
!
!   Esempio:
!       Cubatura 4D su dominio poliedrale convesso in movimento
!
!   Descrizione:
!       Questo esempio dimostra l'utilizzo del metodo
!       OptimalPolyCuba4D Tensor parallelizzato con OpenMP per
!       l'integrazione di una funzione su un dominio spazio-temporale
!       (x, y, z, tau), ottenuto dalla deformazione lineare, per
!       tau in [0,1], di un poliedro convesso rappresentato da una
!       mesh superficiale triangolare.
!
!       La configurazione iniziale (tau = 0) e quella finale (tau = 1)
!       sono lette rispettivamente da concave_vertex.dat e
!       concave_vertex_new.dat, con la stessa connettivita' delle facce.
!
!***********************************************************************

!***********************************************************************
! Dichiarazione delle variabili
!***********************************************************************

INTEGER                         :: ade
INTEGER                         :: n_vertici
INTEGER                         :: n_facce
INTEGER                         :: n_threads
INTEGER                         :: i

REAL(dp)                        :: Integrale
REAL(dp)                        :: elapsedTime
REAL(dp)                        :: t_start
REAL(dp)                        :: t_end

REAL(dp), ALLOCATABLE           :: XYZT(:,:)
REAL(dp), ALLOCATABLE           :: W(:)

TYPE(t_polyhedron)              :: poly

CHARACTER(LEN=256)              :: vertici_file
CHARACTER(LEN=256)              :: vertici_new_file
CHARACTER(LEN=256)              :: tri_file

!***********************************************************************
! Parametri e caricamento mesh
!***********************************************************************

ade = 10

! Numero massimo di thread disponibili/utilizzabili.
n_threads = omp_get_max_threads()

WRITE(*,'(A)')
WRITE(*,'(A)') '**************************************************************'
WRITE(*,'(A)') '                  OPC4D Parallel Tensor'
WRITE(*,'(A)') '    Cubatura su Dominio Poliedrale Convesso in Movimento'
WRITE(*,'(A)') '**************************************************************'
WRITE(*,'(A)')

WRITE(*,'(A,I0)') 'Ade:                   ', ade
WRITE(*,'(A,I0)') 'Numero di thread:      ', n_threads

! Percorsi relativi alla cartella examples.
vertici_file     = 'concave_vertex.dat'
vertici_new_file = 'concave_vertex_new.dat'
tri_file         = 'concave_tri.dat'

! Lettura della mesh superficiale triangolare (configurazione iniziale
! e finale) e della connettivita' delle facce.
CALL MeshReader(vertici_file, vertici_new_file, tri_file, poly)

n_vertici = SIZE(poly%vertici,1)
n_facce   = SIZE(poly%facce,1)

WRITE(*,'(A,I0)') 'Numero di vertici:     ', n_vertici
WRITE(*,'(A,I0)') 'Numero di facce:       ', n_facce

!***********************************************************************
! Definizione funzione integranda
!***********************************************************************

WRITE(*,'(A)') 'Funzione integranda:   f(x,y,z,tau) = 1'

!***********************************************************************
! Inizio regola di cubatura
!***********************************************************************

WRITE(*,'(A)')
WRITE(*,'(A)') 'Inizio Cubatura...'

! Misura del tempo reale di esecuzione.
t_start = omp_get_wtime()

CALL OPC4D_Parallel_Tensor(ade, poly%vertici, poly%vertici_new, poly%facce, 'GJCC', XYZT, W)

t_end = omp_get_wtime()

elapsedTime = t_end - t_start

!***********************************************************************
! Calcolo dell'integrale
!***********************************************************************

Integrale = 0.0_dp

DO i = 1, SIZE(W)
    Integrale = Integrale + &
                W(i) * f(XYZT(i,1), XYZT(i,2), XYZT(i,3), XYZT(i,4))
END DO

!***********************************************************************
! Visualizzazione risultati
!***********************************************************************

WRITE(*,'(A)') 'Fine Cubatura...'

WRITE(*,'(A)')
WRITE(*,'(A,I0)')       'Numero di nodi 4D:      ', SIZE(W)
WRITE(*,'(A,F0.15)')    'Integrale numerico:     ', Integrale
WRITE(*,'(A,ES14.6,A)') 'Tempo di calcolo:       ', elapsedTime, ' s'

CONTAINS

!***********************************************************************
! Funzione integranda
!***********************************************************************

REAL(dp) FUNCTION f(x,y,z,tau)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: x
    REAL(dp), INTENT(IN) :: y
    REAL(dp), INTENT(IN) :: z
    REAL(dp), INTENT(IN) :: tau

    f = 1.0_dp

END FUNCTION f

END PROGRAM example_concave