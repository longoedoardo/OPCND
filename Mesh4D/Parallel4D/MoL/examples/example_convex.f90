PROGRAM example_convex

USE TypesDef, ONLY: dp
USE PolyhedronMesh, ONLY: t_polyhedron, MeshReader
USE OPC4D_Parallel_MoL_Module, ONLY: OPC4D_Parallel_MoL
USE OMP_LIB, ONLY: omp_get_wtime, omp_set_num_threads, omp_set_dynamic

IMPLICIT NONE

!***********************************************************************
!
!   Esempio:
!       Cubatura 4D su dominio poliedrale convesso in movimento.
!       Confronto fra versione seriale e parallela.
!
!   Descrizione:
!       Questo esempio dimostra l'utilizzo del metodo
!       OptimalPolyCuba4D MoL per l'integrazione di una funzione su
!       un dominio spazio-temporale (x, y, z, tau), ottenuto dalla
!       deformazione lineare, per tau in [0,1], di un poliedro convesso
!       rappresentato da una mesh superficiale triangolare.
!
!       La configurazione iniziale (tau = 0) e quella finale (tau = 1)
!       sono lette rispettivamente da convex_vertex.dat e
!       convex_vertex_new.dat, con la stessa connettivita' delle facce.
!
!***********************************************************************

!***********************************************************************
! Dichiarazione delle variabili
!***********************************************************************

INTEGER                         :: ade
INTEGER                         :: n_tau
INTEGER                         :: n_vertici, n_facce
REAL(dp)                        :: Integrale
INTEGER                         :: i

REAL(dp), ALLOCATABLE           :: XYZT(:,:)
REAL(dp), ALLOCATABLE           :: W(:)

REAL(dp)                        :: elapsedTime
REAL(dp)                        :: t_start, t_end

TYPE(t_polyhedron)              :: poly

CHARACTER(LEN=256)              :: vertici_file
CHARACTER(LEN=256)              :: vertici_new_file
CHARACTER(LEN=256)              :: tri_file

!***********************************************************************
! Parametri e caricamento mesh
!***********************************************************************

ade   = 10
n_tau = 100

CALL omp_set_dynamic(.FALSE.)

WRITE(*,'(A)')
WRITE(*,'(A)') '**************************************************************'
WRITE(*,'(A)') '                       OPTIMALPOLYCUBA4D'
WRITE(*,'(A)') '       Cubatura su Dominio Poliedrale Convesso in Movimento'
WRITE(*,'(A)') '**************************************************************'
WRITE(*,'(A)')

WRITE(*,'(A,I0)') 'Ade:                   ', ade
WRITE(*,'(A,I0)') 'Numero di nodi tau:    ', n_tau

! Percorsi relativi alla cartella examples.
vertici_file     = 'convex_vertex.dat'
vertici_new_file = 'convex_vertex_new.dat'
tri_file         = 'convex_tri.dat'

! Lettura della mesh superficiale triangolare nelle configurazioni
! iniziale e finale e della connettivita' delle facce.
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
! Esecuzione seriale
!***********************************************************************

WRITE(*,'(A)')
WRITE(*,'(A)') '**************************************************************'
WRITE(*,'(A)') '                     ESECUZIONE SERIALE'
WRITE(*,'(A)') '                         1 THREAD'
WRITE(*,'(A)') '**************************************************************'
WRITE(*,'(A)')

CALL omp_set_num_threads(1)

WRITE(*,'(A)') 'Inizio Cubatura...'

t_start = omp_get_wtime()

CALL OPC4D_Parallel_MoL( &
     ade, n_tau, &
     poly%vertici, poly%vertici_new, poly%facce, &
     'GJ', XYZT, W)

! CALL OPC4D_Parallel_MoL( &
!      ade, n_tau, &
!      poly%vertici, poly%vertici_new, poly%facce, &
!      'D', XYZT, W)

t_end = omp_get_wtime()

elapsedTime = t_end - t_start

Integrale = 0.0_dp

DO i = 1, SIZE(W)
    Integrale = Integrale + W(i) * &
                f(XYZT(i,1), XYZT(i,2), XYZT(i,3), XYZT(i,4))
END DO

WRITE(*,'(A)') 'Fine Cubatura...'

WRITE(*,'(A)')
WRITE(*,'(A,I0)')       'Numero di nodi 4D:      ', SIZE(W)
WRITE(*,'(A,F0.15)')    'Integrale numerico:     ', Integrale
WRITE(*,'(A,F12.6,A)')  'Tempo di calcolo:       ', elapsedTime, ' s'

!***********************************************************************
! Esecuzione parallela con 8 thread
!***********************************************************************

WRITE(*,'(A)')
WRITE(*,'(A)') '**************************************************************'
WRITE(*,'(A)') '                  ESECUZIONE PARALLELA'
WRITE(*,'(A)') '                         8 THREADS'
WRITE(*,'(A)') '**************************************************************'
WRITE(*,'(A)')

CALL omp_set_num_threads(8)

WRITE(*,'(A)') 'Inizio Cubatura...'

t_start = omp_get_wtime()

CALL OPC4D_Parallel_MoL( &
     ade, n_tau, &
     poly%vertici, poly%vertici_new, poly%facce, &
     'GJ', XYZT, W)

! CALL OPC4D_Parallel_MoL( &
!      ade, n_tau, &
!      poly%vertici, poly%vertici_new, poly%facce, &
!      'D', XYZT, W)

t_end = omp_get_wtime()

elapsedTime = t_end - t_start

Integrale = 0.0_dp

DO i = 1, SIZE(W)
    Integrale = Integrale + W(i) * &
                f(XYZT(i,1), XYZT(i,2), XYZT(i,3), XYZT(i,4))
END DO

WRITE(*,'(A)') 'Fine Cubatura...'

WRITE(*,'(A)')
WRITE(*,'(A,I0)')       'Numero di nodi 4D:      ', SIZE(W)
WRITE(*,'(A,F0.15)')    'Integrale numerico:     ', Integrale
WRITE(*,'(A,F12.6,A)')  'Tempo di calcolo:       ', elapsedTime, ' s'

CONTAINS

!***********************************************************************
! Funzione integranda
!
! Modificare esclusivamente questa funzione per cambiare
! l'integranda dell'esempio.
!***********************************************************************

REAL(dp) FUNCTION f(x,y,z,tau)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: x
    REAL(dp), INTENT(IN) :: y
    REAL(dp), INTENT(IN) :: z
    REAL(dp), INTENT(IN) :: tau

    f = 1.0_dp

END FUNCTION f

END PROGRAM example_convex