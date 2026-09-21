PROGRAM example_convex

    USE TypesDef, ONLY: dp
    USE PolyhedronMesh, ONLY: t_polyhedron, MeshReader
    USE OPC4D_Tensor_Module, ONLY: OPC4D_Tensor

    IMPLICIT NONE

    !***********************************************************************
    !
    !   Esempio:
    !       Cubatura 4D su dominio poliedrale convesso in movimento
    !
    !   Descrizione:
    !       Questo esempio dimostra l'utilizzo del metodo
    !       OptimalPolyCuba4D Tensor per l'integrazione di una funzione su
    !       un dominio spazio-temporale (x, y, z, tau), ottenuto dalla
    !       deformazione lineare, per tau in [0,1], di un poliedro convesso
    !       rappresentato da una mesh superficiale triangolare.
    !       La configurazione iniziale (tau = 0) e quella finale (tau = 1)
    !       sono lette rispettivamente da convex_vertex.dat e
    !       convex_vertex_new.dat, con la stessa connettivita' delle facce.
    !
    !***********************************************************************

    !***********************************************************************
    ! Dichiarazione delle variabili
    !***********************************************************************

    INTEGER                         :: ade
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

    ade = 10

    WRITE(*,'(A)')
    WRITE(*,'(A)') '**************************************************************'
    WRITE(*,'(A)') '                       OPC4D Tensor'
    WRITE(*,'(A)') '    Cubatura su Dominio Poliedrale Convesso in Movimento'
    WRITE(*,'(A)') '**************************************************************'
    WRITE(*,'(A)')

    WRITE(*,'(A,I0)') 'Ade:                   ', ade

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

    CALL CPU_TIME(t_start)

    CALL OPC4D_Tensor(ade, poly%vertici, poly%vertici_new, poly%facce, 'GJCC', XYZT, W)
    ! CALL OPC4D_Tensor(ade, poly%vertici, poly%vertici_new, poly%facce, 'GJL', XYZT, W)
    ! CALL OPC4D_Tensor(ade, poly%vertici, poly%vertici_new, poly%facce, 'DCC', XYZT, W)
    ! CALL OPC4D_Tensor(ade, poly%vertici, poly%vertici_new, poly%facce, 'DGL', XYZT, W)

    CALL CPU_TIME(t_end)

    elapsedTime = t_end - t_start

    Integrale = 0.0_dp
    DO i = 1, SIZE(W)
        Integrale = Integrale + W(i) * f(XYZT(i,1), XYZT(i,2), XYZT(i,3), XYZT(i,4))
    END DO

    !***********************************************************************
    ! Visualizzazione risultati e punti di cubatura
    !***********************************************************************

    WRITE(*,'(A)') 'Fine Cubatura...'

    WRITE(*,'(A)')
    WRITE(*,'(A,I0)') 'Numero di nodi 4D:      ', SIZE(W)
    WRITE(*,'(A,F0.15)') 'Integrale numerico:       ', Integrale
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

END PROGRAM example_convex