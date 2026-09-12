PROGRAM example_concave

    USE TypesDef, ONLY: dp
    USE PolyhedronMesh, ONLY: t_polyhedron, MeshReader
    USE OptimalPolyCuba3D_Module, ONLY: OptimalPolyCuba3D

    IMPLICIT NONE

    !***********************************************************************
    !
    !   Esempio:
    !       Cubatura su dominio poliedrale concavo
    !
    !   Descrizione:
    !       Questo esempio dimostra l'utilizzo del metodo
    !       OptimalPolyCuba3D per l'integrazione di una funzione su un
    !       dominio poliedrale concavo, rappresentato da una mesh
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

    ade = 1

    WRITE(*,'(A)')
    WRITE(*,'(A)') '**************************************************************'
    WRITE(*,'(A)') '                 OPTIMALPOLYCUBA3D'
    WRITE(*,'(A)') '        Cubatura su Dominio Poliedrale concavo'
    WRITE(*,'(A)') '**************************************************************'
    WRITE(*,'(A)')

    WRITE(*,'(A,I0)') 'Ade:                   ', ade

    ! Percorsi relativi alla cartella examples.
    vertici_file = 'concave_vertex.dat'
    tri_file     = 'concave_tri.dat'

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
    !   Inizio regola di cubatura
    !***********************************************************************

    WRITE(*,'(A)')
    WRITE(*,'(A)') 'Inizio Cubatura...'

    CALL CPU_TIME(t_start)

    CALL OptimalPolyCuba3D(ade, poly%vertici, poly%facce, XYZ, W)

    CALL CPU_TIME(t_end)

    elapsedTime = t_end - t_start

    Integrale = 0.0_dp
    DO i = 1, SIZE(W)
        Integrale = Integrale + W(i) * f(XYZ(i,1), XYZ(i,2), XYZ(i,3))
    END DO

    !***********************************************************************
    !   Visualizzazione risultati e punti di cubatura
    !***********************************************************************

    WRITE(*,'(A)') 'Fine Cubatura...'

    WRITE(*,'(A)')
    WRITE(*,'(A,F0.15)') 'Integrale numerico:       ', Integrale
    WRITE(*,'(A,ES14.6,A)') 'Tempo di calcolo:       ', elapsedTime, ' s'

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

END PROGRAM example_concave