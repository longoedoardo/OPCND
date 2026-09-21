PROGRAM example_bunny

USE TypesDef, ONLY: dp
USE PolyhedronMesh, ONLY: t_polyhedron, MeshReader
USE OPC3D_Module, ONLY: OPC3D

    IMPLICIT NONE

    !***********************************************************************
    !
    !   Esempio:
    !       Cubatura su dominio poliedrale generico
    !
    !   Descrizione:
    !       Questo esempio dimostra l'utilizzo del metodo
    !       OptimalPolyCuba3D per l'integrazione di una funzione su un
    !       dominio poliedrale generico, rappresentato da una mesh
    !       superficiale triangolare.
    !
    !***********************************************************************

    !***********************************************************************
    ! Dichiarazione delle variabili
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
    ! Parametri e caricamento mesh
    !***********************************************************************

    ade = 2

    WRITE(*,'(A)')
    WRITE(*,'(A)') '**************************************************************'
    WRITE(*,'(A)') '                            OPC3D'
    WRITE(*,'(A)') '        Cubatura su Dominio Poliedrale "Stanford Bunny"'
    WRITE(*,'(A)') '**************************************************************'
    WRITE(*,'(A)')

    WRITE(*,'(A,I0)') 'Ade:                   ', ade

    ! Percorsi relativi alla cartella examples.
    vertici_file = 'bunny_vertex.dat'
    tri_file     = 'bunny_tri.dat'

    ! Lettura della mesh superficiale triangolare.
    CALL MeshReader(vertici_file, tri_file, poly)

    n_vertici = SIZE(poly%vertici,1)
    n_facce   = SIZE(poly%facce,1)

    WRITE(*,'(A,I0)') 'Numero di vertici:     ', n_vertici
    WRITE(*,'(A,I0)') 'Numero di facce:       ', n_facce

    !***********************************************************************
    !   Definizione funzione integranda
    !***********************************************************************

    WRITE(*,'(A)') 'Funzione integranda:   f(x,y,z) = x^2+y^2+z^2'

    !***********************************************************************
    !   Inizio regola di cubatura
    !***********************************************************************

    WRITE(*,'(A)')
    WRITE(*,'(A)') 'Inizio Cubatura...'

    CALL CPU_TIME(t_start)

    CALL OPC3D(ade, poly%vertici, poly%facce, 'D', XYZ, W)
    ! CALL OPC3D(ade, poly%vertici, poly%facce, 'D', XYZ, W)

    CALL CPU_TIME(t_end)

    elapsedTime = t_end - t_start

    Integrale = 0.0_dp
    DO i = 1, SIZE(W)
        Integrale = Integrale + W(i) * f(XYZ(i,1), XYZ(i,2), XYZ(i,3))
    END DO

    !***********************************************************************
    ! Visualizzazione risultati e punti di cubatura
    !***********************************************************************

    WRITE(*,'(A)') 'Fine Cubatura...'

    WRITE(*,'(A)')
    WRITE(*,'(A,F0.15)') 'Integrale numerico:       ', Integrale
    WRITE(*,'(A,ES14.6,A)') 'Tempo di calcolo:       ', elapsedTime, ' s'

    CONTAINS

    !***********************************************************************
    ! Funzione integranda
    !***********************************************************************

    REAL(dp) FUNCTION f(x,y,z)

        IMPLICIT NONE

        REAL(dp), INTENT(IN) :: x
        REAL(dp), INTENT(IN) :: y
        REAL(dp), INTENT(IN) :: z

        f = x**2+y**2+z**2

    END FUNCTION f

END PROGRAM example_bunny