MODULE PolyhedronMesh
    USE TypesDef
    IMPLICIT NONE

    ! Definizione della struttura dei dati
    TYPE :: t_polyhedron
        REAL(dp), ALLOCATABLE                       :: vertici(:,:)
        REAL(dp), ALLOCATABLE                       :: vertici_new(:,:)
        INTEGER, ALLOCATABLE                        :: facce(:,:)
    END TYPE t_polyhedron

    CONTAINS

    SUBROUTINE MeshReader(vertici_file, vertici_new_file, tri_file, poly)

    IMPLICIT NONE

    !***********************************************************************
    ! Argomenti
    !***********************************************************************
    CHARACTER(LEN=*), INTENT(IN)                    :: vertici_file
    CHARACTER(LEN=*), INTENT(IN)                    :: vertici_new_file
    CHARACTER(LEN=*), INTENT(IN)                    :: tri_file

    TYPE(t_polyhedron), INTENT(OUT) :: poly
    !***********************************************************************
    ! Variabili locali
    !***********************************************************************
    INTEGER :: u_vert
    INTEGER :: u_vert_new
    INTEGER :: u_tri
    INTEGER :: ierr

    INTEGER :: n_vertici
    INTEGER :: n_vertici_new
    INTEGER :: n_facce

    INTEGER :: i
    INTEGER :: j
    !***********************************************************************

    ! Inizio lettura configurazione iniziale
    OPEN(NEWUNIT=u_vert, FILE=vertici_file, STATUS='OLD', &
         ACTION='READ', IOSTAT=ierr)

    IF (ierr /= 0) THEN
        ERROR STOP 'ERRORE: impossibile aprire il file dei vertici iniziali.'
    END IF

    ! Primo valore: numero di vertici
    READ(u_vert, *, IOSTAT=ierr) n_vertici

    IF (ierr /= 0) THEN
        CLOSE(u_vert)
        ERROR STOP 'ERRORE: impossibile leggere il numero di vertici iniziali.'
    END IF

    IF (n_vertici <= 0) THEN
        CLOSE(u_vert)
        ERROR STOP 'ERRORE: numero di vertici iniziali non valido.'
    END IF

    ALLOCATE(poly%vertici(n_vertici, 3))

    ! Coordinate dei vertici
    READ(u_vert, *, IOSTAT=ierr) &
        ((poly%vertici(i,j), j=1,3), i=1,n_vertici)

    IF (ierr /= 0) THEN
        CLOSE(u_vert)
        ERROR STOP 'ERRORE: impossibile leggere i vertici iniziali.'
    END IF

    CLOSE(u_vert)


    ! Inizio lettura configurazione finale

    OPEN(NEWUNIT=u_vert_new, FILE=vertici_new_file, STATUS='OLD', &
         ACTION='READ', IOSTAT=ierr)

    IF (ierr /= 0) THEN
        ERROR STOP 'ERRORE: impossibile aprire il file dei vertici finali.'
    END IF

    ! Primo valore: numero di vertici
    READ(u_vert_new, *, IOSTAT=ierr) n_vertici_new

    IF (ierr /= 0) THEN
        CLOSE(u_vert_new)
        ERROR STOP 'ERRORE: impossibile leggere il numero di vertici finali.'
    END IF

    ALLOCATE(poly%vertici_new(n_vertici_new, 3))

    ! Coordinate dei vertici finali
    READ(u_vert_new, *, IOSTAT=ierr) &
        ((poly%vertici_new(i,j), j=1,3), i=1,n_vertici_new)

    IF (ierr /= 0) THEN
        CLOSE(u_vert_new)
        ERROR STOP 'ERRORE: impossibile leggere i vertici finali.'
    END IF

    CLOSE(u_vert_new)


    ! Lettura facce triangolari
    OPEN(NEWUNIT=u_tri, FILE=tri_file, STATUS='OLD', &
         ACTION='READ', IOSTAT=ierr)

    IF (ierr /= 0) THEN
        ERROR STOP 'ERRORE: impossibile aprire il file delle facce.'
    END IF

    READ(u_tri, *, IOSTAT=ierr) n_facce

    IF (ierr /= 0) THEN
        CLOSE(u_tri)
        ERROR STOP 'ERRORE: impossibile leggere il numero di facce.'
    END IF

    ALLOCATE(poly%facce(n_facce, 3))

    ! Connettività delle facce
    READ(u_tri, *, IOSTAT=ierr) &
        ((poly%facce(i,j), j=1,3), i=1,n_facce)

    IF (ierr /= 0) THEN
        CLOSE(u_tri)
        ERROR STOP 'ERRORE: impossibile leggere le facce triangolari.'
    END IF

    CLOSE(u_tri)

END SUBROUTINE MeshReader

END MODULE PolyhedronMesh