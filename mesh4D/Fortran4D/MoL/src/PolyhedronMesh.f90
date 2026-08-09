MODULE PolyhedronMesh
    USE TypesDef
    IMPLICIT NONE

    ! Definizione della struttura
    TYPE :: t_polyhedron
        REAL, ALLOCATABLE    :: vertici(:,:)      ! Vertici a tau = 0, dimensione n_vertici x 3
        REAL, ALLOCATABLE    :: vertici_new(:,:)  ! Vertici a tau = 1, dimensione n_vertici x 3
        INTEGER, ALLOCATABLE :: facce(:,:)        ! Facce triangolari, dimensione n_facce x 3
    END TYPE t_polyhedron

CONTAINS

    SUBROUTINE MeshReader(vertici_file, vertici_new_file, tri_file, poly)
        CHARACTER(LEN=*), INTENT(IN)    :: vertici_file
        CHARACTER(LEN=*), INTENT(IN)    :: vertici_new_file
        CHARACTER(LEN=*), INTENT(IN)    :: tri_file
        TYPE(t_polyhedron), INTENT(OUT) :: poly

        INTEGER :: u_vert, u_vert_new, u_tri
        INTEGER :: ierr
        INTEGER :: n_vertici, n_vertici_new, n_facce
        INTEGER :: i, j

        ! Lettura dei vertici iniziali
        OPEN(NEWUNIT=u_vert, FILE=vertici_file, STATUS='OLD', ACTION='READ', IOSTAT=ierr)
        IF (ierr /= 0) STOP "ERRORE: Impossibile aprire il file dei vertici iniziali."

        READ(u_vert, *, IOSTAT=ierr) n_vertici
        IF (ierr /= 0) STOP "ERRORE: Impossibile leggere il numero dei vertici iniziali."

        ALLOCATE(poly%vertici(n_vertici, 3))

        READ(u_vert, *, IOSTAT=ierr) ((poly%vertici(i, j), j = 1, 3), i = 1, n_vertici)
        IF (ierr /= 0) STOP "ERRORE: Impossibile leggere i vertici iniziali."

        CLOSE(u_vert)

        ! Lettura dei vertici finali
        OPEN(NEWUNIT=u_vert_new, FILE=vertici_new_file, STATUS='OLD', ACTION='READ', IOSTAT=ierr)
        IF (ierr /= 0) STOP "ERRORE: Impossibile aprire il file dei vertici finali."

        READ(u_vert_new, *, IOSTAT=ierr) n_vertici_new
        IF (ierr /= 0) STOP "ERRORE: Impossibile leggere il numero dei vertici finali."

        IF (n_vertici_new /= n_vertici) THEN
            STOP "ERRORE: Il numero di vertici iniziali e finali non coincide."
        END IF

        ALLOCATE(poly%vertici_new(n_vertici, 3))

        READ(u_vert_new, *, IOSTAT=ierr) ((poly%vertici_new(i, j), j = 1, 3), i = 1, n_vertici)
        IF (ierr /= 0) STOP "ERRORE: Impossibile leggere i vertici finali."

        CLOSE(u_vert_new)

        ! Lettura delle facce triangolari
        OPEN(NEWUNIT=u_tri, FILE=tri_file, STATUS='OLD', ACTION='READ', IOSTAT=ierr)
        IF (ierr /= 0) STOP "ERRORE: Impossibile aprire il file delle facce."

        READ(u_tri, *, IOSTAT=ierr) n_facce
        IF (ierr /= 0) STOP "ERRORE: Impossibile leggere il numero delle facce."

        ALLOCATE(poly%facce(n_facce, 3))

        READ(u_tri, *, IOSTAT=ierr) ((poly%facce(i, j), j = 1, 3), i = 1, n_facce)
        IF (ierr /= 0) STOP "ERRORE: Impossibile leggere le facce."

        CLOSE(u_tri)

    END SUBROUTINE MeshReader

END MODULE PolyhedronMesh