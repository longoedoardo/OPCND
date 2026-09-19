MODULE PolyhedronMesh
    USE TypesDef
IMPLICIT NONE

! Definizione della struttura

    TYPE :: t_hyperfacet
        INTEGER, ALLOCATABLE :: Vertices_ID(:)
    END TYPE t_hyperfacet
    
    TYPE :: t_polyhedron_4D
        REAL(dp), ALLOCATABLE           :: vertici_4D(:,:)   ! Vertici 4D [x,y,z,tau], dimensione (2*n_vertici) x 4
        TYPE(t_hyperfacet), ALLOCATABLE :: Hyperfacets(:)    ! Iperfacce: base, tetto, prismi laterali
        REAL(dp)                        :: bbox_4D(2,4)      ! bbox_4D(1,:) = min, bbox_4D(2,:) = max
    END TYPE t_polyhedron_4D

CONTAINS

SUBROUTINE MeshReader(vertici_file, vertici_new_file, tri_file, poly)
    CHARACTER(LEN=*), INTENT(IN)       :: vertici_file
    CHARACTER(LEN=*), INTENT(IN)       :: vertici_new_file
    CHARACTER(LEN=*), INTENT(IN)       :: tri_file
    TYPE(t_polyhedron_4D), INTENT(OUT) :: poly

    INTEGER :: u_vert, u_vert_new, u_tri
    INTEGER :: ierr
    INTEGER :: n_vertici, n_vertici_new, n_facce
    INTEGER :: i, j, m
    REAL(dp), ALLOCATABLE :: vertici_0(:,:), vertici_1(:,:)
    INTEGER,  ALLOCATABLE :: facce(:,:)

    ! Lettura dei vertici iniziali (tau = 0)
    OPEN(NEWUNIT=u_vert, FILE=vertici_file, STATUS='OLD', ACTION='READ', IOSTAT=ierr)
    IF (ierr /= 0) STOP "ERRORE: Impossibile aprire il file dei vertici iniziali."

    READ(u_vert, *, IOSTAT=ierr) n_vertici
    IF (ierr /= 0) STOP "ERRORE: Impossibile leggere il numero dei vertici iniziali."

    ALLOCATE(vertici_0(n_vertici, 3))

    READ(u_vert, *, IOSTAT=ierr) ((vertici_0(i, j), j = 1, 3), i = 1, n_vertici)
    IF (ierr /= 0) STOP "ERRORE: Impossibile leggere i vertici iniziali."

    CLOSE(u_vert)

    ! Lettura dei vertici finali (tau = 1)
    OPEN(NEWUNIT=u_vert_new, FILE=vertici_new_file, STATUS='OLD', ACTION='READ', IOSTAT=ierr)
    IF (ierr /= 0) STOP "ERRORE: Impossibile aprire il file dei vertici finali."

    READ(u_vert_new, *, IOSTAT=ierr) n_vertici_new
    IF (ierr /= 0) STOP "ERRORE: Impossibile leggere il numero dei vertici finali."

    IF (n_vertici_new /= n_vertici) THEN
    STOP "ERRORE: Il numero di vertici iniziali e finali non coincide."
    END IF

    ALLOCATE(vertici_1(n_vertici, 3))

    READ(u_vert_new, *, IOSTAT=ierr) ((vertici_1(i, j), j = 1, 3), i = 1, n_vertici)
    IF (ierr /= 0) STOP "ERRORE: Impossibile leggere i vertici finali."

    CLOSE(u_vert_new)

    ! Lettura delle facce triangolari
    OPEN(NEWUNIT=u_tri, FILE=tri_file, STATUS='OLD', ACTION='READ', IOSTAT=ierr)
    IF (ierr /= 0) STOP "ERRORE: Impossibile aprire il file delle facce."

    READ(u_tri, *, IOSTAT=ierr) n_facce
    IF (ierr /= 0) STOP "ERRORE: Impossibile leggere il numero delle facce."

    ALLOCATE(facce(n_facce, 3))

    READ(u_tri, *, IOSTAT=ierr) ((facce(i, j), j = 1, 3), i = 1, n_facce)
    IF (ierr /= 0) STOP "ERRORE: Impossibile leggere le facce."

    CLOSE(u_tri)

    !**********************************************************************
    ! Costruzione di vertici_4D: [x,y,z,tau], prima lo strato tau=0, poi tau=1
    !**********************************************************************
    ALLOCATE(poly%vertici_4D(2*n_vertici, 4))

    poly%vertici_4D(1:n_vertici, 1:3)            = vertici_0
    poly%vertici_4D(1:n_vertici, 4)              = 0.0_dp
    poly%vertici_4D(n_vertici+1:2*n_vertici, 1:3) = vertici_1
    poly%vertici_4D(n_vertici+1:2*n_vertici, 4)   = 1.0_dp

    !**********************************************************************
    ! Costruzione delle Hyperfacets: base (tau=0), tetto (tau=1), prismi laterali
    ! Le prime due non contribuiscono all'integrale (n_x = 0 per costruzione)
    ! e vengono saltate da chebyshev_moments_polyhedron_4D (ciclo da k=3).
    !**********************************************************************
    ALLOCATE(poly%Hyperfacets(n_facce + 2))

    ! Hyperfacets(1) = base a tau = 0
    ALLOCATE(poly%Hyperfacets(1)%Vertices_ID(n_vertici))
    poly%Hyperfacets(1)%Vertices_ID = [(i, i = 1, n_vertici)]

    ! Hyperfacets(2) = tetto a tau = 1
    ALLOCATE(poly%Hyperfacets(2)%Vertices_ID(n_vertici))
    poly%Hyperfacets(2)%Vertices_ID = [(i, i = n_vertici + 1, 2*n_vertici)]

    ! Hyperfacets(2+m) = prisma laterale associato alla faccia triangolare m,
    ! con vertici ordinati [A0,B0,C0,A1,B1,C1] come richiesto da PrismQuad4D
    DO m = 1, n_facce
    ALLOCATE(poly%Hyperfacets(2 + m)%Vertices_ID(6))
    poly%Hyperfacets(2 + m)%Vertices_ID(1) = facce(m, 1)
    poly%Hyperfacets(2 + m)%Vertices_ID(2) = facce(m, 2)
    poly%Hyperfacets(2 + m)%Vertices_ID(3) = facce(m, 3)
    poly%Hyperfacets(2 + m)%Vertices_ID(4) = facce(m, 1) + n_vertici
    poly%Hyperfacets(2 + m)%Vertices_ID(5) = facce(m, 2) + n_vertici
    poly%Hyperfacets(2 + m)%Vertices_ID(6) = facce(m, 3) + n_vertici
    END DO

    !**********************************************************************
    ! Bounding box 4D: inviluppo spaziale di entrambi gli strati, tempo in [0,1]
    !**********************************************************************
    poly%bbox_4D(1,1) = MIN(MINVAL(vertici_0(:,1)), MINVAL(vertici_1(:,1)))
    poly%bbox_4D(1,2) = MIN(MINVAL(vertici_0(:,2)), MINVAL(vertici_1(:,2)))
    poly%bbox_4D(1,3) = MIN(MINVAL(vertici_0(:,3)), MINVAL(vertici_1(:,3)))
    poly%bbox_4D(1,4) = 0.0_dp

    poly%bbox_4D(2,1) = MAX(MAXVAL(vertici_0(:,1)), MAXVAL(vertici_1(:,1)))
    poly%bbox_4D(2,2) = MAX(MAXVAL(vertici_0(:,2)), MAXVAL(vertici_1(:,2)))
    poly%bbox_4D(2,3) = MAX(MAXVAL(vertici_0(:,3)), MAXVAL(vertici_1(:,3)))
    poly%bbox_4D(2,4) = 1.0_dp

    DEALLOCATE(vertici_0, vertici_1, facce)

    END SUBROUTINE MeshReader

END MODULE PolyhedronMesh