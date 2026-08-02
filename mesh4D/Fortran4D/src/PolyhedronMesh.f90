MODULE PolyhedronMesh
    USE TypesDef
    IMPLICIT NONE

    ! Struttura per gestire l'iperfaccia spaziotemporale 4D
    ! NB: la normale NON viene piu' precalcolata qui: per una mesh che si
    ! deforma tra t=0 e t=1 (vertici_1 /= vertici_0) la normale locale varia
    ! punto per punto sul prisma, quindi va calcolata dentro PrismQuad4D
    ! tramite il prodotto vettoriale generalizzato a 4D dei tre vettori
    ! tangenti (g_xi, g_eta, g_tau), non da una costante valutata a t=0.
    TYPE :: t_hyperfacet
        INTEGER, ALLOCATABLE                :: Vertices_ID(:)
    END TYPE t_hyperfacet

    ! Definizione della struttura dei dati aggiornata al 4D spazio-temporale
    TYPE :: t_polyhedron_4D
        REAL(dp), ALLOCATABLE               :: vertici_4D(:,:) ! Matrice complessiva (2*n_vertici x 4)
        TYPE(t_hyperfacet), ALLOCATABLE     :: Hyperfacets(:)  ! Array delle facce 4D (prismi spaziotemporali)
        REAL(dp)                            :: bbox_4D(2, 4)   ! Bounding box 4D (min/max per ciascuna delle 4 coordinate)
    END TYPE t_polyhedron_4D

CONTAINS

    SUBROUTINE MeshReader4D(vertici_file, vertici_new_file, tri_file, poly)
        CHARACTER(LEN=*), INTENT(IN)            :: vertici_file
        CHARACTER(LEN=*), INTENT(IN)            :: vertici_new_file
        CHARACTER(LEN=*), INTENT(IN)            :: tri_file
        TYPE(t_polyhedron_4D), INTENT(OUT)      :: poly

        INTEGER                                 :: u_v0, u_v1, u_tri, ierr, n_vertici, n_facce, i, j, k
        REAL(dp), ALLOCATABLE                   :: vertici_0(:,:), vertici_1(:,:)
        INTEGER, ALLOCATABLE                    :: facce(:,:)

        ! 1. Lettura vertici a t=0
        OPEN(NEWUNIT=u_v0, FILE=vertici_file, STATUS='OLD', ACTION='READ', IOSTAT=ierr)
        IF (ierr /= 0) STOP "ERRORE: Impossibile aprire il file dei vertici a t=0."
        READ(u_v0, *) n_vertici
        ALLOCATE(vertici_0(n_vertici, 3))
        READ(u_v0, *) ((vertici_0(i, j), j=1,3), i=1,n_vertici)
        CLOSE(u_v0)

        ! 2. Lettura vertici a t=1
        OPEN(NEWUNIT=u_v1, FILE=vertici_new_file, STATUS='OLD', ACTION='READ', IOSTAT=ierr)
        IF (ierr /= 0) STOP "ERRORE: Impossibile aprire il file dei vertici a t=1."
        READ(u_v1, *) j ! Legge il numero di vertici per controllo opzionale
        IF (j /= n_vertici) STOP "ERRORE: Numero di vertici a t=0 e t=1 non coincidente."
        ALLOCATE(vertici_1(n_vertici, 3))
        READ(u_v1, *) ((vertici_1(i, j), j=1,3), i=1,n_vertici)
        CLOSE(u_v1)

        ! 3. Lettura delle facce (triangoli spaziali di base)
        OPEN(NEWUNIT=u_tri, FILE=tri_file, STATUS='OLD', ACTION='READ', IOSTAT=ierr)
        IF (ierr /= 0) STOP "ERRORE: Impossibile aprire il file delle facce."
        READ(u_tri, *) n_facce
        ALLOCATE(facce(n_facce, 3))
        READ(u_tri, *) ((facce(i, j), j=1,3), i=1,n_facce)
        CLOSE(u_tri)

        ! 4. Assemblaggio dei vertici 4D [x, y, z, tau]
        ! I primi n_vertici sono a tau=0, i successivi n_vertici sono a tau=1
        ALLOCATE(poly%vertici_4D(2 * n_vertici, 4))
        DO i = 1, n_vertici
            poly%vertici_4D(i, 1:3) = vertici_0(i, :)
            poly%vertici_4D(i, 4)   = 0.0_dp

            poly%vertici_4D(n_vertici + i, 1:3) = vertici_1(i, :)
            poly%vertici_4D(n_vertici + i, 4)   = 1.0_dp
        END DO

        ! 5. Creazione delle iperfacce 4D (prismi spaziotemporali per ogni faccia triangolare)
        ! Ogni faccia triangolare genera un prisma spaziotemporale definito dai
        ! 6 vertici: 3 nodi a t=0 e i corrispondenti 3 nodi a t=1.
        ! La normale locale (necessaria al teorema della divergenza) viene
        ! calcolata direttamente in PrismQuad4D nei punti di quadratura, non qui.
        ALLOCATE(poly%Hyperfacets(n_facce))

        DO k = 1, n_facce
            ALLOCATE(poly%Hyperfacets(k)%Vertices_ID(6))
            ! I 6 nodi del prisma spaziotemporale associato al triangolo k
            poly%Hyperfacets(k)%Vertices_ID(1) = facce(k, 1)
            poly%Hyperfacets(k)%Vertices_ID(2) = facce(k, 2)
            poly%Hyperfacets(k)%Vertices_ID(3) = facce(k, 3)
            poly%Hyperfacets(k)%Vertices_ID(4) = facce(k, 1) + n_vertici
            poly%Hyperfacets(k)%Vertices_ID(5) = facce(k, 2) + n_vertici
            poly%Hyperfacets(k)%Vertices_ID(6) = facce(k, 3) + n_vertici
        END DO

        ! 6. Calcolo della bounding box 4D (min/max per ogni coordinata da 1 a 4)
        DO j = 1, 4
            poly%bbox_4D(1, j) = MINVAL(poly%vertici_4D(:, j))
            poly%bbox_4D(2, j) = MAXVAL(poly%vertici_4D(:, j))
        END DO

        DEALLOCATE(vertici_0, vertici_1, facce)

        WRITE(*,*) "Geometria 4D caricata e processata con successo!"
        WRITE(*,*) "**********************************************************************"
    END SUBROUTINE MeshReader4D

END MODULE PolyhedronMesh