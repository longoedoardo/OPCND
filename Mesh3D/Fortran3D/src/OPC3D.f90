MODULE OPC3D_Module

USE TypesDef

IMPLICIT NONE

CONTAINS

SUBROUTINE OPC3D(ade, vertices, facets, method, XYZ, W)

!***********************************************************************
!
! SUBROUTINE OPC3D(ade, vertices, facets, method, XYZ, W)
!
! Calcola un'approssimazione numerica di un integrale di volume
! tridimensionale su un dominio poliedrico rappresentato mediante
! una mesh superficiale triangolare chiusa e orientata.
!
! La regola di cubatura è costruita combinando una base di Chebyshev
! tensoriale shape-independent con i momenti shape-dependent del dominio,
! calcolati a partire dalla rappresentazione della superficie mediante
! il teorema della divergenza.
!
! Il metodo non richiede una decomposizione volumetrica del dominio
! mediante tetraedri e utilizza esclusivamente la triangolazione della
! superficie del poliedro.
!
! La costruzione della regola è suddivisa in due componenti:
!
!   1. PARTE SHAPE-INDEPENDENT
!      Costruzione della griglia tensoriale di Gauss-Chebyshev sul cubo
!      di riferimento [-1,1]^3 e della corrispondente matrice di
!      Vandermonde-Chebyshev.
!
!   2. PARTE SHAPE-DEPENDENT
!      Calcolo dei momenti della base di Chebyshev sul dominio fisico
!      e successivo riscalamento della regola sulla bounding box della
!      mesh.
!
!***********************************************************************
!
! INPUTS:
!
!   - ade:
!       Grado polinomiale totale massimo della regola di cubatura.
!       Il numero di momenti utilizzati è N_mom = (ade+1)(ade+2)(ade+3)/6.
!
!   - vertices:
!       Matrice N x 3 contenente le coordinate cartesiane dei vertici
!       della mesh superficiale. Ogni riga rappresenta un vertice
!       [x_i, y_i, z_i].
!
!   - facets:
!       Matrice M x 3 contenente la connettività della mesh superficiale
!       triangolare. Ogni riga contiene gli indici dei tre vertici che
!       definiscono una faccia triangolare. Le facce devono rappresentare 
!       una superficie chiusa e orientata coerentemente, in modo che le 
!       normali associate risultino compatibili con l'orientazione del 
!       bordo del dominio.
!
!   - method:
!       Metodo utilizzato per la quadratura delle facce triangolari.
!       'D' : quadratura di Dunavant. Disponibile fino al grado 20.
!       'G' : quadratura di Gauss-Jacobi sul triangolo di riferimento.
!
!***********************************************************************
!
! OUTPUTS:
!
!   - XYZ:
!       Matrice N x 3 contenente le coordinate cartesiane dei nodi
!       di cubatura nel dominio fisico. Ogni riga rappresenta un nodo
!       [x_k, y_k, z_k].
!
!   - W:
!       Vettore N x 1 contenente i pesi associati ai nodi di cubatura.
!
!***********************************************************************
!
!   Autore:
!       Edoardo Longo, Università degli Studi di Verona
!
!   Data:
!       Settembre 2026
!
!***********************************************************************

USE ReferenceFunctions
USE CubatureFunctions

IMPLICIT NONE

!**********************************************************************
! Argomenti
!**********************************************************************
INTEGER,                                    INTENT(IN)  :: ade
CHARACTER(LEN=*),                           INTENT(IN)  :: method
REAL(dp),                                   INTENT(IN)  :: vertices(:,:)
INTEGER,                                    INTENT(IN)  :: facets(:,:)
REAL(dp), ALLOCATABLE,                      INTENT(OUT) :: XYZ(:,:)
REAL(dp), ALLOCATABLE,                      INTENT(OUT) :: W(:)
!**********************************************************************
! Variabili locali
!**********************************************************************
INTEGER :: N, N_mom
INTEGER :: k
INTEGER :: ind_curr(3)
REAL(dp), ALLOCATABLE                                   :: XYZW_tens_ref(:,:)
REAL(dp), ALLOCATABLE                                   :: X(:,:)
REAL(dp), ALLOCATABLE                                   :: XYZW_tens(:,:)
REAL(dp), ALLOCATABLE                                   :: V_ref(:,:)
REAL(dp), ALLOCATABLE                                   :: coeffs(:)

INTEGER, ALLOCATABLE                                    :: chebyshev_indices(:,:)

REAL(dp), ALLOCATABLE                                   :: bbox(:)
REAL(dp), ALLOCATABLE                                   :: moments_ch(:)
REAL(dp), ALLOCATABLE                                   :: alpha(:)
!***********************************************************************

! Controllo della disponibilità delle regole di Dunavant
IF (TRIM(method) == 'D' .AND. ade > 19) THEN
    WRITE(*,'(A)') 'WARNING: Le regole di Dunavant sono disponibili solo fino al grado 20.'
    RETURN
END IF

!**********************************************************************
! INIZIO PARTE SHAPE-INDEPENDENT
!**********************************************************************

! Griglia tensoriale di Gauss-Chebyshev sul cubo di riferimento [-1,1]^3
CALL cub_gausscheb_tens3D(2*ade, XYZW_tens_ref)

! Dimensione dello spazio polinomiale di grado totale <= ade
N = SIZE(XYZW_tens_ref, 1)
N_mom = ((ade + 1) * (ade + 2) * (ade + 3)) / 6

! Indici (i,j,k) dei monomi tensoriali di Chebyshev, ordinamento GRLEX
IF (ALLOCATED(chebyshev_indices)) THEN
    DEALLOCATE(chebyshev_indices)
END IF
ALLOCATE(chebyshev_indices(N_mom,3))
chebyshev_indices(1,:) = 0
DO k = 2, N_mom
    ind_curr = chebyshev_indices(k-1,:)
    CALL mono_next_grlex(3, ind_curr)
    chebyshev_indices(k,:) = ind_curr
END DO

IF (ALLOCATED(X)) DEALLOCATE(X)
ALLOCATE(X(N,3))

X(:,1) = XYZW_tens_ref(:,1)
X(:,2) = XYZW_tens_ref(:,2)
X(:,3) = XYZW_tens_ref(:,3)

! Matrice di Vandermonde-Chebyshev 3D sui nodi della griglia di riferimento
CALL dCHEBVAND(ade, X, chebyshev_indices, V_ref)

! Coefficienti di normalizzazione della base di Chebyshev tensoriale
CALL tenscheb_norm2sq(chebyshev_indices, coeffs)

!**********************************************************************
! INIZIO PARTE SHAPE-DEPENDENT
!**********************************************************************

IF (ALLOCATED(bbox)) DEALLOCATE(bbox)
ALLOCATE(bbox(6))

! Bounding box: [xmin, xmax, ymin, ymax, zmin, zmax]
bbox = [MINVAL(vertices(:,1)), MAXVAL(vertices(:,1)), &
        MINVAL(vertices(:,2)), MAXVAL(vertices(:,2)), &
        MINVAL(vertices(:,3)), MAXVAL(vertices(:,3))]

! % Momenti della base sul poliedro (teorema della divergenza sulle facce)
CALL chebyshev_moments_polyhedron(vertices, facets, ade, chebyshev_indices, bbox, method, moments_ch)

! % Riscalamento della griglia di riferimento sulla bounding box reale
CALL scale_rule(XYZW_tens_ref, bbox, XYZW_tens)

! Pesi di cubatura ottenuti ricombinando la regola di riferimento
IF (ALLOCATED(W)) DEALLOCATE(W)
ALLOCATE(W(N))
ALLOCATE(alpha(N_mom))

alpha = moments_ch / coeffs
W = XYZW_tens(:,4) * MATMUL(V_ref, alpha)

IF (ALLOCATED(XYZ)) DEALLOCATE(XYZ)
ALLOCATE(XYZ(N,3))

XYZ(:,1) = XYZW_tens(:,1)
XYZ(:,2) = XYZW_tens(:,2)
XYZ(:,3) = XYZW_tens(:,3)

END SUBROUTINE OPC3D

END MODULE OPC3D_Module