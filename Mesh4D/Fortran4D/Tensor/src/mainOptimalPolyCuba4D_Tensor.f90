PROGRAM mainOptimalPolyCuba4D_Tensor
    USE TypesDef
    USE PolyhedronMesh
    USE PrepCheap
    USE CubaCheap
    USE PrismQuadratureGJ
    
    IMPLICIT NONE
 
    !**********************************************************************
    ! Variabili locali
    !**********************************************************************
    TYPE(t_polyhedron_4D)       :: poly
    INTEGER                     :: k
 
    !**********************************************************************
    ! Setup iniziale
    !**********************************************************************

    WRITE(*,*) ".................................................."
    WRITE(*,*) "Cubatura con OptimalPolyCuba4D (Metodo Unificato)"
    WRITE(*,*) ".................................................."

    ade = 5

    ! Lettura corretta con i file vertex.dat (t=0), vertex_new.dat (t=1) e tri.dat
    CALL MeshReader('data/vertex.dat', 'data/vertex_new.dat', 'data/tri.dat', poly)

    !**********************************************************************
    ! INIZIO PARTE SHAPE-INDEPENDENT (4D)
    !**********************************************************************
    
    ! Griglia tensoriale Gauss-Chebyshev 4D: input ade, output XYZWT_tens_ref (Mx5)
    CALL cub_gausscheb_tens4D(2*ade, XYZWT_tens_ref)
    
    N = SIZE(XYZWT_tens_ref, 1) ! N = numero di punti della griglia tensore 4D

    ! Calcolo del numero di monomi N_mom per grado totale <= ade in 4D: (ade+4)! / (ade! * 4!)
    N_mom = ((ade + 1) * (ade + 2) * (ade + 3) * (ade + 4)) / 24

    ! Indici multi-indice (i,j,k,l) in ordine grlex dimensionati su N_mom x 4
    IF (ALLOCATED(chebyshev_indices)) DEALLOCATE(chebyshev_indices)
    ALLOCATE(chebyshev_indices(N_mom, 4))
    
    ! Inizializza la prima riga a 0 (Monomio di grado 0: [0, 0, 0, 0])
    chebyshev_indices(1,:) = 0
    
    DO k = 2, N_mom
        ind_curr = chebyshev_indices(k-1,:) 
        CALL mono_next_grlex(4, ind_curr)
        chebyshev_indices(k,:) = ind_curr
    END DO

    ! Coordinate dei punti (colonne 1:4 di XYZWT_tens_ref)
    IF (ALLOCATED(X)) DEALLOCATE(X)
    ALLOCATE(X(N, 4))
    X(:,1) = XYZWT_tens_ref(:,1)
    X(:,2) = XYZWT_tens_ref(:,2)
    X(:,3) = XYZWT_tens_ref(:,3)
    X(:,4) = XYZWT_tens_ref(:,4)
 
    !**********************************************************************
    ! INIZIO PARTE SHAPE-DEPENDENT
    !**********************************************************************

    ! Matrice di Vandermonde di Chebyshev 4D: V_ref(N, N_mom)
    CALL dCHEBVAND(ade, X, chebyshev_indices, V_ref)

    ! Coefficienti di norma^2 dei polinomi di Chebyshev in 4D
    CALL tenscheb_norm2sq(chebyshev_indices, coeffs)

    ! Momenti di Chebyshev sul poliedro 4D tramite teorema della divergenza
    CALL chebyshev_moments_polyhedron_4D(poly%vertici_4D, poly%Hyperfacets, ade, chebyshev_indices, poly%bbox_4D, moments_ch)

    ! Riscalamento della regola 4D dal cubo di riferimento al bbox 4D del poliedro
    CALL scale_rule(XYZWT_tens_ref, poly%bbox_4D, XYZWT_tens)

    !**********************************************************************
    ! PESI FINALI E VALUTAZIONE (4D)
    !**********************************************************************

    IF (ALLOCATED(W)) DEALLOCATE(W)
    ALLOCATE(W(N))
    DO k = 1, N
        W(k) = XYZWT_tens(k, 5) * DOT_PRODUCT(V_ref(k,:), moments_ch / coeffs)
    END DO
 
    ! Regola finale XYZWT = [XYZWT_tens(:,1:4), W]
    IF (ALLOCATED(XYZWT)) DEALLOCATE(XYZWT)
    ALLOCATE(XYZWT(N, 5))
    XYZWT(:,1) = XYZWT_tens(:,1)
    XYZWT(:,2) = XYZWT_tens(:,2)
    XYZWT(:,3) = XYZWT_tens(:,3)
    XYZWT(:,4) = XYZWT_tens(:,4)
    XYZWT(:,5) = W
 
    ! Valutazione di f nei punti di cubatura 4D e calcolo integrale
    IF (ALLOCATED(fXYZWT)) DEALLOCATE(fXYZWT)
    ALLOCATE(fXYZWT(N))
    CALL feval(XYZWT(:,1), XYZWT(:,2), XYZWT(:,3), XYZWT(:,4), N, fXYZWT)
 
    I_4D = SUM(XYZWT(:, 5) * fXYZWT)

    WRITE(*,*) ".................................................."
    WRITE(*, *) "I = ", I_4D
    WRITE(*,*) ".................................................."

    CONTAINS
 
    !**********************************************************************
    ! Funzione integranda f(x,y,z,tau)
    !**********************************************************************
    SUBROUTINE feval(x_in, y_in, z_in, tau_in, n_in, fval_out)
        IMPLICIT NONE
        
        INTEGER, INTENT(IN)   :: n_in
        REAL(dp), INTENT(IN)  :: x_in(n_in), y_in(n_in), z_in(n_in), tau_in(n_in)
        REAL(dp), INTENT(OUT) :: fval_out(n_in)
        
        fval_out = 1.0_dp
        
    END SUBROUTINE feval
 
END PROGRAM mainOptimalPolyCuba4D_Tensor