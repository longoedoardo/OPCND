PROGRAM mainOptimalPolyCuba3D
    USE TypesDef
    USE PolyhedronMesh
    USE PrepCheap
    USE CubaCheap
    IMPLICIT NONE
 
    !**********************************************************************
    ! Variabili locali attive
    !**********************************************************************
    TYPE(t_polyhedron)          :: poly
    INTEGER                     :: k
 
    !**********************************************************************
    ! Setup iniziale
    !**********************************************************************

    ! Prima di procedere alla compilazione e all'esecuzione, verificare 
    ! e adattare i seguenti parametri in base alle specifiche esigenze:
    !
    !   - ade (Grado di Esattezza Algebrica):
    !      Rappresenta il grado massimo dei polinomi integrati in modo esatto.
    !      Per funzioni non polinomiali, l'incremento di ade riduce l'errore.
    !
    !   - Geometria (Mesh Files):
    !      I file indicati nella chiamata a MeshReader ('data/vertex.dat' 
    !      e 'data/tri.dat') definiscono i vertici e le facce del poliedro.
    !
    !   - Funzione Integranda (feval):
    !      Espressione analitica di f(x, y, z). Eventuali argomenti di input 
    !      non utilizzati genereranno un warning di compilazione del tipo 
    !      "Unused dummy argument", che non pregiudica l'esecuzione del codice.

    ade = 4

    WRITE(*,*) "Inizio caricamento geometria..."
    ! Legge vertices.dat e tri.dat, calcola bbox, salva tutto in poly
    CALL MeshReader('data/vertex.dat', 'data/tri.dat', poly)

    !**********************************************************************
    ! INIZIO PARTE SHAPE-INDEPENDENT
    !**********************************************************************

    WRITE(*,*) "Inizio parte shape-indipendent..."
    ! Griglia tensoriale Gauss-Chebyshev: input 2*ade, output XYZW_tens_ref
    CALL cub_gausscheb_tens3D(2*ade, XYZW_tens_ref)
    
    N = SIZE(XYZW_tens_ref, 1) ! N = numero di punti della griglia tensore

    ! Calcolo del numero di monomi N_mom per grado totale <= ade in 3D
    N_mom = ((ade + 1) * (ade + 2) * (ade + 3)) / 6

    ! Indici multi-indice (i,j,k) in ordine grlex dimensionati su N_mom
    IF (ALLOCATED(chebyshev_indices)) DEALLOCATE(chebyshev_indices)
    ALLOCATE(chebyshev_indices(N_mom, 3))
    
    ! Inizializza la prima riga a 0 (Monomio di grado 0: [0, 0, 0])
    chebyshev_indices(1,:) = 0
    
    DO k = 2, N_mom
        ind_curr = chebyshev_indices(k-1,:) 
        CALL mono_next_grlex(3, ind_curr)
        chebyshev_indices(k,:) = ind_curr
    END DO

    ! Coordinate dei punti (colonne 1:3 di XYZW_tens_ref)
    IF (ALLOCATED(X)) DEALLOCATE(X)
    ALLOCATE(X(N, 3))
    X(:,1) = XYZW_tens_ref(:,1)
    X(:,2) = XYZW_tens_ref(:,2)
    X(:,3) = XYZW_tens_ref(:,3)

    WRITE(*,*) "Fine parte shape-indipendent!"
    WRITE(*,*) "**********************************************************************"
 
    !**********************************************************************
    ! INIZIO PARTE SHAPE-DEPENDENT
    !**********************************************************************

    ! Matrice di Vandermonde di Chebyshev: V_ref(N, N_mom)
    CALL dCHEBVAND(ade, X, chebyshev_indices, V_ref)

    ! Coefficienti di norma^2 dei polinomi di Chebyshev
    CALL tenscheb_norm2sq(chebyshev_indices, coeffs)

    ! Momenti di Chebyshev sul poliedro
    CALL chebyshev_moments_polyhedron(poly%vertici, poly%facce, ade, chebyshev_indices, poly%bbox, moments_ch)

    ! Riscalamento della regola dal cubo di riferimento al bbox del poliedro
    CALL scale_rule(XYZW_tens_ref, poly%bbox, XYZW_tens)

    !**********************************************************************
    ! PESI FINALI E VALUTAZIONE
    !**********************************************************************

    IF (ALLOCATED(W)) DEALLOCATE(W)
    ALLOCATE(W(N))
    DO k = 1, N
        W(k) = XYZW_tens(k,4) * DOT_PRODUCT(V_ref(k,:), moments_ch / coeffs)
    END DO
 
    ! Regola finale XYZW = [XYZW_tens(:,1:3), W]
    IF (ALLOCATED(XYZW)) DEALLOCATE(XYZW)
    ALLOCATE(XYZW(N, 4))
    XYZW(:,1) = XYZW_tens(:,1)
    XYZW(:,2) = XYZW_tens(:,2)
    XYZW(:,3) = XYZW_tens(:,3)
    XYZW(:,4) = W
 
    ! Valutazione di f nei punti di cubatura e calcolo integrale
    IF (ALLOCATED(fXYZW)) DEALLOCATE(fXYZW)
    ALLOCATE(fXYZW(N))
    CALL feval(XYZW(:,1), XYZW(:,2), XYZW(:,3), N, fXYZW)
 
    ! I = W' * fXYZW
    I_risultato = DOT_PRODUCT(W, fXYZW)
    WRITE(*,*) "**********************************************************************"
    WRITE(*, '("I = ", ES22.15)') I_risultato
    WRITE(*,*) "**********************************************************************"
 
    CONTAINS
 
    !**********************************************************************
    ! Funzione integranda f(x,y,z)
    !**********************************************************************
    SUBROUTINE feval(x_in, y_in, z_in, n_in, fval_out)
        IMPLICIT NONE
        
        INTEGER, INTENT(IN)  :: n_in
        REAL(dp), INTENT(IN)  :: x_in(n_in), y_in(n_in), z_in(n_in)
        REAL(dp), INTENT(OUT) :: fval_out(n_in)
        
       fval_out = 1.0_dp
        
    END SUBROUTINE feval

END PROGRAM mainOptimalPolyCuba3D