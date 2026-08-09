PROGRAM mainOptimalPolyCuba4D_MoL
    USE TypesDef
    USE PolyhedronMesh
    USE TimeDiscretization
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
    !      I file indicati nella chiamata a MeshReader definiscono i vertici 
    !      iniziali, vertici finali e le facce del poliedro.
    !
    !   - Funzione Integranda (feval):
    !      Espressione analitica di f(x, y, z, tau). Eventuali argomenti di input 
    !      non utilizzati genereranno un warning di compilazione del tipo 
    !      "Unused dummy argument", che non pregiudica l'esecuzione del codice.

    WRITE(*,*) "........................"
    WRITE(*,*) "Cubatura con OptimalPolyCuba4D (Metodo delle linee)"
    WRITE(*,*) "........................"

    ade = 4
    n_tau = 100

    ! Legge vertices.dat e tri.dat, calcola bbox, salva tutto in poly
    CALL MeshReader('data/vertex.dat', 'data/vertex_new.dat', 'data/tri.dat', poly)

    CALL ComputeClenshawCurtis(n_tau, tau_nodes, tau_weights)

    !**********************************************************************
    ! INIZIO PARTE SHAPE-INDEPENDENT
    !**********************************************************************

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

    ! Matrice di Vandermonde di Chebyshev: V_ref(N, N_mom)
    CALL dCHEBVAND(ade, X, chebyshev_indices, V_ref)

    ! Coefficienti di norma^2 dei polinomi di Chebyshev
    CALL tenscheb_norm2sq(chebyshev_indices, coeffs)
 
    !**********************************************************************
    ! INIZIO PARTE SHAPE-DEPENDENT
    !**********************************************************************

    I_4D = 0.0

    IF (ALLOCATED(vertices_tau)) DEALLOCATE(vertices_tau)
    ALLOCATE(vertices_tau(SIZE(poly%vertici, 1), 3))

    IF (ALLOCATED(tau_array)) DEALLOCATE(tau_array)
    ALLOCATE(tau_array(N))

    IF (ALLOCATED(W)) DEALLOCATE(W)
    ALLOCATE(W(N))

    IF (ALLOCATED(fXYZW)) DEALLOCATE(fXYZW)
    ALLOCATE(fXYZW(N))

    DO k = 1, n_tau

        current_tau = tau_nodes(k)

        ! Vertici del poliedro al tempo corrente (interpolazione lineare)
        vertices_tau = (1.0_dp - current_tau) * poly%vertici + current_tau * poly%vertici_new

        ! Bounding box dinamica al tempo corrente
        bbox_tau(1) = MINVAL(vertices_tau(:,1))
        bbox_tau(2) = MAXVAL(vertices_tau(:,1))
        bbox_tau(3) = MINVAL(vertices_tau(:,2))
        bbox_tau(4) = MAXVAL(vertices_tau(:,2))
        bbox_tau(5) = MINVAL(vertices_tau(:,3))
        bbox_tau(6) = MAXVAL(vertices_tau(:,3))

        ! Momenti della base sul poliedro (teorema della divergenza sulle facce)
        CALL chebyshev_moments_polyhedron(vertices_tau, poly%facce, ade, chebyshev_indices, bbox_tau, moments_ch)

        ! Riscalamento della griglia di riferimento sulla bbox reale corrente
        CALL scale_rule(XYZW_tens_ref, bbox_tau, XYZW_tens)

        ! Pesi di cubatura ottenuti ricombinando la regola di riferimento
        W = XYZW_tens(:,4) * MATMUL(V_ref, moments_ch / coeffs)

        ! Valutazione della funzione integranda al tempo corrente
        tau_array = current_tau
        CALL feval(XYZW_tens(:,1), XYZW_tens(:,2), XYZW_tens(:,3), tau_array, N, fXYZW)

        ! Accumulo dell'integrale 4D
        I_4D = I_4D + tau_weights(k) * SUM(W * fXYZW)

    END DO

    WRITE(*,*) "**********************************************************************"
    WRITE(*, '("I = ", ES22.15)') I_4D
    WRITE(*,*) "**********************************************************************"
 
    CONTAINS
 
    !**********************************************************************
    ! Funzione integranda f(x,y,z,tau)
    !**********************************************************************
    SUBROUTINE feval(x_in, y_in, z_in, tau_in, n_in, fval_out)
        IMPLICIT NONE
        
        INTEGER, INTENT(IN)  :: n_in
        REAL(dp), INTENT(IN)  :: x_in(n_in), y_in(n_in), z_in(n_in), tau_in(n_in)
        REAL(dp), INTENT(OUT) :: fval_out(n_in)
        
       fval_out = 1.0_dp
        
    END SUBROUTINE feval

END PROGRAM mainOptimalPolyCuba4D_MoL