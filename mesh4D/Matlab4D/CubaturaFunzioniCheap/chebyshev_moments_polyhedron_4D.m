function moments = chebyshev_moments_polyhedron_4D(vertici_4D, Hyperfacets, ade, chebyshev_indices, dbox)
%**************************************************************************
% Calcolo dei momenti di Chebyshev su un iperpoliedro 4D (Versione Punto-per-Punto)
%**************************************************************************
num_indici = size(chebyshev_indices, 1); 
n_facce = length(Hyperfacets);
chebyshev_moms = zeros(num_indici, n_facce); 

for k = 1:n_facce
    % Recupero dei 6 vertici 4D del prisma spazio-temporale k
    nodes_id = Hyperfacets{k}.Vertices_ID;
    V = vertici_4D(nodes_id, :); % Matrice 6x4 [v1; v2; v3; v1_new; v2_new; v3_new]

    % Mappatura sul prisma 4D tramite PrismQuad4D (restituisce punti e pesi locali col Gramiano)
    [XI_ref, ETA_ref, T_ref, W_ref] = reference_prism_quadrature(ade, dbox);
    [XYZTW, WV_CUB] = PrismQuad4D(V, XI_ref, ETA_ref, T_ref, W_ref);

    % Calcola i momenti superficiali per la iperfaccia 4D (i pesi WV_CUB includono già la metrica locale)
    moms_facet_raw = cubature_tens_chebyshev_facet_4D(XYZTW, WV_CUB, ade, chebyshev_indices, dbox);

    % Accumulo dei momenti per la faccia corrente (senza più moltiplicare per nx esterno costante)
    chebyshev_moms(:, k) = moms_facet_raw;
end

moments = sum(chebyshev_moms, 2);
end