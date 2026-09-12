function moments = chebyshev_moments_polyhedron(vertices, facets, ade, chebyshev_indices, dbox)
%**************************************************************************
%
% function moments = chebyshev_moments_polyhedron(vertices, facets, ade, chebyshev_indices, dbox)
%
% Calcolo dei momenti di una base di polinomi di Chebyshev tensoriali:
%                         P(x,y,z) = T_i(x) * T_j(y) * T_k(z)
% con grado totale (i+j+k) <= ade, integrati sul volume di un poliedro.
%
% L'integrale di volume viene ricondotto a un integrale di superficie sulle
% facce del poliedro tramite il Teorema della Divergenza. Si utilizza un
% campo vettoriale V = (f, 0, 0), dove f è la primitiva rispetto a x del
% polinomio di Chebyshev considerato.
%
% L'integrazione su ciascuna faccia triangolare avviene direttamente in 3D
% tramite quadratura di Gauss-Jacobi sul triangolo (TriangleQuadraturePoints
% + mapTrianglePoints). 
%**************************************************************************
% INPUT:
%
% vertices:          Matrice [m x 3] delle coordinate dei vertici.
%
% facets:            Matrice [n_facce x 3] contenente gli indici dei vertici
%                    per ogni faccia triangolare (orientamento antiorario
%                    rispetto alla normale esterna).
%
% ade:               Grado algebrico massimo dei polinomi.
%
% chebyshev_indices: Matrice degli indici (i, j, k) che definisce l'ordine
%                    dei momenti da calcolare.
%
% dbox:              Iper-rettangolo [min_x, max_x, min_y, max_y, min_z, max_z]
%                    che racchiude il dominio per la base di Chebyshev.
%
%**************************************************************************
% OUTPUT:
%
% moments:           Vettore colonna contenente i momenti calcolati per
%                    ogni tripla di indici in chebyshev_indices.
%
%**************************************************************************
% Riferimento bibliografico:
% [1] E.B. Chin, J.B. Lasserre, N. Sukumar: "Numerical integration of
% homogeneous function on convex and nonconvex polygons and polyhedra".
% Computational Mechanics, Vol. 56, No. 6, pp 967-981.
%**************************************************************************

num_indici = size(chebyshev_indices, 1);
n_facce = size(facets, 1);

% Pre-allocazione della matrice dei momenti
chebyshev_moms = zeros(num_indici, n_facce);

v1 = vertices(facets(:, 1), :);
v2 = vertices(facets(:, 2), :);
v3 = vertices(facets(:, 3), :);

% Regola di quadratura sul triangolo di riferimento: dipende solo dal
% grado richiesto, quindi va calcolata una sola volta e riusata per tutte
% le facce.
% Grado di precisione richiesto: ade+1,
% quindi servono almeno nGP punti per dimensione tali che 2*nGP-1 >=
% ade+1.

nGP = ceil((ade + 2) / 2);
[nodi_rif, pesi_rif] = TriangleQuadraturePoints(nGP);

for k = 1:n_facce

    Vk = [v1(k,:); v2(k,:); v3(k,:)];

    % Nodi fisici, pesi fisici e normale esterna sulla faccia triangolare
    [nodes_XYZW, WV_CUB, norm_ext] = mapTrianglePoints(nodi_rif, pesi_rif, Vk);

    if all(norm_ext == 0)
        % Faccia degenere: nessun contributo
        continue;
    end

    % Calcolo dei momenti superficiali grezzi (integrale della primitiva
    % rispetto a x, valutata sui nodi della faccia)
    moms_facet = cubature_tens_chebyshev_facet_V(nodes_XYZW, WV_CUB, chebyshev_indices, dbox);

    % Componente x della normale esterna, come richiesto dal teorema della
    % divergenza applicato al campo V = (f, 0, 0)
    chebyshev_moms(:, k) = norm_ext(1) * moms_facet;
end

% Somma finale dei contributi di tutte le facce
moments = sum(chebyshev_moms, 2);

end