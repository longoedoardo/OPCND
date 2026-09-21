function moments = chebyshev_moments_polyhedron(vertices, facets, ade, chebyshev_indices, dbox, method)

%**************************************************************************
%
% function moments = chebyshev_moments_polyhedron(vertices, facets, ade,
%                                           chebyshev_indices, dbox, method)
%
% Calcolo dei momenti di una base di polinomi di Chebyshev tensoriali:
%                    P(x,y,z) = T_i(x) * T_j(y) * T_k(z)
% con grado totale (i+j+k) <= ade, integrati sul volume di un poliedro.
%
% L'integrale di volume viene ricondotto a un integrale di superficie sulle
% facce del poliedro tramite il Teorema della Divergenza. Si utilizza un
% campo vettoriale V = (f, 0, 0), dove f è la primitiva rispetto a x del
% polinomio di Chebyshev considerato.
%
% La quadratura sulle facce triangolari viene costruita una sola volta
% sul triangolo di riferimento e successivamente trasformata sul triangolo
% fisico per ciascuna faccia.
% Sono disponibili i metodi di Gauss-Jacobi e Dunavant.
%
%**************************************************************************
%
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
% method:
%                    Metodo di integrazione sulle facce triangolari:
%                    * "D": Metodo simmetrico di Dunavant
%                    * "GJ": Metodo di Gauss-Jacobi
%
%**************************************************************************
%
% OUTPUT:
%
% moments:           Vettore colonna contenente i momenti calcolati per
%                    ogni tripla di indici in chebyshev_indices.
%
%**************************************************************************
%
% Riferimento bibliografico:
% [1] E.B. Chin, J.B. Lasserre, N. Sukumar: "Numerical integration of
% homogeneous function on convex and nonconvex polygons and polyhedra".
% Computational Mechanics, Vol. 56, No. 6, pp 967-981.
%
%**************************************************************************

num_indici = size(chebyshev_indices, 1);
n_facce = size(facets, 1);

% Pre-allocazione dei momenti
moments = zeros(num_indici, 1);

v1 = vertices(facets(:, 1), :);
v2 = vertices(facets(:, 2), :);
v3 = vertices(facets(:, 3), :);

% Costruzione della regola di quadratura sul triangolo di riferimento.
[nodes_ref, weights_ref] = TriangleQuadrature(method, ade);

for k = 1:n_facce

    % Vertici della faccia triangolare
    V = [v1(k,:); v2(k,:); v3(k,:)];

    % Trasformazione dei punti e dei pesi dal triangolo di riferimento
    % al triangolo fisico.
    [nodes_XYZW, WV_CUB] = ShiftingTriangleQuadrature( ...
        V, nodes_ref, weights_ref);

    % Calcolo della normale alla faccia
    A = V(2,:) - V(1,:);
    B = V(3,:) - V(1,:);

    % Prodotto vettoriale
    cp = cross(A, B);

    % Norma del prodotto vettoriale
    area2 = norm(cp);

    if area2 <= 1e-14
        % Faccia degenere, nessun contributo
        continue;
    end

    % Normale esterna unitaria
    norm_ext = cp / area2;

    % Calcolo dei momenti superficiali grezzi
    moms_facet = cubature_tens_chebyshev_facet_V( ...
        nodes_XYZW, WV_CUB, chebyshev_indices, dbox);

    % Componente x della normale esterna, come richiesto dal
    % teorema della divergenza applicato al campo V = (f, 0, 0)
    moments = moments + norm_ext(1) * moms_facet;

end

end