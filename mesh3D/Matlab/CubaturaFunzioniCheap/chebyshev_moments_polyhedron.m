function moments=chebyshev_moments_polyhedron(vertices,facets,ade,chebyshev_indices, dbox)
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

% Ciclo sulle facce
for k = 1:n_facce

    % Costruzione vettoriale dei vertici della faccia corrente
    XV = [v1(k, 1); v2(k, 1); v3(k, 1)];
    YV = [v1(k, 2); v2(k, 2); v3(k, 2)];
    ZV = [v1(k, 3); v2(k, 3); v3(k, 3)];

    % Proiezione 3D -> 2D
    [xv0, yv0, zv0, R10, RFM] = maptopolygon2(XV, YV, ZV);

    % Generazione nodi e pesi di cubatura sul poligono 2D piano
    [xyw, ~, ~, ~, ~] = polygauss_2018(ade+1, xv0, yv0);

    % Preparazione coordinate 3D dei nodi nel sistema locale
    n_nodi = size(xyw, 1);
    xv_cub = xyw(:, 1);
    yv_cub = xyw(:, 2);
    zw_cub = zv0(1) * ones(n_nodi, 1);

    % Mapping inverso 2D -> 3D e calcolo normale esterna
    [XV_CUB, YV_CUB, ZV_CUB, norm_ext] = maptopolygon3(xv_cub, yv_cub, zw_cub, R10, RFM);

    % Raggruppamento nodi 3D
    nodes_XYZW = [XV_CUB, YV_CUB, ZV_CUB];
    WV_CUB = xyw(:, 3);

    % Calcolo dei momenti superficiali grezzi
    moms_facet = cubature_tens_chebyshev_facet_V(nodes_XYZW, WV_CUB, chebyshev_indices, dbox);

    chebyshev_moms(:, k) = norm_ext(1) * moms_facet;
end

% Somma finale dei contributi di tutte le facce
moments = sum(chebyshev_moms, 2);
end