function chebyshev_moms = cubature_tens_chebyshev_facet_V(nodes, weights, chebyshev_indices, dbox)
%**************************************************************************
%
% moments = chebyshev_moments_polyhedron(vertices, facets, ade, chebyshev_indices, dbox)
%
% Calcola gli integrali di volume sul poliedro definito da 'vertices' e 'facets'
% per una base di polinomi tensoriali di Chebyshev P(x,y,z) = T_i(x)*T_j(y)*T_k(z),
% sfruttando il Teorema della Divergenza per convertire l'integrale di volume
% in un integrale di superficie sulle facce.
%
%**************************************************************************
%
% INPUT:
%   vertices            Matrice [m x 3] delle coordinate dei vertici.
%   facets              Matrice [n x 3] degli indici dei vertici per ogni faccia.
%   ade                 Grado algebrico massimo (totale o di riferimento).
%   chebyshev_indices   Matrice degli indici (i, j, k) della base tensoriale.
%   dbox                Iper-rettangolo [min_x, max_x, ..., max_z] del dominio.
%
% OUTPUT:
%   moments             Vettore colonna con i momenti integrati per ciascuna tripla.
%
%**************************************************************************

% Estrazione coordinate normalizzate e fattori di scala
inv_dx = 2.0 / (dbox(2) - dbox(1));
inv_dy = 2.0 / (dbox(4) - dbox(3));
inv_dz = 2.0 / (dbox(6) - dbox(5));

mid_x  = (dbox(1) + dbox(2)) * 0.5;
mid_y  = (dbox(3) + dbox(4)) * 0.5;
mid_z  = (dbox(5) + dbox(6)) * 0.5;

XN = (nodes(:, 1) - mid_x) * inv_dx;
YN = (nodes(:, 2) - mid_y) * inv_dy;
ZN = (nodes(:, 3) - mid_z) * inv_dz;

B1 = 1.0 / inv_dx; % Fattore di scala originale per x

% Grado massimo per asse
max_i = max(chebyshev_indices(:, 1));
max_j = max(chebyshev_indices(:, 2));
max_k = max(chebyshev_indices(:, 3));

% Calcolo polinomi di Chebyshev solo fino al grado necessario per asse
TX = chebpolys(max_i + 1, XN);
TY = chebpolys(max_j, YN);
TZ = chebpolys(max_k, ZN);

% Pre-calcolo della matrice delle primitive rispetto a X (IntX)
n = size(TX, 1);
IntX = zeros(n, max_i + 1);
IntX(:, 1) = XN;                  % i = 0

if max_i >= 1
    IntX(:, 2) = 0.5 * (XN .* XN); % i = 1
end

if max_i >= 2
    i_vec = 2:max_i;
    denom1 = i_vec + 1;
    denom2 = i_vec - 1;
    % Formula analitica vettorializzata
    IntX(:, i_vec + 1) = (i_vec .* TX(:, i_vec + 2) ./ (i_vec.^2 - 1)) - ...
        (XN .* TX(:, i_vec + 1) ./ denom2);
end

idx_i = chebyshev_indices(:, 1) + 1;
idx_j = chebyshev_indices(:, 2) + 1;
idx_k = chebyshev_indices(:, 3) + 1;

% Estrazione colonne selettive
IntX_cols = IntX(:, idx_i);
TY_cols   = TY(:,   idx_j);
TZ_cols   = TZ(:,   idx_k);

% Combinazione finale con i pesi tramite moltiplicazione matrice-vettore
w = weights(:);
chebyshev_moms = B1 * (((IntX_cols .* TY_cols .* TZ_cols)' * w));

end