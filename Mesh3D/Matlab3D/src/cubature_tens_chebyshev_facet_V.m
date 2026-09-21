function chebyshev_moms = cubature_tens_chebyshev_facet_V(nodes, weights, chebyshev_indices, dbox)
%**************************************************************************
%
% moments = chebyshev_moments_polyhedron(vertices, facets, ade, chebyshev_indices, dbox)
%
% Calcola gli integrali di volume sul poliedro per una base di polinomi 
% tensoriali di Chebyshev P(x,y,z) = T_i(x)*T_j(y)*T_k(z), sfruttando il 
% Teorema della Divergenza per convertire l'integrale di volume
% in un integrale di superficie sulle facce.
%
%**************************************************************************
%
% INPUT:
%
%   vertices:               Matrice [m x 3] delle coordinate dei vertici
%
%   facets:                 Matrice [n x 3] degli indici dei vertici 
%                           per ogni faccia
%
%   ade:                    Grado algebrico massimo
%
%   chebyshev_indices:      Matrice degli indici (i, j, k) della base tensoriale
%
%   dbox:                   Iper-rettangolo [min_x, max_x, ..., max_z] del dominio
%
%**************************************************************************
%
% OUTPUT:
%   moments:                Vettore colonna con i momenti integrati per ciascuna tripla.
%
%**************************************************************************

% Numero di punti di quadratura e numero di momenti richiesti
num_pts = size(nodes, 1);
num_moments = size(chebyshev_indices, 1);

X_ref = (2.0 * nodes - (dbox(1,:) + dbox(2,:))) ./ ...
    (dbox(2,:) - dbox(1,:));

% Grado massimo per asse
max_i = max(chebyshev_indices(:, 1));
max_j = max(chebyshev_indices(:, 2));
max_k = max(chebyshev_indices(:, 3));

idx_i = chebyshev_indices(:,1) + 1;
idx_j = chebyshev_indices(:,2) + 1;
idx_k = chebyshev_indices(:,3) + 1;

% Calcolo delle matrici di Chebyshev.
%
% La direzione x richiede un grado aggiuntivo poiche' la primitiva di
% T_i(x), per i >= 2, contiene T_(i+1)(x). Nelle altre direzioni e'
% sufficiente il grado massimo effettivamente richiesto.

x = X_ref(:,1);
TX = chebpolys(max_i + 1, x);
TY = chebpolys(max_j, X_ref(:,2));
TZ = chebpolys(max_k, X_ref(:,3));

% Costruzione diretta delle primitive nella direzione x
% Viene allocata solamente la matrice necessaria per i momenti richiesti
PhiX = zeros(num_pts, num_moments);

% Gli indici nella direzione x vengono separati in tre casi per evitare
% integrazione numerica e utilizzare le primitive analitiche
is0 = (chebyshev_indices(:,1) == 0);
is1 = (chebyshev_indices(:,1) == 1);
isN = (chebyshev_indices(:,1) >= 2);

% Primitiva di T_0(x) = 1
PhiX(:,is0) = repmat(x, 1, nnz(is0));

% Primitiva di T_1(x) = x
PhiX(:,is1) = repmat(0.5 .* x.^2, 1, nnz(is1));

% Primitiva analitica di T_i(x), per i >= 2
if any(isN)
    i = chebyshev_indices(isN,1).';
    PhiX(:,isN) = ...
        TX(:,i + 2) .* (i ./ (i.^2 - 1)) - ...
        x .* TX(:,i + 1) ./ (i - 1);
end

% Estrazione delle sole colonne necessarie nelle altre direzioni
TY_cols = TY(:,idx_j);
TZ_cols = TZ(:,idx_k);

% Costruzione vettorializzata dell'integrando
% Il peso della quadratura viene applicato direttamente alla primitiva
% nella direzione x
P = PhiX .* weights(:);
P = P .* TY_cols;
P = P .* TZ_cols;

% La trasformazione dalla coordinata fisica x alla coordinata di
% riferimento introduce il fattore dx = (xmax - xmin)/2 dxi
B1 = 0.5 * (dbox(2,1) - dbox(1,1));

% Somma dei contributi dei punti di quadratura
chebyshev_moms = B1 .* sum(P, 1).';
end