function V = dCHEBVAND(deg, X, chebyshev_indices)
%**************************************************************************
%
% [V,dbox] = dCHEBVAND(deg,X,dbox,chebyshev_indices)
%
% La funzione calcola la matrice di Chebyshev-Vandermonde per il grado "deg"
% sui punti "X".
%
% La base di Chebyshev utilizzata è la base tensoriale di grado totale "deg",
% traslata (shifted) sull'iper-rettangolo definito da "dbox".
%
%**************************************************************************

% default box
dbox = [-1 -1 -1; 1 1 1];
[m, d] = size(X); % m numero di punti e d dimensione spaziale
N = size(chebyshev_indices, 1); % numero di monomi, serve per le colonne
% della matrice di Vandermonde da costruire

% Mappatura punti in [-1,1] per ciascuna dimensione
mins = dbox(1, 1:d);
maxs = dbox(2, 1:d);
map = (2*X - maxs - mins) ./ (maxs - mins); % mappatura di ogni coordinata
% di X nell'intervallo [-1,1] per la dimensione che sto considerando

% Prealloca V con 1 (prodotto neutro)
V = ones(m, N);
for k = 1:d 
    xk = map(:, k); 
    T = chebpolys(deg, xk);
    idx = chebyshev_indices(:, k) + 1; 
    V = V .* T(:, idx); 
end
end