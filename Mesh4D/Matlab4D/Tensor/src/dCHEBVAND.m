function V = dCHEBVAND(deg, X, chebyshev_indices)
%**************************************************************************
% 
% V = dCHEBVAND(deg, X, chebyshev_indices)
%
% Calcola la matrice di Chebyshev-Vandermonde per il grado "deg"
% sui punti "X"
%
%**************************************************************************
%
% INPUT:
%
%   deg:
%               Grado massimo dei polinomi di Chebyshev da calcolare.
%
%   X:
%               Matrice m x d contenente le coordinate dei punti in cui
%               valutare la base. Ogni riga rappresenta un punto e ogni
%               colonna una coordinata.
%
%   chebyshev_indices:
%               Matrice N x d contenente gli indici della base tensoriale
%               di Chebyshev. 
%
%**************************************************************************
%
% OUTPUT:
%
%   V:
%               Matrice m x N di Chebyshev--Vandermonde, nella quale
%               l'elemento V(i,j) contiene il valore del j-esimo
%               polinomio della base nel punto X(i,:).
%
%**************************************************************************

[m, d] = size(X); % m = numero di punti
N = size(chebyshev_indices, 1); % Numero di monomi (colonne di V)

% Definizione automatica del box di riferimento [-1, 1] per le 'd' dimensioni
dbox = [-ones(1, d); ones(1, d)];

mins = dbox(1, :);
maxs = dbox(2, :);

% Mappatura vettorializzata di tutte le coordinate nell'intervallo [-1, 1]
map = (2.0 * X - maxs - mins) ./ (maxs - mins); 

% Prealloca la matrice di Vandermonde con 1 (elemento neutro moltiplicativo)
V = ones(m, N); 

% Ciclo sulle dimensioni (x, y, z, tau)
for k = 1 : d
    xk = map(:, k); % Coordinate della k-esima dimensione
    T = chebpolys(deg, xk); % Calcola tutti i polinomi di Chebyshev fino a 'deg'
    idx = chebyshev_indices(:, k) + 1; % Estrae gli indici degli esponenti
    V = V .* T(:, idx); % Moltiplicazione vettoriale sulle colonne della base
end

end
