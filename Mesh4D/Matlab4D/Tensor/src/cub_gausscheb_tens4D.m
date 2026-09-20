function XYZTW = cub_gausscheb_tens4D(deg)
%**************************************************************************
%
% function XYZW = cub_gausscheb_tens4D(deg)
%
% Calcola i punti e i pesi per l'integrazione numerica (cubatura) in 4D.
% Utilizza una griglia a prodotto tensoriale basata sui nodi di Gauss-Chebyshev.
%
%**************************************************************************
%
% INPUT:
%
%   deg:
%               Grado massimo della base polinomiale per la quale si
%               desidera costruire la griglia di cubatura.
%
%               Il numero di nodi nella singola direzione e' scelto come
%
%                   n = ceil((deg + 1) / 2),
%
%**************************************************************************
%
% OUTPUT:
%
%   XYZTW:
%               Matrice Nq x 5 contenente i nodi e i pesi della quadratura
%               tensoriale 4D. Le prime quattro colonne contengono le 
%               coordinate [x, y, z, tau],
%               mentre la quinta colonna contiene il peso associato a
%               ciascun nodo. Il peso tensoriale
%               e' costante ed uguale a w = (pi / n)^4.
%
%**************************************************************************

n = ceil((deg + 1) / 2);

% Nodi di Chebyshev e pesi 1D

k = (1:n).';
angoli = ((2*k - 1) * pi) / (2 * n);
x = cos(angoli);       % n x 1
w = (pi / n) * ones(n,1); % n x 1 (costante)

% Numero totale di punti
N4 = n^4;

% Coordinate della griglia tensoriale
X = kron(ones(n^3,1), x);
Y = kron(ones(n^2,1), kron(x, ones(n,1)));
Z = kron(ones(n,1), kron(x, ones(n^2,1)));
tau = kron(x, ones(n^3,1));

% Peso tensoriale costante
w4 = (pi / n)^4;

% Costruzione dell'output
XYZTW = zeros(N4, 5);

XYZTW(:,1) = X;
XYZTW(:,2) = Y;
XYZTW(:,3) = Z;
XYZTW(:,4) = tau;
XYZTW(:,5) = w4;

end
