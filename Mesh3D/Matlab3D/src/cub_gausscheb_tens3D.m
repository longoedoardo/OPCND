function XYZW = cub_gausscheb_tens3D(deg)

%**************************************************************************
%
% function XYZW = cub_gausscheb_tens3D(deg)
%
% Calcola i punti e i pesi per l'integrazione numerica in 3D.
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
%                   n = ceil((deg + 1) / 2)
%
%**************************************************************************
%
% OUTPUT:
%
%   XYZW:
%               Matrice Nq x 4 contenente i nodi e i pesi della quadratura
%               tensoriale 3D. Il numero totale di punti e' Nq = n^3 e il 
%               peso tensoriale e' costante ed uguale a w = (pi / n)^3.
%
%**************************************************************************

n = ceil((deg + 1) / 2);

% Nodi 1D di Gauss-Chebyshev
i = (1:n).';
x = cos((2*i - 1) * pi / (2*n));

% Griglia tensoriale 3D tramite ndgrid
[X, Y, Z] = ndgrid(x, x, x);

% Peso 3D costante
w3 = (pi / n)^3;

XYZW = [X(:), Y(:), Z(:), w3*ones(n^3,1)];