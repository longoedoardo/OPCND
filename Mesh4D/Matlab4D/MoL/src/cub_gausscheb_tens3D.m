function XYZW = cub_gausscheb_tens3D(deg)

%**************************************************************************

% function XYZW = cub_gausscheb_tens3D(deg)

% Calcola i punti e i pesi per l'integrazione numerica (cubatura) in 3D.
% Utilizza una griglia a prodotto tensoriale basata sui nodi di Gauss-Chebyshev.

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