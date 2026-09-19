function XYZW = cub_gausscheb_tens4D(deg)
%**************************************************************************

% function XYZW = cub_gausscheb_tens4D(deg)

% Calcola i punti e i pesi per l'integrazione numerica (cubatura) in 4D.
% Utilizza una griglia a prodotto tensoriale basata sui nodi di Gauss-Chebyshev.

%**************************************************************************

n = ceil((deg + 1) / 2);
% Nodi di Chebyshev e pesi 1D

k = (1:n).';
angoli = ((2*k - 1) * pi) / (2 * n);
x = cos(angoli);       % n x 1
w = (pi / n) * ones(n,1); % n x 1 (costante)

% Griglia 3D e pesi 3D
[X, Y, Z, tau] = ndgrid(x, x, x, x);      % n x n x n
[Wx, Wy, Wz, Wtau] = ndgrid(w, w, w, w);   % n x n x n

W4 = Wx .* Wy .* Wz .* Wtau; % peso 3D

N4 = n^4; % numero totale di punti nella griglia 4D
XYZW = zeros(N4, 5);

XYZW(:,1) = reshape(X, N4, 1);
XYZW(:,2) = reshape(Y, N4, 1);
XYZW(:,3) = reshape(Z, N4, 1);
XYZW(:,4) = reshape(tau, N4, 1);
XYZW(:,5) = reshape(W4, N4, 1);
end
