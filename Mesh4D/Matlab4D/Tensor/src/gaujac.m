function [x, w] = gaujac(n, alpha, beta)
%**************************************************************************
%
% function [x, w] = gaujac(n, alpha, beta)
% 
% La funzione calcola i nodi (x) e i pesi (w) della quadratura di 
% Gauss-Jacobi nell'intervallo di definizione [-1, 1]. 
%
% Viene applicato l'algoritmo di Golub-Welsch: nodi e pesi calcolati 
% tramite diagonalizzazione della matrice di Jacobi
%
%**************************************************************************

if n == 0, x = []; w = []; return; end
if n == 1
    x = (beta - alpha) / (alpha + beta + 2);
    w = 2^(alpha + beta + 1) * gamma(alpha + 1) * gamma(beta + 1) / gamma(alpha + beta + 2);
    return;
end

% Indici
i = (1:n-1)';
abi = alpha + beta + 2*i;

% Coefficienti diagonali (aa)
aa = zeros(n, 1);
aa(1) = (beta - alpha) / (alpha + beta + 2);
aa(2:n) = (beta^2 - alpha^2) ./ (abi .* (abi + 2));

% Coefficienti sub-diagonali (bb)
bb = zeros(n-1, 1);
bb(1) = 2 * sqrt(1 * (1 + alpha) * (1 + beta) / ((alpha + beta + 2)^2 * (alpha + beta + 3)));
bb(2:n-1) = 2 ./ (abi(1:end-1) + 2) .* sqrt(i(2:end) .* (i(2:end) + alpha) .* ...
    (i(2:end) + beta) .* (i(2:end) + alpha + beta) ./ ...
    ((abi(1:end-1) + 1) .* (abi(1:end-1) + 3)));

% Costruzione matrice di Jacobi (n x n)
J = diag(aa) + diag(bb, 1) + diag(bb, -1);

% Autovalori e pesi
[V, D] = eig(J);
[x, idx] = sort(diag(D));
V = V(:, idx);

% Fattore di normalizzazione
factor = 2^(alpha + beta + 1) * gamma(alpha + 1) * gamma(beta + 1) / gamma(alpha + beta + 2);
w = factor * (V(1, :)').^2;
end