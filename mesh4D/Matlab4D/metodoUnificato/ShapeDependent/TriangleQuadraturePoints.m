function [IntGaussP, IntGaussW] = TriangleQuadraturePoints(nGP)
%**************************************************************************
%
% function [IntGaussP, IntGaussW] = TriangleQuadraturePoints(nGP)
%
% INPUT:
% - nGP       : Scalare, numero di punti di Gauss per ogni dimensione lineare
%            (nIntGP = nGP^2 punti totali).
%
% OUTPUT:
% - IntGaussP : Matrice (nGP^2 x 2) contenente le coordinate (r,s) sul
%               triangolo di riferimento di vertici (0,0), (1,0), (0,1).
% - IntGaussW : Vettore colonna (nGP^2 x 1) dei pesi di quadratura.
%               sum(IntGaussW) = 0.5 (area del triangolo di riferimento).
%
% Mappa il quadrato [0,1]^2 nel triangolo di riferimento tramite la
% trasformazione: r = mu1, s = mu2 * (1 - mu1).
% Lo Jacobiano di questa mappa è (1 - mu1); per compensarlo si integra
% mu1 con la quadratura di Gauss-Jacobi di peso (1-t)^1 (alpha=1, beta=0)
% e mu2 con Gauss-Legendre standard (alpha=0, beta=0).
%
% Algoritmo di Golub-Welsch: nodi e pesi calcolati tramite diagonalizzazione
% della matrice di Jacobi (funzione locale gaujac).
%
%**************************************************************************

nIntGP = nGP^2;
IntGaussP = zeros(nIntGP, 2);
IntGaussW = zeros(nIntGP, 1);

% Ottenimento posizioni (mu) e pesi (A) tramite Gauss-Jacobi
[mu1, A1] = gaujac(nGP, 1.0, 0.0);
[mu2, A2] = gaujac(nGP, 0.0, 0.0);

% Shift e rescale dall'intervallo [-1, 1] a [0, 1]
mu1 = 0.5 * mu1 + 0.5;  A1 = (0.5^2) * A1;
mu2 = 0.5 * mu2 + 0.5;  A2 = (0.5^1) * A2;

% Generazione dei punti tramite trasformazione di Duffy
iIntGP = 1;
for i = 1:nGP
    for j = 1:nGP
        m1 = mu1(i);
        m2 = mu2(j);

        IntGaussP(iIntGP, 1) = m1;
        IntGaussP(iIntGP, 2) = m2 * (1.0 - m1);

        IntGaussW(iIntGP) = A1(i) * A2(j);

        iIntGP = iIntGP + 1;
    end
end

end

function [x, w] = gaujac(n, alpha, beta)
% Calcola nodi (x) e pesi (w) di Gauss-Jacobi (formato colonna)
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