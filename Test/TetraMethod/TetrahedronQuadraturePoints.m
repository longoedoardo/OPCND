function [IntGaussP, IntGaussW] = TetrahedronQuadraturePoints(nIntGP, nGP)

%**************************************************************************
%
% function [IntGaussP, IntGaussW] = TetrahedronQuadraturePoints(nIntGP, nGP)
%
% INPUT:
% - nIntGP : Scalare, numero totale di punti di integrazione (nGP^3).
% - nGP    : Scalare, numero di punti di Gauss per ogni dimensione lineare.
%
% OUTPUT:
% - IntGaussP : Matrice (nIntGP x 3) contenente le coordinate (x,y,z).
% - IntGaussW : Vettore colonna (nIntGP x 1) dei pesi di quadratura.
%
% Mappa il dominio cubico [0,1]^3 nel dominio tetraedrico per gestire i 
% limiti di integrazione variabili: x = mu1, y = mu2 * (1 - mu1), 
% z = mu3 * (1 - mu2) * (1 - mu1). Quadratura di Gauss-Jacobi: Per 
% compensare lo Jacobiano della mappatura, vengono utilizzati polinomi di 
% Jacobi P^(alpha, beta) con parametri:
%         Dim 1: alpha=2, beta=0 (mu1)
%         Dim 2: alpha=1, beta=0 (mu2)
%         Dim 3: alpha=0, beta=0 (mu3, equivalente a Gauss-Legendre)
%
% Algoritmo di Golub-Welsch: I nodi e i pesi sono calcolati tramite la
% diagonalizzazione della matrice di Jacobi.
%
%**************************************************************************

% Pre-allocazione (nIntGP x 3 per i nodi, nIntGP x 1 per i pesi)
IntGaussP = zeros(nIntGP, 3);
IntGaussW = zeros(nIntGP, 1);

tol = 1e-10; % Tolleranza per verifica coerenza volume

% Ottenimento posizioni (mu) e pesi (A) tramite Gauss-Jacobi
[mu1, A1] = gaujac(nGP, 2.0, 0.0);
[mu2, A2] = gaujac(nGP, 1.0, 0.0);
[mu3, A3] = gaujac(nGP, 0.0, 0.0);

% Shift e rescale dall'intervallo [-1, 1] a [0, 1]
mu1 = 0.5 * mu1 + 0.5;  A1 = (0.5^3) * A1;
mu2 = 0.5 * mu2 + 0.5;  A2 = (0.5^2) * A2;
mu3 = 0.5 * mu3 + 0.5;  A3 = (0.5^1) * A3;

% Generazione dei punti tramite trasformazione di Duffy
iIntGP = 1;
for i = 1:nGP
    for j = 1:nGP
        for k = 1:nGP
            % Coordinate (x, y, z) - Formato riga
            m1 = mu1(i);
            m2 = mu2(j);
            m3 = mu3(k);
            
            IntGaussP(iIntGP, 1) = m1;
            IntGaussP(iIntGP, 2) = m2 * (1.0 - m1);
            IntGaussP(iIntGP, 3) = m3 * (1.0 - m2) * (1.0 - m1);
            
            % Peso associato
            IntGaussW(iIntGP) = A1(i) * A2(j) * A3(k);
            
            iIntGP = iIntGP + 1;
        end
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

% Coefficienti diagonali (aa) - Assicuriamoci che sia un vettore colonna n x 1
aa = zeros(n, 1);
aa(1) = (beta - alpha) / (alpha + beta + 2);
aa(2:n) = (beta^2 - alpha^2) ./ (abi .* (abi + 2));

% Coefficienti sub-diagonali (bb) - Vettore (n-1) x 1
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