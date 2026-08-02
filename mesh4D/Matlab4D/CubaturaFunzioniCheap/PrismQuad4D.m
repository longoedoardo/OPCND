function [XYZTW, WV_CUB] = PrismQuad4D(V, XI_ref, ETA_ref, T_ref, W_ref)
%**************************************************************************
% function [XYZTW, WV_CUB] = PrismQuad4D(V, XI_ref, ETA_ref, T_ref, W_ref)
%
% DESCRIZIONE:
%   Trasforma una regola di cubatura definita su un prisma di riferimento
%   in un prisma 4D deformato (spazio-tempo). Calcola direttamente il
%   flusso orientato locale (componente N1 del prodotto vettoriale 
%   generalizzato 4D) punto per punto, allineato alla versione Fortran.
%
% INPUT:
%   V          - (6x4) Matrice vertici [x, y, z, tau] del prisma fisico.
%   XI_ref, ETA_ref - (Np x 1) Coordinate spaziali di riferimento (triangolo).
%   T_ref      - (Np x 1) Coordinata temporale di riferimento [0, 1].
%   W_ref      - (Np x 1) Pesi della cubatura di riferimento.
%
% OUTPUT:
%   XYZTW      - (Np x 4) Punti mappati nello spazio fisico 4D.
%   WV_CUB     - (Np x 1) Pesi fisici di flusso orientato locale (N1 * W_ref).
%**************************************************************************
num_pts = length(W_ref);

% Calcoliamo la base inferiore e superiore del prisma in modo vettoriale
X_base = zeros(num_pts, 4);
X_top  = zeros(num_pts, 4);
for i = 1:4
    X_base(:, i) = V(1, i) + XI_ref * (V(2, i) - V(1, i)) + ETA_ref * (V(3, i) - V(1, i));
    X_top(:, i)  = V(4, i) + XI_ref * (V(5, i) - V(4, i)) + ETA_ref * (V(6, i) - V(4, i));
end

% Mappatura 4D finale: media pesata tra base e tetto in funzione del tempo T
XYZTW = zeros(num_pts, 4);
for i = 1:4
    XYZTW(:, i) = (1.0 - T_ref) .* X_base(:, i) + T_ref .* X_top(:, i);
end

% Vettori tangenti locali (dipendono da T_ref: variano punto per punto se la mesh si deforma)
g_xi  = zeros(num_pts, 4);
g_eta = zeros(num_pts, 4);
g_tau = zeros(num_pts, 4);
for i = 1:4
    g_xi(:, i)  = (1.0 - T_ref) .* (V(2, i) - V(1, i)) + T_ref .* (V(5, i) - V(4, i));
    g_eta(:, i) = (1.0 - T_ref) .* (V(3, i) - V(1, i)) + T_ref .* (V(6, i) - V(4, i));
    g_tau(:, i) = X_top(:, i) - X_base(:, i);
end

% Componente x (indice 1) del prodotto vettoriale generalizzato 4D di (g_xi, g_eta, g_tau)
N1 = g_xi(:,2) .* (g_eta(:,3).*g_tau(:,4) - g_eta(:,4).*g_tau(:,3)) ...
    - g_xi(:,3) .* (g_eta(:,2).*g_tau(:,4) - g_eta(:,4).*g_tau(:,2)) ...
    + g_xi(:,4) .* (g_eta(:,2).*g_tau(:,3) - g_eta(:,3).*g_tau(:,2));

% Pesi fisici di flusso localmente orientati
WV_CUB = W_ref .* N1;
end