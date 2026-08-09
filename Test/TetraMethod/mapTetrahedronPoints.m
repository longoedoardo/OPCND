function [nodifisici, pesifisici] = mapTetrahedronPoints(nodi_rif, pesi_rif, V)
%**************************************************************************
%
% function [nodi_fisici, pesi_fisici] = mapTetrahedronPoints(nodi_rif, pesi_rif, V)
%
% esegue una trasformazione affine per trasportare i punti di integrazione
% definiti sul tetraedro canonico (vertici (0,0,0), (1,0,0), (0,1,0), (0,0,1))
% su un tetraedro generico definito dai vertici in V.
%
% INPUT:
%     nodi_rif  : Matrice (N x 3) dei nodi nel dominio di riferimento.
%                   Ogni riga rappresenta le coordinate baricentriche (r, s, t).
%     pesi_rif  : Vettore colonna (N x 1) dei pesi di quadratura associati.
%                   Si assume che sum(pesi_rif) = 1/6 (volume del rif.).
%     V         : Matrice (4 x 3) dei vertici del tetraedro target.
%                   V(1,:) è l'origine locale (A), le altre righe sono (B, C, D).
%
% OUTPUT:
%     nodi_fisici : Matrice (N x 3) dei nodi mappati nello spazio fisico (x, y, z).
%     pesi_fisici : Vettore (N x 1) dei pesi riscalati per il volume del target.
%
%**************************************************************************
%
% La trasformazione è definita dalla mappa affine: x = phi(xi) = A + J * xi
% dove J è lo Jacobiano [v1, v2, v3] composto dai vettori spigolo che partono
% dal primo vertice V(1,:). Pesi: I pesi fisici incorporano il determinante
% dello Jacobiano (detJ), che rappresenta il rapporto tra il volume fisico
% e quello di riferimento: pesi_fisici = pesi_rif * |det(J)|.
%
%**************************************************************************

% Calcolo vettori spigolo direttamente dalla matrice V
v1 = V(2,:) - V(1,:);
v2 = V(3,:) - V(1,:);
v3 = V(4,:) - V(1,:);

% Determinante rapido (prodotto misto)
detJ = v1(1)*(v2(2)*v3(3) - v2(3)*v3(2)) - ...
    v1(2)*(v2(1)*v3(3) - v2(3)*v3(1)) + ...
    v1(3)*(v2(1)*v3(2) - v2(2)*v3(1));

if abs(detJ) < eps
    warning('Tetraedro degenere o quasi-piatto.');
end

% Trasformazione affine vettorizzata: x = J*r + A
J_t = [v1; v2; v3]; % Questa è la trasposta della tua J
nodifisici = nodi_rif * J_t + V(1,:);

% I pesi scalano linearmente con il volume
pesifisici = pesi_rif * abs(detJ);
end