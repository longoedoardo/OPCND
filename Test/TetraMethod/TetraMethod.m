function [P, w] = TetraMethod(V, F, nGP)

%**************************************************************************
%
% function [P, w] = TetraMethod(V, F, nGP)
%
% La funzione calcola le coordinate spaziali e i pesi di integrazione 
% per un oggetto solido definito dalla sua mesh superficiale. La funzione 
% scompone il volume interno in tetraedri e distribuisce i punti di Gauss 
% per permettere l'integrazione numerica di funzioni 3D.
%
% INPUT:
% - V: Matrice (N x 3) dei vertici della superficie originale.
% - F: Matrice (M x 3) delle facce triangolari (connettività).
% - nGP: Scalare, numero di punti di Gauss per dimensione lineare. Il 
% numero totale di punti per ogni tetraedro e' solitamente nGP^3.
%
% OUTPUT:
% - P: Matrice (TotPunti x 3) delle coordinate (x,y,z) di tutti i punti di Gauss.
% - w: Vettore colonna (TotPunti x 1) dei pesi di quadratura riscalati.
%
% Inizialmente otteniamo i nodi e i pesi sul tetraedro canonico tramite la 
% trasformazione di Duffy e polinomi di Gauss-Jacobi. Scompone il volume 
% interno collegando ogni faccia (F) al baricentro geometrico della mesh 
% (tramite 'TensMeshTetrahedron'). Trasferisce i punti dal dominio di 
% riferimento [0,1]^3 a ogni singolo tetraedro fisico calcolando lo 
% Jacobiano della trasformazione. Accumula i punti e i pesi in matrici 
% globali pre-allocate per massimizzare l'efficienza computazionale.
%
%**************************************************************************

nPointsPerTetra = nGP^3;
nTetra = size(F, 1);
totalPoints = nTetra * nPointsPerTetra;

% Calcolo nodi e pesi sul tetraedro di riferimento
[IntGaussP, IntGaussW] = TetrahedronQuadraturePoints(nPointsPerTetra, nGP);

% Scomposizione della mesh in tetraedri
V_tetr = TensMeshTetrahedron(V, F);

% Pre-allocazione matrice punti + colonna pesi
P = zeros(totalPoints, 3);
w = zeros(totalPoints, 1);

% Mappatura dei punti di quadratura
for k = 1:nTetra
% Mappa i punti dal tetraedro di riferimento a quello fisico corrente
[phys_nodes, phys_weights] = mapTetrahedronPoints(IntGaussP, IntGaussW, V_tetr(:,:,k));

% Calcolo indici di inserimento
idx_start = (k-1) * nPointsPerTetra + 1;
idx_end   = k * nPointsPerTetra;

% Memorizzazione risultati
P(idx_start:idx_end, :) = phys_nodes;
w(idx_start:idx_end)    = phys_weights;
end

end

