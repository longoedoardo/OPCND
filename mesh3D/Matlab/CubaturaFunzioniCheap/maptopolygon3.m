function [XV, YV, ZV, n] = maptopolygon3(xv, yv, zv, R10, RFM)
%**************************************************************************
%
% function [XV, YV, ZV, n] = maptopolygon3(xv, yv, zv, R10, RFM)
%
% Questa routine prende un poligono piano nello spazio bidimensionale (xv, yv)
% e lo riporta nella sua posizione originale nello spazio tridimensionale 
% (XV, YV, ZV) mediante un "moto rigido" inverso.
%
% Il vettore "zv" e la matrice di rotazione "R10" permettono di riportare
% indietro i punti dal piano locale xy alla giacitura 3D originale.
%**************************************************************************
% INPUT:
%
% xv, yv: Coordinate dei vertici di un poligono nello spazio bidimensionale
%   ottenuto tramite un "moto rigido" di un poligono con vertici (XV,YV,ZV);
%   l'ultimo vertice non deve essere ripetuto; sono vettori colonna.
%
% R10, RFM: Matrici di rotazione e riflessione coinvolte nel processo.
%   Risulta "P = R10' * (RFM * P2 + PS)" dove "PS" è uno specifico spostamento,
%   "P2" è il punto sul piano xy e "P" è il punto risultante nello spazio 3D.
%
%**************************************************************************
% OUTPUT:
%
% XV, YV, ZV: Coordinate dei vertici di un poligono nello spazio 
%   tridimensionale (es. i vertici di una faccia di un tetraedro); l'ultimo
%   vertice non deve essere ripetuto; sono vettori colonna.
%
% n: Versore della normale esterna al piano contenente [XV YV ZV].
%
%**************************************************************************

% Preparazione delle coordinate locali in una matrice [N x 3]
% Nota: la coordinata z locale è data da zv (che rappresenta l'altezza del piano)
n_punti = length(xv);
P2 = [xv(:), yv(:), zeros(n_punti, 1)];

% Applicazione della riflessione in blocco (P_ref = P2 * RFM')
P_ref = P2 * RFM';

% Aggiunta dello spostamento PS (la coordinata Z viene traslata di zv)
P_shifted = P_ref;
P_shifted(:, 3) = P_shifted(:, 3) + zv(:);

% Rotazione inversa e ritorno allo spazio 3D (P_final = P_shifted * R10)
% Poiché cerchiamo P = R10' * P_shifted, trasposto diventa P_final' = P_shifted_matrix * R10
P_final = P_shifted * R10;

XV = P_final(:, 1);
YV = P_final(:, 2);
ZV = P_final(:, 3);

% La normale originale ee = [0, 0, 1] trasformata: n = R10' * (RFM * ee)
% RFM * [0; 0; 1] corrisponde alla terza colonna di RFM
re = RFM(:, 3)';

% Moltiplicazione per R10' in termini vettoriali (re * R10)
n_raw = re * R10;

% Normalizzazione del vettore (norma Euclidea)
norma = norm(n_raw);
n = n_raw / norma;

end