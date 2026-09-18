function [P3D, w_phys, norm_ext] = mapTriangleDunavantPoints( ...
    nodi_rif, pesi_rif, A, B, C)

%**************************************************************************
%
% function [P3D, w_phys, norm_ext] = mapTriangleDunavantPoints( ...
%                                      nodi_rif, pesi_rif, A, B, C)
%
% Esegue la trasformazione affine per trasportare i punti di integrazione
% dal triangolo di riferimento al triangolo fisico definito dai vertici
% A, B e C.
%
%**************************************************************************
% INPUT:
%   nodi_rif   - (Np x 2) Nodi di quadratura sul triangolo di riferimento.
%   pesi_rif   - (Np x 1) Pesi di quadratura sul triangolo di riferimento.
%   A,B,C      - (1 x 3) Coordinate dei vertici del triangolo fisico.
%
% OUTPUT:
%   P3D        - (Np x 3) Coordinate dei punti di quadratura nello spazio.
%   w_phys     - (Np x 1) Pesi fisici corretti per l'area del triangolo.
%   norm_ext   - (1 x 3) Normale unitaria esterna del triangolo.
%
%**************************************************************************

% Calcolo della geometria del triangolo

vec1 = B - A;
vec2 = C - A;

% Vettore normale non normalizzato

cp = cross(vec1, vec2);
area2 = norm(cp);

% Coordinate dei punti di riferimento

xrif = nodi_rif(:,1);
yrif = nodi_rif(:,2);

% Costruzione delle coordinate baricentriche

L1 = 1 - xrif - yrif;
L2 = xrif;
L3 = yrif;

L = [L1, L2, L3];

% Mappatura dei punti dal triangolo di riferimento
% al triangolo fisico

P3D = L(:,1) * A + L(:,2) * B + L(:,3) * C;

% Calcolo dei pesi fisici dai pesi di riferimento

w_phys = pesi_rif * (0.5 * area2);

% Calcolo della normale unitaria esterna

norm_ext = cp / area2;

end