function [P3D, w_phys, norm_ext] = TriangleDunavantQuadraturePoints(rule, A, B, C)

%*************************************************************************
%
% [P3D, w_phys, norm_ext] = DunavantQuad(rule, A, B, C)
%
% Mappatura di una regola di cubatura di Dunavant su un triangolo.
%
% INPUT:
%   rule  - (int) Grado della regola di quadratura Dunavant.
%   A,B,C - (1x3) Coordinate dei vertici del triangolo.
%
% OUTPUT:
%   P3D      - (Np x 3) Coordinate dei punti di quadratura nello spazio.
%   w_phys   - (Np x 1) Pesi fisici corretti per l'area del triangolo.
%   norm_ext - (1 x 3) Normale unitaria esterna del triangolo.
%
% Per poter utilizzare la seguente funzione, sono necessarie gli script
% contenuti nella cartella "Dunavant"
%*************************************************************************

addpath("Dunavant/")

% Calcolo della geometria del triangolo
% vettori dei due lati
vec1 = B - A; % 1x3
vec2 = C - A;
% vettore normale non normalizzato
cp = cross(vec1, vec2); % 1x3
area2 = norm(cp); % = 2 * area

order_num = dunavant_order_num(rule); % Calcolo indice della regola
[xyrif, wrif] = dunavant_rule(rule, order_num); % ottengo le coordinate di riferimento
% e i pesi di riferimento sul triangolo di riferimento

% Formatto le colonne, converto tutto
xrif = xyrif(1,:).';
yrif = xyrif(2,:).';
if isrow(wrif); wrif = wrif.'; end

% Costruzione delle coordinate baricentrice
L1 = 1 - xrif - yrif;
L2 = xrif;
L3 = yrif;
L = [L1, L2, L3];

P3D = L(:,1) * A + L(:,2) * B + L(:,3) * C; % produce i punti di quadratura
% Sfrutto la combinazione lineare baricentrica. I punti di riferimento
% vengono mappati sul triangolo nello spazio

% Calcolo dei pesi fisici dai pesi di riferimento
w_phys = wrif * (0.5 * area2);   % Np x 1

% Calcolo della normale unitaria esterna
norm_ext = cp / area2;  % 1x3

end
