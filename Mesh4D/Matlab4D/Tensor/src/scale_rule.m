function XYZW_tens = scale_rule(XYZTW_tens_ref, dbox)

%**************************************************************************
%
% function XYZW_tens = scale_rule(XYZTW_tens_ref, dbox)
%
% Scala i nodi di riferimento dall'ipercubo [-1,1]^4 alla Bounding Box 4D.
%
%**************************************************************************
%
% INPUT:
%
%   XYZW_tens_ref:
%               Matrice Nq x 4 oppure Nq x 5 contenente i nodi della
%               regola di quadratura nel dominio di riferimento
%
%   dbox:
%               Bounding box 4D
%
%**************************************************************************
%
% OUTPUT:
%
%   XYZW_tens:
%               Matrice Nq x 4 oppure Nq x 5 contenente i nodi trasformati
%               nella bounding box 4D.
%
%               Le prime quattro colonne contengono le coordinate
%               trasformate [x, y, z, tau]. Se presente nella matrice di
%               ingresso, la quinta colonna contiene invariato il peso
%               associato a ciascun nodo.
%
%**************************************************************************

% Preallocazione della matrice finale
XYZW_tens = zeros(size(XYZTW_tens_ref));

% dbox è una matrice 2 x 4:
% Riga 1: [min_x, min_y, min_z, min_tau]
% Riga 2: [max_x, max_y, max_z, max_tau]
min_val = dbox(1, :); % Vettore riga 1x4
max_val = dbox(2, :); % Vettore riga 1x4

% Calcolo vettoriale del centro e della semi-ampiezza per i 4 assi
centro = (min_val + max_val) / 2; % Spostamento (1 x 4)
semi_ampiezza = (max_val - min_val) / 2;  % Dilatazione (1 x 4)

XYZW_tens(:, 1:4) = centro + semi_ampiezza .* XYZTW_tens_ref(:, 1:4);

if size(XYZTW_tens_ref, 2) == 5
    XYZW_tens(:, 5) = XYZTW_tens_ref(:, 5);
end

end