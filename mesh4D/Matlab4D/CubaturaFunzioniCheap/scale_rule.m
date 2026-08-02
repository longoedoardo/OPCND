function XYZW_tens = scale_rule(XYZTW_tens_ref, dbox)

%**************************************************************************
%
% function XYZW_tens = scale_rule(XYZTW_tens_ref, dbox)
%
% Scala i nodi di riferimento dall'ipercubo [-1,1]^4 alla Bounding Box 4D.
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

% Moltiplichiamo e sommiamo le prime 4 colonne in un colpo solo.
% MATLAB applica automaticamente centro (1x4) e semi_ampiezza (1x4) 
% a tutte le righe della matrice XYZTW_tens_ref (Mx4).
XYZW_tens(:, 1:4) = centro + semi_ampiezza .* XYZTW_tens_ref(:, 1:4);

% Se la matrice di riferimento ha 5 colonne, copiamo i pesi originali
% nella quinta colonna dell'output, senza modificarli (come da teoria).
if size(XYZTW_tens_ref, 2) == 5
    XYZW_tens(:, 5) = XYZTW_tens_ref(:, 5);
end

end