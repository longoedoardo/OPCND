function [xv,yv,zv,R10,RFM]=maptopolygon2(XV,YV,ZV)
%**************************************************************************
%
% function [xv,yv,zv,R10,RFM]=maptopolygon2(XV,YV,ZV)
%
% Questa routine proietta un poligono 3D con vertici (XV, YV, ZV) su un 
% piano 2D (xy) tramite una trasformazione rigida (rotazione e, se 
% necessario, riflessione).
%**************************************************************************
% INPUT:
% XV, YV, ZV: Vettori delle coordinate dei vertici (non ripetere l'ultimo).
%
% OUTPUT:
% xv, yv, zv: Coordinate trasformate.
% R10: Matrice di rotazione 3x3.
% RFM: Matrice di riflessione 3x3 (identità se orientamento già corretto).
%**************************************************************************

% Estrazione delle coordinate dei primi tre vertici (A, B, C)
Ax = XV(1); Ay = YV(1); Az = ZV(1);
Bx = XV(2); By = YV(2); Bz = ZV(2);
Cx = XV(3); Cy = YV(3); Cz = ZV(3);

% Coefficienti del piano a*x + b*y + c*z + d = 0
a = (By - Ay) * (Cz - Az) - (Bz - Az) * (Cy - Ay);
b = (Bz - Az) * (Cx - Ax) - (Bx - Ax) * (Cz - Az);
c = (Bx - Ax) * (Cy - Ay) - (By - Ay) * (Cx - Ax);
d = -(Ax * a + Ay * b + Az * c);
% Calcolo della matrice di rotazione R10 (3x3)
ab_sq = a^2 + b^2;
R10 = eye(3); % Pre-impostata come identità

if ab_sq > 0
    den = sqrt(ab_sq + c^2);
    cos_theta = c / den;
    sin_theta = sqrt(ab_sq) / den;

    u1 = b / sqrt(ab_sq);
    u2 = -a / sqrt(ab_sq);
    term_cos = 1 - cos_theta;

    u1_sq = u1^2;
    u2_sq = u2^2;
    u1_u2_term = u1 * u2 * term_cos;
    u_sin1 = u2 * sin_theta;
    u_sin2 = u1 * sin_theta;

    R10(1,1) = cos_theta + u1_sq * term_cos;
    R10(1,2) = u1_u2_term;
    R10(1,3) = u_sin1;

    R10(2,1) = u1_u2_term;
    R10(2,2) = cos_theta + u2_sq * term_cos;
    R10(2,3) = -u_sin2;

    R10(3,1) = -u_sin1;
    R10(3,2) = u_sin2;
    R10(3,3) = cos_theta;
end

% Vettorializzazione totale della rotazione (nessun ciclo for)
coords = [XV(:), YV(:), ZV(:)] * R10';
xv = coords(:, 1);
yv = coords(:, 2);
zv = coords(:, 3);

% Calcolo dell'area tramite Shoelace vettorializzato (Gauss)
% Shift circolare degli indici per evitare cicli
x_next = xv([2:end, 1]);
y_next = yv([2:end, 1]);
signed_area = 0.5 * sum(xv .* y_next - x_next .* yv);

% Gestione della riflessione e orientamento
RFM = diag([1.0, 1.0, 1.0]);
if signed_area < 0
    RFM(2,2) = -1.0;
    yv = -yv; % Vettorizzato anziché iterato
end

end