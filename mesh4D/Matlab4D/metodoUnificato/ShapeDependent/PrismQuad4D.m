function [XYZTW, WV_X] = PrismQuad4D(V, XI_ref, ETA_ref, T_ref, W_ref, baricentro)

%**************************************************************************
%
% [XYZTW, WV_X] = PrismQuad4D(V, XI_ref, ETA_ref, T_ref, W_ref, baricentro)
%
% Costruisce i nodi fisici e i pesi orientati per integrare su una
% iperfaccia laterale spazio-temporale 3D immersa in R^4.
%
% Il prisma di riferimento e':
%
% K_ref = {(xi,eta,t) : xi >= 0, eta >= 0, xi + eta <= 1, 0 <= t <= 1}.
%
% La mappa verso R^4 e' X(xi,eta,t) = (1-t) X0(xi,eta) + t X1(xi,eta),
%
% con:
%
%   X0(xi,eta) = A0 + xi (B0-A0) + eta (C0-A0),
%   X1(xi,eta) = A1 + xi (B1-A1) + eta (C1-A1).
%
% Il vettore superficie orientato locale e' N(xi,eta,t) = n(xi,eta,t) dS,
% ortogonale ai tre vettori tangenti dX/dxi, dX/deta, dX/dt. Per il teorema 
% della divergenza applicato al campo F = (Phi_x, 0, 0, 0), serve integrare:
% Phi_x * n_x dS. Quindi la funzione restituisce direttamente: WV_X(q) = 
% W_ref(q) * N_x(q).
%
% INPUT:
%   V           matrice 6 x 4 dei vertici del prisma spazio-temporale:
%                 V(1,:) = A0, vertice 1 della faccia a tau = 0
%                 V(2,:) = B0, vertice 2 della faccia a tau = 0
%                 V(3,:) = C0, vertice 3 della faccia a tau = 0
%                 V(4,:) = A1, vertice 1 della faccia a tau = 1
%                 V(5,:) = B1, vertice 2 della faccia a tau = 1
%                 V(6,:) = C1, vertice 3 della faccia a tau = 1
%               Ogni riga e' [x, y, z, tau].
%
%   XI_ref      vettore Nq x 1 delle coordinate xi dei nodi di quadratura
%   ETA_ref     vettore Nq x 1 delle coordinate eta dei nodi di quadratura
%   T_ref       vettore Nq x 1 delle coordinate temporali in [0,1]
%   W_ref       vettore Nq x 1 dei pesi sul prisma di riferimento
%   baricentro  vettore 1 x 4 interno al dominio 4D, usato per orientare
%               N verso l'esterno
%
% OUTPUT:
%   XYZT W      matrice Nq x 4 dei nodi fisici 4D [x, y, z, tau]
%
%   WV_X        vettore Nq x 1 dei pesi orientati:
%                 WV_X = W_ref .* N_x
%
%**************************************************************************

% Vertici della faccia triangolare al tempo iniziale tau = 0
A0 = V(1, :);
B0 = V(2, :);
C0 = V(3, :);

% Vertici corrispondenti della faccia triangolare al tempo finale tau = 1
A1 = V(4, :);
B1 = V(5, :);
C1 = V(6, :);

% Vettori lato della faccia iniziale
lato1_base = B0 - A0;
lato2_base = C0 - A0;

% Vettori lato della faccia finale
lato1_top = B1 - A1;
lato2_top = C1 - A1;

% Nodi sulla faccia iniziale X0(xi,eta)
nodi_base_4D = A0 + XI_ref * lato1_base + ETA_ref * lato2_base;

% Nodi sulla faccia finale X1(xi,eta)
nodi_top_4D  = A1 + XI_ref * lato1_top  + ETA_ref * lato2_top;

% Nodi fisici sulla faccia laterale spazio-temporale:
% X(xi,eta,t) = (1-t) X0(xi,eta) + t X1(xi,eta)
XYZTW = (1.0 - T_ref) .* nodi_base_4D + T_ref .* nodi_top_4D;

% Vettori tangenti rispetto a xi, era
dXYZTW_dxi  = (1.0 - T_ref) .* lato1_base + T_ref .* lato1_top;
dXYZTW_deta = (1.0 - T_ref) .* lato2_base + T_ref .* lato2_top;

% Vettore tangente rispetto al tempo
% Dipende da xi ed eta, per questo la normale, in generale, non e' costante.
dXYZTW_dt   = nodi_top_4D - nodi_base_4D;

% c = dX/dxi, r = dX/deta, z = dX/dt
c = dXYZTW_dxi;
r = dXYZTW_deta;
z = dXYZTW_dt;

c1 = c(:,1); c2 = c(:,2); c3 = c(:,3); c4 = c(:,4);
r1 = r(:,1); r2 = r(:,2); r3 = r(:,3); r4 = r(:,4);
z1 = z(:,1); z2 = z(:,2); z3 = z(:,3); z4 = z(:,4);


% Minori 3 x 3 ottenuti eliminando, rispettivamente, la colonna
% x, y, z, tau dalla matrice jacobiana 3 x 4:
%
%   J = [c; r; z]
%
M1 = c2.*(r3.*z4 - r4.*z3) ...
    - c3.*(r2.*z4 - r4.*z2) ...
    + c4.*(r2.*z3 - r3.*z2);

M2 = c1.*(r3.*z4 - r4.*z3) ...
    - c3.*(r1.*z4 - r4.*z1) ...
    + c4.*(r1.*z3 - r3.*z1);

M3 = c1.*(r2.*z4 - r4.*z2) ...
    - c2.*(r1.*z4 - r4.*z1) ...
    + c4.*(r1.*z2 - r2.*z1);

M4 = c1.*(r2.*z3 - r3.*z2) ...
    - c2.*(r1.*z3 - r3.*z1) ...
    + c3.*(r1.*z2 - r2.*z1);

% Vettore superficie 4D con i segni alternati:
%
%   N = [M1, -M2, M3, -M4] = n dS
%
% Ogni riga di Nvec e' ortogonale a dX/dxi, dX/deta e dX/dt

Nvec = [M1, -M2, M3, -M4];

% Orientamento uscente, se N punta verso il baricentro, viene ribaltato
direzione_uscente = XYZTW - baricentro;
flip = sum(Nvec .* direzione_uscente, 2) < 0;
Nvec(flip, :) = -Nvec(flip, :);

% Peso orientato richiesto dalla formula di divergenza con
% F = (Phi_x, 0, 0, 0), F dot n dS = Phi_x * N_x

WV_X = W_ref .* Nvec(:, 1);

end