function [XYZTW, WV_X] = ShiftingPrismQuadrature(V, XI_ref, ETA_ref, T_ref, W_ref)

%**************************************************************************
%
% [XYZTW, WV_X] = ShiftingPrismQuadrature(V, XI_ref, ETA_ref, T_ref, W_ref)
%
% Costruisce i nodi fisici e i pesi orientati per integrare su una
% iperfaccia laterale spazio-temporale 3D immersa in R^4.
%
% Il prisma di riferimento e':
%
% P = {(xi,eta,t) : xi >= 0, eta >= 0, xi + eta <= 1, 0 <= t <= 1}.
% L'iperfaccia e' ottenuta collegando linearmente due configurazioni
% triangolari corrispondenti mediante la parametrizzazione
%
%       X(xi,eta,tau) = (1-tau) X_0(xi,eta) + tau X_1(xi,eta),
%
% dove X_0 e X_1 sono le parametrizzazioni affini delle facce triangolari
% alle configurazioni iniziale e finale.
% La componente x della normale generalizzata all'iperfaccia e' calcolata
% mediante il minore 3x3 associato alla prima coordinata. Poiche' la
% coordinata temporale della mappa coincide con tau, si ottiene
%
%       N_x = y_xi z_eta - z_xi y_eta.
%
% I pesi restituiti includono il fattore geometrico orientato associato
% alla componente N_x della normale:
%
%       WV_X = W_ref .* N_x
%
%**************************************************************************
%
% INPUT:
%
%   V:
%               Matrice 6 x 4 contenente i vertici delle due facce
%               triangolari del prisma spazio-temporale.
%
%               Le prime tre righe rappresentano la configurazione
%               iniziale a tau = 0:
%
%                   V(1,:) = A_0
%                   V(2,:) = B_0
%                   V(3,:) = C_0
%
%               Le ultime tre righe rappresentano la configurazione
%               finale a tau = 1:
%
%                   V(4,:) = A_1
%                   V(5,:) = B_1
%                   V(6,:) = C_1
%
%               Ogni riga ha la forma [x,y,z,tau].
%
%   XI_ref:
%               Vettore Nq x 1 contenente le coordinate xi dei nodi di
%               quadratura sul triangolo di riferimento.
%
%   ETA_ref:
%               Vettore Nq x 1 contenente le coordinate eta dei nodi di
%               quadratura sul triangolo di riferimento.
%
%   T_ref:
%               Vettore Nq x 1 contenente le coordinate temporali tau dei
%               nodi di quadratura, con tau appartenente a [0,1].
%
%   W_ref:
%               Vettore Nq x 1 contenente i pesi della regola di quadratura
%               sul prisma di riferimento.
%
%**************************************************************************
%
% OUTPUT:
%
%   XYZTW:
%               Matrice Nq x 4 contenente i nodi fisici dell'iperfaccia
%               spazio-temporale nella forma
%
%                   [x,y,z,tau].
%
%   WV_X:
%               Vettore Nq x 1 contenente i pesi orientati associati alla
%               componente x della normale generalizzata:
%
%                   WV_X = W_ref .* N_x.
%
%**************************************************************************

% Le due facce triangolari sono le sezioni del prisma per tau = 0 e tau = 1
% Esse rappresentano la stessa faccia osservata nelle due configurazioni
% estreme del modo. La parametrizzazione spazio-temporale verra' costruita
% collegando punti corrispondenti delle due facce. 

% Vertici della faccia triangolare al tempo iniziale tau = 0

A0 = V(1, :);
B0 = V(2, :);
C0 = V(3, :);

% Vertici corrispondenti della faccia triangolare al tempo finale tau = 1
A1 = V(4, :);
B1 = V(5, :);
C1 = V(6, :);

% Per il triangolo iniziale definiamo 
%
%       e_1^0 = B_0 - A_0,
%       e_2^0 = C_0 - A_0,
%
% mentre per il triangolo finale
%
%       e_1^1 = B_1 - A_1,
%       e_2^1 = C_1 - A_1,
%
% Questi vettori costituiscono le due direzioni tangenti associate alle
% coordinate locali xi ed eta.

% Vettori lato della faccia iniziale
lato1_base = B0 - A0;
lato2_base = C0 - A0;

% Vettori lato della faccia finale
lato1_top = B1 - A1;
lato2_top = C1 - A1;


% Sul triangolo di riferimento
%
%       T_ref = {(xi,eta) : xi >= 0, eta >= 0, xi+eta <= 1},
%
% la mappa affine standard e'
%
%       X_0(xi,eta) = A_0 + xi (B_0-A_0) + eta (C_0-A_0).
%
% Per xi = 0, eta = 0 si ottiene A_0.
%
% Per xi = 1, eta = 0 si ottiene B_0.
%
% Per xi = 0, eta = 1 si ottiene C_0.
%
% Questa e' l'analoga 4D della parametrizzazione affine di un triangolo
% utilizzata nella funzione mapTrianglePoints nella versione 3D

% Nodi sulla faccia iniziale X0(xi,eta)
nodi_base_4D = A0 + XI_ref * lato1_base + ETA_ref * lato2_base;


% In maniera completamente analoga si costruisce
%
%       X_1(xi,eta) = A_1 + xi (B_1-A_1) + eta (C_1-A_1).
%
% La stessa coppia (xi,eta) identifica quindi due punti corrispondenti: 
% X_0(xi,eta) e X_1(xi,eta). Questi due punti verranno successivamente 
% collegati lungo la direzione temporale.

% Nodi sulla faccia finale X1(xi,eta)
nodi_top_4D  = A1 + XI_ref * lato1_top  + ETA_ref * lato2_top;

% Facciamo dunque ora una interpolazione lineare spazio-temporale. La 
% parametrizzazione completa dell'iperfaccia laterale e'
%
%       X(xi,eta,tau) = (1-tau) X_0(xi,eta) + tau X_1(xi,eta).
%
% Questa formula puo' essere interpretata come un prodotto tensoriale
% geometrico tra:
%   
%   - la parametrizzazione triangolare nelle variabili xi, eta;
%   - la base lineare {1-tau, tau} nella variabile temporale.
%
% Infatti
%
%       X(xi,eta,tau) = (1-tau) X_0(xi,eta) + tau X_1(xi,eta)
%
% e' precisamente l'interpolazione lineare di X_0 e X_1 lungo tau.
%
% Per tau = 0 otteniamo X(xi,eta,0) = X_0(xi,eta), mentre per tau = 1 si 
% ottiene X(xi,eta,1) = X_1(xi,eta).
%
% In questo modo il triangolo iniziale e quello finale vengono collegati
% generando una iperfaccia tridimensionale in R^4.
%
% L'operazione e' completamente vettorializzata rispetto ai Nq nodi
% della quadratura.

XYZTW = (1.0 - T_ref) .* nodi_base_4D + T_ref .* nodi_top_4D;

% Inizio del calcolo dei vettori tangenti rispetto a xi ed eta.
% Partendo da
%
%       X(xi,eta,tau) = (1-tau) X_0(xi,eta) + tau X_1(xi,eta),
%
% e ricordando che X_0 e X_1 sono affine rispetto a xi ed eta, otteniamo
%
%       dX/dxi = (1-tau) dX_0/dxi + tau dX_1/dxi
%
% e analogamente
%
%       dX/deta = (1-tau) dX_0/deta + tau dX_1/deta.
%
% Poiche'
%
%       dX_0/dxi  = B_0-A_0,
%       dX_1/dxi  = B_1-A_1,
%       dX_0/deta = C_0-A_0,
%       dX_1/deta = C_1-A_1,
%
% i due vettori tangenti sono semplicemente interpolazioni lineari dei
% corrispondenti lati iniziale e finale.

% Vettori tangenti rispetto a xi, eta
dXYZTW_dxi  = (1.0 - T_ref) .* lato1_base + T_ref .* lato1_top;
dXYZTW_deta = (1.0 - T_ref) .* lato2_base + T_ref .* lato2_top;

% Calcoliamo ora il vettore tangente associato alla variabile temporale. 
% Derivando
%
%       X(xi,eta,tau) = (1-tau) X_0(xi,eta) + tau X_1(xi,eta)
%
% rispetto a tau si ottiene
%
%       dX/dtau = X_1(xi,eta) - X_0(xi,eta).
%
% A differenza dei due vettori precedenti, questo vettore dipende
% esplicitamente da xi ed eta attraverso X_0 e X_1.
%
% Geometricamente rappresenta la direzione di moto del punto della
% superficie identificato dalla coppia (xi,eta).
%
% Se tutti i punti si muovessero con la stessa traslazione, questo
% vettore sarebbe costante.
%
% In presenza di rotazione, deformazione o moto non uniforme, esso
% dipende invece dalla posizione sul triangolo.

% Vettore tangente rispetto al tempo
dXYZTW_dt   = nodi_top_4D - nodi_base_4D;

% In R4, l'iperfaccia e' tridimensionale e possiede tre vettori tangenti.
% La normale generalizzata e' il vettore ortogonale a tutti e tre i vettori
% tangenti. Le seu componenti possono essere costruite tramite i minori 3x3
% della matrice jacobiana delle tangenti. I vertici della configurazione
% iniziale hanno tau =0, mentre quelli della configurazione finale hanno
% tau = 1. Di conseguenza, per ogni punto dell'iperfaccia X(xi, eta, tau) =
% (x(xi, eta, tau), y(xi, eta, tau), z(xi, eta, tau), tau). Quindi l'ultima
% coordinata della mappa e' il parametro tau. Pertanto 
%
%       dtau/dxi = 0, dtau/deta = 0, dtau/dtau = 1
%
% Scrivendo quindi X_xi = (x_xi, y_xi, z_xi, 0)^T ecc per X_eta, X_tau, la
% prima componente della normale generalizzata e' il minore ottenuto
% eliminando la prima coordinata, ossia 
%
%            ( y_xi,    z_xi,   0)
%   N_x = det( y_eta,   z_eta,  0)
%            ( y_tau,   z_tau,  1)
%
% Da cui otteniamo N_x = y_xi z_eta - z_xi y_eta
%
% La formula data sopra di N_x misura quanto l'elemento tridimensionale
% dell'iperfaccia contribuisce alla sua proiezione nella direzione x.

T_ref = T_ref(:);
W_ref = W_ref(:);

% Derivate rispetto a xi
dy_dxi = (1.0 - T_ref) .* lato1_base(2) + ...
         T_ref .* lato1_top(2);

dz_dxi = (1.0 - T_ref) .* lato1_base(3) + ...
         T_ref .* lato1_top(3);

% Derivate rispetto a eta
dy_deta = (1.0 - T_ref) .* lato2_base(2) + ...
          T_ref .* lato2_top(2);

dz_deta = (1.0 - T_ref) .* lato2_base(3) + ...
          T_ref .* lato2_top(3);

% Prima componente della normale generalizzata
Nx = dy_dxi .* dz_deta - dz_dxi .* dy_deta;

% Costruzione dei pesi orientati nella direzione x
WV_X = W_ref .* Nx;

end