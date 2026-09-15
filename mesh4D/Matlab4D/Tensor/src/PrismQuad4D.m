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

% Indichiamo con
%
%       c = dX/dxi,
%       r = dX/deta,
%       z = dX/dtau.
%
% Ciascuno di essi e' una matrice [Nq x 4], poiche' esiste un vettore
% tangente distinto in corrispondenza di ciascun nodo di quadratura.
%
% Successivamente le quattro componenti vengono estratte separatamente
% per poter costruire i minori 3 x 3 necessari alla normale generalizzata.

c = dXYZTW_dxi;
r = dXYZTW_deta;
z = dXYZTW_dt;

c1 = c(:,1); c2 = c(:,2); c3 = c(:,3); c4 = c(:,4);
r1 = r(:,1); r2 = r(:,2); r3 = r(:,3); r4 = r(:,4);
z1 = z(:,1); z2 = z(:,2); z3 = z(:,3); z4 = z(:,4);

% Costruzione del vettore normale generalizzato in R^4. In R^4 una 
% ipersuperficie tridimensionale possiede tre vettori tangenti:
%
%       c = dX/dxi,
%       r = dX/deta,
%       z = dX/dtau.
%
% Occorre quindi costruire un vettore in R^4 ortogonale simultaneamente
% a tutti e tre.
%
% La generalizzazione del prodotto vettoriale e' ottenuta mediante i
% cofattori della matrice
%
%       
%
%                    J = [ c1 c2 c3 c4
%                         r1 r2 r3 r4
%                         z1 z2 z3 z4 ].
%
% Per ogni coordinata si elimina la colonna corrispondente e si calcola
% il determinante della matrice 3 x 3 rimanente.
%
% Definiamo quindi
%
%       M1 = det([c2 c3 c4;
%                 r2 r3 r4;
%                 z2 z3 z4]),
%
%       M2 = det([c1 c3 c4;
%                 r1 r3 r4;
%                 z1 z3 z4]),
%
%       M3 = det([c1 c2 c4;
%                 r1 r2 r4;
%                 z1 z2 z4]),
%
%       M4 = det([c1 c2 c3;
%                 r1 r2 r3;
%                 z1 z2 z3]).
%
% Il vettore dei cofattori e'
%
%       N = [M1, -M2, M3, -M4].
%
% I segni alternati derivano dalla struttura dei cofattori (+,-,+,-)
% associati allo sviluppo di un determinante lungo una riga o una
% colonna. Per costruzione si ha
%
%       N . c = 0,
%       N . r = 0,
%       N . z = 0.
%
% Quindi N e' ortogonale allo spazio tangente tridimensionale
% dell'iperfaccia. Inoltre la sua norma incorpora il fattore metrico 
% associato alla parametrizzazione N = n dS

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

Nvec = [M1, -M2, M3, -M4];


% Il vettore N costruito mediante i cofattori determina la direzione
% normale, ma il suo verso dipende dall'ordinamento dei vertici e dalla
% parametrizzazione scelta. Per il Teorema della Divergenza e' invece 
% necessario utilizzare la normale uscente. Per ogni punto di quadratura 
% consideriamo quindi il vettore
%
%       d = X - C,
%
% dove C e' un punto interno al dominio, qui rappresentato dal
% baricentro. Se la normale N punta verso l'esterno, ci aspettiamo che il 
% prodotto sia > 0. Se invece e' < 0, la normale punta verso l'interno e 
% deve essere invertita. 

direzione_uscente = XYZTW - baricentro;
flip = sum(Nvec .* direzione_uscente, 2) < 0;
Nvec(flip, :) = -Nvec(flip, :);

% Costruzione dei pesi orientati nella direzione x. Nel calcolo dei momenti 
% di Chebyshev 4D viene utilizzato il campo F = (Phi_x, 0, 0, 0). Il 
% prodotto scalare con l'elemento di superficie orientato e'
%
%       F * N = Phi_x N_x.
%
% Pertanto, sulla singola iperfaccia, la componente della normale che
% interessa il calcolo dei momenti e' esclusivamente la prima:
%
%       N_x = Nvec(:,1).
%
% I pesi della quadratura sul prisma di riferimento vengono moltiplicati
% per questa componente, ottenendo direttamente
%
%       WV_X = W_ref .* N_x.

WV_X = W_ref .* Nvec(:, 1);

end