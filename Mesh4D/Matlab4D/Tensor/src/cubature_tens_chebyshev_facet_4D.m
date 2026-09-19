function moms_facet_raw = cubature_tens_chebyshev_facet_4D(XYZTW, WV_CUB, chebyshev_indices, dbox)

%**************************************************************************
%
% function moms_facet_raw = cubature_tens_chebyshev_facet_4D(XYZTW, WV_CUB, 
%                                                   chebyshev_indices, dbox)
%
% Momenti superficiali di Chebyshev sulla singola iperfaccia di un
% poliedro quadridimensionale.
%
% Dato un insieme di nodi e pesi di cubatura su una singola iperfaccia
% tridimensionale di un poliedro 4D, calcola il contributo superficiale
% ai momenti della base tensoriale di Chebyshev
%
%       T_i(x) * T_j(y) * T_k(z) * T_l(t).
%
% La primitiva viene calcolata analiticamente lungo la prima coordinata x.
% Questo permette di applicare il Teorema della Divergenza scegliendo
% un campo vettoriale con derivata rispetto a x uguale al prodotto
% tensoriale dei polinomi di Chebyshev.
%
% La costruzione e' completamente vettorializzata rispetto ai punti di
% quadratura e agli indici della base, evitando un ciclo sui momenti.
% In questo modo il prodotto tensoriale dei polinomi viene realizzato
% mediante operazioni matriciali e prodotti elemento-per-elemento.
%
%**************************************************************************
%
% INPUT:
%   XYZTW              Matrice [N x 4] dei nodi di quadratura
%                      [x, y, z, tau].
%
%   WV_CUB             Vettore [N x 1] dei pesi di quadratura comprensivi
%                      del fattore geometrico e della componente orientata
%                      della normale nella direzione x.
%
%   chebyshev_indices  Matrice [M x 4] contenente gli indici
%                      (i,j,k,l) della base tensoriale di Chebyshev.
%
%   dbox               Bounding box 4D nel formato
%                      [xmin ymin zmin taumin;
%                       xmax ymax zmax taumax].
%
% OUTPUT:
%   moms_facet_raw     Vettore [M x 1] contenente i contributi superficiali
%                      ai momenti di Chebyshev.
%
%**************************************************************************

% Numero di punti di quadratura e numero di momenti richiesti
num_pts = size(XYZTW, 1);

% Trasformazione affine dal dominio fisico al dominio di riferimento [-1,1]^4
X_ref = (2.0 * XYZTW - (dbox(1,:) + dbox(2,:))) ./ (dbox(2,:) - dbox(1,:));


% Ogni riga di chebyshev_indices identifica un elemento della base T_i(x)
% T_j(y) T_k(z) T_l(tau). Per costruire le matrici dei polinomi e'
% sufficiente conoscere il massimo grado richiesto separatamente in
% ciascuna direzione. Questa operazione permette di evitare la costruzione
% di polinomi di grado superiore a quello realmente necessario

% Gradi massimi richiesti nelle quattro direzioni
max_i = max(chebyshev_indices(:,1));
max_j = max(chebyshev_indices(:,2));
max_k = max(chebyshev_indices(:,3));
max_l = max(chebyshev_indices(:,4));

% Poiche' gli indici dei polinomi di Chebyshev partono da 0, l'indice
% matematico i corrpisponde alla colonna i+1 della matrice restituita da
% chebpolys, in MATLAB gli indici di array partono da 1

idx_i = chebyshev_indices(:,1) + 1;
idx_j = chebyshev_indices(:,2) + 1;
idx_k = chebyshev_indices(:,3) + 1;
idx_l = chebyshev_indices(:,4) + 1;

% Calcolo delle matrici di Chebyshev
%
% La direzione x richiede un grado aggiuntivo. Per la componente x del 
% campo vettoriale occorre una primitiva di T_i(x). La formula della 
% primitiva per i >= 2 contiene infatti il polinomio T_(i+1)(x). 
% Di conseguenza, per costruire la primitiva del termine di grado massimo 
% T_max_i e' necessario conoscere anche T_(max_i+1).

TX = chebpolys(max_i + 1, X_ref(:,1));
TY = chebpolys(max_j,     X_ref(:,2));
TZ = chebpolys(max_k,     X_ref(:,3));
TT = chebpolys(max_l,     X_ref(:,4));

% Costruzione vettorializzata delle primitive nella direzione x

IntX = zeros(num_pts, max_i + 1);

% Primitiva di T_0(x) = 1
IntX(:,1) = X_ref(:,1);

% Primitiva di T_1(x) = x
if max_i >= 1
    IntX(:,2) = 0.5 * (X_ref(:,1) .* X_ref(:,1));
end

% Primitiva di T_i(x), i >= 2
if max_i >= 2

    i_vec = 2:max_i;

    IntX(:,i_vec + 1) = ...
        i_vec .* TX(:,i_vec + 2) ./ (i_vec.^2 - 1) - ...
        X_ref(:,1) .* TX(:,i_vec + 1) ./ (i_vec - 1);

end

% Fattore di scala associato alla trasformazione della coordinata x

% Nel caso 3D si definiva inv_dx = 2/(x_max-x_min),
% e quindi ottenevamo B1 = 1/inv_dx.
%
% La stessa struttura viene mantenuta qui per garantire coerenza tra
% l'implementazione 3D e quella 4D.

inv_dx = 2.0 / (dbox(2) - dbox(1));
B1 = 1.0 / inv_dx;

% Estrazione delle colonne richieste dalla base di Chebyshev
IntX_cols = IntX(:,idx_i);
TY_cols   = TY(:,idx_j);
TZ_cols   = TZ(:,idx_k);
TT_cols   = TT(:,idx_l);

% Per ogni punto di quadratura e per ogni indice (i,j,k,l) si costruisce
% il prodotto
%
%       Phi_i(x)* T_j(y)* T_k(z)* T_l(tau).
%
% Le quattro matrici [N x M] vengono quindi moltiplicate elemento per
% elemento.
%
% Il risultato e' una matrice [N x M], nella quale:
%
%   - ogni riga corrisponde a un punto di quadratura;
%   - ogni colonna corrisponde a un momento della base.
%
% La trasposizione produce una matrice [M x N], che viene poi moltiplicata
% per il vettore dei pesi WV_CUB [N x 1].
%
% Si ottiene quindi direttamente un vettore [M x 1]:
%
%       moms_facet_raw(m)
%
% rappresenta il contributo della singola iperfaccia al momento
% corrispondente alla m-esima riga di chebyshev_indices.
%
% Il fattore B1 converte infine la primitiva dalla coordinata di
% riferimento alla coordinata fisica.
%
% La struttura e' esattamente analoga al caso 3D:
%
%   3D:
%       (IntX_cols .* TY_cols .* TZ_cols)' * w
%
%   4D:
%       (IntX_cols .* TY_cols .* TZ_cols .* TT_cols)' * WV_CUB
%
% L'unica differenza e' quindi l'aggiunta del quarto fattore
% T_l(tau), corrispondente alla nuova coordinata spazio-temporale.

moms_facet_raw = B1*((IntX_cols .* TY_cols .* TZ_cols .* TT_cols)' * WV_CUB(:));

end