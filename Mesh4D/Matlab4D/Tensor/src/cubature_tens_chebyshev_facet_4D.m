function moms_facet_raw = cubature_tens_chebyshev_facet_4D(XYZTW, WV_CUB, chebyshev_indices, dbox)

%**************************************************************************
%
% function moms_facet_raw = cubature_tens_chebyshev_facet_4D(XYZTW, WV_CUB, 
%                                                   chebyshev_indices, dbox)
%
% Momenti superficiali di Chebyshev sulla singola iperfaccia di un
% poliedro quadridimensionale.
% Per ogni elemento della base tensoriale viene costruita analiticamente la
% primitiva rispetto a x. 
% La costruzione dei polinomi e' completamente vettorializzata rispetto
% ai punti di quadratura e agli indici della base. Nella direzione x viene 
% calcolata una funzione di grado aggiuntivo per costruire analiticamente 
% la primitiva dei polinomi di Chebyshev.
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
%**************************************************************************
%
% OUTPUT:
%   moms_facet_raw     Vettore [M x 1] contenente i contributi superficiali
%                      ai momenti di Chebyshev.
%
%**************************************************************************

% Numero di punti di quadratura e numero di momenti richiesti
num_pts = size(XYZTW, 1);
num_moments = size(chebyshev_indices, 1);

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
% primitiva per i >= 2 contiene infatti il polinomio T_(i+1)(x)
% Di conseguenza, per costruire la primitiva del termine di grado massimo 
% T_max_i e' necessario conoscere anche T_(max_i+1)
% Nelle altre tre direzioni non e' necessaria alcuna primitiva e quindi
% e' sufficiente arrivare al grado massimo richiesto

x  = X_ref(:,1);
TX = chebpolys(max_i + 1, x);
TY = chebpolys(max_j, X_ref(:,2));
TZ = chebpolys(max_k, X_ref(:,3));
TT = chebpolys(max_l, X_ref(:,4));

PhiX = zeros(num_pts, num_moments);

% Gli indici nella direzione x vengono separati in tre casi. Questo
% permette di evitare integrazione numerica lungo x

is0 = (chebyshev_indices(:,1) == 0);
is1 = (chebyshev_indices(:,1) == 1);
isN = (chebyshev_indices(:,1) >= 2);

PhiX(:,is0) = repmat(x, 1, nnz(is0));
PhiX(:,is1) = repmat(0.5 .* x.^2, 1, nnz(is1));

if any(isN)
    i = chebyshev_indices(isN,1).';

PhiX(:,isN) = ...
    TX(:,i + 2) .* (i ./ (i.^2 - 1)) - x .* TX(:,i + 1) ./ (i - 1);
end

% Costruzione della funzione primitiva in modo tensoriale. Per ogni
% multi-indice si costruisce Phi_x = [int T_i(x) dx] T_j(y) T_k(z) T_l(tau).
% Il vettore WV_CUB contiene gia' il peso della quadratura e il fattore
% geometrico orientato n_x dS. Il prodotto elemento per elemento con
% PhiX restituisce l'integrando della formula di superficie

P = PhiX .* WV_CUB(:);
P = P .* TY(:,idx_j);
P = P .* TZ(:,idx_k);
P = P .* TT(:,idx_l);

% La trasformazione dalla coordinata fisica x alla coordinata di
% riferimento xi introduce il fattore dx = (xmax - xmin)/2 dxi

B1 = 0.5 * (dbox(2,1) - dbox(1,1));

moms_facet_raw = B1 .* sum(P, 1).';

end