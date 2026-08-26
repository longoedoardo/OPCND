function chebyshev_moms = cubature_tens_chebyshev_facet_V(nodes, weights, chebyshev_indices, dbox)
%CUBATURE_TENS_CHEBYSHEV_FACET_V  Momenti superficiali per la base tensoriale
%   di Chebyshev, sui nodi di quadratura di una singola faccia triangolare.
%
%--------------------------------------------------------------------------
% Oggetto:
%   Funzione ausiliaria di CHEBYSHEV_MOMENTS_POLYHEDRON. Dato un insieme di
%   nodi e pesi di quadratura su una faccia triangolare del poliedro,
%   calcola il contributo superficiale
%
%       Q_{ijk} = sum_q  w_q * Int_i(x_q) * T_j(y_q) * T_k(z_q)
%
%   per ogni tripla (i,j,k) elencata in chebyshev_indices, dove:
%     - T_p denota il polinomio di Chebyshev di grado p sull'intervallo
%       di riferimento [-1,1];
%     - Int_i(x) = Int_{-1}^{x} T_i(t) dt e' la primitiva di T_i rispetto
%       alla variabile x, necessaria per l'applicazione del Teorema della
%       Divergenza (campo vettoriale F = (Q,0,0) con dQ/dx = T_i*T_j*T_k).
%
%   Il risultato va poi moltiplicato, a livello della funzione chiamante,
%   per la componente x della normale esterna alla faccia, e sommato su
%   tutte le facce del poliedro per ottenere il momento di volume.
%
%   Questa versione e' ottimizzata per essere chiamata migliaia di volte
%   (una per faccia) su mesh grandi: rispetto a un'implementazione naive,
%   evita matrici temporanee superflue e trasposizioni esplicite, pur
%   mantenendo esattamente la stessa aritmetica (stesso ordine delle
%   operazioni in virgola mobile) della formulazione originale.
%
%--------------------------------------------------------------------------
% Input:
%   nodes             Matrice [n x 3] delle coordinate fisiche (x,y,z) dei
%                     nodi di quadratura sulla faccia.
%   weights           Vettore [n x 1] (o [1 x n]) dei pesi di quadratura
%                     fisici associati ai nodi.
%   chebyshev_indices Matrice [num_indici x 3] delle triple (i,j,k) della
%                     base tensoriale T_i(x)*T_j(y)*T_k(z).
%   dbox              Iper-rettangolo [xmin xmax ymin ymax zmin zmax] che
%                     definisce il dominio di normalizzazione per la base
%                     di Chebyshev (mappa affine verso [-1,1]^3).
%
% Output:
%   chebyshev_moms    Vettore colonna [num_indici x 1] dei contributi
%                     superficiali, uno per ogni tripla (i,j,k), nello
%                     stesso ordine di chebyshev_indices.
%
%--------------------------------------------------------------------------
% Note implementative:
%   - I polinomi T_i(x) e le rispettive primitive Int_i(x) vengono
%     costruiti in un'unica ricorrenza, senza calcolare o conservare
%     colonne di T_i che non servono alla formula della primitiva.
%   - Le colonne relative a T_j(y) e T_k(z) sono ottenute tramite la
%     funzione locale CHEBPOLYS (ricorrenza standard a tre termini).
%   - Il prodotto finale pesato viene calcolato applicando i pesi a una
%     sola delle tre matrici coinvolte e sommando poi per colonna, il che
%     evita sia matrici intermedie ridondanti sia una trasposizione
%     esplicita di una matrice di grandi dimensioni.
%
%--------------------------------------------------------------------------

%--- Cambio di variabile verso il cubo di riferimento [-1,1]^3 -----------
inv_dx = 2.0 / (dbox(2) - dbox(1));
inv_dy = 2.0 / (dbox(4) - dbox(3));
inv_dz = 2.0 / (dbox(6) - dbox(5));

mid_x = (dbox(1) + dbox(2)) * 0.5;
mid_y = (dbox(3) + dbox(4)) * 0.5;
mid_z = (dbox(5) + dbox(6)) * 0.5;

XN = (nodes(:,1) - mid_x) * inv_dx;
YN = (nodes(:,2) - mid_y) * inv_dy;
ZN = (nodes(:,3) - mid_z) * inv_dz;

B1 = 1.0 / inv_dx;   % fattore di scala Jacobiano per la sola direzione x
n  = size(nodes, 1);

% Gradi massimi richiesti per asse e indicizzazione 1-based 
max_i = max(chebyshev_indices(:,1));
max_j = max(chebyshev_indices(:,2));
max_k = max(chebyshev_indices(:,3));

idx_i = chebyshev_indices(:,1) + 1;
idx_j = chebyshev_indices(:,2) + 1;
idx_k = chebyshev_indices(:,3) + 1;

w = weights(:);

%--- Costruzione delle primitive Int_i(x), i = 0,...,max_i ---------------
%
% Formule utilizzate:
%   Int_0(x) = x
%   Int_1(x) = x^2 / 2
%   Int_i(x) = i*T_{i+1}(x) / (i^2-1)  -  x*T_i(x) / (i-1),   per i >= 2
%
% La ricorrenza a tre termini di Chebyshev, T_p = 2x*T_{p-1} - T_{p-2},
% viene applicata tenendo in memoria solo gli ultimi due
% polinomi necessari, senza costruire l'intera matrice T_0,...,T_{max_i+1}.

IntX = zeros(n, max_i + 1);
IntX(:,1) = XN;                       % i = 0

if max_i >= 1
    IntX(:,2) = 0.5 * (XN .* XN);      % i = 1
end

if max_i >= 2
    Tim1 = XN;             % T_1
    Tim2 = ones(n,1);      % T_0
    for p = 2:max_i
        Tp   = 2*XN.*Tim1 - Tim2;      % T_p
        Tpp1 = 2*XN.*Tp   - Tim1;      % T_{p+1}
        IntX(:,p+1) = (p * Tpp1) / (p^2 - 1) - (XN .* Tp) / (p - 1);
        Tim2 = Tim1;
        Tim1 = Tp;
    end
end

%--- Polinomi di Chebyshev lungo y e z (nessuna primitiva richiesta) -----
TY = chebpolys(max_j, YN);
TZ = chebpolys(max_k, ZN);

%--- Assemblaggio pesato finale -------------------------------------------
%
% chebyshev_moms(m) = sum_q  w(q) * Int_{i_m}(x_q) * T_{j_m}(y_q) * T_{k_m}(z_q)
%
% Il peso w viene applicato una sola volta (su IntX_cols), poi il prodotto
% elementwise con le altre due matrici viene sommato per colonna: questo
% evita sia la costruzione di tre matrici intermedie complete, sia una
% trasposizione esplicita per il prodotto matrice-vettore.

IntX_cols_w = IntX(:, idx_i) .* w;
TY_cols     = TY(:, idx_j);
TZ_cols     = TZ(:, idx_k);

chebyshev_moms = B1 * sum(IntX_cols_w .* TY_cols .* TZ_cols, 1).';

end


function T = chebpolys(deg, x)
%CHEBPOLYS  Matrice di Vandermonde-Chebyshev sulla retta reale, per ricorrenza.
%--------------------------------------------------------------------------
% Oggetto:
%   Calcola la matrice di Vandermonde-Chebyshev per ricorrenza.
%--------------------------------------------------------------------------
% Input:
%   deg: grado polinomiale massimo.
%   x:   vettore colonna di ascisse.
%--------------------------------------------------------------------------
% Output:
%   T: matrice di Vandermonde-Chebyshev in x, T(i,j+1) = T_j(x_i),
%      j = 0,...,deg.
%--------------------------------------------------------------------------
% I polinomi sono generati dalla ricorrenza a tre termini
%   T_0 = 1,  T_1 = x,  T_p = 2x*T_{p-1} - T_{p-2},
% e sono ortogonali rispetto al peso 1/sqrt(1-x^2) sull'intervallo [-1,1].
%--------------------------------------------------------------------------
% Autori originali:
%   Alvise Sommariva e Marco Vianello, Universita' di Padova, 15/12/2017
%--------------------------------------------------------------------------

n = length(x);
T = zeros(n, deg+1);

T(:,1) = 1;
if deg == 0
    return;
end
T(:,2) = x;

for j = 2:deg
    T(:,j+1) = 2*x.*T(:,j) - T(:,j-1);
end

end