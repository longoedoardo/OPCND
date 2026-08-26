function moms_facet_raw = cubature_tens_chebyshev_facet_4D(XYZTW, WV_CUB, chebyshev_indices, dbox)
%**************************************************************************
%
%   function moms_facet_raw = cubature_tens_chebyshev_facet_4D(XYZTW, 
%           WV_CUB, chebyshev_indices, dbox)
%
%   Momenti superficiali di Chebyshev
%   sulla singola iperfaccia di un poliedro quadridimensionale.
%
%   Dato un insieme di nodi e pesi di cubatura su una singola iperfaccia
%   tridimensionale di un poliedro 4D, calcola il contributo superficiale
%   ai momenti della base tensoriale di Chebyshev
%
%       T_i(x) * T_j(y) * T_k(z) * T_l(t)
%
%   per ogni quadrupla (i,j,k,l) contenuta in chebyshev_indices.
%
%   La primitiva viene calcolata analiticamente lungo la prima
%   coordinata x. Questo permette di applicare il Teorema della
%   Divergenza scegliendo un campo vettoriale con derivata rispetto
%   a x uguale al prodotto tensoriale di Chebyshev.
%
%   Il risultato contiene quindi i contributi superficiali associati
%   alla normale nella direzione x; il fattore relativo alla normale
%   e all'orientazione della faccia viene gestito dalla funzione
%   chiamante.
%
%   Rispetto a una formulazione che costruisce esplicitamente la matrice
%   completa dei prodotti tensoriali, questa versione evita di creare
%   la matrice temporanea V_poly, che può diventare molto grande quando
%   il numero di nodi o il numero di momenti è elevato.
%
%**************************************************************************

% Numero di punti di quadratura sulla faccia e numero di momenti richiesti.
num_pts = size(XYZTW, 1);
num_indici = size(chebyshev_indices, 1);

% Calcoliamo centro e ampiezza di ciascuna direzione per effettuare
% la trasformazione affine dal dominio fisico al cubo di riferimento
% [-1,1]^4.

box_sum  = dbox(1,:) + dbox(2,:);
box_diff = dbox(2,:) - dbox(1,:);

% Trasformazione affine delle coordinate fisiche: x_ref = (2*x - (xmin+xmax)) / (xmax-xmin).

X_ref = (2.0 * XYZTW - box_sum) ./ box_diff;

% Gradi massimi richiesti nelle quattro direzioni.
% Gli indici contenuti in chebyshev_indices partono da 0, mentre
% le colonne delle matrici MATLAB partono da 1.

max_i = max(chebyshev_indices(:,1));
max_j = max(chebyshev_indices(:,2));
max_k = max(chebyshev_indices(:,3));
max_l = max(chebyshev_indices(:,4));

idx_i = chebyshev_indices(:,1) + 1;
idx_j = chebyshev_indices(:,2) + 1;
idx_k = chebyshev_indices(:,3) + 1;
idx_l = chebyshev_indices(:,4) + 1;

% Per la direzione x serve anche il polinomio T_(i+1), perché
% la primitiva di T_i contiene T_(i+1).
%
% Nelle altre tre direzioni sono sufficienti i gradi effettivamente
% richiesti dalla base tensoriale.

Tx = chebpolys(max_i + 1, X_ref(:,1));
Ty = chebpolys(max_j,     X_ref(:,2));
Tz = chebpolys(max_k,     X_ref(:,3));
Tt = chebpolys(max_l,     X_ref(:,4));

% Costruzione delle primitive dei polinomi di Chebyshev nella
% direzione x.
%
% Si utilizzano le formule
%
%   Int T_0(x) dx = x,
%
%   Int T_1(x) dx = x^2/2,
%
% e, per i >= 2,
%
%   Int T_i(x) dx =
%       i/(i^2-1) * T_(i+1)(x)
%       - x/(i-1) * T_i(x).

IntTx = zeros(num_pts, max_i + 1);

% Primitiva di T_0(x) = 1.

IntTx(:,1) = X_ref(:,1);

% Primitiva di T_1(x) = x.

if max_i >= 1
    IntTx(:,2) = 0.5 * (X_ref(:,1) .* X_ref(:,1));
end

% Per i >= 2 utilizziamo direttamente i polinomi di Chebyshev
% già calcolati, evitando di ricostruirli.

if max_i >= 2
    for i = 2:max_i

        IntTx(:,i+1) = ...
            (i * Tx(:,i+2)) / (i^2 - 1) ...
            - (X_ref(:,1) .* Tx(:,i+1)) / (i - 1);

    end
end


% La trasformazione dalla coordinata fisica x alla coordinata
% normalizzata x_ref introduce il fattore
%
%       dx = (xmax-xmin)/2 * dx_ref.
%
% Per questo motivo la primitiva rispetto alla coordinata fisica
% deve essere moltiplicata per la semiampiezza della bounding box
% nella direzione x.

scale_x = box_diff(1) / 2;
IntTx = IntTx * scale_x;

% I pesi di cubatura vengono applicati direttamente alla matrice
% delle primitive.

IntTx_weighted = IntTx .* WV_CUB(:);

% Calcolo dei momenti.
%
% Per ogni quadrupla (i,j,k,l) dobbiamo calcolare
%
%   sum_q WV_CUB(q) * IntT_i(x_q) * T_j(y_q) * T_k(z_q) * T_l(t_q).

moms_facet_raw = zeros(num_indici, 1);

for ii = 1:num_indici

    % Selezioniamo il polinomio richiesto in ciascuna direzione
    % dalla quadrupla corrente (i,j,k,l).

    IntTx_i = IntTx_weighted(:,idx_i(ii));
    Ty_j    = Ty(:,idx_j(ii));
    Tz_k    = Tz(:,idx_k(ii));
    Tt_l    = Tt(:,idx_l(ii));

    % Prodotto tensoriale valutato nei nodi di quadratura e
    % somma sui punti della faccia.

    moms_facet_raw(ii) = sum(IntTx_i .* Ty_j .* Tz_k .* Tt_l);

end

end

function T=chebpolys(deg,x)

%--------------------------------------------------------------------------
% Object:
% This routine computes the Chebyshev-Vandermonde matrix on the real line
% by recurrence.
%--------------------------------------------------------------------------
% Input:
% deg: maximum polynomial degree
% x: 1-column array of abscissas
%--------------------------------------------------------------------------
% Output:
% T: Chebyshev-Vandermonde matrix at x, T(i,j+1)=T_j(x_i), j=0,...,deg.
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

T = zeros(length(x), deg+1);
T(:,1) = 1;

if deg >= 1
    T(:,2) = x;
end

for j = 2:deg
    T(:,j+1) = 2*x.*T(:,j) - T(:,j-1);
end
end