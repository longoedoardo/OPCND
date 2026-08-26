%**************************************************************************
%
%                         MAIN OptimalPolyCuba4D_MoL
%
%**************************************************************************
%
% Metodo che calcola l'integrale di una funzione f(x,y,z) su un poliedro 
% nel tempo [0,1] con il metodo delle linee, usando grado algebrico di 
% esattezza "ade".
% 
% INPUT:
%   vertices : Matrice N x 3 dei vertici [x, y, z]
%   facets   : Struttura o cell array delle facce
%   ade      : Grado di esattezza algebrica
%   f        : Function handle della funzione integranda
%
%**************************************************************************

addpath('ShapeIndependent/');
addpath('ShapeDependent/');

fprintf('........................\n');
fprintf('Cubatura con OptimalPolyCuba4D (Metodo delle linee) \n');
fprintf('........................\n');

ade = 4; 
n_tau = 100;

fprintf('ade: %-3.0f\n', ade);

f = @(x,y,z,tau) x.*tau.*y.^3;

%**************************************************************************
%
% GEOMETRIA POLIEDRO INIZIALE E FINALE
%
%**************************************************************************

vertici_iniziali = load('vertex.dat');    % Configurazione a tau = 0
facets           = load('tri.dat');       % Connettivita' delle facce (costante)
vertici_finali   = load('vertex_new.dat');% Configurazione a tau = 1

%**************************************************************************
%
% DISCRETIZZAZIONE TEMPORALE
%
%**************************************************************************

% La formula di Clenshaw-Curtis utilizza i nodi di Chebyshev-Lobatto 
% 
% x_i = cos(i*pi/N)
%
% con i = 0,...,N, sull'intervallo [-1,1]. I nodi vengono successivamente 
% trasformati nell'intervallo [0,1].

% Se n_tau è il numero totale di nodi, allora il numero di
% sottointervalli è N = n_tau - 1.
num_intervalli = n_tau - 1;

theta = pi * (0:num_intervalli)' / num_intervalli;
nodiChebLob = cos(theta);


% Trasformazione affine da [-1,1] a [0,1]: map di x in (x+1)/2. Uso flipud 
% per riordinare i nodi in ordine crescente

tau_nodi = flipud((nodiChebLob + 1) / 2);

% Nella formula dei pesi compaiono i termini pari, cos(2*j*theta). Sono 
% quindi necessari gli indici j = 1,...,floor(N/2).

indici_serie = (1:floor(num_intervalli/2))';

% Coefficienti che moltiplicano i termini della serie.

coefficienti_serie = 2 * ones(size(indici_serie));

% Se N è pari, l'ultimo termine della serie corrisponde a j = N/2.
% In questo caso il suo coefficiente deve essere 1 anziché 2.

if mod(num_intervalli, 2) == 0
    coefficienti_serie(end) = 1;
end

% Costruiamo una matrice in cui ogni riga corrisponde a un nodo theta_i;
% e ogni colonna corrisponde a un indice j. L'elemento (i,j) è
% [cos(2*j*theta_i)] / [1 - 4*j^2] che compare nella formula dei 
% pesi di Clenshaw-Curtis.

termini_serie = cos(2 * theta * indici_serie')./ (1 - 4 * indici_serie'.^2);

% Sommiamo i termini della serie pesandoli con i coefficienti
% precedentemente definiti.

somma_serie = termini_serie * coefficienti_serie;

% Nella formula dei pesi, i nodi agli estremi theta = 0 e theta = pi
% hanno un fattore 1/2. Per tutti gli altri nodi il fattore vale 1.

fattore_estremi = ones(n_tau, 1);

fattore_estremi([1, end]) = 0.5;

% Formula dei pesi sull'intervallo [-1,1]: w_i = 1/N * fattore_i * (1 + somma_serie_i).
% Poiché i nodi sono stati trasformati nell'intervallo [0,1],
% la stessa formula viene utilizzata con la corrispondente
% normalizzazione.

tau_pesi = ...
    (1 / num_intervalli) ...
    * fattore_estremi ...
    .* (1 + somma_serie);

% I nodi di Chebyshev erano inizialmente ordinati da 1 a 0.
% Dopo flipud, i nodi tau_nodi sono invece ordinati da 0 a 1. 
% I pesi devono avere lo stesso ordinamento dei nodi.

tau_pesi = flipud(tau_pesi);

%**************************************************************************
%
% INIZIO PARTE SHAPE-INDEPENDENT
%
%**************************************************************************

% Griglia tensoriale di Gauss-Chebyshev sul cubo di riferimento [-1,1]^3
XYZW_tens_ref = cub_gausscheb_tens3D(2*ade);

% Dimensione dello spazio polinomiale di grado totale <= ade
N_mom = ((ade+3)*(ade+2)*(ade+1))/6;

% Indici (i,j,k) dei monomi tensoriali di Chebyshev, ordinamento GRLEX
chebyshev_indices = zeros(N_mom,3);
for i = 2:N_mom
    chebyshev_indices(i,:) = mono_next_grlex(3, chebyshev_indices(i-1,:));
end

% Matrice di Vandermonde-Chebyshev 3D sui nodi della griglia di riferimento
X = XYZW_tens_ref(:,1:3);

V_ref = dCHEBVAND(ade, X, chebyshev_indices);

% Coefficienti di normalizzazione della base di Chebyshev tensoriale
coeffs = tenscheb_norm2sq(chebyshev_indices);

%**************************************************************************
%
% INIZIO PARTE SHAPE-DEPENDENT
%
%**************************************************************************
I_4D = 0;
for k = 1:n_tau
    current_tau = tau_nodi(k); 

    vertices_tau = (1 - current_tau) * vertici_iniziali + current_tau * vertici_finali;

    % Bounding box dinamica al tempo corrente
    bbox_min = min(vertices_tau, [], 1);
    bbox_max = max(vertices_tau, [], 1);
    bbox_tau = [bbox_min(1), bbox_max(1), bbox_min(2), bbox_max(2), ...
            bbox_min(3), bbox_max(3)];

    % Momenti della base sul poliedro (teorema della divergenza sulle facce)
    moments_ch = chebyshev_moments_polyhedron(vertices_tau, facets, ade, ...
    chebyshev_indices, bbox_tau);

    % Riscalamento della griglia di riferimento sulla bounding box reale corrente
    XYZW_tens = scale_rule(XYZW_tens_ref, bbox_tau);

    % Pesi di cubatura ottenuti ricombinando la regola di riferimento
    W = (XYZW_tens(:,4)) .* V_ref * (moments_ch ./ coeffs);
    XYZW = [XYZW_tens(:,1:3) W];

    feval = f(XYZW_tens(:,1), XYZW_tens(:,2), XYZW_tens(:,3), current_tau * ones(size(XYZW_tens(:,1))));

    W = XYZW(:,4);
    I_4D = I_4D + sum((tau_pesi(k) * W) .* feval);
end

fprintf('Integrale: %1.15e\n', I_4D);
fprintf('........................\n');