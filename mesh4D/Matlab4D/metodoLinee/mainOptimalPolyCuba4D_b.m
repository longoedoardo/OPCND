%**************************************************************************
%
%                         MAIN OptimalPolyCuba4D
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

N_cc = n_tau - 1;

theta      = pi * (0:N_cc)' / N_cc;
nodes_cheb = cos(theta);
tau_nodes  = flipud((nodes_cheb + 1) / 2);

j = (1:floor(N_cc/2))';           % indici j, vettore colonna
fact = 2 * ones(size(j));
if mod(N_cc,2) == 0
    fact(end) = 1;                 % j = N_cc/2 -> fattore 1
end

% Matrice (n_tau x length(j)): termine cos(2*j*theta_i)/(1-4*j^2)
% per ogni combinazione di nodo i e indice j
M = cos(2 * theta * j') ./ (1 - 4 * j'.^2);   % broadcasting: theta (n_tau x1), j' (1 x nj)
sum_w = M * fact;                              % somma su j per ogni i (n_tau x 1)

g_cc = ones(n_tau,1);
g_cc([1 end]) = 0.5;

tau_weights = (1/N_cc) * g_cc .* (1 + sum_w);
tau_weights = flipud(tau_weights);

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
    current_tau = tau_nodes(k); 

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
    I_4D = I_4D + sum((tau_weights(k) * W) .* feval);
end

fprintf('Integrale: %1.15e\n', I_4D);
fprintf('........................\n');