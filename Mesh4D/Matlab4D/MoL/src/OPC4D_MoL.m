function [XYZtau, W] = OPC4D_MoL(ade, n_tau, vertici_iniziali, vertici_finali, facets, method)

%**************************************************************************
%
% function [XYZtau, W] = OptimalPolyCuba4D_MoL(ade, n_tau, ...
%                               vertici_iniziali, vertici_finali, facets)
%
% Costruisce una regola di cubatura quadridimensionale per un dominio
% poliedrico in movimento per t in [0,1], rappresentato nello spazio
% mediante una mesh superficiale triangolare chiusa.
%
% La regola di cubatura è costruita combinando una base di Chebyshev
% tensoriale shape-independent con i momenti shape-dependent calcolati
% a partire dalla rappresentazione superficiale del dominio poliedrico.
% La discretizzazione temporale è ottenuta mediante una formula di
% quadratura di Clenshaw-Curtis basata sui nodi di Chebyshev-Lobatto.
%
% La funzione restituisce esclusivamente i nodi e i pesi della regola
% di cubatura 4D, senza effettuare la valutazione della funzione integranda
% né il calcolo dell'integrale.
%
%**************************************************************************%
%   INPUTS:
%   - ade:
%       Grado polinomiale totale massimo della regola di cubatura spaziale.
%
%   - n_tau:
%       Numero di nodi della discretizzazione temporale ottenuta mediante
%       la formula di quadratura di Clenshaw-Curtis.
%
%   - vertici_iniziali:
%       Matrice N x 3 contenente le coordinate cartesiane dei vertici
%       della mesh nella configurazione iniziale. Ogni riga rappresenta
%       un vertice [x_i, y_i, z_i].
%
%   - vertici_finali:
%       Matrice N x 3 contenente le coordinate cartesiane dei vertici
%       della mesh nella configurazione finale. Ogni riga rappresenta
%       un vertice [x_i, y_i, z_i].
%
%   - facets:
%       Matrice M x 3 contenente la connettività della mesh superficiale
%       triangolare. Ogni riga contiene gli indici dei tre vertici che
%       definiscono una faccia triangolare.
%   - method:
%       Metodo di integrazione sulle facce triangolari:
%       * "D": Metodo simmetrico di Dunavant
%       * "GJ": Metodo di Gauss-Jacobi
%
%   OUTPUT:
%   - XYZtau:
%       Matrice (N x 4) contenente i nodi della regola di cubatura
%       quadridimensionale. Ogni riga ha la forma [x_i, y_i, z_i, tau_i],
%       dove tau_i appartiene all'intervallo temporale [0,1].
%
%   - W:
%       Vettore colonna (N x 1) contenente i pesi della regola di
%       cubatura quadridimensionale associati ai nodi in XYZtau.
%
%**************************************************************************
%
%   Autore:
%       Edoardo Longo, Università degli Studi di Verona
%
%   Data:
%       Settembre 2026
%
%**************************************************************************

% Percorso della directory principale del progetto
projectRoot = fileparts(fileparts(mfilename('fullpath')));

% Aggiunta della cartella contenente le funzioni ausiliarie
addpath(fullfile(projectRoot, 'src'));

% Controllo della disponibilità delle regole di Dunavant
if strcmp(method, "D") && ade > 20
    warning(['Le regole di Dunavant sono disponibili solo fino al ', ...
        'grado 20. Cambiare metodo o modificare ade']);
    return;
end


%**************************************************************************
%
% INIZIO PARTE SHAPE-INDEPENDENT
%
%**************************************************************************

% Discretizzazione temporale
[tau_nodi, tau_pesi] = ClenshawCurtisTime(n_tau);

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

% Numero di nodi della quadratura spaziale
N_spaz = size(XYZW_tens_ref, 1);

% Preallocazione della matrice dei nodi 4D e del vettore dei pesi
XYZtau = zeros(N_spaz * n_tau, 4);
W = zeros(N_spaz * n_tau, 1);

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
    chebyshev_indices, bbox_tau, method);

    % Riscalamento della griglia di riferimento sulla bounding box reale corrente
    XYZW_tens = scale_rule(XYZW_tens_ref, bbox_tau);

    % Pesi di cubatura ottenuti ricombinando la regola di riferimento, pesi
    % della quadratura spaziale
    
    W_spaziali = (XYZW_tens(:,4)) .* V_ref * (moments_ch ./ coeffs);

    % Indici dei nodi relativi all'istante temporale corrente
    indici = (k-1)*N_spaz + (1:N_spaz);

    % Nodi della quadratura 4D
    XYZtau(indici,1:3) = XYZW_tens(:,1:3);
    XYZtau(indici,4) = current_tau;

    % Pesi della quadratura 4D
    W(indici) = tau_pesi(k) * W_spaziali;
end