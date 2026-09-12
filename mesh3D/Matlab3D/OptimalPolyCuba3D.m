function [XYZ, W] = OptimalPolyCuba3D(ade, vertices, facets)

%**************************************************************************
%
% function I = OptimalPolyCuba3D(ade, vertices, facets, f)
%
% Calcola un'approssimazione numerica di un integrale di volume
% tridimensionale su un dominio poliedrico rappresentato mediante
% una mesh superficiale triangolare chiusa.
%
% La regola di cubatura è costruita combinando una base di Chebyshev
% tensoriale shape-independent con i momenti shape-dependent calcolati
% a partire dalla rappresentazione superficiale del dominio poliedrico.
%
%**************************************************************************
%
%   INPUTS:
%   - ade:
%       Grado polinomiale totale massimo della regola di cubatura.
%
%   - vertices:
%       Matrice N x 3 contenente le coordinate cartesiane dei vertici
%       della mesh. Ogni riga rappresenta un vertice [x_i, y_i, z_i].
%
%   - facets:
%       Matrice M x 3 contenente la connettività della mesh superficiale
%       triangolare. Ogni riga contiene gli indici dei tre vertici che
%       definiscono una faccia triangolare.
%
%   - f:
%       Function handle che rappresenta la funzione integranda f(x,y,z).
%
%   OUTPUT:
%   - I:
%       Approssimazione numerica dell'integrale di volume
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
projectRoot = fileparts(mfilename('fullpath'));

% Aggiunta della cartella contenente le funzioni ausiliarie
addpath(fullfile(projectRoot, 'src'));

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

bbox_min = min(vertices, [], 1);
bbox_max = max(vertices, [], 1);
bbox = [bbox_min; bbox_max]; 

% Momenti della base sul poliedro (teorema della divergenza sulle facce)
moments_ch = chebyshev_moments_polyhedron(vertices, facets, ade, ...
    chebyshev_indices, bbox);

% Riscalamento della griglia di riferimento sulla bounding box reale
XYZW_tens = scale_rule(XYZW_tens_ref, bbox);

% Pesi di cubatura ottenuti ricombinando la regola di riferimento
W = (XYZW_tens(:,4)) .* V_ref * (moments_ch ./ coeffs);
XYZW = [XYZW_tens(:,1:3) W];

XYZ = XYZW(:,1:3);
W = XYZW(:,4);