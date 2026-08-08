%**************************************************************************
%
%                         MAIN OptimalPolyCuba3D
%
%**************************************************************************
%
% Metodo che calcola l'integrale di una funzione f(x,y,z) su un poliedro 
% definito, usando grado algebrico di esattezza "ade".
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
fprintf('Cubatura con OptimalPolyCuba3D \n');
fprintf('........................\n');

ade = 1; 
fprintf('ade: %-3.0f\n', ade);

vertices = load('vertex.dat');
facets   = load('tri.dat');
f = @(x,y,z) ones(size(x));

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

fXYZW = feval(f, XYZW(:,1), XYZW(:,2), XYZW(:,3));
W = XYZW(:,4);
I = W' * fXYZW;

fprintf('Integrale: %1.15e\n', I);
fprintf('........................\n');