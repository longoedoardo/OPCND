
%**************************************************************************
%
%                               MAIN LALO 3D
%
%**************************************************************************
%
% Metodo che calcola l'integrale di una funzione f(x,y,z) su un poliedro 
% definito, usando grado algebrico di esattezza "ade".
% 
% INPUT:
%   vertices : Matrice Nv x 3 dei vertici [x, y, z]
%   facets   : Struttura o cell array delle facce (richiesto da chebyshev_moments_polyhedron)
%   ade      : Grado di esattezza algebrica
%   f        : Function handle della funzione integranda es. @(x,y,z) x.^3 + y + z
%
%**************************************************************************

addpath('GeometryCheap/');
addpath('PrepCheap/');
addpath('CubaturaFunzioniCheap/');

fprintf('........................\n');
fprintf('Cubatura con LALOcheap \n');
fprintf('........................\n');
fprintf('ade: %-3.0f\n', ade);

vertices = load('vertex.dat');
facets   = load('tri.dat');

ade = 4; 

f = @(x,y,z) x.^4 + y + z;

%**************************************************************************
% Fase indipendente dalla geometria: griglia e base polinomiale di
% riferimento, generate una sola volta in funzione del solo grado "ade"
%**************************************************************************
% Griglia tensoriale di Gauss-Chebyshev sul cubo di riferimento [-1,1]^3
XYZW_tens_ref = cub_gausscheb_tens3D(2*ade);

% Dimensione dello spazio polinomiale di grado totale <= ade
N_mom = ((ade+3)*(ade+2)*(ade+1))/6;
fprintf('Dimensione spazio polinomiale (N_mom): %-3.0f\n', N_mom);

% Indici (i,j,k) dei monomi tensoriali di Chebyshev, ordinamento GRLEX
chebyshev_indices = zeros(N_mom,3);
for i = 2:N_mom
    chebyshev_indices(i,:) = mono_next_grlex(3, chebyshev_indices(i-1,:));
end

% Matrice di Vandermonde-Chebyshev 3D sui nodi della griglia di riferimento
X = XYZW_tens_ref(:,1:3);
V_ref = dCHEBVAND(ade, X, chebyshev_indices);

%**************************************************************************
% Bounding box del poliedro (min/max per coordinata)
% IMPORTANTE: Formattato identico al blocco Fortran per evitare scambi di assi
%**************************************************************************
bbox_min = min(vertices, [], 1);
bbox_max = max(vertices, [], 1);
bbox = [bbox_min; bbox_max]; 

%**************************************************************************
% Calcolo dei momenti sul poliedro e assemblaggio dei pesi di cubatura
%**************************************************************************
% Coefficienti di normalizzazione della base di Chebyshev tensoriale
coeffs = tenscheb_norm2sq(chebyshev_indices);

% Momenti della base sul poliedro (teorema della divergenza sulle facce)
moments_ch = chebyshev_moments_polyhedron(vertices, facets, ade, ...
    chebyshev_indices, bbox);

% Riscalamento della griglia di riferimento sulla bounding box reale
XYZW_tens = scale_rule(XYZW_tens_ref, bbox);

% Pesi di cubatura ottenuti ricombinando la regola di riferimento
W = (XYZW_tens(:,4)) .* V_ref * (moments_ch ./ coeffs);
XYZW = [XYZW_tens(:,1:3) W];

%**************************************************************************
% Valutazione della cubatura
%**************************************************************************
fXYZW = feval(f, XYZW(:,1), XYZW(:,2), XYZW(:,3));
W = XYZW(:,4);
I = W' * fXYZW;

fprintf('Integrale: %1.15e\n', I);
fprintf('........................\n');
