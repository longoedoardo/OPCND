projectRoot = fileparts(mfilename('fullpath'));

addpath(fullfile(projectRoot, 'Matlab3D/src'));
addpath(fullfile(projectRoot, 'TetraMethod'));
addpath(fullfile(projectRoot, 'examples'));

%**************************************************************************
%
%   Esempio: 
%       Cubatura su dominio generico
%
%   Descrizione:
%       Questo esempio confronta le prestazioni dei metodi analizzati con e
%       senza tetraedralizzazione del dominio
%
%**************************************************************************

%**************************************************************************
%  Parametri e caricamento mesh
%**************************************************************************

ade = 1;

mesh = readSurfaceMesh('bun_zipper.ply');

vertices = mesh.Vertices;
facets   = mesh.Faces;

fprintf('Numero di vertici:     %d\n', size(vertices,1));
fprintf('Numero di facce:       %d\n', size(facets,1));

%**************************************************************************
% Definizione funzione integranda
%**************************************************************************

f = @(x,y,z) ones(size(x));

%**************************************************************************
% Plot dominio poliedrale
%**************************************************************************

% Riduzione della mesh esclusivamente per la visualizzazione
% 1% dei triangoli utilizzati nella cubatura
[facets_plot, vertices_plot] = reducepatch(facets, vertices, 0.03);

% Rotazione della mesh solo per la visualizzazione
vertices_plot = vertices_plot * [1  0  0;
                                  0  0  1;
                                  0 -1  0];

figure('Color','w');

patch('Vertices', vertices_plot, ...
      'Faces', facets_plot, ...
      'FaceColor', [0.75 0.75 0.75], ...
      'EdgeColor', [0.15 0.15 0.15], ...
      'LineWidth', 0.5);

axis equal;
grid on;
box on;
view(3);

xlabel('$x$', 'Interpreter','latex');
ylabel('$y$', 'Interpreter','latex');
zlabel('$z$', 'Interpreter','latex');

title('Stanford Bunny', 'Interpreter','latex');

%**************************************************************************
% Inizio regola di cubatura con Gauss-Jacobi
%**************************************************************************

[XYZ_GJ, W_GJ] = OPC3D(ade, vertices, facets, "GJ");

elapsedTime_GJ = timeit(@() OPC3D(ade, vertices, facets, "GJ"));

I_GJ = W_GJ' * f(XYZ_GJ(:,1), XYZ_GJ(:,2), XYZ_GJ(:,3));

%**************************************************************************
% Inizio regola di cubatura Thetra
%**************************************************************************

nGP = ceil((ade+1)/2);

[P_T, w_T] = TetraMethod(vertices, facets, nGP);

elapsedTime_T = timeit(@() TetraMethod(vertices, facets, nGP));

I_T = w_T' * f(P_T(:,1), P_T(:,2), P_T(:,3));

%**************************************************************************
% Inizio regola di cubatura con Dunavant
%**************************************************************************

[XYZ_D, W_D] = OPC3D(ade, vertices, facets, "D");

elapsedTime_D = timeit(@() OPC3D(ade, vertices, facets, "D"));

I_D = W_D' * f(XYZ_D(:,1), XYZ_D(:,2), XYZ_D(:,3));

%**************************************************************************
% Visualizzazione risultati e punti di cubatura
%**************************************************************************

% fprintf('Fine Cubatura...\n');
% 
% fprintf('\n');
% fprintf('Integrale numerico:    %.15e\n', I);
% fprintf('Tempo di calcolo:      %.6e s\n', elapsedTime);


% =========================================================================
% RISULTATI
% =========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('                    RISULTATI BUNNY\n');
fprintf('============================================================\n');

fprintf('\n');
fprintf('OptimalPolyCuba3D - Gauss-Jacobi\n');
fprintf('  Nodi     : %d\n', length(W_GJ));
fprintf('  Integrale: %.15e\n', I_GJ);
fprintf('  Tempo    : %.6e s\n', elapsedTime_GJ);

fprintf('\n');
fprintf('OptimalPolyCuba3D - Dunavant\n');
fprintf('  Nodi     : %d\n', length(W_D));
fprintf('  Integrale: %.15e\n', I_D);
fprintf('  Tempo    : %.6e s\n', elapsedTime_D);

fprintf('\n');
fprintf('TetraMethod\n');
fprintf('  Nodi     : %d\n', length(w_T));
fprintf('  Integrale: %.15e\n', I_T);
fprintf('  Tempo    : %.6e s\n', elapsedTime_T);

fprintf('\n');
fprintf('============================================================\n');