clear;
clc;
close all;

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
addpath(fullfile(projectRoot, 'src'));

%**************************************************************************
%
%   Esempio: 
%       Cubatura su dominio generico
%
%   Descrizione:
%       Questo esempio dimostra l'utilizzo del metodo OptimalPolyCuba3D per
%       l'integrazione di una funzione su un dominio poliedrale generico,
%       rappresentato da una mesh superficiale triangolare. 
%
%**************************************************************************

%**************************************************************************
%  Parametri e caricamento mesh
%**************************************************************************

ade = 133;

fprintf('\n');
fprintf('**************************************************************\n');
fprintf('                 OPTIMALPOLYCUBA3D\n');
fprintf('        Cubatura su Dominio Poliedrale "Stanford Bunny" \n');
fprintf('**************************************************************\n');
fprintf('\n');

fprintf('Ade:                   %d\n', ade);

mesh = readSurfaceMesh('bun_zipper.ply');

vertices = mesh.Vertices;
facets   = mesh.Faces;

fprintf('Numero di vertici:     %d\n', size(vertices,1));
fprintf('Numero di facce:       %d\n', size(facets,1));

%**************************************************************************
% Definizione funzione integranda
%**************************************************************************

f = @(x,y,z) ones(size(x));
f_string = 'f(x,y,z) = 1';
fprintf('Funzione integranda:   %s\n', f_string);

%**************************************************************************
% Plot dominio poliedrale
%**************************************************************************

% Riduzione della mesh esclusivamente per la visualizzazione
% 1% dei triangoli utilizzati nella cubatura
[facets_plot, vertices_plot] = reducepatch(facets, vertices, 0.01);

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
% Inizio regola di cubatura
%**************************************************************************

fprintf('\n');
fprintf('Inizio Cubatura...\n');

tic;
[XYZ, W] = OPC3D(ade, vertices, facets, "D");
%[XYZ, W] = OPC3D(ade, vertices, facets, "GJ");
elapsedTime = toc;

I = W' * f(XYZ(:,1), XYZ(:,2), XYZ(:,3));

%**************************************************************************
% Visualizzazione risultati e punti di cubatura
%**************************************************************************

fprintf('Fine Cubatura...\n');

fprintf('\n');
fprintf('Integrale numerico:    %.15e\n', I);
fprintf('Tempo di calcolo:      %.6e s\n', elapsedTime);