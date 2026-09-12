clear;
clc;
close all;

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
addpath(fullfile(projectRoot, 'src'));

%**************************************************************************
%
%   Esempio: 
%       Cubatura su dominio concavo
%
%   Descrizione:
%       Questo esempio dimostra l'utilizzo del metodo OptimalPolyCuba3D per
%       l'integrazione di una funzione su un dominio poliedrale convessa,
%       rappresentato da una mesh superficiale triangolare. 
%
%**************************************************************************

%**************************************************************************
%  Parametri e caricamento mesh
%**************************************************************************

ade = 1;

fprintf('\n');
fprintf('**************************************************************\n');
fprintf('                 OPTIMALPOLYCUBA3D\n');
fprintf('        Cubatura su Dominio Poliedrale Convesso\n');
fprintf('**************************************************************\n');
fprintf('\n');

fprintf('Ade:                   %d\n', ade);

vertices = load('convex_vertex.dat');

facets = load('convex_tri.dat');

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

figure('Color','w');

patch('Vertices', vertices, ...
      'Faces', facets, ...
      'FaceColor', [0.7 0.7 0.7], ...
      'EdgeColor', [0.2 0.2 0.2], ...
      'FaceAlpha', 0.65);

hold on;

axis equal;
grid on;
box on;

xlabel('$x$', 'Interpreter','latex');
ylabel('$y$', 'Interpreter','latex');
zlabel('$z$', 'Interpreter','latex');

title('Dominio Poliedrale Convesso','Interpreter','latex');

view(3);

hold off;

%**************************************************************************
% Inizio regola di cubatura
%**************************************************************************

fprintf('\n');
fprintf('Inizio Cubatura...\n');

tic;
[XYZ, W] = OptimalPolyCuba3D(ade, vertices, facets);
elapsedTime = toc;

I = W' * f(XYZ(:,1), XYZ(:,2), XYZ(:,3));

%**************************************************************************
% Visualizzazione risultati e punti di cubatura
%**************************************************************************

fprintf('Fine Cubatura...\n');

fprintf('\n');
fprintf('Integrale numerico:    %.15e\n', I);
fprintf('Tempo di calcolo:      %.6e s\n', elapsedTime);
