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
%       Questo esempio dimostra l'utilizzo del metodo OPC3D per
%       l'integrazione di una funzione su un dominio poliedrale concavo,
%       rappresentato da una mesh superficiale triangolare. 
%
%**************************************************************************

%**************************************************************************
%  Parametri e caricamento mesh
%**************************************************************************

ade = 1;

fprintf('\n');
fprintf('**************************************************************\n');
fprintf('                          OPC3D\n');
fprintf('        Cubatura su Dominio Poliedrale Concavo\n');
fprintf('**************************************************************\n');
fprintf('\n');

fprintf('Ade:                   %d\n', ade);

vertices = load('concave_vertex.dat');

facets = load('concave_tri.dat');

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

title('Dominio Poliedrale Concavo','Interpreter','latex');

view(3);

hold off;

%**************************************************************************
% Inizio regola di cubatura
%**************************************************************************

fprintf('\n');
fprintf('Inizio Cubatura...\n');

[XYZ, W] = OPC3D(ade, vertices, facets, "D");
%[XYZ, W] = OPC3D(ade, vertices, facets, "GJ");
elapsedTime = timeit(@() OPC3D(ade, vertices, facets, "D"));

I = W' * f(XYZ(:,1), XYZ(:,2), XYZ(:,3));

%**************************************************************************
% Visualizzazione risultati e punti di cubatura
%**************************************************************************

fprintf('Fine Cubatura...\n');

fprintf('\n');
fprintf('Integrale numerico:    %.15e\n', I);
fprintf('Tempo di calcolo:      %.6e s\n', elapsedTime);
