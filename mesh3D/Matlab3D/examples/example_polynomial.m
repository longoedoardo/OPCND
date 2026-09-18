clear;
clc;
close all;

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
addpath(fullfile(projectRoot, 'src'));

%**************************************************************************
%
%   Esempio:
%       Cubatura su dominio cubico
%
%   Descrizione:
%       Questo esempio dimostra l'utilizzo del metodo OptimalPolyCuba3D
%       per l'integrazione di funzioni polinomiali su un cubo, rappresentato
%       mediante una mesh superficiale triangolare chiusa e orientata
%       secondo le normali esterne.
%
%       Vengono considerate tre funzioni integrande di diverso grado
%       polinomiale e i risultati numerici vengono confrontati con
%       i corrispondenti valori analitici.
%
%**************************************************************************

%**************************************************************************
%  Parametri e caricamento mesh
%**************************************************************************

fprintf('\n');
fprintf('**************************************************************\n');
fprintf('                 OPTIMALPOLYCUBA3D\n');
fprintf('                 Cubatura sul Cubo\n');
fprintf('**************************************************************\n');
fprintf('\n');

% Vertici del cubo [-1,1]^3
vertices = [
    -1 -1 -1;
     1 -1 -1;
     1  1 -1;
    -1  1 -1;
    -1 -1  1;
     1 -1  1;
     1  1  1;
    -1  1  1
];

% Facce triangolari della mesh superficiale.
% I vertici di ogni triangolo sono ordinati in senso antiorario
% se osservati dall'esterno del dominio, garantendo normali esterne.

facets = [
    1 3 2;
    1 4 3;

    5 6 7;
    5 7 8;

    1 2 6;
    1 6 5;

    4 8 7;
    4 7 3;

    1 5 8;
    1 8 4;

    2 3 7;
    2 7 6
];

fprintf('Numero di vertici:     %d\n', size(vertices,1));
fprintf('Numero di facce:       %d\n', size(facets,1));

%**************************************************************************
% Plot dominio poliedrale
%**************************************************************************

figure('Color','w');

patch('Vertices', vertices, ...
      'Faces', facets, ...
      'FaceColor', [0.7 0.7 0.7], ...
      'EdgeColor', [0.2 0.2 0.2], ...
      'FaceAlpha', 0.65);

axis equal;
grid on;
box on;

xlabel('$x$', 'Interpreter','latex');
ylabel('$y$', 'Interpreter','latex');
zlabel('$z$', 'Interpreter','latex');

title('Dominio Poliedrale', 'Interpreter','latex');

view(3);

%**************************************************************************
% Definizione funzione integranda f1
%**************************************************************************

ade = 1;

f1 = @(x,y,z) ones(size(x));
f1_string = 'f_1(x,y,z) = 1';

I_exact = 8;

fprintf('\n');
fprintf('--------------------------------------------------------------\n');
fprintf('Funzione integranda 1\n');
fprintf('--------------------------------------------------------------\n');
fprintf('Ade:                   %d\n', ade);
fprintf('Funzione integranda:   %s\n', f1_string);
fprintf('Integrale esatto:      %.15e\n', I_exact);

%**************************************************************************
% Inizio regola di cubatura
%**************************************************************************

fprintf('\n');
fprintf('Inizio Cubatura...\n');

tic;
[XYZ, W] = OPC3D(ade, vertices, facets, "D");
%[XYZ, W] = OPC3D(ade, vertices, facets, "GJ");
elapsedTime = toc;

I = W' * f1(XYZ(:,1), XYZ(:,2), XYZ(:,3));

fprintf('Fine Cubatura...\n');

error_abs = abs(I - I_exact);

fprintf('\n');
fprintf('Integrale numerico:    %.15e\n', I);
fprintf('Errore assoluto:       %.6e\n', error_abs);
fprintf('Tempo di calcolo:      %.6e s\n', elapsedTime);

%**************************************************************************
% Definizione funzione integranda f2
%**************************************************************************

ade = 2;

f2 = @(x,y,z) x.^2 + y.^2 + z.^2;
f2_string = 'f_2(x,y,z) = x^2 + y^2 + z^2';

I_exact = 8;

fprintf('\n');
fprintf('--------------------------------------------------------------\n');
fprintf('Funzione integranda 2\n');
fprintf('--------------------------------------------------------------\n');
fprintf('Ade:                   %d\n', ade);
fprintf('Funzione integranda:   %s\n', f2_string);
fprintf('Integrale esatto:      %.15e\n', I_exact);

%**************************************************************************
% Inizio regola di cubatura
%**************************************************************************

fprintf('\n');
fprintf('Inizio Cubatura...\n');

tic;
[XYZ, W] = OPC3D(ade, vertices, facets, "D");
%[XYZ, W] = OPC3D(ade, vertices, facets, "GJ");
elapsedTime = toc;


I = W' * f2(XYZ(:,1), XYZ(:,2), XYZ(:,3));

fprintf('Fine Cubatura...\n');

error_abs = abs(I - I_exact);

fprintf('\n');
fprintf('Integrale numerico:    %.15e\n', I);
fprintf('Errore assoluto:       %.6e\n', error_abs);
fprintf('Tempo di calcolo:      %.6e s\n', elapsedTime);

%**************************************************************************
% Definizione funzione integranda f3
%**************************************************************************

ade = 6;

f3 = @(x,y,z) x.^2 .* y.^2 .* z.^2;
f3_string = 'f_3(x,y,z) = x^2*y^2*z^2';

I_exact = 8/27;

fprintf('\n');
fprintf('--------------------------------------------------------------\n');
fprintf('Funzione integranda 3\n');
fprintf('--------------------------------------------------------------\n');
fprintf('Ade:                   %d\n', ade);
fprintf('Funzione integranda:   %s\n', f3_string);
fprintf('Integrale esatto:      %.15e\n', I_exact);

%**************************************************************************
% Inizio regola di cubatura
%**************************************************************************

fprintf('\n');
fprintf('Inizio Cubatura...\n');

tic;
[XYZ, W] = OPC3D(ade, vertices, facets, "D");
%[XYZ, W] = OPC3D(ade, vertices, facets, "GJ");
elapsedTime = toc;

I = W' * f3(XYZ(:,1), XYZ(:,2), XYZ(:,3));

fprintf('Fine Cubatura...\n');

error_abs = abs(I - I_exact);

fprintf('\n');
fprintf('Integrale numerico:    %.15e\n', I);
fprintf('Errore assoluto:       %.6e\n', error_abs);
fprintf('Tempo di calcolo:      %.6e s\n', elapsedTime);