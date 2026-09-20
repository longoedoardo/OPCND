clear;
clc;
close all;

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
addpath(fullfile(projectRoot, 'src'));

%**************************************************************************
%
%   Esempio:
%       Cubatura 4D su dominio poliedrale (Cubo) in rototraslazione rigida
%
%   Descrizione:
%       Questo esempio dimostra l'utilizzo del metodo OPC4D_MoL
%       per l'integrazione di funzioni su un dominio poliedrale in movimento
%       per tau in [0,1], rappresentato mediante una mesh superficiale
%       triangolare chiusa e orientata secondo le normali esterne.
%
%       Il dominio è sottoposto a una rototraslazione rigida nel tempo:
%       una rotazione attorno all'asse z combinata con una traslazione.
%
%       La rotazione e la traslazione non modificano la forma né il volume
%       del cubo durante il movimento.
%
%       Vengono considerate funzioni integrande di diverso grado
%       polinomiale e i risultati numerici vengono confrontati con
%       i corrispondenti valori analitici.
%
%**************************************************************************

%**************************************************************************
% Parametri e caricamento mesh
%**************************************************************************

fprintf('\n');
fprintf('**************************************************************\n');
fprintf('                     OPC4D_MoL\n');
fprintf('                Metodo delle Linee\n');
fprintf('       Cubatura su Cubo in Rototraslazione\n');
fprintf('**************************************************************\n');
fprintf('\n');

n_tau = 10;

% Vertici del cubo nella configurazione iniziale
vertici_iniziali = [
-1 -1 -1;
 1 -1 -1;
 1  1 -1;
-1  1 -1;
-1 -1  1;
 1 -1  1;
 1  1  1;
-1  1  1
];

% La configurazione finale è ottenuta applicando una rotazione di pi/4
% attorno all'asse z e una traslazione [0.20, 0.10, 0.15].

theta_finale = pi/4;

R_finale = [
cos(theta_finale), -sin(theta_finale), 0;
sin(theta_finale),  cos(theta_finale), 0;
0,                  0,                 1
];

traslazione_finale = [0.20, 0.10, 0.15];

vertici_finali = (R_finale * vertici_iniziali')' + traslazione_finale;

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

fprintf('Numero di vertici:     %d\n', size(vertici_iniziali,1));
fprintf('Numero di facce:       %d\n', size(facets,1));
fprintf('Numero di nodi tau:    %d\n', n_tau);

%**************************************************************************
% Plot dominio poliedrale iniziale e finale
%**************************************************************************

figure('Color','w');

h = patch('Vertices', vertici_iniziali, ...
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

view(3);

% Limiti degli assi fissati sull'intero movimento
vertici_complessivi = [vertici_iniziali; vertici_finali];

xlim([min(vertici_complessivi(:,1)), max(vertici_complessivi(:,1))]);
ylim([min(vertici_complessivi(:,2)), max(vertici_complessivi(:,2))]);
zlim([min(vertici_complessivi(:,3)), max(vertici_complessivi(:,3))]);

% Animazione del dominio
n_frame = 100;

for k = 1:n_frame

    tau = (k-1)/(n_frame-1);

    % Rotazione progressiva attorno all'asse z
    theta = (pi/4) * tau;

    R = [
        cos(theta), -sin(theta), 0;
        sin(theta),  cos(theta), 0;
        0,           0,          1
    ];

    % Traslazione progressiva
    translation = [
        0.20 * tau, ...
        0.10 * tau, ...
        0.15 * tau
    ];

    % Rototraslazione rigida
    vertici_tau = (R * vertici_iniziali')' + translation;

    % Aggiornamento della geometria
    set(h, 'Vertices', vertici_tau);

    % Aggiornamento del titolo
    title(sprintf( ...
        'Cubo in Rototraslazione Rigida, $\\tau = %.2f$', tau), ...
        'Interpreter','latex');

    drawnow;

    pause(0.01);

end

%**************************************************************************
% Definizione funzione integranda f1
%**************************************************************************

ade = 1;

f1 = @(x,y,z,tau) ones(size(x));
f1_string = 'f_1(x,y,z,\tau) = 1';

I_exact = 16/3 + 4*sqrt(2)/3;

fprintf('\n');
fprintf('--------------------------------------------------------------\n');
fprintf('Funzione integranda 1\n');
fprintf('--------------------------------------------------------------\n');
fprintf('Ade:                   %d\n', ade);
fprintf('Numero di nodi tau:    %d\n', n_tau);
fprintf('Funzione integranda:   %s\n', f1_string);
fprintf('Integrale esatto:      %.15e\n', I_exact);

%**************************************************************************
% Inizio regola di cubatura
%**************************************************************************

fprintf('\n');
fprintf('Inizio Cubatura...\n');

tic;

[XYZtau, W] = OPC4D_MoL(ade, n_tau, vertici_iniziali, vertici_finali, facets, 'GJ');
% [XYZtau, W] = OPC4D_MoL(ade, n_tau, vertici_iniziali, vertici_finali, facets, 'D');

elapsedTime = toc;

I = W' * f1(XYZtau(:,1), XYZtau(:,2), XYZtau(:,3), XYZtau(:,4));

fprintf('Fine Cubatura...\n');

error_abs = abs(I - I_exact);

fprintf('\n');
fprintf('Numero di nodi 4D:     %d\n', size(XYZtau,1));
fprintf('Integrale numerico:    %.15e\n', I);
fprintf('Errore assoluto:       %.6e\n', error_abs);
fprintf('Tempo di calcolo:      %.6e s\n', elapsedTime);

%**************************************************************************
% Definizione funzione integranda f2
%**************************************************************************

ade = 2;

f2 = @(x,y,z,tau) tau.^2;
f2_string = 'f_2(x,y,z,\tau) = \tau^2';

I_exact = 28/15 + 2*sqrt(2)/5;

fprintf('\n');
fprintf('--------------------------------------------------------------\n');
fprintf('Funzione integranda 2\n');
fprintf('--------------------------------------------------------------\n');
fprintf('Ade:                   %d\n', ade);
fprintf('Numero di nodi tau:    %d\n', n_tau);
fprintf('Funzione integranda:   %s\n', f2_string);
fprintf('Integrale esatto:      %.15e\n', I_exact);

%**************************************************************************
% Inizio regola di cubatura
%**************************************************************************

fprintf('\n');
fprintf('Inizio Cubatura...\n');

tic;

[XYZtau, W] = OPC4D_MoL(ade, n_tau, vertici_iniziali, vertici_finali, facets, 'GJ');
% [XYZtau, W] = OPC4D_MoL(ade, n_tau, vertici_iniziali, vertici_finali, facets, 'D');

elapsedTime = toc;

I = W' * f2(XYZtau(:,1), XYZtau(:,2), XYZtau(:,3), XYZtau(:,4));

fprintf('Fine Cubatura...\n');

error_abs = abs(I - I_exact);

fprintf('\n');
fprintf('Numero di nodi 4D:     %d\n', size(XYZtau,1));
fprintf('Integrale numerico:    %.15e\n', I);
fprintf('Errore assoluto:       %.6e\n', error_abs);
fprintf('Tempo di calcolo:      %.6e s\n', elapsedTime);