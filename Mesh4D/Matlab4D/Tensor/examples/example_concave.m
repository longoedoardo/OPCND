clear;
clc;
close all;

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
addpath(fullfile(projectRoot, 'src'));

%**************************************************************************
%
%   Esempio:
%       Cubatura 4D su dominio poliedrale concavo in movimento
%
%   Descrizione:
%       Questo esempio dimostra l'utilizzo del metodo
%       OPC4D_Tensor per l'integrazione di una funzione su un
%       dominio poliedrale concavo in movimento per tau in [0,1],
%       rappresentato mediante una mesh superficiale triangolare.
%
%       Il test considera la funzione costante f = 1. L'integrale 4D esatto
%       viene calcolato integrando analiticamente il volume istantaneo del
%       dominio lungo l'intervallo temporale [0,1].
%
%**************************************************************************

%**************************************************************************
% Parametri e caricamento mesh
%**************************************************************************

ade = 1;

fprintf('\n');
fprintf('**************************************************************\n');
fprintf('                   OPC4D_Tensor\n');
fprintf('              - Metodo Tensoriale -\n');
fprintf('    Cubatura su Dominio Poliedrale Concavo in Movimento\n');
fprintf('**************************************************************\n');
fprintf('\n');

fprintf('Ade:                   %d\n', ade);

vertici_iniziali = load('concave_vertex.dat');
facets           = load('concave_tri.dat');
vertici_finali   = load('concave_vertex_new.dat');

fprintf('Numero di vertici:     %d\n', size(vertici_iniziali,1));
fprintf('Numero di facce:       %d\n', size(facets,1));

%**************************************************************************
% Definizione funzione integranda
%**************************************************************************

f = @(x,y,z,tau) ones(size(x));
f_string = 'f(x,y,z,tau) = 1';

fprintf('Funzione integranda:   %s\n', f_string);

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

%**************************************************************************
% Animazione del dominio
%**************************************************************************

n_frame = 100;

for k = 1:n_frame

    tau = (k-1)/(n_frame-1);

    % Interpolazione convessa tra configurazione iniziale e finale
    vertici_tau = (1-tau) * vertici_iniziali + ...
                  tau * vertici_finali;

    % Aggiornamento della geometria
    set(h, 'Vertices', vertici_tau);

    % Aggiornamento del titolo
    title(sprintf( ...
        'Dominio Poliedrale Concavo in Movimento, $\\tau = %.2f$', tau), ...
        'Interpreter','latex');

    drawnow;

    % Controllo della velocita' dell'animazione
    pause(0.01);

end


%**************************************************************************
% Inizio regola di cubatura
%**************************************************************************

fprintf('\n');
fprintf('Inizio Cubatura...\n');

tic;
[XYZT, W] = OPC4D_Tensor(ade, vertici_iniziali,vertici_finali,facets, 'GJCC');
% [XYZT, W] = OPC4D_Tensor(ade, vertici_iniziali,vertici_finali,facets, 'GJL');
% [XYZT, W] = OPC4D_Tensor(ade, vertici_iniziali,vertici_finali,facets, 'DGL');
% [XYZT, W] = OPC4D_Tensor(ade, vertici_iniziali,vertici_finali,facets, 'DCC');
elapsedTime = toc;

I = W' * f(XYZT(:,1), XYZT(:,2), XYZT(:,3), XYZT(:,4));


%**************************************************************************
% Visualizzazione risultati e punti di cubatura
%**************************************************************************

fprintf('Fine Cubatura...\n');
fprintf('\n');

fprintf('Numero di nodi 4D:    %d\n', size(XYZT,1));
fprintf('Integrale numerico:   %.15e\n', I);
fprintf('Tempo di calcolo:     %.6e s\n', elapsedTime);
