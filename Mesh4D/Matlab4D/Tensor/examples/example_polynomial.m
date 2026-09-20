clear;
clc;
close all;

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
addpath(fullfile(projectRoot, 'src'));

%**************************************************************************
%
%   Esempio:
%       Cubatura 4D su dominio poliedrale (Cubo) in movimento (traslazione
%       e deformazione)
%
%   Descrizione:
%       Questo esempio dimostra l'utilizzo del metodo
%       OPC4D_Tensor per l'integrazione di funzioni su un 
%       dominio poliedrale in movimento per tau in [0,1], rappresentato 
%       mediante una mesh superficiale triangolare chiusa e orientata 
%       secondo le normali esterne.
%
%       Vengono considerate tre funzioni integrande di diverso grado
%       polinomiale e i risultati numerici vengono confrontati con
%       i corrispondenti valori analitici.
%
%       Il dominio viene deformato nel tempo mediante un'interpolazione
%       lineare tra una configurazione iniziale e una configurazione finale.
%
%**************************************************************************

%**************************************************************************
% Parametri e caricamento mesh
%**************************************************************************
fprintf('\n');
fprintf('**************************************************************\n');
fprintf('                   OPC4D_Tensor\n');
fprintf('                 Metodo Tensoriale\n');
fprintf('            Cubatura su Cubo in Movimento\n');
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

% Vertici del cubo nella configurazione finale.
% Il dominio viene traslato e deformato lungo le tre direzioni.
vertici_finali = [
-0.90 -0.95 -0.90;
 1.10 -0.95 -0.90;
 1.10  1.05 -0.90;
-0.90  1.05 -0.90;
-0.90 -0.95  1.10;
 1.10 -0.95  1.10;
 1.10  1.05  1.10;
-0.90  1.05  1.10
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

    % Interpolazione convessa tra configurazione iniziale e finale
    vertici_tau = (1-tau) * vertici_iniziali + ...
                  tau * vertici_finali;

    % Aggiornamento della geometria
    set(h, 'Vertices', vertici_tau);

    % Aggiornamento del titolo
    title(sprintf('Dominio Poliedrale in Movimento, $\\tau = %.2f$', tau), ...
          'Interpreter','latex');

    drawnow;

    % Controllo della velocita' dell'animazione
    pause(0.01);

end


%**************************************************************************
% Definizione funzione integranda f1
%**************************************************************************

ade = 1;

f1 = @(x,y,z,tau) ones(size(x));
f1_string = 'f_1(x,y,z,\tau) = 1';

% Volume iniziale = 8
% Volume finale   = 2*2*2 = 8
% Il volume rimane quindi costante durante il movimento.

I_exact = 8;

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
[XYZT, W] = OPC4D_Tensor(ade, vertici_iniziali,vertici_finali,facets, 'GJCC');
% [XYZT, W] = OPC4D_Tensor(ade, vertici_iniziali,vertici_finali,facets, 'GJL');
% [XYZT, W] = OPC4D_Tensor(ade, vertici_iniziali,vertici_finali,facets, 'DGL');
% [XYZT, W] = OPC4D_Tensor(ade, vertici_iniziali,vertici_finali,facets, 'DCC');
elapsedTime = toc;

I = W' * f1(XYZT(:,1), XYZT(:,2), XYZT(:,3), XYZT(:,4));

fprintf('Fine Cubatura...\n');

error_abs = abs(I - I_exact);

fprintf('\n');
fprintf('Numero di nodi 4D:    %d\n', size(XYZT,1));
fprintf('Integrale numerico:    %.15e\n', I);
fprintf('Errore assoluto:       %.6e\n', error_abs);
fprintf('Tempo di calcolo:      %.6e s\n', elapsedTime);

%**************************************************************************
% Definizione funzione integranda f2
%**************************************************************************

ade = 2;

f2 = @(x,y,z,tau) x.^2 + y.^2 + z.^2 + tau.^2;
f2_string = 'f_2(x,y,z,\tau) = x^2 + y^2 + z^2 + \tau^2';

I_exact = 8 + 8.18/3;

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
[XYZT, W] = OPC4D_Tensor(ade, vertici_iniziali,vertici_finali,facets, 'GJCC');
% [XYZT, W] = OPC4D_Tensor(ade, vertici_iniziali,vertici_finali,facets, 'GJL');
% [XYZT, W] = OPC4D_Tensor(ade, vertici_iniziali,vertici_finali,facets, 'DGL');
% [XYZT, W] = OPC4D_Tensor(ade, vertici_iniziali,vertici_finali,facets, 'DCC');
elapsedTime = toc;

I = W' * f2(XYZT(:,1), XYZT(:,2), XYZT(:,3), XYZT(:,4));

fprintf('Fine Cubatura...\n');

error_abs = abs(I - I_exact);

fprintf('\n');
fprintf('Numero di nodi 4D:    %d\n', size(XYZT,1));
fprintf('Integrale numerico:    %.15e\n', I);
fprintf('Errore assoluto:       %.6e\n', error_abs);
fprintf('Tempo di calcolo:      %.6e s\n', elapsedTime);

%**************************************************************************
% Definizione funzione integranda f3
%**************************************************************************

ade = 8;

f3 = @(x,y,z,tau) x.^2 .* y.^2 .* z.^2 .* tau.^2;
f3_string = 'f_3(x,y,z,\tau) = x^2*y^2*z^2*\tau^2';

I_exact = 29150263/283500000;

fprintf('\n');
fprintf('--------------------------------------------------------------\n');
fprintf('Funzione integranda 3\n');
fprintf('--------------------------------------------------------------\n');
fprintf('Ade:                   %d\n', ade);
fprintf('Numero di nodi tau:    %d\n', n_tau);
fprintf('Funzione integranda:   %s\n', f3_string);
fprintf('Integrale esatto:      %.15e\n', I_exact);

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

I = W' * f3(XYZT(:,1), XYZT(:,2), XYZT(:,3), XYZT(:,4));

fprintf('Fine Cubatura...\n');

error_abs = abs(I - I_exact);

fprintf('\n');
fprintf('Numero di nodi 4D:    %d\n', size(XYZT,1));
fprintf('Integrale numerico:    %.15e\n', I);
fprintf('Errore assoluto:       %.6e\n', error_abs);
fprintf('Tempo di calcolo:      %.6e s\n', elapsedTime);