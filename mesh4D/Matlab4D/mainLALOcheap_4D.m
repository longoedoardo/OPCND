
%**************************************************************************
%
%                               MAIN LALO
%
%**************************************************************************
% Script principale per il calcolo cubatura cheap
% Data: 29/05/2026
%**************************************************************************

fprintf('........................\n ');

% Parametri modificabili di cubatura
ade = 5; % Grado di esattezza algebrica

fprintf('Cubatura con LALOcheap \n');
fprintf('........................\n ');
fprintf('ade: %-3.0f\n ', ade);

f = @(x,y,z,tau) ones(size(x)); % Funzione da inum_triegrare

%**************************************************************************

V0 = load('vertex.dat');
V1 = load('vertex_new.dat');
tri = load('tri.dat');

num_vertici = size(V0, 1);
num_tri = size(tri, 1);

% Generazione Matrice Vertici 4D [x, y, z, tau]
vertici_4D = zeros(2*num_vertici, 4);

% Blocco tau = 0
vertici_4D(1:num_vertici, 1:3) = V0;
vertici_4D(1:num_vertici, 4) = 0.0;

% Blocco tau = 1
vertici_4D(num_vertici+1:2*num_vertici, 1:3) = V1;
vertici_4D(num_vertici+1:2*num_vertici, 4) = 1.0;

% Allocazione Array delle Iperfacce 3D
Hyperfacets = cell(2 + num_tri, 1);

% IPERFACCIA 1: Base (tau = 0)
Hyperfacets{1}.Vertices_ID = 1:num_vertici;
Hyperfacets{1}.Sub_Topology = tri;
Hyperfacets{1}.Normal_4D = [0, 0, 0, -1];

% IPERFACCIA 2: Tetto (tau = 1)
Hyperfacets{2}.Vertices_ID = (num_vertici+1):(2*num_vertici);
Hyperfacets{2}.Sub_Topology = tri + num_vertici; % Shift degli indici per mappare il tetto
Hyperfacets{2}.Normal_4D = [0, 0, 0, 1];

% Calcolo baricentro
baricentro = mean(vertici_4D, 1);

% IPERFACCE LATERALI: Prismi spaziotemporali
for k = 1:num_tri
    n1 = tri(k, 1);
    n2 = tri(k, 2);
    n3 = tri(k, 3);

    n1_new = n1 + num_vertici;
    n2_new = n2 + num_vertici;
    n3_new = n3 + num_vertici;

    idx = 2 + k;
    Hyperfacets{idx}.Vertices_ID = [n1, n2, n3, n1_new, n2_new, n3_new];

    % Sotto-connettività 2D del prisma
    Hyperfacets{idx}.Sub_Topology = [
        n1, n3, n2;
        n1_new, n2_new, n3_new;
        n1, n2, n2_new;
        n1, n2_new, n1_new;
        n2, n3, n3_new;
        n2, n3_new, n2_new;
        n3, n1, n1_new;
        n3, n1_new, n3_new
        ];

    % Calcolo analitico della normale 4D ortogonale alle iperfacce laterali
    v1 = vertici_4D(n2, :) - vertici_4D(n1, :);
    v2 = vertici_4D(n3, :) - vertici_4D(n1, :);
    v3 = vertici_4D(n1_new, :) - vertici_4D(n1, :);

    % Calcolo dello spazio nullo (vettore ortogonale a v1, v2, v3 in 4D)
    null_vector = null([v1; v2; v3]);

    baricentro_hyperfacet = mean(vertici_4D(Hyperfacets{idx}.Vertices_ID, :), 1);
    v_out = baricentro_hyperfacet - baricentro;

    if dot(null_vector, v_out) < 0
        null_vector = -null_vector;
    end

    Hyperfacets{idx}.Normal_4D = null_vector' / norm(null_vector);
end

% Estrae tutte le normali 4D in una matrice (N_iperfacce x 4)
N_iperfacce = length(Hyperfacets);
matrice_normali = zeros(N_iperfacce, 4);

for j = 1:N_iperfacce
    matrice_normali(j, :) = Hyperfacets{j}.Normal_4D;
end

%**************************************************************************

% Preparazione iniziale, questa parte è indipendente dal poliedro specifico
% e dipende solo dal grado scelto

addpath('PrepCheap/')

% Calcolo punum_trii e pesi Gauss-Chebyshev, generazione griglia tensoriale di
% punum_trii di Gauss-Chebyshev

XYZTW_tens_ref=cub_gausscheb_tens4D(2*ade);

N = ((ade+4)*(ade+3)*(ade+2)*(ade+1))/24; % calcolo quanum_trii termini
% ci sono in un polinomio 4d di grado ade (dimensione spazio polinomiale)

fprintf('Dimensione spazio polinomiale: %-3.0f\n ', N);

chebyshev_indices = zeros(N,4); % terna degli esponenum_trii (i,j,k,l) per (x,y,z,tau)

% Riempio la tabella seguendo l'ordine GRLEX, inum_verticiece di generare
% combinazioni a caso prima li metto in ordine per somma degli esponenum_trii e
% poi in ordine alfabetico/numerico

for i=2:N
    chebyshev_indices(i,:) = mono_next_grlex(4,chebyshev_indices(i-1,:));
end

% Prendo solo le coordinate e scarto i pesi per il momenum_trio
X=XYZTW_tens_ref(:,1:4);
% Costruzione della matrice V, dove ogni colonna rappresenum_tria un polinomio
% di Chebyshev 3d valutato nei punum_trii della griglia, tutto su cubo [-1,1]^4
V_ref = dCHEBVAND(ade,X,chebyshev_indices);
W = XYZTW_tens_ref(:,5); % Estrae i pesi della griglia 4D

%**************************************************************************
%
% Calcolo i limiti minimi e massimi dell'iperpoliedro 4D
% Calcolo una sola volta la Bounding Box 4D (Iper-rettangolo spaziotemporale).
% Serve per scalare la griglia di punum_trii dal cubo di riferimenum_trio [-1,1]^4
% al volume/tempo reale del dominio.
%
%**************************************************************************


limiti_min = min(vertici_4D, [], 1); % Restituisce: [min_x, min_y, min_z, min_tau]
limiti_max = max(vertici_4D, [], 1); % Restituisce: [max_x, max_y, max_z, max_tau]

% Salvataggio nella matrice bbox (Formato 2 x 4)
% Riga 1: Valori Minimi [x, y, z, tau]
% Riga 2: Valori Massimi [x, y, z, tau]
bbox = [limiti_min; limiti_max];

%**************************************************************************

addpath('CubaturaFunzioniCheap/')

% Calcolo delle funzioni phi_i, phi_k, per la normalizzzione
[coeffs]=tenscheb_norm2sq(chebyshev_indices);

moments_ch = chebyshev_moments_polyhedron_4D(vertici_4D, Hyperfacets, ade, chebyshev_indices, bbox);

XYZTW_tens=scale_rule(XYZTW_tens_ref,bbox);

W=(XYZTW_tens(:,5)).*V_ref*(moments_ch./(coeffs)); 

XYZTW=[XYZTW_tens(:,1:4) W];

fXYZTW=feval(f,XYZTW(:,1),XYZTW(:,2),XYZTW(:,3), XYZTW(:,4));
W=XYZTW(:,5); % pesi
I= W'*fXYZTW; % Risultato finale cubatura

fprintf('Integrale: %1.15e\n',I)

fprintf('........................\n ');

%**************************************************************************
%
% CONFRONTO POLIEDRO INIZIALE VS FINALE
% 
%**************************************************************************

figure;
hold on; grid on; box on;

% Palette coordinata
col_t0 = [0.20 0.45 0.95];
col_t1 = [0.95 0.30 0.35];

% 1. Poliedri Ghost
patch('Vertices', V0, 'Faces', tri, 'FaceColor', col_t0, 'FaceAlpha', 0.05, ...
      'EdgeColor', col_t0, 'LineWidth', 1.0, 'DisplayName', 'Poliedro T=0');
patch('Vertices', V1, 'Faces', tri, 'FaceColor', col_t1, 'FaceAlpha', 0.05, ...
      'EdgeColor', col_t1, 'LineWidth', 1.0, 'DisplayName', 'Poliedro T=1');

% 2. Linee di spostamenum_trio (HandleVisibility off per leggenda pulita)
for i = 1:num_vertici
    plot3([V0(i,1) V1(i,1)], [V0(i,2) V1(i,2)], [V0(i,3) V1(i,3)], ...
          'k--', 'LineWidth', 0.6, 'HandleVisibility', 'off'); 
end

% Vertici evidenziati
plot3(V0(:,1), V0(:,2), V0(:,3), 'o', 'MarkerSize', 4, 'MarkerFaceColor', col_t0, 'MarkerEdgeColor', 'w', 'HandleVisibility', 'off');
plot3(V1(:,1), V1(:,2), V1(:,3), 'o', 'MarkerSize', 4, 'MarkerFaceColor', col_t1, 'MarkerEdgeColor', 'w', 'HandleVisibility', 'off');

% Estetica
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Evoluzione geometrica: Poliedri 4D');
view(3); camlight; lighting gouraud; axis equal;
legend('Location', 'northeast', 'Box', 'off');

%**************************************************************************
%
% PLOT PRISMA 4D
%
%**************************************************************************
% 
% figure
% hold on; grid on; box on;
% 
% % Palette coordinata
% col_t0 = [0.20 0.45 0.95];
% col_t1 = [0.95 0.30 0.35];
% col_prisma_0 = [0.6, 0.85, 0.90];
% col_prisma_1 = [0.90, 0.60, 0.75];
% 
% % Offset visuale
% offset = [1.8 0 0];
% V1_vis = V1 + offset;
% 
% % Poliedri
% patch('Vertices', V0, 'Faces', tri, 'FaceColor', col_t0, 'FaceAlpha', 0.05, ...
%       'EdgeColor', col_t0, 'LineStyle', ':', 'LineWidth', 1.1, 'DisplayName', 'Poliedro T=0');
% patch('Vertices', V1_vis, 'Faces', tri, 'FaceColor', col_t1, 'FaceAlpha', 0.05, ...
%       'EdgeColor', col_t1, 'LineStyle', ':', 'LineWidth', 1.1, 'DisplayName', 'Poliedro T=1');
% 
% % Estrazione Prisma Automatico
% [~, k_prisma] = min(vecnorm(mean(V0(tri,:), 2) - mean(V0, 1), 2, 2));
% nodes = Hyperfacets{2+k_prisma}.Vertices_ID;
% V_p = vertici_4D(nodes, :);
% V_prisma = [V_p(1:3, 1:3); V_p(4:6, 1:3) + offset];
% 
% % Disegno Prisma
% F_quad = [1 2 5 4; 2 3 6 5; 3 1 4 6];
% patch('Vertices', V_prisma, 'Faces', F_quad, 'FaceColor', col_prisma_0, 'FaceAlpha', 0.82, ...
%       'EdgeColor', [0.1 0.1 0.1], 'LineWidth', 1.7, 'DisplayName', 'Prisma 4D');
% patch('Vertices', V_prisma, 'Faces', [1 2 3], 'FaceColor', col_prisma_0, 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'HandleVisibility', 'off');
% patch('Vertices', V_prisma, 'Faces', [4 5 6], 'FaceColor', col_prisma_1, 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'HandleVisibility', 'off');
% 
% % Connessioni e Vertici
% for i = 1:3
%     plot3([V_prisma(i,1) V_prisma(i+3,1)], [V_prisma(i,2) V_prisma(i+3,2)], [V_prisma(i,3) V_prisma(i+3,3)], ...
%           'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
% end
% scatter3(V_prisma(1:3,1), V_prisma(1:3,2), V_prisma(1:3,3), 4, col_t0, 'filled', 'MarkerEdgeColor', 'w');
% scatter3(V_prisma(4:6,1), V_prisma(4:6,2), V_prisma(4:6,3), 4, col_t1, 'filled', 'MarkerEdgeColor', 'w');
% 
% title('Prisma Spaziotemporale 4D');
% xlabel('X'); ylabel('Y'); zlabel('Z');
% axis equal; view(38,24); camlight; lighting gouraud; material shiny;
% legend('Location', 'bestoutside', 'Box', 'off');