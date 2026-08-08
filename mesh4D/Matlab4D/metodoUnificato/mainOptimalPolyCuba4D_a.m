
%**************************************************************************
%
%                         MAIN OptimalPolyCuba4D
%
%**************************************************************************
%
% Metodo che calcola l'integrale di una funzione f(x,y,z) su un poliedro 
% nel tempo [0,1] con approccio unificato, usando grado algebrico di esattezza "ade".
% 
% INPUT:
%   vertices : Matrice N x 3 dei vertici [x, y, z]
%   facets   : Struttura o cell array delle facce
%   ade      : Grado di esattezza algebrica
%   f        : Function handle della funzione integranda
%
%**************************************************************************

addpath('ShapeIndependent/');
addpath('ShapeDependent/');

fprintf('........................\n');
fprintf('Cubatura con OptimalPolyCuba4D (Metodo unificato) \n');
fprintf('........................\n');

ade = 4; % Grado di esattezza algebrica (Algebraic Degree of Exactness)
f = @(x,y,z,tau) x.*tau.*y.^3; % Funzione integranda f(x, y, z, tau)

fprintf('ade: %-3.0f\n', ade);

V0 = load('vertex.dat'); % Vertici 3D a istante iniziale tau = 0
V1 = load('vertex_new.dat'); % Vertici 3D a istante finale tau = 1
tri = load('tri.dat'); % Connettività delle facce triangolari 3D

num_vertici = size(V0, 1);
num_tri = size(tri, 1);

%**************************************************************************
%
% Costruzione dominio 4D
%
%**************************************************************************

vertici_4D = zeros(2 * num_vertici, 4);

% Configurazione iniziale (tau = 0)
vertici_4D(1:num_vertici, 1:3) = V0;
vertici_4D(1:num_vertici, 4) = 0.0;

% Configurazione finale (tau = 1)
vertici_4D(num_vertici + 1:2 * num_vertici, 1:3) = V1;
vertici_4D(num_vertici + 1:2 * num_vertici, 4) = 1.0;

% Costruzione delle iperfacce 3D del poliedro spaziotemporale
% Il poliedro 4D è racchiuso da: 1 Base (tau=0), 1 Tetto (tau=1) e 
% N_tri facce laterali (prismi che descrivono il movimento dei triangoli).
Hyperfacets = cell(2 + num_tri, 1);

% Iperfaccia 1: Base inferiore (tau = 0)
Hyperfacets{1}.Vertices_ID = 1:num_vertici;
Hyperfacets{1}.Sub_Topology = tri;
Hyperfacets{1}.Normal_4D = [0, 0, 0, -1];  % Normale uscente orientata nel passato

% Iperfaccia 2: Tetto superiore (tau = 1)
Hyperfacets{2}.Vertices_ID = (num_vertici + 1):(2 * num_vertici);
Hyperfacets{2}.Sub_Topology = tri + num_vertici; % Shift degli indici
Hyperfacets{2}.Normal_4D = [0, 0, 0, 1];   % Normale uscente orientata nel futuro

% Baricentro globale 4D per determinare l'orientamento corretto delle normali
baricentro = mean(vertici_4D, 1);

% Iperfacce laterali (Prismi spazio-temporali per ciascun triangolo)
for k = 1:num_tri
    n1 = tri(k, 1);
    n2 = tri(k, 2);
    n3 = tri(k, 3);

    n1_new = n1 + num_vertici;
    n2_new = n2 + num_vertici;
    n3_new = n3 + num_vertici;

    idx = 2 + k;
    Hyperfacets{idx}.Vertices_ID = [n1, n2, n3, n1_new, n2_new, n3_new];

    % Sotto-connettività 2D delle facce del prisma laterale
    Hyperfacets{idx}.Sub_Topology = [
        n1,     n3,     n2;     
        n1_new, n2_new, n3_new; 
        n1,     n2,     n2_new; 
        n1,     n2_new, n1_new;
        n2,     n3,     n3_new;
        n2,     n3_new, n2_new;
        n3,     n1,     n1_new;
        n3,     n1_new, n3_new
    ];

    % Calcolo analitico del vettore normale ortogonale in 4D
    v1 = vertici_4D(n2, :) - vertici_4D(n1, :);
    v2 = vertici_4D(n3, :) - vertici_4D(n1, :);
    v3 = vertici_4D(n1_new, :) - vertici_4D(n1, :);

    null_vector = null([v1; v2; v3]); % Spazio nullo ortogonale ai 3 vettori guida

    % Controllo orientamento verso l'esterno rispetto al baricentro
    baricentro_hyperfacet = mean(vertici_4D(Hyperfacets{idx}.Vertices_ID, :), 1);
    v_out = baricentro_hyperfacet - baricentro;

    if dot(null_vector, v_out) < 0
        null_vector = -null_vector;
    end

    Hyperfacets{idx}.Normal_4D = null_vector' / norm(null_vector);
end

%**************************************************************************
%
% INIZIO PARTE SHAPE-INDEPENDENT
%
%**************************************************************************

% Generazione griglia tensoriale di Gauss-Chebyshev sul dominio di riferimento [-1, 1]^4
XYZTW_tens_ref = cub_gausscheb_tens4D(2 * ade);

% Calcolo della dimensione dello spazio polinomiale 4D per il grado 'ade'
N = ((ade + 4) * (ade + 3) * (ade + 2) * (ade + 1)) / 24;
fprintf('Dimensione spazio polinomiale: %d\n', N);

% Inizializzazione matrice degli esponenti in ordine Grlex
chebyshev_indices = zeros(N, 4);
for i = 2:N
    chebyshev_indices(i, :) = mono_next_grlex(4, chebyshev_indices(i - 1, :));
end

% Estrazione coordinate e calcolo della matrice di Vandermonde di riferimento
X = XYZTW_tens_ref(:, 1:4);
V_ref = dCHEBVAND(ade, X, chebyshev_indices); 

% Calcolo dei coefficienti di normalizzazione al quadrato dei polinomi di Chebyshev
[coeffs] = tenscheb_norm2sq(chebyshev_indices);

%**************************************************************************
%
% INIZIO PARTE SHAPE-DEPENDENT
%
%**************************************************************************

% Calcolo dell'iper-rettangolo che racchiude l'intero dominio spazio-temporale
limiti_min = min(vertici_4D, [], 1);
limiti_max = max(vertici_4D, [], 1);
bbox = [limiti_min; limiti_max];

% Calcolo dei momenti di Chebyshev sul poliedro 4D tramite il teorema della divergenza
moments_ch = chebyshev_moments_polyhedron_4D(vertici_4D, Hyperfacets, ade, chebyshev_indices, bbox);

% Riscalamento della regola di cubatura sulla bounding box reale 4D
XYZTW_tens = scale_rule(XYZTW_tens_ref, bbox);

% Ricombinazione dei pesi di cubatura basata sui momenti del poliedro
W = (XYZTW_tens(:, 5)) .* V_ref * (moments_ch ./ coeffs);
XYZTW = [XYZTW_tens(:, 1:4), W];

% Valutazione della funzione e integrazione numerica finale
fXYZTW = f(XYZTW(:, 1), XYZTW(:, 2), XYZTW(:, 3), XYZTW(:, 4));
I = XYZTW(:, 5)' * fXYZTW;

fprintf('Integrale: %1.15e\n', I);
fprintf('........................\n');

%**************************************************************************
%
% Visualizzazione grafica
%
%**************************************************************************

figure('Name', 'Evoluzione Geometrica 3D', 'Color', 'w');
hold on; grid on; box on;

col_t0 = [0.20 0.45 0.95]; % Blu per t = 0
col_t1 = [0.95 0.30 0.35]; % Rosso per t = 1

patch('Vertices', V0, 'Faces', tri, 'FaceColor', col_t0, 'FaceAlpha', 0.05, ...
      'EdgeColor', col_t0, 'LineWidth', 1.0, 'DisplayName', 'Poliedro T=0');
patch('Vertices', V1, 'Faces', tri, 'FaceColor', col_t1, 'FaceAlpha', 0.05, ...
      'EdgeColor', col_t1, 'LineWidth', 1.0, 'DisplayName', 'Poliedro T=1');

for i = 1:num_vertici
    plot3([V0(i, 1) V1(i, 1)], [V0(i, 2) V1(i, 2)], [V0(i, 3) V1(i, 3)], ...
          'k--', 'LineWidth', 0.6, 'HandleVisibility', 'off'); 
end

plot3(V0(:, 1), V0(:, 2), V0(:, 3), 'o', 'MarkerSize', 4, 'MarkerFaceColor', col_t0, 'MarkerEdgeColor', 'w', 'HandleVisibility', 'off');
plot3(V1(:, 1), V1(:, 2), V1(:, 3), 'o', 'MarkerSize', 4, 'MarkerFaceColor', col_t1, 'MarkerEdgeColor', 'w', 'HandleVisibility', 'off');

xlabel('X'); ylabel('Y'); zlabel('Z');
title('Evoluzione geometrica: Poliedro tra T=0 e T=1');
view(3); camlight; lighting gouraud; axis equal;
legend('Location', 'northeast', 'Box', 'off');

figure('Name', 'Prisma Spaziotemporale 4D', 'Color', 'w');
hold on; grid on; box on;

col_prisma_0 = [0.6, 0.85, 0.90];
col_prisma_1 = [0.90, 0.60, 0.75];

offset = [1.8 0 0];
V1_vis = V1 + offset;

patch('Vertices', V0, 'Faces', tri, 'FaceColor', col_t0, 'FaceAlpha', 0.05, ...
      'EdgeColor', col_t0, 'LineStyle', ':', 'LineWidth', 1.1, 'DisplayName', 'Poliedro T=0');
patch('Vertices', V1_vis, 'Faces', tri, 'FaceColor', col_t1, 'FaceAlpha', 0.05, ...
      'EdgeColor', col_t1, 'LineStyle', ':', 'LineWidth', 1.1, 'DisplayName', 'Poliedro T=1');

% Estrazione di un singolo prisma rappresentativo dalla mesh
[~, k_prisma] = min(vecnorm(mean(V0(tri, :), 2) - mean(V0, 1), 2, 2));
nodes = Hyperfacets{2 + k_prisma}.Vertices_ID;
V_p = vertici_4D(nodes, :);
V_prisma = [V_p(1:3, 1:3); V_p(4:6, 1:3) + offset];

F_quad = [1 2 5 4; 2 3 6 5; 3 1 4 6];
patch('Vertices', V_prisma, 'Faces', F_quad, 'FaceColor', col_prisma_0, 'FaceAlpha', 0.82, ...
      'EdgeColor', [0.1 0.1 0.1], 'LineWidth', 1.7, 'DisplayName', 'Prisma 4D');
patch('Vertices', V_prisma, 'Faces', [1 2 3], 'FaceColor', col_prisma_0, 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'HandleVisibility', 'off');
patch('Vertices', V_prisma, 'Faces', [4 5 6], 'FaceColor', col_prisma_1, 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'HandleVisibility', 'off');

for i = 1:3
    plot3([V_prisma(i, 1) V_prisma(i+3, 1)], [V_prisma(i, 2) V_prisma(i+3, 2)], [V_prisma(i, 3) V_prisma(i+3, 3)], ...
          'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
end
scatter3(V_prisma(1:3, 1), V_prisma(1:3, 2), V_prisma(1:3, 3), 4, col_t0, 'filled', 'MarkerEdgeColor', 'w');
scatter3(V_prisma(4:6, 1), V_prisma(4:6, 2), V_prisma(4:6, 3), 4, col_t1, 'filled', 'MarkerEdgeColor', 'w');

title('Sezione di Prisma Spaziotemporale 4D');
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal; view(38, 24); camlight; lighting gouraud; material shiny;
legend('Location', 'bestoutside', 'Box', 'off');