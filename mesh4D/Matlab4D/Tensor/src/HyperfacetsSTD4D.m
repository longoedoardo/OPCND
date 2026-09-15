function [vertici_4D, Hyperfacets] = HyperfacetsSTD4D(vertici_iniziali, vertici_finali, facets)

%**************************************************************************
%
% function [vertici_4D, Hyperfacets] = HyperfacetsSTD4D(vertici_iniziali, ...
%                                               vertici_finali, facets)
%
% Costruisce la rappresentazione geometrica del dominio spazio-temporale
% 4D associato all'evoluzione di un dominio poliedrico tridimensionale
% dalla configurazione iniziale, corrispondente a tau = 0, alla
% configurazione finale, corrispondente a tau = 1.
%
% Il dominio spazio-temporale è rappresentato in R^4 mediante le
% coordinate (x,y,z,tau). Il suo bordo è costituito da due iperfacce
% tridimensionali, corrispondenti alle configurazioni iniziale e finale,
% e da una collezione di iperfacce laterali associate alle facce
% triangolari della mesh superficiale.
%
% Ogni iperfaccia laterale è ottenuta collegando una faccia triangolare
% della configurazione iniziale con la corrispondente faccia della
% configurazione finale, formando un prisma triangolare nello spazio-
% tempo.
%
% Per ciascuna iperfaccia laterale viene inoltre calcolata la normale
% esterna unitaria in R^4. Il verso della normale viene determinato
% confrontandola con il baricentro globale del dominio spazio-temporale.
%
% La funzione viene utilizzata nella costruzione della regola di
% cubatura quadridimensionale basata sui momenti di Chebyshev.
%
%**************************************************************************

%
% INPUT:
%
%   vertici_iniziali:
%       Matrice N_v x 3 contenente le coordinate cartesiane dei vertici
%       della mesh nella configurazione iniziale, corrispondente a
%       tau = 0.
%
%   vertici_finali:
%       Matrice N_v x 3 contenente le coordinate cartesiane dei vertici
%       della mesh nella configurazione finale, corrispondente a
%       tau = 1.
%
%       La connettività della mesh deve essere la stessa utilizzata per
%       vertici_iniziali.
%
%   facets:
%       Matrice N_f x 3 contenente la connettività della mesh superficiale
%       triangolare. Ogni riga contiene gli indici dei tre vertici che
%       definiscono una faccia triangolare.
%
%**************************************************************************

%
% OUTPUT:
%
%   vertici_4D:
%       Matrice (2*N_v) x 4 contenente i vertici del dominio
%       spazio-temporale.
%
%       I primi N_v vertici appartengono alla configurazione iniziale
%       e hanno coordinata tau = 0; i successivi N_v vertici appartengono
%       alla configurazione finale e hanno coordinata tau = 1.
%
%   Hyperfacets:
%       Cella contenente la rappresentazione delle iperfacce del bordo
%       del dominio spazio-temporale 4D.
%
%       Le prime due celle corrispondono rispettivamente alle sezioni
%       tau = 0 e tau = 1. Le celle successive rappresentano le
%       iperfacce laterali associate alle facce triangolari della mesh.
%
%       Per ciascuna iperfaccia vengono memorizzati:
%
%           Vertices_ID
%               Indici dei vertici appartenenti all'iperfaccia.
%
%           Sub_Topology
%               Connettività triangolare utilizzata per rappresentare
%               la superficie dell'iperfaccia.
%
%           Normal_4D
%               Normale esterna unitaria all'iperfaccia in R^4.
%
%**************************************************************************

% Inizio costruzione dominio spazio-temporale 4D

% Numero di vertici della mesh tridimensionale
num_vertici = size(vertici_iniziali, 1);

% Numero di facce triangolari
num_tri = size(facets, 1);

% Il dominio spazio-temporale e' delimitato dalle due configurazioni
% tridimensionali agli istanti tau = 0, tau = 1.
% I primi num_vertici corrispondono alla configurazione iniziale, i
% successivi corrispondono alla configurazione finale

vertici_4D = zeros(2 * num_vertici, 4);

% Configurazione iniziale (tau = 0)
vertici_4D(1:num_vertici, 1:3) = vertici_iniziali;
vertici_4D(1:num_vertici, 4) = 0.0;

% Configurazione finale (tau = 1)
vertici_4D(num_vertici + 1:2 * num_vertici, 1:3) = vertici_finali;
vertici_4D(num_vertici + 1:2 * num_vertici, 4) = 1.0;


% Inizio costruzione delle iperfacce del bordo 4D

% Il bordo del dominio spaziotemporale e' costituito da una iperfaccia
% inferiore in tau = 0, una iperfaccia superiore in tau = 1 e una iperfaccia
% laterale per ciascuna faccia triangolare della mesh. Le iperfacce
% laterali sono prismi triangolari nello spaziotempo e rappresentano
% l'evoluzione delle facce della mesh tra tau = 0 e tau = 1.

Hyperfacets = cell(2 + num_tri, 1);

% Iperfaccia 1: Base inferiore (tau = 0)
Hyperfacets{1}.Vertices_ID = 1:num_vertici;
Hyperfacets{1}.Sub_Topology = facets;
Hyperfacets{1}.Normal_4D = [0, 0, 0, -1];  % Normale uscente orientata nel passato

% Iperfaccia 2: Tetto superiore (tau = 1)
Hyperfacets{2}.Vertices_ID = (num_vertici + 1):(2 * num_vertici);

Hyperfacets{2}.Sub_Topology = facets + num_vertici; % Shift degli indici,
% gli indici devono essere traslati di num_vertici per riferirsi ai
% vertici della configurazione finale

Hyperfacets{2}.Normal_4D = [0, 0, 0, 1];   % Normale uscente orientata nel futuro

% Baricentro globale 4D per determinare l'orientamento corretto delle normali
baricentro = mean(vertici_4D, 1);

% Iperfacce laterali (Prismi spazio-temporali per ciascun triangolo)
for k = 1:num_tri

    % Vertici della faccia triangolare nella configurazione iniziale
    n1 = facets(k, 1);
    n2 = facets(k, 2);
    n3 = facets(k, 3);

    % Vertici della faccia triangolare nella configurazione finale
    n1_new = n1 + num_vertici;
    n2_new = n2 + num_vertici;
    n3_new = n3 + num_vertici;

    % Indice dell'iperfaccia laterale nella struttura Hyperfacets.
    idx = 2 + k;

    % Vertici che costituiscono il prisma spazio-temporale associato alla
    % faccia triangolare k-esima
    Hyperfacets{idx}.Vertices_ID = [n1, n2, n3, n1_new, n2_new, n3_new];

    % Sotto-connettività 2D delle facce del prisma laterale. Il prisma
    % triangolare e' costituito da una facia triangolare iniziale, una
    % faccia triangolare finale e tre facce laterali quadrilatere, ciascuna
    % suddivisa in due triangoli per ottenere una rappresentazione
    % triangolare della superficie

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
    % La normale viene calcolata determinando un vettore appartenente allo
    % spazio nullo dei tre vettori tangenti linearmente indipendenti

    v1 = vertici_4D(n2, :) - vertici_4D(n1, :);
    v2 = vertici_4D(n3, :) - vertici_4D(n1, :);
    v3 = vertici_4D(n1_new, :) - vertici_4D(n1, :);

    null_vector = null([v1; v2; v3]); % Spazio nullo ortogonale ai 3 vettori tangenti

    % Controllo orientamento verso l'esterno rispetto al baricentro
    baricentro_hyperfacet = mean(vertici_4D(Hyperfacets{idx}.Vertices_ID, :), 1);

    % Vettore che punta dal baricentro globale verso l'iperfaccia.
    % Il verso della normale viene scelto in modo che punti verso l'esterno
    v_out = baricentro_hyperfacet - baricentro;

    if dot(null_vector, v_out) < 0
        null_vector = -null_vector;
    end

    % Normalizzazione della normale 4D
    Hyperfacets{idx}.Normal_4D = null_vector' / norm(null_vector);
end


end