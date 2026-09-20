function [XYZT, W] = OPC4D_Tensor(ade, vertici_iniziali, vertici_finali, facets, method)

%**************************************************************************
%
% function [XYZT, W] = OPC4D_Tensor(ade, vertici_iniziali, ...
%                                   vertici_finali, facets)
%
% Costruisce una regola di cubatura quadridimensionale per un dominio
% poliedrico spazio-temporale in (x,y,z,tau), ottenuto dall'evoluzione
% di un dominio poliedrico tridimensionale tra le configurazioni iniziale
% e finale corrispondenti rispettivamente a tau = 0 e tau = 1.
%
% Il dominio 4D è rappresentato mediante una mesh superficiale composta
% da due sezioni tridimensionali, poste agli istanti tau = 0 e tau = 1,
% e dalle corrispondenti iperfacce laterali che collegano le facce
% triangolari della mesh iniziale con quelle della mesh finale.
%
% La regola di cubatura è costruita mediante:
%
%   - una base tensoriale di polinomi di Chebyshev in quattro variabili;
%   - una regola di Gauss-Chebyshev tensoriale sul dominio di riferimento
%     [-1,1]^4;
%   - il calcolo dei momenti shape-dependent della base sul dominio
%     spazio-temporale mediante il Teorema della Divergenza in 4D;
%   - il riscalamento della regola di riferimento sulla bounding box
%     del dominio fisico;
%   - la ricostruzione dei pesi della regola mediante la matrice di
%     Vandermonde di Chebyshev.
%
% La funzione restituisce esclusivamente i nodi e i pesi della regola
% di cubatura 4D. In particolare, per una funzione f(x,y,z,tau), 
% l'integrale può essere calcolato esternamente mediante
%
%       f_val = f(XYZT(:,1), XYZT(:,2), XYZT(:,3), XYZT(:,4));
%       I = sum(W .* f_val);
%
%**************************************************************************
%
% INPUT:
%
%   ade:
%       Grado polinomiale totale massimo della base di Chebyshev
%       utilizzata per la costruzione della regola di cubatura 4D.
%
%       La dimensione dello spazio polinomiale di grado totale <= ade
%       in quattro variabili è
%
%           N = (ade+1)(ade+2)(ade+3)(ade+4) / 24.
%
%   vertici_iniziali:
%       Matrice N_v x 3 contenente le coordinate cartesiane dei vertici
%       della mesh nella configurazione iniziale, corrispondente a
%       tau = 0.
%
%       Ogni riga rappresenta un vertice nella forma
%
%           [x_i, y_i, z_i].
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
%       triangolare.
%
%   method:
%       Stringa che specifica la regola di quadratura utilizzata per
%       l'integrazione sulle iperfacce laterali del dominio 4D.
%
%       Sono disponibili i seguenti metodi:
%
%           'DCC'  : Dunavant--Clenshaw--Curtis
%           'DGL'  : Dunavant--Gauss--Legendre
%           'GJCC' : Gauss--Jacobi--Clenshaw--Curtis
%           'GJL'  : Gauss--Jacobi--Gauss--Legendre
%
%       Le regole sono costruite come prodotto tensoriale tra una
%       quadratura sul triangolo di riferimento e una quadratura
%       unidimensionale sull'intervallo temporale [0,1].
%
%**************************************************************************
%
% OUTPUT:
%
%   XYZT:
%       Matrice N_q x 4 contenente i nodi della regola di cubatura
%       quadridimensionale.
%
%       Ogni riga contiene un nodo nella forma
%
%           [x_i, y_i, z_i, tau_i],
%
%       con tau_i appartenente all'intervallo [0,1].
%
%   W:
%       Vettore colonna N_q x 1 contenente i pesi della regola di
%       cubatura quadridimensionale associati ai nodi contenuti in XYZT.
%
%**************************************************************************
%
% NOTE:
%
%   Il dominio spazio-temporale 4D è costruito a partire dalle
%   configurazioni iniziale e finale della mesh. Le iperfacce laterali
%   sono ottenute collegando ciascuna faccia triangolare iniziale con la
%   corrispondente faccia triangolare finale, formando un prisma
%   spazio-temporale tridimensionale.
%
%   Le normali alle iperfacce del bordo 4D vengono orientate verso
%   l'esterno del dominio mediante il confronto con il baricentro
%   globale del dominio spazio-temporale.
%
%**************************************************************************
%
% Autore:
%       Edoardo Longo
%       Università degli Studi di Verona
%
% Data:
%       Settembre 2026
%
%**************************************************************************

% Percorso della directory principale del progetto
projectRoot = fileparts(fileparts(mfilename('fullpath')));

% Aggiunta della cartella contenente le funzioni ausiliarie
addpath(fullfile(projectRoot, 'src'));

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

% Numero di vertici della mesh tridimensionale
num_vertici = size(vertici_iniziali, 1);

% Inizializzazione dei vertici spazio-temporali
vertici_4D = zeros(2 * num_vertici, 4);

% Configurazione iniziale (tau = 0)
vertici_4D(1:num_vertici, 1:3) = vertici_iniziali;
vertici_4D(1:num_vertici, 4) = 0.0;

% Configurazione finale (tau = 1)
vertici_4D(num_vertici + 1:2*num_vertici, 1:3) = vertici_finali;
vertici_4D(num_vertici + 1:2*num_vertici, 4) = 1.0;


% Calcolo dell'iper-rettangolo che racchiude l'intero dominio spazio-temporale
limiti_min = min(vertici_4D, [], 1);
limiti_max = max(vertici_4D, [], 1);
bbox = [limiti_min; limiti_max];

% Calcolo dei momenti di Chebyshev sul poliedro 4D tramite il teorema della divergenza
moments_ch = chebyshev_moments_polyhedron_4D(vertici_4D, facets, ade, chebyshev_indices, bbox, method);

% Riscalamento della regola di cubatura sulla bounding box reale 4D
XYZTW_tens = scale_rule(XYZTW_tens_ref, bbox);

% Ricombinazione dei pesi di cubatura basata sui momenti del poliedro
W = (XYZTW_tens(:, 5)) .* V_ref * (moments_ch ./ coeffs);
XYZT = XYZTW_tens(:, 1:4);
