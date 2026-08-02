function V = dCHEBVAND(deg, X, duples)
%**************************************************************************
% 
% V = dCHEBVAND(deg, X, duples)
%
% La funzione calcola la matrice di Chebyshev-Vandermonde per il grado "deg"
% sui punti "X". Modificata e ottimizzata per 4D (vettorializzata).
%
%**************************************************************************

[m, d] = size(X); % d sarà 4, m è il numero totale di nodi spaziotemporali
N = size(duples, 1); % Dimensione dello spazio polinomiale 4D

% Crea automaticamente un box [-1, 1] per ogni dimensione presente (4D).
dbox = [-ones(1, d); ones(1, d)]; 

% Evitiamo il ciclo riga per riga. Lavoriamo direttamente sull'intera matrice X.
map = zeros(m, d);
for i = 1 : d
    min_val = dbox(1, i);
    max_val = dbox(2, i);
    % Vettorializzazione: calcola il map per tutti gli m punti in un colpo solo
    map(:, i) = (2.0 * X(:, i) - max_val - min_val) / (max_val - min_val);
end

V = ones(m, N); % Inizializzazione a 1 (neutro moltiplicativo)

% Ciclo sulle dimensioni (x, y, z, tau)
for k = 1 : d
    
    % Calcolo dei Polinomi di Chebyshev 1D per l'intera colonna k-esima
    % Invece di calcolarlo punto per punto, calcoliamo i vettori colonna
    T_dim = zeros(m, deg + 1);
    
    x_val = map(:, k); % Vettore colonna di tutti gli m punti per la dim k
    
    % Grado 0: T0(x) = 1 per tutti i punti
    T_dim(:, 1) = 1.0;
    
    if deg >= 1
        % Grado 1: T1(x) = x per tutti i punti
        T_dim(:, 2) = x_val;
        
        % Gradi successivi con formula di ricorrenza T_n+1 = 2*x*T_n - T_n-1
        for g = 2 : deg
            % Moltiplicazione elemento per elemento (.*) sui vettori!
            T_dim(:, g+1) = 2.0 * x_val .* T_dim(:, g) - T_dim(:, g-1);
        end
    end
    % Moltiplichiamo il contributo della dimensione k ai vari monomi multi-indice
    for col = 1 : N
        grado_richiesto = duples(col, k);
        
        % Invece del ciclo for riga = 1:m, moltiplichiamo l'intera colonna
        % V(:, col) in una singola operazione vettoriale.
        V(:, col) = V(:, col) .* T_dim(:, grado_richiesto + 1);
    end
    
    % In 4 dimensioni, il termine al denominatore e' pi^4.
    % V(:,k) = V(:,k) * sqrt( 2^(sum(duples(k,:) > 0)) / (pi^d) ); 
end

end