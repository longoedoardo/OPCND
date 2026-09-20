function [XI_ref, ETA_ref, T_ref, W_ref] = PrismQuadrature(ade, method)

%**************************************************************************
%
% [XI_ref, ETA_ref, T_ref, W_ref] = PrismQuadrature(ade, method)
%
% Costruisce la regola di quadratura di riferimento sul prisma
% triangolare spazio-temporale, ottenuta come prodotto tensoriale.
%
%**************************************************************************
%
% INPUT:
%
%   ade:
%               Grado massimo della base polinomiale di Chebyshev per la
%               quale vengono calcolati i momenti.
%
%   method:
%               Metodo di quadratura utilizzato sul prisma di riferimento:
%
%                   'DCC'  = Dunavant--Clenshaw--Curtis
%                   'DGL'  = Dunavant--Gauss--Legendre
%                   'GJCC' = Gauss--Jacobi--Clenshaw--Curtis
%                   'GJL'  = Gauss--Jacobi--Gauss--Legendre
%
%**************************************************************************
%
% OUTPUT:
%
%   XI_ref, ETA_ref:
%               Coordinate (Nq x 1) dei nodi sul triangolo di riferimento.
%
%   T_ref:
%               Coordinata temporale tau (Nq x 1) dei nodi nella direzione
%               temporale, con tau appartenente all'intervallo [0,1].
%
%   W_ref:
%               Pesi complessivi (Nq x 1) della regola di quadratura sul
%               prisma di riferimento.
%
%**************************************************************************


if strcmp(method, "DCC")

    % Regola di quadratura sul prisma di riferimento
    [XI_ref, ETA_ref, T_ref, W_ref] = PrismaDCC(ade);

elseif strcmp(method, "DGL")

    % Regola di quadratura sul prisma di riferimento
    [XI_ref, ETA_ref, T_ref, W_ref] = PrismaDGL(ade);

elseif strcmp(method, "GJCC")

    % Regola di quadratura sul prisma di riferimento
    [XI_ref, ETA_ref, T_ref, W_ref] = PrismaGJCC(ade);

elseif strcmp(method, "GJL")

    % Regola di quadratura sul prisma di riferimento
    [XI_ref, ETA_ref, T_ref, W_ref] = PrismaGJL(ade);

else

    error('Metodo scelto non valido. Usare "GJL", "GJCC", "DGL" o "DCC".');

end


end

function [XI_ref, ETA_ref, T_ref, W_ref] = PrismaDCC(ade)

%**************************************************************************
%
% [XI_ref, ETA_ref, T_ref, W_ref] = PrismaDCC(ade)
%
% Costruisce la regola di quadratura tensoriale Dunavant--Clenshaw--Curtis
% (DCC) sul prisma triangolare di riferimento
%
%       P_ref = { (xi,eta,tau) : xi>=0, eta>=0, xi+eta<=1, 0<=tau<=1 },
%
% ottenuta come prodotto tensoriale tra:
%
%   - la regola simmetrica di Dunavant sul triangolo di riferimento
%
%   - la regola di Clenshaw--Curtis sull'intervallo temporale [0,1]
%
%**************************************************************************
%
%   INPUT:
%       ade      : (scalare intero) grado di esattezza algebrica target che si
%                   vuole raggiungere sul prisma.
%
%   OUTPUT:
%       XI_ref, ETA_ref  : (Nq x 1) coordinate baricentriche (xi,eta) dei nodi
%                           di quadratura sul triangolo di riferimento,
%                           ripetute per ciascun nodo temporale.
%       T_ref            : (Nq x 1) coordinata temporale tau in [0,1] di
%                           ciascun nodo.
%       W_ref            : (Nq x 1) peso complessivo di ciascun nodo sul
%                           prisma di riferimento, dato dal prodotto tra il
%                           peso spaziale e il peso temporale.
%
%**************************************************************************

%**************************************************************************
%
% Inizio Quadratura Spaziale
%
%**************************************************************************

rule_tri = ade+1;
n_tri = dunavant_order_num(rule_tri);
[xyrif, wrif] = dunavant_rule(rule_tri, n_tri);

r = xyrif(1, :).';
s = xyrif(2, :).';

if isrow(wrif)
    wrif = wrif.';
end
tri_weights = 0.5 * wrif;

%**************************************************************************
%
% Inizio Quadratura Temporale
%
%**************************************************************************

n1D = ade+4;
[tau_nodes, tau_weights] = clenshaw_curtis(n1D, 0, 1);



n_tau = numel(tau_nodes);

XI_ref  = kron(ones(n_tau, 1), r);
ETA_ref = kron(ones(n_tau, 1), s);
T_ref   = kron(tau_nodes, ones(n_tri, 1));
W_ref   = kron(tau_weights, tri_weights);

end

function [XI_ref, ETA_ref, T_ref, W_ref] = PrismaDGL(ade)

%**************************************************************************
%
% function [XI_ref, ETA_ref, T_ref, W_ref] = DGL(ade)
%
% Costruisce la regola di quadratura tensoriale Dunavant--Gauss--Legendre
% (DGL) sul prisma triangolare di riferimento
%
%       P_ref = { (xi,eta,tau) : xi>=0, eta>=0, xi+eta<=1, 0<=tau<=1 },
%
% ottenuta come prodotto tensoriale tra:
%
%   - una regola simmetrica di Dunavant sul triangolo di riferimento
%     (funzioni dunavant_order_num e dunavant_rule), che sfrutta le
%     simmetrie del gruppo diedrale D3;
%
%   - la stessa regola di Gauss--Legendre sull'intervallo [0,1] impiegata
%     in GJL (funzione lgwt).
%
%**************************************************************************
%
%   INPUT:
%       ade      : (scalare intero) grado di esattezza algebrica target che si
%                   vuole raggiungere sul prisma.
%
%   OUTPUT:
%       XI_ref, ETA_ref  : (Nq x 1) coordinate baricentriche (xi,eta) dei nodi
%                           di quadratura sul triangolo di riferimento,
%                           ripetute per ciascun nodo temporale.
%       T_ref            : (Nq x 1) coordinata temporale tau in [0,1] di
%                           ciascun nodo.
%       W_ref            : (Nq x 1) peso complessivo di ciascun nodo sul
%                           prisma di riferimento, dato dal prodotto tra il
%                           peso spaziale e il peso temporale.
%
%**************************************************************************

%**************************************************************************
%
% Inizio Quadratura Spaziale
%
%**************************************************************************

rule_tri = ade+1;

% dunavant_order_num restituisce il numero totale di punti fisici n_tri
% che la regola di grado rule_tri produrrà una volta espansa (somma delle
% molteplicità 1, 3 o 6 di ciascuna orbita di simmetria).
n_tri = dunavant_order_num(rule_tri);

% dunavant_rule espande i subordini tabulati nell'elenco completo dei
% nodi e dei pesi, applicando le simmetrie del gruppo D3 (permutazioni
% cicliche e riflessioni delle coordinate baricentriche). Restituisce le
% coordinate come matrice 2 x n_tri e i pesi come vettore di lunghezza n_tri.
[xyrif, wrif] = dunavant_rule(rule_tri, n_tri);

% Portiamo i nodi spaziali in formato colonna
r = xyrif(1, :).';
s = xyrif(2, :).';

if isrow(wrif)
    wrif = wrif.';
end

tri_weights = 0.5 * wrif;

%**************************************************************************
%
% Inizio Quadratura Temporale
%
%**************************************************************************

n1D = ceil((ade + 4) / 2);
[tau_nodes, tau_weights] = lgwt(n1D, 0, 1);

n_tau = numel(tau_nodes);

XI_ref  = kron(ones(n_tau, 1), r);
ETA_ref = kron(ones(n_tau, 1), s);
T_ref   = kron(tau_nodes, ones(n_tri, 1));
W_ref   = kron(tau_weights, tri_weights);

end


function [XI_ref, ETA_ref, T_ref, W_ref] = PrismaGJCC(ade)
%**************************************************************************
%
% function [XI_ref, ETA_ref, T_ref, W_ref] = PrismaGJCC(ade)
%
% Costruisce la regola di quadratura tensoriale Gauss--Jacobi--Clenshaw--
% Curtis (GJCC) sul prisma triangolare di riferimento
%
%       P_ref = { (xi,eta,tau) : xi>=0, eta>=0, xi+eta<=1, 0<=tau<=1 },
%
% ottenuta come prodotto tensoriale tra:
%
%   - la regola di Gauss--Jacobi sul triangolo di riferimento;
%
%   - una regola di Clenshaw--Curtis sull'intervallo temporale [0,1]
%     (funzione clenshaw_curtis), che utilizza i nodi di Čebyšëv--Lobatto.
%
%**************************************************************************
%
%   INPUT:
%       ade      : (scalare intero) grado di esattezza algebrica target che si
%                   vuole raggiungere sul prisma.
%
%   OUTPUT:
%       XI_ref, ETA_ref  : (Nq x 1) coordinate baricentriche (xi,eta) dei nodi
%                           di quadratura sul triangolo di riferimento,
%                           ripetute per ciascun nodo temporale.
%       T_ref            : (Nq x 1) coordinata temporale tau in [0,1] di
%                           ciascun nodo.
%       W_ref            : (Nq x 1) peso complessivo di ciascun nodo sul
%                           prisma di riferimento, dato dal prodotto tra il
%                           peso spaziale e il peso temporale.
%
%**************************************************************************

%**************************************************************************
%
% Inizio Quadratura Spaziale
%
%**************************************************************************

nGP_tri = ceil((ade + 2) / 2);
[tri_points, tri_weights] = TriangleGJQuadraturePoints(nGP_tri);
r = tri_points(:, 1);
s = tri_points(:, 2);

%**************************************************************************
%
% Inizio Quadratura Temporale
%
%**************************************************************************

n1D = ade + 4;

[tau_nodes, tau_weights] = clenshaw_curtis(n1D, 0, 1);

n_tri = numel(r);
n_tau = numel(tau_nodes);

XI_ref  = kron(ones(n_tau, 1), r);
ETA_ref = kron(ones(n_tau, 1), s);
T_ref   = kron(tau_nodes, ones(n_tri, 1));
W_ref   = kron(tau_weights, tri_weights);

end


function [XI_ref, ETA_ref, T_ref, W_ref] = PrismaGJL(ade)
%**************************************************************************
%
% [XI_ref, ETA_ref, T_ref, W_ref] = GJL(ade)
%
% Costruisce la regola di quadratura tensoriale Gauss--Jacobi--Legendre
% (GJL) sul prisma triangolare di riferimento
%
%       P_ref = { (xi,eta,tau) : xi>=0, eta>=0, xi+eta<=1, 0<=tau<=1 },
%
% ottenuta come prodotto tensoriale tra:
%
%   - una regola di Gauss--Jacobi sul triangolo di riferimento, costruita
%     tramite la trasformazione di Duffy (funzione TriangleQuadraturePoints,
%     che internamente richiama l'algoritmo di Golub--Welsch gaujac);
%
%   - una regola di Gauss--Legendre sull'intervallo temporale [0,1]
%     (funzione lgwt), che è il caso particolare alpha=beta=0 della
%     famiglia di Gauss--Jacobi.
%
%**************************************************************************
%
%   INPUT:
%       ade      : (scalare intero) grado di esattezza algebrica target che si
%                   vuole raggiungere sul prisma.
%
%   OUTPUT:
%       XI_ref, ETA_ref  : (Nq x 1) coordinate baricentriche (xi,eta) dei nodi
%                           di quadratura sul triangolo di riferimento,
%                           ripetute per ciascun nodo temporale.
%       T_ref            : (Nq x 1) coordinata temporale tau in [0,1] di
%                           ciascun nodo.
%       W_ref            : (Nq x 1) peso complessivo di ciascun nodo sul
%                           prisma di riferimento, dato dal prodotto tra il
%                           peso spaziale e il peso temporale.
%
%**************************************************************************

%**************************************************************************
%
% Inizio Quadratura Spaziale
%
%**************************************************************************

nGP_tri = ceil((ade + 2) / 2);

[tri_points, tri_weights] = TriangleGJQuadraturePoints(nGP_tri);
r = tri_points(:, 1);   % coordinata xi dei nodi spaziali
s = tri_points(:, 2);   % coordinata eta dei nodi spaziali

%**************************************************************************
%
% Inizio Quadratura Temporale
%
%**************************************************************************

n1D = ceil((ade + 4) / 2);

[tau_nodes, tau_weights] = lgwt(n1D, 0, 1);



n_tri = numel(r);
n_tau = numel(tau_nodes);

XI_ref  = kron(ones(n_tau, 1), r);            % (n_tri*n_tau x 1)
ETA_ref = kron(ones(n_tau, 1), s);            % (n_tri*n_tau x 1)
T_ref   = kron(tau_nodes, ones(n_tri, 1));    % (n_tri*n_tau x 1)
W_ref   = kron(tau_weights, tri_weights);     % (n_tri*n_tau x 1)

end

function [IntGaussP, IntGaussW] = TriangleGJQuadraturePoints(nGP)
%**************************************************************************
%
% function [IntGaussP, IntGaussW] = TriangleGJQuadraturePoints(nGP)
%
% INPUT:
% - nGP       : Scalare, numero di punti di Gauss per ogni dimensione lineare
%            (nIntGP = nGP^2 punti totali).
%
% OUTPUT:
% - IntGaussP : Matrice (nGP^2 x 2) contenente le coordinate (r,s) sul
%               triangolo di riferimento di vertici (0,0), (1,0), (0,1).
% - IntGaussW : Vettore colonna (nGP^2 x 1) dei pesi di quadratura.
%               sum(IntGaussW) = 0.5 (area del triangolo di riferimento).
%
% Mappa il quadrato [0,1]^2 nel triangolo di riferimento tramite la
% trasformazione: r = mu1, s = mu2 * (1 - mu1).
% Lo Jacobiano di questa mappa è (1 - mu1); per compensarlo si integra
% mu1 con la quadratura di Gauss-Jacobi di peso (1-t)^1 (alpha=1, beta=0)
% e mu2 con Gauss-Legendre standard (alpha=0, beta=0).
%
% Algoritmo di Golub-Welsch: nodi e pesi calcolati tramite diagonalizzazione
% della matrice di Jacobi (funzione locale gaujac, identica a quella usata
% in TetrahedronQuadraturePoints).
%
%**************************************************************************

nIntGP = nGP^2;
IntGaussP = zeros(nIntGP, 2);
IntGaussW = zeros(nIntGP, 1);

% Ottenimento posizioni (mu) e pesi (A) tramite Gauss-Jacobi
[mu1, A1] = gaujac(nGP, 1.0, 0.0);
[mu2, A2] = gaujac(nGP, 0.0, 0.0);

% Shift e rescale dall'intervallo [-1, 1] a [0, 1]
mu1 = 0.5 * mu1 + 0.5;  A1 = (0.5^2) * A1;
mu2 = 0.5 * mu2 + 0.5;  A2 = (0.5^1) * A2;

% Generazione dei punti tramite trasformazione di Duffy
iIntGP = 1;
for i = 1:nGP
    for j = 1:nGP
        m1 = mu1(i);
        m2 = mu2(j);

        IntGaussP(iIntGP, 1) = m1;
        IntGaussP(iIntGP, 2) = m2 * (1.0 - m1);

        IntGaussW(iIntGP) = A1(i) * A2(j);

        iIntGP = iIntGP + 1;
    end
end

end

function [x, w] = gaujac(n, alpha, beta)
% Calcola nodi (x) e pesi (w) di Gauss-Jacobi (formato colonna)
if n == 0, x = []; w = []; return; end
if n == 1
    x = (beta - alpha) / (alpha + beta + 2);
    w = 2^(alpha + beta + 1) * gamma(alpha + 1) * gamma(beta + 1) / gamma(alpha + beta + 2);
    return;
end

% Indici
i = (1:n-1)';
abi = alpha + beta + 2*i;

% Coefficienti diagonali (aa)
aa = zeros(n, 1);
aa(1) = (beta - alpha) / (alpha + beta + 2);
aa(2:n) = (beta^2 - alpha^2) ./ (abi .* (abi + 2));

% Coefficienti sub-diagonali (bb)
bb = zeros(n-1, 1);
bb(1) = 2 * sqrt(1 * (1 + alpha) * (1 + beta) / ((alpha + beta + 2)^2 * (alpha + beta + 3)));
bb(2:n-1) = 2 ./ (abi(1:end-1) + 2) .* sqrt(i(2:end) .* (i(2:end) + alpha) .* ...
    (i(2:end) + beta) .* (i(2:end) + alpha + beta) ./ ...
    ((abi(1:end-1) + 1) .* (abi(1:end-1) + 3)));

% Costruzione matrice di Jacobi (n x n)
J = diag(aa) + diag(bb, 1) + diag(bb, -1);

% Autovalori e pesi
[V, D] = eig(J);
[x, idx] = sort(diag(D));
V = V(:, idx);

% Fattore di normalizzazione
factor = 2^(alpha + beta + 1) * gamma(alpha + 1) * gamma(beta + 1) / gamma(alpha + beta + 2);
w = factor * (V(1, :)').^2;
end