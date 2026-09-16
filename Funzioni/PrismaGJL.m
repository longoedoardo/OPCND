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

nGP_tri = ceil((ade + 1) / 2);

[tri_points, tri_weights] = TriangleQuadraturePoints(nGP_tri);
r = tri_points(:, 1);   % coordinata xi dei nodi spaziali
s = tri_points(:, 2);   % coordinata eta dei nodi spaziali

%**************************************************************************
%
% Inizio Quadratura Temporale
%
%**************************************************************************

n1D = ceil((ade + 1) / 2);

[tau_nodes, tau_weights] = lgwt(n1D, 0, 1);



n_tri = numel(r);
n_tau = numel(tau_nodes);

XI_ref  = kron(ones(n_tau, 1), r);            % (n_tri*n_tau x 1)
ETA_ref = kron(ones(n_tau, 1), s);            % (n_tri*n_tau x 1)
T_ref   = kron(tau_nodes, ones(n_tri, 1));    % (n_tri*n_tau x 1)
W_ref   = kron(tau_weights, tri_weights);     % (n_tri*n_tau x 1)

end