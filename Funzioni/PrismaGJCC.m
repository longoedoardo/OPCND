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

nGP_tri = ceil((ade + 1) / 2);
[tri_points, tri_weights] = TriangleQuadraturePoints(nGP_tri);
r = tri_points(:, 1);
s = tri_points(:, 2);

%**************************************************************************
%
% Inizio Quadratura Spaziale
%
%**************************************************************************

n1D = ade+1;

[tau_nodes, tau_weights] = clenshaw_curtis(n1D, 0, 1);

n_tri = numel(r);
n_tau = numel(tau_nodes);

XI_ref  = kron(ones(n_tau, 1), r);
ETA_ref = kron(ones(n_tau, 1), s);
T_ref   = kron(tau_nodes, ones(n_tri, 1));
W_ref   = kron(tau_weights, tri_weights);

end