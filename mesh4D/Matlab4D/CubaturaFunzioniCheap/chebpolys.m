function T=chebpolys(deg,x)
%--------------------------------------------------------------------------
% Object:
% This routine computes the Chebyshev-Vandermonde matrix on the real line
% by recurrence.
%--------------------------------------------------------------------------
T=zeros(length(x),deg+1);
t0=ones(length(x),1); T(:,1)=t0;

if deg > 0 % <- FIX DI ROBUSTEZZA: Evita il crash se si richiede il grado 0
    t1=x; T(:,2)=t1;
    % costruisce la matrice di vandermonde usando come base polinomiale quella
    % di chebyshev data dalla formula 2xT_n -T_{n-1}.
    for j=2:deg
        t2=2*x.*t1-t0;
        T(:,j+1)=t2;
        t0=t1;
        t1=t2;
    end
end
end