function T=chebpolys(deg,x)

%--------------------------------------------------------------------------
% Object:
% This routine computes the Chebyshev-Vandermonde matrix on the real line
% by recurrence.
%--------------------------------------------------------------------------
% Input:
% deg: maximum polynomial degree
% x: 1-column array of abscissas
%--------------------------------------------------------------------------
% Output:
% T: Chebyshev-Vandermonde matrix at x, T(i,j+1)=T_j(x_i), j=0,...,deg.
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

T = zeros(length(x), deg+1);
T(:,1) = 1;

if deg >= 1
    T(:,2) = x;
end

for j = 2:deg
    T(:,j+1) = 2*x.*T(:,j) - T(:,j-1);
end