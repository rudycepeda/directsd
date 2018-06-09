function [a,n] = striplz (a,tol)
%STRIPLZ Strip leading zeros of a vector.
%
%   A = STRIPLZ(A,tol) removes all the leading elements of vector A whose
%   magnitude is less than tol. If no tolerance is given, its value is set
%   to eps by default.
%
%   [A,N] = STRIPLZ(A,tol) also returns N, the number of leading elements
%   removed.
%
%   See also DEG.

%------------------------------------------------------
% Copyright 1999-2014
% Originally written by K. Yu. Polyakov
% Revised by Rudy Cepeda-Gomez
% Revision: 4.00
% Date: 14-Jul-2014
%------------------------------------------------------

%Check data
if ~isvector(a)
    error('A must be a vector')
end
if ~exist('tol','var')
    tol = eps;
end
%Making sure a is a row vector.
a=a(:)';

%Special case - zero polynomial
if norm(a) < max(tol*tol,100*eps)
    n = size(a,2);
    a = zeros(size(a,1),1);
    return;
end

%General case
n = 0;
c = length(a);
while abs(a(1))<tol*norm(a(2:end))
    a=a(2:end);
    n=n+1;
    if n == c
        break
    end
end
