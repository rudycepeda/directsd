function [X,Y,err] = dioph2 ( PA, PB, PC, degX )
%DIOPH2 Solves a special polynomial equations X*A + X~*B + Y*C = 0.
%
%     [X,Y,ERR] = DIOPH2 ( A, B, C, DEGX )
%
%   Inputs:
%     A, B, C - polynomials
%     DEGX - desired degree of X (optional)
%
%   Outputs::
%     X, Y - polynomial solution
%     ERR  - solution error
%
%   See also DIOPH, DIOPHSYS, DIOPHSYS2.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Scaling
%------------------------------------------------------
polnMode = isa(PA,'poln');
%------------------------------------------------------
%       Transform to coefficient form
%------------------------------------------------------
if polnMode
    [PA,PB,PC] = compat(PA, PB, PC);
    A = pol2coef(PA);
    B = pol2coef(PB);
    C = pol2coef(PC);
else
    A = PA; B = PB; C = PC;
end;
%------------------------------------------------------
%       Scaling
%------------------------------------------------------
%norms = [norm(A) norm(B) norm(C)];
%k0 = min( norm( find(norms > 0) ) );
k0 = 1;
A = A / k0;
B = B / k0;
C = C / k0;
%------------------------------------------------------
%       Determine degrees
%------------------------------------------------------
if ~exist('degX','var'), degX = deg(C); end;
degY = max(deg(A),deg(B));
degAll = max(degX + max(deg(A),deg(B)), degY+deg(C));
%------------------------------------------------------
%       Form complete system of equations
%------------------------------------------------------
AmX = toep(A,degAll+1,degX+1) + fliplr(toep(B,degAll+1,degX+1));
AmY = toep(C,degAll+1,degY+1);
Am = [AmY AmX];
%-----------------------------------------------
%   The highest coefficient of X equals 1
%-----------------------------------------------
Cm = - Am(:,end);
Am(:,end) = [];
%------------------------------------------------------
%       Solve equations
%------------------------------------------------------
%condA = cond(Am);
Bm = linsys ( Am, Cm, 'svd', 'refine' );
if isempty(Bm) %-- It cannot be worse!
    ws = warning; warning off;
    Bm = Am \ Cm;
    warning(ws);
end;
Bm = findzero(Bm);
%------------------------------------------------------
%       Extract results
%------------------------------------------------------
Y = fliplr(Bm(1:degY+1).');
X = fliplr([Bm(degY+2:end).' 1]);
if isempty(X), X = 0; end;
if isempty(Y), Y = 0; end;
Y = striplz(Y);
%------------------------------------------------------
%       Transform to polynomial form
%------------------------------------------------------
if polnMode
    X = poln(X, PA.var);
    Y = poln(Y, PA.var);
    p = PA*X + PB*recip(X) + PC*Y;
    err = norm(p);
else
    p = sumpol( conv(A,X), conv(B,fliplr(X)) );
    p = sumpol( p, conv(C,Y) );
    err = norm(p);
end;

%------- End of DIOPH2.M --------- KYuP ----------

