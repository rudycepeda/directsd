function [X,Y1,Y2,err] = diophsys2 ( PA1, PB1, PC1, PA2, PB2, PC2, degX )
%DIOPHSYS2 Solves a system of two special polynomial equations.
%
%     [X,Y1,Y2,ERR] = DIOPHSYS2 ( A1, B1, C1, A2, B2, C2, DEGX )
%
%   Inputs:
%     A1, B1, C1, A2, B2, C2 - polynomials
%     DEGX - desired degree of X (optional)
%
%   Outputs::
%     X, Y1, Y2 - solution, with X of minimal degree of the system
%              X*A1 + X~*B1 + Y1*C1 = 0
%              X*A2 + X~*B2 + Y2*C2 = 0
%     ERR - solution error
%
%   See also DIOPH2, DIOPHSYS.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Preliminaries
%------------------------------------------------------
        [PA1,PB1,PC1,PA2,PB2,PC2] = compat(PA1,PB1,PC1,PA2,PB2,PC2);
        polnMode = isa(PA1,'poln');
%------------------------------------------------------
%       Transform to coefficient form
%------------------------------------------------------
        if polnMode
           A1 = pol2coef(PA1);
           B1 = pol2coef(PB1);
           C1 = pol2coef(PC1);        
           A2 = pol2coef(PA2);
           B2 = pol2coef(PB2);
           C2 = pol2coef(PC2);
        else
           A1 = PA1; B1 = PB1; C1 = PC1; 
           A2 = PA2; B2 = PB2; C2 = PC2; 
        end;
%------------------------------------------------------
%       Scaling
%------------------------------------------------------
        norms = [norm(A1) norm(B1) norm(C1) norm(A2) norm(B2) norm(C2)];
        k0 = min( norm( find(norms > 0) ) );
        A1 = A1 / k0;
        B1 = B1 / k0;
        C1 = C1 / k0;
        A2 = A2 / k0;
        B2 = B2 / k0;
        C2 = C2 / k0;
%------------------------------------------------------
%       Determine degrees
%------------------------------------------------------
        if ~exist('degX','var'), 
           degX = max(deg(C1),deg(C2)); 
        end;
        deg1 = max(degX + max(deg(A1),deg(B1)), deg(C1));
        deg2 = max(degX + max(deg(A2),deg(B2)), deg(C2));
        degY1 = deg1 - deg(C1);
        degY2 = deg2 - deg(C2);        
%------------------------------------------------------
%       Form complete system of equations
%------------------------------------------------------       
        AmX1 = toep(A1,deg1+1,degX+1) + fliplr(toep(B1,deg1+1,degX+1));
        AmY1 = toep(C1,deg1+1,degY1+1);
        AmX2 = toep(A2,deg2+1,degX+1) + fliplr(toep(B2,deg2+1,degX+1));
        AmY2 = toep(C2,deg2+1,degY2+1);        
        Am = [fliplr(AmX1) fliplr(AmY1) zeros(deg1+1,degY2+1)
              fliplr(AmX2) zeros(deg2+1,degY1+1) fliplr(AmY2)];
        Cm = - Am(:,1);        
        Am(:,1) = [];
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
%------------------------------------------------------
%       Extract results
%------------------------------------------------------
        X  = [1 Bm(1:degX)'];
        Y1 = Bm(degX+1:degX+degY1+1)';
        Y2 = Bm(degX+degY1+2:end)';
        if isempty(X),  X = 0; end;
        if isempty(Y1), Y1 = 0; end;
        if isempty(Y2), Y2 = 0; end;       
%------------------------------------------------------
%       Transform to polynomial form
%------------------------------------------------------
        if polnMode
           X  = poln(X, PA1.var);
           Y1 = poln(Y1, PA1.var);
           Y2 = poln(Y2, PA1.var);
           p1 = PA1*X + PB1*recip(X) + PC1*Y1;
           p2 = PA2*X + PB2*recip(X) + PC2*Y2;
           err = norm(p1) + norm(p2);
        else
           p1 = sumpol( conv(A1,X), conv(B1,fliplr(X)) );
           p1 = sumpol( p1, conv(C1,Y1) );
           p2 = sumpol( conv(A2,X), conv(B2,fliplr(X)) );
           p2 = sumpol( p2, conv(C2,Y2) );
           err = norm(p1) + norm(p2);
        end;

%------- End of DIOPHSYS2.M --------- KYuP ----------

