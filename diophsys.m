function [X,Y1,Y2,err,condA] = diophsys ( PA1, PB1, PC1, PA2, PB2, PC2 )
%DIOPHSYS Solves a system of two Diophantine polynomial equations.
%
%     [X,Y1,Y2,ERR] = DIOPHSYS ( A1, B1, C1, A2, B2, C2 )
%
%   Inputs:
%     A1, B1, C1, A2, B2, C2 - polynomials
%
%   Note: Polynomial A2*C1-A1*C2 should be divisible by GCD(B1,B2)
%
%   Outputs:
%     X, Y1, Y2 - solution, with X of minimal degree of the system
%                   X*A1 + Y1*B1 = C1
%                   X*A2 + Y2*B2 = C2
%     ERR - solution error
%
%   See also DIOPH, DIOPHSYS2.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Scaling
%------------------------------------------------------
        [PA1,PB1,PC1,PA2,PB2,PC2] = compat(PA1,PB1,PC1,PA2,PB2,PC2);
        polnMode = isa(PA1,'poln');
%------------------------------------------------------
%       Transform to coefficient form
%------------------------------------------------------
        [B1r,B2r,V] = coprime(PB1, PB2);
        if polnMode
           A1 = pol2coef(PA1);
           B1 = pol2coef(PB1);
           C1 = pol2coef(PC1);        
           A2 = pol2coef(PA2);
           B2 = pol2coef(PB2);
           C2 = pol2coef(PC2);
           V  = pol2coef(V);
           B  = pol2coef(PB1*B2r);            
        else
           A1 = PA1; B1 = PB1; C1 = PC1; 
           A2 = PA2; B2 = PB2; C2 = PC2; 
           B  = conv(PB1, B2r);            
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
%       Check solvability
%------------------------------------------------------        
        P = sumpol(conv(A2,C1), -conv(A1,C2));
        if norm(P) > sqrt(eps)
           [P0,rem] = deconv(P, V);
           P0 = striplz(P0);
           err = norm(rem) / norm(P);
           if err > 1.e-6*norm(P)
              warning('Polynomial A2*C1-A1*C2 should be divisible by V (err=%g)',err);
           end;
        end;
%------------------------------------------------------
%       Find polynomial degrees
%------------------------------------------------------
        degX = deg(B) - 1;
        deg1 = max(degX + max(deg(A1),deg(B1)), deg(C1));
        deg2 = max(degX + max(deg(A2),deg(B2)), deg(C2));
        degY1 = deg(A1) + deg(B2r) - 1;
        degY2 = deg(A2) + deg(B1r) - 1;   
        degC12 = deg(P0);
        if degY1+deg(A2)+deg(B1) < degC12
           degY2 = degC12 - deg(A1) - deg(B2r);
        end;
        if degX+deg(A1) < deg(C1)
           degY1 = deg(C1) - deg(B1);
        end;
        if degX+deg(A2) < deg(C2)
           degY2 = deg(C2) - deg(B2);
        end;
%------------------------------------------------------
%       Solve complete system of equations
%------------------------------------------------------       
        m  = degX + 1;
        n1 = degY1 + 1;
        n2 = degY2 + 1;
        Am1X  = toep(A1, deg1+1, m);
        Am1Y1 = toep(B1, deg1+1, n1);
        Am2X  = toep(A2, deg2+1, m);
        Am2Y2 = toep(B2, deg2+1, n2);
        Am = [Am1X Am1Y1 zeros(deg1+1,n2)
              Am2X zeros(deg2+1,n1) Am2Y2];
        Cm = [toep(C1,deg1+1,1); toep(C2,deg2+1,1)];        
%------------------------------------------------------
%       Solve equations
%------------------------------------------------------
        condA = cond(Am);
        Bm = linsys ( Am, Cm, 'svd', 'refine' );       
        if isempty(Bm) %-- It cannot be worse!
           ws = warning; warning off;
           Bm = Am \ Cm;
           warning(ws); 
        end;
%------------------------------------------------------
%       Extract results
%------------------------------------------------------
        X  = fliplr ( Bm(1:m)' );
        Y1 = fliplr ( Bm(m+1:m+n1)' );
        Y2 = fliplr ( Bm(m+n1+1:m+n1+n2)' );
        if isempty(X),  X = 0; end;
        if isempty(Y1), Y1 = 0; end;
        if isempty(Y2), Y2 = 0; end;       
%------------------------------------------------------
%       Transform to polynomial form
%------------------------------------------------------
        if polnMode
           X = poln(X, PA1.var);
           Y1 = poln(Y1, PA1.var);
           Y2 = poln(Y2, PA1.var);
           p1 = PA1*X + PB1*Y1;
           p2 = PA2*X + PB2*Y2;
           err = norm(p1-PC1) + norm(p2-PC2);
        else
           p1 = sumpol( conv(A1,X), conv(B1,Y1) );
           p2 = sumpol( conv(A2,X), conv(B2,Y2) );
           err = norm( sumpol(p1, -C1) ) + norm( sumpol(p2, -C2) );
        end;

%------- End of DIOPHSYS.M --------- KYuP ----------

