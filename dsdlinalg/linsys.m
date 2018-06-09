function [flag,x] = linsys ( a, b, method, refine )
%LINSYS Solves systems of linear equations A*X=B.
%
%     X = LINSYS ( A, B )
%     X = LINSYS ( A, B, TYPE )
%     X = LINSYS ( A, B, TYPE, REFINE )
%     [FLAG,X] = LINSYS ( A, B )
%     [FLAG,X] = LINSYS ( A, B, TYPE )
%     [FLAG,X] = LINSYS ( A, B, TYPE, REFINE )
%
%   Inputs:
%     A      - rectangular system matrix
%     B      - rectangular coefficient matrix 
%     METHOD - method (optional) of
%           'qr'  - using QR-decomposition (default)
%           'svd' - using SVD-decomposition
%     REFINE - refinement (optional) of
%           'basic'  - no refinement (default)
%           'refine' - with iterative refinement
%
%   Outputs:
%     FLAG - indicator (0 - no solution, 1 - single solution, 2 - multiple solutions)
%     X    - a solution if available, or empty matrix
%
%   See also CAXBY.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check matrix dimensions
%------------------------------------------------------
        [ra,ca] = size(a);
        [rb,cb] = size(b);
        if ra ~= rb, error('Matrices A and B must have the same number of rows'); end;
%------------------------------------------------------
%       Check type
%------------------------------------------------------
        if ~exist('method','var'), 
            method = 'qr'; 
        elseif strcmp(method,'basic')  ||  strcmp(method,'refine')
            refine = method; method = 'qr';
        end;        
        if ~strcmp(method,'qr') && ~strcmp(method,'svd')
            error('Unknown solution method');
        end;       
        if ~exist('refine','var'), refine = 'basic'; end;
        if ~strcmp(refine,'basic') && ~strcmp(refine,'refine')
            error('Incorrect type of refinement');
        end;       
%########################################################        
%   Solution using QR-decomposition
%########################################################        
    tol = sqrt(eps);
    if strcmp(method, 'qr')
       %------------------------------------------------------
       %       Check solvability
       %------------------------------------------------------
        [u,ac,r] = compress(a, 'row');
        ub = u*b;
        if r < ra
           err = norm(ub(r+1:end,:));
           if err > sqrt(eps)*norm(ub)
             flag = 0;
             x = [];
             if nargout < 2, flag = x; end;
             return;
           end;
        end;    
        flag = 1;
       %------------------------------------------------------
       %   Column compression to find multiple solutions
       %------------------------------------------------------
        a0 = ac(1:r,:); 
        b0 = ub(1:r,:);
        if r < ca
           flag = 2;
           [v,a0,r] = compress(a0, 'col');
           v = fliplr(v);
           a0 = fliplr(a0);
           a0 = a0(:,1:r);
        end;
       %------------------------------------------------------
       %   Solution using QR-decomposition
       %------------------------------------------------------
        [q,R] = qr(a0);
        x = R\(R'\(a0'*b0));
        err = b0 - a0*x;   
        while 1
           dx = R\(R'\(a0'*err));
           x1 = x + dx;     
           err1 = b0 - a0*x1;
           if norm(err1) >= norm(err)  ||  norm(err1) < eps, break; end;
           x = x1;
           err = err1;
        end;
       %------------------------------------------------------
       %   Take multiple solution into account
       %------------------------------------------------------
        if flag == 2
           x = v*[x; zeros(ca-r,1)];
        end;
%########################################################        
%   Solution using SVD-decomposition
%########################################################        
    else
        [u,s,v] = svd(a);
        r = rank(s);
        ub = u' * b;
        err = norm(ub(r+1:end,:));
        if (err > tol*norm(ub))  &&  ~strcmp(refine,'refine') 
           flag = 0;
           x = [];
        else
           flag = 1;
           x = s(1:r,1:r) \ ub(1:r,:);
           if ca > r
              flag = 2;
              x = [x; zeros(ca-r,cb)];
           end;
           x = v * x;
        end;
    end;
%########################################################        
%   Refinement
%########################################################        
    if flag == 0  ||  ~strcmp(refine,'refine'), 
       if nargout < 2, flag = x; end;
       return; 
    end;
    e = b - a*x;
    while 1
       [flagRef,dx] = linsys(a, e, 'basic');
       if flagRef == 0, break; end;
       x1 = x + dx;
       e1 = b - a*x1;
       if norm(e1) >= norm(e), break; end;
       x = x1;
       e = e1;
    end;
    if norm(e) > tol*norm(b)
       flag = 0;
       x = [];
    end;
    if nargout < 2, flag = x; end;

%------- End of LINSYS.M --------- KYuP ----------           
    
