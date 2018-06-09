function [w,r] = rectchol ( a, m )
%RECTCHOL Cholesky-like rectangular form of a square matrix.
%         The matrix may have incomplete rank.
%
%     W = RECTCHOL ( A, M )
%     W = RECTCHOL ( A, TOL )
%
%   Inputs:
%  	  A - square N x N nonnegative definite matrix
%  	  M - row dimension of the desired matrix factor 
%  	  TOL - tolerance for rank evaluation
%
%   Outputs:
%	  W - an M x N matrix such that W' * W = A
%
%   See also CHOL.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 05-Jan-2007 $
%------------------------------------------------------
%       Check input data
%------------------------------------------------------
    	[n,c] = size(a);
    	if n ~= c, error('Matrix is not square');end;
        tol = eps*norm(a);
        if exist('m','var')
           if ~isint(m) 
              tol = m; m = -1; 
           end;
        else
           m = -1;
        end;
%------------------------------------------------------
%       Rank evaluation
%------------------------------------------------------
        r = rank(a, tol);
        if m > 0, r = min(m,r); end; % need truncation
%------------------------------------------------------
%       Symmetrize using SQRTM
%------------------------------------------------------
        ws = warning;    
        warning off MATLAB:sqrtm:SingularMatrix;
        a = real(sqrtm(a*a')); % --- symmetrize
        warning(ws);    
%	if norm(a-a') > 1.e-4*norm(a)
%	   warning('Matrix is not symmetric');
%	end;
%	a = (a + a') / 2;       % regularization    
%------------------------------------------------------
%       Find symmetric SVD decomposition
%------------------------------------------------------
        if any(isnan(a)), keyboard; end;
    	[~,S0,V] = svd ( a );
       %-----------------------------------------------
       %       Check positive semidefiniteness
       %-----------------------------------------------
        %if norm(V'*a*V-S0) > sqrt(eps)*norm(S0), 
        %  warning('Matrix might be not positive semidefinite'); 
        %end;
%------------------------------------------------------
%       Rank evaluation
%------------------------------------------------------
        S = S0(1:r,:);
        S(1:r,1:r) = sqrtm(S(1:r,1:r));
        w = S * V';
        w = w';
%------------------------------------------------------
%       Expand matrix is necessary
%------------------------------------------------------
        if m > r, w = [w zeros(n,m-r)]; end;
%------------------------------------------------------
%   	Check correctness of solution
%------------------------------------------------------
        %if norm(w*w'-a) > sqrt(eps)*norm(S0), 
        % warning('Matrix might be not positive semidefinite'); 
        %end;
    	w = w';

%------- End of RECTCHOL.M --------- KYuP ----------
