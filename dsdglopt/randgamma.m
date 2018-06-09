function rg = randgamma ( a, rows, cols )
%RANDGAMMA Gamma-distributed random numbers.
%
%     X = RANDGAMMA ( A, R, C )
%
%   Inputs: 
%     A - parameter
%     R - rows
%     C - columns
%
%   Outputs:
%     X - array of gamma-distributed values
%
%   See also RANDBETA.

% [1] Marsaglia, G. and Tsang, W.W. (2000) "A Simple Method for Generating
%     Gamma Variables", ACM Trans. Math. Soft. 26(3):363-372.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%%       Check data
%------------------------------------------------------
        if ~exist('rows','var'), rows = 1; end;
        if ~exist('cols','var'), cols = rows; end;
        rg = zeros(rows, cols);
%------------------------------------------------------
%%       Case a < 1
%------------------------------------------------------
        if a == 0, return; end;
        if a <= 1,
          rg = randgamma(a+1,rows,cols).*rand(rows,cols).^(1/a);
          return;  
        end;
%------------------------------------------------------
%%       Initialization
%------------------------------------------------------
        d = a - 1./3.; 
        c = 1. / sqrt(9.*d);
%------------------------------------------------------
%%       Main algorithm
%------------------------------------------------------
        for i=1:rows*cols
          while 1
            while 1 
              x = randn;     % normal distribution
              v = 1. + c*x;
              if v > 0, break; end;
            end;
            v = v*v*v; 
            u = rand;        % uniform distribution 
            rg(i) = d*v;
            if u < 1-0.0331*(x*x)*(x*x), break; end;
            if log(u) < 0.5*x*x+d*(1-v+log(v)), break; end;
           end; 
        end;

%------- End of RANDGAMMA.M --------- KYuP ----------
