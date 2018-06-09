function c = sumpol2 (a, b, varargin)
%SUMPOL2 Find sum of two symmetric quasipolynomials.
%
%     C = SUMPOL2 ( A, B )
%
%   Inputs:
%     A, B - symmetric polynomials as coefficient arrays
%
%   Outputs:
%     C - the sum as an array
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Multiple polynomials version
%------------------------------------------------------
        if nargin > 2, 
           c = sumpol2(a, b);
           for i=3:nargin
             c = sumpol2(c, varargin{i-2});  
           end;
           return;
        end;
%------------------------------------------------------
%       Simplest version for two quasipolynomials
%------------------------------------------------------
        la = length(a);
        lb = length(b);        
        if mod(la,2)==0  ||  mod(lb,2)==0
           error('Incorrect symmetric quasipolynomials');
        end;        
        n = (la - lb) / 2;
        da = zeros(1,-n);
        db = zeros(1,n);
        c = [da a da] + [db b db];

%------- End of SUMPOL2.M --------- KYuP ----------
