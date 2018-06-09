function G = gcd ( A, B, varargin )
%GCD Monic greatest common divisor of two or more polynomials.
%       
%     G = GCD ( A, B, ... )
%
%   Inputs:
%     A, B, ... - polynomials
%
%   Outputs:
%     G - greatest common divisor of all polynomials
%
%   See also COPRIME, TRIPLE. 

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check syntax
%------------------------------------------------------
        [A,B] = compat ( A, B );
        [xx,xx,G] = coprime ( A, B );
%------------------------------------------------------
%       Multiple polynomial case
%------------------------------------------------------
        for i=1:nargin-2
          [xx,C] = compat(G, varargin{i});  
          [xx,xx,G] = coprime ( G, C );  
        end;
        
%------- End of GCD.M --------- KYuP ----------           

