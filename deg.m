function n = deg ( p )
%DEG Degree of a polynomial or transfer function.
%
%     N = DEG ( P )
%
%   Inputs:
%     P - a polynomial or transfer function
%
%   Outputs:
%     N - degree of the polynomial or order of transfer function 
%
%   See also STRIPLZ.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Denominator degree of a transfer function
%------------------------------------------------------
        if isa(p,'lti')
           p = zpk(p);     
           [zero,pole] = zpkdata(p, 'v');
           n = max(length(zero),length(pole));
           return;     
        end;       
%------------------------------------------------------
%       Polynomial degree
%------------------------------------------------------
        n = length(p) - 1;

%------- End of DEG.M --------- KYuP ----------
