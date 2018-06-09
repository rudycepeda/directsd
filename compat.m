function varargout = compat ( varargin )
%COMPAT Make compatible polynomials with the same variable.
%
%       [A1,B1,...] = COMPAT ( A, B, ... )
%
%   Inputs:
%     A, B, ... - polynomials or vectors of coefficients
%
%   Outputs:
%     A1, B1, ... - polynomials with the same variable
%
%   See also S, P, D, Z, Q.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Equalize variable
%------------------------------------------------------
var = [];
for i=1:nargin
    if isa(varargin{i},'poln'),
        var = varargin{i}.var;
        break;
    end;
end;
%------------------------------------------------------
%       If there is at least one polynomial, ...
%------------------------------------------------------
varargout = varargin;
for i=1:nargin
        varargout{i} = poln(varargin{i},var);
end
%------- End of COMPAT.M --------- KYuP ----------
