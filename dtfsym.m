function dF = dtfsym ( F, T, t, Hs, Hz )
%DTFSYM Symbolic discretization of transfer function with
%       a preceding hold
%
%       D = DTFSYM ( F, T[, Hs[, Hz]] )
%
%  Inputs:
%	  F(s) - scalar transfer function (LTI-model)
%     T    - sampling period (real value)
%     t    - time parameter (0 <= t < T)
%     H(s) - hold transfer function (symbolic) or special cases:
%            'zoh' - zero-order hold (default)
%            'foh' - first-order hold
%
%  Outputs:
%     D - discretized model as LTI-model
%
%   See also ZTM, DTFM, DTFM2.

%------------------------------------------------------
% Copyright 1999-2008 K. Polyakov 
% $Revision: 3.00 $    $Date: 09-Jan-2009 $
%------------------------------------------------------
%   Handle hold transfer function    
%------------------------------------------------------
    syms r;
    varZ = sym('z');
    if ~exist('Hs','var'), Hs = 'zoh'; end;
    if isequal(Hs,'zoh')
      Hz = (varZ-1)/varZ;
      Hs = 1/r;
    elseif isequal(Hs,'foh')    
      Hz = ((varZ-1)/varZ)^2;
      Hs = (1+r*T)/T/r^2;
    end  
    var = findsym(Hs);
    if var ~= 'r', Hs = subs(Hs, var, r); end;
%------------------------------------------------------
%   Transformation of LTI-data to symbolic form
%------------------------------------------------------
    [n,d] = tfdata(tf(F),'v');
    F = poly2sym(n,r)/poly2sym(d,r);
    FH = F*Hs;
%------------------------------------------------------
%   Find inverse Laplace transform and lattice function
%------------------------------------------------------
    tFH = ilaplace(FH);
    tFH = subs(tFH, sym('t'), T*sym('n')+t);
%------------------------------------------------------
%   Find z-transform
%------------------------------------------------------
    dF = simplify(ztrans(tFH)*Hz);
%------------------------------------------------------
%   Form final LTI-model
%------------------------------------------------------
    [n,d] = numden(dF);
    n = sym2poly(n);
    d = sym2poly(d);
    dF = zpk(minreal(tf(n,d,T)));

%------- End of DTFSYM.M --------- KYuP ----------
