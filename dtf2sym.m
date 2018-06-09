function dF = dtf2sym ( F, T, Hs, Hz )
%DTF2SYM Symbolic "double" discretization of transfer 
%        function with a preceding hold
%
%       D = DTF2SYM ( F, T[, Hs[, Hz]] )
%
%  Inputs:
%	  F(s) - scalar transfer function (LTI-model)
%     T    - sampling period (real value)
%     Hs   - s-part of the hold transfer function (symbolic) 
%            or special cases:
%             'zoh' - zero-order hold (default)
%             'foh' - first-order hold
%     Hz   - z-part of the hold transfer function (symbolic)
%
%  Outputs:
%     D - discretized model as LTI-model
%
%   See also DTFSYM, ZTM, DTFM, DTFM2.

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
    end; 
    var = findsym(Hs);
    if var ~= 'r', Hs = subs(Hs, var, r); end;
%------------------------------------------------------
%   Transformation of LTI-data to symbolic form
%------------------------------------------------------
    [n,d] = tfdata(tf(F),'v');
    F = poly2sym(n,r)/poly2sym(d,r);
    FH = F*Hs;
    iFH = subs(FH, r, -r);
    iHz = subs(Hz, sym('z'), sym('1/z'));
%------------------------------------------------------
%   Find inverse Laplace transform and lattice function
%------------------------------------------------------
    tFH2 = ilaplace(FH*iFH);
    tFH2 = subs(tFH2, sym('t'), T*sym('n'));
%------------------------------------------------------
%   Find z-transform
%------------------------------------------------------
    dF = simplify(ztrans(tFH2)*Hz*iHz);
%------------------------------------------------------
%   Form final LTI-model
%------------------------------------------------------
    [n,d] = numden(dF);
    n = sym2poly(n);
    d = sym2poly(d);
    dF = zpk(minreal(tf(n,d,T)));

%------- End of DTF2SYM.M --------- KYuP ----------
