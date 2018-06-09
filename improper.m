function [P,R0] = improper ( R, type )
%IMPROPER Separate improper (polynomial) part of a rational matrix.
%
%     [P,R0] = IMPROPER ( R, TYPE )
%
%   Inputs:
%     R(x) - rational matrix in the form of LTI model
%     TYPE - properness type of
%        'p'  - R0 is proper but may contain a direct term, 
%               P is a polynomial without a free term
%        'sp' - R0 is strictly proper,  
%               P is a polynomial with a free term (default)
%        'symm' - the same as 'sp', for symmetric spectral densities
%
%   Outputs:
%     P(x)  - polynomial (improper) part
%     R0(x) - proper or strictly proper part
%
%   See also TF2ND.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Determine type of result
%------------------------------------------------------       
        if ~exist('type','var'), type = 'sp'; end;
        switch type
          case {'p','sp','symm'},
          otherwise, error('Unknown properness type ''%s''',type);
        end;
%------------------------------------------------------
%       Special case: constant matrix
%------------------------------------------------------       
        if isnumeric(R)
           if isequal(type,'p')         
                P = R*0; R0 = zpk(R); 
           else P = R;   R0 = zpk(R*0); 
           end;
           return;
        end;
        cls = class(R);
        if isequal(type,'symm') && isdt(R)
          warning('Symmetrization is senseless for discrete-time models');  
          type = 'sp';
        end;
%------------------------------------------------------
%       Handle MIMO case
%------------------------------------------------------       
        [rows,cols] = size(R);
        if rows > 1  ||  cols > 1
          P = cell(rows,cols);
          R0 = tf(zeros(rows,cols)); R0.Ts = R.Ts;
          if isequal(cls,'ss'),  R0 = ss(R0); end; 
          if isequal(cls,'zpk'), R0 = zpk(R0); end; 
          for r=1:rows
            for c=1:cols
              [Pi,R0i] = improper(R(r,c), type);  
              P{r,c} = Pi;
              R0(r,c) = R0i;
            end;
          end;
          return;  
        end;
%------------------------------------------------------
%       Special case: state-space
%------------------------------------------------------        
        if isequal(cls,'ss')
          switch type
            case 'p',  P = 0; R0 = R; return;
            case 'sp', P = R.d; R0 = R; R0.d = R0.d*0; return;
            otherwise, R = zpk(R); % symmetrization follow
          end;
        end;
%------------------------------------------------------
%       General case: improper function
%------------------------------------------------------       
        P = 0;
        var = R.var;
        [n0,d0] = tfdata(tf(R), 'v');
        [z,p,xx,T] = zpkdata(zpk(R), 'v');
%------------------------------------------------------
%       Special case: no changes need
%------------------------------------------------------
        R0 = R;
        switch type
          case 'p',  if length(z) <= length(p), return; end;
          case 'sp', if length(z) < length(p), return; end;
        end;
%------------------------------------------------------
%       Divide numerator by the denominator using
%       a polynomial equation
%------------------------------------------------------
        if norm(n0) < eps
           n = 0;
        else   
           n0 = striplz(n0);
           d0 = striplz(d0);
           [n,P] = dioph(1, d0, n0);                     
           if isequal(type,'symm'), 
             n(end-1:-2:1) = 0;
             n = striplz(n);
           end;
           n = findzero(n);
        end;
        R0 = zpk( roots(n), p, n(1)/d0(1), T );
        R0.var = var;
%------------------------------------------------------
%       Restore initial form
%------------------------------------------------------          
        if isequal(type,'p')
           k0 = P(end); P(end) = 0;
           R0 = minreal(R0 + k0);
        end;
        if isequal(cls, 'tf'), R0 = tf(R0); end;
        if isequal(cls, 'ss'), R0 = ss(R0); end;
        P = poln(P, var);

%------- End of IMPROPER.M --------- KYuP ----------
