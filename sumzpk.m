function s = sumzpk ( varargin )
%SUMZPK Reliable summation of zpk-models with common poles.
%
%     S = SUMZPK ( F1, F2, ... )
%
%   Inputs:
%     F1, F1, ... - LTI-models or constants
%
%   Outputs:
%     S - the sum of these operands
%
%   See also PLUS.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Preliminaries
%------------------------------------------------------
        i = 1;
        s = [];
        tol = sqrt(eps);        
        T = -1;
%------------------------------------------------------
%       Main loop
%------------------------------------------------------
        while i <= nargin
          f = varargin{i};  
          if T < 0  &&  isa(f,'lti')
              T = f.Ts; 
          end;
          if ~isa(f,'lti');%%Added by RCG to make sure the variable is kept
              f = sdzpk([],[],f,T);
          end          
          f = stripzpk(sdzpk(f));          
          if isempty(s), 
             [rows,cols] = size(f);
             s = sdzpk(zeros(rows,cols));
             s.Ts = T;
          end;
          if any(size(f) ~= [rows cols])
             error('All the summands must be of the same size'); 
          end;
          for r=1:rows
             for c=1:cols          
               [sz,sp] = zpkdata(s(r,c), 'v');
               ind = find(abs(imag(sz)) < tol);
               sz(ind) = real(sz(ind));
               ind = find(abs(imag(sp)) < tol);
               sp(ind) = real(sp(ind));

               f(r,c) = minreal(f(r,c));
               [fz,fp] = zpkdata(f(r,c), 'v');
               ind = find(abs(imag(fz)) < tol);
               fz(ind) = real(fz(ind));
               ind = find(abs(imag(fp)) < tol);
               fp(ind) = real(fp(ind));

               [fz,zcommon] = others(fz, sz, tol);
               sz = others(sz, zcommon, tol);
               zcommon = roots(real(poly(zcommon)));
               s.z{r,c} = sz;
               f.z{r,c} = fz;

               [fp,pcommon] = others(fp, sp, tol);
               sp = others(sp, pcommon, tol);
               s.p{r,c} = sp;
               f.p{r,c} = fp;

               %s(r,c) = s(r,c) + f(r,c);
               [ns,ds] = tf2nd ( s(r,c) );
               [nf,df] = tf2nd ( f(r,c) );
               n = ns*df + nf*ds;
               d = ds * df;
               s(r,c) = sdzpk(n.z, d.z, n.coef(1), T);
               
               s.z{r,c} = [s.z{r,c}; zcommon];
               s.p{r,c} = [s.p{r,c}; pcommon];
             end;
          end;         
          s = minreal(s);
          i = i + 1;
        end;
%------------------------------------------------------
%       Remove zeros near infinity by truncating
%       leading coefficients of the numerator
%------------------------------------------------------
        s = stripzpk(s);

%------- End of SUMZPK.M --------- KYuP ----------
