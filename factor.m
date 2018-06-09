function [fs,fu,f0] = factor ( f, type )
%FACTOR Stable-unstable factorization of a polynomial.
%
%    [FS,FU,F0] = FACTOR ( F, TYPE )
%
%  Inputs:
%    S - a continuous-time or discrete LTI system 
%    TYPE - factorization type of 
%        's' (default for continuous-time models)
%        'd' (default for discrete-time models)
%        'z' for z-plane
%
%  Outputs:
%    FS  - stable factor
%    FU  - unstable factor
%    F0  - neutral factor
%
%   See also SFACTOR, SFACTFFT. 

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Special case: array of coefficients
%------------------------------------------------------
        if isnumeric(f)
          if ~exist('type','var'), 
               [fs,fu,f0] = factor ( poln(f) );
          else [fs,fu,f0] = factor ( poln(f), type ); 
          end;
          fs = fs.coef;
          fu = fu.coef;
          f0 = f0.coef;
          return;  
        end;
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if ~exist('type','var'), 
           if isdt(f), type = 'd'; else type = 's'; end;
        end;
        if ~strcmp(type,'s')  &&  ~strcmp(type,'z')  &&  ~strcmp(type,'d')
            error('Unknown factorization type: ''%s''', type);
        end;         
%------------------------------------------------------
%       Check for zero function
%------------------------------------------------------
        if abs(f.k) < eps, 
          fs = poln(0, f.var);
          fu = poln(1, f.var);
          f0 = fu;
          return;
        end;
%------------------------------------------------------
%       Factorize a polynomial
%------------------------------------------------------
        f0 = f.z; 
        switch type
          case 's',   
            inds = find(real(f0) < 0);  
            indu = find(real(f0) > 0);  
          case 'z',   
            inds = find(abs(f0) < 1);  
            indu = find(abs(f0) > 1);  
          case 'd',   
            inds = find(abs(f0) > 1);  
            indu = find(abs(f0) < 1);  
        end;
        fs = f0(inds);
        fu = f0(indu);
        f0([inds; indu]) = [];
%------------------------------------------------------
%       Form results
%------------------------------------------------------
        fs = f.k * poln ( fs, ['r' f.var] );
        fu = poln ( fu, ['r' f.var], f.shift );
        f0 = poln ( f0, ['r' f.var] );
        
%------- End of FACTOR.M --------- KYuP ----------
