function [r,err] = remove (r,x,tol)
%REMOVE remove nearest item(s) from an array
%
%     [R0,ERR] = REMOVE ( R, X )  
%     [R0,ERR] = REMOVE ( R, X, TOL )  
%
%   Inputs:
%     R - initial array
%     X - vector of terms to be removed
%     TOL - tolerance for removing terms (default 1e-8)
%
%   Outputs:
%     R0  - reduced array
%     ERR - error between given and removed terms
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if ~exist('tol','var'), tol = 1e-8; end;
%------------------------------------------------------
%       Main loop
%------------------------------------------------------
        err = [];
        for i=1:length(x)
           if length(r) < 1, 
             error('No more elements to remove');  
           end;
           [xx,ind] = sort(abs(r-x(i)));
           k = ind(1);
           if abs(xx(1)) > tol, 
              error('Unable to find a term to remove: err=%g',abs(xx(1))); 
           end;
           if imag(x(i)) == 0  &&  imag(r(k)) ~= 0
             [xx,cc] = sort(abs(r-conj(r(k)))); 
             r(cc) = real(r(cc));
           end;
           r(k) = [];
           err = [err; xx(1)];
        end;
        
%------- End of REMOVE.M --------- KYuP ---------------
