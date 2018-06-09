function [zs,z,z0] = extrpair ( z, type, tol )
%EXTRPAIR Extract a pair of symmetric zeros.
%
%     [ZS,R,Z0] = EXTRPAIR ( Z, TYPE, TOL )
%
%   Inputs:
%     Z    - initial array
%     TYPE - factorization type of 's' (default), 'z' and 'd'
%     TOL  - relative tolerance (default 1e-4)
%
%   Outputs:
%     ZS - factorized stable zeros
%     R  - unpaired nonzero zeros
%     Z0 - unpaired zeros at origin
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if ~exist('type','var'), type = 's'; end;
        if isempty(find(type == 'szd'))
          error('Unknown factorization type ''%s''', type);  
        end;
        if ~exist('tol','var'), tol = 1e-4; end;
%------------------------------------------------------
%       Start from the terms with maximal absolute values
%------------------------------------------------------
        zs = []; z = z(:);
        [xx,ind] = sort(abs(z));
        z = z(flipud(ind));
%------------------------------------------------------
%       Iteration
%------------------------------------------------------
        while 1 
          len = length(z);
          if len < 2, break; end;
          minDist = Inf; 
          iMin = 0; jMin = 0;
          for i=1:len
            for j=1:len
              if i ~= j  
                switch type   
                  case 's',        dij = abs(z(i) + z(j));  
                  case {'z', 'd'}, dij = abs(z(i)*z(j)-1);  
                end;
                if dij < minDist && ...
                  (abs(z(i)-1) < tol  ||  ...     
                   ~(abs(z(i)) == 1  &&  abs(imag(z(i))) > eps)) 
                  iMin = i; jMin = j; minDist = dij; 
                end;
              end;
            end;
          end;
          if minDist > tol  &&  minDist > tol*abs(z(iMin)), break; end;
        %------------------------------------------------------
        %       Extract pair
        %------------------------------------------------------
          if isequal(type,'s')
             if z(iMin) < 0, zs1 = z(iMin); else zs1 = z(jMin); end; 
          elseif isequal(type,'z')
             if abs(z(iMin)) < 1, zs1 = z(iMin); else zs1 = z(jMin); end; 
          else
             if abs(z(iMin)) > 1, zs1 = z(iMin); else zs1 = z(jMin); end; 
          end;
          %z([iMin jMin])
          zs = [zs; zs1];
          z([iMin jMin]) = [];
        end;
%------------------------------------------------------
%       Find unpaired complex zeros
%------------------------------------------------------
        indc = find(imag(zs) ~= 0).';
        if ~isempty(indc)
          for i = indc
            [zz,ii] = sort(abs(zs - conj(zs(i)))); 
            if ii(1) == i, 
              if imag(zs(i)) > 1e-6,    
                warning('Unpaired complex conjugate zeros');
              end;
              zs(i) = real(zs(i));  
            end;
          end;
        end;
%------------------------------------------------------
%       Find zero zeros
%------------------------------------------------------
        z0 = 0;
        ind0 = find(abs(z) < tol);
        if ~isempty(ind0)
           z0 = length(ind0); 
           z(ind0) = []; 
        end;
        
%------- End of EXTRPAIR.M --------- KYuP ----------           

