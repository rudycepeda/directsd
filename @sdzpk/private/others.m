function [a,x] = others ( a, b, tol )
%OTHERS Find all members of A not included in B.
%
%     [C,X] = OTHERS ( A, B )
%     [C,X] = OTHERS ( A, B, TOL )
%
%   Inputs:
%     A, B - arrays
%
%   Outputs:
%     C - members of A that are not members of B
%     X - common members of A and B
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Preliminaries
%------------------------------------------------------
        if ~exist('tol','var'), tol = sqrt(eps); end;
%------------------------------------------------------
%       Remove common elements
%------------------------------------------------------
        x = [];       
        while 1
          ii = -1;      
          minDist = tol;          
         %--------------------------------------------- 
         %   Find a pair with minimal distance < tol
         %--------------------------------------------- 
          for i=1:length(a)
            for j=1:length(b)
              if abs(a(i)-b(j)) < minDist
                 ii = i;
                 jj = j;
                 minDist = abs(a(i)-b(j));
              end;
            end;
          end;
          if ii < 0, break; end;
         %--------------------------------------------- 
         %   Remove an element
         %--------------------------------------------- 
          x = [x; a(ii)]; 
          a(ii) = [];
          b(jj) = [];          
        end

%------- End of OTHERS.M --------- KYuP ----------

