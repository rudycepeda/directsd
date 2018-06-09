function V = vec ( P,  type )
%VEC Vectorize coefficients of a (quasi)polynomial.
%
%     V = VEC ( P, TYPE );
%       
%   Inputs:  
%     P - a polynomial or a quasipolynomial
%     TYPE - type of vectorization of
%         'all' - default
%          '+'  - by ascending non-negative powers
%          '-'  - by ascending non-positive powers
%   Outputs: 
%     V - coefficient vector
%
%   See also TOEP, HANK.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if ~exist('type','var'), type = 'all'; end;
        if ~isequal(type,'all')  &&  ~isequal(type,'-')  &&  ~isequal(type,'+')
           error('Unknown type of vector ''%s''',type);  
        end;
        if isobject(P)  
           if isa(P,'poln'), 
                C = P.coef;
                shift = P.shift;
           else error('First argument must be a polynomial'); 
           end;
        else
           C = P(:)';
           shift = floor(deg(C)/2);
        end;        
%------------------------------------------------------
%       Form required vector
%------------------------------------------------------
        V = fliplr(C);
        switch type,
          case 'all',
          case '-',
              V = V(1:shift+1);
          case '+',
              V = V(shift+1:end);
        end;
        V = V(:);
        
%------- End of VEC.M --------- KYuP ----------           
