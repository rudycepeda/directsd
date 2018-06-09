function Kz = c2dpz ( Kc, T, type )
%C2DPZ discretization using pole-zero mapping
%
%     K = C2ZPZ ( KC, T )
%     K = C2ZPZ ( KC, T, TYPE )
%
%   Inputs:
%     KC - analog system
%     T  - sampling period
%     TYPE - discretization type of
%             'p'  - proper controller (default) 
%             'sp' - strictly proper controller
%   Outputs:
%     K - discretized controller
%

%------------------------------------------------------
%       K. Yu. Polyakov         28 Nov 2006
%                               2006
%------------------------------------------------------
%       Check syntax
%------------------------------------------------------
        if ~exist('type','var'), type = 'p'; end;
        if ~isequal(type,'p')  &  ~isequal(type,'sp')
          error('Incorrect type: ''%s''', type);  
        end;
%------------------------------------------------------
%       Find pole excess
%------------------------------------------------------
        [z,p,k] = zpkdata ( Kc, 'v' );
        q = length(p) - length(z);
        if isequal(type,'sp')
          q = q - 1;
          if q < 0, 
             error('Unable to discretize biproper controller');
          end;
        end
%------------------------------------------------------
%       Discretization
%------------------------------------------------------
        pd = exp(p*T);
        zd = exp(z*T);
        kd = prod(-z)*prod(1-pd)/prod(-p)/prod(1-zd)/2^q/k;
        Kz = zpk([zd; -ones(q,1)], pd, kd, T);

%------- End of C2DPZ.M --------- KYuP ----------           
