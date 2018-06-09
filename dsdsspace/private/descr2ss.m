function [A,B,C,D,flag] = descr2ss ( AD, BD, CD, DD, ED )
%DESCR2SS Transform descriptor state space to state space.
%
%     [A1,B1,C1,D1] = DESCR2SS ( A, B, C, D, E, MODE )
%     [A1,B1,C1,D1,FLAG] = DESCR2SS ( A, B, C, D, E, MODE )
%
%   Inputs: 
%     A, B, C, D, E - matrices of compatible dimensions
%
%   Outputs:
%     A1,B1,C1,D1 - state space matrices
%     FLAG - flag indicating (by the value of 0) unsuccessful transformation
%            (transfer function is improper);
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check input data
%------------------------------------------------------
        [n,n1] = size(AD);
        if n ~= n1, error('Matrix A must be square'); end;
        if any(size(ED) ~= n), error('Matrices A and E must have the same dimensions'); end;        
        nin  = size(BD,2);
        nout = size(CD,1);
%------------------------------------------------------
%       Separate regular and irregular parts of ED
%------------------------------------------------------
        [A,E22,Q,Z] = weierstr ( AD, ED );     
        n1 = max(size(A));
        n2 = max(size(E22));        
        CZ = CD*Z;
        QB = Q*BD;
%------------------------------------------------------
%       Form state space matrices
%------------------------------------------------------
        flag = 1;
        if n1 == 0
           B = zeros(0, nin);
           C = zeros(nout, 0);
        else
           B = QB(1:n1,:);
           C = CZ(:,1:n1);
        end;        
        if n2 == 0, 
           D = DD;     
        else           
           k = 0;
           temp = E22;
           while norm(temp) > eps, % -- determine degree of improper part
             temp = temp * E22;
             k = k + 1;
           end;
           if k > 0  
              flag = 0;
              if nargin < 5,
                 error('Transfer matrix is not proper, transition to state space impossible');  
              end;
           end           
           B2 = QB((n1+1):end,:);
           C2 = CZ(:,(n1+1):end);
           D = C2*eye(n2)*B2 + DD;   % --- for proper systems only
        end

%------- End of DESCR2SS.M --------- KYuP ----------           

