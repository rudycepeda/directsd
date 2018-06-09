function M = blckdiag ( varargin )
%BLCKDIAG Build block diagonal matrix
%
%         M = BLCKDIAG ( A1, A2, ... )
%
%       Any Ai can be numeric matrix or a cell 
%       containing an integer (in this case the identity
%       matrix of the corresponding size is inserted)
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
        n = 0;
        m = 0;
        M = [];

        for i=1:length(varargin)
          A = varargin(i);
          A = A{1};
          if iscell(A)
            A = eye(A{1});
          end;
          [n1,m1] = size(A);
          M = [M zeros(n,m1); zeros(n1,m) A];
          n = n + n1;
          m = m + m1;
        end;
        
%------- End of BLCKDIAG.M --------- KYuP ----------           
