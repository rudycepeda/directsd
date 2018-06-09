function sys = mtimes ( sys1, sys2 )
%MTIMES Multiplication for STPBC C = A*B.
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
    n1 = size(sys1.a,1);
    n2 = size(sys2.a,1);
    sys = stpbc();
    sys.a = [sys1.a  mul(sys1.b,sys2.c)
             zeros(n2,n1) sys2.a];
    sys.om = blckdiag ( sys1.om, sys2.om );
    sys.ups = blckdiag ( sys1.ups, sys2.ups );
    sys.b = [mul(sys1.b,sys2.d); sys2.b];
    sys.c = [sys1.c mul(sys1.d,sys2.c)];
    sys.d = mul(sys1.d,sys2.d);
   
%------- End of MTIMES.M --------- KYuP ----------

function X = mul(A, B)
%MUL Symbolic matrix multiplication.
%      
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Jan-2007 $
%------------------------------------------------------
%   Numerical case
%------------------------------------------------------
    try 
       AN = double(A); 
       if isnumeric(AN), A = AN; end;
    catch end;
    try 
       BN = double(B); 
       if isnumeric(BN), B = BN; end;
    catch end;
    if isnumeric(A)  &&  isnumeric(B)
       X = A*B;
       return;
    end;
%------------------------------------------------------
%   Change conj for transp
%------------------------------------------------------
    A = subs(sym(A), 'conj', 'transp');
    B = subs(sym(B), 'conj', 'transp');
    X = char(simplify(A*B));
%------------------------------------------------------
%   Back substitution
%------------------------------------------------------
    pos = regexp(X,'transp');
    for i=length(pos):-1:1
      X = [X(1:pos(i)-1) 'conj' X(pos(i)+6:end)];  
    end;
    X = sym(X);
    
%------- End of MUL.M --------- KYuP ----------

