function [Q,eAT,Xi] = sampled ( sys, muo, mui, T )
%SAMPLED Sampling for STPBC.
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Determine system sizes
%------------------------------------------------------
    n = size(sys.a,1);
    m = size(sys.b,2);    
    p = size(sys.c,1);    
    isSym = isa(sys.a,'sym');
%------------------------------------------------------
%       Xi matrix
%------------------------------------------------------
    [Xi,eAT,invXi] = xi ( sys );
%------------------------------------------------------
%       (0,T) and (T,0) cases
%------------------------------------------------------
    if muo == 1  &&  mui == 0 
      Q = sys.c*eAT*invXi*sys.om*sys.b;      
      %if isSym, Q = maple('map','radsimp',Q); end;
      return;  
    end     
    if muo == 0  &&  mui == 1 
      Q = - sys.c*invXi*sys.ups*sys.b;      
      %if isSym, Q = maple('map','radsimp',Q); end;
      return;  
    end     
%------------------------------------------------------
%     Check condition CB=0 for (0,0) and (T,T) cases
%------------------------------------------------------
    CB = sys.c*sys.b;
    if isSym
      OK = isequal(sym(CB),sym(zeros(p,m)));
    else
      OK = norm(CB) < 1e-6;
    end;
    if ~OK
      error('CB <> 0, kernel is not continuous.'); 
    end;
%------------------------------------------------------
%       (0,0) and (T,T) cases
%------------------------------------------------------
    if muo == 0  &&  mui == 0 
      Q = sys.c*invXi*sys.om*sys.b;      
      %if isSym, Q = maple('map','radsimp',Q); end;
      return;  
    end     
    if muo == 1  &&  mui == 1 
      Q = - sys.c*eAT*invXi*sys.ups*sys.b;      
      %if isSym, Q = maple('map','radsimp',Q); end;
      return;  
    end     
    error('Incorrect MUO and/or MUI (both should be either 0 or 1)');
    
%------- End of SAMPLED.M --------- KYuP ----------

