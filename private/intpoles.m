function pInt = intpoles ( sys )
%INTPOLES Find internal poles of standard system.
%         These poles are controllable from w and u, and are 
%         observable from eps and y.
%
%     P = INTPOLES ( SYS )
%
%   Inputs:
%     SYS - a standard LTI plant
%
%   Outputs:
%      P - internal poles
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Construct four systems: w->eps, w->y, u->eps, u->y
%------------------------------------------------------
        sys = ss ( sys );
        [nout,nin] = size ( sys );
        o2 = 1; i2 = 1;
        i1 = nin - i2;
        o1 = nout - o2;
        tol = sqrt(eps);
        sys1 = sys(1:o1,1:i1);
        sys2 = sys(1:o1, end);
        sys3 = sys(end, 1:i1);
        sys4 = sys(end, end);
        
        sys1 = minreal (sys1.a,sys1.b,sys1.c,sys1.d,tol);
        sys2 = minreal (sys2.a,sys2.b,sys2.c,sys2.d,tol);
        sys3 = minreal (sys3.a,sys3.b,sys3.c,sys3.d,tol);
        sys4 = minreal (sys4.a,sys4.b,sys4.c,sys4.d,tol);
%------------------------------------------------------
%       Find common modes of all these systems
%------------------------------------------------------
        pInt = allpoles(sdzpk(sys1));
        [~,pInt] = others(pInt, allpoles(sdzpk(sys2)));
        [~,pInt] = others(pInt, allpoles(sdzpk(sys3)));
        [~,pInt] = others(pInt, allpoles(sdzpk(sys4)));

%------- End of INTPOLES.M --------- KYuP ----------           
        
