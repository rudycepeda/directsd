function display ( sys )

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
    iname = inputname(1);
    disp([iname,': ','[0,T]-system with two-point boundary conditions']);
    A = sys.a
    B = sys.b
    C = sys.c
    D = sys.d
    Om = sys.om
    Ups = sys.ups

%------- End of DISPLAY.M --------- KYuP ----------
   
