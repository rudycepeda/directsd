function [valid,report] = syscheck ( sys, o2, i2 )
%SYSCHECK Check standard assumtions for H2- and Hinf-problems.
%
%     [VALID,REPORT] = SYSCHECK ( SYS, NMEAS, NCON )
%
%   Inputs:
%     SYS - continuous-time standard system 
%               | z | = | P11 P12 | | w |
%               | y |   | P21 P22 | | u |
%     NMEAS - dimension of the vector 'y' 
%     NCON  - dimension of the vector 'u' 
%
%   Outputs:
%     VALID - validity flag 
%        -2 - the system cannot be optimized in the state space
%             (empty blocks, not stabilizable and/or detectable)
%        -1 - assumptions dim(z)>=dim(u) or dim(w)>=dim(y) do not hold,
%             it is required to increase number of inputs or outputs  
%         0 - the system can be regularized (rank assumptions on
%             D12 and D21 do not hold) 
%         1 - the system is valid for Riccati-based optimization
%     REPORT - string error report
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 21-Sep-2006 $
%------------------------------------------------------
%	Continuous-time state space
%------------------------------------------------------
    sys = minreal(ss(sys), sqrt(eps), 0);
    [a,b,c,d] = ssdata(sys);
    [nout,nin] = size(sys);        
    n = size(a, 1);
%------------------------------------------------------
%	Check validity
%------------------------------------------------------
    valid = -2;
    report = '';
	i1 = nin - i2;
	o1 = nout - o2;
	if i1 < 1, report = 'Negative or zero number of reference inputs';  return; end;
	if o1 < 1, report = 'Negative or zero number of working outputs';  return; end;
	if i1 >= nin, report = 'No control signals specified'; return; end;
	if o1 >= nout, report = 'No measured signals specified';  return;  end;    
%------------------------------------------------------
%	Check stabilizabilty and detectability  of the 2-2 block
%------------------------------------------------------
	tol = 1.e-10;    
	[b1,b2,c1,c2,d11,d12,d21,d22] = blocks4 ( b, c, d, o2, i2 );
    [syss,sysu] = separss ( ss (a, b2, c2, d22, sys.Ts) );
    if ~isempty(sysu.a)
       h = hankelsv ( - sysu.a - eye(size(sysu.a)), sysu.b, sysu.c );
       if min(h) < tol*max(h)
          report = 'The system is not stabilizable and/or detectable';
          return;
       end
    end
%------------------------------------------------------
%   Nondegeneracy assumptions    
%------------------------------------------------------
    valid = -1;
    if i2 > o1
       report = 'Assumption dim(z) >= dim(u) does not hold'; 
       return; 
    end;
    stateFdb = (rank(c2) == n & norm(d21) == 0  & norm(d22) == 0);
    if o2 > i1 && ~stateFdb  % 21.09.06 for state feedback case 
       report = 'Assumption dim(w) >= dim(y) does not hold'; 
       return; 
    end;  
%------------------------------------------------------
%   Full-rank matrices D12 and D21
%------------------------------------------------------
    valid = 0;
	if det(d12'*d12) == 0
	   report = 'Matrix D12 must have full column rank';
	   return;
	end;
	if det(d21*d21') == 0
	   report = 'Matrix D21 must have full row rank';
	   return;
	end;
%------------------------------------------------------
%   Completely OK
%------------------------------------------------------
    valid = 1;

%------- End of SYSCHECK.M --------- KYuP ----------
