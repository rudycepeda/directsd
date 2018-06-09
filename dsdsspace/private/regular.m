function [valid,sys,report] = regular ( sys, o2, i2, tol, method )
%REGULAR Regularize standard system.
%
%     [VALID,SYS1,REPORT] = REGULAR ( SYS, NMEAS, NCON, TOL, METHOD )
%
%   Inputs:
%     SYS - standard 4-block system model 
%             | z | = | P11 P12 | | w |
%             | y |   | P21 P22 | | u |
%     NMEAS - dimension of the vector 'y' 
%     NCON  - dimension of the vector 'u' 
%     TOL - size of regularizing perturbation (default 1.e-6)
%     METHOD - method code:
%            1 - based on diagonal perturbations (default)
%            2 - randomized
%   Outputs:
%	  VALID - validity flag for initial system
%        -1 - the system cannot be optimized in the state space
%         0 - the system has been regularized
%         1 - the system is valid for Riccati-based optimization
%	  SYS1   - regularized 4-block system model
%     REPORT - report message (if specified, errors and warnings are suppressed)
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 21-Sep-2006 $
%------------------------------------------------------
%	Continuous-time state space
%------------------------------------------------------
    report = '';
	sys = minreal(ss(sys), sqrt(eps), 0);
    [a,b,c,d] = ssdata(sys);
    n = size(a,1);
    [nout,nin] = size(sys);
    i1 = nin - i2;
    o1 = nout - o2;
    if ~exist('tol', 'var'),  tol = 1.e-6;  end;       
    if ~exist('method', 'var')
       method = 1;
    elseif method ~= 1  &&  method ~= 2
       report = 'Incorrect regularization method';
       if nargin < 3, error(report); end;
       return;
    end;    
%------------------------------------------------------
%	Check validity
%------------------------------------------------------
	[valid,info] = syscheck(sys, o2, i2);
    if valid == 1, return; end;      
    if valid < -1, 
       warning(info);
       report = 'System cannot be regularized for Riccati-based design';
       if nargin < 3, error(report); end;
       return;
    end;
    report = 'The standard system is singular (regularization performed)';
    if nargin < 3, warning(report); end;
%------------------------------------------------------
%	Extract partial matrices
%------------------------------------------------------
	[b1,b2,c1,c2,d11,d12,d21,d22] = blocks4 ( b, c, d, o2, i2 );
%-----------------------------------------------------------
%	Regularization of d12 so that d12'*d12 is nonsingular
%   (d12 has full column rank)
%------------------------------------------------------
    if det(d12'*d12) == 0
       [r12,c12] = size(d12);     
       if rank(b2) < n  % 21.09.06 for full-rank control case 
         if r12 < c12     % --- singular problem - add a square output block
            %i0 = c12 - r12;
            c1  = [c1;  zeros(c12,n)];
            d11 = [d11; zeros(c12,i1)];
            if method == 1 
                 d12 = [d12; ones(c12)*tol];
            else d12 = [d12; rand(c12)*tol]; 
            end;
         else % if isct(sys)   % 12.09.2002 for Safonov-Chiang design in discrete case
            if method == 1 
                 d12 = d12 + eye(r12,c12)*tol;
            else d12 = d12 + rand(r12,c12)*tol; 
            end;
         end;
       end;                    
	end;
%-----------------------------------------------------------
%	Regularization of d21 so that d21*d21' is nonsingular
%   (d21 has full row rank)
%------------------------------------------------------
	if det(d21*d21') == 0            
       [r21,c21] = size(d21);
       stateFdb = (rank(c2) == n & norm(d21) == 0  & norm(d22) == 0);
       if ~stateFdb  % 21.09.06 for state feedback case 
         if r21 > c21     % --- singular problem - add a square input block
            b1  = [b1 zeros(n,r21)];
            d11 = [d11 zeros(o1,r21)];
            if method == 1 
                 d21 = [d21 eye(r21)*tol];
            else d21 = [d21 rand(r21)*tol]; 
            end;
         else  % if isct(sys)  % 12.09.2002 for Safonov-Chiang design in discrete case
            if method == 1 
                 d21 = d21 + eye(r21,c21)*tol;
            else d21 = d21 + rand(r21,c21)*tol; 
            end;
         end;
       end;       
       
	end;
%-----------------------------------------------------------
%	Construct regularized system
%------------------------------------------------------
    sys = ss ( a, [b1 b2], [c1; c2], [d11 d12; d21 d22], sys.Ts );
  
%------- End of REGULAR.M --------- KYuP ----------
