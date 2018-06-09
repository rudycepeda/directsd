function [A,B,E,kCoef,n,d,A0,A1] = sd2dofcoef ( sys, T, H )
%SD2DOFCOEF Coefficients for optimal 2-DOF sampled-data system design.
%
%     [A,B,E,C,N,D,A0,A1] = SD2DOFCOEF ( SYS, T, H )
%
%   Inputs:
%     SYS - an LTI system
%     T   - sampling period
%     H   - generalized hold as an LTI model 
%
%   Outputs:
%     A, B, E - function coefficients in quadratic functionsl
%     C       - direct transfer coefficient for 2-DOF design 
%     N, D    - numerator and denominator of discrete model D22
%     A0, A1  - functions such that A = A0 * A1 
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if ~exist('H','var'), H = []; end;
        [nout,nin] = size(sys);
        i2 = 1; i1 = nin - i2;
        if i1 ~= 1, error('L2-problem is formulated for a single input'); end;
        o2 = 2; o1 = nout - o2;
        if o1 < 1, error('Incorrect number of output signals'); end;
        P221 = sys(end-1,end);
        if abs(P221.k) > eps,
           error('Incorrect 2-DOF structure (P221 <> 0)'); 
        end;
%------------------------------------------------------
%       Construct separated systems
%------------------------------------------------------
        sys2 = sys; sys2(end,:) = [];               
        P212 = sys (end,1:i1);
        P211 = sys2(end,1:i1);
        kCoef = minreal(P212/P211);
        if ~isstatic(kCoef)  &&  abs(kCoef.k) > eps
          error('Ensure P212=k*P211 in this version');  
        end;
        kCoef = kCoef.k;
%------------------------------------------------------
%       Construct a/Delta
%------------------------------------------------------
        D22 = dtfm(sys(end,end), T, 0, H)';
        [n,d] = tf2nd ( D22 );
        dzpk = zpk(d, 1, T);
%------------------------------------------------------
%       Construct coefficients
%------------------------------------------------------
        [A,B,E,A0,A1] = sdl2coef ( sys2, T, H );        
        A  = minreal(A*dzpk*dzpk');
        A0 = minreal(A0*dzpk*dzpk');
        B = minreal(B*dzpk);

%------- End of SD2DOFCOEF.M --------- KYuP ----------           
