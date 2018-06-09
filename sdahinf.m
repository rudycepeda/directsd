function [K,err,Sigma] = sdahinf ( sys, T, H )
%SDAHINF AHinf-optimal controller for sampled-data systems.
%
%     [K,ERR] = SDAHINF ( SYS, T )
%     [K,ERR] = SDAHINF ( SYS, T, H )
%
%   Inputs:
%     SYS - an LTI system
%     T   - sampling period
%     H   - generalized hold as an LTI model 
%
%   Outputs:
%     K   - AHinf-optimal discrete controller
%     ERR - cost function for the optimal system
%
%   See also SDAHINORM, SDH2HINF, SDTRHINF.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Transform to negative feedback
%------------------------------------------------------
        if ~exist('H','var'), H = []; end;
        sys = zpk(sys);
        sys(end,:) = - sys(end,:);
%------------------------------------------------------
%       Construct D22
%------------------------------------------------------
        P22 = sys(end,end);
        D22 = dtfm(P22, T, 0, H)';
        D22 = setpoles(D22, exp(-allpoles(P22)*T) );
%------------------------------------------------------
%       Construct coefficients
%------------------------------------------------------
        [A,B,E] = sdh2coef ( sys, T, [], H );
%------------------------------------------------------
%       Solve Hinf-type problem
%------------------------------------------------------
        PCancel = intpoles ( sys );  
        switch nargout 
          case {0, 1, 2}, 
            [K,err] = polhinf ( A, B, E, D22, PCancel );
          otherwise, 
            [K,err,Sigma] = polhinf ( A, B, E, D22, PCancel );
        end;
        err = sqrt(T)*err;

%------- End of SDHAHINF.M --------- KYuP ----------           
        
