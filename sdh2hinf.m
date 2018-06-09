function [K,err] = sdh2hinf ( sys, T, rho, o11, i11, H  )
%SDH2HINF Mixed H2/AHinf-optimization of sampled-data systems.
%
%     [K,ERR] = SDH2HINF ( SYS, T, RHO, OH2, IH2 )
%     [K,ERR] = SDH2HINF ( SYS, T, RHO, OH2, IH2, H )
%
%   Inputs:
%     SYS - an LTI system
%     T   - sampling period
%     RHO - weight factor 
%     OH2 - number of outputs in H2-block 
%     IH2 - number of inputs in H2-block 
%     H   - generalized hold as an LTI model (default ZOH)
%
%   Outputs:
%     K   - optimal discrete controller as an LTI system
%     ERR - cost function for the optimal system
%
%   See also SDH2, SDAHINF. 

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Transform to negative feedback
%------------------------------------------------------
        sys = zpk(sys);
        sys(end,:) = - sys(end,:);
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if ~exist('H', 'var'), H = []; end;
        [nout,nin] = size(sys);
        if ~exist('i11','var'), i11 = 1; end;
        if ~exist('o11','var'), o11 = 1; end;
        i12 = nin - 1 - i11;
        o12 = nout - 1 - o11;
        if i11 < 1  ||  o11 < 1, error('Incorrect size of H2-blocks'); end;
        if i12 < 0  ||  o12 < 0, error('Incorrect size of AHinf-blocks'); end;
        indI1 = 1:i11; 
        indO1 = 1:o11;        
        if i12 == 0, indI2 = indI1; else indI2 = i11+1:i11+i12; end;
        if o12 == 0, indO2 = indO1; else indO2 = o11+1:o11+o12; end;
%------------------------------------------------------
%       Construct H2 and Hinf-system
%------------------------------------------------------
        sysH2   = sys([indO1 nout],[indI1 nin]);
        sysHinf = sys([indO2 nout],[indI2 nin]);
%------------------------------------------------------
%       Construct D22
%------------------------------------------------------
        P22 = sys(end,end);
        D22 = dtfm(P22, T, 0, H)';
        D22 = setpoles(D22, exp(-allpoles(P22)*T) );
%------------------------------------------------------
%       Construct coefficients
%------------------------------------------------------
        [A2,B2,E2] = sdh2coef ( sysH2, T, [], H ); 
        [Ainf,Binf,Einf] = sdh2coef ( sysHinf, T, [], H ); 
%------------------------------------------------------
%       Solve Hinf-problem and find Sigma
%------------------------------------------------------
        [Kinf,err,Sigma] = polhinf ( Ainf, Binf, Einf, D22 );
%------------------------------------------------------
%       Solve generalized quadratic problem
%------------------------------------------------------
        A = sumzpk(rho*A2, (1-rho)*Ainf*Sigma);
        B = sumzpk(rho*B2, (1-rho)*Binf*Sigma);
        E = sumzpk(rho*E2, (1-rho)*Einf*Sigma);
        [K,err] = polquad ( A, B, E, D22 );
        err = sqrt(err);

%------- End of SDH2HINF.M --------- KYuP ----------           
        
