function [KR,err,errPrev] = sd2dof ( sys, K, H, method )
%SD2DOF Optimal feedforward controller for 2-DOF sampled-data systems.
%
%     [KR,ERR,PERR] = SD2DOF ( SYS, K, H, METHOD )
%
%   Inputs:
%     SYS - an LTI system
%     K   - known loop controller 
%     H   - generalized hold as an LTI model 
%     METHOD - calculation method of:
%        'pol' polynomial (default)
%        'ssf' state-space based on frequency domain
%        'wh'  classical Wiener-Hopf
%
%   Outputs:
%     KR   - optimal discrete reference controller
%     ERR  - cost function for the optimal system
%     PERR - limiting cost function for the optimal systems
%            as preview interval tends to Inf
%
%   See also SDL2, SD2DOFERR, SPLIT2DOF.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if ~exist('H','var'), H = []; 
        elseif ischar(H)
          method = H; H = [];  
        end;
        if ~exist('method','var'), method = 'pol'; end;         
        if ~isequal(method,'pol')  &&  ~isequal(method,'wh') ...
           &&  ~isequal(method,'ssf')     
          error('Unknown solution method ''%s''',method);  
        end;
%------------------------------------------------------
%       Transform to negative feedback
%------------------------------------------------------
        sys = zpk(sys);
        sys(end-1,:) = -sys(end-1,:);
        K = - K;
%------------------------------------------------------
%       Construct coefficients
%------------------------------------------------------
        T = K.Ts;
        [A,B,E,kCoef,n,d,A0,A1] = sd2dofcoef ( sys, T, H );
%------------------------------------------------------
%       Solve general quadratic problem
%------------------------------------------------------
        D22 = zpk([], [], 0, T); 
        if isequal(method,'pol')
           switch nargout 
            case {0, 1}, 
              PsiX = polquad ( A, B, E, D22 );
            case 2,      
              [PsiX,err] = polquad ( A, B, E, D22 );
            case 3,
              [PsiX,err,err0,err2] = polquad ( A, B, E, D22 );
              errPrev = err0 + err2;   
           end;
        elseif isequal(method,'wh')
           if nargout < 3, 
             [PsiX,err] = whquad  ( A, B, E, D22 ); 
           else
             [PsiX,err,err0,err2] = whquad  ( A, B, E, D22 ); 
             errPrev = err0 + err2;   
           end;
        elseif isequal(method,'ssf')
           errPrev = 0; 
           [PsiX,err] = ssquad ( A0, A1, B, E, D22 );
        end;
        PsiX = PsiX';
%------------------------------------------------------
%       Construct b and Delta
%------------------------------------------------------
        [a,b] = tf2nd ( K' );        
        Delta = a*n + b*d;
%------------------------------------------------------
%       Form reference controller
%------------------------------------------------------
        KR = sumzpk(minreal(PsiX*zpk(Delta,b,T)), -kCoef*K);
        KR = KR';

%------- End of SD2DOF.M --------- KYuP ----------           
        
