function DF = dtfm2 ( F, T, H, type )
%DTFM2 Special symmetric discrete Laplace transform with a hold.
%
%     D = DTFM2 ( F, T )
%     D = DTFM2 ( F, T, H )
%     D = DTFM2 ( F, T, TYPE )
%     D = DTFM2 ( F, T, H, TYPE )
%
%   Inputs:
%     F(s) - transfer matrix of continuous-time plant
%     T    - sampling period
%     TYPE - class of return value: 'zpk' (default), 'ss', or 'tf'
%     H(s) - transfer matrix of the hold:
%              h(t)=K*exp(L*t)*M => H(s)=K*(sI-L)^{-1}*M
%            default hold is the zero-order hold (ZOH) with H(s)=1/s.
%   Outputs:
%     D(z) - discrete transformation D{Gh~*F~*F*Gh}(0), where
%             Gh(s) = int[0..T] h(t) exp(-st) dt 
%                   = K*(sI-L)^{-1}*[I-exp(L*T)*exp(-sT)]*M
%
%   See also ZTRM, DTFM.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 09-Jan-2009 $
%------------------------------------------------------
%       Transform to state space and check correctness
%------------------------------------------------------
        if ~isproper(F), error('Transfer matrix must be proper!'); end;
        F = ss(F);
        cF = size(F,2);
        A = F.a; B = F.b; C = F.c; D = F.d;
%------------------------------------------------------
%       Extrapolator
%------------------------------------------------------
        if ~exist('H','var'), H = [];
        elseif ischar(H)
           type = H; H = [];  
        end;
        if isempty(H), 
             H = ss(0, 1, 1, 0); % scalar ZOH by default
        else H = ss(H); 
        end;     
        if size(H,1)+size(H,2) == 2,  % scalar extrapolator  
           H = H * eye(size(F,2)); 
        end; % scalar extrapolator              
        L = H.a; M = H.b; K = H.c;
        rH = size(H,1);
        if cF ~= rH, error('Dimensions of plant and hold are incompatible'); end;
        n = size(A,1);       
        m = size(L,1);
%------------------------------------------------------
%       Type of return function
%------------------------------------------------------       
        if ~exist('type','var')  ||  isempty(type) 
             type = 'zpk';         
        else type = lower(type); 
        end;
        if ~isequal(type,'zpk')  &&  ~isequal(type,'ss')  &&  ~isequal(type,'tf')
           error('Incorrect type %s',type); 
        end;       
%------------------------------------------------------
%       Special case: first order hold
%       including fractional order hold
%------------------------------------------------------
        if size(H,3) > 1 
          DF = dtfm2foh ( F, T, H );
          return;
    	end;
%------------------------------------------------------
%       Computing matrix integrals
%------------------------------------------------------
        [Gamma,Phi] = intaba ( A, B*K, L, T );
        Theta = inv(Phi');        
        Av = [A B*K; zeros(m,n) L];
        CDK = [C D*K];
        [G,PhiG] = intaba ( -Av', CDK'*CDK, Av, T );
        G = PhiG \ G;
        g11 = G(1:n,1:n);
        g12 = G(1:n,n+1:end);
        g21 = G(n+1:end,1:n);
        g22 = G(n+1:end,n+1:end);        
%------------------------------------------------------
%       State-space matrices
%------------------------------------------------------
%        Ad = [Theta Theta*Theta*g11; zeros(n) Phi];
%        Bd = [Theta*Theta*g12; Gamma]*M;
%        Cd = M'*[-Gamma' g21-Gamma'*Theta*g11];
        Ad = [Theta Theta*g11; zeros(n) Phi];
        Bd = [Theta*g12; Gamma]*M;
        Cd = M'*[-Gamma'*Theta g21-Gamma'*Theta*g11];
        Dd = M'*(g22 - Gamma'*Theta*g12)*M;       
%------------------------------------------------------
%       Final result
%------------------------------------------------------
        DF = ss(Ad, Bd, Cd, Dd, T);       
        zpkDF = sdzpk(DF);
        if norm(zpkDF.k) < eps,
           DF = ss(zeros(size(F)));
           DF.Ts = T;
        end;
        switch(type)
          case 'zpk', DF = sdzpk(DF);        
          case 'tf',  DF = sdtf(DF);        
        end;        
        
%------- End of DTFM2.M --------- KYuP ----------
