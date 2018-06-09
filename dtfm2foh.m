function DF = dtfm2foh ( F, T, H )
%DTFM2FOH Special symmetric discrete Laplace transform with 
%         a first-order holdhold.
%
%     D = DTFM2FOH ( F, T )
%     D = DTFM2FOH ( F, T, H )
%
%   Inputs:
%     F(s) - transfer matrix of continuous-time plant
%     T    - sampling period
%     H    - hold matrix suzh that size(H,3)=2 and beta=norm(H(:,:,2))   
%            beta=0 for ZOH and beta=1 for FOH
%   Outputs:
%     D(z) - discrete transformation D{H~*F~*F*H}(0) in zpk-form
%
%   See also ZTRM, DTFM, DTFMFOH, DTFM2.

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
        if ~exist('H','var'),
             H = ss(zeros(1,1,2));
             H(:,:,1) = ss(tf(1)); 
             H(:,:,2) = H(:,:,1); % scalar FOH by default
        else H = ss(H); 
        end;     
        if size(H,1)+size(H,2) == 2,  % scalar extrapolator  
           scalarH = 1; H0 = H; 
           H = H * eye(size(F,2)); 
        else 
           error('Matrix first order hold is not admissible so far');
           scalarH = 0; 
           H0 = H;
        end; % scalar extrapolator              
%------------------------------------------------------
%       Size compatibility
%------------------------------------------------------       
        rH = size(H,1);
        if cF ~= rH, error('Dimensions of plant and hold are incompatible'); end;
        nx = size(A,1);       
        nu = size(F,2);
        ny = size(F,1);
%------------------------------------------------------
%       Construct discrete plant model with FOH
%       and extract 1/z
%------------------------------------------------------
    	beta = norm(H.d(:,:,2));
    	Av = [  A      B        zeros(nx,nu)
               zeros(nu,nx+nu) -beta/T*eye(nu)
               zeros(nu,nx+nu)    zeros(nu)        ];
    	E = expm(Av*T);
    	Phi = E(1:nx,1:nx);
    	Gamma10 = E(1:nx,nx+1:nx+nu);
    	Gamma1 = E(1:nx,nx+nu+1:end);
    	Gamma0 = Gamma10 - Gamma1;
        Phi = [Phi Gamma1; zeros(nu,nx+nu)];
    	Gamma = [Gamma0; eye(nu)];
        %Theta = inv(Phi');        
%------------------------------------------------------
%       Computing matrix integrals
%------------------------------------------------------
        CD = [C D zeros(ny,nu)];
        [G,PhiG] = intaba ( -Av', CD'*CD, Av, T );
        G = PhiG \ G;
%------------------------------------------------------
%       Transform to xv[k] = [ x[k]; nu[k-1] ]
%------------------------------------------------------
        Inu = eye(nu);
        R = [eye(nx) zeros(nx,2*nu)
             zeros(nu,nx+nu) Inu
             zeros(nu,nx) Inu -Inu ];
        G = R' * G * R;
        n = nx + nu;
        g11 = G(1:n,1:n);
        g12 = G(1:n,n+1:end);
        g21 = G(n+1:end,1:n);
        g22 = G(n+1:end,n+1:end);        
%------------------------------------------------------
%       Equivalent DTF
%------------------------------------------------------
        Ad = Phi;
        Bd = Gamma;
        CD = chol(G);
        Cd = CD(:,1:nx+nu);
        Dd = CD(:,nx+nu+1:end);
        DF = ss(Ad, Bd, Cd, Dd, T);       
        zpkDF = zpk(DF);
%------------------------------------------------------
%       Final result
%------------------------------------------------------
        DF = minreal(zpkDF'*zpkDF, sqrt(eps));
        
%------- End of DTFM2FOH.M --------- KYuP ----------
