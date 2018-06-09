function [D,msg] = dtfmfoh ( F, T, t, H, type )
%DTFMFOH Discrete Laplace transform with a first-order hold.
%
%     [D,MSG] = DTFMFOH ( F, T, t )
%     [D,MSG] = DTFMFOH ( F, T, t, H )
%     [D,MSG] = DTFMFOH ( F, T, t, TYPE )
%     [D,MSG] = DTFMFOH ( F, T, t, H, TYPE )
%
%   Inputs:
%     F(s) - transfer matrix of continuous-time plant
%     T    - sampling period
%     t    - time parameter
%     H    - hold matrix suzh that size(H,3)=2 and beta=norm(H(:,:,2))   
%            beta=0 for ZOH and beta=1 for FOH
%     TYPE - class of return value: 'zpk' (default), 'ss', or 'tf'
%
%   Outputs::
%      D(z) - discrete LTI model
%      MSG  - error or warning message if any
%
%   See also DTFM, ZTRM, DTFM2.

%------------------------------------------------------
% Copyright 1999-2008 K. Polyakov 
% $Revision: 3.10 $    $Date: 09-Jan-2008 $
%------------------------------------------------------
%       Handle possible preview or delay
%------------------------------------------------------
        if ~exist('t','var') | isempty(t), t = 0; end;
        delay = totaldelay ( F );
        comdelay = min(min(delay));
        delay = delay - comdelay;        
        msg = '';
%------------------------------------------------------
%       Transform to state space
%------------------------------------------------------
        if ~isproper(F), error('Transfer matrix must be proper!'); end;
	    F = ss(F);
        A = F.a; B = F.b; C = F.c; D = F.d;
%------------------------------------------------------
%       Extrapolator
%------------------------------------------------------
        if ~exist('H','var'), H = [];
        elseif ischar(H)
           type = H; H = [];  
        end;
        if isempty(H),
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
%       Check data compatibility
%------------------------------------------------------
        if size(F,2) ~= size(H,1), error('Dimensions of plant and hold are incompatible'); end;
%------------------------------------------------------
%       Check type
%------------------------------------------------------
        if ~exist('type','var')  ||  isempty(type) 
             type = 'zpk';         
        else type = lower(type); 
        end;
        if ~isequal(type,'zpk')  &&  ~isequal(type,'ss')  &&  ~isequal(type,'tf')
           error('Incorrect type ''%s''',type); 
        end;        
%------------------------------------------------------
%       First order hold discretization
%       including fractional order hold
%	NB: so far only for delay-free case!
%------------------------------------------------------
 	    if norm(comdelay+delay) > sqrt(eps)
	      error('Unable to handle first-order hold with delay');
        end;
%------------------------------------------------------
%       Build matrices Ad and Bd of discrete state space 
%       model, they are time-independent 
%------------------------------------------------------
        nx = size(A,1);
        nu = size(F,2);
	    ny = size(F,1);
	    beta = norm(H.d(:,:,2));
   	    Ah = [  A      B      zeros(nx,nu)
              zeros(nu,nx+nu) -beta/T*eye(nu)
              zeros(nu,nx+nu)   zeros(nu)     ];
    	E = expm(Ah*T);
    	Phi = E(1:nx,1:nx);
    	Gamma10  = E(1:nx,nx+1:nx+nu);
    	Gamma1 = E(1:nx,nx+nu+1);
    	Gamma0 = Gamma10 - Gamma1;
    	Ad = [Phi Gamma1; zeros(1,nx+nu)];
    	Bd = [Gamma0; eye(nu)];
%------------------------------------------------------
%       Build Cd and Dd matrices depending on t
%------------------------------------------------------        
    	Et = expm(Ah*t);
    	Phit = Et(1:nx,1:nx);
    	Gamma10t  = Et(1:nx,nx+1:nx+nu);
    	Gamma1t = Et(1:nx,nx+nu+1);
    	Gamma0t = Gamma10t - Gamma1t;
    	Cd = [C*Phit C*Gamma1t-D*beta*t/T];
        Dd = C*Gamma0t+D*(1 + beta*t/T);
%------------------------------------------------------
%       Transform to final form
%------------------------------------------------------
        D = ss ( Ad, Bd, Cd, Dd, T );       
        switch ( type )
          case 'zpk', D = zpk(D);  
          case 'tf',  D = tf(D);  
        end;
        if isequal(class(D), 'ss')
             D = minreal(D, sqrt(eps), 0); 
        else D = minreal(D, sqrt(eps)); 
        end;
        
%------- End of DTFMFOH.M --------- KYuP ----------           
