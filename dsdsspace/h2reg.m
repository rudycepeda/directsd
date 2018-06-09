function [K,h2norm,clSys] = h2reg ( sys0, o2, i2, options )
%H2REG H2-optimal controller for an LTI-system.
%
%     [C,N,CL] = H2REG ( SYS, NMEAS, NCON, OPTIONS )
%
%   Inputs:
%     SYS - an LTI system
%     NMEAS - dimension of the vector 'y' (default 1)
%     NCON  - dimension of the vector 'u' (default 1) 
%     OPTIONS - structure with possible fields
%        .tol - tolerance for relative error (default 1.e-4)
%        .method  - solution method of
%            'sa' - method by Safonov & Chiang - STRICTLY proper controllers 
%                   (default and the only method for continuous systems)
%            'ch' - formulas by Chen and Francis - proper discrete controllers
%                   (default method for discrete systems)   
%   Outputs:
%     K  - H2-optimal controller as an LTI system
%     N  - H2-norm for the optimal closed-loop system
%     CL - transfer function the optimal closed-loop system
%
%   See also H2REG, SDH2, CH2.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%	Check syntax
%------------------------------------------------------
    if ~exist('o2','var'), o2 = 1; end;
    if ~exist('i2','var'), i2 = 1; end;
    if exist('options','var'),
       if isfield(options,'tol'),     tol = options.tol; end;
       if isfield(options,'method'),  method = options.method; end;
    end;
    if ~exist('gammaTol', 'var'), tol = 1.e-4; end;
    if ~exist('method', 'var'),   
       if isct(sys0), method = 'sa'; 
       else           method = 'ch'; end;
    end;
    if ~isequal(method,'sa')  &&  ~isequal(method,'ch')
       error('Unknown solution method ''%s''', method); 
    end;    
    if isequal(method,'ch')  &&  isct(sys0)
       warning('Only Safonov & Chiang method is allowed for continuous-time systems');        
       method = 'sa';
    end;
%------------------------------------------------------
%	Continuous-time state space
%------------------------------------------------------    
	sys0 = minreal(ss(sys0), sqrt(eps), 0);  
    sys = sys0;
    [valid,regsys] = regular ( sys, o2, i2, tol );
    if isequal(method,'sa') ||  valid == -1    % --- method of Safonov and Chiang requires regular system
       sys = regsys;
    end;
    [a,b,c,d] = ssdata(sys);
    n = size(a,1);
    %[nout,nin] = size(sys);
    %i1 = nin - i2;
	%o1 = nout - o2;
%------------------------------------------------------
%	Extract state-space matrices
%------------------------------------------------------
    [b1,b2,c1,c2,d11,d12,d21,d22] = blocks4 ( b, c, d, o2, i2 );
%------------------------------------------------------
%	Assumption on D11 for Safonov and Chiang technique
%------------------------------------------------------
	if isequal(method,'sa')  &&  norm(d11) > 1.e-10
	   warning('Matrix D11 must be (assumed to be) zero');
	   d11 = d11 * 0;
	end;
%------------------------------------------------------
%	Check by means of Robust Toolbox
%------------------------------------------------------
	if 0  &&  exist('mksys.m','file')   
       TSS_ = mksys (a, b1, b2, c1, c2, d11, d12, d21, d22,'tss');
	   if isct(sys)
            SS_CP = h2lqg (TSS_);    %[SS_CP,SS_CL] = h2lqg (TSS_);
       else SS_CP = dh2lqg (TSS_);   %[SS_CP,SS_CL] = dh2lqg (TSS_);
       end;	
       disp(' ')
	   disp ('Optimal controller as given by Robust Toolbox:');
       disp ('-----------------------------------------------');
	   [ac,bc,cc,dc] = branch (SS_CP,'a,b,c,d');
	   K = - zpk(ss(ac,bc,cc,dc,sys.Ts));
	   % disp ('Closed-loop system:');
	   % [ac,bc,cc,dc] = branch (SS_CL,'a,b,c,d');
	   % CL = zpk(ss(ac,bc,cc,dc,sys.Ts));
	end;    
%------------------------------------------------------
%	Check by means of Mu-Tools Toolbox
%------------------------------------------------------
	if 0  &&  isct(sys)  &&  exist('h2syn.m','file')       
       plant = pck (sys.a, sys.b, sys.c, [d11 d12; d21 d22] );
       [k,g,norms,kfi,gfi,hamx,hamy] = h2syn(plant, 1, 1, 2, 0);
       disp(' ')
       disp ('Optimal controller as given by Mu-Toolbox:');
       disp ('-------------------------------------------');
	   [ac,bc,cc,dc] = unpck (k);
	   K = zpk(ss(ac,bc,cc,dc,sys.Ts));
   end;        
%------------------------------------------------------
%	Normalization so that D12'*D12=I, D12*D12'=I
%------------------------------------------------------
    if isequal(method,'sa')
	  tt1 = sqrtm(d12'*d12);
	  t1  = inv(tt1);
	  tt2 = sqrtm(d21*d21');
	  t2  = inv(tt2);
	  b2  = b2 * t1;
	  c2  = t2 * c2;
	  d12 = d12 * t1;
	  d21 = t2 * d21;
	  d22 = t2 * d22 * t1;
    end;	
%------------------------------------------------------
%	Special case: state feedback
%------------------------------------------------------
    stateFdb = (isequal(c2,eye(n)) & norm(d21) == 0  & norm(d22) == 0);
%------------------------------------------------------
%	Balanced realization
%------------------------------------------------------
    temp = ss ( a, [b1 b2], [c1; c2], [d11 d12; d21 d22], sys.Ts );
    [a,b,c,d] = ssdata ( ssbal ( temp ) );   
    [b1,b2,c1,c2,d11,d12,d21,d22] = blocks4 ( b, c, d, o2, i2 );
%------------------------------------------------------
%	Solution to LQR-problem with state feedback
%------------------------------------------------------
	Q = c1'*c1;
	R = d12'*d12;   % must be = eye(i2) for method == 1
	N = c1'*d12;  
	if isct(sys)
       Kc = lqr1 (a, b2, Q, R, N);
    else
       if isequal(method,'sa')
          Kc = dlqr1 (a, b2, Q, R, N);
       else
          %---- method by Chen and Francis (for general form of D12 and D21)
          %     original formulas by Chen and Francis, for nonsingular MM
          %  MM = d12'*d12;
          %  if rank(MM) < size(MM,1)
          %     error('Optimal control problem is singular'); 
          %  end;    
          %  ax = a - (b2 / MM) * d12'*c1;
          %  Qx = c1' * (eye(o1) - (d12/MM)*d12') * c1;
          %  Rx = MM;
          %  [dummy,X,cerr1] = dlqr1 ( ax, b2, Qx, Rx );
          %  F = - (MM + b2'*X*b2) \ (b2'*X*a + d12'*c1);
          %  F0 = - (MM + b2'*X*b2) \ (b2'*X*b1 + d12'*d11);
          X = dare1 ( a, b2, Q, R, N);
          F = - (R + b2'*X*b2) \ (b2'*X*a + d12'*c1);
          F0 = - (R + b2'*X*b2) \ (b2'*X*b1 + d12'*d11);
      end;
   end;
%------------------------------------------------------
%	Solution to optimal filtering problem
%------------------------------------------------------
	Q = b1*b1';
    R = d21*d21'; % must be = eye(o2) for method == 1
	N = b1*d21';
    if isct(sys)
        Kf = lqr1 (a', c2', Q, R, N );
        Kf = Kf';
    else
       if isequal(method,'sa')
          Kf = dlqr1 (a', c2', Q, R, N );
          Kf = Kf';
       else 
          %---- method by Chen and Francis (for general form of D12 and D21)
          % --- original formulas by Chen and Francis, for nonsingular NN  
          %  NN = d21*d21';
          %  if rank(NN) < size(NN,1)
          %     error('Optimal filtering problem is singular'); 
          %  end;    
          %  ay = a - (b1*d21' / NN) * c2;
          %  Qy = b1 * (eye(i1) - (d21'/NN)*d21) * b1';
          %  Ry = NN;
          % [dummy,Y,ferr1] = dlqr1 ( ay', c2', Qy, Ry );
          % L = - (a*Y*c2' + b1*d21') / (NN + c2*Y*c2');
          % L0 = (F*Y*c2' + F0*d21') / (NN + c2*Y*c2');          
          if stateFdb
               Y = eye(n); % any value of compatible dimensions!
          else Y = dare1 ( a', c2', Q, R, N); end;
          L = - (a*Y*c2' + b1*d21') / (R + c2*Y*c2');
          L0 = (F*Y*c2' + F0*d21') / (R + c2*Y*c2');
      end;
    end;
%------------------------------------------------------
%	Final H2-optimal controller
%------------------------------------------------------
    if isequal(method,'sa')
	   ac = a - Kf*c2 - b2*Kc + Kf*d22*Kc;
	   bc = Kf;
	   cc = - Kc;
	   dc = zeros( size(cc,1), size(bc,2) );
    else
	   ac = a + b2*F + L*c2 - b2*L0*c2;
	   bc = L - b2*L0;
	   cc = L0*c2 - F;
	   dc = L0;
    end;   
%------------------------------------------------------
%	Scale the controller back
%------------------------------------------------------
    if isequal(method,'sa')
       bc = bc * t2;
	   cc = t1 * cc;
	   dc = t1 * dc * t2;
    end;    
	K = minreal ( ss ( ac, bc, cc, dc, sys.Ts ), sqrt(eps), 0 );		
%------------------------------------------------------
%	Regularized case - not strictly proper controller
%------------------------------------------------------
	if valid <1  &&  isct(sys)
       [dr,nr] = tfdata ( 1/zpk(K), 'v' );
	   dr = striplz(dr, 1.e-6*norm(dr));
	   nr = striplz(nr, 1.e-6*norm(nr));
	   K = zpk(tf(nr, dr));
    elseif i2 == 1  &&  o2 == 1 % cancellation for SISO
       K = minreal(zpk(K));
       K = ss(K);
	end;
%------------------------------------------------------
%	Additional outputs
%------------------------------------------------------
    if nargout > 1
       if isa(K,'ss') 
            clSys = minreal(lft(sys0,K), sqrt(eps), 0);  
       else clSys = minreal(lft(sys0,K)); 
       end;
       h2norm = norm(clSys);
    end; 

%------- End of H2REG.M --------- KYuP ----------



