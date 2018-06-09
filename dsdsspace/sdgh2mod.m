function [dsys,gamma,dsysDelta,DP11,DP12,DP21,DP22,DP0] = sdgh2mod ( csys, T, o2, i2, t )
%SDGH2MOD Generalized H2-equivalent model for sampled-data system.
%
%     [DSYS,GAMMA,DSYSDELTA,DP11,DP12,DP21,DP22,DP0] = SDGH2MOD ( SYS, T, NMEAS, NCON, THETA )
%
%   Inputs:
%     SYS - continuous-time standard system 
%               | z | = | P11 P12 | | w |
%               | y |   | P21 P22 | | u |
%     T   - sampling period
%     NMEAS - dimension of the vector 'y' (default 1)
%     NCON  - dimension of the vector 'u' (default 1) 
%     THETA - time instant for point-wise stochastic optimization
%
%   Outputs:
%	  DSYS  - discrete-time H2-equivalent of the system
%               | zd | = | DP11 DP12 | | dw |
%               | yd |   | DP21 DP22 | | du |
%     GAMMA - additional term for variance calculation	
%     DP0   - discretized free term P11 (independent of controller)
%
%   See also SDFAST, SDH2SIMPLE.

% References:
%   [1] T. Hagiwara, M. Araki, FR-operator approach to the H2-analysis
%       and synthesis of sampled-data systems, IEEE Trans. Automat. Control,
%       vol. AC-40, no. 8, 1995.
%   [2] T. Chen, B. Francis, Optimal Sampled-Data Control Systems, 
%       Springer-Verlag, Berlin etc., 1995.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%	Continuous-time state space
%------------------------------------------------------
	    sys = minreal(ss(csys), sqrt(eps), 0);
	    if ~isct(sys)
	       error('Initial system must live in continuous time');
	    end;
        [a,b,c,d] = ssdata(sys);
        n = size(a, 1);
        [nout,nin] = size(sys);        
%------------------------------------------------------
%	Check validity
%------------------------------------------------------
        if ~exist('i2','var'), i2 = 1; end;
        if ~exist('o2','var'), o2 = 1; end;
        if ~exist('t','var'),  t = []; end;
	    i1 = nin - i2;
	    o1 = nout - o2;
	    if i1 < 1, error('Negative or zero number of reference inputs'); end;
	    if o1 < 1, error('Negative or zero number of working outputs');	end;
	    if i1 >= nin, error('No control signals specified'); end;
	    if o1 >= nout, error('No measured signals specified'); end;
%------------------------------------------------------
%	Extract partial matrices
%	ans check assumptions
%------------------------------------------------------
	    [b1,b2,c1,c2,d11,d12,d21,d22] = blocks4 ( b, c, d, o2, i2 );
	    if norm(d11) > eps, warning('D11 must be (assumed) zero');	end;
	    %d11 = 0*d11;
	    %if norm(d21) > eps, warning('D21 must be (assumed) zero');	end;
	    %d21 = 0*d21;
	    if norm(d22) > eps, warning('D22 must be (assumed) zero');	end;
	    %d22 = 0*d22;
%------------------------------------------------------
%	Input discretization
%------------------------------------------------------
        BB1d = intaxa(a,b1*b1',T);
        B1d = rectchol(BB1d,size(BB1d,1))'; % sqrtm(BB1d);        
%------------------------------------------------------
%	Loop discretization
%------------------------------------------------------
        Abar = [   a   b2
                zeros(i2,n+i2)];    
        E = expm(Abar*T);
        Ad  = E(1:n,1:n);
        B2d = E(1:n,n+1:end);        
%------------------------------------------------------
%	Output discretization
%------------------------------------------------------
        c1d12 = [c1 d12];
        Qc = c1d12'*c1d12;
        if isempty(t)
          Qd = (1/T)*intaxa(Abar', Qc, T);
        else % for special t in [0,T]
          Psi = expm(Abar*t);      
          Qd = Psi'*Qc*Psi;  
        end;
        C1dD12d = rectchol(Qd);
        C1d  = C1dD12d(:,1:n);
        D12d = C1dD12d(:,n+1:end);
%------------------------------------------------------
%	Input and outpus size of equivalent system
%------------------------------------------------------
        i1d = size(B1d,2);
        o1d = size(C1d,1);
%------------------------------------------------------
%	Additional \gamma-term in H2-norm
%------------------------------------------------------
        if isempty(t) 
          gamma = (1/T)*trace(c1*int2axa(a,b1*b1',T)*c1');
          %gamma = (1/T)*trace(b1'*int2axa(a',c1'*c1,T)*b1);
%          gamma1 = (1/T) * trace( b1' * matint(a',c1'*c1,a,T,'double') * b1 );
%   	      if isempty ( find( abs(real(eig(a))) < 1.e-6 ) )
%             Vinf = lyap ( a', c1'*c1 );
%             gamma2 = sum(diag( Vinf*(b1*b1' - BB1d/T) ));
%          else gamma2 = gamma; end;    
%          tol = 1.e-6*norm(gamma);
%          if (norm(gamma-gamma1) > tol)  ||  (norm(gamma-gamma2) > tol)
%              warning('Two procedures for calculating gamma gave different values');
%          end;
        else
          gamma = trace(c1*intaxa(a,b1*b1',T)*c1');
        end;
%------------------------------------------------------
%	D11d is chosen such that trace D11d*D11d' = gamma
%------------------------------------------------------
        D11d = zeros(o1d,i1d);
        D11d(1,1) = sqrt(gamma);
%------------------------------------------------------
%	Final discrete model with the same H2-optimal controller
%------------------------------------------------------   
        Cd = [C1d; c2];
        if norm(d21) == 0, 
          Bd = [B1d B2d];
          Dd = [ D11d  D12d  
                 zeros(o2,i1d+i2) ];
       %-----------------------------------------------   
       %  Discrete white noise case
       %-----------------------------------------------   
        else % discrete white noise case !!!
          Bd = [B1d zeros(n,i1) B2d];
          Cd = [C1d; c2];
          Dd = [ zeros(o1d,i1d)  zeros(o1d,size(d21,2))      D12d;  
                  zeros(o2,i1d)       d21                 zeros(o2,i2) ];
        end;
        dsys = minreal ( ss ( Ad, Bd, Cd, Dd, T ), sqrt(eps), 0 );       
%------------------------------------------------------
%	System with Delta-block (see [1])
%------------------------------------------------------   
        if nargout > 2
          %----------------------------------------------
          %  Construct Delta-block
          %----------------------------------------------
           ah1 = [   a      b1*b1'
 	          zeros(n)   -a'];
   	       ah2 = [   a       b2
 	           zeros(i2,n+i2) ];
 	       bh  = [    zeros(n,n+i2)
	              -c1'*c1   -c1'*d12 ];
   	       H = [   ah1           bh
	            zeros(n+i2,2*n)  ah2 ];
	       expH = expm(H * T);
           ind = n+1:2*n;
           PhiInvT = expH(ind,ind);
	       q = expH (1:n, 2*n+1:3*n+i2) + BB1d*PhiInvT*C1d'*C1dD12d;
%        Delta = frdelta ( B1d, C1dD12d, q ); 
           Delta = (pinv(B1d)*q*pinv(C1dD12d')')'; 
           Delta = Delta / sqrt(T);
          %----------------------------------------------
          %  Take into account discrete white noise
          %----------------------------------------------
           if norm(d21) == 0, 
             DdDelta = [ Delta     D12d
                      zeros(o2,i1d+i2) ];
           else
             DdDelta = [     Delta     zeros(o1d,i1)      D12d  
                        zeros(o2,i1d)       d21     zeros(o2,i2) ];
           end;
          %----------------------------------------------
           dsysDelta = minreal ( ss ( Ad, Bd, Cd, DdDelta, T ), sqrt(eps), 0 );       
        end;
%------------------------------------------------------
%	Additional data if necessary
%------------------------------------------------------   
        if nargout > 3
           DP11 = minreal( ss(Ad, B1d, C1d, Delta, T), sqrt(eps), 0 );
        end;    
        if nargout > 4
           DP12 = minreal( ss(Ad, B2d, C1d, D12d, T), sqrt(eps), 0 );
        end;
        if nargout > 5
           DP21 = minreal( ss(Ad, B1d, c2, zeros(o2,i1d), T), sqrt(eps), 0 );
        end;
        if nargout > 6
           DP22 = minreal( ss(Ad, B2d, c2, zeros(o2,i2), T), sqrt(eps), 0 );            
        end;      
        if nargout > 7
           DP0  = minreal( ss(Ad, b1, C1d, zeros(o1d,i1), T), sqrt(eps), 0 );
        end;   
 
%------- End of SDGH2MOD.M --------- KYuP ----------

%###########################################################
function d = frdelta ( w, v, q )
%FRDELTA Delta-matrix in FR-operator method.
%
%     D = FRDELTA ( W, V, Q )
%
%   Inputs:
%     W, V, Q - known constant matrices matrix
%
%   Outputs:
%     D - a matrix such that   W * D' * V = Q
%

%------------------------------------------------------
%       K. Yu. Polyakov         11 Oct 2001
%                               06 Feb 2005
%------------------------------------------------------
%	Check input data
%------------------------------------------------------
	[n,m] = size(q);
	if n ~= size(w,1)  ||  m ~= size(v,2)
	   error('Incompatible matrix dimensions');
	end;
%------------------------------------------------------
%	Find the SVD-form for matrices W and V
%------------------------------------------------------
	[lw,sw,rw] = svd ( w ); rw = rw';
	[lv,sv,rv] = svd ( v );	rv = rv';
%------------------------------------------------------
%	Form matrices Q1, Sw1, Sv1
%------------------------------------------------------
	rankW = rank ( sw );	indW = 1:rankW;
	rankV = rank ( sv );	indV = 1:rankV;
    q1 = lw' * q * rv';
    sw1 = sw ( indW, indW );
    sv1 = sv ( indV, indV );
%------------------------------------------------------
%	Check solvability
%------------------------------------------------------
	q1v = q1 ( indW, indV );
	q1 ( indW, indV ) = 0;
	if norm(q1) > 1.e-6
	   warning('Solution is badly conditioned');
	end;
%------------------------------------------------------
%	Final result
%------------------------------------------------------
	d1 = zeros ( size(w,2), size(v,1) ); 
	d1 ( indW, indV ) = sw1 \ q1v / sv1;
	d = lv * d1' * rw;

%------- End of FRDELTA.M --------- KYuP ----------
