function [G,GK,a0,b0] = psigain ( sys, T, Delta )
%PSIGAIN DC-gain of Psi-parameter for tracking system.
%
%     [G,GK,A0,B0] = PSIGAIN ( SYS, T, DELTA )
%
%   Inputs:
%     SYS   - an LTI system
%     T     - sampling period 
%     DELTA - basic characteristic polynomial such that
%               A0*n + B0*d = Delta
%             where D_{P22H}(T,\z,0) = n(\z) / d(\z). 
%             By default, Delta = 1.
%   Outputs:
%     G  - DC-gain of Psi parameters or NaN if arbitrary
%     GK - DC-gain of controller
%     A0,B0 - numerator and denominator of the basic controller
%
%   See also MODSDL2.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Transform to negative feedback
%------------------------------------------------------
        sys = zpk(sys);
        sys(end,:) = - sys(end,:);
%------------------------------------------------------
%       Default polynomial Delta = 1
%------------------------------------------------------
        if ~exist('Delta','var'), Delta = 1; end;
        if any(abs(roots(Delta)-1) < 1e-5),
          error('Polynomial Delta is marginally stable');  
        end;
%------------------------------------------------------
%       Extract matrices Pij
%------------------------------------------------------
        P11 = sys(1:end-1,1:end-1);
        P12 = sys(1:end-1,end);
        P21 = sys(end,1:end-1);
        P22 = sys(end,end);
%------------------------------------------------------
%       Construct D22
%------------------------------------------------------
        D22 = dtfm(P22, T, 0)';
        p22 = allpoles(P22);
        D22 = setpoles(D22, exp(-p22*T) );
        [n,d] = tf2nd ( D22 ); 
%------------------------------------------------------
%       Basic controller
%------------------------------------------------------        
        [a0,b0] = dioph ( n, d, poln(Delta,n.var) );
%------------------------------------------------------
%       Assume that the input is 1/s 
%------------------------------------------------------
        sZpk = zpk(0, [], 1);
        P11 = minreal(P11*sZpk);
        P21 = minreal(P21*sZpk);
        Px = minreal(P12*P21);
        Py = minreal(Px/P22);
        Pa = sumzpk(P11,-Py);
%------------------------------------------------------
%       Number of poles at s=0 or \z=1
%------------------------------------------------------        
        pa = allpoles(Pa);
        za = length(find(abs(pa) < 1e-6));
        z22 = length(find(abs(p22) < 1e-6));
        px  = allpoles(Px);
        zx  = length(find(abs(px) < 1e-6));
        PxF = minreal(Px*zpk(zeros(zx,1),[],1));          
        py  = allpoles(Py);
        zy  = length(find(abs(py) < 1e-6));
        PyF = minreal(Py*zpk(zeros(zy,1),[],1));          
%------------------------------------------------------
%       Number of poles of V2=P12*P21*d^2 at s=0 or \z=1
%------------------------------------------------------        
        z2 = zx - 2*z22;
        if z2 < 0, 
          V2 = zeros(size(Px));
          z2 = 0;
          scale = 1;
        else
          d2R = (-T)^(2*z22)*d*d / poln(ones(2*z22,1),'rz');
          V2 = dcgain(PxF) * polyval(d2R,1) / polyval(Delta,1); 
          scale = polyval(d2R,1);
        end;    
%------------------------------------------------------
%       Number of poles at s=0 or \z=1 of 
%         V1a=(P11-P12*P21/P22) + Py*b0*d 
%------------------------------------------------------        
        z1 = max(za, zy-z22);
        Pa = minreal(Pa*sZpk^z1);
        z0 = min(z22,zy);
        dR = (-T)^z0*d / poln(ones(z0,1),'rz');
        Py2 = minreal(PyF*sZpk^(z1+z22-zy));
        b0Val1 = polyval(b0,1);
        V1 = dcgain(Pa) + ...
            dcgain(Py2)*b0Val1*polyval(dR,1)/polyval(Delta,1); 
%------------------------------------------------------
%       Scaling with account for d(1), which can be near zero
%------------------------------------------------------      
        V1 = V1 / scale;
        V2 = V2 / scale;
%------------------------------------------------------
%       Find DC-gain of Psi
%------------------------------------------------------        
        msgImposs = 'Asymptotic tracking is impossible';
        if z1 == z2
          if norm(V2) > 1e-6
             if norm(V1) < 1e-6, G = NaN;  
             else                G = V2 \ V1; end;
             if norm(V1-V2*G) > 1e-6, error(msgImposs);  end;          
          else
             if norm(V1) < 1e-6, G = NaN;  
             else error(msgImposs);  end;          
          end;
        else            
          error('Special case of asymptotioc tracking');  
        end
%------------------------------------------------------
%       Find DC-gain of controller
%------------------------------------------------------        
        if isnan(G),
          GK = NaN;
        else
          G = real(G);
          num = a0 + d*G;
          den = b0 - n*G;
          [num,den] = coprime(num, den);
          den = polyval(den,1);
          if abs(den) < 1e-10
            GK = Inf;  
          else
            GK = polyval(num,1) / den;
          end;
        end;

%------- End of PSIGAIN.M --------- KYuP ----------           
        
