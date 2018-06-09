function [N,J0] = sdnorm ( csys, K, type, tau )
%SDNORM Compute a norm for sampled-data system.
%
%     N = SDNORM ( CSYS, K, TYPE, TAU )
%
%   Inputs:
%     CSYS - continuous-time standard system 
%	  K    - discrete controller (positive feedback)
%	  TYPE - type of norm:
%         'gh2'  - generalized H2-norm (default)
%         'sh2'  - simple H2-norm
%         'inf'  - classical Hinf-norm (L2-induced norm)
%         'ainf' - associated Hinf-norm
%     TAU  - input time-delay for simple H2-problem (default 0), 
%            possibly a vector 
%
%   Outputs:
%	  N - required norm of the closed-loop system
%
%   See also SDH2NORM, SDHINORM, SDAHINORM.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 29-Oct-2006 $
%------------------------------------------------------
%       Check syntax
%------------------------------------------------------
        if ~exist('type','var'), type = 'gh2'; end;
        switch type
          case {'gh2', 'h22', 'sh2', 'inf', 'ainf', 'ainfx'},  
          otherwise, error('Unknown norm ''%s''',type);
        end;
%------------------------------------------------------
%       Preliminaries
%------------------------------------------------------
        tol = 1e-3;
        T = K.Ts;    
        [nCon,nMeas] = size(K);
        i1 = size(csys,2) - nCon;  
        o1 = size(csys,1) - nMeas;  
        csys = ss ( csys );
        B1 = csys.b(:,1:i1);   
        C1 = csys.c(1:o1,:);   
%------------------------------------------------------
%       Special case: simple H2-norm
%------------------------------------------------------
        if isequal(type,'sh2')
           if ~exist('tau','var'), tau = 0; end;
           N = zeros(size(tau)); J0 = N;
           for i=1:length(tau)
             if tau(i) > T, 
                warning('Tau > T, assuming tau=0.'); 
                tau(i) = 0;
             end;
            %------------------------------------------- 
            %  Equivalent discrete model
            %------------------------------------------- 
             dsys = sdh2simple ( csys, T, nMeas, nCon, tau(i) );
            %------------------------------------------- 
            %  Find output norm of closed-loop system
            %------------------------------------------- 
             tol = 1e-8;
             dcl = lft( dsys, K );	    
             while 1
                dcl = minreal ( dcl, tol, 0 );	    
                if all(abs(eig(dcl.a)) < 1), break; end;
                tol = tol * 10;
                if tol > 1e-2, break; end;
             end;
             N(i) = norm(dcl);
           end;
           return;
        end;
%------------------------------------------------------
%       Special case: classical Hinf-norm
%------------------------------------------------------
        if isequal(type,'inf')
           N = sdhinorm ( csys, K );
           return;
        end;
%------------------------------------------------------
%       Construct H2-equivalent discrete model
%       for closed-loop system
%------------------------------------------------------
        [dsys,addVar,dsysDelta,G11,G12,G21,G22,G0] = sdgh2mod ( csys, T );        
%------------------------------------------------------
%       Special case: H2-norm 
%------------------------------------------------------
        if isequal(type,'gh2')
           %dclDelta = minreal ( lft( dsysDelta, K ), tol, 0 );	    
           %N = sqrt(norm(dclDelta)^2+T*norm(G0)^2-norm(G11)^2);
           tol = 1e-8;
           dcl = lft( dsys, K );	    
           while 1
              dcl = minreal ( dcl, tol, 0 );	    
              if all(abs(eig(dcl.a)) < 1), break; end;
              tol = tol * 10;
              if tol > 1e-2, break; end;
           end;
           N = norm(dcl);
           return;
        end;
%------------------------------------------------------
%       Check critical poles of G0 and G11
%------------------------------------------------------        
        dclDelta = minreal ( lft( dsysDelta, K ), tol, 0 );	    
        eig0 = eig(G0.a);
        critG0 = any(abs(abs(eig0)-1) < tol);
%------------------------------------------------------
%       Compute T*trace(G0'*G0) - trace(G11'*G11)
%           or  T*G0'*G0 - G11'*G11 for 'ahinfx' (if possible)
%------------------------------------------------------        
        if isequal(type,'ainfx')
          try
             PP = minreal(T*G0'*G0, -G11'*G11, sqrt(eps), 0);
          catch
             G0 = compadibleDim ( G0, G11 ); 
             G0 = zpk(G0); G11 = zpk(G11);  
             PP = sumzpk(T*G0'*G0, -G11'*G11);
          end;
        else %if ~isequal(type,'ainfx')
          try
             PP = minreal(T*trace(G0'*G0), -trace(G11'*G11), sqrt(eps), 0);
          catch
             G0 = zpk(G0); G11 = zpk(G11);  
             PP = sumzpk(T*trace(G0'*G0), -trace(G11'*G11));
          end;
        end;
%------------------------------------------------------
%       Remove critical poles
%------------------------------------------------------        
        if critG0
           PP = removeCritical ( PP, 1e-6 ); 
           try   PP = ss ( PP );
           catch dclDelta = zpk(dclDelta); end;
        end;
%------------------------------------------------------
%       AHinfX-norm is the square root of ||X||_\infty
%------------------------------------------------------
        if isequal(type,'ainfx')
          try  
            X = minreal(dclDelta'*dclDelta + PP, 1e-6, 0 );            
          catch
            PP = zpk(PP); dclDelta = zpk(dclDelta); 
            X = minreal(dclDelta'*dclDelta + PP);            
          end;
	      N = sqrt(norm(X,'inf')) * sqrt(T);
        end;
%------------------------------------------------------
%       AHinf-norm is the square root of ||trace X||_\infty
%------------------------------------------------------
        if isequal(type,'ainf')
          try  
            X = minreal(trace(dclDelta'*dclDelta) + PP, 1e-6, 0 );            
          catch
            PP = zpk(PP); dclDelta = zpk(dclDelta); 
            X = minreal(trace(dclDelta'*dclDelta) + PP);            
          end;
	      N = sqrt(norm(X,'inf')) * sqrt(T);
        end;
%------------------------------------------------------
%       Special case: H2-norm, type II formulas 
%       (unreliable for models with unit-circle poles)
%------------------------------------------------------
        if isequal(type,'h22')
           if ~critG0
              N = sqrt(norm(dclDelta)^2 + T*norm(G0)^2 - norm(G11)^2);
           else
              FP = minreal(sfactz(PP), 1e-4);
              N = sqrt(norm(dclDelta)^2 + norm(FP)^2);
           end;
        end;

%------- End of SDNORM.M --------- KYuP ----------

%#################################################
function G0 = compadibleDim ( G0, G11 )
%------------------------------------------------------        
%
% COMPATIBLEDIM add extra zero column to G0
%
%------------------------------------------------------                
        c0  = size(G0,2);
        c11 = size(G11,2);
        diff = c11 - c0;
        if diff > 0            
           [a,b,c,d,Ts] = ssdata(G0);
           b = [b zeros(size(b,1),diff)];
           d = [d zeros(size(d,1),diff)];
           G0 = ss(a, b, c, d, Ts);
        end;

%------- End of COMPATIBLEDIM.M --------- KYuP ----------

%#################################################
function F = removeCritical ( F0, tol )
%------------------------------------------------------        
%
% REMOVECRITICAL remove unit-circle poles by
%                approximate cancellation
%------------------------------------------------------                
        F = minreal(zpk(F0), tol);
%------------------------------------------------------                
%       MIMO-case - for each term separately        
%------------------------------------------------------                
        [rF,cF] = size(F);
        if any([rF cF] > 1),
           for i=1:rF
             for j=1:cF
               F(i,j) = removeCritical ( F(i,j), tol );
             end;
           end;
           return;
        end;
%------------------------------------------------------                
%       SISO-case
%------------------------------------------------------                
        [zz,pp,K,T] = zpkdata(F, 'v');
        critNo = find(abs(abs(pp)-1) < tol)';
        if isempty(critNo)
           F = F0;
        else
           cancErr = 0;
           for i=critNo
             [xx,ind,ee] = nearest ( pp(i), zz );
             if ~isempty(ind)
                zz(ind) = []; 
                cancErr = cancErr + ee;
             end;
           end;
           pp(critNo) = [];
           F = zpk(zz, pp, K, T);
           %if cancErr > 1e-3, 
           %   warning(sprintf('Incorrect cancellation %g',cancErr)); 
           %end;
        end;
           
%------- End of REMOVECRITICAL.M --------- KYuP ----------

%#################################################
function [y0,ind,ee] = nearest ( x, y )
%-------------------------------------------------
%
% NEAREST find the nearest element in Y to X
%         and its index in Y
%
%-------------------------------------------------
        err = abs(y - x);
        [ee,ind] = min(err);
        if ~isempty(ind), 
             y0 = y(ind);
        else y0 = []; 
        end;
        
%------- End of NEAREST.M --------- KYuP ----------
