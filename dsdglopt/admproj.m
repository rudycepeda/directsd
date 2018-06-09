function p = admproj ( p, alpha, beta, dom, T )
%ADMPROJ Proection of poles onto an admissible region.
%        Given degree of stability and damping factor.
%
%	  PA = ADMPROJ ( P, ALPHA, BETA, DOM, T )
%
%   Inputs:
%     P     - initial poles
%     ALPHA - limiting degree of stability
%     BETA  - limiting damping factor
%             BETA>0 for truncated sector and
%             BETA<0 for shifted sector
%     DOM   - domain of 's' (default), 'z', or 'd'
%     T     - sampling period (default 1 for discrete-time systems)
%
%   Outputs:
%     PR    - projected poles 
%
%   See also MODSDH2, MODSDL2.
	
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if ~exist('alpha','var'),  alpha = 0; end;
        if ~exist('beta','var'),   beta = Inf; end;
        if ~exist('dom','var'),    dom = 's'; end; 
        if ~exist('T','var'), T = 1; end;
        Shifted = (beta < 0);
        beta = abs(beta);
%------------------------------------------------------
%       Transform to s-plane
%------------------------------------------------------
        if dom == 'z'       
          ind = find(abs(p) < 1e-12);
          pZero = zeros(length(ind),1);
          p(ind) = [];
        end;
        switch dom
          case 's', 
          case 'z', p = log(p) / T;
          case 'd', p = - log(p) / T;
          otherwise, error('Incorrect domain ''%s''',dom);
        end;
%------------------------------------------------------
%       Project poles
%------------------------------------------------------
        p = p + Shifted*alpha;
        alpha0 = (1-Shifted)*alpha; 
        for i=1:length(p)
          re = real(p(i)); im = imag(p(i));   
          re = min(re, -alpha0);
          betaX = abs(im/re);  
          if betaX > beta
            if dom ~= 's'  &&  im == pi  % real poles in 'z' or 'zeta'
              re = - pi/beta;  
            else
              im = im*beta/betaX;
            end; 
          end;
          p(i) = re+sqrt(-1)*im;
        end;
%------------------------------------------------------
%       Finalization
%------------------------------------------------------
        p = p - Shifted*alpha;
        if dom ~= 's'
          ind = find(imag(p) == pi);  
          switch dom
            case 'z', p = exp(p*T);
            case 'd', p = exp(-p*T);
          end;
          p(ind) = real(p(ind));
          if dom == 'z'       
            p = [p; pZero];  
          end;
        end;        

%------- End of ADMPROJ.M --------- KYuP ----------
