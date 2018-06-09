function [xBest,zBest,iter,xTrace] = infglob ( funf, opt )
%INFGLOB Global optimization using information algorithm.
%
%	  [X,Y,ITER] = INFGLOB ( F, OPTIONS )
%
%   Inputs:
%	  F - function to be minimized as a string
%         F(x) assumes that 0 <= x <= 1 
%     OPTIONS - structure with the following fields: 
%            .display  - verbosity mode: 'on', 'off' or 'final'
%            .tol      - tolerance (default 1e-4)
%            .maxIter  - maximal number of iterations (default 500)
%            .dispIter - step of dislaying results (default 10)
%            .dim      - number of unknown (default 1)
%            .r        - optimization parameter (default 2) 
%            .ksi      - optimization parameter (default 1)
%
%   Outputs:
%	  X    - minimizing solution
%     Y    - minimal function value 
%     ITER - number of iterations
%
%   See also INFGLOBC.
	
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Default values
%------------------------------------------------------
        options = struct('display',    'on', ...
                         'tol',        1e-4, ...
                         'maxIter',    500, ...
                         'dispIter',   10, ...                        
                         'dim',        1, ...
                         'r',          2, ...
                         'ksi',        1, ...
                         'guess',      [] ...
                         );
        if nargin==1 && nargout <= 1 && isequal(funf,'defaults')
           xBest = options;
           return
        end;        
%------------------------------------------------------
%	    Load options
%------------------------------------------------------
        if exist('opt','var')  &&  ~isempty(opt)
           options = updateopt ( options, opt );
        end;
        display  = options.display; 
        tol      = options.tol;
        maxIter  = options.maxIter;
        dispIter = options.dispIter;
        n        = options.dim;
        r        = options.r;
        ksi      = options.ksi;
        guess    = options.guess;
%------------------------------------------------------
%	    Determine verbosity level
%------------------------------------------------------
        switch display
          case 'off',   verbosity = 0;
          case 'final', verbosity = 1;
          case 'on',    verbosity = 2;
          otherwise,    error('Unknown verbosity level ''%s''', display);
        end        
        if verbosity > 1, 
           disp('# Information global search:');
           disp('      Iter       xBest        FunVal');
           disp('--------------------------------------');
        end;        
%------------------------------------------------------
%       Initialization
%------------------------------------------------------
        x = (0:0.1:1);
        if ~isempty(guess), x = sort([x guess]); end;
        z = x;
        k = length(x);    
        xTrace = x;
        for i=1:k, z(i) = feval(funf,x(i)); end;
        xOld = Inf;
%------------------------------------------------------
%       Initialize loop variables
%------------------------------------------------------
        iter = 0;
        lam = zeros(k+maxIter,1);
        [zBest,ind] = min(z);
        xBest = x(ind);
%######################################################
%       Main iteration loop
%######################################################
        while 1 
          iter = iter + 1;   
          if iter > maxIter, break; end;
       %-----------------------------------------------
       %   Step 3. Sorting
       %-----------------------------------------------
          [x,ind] = sort(x);  
          z = z(ind);
       %-----------------------------------------------
       %   Step 4. Find X, mu, gamma, lam
       %-----------------------------------------------
          dxn = (x(2:end) - x(1:end-1)).^(1/n);
          dz  = z(2:end) - z(1:end-1);
          adz = abs(dz);
          zDx = adz ./ dxn;
          X = max(dxn);
          mu = max(zDx);
          gamma = mu*dxn/X;
          zDx = [zDx(1) zDx zDx(end)];
          for i=1:k-1
            lam(i) = max(zDx(i:i+2));  
          end;
       %-----------------------------------------------
       %   Step 5. Estimate the local Lipschitz constants
       %-----------------------------------------------
          hatL = max( [gamma; lam(1:k-1)'; ksi*ones(1,k-1)] );
          rLx = r*hatL.*dxn;
       %-----------------------------------------------
       %   Step 6. Find R
       %-----------------------------------------------
          R = rLx - 2*(z(2:end) + z(1:end-1)) + dz.^2 ./ rLx;
%          if any(isinf(R)), keyboard; end;
       %-----------------------------------------------
       %   Step 7. Find the closest interval
       %-----------------------------------------------
          [Rmax,t] = max(R);
       %-----------------------------------------------
       %   Step 8. Find the next point
       %-----------------------------------------------
          xNew = (x(t+1) + x(t))/2 - ...
                 sign(dz(t)) * (adz(t)/hatL(t))^n / (2*r);
          zNew = feval(funf,xNew);
          %disp(zNew)
       %-----------------------------------------------
       %   Store the best solution
       %-----------------------------------------------
          if zNew < zBest
             xBest = xNew;
             zBest = zNew;
          end;
       %-----------------------------------------------
       %   Display the best result
       %-----------------------------------------------
          if verbosity > 1  &&  mod(iter,dispIter) == 0
            fprintf('%10d   %10.5f   %10.5g\n', iter, xBest, zBest);  
          end;
       %-----------------------------------------------
       %   Check finalization condition
       %-----------------------------------------------
          if abs(x(t+1)-x(t)) < tol, break; end;
%          if abs(xNew - xOld) < tol, break; end;
          xOld = xNew;
       %-----------------------------------------------
       %   Augment arrays
       %-----------------------------------------------
          x = [x xNew]; 
          xTrace = [xTrace xNew];
          z = [z zNew];
          k = k + 1;
        end;
%######################################################
%       End of Main iteration loop
%######################################################
        if verbosity > 0
           if iter > maxIter,
                disp(sprintf('Solution not found in %d iterations.', maxIter)); 
           else disp('Optimization completed successfully.'); 
           end;
        end;

%------- End of INFGLOB.M --------- KYuP ----------
