function [xBest,zBest,iter,xTrace] = infglobc ( func, options )
%INFGLOBC Constrained global optimization using information algorithm. 
%
%	  [X,Y,ITER] = INFGLOBC ( F, OPTIONS )
%
%   Inputs:
%	  F - function to be minimized as a string
%         F(x) assumes that 0 <= x <= 1 
%     OPTIONS - structure with the following fields: 
%          .dim - number of unknown (default 1)
%          .nConstr - number of constraints (default 1)
%          .tol - tolerance (default 1e-4)
%          .r   - optimization parameter (default 2) 
%          .maxIter  - maximal number of iterations (default 500)
%          .dispIter - step of dislaying results (default 50)
%          .display  - flag for dislaying results ('on' or 'off', default 'on')
%
%   Outputs:
%     X    - minimizing solution
%     Y    - minimal function value 
%     ITER - number of iterations
%
%   See also INFGLOB.
	
% [1] R.G.Strongin and Ya. D. Sergeev, Global optimization: fractal approach
%     and non-redundant parallelism, Journal of Global Optimization, 2003, 27(1), 25-50.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Default values
%------------------------------------------------------
        n = 1; tol = 1e-4; r = 2; ksi = 1;
        nConstr = 1;
        maxIter = 500;
        dispIter = 50;
        display = 'off';
%------------------------------------------------------
%       Check syntax
%------------------------------------------------------
        if exist('options','var')
           if isfield(options,'dim'), n = options.dim; end;
           if isfield(options,'tol'), tol = options.tol; end;
           if isfield(options,'r'),   r = options.r; end;
           if isfield(options,'nConstr'), nConstr = options.nConstr; end;
           if isfield(options,'maxIter'), maxIter = options.maxIter; end;
           if isfield(options,'dispIter'), dispIter = options.dispIter; end;
           if isfield(options,'display'), display = options.display; end;
        end
%------------------------------------------------------
%       Initialization
%------------------------------------------------------
        x = [0 1];
        xTrace = x;
        g  = feval(func,0); nu = length(g);       
        g1 = feval(func,1); nu1 = length(g1);
        nu = [nu nu1];
        z = [g(end) g1(end)];
        g(1:nu1,2) = g1;
        epsR = 0.1*ones(nConstr,1);
        mu = ones(nConstr+1,1);
        k = 2;    
        xOld = Inf;
%------------------------------------------------------
%       The best known solution       
%------------------------------------------------------   
        indFeas = find(nu == nConstr+1);
        if isempty(indFeas)
           zBest = Inf; xBest = [];   
        else            
           xFeas = x(indFeas);
           [zBest,ind] = min(z(indFeas));
           xBest = xFeas(ind);
        end;
%------------------------------------------------------
%       Main iteration loop
%------------------------------------------------------   
        iter = 0;
        while 1 
          iter = iter + 1;   
       %-----------------------------------------------
       %   Step 1. Sorting
       %-----------------------------------------------
          [x,ind] = sort(x);  
          g = g(:,ind);
          z = z(ind);
          nu = nu(ind);
       %-----------------------------------------------
       %   Step 2-3. Find I_nu, mu
       %-----------------------------------------------
          IvNext = 1:k;          
          for v=1:nConstr+1
            Iv = IvNext;
          %-----------------------------------------------
          %   Step 3. Find mu(v)
          %-----------------------------------------------
            if length(Iv) < 2, mu(v) = 1;
            else
              muV = 0;  
              for j=1:length(Iv)-1  
                jj = Iv(j);  
                for i=j+1:length(Iv)  
                  ii = Iv(i);  
                  muV = max(muV, abs( (g(v,ii)-g(v,jj)) / (x(ii) - x(jj))^(1/n)) );  
                end;
              end;
              if muV == 0, mu(v) = 1;
              else         mu(v) = muV; end;
            end;
          %-----------------------------------------------
          %   Step 4. Find zAst(v)
          %-----------------------------------------------
            IvNext = find(nu >= v+1);
            if isempty(IvNext)
                 zAst(v) = min(z(Iv)); break; 
            else zAst(v) = - epsR(v); 
            end;
          end;
       %-----------------------------------------------
       %   Step 5. Find R
       %-----------------------------------------------
          Delta = (x(2:end) - x(1:end-1)).^(1/n);
          dz  = z(2:end) - z(1:end-1);
          for i=1:k-1      
            v = max(nu(i), nu(i+1));  
            rMu = r * mu(v);
            rMuDelta = rMu * Delta(i);
            if nu(i) == nu(i+1)
               R(i) = rMuDelta + dz(i)^2/rMuDelta - 2*(z(i+1)+z(i)-2*zAst(v)); 
            elseif nu(i) < nu(i+1)
               R(i) = 2*Delta(i) - 4*(z(i+1) - zAst(v))/rMu; 
            else
               R(i) = 2*Delta(i) - 4*(z(i) - zAst(v))/rMu; 
            end;
          end;
       %-----------------------------------------------
       %   Step 6. Find the closest interval
       %-----------------------------------------------
          [Rmax,t] = max(R);
       %-----------------------------------------------
       %   Step 8. Find the next point
       %-----------------------------------------------
          xNew = (x(t+1) + x(t))/2;
          v = nu(t);
          if v == nu(t+1)
             xNew = xNew - sign(dz(t)) * (abs(dz(t))/mu(v))^(1/n) / (2*r); 
          end;
          zNew = feval(func,xNew);   
          nuNew = length(zNew);
       %-----------------------------------------------
       %   Store the best solution
       %-----------------------------------------------
          if nuNew == nConstr+1  % all constraints hold
             if zNew(end) < zBest
                xBest = xNew;
                zBest = zNew(end);
             end;
          end;
       %-----------------------------------------------
       %   Display the best result
       %-----------------------------------------------
          if isequal(display,'on')  &&  mod(iter,dispIter) == 0
            if isempty(xBest)  
               fprintf('x = ???\n'); 
            else    
               fprintf(sprintf('x = %9.6f   f = %g\n', xBest, zBest));  
            end;
          end;
       %-----------------------------------------------
       %   Check finalization condition
       %-----------------------------------------------
          if iter >= maxIter
             warning('Solution not found in %g iterations',maxIter); 
             break; 
          end;
          if abs(x(t+1)-x(t)) < tol, break; end;
%          if ~isempty(xBest)
%             if abs(xNew - xOld) < tol, break; end;
%          end;
          xOld = xNew;
       %-----------------------------------------------
       %   Augment arrays
       %-----------------------------------------------
          x = [x xNew]; 
          xTrace = [xTrace xNew];
          k = k + 1;
          g(1:nuNew,k) = zNew;
          z = [z zNew(end)];
          nu = [nu length(zNew)];
        end;
%------- End of INFGLOBC.M --------- KYuP ----------
