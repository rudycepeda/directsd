function [xBest,yBest,func_evals] = sasimplex (func, x0, opt)
%SASIMPLEX Simulated annealing using Nelder-Mead method. 
%
%     [X,Y,N] = SASIMPLEX ( FUNC, X0, OPTIONS )
%
%   Inputs:
%     FUNC - function to minimize
%     X0   - initial guess
%     OPTIONS - options, structure with the following possible fields
%             .display     - verbosity mode: 'on', 'off' or 'final'
%             .tol         - parameter tolerance (1e-4)
%             .maxFunEvals - maximal number of function call (default 1000)
%             .dispIter    - step of dislaying results (default 10)
%             .multiStep   - number of iterations with the same temperature (10)
%             .startTemp   - starting Boltsman temperature (100)
%             .tempDecRate - rate of temperature decrease (0.01)  
%   Outputs:
%     K   - modal H2-optimal discrete controller of order N
%     ERR - cost function for the optimal system
%
%   See also SIMANNEAL.

% References:
% [1] Numerical Recipes, Section 10.9.
%
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%%       Initialize parameters
%------------------------------------------------------
        options = struct('display',     'on', ...
                         'tol',         1e-5, ...
                         'maxFunEvals', 10000, ...
                         'dispIter',    100, ...                        
                         'multiStep',   20, ...                        
                         'startTemp',   100, ...
                         'tempDecRate', [] ...
                         );
        func_evals = 0;
        iter = 0;
        rho = 1; chi = 2; psi = 0.5; sigma = 0.5;
%------------------------------------------------------
%	    Load options
%------------------------------------------------------
        if exist('opt','var')  &&  ~isempty(opt)
           options = updateopt ( options, opt );
        end;
%------------------------------------------------------       
        display   = options.display;
        tol       = options.tol;
        maxFunEvals = options.maxFunEvals;
        dispIter  = options.dispIter;        
        M         = options.multiStep;
        T         = - options.startTemp;
        alpha     = options.tempDecRate;
        T0 = T;
        if isempty(alpha), alpha = - log(0.1)/200; end;
%------------------------------------------------------
%%	    Determine verbosity level
%------------------------------------------------------
        switch display
          case 'off',   verbosity = 0;
          case 'final', verbosity = 1;
          case 'on',    verbosity = 2;
          otherwise,    error('Unknown verbosity level ''%s''', display);
        end        
        if verbosity > 1, 
           disp('# Simulated annealing procedure:');
           disp('      Iter       FunVal         T   ');
           disp('--------------------------------------');
        end;        
%------------------------------------------------------
%%       Initialize arrays
%------------------------------------------------------
        x0 = x0(:);
        n = length(x0);
        v = zeros(n,n+1); 
        y = zeros(1,n+1);
%------------------------------------------------------
%%      Generate a simplex near the initial guess
%------------------------------------------------------
        x(:,1) = x0;
        y(1) = feval(func, x0);
%------------------------------------------------------
        for j = 1:n
          z = x0;
          if z(j) ~= 0
               z(j) = 1.05*z(j);    % 5% shift
          else z(j) = 0.00025; 
          end;
          x(:,j+1) = z;
          y(j+1) = feval(func, y);
        end
%------------------------------------------------------
        xBest = x0;
        yBest = y(1);
%######################################################
%   Main iteration loop
%######################################################
        while max(max(abs(x(:,2:n+1)-x(:,ones(1,n))))) > tol
         %---------------------------------------------   
         %  Stochastic fluctuation for all y's except for the best point
         %---------------------------------------------   
          [y,ind] = sort ( y );
          x = x(:,ind);
          yFlu = [y(1) y(2:end)+T*log(rand(1,n))];
         %---------------------------------------------   
         %  Sort for increasing function value
         %---------------------------------------------   
          [yFlu,ind] = sort(yFlu);
          y = y(ind);
          x = x(:,ind);
         %---------------------------------------------   
         %  Calculate center using n (NOT n+1) best points
         %---------------------------------------------   
          xbar = sum(x(:,1:n),2)/n;
         %---------------------------------------------   
         %  Perform reflection of the worst point
         %---------------------------------------------   
          xr = (1 + rho)*xbar - rho*x(:,end);
          fxr = feval(func, xr);
          fxrFlu = fxr - T*log(rand);
          func_evals = func_evals+1;     
         %---------------------------------------------   
         %  If reflected point is the best of all ...
         %---------------------------------------------   
          if fxrFlu < yFlu(1) 
           %---------------------------------------------   
           %  Try expansion (extrapolation)
           %---------------------------------------------   
            xe = (1 + rho*chi)*xbar - rho*chi*x(:,end);
            fxe = feval(func, xe);
            fxeFlu = fxe - T*log(rand);
            func_evals = func_evals+1;
            if fxeFlu < fxrFlu
                 x(:,end) = xe; y(end) = fxe; how = 'expand';
            else x(:,end) = xr; y(end) = fxr; how = 'reflect'; 
            end;
         %---------------------------------------------   
         %  If reflected point is NOT the best of all ...
         %---------------------------------------------   
          else % y(1) <= fxr
            if fxrFlu < yFlu(n)
              x(:,end) = xr; y(end) = fxr; how = 'reflect';
           %---------------------------------------------   
           %  If reflected point is between y(n) and y(n+1) ...
           %---------------------------------------------   
            else % fxr >= y(n)
             %-------------------------------------------   
             % Perform contraction
             %-------------------------------------------   
              if fxrFlu < yFlu(end)
               %-------------------------------------------   
               % Outside contraction
               %-------------------------------------------   
                xc = (1 + psi*rho)*xbar - psi*rho*x(:,end);
                fxc = feval(func, xc);
                fxcFlu = fxc - T*log(rand);
                func_evals = func_evals+1;    
                if fxcFlu <= fxrFlu
                     x(:,end) = xc; y(end) = fxc; 
                     how = 'contract outside';
                else how = 'shrink'; 
                end;
              else
               %-------------------------------------------   
               % Inside contraction
               %-------------------------------------------   
                xcc = (1-psi)*xbar + psi*x(:,end);
                fxcc = feval(func, xcc);
                fxccFlu = fxcc - T*log(rand);
                func_evals = func_evals+1;   
                if fxccFlu < yFlu(end)
                     x(:,end) = xcc; y(end) = fxcc; 
                     how = 'contract inside';
                else how = 'shrink'; 
                end;
              end
            %---------------------------------------------   
            %  Shrinking ...
            %---------------------------------------------   
              if isequal(how,'shrink')
                for j=2:n+1
                  x(:,j) = x(:,1)+sigma*(x(:,j) - x(:,1));
                  y(j) = feval(func, x(:,j));
                end
                func_evals = func_evals + n;
              end
            end
          end
         %------------------------------------------------
         %   Change the best solution
         %------------------------------------------------
          [ymin,ind] = min(y);
          if ymin < yBest, 
             yBest = ymin;
             xBest = x(:,ind);
          end;
         %------------------------------------------------
         %   Iteration count
         %------------------------------------------------
          iter = iter + 1;          
          if func_evals >= maxFunEvals, 
             warning('Maximal number of function evaluations exceeded'); 
             break; 
          end;
         %--------------------------------------------
         % [1] L.Ingber, Very Fast Simulated Re-Annealing,
         %     J. Mathl. Comput. Modelling, vol. 12, pp.967-973,
         %     1989.
         %--------------------------------------------
           T = T0 * exp(-alpha*iter); 
         %--------------------------------------------
         %  Change temperature after each M steps
         %--------------------------------------------          
          %if mod(iter,M) == 0, T = T * alpha; end;
          if mod(iter,dispIter) == 0,
             showinfo ( iter, yBest, T ) 
          end;
        end;
%######################################################
%   END of Main iteration loop
%######################################################
        
%------- End of SASIMPLEX.M --------- KYuP ----------           

%####################################################
function showinfo ( iter, fval, T )
      fprintf('%10d   %10.5f   %10.5f\n', iter, fval, T );

%------- End of SHOWINFO.M --------- KYuP ----------
