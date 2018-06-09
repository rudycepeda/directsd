function [x,y,func_evals] = neldermead (func, x0, tolx)
%NELDERMEAD Simplex minimization method by Nelder and Mead.
%
%     [X,Y,N] = NELDERMEAD ( FUNC, X0, TOLX )
%
%   Inputs:
%     FUNC - function to minimize
%     X0   - initial guess
%     TOLX - tolerance with respect to x
%
%   Outputs:
%     K   - modal H2-optimal discrete controller of order N
%     ERR - cost function for the optimal system
%     N   - number of function evaluation 
%
%   See also MODSDH2, MODSDL2.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%%       Initialize parameters
%------------------------------------------------------
        if ~exist('tolx','var'), tolx = 1e-5; end; 
        MAXITER = 1000;
        MAXFUNEVAL = 5000;
        func_evals = 0;
        itercount = 0;
        rho = 1; chi = 2; psi = 0.5; sigma = 0.5;
%------------------------------------------------------
%%       Initialize arrays
%------------------------------------------------------
        x0 = x0(:);
        n = length(x0);
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
          y(j+1) = feval(func, z);
        end
%######################################################
%   Main iteration loop
%######################################################
        while max(max(abs(x(:,2:n+1)-x(:,ones(1,n))))) > tolx
         %---------------------------------------------   
         %  Sort for increasing function value
         %---------------------------------------------   
          [y,ind] = sort(y);
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
          func_evals = func_evals+1;     
         %---------------------------------------------   
         %  If reflected point is the best of all ...
         %---------------------------------------------   
          if fxr < y(1)
           %---------------------------------------------   
           %  Try expansion (extrapolation)
           %---------------------------------------------   
            xe = (1 + rho*chi)*xbar - rho*chi*x(:,end);
            fxe = feval(func, xe);
            func_evals = func_evals+1;
            if fxe < fxr
                 x(:,end) = xe; y(end) = fxe; how = 'expand';
            else x(:,end) = xr; y(end) = fxr; how = 'reflect'; 
            end;
         %---------------------------------------------   
         %  If reflected point is NOT the best of all ...
         %---------------------------------------------   
          else % y(1) <= fxr
            if fxr < y(n)
              x(:,end) = xr; y(end) = fxr; how = 'reflect';
           %---------------------------------------------   
           %  If reflected point is between y(n) and y(n+1) ...
           %---------------------------------------------   
            else % fxr >= y(n)
             %-------------------------------------------   
             % Perform contraction
             %-------------------------------------------   
              if fxr < y(end)
               %-------------------------------------------   
               % Outside contraction
               %-------------------------------------------   
                xc = (1 + psi*rho)*xbar - psi*rho*x(:,end);
                fxc = feval(func, xc);
                func_evals = func_evals+1;    
                if fxc <= fxr
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
                func_evals = func_evals+1;   
                if fxcc < y(end)
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
          itercount = itercount + 1;
          if itercount >= MAXITER, break; end;
          if func_evals >= MAXFUNEVAL, break; end;
        end;
%######################################################
%   END of Main iteration loop
%######################################################
        x = x(:,1);
        y = min(y);

%------- End of NELDERMEAD.M --------- KYuP ----------           
