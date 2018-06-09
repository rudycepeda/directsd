function [xBest,valBest,func_evals,xTrace] = simanneal ( funf, x, opt, fung, funproj )
%SIMANNEAL Direct search using simulated annealing technique.
%
%     [X,FVAL,N] = SIMANNEAL ( FUNF, X0, OPTIONS, FUNG, FUNPROJ )
%
%   Inputs: 
%     FUNF - cost function to be minimized
%     X0   - initial guess for parameter values
%     OPTIONS - options, structure with the following possible fields
%            .display     - verbosity mode: 'on', 'off' or 'final'
%            .tol         - parameter tolerance (1e-4)
%            .maxFunEvals - maximal number of function call (default 1000)
%            .dispIter    - step of dislaying results (default 10)
%            .iniStep     - length of initial step in n-dimensional space (0.1)
%            .multiStep   - number of iterations with the same temperature (10)
%            .maxFail     - number of unsuccessful iterations before step decrease
%            .decStepBy   - step size multiplier if maxFail steps were unsuccessful (0.5)
%            .decFailRate - rate of decreasing maxFail when step decreases (1.1)
%            .adaptRate   - adaptation (self-learning) rate (0.1) 
%            .startTemp   - starting Boltsman temperature (100)
%            .tempDecRate - rate of temperature decrease (0.95)  
%     FUNG    - constraints function, returns constraint values that
%               should be negative
%     FUNPROJ - function that projects point onto admissible region
%
%   Outputs:
%     X    - optimal parameter values
%     FVAL - minimal value of the cost function
%     ITER - number of iteration performed
%
%   See also SASIMPLEX.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%	    Default options
%------------------------------------------------------
        options = struct('display',     'on', ...
                         'tol',         1e-4, ...
                         'maxFunEvals', 1000, ...
                         'dispIter',    100, ...                        
                         'iniStep',     0.1, ...                        
                         'multiStep',   20, ...
                         'maxFail',     100, ...
                         'decStepBy',   0.5, ...
                         'decFailRate', 1.1, ...
                         'adaptRate',   0.1, ...
                         'startTemp',   100, ...
                         'tempDecRate', [] ...
                         );
        if nargin==1 && nargout <= 1 && isequal(funf,'defaults')
           xBest = options;
           return
        end;       
%------------------------------------------------------
%	    Parameters
%------------------------------------------------------
        x = x(:);
        len = length(x);
        funf = fcnchk ( funf );
%------------------------------------------------------       
        withConstr = exist('fung','var');
        if withConstr, withConstr = ~isempty(fung); end;
        if withConstr, fung = fcnchk (fung); end;
%------------------------------------------------------
        withProj   = exist('funproj','var');
        if withProj, withProj = ~isempty(funproj); end;
        if withProj, funproj = fcnchk (funproj); end;
%------------------------------------------------------
%	    Load options
%------------------------------------------------------
        if exist('opt','var')  &&  ~isempty(opt)
           options = updateopt ( options, opt );
        end;
%------------------------------------------------------       
        display     = options.display;
        tol         = options.tol;
        maxFunEvals = options.maxFunEvals;
        dispIter  = options.dispIter;        
        step      = options.iniStep;       
        M         = options.multiStep;
        maxFail   = options.maxFail;
        decStepBy = options.decStepBy;
        decFail   = options.decFailRate;
        adaptRate = options.adaptRate;
        T         = options.startTemp;
        alpha     = options.tempDecRate;
%------------------------------------------------------
%	    Very fast SA
%------------------------------------------------------
        T0 = T;
        if isempty(alpha), alpha = - log(0.1)/200; end;
        beta = 350; % standard p = exp(-3.5*delta/val)
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
           disp('# Simulated annealing procedure:');
           disp('      Iter       FunVal         Step');
           disp('--------------------------------------');
        end;        
%------------------------------------------------------
%	    Initial values
%------------------------------------------------------
        val = feval ( funf, x );
        func_evals = 1;
        xBest = x;
        valBest = val;        
%------------------------------------------------------
%	    Prepare adaptation arrays
%------------------------------------------------------
        memLen = 10;
        dJ = zeros(1, memLen);
        dxAdapt = zeros(len, memLen);
        adaptDir = zeros(len, 1);       
        adaptRate = 0;
%------------------------------------------------------
%	    Initialize loop data
%------------------------------------------------------
        failCount = 0;
        iter = 0;
        if verbosity > 1, 
           showinfo(iter, valBest, step, '  '); 
        end;    
        xTrace = [];
%######################################################
%	    Main loop
%######################################################
        while step > tol            
           iter = iter + 1;
           if func_evals > maxFunEvals, break; end;
           bShown = 0;
          %--------------------------------------------
          %  Choose random direction
          %--------------------------------------------          
           dirX = adaptDir + rand(len,1) - 0.5;
           dirX = dirX / norm(dirX);
           dx = dirX * step * (0.5+rand);
           xNew = x + dx;           
          %--------------------------------------------
          %  Check constraints
          %--------------------------------------------          
           if withProj, xNew = feval ( funproj, xNew ); end;
           if withConstr
              if any(feval(fung, xNew) > 0),
                 failCount = failCount + 1; 
                 if failCount >= maxFail
                    step = step * decStepBy; 
                    if verbosity > 1, 
                       showinfo(iter, valBest, step, ' >');
                       bShown = 1;
                    end;
                    failCount = 0;
                 end;
                 dJ = zeros(1, memLen);
                 dxAdapt = zeros(len, memLen);
                 adaptDir = zeros(len, 1);
                 continue;
              end;
           end;
          %--------------------------------------------
          %  Compute function value
          %--------------------------------------------          
           dx = xNew - x;
           valNew = feval ( funf, xNew ); 
           func_evals = func_evals + 1;
           %fprintf('x = %g,  f = %g\n', xNew, valNew);
           %fprintf('f = %g\n', valNew);pause
           xTrace = [xTrace; xNew.' 1];
          %--------------------------------------------
          %  Adaptation
          %--------------------------------------------          
           dxAdapt = [dx dxAdapt];
           dxAdapt(:,end) = [];
           dJ = [valNew-val dJ];
           dJ(end) = [];           
           adaptDir = dxAdapt * diag(dJ) * diag(memLen:-1:1);
           adaptDir = sum(adaptDir,2);
           if norm(adaptDir) > sqrt(eps)
              adaptDir = adaptRate * adaptDir / norm(adaptDir);
           end;           
          %--------------------------------------------
          %  Accept step or not?
          %--------------------------------------------          
           failCount = failCount + 1;  
           if valNew < val,     % accept surely, downhill
              if valNew < valBest
                 xBest = xNew;
                 valBest = valNew;
                 failCount = 0;
                 if verbosity > 1, showinfo(iter, valBest, step, '> '); bShown = 1; end;
              end;   
              x = xNew;
              val = valNew;
            else                    % accept probablistically
              prob = exp(-beta*(valNew - val)/T);
              if rand < prob, 
                 x = xNew;
                 val = valNew;
              else 
                 xTrace(end,len+1) = 0;
              end;    
              if failCount >= maxFail
                 x = xBest;
                 val = valBest;
                 step = step * decStepBy; 
                 maxFail = max(20, floor(maxFail/decFail));
                 if verbosity > 1, showinfo(iter, valBest, step, ' >'); bShown = 1; end;
                 failCount = 0;
              end;
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
          %----------------------------------------------
          %  Show current results
          %----------------------------------------------
           if verbosity > 1  &&  mod(iter,dispIter) == 0  &&  ~bShown
             showinfo(iter, valBest, step, '  ');   
           end;
        end;
%######################################################
%	    End of Main loop
%######################################################
        if verbosity > 0
           if func_evals > maxFunEvals,
                disp(sprintf('Solution not found in %d function call.', maxFunEvals)); 
           else disp('Optimization completed successfully.'); 
           end;
        end;
        
%------- End of SIMANNEAL.M --------- KYuP ----------

%####################################################
function showinfo ( iter, fval, step, marks )
      fprintf('%10d   %s%10.5g   %s%10.5f\n', iter, marks(1), fval, marks(2), step );

%------- End of SHOWINFO.M --------- KYuP ----------
