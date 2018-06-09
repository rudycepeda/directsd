function [x,valBest,func_evals] = arandsearch ( funf, x, opt, fung, funproj )
%ARANDSEARCH Accelerated random search optimization.
%          
%    [X,VAL,ITER] = ARANDSEARCH ( FUNF, X0, OPTIONS, FUNG, FUNPROJ )
%                                       
%  Inputs:
%    FUNF    - name of a fitness function
%    X0      - initial guess
%    OPTIONS - options, structure with the following possible fields
%           .display     - verbosity mode: 'on', 'off' or 'final'
%           .tol         - parameter tolerance (1e-4)
%           .maxFunEvals - maximal number of function call (default 1000)
%           .dispIter    - step of dislaying results (default 10)
%           .iniStep     - length of initial step in n-dimensional space (0.1)
%           .multiStep   - number of trial points in an iteration (10)
%           .maxFail     - number of unsuccessful iterations before step change (10)
%           .decStepBy   - step size multiplier if maxFail steps were unsuccessful (0.5)
%           .maxSuccess  - number of successful iterations before step increase
%           .incStepBy   - step size multiplier if the current step is successful (1.05)
%           .adaptRate   - adaptation (self-learning) rate (0.1)
%    FUNG    - name of a constraints function
%    FUNPROJ - function that projects point onto admissible region
%
%  Outputs:
%    X   - solution found
%    VAL - minimal attainable value
%    N   - number of function evaluation 
%
%   See also RANDSEARCH, CRANDSEARCH.

% [1] M.J. Appel, R. Labarre, D. Radulovic, On Accelerated 
%     Random Search, SIAM Journal of Optimization, 
%     2004, vol. 14; part 3, pages 708-731.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%	    Default options
%------------------------------------------------------
        options = struct('display',    'on', ...
                         'tol',        1e-4, ...
                         'maxFunEvals', 1000, ...
                         'dispIter',   100, ...                        
                         'iniStep',    0.1, ...
                         'multiStep',  1, ...
                         'maxFail',    10, ...
                         'decStepBy',  0.5, ...
                         'maxSuccess', 2, ...
                         'incStepBy',  2, ...
                         'adaptRate',  0.1 ...
                         );
        if nargin==1 && nargout <= 1 && isequal(funf,'defaults')
           xBest = options;
           return
        end;        
%------------------------------------------------------
%	    Parameters
%------------------------------------------------------
        len = length(x);
        x = reshape(x, len, 1, 1);
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
        dispIter    = options.dispIter;
        step        = options.iniStep;      
        multiStep   = options.multiStep;
        maxFail     = options.maxFail;
        decStepBy   = options.decStepBy;       
        maxSuccess  = options.maxSuccess;
        adaptRate   = options.adaptRate;
        incStepBy   = options.incStepBy;
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
           disp('# Random search procedure:');
           disp('      Iter       FunVal         Step');
           disp('--------------------------------------');
        end;        
%------------------------------------------------------
%	    Initial values
%------------------------------------------------------
        xBest = x;
        valBest = feval ( funf, x );
        func_evals = 1;
%------------------------------------------------------
        failCount = 0;
        successCount = 0;
        iter = 0;
        if verbosity > 1, showinfo(iter, valBest, step, '  '); end;
%------------------------------------------------------       
        valMult = zeros(multiStep, 1);
        dxMult  = zeros(len, multiStep);
        adaptDir = zeros(len, 1);
        h = step;
        hSuccess = ones(1,5)*step/2;
        ySuccess = (5:-1:1)*(1 + valBest);
%######################################################
%	    Main loop
%######################################################
        while step > tol
           iter = iter + 1;
           if func_evals > maxFunEvals, break; end;
           bShown = 0; % display flag for this step
          %----------------------------------------------
          %  Multiple attempts
          %----------------------------------------------
           for k=1:multiStep
                %---  Choose random direction
              dirX = adaptDir + rand(len,1) - 0.5;
              dirX = dirX / norm(dirX);
              dx = dirX * h;
              xNew = x + dx;
                %---  Check constraints
              if withProj, xNew = feval ( funproj, xNew ); end;
              if withConstr && any(feval(fung,xNew) > 0), 
                   valMult(k) = Inf;
              else valMult(k) = feval ( funf, xNew ); 
              end;
              dxMult(:,k) = xNew - x;
              func_evals = func_evals + 1;
           end;   
           [valMin,ind] = min(valMult);
            dx = dxMult(:,ind);           
          %----------------------------------------------
          %  Adaptation
          %----------------------------------------------
           if ~isinf(valMin), 
              if ~isinf(valMin)
                   adaptDir = (1-adaptRate)* adaptDir + adaptRate*dx*sign(valBest-valMin);
              else adaptDir = adaptDir*0; 
              end;
              if norm(adaptDir) > eps,
                 adaptDir = adaptRate * adaptDir / norm(adaptDir);
              end;
           end;
          %----------------------------------------------
          %  Choose the best attempt
          %----------------------------------------------
           if valMin < valBest
              valBest = valMin;
              x = x + dx;
              successCount = successCount + 1;
              hSuccess = [hSuccess(2:end) h];
             %----------------------------------------------
             %  Switch to local optimization
             %----------------------------------------------
              ySuccess = [ySuccess(2:end) valBest];              
              dyAll = max(ySuccess) - min(ySuccess);
              if dyAll < 1e-4*valBest  ||  dyAll < tol, break; end;
             %----------------------------------------------
             %  Increase step on several successful sampling
             %----------------------------------------------
              if successCount >= maxSuccess,
                 step = step * incStepBy;
                 if verbosity > 1, showinfo(iter, valBest, step, 'v^'); bShown = 1; end;
                 successCount = 0;
              elseif verbosity > 1, showinfo(iter, valBest, step, 'v '); bShown = 1; 
              end;
              failCount = 0;
              h = step;
           else
              successCount = 0;
              h = h * decStepBy;
              if h < tol,   
                failCount = failCount + 1;
                if failCount >= maxFail, 
                  step = step * decStepBy; 
                  if verbosity > 1, showinfo(iter, valBest, step, ' v'); bShown = 1; end;
                  failCount = 0;
                end;
                h = step;
              end;
           end; 
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
%------------------------------------------------------
%	    Local optimization
%------------------------------------------------------
        %[x,valBest,func_loc] = neldermead(funf, x, tol); 
        %func_evals = func_evals + func_loc;
        
%------- End of ARANDSEARCH.M --------- KYuP ----------

%####################################################
function showinfo ( iter, fval, step, marks )

        fprintf('%10d   %s%10.5g   %s%10.5f\n', iter, marks(1), fval, marks(2), step );

%------- End of SHOWINFO.M --------- KYuP ----------
