function [x,valBest,func_evals,xTrace] = crandsearch ( funf, x0, opt, fung, funproj )
%CRANDSEARCH Controlled random search (CRS) optimization.
%          
%    [X,VAL,N] = CRANDSEARCH ( FUNF, X0, OPTIONS, FUNG, FUNPROJ )
%                                       
%  Inputs:
%    FUNF    - name of a fitness function
%    X0      - initial guess
%    OPTIONS - options, structure with the following possible fields
%           .display     - verbosity mode: 'on', 'off' or 'final'
%           .tol         - parameter tolerance (1e-4)
%           .maxFunEvals - maximal number of function call (default 1000)
%           .dispIter    - step of dislaying results (default 11)
%           .iniStep     - length of initial step in n-dimensional space (0.1)
%           .multiStep   - number of trial points in an iteration (20)
%           .maxFail     - number of unsuccessful iterations before step decrease (10)
%           .decStepBy   - step size multiplier if maxFail steps were unsuccessful (0.75)
%           .popSize     - initial population size (default 10*len)
%           .addSize     - number of new points at each iteration (2)
%    FUNG    - constraints function (should be gi(x) < 0)
%    FUNPROJ - function that projects point onto admissible region
%
%  Outputs:
%    X   - solution found
%    VAL - minimal attainable value
%    N   - number of function evaluation 
%
%   See also RANDSEARCH, ARANDSEARCH.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%	    Default options
%------------------------------------------------------
        options = struct('display',    'on', ...
                         'tol',        1e-4, ...
                         'maxFunEvals', 1000, ...
                         'dispIter',   1, ...                        
                         'iniStep',    0.1, ...
                         'multiStep',  20, ...
                         'maxFail',    10, ...
                         'decStepBy',  0.75, ...
                         'popSize',    [], ...
                         'addSize',    20 ...
                         );
        if nargin==1 && nargout <= 1 && isequal(funf,'defaults')
           xBest = options;
           return
        end;        
%------------------------------------------------------
%	    Parameters
%------------------------------------------------------
        x0 = x0(:);
        len = length(x0);
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
        N           = options.popSize;       
        M           = options.addSize;       
        if isempty(N), N = 10*len; end;
        N = max(N, 2*len+1);
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
           disp('# Controlled random search procedure:');
           disp('      Iter       FunVal         Step');
           disp('--------------------------------------');
        end;        
%------------------------------------------------------
%	    Initialization
%------------------------------------------------------
        nAdd = N;
        indN = 1:N+1;
        x = x0; 
        y = feval(funf, x0);
        func_evals = 1;
        valBest = y;
%------------------------------------------------------
        failCount = 0;
        iter = 0;
        if verbosity > 1, showinfo(iter, valBest, step, '  '); end;
        xTrace = [x0' 1];
        ySuccess = (5:-1:1)*(1 + valBest);
%######################################################
%	    Main loop
%######################################################
        while 1
           iter = iter + 1;
           if func_evals > maxFunEvals, break; end;
           bShown = 0; % display flag for this step
          %----------------------------------------------
          %  Add N new points around the best point 
          %----------------------------------------------
           while 1 
             bFound = 0;
             valWorst = y(end);
             for i=1:nAdd
               j = 1 + floor(rand*length(y));  
               dirX = rand(len,1) - 0.5;
               xi = x(:,j) + (0.5 + rand)*step*dirX/norm(dirX);
               yi = feval(funf, xi);
               x = [x xi];
               y = [y yi];
               xTrace = [xTrace; xi' 1];
               if yi < valWorst
                    bFound = 1;
               else xTrace(end,end) = 0; 
               end;             
             end;    
             func_evals = func_evals + nAdd; 
            %----------------------------------------------
            %  Decrease step until improvement
            %----------------------------------------------
             if bFound, break;
             else
               step = step * decStepBy; 
               decStepBy = 0.75 * decStepBy;
               if verbosity > 1, showinfo(iter, valBest, step, ' v'); bShown = 1; end;
             end;
           end;
           nAdd = M; % for all iterations except for the first one
          %----------------------------------------------
          %  Sort by increasing function values
          %----------------------------------------------
           [y,ind] = sort(y);
           x = x(:,ind);
          %----------------------------------------------
          %  Truncate to N+1 best points
          %----------------------------------------------
            y = y(indN);
            x = x(:,indN);
            if y(end)-y(1) < tol && step < tol, break; end;
          %----------------------------------------------
          %  Multiple attempts
          %----------------------------------------------
           for k=1:multiStep
            %----------------------------------------------
            %  Reflection of a randomly chosen simplex
            %----------------------------------------------
             xs = [x(:,1) sampling(x,len)];  
             xbar = sum(xs,2) / (len + 1);
             xNew = 2*xbar - xs(:,end);
             yNew = feval(funf, xNew);
             xTrace = [xTrace; xNew' 1];
             if yNew < y(end)
               x(:,end) = xNew;
               y(end) = yNew;  
               failCount = 0;
             else
               failCount = failCount + 1;
               xTrace(end,end) = 0;
             end;             
            %----------------------------------------------
            %  Sort by values of target function
            %----------------------------------------------
             [y,ind] = sort(y);
             x = x(:,ind);
             valBest = y(1);
             if y(end)-y(1) < tol, break; end;
           end;   
           func_evals = func_evals + k; 
          %----------------------------------------------
          %  Termination condition
          %----------------------------------------------
           if y(end)-y(1) < tol && step < tol, break; end;
          %----------------------------------------------
          %  Shrink simplex if necessary
          %----------------------------------------------
           if failCount >= maxFail, 
             for j=2:N
               x(:,j) = x(:,1)+decStepBy*(x(:,j) - x(:,1));
               y(j) = feval(funf, x(:,j));
             end               
             func_evals = func_evals + N - 1; 
             valBest = min(y);
             failCount = 0;
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
        [valBest,i] = min(y);
        x = x(:,i);
        if verbosity > 0
           if func_evals > maxFunEvals,
                disp(sprintf('Solution not found in %d function call.', maxFunEvals)); 
           else disp('Optimization completed successfully.'); 
           end;
        end;
%------------------------------------------------------
%	    Local optimization
%------------------------------------------------------
        [x,valBest,func_loc] = neldermead(funf, x, tol); 
        func_evals = func_evals + func_loc;
        
%------- End of CRANDSEARCH.M --------- KYuP ----------

%####################################################
function x = sampling ( x, n )
%
% SAMPLING choose n of all columns of x
%          not touching the first one
%
%------------------------------------------------------
%       K. Yu. Polyakov         01 Jan 2005
%                               
%------------------------------------------------------
        N = size(x,2);
        ind = 2:N;
        for i=N-1:-1:n+1
          j = 1 + floor(rand*i);
          ind(j) = [];
        end;
        x = x(:,ind);
        return;   

%------- End of SAMPLING.M --------- KYuP ----------

%####################################################
function showinfo ( iter, fval, step, marks )

        fprintf('%10d   %s%10.5g   %s%10.5f\n', iter, marks(1), fval, marks(2), step );

%------- End of SHOWINFO.M --------- KYuP ----------
