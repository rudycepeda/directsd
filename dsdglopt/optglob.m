function [xBest,zBest,coef,iter,xTrace] = optglob ( func, options, f_func )
%OPTGLOB Global optimization using information algorithm.
%
%	  [X,Y,COEF,ITER] = OPTGLOB ( F, OPTIONS )
%
%   Inputs:
%	  F - function to be minimized as a string
%         F(x) assumes that 0 <= x <= 1 
%     OPTIONS - structure with the following fields: 
%          .dim - number of unknown (default 1)
%          .tol - tolerance (default 1e-4)
%          .r   - optimization parameter (default 2) 
%          .ksi - optimization parameter (default 1)
%          .maxIter  - maximal number of iterations (default 100)
%          .maxLoop  - maximal number of loops (default 20)
%          .dispIter - step of dislaying results (defau;t 50)
%          .display  - flag for dislaying results ('on' or 'off', default 'on')
%          .bounds   - whether to check initial bounds (default 0)
%          .decLim   - decrease ratio for limiting bounds (default 2.5);
%          .incR     - increase ratio of parameter 'r' (default 1.2);
%
%   Outputs:
%	  X    - minimizing solution
%     Y    - minimal function value 
%     COEF - vector of unknown coefficients
%     ITER - number of iterations
%
%   See also OPTGLOBC.
	
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Global variables
%------------------------------------------------------
        global gloptData;
%------------------------------------------------------
%       Check syntax
%------------------------------------------------------
        tol = 1e-4;
        tolLim = 0.01;
        maxLoop = 20;
        display = 'off';
        bounds = 0;
        decLim = 2.5;
        incR = 1.2;
        if exist('options','var')
           if isfield(options,'tol'),     tol = options.tol; end;
           if isfield(options,'tolLim'),  tol = options.tolLim; end;
           if isfield(options,'display'), display = options.display; end;
           if isfield(options,'maxLoop'), maxLoop = options.maxLoop; end;
           if isfield(options,'bounds'),  bounds = options.bounds; end;
           if isfield(options,'decLim'),  decLim = options.decLim; end;
           if isfield(options,'incR'),    incR = options.incR; end;
        else
           options.display = 'on';
           options.maxIter = 100;
        end
%------------------------------------------------------
%       Initialization
%------------------------------------------------------
        ws = warning; warning off;
        coefSup = [];
        zSup = Inf;
        k = 0;     
        Lim  = gloptData.Lim;
        Lim0 = gloptData.Lim;
        Ncoef = size(Lim,1);
        w = zeros(Ncoef,1);
        if isequal(display,'on') 
          disp(' ')  
          disp('Global optimization:');  
        end;
%------------------------------------------------------
%       Main loop
%------------------------------------------------------
        while 1
          k = k + 1;
          if isequal(display,'on') 
            disp(' ')  
            disp(sprintf('Run %d', k));  
            disp('------------------------------');  
          end;
          [xBest,zBest,iter,xTrace] = infglob ( func, options ); 
         %---------------------------------------------- 
         %   Solution found
         %---------------------------------------------- 
          if ~isempty(xBest)
             [xx,coef] = feval(func, xBest);
             if zBest < zSup,
                zSup = zBest; 
                coefSup = coef;
             end;
          end;
         %---------------------------------------------- 
         %   Check loop counter
         %---------------------------------------------- 
          if k >= maxLoop, break; end;
         %---------------------------------------------- 
         %   Change limits
         %---------------------------------------------- 
          if ~isempty(xBest)
             if abs(Lim(:,2)-Lim(:,1)) < tolLim, break; end;
             try
               for i=1:Ncoef
                 w(i) = (Lim(i,2) - Lim(i,1)) / decLim;  
                 Lim(i,1) = coefSup(i) - w(i);
                 Lim(i,2) = coefSup(i) + w(i);
               end;
             catch keyboard; 
             end;
         %---------------------------------------------- 
         %   No solution found
         %---------------------------------------------- 
          else
             for i=1:Ncoef
               ci = (Lim(i,2) + Lim(i,1)) / 2;  
               w(i) = (Lim(i,2) - Lim(i,1)) / decLim;  
               Lim(i,1) = ci - w(i);
               Lim(i,2) = ci + w(i);
             end;
          end;   
         %---------------------------------------------- 
         %   Check bounds
         %---------------------------------------------- 
          if bounds
             for i=1:Ncoef
               wi = Lim(i,2) - Lim(i,1);  
               if Lim(i,1) < Lim0(i,1)
                  Lim(i,1) = Lim0(i,1); 
                  Lim(i,2) = Lim0(i,1) + 2*w(i); 
               end;   
               if Lim(i,2) > Lim0(i,2)
                  Lim(i,2) = Lim0(i,2); 
                  Lim(i,1) = Lim0(i,2) - 2*w(i); 
               end;   
             end;
          end;
          options.r = min(10, incR*options.r);
        end;
%------------------------------------------------------
%       Finalization
%------------------------------------------------------       
        zBest = zSup;
        coef = coefSup;
        if nargin > 2
           if isequal(display,'on')
              disp(' ') 
              disp(' Nonlinear local optimization ... ') 
              disp('-----------------------------------') 
           end;
           if bounds
              ops = optimset('MaxFunEvals', 1000*length(coef));
              [coef,zBest] = fmincon(f_func,coef,[],[],[],[],Lim0(:,1),Lim0(:,2),[],ops);
           else 
              [coef,zBest] = fminsearch(f_func, coef);
           end;
           xBest = [];            
        end;
        warning(ws);
        
%------- End of OPTGLOB.M --------- KYuP ----------
