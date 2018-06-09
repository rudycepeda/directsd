function [K,err,coef] = modalopt ( sys, T, ordK, alpha, beta, opt )
%MODALOPT Modal reduced-order optimal controller design.
%
%     [K,ERR,PAR] = MODALOPT ( SYS, T, N, ALPHA, BETA, OPTIONS )
%
%   Inputs:
%     SYS - an LTI system
%     T   - sampling period
%     N   - controller order (default = plant order - 1)
%     ALPHA - degree of stability (default 0)
%     BETA  - oscillation degree (default Inf)
%     OPTIONS - options, structure with the following possible fields
%            .go_func  - target function for global optimization
%            .f_func   - target function for direct optimization
%            .dcgain   - flag (1/0) indicating that static gain 
%                        should be fixed
%            .ssmodel  - command line for generating state-space model
%            .ABE_func - function for calculating coefficients in
%                        functional
%            .method   - optimization method: 'infglob' (default),
%                       'randsearch', 'simanneal' or 'local'
%            for other fields, see help on 'infglob', 'randsearch',
%            'simanneal' or 'fmincon', respectively  
%
%   Outputs:
%     K   - modal H2-optimal discrete controller of order N
%     ERR - cost function for the optimal system
%     PAR - vector of parameters in [0,1]
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%%       Check data
%------------------------------------------------------
        if ~exist('alpha','var'), alpha = 0; end; 
        alpha = max(1e-3, alpha)*T;
        if ~exist('beta','var'), beta = Inf; end; 
        if ~exist('opt','var'), opt = struct; end; 
%------------------------------------------------------
%%       Loop plant model
%------------------------------------------------------
        sys = zpk(sys);
        P22 = minreal(sys(end,end));
        p22 = allpoles(P22);
        Pz  = - dtfm ( P22, T, 0 );
        Pz  = setpoles(Pz, exp(p22*T));
%------------------------------------------------------
%%      Global data
%------------------------------------------------------
        global gloptData;
        gloptData.sys = sys;
        gloptData.plant = Pz; 
        gloptData.dK0 = 1;
        gloptData.alpha = alpha;
        gloptData.beta = beta;        
        gloptData = updateopt ( gloptData, opt );
%------------------------------------------------------
%%       Auxiliary data for maintaining static gain
%------------------------------------------------------
        if isfield(opt,'dcgain')  &&  opt.dcgain 
             [gloptData.gPsi,gloptData.gK,...
              gloptData.a0,gloptData.b0] = psigain ( sys, T );     
        else gloptData.gPsi = NaN; 
        end;
%------------------------------------------------------
%%       Check controller order
%------------------------------------------------------
        minOrdK = deg(Pz) - 1 + deg(gloptData.dK0);
        if ~isnan(gloptData.gPsi), minOrdK = minOrdK + 1; end;
        if ordK < minOrdK,
          disp('Warning: Controller order is too low for such a plant');  
          disp(sprintf('         Order %d is taken.',minOrdK));
          ordK = minOrdK;
        end;
        gloptData.ordK = ordK;
        gloptData.dim = deg(Pz) + ordK;
        gloptData.Lim = kron([0 1],ones(gloptData.dim,1));
%------------------------------------------------------
%%      Maximal possible precision
%------------------------------------------------------
        gloptData.bits = floor(52 / gloptData.dim);
%------------------------------------------------------
%%      Correct number of possible complex poles
%------------------------------------------------------
        maxPairs = floor(gloptData.dim/2);
        if ~isfield(gloptData,'nPairs')
          gloptData.nPairs = maxPairs;
        end;
        gloptData.nPairs = max(0, min(gloptData.nPairs, maxPairs));
%------------------------------------------------------
%%      Initial guess: poles of optimal system when given
%------------------------------------------------------
        if isfield(opt,'optPoles')
          poles = admproj(opt.optPoles, alpha, beta, 'z');  
          p = guesspoles ( poles, gloptData.dim );
          DeltaZ = poln(p, 'rz');
          paramGuess = cp2par ( DeltaZ, alpha, beta, gloptData.nPairs );
        else
          paramGuess = 0.5*ones(gloptData.dim,1);  
        end;        
        if isfield(opt,'paramGuess')
          opt.paramGuess = opt.paramGuess(:);  
          if length(opt.paramGuess) ~= gloptData.dim 
            warning('Incorrect size of options.paramGuess. Ingnored...');  
          else
            paramGuess = opt.paramGuess(:);          
          end;
        end;
%------------------------------------------------------
%%       Possible use of state-space analysis
%------------------------------------------------------
        gloptData.ss = (norm(totaldelay(sys))==0);
        if gloptData.ss
          eval(opt.ssmodel);             
        end;
%------------------------------------------------------
%%       Find A, B and E
%------------------------------------------------------
        if gloptData.ss  ||  ordK > minOrdK
          sys(end,:) = - sys(end,:);
          if gloptData.ss          
               [gloptData.A,gloptData.B] = feval(opt.ABE_func, sys, T );
          else [gloptData.A,gloptData.B,gloptData.E] = feval(opt.ABE_func, sys, T );
          end;
          sys(end,:) = - sys(end,:);
        end;
%------------------------------------------------------
%%   Optimization options
%------------------------------------------------------
        if isfield(opt,'method')
             options.method = opt.method; 
        else options.method = 'infglob'; 
        end;
        switch options.method
          case 'infglob',
            options.dim = gloptData.dim;
            options.r = 2;
            options.tol = 1e-4;
            options.maxIter  = 150;
            options.display  = 'on';
            options.dispIter = 10;
            options.bounds   = 1;
            options.maxLoop  = 1;
            options.decLim   = 0.5;
            options.incR     = 1.1;
          case 'randsearch',
            options.display = 'on';
            options.tol = 1e-4;
            options.maxFunEvals = 1000;
            options.dispIter = 100;
            options.iniStep = 1;
            options.multiStep = 1;
            options.decStepBy = 0.5;
            options.maxSuccess = 2;
            options.incStepBy = 2;
            options.adaptRate = 0.1;
          case 'simanneal',
            options.display = 'on';
            options.tol = 1e-4;
            options.maxFunEvals = 1000;
            options.dispIter = 100;
            options.iniStep = 1;
            options.multiStep = 20;
            options.decStepBy = 0.5;
            options.adaptRate = 0.1;
            options.decFailRate = 1.1;
            options.startTemp = 100;
            options.tempDecRate = [];
          case 'local',
            options.display = 'on';
          otherwise, error('Unknown optimization method ''%s''.',options.method);    
        end;    
        options = updateopt ( options, opt );
        if isfield(opt,'optPoles')
          options.guess = coord2hilb (paramGuess, gloptData.bits);
        end;
%------------------------------------------------------
%%   Perform optimization
%------------------------------------------------------
        switch options.method
          case 'infglob',
             [x,err,coef,iter] = optglob ( opt.go_func, options ); 
             [err,coef,K] = feval( opt.go_func, x );
          case 'randsearch',
             [coef,err,iter] = arandsearch ( opt.f_func, ...
                               paramGuess, options, [], 'uniproj' ); 
             [err,K] = feval( opt.f_func, coef );
          case 'simanneal',
             [coef,err,iter] = simanneal ( opt.f_func, ... 
                               paramGuess, options, [], 'uniproj' ); 
             [err,K] = feval( opt.f_func, coef );
          case 'local',             
             if isequal(options.display,'on')
                disp(' ') 
                disp(' Nonlinear local optimization ... ') 
                disp('-----------------------------------') 
             end;
             ops = optimset('MaxFunEvals', 1000*length(paramGuess));
             options = updateopt ( ops, options );
             [coef,zBest] = fmincon(opt.f_func, paramGuess, [], [], [], [], ...
                            gloptData.Lim(:,1), gloptData.Lim(:,2),[],options);
             [err,K] = feval( opt.f_func, coef );
        end;

%------- End of MODALOPT.M --------- KYuP ----------           
