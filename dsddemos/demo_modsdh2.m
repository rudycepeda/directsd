%DEMO_MODSDH2 Quasioptimal reduced-order H2-controller. 
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%%	Preliminaries
%------------------------------------------------------
        clc
        clear all;
        close all;
        format short;
        format compact;
        paused = 0;
        disp('#############################################################')
        disp(' DirectSD Toolbox Example:')
        disp('   Reduced-order H2-optimal controller.')
        disp('#############################################################')
%------------------------------------------------------
%%   System model
%------------------------------------------------------
        F = tf ( 0.051, [25 1 0] );
        rho = sqrt(0.1);
        Sw = tf ( 0.0757, [1 0 2.489 0 1.848] );
        Fw = sfactor ( Sw );
        T  = 1;
        if paused
          disp(' ')
          disp('# Plant:'); F
          disp('# Disturbance spectrum:'); Fw
          disp('# Control weight:'); rho
          disp('# Sampling period:'); T
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       Form standard system
%------------------------------------------------------
        P11 = [minreal(F*Fw); 0];
        P12 = [F; rho];
        P21 = -minreal(F*Fw);
        P22 = -F;
        sys = [P11 P12;P21 P22];
        if paused
          clc
          disp('# Standard system:'); sys
          disp('Press any key ...')
          pause
          clc
        end;
%------------------------------------------------------
%%   Full-order optimal controller 
%------------------------------------------------------   
        disp(' ')
        disp('=========================================')
        disp('     Full-order optimal controller')
        disp('=========================================')
        [KOpt,errOpt] = sdh2 ( sys, T )
        disp(sprintf('H2-norm   %g', errOpt));
        disp('Closed-loop poles:')
        rOpt = roots(charpol(sys, KOpt))
        if paused
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%   Reduced-order controller
%------------------------------------------------------
        disp(' ')
        disp('# Global optimization ...');    
        disp(' ')
        ordK = 1;
        alpha = 0;
        beta = Inf;
          %options.method = 'infglob';        
          %options.method = 'simanneal';        
          %options.method = 'randsearch';        
        options.optPoles = rOpt;
        [KRaw,err,x] = modsdh2 ( sys, T, ordK, alpha, beta, options );
        disp(' ')
        disp('=========================================')
        disp('       Reduced-order modal controller')
        disp('=========================================')
        KRaw
        disp(sprintf('H2-norm   %g', err));
        disp('Closed-loop poles:')
        r = roots(charpol(sys, KRaw))
        if paused
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%   Reduced-order controller (further local optimization)
%------------------------------------------------------
        disp(' ')
        disp('# Local optimization ...');    
        disp(' ')
        options.paramGuess = x;
        options.method = 'local';
        [K,err] = modsdh2 ( sys, T, ordK, alpha, beta, options );
        disp(' ')
        disp('=========================================')
        disp('    Optimized reduced-order controller')
        disp('=========================================')
        K
        disp(sprintf('Integral error   %g', err));
        disp('Closed-loop poles:')
        r = roots(charpol(sys, K))
    
%------- End of DEMO_MODSDH2.M --------- KYuP ----------    
