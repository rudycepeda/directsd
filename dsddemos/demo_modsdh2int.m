%DEMO_MODSDH2INT Reduced-order H2-optimal controller with an integrator.
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
        disp('   Reduced-order H2-optimal controller with an integrator.')
        disp('#############################################################')
%------------------------------------------------------
%%   System model
%------------------------------------------------------
        F1 = tf(0.0694, [18.22 1]);
        F2 = tf ( 1, [1 0] );
        lam = 0.3;
        w0 = 0.3;
        sigma = 7.25;
        Fw = tf ( [2*lam*w0*sigma 0], [1 2*lam*w0 w0^2] );
        rho = 2;
        T  = 2;
        if paused
          disp(' ')
          disp('# Plant:'); F1, F2
          disp('# Disturbance spectrum:'); Fw
          disp('# Control weight:'); rho
          disp('# Sampling period:'); T
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       Form standard system
%------------------------------------------------------
        P11 = [minreal(F2*Fw); 0];
        P12 = [minreal(F2*F1); rho];
        P21 = -minreal(F2*Fw);
        P22 = -minreal(F2*F1);
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
        alpha = 0.02; 
        beta = 2;
        ordK = 2;
        options.r = 10;
        options.dK0 = [1 -1];
        options.method = 'infglob';        
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
%------------------------------------------------------
%%   Degree of stability and oscillation
%------------------------------------------------------
        [alphaR,betaR] = sector ( r, 'z' );
        alphaR = alphaR / T
        betaR
        if paused
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%   Reduced-order controller (local optimization)
%------------------------------------------------------
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
%------------------------------------------------------
%%   Degree of stability and oscillation
%------------------------------------------------------
        [alphaR,betaR] = sector ( r, 'z' );
        alphaR = alphaR / T
        betaR
    
%------- End of DEMO_MODSDH2INT.M --------- KYuP ----------    
