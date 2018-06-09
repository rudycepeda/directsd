%DEMO_AT96 Optimal ship course stabilization.
%

% References:
% [1]  Y.N. Rosenwasser, K.Y. Polyakov, and B.P. Lampe,
%      Entwurf optimaler digitaler Kursregler mit Hilfe von
%      Parametrischen Uebertragungsfunktionen, Automatisierungstechnik,
%      1996, vol. 44, no. 10, pp. 487-495.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%%      Preliminaries
%------------------------------------------------------
        clc
        clear all;
        close all;
        format short;
        format compact;
        paused = 1;
        disp('#############################################################')
        disp(' DirectSD Toolbox Example:')
        disp('   Stochastic optimization for "Kazbek" type tanker.')
        disp('#############################################################')
%------------------------------------------------------
%%       Initial data            
%------------------------------------------------------
        F = tf ( 0.051, [25 1 0] );
        rho = sqrt(0.1);
        Sw = tf ( 0.0757, [1 0 2.489 0 1.848] );
        Fw = sfactor ( Sw );
        T  = 1;
        if paused
          disp(' ')
          disp('# Plant:'); F
          disp('# Wave spectrum:'); Fw
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
%%       H2-optimization
%------------------------------------------------------                             
        disp(' ')
        disp('=========================================')
        disp('          H2-optimal controller')
        disp('=========================================')
        [K,errOpt] = sdh2 ( sys, T );
        K
        fprintf('Optimal H2-cost:    %g\n',errOpt);    
        err = sdh2norm ( sys, K );
        fprintf('Direct calculation: %g\n',err);    
        if paused
          disp('Press any key ...')
          pause
          clc
        end;
%------------------------------------------------------
%%       Variate rho
%------------------------------------------------------                             
        rr = [0.0001 0.0003 0.001:0.001:0.005 0.007 0.01:0.01:0.1 0.2:0.1:1 2:1:10];
        sigma_psi = rr*0;
        sigma_u = rr*0;
        sys_psi   = [1 0 0; 0 0 1] * sys;
        sys(2,2)  = 1;
        sys_u = [0 1 0;0 0 1] * sys;
        disp(' ')
        disp('# Computing curve sigma_u(sigma_psi) ...');    
        disp(' ')
        for i=1:length(rr) 
           r = rr(i); disp(sprintf('r = %g',r));
           P12(2,1) = sqrt(r);
           sys = [P11 P12;P21 P22];           
           K = sdh2 ( sys, T );
           sigma_psi(i) = sdh2norm ( sys_psi, K );
           sigma_u(i) = sdh2norm ( sys_u, K );
        end;
        sigma_psi   = sigma_psi * 180 / pi;
        sigma_u = sigma_u * 180 / pi;
%------------------------------------------------------
%%       Plot curve
%------------------------------------------------------                             
        disp(' ')
        disp('# Plotting curve ...');    
        disp(' ')
        plot(sigma_psi, sigma_u);
        set(gca, 'yLim', [0 10], 'FontSize', 12);
        v = version;
        if v(1) > '6'
           h = xlabel('$\sigma_\psi, ^o$', 'Interpreter', 'LaTeX');
        else
           h = xlabel('\sigma_\psi, ^o');
        end;
        pos = get(h, 'Position');
        set ( h, 'FontSize', 14, 'Position', [pos(1) pos(2)+0.5]);
        if v(1) > '6'
           h = ylabel('$\sigma_u, ^o$', 'Rotation', 0, 'Interpreter', 'LaTeX' );
        else
           h = ylabel('\sigma_u, ^o', 'Rotation', 0 );
        end;
        pos = get(h, 'Position');
        set ( h, 'FontSize', 14, 'Position', [pos(1) 9.8]);
        disp('OK');    
                
%------- End of DEMO_AT96.M --------- KYuP ----------
        





