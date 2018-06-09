%DEMO_MODSDL2B Reduced-order L2-optimal controller with fixed static gain.
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
        disp('   Reduced-order L2-optimal controller with fixed')
        disp('   static gain of closed-loop system.');
        disp('#############################################################')
%------------------------------------------------------
%%       Initial data            
%------------------------------------------------------
        F = tf(1, [1 1 1]);
        R = tf(1, [1 0]);
        Q = tf(1, [5 1]);
        Qu = tf(1);
        rho = 0;
        T = 0.2;       
        if paused
          disp(' ')
          disp('# Plant:'); F
          disp('# Ideal operator:'); Q
          disp('# Control weight:'); rho
          disp('# Sampling period:'); T
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       Form standard system
%------------------------------------------------------
        sys = [    -Q*R     F
                -rho*Qu*R   rho
                     R     -F];          
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
        disp('      Full-order optimal controller')
        disp('=========================================')
        [KOpt,errOpt] = sdl2 ( sys, T )
        disp(sprintf('Integral error   %g', errOpt));
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
        options.optPoles = rOpt;
        options.method = 'infglob';
          %options.method = 'simanneal';
          %options.method = 'randsearch';
        [KRaw,err,rho] = modsdl2 ( sys, T, ordK, alpha, beta, options );
        disp(' ')
        disp('=========================================')
        disp('       Reduced-order modal controller')
        disp('=========================================')
        KRaw
        disp(sprintf('Integral error   %g', err));
        disp('Closed-loop poles:')
        r = roots(charpol(sys, KRaw))
        if paused
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%   Reduced-order controller (local optimization)
%------------------------------------------------------
        disp(' ')
        disp('# Local optimization ...');    
        disp(' ')
        options.paramGuess = rho;
        options.method = 'local';
        [K,err] = modsdl2 ( sys, T, ordK, alpha, beta, options );
        disp(' ')
        disp('=========================================')
        disp('    Optimized reduced-order controller')
        disp('=========================================')
        K
        disp(sprintf('Integral error   %g', err));
        disp('Closed-loop poles:')
        r = roots(charpol(sys, K))
        if paused
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       Simulation
%------------------------------------------------------                             
        disp(' ')
        disp('# Simulation ...');    
        disp(' ')
        Fo = 1;
        G  = 1;
        [nR,dR]   = tfdata ( tf(R), 'v' );
        [nFo,dFo] = tfdata ( tf(Fo), 'v' );
        [nF,dF]   = tfdata ( tf(F), 'v' );
        [nG,dG]   = tfdata ( tf(G), 'v' );
        [nQ,dQ]   = tfdata ( tf(Q), 'v' );
        [nQu,dQu] = tfdata ( tf(Qu), 'v' );
        [nK,dK]   = tfdata ( tf(KOpt), 'v' );        
        Tmax = 20;
        tauR  = R.iodelay;
        tauF  = F.iodelay;
        tauG  = 0; %G.iodelay;
        tauQ  = 0; %Q.iodelay;
        tauQu = 0; %Qu.iodelay;
        
        opt = simset('SrcWorkspace', 'current');
        warning off;
        sim ( 'sim_track', [], opt );
        yOpt = y; uOpt = u;
        [nK,dK] = tfdata ( tf(K), 'v' );        
        sim ( 'sim_track', [], opt );
        warning on;
%------------------------------------------------------
%%       Plot curves
%------------------------------------------------------
        close all;
        figure;
        subplot(2, 1, 1);
        plot(y(:,1),y(:,3),yOpt(:,1),yOpt(:,2),y(:,1),y(:,2));
        title('Output transients in ideal and actual systems');
        xlabel('Time, t');
        ylabel('Output');
        legend('Ideal respose', 'Optimal response', 'Suboptimal respose');

        subplot(2, 1, 2);      
        plot(u(:,1),u(:,3),uOpt(:,1),uOpt(:,2),u(:,1),u(:,2));
        title('Control transients in ideal and actual systems');
        xlabel('Time, t');
        ylabel('Control');
        legend('Ideal respose', 'Optimal response', 'Suboptimal respose');
    
%------- End of DEMO_MODSDL2B.M --------- KYuP ----------    
