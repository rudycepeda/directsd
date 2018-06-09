%DEMO_MODSDL2A Quasioptimal reduced-order L2-redesign.
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
        paused = 1;
        disp('#############################################################')
        disp(' DirectSD Toolbox Example:')
        disp('   L2-optimal redesign problem with fixed')
        disp('   static gain of closed-loop system.');
        disp('#############################################################')
%------------------------------------------------------
%%       Initial data            
%------------------------------------------------------
        F = tf(1, [1 -1 0]);
        Kc = tf([5 1],[1 3]);
        R = tf(1, [1 0]);   
        Q = feedback(F,Kc);
        T = 0.5;
        if paused
          disp(' ')
          disp('# Plant:'); F
          disp('# Continuos controller:'); Kc
          disp('# Sampling period:'); T
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       Form standard system
%------------------------------------------------------
        sys = [(F-Q)*R  F
                -F*R   -F];
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
        ordK = 2;
        alpha = 0.001;
        beta = Inf;
        options.optPoles = rOpt;
        options.method = 'infglob';
          %options.method = 'simanneal';
          %options.method = 'randsearch';
        options.maxIter = 100;
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
        G  = 1;
        Fo = 1;
        Qu = 0;
        [nR,dR]   = tfdata ( tf(R), 'v' );
        [nFo,dFo] = tfdata ( tf(Fo), 'v' );
        [nF,dF]   = tfdata ( tf(F), 'v' );
        [nG,dG]   = tfdata ( tf(G), 'v' );
        [nQ,dQ]   = tfdata ( tf(Q), 'v' );
        [nQu,dQu] = tfdata ( tf(Qu), 'v' );
        [nK,dK]   = tfdata ( tf(K), 'v' );        
        Tmax = 10;
        tauR  = R.iodelay;
        tauF  = F.iodelay;
        tauG  = 0; %G.iodelay;
        tauQ  = 0; %Q.iodelay;
        tauQu = 0; %Qu.iodelay;
        
        opt = simset('SrcWorkspace', 'current');
        warning off;
        sim ( 'sim_track2', [], opt );
        warning on;
%------------------------------------------------------
%%       Plot curves
%------------------------------------------------------
        close all;
        figure;
        %subplot(2, 1, 1);
        plot(y(:,1),y(:,3),'b--',y(:,1),y(:,3),'g-');
        title('Output transients in ideal and actual systems');
        xlabel('Time, t');
        ylabel('Output');
        legend('Reference', 'Step response');
        if 0
          subplot(2, 1, 2);      
          plot(u(:,1),u(:,3),'b--',u(:,1),u(:,2),'g-');
          title('Control transients in ideal and actual systems');
          xlabel('Time, t');
          ylabel('Control');
          legend('Reference', 'Step response');
        end;
    
%------- End of DEMO_MODSDL2A.M --------- KYuP ----------    
