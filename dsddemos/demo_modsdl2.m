%DEMO_MODSDL2 Quasioptimal reduced-order L2-controller.
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
        disp('   Reduced-order L2-optimal controller.')
        disp('#############################################################')
%------------------------------------------------------
%%       Initial data            
%------------------------------------------------------
        F  = tf ( 10, [2 1 0] );
        R  = tf ( 1, [1 0] );
        Q  = tf(1, [1 2 1]);
        T  = 0.2;
        if paused
          disp(' ')
          disp('# Plant:'); F
          disp('# Ideal operator:'); Q
          disp('# Sampling period:'); T
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       Form standard system
%------------------------------------------------------
        P11 = -minreal(Q*R);
        P12 = F;
        P21 = R;
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
        alpha = 0.1;
        beta = Inf;
        options.optPoles = rOpt;
        options.method = 'infglob';
          %options.method = 'simanneal';
          %options.method = 'randsearch';
        [KRaw,err,x] = modsdl2 ( sys, T, ordK, alpha, beta, options );
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
        options.paramGuess = x;
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
        sim ( 'sim_track', [], opt );
        warning on;
%------------------------------------------------------
%%       Plot curves
%------------------------------------------------------
        close all;
        figure;
        subplot(2, 1, 1);
        plot(y(:,1),y(:,2),y(:,1),y(:,3));
        title('Output transients in ideal and actual systems');
        xlabel('Time, t');
        ylabel('Output');
        legend('Actual response', 'Ideal respose');

        subplot(2, 1, 2);      
        plot(u(:,1),u(:,2),u(:,1),u(:,3));
        title('Control transients in ideal and actual systems');
        xlabel('Time, t');
        ylabel('Control');
        legend('Actual response', 'Ideal respose');               
    
%------- End of DEMO_MODSDL2.M --------- KYuP ----------    
