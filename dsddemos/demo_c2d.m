%DEMO_C2D Optimal digital redesign of a continuous-time controller.
%

% References:
% [1]  Rattan K., Compensating for computational delay in digital equivalent 
%      of continuous control systems // IEEE Trans. Automat. Contr., vol. AC-34,
%      pp. 895-899, 1989.
% [2]  E. Rosenwasser, K. Polyakov, and B. Lampe, Application of Laplace
%      transformation for digital redesign of continuous control systems 
%      // IEEE Trans. Automat. Contr., vol. AC-44, no. 4, pp. 883-886, 1999.
% [3]  K. Yu. Polyakov, Polynomial design of optimal sampled-data 
%      tracking systems. I & II. Automation & Remote Control,
%      2001; vol. 62, no. 2, pp. 304-316; no. 3, pp. 430-442. 
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
        paused = 0;
        disp('#############################################################')
        disp(' DirectSD Toolbox Example:')
        disp('   Optimal digital redesign (Rattan''s example).')
        disp('#############################################################')
%------------------------------------------------------
%%       Initial data            
%------------------------------------------------------
        F  = tf ( 10, [1 1 0] );
        Kc = tf ( [0.416 1], [0.139 1] );
        Q = minreal(feedback(F*Kc,1));
        Qu = minreal(feedback(Kc,F));
        F.iodelay = 0.01;
%        Q  = tf ( [29.9281 71.9424], ...
%                  [1.0000 8.1942 37.1223 71.9424] );
        R  = tf ( 1, [1 0] );
        T  = 0.04;
        if paused
          disp(' ')
          disp('# Plant:'); F
          disp('# Continuous-time controller:'); Kc
          disp('# Sampling period:'); T
          disp('Press any key ...')
          pause
        end;       
%------------------------------------------------------
%%       Form standard system
%------------------------------------------------------
        P11 = minreal(Q*R);
        P12 = -F;
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
%%       L2-optimization
%------------------------------------------------------                             
        disp(' ')
        disp('=========================================')
        disp('           L2-optimization')
        disp('=========================================')
        disp(' ')
        [K,errOpt] = sdl2 ( sys, T );
        disp('Optimal controller');    
        K
        disp('Closed-loop poles');    
        disp(roots(charpol(sys,K)));
        fprintf('Optimal L2 cost:    %g\n',errOpt);    
        err = sdl2err ( sys, K );
        fprintf('Direct calculation: %g\n',err);    
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
        [nR,dR]   = tfdata ( tf(R), 'v' );
        [nFo,dFo] = tfdata ( tf(1), 'v' );
        [nF,dF]   = tfdata ( tf(F), 'v' );
        [nG,dG]   = tfdata ( tf(1), 'v' );
        [nQ,dQ]   = tfdata ( tf(Q), 'v' );
        [nQu,dQu] = tfdata ( tf(Qu), 'v' );
        [nK,dK] = tfdata ( tf(K), 'v' );        
        Tmax = 2;
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
        plot(y(:,1),y(:,3),y(:,1),y(:,2));
        title('Output transients in ideal and actual systems');
        xlabel('Time, t');
        ylabel('Output');
        legend('Original', 'Redesigned');

        subplot(2, 1, 2);      
        plot(u(:,1),u(:,3),u(:,1),u(:,2));
        title('Control transients in ideal and actual systems');
        xlabel('Time, t');
        ylabel('Control');
        legend('Original', 'Redesigned');
        
%------- End of DEMO_C2D.M --------- KYuP ----------
        





