%DEMO_CF2 Examples 6.6.1, 8.4.2, 12.1.1 from the book by T. Chen and B. Francis.

% References:
% [1]  T. Chen & B.A. Francis, Optimal Sampled-Data Control 
%      Systems, Springer-Verlag, Berlin, 1995. 
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
        disp('   Example 6.6.1, 8.4.2, 12.1.1 from the book "Optimal')
        disp('   Sampled-Data Control Systems" by T. Chen and B. Francis.')
        disp('#############################################################')
%------------------------------------------------------
%%       Initial data            
%------------------------------------------------------
        F  = tf ( 1, conv([10 1], [25 1]) );
        Q  = 1;
        R  = tf ( 1, [1 0] );
        T = 1;
        if paused
          disp(' ')
          disp('# Plant:'); F
          disp('# Sampling period:'); T
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       Form standard system
%------------------------------------------------------
        P11 = minreal(-Q*R);
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
%%       L2-optimization
%------------------------------------------------------                             
        disp(' ')
        disp('=========================================')
        disp('            L2-optimization')
        disp('=========================================')
        [K,errOpt] = sdl2 ( sys, T );
        K
        fprintf('L2-cost:            %g\n',errOpt);    
        err = sdl2err ( sys, K );
        fprintf('Direct computation  %g\n',err);                       
%------------------------------------------------------
%%       Simulation
%------------------------------------------------------                                     
        disp(' ')
        disp('# Simulation ...');    
        disp(' ')
        Fo = 1; G = 1; Qu = 1;
        [nR,dR]   = tfdata ( tf(R), 'v' );
        [nFo,dFo] = tfdata ( tf(Fo), 'v' );
        [nF,dF]   = tfdata ( tf(F), 'v' );
        [nG,dG]   = tfdata ( tf(G), 'v' );
        [nQ,dQ]   = tfdata ( tf(Q), 'v' );
        [nQu,dQu] = tfdata ( tf(Qu), 'v' );
        Tmax = 10;
        tauR  = R.iodelay;
        tauF  = F.iodelay;
        tauG  = 0; %G.iodelay;
        tauQ  = 0; %Q.iodelay;
        tauQu = 0; %Qu.iodelay;

        [nK,dK] = tfdata ( K, 'v' );
        opt = simset('SrcWorkspace', 'current');
        warning off;
        sim ( 'sim_track', [], opt );
        warning on;                        
%------------------------------------------------------
%%       Plot curves
%------------------------------------------------------
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
        
%------- End of DEMO_CF2.M --------- KYuP ----------           
