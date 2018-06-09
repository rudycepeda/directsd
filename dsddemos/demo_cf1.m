%DEMO_CF1 Example 12.4.2 from the book by T. Chen and B. Francis.

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
        paused = 1;
        disp('#############################################################')
        disp(' DirectSD Toolbox Example:')
        disp('   Example 12.4.2 from the book by T. Chen and B. Francis.')
        disp('#############################################################')
%------------------------------------------------------
%%       Initial data            
%------------------------------------------------------
        pd = [1/12 1/2 1];
        Gh = tf ( 2, conv(pd,pd) );
        Gm = tf (1, [1 0]);
        Q  = 1;
        T  = 1;
        if paused
          disp(' ')
          disp('# Plant1:'); Gh
          disp('# Plant2:'); Gm
          disp('# Sampling period:'); T
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       Form standard system
%------------------------------------------------------
        P11 = minreal((Gm-Q)*Gh);
        P12 = Gm;
        P21 = -minreal(Gm*Gh);
        P22 = -Gm;
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
%%       H2-optimization
%------------------------------------------------------                             
        disp(' ')
        disp('=========================================')
        disp('            H2-optimization')
        disp('=========================================')
        [KH2,errOpt] = sdh2 ( sys, T );
        KH2
        fprintf('H2-cost:            %g\n',errOpt);    
        err = sdl2err ( sys, KH2 );
        fprintf('L2-cost             %g\n',err);            
%------------------------------------------------------
%%       Simulation
%------------------------------------------------------                             
        disp(' ')
        disp('# Simulation ...');    
        disp(' ')
        [nGh,dGh] = tfdata ( Gh, 'v' );
        [nGm,dGm] = tfdata ( Gm, 'v' );
        tau = 0.8;

        [nK,dK] = tfdata ( KH2, 'v' );
        opt = simset('SrcWorkspace', 'current');
        warning off;
        sim ( 'sim_cf95', [], opt );
        warning on;                
        yH2 = y; yAppH2 = yApp;

        [nK,dK] = tfdata ( K, 'v' );
        warning off;
        sim ( 'sim_cf95', [], opt );
        warning on;                      
%------------------------------------------------------
%%       Plot curves
%------------------------------------------------------
        figure;
        subplot(2, 1, 1);
        plot(y(:,1),y(:,2),'r',y(:,1),y(:,3),'b',yH2(:,1),yH2(:,3),'g');
        title('Output transients in ideal and actual systems');
        xlabel('Time, t');
        ylabel('Output');
        legend('Ideal response', 'L2-optimal respose', 'H2-optimal response');

        subplot(2, 1, 2);      
        plot(yApp(:,1),yApp(:,2),'r',yApp(:,1),yApp(:,3),'b',yAppH2(:,1),yAppH2(:,3),'g--');
        title('Control transients in ideal and actual systems');
        xlabel('Time, t');
        ylabel('Control');
        legend('Ideal response', 'L2-optimal respose', 'H2-optimal response');
        
%------- End of DEMO_CF1.M --------- KYuP ----------           
