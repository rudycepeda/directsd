%DEMO_AIT01B AHinf-optimization for tracking system.
%

% References:
% [1]  K. Yu. Polyakov, Polynomial design of optimal sampled-data 
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
        disp('   AHinf-optimization for tracking system.')
        disp('#############################################################')
%------------------------------------------------------
%%       Initial data            
%------------------------------------------------------
        F  = tf ( 1, [1 1] );
        G  = 1;
        Fo = 1;
        Q  = 1;
        Qu = 0;
        R  = tf ( 1, [1 0] );
        T  = 0.2;
        Ve = 1;
        Vu = 0;
        if paused
          disp(' ')
          disp('# Plant:'); F
          disp('# Control weight:'); Vu
          disp('# Sampling period:'); T
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       Form standard system
%------------------------------------------------------
        P11 = [-minreal(Ve*Q*R); -minreal(Vu*Qu*R)];
        P12 = [minreal(Ve*F); Vu];
        P21 = minreal(Fo*R);
        P22 = -minreal(Fo*G*F);
        sys = [P11 P12;P21 P22];
        if paused
          clc
          disp('# Standard system:'); sys
          disp('Press any key ...')
          pause
          clc
        end;
%------------------------------------------------------
%%       L2-optimal controller
%------------------------------------------------------
        disp(' ')
        disp('=========================================')
        disp('           L2-optimization')
        disp('=========================================')
        disp(' ')
        [KL2,err] = sdl2 ( sys, T );
        KL2
        disp('Closed-loop poles');    
        disp(roots(charpol(sys,KL2)));
        fprintf('L2 cost function:   %g\n',err);    
        lamL2 = sdtrhinferr ( sys, KL2 );
        fprintf('AHinf-cost:         %g\n\n', lamL2);    
%------------------------------------------------------
%%       Hinf-optimal controller
%------------------------------------------------------
        disp(' ')
        disp('=========================================')
        disp('           AHinf-optimization')
        disp('=========================================')
        disp(' ')
        [K,lam] = sdtrhinf ( sys, T );
        K=K{1}
        disp('Closed-loop poles');    
        disp(roots(charpol(sys,K)));
        err = sdl2err ( sys, K );
        fprintf('L2 cost function:   %g\n',err);    
        fprintf('Optimal AHinf-cost: %g\n', lam);    
        lam1 = sdtrhinferr ( sys, K );
        fprintf('AHinf-cost:         %g\n', lam1);    
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
        [nFo,dFo] = tfdata ( tf(Fo), 'v' );
        [nF,dF]   = tfdata ( tf(F), 'v' );
        [nG,dG]   = tfdata ( tf(G), 'v' );
        [nQ,dQ]   = tfdata ( tf(Q), 'v' );
        [nQu,dQu] = tfdata ( tf(Qu), 'v' );
        Tmax = 1;
        tauR  = R.iodelay;
        tauF  = F.iodelay;
        tauG  = 0; %G.iodelay;
        tauQ  = 0; %Q.iodelay;
        tauQu = 0; %Qu.iodelay;
        
        opt = simset('SrcWorkspace', 'current');
        warning off;
        [nK,dK] = tfdata ( tf(KL2), 'v' );        
        sim ( 'sim_track', [], opt ); 
        yL2 = y; uL2 = u;
        [nK,dK] = tfdata ( tf(K), 'v' );        
        sim ( 'sim_track', [], opt );
        warning on;
%------------------------------------------------------
%%       Plot curves
%------------------------------------------------------
        close all;
        figure;
        subplot(2, 1, 1);
        plot(y(:,1),y(:,2),yL2(:,1),yL2(:,2),'--');
        title('Output transients');
        xlabel('Time, t');
        ylabel('Output');
        legend('AHinf-optimal', 'L2-optimal');

        subplot(2, 1, 2);      
        plot(u(:,1),u(:,2),uL2(:,1),uL2(:,2),'--');
        title('Control transients');
        xlabel('Time, t');
        ylabel('Control');
        legend('AHinf-optimal', 'L2-optimal');
        
%------- End of DEMO_AIT01B.M --------- KYuP ----------
        





