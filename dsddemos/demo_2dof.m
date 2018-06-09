%DEMO_2DOF Optimal 2-DOF tracking system.
%

% References:
% [1]  K. Yu. Polyakov, Design of optimal 2-DOF sampled-data systems,
%      Automation & Remote Control, vol. 62, 2001; no. 6, pp. 934-942. 
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
        disp('   Optimal digital 2-DOF controller.')
        disp('#############################################################')
%------------------------------------------------------
%%       Initial data            
%------------------------------------------------------
        F = tf ( 1, [1 -1] );
        R = tf ( 1, [1 0] );
        Q  = tf ( 1, [1 2] );
        T  = 0.5;
        k  = 0;
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
        P11 = minreal(Q*R);
        P12 = -F;
        P21 = [R; k*R];
        P22 = [0; -F];
        sys = [P11 P12;P21 P22];
        sys1DOF = sys; 
        sys1DOF(end-1,:) = [];
        sys1DOF(end,1) = P21(1);
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
        disp('=========================================')
        disp('          1-DOF L2-optimization')
        disp('=========================================')
        disp(' ')
        [K,errOpt] = sdl2 ( sys1DOF, T, [], 'pol' );
        K
        fprintf('L2-cost:            %g\n',errOpt);    
        err = sdl2err ( sys1DOF, K, [], 'pol' );
        fprintf('Direct computation  %g\n',err);            
%------------------------------------------------------
%%       2-DOF L2-optimization
%------------------------------------------------------                             
        disp(' ')
        disp('=========================================')
        disp('           2-DOF optimization')
        disp('=========================================')
        disp(' ')
        [KR,errOpt] = sd2dof ( sys, K, [], 'pol' );
        disp('Optimal reference controller');    
        KR
        fprintf('L2-cost:             %g\n',errOpt);    
        err = sd2doferr ( sys, K, KR );
        fprintf('Direct calculation:  %g\n',err);    
        if paused
          disp('Press any key ...')
          pause
        end;
%---------------------------------------------------
%%       Simulate 1-DOF system
%---------------------------------------------------
        disp(' ')
        disp('# Simulation 1-DOF system ...');    
        disp(' ')
        G = 1;
        Qu = - 0.5;
        [nR,dR]   = tfdata ( tf(R), 'v' );
        [nFo,dFo] = tfdata ( tf(1), 'v' );
        [nF,dF]   = tfdata ( tf(F), 'v' );
        [nG,dG]   = tfdata ( tf(G), 'v' );
        [nQ,dQ]   = tfdata ( tf(Q), 'v' );
        [nQu,dQu] = tfdata ( tf(Qu), 'v' );
        [nK,dK] = tfdata ( tf(K), 'v' );        
        Tmax = 5;
        tauR  = R.iodelay;
        tauF  = F.iodelay;
        tauG  = 0; %G.iodelay;
        tauQ  = 0; %Q.iodelay;
        tauQu = 0; %Qu.iodelay;

        opt = simset('SrcWorkspace', 'current');
        warning off;
        sim ( 'sim_track', [], opt );
        warning on;        
        y1dof = y; u1dof = u;       
%---------------------------------------------------
%%       Simulate  2-DOF system
%---------------------------------------------------
        disp(' ')
        disp('# Simulation 2-DOF system ...');    
        disp(' ')
        [Kx,KRx,KC] = split2dof ( K, KR );
        [nK,dK]   = tfdata(Kx, 'v');
        [nKR,dKR] = tfdata(KRx, 'v');
        [nKC,dKC] = tfdata(KC, 'v');        
        warning off;
        sim ( 'sim_2dof', [], opt );
        warning on;      
        disp('OK')
%------------------------------------------------------
%%       Plot curves
%------------------------------------------------------
        close all;
        figure;
        subplot(2, 1, 1);
        plot(y(:,1),y(:,3),'r',y1dof(:,1),y1dof(:,2),'m',y(:,1),y(:,2),'b');
        title('Output transients in ideal and actual systems');
        xlabel('Time, t');
        ylabel('Output');
        legend('Reference', '1-DOF', '2-DOF');

        subplot(2, 1, 2);      
        plot(u1dof(:,1),u1dof(:,2),'m',u(:,1),u(:,2),'b');
        title('Control transients in ideal and actual systems');
        xlabel('Time, t');
        ylabel('Control');
        legend('1-DOF', '2-DOF');
        
%------- End of DEMO_2DOF.M --------- KYuP ----------
        





