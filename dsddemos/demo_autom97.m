%DEMO_AUTOM97 H2-optimal control for delayed plants.
%

% References:
% [1]  Y.N. Rosenwasser, K.Y. Polyakov, and B.P. Lampe,
%      Frequency-domain method for H2-optimization of time-delayed 
%      sampled-data systems, Automatica, 1997, vol. 32, no. 7, 
%      pp. 1387-1392
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
        disp('   Stochastic optimization with time-delay.')
        disp('#############################################################')
%------------------------------------------------------
%%       Initial data            
%------------------------------------------------------
        F1 = tf ( 1 );
        F1.iodelay = 0.093;
        F = tf ( [1 0.5], [1 -1 0] );
        G  = 1;
        Ve = 1;
        Vu = 0;
        Fw = tf( 2, [1 2] );
        T  = 0.2;
        if paused
          disp(' ')
          disp('# Plant:'); F
          disp('# Disturbance spectrum:'); Fw
          disp('# Control delay:'); F1.iodelay
          disp('# Sampling period:'); T
          disp('# Control weight:'); Vu
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       Form standard system
%------------------------------------------------------
        P11 = [minreal(Ve*F*Fw); 0];
        P12 = [minreal(Ve*F*F1); Vu];
        P21 = -minreal(G*F*Fw);
        P22 = -minreal(G*F*F1);
        sys = [P11 P12;P21 P22];
        %sys = [ F*Fw  F*F1     
        %       -F*Fw -F*F1]; 
        if paused
          clc
          disp('# Standard system:'); sys
          disp('Press any key ...')
          pause
          clc
        end;
%------------------------------------------------------
%%       Sampled-data H2-optimization
%------------------------------------------------------                             
        disp(' ')
        disp('=========================================')
        disp('       Sampled-data H2-optimization')
        disp('=========================================')
        disp(' ')
        [K,errOpt] = sdh2 ( sys, T );
        disp('H2-optimal controller');    
        K
        fprintf('Optimal H2-cost:    %g\n',errOpt);    
        err = sdh2norm ( sys, K, [], [], 'pol' );
        fprintf('Average variance:   %g\n',err);    
        err0 = sdh2norm ( sys, K, 0, [], 'pol' );
        fprintf('Variance at t=0:    %g\n',err0);    
        disp('Computing variance between sampling instants ...');
        t = (0:0.025:1)*T;
        err_t = sdh2norm ( sys, K, t, [], 'pol' );
        disp('OK')
        if paused
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       Discrete H2-optimization for t=0
%------------------------------------------------------                             
        disp(' ')
        disp('=========================================')
        disp('         H2-optimization for t=0')
        disp('=========================================')
        disp(' ')
        [K0,errOpt] = sdh2 ( sys, T, 0 );
        disp('DH2-optimal controller');    
        K0
        fprintf('Variance at t=0:    %g\n',errOpt);    
        err0 = sdh2norm ( sys, K0, 0, [], 'pol' );
        fprintf('Direct computation  %g\n',err0);    
        err = sdh2norm ( sys, K0,[], [], 'pol' );
        fprintf('Average variance:   %g\n',err);    
        disp('Computing variance between sampling instants ...');
        err_t0 = sdh2norm ( sys, K0, t, [], 'pol' );
        disp('OK')
        if paused
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       Plot variance
%------------------------------------------------------                             
        disp(' ')
        disp('# Plotting curve ...');    
        disp(' ')
        plot ( t, err_t, t, err_t0, '--' );
        xlabel('Time'); ylabel('Variance');
        title('Variation of variance between sampling instants');
        legend('CMV-optimal', 'MV-optimal');
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
        R = tf ( 1, [1 0] );
        [nR,dR] = tfdata ( tf(R), 'v' );
        [nF,dF] = tfdata ( tf(F), 'v' );
        [nK,dK] = tfdata ( tf(K0), 'v' );        
        Tmax = 5;
        tauF  = F1.iodelay;       
        opt = simset('SrcWorkspace', 'current');
        warning off;
        sim ( 'sim_autom97', [], opt );
        y0 = y; u0 = u;
        [nK,dK] = tfdata ( tf(K), 'v' );        
        sim ( 'sim_autom97', [], opt );
        warning on;       
%------------------------------------------------------
%%       Generate arrays for the sampling instants
%------------------------------------------------------
        yT = (0:T:Tmax)';
        uT = yT;
        for i = 1:size(yT,1)
          ti = yT(i,1);  
          [xx,ind] = sort(abs(y0(:,1)-ti));  
          yT(i,2) = y0(ind(1),2);
          [xx,ind] = sort(abs(y(:,1)-ti));  
          yT(i,3) = y(ind(1),2);         
          [xx,ind] = sort(abs(u0(:,1)-ti));  
          uT(i,2) = u0(ind(1),2);         
          [xx,ind] = sort(abs(u(:,1)-ti));  
          uT(i,3) = u(ind(1),2);         
        end
%------------------------------------------------------
%%       Plot curves
%------------------------------------------------------
        disp(' ')
        disp('# Plotting curves ...');    
        disp(' ')
        figure(2);
        subplot(2, 1, 1);
        plot(y(:,1),y(:,2),y0(:,1),y0(:,2),'--',yT(:,1),yT(:,2),'o',yT(:,1),yT(:,3),'s');
        set(gca,'FontName','Arial','FontSize',18);
        xlabel('Time, t');
        ylabel('Output');
        legend('CMV-control', 'MV-control');

        subplot(2, 1, 2);      
        plot(u(:,1),u(:,2),u0(:,1),u0(:,2),'--');
        set(gca,'FontName','Arial','FontSize',18);
        xlabel('Time, t');
        ylabel('Control');
        legend('CMV-control', 'MV-control');
        
%------- End of DEMO_AUTOM97.M --------- KYuP ----------
        





