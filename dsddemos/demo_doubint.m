%DEMO_DOUBINT Optimal control for double integrator.
%

% References:
% [1]  Polyakov K. Y., Rosenwasser E. N., and Lampe B., 
%      Optimal digital controllers for double integrator: 
%      Comparison of four methods // Proc. CCA and
%      CACSD, Glasgow, UK, CACSDREG 1026, September 2002.
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
        disp('   Optimal digital control of double integrator.')
        disp('#############################################################')
%------------------------------------------------------
%%       Initial data            
%------------------------------------------------------
        F  = tf ( 1, [1 0 0] );
        T  = 0.1;
        Ve = 1;
        Vu = 1;
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
        P11 = [-minreal(Ve*F); 0];
        P12 = [-minreal(Ve*F); Vu];
        P21 = -F;
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
%%       Sampled-data H2-optimization
%------------------------------------------------------                             
        disp(' ')
        disp('=========================================')
        disp('       Sampled-data H2-optimization')
        disp('=========================================')
        disp(' ')
        Kopt = sdh2(sys, T);
        disp('H2-optimal controller');    
        Kopt
        disp('Closed-loop poles');    
        disp(roots(charpol(sys,Kopt)));
        err = sdh2norm ( sys, Kopt );
        fprintf('Sampled-data H2 cost:    %g\n',err);    
        if paused
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       Contrinuous-time H2-optimization
%------------------------------------------------------                             
        disp(' ')
        disp('=========================================')
        disp('     Continuous-time H2-optimization')
        disp('=========================================')
        disp(' ')
        Kc = ch2 ( sys );
        disp('Optimal continuous-time controller');    
        Kc
        Kc2d = bilintr ( Kc, 'tustin' );
        disp('Discretized controller');    
        Kc2d        
        disp('Closed-loop poles');    
        disp(roots(charpol(sys,Kc2d)));
        err = sdh2norm ( sys, Kc2d );
        fprintf('Sampled-data H2 cost:    %g\n',err);    
        if paused
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       H2-optimization for discretized model
%------------------------------------------------------                             
        disp(' ')
        disp('=========================================')
        disp('  H2-optimization for discretized model')
        disp('=========================================')
        disp(' ')
        dsys = c2d ( sys, T );
        Kd = minreal(zpk(h2reg(dsys)));
        disp('Optimal discrete-time controller');    
        Kd
        disp('Closed-loop poles');    
        disp(roots(charpol(sys,Kd)));
        err = sdh2norm ( sys, Kd );
        fprintf('Sampled-data H2 cost:    %g\n',err);    
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
        Tmax = 10;
        opt = simset('SrcWorkspace', 'current');
        warning off;
        [nK,dK] = tfdata ( tf(Kopt), 'v' );        
        sim ( 'sim_doubint', [], opt );
        Out{1} = out;
        [nK,dK] = tfdata ( tf(Kc2d), 'v' );        
        sim ( 'sim_doubint', [], opt );
        Out{2} = out;
        [nK,dK] = tfdata ( tf(Kd), 'v' );        
        sim ( 'sim_doubint', [], opt );
        Out{3} = out;
        warning on;
%------------------------------------------------------
%%       Plot curves
%------------------------------------------------------
        close all;
        titles = {'Sampled-data design', ...
                  'Redesign technique', ...
                  'Discrete-time design' };
        figure;
        for i=1:3
          y = Out{i};  
          subplot(3, 2, 2*(i-1)+1);
          plot(y(:,1),y(:,2));
          xlabel('Time, t');
          ylabel('Output');
          title(titles{i});
          subplot(3, 2, 2*(i-1)+2);
          plot(y(:,1),y(:,3));
          xlabel('Time, t');
          ylabel('Control');
        end;
        
%------- End of DEMO_DOUBINT.M --------- KYuP ----------
        





