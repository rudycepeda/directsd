%DEMO_L2 Design of L2-optimal controller.
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
        disp('   Sampled-data L2-optimization .')
        disp('#############################################################')
%------------------------------------------------------
%%       Initial data            
%------------------------------------------------------
        F = tf ( 1, [5 1 0] );
        Q = tf ( 1, [1 1] );
        R = tf ( 1, [1 0] );
        T = 0.2;
        rho = 0;
        if paused
          disp(' ')
          disp('# Plant:'); F
          disp('# Step reference:'); R
          disp('# Sampling period:'); T
          disp('# Control weight:'); rho
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       Form standard system
%------------------------------------------------------
        P11 = [minreal(Q*R); 0];
        P12 = [-F; rho];
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
%%       Sampled-data L2-optimization
%------------------------------------------------------                             
        disp(' ')
        disp('=========================================')
        disp('       Sampled-data L2-optimization')
        disp('=========================================')
        [K,errOpt] = sdl2 ( sys, T );
        K
        disp('Closed-loop poles');    
        disp(roots(charpol(sys,K)));
        fprintf('Optimal L2 cost:    %g\n',errOpt);    
        err = sdl2err ( sys, K );
        fprintf('Direct calculation: %g\n',err); 
        errdss = sdl2err ( sys, K, [], 'ss' );
        fprintf('L2 cost (ss):       %g\n',errdss);    
        [KWH,errOpt] = sdl2 ( sys, T, [], 'wh' )
        if paused
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       Sampled-data L2-optimization using lifting
%------------------------------------------------------                             
        disp(' ')
        disp('=========================================')
        disp('  Sampled-data L2-optimization (lifting)')
        disp('=========================================')
        sysH2s = sdh2simple ( sys, T );
        Ks = zpk( h2reg(sysH2s) );
        Ks
        disp('Closed-loop poles');    
        disp(roots(charpol(sys,Ks)));
        fprintf('Optimal L2 cost:    %g\n',errOpt);    
        err = sdl2err ( sys, Ks );
        fprintf('Direct calculation: %g\n',err); 
        errdss = sdl2err ( sys, Ks, [], 'ss' );
        fprintf('L2 cost (ss):       %g\n',errdss);    
        if paused
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       State-space sampled-data simple H2-optimization
%------------------------------------------------------                             
        disp(' ')
        disp('=========================================')
        disp('   Sampled-data simple H2-optimization')
        disp('=========================================')
        dsysL2 = sdh2simple ( sys, T );
        Kss = zpk(h2reg(dsysL2, 1, 1))
        errss = sdl2err ( sys, Kss );
        fprintf('Optimal L2 cost:    %g\n',errss);    
        if paused
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       Optimal redesign approach
%------------------------------------------------------                             
        disp(' ')
        disp('=========================================')
        disp('        Optimal redesign approach')
        disp('=========================================')
        Kc = ch2 ( sys );
        Kcd = bilintr(Kc, 'tustin', T);
        disp('Discretized C-optimal controller');    
        Kcd
        disp('Closed-loop poles');    
        rtsC = roots(charpol(sys,Kcd));
        disp(rtsC);
        errc = sdl2err ( sys, Kcd );
        fprintf('L2 cost:       %g\n',errc);                    
        errcss = sdl2err ( sys, Kcd );
        fprintf('L2 cost (ss):  %g\n',errcss);    
        if paused
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       Discrete-time L2-optimization
%------------------------------------------------------                             
        disp(' ')
        disp('=========================================')
        disp('        Discrete-time optimization')
        disp('=========================================')
        disp(' ')
        dsys = c2d ( sys, T );
        Kd = h2reg ( dsys, 1, 1);
        disp('Optimization based on discrete-time model');    
        Kd = zpk(Kd)
        disp('Closed-loop poles');    
        disp(roots(charpol(sys,Kd)));
        errd = sdl2err ( sys, Kd );
        fprintf('L2 cost:       %g\n',errd);            
        errdss = sdl2err ( sys, Kd );
        fprintf('L2 cost (ss):  %g\n',errdss);    
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
        G = 1;
        Fo = 1;
        Qu = 0;
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
        
        opt = simset('SrcWorkspace', 'current');
        warning off;
        [nK,dK] = tfdata ( tf(K), 'v' );        
        sim ( 'sim_track', [], opt );
        Out{1} = y; U{1} = u;

        [nK,dK] = tfdata ( tf(Kcd), 'v' );        
        sim ( 'sim_track', [], opt );
        Out{2} = y; U{2} = u;
        
        [nK,dK] = tfdata ( tf(Kd), 'v' );        
        sim ( 'sim_track', [], opt );
        Out{3} = y; U{3} = u;

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
          u = U{i};  
          subplot(3, 2, 2*(i-1)+1);
          plot(y(:,1),y(:,2));
          set(gca,'yLim',[0 1.2]);
          if i == 3, xlabel('Time, t'); end;
          ylabel('Output');
          title(titles{i});
          subplot(3, 2, 2*(i-1)+2);
          plot(u(:,1),u(:,2));
          if i == 3, xlabel('Time, t'); end;
          ylabel('Control');
        end;
        
%------- End of DEMO_L2.M --------- KYuP ----------
        





