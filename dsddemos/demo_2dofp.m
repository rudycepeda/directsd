%DEMO_2DOFP 2-DOF optimal control with preview.

% References:
% [1] K. Y. Polyakov, E. N. Rosenwasser, and B.P. Lampe,
%     Optimal design of 2-DOF digital controller for sampled-data tracking 
%	  systems with preview, Proc. 43th CDC, Bahama Isl., USA,
%     pp. 2352--2357, 2004.
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
        disp('   2-DOF optimal controller with preview')
        disp('#############################################################')
%------------------------------------------------------
%%       Initial data            
%------------------------------------------------------
        prevMax = 4; 
        F  = tf ( 1, [5 -1] );
        Q  = tf ( 1, [0.1 1] );
        R  = tf ( 1, [1 1 0] );
        T  = 1;
        F.iodelay = 1.5;
        preview  = 2;
        if paused
          disp(' ')
          disp('# Plant:'); F
          disp('# Signal model:'); R
          disp('# Ideal operator:'); Q          
          disp('# Preview:'); preview
          disp('# Sampling period:'); T
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%      Construct equivalent system with delays only
%------------------------------------------------------
        sigma = ceil(preview/T);
        theta = sigma*T - preview;
        R.iodelay = theta;
        Q.iodelay = preview;
%------------------------------------------------------
%%       Standard system
%------------------------------------------------------
        sys1DOF = [Q*R  -F
                     R  -F];
        k = 0;
        sys = [Q*R -F
                 R  0
               k*R -F];
        if paused
          clc
          disp('# Standard 2-DOF system:'); sys
          disp('Press any key ...')
          pause
          clc
        end;
%------------------------------------------------------
%%       Optimal 1-DOF controller
%------------------------------------------------------
        disp('=========================================')
        disp('          1-DOF L2-optimization')
        disp('=========================================')
        disp(' ')
        [K,errOpt] = sdl2 ( sys1DOF, T );
        disp('Optimal controller');    
        K
        disp('Closed-loop poles');    
        disp(roots(charpol(sys1DOF,K)));
        fprintf('Optimal L2 cost:    %g\n',errOpt);    
        err = sdl2err ( sys1DOF, K );
        fprintf('Direct calculation: %g\n',err); 
        if paused
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       Optimal 2-DOF controller
%------------------------------------------------------
        disp(' ')
        disp('=========================================')
        disp('           2-DOF optimization')
        disp('=========================================')
        disp(' ')
        [KR,errOpt] = sd2dof ( sys, K );
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
disp('    <<<  Simulating the transients of both systems. Please wait ... >>>')
        disp(' ')
        [nR,dR]   = tfdata ( tf(R), 'v' );
        [nFo,dFo] = tfdata ( tf(1), 'v' );
        [nF,dF]   = tfdata ( tf(F), 'v' );
        [nG,dG]   = tfdata ( tf(1), 'v' );
        [nQ,dQ]   = tfdata ( tf(Q), 'v' );
        [nQu,dQu] = tfdata ( tf(0), 'v' );
        [nK,dK]   = tfdata ( tf(K), 'v' );        
        Tmax = 20;
        tauR  = R.iodelay;
        tauF  = F.iodelay;
        tauG  = 0; %G.iodelay;
        tauQ  = Q.iodelay;
        tauQu = 0; %Qu.iodelay;

        opt = simset('SrcWorkspace', 'current');
        warning off;
        sim ( 'sim_track', [], opt );
        warning on;        
        y1dof = y; u1dof = u;       
%---------------------------------------------------
%%       Simulate  2-DOF system
%---------------------------------------------------
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
        plot(y(:,1),y(:,3),'b--',y1dof(:,1),y1dof(:,2),'r-',y(:,1),y(:,2),'g-');
        xlabel('Time, t');
        ylabel('Output');
        legend('Ideal response', '1-DOF control', '2-DOF control');
        if paused
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       Preview time array 
%------------------------------------------------------
        hPrev = 0.1;
        Nh = T / hPrev;
        arPrev = sort ( [(0:hPrev:prevMax) (0.9999:1:prevMax)] );
        Jopt = [];
        disp(' ')
        disp('Calculating series of L2-optimal predictive controllers ...')
        disp(' ')
        for preview=arPrev
          sigma = ceil(preview/T);
          theta = sigma*T - preview;
          R.iodelay = theta;
          Q.iodelay = preview;
          sys = [-Q*R  F
                   R   0
                  k*R -F];
          [KRt,Jt,Jprev] = sd2dof ( sys, K );
          fprintf('tau = %3.1f: Jopt = %g\n', preview, Jt );
          Jopt = [Jopt; [Jt Jprev]];
        end;                
        ind = 1:Nh+1:(prevMax/hPrev+ceil(prevMax)+1);
%------------------------------------------------------
%%       Plot curves
%------------------------------------------------------
        close all;
        figure;
        plot(arPrev,Jopt(:,1),'-',arPrev,Jopt(:,2),'--',arPrev(ind),Jopt(ind,1),'o');
        ylabel('Optimal cost J_{min}(\pi)');
        xlabel('Preview interval, \pi');
        
%------- End of DEMO_2DOFP.M --------- KYuP ----------
