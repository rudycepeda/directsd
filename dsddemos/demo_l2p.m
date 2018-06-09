%DEMO_L2P L2-optimal control with preview.

% References:
% [1] K. Y. Polyakov, E. N. Rosenwasser, and B.P. Lampe,
%     Design of optimal sampled-data tracking systems with preview,
%     Proc. 4th IFAC Workshops on Time-Delay Systems, 
%     Rocquencourt, France, September 8-10, 2003.  
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
        disp('   L2-optimal controller with preview')
        disp('#############################################################')
%------------------------------------------------------
%%       Initial data            
%------------------------------------------------------
        exNo = 2;
        if exNo == 1
          prevMax = 4; 
          F  = tf ( 1, [5 1] );
          F.iodelay = 1.5;
          ymax = 1.2;
        else
          prevMax = 12; 
          F = tf ( 1, [5 -1] );
          F.iodelay = 1.5;
          ymax = 2.2;
        end;
        Q = tf ( 1, [0.1 1] );
        R = tf ( 1, [1 1 0] );
        T = 1;
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
        sys = [Q*R  -F
                 R  -F];
        if paused
          clc
          disp('# Standard system:'); sys
          disp('Press any key ...')
          pause
          clc
        end;
%------------------------------------------------------
%%       Optimal controller
%------------------------------------------------------
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
        yOpt = y; uOpt = u;
        warning on;
%------------------------------------------------------
%%       Plot curves
%------------------------------------------------------
        close all;
        figure;
        plot(y(:,1),y(:,3),'--',yOpt(:,1),yOpt(:,2),'-');
        ylabel('Output');
        xlabel('Time');
        if paused
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       Preview time array 
%------------------------------------------------------
        hPrev = 0.1;
        Nh = T / hPrev;
        arPrev = sort ( [0:hPrev:prevMax 0.9999:1:prevMax] );
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
                   R  -F];
          [Kt,Jt,Jprev] = sdl2 ( sys, T );
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
        
%------- End of DEMO_L2P.M --------- KYuP ----------
