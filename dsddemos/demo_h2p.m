%DEMO_H2P H2-optimal control with preview.

% References: 
% [1] K. Y. Polyakov, E. N. Rosenwasser, and B.P. Lampe,
%     Optimal stochastic sampled-data control with preview,
%     submitted to IEEE AC, 2002.  
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
        disp('   H2-optimal controller with preview')
        disp('#############################################################')
%------------------------------------------------------
%%       Initial data            
%------------------------------------------------------
        exNo = 2;
        if exNo == 1
          prevMax = 4; 
          F  = sdtf ( 1, [1 1] );
          F.iodelay = 1.5;
          yLim = [0 0.07];
        else
          prevMax = 10; 
          F = sdtf ( 1, [1 -1] );
          F.iodelay = 1.5;
          yLim = [11.5 13.2];
        end;
        Fr = sdtf ( 1, [5 1] );
        Fn = 0.2;
        T  = 1;
        preview  = 2;
        Q  = sdtf ( 1 );
        Q.iodelay = preview;
        if paused
          disp(' ')
          disp('# Plant:'); F
          disp('# Signal spectrum:'); Fr
          disp('# Noise spectrum:'); Fn
          disp('# Ideal operator:'); Q
          disp('# Preview:'); preview
          disp('# Sampling period:'); T
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       Standard system
%------------------------------------------------------
        sys = [Q*Fr  0  -F
                 Fr Fn  -F];
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
        [K,errOpt] = sdh2 ( sys, T );
        disp('Optimal controller');    
        K
        disp('Closed-loop poles');    
        disp(roots(charpol(sys,K)));
        fprintf('Optimal H2 cost:    %g\n',errOpt^2);    
        err = sdh2norm ( sys, K );
        fprintf('Direct calculation: %g\n',err^2); 
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
        disp('Calculating series of H2-optimal predictive controllers ...')
        disp(' ')
        for preview=arPrev
          Q.iodelay = preview;
          sys = [Q*Fr  0 -F
                   Fr Fn -F];
          [Kt,Jt,Jprev] = sdh2 ( sys, T );
          fprintf('tau = %3.1f: Jopt = %g\n', preview, Jt^2 );
          Jopt = [Jopt; [Jt^2 Jprev^2]];
        end;                
        ind = 1:Nh+1:(prevMax/hPrev+ceil(prevMax)+1);
%------------------------------------------------------
%%       Plot curves
%------------------------------------------------------
        close all;
        figure;
        plot(arPrev,Jopt(:,1),'-',arPrev,Jopt(:,2),'--',arPrev(ind),Jopt(ind,1),'o');
        set(gca,'FontName','Arial','FontSize',18,...
            'xLim', [0 prevMax], 'yLim', yLim);
        ylabel('Optimal cost J_{min}(\pi)');
        xlabel('Preview interval, \pi');
        h = get(gca,'children');
        set(h(1),'LineWidth',1.5);
        set(h(2),'LineWidth',1.5);        
        set(h(3),'LineWidth',1.5);        
        
%------- End of DEMO_H2P.M --------- KYuP ----------
