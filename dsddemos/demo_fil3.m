%DEMO_FIL3 Optimal sampled-data filtering for 2nd order plant.
%

% References:
% [1]  Y.N. Rosenwasser, K.Y. Polyakov, and B.P. Lampe,
%      Optimal discrete filtering for time-delayed 
%      systems with respect to mean-square continuous-time 
%      error criterion, Int. J. Adapt. Contr. Signal Proc.,
%      1998, vol. 12, pp. 389-406. 
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
        disp('   Optimal sampled-data filtering for 2nd order plant.')
        disp('#############################################################')
%------------------------------------------------------
%%       Initial data            
%------------------------------------------------------
        F = tf ( 1, [1 1  1] );
        Fw = tf ( 2, [1 2] );
        Fn = 1;
        Q  = 1;
        T  = 0.1;
        if paused
          disp(' ')
          disp('# Plant:'); F
          disp('# Signal spectrum:'); Fw
          disp('# Noise spectrum:'); Fn
          disp('# Sampling period:'); T
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       Form standard system
%------------------------------------------------------
        P11 = [-minreal(Q*Fw) 0];
        P12 = F;
        P21 = [Fw Fn];
        P22 = 0;
        sys = [P11 P12;P21 P22];
        if paused
          clc
          disp('# Standard system:'); sys
          disp('Press any key ...')
          pause
          clc
        end;
%------------------------------------------------------
%%       H2-optimization
%------------------------------------------------------                             
        disp(' ')
        disp('=========================================')
        disp('           H2-optimization')
        disp('=========================================')
        [K,errOpt] = sdh2 ( sys, T );
        K
        fprintf('Optimal H2-cost:    %g\n',errOpt);    
        err = sdh2norm ( sys, K, [], [], 'pol' );
        fprintf('Average variance:   %g\n',err);    
        err0 = sdh2norm ( sys, K, 0, [], 'pol');
        fprintf('Variance at t=0:    %g\n',err0);    
        disp('Computing variance between sampling instants ...');
        t = (0:0.025:1)*T;
        err_t = sdh2norm ( sys, K, t, [], 'pol' );
        disp('OK')
%------------------------------------------------------
%%       Discrete H2-optimization for t=0
%------------------------------------------------------                             
        disp(' ')
        disp('=========================================')
        disp('         H2-optimization for t=0')
        disp('=========================================')
        [K0,errOpt] = sdh2 ( sys, T, 0 );
        K0
        fprintf('Variance at t=0:    %g\n',errOpt);    
        err0 = sdh2norm ( sys, K0, 0, [], 'pol' );
        fprintf('Direct computation  %g\n',err0);    
        err = sdh2norm ( sys, K0, [], [], 'pol' );
        fprintf('Average variance:   %g\n',err);    
        disp('Computing variance between sampling instants ...');
        err_t0 = sdh2norm ( sys, K0, t, [], 'pol' );
        disp('OK')
%------------------------------------------------------
%%       Plotting curves
%------------------------------------------------------                             
        plot ( t, err_t, t, err_t0, '--' );
        xlabel('Time'); ylabel('Variance');
        title('Variation of variance between sampling instants');
        legend('CMV-optimal','MV-optimal')
                
%------- End of DEMO_FIL3.M --------- KYuP ----------
        





