%DEMO_DHINF Discrete-time polynomial Hinf-optimization.
%

% References:
% [1]  M. Grimble, Robust Industrial Control: Optimal Design 
%      Approach for Polynomial Systems. UK: Prentice-Hall, 
%      Hemel Hempstead, 1994.
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
        disp('   Polynomial Hinf-design for discrete-time systems.')
        disp('#############################################################')
%------------------------------------------------------
%%       Example 1
%------------------------------------------------------                             
        disp(' ')
        disp('=========================================')
        disp('           Example 1')
        disp('=========================================')
        T = 1;
        F = tf ( [1 0],  [-2 1], T );
        Fw = tf ( [1 -2], [-2 1], T );
        V1 = 1;
        V2 = 1;
        if paused
          disp(' ')
          disp('# Plant:'); F
          disp('# Disturbance:'); Fw
          disp('# Error weight:'); V1
          disp('# Control weight:'); V2
          disp('Press any key ...')
          pause
        end;
        sys = [V1*Fw V1*F
                 0    V2
                -Fw   -F];
        sys = z2zeta(sdtf(sys));
        if paused
          disp(' ')
          disp('# Standard system:'); sys
          disp('Press any key ...')
          pause
        end;
       %-----------------------------------------------
       %   Polynomial design
       %-----------------------------------------------
        [K,errOpt] = dhinf ( sys );
        disp(' ')
        disp('# Hinf-optimal controller');    
        K
        fprintf('Optimal Hinf-cost:    %g\n',errOpt);            
        err = dahinorm ( sys, K );
        fprintf('Direct calculation:   %g\n',err);            
        fprintf('Closed-loop poles\n');    
        disp(roots(charpol(sys,K)))
        if paused
          disp('Press any key ...')
          pause
          clc;
        end;
%------------------------------------------------------
%%       Example 2
%------------------------------------------------------                             
        disp(' ')
        disp('=========================================')
        disp('           Example 2 (generic)')
        disp('=========================================')
        T = 1;
        F2 = tf ( 1, [-1 -2.1 1], T );
        F1 = tf ( [2 1 0], 1, T );
        Fw = tf ( [0.3 1], 1, T );
        V1 = tf ( 1, 1, T );
        V2 = tf ( 1, 1, T );
        if paused
          disp(' ')
          disp('# Plant:'); F1, F2
          disp('# Disturbance:'); Fw
          disp('# Error weight:'); V1
          disp('# Control weight:'); V2
          disp('Press any key ...')
          pause
        end;
        sys = [V1*F2*Fw V1*F2*F1
                   0       V2
                -F2*Fw   -F2*F1];
        sys = z2zeta(sdtf(sys));
        if paused
          disp(' ')
          disp('# Standard system:'); sys
          disp('Press any key ...')
          pause
        end;
       %-----------------------------------------------
       %   Polynomial design
       %-----------------------------------------------
        [K,errOpt] = dhinf ( sys );
        disp(' ')
        disp('# Hinf-optimal controller');    
        K
        fprintf('Optimal Hinf-cost:    %g\n',errOpt);            
        err = dahinorm ( sys, K );
        fprintf('Direct calculation:   %g\n',err);    
        fprintf('Closed-loop poles\n');    
        disp(roots(charpol(sys,K)))
        if paused
          disp('Press any key ...')
          pause
          clc;
        end;               
%------------------------------------------------------
%%       Example 3 (Non-generic case)
%------------------------------------------------------                             
        disp(' ')
        disp('=========================================')
        disp('         Example 3 (non-generic)')
        disp('=========================================')
        T = 1;
        F2 = tf ( 1, [0.32 -1.2 1], T );
        F1 = tf ( [-0.8 1 0], 1, T );
        Fw = tf ( [-0.4 1], 1, T );
        V1 = tf ( 1, 1, T );
        V2 = tf ( 1, 1, T );
        if paused
          disp(' ')
          disp('# Plant:'); F1, F2
          disp('# Disturbance:'); Fw
          disp('# Error weight:'); V1
          disp('# Control weight:'); V2
          disp('Press any key ...')
          pause
        end;
        sys = [V1*F2*Fw V1*F2*F1
                   0       V2
                -F2*Fw   -F2*F1];
        sys = z2zeta(sdtf(sys));
        if paused
          disp(' ')
          disp('# Standard system:'); sys
          disp('Press any key ...')
          pause
        end;
       %-----------------------------------------------
       %   Polynomial design
       %-----------------------------------------------
        [K,errOpt] = dhinf ( sys );
        nonGen = iscell(K);
        if nonGen, 
           K1 = K{2}; K = K{1}; 
           disp(' ')
           disp('# Nongeneric solution found!');
        end;
        disp(' ')
        disp('# Hinf-optimal controller');    
        K
        fprintf('Optimal Hinf-cost:    %g\n',errOpt);            
        err = dahinorm ( sys, K );
        fprintf('Direct calculation:   %g\n',err);    
        fprintf('Closed-loop poles\n');    
        disp(roots(charpol(sys,K)))
        if nonGen, 
          disp(' ')  
          disp('# Second Hinf-optimal controller');    
          K1
          fprintf('Optimal Hinf-cost:    %g\n',errOpt);            
          err = dahinorm ( sys, K1 );
          fprintf('Direct calculation:   %g\n',err);    
          fprintf('Closed-loop poles\n');    
          disp(roots(charpol(sys,K1)))
        end;
        
%------- End of DEMO_DHINF.M --------- KYuP ----------
        
