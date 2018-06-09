%DEMO_AIT98 H2- and AHinf-optimization for sampled-data system.

% References:
% [1]  K. Yu. Polyakov, Polynomial design of sampled-data 
%      systems. I & II. Automation & Remote Control,
%      1998; no. 10, pp. 74-89; no. 12, pp. 94-108. 
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
        disp('   H2 and AHinf-optimization of sampled-data system.')
        disp('#############################################################')
%------------------------------------------------------
%%       Initial data            
%------------------------------------------------------
        F1 = tf(1);
        F2 = tf ( 1, [1 -1] );
        G  = 1;
        Ve = 1;
        Vu = 0;
        Fw = 1;
        T  = 1;
        if paused
          disp(' ')
          disp('# Plant:'); F2
          disp('# Wave spectrum:'); Fw
          disp('# Sampling period:'); T
          disp('# Control weight:'); Vu
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       Form standard system
%------------------------------------------------------
        P11 = [minreal(Ve*F2*Fw); 0];
        P12 = [minreal(Ve*F2*F1); Vu];
        P21 = -minreal(G*F2*Fw);
        P22 = -minreal(G*F2*F1);
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
        disp('=========================================')
        disp('           H2-optimization')
        disp('=========================================')
        disp(' ')
        [K,errOpt] = sdh2 ( sys, T );
        disp('H2-optimal controller');    
        K
        fprintf('H2-cost:            %g\n',errOpt);    
        err = sdh2norm ( sys, K );
        fprintf('Direct computation  %g\n',err);    
        lam = sdahinorm ( sys, K );
        fprintf('AHinf-cost:         %g\n',lam);    
        lam1 = sdnorm ( sys, K, 'inf' );
        fprintf('Hinf-norm:          %g\n',lam1);    
%------------------------------------------------------
%%       Hinf-optimization
%------------------------------------------------------                             
        disp(' ')
        disp('=========================================')
        disp('           AHinf-optimization')
        disp('=========================================')
        disp(' ')
        [Kinf,lam] = sdahinf ( sys, T );
        disp('Hinf-optimal controller');    
        Kinf
        err  = sdh2norm ( sys, Kinf );
        fprintf('H2-cost:             %g\n',err);    
        fprintf('Optimal AHinf-cost:  %g\n',lam);    
        lam1 = sdahinorm ( sys, Kinf );
        fprintf('Direct computation:  %g\n',lam1);    
        lam2 = sdnorm ( sys, Kinf, 'inf' );
        fprintf('Hinf-norm:           %g\n',lam2);    
        
%------- End of DEMO_AIT98.M --------- KYuP ----------
        





