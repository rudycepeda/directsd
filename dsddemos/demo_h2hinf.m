%DEMO_H2HINF Mixed H2/AHinf-optimization.
%

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
        disp('   Mixed H2/Hinf optimization.')
        disp('#############################################################')
%------------------------------------------------------
%%       Initial data            
%------------------------------------------------------
        F = tf ( 1, [5 1 0] );
        rho = 1;
        Fw = 1;
        T  = 1;
        if paused
          disp(' ')
          disp('# Plant:'); F
          disp('# Wave spectrum:'); Fw
          disp('# Control weight:'); rho
          disp('# Sampling period:'); T
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       Form standard system (H2 and Hinf channels coincide)
%------------------------------------------------------
        P11 = [-F*Fw; 0];
        P12 = [-F; rho];
        P21 = - F*Fw;
        P22 = - F;
        sys = [P11 P12 
               P21 P22];
        if paused
          clc
          disp('# Standard system:'); sys
          disp('Press any key ...')
          pause
          clc
        end;
        sys = minreal(sys);
%------------------------------------------------------
%%       H2-optimization
%------------------------------------------------------                             
        disp('=========================================')
        disp('            H2-optimization')
        disp('=========================================')
        [KH2,errOpt] = sdh2 ( sys, T );
        KH2
        fprintf('H2-cost:            %g\n',errOpt);    
        err = sdh2norm ( sys, KH2 );
        fprintf('Direct computation: %g\n',err);    
        lamH2 = sdahinorm ( sys, KH2 );
        fprintf('AHinf cost:         %g\n',lamH2);    
        lamH2x = sdnorm ( sys, KH2, 'inf' );
        fprintf('Hinf norm:          %g\n',lamH2x);    
%------------------------------------------------------
%%       Hinf-optimization
%------------------------------------------------------                             
        disp(' ')
        disp('=========================================')
        disp('           Hinf-optimization')
        disp('=========================================')
        [Kinf,lam] = sdahinf ( sys, T );
        Kinf
        err = sdh2norm ( sys, Kinf );
        fprintf('H2-cost:              %g\n',err);    
        fprintf('Optimal AHinf-cost:   %g\n',lam);    
        lam1 = sdahinorm ( sys, Kinf );
        fprintf('Direct computation:   %g\n',lam1);    
        lam2 = sdnorm ( sys, Kinf, 'inf' );
        fprintf('Hinf norm:            %g\n',lam2);    
        if paused
          disp('Press any key ...')
          pause
        end;
%------------------------------------------------------
%%       H2/Hinf-optimization
%------------------------------------------------------                             
        disp(' ')
        disp('=========================================')
        disp('          H2/Hinf-optimization')
        disp('=========================================')
        [Kmix,errOpt] = sdh2hinf ( sys, T, 0.5, 2, 1 );
        Kmix
        errMix = sdh2norm ( sys, Kmix );
        fprintf('H2-cost:              %g\n',errMix);    
        lamMix = sdahinorm ( sys, Kmix );
        fprintf('AHinf cost:           %g\n',lamMix);    
        lam1 = sdnorm ( sys, Kmix, 'inf' );
        fprintf('Hinf norm:            %g\n',lam1);    

%------------------------------------------------------
%%       Sampled-data H2/AHinf-optimization
%------------------------------------------------------                             
        disp(' ')
        disp('=========================================')
        disp('     Sampled-data H2/AHinf-optimization')
        disp('=========================================')
        disp(' ')
        arRho = (0:0.05:1);
        h2err = []; hinferr = [];        
        for rho=arRho     
          rho  
          Kmix = sdh2hinf(sys, T, rho, 2, 1);
          disp('H2/AHinf-optimal controller');    
          Kmix
          disp('Closed-loop poles');    
          disp(roots(charpol(sys,Kmix)));        
          err = sdh2norm ( sys, Kmix );
          fprintf('Sampled-data H2 cost:    %g\n',err);    
          h2err = [h2err; err];
          err = sdahinorm ( sys, Kmix );
          fprintf('Sampled-data AHinf norm: %g\n',err);    
          hinferr = [hinferr; err];
          err = sdhinorm ( sys, Kmix );
          fprintf('L2-induced norm:         %g\n',err);    
          disp(' ')
        end;
%------------------------------------------------------
%%       Draw curves
%------------------------------------------------------                             
        figure(1);
        plot(arRho, h2err, 'b-', arRho, hinferr, 'g--');
        set(gca,'FontName','Arial','FontSize',18);
        xlabel('$\rho$','Interpreter', 'latex');
        h = legend('${\mathcal H}_2$-norm','${\mathcal{AH}}_\infty$-norm');
        set(h,'Interpreter', 'latex');
        h = get(gca,'children');
        set(h(1),'LineWidth',1.5);
        set(h(2),'LineWidth',1.5);        
        figure(2);
        plot(h2err, hinferr);
        set(gca,'FontName','Arial','FontSize',18);
        xlabel('${\mathcal H}_2$-norm','Interpreter', 'latex');
        ylabel('${\mathcal{AH}}_\infty$-norm','Interpreter', 'latex');
        h = get(gca,'children');
        set(h(1),'LineWidth',1.5);
        if paused
          disp('Press any key ...')
          %pause
        end;
        
%------- End of DEMO_H2HINF.M --------- KYuP ----------
        





