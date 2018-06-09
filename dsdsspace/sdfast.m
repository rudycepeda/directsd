function [dsys,DP11,DP12,DP21,DP22] = sdfast ( csys, T, Nfast, o2, i2, type, H, Hf )
%SDFAST Fast discretization for sampled-data system. 
%
%     [DSYS,DP11,DP12,DP21,DP22] = SDFAST ( SYS, T, NF, NMEAS, NCON, TYPE )
%     [DSYS,DP11,DP12,DP21,DP22] = SDFAST ( SYS, T, NF, NMEAS, NCON, TYPE, H, HF )
%
%   Inputs:
%     SYS - continuous-time standard system 
%             | z | = | P11 P12 | | w |
%             | y |   | P21 P22 | | u |
%  	  T     - sampling period
%     NF    - frequency multiplication factor for fast discretization (default 5)
%     NMEAS - dimension of the vector 'y' (default 1)
%     NCON  - dimension of the vector 'u' (default 1)
%     TYPE - discretization method of
%         'std' - standard technique, discretize input and 
%                 output [1] (default)
%         'l2'  - discretize input and preserve L2-norm of the 
%                 output [2]
%     H  - hold device for main loop as an LTI system (default ZOH)
%     HF - hold device for fast sampling as an LTI system (default ZOH)
%
%   Outputs:
%	  DSYS  - discrete-time fast sampled model 
%             | zd | = | DP11 DP12 | | dw |
%             | yd |   | DP21 DP22 | | du |
%
%   See also SDH2SIMPLE, SDGH2MOD, SDH2REG.

% References:
% [1] T. Chen, B. Francis, Optimal Sampled-Data Control Systems, 
%     Springer-Verlag, Berlin etc., 1995.
% [2] S. Hara, H. Fujioka, P. Khargonekar, Y. Yamamoto, Computational
%     aspects of gain-frequency response for sampled-data systems //
%     Proc 34th CDC, New Orleans, LA, pp. 1784-1789, 1995.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%	Continuous-time state space
%------------------------------------------------------
	    sys = minreal(ss(csys), sqrt(eps), 0);
	    if ~isct(sys)
	       error('Initial system must be continuous time');
	    end;
        [a,b,c,d] = ssdata(sys);
        n = size(a,1);
        [nout,nin] = size(sys);        
        if ~exist('type', 'var'), type = 'std'; end;
        if ~isequal(type, 'std')  &&  ~isequal(type, 'l2')
           error('Unknown type of discretized model ''%s''', type);  
        end;        
%------------------------------------------------------
%	Determine inputs and outputs
%------------------------------------------------------
        if ~exist('i2','var'), i2 = 1; end;
        if ~exist('o2','var'), o2 = 1; end;
%------------------------------------------------------
%	Define extrapolator for the main loop
%------------------------------------------------------
        HZOH = 0;
        if ~exist('H','var'),    % ZOH 
          H = ss ( 0, 1, 1, 0 ); 
          HZOH = 1;
        end; 
        H = ss ( H );
        if isscalar(H)
           if ~HZOH, HZOH = isequal(H,ss(0,1,1,0)); end;
           H = H * eye(i2); 
        else
           if ~HZOH, HZOH = isequal(H,eye(i2)*ss(0,1,1,0)); end;
           [i2,iH] = size(H,1);
        end;
        [ah,bh,ch] = ssdata ( H );
        nh = size(ah,1);
%------------------------------------------------------
%	Check validity of inputs and outputs
%------------------------------------------------------
        canBeZero = 1;
        [o1,i1] = checkio ( nout, nin, o2, i2, canBeZero);
%------------------------------------------------------
%	Define extrapolator for fast sampling
%------------------------------------------------------
        HfZOH = 0;
        if ~exist('Hf','var');     % ZOH 
          Hf = ss ( 0, 1, 1, 0 ); 
          HfZOH = 1;
        end; 
        Hf = ss ( Hf );
        if isscalar(Hf)
           if ~HfZOH, HfZOH = isequal(Hf,ss(0,1,1,0)); end;
           Hf = Hf * eye(i1); % input sampling hold
        else
           if ~HfZOH, HfZOH = isequal(Hf,eye(i1)*ss(0,1,1,0)); end;
           [oHf,iHf] = size(Hf,1); 
           if oHf ~= i1
              error('Fast sampling hold improperly dimensioned');
           end;
        end;
        [ahf,bhf,chf] = ssdata ( Hf );
        nhf = size(ahf,1);
%------------------------------------------------------
%	Extract partial matrices
%	ans check assumptions
%------------------------------------------------------
	    [b1,b2,c1,c2,d11,d12,d21,d22] = blocks4 ( b, c, d, o2, i2 );
%------------------------------------------------------
%   Preliminary slow rate discretization (generalzied hold)
%------------------------------------------------------
        if ~exist('Nfast','var'), Nfast = 5; end;
        Tf = T / Nfast;
        AB = [a b2*ch
              zeros(nh,n) ah];
        E = expm(AB*T);
        Ad = E(1:n, 1:n);
        B2d = E(1:n, n+1:end) * bh; 
        %[Ad,B2d] = c2d ( a, b2, T ); % ZOH
%------------------------------------------------------
%   Fast sampling always with generalized hold: Af, B1f, B2f
%------------------------------------------------------
        if HZOH  &&  HfZOH
          [Af,Bf] = c2d ( a, [b1 b2], Tf ); 
          B1f = Bf(:, 1:i1);
          B2f = Bf(:, (i1+1):nin);
        else
          if HfZOH
            [Af,B1f]  = c2d ( a, b1, Tf ); 
          else
            AB1f = [a b1*chf
                    zeros(nhf,n) ahf];
            E1 = expm(AB1f*Tf);
            Af = E1(1:n, 1:n);
            B1f = E1(1:n, n+1:end) * bhf; 
          end;
          if HZOH
             [Af,B2f]  = c2d ( a, b2, Tf ); 
          else
             AB2f = [a b2*ch
                     zeros(nh,n) ah];          
          end;
        end;
%------------------------------------------------------
%   Initialization for both the methods
%------------------------------------------------------      
        B1d = [];
        C1d = [];        
        D21d = [d21 zeros(o2,(Nfast-1)*i1)];
        D12d = d12*ch*bh;     % h(0)=Ch*Bh
        D12i = zeros(size(d12));         
%------------------------------------------------------
%   Compute system matrices for the standard method
%   by Chen and Francis
%------------------------------------------------------      
        if isequal(type,'std')
           D11d = d11*chf*bhf;  % hf(0)=Chf*Bhf
           D11i = D11d;
           Ai = eye(n);
           for i=0:Nfast-1 %--- i is the degree of Ai = af^i
             if i < Nfast-1
                if HZOH
                  D12i = c1*Ai*B2f + D12i;  
                else
                  iTf = (i+1)*Tf;
                  E2 = expm(AB2f*iTf);          
                  X = c1*E2(1:n,n+1:end)*bh;
                  D12i = X + d12*ch*expm(ah*iTf)*bh;
                end;
                D12d = [D12d; D12i];
                D11i = [c1*Ai*B1f D11i];
                D11d = [D11d zeros(o1*(i+1),i1)
                            D11i  ];
             end;
             B1d = [Ai*B1f B1d]; 
             C1d = [C1d; c1*Ai];              
             Ai = Ai * Af;  
           end;
        end;
%------------------------------------------------------
%   Compute system matrices using the formulas
%   by Hara et al. - L2-norm preserving output discretization
%------------------------------------------------------      
        if isequal(type,'l2')
           ne = n + i2;
           nbar = n + i2 + i1;
          %-------------------------------------------- 
          % Compute Lam1 and Lam2 using matrix exponential
          %-------------------------------------------- 
           if HfZOH
             bhf = bhf / sqrt(Tf);
             B1f = B1f / sqrt(Tf);
           else
             Y = intaxa ( ahf', chf'*chf, Tf );
             Y = bhf'*Y*bhf;
             F = rectchol(Y);
             bhf = bhf*pinv(F);
             B1f = B1f*pinv(F);
           end;
          %-------------------------------------------- 
          % Form Ae, Be and Abar
          %-------------------------------------------- 
           Ae = [a b2*ch; zeros(i2,ne)];
           Be = [b1*chf; zeros(i2,i1)];
           Abar = [Ae           Be; 
                   zeros(i1,ne) ahf];
          %-------------------------------------------- 
          % Compute Lam1 and Lam2 using matrix exponbarntial
          %-------------------------------------------- 
           CDe = [c1 d12*ch d11*chf];
           Q = CDe'*CDe;
           [LL,E,R11,R22] = intaxa ( Abar', Q, Tf );
           Lam12 = rectchol(LL); 
           Lam1 = Lam12(:,1:ne);
           Lam2 = Lam12(:,ne+1:end);
           Phie = R22(1:ne,1:ne);
           Gammae = R22(1:ne,ne+1:end);
          %-------------------------------------------- 
          % Generalized holds
          %-------------------------------------------- 
           vBH = blckdiag ( eye(n), bh );
           Lam2f = Lam2 * bhf;
           Gammaef = Gammae*bhf;
          %-------------------------------------------- 
          % Compute B1d, C1d, D12d and D11d
          %-------------------------------------------- 
           D11d = Lam2f;
           D11i = Lam2f;
           o1d = size(Lam1,1);
           C1D12d = [];        
           Ai = eye(n);
           Phii = eye(ne);
           for i=0:Nfast-1 %--- i is the degree of Ai = af^i
             if i < Nfast-1
                D11i = [Lam1*Phii*Gammaef D11i];
                D11d = [D11d zeros(o1d*(i+1),i1)
                            D11i  ];
             end;
             B1d = [Ai*B1f B1d]; 
             C1D12d = [C1D12d; Lam1*Phii*vBH];              
             Phii = Phii * Phie;  
             Ai = Ai * Af;  
           end;
           C1d = C1D12d(:,1:n);
           D12d = C1D12d(:,n+1:end);
        end;
%------------------------------------------------------
%	Final fast discretized discrete model
%------------------------------------------------------   
        Bd = [B1d B2d];
        Cd = [C1d; c2];
        Dd = [D11d D12d  
              D21d d22];
        dsys = minreal ( ss ( Ad, Bd, Cd, Dd, T ), sqrt(eps), 0 );       
%------------------------------------------------------
%	Additional data if necessary
%------------------------------------------------------   
        if nargout > 1
           DP11 = minreal( ss(Ad, B1d, C1d, D11d, T), sqrt(eps), 0 );
        end;    
        if nargout > 2
           DP12 = minreal( ss(Ad, B2d, C1d, D12d, T), sqrt(eps), 0 );
        end;
        if nargout > 3
           DP21 = minreal( ss(Ad, B1d, c2, D12d, T), sqrt(eps), 0 );
        end;
        if nargout > 4
           DP22 = minreal( ss(Ad, B2d, c2, d22, T), sqrt(eps), 0 );            
        end;      
    
%------- End of SDFAST.M --------- KYuP ----------
