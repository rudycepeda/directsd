function [D,msg] = dtfm ( F, T, t, H, type, recFlag )
%DTFM Discrete Laplace transform with a hold.
%
%     [D,MSG] = DTFM ( F, T, t )
%     [D,MSG] = DTFM ( F, T, t, H )
%     [D,MSG] = DTFM ( F, T, t, TYPE )
%     [D,MSG] = DTFM ( F, T, t, H, TYPE )
%
%   Inputs:
%     F(s) - transfer matrix of continuous-time plant
%     T    - sampling period
%     t    - time parameter
%     H(s) - transfer matrix of the hold:
%              h(t)=K*exp(L*t)*M => H(s)=K*(sI-L)^{-1}*M
%            default hold is the scalar zero-order hold (ZOH) with H(s)=1/s.
%     TYPE - class of return value: 'zpk' (default), 'ss', or 'tf'
%
%   Outputs::
%      D(z) - discrete LTI model
%      MSG  - error or warning message if any
%
%   See also ZTRM, DTFM2.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 04-Dec-2007 $
% $Revision: 3.10 $    $Date: 09-Jan-2008 $
%------------------------------------------------------
%       Handle possible preview or delay
%------------------------------------------------------
        if ~exist('t','var') | isempty(t), t = 0; end;
        delay = totaldelay ( F );
        comdelay = min(min(delay));
        delay = delay - comdelay;        
%------------------------------------------------------
%       Transform to state space
%------------------------------------------------------
        if ~isproper(F), error('Transfer matrix must be proper!'); end;
	    F = ss(F);
        A = F.a; B = F.b; C = F.c; D = F.d;
%------------------------------------------------------
%       Extrapolator
%------------------------------------------------------
        if ~exist('H','var'), H = [];
        elseif ischar(H)
           type = H; H = [];  
        end;
        if isempty(H), 
             H = ss(0, 1, 1, 0); % scalar ZOH by default
        else H = ss(H); 
        end;     
        if size(H,1)+size(H,2) == 2,  % scalar extrapolator  
           scalarH = 1; H0 = H; 
           H = H * eye(size(F,2)); 
        else
           scalarH = 0; 
           H0 = H;
        end; % scalar extrapolator              
%------------------------------------------------------
%       Check input arguments
%------------------------------------------------------
        if ~exist('type','var')  ||  isempty(type) 
             type = 'zpk';         
        else type = lower(type); 
        end;
        if ~isequal(type,'zpk')  &&  ~isequal(type,'ss')  &&  ~isequal(type,'tf')
           error('Incorrect type ''%s''',type); 
        end;        
%------------------------------------------------------
%       Special case: first order hold
%       including fractional order hold
%	NB: so far only for delay-free case!
%------------------------------------------------------
        if size(H,3) > 1 
          if norm(comdelay+delay) > sqrt(eps)
	         error('Unable to handle first-order hold with delay');
          end;   
          [D,msg] = dtfmfoh ( F, T, t, H, type );
          return;
    	end;
%------------------------------------------------------
%       Check data compatibility
%------------------------------------------------------
        if size(F,2) ~= size(H,1), error('Dimensions of plant and hold are incompatible'); end;
        L = H.a; M = H.b; K = H.c;
        if ~exist('recFlag','var'), recFlag = 0; end;
        
%####################################################################        
%       Time-delayed case: if delays are different, use
%       row-wise or even element-wise transformations
%####################################################################
        
        if norm(delay) > sqrt(eps)
           [n,ell] = size(F); m = size(H,2);
           switch type
             case 'ss',         D = ss(zeros(n,m));  type0 = 'ss';
             case {'zpk','tf'}, D = zpk(zeros(n,m)); type0 = 'zpk';
           end;
           D.Ts = T;
           msg = '';
          %------------------------------------------------------
          %     Multirow case  
          %------------------------------------------------------         
           recFlag0 = 1;
           if n > 1
              for i=1:n
                 [Di,msg0] = dtfm ( F(i,:), T, t, type0, H0, recFlag0 );
                 if ~isempty(msg0), D = zpk(D); type0 = 'zpk'; msg = msg0; end;
                 D(i,:) = Di;
              end;              
          %------------------------------------------------------
          %     Single row case  
          %------------------------------------------------------         
           else         
              for j=1:m
               %------------------------------------------------------
               %     Scalar extrapolator: D(i,j) = D{F(i,j)*H}
               %------------------------------------------------------         
                if scalarH
                   [Dj,msg0] = dtfm ( F(1,j), T, t, type0, H0, recFlag0 );
                   if ~isempty(msg0), msg = msg0; end;      
               %------------------------------------------------------
               %     Matrix extrapolator: D(i,j) = sum_1^ell D{F(i,k)*H(k,j)}
               %------------------------------------------------------         
                else   
                   switch type0
                     case 'ss',  Dj = ss(0); 
                     case 'zpk', Dj = zpk(0); 
                   end;
                   Dj.Ts = T;
                   for k=1:ell
                     [Djk,msg0] = dtfm ( F(1,k), T, t, type0, H(k,j), recFlag0 );
                     if ~isempty(msg0), Dj = zpk(Dj); type0 = 'zpk'; msg = msg0; end;
                     Djk = delay2z(Djk);
                     if isequal(class(Dj),'ss')
                          Dj = minreal(Dj + Djk, sqrt(eps), 0); 
                     else Dj = minreal(sumzpk(Dj,Djk)); 
                     end;
                   end;    
                end;
                if ~isempty(msg0), D = zpk(D); type0 = 'zpk'; end;
                D(1,j) = Dj;
              end;
           end;          
          %------------------------------------------------------
          %     Final transformations
          %------------------------------------------------------          
           switch ( type )
             case 'zpk', D = zpk(D);  
             case 'tf',  D = tf(D);  
           end;
           if ~recFlag, 
              D = delay2z(D); 
              if isequal(class(D), 'ss')
                 D = minreal(D, sqrt(eps), 0); 
              end;
           end;                   
           return; 
        end;
        
%####################################################################        
%                      BASIC ALGORITHM                              
%####################################################################        
%------------------------------------------------------
%       Handle time
%------------------------------------------------------
        t = t - comdelay;
        tInt = floor ( t/T );
        m = t - tInt*T;       
        if tInt > 0, 
           zNum = tInt;
           zDen = 0;
        else % to avoid problems with static gain !
           zNum = 0;
           zDen = - tInt;
        end;
%------------------------------------------------------
%       Build matrices of discrete state space model 
%------------------------------------------------------
        if isstatic(F)
           Ad = []; Bd = []; Cd = [];
           Dd = D*K*expm(L*t)*M;
        else
           [Gammam,Phi1m,Phi2m] = intaba ( A, B*K, L, m);
           [Gamma,Phi]   = intaba ( A, B*K, L, T);
           Ad = Phi;
           Cd = C*Phi1m;
           Bd = Gamma*M;
           Dd = (C*Gammam + D*K*Phi2m)*M;
        end;
        D = ss ( Ad, Bd, Cd, Dd, T );       
%------------------------------------------------------
%       Use additional shift z^k or 1/z^k
%------------------------------------------------------
        msg = '';
        if zNum > 0
           add = zpk( zeros(1,zNum), [], 1 ); 
           if isequal(type,'tf'), D = tf(D); add = tf(add); end;
           if isequal(type,'ss') || isequal(type,'zpk'), D = zpk(D);    end;
           add.Ts = T;
           D = minreal(D * add);        
           if isequal(type, 'ss') &&  isa(D, 'zpk')
              if isproper(D)
                   D = ss(D); 
              else msg = 'Improper functions cannot be converted to state space'; 
              end;
              if nargout < 2, warning(msg); end;
           end;
        else
           switch ( type )
             case 'zpk', D = zpk(D);  
             case 'tf',  D = tf(D);  
           end;
           if zDen > 0, D.ioDelay = zDen * ones(size(D)); end;        
        end;
%------------------------------------------------------
%       Transform discrete delay to state space
%------------------------------------------------------
        if ~recFlag, 
           D = delay2z(D); 
           if isequal(class(D), 'ss')
                D = minreal(D, sqrt(eps), 0); 
           else D = minreal(D, sqrt(eps)); 
           end;
        end;        
        
%------- End of DTFM.M --------- KYuP ----------