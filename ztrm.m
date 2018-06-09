function [DF,msg] = ztrm ( F, T, t, type, recFlag )
%ZTRM Modified Z-transform for a transfer matrix.
%
%     D = ZTRM ( F, T, t, TYPE )
%
%   Inputs:
%     F(s) - transfer matrix of continuous-time plant
%     T    - sampling period
%     t    - time parameter
%     TYPE - class of return value: 'zpk' (default), 'ss', or 'tf'
%
%   Outputs::
%     D(z) - modified Z-transform in the zpk-form       
%
%   See also DTFM, DTFM2.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 04-Dec-2007 $
%------------------------------------------------------
%       Check syntax
%------------------------------------------------------
        if ~exist('recFlag','var'), recFlag = 0; end;
        if ~exist('type','var'), type = 'zpk'; end;
        if ~exist('t','var')  |  isempty(t), t = 0; end;
        type = lower(type);
        if ~isequal(type,'zpk')  &&  ~isequal(type,'ss')  &&  ~isequal(type,'tf')
           error('Unknown type ''%s''',type); 
        end;       
        [rows,cols] = size(F);
%####################################################################        
%       Time-delayed case: if delays are different, use
%       row-wise or even element-wise transformations
%####################################################################        
        delay = totaldelay ( F );
        comdelay = min(min(delay));
        delay = delay - comdelay;
        
        if norm(delay) > sqrt(eps)
           switch type
             case 'ss',         DF = ss(zeros(n,m));  type0 = 'ss';
             case {'zpk','tf'}, DF = zpk(zeros(n,m)); type0 = 'zpk';
           end;
           DF.Ts = T;
           msg = '';
          %------------------------------------------------------
          %     Multirow case  
          %------------------------------------------------------         
           recFlag0 = 1;
           if rows > 1
              for i=1:rows
                 [Di,msg0] = ztrm ( F(i,:), T, t, type0, recFlag0 );
                 if ~isempty(msg0), DF = zpk(DF); type0 = 'zpk'; msg = msg0; end;
                 DF(i,:) = Di;
              end;              
          %------------------------------------------------------
          %     Single row case  
          %------------------------------------------------------         
           else         
              for j=1:cols
                [Dj,msg0] = ztrm ( F(1,j), T, t, type0, recFlag0 );
                if ~isempty(msg0), msg = msg0; end;      
                if ~isempty(msg0), DF = zpk(DF); type0 = 'zpk'; end;
                DF(1,j) = Dj;
              end;
           end;          
          %------------------------------------------------------
          %     Final transformations
          %------------------------------------------------------          
           switch ( type )
             case 'zpk', DF = zpk(DF);  
             case 'tf',  DF = tf(DF);  
           end;
           if ~recFlag, 
              DF = delay2z(DF); 
              if isequal(class(DF), 'ss')
                 DF = minreal(DF, sqrt(eps), 0); 
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
        t = t - tInt*T;
        add = zpk( zeros(1,abs(tInt)), [], 1 ); 
        if tInt > 0, 
           zNum = tInt;
           zDen = 0;
           tInt = 0; 
        else % to avoid problems with static gain !
           zNum = 0;
           zDen = - tInt;
           tInt = 0; 
        end;
%------------------------------------------------------
%       Special case: discrete white noise
%------------------------------------------------------
        ssF = minreal( ss(F), sqrt(eps), 0 );
        if size(ssF.a) == 0
           DF = zpk(ssF.d);
           if tInt > 0, DF = minreal(DF*zAdd); end;
           if tInt < 0, DF = minreal(DF/zAdd); end;
           switch type 
             case 'ss', if tInt <= 0, DF = ss(DF); end;
             case 'tf', DF = tf(DF);
           end;
           DF.Ts = T;
           return;
        end;
%------------------------------------------------------
%       Build discrete Laplace transform
%------------------------------------------------------
        if norm(ssF.d) > eps, 
           % Is admissible, the case of discrete white noise!!! 
           %error('Function is apllicable to strictly proper Laplace images only'); 
        end;
        Phi = expm(ssF.a*T);
        DF = ss(ssF.a, ssF.b, ssF.c, ssF.d);
        DF.a = Phi;
        DF.b = Phi*ssF.b;
        DF.c = ssF.c*expm(ssF.a*t);
        DF.d = DF.c*ssF.b + ssF.d;
        DF.Ts = T;
%------------------------------------------------------
%       Mandatory zero at z=0 for t > 0
%------------------------------------------------------
        zpkDF = zpk(DF);
        for r=1:rows 
          for c=1:cols 
            z = zpkDF.z{r,c};
            if length(z) > 0
              [xx,ind] = sort(abs(z));
              z(ind(1)) = 0; 
              zpkDF.z{r,c} = z;
            end;
          end;
        end;
        if norm(zpkDF.k) < eps,
           DF = zpk ( zeros(size(F)) );
           DF.Ts = T;
        end;       
%------------------------------------------------------
%       Use additional shift z^k or 1/z^k
%------------------------------------------------------
        msg = '';
        if zNum > 0
           if isequal(type,'tf'), DF = tf(DF); add = tf(add); end;
           if isequal(type,'ss') || isequal(type,'zpk'), DF = zpk(DF);    end;
           add.Ts = T;
           DF = minreal(DF * add);        
           if isequal(type, 'ss') &&  isa(DF, 'zpk')
              if isproper(DF)
                   DF = ss(DF); 
              else msg = 'Improper functions cannot be converted to state space'; 
              end;
              if nargout < 2, warning(msg); end;
           end;
        else
           switch ( type )
             case 'zpk', DF = sdzpk(DF);  
             case 'tf',  DF = sdtf(DF);  
           end;
           if zDen > 0, DF.ioDelay = zDen * ones(size(DF)); end;        
        end;
%------------------------------------------------------
%       Transform discrete delay to state space
%------------------------------------------------------
        if ~recFlag, 
           DF = delay2z(DF); 
           if isequal(class(DF), 'ss')
                DF = minreal(DF, sqrt(eps), 0); 
           else
               DF = minreal(DF, sqrt(eps));
           end
        end      
        
%------- End of ZTRM.M --------- KYuP ----------
