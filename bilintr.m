function B = bilintr ( F, type, T )
%BILINTR Bilinear transformation for SISO system.
%
%     B = BILINTR ( F, TYPE );
%     B = BILINTR ( F, TYPE, T );
%     B = BILINTR ( F, COEFS );
%       
%   Inputs:  
%     F    - an LTI system model
%     TYPE - transformation type:
%        's2z' from 's' to 'z', maps LHP into the unit disk
%        'z2s' from 'z' to 's' (default), maps the unit disk into LHP
%        's2d' from 's' to d=1/z, maps RHP into the unit disk
%        'd2s' from 'd' to 's', maps the unit disk into RHP
%        'tustin' Tustin (and inverse Tustin) transform
%     COEFS - coefficients of a nonstandard bilinear transformation 
%             as an array [a, b, c, d] so that
%                y = (a*x + b)/(c*x + d)    
%             If only one coefficient is given, it is assumed that 
%             it is 'b', and y = (x + b)/(x - b)
%     T - sampling period for Tustin transform
%
%   Outputs: 
%     B - transformed LTI system
%
%   See also BILINSS.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if ~isobject(F)
           error('First argument must be an LTI system');
        end;        
        [rows,cols] = size(F);
        if ~exist('type','var') 
           if isct(F)  
                type = 's2z'; 
           else type = 'z2s'; 
           end;
        end;
        if ischar(type), type = lower(type); end;
        if ~exist('T','var'), T = 1; end;
%------------------------------------------------------
%       MIMO case - termwise
%------------------------------------------------------
        if rows ~= 1  ||  cols ~= 1
           B = zpk(zeros(rows,cols));
           if ~isdt(F), B.Ts = 1; end;
           if isa(F, 'tf'), B = tf(B); end;
           if isa(F, 'ss'), B = ss(B); end;
           for r=1:rows
             for c=1:cols
               B(r,c) = bilintr(F(r,c), type, T);  
             end;
           end;
           return;
        end;
%------------------------------------------------------
%       Standard coefficients 
%------------------------------------------------------
        a0 = 1;   b0 = 1;      % z = (1 + s)/(1 - s)
        c0 = - 1; d0 = 1;      % s = (z - 1)/(z + 1)
                               % s = (-d + 1)/(d + 1)
                               % d = (1 - s)/(1 + s)
%------------------------------------------------------
%       Special bilinear transformation 
%------------------------------------------------------
        if isnumeric(type)
           if length(type) == 1
              type = [1 type 1 -type];
           end;
           if length(type) ~= 4
              error('Incorrect number of coefficients');
           end;
           a0 = type(1);
           b0 = type(2);
           c0 = type(3);
           d0 = type(4);
           type = 'special';
        end;
%------------------------------------------------------
%       Select transformation type 
%------------------------------------------------------
%        z2s: z = (as + b)/(cs + d)
%        s2z: s = (-dz + b)/(cz - a)
%        d2s: d = (cs + d)/(as + b)
%        s2d: s = (bz - d)/(-az + c)
%------------------------------------------------------
        if strcmp(type,'z2s')
             a = a0; b = b0; c = c0; d = d0;
        elseif strcmp(type,'s2z')
             a = -d0; b = b0; c = c0; d = -a0;
        elseif strcmp(type,'d2s')
             a = c0; b = d0; c = a0; d = b0;
        elseif strcmp(type,'s2d')
             a = b0; b = -d0; c = -a0; d = c0;
        elseif strcmp(type,'special')
             a = a0; b = b0; c = c0; d = d0;
        elseif strcmp(type,'tustin')
          mu = 2/T;  
          if isct(F)
               a = mu; b = -mu; c = 1; d = 1;  
          else a = 1; b = mu; c = -1; d = mu; 
          end;
        else error('Incorrect bilinear transformation type');
        end;
%------------------------------------------------------
%       Perform bilinear transformation 
%------------------------------------------------------
        [z,p,K,Ts] = zpkdata(zpk(F), 'v');
        nz = length(z);
        np = length(p);
       %-----------------------------------------------
       %       Zeros
       %-----------------------------------------------
        for i=1:nz
           if abs(a - z(i)*c) > eps
              K = K * (a/c - z(i));
              z(i) = (z(i)*d - b) / (a - z(i)*c);
           else
              K = K * (b - z(i)*d) / c;
              z(i) = Inf;
           end;
        end;
       %-----------------------------------------------
       %       Poles
       %-----------------------------------------------
        for i=1:np
           if abs(a - p(i)*c) > eps
              K = K / (a/c - p(i)) ;
              p(i) = (p(i)*d - b) / (a - p(i)*c);
           else
              K = K / ((b - p(i)*d) / c);
              p(i) = Inf;
           end;
        end;
       %-----------------------------------------------
       %       Check poles and zeros at infinity
       %-----------------------------------------------
        z(z == inf) = [];
        p(p == inf) = [];        
       %-----------------------------------------------
       %  If number of zeros and poles are not the same
       %-----------------------------------------------
        if nz > np
          p = [p; ones(nz-np,1)*(-d/c)];
        else
          z = [z; ones(np-nz,1)*(-d/c)];
        end;                
%------------------------------------------------------
%       Finalization
%------------------------------------------------------
        K = real(K);
        if Ts == 0
             B = zpk(z, p, K, T);
        else B = zpk(z, p, K); 
        end;
        if isa(F, 'tf'), B = tf(B); end;
        if isa(F, 'ss'), B = ss(B); end;
        
%------- End of BILINTR.M --------- KYuP ----------           
