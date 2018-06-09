function [Fs,Fu,F0] = separtf ( F, type, stype )
%SEPARTF Proper separation (polynomial equations technique).
%            F = FS + FU + F0
%
%     [FS,FU,F0] = SEPARTF ( F, TYPE, STYPE )
%
%   Inputs:   
%     F    - an LTI system model
%     TYPE - type of separation
%         's'  - stable poles are in the left half-plane  
%                (default for continuous-time systems) 
%         'z'  - stable poles are inside the unit disk
%                (default for discrete-time systems) 
%         'd'  - stable poles are outside the unit disk
%         'T'  - stable poles are inside the delta-disk 
%                with center at (-1/T,0) with radius 1/T
%         's0', 'z0', 'd0' - include neutral poles in the
%                stable part
%     STYPE - type of stable part
%         'infs' - stable part include polynomial (default) 
%         'infu' - unstable part include polynomial
%
%   Outputs:  
%     FS - stable term
%     FU - strictly proper strictly antistable term
%     F0 - neutral term
%
%   See also SFACTOR, SEPARSS.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if ~isobject(F)
           error('First argument must be an LTI system');
        end;
%------------------------------------------------------
        if ~exist('type', 'var'),
           if isdt(F),
                type = 'z'; 
           else type = 's'; 
           end;
        end;        
        if ~strcmp(type,'s')  &&  ~strcmp(type,'z')  &&  ~strcmp(type,'d') && ...
           ~strcmp(type,'s0')  &&  ~strcmp(type,'z0')  &&  ~strcmp(type,'d0')  &&  ~strcmp(type,'T')
           error('Incorrect separation type');
        end;        
%------------------------------------------------------
        if ~exist('stype', 'var'),
           stype = 'infs'; 
        end;
        if ~strcmp(stype,'infs')  &&  ~strcmp(stype,'infu')
           error('Incorrect properness type');
        end;
%------------------------------------------------------
        alpha = 0;
        dalpha = exp(-alpha*F.Ts);       
%------------------------------------------------------
%       Find stable poles
%------------------------------------------------------
        [z,p,K,T] = zpkdata(zpk(F), 'v');        
        p0 = [];
        switch type
          case 's',
             p0 = find(real(p) == - alpha);   
             ps = find(real(p) < - alpha);
          case 'z',
             p0 = find(abs(p) == dalpha);
             ps = find(abs(p) < dalpha);
          case 'd',
             p0 = find(abs(p) == dalpha);
             ps = find(abs(p) > dalpha);
          case 'T',
             T = F.Ts;
             if (F.var ~= 'q')  ||  (T <= 0) 
                error('Transfer function must be in delta-domain with Ts > 0'); 
             end;
             p0 = find(abs(p+1/T) == 1/T);   
             ps = find(abs(p+1/T) < 1/T);   
          case 's0',
             ps = find(real(p) <= - alpha);
          case 'z0',
             ps = find(abs(p) <= dalpha);
          case 'd0',
             ps = find(abs(p) >= dalpha);
        end;        
%------------------------------------------------------
%       Separation via polynomial equation
%------------------------------------------------------
        ps = p(ps);
        pu = others(p, ps);
        ds = poly(ps);
        du = poly(pu);
        n = K*poly(z);       
        if isequal(stype, 'infs')
             [nu,ns] = dioph (ds, du, n);
        else [ns,nu] = dioph (du, ds, n); 
        end;        
        Fs = minreal ( zpk ( roots(ns), ps, ns(1), T ) );
        Fu = minreal ( zpk ( roots(nu), pu, nu(1), T ) );
%------------------------------------------------------
%       Handle neutral poles
%------------------------------------------------------
        if ~isempty(p0)
          %error('Function has poles at the stability boundary');
          p0 = p(p0);
          pu = others(pu, p0);
          d0 = poly(p0);
          du = poly(pu);
          n = tfdata ( Fu, 'v');
          if isequal(stype, 'infs')
               [nu,n0] = dioph (d0, du, n);
          else [n0,nu] = dioph (du, d0, n); 
          end;        
          F0 = minreal ( zpk ( roots(n0), p0, n0(1), T ) );
          Fu = minreal ( zpk ( roots(nu), pu, nu(1), T ) );          
        else
          F0 = zpk([], [], 0, T);  
        end;       
%------------------------------------------------------
%       Tranform to required form
%------------------------------------------------------
        if isa(F, 'tf'), 
           Fs = tf(Fs);
           Fu = tf(Fu);	  
           F0 = tf(F0);	  
        end;       
        if isequal(type,'T')  &&  (isa(F, 'tf')  ||  isa(F, 'zpk'))  
           Fs.var = 'q';
           Fu.var = 'q';	  
           F0.var = 'q';	  
        end;
        if isa(F, 'ss'),
           Fs = ss(Fs);
           Fu = ss(Fu);
           F0 = ss(F0);
        end;
        
%------- End of SEPARTF.M --------- KYuP ----------
