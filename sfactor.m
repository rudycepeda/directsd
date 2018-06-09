function [fs,fs0] = sfactor ( s, type )
%SFACTOR Spectral factorization for polynomials and transfer functions.
%        Neutral zeros and poles are divided equally between stable
%        and unstable parts.
%
%    [FS,FS0] = SFACTOR ( S, TYPE )
%
%  Inputs:
%    S - a continuous-time or discrete LTI system 
%    TYPE - factorization type of 
%        's' (default for continuous-time models)
%        'd' (default for discrete-time models)
%        'z' for z-plane
%  Outputs:
%    FS  - minimal realization of spectral factor
%    FS0 - without cancellations
%
%   See also FACTOR, SFACTFFT.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 29-Oct-2006 $
%------------------------------------------------------    
%       Check data
%------------------------------------------------------           
        if ~exist('type','var'),
           if isnumeric(s), s = poln(s); end; 
           if isdt(s), type = 'd'; else type = 's'; end;
        end;
        if ~strcmp(type,'s') &&  ~strcmp(type,'z')  &&  ~strcmp(type,'d')
            error('Unknown factorization type: ''%s''', type);
        end;         
        if isnumeric(s), s = poln(s, type); end;
        if isa(s,'poln')
           fs = sfactor ( s, type );
           fs = fs.coef;
           fs0 = fs;
           return;
        end;
        if ~isa(s,'lti')
            error('This factorization is applicable only to LTI models');
        end;         
%------------------------------------------------------
%       Check for zero function
%------------------------------------------------------
        s = zpk(s); 
        if s.K == 0, 
          fs = zpk(0); fs.Ts = s.Ts;
          fs0 = fs; 
          return;
        end;
        errmsg = 'Exact Hermitian factorization is impossible';
%------------------------------------------------------
%       Symmetrize zeros
%------------------------------------------------------
        tol = 1e-2;   % tolerance 1%
        zz = s.z{1}; zCorr = [];
        
        i = 1;          
        while i < length(zz)
         %----------------------------------------------- 
         % Find the nearest to counterpart of z0
         %----------------------------------------------- 
          z0 = zz(i); zz(i) = [];  
          if type == 's', 
            z0H = - z0; 
          elseif abs(z0) > sqrt(eps), 
            z0H = 1/z0; 
          else
            zCorr = [zCorr; z0];           
            continue;  
          end;
          diff = abs(zz - z0H); 
          [diff,ind] = sort(diff);  
          if diff(1) < tol*(abs(z0)+abs(z0H))/2  ||  diff(1) < sqrt(eps)
            z1 = zz(ind(1));
            im0 = imag(z0); im1 = imag(z1);
           %----------------------------------------------- 
           % Consider real-complex pair, make both real
           %----------------------------------------------- 
            if (im0 == 0) ~= (im1 == 0),
              z0 = real(z0); z1 = real(z1);  
              err = abs(im0+im1);  
              if err > 1e-3,
                warning('Symmetrization error %g',err);
              end;
            end;
           %----------------------------------------------- 
           % Symmetrization
           %----------------------------------------------- 
            if type == 's', 
              z0 = (z0 - z1) / 2;  
              zCorr = [zCorr; z0; -z0];
            else
              z0 = (1/z0 + z1) / 2;
              zCorr = [zCorr; z0; 1/z0];
            end;
            zz(ind(1)) = [];
          else
            zCorr = [zCorr; z0];  %error(errmsg); 
          end;
        end;        
%------------------------------------------------------
%       Unfactored zeros
%------------------------------------------------------
        if ~isempty(zz)
           zCorr = [zCorr; zz];
        end;
%------------------------------------------------------
%       Factorize a rational function
%------------------------------------------------------
        [zs,zRem,z0] = extrpair ( zCorr, type );
        [ps,pRem,p0] = extrpair ( s.p{1}, type );       
        if ~isempty(zRem)  ||  ~isempty(pRem), 
          if length(zRem) == length(pRem)
            rem = others(zRem,pRem,1e-2);
            if ~isempty(rem), error(errmsg); end;
          else warning(errmsg); %Changed error to warning RCG 31.07.2014 
          end;
        end;
%------------------------------------------------------
%       Check conjugation of complex zeros
%------------------------------------------------------
        ind = find(imag(zs) ~= 0);
        for i = 1:length(ind);  
          zc = zs(ind(i));  
          if isempty(find(zs == conj(zc)))   
            zs(ind(i)) = real(zc);   
          end;
        end;
%------------------------------------------------------
%       Check gain
%------------------------------------------------------
        if type == 's'
           K = real(s.K);
           if mod(length(ps)+length(zs),2) == 1, K = - K; end;
        else
           if length(zs)+z0 ~= length(ps)+p0, warning(errmsg); end;       %Changed error to warning RCG 31.07.2014
           K = real(s.K*prod(-ps) / prod(-zs));
        end;
        errmsg = 'Function is negative definite';
        if K < 0 
            warning(errmsg); %Changed by RCG 28.07.2013
            %return; 
        end        
%------------------------------------------------------
%       Form results
%------------------------------------------------------
        fs0 = zpk ( zs, ps, sqrt(K), s.Ts );
        fs = minreal(fs0);           
        
%------- End of SFACTOR.M --------- KYuP ----------
