function [fp,fm] = sfactfft ( p0, type, RN )
%SFACTFFT Polynomial spectral factorization using FFT.
%
%     [FS,FU] = SFACTFFT ( P, TYPE, N )
%
%   Inputs:
%     P    - initial symmetric quesipolynomial
%     TYPE - factorization type of 'd' (default) or 'z' 
%     N    - number of points is N*deg(P) (default N=10)
%
%   Outputs:
%     FS - stable spectral factor
%     FU - unstable spectral factor
%
%   See also FACTOR, SFACTOR. 

% Reference:
% [1] Hromcik M., Jezek J., and Sebek M., New algorithm for spectral
%     factorization and its practical application //  in Proc. European
%     Control Conference ECC'2001, Porto, Portugal, September 1-5, 2001.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if ~exist('RN','var'), RN = 10; end;
        polnMode = isa(p0,'poln');
        if polnMode, 
          var = p0.var;
          if ~isequal(var,'z') && ~isequal(var,'d')  &&  ~isequal(var,'q')  
            error('FFT factorization not applicable to continuous-time polynomials');  
          end;
          p0 = p0.coef; 
        end;
        if ~exist('type','var'), 
          type = 'd'; 
        end;
        if ~isequal(type,'z')  &&  ~isequal(type,'d')  
           error('Unknown factorization type ''%s''', type); 
        end;
%------------------------------------------------------
%       Find R
%------------------------------------------------------
        dg = floor( (length(p0)-1)/2 );
        R = RN*length(p0);
%------------------------------------------------------
%       FFT-I
%------------------------------------------------------
        nZeros = 2*R+1 - length(p0);
        p = [fliplr(p0(1:dg+1)) zeros(1,nZeros) p0(1:dg)];
        P = fft(p);
%------------------------------------------------------
%       Logarithmization
%------------------------------------------------------
        N = log(P);
%------------------------------------------------------
%       IFFT-I
%       For symmetric spectral factorzation we need not
%       find fm specially
%------------------------------------------------------
        n = ifft(N);
        xp = n(1:R+1); xp(1) = n(1)/2;
        %xm = [n(1)/2 fliplr(n(R+2:end))];
%------------------------------------------------------
%       FFT-II
%------------------------------------------------------
        Xp = fft(xp);
        %Xm = fft(xm);
%------------------------------------------------------
%       Exponential for each part
%------------------------------------------------------
        Pvp = exp(Xp);
        %Pvm = exp(Xm);
%------------------------------------------------------
%       IFFT-II
%------------------------------------------------------
        pvp = ifft(Pvp);
        %pvm = ifft(Pvm);
%------------------------------------------------------
%       Extract results
%------------------------------------------------------
        fp = real(pvp(1:dg+1));
        %fm = pvm(1:dg+1);
        if type(1) == 'd', fp = fliplr(fp); end;
%------------------------------------------------------
%       Optimization
%------------------------------------------------------
        if length(type) > 1  &&  exist('fminunc','file')
           options = optimset ( 'display', 'off' );
           ws = warning; warning off;           
           fp = fminunc ( @f_sfact, fp, options, p0);
           warning(ws);        
        end;
        fm = fliplr(fp);
%------------------------------------------------------
%       Restore polynomial form
%------------------------------------------------------
        if polnMode, 
           fp = poln(fp, var);  
           fm = poln(fm, var);  
        end;
                
%------- End of SFACTFFT.M --------- KYuP ----------           

%######################################################
function err = f_sfact ( x, p2 )
%
% F_SFACT auxiliary function for refining spectral factorization
%

%------------------------------------------------------
%       K. Yu. Polyakov         28 Oct 2005
%                               
%------------------------------------------------------
        x = x(:)';
        p1 = conv(x, fliplr(x));
        e = sumpol2 ( p1, -p2 );
        err = norm(e);
 %------- End of F_SFACT.M --------- KYuP ----------           

