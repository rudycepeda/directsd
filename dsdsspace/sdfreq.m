function  [R,w] = sdfreq ( csys, K, w, type )
%SDFREQ Frequency response gain for sampled-adta system.
%
%     R = SDFREQ ( CSYS, K, W, TYPE )
%
%   Inputs:
%     CSYS - continuous-time LTI system
%     K    - discrete-time controller
%     TYPE - type of frequency responce of
%         'std'  - standard (default)
%         'conj' - Hermitian conjugate
%         'spec' - spectral (for W'(-s,t)*W(s,t))
%         'sing' - sqrt(eigenvalues) for (1/T) int_0^T W'(-jw,t)*W(jw,t) dt
%         'dspec' - discretized spectrum
%         'dsing' - discretized singular values
%   Outputs:
%     R(:,:,i) - average frequency response for W(i)
%          R = (1/T) int_{0}^T W(s,t) dt =
%            = P_{11} + (1/T)P_{12}*H*M*P_{21},
%          where M = K*(I - D_{P_{22}*H}(0)K)^{-1} 
%
%   See also SDHINORM.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%	Check syntax
%------------------------------------------------------
    if nargin == 3  &&  ischar(w),
       type = w;
       w = [];
    end;
    if ~exist('type','var'), type = 'std'; end;
    if ~isequal(type,'std') && ~isequal(type,'conj') && ~isequal(type,'spec') ...
       && ~isequal(type,'dspec')  && ~isequal(type,'sing') && ~isequal(type,'dsing')   
       error('Unknown type of frequency responce ''%s''',type); 
    end;
    if ~isct(csys), error('Plant must be a continuous-time system'); end;
    if ~isdt(K),    error('Controller must be a discrete-time system'); end;
    [nout,nin] = size(csys);
    [i2,o2] = size(K);
    i1 = nin - i2;
	o1 = nout - o2;
    if i1 < 1, error('No input signals detected'); end;
    if o1 < 1, error('No output signals detected'); end;
    T = K.Ts;
    if ~exist('w','var')  ||  isempty(w)
       ws = 2*pi / T; 
       w = linspace(0.001*ws, ws*0.999, 50);
    end;    
    n = length(w);
%------------------------------------------------------
%	Extract matrices
%------------------------------------------------------
    P11 = csys(1:o1,1:i1);
    P12 = csys(1:o1,i1+1:end);
    P21 = csys(o1+1:end,1:i1);
    P22 = csys(o1+1:end,i1+1:end);
%------------------------------------------------------
%   Construct M = K*(I - D22*K)^{-1} 
%------------------------------------------------------
    D22 = c2d ( P22, T );
    M = feedback ( K, - D22 );
    rM = freqresp(M, w);
%------------------------------------------------------
%   Frequency responces of continuous parts
%------------------------------------------------------
    j = sqrt(-1);
    rH = (1-exp(-j*w*T)) ./ (j*w);
    rH = reshape(rH, 1, 1, length(rH));
    rP11 = freqresp(P11, w);
    rP12 = freqresp(P12, w);
    rP21 = freqresp(P21, w);
%------------------------------------------------------
%   Standard averaged frequency responce
%------------------------------------------------------
    RPMP = mulfreq(rP12,mulfreq(rH,mulfreq(rM,rP21))) / T;
    if isequal(type,'std')
       R = rP11 + RPMP;    
       return;
    end;
%------------------------------------------------------
%   Complex conjugate transposed averaged frequency responce
%------------------------------------------------------
    if isequal(type,'conj')
       R = conjfreq(R);
       return;
    end;
%------------------------------------------------------
%   Spectrum of averaged frequency responce
%------------------------------------------------------
    if isequal(type,'spec')  ||  isequal(type,'sing')
       SD12 = dtf2 ( P12, T, 0 );       
       rSD12 = freqresp(SD12, w) / T;
       rP11conj = conjfreq(rP11);
       R11 = mulfreq(rP11conj,rP11);
       R21 = mulfreq(rM, rP21);
       RA = mulfreq(conjfreq(R21), mulfreq(rSD12, R21));
       RB = mulfreq(rP11conj, RPMP);
       R = real(R11 + RA + RB + conjfreq(RB));
       if isequal(type,'sing')
          S = zeros(1,n); 
          for i=1:n
             S(i) = sqrt(max(svd(R(:,:,i)))); 
          end;
          R = S;
       end;
    end;      
%------------------------------------------------------
%   Discretized spectrum of averaged frequency responce
%------------------------------------------------------
    if isequal(type,'dspec')  ||  isequal(type,'dsing')
       SD12 = dtf2 ( P12, T, 0 ) / T;       
       rSD12 = freqresp(SD12, w);
       
       SD11 = ztr(P11'*P11, T, 0);
       R11 = freqresp(SD11, w);
       
       RMAM = mulfreq(conjfreq(rM), mulfreq(rSD12,rM));
       SD21 = ztr(P21*P21', T, 0);
       rSD21 = freqresp(SD21, w);
       RA = mulfreq(RMAM, rSD21);
       
       B = dtfn (P21*P11'*P12, T, 0) / T;
       rB = freqresp(B, w);
       RBM = mulfreq(rB, rM);
       
       R = T * real(trfreq(R11) + RA + RBM + conjfreq(RBM));
       if isequal(type,'dsing')
          S = zeros(1,n); 
          for i=1:n
             S(i) = sqrt(max(svd(R(:,:,i)))); 
          end;
          R = S;
       end;
    end;
       
%------- End of SDFREQ.M --------- KYuP ----------

%###########################################################
function trR = trfreq ( R )
%TRFREQ Trace of a matrix frequency responce.
%
%     T = TRFREQ ( R )
%
%------------------------------------------------------
%       K. Yu. Polyakov         13 Feb 2005
%                               
%------------------------------------------------------
       n = size(R, 3);
       trR = zeros(1,1,n);
       for i=1:n, 
         trR(:,:,i) = trace(R(:,:,i)); 
       end;

%------- End of TRFREQ.M --------- KYuP ----------

%###########################################################
function Rconj = conjfreq ( R )
%CONJFREQ Complex conjugate transpose for a matrix frequency responces.
%
%     RCONJ = CONJFREQ ( R )
%
%------------------------------------------------------
%       K. Yu. Polyakov         11 Feb 2005
%                               
%------------------------------------------------------
       [r,c,n] = size(R);
       Rconj = zeros(c,r,n);
       for i=1:n, 
         Rconj(:,:,i) = R(:,:,i)'; 
       end;

%------- End of CONJFREQ.M --------- KYuP ----------

%###########################################################
function  R = mulfreq ( R1, R2 )
%MULFREQ Multiply matrix frequency responces.
%
%     R = MULFREQ ( R1, R2 )
%
%------------------------------------------------------
%       K. Yu. Polyakov         11 Feb 2005
%                               
%------------------------------------------------------
    N = size(R1,3);
    R = zeros(size(R1,1), size(R2,2), N);    
    for i=1:N
      R(:,:,i) = R1(:,:,i)*R2(:,:,i);  
    end;

%------- End of MULFREQ.M --------- KYuP ----------
