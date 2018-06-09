function tsys = ctranspose(sys)
%CTRANSPOSE  Pertransposition of zero-pole-gain models.
%
%   TSYS = CTRANSPOSE(SYS) is invoked by TSYS = SYS'
%
%   If SYS represents the continuous-time transfer function
%   H(s), TSYS represents its pertranspose H(-s).' .   In
%   discrete time, TSYS represents H(1/z).' if SYS represents
%   H(z).
%
%   See also TRANSPOSE, ZPK, LTIMODELS.

%   Author(s): A. Potvin, P. Gahinet
%   Copyright 1986-2002 The MathWorks, Inc.
%   $Revision: 1.15 $  $Date: 07-Apr-2006 $

%   Revised by K. Polyakov, 10 Apr 2003
%                           05 May 2003

% Extract data
[z,p,k,Ts] = zpkdata(sys);
if isempty(k),
    tsys = sys.';
    return
end

sizes = size(k);
sizes([1 2]) = sizes([2 1]);

if Ts==0,
    % Continuous-time case: replace s by -s
    for j=1:prod(sizes),
        z{j} = -sort(conj(z{j}));
        p{j} = -sort(conj(p{j}));
        dl = length(p{j}) - length(z{j});
        if mod(dl,2),
            k(j) = -k(j);
        end
        k(j) = conj(k(j));
    end
else
    % Discrete-time case: replace z by z^-1
    for j=1:prod(sizes),
        zj = sort(conj(z{j}));
        pj = sort(conj(p{j}));
        
        %-------------------------------------------
        % Original version does not work because
        % NOT operation cannot be applied to complex
        % values starting from version 6.5
        %	idz = find(~zj);   zj(idz) = [];
        %	idp = find(~pj);   pj(idp) = [];
        % Changed by K. Polyakov, 10 Apr 2003
        %-------------------------------------------
        
        idz = find(~abs(zj));   zj(idz) = [];
        idp = find(~abs(pj));   pj(idp) = [];
        
        %-------------------------------------------
        
        k(j) = conj(k(j)) * prod(-zj) / prod(-pj);
        zj = 1./zj;
        pj = 1./pj;
        zpow = length(idp) + length(pj) - (length(idz) + length(zj));
        z{j} = [zj ; zeros(zpow,1)];
        p{j} = [pj ; zeros(-zpow,1)];
    end
end

% Create result
% tsys = sys;
% tsys.z = permute(z,[2 1 3:length(sizes)]);
% tsys.p = permute(p,[2 1 3:length(sizes)]);
% tsys.k = permute(k,[2 1 3:length(sizes)]);

%Modified by RCG on 29.07.2014

zz = permute(z,[2 1 3:length(sizes)]);
pp = permute(p,[2 1 3:length(sizes)]);
kk = permute(k,[2 1 3:length(sizes)]);

args ={'DisplayFormat'    sys.DisplayFormat...
    'Variable'         sys.Variable...
    'ioDelay'          sys.ioDelay'...
    'InputDelay'       sys.OutputDelay...
    'OutputDelay'      sys.InputDelay...
    'TimeUnit'         sys.TimeUnit...
    'InputName'        sys.InputName...
    'InputUnit'        sys.InputUnit...
    'InputGroup'       sys.InputGroup...
    'OutputName'       sys.OutputName...
    'OutputUnit'       sys.OutputUnit...
    'OutputGroup'      sys.OutputGroup...
    'Name'             sys.Name...
    'Notes'            sys.Notes...
    'UserData'         sys.UserData};



%-------------------------------------------
% Original version does not work because of
%  a misprint
%	sys.lti = (sys.lti)';
% Changed by K. Polyakov, 04 May 2003
%-------------------------------------------

%tsys.lti = (sys.lti)';

%-------------------------------------------
% Original version does could return complex
% coefficient
% Changed by K. Polyakov, 04 May 2003
%-------------------------------------------

%tsys.k = real(tsys.k);

kk = real(kk);
tsys = sdzpk(zz,pp,kk,sys.Ts,args{:});
