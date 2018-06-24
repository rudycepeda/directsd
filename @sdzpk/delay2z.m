function sys = delay2z(sys)
%DELAY2Z  Replaces delays by poles at z=0 or FRD phase shift.  
%
%   For discrete-time TF, ZPK, or SS models SYS,
%      SYSND = DELAY2Z(SYS) 
%   maps all time delays to poles at z=0.  Specifically, a 
%   delay of k sampling periods is replaced by (1/z)^k.
%
%   For state-space models,
%      [SYSND,G] = DELAY2Z(SYS)
%   also returns the matrix G mapping the initial state x0
%   of SYS to the corresponding initial state G*x0 for SYSND.
%   
%   For FRD models, DELAY2Z absorbs all time delays into the 
%   frequency response data, and is applicable to both 
%   continuous- and discrete-time FRDs.
%
%   See also HASDELAY, PADE, LTIMODELS.

%   Author(s):  P. Gahinet 
%   Copyright 1986-2002 The MathWorks, Inc. 
%   $Revision: 1.9 $  $Date: 07-Apr-2006 $

%   ### Revised by K. Polyakov, 13 Dec 2003
%   ### Allow for negative delays                          

narginchk(1,1);

Td = totaldelay(sys);
if sys.Ts==0
   error('LTI model SYS must be discrete.')
elseif ~any(Td(:)),
   return
end

% Map the I/O delays to poles at zero
sizes = size(sys.K);
if ndims(Td)<length(sizes), 
   Td = repmat(Td,[1 1 sizes(3:end)]);
end

for j=find(Td(:))',
   Tdj = Td(j); 
   if Tdj > 0					% Revision by K. Polyakov
        sys.P{j} = [sys.P{j} ; zeros(Tdj,1)];
   else sys.Z{j} = [sys.Z{j} ; zeros(-Tdj,1)]; 
   end;
end

% Set I/O delays to zero
sys = pvset(sys,'InputDelay',zeros(sizes(2),1),...
   'OutputDelay',zeros(sizes(1),1),...
   'ioDelay',zeros(sizes(1:2)));

sys = minreal(sys);
