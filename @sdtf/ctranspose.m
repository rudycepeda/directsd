function L = ctranspose(L)
%CTRANSPOSE   Manages LTI properties in pertransposition.
%
%   SYS.LTI = (SYS.LTI)'

%       Author(s): P. Gahinet, 5-28-96
%       Copyright 1986-2002 The MathWorks, Inc. 
%       $Revision: 1.11 $

%	### Revised by K. Polyakov, 13.12.2003 ###
%   ### Allow for negative delays          ###

%if hasdelay(L),
%   error('Time delays makes pertransposed model non causal.')
%end

% Get I/O dims
nu = length(L.InputName);
ny = length(L.OutputName);

% Delete I/O names and groups
EmptyStr = {''};
L.InputName = EmptyStr(ones(ny,1),1);
L.OutputName = EmptyStr(ones(nu,1),1);
L.InputGroup = cell(0,2);
L.OutputGroup = cell(0,2);


    % Set delay times to zero
    %L.ioDelay = zeros(nu,ny);
    %L.InputDelay = zeros(ny,1);
    %L.OutputDelay = zeros(nu,1);

% Convert delay times
L.ioDelay = - L.ioDelay';
temp = L.InputDelay;
L.InputDelay = - L.OutputDelay;
L.OutputDelay = - temp;
