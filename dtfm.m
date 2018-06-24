function sysd = dtfm(sys,T,t,H,type)
%DTFM Discrete transfer function
%   SYSD = DTFM(SYSC,T,t) is the Discrete Transfer Function of the
%   continuous system SYSC, obtained with a sample time T at a given instant
%   t. SYSC can be given in any of the standard Matlab representations for
%   LTI systems. If the argument t is not given, it is assumed as t = 0. A
%   zero order hold with m(t) = 1 for 0 <= t <= T is assumed.
%
%   SYSD = DTFM(SYSC,T,t,H) where H is an LTI model, uses a ZOH with
%   transfer function equivalent to H.
%
%   SYSD = DTFM(SYSC,T,t,H,TYPE) allows the selection of the type of LTI
%   model for SYSD . TYPE can be 'tf', 'zpk' (default) or 'ss'.
%
%See also DTFM2 ZTRM

%------------------------------------------------------
% Copyright 1999-2014
% Originally written by K. Polyakov
% Revised by R. Cepeda-Gomez
% $Revision: 4.00$    $Date: 6-Nov-2014$
%------------------------------------------------------

%Check inputs
%System
if isa(sys,'lti')
    if ~isproper(sys)
        error('The system must be proper');
    else
        sys = ss(sys);
        A = sys.a;
        B = sys.b;
        C = sys.c;
        D = sys.D;
    end
else
    error('SYS must be an LTI system');
end
%Time instant
if ~exist('t','var')
    t = 0;
elseif isobject(t)
    H = t;
    t = 0;
elseif ischar(t)
    type = t;
    H = [];
    t = 0;
end
%Extrapolator
if ~exist('H','var')
    H = [];
elseif ischar(H)
    type = H;
    H = [];
end
if isempty(H)
    H = ss(0,1,1,0); % ZOH with m=1 by default
else
    H = ss(H);
end
if isscalar(H),  % scalar extrapolator
    H = H*eye(size(sys,2));
end
%Check dimensions of the hold.
%Compatibility for FOH was removed by RCG on 16.10.2014
%Revision pending
if size(H,3) == 1
    zoh = 1;
else
    zoh = 0;
    %error('Only zero order hold is supported')
end
%Check compatibility of the dimensions
if size(sys,2) ~= size(H,1)
    error('Dimensions of plant and hold are incompatible');
end

%Type
if ~exist('type','var')  ||  isempty(type)
    type = 'zpk';
else
    type = lower(type);
    if ~isequal(type,'zpk') && ~isequal(type,'ss') && ~isequal(type,'tf')
        error(['Incorrect type ' type]);
    end
end

%Find the minimum delay
delay = totaldelay(sys);
comdelay = min(min(delay));
delay = delay - comdelay;

if norm(delay) < eps %%All delays are equal, no need for special cases
    %Determine the location within the period
    t = t - comdelay;   %Compensating time delay
    m = floor(t/T);     %How many periods have passed?
    h = t - m*T;        %Where within the period is the desired time?
    if m > 0,
        zNum = m;
        zDen = 0;
    else % to avoid problems with static gain !
        zNum = 0;
        zDen = -m;
    end
    if zoh
        %Extract information of the hold device
        L = H.a;
        M = H.b;
        K = H.c;
        if ~isstatic(sys)
            %Perform the integration
            [Gammah,Phi1h,Phi2h] = intaba(A,B*K,L,h);
            [Gamma,Phi] = intaba(A,B*K,L,T);
            %Create the SS matrices
            Ad = Phi;
            Cd = C*Phi1h;
            Bd = Gamma*M;
            Dd = (C*Gammah + D*K*Phi2h)*M;
        else %Static gain case
            Ad = [];
            Bd = [];
            Cd = [];
            Dd = D*K*expm(L*t)*M;
        end
    else
        nx = size(A,1);
        nu = size(H,2);
        %ny = size(H,1);
        beta = norm(H.d(:,:,2));
        %Build matrices Ad and Bd of discrete state space
        %model, they are time-independent
        Ah = [  A      B      zeros(nx,nu)
            zeros(nu,nx+nu) -beta/T*eye(nu)
            zeros(nu,nx+nu)   zeros(nu)     ];
        E = expm(Ah*T);
        Phi = E(1:nx,1:nx);
        Gamma10  = E(1:nx,nx+1:nx+nu);
        Gamma1 = E(1:nx,nx+nu+1);
        Gamma0 = Gamma10 - Gamma1;
        Ad = [Phi Gamma1; zeros(1,nx+nu)];
        Bd = [Gamma0; eye(nu)];
        %Build Cd and Dd matrices depending on t
        Et = expm(Ah*h);
        Phit = Et(1:nx,1:nx);
        Gamma10t  = Et(1:nx,nx+1:nx+nu);
        Gamma1t = Et(1:nx,nx+nu+1);
        Gamma0t = Gamma10t - Gamma1t;
        Cd = [C*Phit C*Gamma1t-D*beta*h/T];
        Dd = C*Gamma0t+D*(1 + beta*h/T);
    end
    sysd = ss(Ad,Bd,Cd,Dd,T);
    %Add zeros or poles according to the ammount of periods passed
    if zNum > 0 %Add zeros
        add = zpk(zeros(1,zNum),[],1,T);
        if isequal(type,'tf')
            sysd = tf(sysd);
            add = tf(add);
        end
        if isequal(type,'ss') || isequal(type,'zpk')
            sysd = zpk(sysd);
        end
        sysd = minreal(sysd*add);
        if isequal(type,'ss') && ~isproper(sysd,'zpk')
            m = 'Improper functions cannot be converted to state space, ';
            m = [m 'reporting result in zpk form.'];
            warning(m);
            type = 'zpk';
        else
            sysd=ss(sysd);
        end
    elseif zDen > 0
        sysd.ioDelay = zDen*ones(size(sysd));
        sysd = minreal(absorbDelay(sysd));
        switch (type)
            case 'zpk'
                sysd = zpk(sysd);
            case 'tf'
                sysd = tf(sysd);
        end
    end
else %Special cases, different delays for each input-output pair
    [no,ni] = size(sys); %Number of outputs and inputs of the plant
    tni = size(H,2); %Total number of inputs
    sysd = feval(type,zeros(no,tni));%Creating space for the output
    sysd.Ts = T;
    if no > 1
        for i = 1:no %Recursive call row by row
            sysd(i,:) = dtfm(sys(i,:),T,t,H,type);
        end
    else
        for j = 1:tni %recursive call element by element inside a row
            for k = 1:ni %input to the hold affects all inputs to the plant
                temp = dtfm(sys(1,j),T,t,H(k,j),type);
                sysd(1,j) = minreal(sysd(1,j)+temp);
            end
        end
    end
end

switch type
    case 'zpk'
        sysd = sdzpk(minreal(sysd));
    case 'tf'
        sysd = sdtf(minreal(sysd));
end

