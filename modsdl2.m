function [K,err,coef] = modsdl2 ( sys, T, ordK, alpha, beta, opt )
%MODSDL2 L2-optimal reduced-order modal digital controller.
%
%     [K,ERR,PAR] = MODSDL2 ( SYS, T, N, ALPHA, BETA, OPTIONS )
%
%   Inputs:
%     SYS - an LTI system
%     T   - sampling period
%     N   - controller order (default = plant order - 1)
%     ALPHA - degree of stability (default 0)
%     BETA  - oscillation degree (default Inf)
%     OPTIONS - options, structure with the following possible fields
%            .method - optimization method: 'global' (default),
%                      'randsearch', 'simanneal' or 'local'
%            for other fields, see help on 'infglob', 'randsearch',
%            'simanneal' or 'fmincon', respectively  
%   Outputs:
%     K   - modal L2-optimal discrete controller of order N
%     ERR - cost function for the optimal system
%     PAR - vector of parameters in [0,1]
%
%   See also MODSDH2, SDL2, SECTOR.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%%       Check data
%------------------------------------------------------
        if ~exist('ordK','var'),  ordK = []; end; 
        if ~exist('alpha','var'), alpha = 0; end; 
        if ~exist('beta','var'),  beta = Inf; end; 
        if ~exist('opt','var'),   opt = struct; end; 
%------------------------------------------------------
%%       Specific function names
%------------------------------------------------------
        opt.go_func  = 'go_sdl2p';
        opt.f_func   = 'f_sdl2p';
        if ~isfield(opt,'dcgain'), opt.dcgain   = 1; end;
        opt.ssmodel  = 'gloptData.dsys = sdh2simple ( sys, T );';
        opt.ABE_func = 'sdl2coef';
%------------------------------------------------------
%%       Call function 'modalopt'
%------------------------------------------------------
        [K,err,coef] = modalopt ( sys, T, ordK, alpha, beta, opt );

%------- End of MODSDL2.M --------- KYuP ----------           

