function [K,err,coef] = modsdh2 ( sys, T, ordK, alpha, beta, opt )
%MODSDH2 H2-optimal reduced-order modal digital controller.
%
%     [K,ERR,PAR] = MODSDH2 ( SYS, T, N, ALPHA, BETA, OPTIONS )
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
%     K   - modal H2-optimal discrete controller of order N
%     ERR - cost function for the optimal system
%     PAR - vector of parameters in [0,1]
%
%   See also MODSDL2, SDH2, SECTOR.

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
        opt.go_func  = 'go_sdh2p';
        opt.f_func   = 'f_sdh2p';
        opt.dcgain   = 0;
        opt.ssmodel  = '[gloptData.dsys,gloptData.gamma] = sdgh2mod ( sys, T );';
        opt.ABE_func = 'sdh2coef';
%------------------------------------------------------
%%       Call function 'modalopt'
%------------------------------------------------------
        [K,err,coef] = modalopt ( sys, T, ordK, alpha, beta, opt );

%------- End of MODSDH2.M --------- KYuP ----------           
