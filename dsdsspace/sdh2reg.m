function [Kopt,gamma,Poles] = sdh2reg ( csys, T, o2, i2, options )
%SDH2REG H2-optimal controller for sampled-data system.
%	     Equivalent discrete model approach.	
%
%     [K,GAMMA,POLES] = SDH2REG ( SYS, T, NMEAS, NCON, OPTIONS )
%
%   Inputs:
%     SYS - an LTI system
%     T   - sampling period
%     NMEAS - dimension of the vector 'y' (default 1)
%     NCON  - dimension of the vector 'u' (default 1)
%     OPTIONS - structure with possible fields
%         .method  - solution method of
%           'gh2'  - generalized H2-measure (default)
%           'sh2'  - simple H2-problem
%           'fast' - using fast discretization
%         .Nfast   - number of intervals for fast discretization
%                    (default value is 5) 
%   Outputs:
%     K     - H2-optimal controller as an LTI model
%	  GAMMA - optimal value of gamma = ||Tzw||_2
%	  POLES - poles of the closed-loop system in 'z'-plane
%
%   See also SDHINFREG, SDH2, SDH2NORM, SDL2.

% References:
%   [1] T. Chen, B. Francis, Optimal Sampled-Data Control Systems, 
%       Springer-Verlag, Berlin etc., 1995.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%	Check syntax
%------------------------------------------------------
    if exist('options','var'),
       if isfield(options,'Nfast'),  Nfast = 5; end;
       if isfield(options,'method'), method = options.method; end;
    end;
    if ~exist('Nfast', 'var'),  Nfast = 5; end;
    if ~exist('method', 'var'), method = 'gh2'; end;
    if ~isequal(method,'gh2')  &&  ~isequal(method,'sh2')  &&  ~isequal(method,'fast')
       error('Unknown solution method ''%s''', method); 
    end;    
%------------------------------------------------------
%	Continuous-time state space
%------------------------------------------------------
    csys = minreal(ss(csys), sqrt(eps), 0);
  	[nout,nin] = size(csys);    
%------------------------------------------------------
%	Check inputs and outputs
%------------------------------------------------------
    if ~exist('i2','var'), i2 = 1; end;
    if ~exist('o2','var'), o2 = 1; end;
	i1 = nin - i2;
	o1 = nout - o2;
	if i1 < 1, error('Negative or zero number of reference inputs'); end;
	if o1 < 1, error('Negative or zero number of working outputs');	end;
	if i1 >= nin, error('No control signals specified'); end;
	if o1 >= nout, error('No measured signals specified'); end;      
%------------------------------------------------------
%   Construct equivalent discrete-time model
%------------------------------------------------------
    gamma = 0;
    switch method
      case 'gh2',  [dsys,gamma] = sdgh2mod (csys, T, o2, i2);
      case 'sh2',  dsys = sdh2simple (csys, T, o2, i2);
      case 'fast', dsys = sdfast (csys, T, Nfast, o2, i2, 'l2');
    end;    
%--------------------------------------------------------------
%   Using discrete-time H2-optimal design
%--------------------------------------------------------------
    Kopt = h2reg ( dsys, o2, i2 );
%--------------------------------------------------------------
%   Constructing closed-loop system model
%--------------------------------------------------------------
    if nargout > 1
       ClSys = minreal ( lft ( dsys, Kopt ), sqrt(eps), 0 );
       gamma = sqrt(norm(ClSys)^2 + gamma);
       if nargout > 2
          Poles = eig ( ClSys.a );
       end;
    end;

%------- End of SDH2REG.M --------- KYuP ----------
