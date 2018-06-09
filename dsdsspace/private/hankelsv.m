function [h,ps,qs,pu,qu] = hankelsv ( a, b, c, type )
%HANKELSV Hankel singular values for an LTI system.
%         For unstable part of the system, Hankel singular 
%         values of the conjugate system are used.
%
% 	  [H,PS,QS,PU,QU] = HANKELSV ( SYS ) 
%     [H,PS,QS,PU,QU] = HANKELSV ( A, B, C ) 
%     [H,PS,QS,PU,QU] = HANKELSV ( A, B, C, 'z' ) 
%
%   Inputs:
%	  SYS - an LTI system model
%	  A, B, C - system state space matrices
%
%   Outputs:
%	  H   - Hankel singular values	
%	  PS, PU - reachability grammians of stable and antistable parts
%	  QS, QU - observability grammians of stable and antistable parts
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check parameters
%------------------------------------------------------
        if nargin == 1
           sys = ss(a);  	
	    else
	       sys = ss(a, b, c, zeros( size(c,1), size(b,2) ));
           if exist('type','var')
              if isequal(type,'z'), sys.Ts = 1; 
              elseif ~isequal(type,'s')
                 error('Unknown type of system: ''%s''',type); 
              end;
           end;
        end;
        if isdt(sys), type = 'z'; else type = 's'; end;	
%------------------------------------------------------
%       Separate stable and antistable parts
%------------------------------------------------------
   	    [syss,sysu,ms,mu] = separss ( sys, type );
        n = size(sys.a);
   	  %-------------------------------------------
   	  %  Stable part
      %-------------------------------------------
        if ms > 0           
   	       ps = gram ( syss, 'c' );
   	       qs = gram ( syss, 'o' );
   	       hs = sqrt ( eig(ps*qs));
   	       hs = sort ( real(hs) );
   	       hs = hs ( ms:-1:1, : );
        else
           hs = []; ps = []; qs = []; 
        end;
   	  %-------------------------------------------
   	  %  Anti-stable part
      %-------------------------------------------
        if mu > 0
	       sysu = sysu';
           pu = gram ( sysu, 'c' );
	       qu = gram ( sysu, 'o' );
	       hu = sqrt ( eig(pu*qu) );
	       hu = sort ( real(hu) );
	       hu = hu ( n-ms:-1:1, : );
        else
           hu = []; pu = []; qu = [];
        end;
	
	    h = [hs; hu];

%------- End of HANKELSV.M --------- KYuP ----------
