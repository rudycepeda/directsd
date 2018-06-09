function [syss,sysu,ms,mu] = separssnc ( sys, type, stype, alpha )
%SEPARSSNC Stable/antistable separation of an LTI model: SYS = SYSS + SYSU.
%          Without minimal realization.
%
%     [SYSS,SYSU,MS,MU] = SEPARSSNC ( SYS, TYPE, STYPE, ALPHA )
%     [SYSS,SYSU,MS,MU] = SEPARSSNC ( SYS, E )
%
%   Inputs:
%     SYS - an LTI system
%     TYPE - type of separation
%         's'  - stable poles are in the left half-plane  
%                (default for continuous-time systems) 
%         'z'  - stable poles are inside the unit disk
%                (default for discrete-time systems) 
%         'd'  - stable poles are outside the unit disk
%         'T'  - stable poles are inside the delta-disk 
%                with center at (-1/T,0) with radius 1/T
%         's0', 'z0', 'd0' - include neutral poles in the
%                stable part
%     STYPE - type of stable part
%         'infs' - stable part include DC gain (default)
%         'infu' - unstable part include DC gain
%         'symm' - symmetric separation of DC gain
%	  ALPHA - desired degree of stability of "stable" part (in 's')	
%     E     - eigenvalues that must be in the first diagonal block
%             (special syntax)  
%
%   Outputs:
%     SYSS - stable part of the LTI system 
%     SYSU - antistable part of the LTI system (including direct gain)
%     MS   - number of states of the stable part
%     MU   - number of states of the antistable part
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%	Determine type of separation
%------------------------------------------------------   
        if ~isobject(sys)
           error('First argument must be an LTI system');
        end;
        if ~exist('type', 'var'),
           if isdt(sys),
                type = 'z'; 
           else type = 's'; 
           end;
        else 
           if isnumeric(type), 
              eig1 = type; type = 'e'; 
              if nargin > 3, error('Too many arguments for this syntax'); end;
           end; 
        end;        
        if ~strcmp(type,'s')  &&  ~strcmp(type,'z')  &&  ~strcmp(type,'d')  &&  ~strcmp(type,'e') && ...
           ~strcmp(type,'s0')  &&  ~strcmp(type,'z0')  &&  ~strcmp(type,'d0')  &&  ~strcmp(type,'T')
           error('Incorrect separation type ''%s''', type);
        end;        
        if ~exist('stype', 'var'), stype = 'infs'; end;
        if ~strcmp(stype,'infs')  &&  ~strcmp(stype,'infu')  &&  ~strcmp(stype,'symm')
           error('Incorrect properness type');
        end;        
        if ~exist('alpha', 'var'), alpha = 0; end;
        dalpha = exp(-alpha*sys.Ts);        
%------------------------------------------------------
%	Extract system matrices
%------------------------------------------------------
        [a,b,c,d] = ssdata ( ss(sys) );
        n = size(a,1);
        nin  = size(b,2);
        nout = size(c,1);
%------------------------------------------------------
%	Complex Schur decompisition
%	    A = U * T * U'
%------------------------------------------------------
        [u,t] = schur ( a );       
        [u,t] = rsf2csf ( u, t );
        eigA = diag(t);
%------------------------------------------------------
%	Handle poles at the stability region
%------------------------------------------------------
        if isct(sys)
           ind = find( abs(eigA) < sqrt(eps) );
           eigA(ind) = 0;
        else
           ind = find( abs(eigA-1) < sqrt(eps) );
           eigA(ind) = 1;
        end;
%------------------------------------------------------
%	Determine order of eigenvalues in Schur decompisition
%------------------------------------------------------
        switch type
          case 's',
             ms = length ( find(real(eigA) < - alpha) );   
             ord = sign ( real(eigA) + alpha );
          case 'z',
             ms = length ( find(abs(eigA) < dalpha) );   
             ord = sign ( abs(eigA) - dalpha );
          case 'd',
             ms = length ( find(abs(eigA) > dalpha) );   
             ord = - sign ( abs(eigA) - dalpha );
          case 'T',
             T = sys.Ts;
             if (sys.var ~= 'q')  ||  (T <= 0), 
                error('Transfer function must be in delta-domain with Ts > 0'); 
             end;
             ms = length ( find(abs(eigA+1/T) < 1/T) );   
             ord = sign ( abs(eigA+1/T) - 1/T );
          case 's0',
             ms = length ( find(real(eigA) <= - alpha) );   
             ord = sign ( real(eigA) + alpha - eps );
          case 'z0',
             ms = length ( find(abs(eigA) <= dalpha) );   
             ord = sign ( abs(eigA) - dalpha - eps );
          case 'd0',
             ms = length ( find(abs(eigA) >= dalpha) );   
             ord = - sign ( abs(eigA) - dalpha + eps );
          case 'e', % if isequal(type,'e')
             remEigs = eigA;
             ms = length ( eig1 );   
             ord = ones(1,n); 
             ind = 1:n;
             while length(eig1) > 0
               if isempty(remEigs)
                  warning('Too many eigs specified for the first block'); 
                  break; 
               end;
               diff = abs(remEigs - eig1(1));
               [err,i] = min(diff);
               if err > 1e-6, warning('Too large eig error %g',err); end;    
               ord(ind(i)) = 0;
               remEigs(i) = []; ind(i) = []; eig1(1) = [];
             end;
        end;

        mu = n - ms;

%------------------------------------------------------
%	Order eigenvalues in Schur decompisition
%   using complex Givens rotation
%------------------------------------------------------
        [u,t] = cschurord ( u, t, ord );
%------------------------------------------------------
%	Transformation from complex to real ordered 
%   block Schur form such that 
%	    A = U * ABLOCK * U'
%------------------------------------------------------
        [u,ablock] = csf2rsf ( u, t, 1e-5 );       
%------------------------------------------------------
%	State-space separation
%------------------------------------------------------     	
    %--------- No stable part ---------------------  
	if ms == 0   
	  as = zeros(n); bs = zeros(n,nin); cs = zeros(nout,n);
  	  au = a; bu = b; cu = c; 
      switch stype
         case 'symm',
           error('Incorrect symmetric separation');
         case 'infs',
           ds = d; 
           du = zeros(nout,nin);
         case 'infu',
           du = d; 
           ds = zeros(nout,nin);
      end;
    %--------- No anti-stable part ---------------------  
	elseif mu == 0 
  	  as = a; bs = b; cs = c;
	  au = zeros(n); bu = zeros(n,nin); cu = zeros(nout,n);
      switch stype
         case 'symm',
           error('Incorrect symmetric separation');
         case 'infs',
           ds = d; 
           du = zeros(nout,nin);
         case 'infu',
           du = d; 
           ds = zeros(nout,nin);
       end;
    %--------- General case ---------------------  
    else   
  	   as = ablock (1:ms, 1:ms);
   	   au = ablock (ms+1:n, ms+1:n);
   	   a12 = ablock (1:ms, ms+1:n);

   	   v1 = u (:, 1:ms); 
   	   v2 = u (:, ms+1:n);

   	   x = lyap ( as, -au, a12 );
    
   	   t1 = v1;
   	   t2 = v1 * x + v2;
   	   s1 = v1' - x * v2';
   	   s2 = v2';
   	   bs = s1 * b;    bu = s2 * b;
   	   cs = c * t1;    cu = c * t2;

       switch stype
         case 'symm',
           ds = d;
           for i=1:nin, ds(i,i) = ds(i,i) / 2; end;
           for i=2:nin, for j=1:i-1, ds(i,j) = 0; end; end;
           du = ds';
         case 'infs',
           ds = d; 
           du = zeros(nout,nin);
         case 'infu',
           du = d; 
           ds = zeros(nout,nin);
       end;
	end
%------------------------------------------------------
%	Final result
%------------------------------------------------------
	syss = ss ( as, bs, cs, ds, sys.Ts );
	sysu = ss ( au, bu, cu, du, sys.Ts );
%------------------------------------------------------
%	Restore initial form
%------------------------------------------------------
	if isa(sys,'tf'),  syss = tf(syss); sysu = tf(sysu);	end;
	if isa(sys,'zpk'), syss = zpk(syss); sysu = zpk(sysu);	end;

%------- End of SEPARSSNC.M --------- KYuP ----------


