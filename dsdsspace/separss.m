function [syss,sysu,ms,mu] = separss ( sys0, type, stype, alpha )
%SEPARSS Stable/antistable separation of an LTI model.
%        State-space technique.
%             SYS = SYSS + SYSU
%
%     [SYSS,SYSU,MS,MU] = SEPARSS ( SYS, TYPE, STYPE, ALPHA )
%     [SYSS,SYSU,MS,MU] = SEPARSS ( SYS, E )
%
%   Inputs:
%     SYS - an LTI system
%     TYPE - type of separation
%          's'  - stable poles are in the left half-plane  
%                 (default for continuous-time systems) 
%          'z'  - stable poles are inside the unit disk
%                 (default for discrete-time systems) 
%          'd'  - stable poles are outside the unit disk
%          'T'  - stable poles are inside the delta-disk 
%                 with center at (-1/T,0) with radius 1/T
%          's0', 'z0', 'd0' - include neutral poles in the
%                 stable part
%      STYPE - type of stable part
%          'infs' - stable part include DC gain (default)
%          'infu' - unstable part include DC gain
%          'symm' - symmetric separation of DC gain
%	   ALPHA - desired degree of stability of "stable" part (in 's')	
%      E     - eigenvalues that must be in the first diagonal block
%              (special syntax)  
%
%   Outputs:
%      SYSS - stable part of the LTI system 
%      SYSU - antistable part of the LTI system (including direct gain)
%      MS   - number of states of the stable part
%      MU   - number of states of the antistable part
%
%   See also SEPARTF.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-May-2006 $
%------------------------------------------------------
%	    Call separation routine without taking minimal realization
%------------------------------------------------------   
        if ~exist('stype', 'var'), stype = 'infs'; end;
        if ~strcmp(stype,'infs')  &&  ~strcmp(stype,'infu')  &&  ~strcmp(stype,'symm')
           error('Incorrect properness type ''%s''',stype);
        end;
        cls = class(sys0);
        if isequal(cls,'ss')
          sys = sys0;  
          imp = num2cell ( zeros(size(sys0)) );            
        else
          [imp,sys] = improper(sys0);
        end;
        sys = ss(sys);        
        switch nargin
          case 1, [syss,sysu,ms,mu] = separssnc ( sys );
          case 2, [syss,sysu,ms,mu] = separssnc ( sys, type );
          case 3, [syss,sysu,ms,mu] = separssnc ( sys, type, stype );
          case 4, [syss,sysu,ms,mu] = separssnc ( sys, type, stype, alpha );
        end;        
%------------------------------------------------------
%	    Transform to required type
%------------------------------------------------------
	    syss = minreal ( syss, sqrt(eps), 0 ); 
        ms = size(syss.a,1);
	    sysu = minreal ( sysu, sqrt(eps), 0 ); 
        mu = size(sysu.a,1);
	    if isequal(cls,'tf')  
           syss = tf(syss);  
           sysu = tf(sysu); 
        end;
    	if isequal(cls,'zpk') 
           syss = zpk(syss); 
           sysu = zpk(sysu); 
        end;
%------------------------------------------------------
%	    Take improper part
%------------------------------------------------------
        [rows,cols] = size(sys);
        if ~iscell(imp), imp = {imp}; end;
        for r=1:rows
          for c=1:cols
            if norm(imp{r,c}) > eps
              add = zpk(imp{r,c}, 1, sys.Ts);  
              if isa(sys0,'tf'), add = tf(add); end; 
              switch stype
                case 'infs',
                  if isa(syss,'zpk'), 
                       syss = sumzpk(syss, add);
                  else syss = minreal(syss + add); 
                  end;
                case 'infu',
                  if isa(sysu,'zpk'), 
                       sysu = sumzpk(sysu, add);
                  else sysu = minreal(sysu + add); 
                  end;
                case 'symm',
                  if isa(syss,'zpk'), 
                    syss = sumzpk(syss, add/2);
                    sysu = sumzpk(sysu, add/2);
                  else
                    syss = minreal(syss + add/2);
                    sysu = minreal(sysu + add/2);
                  end;
              end;                           
            end;
          end;
        end;

%------- End of SEPARSS.M --------- KYuP ----------


