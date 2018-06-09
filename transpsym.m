function y = transpsym ( x )
%TRANSPSYM Transpose for symbolic expression
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
    if ~isa(x,'sym'), y = x'; return; end;
%-------------------------------------------------    
%   MIMO case
%-------------------------------------------------    
    [rr,cc] = size(x);
    y = x.';
    for r=1:rr
      for c=1:cc
        y(c,r) = sym(transp(char(x(r,c))));  
      end;
    end;

%------- End of TRANSPSYM.M --------- KYuP ----------

%###################################################
function s = transp ( s )
%-------------------------------------------------    
%TRANSP Transpose for symbolic expression
%-------------------------------------------------    
    prty = 100; 
    nest = 0;
    pos = -1;
    for i=1:length(s)
      switch s(i)
        case '(', nest = nest + 1;
        case ')', nest = nest - 1;
        case {'+','-'}
           if nest == 0 && prty >= 1, 
              pos = i; prty = 1; 
           end; 
        case {'*','/','\'}
           if nest == 0  &&  prty >= 2, 
              pos = i; prty = 2; 
           end; 
      end;
    end;
%-------------------------------------------------    
%   Operation not found
%-------------------------------------------------        
    if pos < 0 % operation not found
       if isempty(s), return; end; 
       if s(1) == '('  &&  s(end) == ')'
          s = transp(s(2:end-1)); 
          return;
       end;
       pos = regexp(s, 'conj\(.*?\)' );
       if ~isempty(pos) && pos(1) == 1
          s = s(6:end-1); 
       else
          if ~all(isstrprop(s,'digit')) 
             s = ['conj(' s ')']; 
          end;
       end;
       return; 
    end;
%-------------------------------------------------    
%   Operations
%-------------------------------------------------    
    if prty == 1
       s = [transp(s(1:pos-1)) s(pos) transp(s(pos+1:end))]; 
       return; 
    end;
    if s(pos) == '*' 
       s = [transp(s(pos+1:end)) '*' transp(s(1:pos-1))]; 
       return; 
    end;
    if s(pos) == '/' 
       p = regexp(s,'[^01234567890]');
       if ~isequal(p,pos)
         s = [transp(s(pos+1:end)) '\' transp(s(1:pos-1))]; 
       end;
       return; 
    end;
    if s(pos) == '\' 
       s = [transp(s(pos+1:end)) '/' transp(s(1:pos-1))]; 
       return; 
    end;  
    
%------- End of TRANSP.M --------- KYuP ----------
