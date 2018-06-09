function disp(X)
%DISP   Displays a sym as text.
%   DISP(S) displays the scalar or array sym,
%   without printing the the sym name.

%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.12.4.2 $  $Date: 2004/04/16 22:22:16 $

loose = isequal(get(0,'FormatSpacing'),'loose');
margin = '     '; % Added by KYuP 03.01.07
if isempty(X) 
   disp('[ empty sym ]')
elseif all(size(X) == 1)
   s = X.s;
   s(s=='`') = [];
   %disp(s)% Commented by KYuP 03.01.07
   s = deconj(s);% Added by KYuP 03.01.07
   disp([margin s])% Added by KYuP 03.01.07
else
   % Find maximum string length of the elements of a X
   p = size(X);
   d = ndims(X);
   for k = 1:prod(p)
      lengths(k) = length(deconj(X(k).s));
   end;
   len = max(lengths);

   for k = 1:prod(p(3:end))
      if d > 2
         if loose, disp(' '), end
         disp([inputname(1) '(:,:,' int2strnd(k,p(3:end)) ') = '])
         if loose, disp(' '), end
      end
      % Pad each element with the appropriate number of zeros
      for i = 1:p(1)
         str = '[';
         for j = 1:p(2)
            s = X(i,j,k).s;
            s(s=='`') = [];
            s = deconj(s);% Added by KYuP 03.01.07
            str = [str blanks(len-length(s)+1) s ','];
         end
         str(end) = ']';
         %if p(2) == 1; str = str(2:end-1); end % Commented by KYuP 03.01.07
         %disp(str)% Commented by KYuP 03.01.07
         disp([margin str])% Added by KYuP 03.01.07
      end
   end
end
if loose, disp(' '), end

% ------------------------

function s = int2strnd(k,p)
s = '';
k = k-1;
for j = 1:length(p)
   d = mod(k,p(j));
   s = [s int2str(d+1) ','];
   k = (k - d)/p(j);
end
s(end) = [];

%------- End of old DISP.M --------- KYuP ----------

%###################################################
% Added by KYuP 03.01.07
%###################################################
function s = deconj ( s )
%-------------------------------------------------    
%DECONJ Change conj(x) for x'
%-------------------------------------------------    
    [pos,match] = regexp(s, 'conj\(.*?\)', 'start', 'match');
    for i=length(pos):-1:1
      x = match{i};
      x = transp(x(6:end-1));
      s = [s(1:pos(i)-1) x s(pos(i)+length(match{i}):end)];
    end;
%------- End of DECONJ.M --------- KYuP ----------

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
        case {'*','*'}
           if nest == 0  &&  prty >= 2, 
              pos = i; prty = 2; 
           end; 
      end;
    end;
%-------------------------------------------------    
%   Operation not found
%-------------------------------------------------        
    if pos < 0 % operation not found
       if s(1) == '('  &&  s(end) == ')'
          s = transp(s(2:end-1)); 
       else
          s = [s '''']; 
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
       s = [transp(s(pos+1:end)) '\' transp(s(1:pos-1))]; 
       return; 
    end;
    if s(pos) == '\' 
       s = [transp(s(pos+1:end)) '/' transp(s(1:pos-1))]; 
       return; 
    end;  
    
%------- End of TRANSP.M --------- KYuP ----------
