function Mi = invmat ( M, recFlag )    
%INVMAT Inversion of symbolic block matrix.
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 03-Jan-2007 $
%------------------------------------------------------
    if ~exist('recFlag','var'), recFlag = 0; end;
    try
      M = double(M);
      if isnumeric(M)
        Mi = sym(inv(M));
        return;
      end;
    catch
    end;
%----------------------------------------------
% 1x1 block case
%----------------------------------------------
    n = size(M,1);
    if n == 1, 
       M = char(M);
       flag = (M(1) == '-');
       if flag, M = M(2:end); end; 
       Mi = sym(sprintf('inv(%s)',M));
       Mi = simprod(Mi);
       if flag, Mi = - Mi; end;           
       return;
    end;
%----------------------------------------------
%   Find nonsingular diagonal A-block
%----------------------------------------------
    for i=1:n-1
      A = M(1:i,1:i);
      if det(A) ~= 0, 
        B = M(1:i,i+1:end);
        C = M(i+1:end,1:i);
        D = M(i+1:end,i+1:end);
        Ai = invmat(A);
        W = invmat(D - C*Ai*B);
        X = Ai + Ai*B*W*C*Ai;
        Mi = [X -Ai*B*W
              -W*C*Ai W];
        Mi = simpinv(Mi);
        return;  
      end;
    end;
%----------------------------------------------
% Find nonsingular diagonal D-block
%----------------------------------------------
    for i=n:-1:2
      D = M(i:n,i:n);
      if det(D) ~= 0, 
        A = M(1:i-1,1:i-1);
        B = M(1:i-1,i:n);
        C = M(i:n,1:i-1);
        Di = invmat(D);
        X = invmat(A - B*Di*C);
        W = Di + Di*C*X*B*Di;
        Mi = [X -X*B*Di
              -Di*C*X W];
        Mi = simpinv(Mi);
        return;  
      end;
    end;   
%----------------------------------------------
% Try block-wise inversion for unitary transform
%----------------------------------------------
    if ~recFlag
      Mi = invmat(flipud(M), 1);
      Mi = flipud(Mi);
      return;
    end;
%----------------------------------------------
% No simplifications...
%----------------------------------------------
     Mi = inv(M);
    
%------- End of INVMAT.M --------- KYuP ----------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = simpinv ( M )   
%SIMPINV Simplify product after inversion
%
%  inv(A)*A = A*inv(A) = I
%  inv(A)^2*A = A*inv(A)^2 = inv(A)
%
%----------------------------------------------
%   Matrix case
%----------------------------------------------   
    [rr,cc] = size(M);
    if any([rr cc] > 1)
      for r=1:rr
        for c=1:cc
          M(r,c) = simpinv(M(r,c));  
        end;
      end;
      return; 
    end;
%----------------------------------------------
%   Scalar case
%----------------------------------------------   
    X = char(M);
    [pos,Xinv] = regexp(X,'inv\(.*?\)','start','match');
    if isempty(pos), return; end;
    for i=1:length(pos)
      Q = Xinv{i}(5:end-1);        
      for k=1:4
        switch k
          case 1, rexp = [Q '\*inv\(' Q '\)']; what = '1';      
          case 2, rexp = ['inv\(' Q '\)\*' Q]; what = '1';  
          case 3, rexp = [Q '\*inv\(' Q '\)\^2']; what = ['inv(' Q ')'];  
          case 4, rexp = ['inv\(' Q '\)\^2\*' Q]; what = ['inv(' Q ')'];  
        end;
        [p,Xi] = regexp(X, rexp, 'start', 'match');      
        if ~isempty(p)
          p = p(1); 
          X = [X(1:p-1) what X(p+length(Xi{1}):end)];
          M = simpinv(simplify(sym(X)));
          return;
        end;
      end;
    end;

%------- End of SIMPINV.M --------- KYuP ----------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = simprod ( M, strip )   
%SIMPROD Simplify inversion of product
%
%  inv(A*(inv(B)*C) = inv(C)*B*inv(A)
%
%----------------------------------------------
%   Default values
%----------------------------------------------
    if ~exist('strip','var'), strip = 1; end;
%----------------------------------------------
%   Default values
%----------------------------------------------
    X = char(M);
    if strip, X = X(5:end-1); end;
    [pos,Xinv] = regexp(X,'inv\(.*?\)','start','match');
    if ~isempty(pos)
      if isequal(pos(1),1)  &&  isequal(Xinv{1},X)
        M = sym(X(5:end-1));  
        return;  
      end;
    end;
%----------------------------------------------
%   Is it a product?
%----------------------------------------------
    nest = X*0;
    n = length(X);
    nest = 1; pos = -1;
    for i=1:n
      switch X(i)  
        case '(', nest = nest+1;    
        case ')', nest = nest-1;    
        case {'*'}, 
          if nest == 0, pos = i; end;  % a product 
        case {'+','-'}, 
          if nest == 0, return; end;  % not a product 
      end;
    end;
%----------------------------------------------
%   Handle a product recursively
%----------------------------------------------
    if pos > 0
      M = ['inv(' X(pos+1:end) ')*'];
      M = [M char(simprod(X(1:pos-1),0))];
      try
        M = sym(M);
      catch
        keyboard
      end;
      return;
    end;
%------- End of SIMPROD.M --------- KYuP ----------
    
    

