function sys = horzcat(varargin)
%HORZCAT Horizontal concatenation for STPBC: C = [A B] 
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
        sys = varargin{1};
        [rows,cols] = size(sys);
        if isnumeric(sys), sys = stpbc(sys); end;       
%------------------------------------------------------
%	Horizontal concatenation
%------------------------------------------------------
        for i=2:nargin
           X = varargin{i};
           [rowsX,colsX] = size(X);
           if min(rowsX,colsX) > 0
              if rows == 0,
                 sys = stpbc(zeros(rowsX,0));
                 rows = rowsX;
              end;              
              if isnumeric(X), X = stpbc(X); end;             
              if rowsX ~= rows, error('Incompatible dimensions'); end;
              sys.a = blckdiag(sys.a, X.a);
              sys.b = blckdiag(sys.b, X.b);
              sys.c = [sys.c X.c];
              sys.d = [sys.d X.d];
              sys.om = blckdiag(sys.om, X.om);
              sys.ups = blckdiag(sys.ups, X.ups);
           end           
        end;
%------------------------------------------------------
%	Minimal realization
%------------------------------------------------------
    sys = minreal(sys);
    
%------- End of HORZCAT.M --------- KYuP ----------

