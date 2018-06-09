function sys = vertcat(varargin)
%VERTCAT Vertical concatenation for STPBC: C = [A; B] 
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
              if cols == 0,
                 sys = stpbc(zeros(0,colsX));
                 cols = colsX;
              end;              
              if isnumeric(X), X = stpbc(X); end;             
              if colsX ~= cols, error('Incompatible dimensions'); end;
              sys.a = blckdiag(sys.a, X.a);
              sys.b = [sys.b; X.b];
              sys.c = blckdiag(sys.c, X.c);
              sys.d = [sys.d; X.d];
              sys.om = blckdiag(sys.om, X.om);
              sys.ups = blckdiag(sys.ups, X.ups);
           end           
        end;
%------------------------------------------------------
%	Minimal realization
%------------------------------------------------------
    try
      sys = minreal(sys);
    catch end;
    
%------- End of VERCAT.M --------- KYuP ----------
