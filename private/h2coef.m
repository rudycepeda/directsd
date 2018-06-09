function [A,B,E] = h2coef ( sys )
%H2COEF Coefficients for H2-optimal LTI system design.
%
%     [A,B,E] = H2COEF ( SYS )
%
%   Inputs:
%     SYS - an LTI system
%
%   Outputs:
%     A, B, E - function coefficients in quadratic functional
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
[nout,nin] = size(sys);
o2 = 1; i2 = 1;
i1 = nin - i2;
o1 = nout - o2;
if i1 < 1, error('Incorrect number of inputs'); end;
if o1 < 1, error('Incorrect number of outputs'); end;
%------------------------------------------------------
%       Extract transfer matrices
%------------------------------------------------------
P11 = sys(1:o1,1:i1);
P12 = sys(1:o1,i1+1:nin);
P21 = sys(o1+1:nout,1:i1);
%P22 = sys(o1+1:nout,i1+1:nin);
%------------------------------------------------------
%       Construct A
%------------------------------------------------------
A0 = minreal(P12'*P12);
A1 = minreal(P21*P21');
A = minreal(A0*A1);
%------------------------------------------------------
%       Construct B
%------------------------------------------------------
%B = minreal(P21*P11H*P12);
P11H = P11';
B = zpk(0); B.Ts = sys.Ts;
for i=1:i1
    Bi = zpk(0); Bi.Ts = sys.Ts;
    for j=1:o1
        Bi = sumzpk(Bi,minreal(P11H(i,j)*P12(j)));
    end;
    B = sumzpk(B, minreal(P21(i)*Bi));
end;
%------------------------------------------------------
%       Construct E
%------------------------------------------------------
if nargout > 2
    PT = P11H*P11;
    n = size(PT,1);
    E = PT(1,1);
    for jj = 2:n
        E = E+PT(jj);
    end    
end;

%------- End of SDH2COEF.M --------- KYuP ----------
