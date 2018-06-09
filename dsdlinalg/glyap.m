function [X,E] = glyap ( A1, B1, A2, B2, C )
%GLYAP Solve generalized Lyapunov equation A1*X*B1 + A2*X*B2 = C.
%
%     [X,E] = GLYAP ( A1, B1, A2, B2, C )
%
%   Inputs:
%       A1, B1, A2, B2 - matrices of compatible dimension
%
%   Outputs:
%       X - solution matrix
%       E - error matrix
%
%   See also LYAP.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        [ra1,ca1] = size(A1);
        [ra2,ca2] = size(A2);
        [rb1,cb1] = size(B1);
        [rb2,cb2] = size(B2);
        [rc,cc]   = size(C);       
%        if ra1 ~= ca1  ||  ra2 ~= ca2  ||  rb1 ~= cb1  ||  rb1 ~= cb2  || rc ~= cc
%           error('All matrices should be square');
%        end;       
        if ra1 ~= ra2  ||  ca1 ~= ca2  ||  rb1 ~= rb2  ||  cb1 ~= cb2  ||  ra1 ~= rc  ||  cb1 ~= cc
           error('Matrices should have compatible dimensions');
        end;
%------------------------------------------------------
%       Use equality
%           vec(A*X*B) = kron(B',A) * vec(X),
%       where vec(Y) = [y1' y2' ... yn']', yi is the i-th column of Y.
%       Then, 
%           vec(A1*X*B1 + A2*X*B2) 
%                = [kron(B1',A1) + kron(B2',A2)] * vec(X) = vec(C)
%
%   [1] K. Zhou, J.C. Doyle, K. Glover, Robust and Optimal Control,
%       Prentice Hall, Upper Saddle River, New Jersey, 1995.
%------------------------------------------------------
        A = kron(B1',A1) + kron(B2',A2);
        vecC = C(:);
        [vecX,xx,flag] = caxby ( A, 1, vecC );       
        if flag == 0
           warning('Generalized Lyapunov equation has no exact solution');
        end;
        X = reshape(vecX, ca1, rb1);
        E = A1*X*B1 + A2*X*B2 - C;

%------- End of GLYAP.M --------- KYuP ----------           
