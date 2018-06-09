function [U,S,invU,n1,n2,E] = schurdiag ( S, type, alpha )
%SCHURDIAG Block diagonalization of the Schur form.
%          Apply ordering of eigenvalues.
%              F = U * S * inv(U)
%
%     [U,S,IU,N1,N2,E] = SCHURDIAG(F, TYPE, ALPHA)
%
%   Input:
%     F    - a constant matrix
%     TYPE - type of sorting
%         's' - ascending real parts (default)
%         'i' - descending real parts
%         'z' - ascending absolute value
%         'd' - descending absolute value
%     ALPHA - parameter for eigenvalue separation
%             (default is 0 for types 's' and 'i', and 
%                         1 for types 'z' and 'd')
%   Outputs:
%     S - block diagonal matrix with ordered complex eigenvalues
%                | T1  0 |
%            S = |  0 T2 | , where T1 and T2 are square blocks such that
%         spec(T1) in Re x < alpha (type = 's')   
%                  in Re x > alpha (type = 'i')   
%                  in  |x| < alpha (type = 'z')   
%                  in  |x} > alpha (type = 'd')   
%     U - nonsingular matrix
%     IU = inv(U)
%     N1, N2 - size of diaginal blocks
%     E - eigenvalues of S
%
%   See also SCHUR, CSCHURORD, CSF2RSF.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check parameters
%------------------------------------------------------
        if exist('type','var')
           if ~isequal(type,'s') && ~isequal(type,'i') && ...
              ~isequal(type,'z') && ~isequal(type,'d')
              error('Incorrect sorting parameter');
           end;
        else type = 's'; 
        end;
        if ~exist('alpha','var')
           if isequal(type,'s') || isequal(type,'i'), 
                alpha = 0; 
           else alpha = 1; 
           end;
        end;
%------------------------------------------------------
%       Construct ordered real Schur form
%       and find block sizes
%------------------------------------------------------
        [Tbal,S] = balance(S);
        [U,S,E] = schurord ( S, type );
        switch type 
           case 's', bl1 = find(real(E) < alpha);  
           case 'i', bl1 = find(real(E) > alpha);  
           case 'z', bl1 = find( abs(E) < alpha);  
           case 'd', bl1 = find( abs(E) > alpha);  
        end;
        n1 = length(bl1);
%------------------------------------------------------
%       Eliminate (1,2) block using Silvester equation
%         T11*X + X*T22 + T12 = 0  
%------------------------------------------------------
        n2 = size(S,1) - n1;
        if n1 == 0  ||  n2 == 0, 
           U = Tbal * U;
           invU = U' / Tbal;
           return; 
        end; 
        ind1 = 1:n1; ind2 = n1+1:n1+n2;
        T11 = S(ind1,ind1);
        T12 = S(ind1,ind2);
        T22 = S(ind2,ind2);
        X = lyap ( T11, -T22, T12 ); 
        Y = [eye(n1) X;zeros(n2,n1) eye(n2)];
        invY = [eye(n1) -X;zeros(n2,n1) eye(n2)];
        invU = invY * U' / Tbal;
        U = Tbal * U * Y;
        S(ind1,ind2) = 0;

%------- End of SCHURDIAG.M --------- KYuP ----------           
