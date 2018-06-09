function [U,S,E] = schurord ( S, type, rctype )
%SCHURORD Complex Schur transformation with ordering eigenvalues.
%              F = U * S * U'
%
%     [U,S,E] = SCHURORD(F, TYPE, RCTYPE)
%
%   Input:
%     F    - a symmetric matrix
%     TYPE - can be 's' - ascending real parts (default)
%                   'i' - descending real parts
%                   'z' - ascending absolute value
%                   'd' - descending absolute value
%     RCTYPE - can be 'real' (default) of 'complex'
%
%   Outputs:
%     S - Schur matrix (see SCHUR) with ordered complex eigenvalues
%     U - Hermitian matrix
%     E - eigenvalues of S
%
%   See also SCHUR, CSCHURORD, CSF2RSF.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check parameters
%------------------------------------------------------
        if ~exist('type','var'), type = 's'; end;
        if ~isequal(type,'s') && ~isequal(type,'i') && ...
           ~isequal(type,'z') && ~isequal(type,'d')
           error('Incorrect sorting parameter ''%s''', type);
        end;
        if ~exist('rctype','var'), rctype = 'real'; end;
        if ~isequal(rctype,'real') && ~isequal(rctype,'complex') && ...
           error('Incorrect result type ''%s''', rctype);
        end;
%------------------------------------------------------
%       Construct complex Schur form
%------------------------------------------------------
        [U,S] = schur ( S );
        [U,S] = rsf2csf ( U, S );
        n = size(S, 1);
%------------------------------------------------------
%       Reorder eigenvalues using Givens rotation
%------------------------------------------------------
        for i = 1:n-1,
           for j = n-1:-1:i,
              Sj = S(j,j); Sj1 = S(j+1,j+1);
              if needswap(Sj, Sj1, type)
                 b = Sj - Sj1; 
                 a = S(j,j+1);
                 aa = abs(a);
                 if aa == 0, 
                    c = 0; 
                    s = 1;
                 else 
                    k = norm([a b]); 
                    c = aa / k; 
                    s = sign(a) * (conj(b) / k);
                 end;
                 G = [c s;-conj(s) c];                 
                 S(j:j+1,:) = G'*S(j:j+1,:); 
                 S(:,j:j+1) = S(:,j:j+1)*G;
                 U(:,j:j+1) = U(:,j:j+1)*G;
                 S(j,j) = Sj1;
                 S(j+1,j+1) = Sj;
              end;
           end;
        end;        
%---------------------------------------------------
%       Restore triangular form
%------------------------------------------------------
        S = triang(S);
        E = diag(S);
        if isequal(rctype,'real') 
           [U,S] = csf2rsf(U, S);        
        end;

%------- End of SCHURORD.M --------- KYuP ----------           

%####################################################
function flag = needswap ( d, e, type )
%
% NEEDSWAP  determines whether two blocks need be swapped
%           according to a prescribed eigenvalue sorting
%
%         F = NEEDSWAP ( d, e, TYPE )
%
%   Inputs:
%       d - eigenvalue of the previous blocks
%       e - eigenvalue of the current block
%       TYPE - can be 's' - ascending real parts
%                     'i' - descending real parts
%                     'z' - ascending absolute value
%                     'd' - ascending absolute value
%
%------------------------------------------------------
        flag = 0;
        if abs(e-d) < 1.e-10, return; end;
        if type == 's'  &&  real(e) < real(d)
           flag = 1;
        end;
        if type == 'i'  &&  real(e) > real(d)
           flag = 1;
        end;
        if type == 'z'  &&  abs(e) < abs(d)
           flag = 1;
        end;
        if type == 'd'  &&  abs(e) > abs(d)
           flag = 1;
        end;
        
%------- End of NEEDSWAP.M --------- KYuP ----------           
