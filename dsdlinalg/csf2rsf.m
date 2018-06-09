function [Q,Z,T,TB,D] = csf2rsf (Q, Z, T, TB, tol)
%CSF2RSF Conversion from complex Schur form to real (block) Schur form.
%        2 x 2 blocks are associated with complex eigenvalues.
%
%     [QR,TR,D] = CSF2RSF ( Q, T );
%     [QR,TR,D] = CSF2RSF ( Q, T, TOL );
%
%   Inputs:  
%     Q - unitary matrix so that Q*Q'=Q'*Q=I
%     T - matrix in a complex Schur form (upper triangular)
%         such that Q'*A *Q = T for some square matrix A 
%
%   Outputs: 
%     QR - unitary matrix
%     TR - matrix in a real (block) Schur form such that QR'*A*QR=TR
%     D  - eigenvalues 
%
% Generalized Schur form for a pencil (A - lam*B)
% ------------------------------------------------
%
%     [QR,ZR,TR,TBR,D] = CSF2RSF ( Q, Z, T, TB )
%     [QR,ZR,TR,TBR,D] = CSF2RSF ( Q, Z, T, TB, TOL )
%
%   Inputs:  
%      Q, Z  - square unitary matrices
%      T, TB - complex upper triangular matrices such that
%              Q'*A*Z=T and Q'*B*Z=TB
%      TOL   - tolerance for detecting real eigenvalue
%
%   Outputs:  
%      QR, ZR  - square unitary matrices
%      TR, TBR - complex upper triangular matrices such that
%                     QR'*A*ZR=TR and QR'*B*ZR=TBR
%                TR is real Schur (upper quasi(block)triangular)
%                TBR is real upper triangular
%      D - pencil eigenvalues
%
%   See also SCHUR, CSCHURORD, QZ, QR.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check syntax
%------------------------------------------------------
        [n,c] = size(Q);
        if n ~= c, error('Matrices should be square'); end;               
        switch nargin
          case {2,3}
              if nargout > 3, error('Too many output values'); end;
              if nargin == 3, tol = T; else tol = 1e-10; end;
              T = Z; Z = Q; TB = eye(n);
          case {4,5}    
              if nargin == 4, tol = 1e-10; end;
          otherwise
           error('Incorrect number of arguments');
        end
        [rz,cz] = size(Z);
        [rt,ct] = size(T);
        [rb,cb] = size(TB);
        if any(n ~= [rz cz rt ct rb cb])
           error('All matrices should have the same size');
        end;       
%------------------------------------------------------
%       Main loop
%------------------------------------------------------
        j = sqrt(-1);
        i = 1;
        D = [];       
        %save cs T TB Q Z;        
        while i <= n,

          %#####################################################   
          %              REAL eigenvalues
          %#####################################################            
          if abs( imag ( T(i,i) * conj(TB(i,i)) ) ) < tol 

             [U,R,xx] = qr ( [real(Q(:,i)) imag(Q(:,i))]' );

            %----------------------------------------------------
            %              SIMPLE REAL eigenvalues
            %----------------------------------------------------
             if rank(R,1.e-6) == 1 %--- generic case
                Q1 = j*U(1,2) + U(2,2);
                Q(:,i) = Q(:,i) * Q1;

                if nargin > 2
                   [U,xx] = qr ( [real(Z(:,i)) imag(Z(:,i))]' );
                   Z1 = j*U(1,2) + U(2,2);
                   Z(:,i) = Z(:,i) * Z1;
                else
                   Z1 = Q1;
                end;

                T(:,i) = T(:,i) * Z1;
                T(i,:) = Q1' * T(i,:);

                if nargin > 2
                   TB(:,i) = TB(:,i) * Z1;
                   TB(i,:) = Q1' * TB(i,:);
                end;
             
                if abs(TB(i,i)) < 1.e-10    
                     D = [D; inf];
                else D = [D; real(T(i,i)/TB(i,i))]; 
                end;
             
                i = i + 1;

            %----------------------------------------------------
            %              REPEATED REAL eigenvalues
            %----------------------------------------------------
             else                 
                k = 2;                 
                while 1
                   ind = (i:i+k-1);
                   [U,R,xx] = qr ( [real(Q(:,ind)) imag(Q(:,ind))]' );
                   if rank(R,1.e-5) == k, break; end;
                   k = k + 1;
                   if i + k - 1 > n
                     msg = 'Unable to transform non-standard eigenvalue blocks';
                     disp(msg); %keyboard  
                     error(msg);
                   end;   
                end;
               
                Q1 = j*U(1:k,k+1:2*k) + U(k+1:2*k,k+1:2*k);
                Q(:,ind) = Q(:,ind) * Q1;
             
                if nargin > 3
                   [U,xx] = qr ( [real(Z(:,ind)) imag(Z(:,ind))]' );
                   Z1 = j*U(1:k,k+1:2*k) + U(k+1:2*k,k+1:2*k);
                   Z(:,ind) = Z(:,ind) * Z1;
                else
                   Z1 = Q1;
                end;

                T(:,ind) = T(:,ind) * Z1;
                T(ind,:) = Q1' * T(ind,:);

                if nargin > 3 % ???
                   TB(:,ind) = TB(:,ind) * Z1;
                   TB(ind,:) = Q1' * TB(ind,:);
                      %--- lower k lines of TB should be zero ???
%                   if norm(TB(i:n,i:n)) > sqrt(eps)
%                     msg = 'Unable to transform non-standard eigenvalue blocks';
%                     disp(msg); keyboard  
%                     error(msg);
%                   end;  
                end;
  
                %--- make T upper triangular by Housholder transformation 
                
                for p=i:i+k-2      
                  rowind = p:i+k-1;
                  v = house ( T(rowind,p) );
                  T(rowind,:)  = rowhouse(T(rowind,:), v);
                  TB(rowind,:) = rowhouse(TB(rowind,:), v);
                  Q(:,rowind)  = colhouse ( Q(:,rowind), v );
                  if nargin < 4
                     T (:,rowind) = colhouse(T(:,rowind), v);
                     TB(:,rowind) = colhouse(TB(:,rowind), v);
                  end;
                end; 
             
                %--- this block of TB must be upper triangular (?) or zero
                
                for ell=0:k-1
                   if abs(TB(i+ell,i+ell)) < sqrt(eps)  
                        D = [D; inf];
                   else D = [D; T(i+ell,i+ell)/TB(i+ell,i+ell)]; 
                   end;
                end;
               
                i = i + k;
             end;   

          %#####################################################   
          %-------- complex eigenvalue pair --------------------
          %#####################################################   
          else
             if i == n
                error('Incorrect complex Schur form');
             end; 
                
             [U,xx,xx] = qr ( [real(Q(:,i:i+1)) imag(Q(:,i:i+1))]' );
             Q1 = j*U(1:2,3:4) + U(3:4,3:4);
             Q(:,i:i+1) = Q(:,i:i+1) * Q1;
             
             if nargin > 2
                [U,xx] = qr ( [real(Z(:,i:i+1)) imag(Z(:,i:i+1))]' );
                Z1 = j*U(1:2,3:4) + U(3:4,3:4);
                Z(:,i:i+1) = Z(:,i:i+1) * Z1;
             else
                Z1 = Q1;
             end;

             T(:,i:i+1) = T(:,i:i+1) * Z1;
             T(i:i+1,:) = Q1' * T(i:i+1,:);

             if nargin > 2
                TB(:,i:i+1) = TB(:,i:i+1) * Z1;
                TB(i:i+1,:) = Q1' * TB(i:i+1,:);

                %--- make TB upper triangular by Givens rotation
                 
                G = cgivens ( TB(i,i), TB(i+1,i) );
                T(i:i+1,:)  = G * T(i:i+1,:);
                TB(i:i+1,:) = G * TB(i:i+1,:);
                Q(:,i:i+1) = Q(:,i:i+1) * G';
             end;
             
             if abs(det(TB(i:i+1,i:i+1))) < 1.e-10
                  D = [D; inf; inf];
             else D = [D; eig(T(i:i+1,i:i+1)/TB(i:i+1,i:i+1))]; 
             end;

             i = i + 2;
          end           
        end 
%------------------------------------------------------
%       Finalization
%------------------------------------------------------
        Q = real(Q);
        Z = real(Z);
        T = hessform(real(T));
        TB = triang(real(TB));        
        if nargin < 4, 
          Z = T; T = D; 
        end;

%------- End of CSF2RSF.M --------- KYuP ----------           
