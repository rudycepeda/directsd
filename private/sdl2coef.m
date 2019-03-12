function [A,B,E,A0,A1] = sdl2coef ( sys, T, H )
%SDL2COEF Coefficients for L2-optimal sampled-data system design.
%
%     [A,B,E,A0,A1] = SDL2COEF ( SYS, T )
%     [A,B,E,A0,A1] = SDL2COEF ( SYS, T, H )
%
%   Inputs:
%     SYS - an LTI system
%     T   - sampling period
%     H   - generalized hold as an LTI model 
%
%   Outputs:
%      A, B, E - function coefficients in quadratic functionsl
%      A0, A1  - functions such that A = A0 * A1 
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
        if i1 ~= 1
            error('L2-problem is formulated for a single input')
        end
        o1 = nout - o2;
        if o1 < 1
            error('Incorrect number of outputs'); 
        end
%------------------------------------------------------
%       Extrapolator
%------------------------------------------------------
        if ~exist('H', 'var') 
           H = ss(0, 1, 1, 0); % scalar ZOH by default
        end
%------------------------------------------------------
%       Extract transfer matrices
%------------------------------------------------------
        P11 = sdzpk(sys(1:o1,1:i1));
        P12 = sdzpk(sys(1:o1,i1+1:nin));
        P21 = sdzpk(sys(o1+1:nout,1:i1));
        %P22 = sys(o1+1:nout,i1+1:nin);
%------------------------------------------------------
%       Construct A0
%------------------------------------------------------
        scale12 = norm(P12.k);
        if scale12 ~= 0
             P12 = P12 / scale12;
        else
            scale12 = 1;
        end
        A0 = dtfm2(P12, T, H);
        p12 = allpoles(P12);
        A0 = scale12^2*setpoles(sdzpk(A0), exp(-[p12; -p12]*T) );
%------------------------------------------------------
%       Construct A1
%------------------------------------------------------
        scale21 = norm(P21.k); 
        if scale21 ~= 0
             P21 = P21 / scale21;
        else
            scale21 = 1;
        end
        A1 = z2zeta(sdzpk(ztrm(P21, T, 0)));
        p21 = allpoles(P21);
        A1 = scale21*setpoles(A1, exp(-p21*T) );
%------------------------------------------------------
%       Construct A
%------------------------------------------------------
        A = minreal(A0*(A1*A1'));        
%------------------------------------------------------
%       Construct B
%------------------------------------------------------
        scale11 = norm(P11.k);
        if scale11 ~= 0
             P11 = P11 / scale11;
        else
            scale11 = 1;
        end
        P11H = sdzpk(P11)';
        try 
          Px = minreal(P11H*P12);
          B1 = z2zeta(dtfm(Px, T, 0, H));
          px = allpoles(Px);
          B1 = setpoles(B1, exp(-px*T) );
        catch
          B1 = 0;  
          for i=1:size(P11H,2)
            Pxi = minreal(P11H(i)*P12(i));  
            B1i = z2zeta(dtfm(Pxi, T, 0, H));
            px = allpoles(Pxi);
            B1i = setpoles(B1i, exp(-px*T) );
            B1 = sumzpk(B1, B1i);
          end
        end
        B = minreal(A1*B1*scale11*scale12);
%------------------------------------------------------
%       Construct spectrum of A1
%------------------------------------------------------
        A1 = minreal(A1*A1');
%------------------------------------------------------
%       Construct E
%------------------------------------------------------
        if nargout > 2
            temp = P11H * P11;
            nt = size (temp,1);
            SP11 = temp(1,1); 
            for jj=1:nt
                SP11 = SP11+temp(jj);
            end
          %SP11 = trace(P11H * P11);
          E = ztrm(SP11, T, 0);
          p11 = allpoles(SP11);
          E = scale11^2*setpoles(sdzpk(E), exp(p11*T) );
	      FE = sfactor(E);
          FET=FE';
          FET.k=real(FET.k);
          E = FE*FET;
        end

%------- End of SDL2COEF.M --------- KYuP ----------           
