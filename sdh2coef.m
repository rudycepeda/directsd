function [A,B,E,A0,A1] = sdh2coef ( sys, T, t, H )
%SDH2COEF Coefficients for H2-optimal sampled-data system design.
%
%     [A,B,E,A0,A1] = SDH2COEF ( SYS, T )
%     [A,B,E,A0,A1] = SDH2COEF ( SYS, T, t )
%     [A,B,E,A0,A1] = SDH2COEF ( SYS, T, [], H )
%     [A,B,E,A0,A1] = SDH2COEF ( SYS, T, t, H )
%
%   Inputs:
%     SYS - an LTI system
%     T   - sampling period
%     t   - array time instants for computing periodic error
%     H   - generalized hold as an LTI model
%
%   Outputs:
%     A, B, E - function coefficients in quadratic functionsl
%     A0, A1  - functions such that A = A0 * A1
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
%       Extrapolator
%------------------------------------------------------
if ~exist('H', 'var'),
    H = ss(0, 1, 1, 0); % scalar ZOH by default
end;
%------------------------------------------------------
%       Extract transfer matrices
%------------------------------------------------------
P11 = sdzpk(sys(1:o1,1:i1));
P12 = sdzpk(sys(1:o1,i1+1:nin));
P21 = sdzpk(sys(o1+1:nout,1:i1));
%P22 = sys(o1+1:nout,i1+1:nin);
P11 = sdzpk(P11);
P11H = P11';
p12 = allpoles(P12);
p21 = allpoles(P21);
%------------------------------------------------------
%       Construct A1
%------------------------------------------------------
if hasdelay(P21)
    P212 = P21;
    P212.iodelay = 0;
    P212.inputdelay = 0;
    P212.outputdelay = 0;
    A1 = ztrm(minreal(P212*P212'), T, 0);
else
    A1 = ztrm(minreal(P21*P21'), T, 0);
end
A1 = setpoles(sdzpk(A1), exp(-[p21; -p21]*T) );
A1 = symmetr(sdzpk(A1), 'd');
%-----------------------------------------------
%       Average variance
%-----------------------------------------------
if ~exist('t','var')  ||  isempty(t)
    %--------------------------------------------
    %       Construct A0
    %--------------------------------------------
    A0 = dtfm2(P12, T, H) / T;
    A0 = setpoles(sdzpk(A0), exp(-[p12; -p12]*T) );
    A0 = symmetr(A0, 'd');
    A = minreals(sdzpk(A0*A1));
    %--------------------------------------------
    %       Construct B
    %--------------------------------------------
    try
        Px = minreal(P21*P11H*P12);
        B = z2zeta(sdzpk(dtfm(Px, T, 0, H)))/T;
        px = allpoles(Px);
        B = setpoles(B, exp(-px*T) );
    catch
        B = 0;
        for i=1:size(P11H,2)
            for j=1:size(P21,2)
                Px = minreal(P21(j)*P11H(j,i)*P12(i));
                Bij = z2zeta(sdzpk(dtfm(Px, T, 0, H))) / T;
                px = allpoles(Px);
                Bij = setpoles(Bij, exp(-px*T) );
                B = sumzpk(B, Bij);
            end;
        end;
    end;
    %-----------------------------------------------
    %       Instantaneous variance
    %-----------------------------------------------
else
    for i=1:length(t)
        if t(i) >= T, t(i) = T - sqrt(eps)*T; end;
        %--------------------------------------------
        %       Construct A0
        %--------------------------------------------
        A0 = z2zeta(sdzpk(dtfm(P12, T, t(i), H)));
        A0 = setpoles(A0, exp(-p12*T) );
        A(i) = minreal(A0'*A0*A1);
        %--------------------------------------------
        %       Construct B
        %--------------------------------------------
        Bs = (minreal(P21*P11H)').';
        pBs = allpoles(Bs);
        Bx = ztrm(Bs, T, t(i));
        Bx = setpoles(sdzpk(Bx), exp(pBs*T) );
        B(i) = zpk([], [], 0, T);
        for k=1:max(size(Bx))
            B(i) = sumzpk(B(i), Bx(k)*A0(k));
        end;
    end;
end;
%------------------------------------------------------
%       Construct E
%------------------------------------------------------
if nargout > 2 %Modified by RCG to avoid using "trace"
    temp = P11H * P11;
    nt = size(temp,1);
    SP11 = temp(1,1);
    for jj = 2:nt
        SP11 = SP11 + temp(jj,jj);
    end
    E = ztrm(SP11, T, 0);
    p11 = allpoles(SP11);
    E = setpoles(sdzpk(E), exp(p11*T) );
    FE = sfactor(sdzpk(E));
    E = FE*FE';
end;

%------- End of SDH2COEF.M --------- KYuP ----------
