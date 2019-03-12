function err = quaderr ( A, B, E, M )
%QUADERR Squared integral error for sampled-data system.
%
%     ERR = QUADERR ( A, B, E, M )
%
%   Inputs:
%     A, B, E - coefficients of the integrand
%     M = K / (1 + D22*K), where K is controller
%
%   Outputs:
%     ERR - integral quadratic error or variance
%
%   See also POLQUAD, POLHINF.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov
% $Revision: 3.00 $    $Date: 08-Oct-2006 $
%------------------------------------------------------
%       Bilinear transformation
%------------------------------------------------------
b = 1;
BL = fbilin ( b );
Mb = bilintr (M, b);
Mb = stripzpk( minreal(Mb, 1e-4) );
Eb = minreal(bilintr(E,b)*BL);
%------------------------------------------------------
%       Loop over all t's
%------------------------------------------------------
err = [];
for i=1:max(size(A))
    Ab = bilintr (A(i), b);
    Bb = bilintr (B(i), b);
    %--------------------------------------------
    %       Construct A*M*M' and B*M
    %--------------------------------------------
    AMM  = minreals(sdzpk(Ab*(Mb*Mb')*BL));
    BM   = minreal (Bb*Mb*BL);
    %--------------------------------------------
    %       Construct full integrand
    %--------------------------------------------
    X1 = sumzpk(AMM, Eb);
    X2 = sumzpk(BM, BM');
    X  = sumzpk(X1, -X2);
    saveK = X.K; X.K = 1;
    X = sumzpk(X, X') / 2;
    X.K = X.K * saveK;
    %--------------------------------------------
    %       Factorization
    %--------------------------------------------
    [~,X] = improper ( X, 'symm' );
    try 
        FX = sfactor ( X );
    catch
        keyboard;
    end
    %[xx,FX] = improper ( FX );
    err = [err; norm(FX)^2]; %#ok<AGROW>
end

%------- End of QUADERR.M --------- KYuP ----------

