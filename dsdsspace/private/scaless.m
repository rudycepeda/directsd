function [b1,b2,c1,c2,d11,d12,d21,d22,su1,su2,sy1,sy2] = scaless (b1,b2,c1,c2,d11,d12,d21,d22)
%SCALESS Scaling of standard model with respect to D12 and D21.
%
%     [B1,B2,C1,C2,D11,D12,D21,D22,SU1,SU2,SY1,SY2] = SCALESS (B1,B2,C1,C2,D11,D12,D21,D22)
%
%   Inputs:
%     B1,B2,C1,C2,D11,D12,D21,D22 - state-space constant matrices of 4-block model
%
%   Outputs:
%     B1,B2,C1,C2,D11,D12,D21,D22 - matrices of the 4-block model
%     SU1,SU2,SY1,SY2 - scale factors
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%-----------------------------------------------------
%   Scaling factors for D12
%---------------------------------------------------------
    i1 = size(b1,2);
    i2 = size(b2,2);
    o1 = size(c1,1);
    o2 = size(c2,1);

    [u12,s12,v12] = svd(d12);

    u1 = u12 (:,1:i2);
    u2 = u12 (:, (i2+1):o1);
    su2 = v12 * inv( s12(1:i2,1:i2) );
    sy1 = [u2';u1'];
%---------------------------------------------------------
%   Scaling factors for D21
%---------------------------------------------------------
    [u21,s21,v21] = svd(d21);
    v1 = v21 (:, 1:o2);
    v2 = v21 (:, (o2+1):i1);
    sy2 = inv ( s21(1:o2,1:o2) ) * u21';
    su1 = [v2 v1];
%---------------------------------------------------------
%   Scale all matrices 
%---------------------------------------------------------
    b1 = b1 * su1;
    b2 = b2 * su2;
    c1 = sy1 * c1;
    c2 = sy2 * c2;
    d11 = sy1 * d11 * su1;
    d12 = sy1 * d12 * su2;
    d21 = sy2 * d21 * su1;

%------- End of SCALESS.M --------- KYuP ----------

