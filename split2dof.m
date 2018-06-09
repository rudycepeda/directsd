function [K,KR,KC] = split2dof ( K, KR, objClass )
%SPLIT2DOF Stable realization of 2-DOF sampled-data systems.
%
%     [KF,KR,KC] = SPLIT2DOF ( K, KR )
%
%   Inputs:
%     K(z)  - full feedback controller
%     KR(z) - full reference controller
%
%   Outputs:
%     KF(z) - feedback controller
%     KR(z) - reference controller
%     KC(z) - common part of both the controllers
%
%   See also SD2DOF, SD2DOFERR.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%---------------------------------------------------
%       Check data
%---------------------------------------------------
        if ~exist('objClass','var'),
           objClass = class ( K );
        end;       
        K = zpk(K);  
        pK = K.p{1};
        KR = zpk(KR); 
        pKR = KR.p{1};
%---------------------------------------------------
%       Separate all common poles of K and KR
%         pK  = pCommon*pKa
%         pKR = pCommon*pKRa
%---------------------------------------------------
        %[pKa,pCommon]  = others(pK,pKR,1e-5);
        [pKRa,pCommon] = others(pKR,pK,1e-5);
%---------------------------------------------------
%       If deg(dCommon) = n, then
%           KC    = z^n / dCommon
%           CFeedback  = nK / (z^n * dKa)
%           CReference = nKR / (z^n * dKRa)
%---------------------------------------------------
        n = length(pCommon);
        KC = zpk(zeros(n,1), pCommon, 1, K.Ts);
        K = minreal( K / KC, 1e-4 );
        KR = minreal( KR / KC, 1e-4 );        
%---------------------------------------------------
%       Finalization
%---------------------------------------------------
        switch (objClass)
          case 'tf', K = tf(K); KR = tf(KR); KC = tf(KC);
          case 'ss', K = ss(K); KR = ss(KR); KC = ss(KC);
        end;            
        if ~isempty(find(abs(pKRa)>=1)),
           error('Reference controler cannot be unstable');
        end;

%------- End of SPLIT2DOF.M --------- KYuP ----------
