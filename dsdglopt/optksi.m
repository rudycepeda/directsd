function [ksi,aDelta,bDelta] = optksi ( data, Delta, degKsi, gKsi, aDelta, bDelta )
%OPTKSI Find optimal 'ksi' polynomial in parameterization.
%
%	  [KSI,ADELTA,BDELTA] = OPTKSI ( DATA, DELTA, DEG, GAIN )
%
%   Inputs:
%     DATA - structure with the following fields
%         .A, .B - coefficients in the functional, in 'zeta'
%         .plant - discrete model of the plant, in 'z' 
%         .dK0   - mandatory part of controller denominator, 
%                  for example [1 -1] to include integrator   
%     DELTA - characteristic polynomial  
%     DEG   - degree of the 'ksi'-polynomial  
%     GAIN  - 'static gain' of the 'ksi'-polynomial (NaN if any) 
%
%   Outputs:
%     KSI - 'ksi'-polynomial
%     ADELTA, BDELTA - minimal controller associated with Delta 
%
%   See also GO_PAR2K.
	
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 26-Dec-2006 $
%------------------------------------------------------
%       Check params
%------------------------------------------------------
    if ~exist('gKsi','var'), gKsi = NaN; end;
    if ~isfield(data,'dK0'), data.dK0 = 1; end;
%------------------------------------------------------
%       Find aDelta and bDelta if not given
%------------------------------------------------------
    if ~exist('aDelta','var')  ||  ~exist('bDelta','var')
      [aDelta,bDelta,data] = abdelta ( data, Delta );  
    end;
%------------------------------------------------------
%       Special case: degKsi < 0, no choice       
%------------------------------------------------------
    if degKsi < 0
      if ~isnan(gKsi), 
           ksi = gKsi;
      else ksi = 0; end;       
      return;
    end;
%------------------------------------------------------
%       Find optimal ksi using polynomial optimization
%------------------------------------------------------
    T = data.plant.Ts;
    D22 = zpk(data.plant)';
    dZpk = zpk(D22.p{1}, [], 1, T);
    dDelta = dZpk * zpk(tf(1,Delta,T));
    d2Delta = minreal(dZpk*dDelta);
    aDeltaZpk = zpk(tf(aDelta, 1, T));
    dK0Zpk = zpk(tf(data.dK0, 1, T));
   %-----------------------------------------------
   %  Find Av abd Bv
   %-----------------------------------------------
    tol = 1e-4;
    Av2 = d2Delta*d2Delta';
    Av3 = dK0Zpk*dK0Zpk';
    Bv1 = minreal(d2Delta*data.B, tol);
    Bv2 = minreal(d2Delta*dDelta'*aDeltaZpk'*data.A, tol);
    Bv = sumzpk(Bv1, -Bv2);
    Bv = minreal(Bv*dK0Zpk, tol);
   %-----------------------------------------------
   %  Find Lam and Lplus
   %-----------------------------------------------
    Lam1 = minreal(sfactor(data.A), tol);
    Lam2 = minreal(sfactor(Av2), tol);
    Lam3 = minreal(sfactor(Av3), tol);
    Lam  = minreal(Lam1*Lam2*Lam3, tol);
   %-----------------------------------------------
   %  Separation
   %-----------------------------------------------
    L = minreal(Bv'/Lam', tol);
   %-----------------------------------------------
   %  Taking account of static gain
   %-----------------------------------------------
    if ~isnan(gKsi), 
      L = sumzpk(L, -gKsi*Lam);  
      Lam = minreal(Lam*zpk(1,[],1,Lam.Ts),tol);  
    end;
   %-----------------------------------------------
   %  Extract stable part in separation
   %-----------------------------------------------
    Lplus = separtf ( L, 'd' );         
%------------------------------------------------------
%  Find optimal polynomial ksi
%------------------------------------------------------
    ksi = polopth2 ( Lplus, Lam, degKsi );
%------------------------------------------------------
%       Preserve static gain
%------------------------------------------------------
    if ~isnan(gKsi), 
      ksi = sumpol(gKsi, conv(ksi,[1 -1]));
    end;
        
%------- End of OPTKSI.M --------- KYuP ----------
