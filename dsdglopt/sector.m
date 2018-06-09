function [alpha,beta] = sector ( p, dom, T )
%SECTOR Computes stability degree abs relative oscillation.
%
%    [ALPHA,BETA] = SECTOR ( P, D, T )
%
%  Inputs:
%    P - roots of the characteristic equation
%    D - domain of 's' (default), 'z', or 'd'
%    T - sampling period (default 1 for discrete-time system)
%
%  Outputs:
%    ALPHA - stability degree
%    BETA  - oscillation degree for each of the poles
%
%   See also MODSDH2, MODSDL2.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check errors
%------------------------------------------------------
        if ~exist('dom','var'), dom = 's'; end; 
        if ~exist('T','var'), T = 1; end;
%------------------------------------------------------
%       Transform poles to s-plane
%------------------------------------------------------
        switch dom
          case 's', 
          case 'z', p = log(p) / T;
          case 'd', p = - log(p) / T;
          otherwise, error('Incorrect domain ''%s''',dom);
        end;
%------------------------------------------------------
%       Find oscillation degree
%------------------------------------------------------
        alpha = - real(p);        
        beta  = abs( imag(p) ./ real(p) );        

%------- End of OSCDEG.M --------- KYuP ---------------
