function  coord = hilb2coord ( x, dimensions, precision )
%HILB2COORD Map an 1D-value in [0,1] to ND-vector using Hilbert curve.
%
%	  ALPHA = HILB2COORD ( X, DIM, BITS )
%
%   Inputs:
%	  X    - scalar real value in [0,1]
%     DIM  - number of dimensions
%     BITS - number of bits in a value 
%
%   Outputs:
%	  ALPHA - vector having DIM elements
%
%   See also COORD2HILB.

% [1] A. R. Butz, Alternative Algorithm for Hilbert's Space-Filling Curve,
%     IEEE Trans. Comp., April, 1971, pp 424-426.
% [2] J.K.Lawder. Calculation of Mappings Between One and n-dimensional 
%     Values Using the Hilbert Space-filling Curve. Research Report 
%     BBKCS-00-01 (formerly JL1/00), Birkbeck College, University of
%     London, August 2000.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Get binary string
%------------------------------------------------------
%        x = bin2val('0', '10011000100010111000');
%        dimensions = 5;
%        precision = 4;
        [xx,r] = val2bin ( x, dimensions*precision );
        r = xor(r-48, 0);
%------------------------------------------------------
%       Find 'ro' by splitting 'r' onto chunks of
%       length 'precision' and adding a zero chunk
%------------------------------------------------------
        rho = reshape(r, dimensions, precision )';
%------------------------------------------------------
%       Find principal position J for each 'ro' 
%------------------------------------------------------
        J = zeros(1,precision);
        for i=1:precision
          rhoi = rho(i,:);  
		  ind = find(rhoi ~= rhoi(end));
          if isempty(ind)
               J(i) = dimensions;
          else J(i) = max(ind); 
          end;
        end;
%------------------------------------------------------
%       Build 'sigma' array: improved method by J.K. Lawder
%       sigma = rho xor (rho >> 1) 
%------------------------------------------------------
        sigma = rho;
        bitmask = eye(dimensions);
        for i=1:precision
          sig = xor([rho(i,:) 0], [0 rho(i,:)]);
          sigma(i,:) = sig(1:dimensions);
        end;        
%------------------------------------------------------
%       Build 'tau' array
%       Complement the last bit and then, if no parity,
%       complement the principal bit
%------------------------------------------------------
        tau = sigma;
        tau(:,dimensions) = ~tau(:,dimensions); 
	    for i=1:precision
          if mod(sum(tau(i,:)),2) == 1, 
             tau(i,J(i)) = ~tau(i,J(i)); 
          end;
        end;
%------------------------------------------------------
%       Build array of shifts: 
%          shift(i) = sum_1^(i-1) (J(i) - 1)
%------------------------------------------------------
        shift = zeros(1,dimensions);
	    for i=2:precision
		  shift(i) = mod(shift(i-1)+J(i-1)-1, dimensions);
        end;
%------------------------------------------------------
%       Build 'tilde sigma' array by cyclic shifing sigma
%       to the right by 'shift' bits
%------------------------------------------------------
        tilde_sigma = sigma;
	    for i=2:precision   % tilde_sigma(1,:) = sigma(1,:)
          tisig = [sigma(i,:) sigma(i,:)];
          k = dimensions - shift(i);
		  tilde_sigma(i,:) = tisig(k+1:k+dimensions);
        end;
%------------------------------------------------------
%       Build 'tilde tau' array by cyclic shifing tau
%       to the right by 'shift' bits
%------------------------------------------------------
	    tilde_tau = tau;
	    for i=2:precision   % tilde_tau(1,:) = tau(1,:)
          titau = [tau(i,:) tau(i,:)];
          k = dimensions - shift(i);
		  tilde_tau(i,:) = titau(k+1:k+dimensions);
        end;        
%------------------------------------------------------
%       Build 'omega' array:
%        omega(i) = omega(i-1) xor tilde_tau(i-1)
%------------------------------------------------------
	    omega = zeros(1,dimensions);
	    for i=2:precision
	       omega(i,:) = xor(omega(i-1,:), tilde_tau(i-1,:));
        end;
%------------------------------------------------------
%       Build 'alpha' array:
%        alpha(i) = omega(i) xor tilde_sigma(i)
%------------------------------------------------------
        alpha = omega;
	    for i=1:precision
		   alpha(i,:) = xor(omega(i,:), tilde_sigma(i,:));        
        end
%------------------------------------------------------
%       Convert 'alpha' array into nD-vector
%------------------------------------------------------
        coord = zeros(dimensions,1);
	    for i = 1:dimensions
		  scoord(i,:) = char(alpha(:,i)'+48);
          coord(i) = bin2val('0',scoord(i,:));
        end
        %rho, J, sigma, tau, tilde_sigma, tilde_tau,
        %omega, alpha, scoord
        %keyboard
        
%------- End of HILB2COORD.M --------- KYuP ----------
