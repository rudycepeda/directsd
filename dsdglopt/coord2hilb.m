function  x = coord2hilb ( coord, precision )
%COORD2HILB Map an ND-vector to 1D-value in [0,1] using Hilbert curve.
%
%	  X = COORD2HILB ( ALPHA, BITS )
%
%   Inputs:
%	  ALPHA - real vector of parameters
%     BITS  - number of bits in a value 
%
%   Outputs:
%	  X    - scalar real value in [0,1]
%
%   See also HILB2COORD.
	
% [1] J.K.Lawder. Calculation of Mappings Between One and n-dimensional 
%     Values Using the Hilbert Space-filling Curve. Research Report 
%     BBKCS-00-01 (formerly JL1/00), Birkbeck College, University of
%     London, August 2000.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Find dimensions
%------------------------------------------------------
        dimensions = length(coord); 
%-------------------------------------------------
%  Transform to binary values
%-------------------------------------------------
        scoord = char(48+zeros(dimensions,precision));
        for i=1:dimensions
          [xx,s] = val2bin(coord(i), precision);  
          scoord(i,:) = s;
        end;
%-------------------------------------------------
%  Construct alpha by transposition
%-------------------------------------------------
        alpha = zeros(precision,dimensions);    
        for i=1:dimensions
          alpha(:,i) = str2num(scoord(i,:)'); 
        end;        
%-------------------------------------------------
%  Initialization
%-------------------------------------------------
        Z = zeros(precision, dimensions);
        omega = Z; rho = Z;  
        sigma = Z; tilde_sigma = Z;  
        tau = Z;   tilde_tau = Z;  
        J = zeros(1,dimensions); shift = J;
%#################################################
%  Main loop
%#################################################
        for i=1:precision
          if i > 1,
            omega(i,:) = xor(omega(i-1,:), tilde_tau(i-1,:));  
          end;
         %-------------------------------------------------
         %  Construct tilde_sigma = alpha xor omega
         %-------------------------------------------------
          tilde_sigma(i,:) = xor(alpha(i,:), omega(i,:));  
         %-------------------------------------------------
         %  Construct sigma by cyclic shifing tilde_sigma
         %  to the left by sum(J-1) bits
         %-------------------------------------------------
          sig = [tilde_sigma(i,:) tilde_sigma(i,:)];
          k = shift(i);
          sigma(i,:) = sig(k+1:k+dimensions);
         %-------------------------------------------------
         %  Construct rho:
         %    rho(1) = sigma(1)
         %    rho(i) = sigma(i) xor rho(i-1), i > 1
         %-------------------------------------------------
          rho(i,1) = sigma(i,1);
          for j=2:dimensions
            rho(i,j) = xor(sigma(i,j), rho(i,j-1));             
          end;
         %-------------------------------------------------
         %  Find J(i) - a principal position
         %-------------------------------------------------
          rhoi = rho(i,:);  
		  ind = find(rhoi ~= rhoi(end));
          if isempty(ind)
               J(i) = dimensions;
          else J(i) = max(ind); 
          end;
          shift(i+1) = mod(shift(i)+J(i)-1, dimensions);
         %-------------------------------------------------
         %  Construct tau:
         %  1) invert the last bit of sigma
         %  2) if odd parity, invert principal bit
         %-------------------------------------------------
          tau(i,:) = sigma(i,:);
          tau(i,dimensions) = ~tau(i,dimensions);
          if mod(sum(tau(i,:)),2) == 1
            tau(i,J(i)) = ~tau(i,J(i));              
          end;
         %-------------------------------------------------
         %  Construct tilde_tau by cyclic shifing tau
         %  to the left by sum(J-1) bits
         %-------------------------------------------------
          titau = [tau(i,:) tau(i,:)];
          k = dimensions - shift(i);
		  tilde_tau(i,:) = titau(k+1:k+dimensions);
        end;
%#################################################
%  End of Main loop
%#################################################
        r = rho';
        r = char(r(:)'+48);
        if length(r) > 52,  % allowable precision is 52 bits
           r = r(1:52);
        end;
        x = bin2val('0', r); 
        
%------- End of COORD2HILB.M --------- KYuP ----------
