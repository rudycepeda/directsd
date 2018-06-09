function [b1,b2,c1,c2,d11,d12,d21,d22] = blocks4 ( b, c, d, o2, i2 )
%BLOCKS4 Extract 4-block model from a state-space data.
%
%     [B1,B2,C1,C2,D11,D12,D21,D22] = BLOCKS4 ( B, C, D, NMEAS, NCON )
%
%   Inputs:
%     B, C, D - state-space constant matrices 
%     NMEAS - dimension of the vector 'y' 
%     NCON  - dimension of the vector 'u' 
%
%   Outputs:
%     B1, B2, C1, C2, D11, D12, D21, D22 - matrices of the 4-block model
%                                w    u  
%                           | A  B1  B2  |
%           P = | A B | = z | C1 D11 D12 |
%               | C D |   y | C2 D21 D22 |
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
        nin = size(b, 2);
        nout = size(c, 1);   
        i1 = nin - i2;
        o1 = nout - o2;
    
        b1 = b ( :, 1:i1 );
    	b2 = b ( :, (i1+1):nin );
    	c1 = c ( 1:o1, : );
    	c2 = c ( (o1+1):nout, : );
    	d11 = d ( 1:o1, 1:i1 );
    	d12 = d ( 1:o1, (i1+1):nin );
    	d21 = d ( (o1+1):nout, 1:i1 );
    	d22 = d ( (o1+1):nout, (i1+1):nin);

%------- End of BLOCKS4.M --------- KYuP ----------
