classdef stpbc
    properties
        a
        b
        c
        d
        om
        ups
    end
    methods
        function sys = stpbc(A,B,C,D,Om,Ups)
            %------------------------------------------------------
            %   Check syntax
            %------------------------------------------------------
            if ~exist('A','var')
                A = []; 
                B = zeros(0,1); 
                C = zeros(1,0); 
                D = 0;
                Om = []; 
                Ups = [];
            end
            if ~exist('B','var')
                D = A;
                A = []; 
                B = zeros(0,1); 
                C = zeros(1,0);
                Om = []; Ups = [];
            end
            n = size(A,1);
            m = size(B,2);
            p = size(C,1);
            if ~exist('D','var')
                D = zeros(p,m); 
            end
            if ~exist('Ups','var')
                Ups = zeros(n); 
            end
            if ~exist('Om','var')
                Om = eye(n); 
            end
            %------------------------------------------------------
            %   Check correctness
            %------------------------------------------------------
            error(abcdchk(A,B,C,D));
            if ~isnumeric(Om)
                error('Omega must be numeric');
            end;
            if ~isnumeric(Ups)
                error('Omega must be numeric');
            end;
            [rr,cc] = size(Om);
            if (rr ~= n)  ||  (cc ~= n),
                error('Dimensions of A and Omega are not equal');
            end;
            [rr,cc] = size(Ups);
            if (rr ~= n)  ||  (cc ~= n),
                error('Dimensions of A and Upsilon are not equal');
            end
            %------------------------------------------------------
            %   Form a new STPBC
            %------------------------------------------------------
            sys.a = A;
            sys.b = B;
            sys.c = C;
            sys.d = D;
            sys.om = Om;
            sys.ups = Ups;
        end
        
    end
end