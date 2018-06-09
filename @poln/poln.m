classdef poln
    %POLN constructs a polynomial object
    %
    %   P = POLN(A) creates a polynomial with the coefficients given in the
    %   vector A.
    %
    %   P = POLN(A,'r') creates a polynomial of which the roots are given in
    %   the vector A.
    %
    %   P = POLN(A,VAR) where VAR is 's', 'p', 'd', 'z', or 'q' creates a
    %   polynomial in the variable VAR with coefficients given in the vector
    %   A. If VAR is a string like 'rs' or 'rz', the roots of the polynomial
    %   are taken from vector A.
    %
    %   P = POLN(A,VAR,SHIFT) allows the construction of discrete
    %   quasi-polynomials. The argument SHIFT must be positive, it indicates
    %   the least negative degree of the quasi-polynomial. This can be used
    %   only with variables d, z, or q.
    
    %--------------------------------------------------------------------------
    % Original methods written by Dr. K. Yu. Polyakov.
    % Revised and updated by Rudy Cepeda-Gomez
    % Revision: 4.0
    % Date: 15-Jul-2014
    %--------------------------------------------------------------------------
    
    properties
        coef; %Coefficients of the Polynomial
        k = 1; %Gain: coefficient of the largest power
        z; %Roots of the Polynomial
        var = 's'; %Variable of the the polynomial
        shift = 0; %Shift for discrete time polynomials
    end
    
    methods (Access = public)
        %Basic Methods: construction, display and refrencing
        function p = poln(varargin)
            %Class constructor for polynomial objects
            if nargin ==  0
                p.coef = [];
                p.k = [];
                p.z =[];
                p.var =[];
                p.shift = [];
                return
            end
            A = varargin{1};%First argument
            RootMode = 0;%Flag to indicate if roots were provided.
            %Check input data
            %Check class of first argument, must be a double or a
            %polynomial if a property is to be reassigned
            if ~isa(A,'double')
                if isa(A,'poln')
                    %If  input is a poln, copy the current properties
                    p.coef = A.coef;
                    p.z = A.z;
                    p.k = A.k;
                    p.var = A.var;
                    p.shift = A.shift;
                else
                    error('Coefficients or roots must be of class double');
                end
            end
            %More input arguments, default variable and shif not to be used
            if nargin > 1
                %Check variable and root mode
                variable = varargin{2};
                switch length(variable)
                    case 1
                        if variable =='r';
                            RootMode = 1;%Roots are provided
                            variable = 's';%Default variable
                        end
                    case 2
                        if variable(1) == 'r'
                            RootMode = 1;%Roots and variable are provided
                            variable = variable(2);
                        else%Unknown first letter for var
                            error(['Unknown option ' variable]);
                        end
                    otherwise
                        error(['Unknown option ' variable]);
                end
                %Check if a valid variable is requested
                if isempty(find(variable=='spdzq',1))
                    error('Variable must be s, p, d, z or q');
                end
                %If variable and shift are provided
                if nargin == 3
                    sh = varargin{3};
                    if sh < 0 %Only positive shifts are allowed
                        error('Parameter ''shift'' must be nonnegative')
                    elseif isempty(find(variable=='dzq',1)) && sh~=0 %Shift is used only in discrete time
                        error('Parameter ''shift'' cannot be used with continuous-time polynomials');
                    end
                end
            end
            if nargin > 3
                error('Too many input arguments')
            end
            %Create the object
            %Roots and coefficients
            if isempty(p.coef)
                if ~RootMode
                    A = striplz(A);
                    p.coef = A;
                    p.z = roots(A);
                    p.k = A(1);
                else
                    p.z = A;
                    p.coef = poly(A);
                end
            end
            %If variable and shift are provided:
            if nargin > 1
                p.var = variable;
                if isct(p) && p.shift~=0
                    warning('Converting to continuous time, ignoring shift')
                    p.shift = 0;
                end
                if nargin ==3
                    p.shift = sh;
                end
            end
            %End of constructor method-------------------------------------
        end
        function display(P)
            %Display method for class poln
            iname = inputname(1);
            if isempty(P.coef)
                disp([iname ' = []']);
                disp('    Empty polynomial object');
                return
            end
            [r,c] = size(P);
            for i = 1:r
                for j = 1:c
                    %Display header
                    if P(i,j).shift == 0
                        sname = 'polynomial';
                    else
                        sname = 'quasipolynomial';
                    end;
                    idx=[num2str(i) ',' num2str(j)];
                    if r*c > 1
                        disp([iname '(' idx '): ',sname]);
                    else
                        disp([iname ': ',sname]);
                    end
                    %Display Polynomial
                    margin = '  ';
                    disp([margin, char(P(i,j))]);
                end
            end
            % End of DISPLAY Method----------------------------------------
        end
        function R = subsref(P,S)
            %SUBSREF Subscripted refrence method for polynomial objects
            for i=1:length(S)
                switch S(i).type
                    case '.'
                        if strcmp(S(i).subs,'shift'),
                            R = P.shift;
                        elseif strcmp(S(i).subs,'deg'),
                            R = deg(P);
                        elseif strcmp(S(i).subs,'var'),
                            R = P.var;
                        elseif strcmp(S(i).subs,'coef'),
                            if i+1 > length(S)
                                R = P.coef;
                            else
                                Q.type = S(i+1).type;
                                Q.subs = S(i+1).subs;
                                if strcmp(Q.type,'()')
                                    R = subsref(P.coef, Q);
                                else
                                    error('Incorrect query');
                                end
                            end
                            return
                        elseif strcmp(S(i).subs,'z')
                            if i+1 > length(S)
                                R = P.z;
                            else
                                Q.type = S(i+1).type;
                                Q.subs = S(i+1).subs;
                                if strcmp(Q.type,'()')
                                    R = subsref(P.z,Q);
                                else
                                    error('Incorrect query')
                                end
                            end
                        elseif strcmp(S(i).subs,'k')
                            R = P.k;
                        else
                            error('Unknown property')
                        end
                    case '()'
                        idx = [S.subs{1:end}];
                        if numel(idx) > 1
                            if numel(idx) == 2
                                if idx(2) == 1
                                    idx = idx(1);
                                else
                                    error('Only column vectors are supported');
                                end
                            else
                                error('Only vectors are supported');
                            end
                        end
                        if idx > numel(P)
                            error('Index out of bounds');
                        end
                        R = P(idx);
                end
            end
            
            %End of SUBSREF method
        end
        function P = subsasgn(P,S,B)
            %SUBSASGN Subscripted assignment method for polynomial objects
            if ~isa(B,'poln')
                if isa(B,'double')
                    B = poln(B,'s');
                    if ~isempty(P)
                        B.var = P(1).var;
                    end
                else
                    msg=('In an assignment A(I)=B, B must be a polynomial');
                    msg=[msg ' or a vector of doubles'];
                    error(msg);
                end
            else
                if ~isempty(P) && (P(1).var ~= B.var)
                    error('Polynomials within the array must share the variable')
                end
            end
            switch S.type
                case '()'
                    idx = [S.subs{1:end}];
                    if length(idx) ~= 1
                        idx = idx(1);
                        warning('Trailing dimensions ignored');
                    end
                    if isempty(P)
                        P = poln(0,B.var);
                        for i = 2:idx(1);
                            P(i,1) = poln(0,B.var);
                        end
                    end
                    P(idx,1) = B;
                otherwise
                    error('Subscripted assingment of properties not supported');
            end
        end
        %Math
        function P = ctranspose(P)
            % ' Complex conjugate transpose
            %
            %   P' is the complex conjugate transpose of polynomial P. If P
            %   is an array of polynomials, P' is an array that contains
            %   the complex conjugate transpose of the elements in P.
            
            %Check size
            if length(P) ==1
                %Polynomials in continuous time
                if isct(P)
                    P.z = - P.z;
                    n = deg(P);
                    if mod(n,2)
                        P.k = -P.k;
                    end
                    c = fliplr(P.coef);
                    c(2:2:n+1) = - c(2:2:n+1);
                    P.coef = fliplr(c);
                    return;
                end
                %Polynomials in discrete time
                indZero = find(abs(P.z) < eps);
                P.z(indZero) = [];
                P.k = real(P.k*prod(-P.z));
                P.z = 1./P.z;
                P.shift = length(P.z) + length(indZero) - P.shift;
                if P.shift < 0
                    P.z = [P.z; zeros(-P.shift,1)];
                    P.shift = 0;
                end;
                P.coef = striplz(fliplr(P.coef),eps);
                add = P.shift + 1 - length(P.coef);
                if add > 0
                    P.coef = [zeros(1,add) P.coef];
                end
            else
                for i = 1:length(P)
                    P(i) = P(i)';
                end
            end
        end
        function N = norm(P)
            %NORM Euclidean norm of polynomial coefficients
            %
            %   N = NORM(P) is the euclidean norm of the vector of
            %   coefficients of the polynomial P.
            N = norm(P.coef);
        end
        function [degP,degM] = deg(P)
            %DEG Degree(s) of a polynomial
            %
            %   [DEGP,DEGM] = DEG(P) returns the largest DEGP, and the
            %   smallest DEGM, powers of polynomial P. DEGM is nonzero only
            %   for discrete quasipolynomials with nonzero shift.
            if norm(P.k) < eps,
                degP = 0;
                degM = 0;
            else
                degP = length(P.z) - P.shift;
                degM = P.shift;
            end;
            %End of DEG method -------------------------------------------
        end
        function Val = polyval(P,x)
            %POLYVAL Evaluates polynomial
            %
            %   VAL = POLYVAL(P,X) returns the value of the polynomial P
            %   when P.var=x.
            if P.shift > 0  &&  x == 0
                warning('Divide by zero');
                Val = inf;
            else
                Val = P.k*prod(x - P.z)/x^P.shift;
                if imag(x) < eps
                    Val = real(Val);
                end
            end
            %End of POLYVAL method ----------------------------------------
        end
        function P = derive(P,n,method)
            %DERIVE Derivative of a polynomial
            %
            %   P = DERIVATIVE(P) is the first derivative of polynomial P
            %   with respect to iits variable.
            %
            %   P = DERIVATIVE(P,n) is the n-th derivative of polynomial P
            %   with respect to its variable.
            %
            %   P = DERIVATIVE(P,n,method) allows the selection of a method
            %   for the computation of the derivative. 'roots', the default
            %   method, tries to preserve the value of the roots of the
            %   original polynomial. 'coef' operates on the coefficients
            %   and may introduce bigger numerical errors. The initial
            %   letter of each method is enough.
            
            %Check data
            if ~exist('n','var')
                n = 1;
            end;
            if ~exist('method','var')
                method = 'roots';
            end
            if ~any([isequal(method,'coef') ~isequal(method,'roots')...
                    isequal(method,'c') ~isequal(method,'r')])
                error(['Unknown method' method]);
            end
            method = method(1);
            if P.shift == 0%Polynomials
                if method == 'r'%Roots method: preserve roots if possible
                    for i=1:n
                        D = poln(0,P.var);
                        l = 1;
                        while l <= length(P.z)
                            zRem = P.z;
                            if imag(zRem(l)) < eps
                                zRem(l) = [];
                                l = l + 1;
                                D = D + poln(zRem, ['r' P.var]);
                            else
                                z0 = real(zRem(i));
                                zRem(l:l+1) = [];
                                l = l + 2;
                                D = D + 2*poln(z0,['r' P.var])*poln(zRem,['r' P.var]);
                            end
                        end
                        P = D*P.k;
                    end
                else %Direct method, via coefficients
                    c = P.coef;
                    for i = 1:n
                        powers = length(c)-1:-1:0;
                        c = c.*powers;
                        c = c(1:end-1);
                        if isempty(c)
                            c = 0;
                            break
                        end
                    end
                    P = poln(striplz(c),P.var);
                end
            else
                %----------------------------------------------------------
                % Quasipolynomial: P(z) = P0(z)*(1/z^m)
                %                 P'(z) = P0'(z)*(1/z^m) - m*P0*(1/z^(m+1))
                %                       = (P0'(z)*z - m*P0)*(1/z^(m+1))
                %----------------------------------------------------------
                m = P.shift;
                P.shift = 0;
                U = poln([1 0],P.var);
                for i=1:n
                    P = derive(P, 1, method)*U - m*P;
                    m = m + 1;
                end
                P.shift = m;
            end
            %End of DERIVE method ---------------------------------------------
        end
        function Val = polyder(P,x,N)
            %POLYDER Derivative of a polynomial at a given point
            %
            %   VAL = POLYDER(P,X) is the value of the first derivative of
            %   polynomial P evaluated at the point P.var = X.
            %
            %   VAL = POLYDER(P,X,N) evaluates the N-th derivative of P at
            %   the point P.var = X.
            
            if ~exist(N,'var')
                N = 1;
            end
            if P.shift > 0
                warning('Division by zero');
                Val = inf;
            else
                Val = polyval(derive(P,N),x);
            end
            %End of POLYDER method ----------------------------------------
        end
        function P = uminus(P)
            %- Unary minus
            %
            %   P = -P Negates the coeffients of P
            P.k=-P.k;
            P.coef=-P.coef;
        end
        function R = roots(P)
            %ROOTS Roots of a polynomial object
            %
            %   R = ROOTS(P) are the roots of the polynomial P. This is
            %   equivalent to the sintax P.z.
            R = P.z;
        end
        function P = recip(P)
            %RECIP Reciprocal polynomial
            %
            %   Q = RECIP(P) is a polynomial with the same degree of P, of
            %   which the coefficients are presented in revres order, i.e.,
            %   if P = a0 + a1x + a2x^2, then Q = a0x^2 + a1x + a2.
            if P.shift > 0
                error('RECIP is not applicable to quasipolynomials');
            else
                tol = sqrt(eps);
                P.z(abs(P.z) < tol) = [];
                P.k = real(P.k*prod(-P.z));
                P.z = 1./P.z;
                P.coef = striplz(fliplr(P.coef));
            end
        end
        function R = comden(P,Q)
            %COMDEN Least common multiple of two polynomials.
            %
            %   R = COMDEN(P,Q) is a polynomial containing all the common
            %   and non common roots of polynomials P and Q. If P and Q are
            %   different denominators in a sum of rational functions, R is
            %   the common denominator.
            
            %Check variables
            if P.var == Q.var
                allPoles = [P.z; others(Q.z,P.z)];
                R = poln(allPoles,['r' P.var]);
            else
                error('Incompatible variables')
            end
        end
        function C = plus(A,B)
            %+ Sum of polynomials
            %
            %   C = A+B is the sum of polynomials A and B. They must have
            %   the same variable.
            
            if ~isa(B,'poln')
                if isnumeric(B)
                    B = poln(B,A.var);
                else
                    error('Unsuported class for B');
                end
            end
            if A.var == B.var
                %Special cases: very small polynomials
                if norm(A) < eps*norm(B)
                    C = B;
                elseif norm(B) < eps*norm(A)
                    C = A;
                elseif (norm(A) < eps)  &&  (norm(B) < eps)
                    C = poln(0,A.var);
                else
                    tol=eps;
                    %Remove common factor
                    [A,B,X] = coprime(A,B);
                    %Obtain maximum and minimum degrees
                    [pA,mA] = deg(A);
                    [pB,mB] = deg(B);
                    p = max(pA,pB);
                    m = max(mA,mB);
                    %Equalize length of vectors of coefficients
                    cA = [zeros(1,p-pA) A.coef zeros(1,m-mA)];
                    cB = [zeros(1,p-pB) B.coef zeros(1,m-mB)];
                    %Add
                    C = cA + cB;
                    C = findzero(C,tol);
                    C = poln(C, A.var);
                    C.z = [C.z; X.z];
                    C.coef = conv(C.coef, X.coef);
                    C.shift = m;
                    C = zeroing(C);
                end
            else
                error('Incompatible variables')
            end
        end
        function C = minus(A,B)
            %- Difference of polynomials
            %
            %   C = A-B is the difference between polynomials A and B.
            %   They must have the same variable.
            C = A + (-B);
        end
        function [A,B,G] = coprime(A,B)
            %COPRIME Make two polynomials coprime
            %
            %   [A,B] = COPRIME(A,B) removes the common factors from
            %   polynomials A and B, returning the coprime parts.
            %
            %   [A,B,G] = COPRIME(A,B) also returns G, the greatest common
            %   divisor of the two polynomials.
            
            %Set tolerance
            tol = sqrt(eps);
            %Check input
            if A.var ~= B.var
                error('Incompatible variables')
            end
            %initialize result
            G = poln(1,A.var);
            %Check special cases
            if (deg(A) < 1)  || (norm(A) < 1e-10)
                return
            end
            if (deg(B) < 1)  || (norm(B) < 1e-10)
                return
            end
            %General case
            %Extract roots
            rtsA = A.z;
            kA = A.k;
            rtsB = B.z;
            kB = B.k;
            rtsG = [];
            %Find common roots
            i = 1;
            modified = 0;
            while i <= length(rtsA)
                R = rtsA(i);
                tolR = max(tol,tol*abs(R));
                if length(rtsB) < 1
                    break
                end;
                errB = sort(abs(rtsB - R));
                if errB(1) < tolR
                    if abs(imag(R)) > eps % complex conjugate pair
                        R = [R; conj(R)];
                    end
                    rtsG = [rtsG; R];
                    rtsA = remove(rtsA,R,tolR );
                    rtsB = remove(rtsB,R,tolR );
                    modified = 1;
                else
                    i = i + 1;
                end
            end
            %Final results
            if modified
                A = kA*poln(rtsA, ['r' A.var], A.shift);
                B = kB*poln(rtsB, ['r' A.var], B.shift);
                G = poln(rtsG, ['r' A.var]);
            end
        end
        function [U,V] = mldivide(B,A)
            %\ Left division of polynomial objects
            %
            %   [U,V] = B \ A returns the quotient U and the remainder V
            %   such that A = U*B+V. B can be either a scalar or a
            %   polynomial, A must be a polynomial.
            [U,V] = mldivide(A,B);
            %End of MLDIVIDE Method ---------------------------------------
        end
        function [U,V] = mrdivide(A,B)
            %/ Right division of polynomial objects
            %
            %   [U,V] = A / B returns the quotient U and the remainder V
            %   such that A = U*B+V. B can be either a scalar or a
            %   polynomial, A must be a polynomial.
            
            %Check input
            if ~exist('tol','var')
                tol = sqrt(eps);
            end
            %If both are polynomials
            if ~(isnumeric(B) && length(B) == 1)
                if B.var == A.var
                    shiftA = A.shift;
                    shiftB = B.shift;
                    A.shift = 0;
                    B.shift = 0;
                    [A,B] = coprime (A,B);
                    if deg(B) > 0,
                        [V,U] = dioph(poln(1,A.var),B,A);
                    else
                        U = A/B.k;
                        V = poln(0,U.var);
                    end
                    U = mtimesxn (U,shiftB-shiftA );
                    V = mtimesxn (V,-shiftA);
                else
                    error('A and B must have the same variable')
                end
            else%B is a scalar
                U = A; U.k = U.k / B;
                U.coef = A.coef / B;
                V = poln(0, A.var);
            end
            %End of MRDIVIDE Method ---------------------------------------
        end
        function B = mpower(P,n)
            %^ Polynomial power
            %
            %   B = P^n is the result of multiplying polynamial P times
            %   itself n times.
            
            %Check power
            if floor(n) - n ~= 0
                error('Only integer powers are supported')
            end
            B = P;
            %Special case z^n
            if B.shift == 0  &&  deg(B) == 1  &&  B.z == 0
                B.z = [];
                if  n < 0
                    if isct(B)
                        error('Negative powers not supported for continuous-time polynomials');
                    end;
                    B.shift = -n;
                elseif  n > 0
                    B.z = zeros(n,1);
                end;
                B.coef = pol2coef(B);
                return;
            end;
            %General case
            if n < 0
                error('Incorrect negative power for a polynomial');
            elseif n == 0
                B = poln(1,P.var);
            else
                for i=2:n
                    B = B * P;
                end
            end
        end
        function C = mtimes(A,B)
            %* Polinomial multiplication
            %
            %   C = A*B when A and B are polynomials is the standard
            %   polynomial product.
            %
            %   C = A*B when either A or B is a polynomial and the other
            %   variable is a scalar, is the product of the polynomial
            %   times the constant.
            %
            %   C = A*B when B is an LTI object, is another LTI object,
            %   resulting from the series connection of a system with
            %   transfer function A with B.
            
            if isa(A,'poln') && isa(B,'poln')
                if A.var == B.var
                    C = A;
                    C.k = C.k * B.k;
                    C.z = [C.z; B.z];
                    C.coef = conv(A.coef, B.coef);
                    C.shift = A.shift + B.shift;
                else
                    error('Incompatible variables')
                end
            elseif isnumeric(A)
                C = B; C.k = C.k * A;
                C.coef = A*B.coef;
            elseif isnumeric(B)
                C = A; C.k = C.k * B;
                C.coef = B*A.coef;
            elseif isa(B,'lti')
                C = zpk(B);
                C = C * zpk(A,1,C.Ts);
            end
            if isa(C,'poln')
                C = zeroing(C);
            end
        end
        function P = mtimesxn(P,n)
            %MTIMESXN multiplication times x^n
            %
            %   C = mtimesxn(P,n) is the result of multiplying polynomial P
            %   times its variable to the power n.
            
            %Check power
            if floor(n) - n ~= 0
                error('Only integer powers are supported')
            end
            if n == 0
                return
            end
            %Discrete-time polynomials
            if isdt(P)
                P.shift = P.shift - n;
                if n > 0
                    if P.shift < 0
                        P.z = [P.z; zeros(-P.shift,1)];
                        P.shift = 0;
                    end
                    P.coef = pol2coef(P);
                end
            else %Continuous-time polynomials
                indZero = find(abs(P.z)<eps);
                if length(indZero) < P.shift
                    error('Multiplication by x^(-n) cannot be applied to continuous polynomials');
                end
                P.z(indZero(1:P.shift)) = [];
                P.shift = 0;
            end
            P.coef = pol2coef(P);
        end
        function fs = sfactor ( s, type )
            %SFACTOR Spectral factorization of a Hermitian polynomial.
            %        Neutral zeros are divided equally between stable
            %        and unstable parts.
            %
            %    FS = SFACTOR ( S, TYPE )
            %
            %  Inputs:
            %    S - a continuous-time or discrete LTI system
            %    TYPE - factorization type of
            %        's' (default for continuous-time models)
            %        'd' (default for discrete-time models)
            %        'z' for z-plane
            %  Outputs:
            %    FS - stable spectral factor
            %
            
            %------------------------------------------------------
            % Copyright 1999-2006 K. Polyakov
            % $Revision: 3.0$    $Date: 07-Apr-2006 $
            %------------------------------------------------------
            if ~exist('type','var'),
                if isdt(s), type = 'd'; else type = 's'; end;
            end;
            if ~strcmp(type,'s')  &&  ~strcmp(type,'z')  &&  ~strcmp(type,'d')
                error('Unknown factorization type: ''%s''', type);
            end;
            %------------------------------------------------------
            %       Check for zero function
            %------------------------------------------------------
            if s.k == 0,
                fs = poln(0, s.var);
                return;
            end;
            %------------------------------------------------------
            %       Factorize a polynomial
            %------------------------------------------------------
            [zs,zRem] = extrpair ( s.z, type );
            errmsg = 'Exact Hermitian factorization is impossible';
            if ~isempty(zRem)
                error(errmsg); 
            end
            % if z0 ~= 0, error(errmsg); end;
            %------------------------------------------------------
            %       Check gain
            %------------------------------------------------------
            if type == 's'
                K = real(s.k);
                if mod(length(zs),2) == 1, K = - K; end;
            else
                K = real(s.k/prod(-zs));
            end;
            if K < 0, error(errmsg); end;
            %------------------------------------------------------
            %       Form results
            %------------------------------------------------------
            fs = sqrt(K) * poln ( zs, ['r' s.var] );
        end
        %Domain Conversion Methods
        function Q = c2d(P,T,var)
            %C2D Discretization of continuous polynomials
            %
            %   Q = C2D(P) is the discrete version of polynomial P, with
            %   sampling period 1, obtained by replacing the continuous
            %   time variable s with the discrete time version exp(-z).
            %
            %   Q = C2D(P,T) is the discrete version of polynomial P, with
            %   sampling period T, obtained by replacing the continuous
            %   time variable s with the discrete time version exp(-T*z).
            %
            %   Q = C2D(P,T,var) where var is either 'z', 'd' or 'q',
            %   replaces the cntinuous time variable with exp(-var*T)
            %
            %   See also C2Z
            if isct(P)
                if ~exist('T','var')
                    T = 1;
                end
                if ~exist('var','var')
                    var = 'z';
                end
                Q = P;
                Q.z = exp(-P.z*T);
                Q.coef = pol2coef(Q);
                Q.var = var;
            else
                error('The polynomial is already in discrete domain');
            end
            %End of C2D method --------------------------------------------
        end
        function Q = c2z(P,T,var)
            %C2Z Discretization of continuous polynomials
            %
            %   Q = C2Z(P) is the discrete version of polynomial P, with
            %   sampling period 1, obtained by replacing the continuous
            %   time variable s with the discrete time version exp(z).
            %
            %   Q = C2Z(P,T) is the discrete version of polynomial P, with
            %   sampling period T, obtained by replacing the continuous
            %   time variable s with the discrete time version exp(z*T).
            %
            %   Q = C2Z(P,T,var) where var is either 'z', 'd' or 'q',
            %   replaces the cntinuous time variable with exp(var*T)
            %
            %   See also C2D
            if isct(P)
                if ~exist(T,'var')
                    T = 1;
                end
                if ~exist(var,'var')
                    var = 'z';
                end
                Q = P;
                Q.z = exp(P.z*T);
                Q.coef = pol2coef(Q);
                Q.var = var;
            else
                error('The polynomial is already in discrete domain');
            end
            %End of C2Z method --------------------------------------------
        end
        % Type Conversion Methods
        function str = char(P)
            %Transforms a polynomial object into a string to be displayed
            c = P.coef;
            n = length(c);
            degree = n - 1 - P.shift; % highest degree
            str = [];
            for i=1:n
                term=[];
                if c(i) ~= 0
                    if (abs(c(i)) ~= 1) || (degree == 0)
                        term = num2str(abs(c(i)));
                    end
                    if degree == 1
                        term = [term P.var];
                    elseif degree ~= 0
                        term = [term P.var '^' num2str(degree)];
                    end
                    if c(i) < 0
                        term = ['- ' term];
                    elseif i > 1
                        term = ['+ ' term];
                    end
                    str = [str term ' '];
                end
                degree=degree-1;
            end
            if isempty(str)
                str = '0 ';
            end
            %Delete the last space:
            str = str(1:end-1);
            %End of CHAR method--------------------------------------------
        end
        function c = pol2coef(P)
            %POL2COEF transforms a polynomial object into a coefficient vector
            %
            %   C = POL2COEF(P) returns the vector of coefficients of the
            %   polynomial object P.
            
            if abs(P.k) < eps
                c = 0;
            else
                c = P.k * poly(P.z);
            end
            %End of POL2COEF method ---------------------------------------
        end
        function c = double(P)
            %DOUBLE Transforms a polynomial object into a double array
            %
            %   C = DOUBLE(P) is a vector of class double which contains
            %   the coeffientes of the polynomial P.
            %
            %   See also POL2COEF
            
            c = pol2coef(P);
            %End of DOUBLE method -----------------------------------------
        end
        function sys = zpk(N,D,T)
            %ZPK Transforms polynomials into zpk objects
            if ~isa(D,'poln');
                if isa(D,'double')
                    D = poln(D,N.var);
                else
                    error('Unsupported type for D');
                end
            end
            
            if ~isa(N,'poln');
                if isa(N,'double')
                    N = poln(N,D.var);
                else
                    error('Unsupported type for N');
                end
            end
            
            if nargin == 1
                D = poln(1,N.var);
            else
                if N.var~=D.var
                    error('Polynomials must have the same variable');
                end
            end
            if nargin < 3
                if isdt(N)
                    T=1;
                else
                    T=0;
                end
            else
                if isct(N)
                    warning('Continuous time polynomials, ignoring T');
                    T=0;
                end
            end
            [N,D] = coprime( N, D );
            zs = N.z;
            ps = D.z;
            zn = D.shift - N.shift;
            if zn > 0
                zs = [zs; zs(zn,1)];
            elseif zn < 0
                ps = [ps; zs(-zn,1)];
            end;
            sys = zpk(zs,ps,N.k/D.k,T);
        end
        function sys = tf(N,D,T)
            %TF Transforms polynomials into TF objects
            if nargin == 1
                D = poln(1,N.var);
            else
                if N.var~=D.var
                    error('Polynomials must have the same variable');
                end
            end
            if nargin < 3
                if isdt(N)
                    T=1;
                else
                    T=0;
                end
            else
                if isct(N)
                    warning('Continuous time polynomials, ignoring T');
                    T=0;
                end
            end
            [N,D] = coprime(N,D);
            n = N.coef;
            d = D.coef;
            zn = D.shift - N.shift;
            if zn > 0
                n = [n  zeros(1,zn)];
            elseif zn < 0
                d = [d zeros(1,-zn)];
            end
            sys = tf(n,d,T);
        end
        function sys = sdzpk(N,D,T)
            %SDZPK Transforms polynomials into sdzpk objects
            sys = zpk(N,D,T);
            sys = sdzpk(sys);
        end
        function sys = sdtf(N,D,T)
            %SDTF Transforms polynomials into TF objects
            sys = tf(N,D,T);
            sys = sdtf(sys);
        end
        %Miscleaneous Methods
        function flag = isdt(P)
            %ISDT True for discrete time polynomials
            %
            %ISDT(P) returns true if P is a discrete time polynimial, i.e.,
            %if P.var is either 'd', 'z', or 'q'.
            
            if length(P) == 1
                ind = find(P.var == 'zdq',1);
                flag = ~isempty(ind);
            else
                flag = isdt(P(1));
            end
            
            %End of ISDT method -------------------------------------------
        end
        function flag = isct(P)
            %ISCT True for continuous time polynomials
            %
            %ISDT(P) returns true if P is a continuous time polynimial,
            %i.e., if P.var is either 's' or 'p'.
            if length(P) == 1
                ind = find(P.var == 'sp',1);
                flag = ~isempty(ind);
            else
                flag = isct(P(1));
            end
            %End of ISCT method -------------------------------------------
        end
        function A = delzero(A,tol)
            %DELZERO Delete trailing zeros
            %
            %   A = DELZERO(A) removes all the unwanted zero coefficients
            %   from polynomial A.
            %
            %   A = DELZERO(A,TOL) uses tolerance TOL to decide if a number
            %   is zero.
            %
            %   [A,nz] = DELZERO(A) also indicates the number of zeros that
            %   have neen removed.
            %
            
            if ~exist('tol','var')
                tol = eps;
            end
            ind = find(abs(A.z) < tol);
            nz = length(ind);
            if nz ~= 0
                A.z(ind) = [];
                A.coef = A.coef(1:end-nz);
            end
        end
        function A = zeroing(A)
            %ZEROING remove zero coefficieents
            
            if abs(A.k) < eps,
                if abs(A.k*sum(abs(A.z))) < eps,
                    A.k = 0;
                    A.z = [];
                    A.coef = 0;
                    A.shift = 0;
                else
                    [~,i] = max(abs(A.z));
                    A.k = - A.k*A.z(i);
                    A.z(i) = [];
                    A.coef(1) = [];
                    A = zeroing(A);
                end
            else
                %       Check for zero negative degrees
                indZero = find(abs(A.z) < eps);
                nZero = min(A.shift, length(indZero));
                if nZero > 0
                    A.shift = A.shift - nZero;
                    A.z(indZero(1:nZero)) = [];
                    A.coef(length(A.coef)-nZero+1:end) = [];
                end
                %       Check for zero positive degrees
                nZero = length(A.coef) - length(A.z) - 1;
                if nZero > 0
                    A.coef(1:nZero) = [];
                end
            end
        end
        %Likely unnecessary, candidates to deletion
        function P = change(P,k,z,coef)
            %CHANGE Change all data in a polynomial
            %
            %   Q = CHANGE(P,K,Z,COEF) Changes the properties of polynomial
            %   P to the gain K, the roots given in vector Z and the
            %   coefficients in vector COEF.
            
            tol = sqrt(eps);
            P.k = k;
            P.z = z;
            P.coef = pol2coef ( P );
            if exist('coef','var')
                if length(P.coef) ~= length(coef)  ||  ...
                        norm(P.coef-coef) > tol*norm(P.coef)
                    error('Incorrect coefficient array');
                end
                P.coef = coef;
            end
            %End of CHANGE method ---------------------------------------------
        end
        function P = fresh(P)
            %fresh coefficient array after changing roots
            P.coef = pol2coef(P);
        end
        function [q,r] = quorem (a,b)
            %QUOREM Polynomial quotient and remainder.
            [q,r] = a/b; %#ok<RHSFN>
        end
    end
    %Protected methods: accessible only to objects of this class or
    %subclasses
    methods (Access = protected)
        %Former private methods
        %This one is also likely unnecessary.
        function [zs,z,z0] = extrpair ( z, type, tol )
            %EXTRPAIR Extract a pair of symmetric zeros.
            %
            %     [ZS,R,Z0] = EXTRPAIR ( Z, TYPE, TOL )
            %
            %   Inputs:
            %     Z    - initial array
            %     TYPE - factorization type of 's' (default), 'z' and 'd'
            %     TOL  - relative tolerance (default 1e-4)
            %
            %   Outputs:
            %     ZS - factorized stable zeros
            %     R  - unpaired nonzero zeros
            %     Z0 - unpaired zeros at origin
            %
            
            %Check data
            if ~exist('type','var')
                type = 's';
            end
            if isempty(find(type == 'szd', 1))
                error('Unknown factorization type ''%s''', type);
            end
            if ~exist('tol','var')
                tol = 1e-4;
            end
            %------------------------------------------------------
            %       Start from the terms with maximal absolute values
            %------------------------------------------------------
            zs = []; z = z(:);
            [~,ind] = sort(abs(z));
            z = z(flipud(ind));
            %------------------------------------------------------
            %       Iteration
            %------------------------------------------------------
            while 1
                len = length(z);
                if len < 2
                    break
                end
                minDist = Inf;
                iMin = 0;
                jMin = 0;
                for i=1:len
                    for j=1:len
                        if i ~= j
                            switch type
                                case 's'
                                    dij = abs(z(i) + z(j));
                                case {'z', 'd'}
                                    dij = abs(z(i)*z(j)-1);
                            end
                            if dij < minDist && ...
                                    (abs(z(i)-1) < tol  ||  ...
                                    ~(abs(z(i)) == 1  &&  abs(imag(z(i))) > eps))
                                iMin = i;
                                jMin = j;
                                minDist = dij;
                            end
                        end
                    end
                end
                if minDist > tol  &&  minDist > tol*abs(z(iMin))
                    break
                end
                %------------------------------------------------------
                %       Extract pair
                %------------------------------------------------------
                if isequal(type,'s')
                    if z(iMin) < 0
                        zs1 = z(iMin);
                    else
                        zs1 = z(jMin);
                    end
                elseif isequal(type,'z')
                    if abs(z(iMin)) < 1
                        zs1 = z(iMin);
                    else
                        zs1 = z(jMin);
                    end
                else
                    if abs(z(iMin)) > 1
                        zs1 = z(iMin);
                    else
                        zs1 = z(jMin);
                    end
                end
                %z([iMin jMin])
                zs = [zs; zs1];
                z([iMin jMin]) = [];
            end
            %------------------------------------------------------
            %       Find unpaired complex zeros
            %------------------------------------------------------
            indc = find(imag(zs) ~= 0).';
            if ~isempty(indc)
                for i = indc
                    [zz,ii] = sort(abs(zs - conj(zs(i))));
                    if ii(1) == i
                        if imag(zs(i)) > 1e-6,
                            warning('Unpaired complex conjugate zeros');
                        end
                        zs(i) = real(zs(i));
                    end
                end
            end
            %------------------------------------------------------
            %       Find zero zeros
            %------------------------------------------------------
            z0 = 0;
            ind0 = find(abs(z) < tol);
            if ~isempty(ind0)
                z0 = length(ind0);
                z(ind0) = [];
            end
        end
    end
end