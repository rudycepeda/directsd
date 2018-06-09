function [x,y] = dioph(A,B,C)
%DIOPH Diophantine equation solver
%   [X,Y] = DIOPH(A,B,C) gives the solution to the Diophantine equation
%   A*X+B*Y=C. The equation can be improper, i.e., deg(C) >= deg(A)+deg(B)

%Original version by Dr. K. Yu. Polyakov.
%Revised and updated by R. Cepeda-Gomez

%Version 4.0 17.07.2014%

%Check if the inputs are polynomials
if isa(A,'poln') && isa(B,'poln') && isa(C,'poln')
    if (A.var == B.var) && (B.var == C.var)
        polnMode = 1;
        %[A,B,C] = triple(A,B,C);
        a = A.coef;
        b = B.coef;
        c = C.coef;
    else
        error('Incompatible variables')
    end
elseif ~(isnumeric(A)&&isnumeric(B)&&isnumeric(C))
    error('DIOPH defined only for input arguments of class poln or numeric')
else
    a=A;
    b=B;
    c=C;
    polnMode = 0;
end
%Making sure the inputs are row vectors
a = a(:)';
b = b(:)';
c = c(:)';
%Changing the order of the coefficients to form the linear system
a = fliplr(a);
b = fliplr(b);
c = fliplr(c);
%Determining the sizes of the array
na = length(a)-1;
nb = length(b)-1;
nc = length(c)-1;
m = max(na,nc-nb+1);
%Forming the array
C = zeros(nb+m,1);
C(1:nc+1) = c';
M = zeros(nb+m,nb+m);
for i = 1:nb
    M(i:i+na,i)= a';
end
for i = 1:m
    M(i:i+nb,i+nb) = b';
end
%Solving the linear system
sol = linsolve(M,C);
%Extracting the results
x = flipud(sol(1:nb));
y = flipud(sol(nb+1:nb+m));
%Formating output
x=striplz(x(:)');
y=striplz(y(:)');
if polnMode
    x = poln(x,A.var);
    y = poln(y,A.var);
end
end