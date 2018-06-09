function flag = isint(n)
%ISINT True for integer numbers
%
%   FLAG = ISINT(N) returns a logical true if N is an integer number
%   (regardles of the actual class of the variable).
%
%   See also ISINTEGER

%Written by R. Cepeda-Gomez.
%29.07.2014

if floor(n)-n == 0
    flag = 1;
else
    flag = 0;
end
flag = boolean(flag);