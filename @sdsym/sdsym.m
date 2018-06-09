classdef sdsym < sym
    methods
        function obj = sdsym(varargin)
            obj = obj@sym(varargin{:});
        end
        disp(X);
        display(X);
        Mi = invmat(M,recFlag);
        y = transpsym(x);
    end
end