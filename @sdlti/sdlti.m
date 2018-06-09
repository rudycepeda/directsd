classdef sdlti < lti
    methods (Access = public)
        function obj = lti(varargin)
            obj = lti@lti(varargin{:});
        end
        L = ctranspose(L);
        sys = pvset(sys,varargin);
        T = trace(sys);
        v = val(sys,x);                    
    end   
end