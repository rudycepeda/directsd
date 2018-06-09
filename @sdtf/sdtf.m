classdef sdtf < tf & sdlti
    methods (Access = public)
        function obj = sdtf(varargin)
            if isa(varargin{1},'tf') 
                old = varargin{1};
                args = {old.num old.den old.Ts ...
                    'Variable'         old.Variable...
                    'ioDelay'          old.ioDelay...
                    'InputDelay'       old.InputDelay...
                    'OutputDelay'      old.OutputDelay...
                    'TimeUnit'         old.TimeUnit...
                    'InputName'        old.InputName...
                    'InputUnit'        old.InputUnit...
                    'InputGroup'       old.InputGroup...
                    'OutputName'       old.OutputName...
                    'OutputUnit'       old.OutputUnit...
                    'OutputGroup'      old.OutputGroup...
                    'Name'             old.Name...
                    'Notes'            old.Notes...
                    'UserData'         old.UserData};
            else
                args = varargin;
            end
            obj = obj@tf(args{:});
        end
        F = z2zeta (F,tol)
        sys = delay2z(sys)            
    end   
end