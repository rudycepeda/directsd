function options = updateopt ( options, opt )
%UPDATEOPT Update options structure.
%
%     UOPTIONS = UPDATEOPT ( OPTIONS, NEWOPT )
%
%   Inputs:
%     OPTIONS - existing options structure
%     NEWOPT  - update
%
%   Outputs:
%     UOPTIONS - updated options structure
%
%   See also FIELDS.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Special case
%------------------------------------------------------
        if ~exist('opt','var') || ~isstruct(opt)
          return;
        end;
%------------------------------------------------------
%       Do update
%------------------------------------------------------
        fOpt = fields ( opt );
        for i=1:length(fOpt)
          if isfield(options, fOpt{i})  
            options = setfield (options, fOpt{i}, ...
                                getfield(opt,fOpt{i}));  
          end;
        end;

%------- End of UPDATEOPT.M --------- KYuP ----------           
