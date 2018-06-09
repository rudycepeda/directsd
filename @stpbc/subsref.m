function R = subsref ( P, S )
%SUBSREF Property of a STPBC.
%       
%       S.a, S.b, S.c, S.d - state-space matrices
%       S.om, S.ups        - boundary condition matrices
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
        for i=1:length(S)
           switch S(i).type
             %-----------------------------------------
             %  Extract properties
             %-----------------------------------------
              case '.',
                if strcmpi(S(i).subs,'a'),
                   R = P.a;
                elseif strcmpi(S(i).subs,'b'),
                   R = P.b;
                elseif strcmpi(S(i).subs,'c'),
                   R = P.c;
                elseif strcmpi(S(i).subs,'d'),
                   R = P.d;
                elseif strcmpi(S(i).subs,'om'),
                   R = P.om;
                elseif strcmp(S(i).subs,'ups'),
                   R = P.ups;
                else 
                   error('Unknown property');
                end;         
           end; % switch
        end; % for        

%------- End of SUBSREF.M --------- KYuP ----------           

