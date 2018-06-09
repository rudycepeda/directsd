function setdate ( dir, rev )
%SETDATE Set current date to M-files.
%
%   SETDATE ( DIR )
%   SETDATE ( DIR, VER )
%
%   Inputs:
%     DIR - directory
%     VER - toolbox version
%
%   Outputs:
%     None
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------        
        if ~exist('dir','var') || isempty(dir)
           dir = '.';
        end;
        if ~exist('rev','var'), 
           rev = []; 
        elseif ~ischar(rev), 
           rev = sprintf('%.2f',rev);
        end;
%------------------------------------------------------
%       Get director list
%------------------------------------------------------
        fList = ls(dir);
        for i=1:size(fList,1)
          fname = deblank(fList(i,:));
          if fname(1) == '.', continue; end;
          fname = [dir '\' fname];
          [path,f,ext] = fileparts(fname);
          if isempty(ext)
             if ~isempty(path), f = [path '\' f]; end; 
             disp ( f );
             setdate ( f ); 
             continue; 
          end;
          if ~isequal(ext,'.m'), continue; end;              
          xname = [f '.$$$'];
          if ~isempty(path), 
             xname = [path '\' xname]; 
          end;
         %---------------------------------------------
         %       Handle file
         %---------------------------------------------
          in = fopen(fname);
          out = fopen(xname, 'wt');
          changed = 0;
          while 1
            s = fgetl(in);
            if ~ischar(s), break, end;
            if ~isempty(strfind(s, 'strfind'))
              fprintf(out, '%s\n', s);
              continue;  
            end;
            n = strfind(s, '$Date: ');
            if ~isempty(n)
              s = [s(1:n+6) datestr(today,1) ' $']; 
              changed = 1;
            end;
            if ~isempty(rev)
              n = strfind(s, '$Revision: ');
              if ~isempty(n)
                s1 = s(1:n+10);
                s = s(n+11:end);
                n1 = strfind(s, '$');
                if ~isempty(n1)
                  s = [s1 rev ' ' s(n1:end)];
                  changed = 1;
                end;
              end;
            end;
            fprintf(out, '%s\n', s);
          end
          fclose(in);
          fclose(out);
          if changed
            delete(fname);
            try 
              movefile(xname, fname);
              disp([fname ': OK']);  
            catch
              disp([fname ': ERROR']);  
            end;
          else
            delete(xname);
            disp(['    ' fname]);  
          end;
        end;

%------- End of SETDATE.M --------- KYuP ----------
