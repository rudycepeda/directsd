function hlp2demo ( dir )
%HLP2DEMO Transform help files to demo HTML-files.
%
%   HLP2DEMO ( DIR )
%
%   Inputs:
%     DIR - directory
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
%------------------------------------------------------
%       Get director list
%------------------------------------------------------
        fList = ls(dir);
        for i=1:size(fList,1)
          fname = deblank(fList(i,:));
          if fname(1) == '.', continue; end;
          fname = [dir '\' fname];
          [path,f,ext] = fileparts(fname);
          if ~isequal(ext,'.html')  &&  ~isequal(ext,'.htm'), 
            continue; 
          end;              
          xname = [f '.$$$'];
          if ~isempty(path), 
             xname = [path '\' xname]; 
          end;
         %---------------------------------------------
         %       Read file
         %---------------------------------------------
          in = fopen(fname);
          out = fopen(xname, 'wt');
         %---------------------------------------------
         %       Read file
         %---------------------------------------------
          sAll = [];
          while 1
            s = fgetl(in);
            if ~ischar(s), break, end;
            if isempty(sAll), 
                 sAll = s;
            else sAll = [sAll char(13) s];
            end;
          end
         %---------------------------------------------
         %       Remove anchors
         %---------------------------------------------
          changed = 0;
          while ( 1 )
            n = strfind(upper(sAll), '<A HREF="..');
            if isempty(n), 
               n = strfind(upper(sAll), '<A HREF="#');
            end;
            if isempty(n), break; end;
            tail = sAll(n+10:end);
            n1 = strfind(tail, '>');
            if isempty(n1), break; end;
            tail = tail(n1+1:end);
            n2 = strfind(upper(tail), '</A>');
            if isempty(n2), break; end;
            s = tail(1:n2-1);
            if isempty(strfind(s,' ')) 
              tail = [tail(1:n2) '/code' tail(n2+3:end)];
              sAll = [sAll(1:n-1) '<code>' tail];;             
            else
              tail(n2:n2+3) = [];
              sAll = [sAll(1:n-1) tail];
            end;
            changed = 1;
          end;          
         %---------------------------------------------
         %       Remove ../
         %---------------------------------------------
          while ( 1 )
            n = strfind(sAll, '../');
            if isempty(n), break; end;
            sAll(n:n+2) = [];             
            changed = 1;
          end;          
         %---------------------------------------------
         %       Remove top references
         %---------------------------------------------
          n = strfind(sAll,'<table border=0 width="100%"');
          if ~isempty(n)
            tail = sAll(n:end);
            n1 = strfind(tail, '</table>');
            if ~isempty(n1)
              tail(1:n1+7) = [];  
              sAll = [sAll(1:n-1) tail];
              changed = 1;
            end;
          end;
         %---------------------------------------------
         %       Remove bottom references
         %---------------------------------------------
          n = strfind(sAll,'<p><table bgcolor="#e7ebf7"');
          if ~isempty(n)
            tail = sAll(n:end);
            n1 = strfind(tail, '</table>');
            if ~isempty(n1)
              tail(1:n1+7) = [];  
              sAll = [sAll(1:n-1) tail];
              changed = 1;
            end;
          end;
         %---------------------------------------------
         %       Save file and finalize
         %---------------------------------------------
          fprintf(out, '%s', sAll);
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

%------- End of HLP2DEMO.M --------- KYuP ----------
