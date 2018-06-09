function fixepsdash ( fname ) 
%FIXEPSDASH Make proper dash length in EPS-files.
%
%     FIXEPSDASH ( FNAME )
%
%   Inputs:
%     FNAME - file name to fix
%

%------------------------------------------------------
% Author: Paul Skoczylas
% Modified by K. Polyakov for Matlab 7.02, 2006 
% $Revision: 1.00 $    $Date: 01-Dec-2006 $
%------------------------------------------------------
%       Open file in read-write mode
%------------------------------------------------------
    fid = fopen(fname, 'r+'); 
%------------------------------------------------------
%       Find overall length    
%------------------------------------------------------
    fseek(fid, 0, 'eof'); 
    iEnd = ftell(fid); 
%------------------------------------------------------
%       Start from the beginning
%------------------------------------------------------
    fseek(fid, 0, 'bof'); 
    i = ftell(fid); 
    L = fgetl(fid); 
    done = 0; 
%------------------------------------------------------
%       Main loop
%------------------------------------------------------
    while (i ~= iEnd & done<3) 
       %fprintf('.%s.\n',L)
       %length(L)
       if length(L) == 40  &  isequal(L(1:3),'/DA') % L=='/DA { [6 dpi2point mul] 0 setdash } bdef' 
          fseek(fid, i, 'bof'); 
          fprintf(fid,'/DA { [2 dpi2point mul] 0 setdash } bdef'); 
          fgetl(fid);
          done=done+1; 
       end; 
       if length(L) == 57  &  isequal(L(1:3),'/DO') % L=='/DO { [.5 dpi2point mul 4 dpi2point mul] 0 setdash } bdef' 
          fseek(fid, i, 'bof'); 
          fprintf(fid,'/DO { [.5 dpi2point mul 2 dpi2point mul] 0 setdash } bdef'); 
          fgetl(fid);
          done=done+1; 
       end; 
       if length(L) == 57  &  isequal(L(1:3),'/DD') % L=='/DD { [.5 dpi2point mul 4 dpi2point mul 6 dpi2point mul 4 dpi2point mul] 0 setdash } bdef' 
          fseek(fid, i, 'bof');                     
          fprintf(fid,'/DD { [.5 dpi2point mul 2 dpi2point mul 3 dpi2point mul 2');                    
          fgetl(fid);
          done=done+1; 
       end; 
       i = ftell(fid);   
       L = fgetl(fid); 
    end; 
%------------------------------------------------------
%       Close file
%------------------------------------------------------
    fid = fclose(fid); 

%------- End of FIXEPSDASH.M --------- KYuP ----------
