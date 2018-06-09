function FP = setpoles ( F, p )
%SETPOLES Set specified poles to an LTI model.
%
%     FP = SETPOLES ( F, P )
%
%   Inputs:
%     F - initial LTI model
%     P - vector of poles it must have
%
%   Outputs:
%     FP - final model which has the desired poles
%          If F has not some poles which are in P, they are added.
%          If F has extra poles (not present in P), an error message
%          is generated.
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
if ~isobject(F) 
    F = sdzpk(F);  
end;
FP = minreal(sdzpk(F));
[rows,cols] = size(FP);
tol = 1e-4;
p0 = p;
%------------------------------------------------------
%       Loop over all items
%------------------------------------------------------
for r=1:rows
    for c=1:cols
        p = p0;
        pF = FP.p{r,c};
        zF = FP.z{r,c};
        %-----------------------------------------------
        %       Correct existing poles
        %-----------------------------------------------
        pFCorr = [];
        for i=length(pF):-1:1
            if isempty(p), break; end;
            [xx,ind] = sort(abs(p-pF(i)));
            if xx(1) < tol
                pFCorr = [pFCorr; p(ind(1))];
                pF(i) = [];
                p(ind(1)) = [];
            end;
        end;
        %-----------------------------------------------
        %       Try to cancel extra poles
        %-----------------------------------------------
        lack = others ( pF, p, tol );
        errmsg = 'Initial model contains an extra pole %g';
        if isdt(F) % do not cancel extra poles at z=0
            ind = find(abs(lack)<eps);
            pFCorr = [pFCorr; lack(ind)];
            lack(ind) = [];
        end;
        for i=1:length(lack)
            if isempty(zF), error(errmsg,lack(i)); end;
            [~,ind] = sort(abs(zF-lack(i)));
            zF(ind(1)) = [];
        end;
        %-----------------------------------------------
        %        Insert additional desired poles
        %-----------------------------------------------
        new = others ( p, pF, tol );
        if size(new,2) > 1, new = new'; end;
        zF = [zF; new];
        pF = [pFCorr; new];
        FP.z{r,c} = zF;
        FP.p{r,c} = pF;
    end;
end;

%------- End of SETPOLES.M --------- KYuP ----------
