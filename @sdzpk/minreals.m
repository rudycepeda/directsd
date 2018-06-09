function F = minreals ( F, tol )
%MINREALS Minimal realization of symmetric LTI-model.
%
%     F0 = MINREALS ( F, TOL )
%
%   Inputs:
%     F  - symmetric LTI-model (spectral density)
%     TOL - tolerance
%
%   Outputs:
%     F0 - the function with cancelled common pole-zero pairs
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov
% $Revision: 3.0$    $Date: 10-Jan-2009 $
%------------------------------------------------------
%         Preliminaries
%------------------------------------------------------
if ~exist('tol','var'), tol = sqrt(eps); end;
if isct(F), type = 's'; else type = 'd'; end;
F = zpk(F);
zz = F.z{1};
pp = F.p{1};
%------------------------------------------------------
%         Extract pairs of roots
%------------------------------------------------------
[zs,zRem,z0] = extrpair ( zz, type, tol );
[ps,pRem,p0] = extrpair ( pp, type, tol );
if ~isempty(zRem) ||  ~isempty(pRem)
    error('The function is not symmetric');
end;
modified = 0;
%------------------------------------------------------
%         Cancel common roots at the origin
%------------------------------------------------------
c0 = min(z0,p0);
if c0 > 0, modified = 1; end;
z0 = z0 - c0; p0 = p0 - c0;
if isct(F)
    if mod(z0,2) ~= 0  ||  mod(p0,2) ~= 0
        error('Incorrect zeros or poles at the origin');
    end;
else
    if length(zs)+z0 ~= length(ps)+p0
        error('Incorrect zeros or poles at the origin');
    end;
end;
%------------------------------------------------------
%         Find common roots
%------------------------------------------------------
i = 1;
ind = find(abs(imag(zs)) < tol); zs(ind)=real(zs(ind));
ind = find(abs(imag(ps)) < tol); ps(ind)=real(ps(ind));
while i <= length(ps)
    A = ps(i);
    tolA = max(tol*abs(A), tol);
    if length(zs) < 1, break; end;
    err = sort(abs(zs - A));
    if err(1) < tolA
        modified = 1;
        if abs(imag(A)) > eps % pair of complex eigenvalues
            A = [A; conj(A)];
        end;
        zs = remove ( zs, A, 100*tolA );
        ps = remove ( ps, A, 100*tolA );
    else
        i = i + 1;
    end;
end
%------------------------------------------------------
%       Final results
%------------------------------------------------------
if modified
    if isct(F)
        zz = [zs; -zs; zeros(z0/2,1)];
        pp = [ps; -ps; zeros(p0/2,1)];
    else
        zz = [zs; 1./zs; zeros(z0,1)];
        pp = [ps; 1./ps; zeros(p0,1)];
    end;
    if ~isempty(zz)
        F.z{1} = zz;
    else
        F.z{1} = [];
    end
    if ~isempty(pp)
        F.p{1} = pp;
    else
        F.p{1} = [];
    end
end

%------- End of MINREALS.M --------- KYuP ----------
