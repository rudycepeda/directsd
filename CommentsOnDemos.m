%Comments on the demo files:

%   demo_fil1
%       Runs. The resulting controllers are the same as listed in the help
%       files. The plot with the variance between sampling instants is not.
%       This does not make sense, the variance should be the same at
%       sampling instants.
%
%       Fixed. There was a problem with the computation of the H2 norm.
%       There are two methods to perform this compuation: State-space and
%       polynomial approaches. The SS method is the default, but it is
%       giving incorrect results. The Polynomial method takes a little bit
%       longer to run, but it works fine. Forcing all the computations to
%       use the polynomial method solved the issue.

%   demo_fil2
%       Runs as expected.

%   demo_fil3
%       Runs. The controller K is slightly different from the one found in
%       the help files. A problem similar to demo_fil1 appears: the plot of
%       variance between sampling instants is different and not symmetric.
%
%       Fixed in the same way as in demo_fil1

%   demo_hold
%       Runs. Again, there is a mismatch in the plot.
%
%       Fixed in the same way.

%   demo_at96
%       The gain of the controller K reported is different than the one in
%       the help file. The optimal cost, however, is the same. 

%   demo_ait98      
%       Runs. Reports static controllers. There is no HTML help for this
%       demo. Not a clear match to the examples in the referenced papers.

%   demo_autom97    
%       Runs. The results match those in the paper. There is no HTML help 
%       for this demo. The plot of the variance between sampling instants
%       is different.

%   demo_doubint    
%       Runs as expected.

%   demo_l2
%       Runs as expected.

%   demo_2dof
%       Runs. There are warnings. The values of the errors are not the same 
%       as in the help. The controllers and plots are correct.

%   demo_c2d
%       Runs as expected.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  The Ahinf part is the one that seems to not be working in the
%%%%%  following two demos
%   demo_ait01b
%       There is no help file. The results are different from those in the
%       refernced paper. However, it is not clear to which examples in the
%       papers the demo corresponds.

%   demo_l2hinf
%       There is no help file. The results are different from those in the
%       refernced paper. However, it is not clear to which examples in the
%       papers the demo corresponds.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  The Ahinf part is the one that seems to not be working.

%   demo_cf1
%       Runs. No help file. I have not checked the book.

%   demo_cf2
%       Runs. No help file. I have not checked the book.


%   demo_dhinf
%       For the first example, the static controller is reported as 3 and
%       not as 1.5 as in the help file. For the other two examples, the
%       results are completely different from the expected results
%       accroding to the help file.

%   demo_h2hinf
%       Runs. No help or corresponding paper.

%   demo_modsdh2
%       Runs. No help file or corresponding paper to compare.

%   demo_modsdh2int
%       Runs. No help file or corresponding paper to compare.

%   demo_modsdl2
%       Runs. No help file or corresponding paper to compare.

%   demo_modsdl2a
%       Does not run

%   demo_modsdl2b
%       Does not run

%   demo_h2p
%       Runs. There are several warnings about hermitian factorization not
%       being possible. (The error messages were changed to warnings).
%       Results are not as expected.

%   demo_l2p
%       Runs. There are several warnings about hermitian factorization not
%       being possible. (The error messages were changed to warnings).
%       Results are not as expected.

%   demo_2dofp
%       Runs. There are several warnings about hermitian factorization not
%       being possible. (The error messages were changed to warnings).
%       Results are not as expected.
