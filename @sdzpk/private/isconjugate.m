function boo = isconjugate(r)
%ISCONJUGATE   Checks if roots from a complex conjugate set.

%      Author: P. Gahinet
%      Copyright 1986-2002 The MathWorks, Inc. 
%      $Revision: 1.3 $  $Date: 07-Apr-2006 $

boo = isequal(sort(r(imag(r)>0)),sort(conj(r(imag(r)<0))));
