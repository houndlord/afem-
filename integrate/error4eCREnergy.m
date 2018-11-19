function val = error4eCREnergy(c4n,n4e,gradExact,uApprox)
    nrElems = size(n4e,1);
    s4e = computeS4e(n4e);
    gradUApprox = zeros(nrElems,2);
    
    % Compute gradient for x
    for elem = 1:size(n4e,1);
        grads = [c4n(n4e(elem,:),:)'; 1 1 1] \ [-2 0; 0 -2; 0 0];
        gradUApprox(elem,:) = uApprox(s4e(elem,:))' * grads([3 1 2],:);
    end
    
    val = sum(integrate(c4n,n4e, ...
        @(n4p, Gpts4p, Gpts4ref) (gradExact(Gpts4p) - gradUApprox).^2, 6),2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2009-2015
% Numerical Analysis Group
% Prof. Dr. Carsten Carstensen
% Humboldt-University
% Departement of Mathematics
% 10099 Berlin
% Germany
%
% This file is part of AFEM.
%
% AFEM is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% AFEM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
