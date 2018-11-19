function val = error4eRT0Energy(c4n, n4e, gradExact, pApprox)
%% error4eRT0energy - energy error for RT0 element
% Calculates the energy error of the RT0 finite element solution with given
% exact gradiant function gradExact.
%
% Input:     c4n       coordinates for the nodes of the mesh
%            n4e       nodes for the elements of the mesh
%            p         first component of the solution
%            u         second component of the solution
%            gradExact exact gradiant of the solution
%
% Output:    val       L2 error of the gradiant on each element.

val = sum(...
    integrate(c4n, n4e, @(n4p,Gpts4p,Gpts4ref)...
    (pApprox(:,1)*[1 0] + pApprox(:,2)*[0 1]...
    + [pApprox(:,3) pApprox(:,3)] .* (Gpts4p - computeMid4e(c4n,n4p)) ...
    - gradExact(Gpts4p)).^2, 6)...
    ,2);
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
