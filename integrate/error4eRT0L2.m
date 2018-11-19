function val = error4eRT0L2(c4n, n4e, uExact, uApprox)
%% error4eRT0energy - energy error for RT0 element
% Calculates the L2 error of the RT0 finite element solution with given
% exact solution uExact.
%
% Input:     c4n       coordinates for the nodes of the mesh
%            n4e       nodes for the elements of the mesh
%            p         first component of the solution
%            u         second component of the solution
%            uExact    exact gradiant of the solution
%
% Output:    val       L2 error of the solution u on each element.

val = integrate(c4n, n4e, @(n4p,Gpts4p,Gpts4ref)...
    (uApprox - uExact(Gpts4p)).^2,2);
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
