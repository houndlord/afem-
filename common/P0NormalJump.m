function jump4s = P0NormalJump(c4n,n4e,n4sDb,n4sNb,sigma4e,g)
%% P0NormalJump - computes jumps in normal direction
%   c4n,n4e,n4sDb,n4sNb: mesh data
%   sigma4e: the values of a P0 function given for each element
%   g: a function giving values for the Neumann boundary
%
%   jump4s: the jumps of v across each side, for Neumann boundary
%           sides, g is the expected value.

    %% Initialisation
    s4e = computeS4e(n4e);
    s4n = computeS4n(n4e);
    normal4e = computeNormal4e(c4n,n4e);
    n4s = computeN4s(n4e);
    mid4s = computeMid4s(c4n, n4s);
    jump4s = zeros(size(n4s,1),1);
    
    %% inner jumps
    for elem = 1 : size(n4e,1)
        jump4s(s4e(elem,:)) = jump4s(s4e(elem,:)) ...
                              + normal4e(:,:,elem)*sigma4e(elem,:)';
    end

    %% Neumann jumps
    for nodes = n4sNb'
        side = s4n(nodes(1),nodes(2));
        jump4s(side) = jump4s(side) - g(mid4s(side,:));
    end

    %% Dirichlet jumps = 0
    jump4s(diag(s4n(n4sDb(:,1),n4sDb(:,2)))) = 0;
    
    jump4s = abs(jump4s);

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
