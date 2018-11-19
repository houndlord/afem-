function jump4s = P0TangentJump(c4n,n4e,n4sDb,n4sNb,sigma4e,u4Db)
%% P0TangentJump - computes jumps in tangential direction
%   c4n,n4e,n4sDb,n4sNb: mesh data
%   sigma4e: the values of a P0 function given for each element
%   u4Db: a function giving values for the Dirichlet boundary
%
%   jump4s: the jumps of v across each side, for Dirichlet boundary
%           sides, u4Db is the expected value.

    %% Initialisation
    s4e = computeS4e(n4e);
    s4n = computeS4n(n4e);
    n4s = computeN4s(n4e);
    tangent4e = computeTangent4e(c4n,n4e);
    length4s = computeLength4s(c4n,n4s);
    jump4s = zeros(size(n4s,1),1);
    
    %% inner jumps
    for elem = 1 : size(n4e,1)
        jump4s(s4e(elem,:)) = jump4s(s4e(elem,:)) ...
                              + tangent4e(:,:,elem)*sigma4e(elem,:)';
    end
    
    %% Neumann jumps
    jump4s(diag(s4n(n4sNb(:,1),n4sNb(:,2)))) = 0;

    %% Dirichlet jumps
    for nodes = n4sDb'
        side = s4n(nodes(1),nodes(2));
        jump4s(side) = jump4s(side) ...
            - (u4Db(c4n(nodes(2),:))-u4Db(c4n(nodes(1),:)))/length4s(side);
    end
    
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
