function [eta4s,n4s] = estimateCREtaSides(f,g,u4Db,x,c4n,n4e,n4sDb,n4sNb)
%% estimateCREtaSides - error estimator for CR element
% Estimate the energy error of the CR finite element solution by the
% jumps of the discrete solution's tangents along the sides.
%
% Input:     f	       right-hand side of the problem definition
%            g         Neumann boundary condition
%            u4Db      Dirichlet boundary condition
%            x         CR basis coefficients of u given by solve
%            c4n       coordinates for the nodes of the mesh
%            n4e       nodes for the elements of the mesh
%            n4sDb	   the nodes of the sides in the Dirichlet boundary
%            n4sNb	   the nodes of the sides in the Neumann boundary
%
% Output:    eta4s     error for each side
%            n4s       nodes of the sides for which the error was computed

    %% Initialisation
    s4e      = computeS4e(n4e);
    n4s      = computeN4s(n4e);
    length4s = computeLength4s(c4n,n4s);
    gradU    = zeros(size(n4e,1),2);
    
    %% Compute gradient.
    for elem = 1:size(n4e,1);
        grads = [c4n(n4e(elem,:),:)'; 1 1 1] \ [-2 0; 0 -2; 0 0];
        gradU(elem,:) = x(s4e(elem,:))' * grads([3 1 2],:);
    end
    
    %% Compute the L2-norm of the jumps and weigh them with length4s.
    eta4sNormal = P0NormalJump(c4n,n4e,n4sDb,n4sNb,gradU,g);
    eta4sTangent = P0TangentJump(c4n,n4e,n4sDb,n4sNb,gradU,u4Db);
    eta4s = ((eta4sNormal+eta4sTangent).*length4s).^2;
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
