function eta4e = estimateP1EtaElements(f,g,u4Db,x,c4n,n4e,n4sDb,n4sNb)
% estimateP1EtaSides.m
% Estimate the energy error of the P1 finite element solution by the
% jumps of the discrete solution over the sides and a volume term.
%
% Input:     f	       right-hand side of the problem definition
%            g         Neumann boundary condition
%            u4Db      Dirichlet boundary condition
%            u         solution u given by solve
%            c4n       coordinates for the nodes of the mesh
%            n4e       nodes for the elements of the mesh
%            n4sDb	   the nodes of the sides in the Dirichlet boundary
%            n4sNb	   the nodes of the sides in the Neumann boundary
%
% Output:    eta4e     error for each element

    %% initialisation
    s4e      = computeS4e(n4e);
    area4e   = computeArea4e(c4n,n4e);
    mid4e    = computeMid4e(c4n,n4e);
    
    %% Compute eta4s.
    eta4s = estimateP1EtaSides(f,g,u4Db,x,c4n,n4e,n4sDb,n4sNb);

    %% Compute eta4e. 
    eta4e = (area4e.*f(mid4e)).^2 + sum(eta4s(s4e),2);
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
