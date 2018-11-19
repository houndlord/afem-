function eta4e = estimateP1AveragingP1(f,g,u4Db,x,c4n,n4e,n4sDb,n4sNb)
%% estimateP1AveragingP1 - averaging error estimator for P1 element
% Estimate the gradient error of the P1 finite element solution by
% comparing the P0-gradient of the discrete solution to the P1 average
% of itself.
%
% Input:     f	       right-hand side of the problem definition
%            g         Neumann boundary condition
%            u4Db      Dirichlet boundary condition
%            x         P1 basis coefficients of u given by solve
%            c4n       coordinates for the nodes of the mesh
%            n4e       nodes for the elements of the mesh
%            n4sDb	   the nodes of the sides in the Dirichlet boundary
%            n4sNb	   the nodes of the sides in the Neumann boundary
%
% Output:    eta4e     error for each element

    %% Initialisation
    nrElems = size(n4e,1);
    grad4e  = zeros(nrElems,2);

    %% Compute the P0-gradient of u
    for elem = 1 : nrElems
        grads = [1,1,1;c4n(n4e(elem,:),:)'] \ [0,0;eye(2)];
        grad4e(elem,:) = x(n4e(elem,:))'*grads;
    end
    
    eta4e = estimateSigmaAveragingP1(c4n,n4e,grad4e);
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
