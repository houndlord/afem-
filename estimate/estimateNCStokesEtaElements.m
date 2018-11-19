function [eta4e,n4s] = estimateNCStokesEtaElements(c4n,n4e,n4sDb,f,Du4Db1,Du4Db2,u,gradU4e)
%% estimateCREtaSides - error estimator for CR solution of Stokes problem 
% Estimate the energy error of the CR finite element solution of the Stokes
% problem by the jumps of the discrete solution's tangents along the sides.
%
% Input:    c4n     coordinates for the nodes of the mesh
%           n4e     nodes for the elements of the mesh
%           n4sDb   the nodes of the sides in the Dirichlet boundary
%           f       right-hand side of the problem definition
%           Du4Db1  first row of the gradient of the Dirichlet boundary 
%                   condition for computation of the normal derivative
%           Du4Db2  second row of the gradient of the Dirichlet boundary 
%                   condition for computation of the normal derivative
%           u       basis coefficients of the numerical solution of the 
%                   velocity field w.r.t. the Crouzeix-Raviart basis
%           gradU4e piecewise gradient of the discrete solution u
%
% Output:   eta4e   error for each element
%           n4s     vector which contains in one row for each side the 
%                   numbers of the corresponding nodes

    %% Initialization
    s4e=computeS4e(n4e);
    n4s=computeN4s(n4e);
    s4n=computeS4n(n4e);
    tangent4e=computeTangent4e(c4n,n4e);
    tangent4s=computeTangent4s(c4n,n4s);
    nrElems=size(n4e,1);
    nrSides=size(n4s,1);
    length4s=computeLength4s(c4n,n4s);
    area4e=computeArea4e(c4n,n4e);
    
    %% Compute gradient.
    if nargin < 8
        gradU4e=zeros(2,2,nrElems);
        for j=1:nrElems
            grads=[c4n(n4e(j,:),:)'; 1 1 1]\[-2 0; 0 -2; 0 0];
            gradU4e(1,:,j)=u(s4e(j,:),1)' * grads([3 1 2],:);
            gradU4e(2,:,j)=u(s4e(j,:),2)' * grads([3 1 2],:);
        end
    end
    
    %% Compute inner jumps
    jump4s=zeros(nrSides,2);    
    for j=1:nrElems
	    jump4s(s4e(j,:),:) = jump4s(s4e(j,:),:) ...
                                    + [tangent4e(:,:,j)*gradU4e(1,:,j)' ...
                                       tangent4e(:,:,j)*gradU4e(2,:,j)'];
    end
    
    %% Dirichlet jumps
    l4DbS=computeLength4s(c4n,n4sDb);
    mean4DbSides1=integrate(c4n,n4sDb,@(x,y,z)(Du4Db1(y)),10)...
                                                        ./[l4DbS l4DbS];
    mean4DbSides2=integrate(c4n,n4sDb,@(x,y,z)(Du4Db2(y)),10)...
                                                        ./[l4DbS l4DbS];
    for j=1:size(n4sDb,1)
        nodes=n4sDb(j,:)';
        side=s4n(nodes(1),nodes(2));
        jump4s(side,:) = jump4s(side,:) ...
                            - [mean4DbSides1(j,:)*tangent4s(side,:)'...
                               mean4DbSides2(j,:)*tangent4s(side,:)'];
    end
    
    %% Compute eta4e
    jump4s_L2normSq=(jump4s(:,1).^2 + jump4s(:,2).^2) .* length4s;
    [~,L2normSq4e]=L2Norm(c4n,n4e,@(x,y,z)f(y),6); 
    eta4e = area4e.*sum(L2normSq4e,2) +...
                                sqrt(area4e).*sum(jump4s_L2normSq(s4e),2);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2009 
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
