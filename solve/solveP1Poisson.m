function [x,nrDof,A,b] = solveP1Poisson(f,g,u4Db,c4n,n4e,n4sDb,n4sNb)
% solveP1Poisson.m
% solves the Poisson problem for given righthand side f and 
% Neumann boundary condition g on the domain given by c4n and n4e.
%
% Input:     f	       right side of the problem definition
%            g         Neumann boundary condition
%            u4Db      Dirichlet boundary condition
%            c4n       coordinates for the nodes of the mesh
%            n4e       nodes for the elements of the mesh
%            n4sDb	   the nodes of the sides in the Dirichlet boundary
%            n4sNb	   the nodes of the sides in the Neumann boundary
%
% Output:    x         solution u for each node corresponding to c4n
%            nrDof     number of degrees of freedom
%            A         the matrix A of the linear system created (optional)
%            b         the right side of the linear system created (optional)

    %% Initialisation
    nrNodes = size(c4n,1);            % number of nodes
    nrElems = size(n4e,1);            % number of elements
    DbNodes = unique(n4sDb);           % Dirichlet boundary nodes
    dof = setdiff(1:nrNodes,DbNodes); % free nodes to be approximated
    nrDof = length(dof);
    Alocal = zeros(3,3,nrElems);      % local stiffness matrices
    b = zeros(nrNodes,1);             % vector for right-hand side

    %% Create the stiffness matrix A and right-hand side b
    area4e = computeArea4e(c4n,n4e);
    mid4e = computeMid4e(c4n,n4e);
    % calculate for each element the gradients of the three nodal
    % basisfunctions, the right-hand side and the local stiffness matrix A 
    for elem = 1 : nrElems
        nodes = n4e(elem,:);   % nodes of the triangle
        coords = c4n(nodes,:); % coordinates of the three nodes
        area = area4e(elem);   % area of the current element
        mid = mid4e(elem,:);   % midpoint of the triangle
        grads = [1,1,1;coords'] \ [0,0;eye(2)];
        b(nodes) = b(nodes)+(1/3) * area * f(mid)*[1;1;1];
        Alocal(:,:,elem) = area * grads * grads';
    end

    % assembly of the global stiffness matrix A
    n4eT = n4e';
    I = [n4eT;n4eT;n4eT];
    J = [n4eT(:),n4eT(:),n4eT(:)]';
    A = sparse(I(:),J(:),Alocal(:));

    %% Neumann boundary conditions
    % involving Neumann boundary for all Neumann edges    
    nrNbSides = size(n4sNb,1);
    length4NbSides = computeLength4s(c4n,n4sNb);
    mid4NbSides = computeMid4s(c4n,n4sNb);
    for NbSide = 1 : nrNbSides
        nodes = n4sNb(NbSide,:);         % nodes of the edge
        len = length4NbSides(NbSide);   % length of the edge
        mid = mid4NbSides(NbSide,:);      % midpoint of the edge
        b(nodes) = b(nodes) + (1/2) * len * g(mid)*[1;1];
    end

    %% Dirichlet boundary conditions
    x = zeros(nrNodes,1); % get the Dirichlet values
    DbCoords = c4n(DbNodes,:); % coordinates of Dirichlet nodes
    x(DbNodes) = u4Db(DbCoords);
    b = b - A * x;  % substract inhomogenous boundary

    %% solve the algebraic equation
    x(dof) = A(dof,dof) \ b(dof);

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
