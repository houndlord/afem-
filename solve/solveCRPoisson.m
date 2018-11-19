function [x,nrDof,A,b] = solveCRPoisson(f,g,u4Db,c4n,n4e,n4sDb,n4sNb)
%% solveCR - solve the Possion problem using the CR element.
% Solves the Poisson problem for given right-hand side f and 
% Neumann boundary condition g on the domain given by c4n and n4e.
%
% Input:     f	       right-hand side of the problem definition
%            g         Neumann boundary condition
%            u4Db      Dirichlet boundary condition
%            c4n       coordinates for the nodes of the mesh
%            n4e       nodes for the elements of the mesh
%            n4sDb	   the nodes of the sides in the Dirichlet boundary
%            n4sNb	   the nodes of the sides in the Neumann boundary
%
% Output:    x         basis coefficients of the numerical solution w.r.t.
%                      the Crouzeix-Raviart basis.
%            nrDof     number of degrees of freedom
%            A         the matrix A of the linear system created
%            b         the right side of the linear system created

    %% Initialisation
    nrElems = size(n4e,1);
    s4e = computeS4e(n4e);
    nrSides = max(max(s4e));
    area4e = computeArea4e(c4n,n4e);
    mid4e = computeMid4e(c4n,n4e);
    s4n = computeS4n(n4e);
    % Dirichlet boundary sides
    DbSides = zeros(1,size(n4sDb,1));
    for i = 1:size(n4sDb,1)
        DbSides(i) = s4n(n4sDb(i,1),n4sDb(i,2));
    end
    % Neumann boundary sides
    NbSides = zeros(1,size(n4sNb,1));
    for i = 1:size(n4sNb,1)
        NbSides(i) = s4n(n4sNb(i,1),n4sNb(i,2));
    end
    % degrees of freedom: one per non-Dirichlet side
    dof = setdiff(1:nrSides,DbSides);
    nrDof = length(dof);
    
    Alocal = zeros(3,3,nrElems);
    b = zeros(nrSides,1);

    %% Create the stiffness matrix A and right-hand side b
    for elem = 1 : nrElems
        nodes = n4e(elem,:);   % nodes of this element
        sides = s4e(elem,:);   % sides of this element
        coords = c4n(nodes,:); % coordinates for the nodes
        area = area4e(elem);   % area of this element
        grads = [coords';1 1 1]\[-2 0; 0 -2; 0 0]; % gradients for CR basis
        grads = grads([3 1 2],:); % reorder to fit DoF numbering
        Alocal(:,:,elem) = area * grads * grads'; % local stiffness matrix
        mid = mid4e(elem,:);     % midpoint of this element
        b(sides) = b(sides) + area*f(mid)*ones(3,1)/3; % right-hand side
    end

    % assembly of the global stiffness matrix A
    s4eT = s4e';
    I = [s4eT;s4eT;s4eT];
    J = [s4eT(:),s4eT(:),s4eT(:)]';
    A = sparse(I(:),J(:),Alocal(:));

    %% Neumann boundary conditions
    length4NbSides = computeLength4s(c4n,n4sNb);
    mid4NbSides = computeMid4s(c4n,n4sNb);
    b(NbSides) = b(NbSides) + length4NbSides .* g(mid4NbSides);
    %% Dirichlet boundary conditions
    x = zeros(nrSides,1);
    x(DbSides) = u4Db(computeMid4s(c4n, n4sDb));
    b = b - A * x;
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
