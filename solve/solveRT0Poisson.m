function [p,u,nrDof] = solveRT0Poisson(f,g,u4Db,c4n,n4e,n4sDb,n4sNb)
% solveRT0Poisson.m
% solves the Poisson problem for given right-hand side f and 
% Neumann boundary condition g on the domain given by c4n and n4e.
% The realised implementation is based on the Lagrange Multiplier technique
%
%
% Input:     f	       right side of the problem definition
%            g         Neumann boundary condition
%            u4Db      Dirichlet boundary condition
%            c4n       coordinates for the nodes of the mesh
%            n4e       nodes for the elements of the mesh
%            n4sDb	   the nodes of the sides in the Dirichlet boundary
%            n4sNb	   the nodes of the sides in the Neumann boundary
%
% Output:    x         solution x = (u,p) from RT0 and P0; 
%                      u includes an vector for each element, 
%                      while p includes an value for each element
%            nrDof     number of degrees of freedom
%            A         the matrix A of the linear system created (optional)
%            b         the right side of the linear system created
%            (optional)

    %% Initialisation
    n4s = computeN4s(n4e); % Nodes for sides
    s4n = computeS4n(n4e); % Sides for nodes
    s4e = computeS4e(n4e); % Sides for elements
    e4s = computeE4s(n4e); % Elements for sides
    area4e = computeArea4e(c4n,n4e); % Area for elements
    mid4e = computeMid4e(c4n,n4e); % Midpoints of elements
    mid4s = computeMid4s(c4n,n4s); % Midpoints of sides
    length4s= computeLength4s(c4n,n4s); % Length of Sides 
    normal4s = computeNormal4s(c4n,n4s); % Normals for sides
    
    nrElems  = size(n4e,1); % Number of elements
    nrISides = length(find(e4s(:,2)~=0)); % Number of inner Sides
    nrNbSides = size(n4sNb,1); % Number of Neumann  boundary sides
    nrDbSides = size(n4sDb,1); % Number of Dirichlet boundary sides

    n4iSides = n4s(e4s(:,2)~=0,:); % Nodes for interiour sides

    % Blocks of the global stiffnes Matrices
    B = sparse( 3*nrElems, 3*nrElems);
    C = sparse( 3*nrElems,   nrElems);
    D = sparse( 3*nrElems,   nrISides);
    F = sparse( 3*nrElems,   nrNbSides);

    b = zeros(4*nrElems+nrISides+nrNbSides,1); % right-hand side (RHS)

    % Assembling the blocks B, C and the RHS
    for curElem = 1 : nrElems
        % Summing up the length (norm) of all sides of curElem
        s = sum(length4s(s4e(curElem,:)).^2);
        % assembling the blocks of global stiffness matrices B and C
        B( 3*curElem-[2,1,0],3*curElem-[2,1,0] ) = area4e(curElem) * diag([1,1,s/36]); 
        C( 3*curElem-[2,1,0],curElem ) = [0;0;2*area4e(curElem)];
        % assembling the RHS b
        b(3*nrElems+curElem) = -area4e(curElem) * f(mid4e(curElem,:));
    end

    %% Assembling the block D (Condition for interior sides)
    for curISide = 1 : nrISides
        side = s4n(n4iSides(curISide,1),n4iSides(curISide,2));
        curNormal = normal4s(side,:)';
        h1 = (c4n(n4iSides(curISide,1),:)-mid4e(e4s(side,1),:)) * curNormal;
        h2 = (c4n(n4iSides(curISide,1),:)-mid4e(e4s(side,2),:)) * curNormal;
        % assembling the blocks of global stiffness matrix D
        D( [3*e4s(side,1)-[2,1,0],3*e4s(side,2)-[2,1,0]],curISide ) ...
           = -length4s(side)*[curNormal;h1;-curNormal;-h2];
    end

    %% Dirichlet conditions; Assembly of b
    for curDbSide = 1 : nrDbSides
        side = s4n(n4sDb(curDbSide,1),n4sDb(curDbSide,2));
        curNormal = normal4s(side,:)';
        curElement = e4s(side,1);
        h = (c4n(n4sDb(curDbSide,1),:)-mid4e(curElement,:)) * curNormal;
        b(3*curElement-[2,1,0]) = b(3*curElement-[2,1,0]) ...
            + u4Db(mid4s(side,:)) *length4s(side)*[curNormal;h];
    end

    %% Neumann conditions; Assembling block F
    for curNbSide = 1 : nrNbSides
        side = s4n(n4sNb(curNbSide,1),n4sNb(curNbSide,2));
        curNormal = normal4s(side,:)';
        curElem = e4s(side,1);
        h = (c4n(n4sNb(curNbSide,1),:)-mid4e(curElem,:)) * curNormal;
        F(3*curElem-[2,1,0],curNbSide) = length4s(side)*[curNormal;h]; 
        b(4*nrElems+nrISides+curNbSide) = length4s(side)*g(mid4s(side,:));
    end

    %% Assembling the global stiffness matrix, including zero block matrices
    O1 = sparse(size(C,2),size(C,2)+size(D,2)+size(F,2));
    O2 = sparse(size(D,2),size(C,2)+size(D,2)+size(F,2));
    O3 = sparse(size(F,2),size(C,2)+size(D,2)+size(F,2));
    A = [B ,  C, D, F; ...
         C', O1      ; ...
         D', O2      ; ...
         F', O3      ];

    nrDof = nrISides + nrNbSides;
    
    %% Solve the linear System of equations
    x = A \ b;
    
    p = reshape(x(1:3*nrElems),3,nrElems)';
    u = x(3*nrElems+1:4*nrElems);
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
