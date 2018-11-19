function [av4e] = estimateL2Averaging(c4n,n4e,p)
%% error4eAveragingL2 - 
% Compute the averaging L2 error for a given finite element solution p of
% stress or flux.
%
% Input:    c4n     coordinates for the nodes of the mesh
%           n4e     nodes for the elements of the mesh
%           p       discrete solution of stres or flux
%
% Output:   av4e    averaging L2 error on each element
%

    %% Initialization
    area4e=computeArea4e(c4n,n4e); 
    mid4e=computeMid4e(c4n,n4e); 
    p=sparse(p);
    
    
    %% Computation of averaging
    stack3=kron(ones(3,1),1:size(n4e,1));
    n4e2Area=sparse(stack3,n4e',[area4e; area4e; area4e]');
    Ap=(n4e2Area'*p(:,[1 2]))./(sum(n4e2Area,1)'*sparse([1 1]));

    % calculate nodal values
    tmp = reshape(p(stack3, [1 2]) - Ap(n4e',:) ...
    + p(stack3,[3 3]).*(c4n(n4e',:)-mid4e(stack3,:)),3,2*size(n4e,1));
    
    % calculate L2-integral
    av4e = sum(sum(repmat(area4e',3,2).*tmp.*tmp,1)...
    + repmat(area4e',1,2).*sum(tmp,1).*sum(tmp,1))/12;
    
    
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
