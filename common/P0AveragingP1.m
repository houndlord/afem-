function val = P0AveragingP1(c4n,n4e,sigma4e)
%% Averaging of a given P0 function to a P1 function
%  c4n, n4e - a triangular mesh
%  sigma4e  - values of a pieceweise P0 function v on the mesh

    %% Initialisation
    d = size(sigma4e,2);
    area4n = computeArea4n(c4n,n4e);
    area4e = computeArea4e(c4n,n4e);
    
    %% Compute node values.
    weightedV = sigma4e.*(area4e*ones(1,d))*[eye(d),eye(d),eye(d)];
    I = n4e(:,[ones(1,d),2*ones(1,d),3*ones(1,d)]);          % node indeces
    J = (ones(size(n4e,1),1)*(1:d))*[eye(d),eye(d),eye(d)]; % component indeces
    val = accumarray([I(:),J(:)],weightedV(:))./(area4n*ones(1,d));
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
