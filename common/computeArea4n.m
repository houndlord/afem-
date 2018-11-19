function area4n = computeArea4n(c4n,n4e)
%% computeArea4n - Area for node patches.
%   A node patch is the union of all elements a node is part of.
%   computeArea4n(c4n, n4e) computes the area of the node patch for each node
%                       of a decomposition where c4n, n4e are as specified
%                       in the documentation.
%
%   See also: computeArea4e

    if isempty(n4e)
        area4n = zeros(0,1);
        return;
    end

    %% Compute area4e.
    % Get the x- and y-coordinates for each node of each element and
    % compute the area of all elements simulateously.
    x1 = c4n(n4e(:,1),1);
    x2 = c4n(n4e(:,2),1);
    x3 = c4n(n4e(:,3),1);
    y1 = c4n(n4e(:,1),2);
    y2 = c4n(n4e(:,2),2);
    y3 = c4n(n4e(:,3),2);
    
    area4e = 1/2 * (x1.*(y2 - y3) + x2.*(y3 - y1) + x3.*(y1 - y2));
    
    %% Compute area4n.
    % Use area4n and n4e to accumulate area4n.
    area4e = area4e * ones(1,3);
    nrNodes = size(c4n,1);
    area4n = accumarray(n4e(:),area4e(:),[nrNodes 1]);
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
