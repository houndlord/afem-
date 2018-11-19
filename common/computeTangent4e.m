function tangent4e = computeTangent4e(c4n,n4e)
%% computeTangent4e - Tangents for elements.
%   computeTangent4e(c4n, n4e) computes the tangent vectors of all sides of all
%                          elements in the decomposition. The tangent
%                          vectors are normed and point in counterclockwise
%                          direction. c4n and n4e are as specified in the
%                          documentation.
%
%   See also: computeTangent4s, computeNormal4e

    if isempty(n4e)
        tangent4e = zeros(0,2);
        return;
    end

    %% Compute tangent4e.
    allSides = [n4e(:,[1 2]); n4e(:,[2 3]); n4e(:,[3 1])];
    c4start  = c4n(allSides(:,1),:);
    c4end    = c4n(allSides(:,2),:);
    lengths  = sqrt(sum((c4end-c4start).^2,2));
    tangents = (c4end - c4start)./[lengths lengths];
    tangent4e(1,:,:) = tangents(1:size(n4e,1),:)';
    tangent4e(2,:,:) = tangents(size(n4e,1)+1:2*size(n4e,1),:)';
    tangent4e(3,:,:) = tangents(2*size(n4e,1)+1:3*size(n4e,1),:)';
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
