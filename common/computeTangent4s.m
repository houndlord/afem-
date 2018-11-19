function tangent4s = computeTangent4s(c4n,n4s)
%% computeTangent4s - Tangents for sides.
%   computeTangent4s(c4n, n4s) computes the tangent vector of each side of the
%                          decomposition. c4n and n4s are as specified in
%                          the documentation.
%
%   See also: computeN4s, computeTangent4e, computeNormal4s
    
    if isempty(n4s)
        tangent4s = zeros(0,2);
        return;
    end

    %% Compute length4s.
    c4start = c4n(n4s(:,1),:);
    c4end = c4n(n4s(:,2),:);
    length4s = sqrt(sum((c4end-c4start).^2,2));
    
    %% Compute tangent4s.              
    tangent4s = (c4end-c4start)./[length4s length4s];
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
