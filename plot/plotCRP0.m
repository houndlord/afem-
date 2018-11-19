function plotCRP0(c4n, n4e, u, p)
%% plotCRP0 - plot a discrete solution for CR-P0 finite elements
%  Input: c4n, n4e - mesh
%         u        - basis coefficients for (CR)^2 basis functions
%         p        - pressure on each element (optional)

    %% Plot p.
    if nargin > 3
        hold on;
        CData = zeros(max(n4e(:)), 3);
        CData(n4e(:, 1), :) = 0.5 + 0.5*(max(p)-p(:))/(max(p)-min(p))*[1 1 1];
        patch('Vertices', c4n, 'Faces', n4e, ...
              'FaceVertexCData', CData, 'FaceColor', 'flat');
    end

    %% Plot u.
    mid4s = computeMid4s(c4n,computeN4s(n4e));
    quiver2(mid4s(:, 1), mid4s(:, 2), u(:, 1), u(:, 2), ...
            'n=', 0.1, 'w=', [1 1]);
    axis equal tight;
    title({'Discrete flux'; [num2str(numel(u) + size(n4e, 1)), ' dofs']});
    drawnow;

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
