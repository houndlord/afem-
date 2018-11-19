function plotRT04e(c4n, n4e, p, OPTtitle4p)
%% plotRT04e plots piecewise constant approximation of Raviart-Thomas flux p
% INPUT:  c4n, n4e      - the triangulation
%         p             - coefficient of the flux with respect to broken
%                         RT basis functions
%         OPTtitle4p    - optional parameter, which defines the title
%                         for the flux p. The default value is
%                         'RT flux plot'.

    if nargin == 4
        title4p = OPTtitle4p;
    else
        title4p = 'RT flux plot';
    end

    %% Compute components of flux
    mid4e = computeMid4e(c4n, n4e);
    u4e = p(:, 1) + p(:, 3).*mid4e(:, 1);
    v4e = p(:, 2) + p(:, 3).*mid4e(:, 2);

    %% Plot
    quiver2(mid4e(:, 1), mid4e(:, 2), u4e, v4e,...
            'n=', 0.1, 'w=', [1 1]);
    axis equal tight;
    title(title4p);
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
